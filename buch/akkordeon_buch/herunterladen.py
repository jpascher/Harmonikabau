#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
herunterladen.py
================
Lädt alle Artikel und Bilder des Wikipedia-Buches
"Das Akkordeon und seine Geschichte" in ein lokales Verzeichnis
und packt das Ganze in ein ZIP, das im Chat hochgeladen werden kann.

Voraussetzungen: NUR Python 3.8+ und Internet — sonst nichts.
                 (kein pandoc, kein LaTeX nötig für diesen Schritt)

Aufruf:
  python3 herunterladen.py
  python3 herunterladen.py --start-bei "Bandoneon"   # Wiederaufnahme
  python3 herunterladen.py --keine-bilder            # nur Texte holen

Ausgabe:
  akkordeon_quellen.zip       <-- diese Datei im Chat hochladen
  akkordeon_quellen/          (das entpackte Arbeitsverzeichnis)
    artikel/                  *.wiki  — Wikitext jedes Artikels
    bilder/                   *.jpg/png — heruntergeladene Bilder
    bilder_index.json         Zuordnung Artikel → Bilddateien
    info.json                 Metadaten zum Lauf
"""

import argparse
import json
import re
import shutil
import sys
import time
import urllib.parse
import urllib.request
import zipfile
from pathlib import Path
from typing import Dict, List, Optional


# ---------------------------------------------------------------------------
# Buchstruktur — exakt nach der Wikipedia-Buchseite
# ---------------------------------------------------------------------------

KAPITEL = [
    # ===== BAND 1: Das Akkordeon — Aufbau, Bauteile, Hersteller =====
    ("Balginstrumente", [
        "Harmonikainstrument",
        "Handzuginstrument",
        "Akkordeon",
    ]),
    ("Funktion und Geschichte der durchschlagenden Stimmzunge", [
        "Zunge (Tonerzeuger)",
        "Durchschlagende Zunge",
    ]),
    ("Bauteile und Mechanik", [
        "Kanzelle",
        "Cassotto",
        "Register (Akkordeon)",
        "Tremolo (Akkordeon)",
        "Stradella-Bass",
        "Konverterbass",
        "Melodiebass",
        "Helikonstimmplatten",
        "Gleichton",
        "Wechseltönig",
        "Wechselbass",
        "Stimmung (Musik)",
    ]),
    ("Theorie und Berufsbezogenes", [
        "Diatonik",
        "Vallotti-Stimmung",
        "Handzuginstrumentenmacher",
        "Akkordeonschule",
    ]),
    ("Erfinder und Pioniere des Akkordeons", [
        "Christian Gottlieb Kratzenstein",
        "Bernhard Eschenbach",
        "Johann Caspar Schlimbach",
        "Carl Friedrich Voit",
        "Friedrich Sturm (Instrumentenbauer)",
        "Anton Haeckl",
        "Anton Reinlein",
        "Christian Friedrich Ludwig Buschmann",
        "Cyrill Demian",
        "Charles Wheatstone",
        "William M. Goodrich",
        "Adolf Müller senior",
        "Johann Wilhelm Rudolph Glier",
        "Paolo Soprani",
    ]),
    ("Bedeutende Harmonikabauer (Helikonbässe)", [
        "Josef Hlaváček",
        "Lubas & Sohn",
        "Josef Fleiß",
        "Anton Mervar",
        "Peter Stachl",
        "Novak Harmonikas",
        "Otto Ludwig (Harmonikabauer)",
        "Roman Gombotz",
        "Othmar Kühn",
    ]),
    ("Produktionszentren", [
        "Geschichte des Akkordeonbaus in Klingenthal",
        "Castelfidardo",
    ]),
    # ===== BAND 2: Akkordeon-Arten und verwandte Instrumente =====
    ("Akkordeon-Vorläufer", [
        "Maultrommel",
        "Sheng (Instrument)",
        "Aeoline (Musikinstrument)",
        "Physharmonika",
        "Symphonium",
        "Schrammelharmonika",
        "Konzertina",
        "Bandoneon",
    ]),
    ("Diatonische Akkordeons", [
        "Diatonisches Akkordeon",
        "Wiener Modell",
        "Organetto",
        "Französisches Akkordeon",
        "Schottisches Akkordeon",
        "Irisches Akkordeon",
        "Heligonka",
        "Steirische Harmonika",
        "Trikitixa",
        "Schwyzerörgeli",
    ]),
    ("Chromatische Akkordeons", [
        "Chromatisches Knopfakkordeon",
        "Garmon",
        "Bajan",
        "Pianoakkordeon",
        "Bassakkordeon",
        "Harmonium",
    ]),
    ("Mechanische Automaten mit durchschlagenden Zungen", [
        "Orchestrion",
        "Panharmonikon",
        "Apollonicon",
        "Belloneon",
        "Mechanischer Trompeter",
        "Schachtürke",
    ]),
    ("Förderer der Durchschlagzungen", [
        "Georg Joseph Vogler",
        "Wolfgang von Kempelen",
        "Johann Nepomuk Mälzel",
        "Leonhard Mälzel",
        "Friedrich Kaufmann (Instrumentenbauer)",
        "Ignaz Kober",
        "Johannes Weinrich (Volkskünstler)",
        "Franz Leppich",
    ]),
    ("Andere Instrumente mit Durchschlagzunge", [
        "Harmonichord",
        "Terpodion",
        "Aerophon",
        "Mundharmonika",
        "Glasharmonika",
        "Äolsharfe",
        "Anemochord",
    ]),
]

# Pfade
ROOT = Path(__file__).resolve().parent
WORK = ROOT / "akkordeon_quellen"
ART_DIR = WORK / "artikel"
BILD_DIR = WORK / "bilder"
INDEX_PFAD = WORK / "bilder_index.json"
INFO_PFAD = WORK / "info.json"
ZIP_PFAD = ROOT / "akkordeon_quellen.zip"

WIKI_API = "https://de.wikipedia.org/w/api.php"
USER_AGENT = ("AkkordeonBuchDownload/1.1 "
              "(https://github.com/jpascher; johann.pascher@gmail.com) "
              "Python-urllib")

# Welche Bildformate herunterladen — SVG, Audio, Video überspringen
GUTE_ENDUNGEN = {".jpg", ".jpeg", ".png", ".gif", ".bmp", ".tif", ".tiff"}

# Maximale Größe pro Bild (verhindert versehentliche Riesen-Downloads)
MAX_BILD_BYTES = 8 * 1024 * 1024   # 8 MB


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def slugify(s: str) -> str:
    """Macht aus einem Titel einen sicheren Dateinamen."""
    s = s.lower()
    s = (s.replace("ä", "ae").replace("ö", "oe").replace("ü", "ue")
          .replace("ß", "ss"))
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s


def http_get(url, retries=4, pause=1.0, max_bytes=None):
    """Robuster HTTP-GET mit User-Agent.
    
    Bei HTTP 429 (rate limited) warten wir deutlich länger als beim normalen
    Wiederholen — Wikimedia weist eingehende Anfragen sonst dauerhaft ab.
    """
    last_err = None
    for versuch in range(retries):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with urllib.request.urlopen(req, timeout=30) as resp:
                if max_bytes is not None:
                    return resp.read(max_bytes + 1)
                return resp.read()
        except urllib.error.HTTPError as e:
            last_err = e
            if e.code == 429:
                # Lange Auszeit gegen Rate-Limit (5, 15, 30, 60 Sekunden)
                wartezeit = [5, 15, 30, 60][min(versuch, 3)]
                print(f"        Rate-Limit (429) — warte {wartezeit}s ...",
                      file=sys.stderr, flush=True)
                time.sleep(wartezeit)
                continue
            # Andere HTTP-Fehler: kürzeres Zurücksetzen
            time.sleep(pause * (versuch + 1))
        except Exception as e:
            last_err = e
            time.sleep(pause * (versuch + 1))
    raise RuntimeError(f"Abruf fehlgeschlagen ({url}): {last_err}")


def fetch_wikitext(titel: str) -> Optional[str]:
    """Holt den Wiki-Quelltext (Wikitext) eines deutschen Wikipedia-Artikels."""
    params = {
        "action": "query",
        "prop": "revisions",
        "rvprop": "content",
        "rvslots": "main",
        "format": "json",
        "formatversion": "2",
        "titles": titel,
        "redirects": "1",
    }
    url = WIKI_API + "?" + urllib.parse.urlencode(params)
    try:
        data = json.loads(http_get(url).decode("utf-8"))
    except Exception as e:
        print(f"      FEHLER (Wikitext): {e}", file=sys.stderr)
        return None
    pages = data.get("query", {}).get("pages", [])
    if not pages or "missing" in pages[0]:
        return None
    revs = pages[0].get("revisions") or []
    if not revs:
        return None
    return revs[0]["slots"]["main"]["content"]


def fetch_bild_liste(titel: str) -> List[str]:
    """Holt die Liste der eingebundenen Bilddateien (Datei:Foo.jpg) für einen Artikel."""
    params = {
        "action": "query",
        "prop": "images",
        "imlimit": "100",
        "format": "json",
        "formatversion": "2",
        "titles": titel,
        "redirects": "1",
    }
    url = WIKI_API + "?" + urllib.parse.urlencode(params)
    try:
        data = json.loads(http_get(url).decode("utf-8"))
    except Exception:
        return []
    pages = data.get("query", {}).get("pages", [])
    if not pages:
        return []
    namen = [img["title"] for img in pages[0].get("images", [])]
    # Aussortieren: Endung gehört zu unseren akzeptierten
    raus: List[str] = []
    for n in namen:
        endung = Path(n).suffix.lower()
        if endung in GUTE_ENDUNGEN:
            raus.append(n)
    return raus


def fetch_bild_url(dateititel: str) -> Optional[str]:
    """Bestimmt die URL einer Wiki-Datei (Datei:Foo.jpg).
    
    Wir bitten Wikimedia gezielt um eine THUMBNAIL-Version mit max. 1200 px
    Breite — das wird in den HTTP-429-Meldungen explizit als richtiger
    Weg für Skript-Downloads genannt. Originalauflösungen sind viel zu
    groß und werden bei automatisierten Zugriffen abgelehnt.
    """
    params = {
        "action": "query",
        "prop": "imageinfo",
        "iiprop": "url|size",
        "iiurlwidth": "1200",       # <-- Thumbnail dieser Breite anfordern
        "format": "json",
        "formatversion": "2",
        "titles": dateititel,
    }
    url = WIKI_API + "?" + urllib.parse.urlencode(params)
    try:
        data = json.loads(http_get(url).decode("utf-8"))
    except Exception:
        return None
    pages = data.get("query", {}).get("pages", [])
    if not pages:
        return None
    info = pages[0].get("imageinfo") or []
    if not info:
        return None
    # Wenn Wikimedia eine thumburl liefert (Bild ist breiter als 1200 px),
    # nehmen wir die. Sonst die normale URL (Bild ist von Haus aus klein).
    return info[0].get("thumburl") or info[0].get("url")


def lade_bild(url: str, ziel: Path) -> int:
    """Lädt ein Bild herunter; gibt die Größe in Bytes zurück, 0 bei Misserfolg."""
    if ziel.exists() and ziel.stat().st_size > 0:
        return ziel.stat().st_size
    try:
        daten = http_get(url, max_bytes=MAX_BILD_BYTES)
        if len(daten) > MAX_BILD_BYTES:
            return 0  # zu groß — überspringen
        ziel.parent.mkdir(parents=True, exist_ok=True)
        ziel.write_bytes(daten)
        return len(daten)
    except Exception as e:
        print(f"        Bild-Download fehlgeschlagen: {e}", file=sys.stderr)
        return 0


def sicherer_dateiname(name: str) -> str:
    """Aus 'Datei:Foo Bar.jpg' wird 'foo_bar.jpg'."""
    base = re.sub(r"^Datei:|^File:|^Bild:", "", name).strip().replace(" ", "_")
    return slugify(Path(base).stem) + Path(base).suffix.lower()


# ---------------------------------------------------------------------------
# Hauptablauf
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--start-bei", type=str, default=None,
                        help="Erst ab diesem Artikeltitel weitermachen (Wiederaufnahme).")
    parser.add_argument("--keine-bilder", action="store_true",
                        help="Bilder NICHT herunterladen, nur Wikitext.")
    parser.add_argument("--kein-zip", action="store_true",
                        help="Am Ende KEIN ZIP packen (für eigene Weiterverarbeitung).")
    args = parser.parse_args()

    ART_DIR.mkdir(parents=True, exist_ok=True)
    BILD_DIR.mkdir(parents=True, exist_ok=True)

    print("=" * 64)
    print(" Wikipedia-Buch: Das Akkordeon und seine Geschichte")
    print(" Quellen-Download (Wikitext + Bilder)")
    print("=" * 64)

    # Bestehenden Bilder-Index laden, falls vorhanden (Wiederaufnahme)
    bilder_index: Dict[str, List[str]] = {}
    if INDEX_PFAD.exists():
        try:
            bilder_index = json.loads(INDEX_PFAD.read_text(encoding="utf-8"))
        except Exception:
            bilder_index = {}

    skip = args.start_bei is not None
    fehlend: List[str] = []
    erfolgreich: List[str] = []
    bilder_gesamt = 0

    for teil_titel, artikel in KAPITEL:
        print(f"\n=== {teil_titel} ===")
        for art in artikel:
            if skip:
                if art == args.start_bei:
                    skip = False
                else:
                    continue

            print(f"  • {art}", flush=True)

            # 1) Wikitext holen
            wikitext = fetch_wikitext(art)
            if wikitext is None:
                print("      NICHT GEFUNDEN — übersprungen.", file=sys.stderr)
                fehlend.append(art)
                time.sleep(0.4)
                continue

            slug = slugify(art)
            wiki_pfad = ART_DIR / f"{slug}.wiki"
            wiki_pfad.write_text(wikitext, encoding="utf-8")
            erfolgreich.append(art)

            # 2) Bilder holen
            if not args.keine_bilder:
                bilder_namen = fetch_bild_liste(art)
                lokale_dateien: List[str] = []
                geladen = 0
                for fn in bilder_namen:
                    sicher = sicherer_dateiname(fn)
                    ziel = BILD_DIR / sicher
                    if ziel.exists() and ziel.stat().st_size > 0:
                        lokale_dateien.append(sicher)
                        continue
                    url = fetch_bild_url(fn)
                    if url is None:
                        continue
                    if lade_bild(url, ziel) > 0:
                        lokale_dateien.append(sicher)
                        geladen += 1
                    time.sleep(0.4)  # höflich gegenüber Wikimedia (Bild-Server)
                bilder_index[slug] = lokale_dateien
                bilder_gesamt += geladen
                if bilder_namen:
                    print(f"      → {geladen} neue Bilder geladen "
                          f"({len(lokale_dateien)} insgesamt für diesen Artikel)",
                          flush=True)
                # Index nach jedem Artikel speichern (Wiederaufnahme robust)
                INDEX_PFAD.write_text(
                    json.dumps(bilder_index, ensure_ascii=False, indent=2),
                    encoding="utf-8",
                )

            time.sleep(0.4)  # höflich gegenüber Wikipedia-API

    # Info-Datei
    INFO_PFAD.write_text(json.dumps({
        "buch": "Das Akkordeon und seine Geschichte",
        "quelle": "Wikipedia:Bücher/Das_Akkordeon_und_seine_Geschichte",
        "lizenz": "CC-BY-SA 4.0",
        "abgerufen_am": time.strftime("%Y-%m-%d"),
        "kapitel_struktur": [
            {"teil": t, "artikel": a} for t, a in KAPITEL
        ],
        "erfolgreich_abgerufen": erfolgreich,
        "fehlend": fehlend,
        "bilder_anzahl": sum(len(v) for v in bilder_index.values()),
    }, ensure_ascii=False, indent=2), encoding="utf-8")

    # Zusammenfassung
    print()
    print("=" * 64)
    print(f" {len(erfolgreich)} Artikel erfolgreich abgerufen")
    if fehlend:
        print(f" {len(fehlend)} Artikel nicht gefunden:")
        for f in fehlend:
            print(f"    - {f}")
    if not args.keine_bilder:
        bilder_dateien = list(BILD_DIR.glob("*"))
        gesamt_kb = sum(p.stat().st_size for p in bilder_dateien) / 1024
        print(f" {len(bilder_dateien)} Bilder gespeichert ({gesamt_kb:.0f} KB)")
    print("=" * 64)

    # ZIP packen
    if not args.kein_zip:
        print(f"\n>>> Packe {ZIP_PFAD.name} ...")
        if ZIP_PFAD.exists():
            ZIP_PFAD.unlink()
        with zipfile.ZipFile(ZIP_PFAD, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as zf:
            for p in sorted(WORK.rglob("*")):
                if p.is_file():
                    zf.write(p, p.relative_to(ROOT))
        groesse_mb = ZIP_PFAD.stat().st_size / (1024 * 1024)
        print(f"    {ZIP_PFAD}  ({groesse_mb:.1f} MB)")
        print()
        print("Diese ZIP-Datei kannst du jetzt im Chat hochladen,")
        print("damit Claude daraus das fertige Buch baut.")

    return 0


if __name__ == "__main__":
    sys.exit(main())
