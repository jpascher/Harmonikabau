#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
verarbeite_quellen.py
=====================
Verarbeitet die mit herunterladen.py erzeugten lokalen Quellen
(akkordeon_quellen/) zu fertigem LaTeX im KDP-6x9-Format.

Voraussetzung: pandoc + LuaLaTeX müssen installiert sein.
KEIN Internet nötig — alle Texte und Bilder sind bereits lokal.

Aufruf:
  python3 verarbeite_quellen.py [--quellen PFAD] [--ausgabe PFAD] [--kompilieren]

Standardpfade:
  --quellen   ./akkordeon_quellen        (Eingabe vom Download-Skript)
  --ausgabe   ./buch_ausgabe             (LaTeX-Dateien + PDF)
"""

import argparse
import json
import re
import shutil
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Buchstruktur — exakt nach der Wikipedia-Buchseite
# (Christian Messner ist nicht mehr als eigenständiger Wikipedia-Artikel
# vorhanden und wird automatisch übersprungen.)
# ---------------------------------------------------------------------------

BUCH_TITEL = "Das Akkordeon und seine Geschichte"
BUCH_AUTOR = "zusammengestellt von Johann Pascher"

KAPITEL: List[Tuple[str, List[str]]] = [
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


# ---------------------------------------------------------------------------
# Hilfsfunktionen
# ---------------------------------------------------------------------------

def slugify(s: str) -> str:
    s = s.lower()
    s = (s.replace("ä", "ae").replace("ö", "oe").replace("ü", "ue")
          .replace("ß", "ss"))
    s = re.sub(r"[^a-z0-9]+", "_", s).strip("_")
    return s


def escape_latex(s: str) -> str:
    repl = {
        "&": r"\&", "%": r"\%", "$": r"\$", "#": r"\#",
        "_": r"\_", "{": r"\{", "}": r"\}",
        "~": r"\textasciitilde{}", "^": r"\textasciicircum{}",
    }
    return "".join(repl.get(c, c) for c in s)


def strip_braces(text: str, open_tok: str, close_tok: str) -> str:
    out = []
    i = 0
    depth = 0
    while i < len(text):
        if text[i:i + len(open_tok)] == open_tok:
            depth += 1
            i += len(open_tok)
        elif text[i:i + len(close_tok)] == close_tok and depth > 0:
            depth -= 1
            i += len(close_tok)
        else:
            if depth == 0:
                out.append(text[i])
            i += 1
    return "".join(out)


# ---------------------------------------------------------------------------
# Wikitext-Vorbereitung
# ---------------------------------------------------------------------------

def vorbereite_wikitext(text: str) -> str:
    """Räumt Wikitext auf, bevor pandoc übersetzt."""
    # 0) <gallery>-Bloecke aufloesen — pandoc kann das mediawiki-Konstrukt
    #    nicht. Jede Zeile hat das Format `Dateiname.ext|Bildunterschrift`.
    #    Wir wandeln das in eine Reihe von [[File:...]]-Bildern um, damit
    #    der spaetere Bilder-Schritt sie wie Einzelbilder behandelt.
    def _gallery_aufloesen(m: re.Match) -> str:
        inhalt = m.group(1)
        zeilen = []
        for raw in inhalt.split("\n"):
            zeile = raw.strip()
            if not zeile:
                continue
            # `Dateiname.ext|Beschreibung` oder nur `Dateiname.ext`
            if "|" in zeile:
                fn, _, _beschr = zeile.partition("|")
                fn = fn.strip()
            else:
                fn = zeile
            if not fn:
                continue
            # Sicherheit: keine eckigen Klammern im Dateinamen
            fn = fn.replace("[", "").replace("]", "")
            zeilen.append(f"[[File:{fn}|thumb]]")
        # Mit Leerzeilen davor/danach, damit pandoc den Block als Bilder nimmt
        return "\n\n" + "\n\n".join(zeilen) + "\n\n"

    text = re.sub(
        r"(?is)<gallery[^>]*?>(.*?)</gallery>",
        _gallery_aufloesen,
        text,
    )

    # 1) Vorlagen {{...}} entfernen (rekursiv)
    text = strip_braces(text, "{{", "}}")
    # 2) Tabellen entfernen
    text = re.sub(r"(?ms)^\{\|.*?^\|\}\s*$", "", text)
    # 3) <ref>-Belege entfernen
    #    WICHTIG: nach der Entfernung darf kein einsames Leerzeichen am
    #    Zeilenanfang bleiben — sonst erkennt unsere spaetere Heuristik den
    #    Text faelschlicherweise als Verbatim-Block. Die Loesung: bei einem
    #    Leerzeichen unmittelbar VOR und nach <ref>...</ref> wird das eine
    #    einbehalten, das andere weg.
    text = re.sub(r"(?is)<ref[^>]*?/>", "", text)
    text = re.sub(r"(?is)<ref[^>]*?>.*?</ref>", "", text)
    # Doppelte Leerzeichen direkt nach der Entfernung normalisieren
    text = re.sub(r"  +", " ", text)
    # Leerzeichen am Zeilenanfang entfernen, falls direkt nach Entfernung
    # entstanden (Pattern: " text..." am Zeilenstart, max 1 Leerzeichen).
    # Das raeumt Reste auf, ohne echte Wiki-Praeformat-Bloecke zu zerstoeren
    # (die haben mehrere Zeilen mit Leerzeichen, oder beginnen direkt nach
    # einem Absatzwechsel — die Heuristik filtert in Schritt 11 dann sauber).
    # Hier wirken wir prophylaktisch: wenn ein Leerzeichen direkt VOR einem
    # alphanumerischen Zeichen am Zeilenanfang steht UND die Zeile davor
    # KEIN Leerzeichen-Block-Anfang ist, dann strippen.
    zeilen_tmp = text.split("\n")
    for j in range(len(zeilen_tmp)):
        if (zeilen_tmp[j].startswith(" ")
                and not zeilen_tmp[j].startswith("    ")
                and len(zeilen_tmp[j]) > 1
                and zeilen_tmp[j][1].isalnum()):
            # Nur wenn die VORHERGEHENDE Zeile NICHT auch mit Leerzeichen anfaengt
            if j == 0 or not zeilen_tmp[j-1].startswith(" "):
                zeilen_tmp[j] = zeilen_tmp[j].lstrip(" ")
    text = "\n".join(zeilen_tmp)
    # 4) HTML-Kommentare
    text = re.sub(r"(?s)<!--.*?-->", "", text)
    # 5) Bilder normalisieren UND Dateinamen sofort auf den slugified Namen
    #    setzen, weil sonst pandoc bei Apostrophen, Kommas, Klammern und
    #    anderen Sonderzeichen im Dateinamen das Bild zerlegt und der
    #    Bildname als Klartext durchrutscht.
    #
    #    WICHTIG: Datei-Direktiven koennen verschachtelte [[Links]] in der
    #    Bildunterschrift haben, z.B.:
    #       [[Datei:Foto.jpg|mini|Eine Zeichnung von [[Karl Frech]]]]
    #    Ein simpler Regex schliesst dann zu frueh am inneren ]], und der
    #    Rest "(wurde abgerissen).]]" rutscht als Klartext durch. Wir nutzen
    #    daher einen balancierten Parser fuer [[...]].
    def _datei_normalisieren_balanced(s: str) -> str:
        out = []
        i = 0
        L = len(s)
        while i < L:
            # Suche naechstes [[
            idx = s.find("[[", i)
            if idx < 0:
                out.append(s[i:])
                break
            out.append(s[i:idx])
            # Steht hier eine Datei-Direktive?
            kopf_end = s.find(":", idx + 2, idx + 30)
            ist_datei = False
            if kopf_end > 0:
                kopf = s[idx + 2:kopf_end]
                if kopf in ("Datei", "Bild", "File"):
                    ist_datei = True
            if not ist_datei:
                # Normales [[...]]: einfach durchreichen, weiterer Code uebersetzt
                out.append("[[")
                i = idx + 2
                continue
            # Balanciert das Ende der Datei-Direktive finden
            depth = 1
            j = idx + 2
            while j < L:
                if s[j:j+2] == "[[":
                    depth += 1
                    j += 2
                elif s[j:j+2] == "]]":
                    depth -= 1
                    if depth == 0:
                        break
                    j += 2
                else:
                    j += 1
            if j >= L:
                # Unbalanciert - sicherheitshalber durchreichen
                out.append("[[")
                i = idx + 2
                continue
            # Wir haben den balancierten Block gefunden: s[idx:j+2]
            innen = s[idx + 2:j]   # ohne aussere [[ ]]
            # Dateiname ist alles bis zum ersten | (top-level)
            tiefe = 0
            sep = -1
            for k in range(len(innen)):
                if innen[k:k+2] == "[[":
                    tiefe += 1
                elif innen[k:k+2] == "]]":
                    tiefe -= 1
                elif innen[k] == "|" and tiefe == 0:
                    sep = k
                    break
            if sep < 0:
                roh = innen
            else:
                roh = innen[:sep]
            # Praefix "Datei:" oder "Bild:" oder "File:" wegnehmen
            if ":" in roh:
                roh = roh.split(":", 1)[1]
            roh = roh.strip()
            if not roh:
                # Defekt — ueberspringen
                i = j + 2
                continue
            slug_name = slugify(Path(roh).stem) + Path(roh).suffix.lower()
            out.append(f"[[File:{slug_name}|thumb]]")
            i = j + 2
        return "".join(out)

    text = _datei_normalisieren_balanced(text)
    # 6) Kategorien
    text = re.sub(r"(?m)^\[\[(?:Kategorie|Category):[^\]]+\]\]\s*$", "", text)
    # 7) Interwiki-Links
    text = re.sub(r"(?m)^\[\[[a-z][a-z\-]{1,11}:[^\]]+\]\]\s*$", "", text)

    # 8) Interne Wiki-Links → Klartext, aber File-Einbindungen schützen
    def _link_zu_text(m: re.Match) -> str:
        inner = m.group(1)
        if inner.startswith(("File:", "Datei:", "Bild:")):
            return m.group(0)
        if "|" in inner:
            return inner.split("|", 1)[1]
        return inner
    text = re.sub(r"\[\[([^\[\]]+?)\]\]", _link_zu_text, text)

    # 9) Externe Links → Anzeigetext
    text = re.sub(r"\[https?://\S+[ \t]+([^\]\n]+)\]", r"\1", text)
    text = re.sub(r"\[https?://[^\s\]]+\]", "", text)

    # 9b) Eigenstaendige Video-/Audio-Verweise entfernen, deren Link nicht
    #     einbettbar ist. "Video eines neueren rocking melodeon." als
    #     einzelne Zeile macht ohne Link keinen Sinn.
    text = re.sub(
        r"(?im)^\s*(?:Video|Hoerbeispiel|Hörbeispiel|Audio|Audiobeispiel)\b[^\n]*\.?\s*$",
        "",
        text,
    )

    # 10) Mehrere Leerzeilen reduzieren
    text = re.sub(r"\n{3,}", "\n\n", text)

    # 11) Fuehrende Leerzeichen am Zeilenanfang behandeln.
    #     Wikipedia-Wikitext rendert Zeilen mit fuehrendem Leerzeichen als
    #     vorformatierten Block (so wie HTML <pre>). Pandoc erkennt das mit
    #     dem mediawiki-Reader nicht zuverlaessig — wir wandeln Bloecke aus
    #     2+ aufeinanderfolgenden Zeilen mit fuehrendem Leerzeichen explizit
    #     in <pre>...</pre> um, was pandoc als Verbatim erkennt.
    #     Einzelne Leerzeichen-Zeilen (Reste der <ref>-Entfernung) werden
    #     entfernt.
    zeilen = text.split("\n")
    ergebnis = []
    i = 0

    def _bricht_lange_zeile(z: str, max_len: int = 70) -> List[str]:
        """Bricht eine zu lange Verbatim-Zeile an gutem Trennpunkt um.
        Nur Zeilen ueber max_len werden bearbeitet — kuerzere bleiben unveraendert.
        """
        if len(z) <= max_len:
            return [z]
        # Bevorzugte Trennstellen — immer am sinnvollsten an Klammerkommentaren
        # wie " (da kein...". Wir suchen sie nur bei zu langen Zeilen.
        for trennzeichen in [" (da ", " ←", " (auch ", " ("]:
            if trennzeichen in z:
                idx = z.find(trennzeichen)
                # Nur trennen, wenn die erste Haelfte mindestens 30 Zeichen
                # hat und die zweite Haelfte auch nicht trivial kurz ist.
                if 30 < idx < max_len:
                    erste = z[:idx]
                    zweite = "  " + z[idx:].lstrip()
                    return [erste, zweite]
        # Notfall: hart umbrechen
        return [z[:max_len], "  " + z[max_len:]]

    while i < len(zeilen):
        z = zeilen[i]
        if z.startswith(" ") and not z.startswith("    "):
            # WICHTIG: nur als Verbatim-Block erkennen, wenn die Zeile
            # DAVOR leer war oder es die erste Zeile ist. Sonst ist die
            # eingerueckte Zeile vermutlich nur eine Pandoc-Wrap-
            # Fortsetzung der vorigen Zeile und KEIN echter Wikipedia-
            # Praeformat-Block.
            ist_block_anfang = (i == 0) or (zeilen[i-1].strip() == "")
            if not ist_block_anfang:
                # Wahrscheinlich Pandoc-Wrap-Fortsetzung — Leerzeichen weg
                # und an die vorige Zeile anhaengen
                if ergebnis and ergebnis[-1]:
                    ergebnis[-1] = ergebnis[-1].rstrip() + " " + z.lstrip(" ")
                else:
                    ergebnis.append(z.lstrip(" "))
                i += 1
                continue
            block = []
            while i < len(zeilen) and zeilen[i].startswith(" ") and not zeilen[i].startswith("    "):
                # Wikipedia-Fett/Kursiv-Marker innerhalb von <pre>-Bloecken
                # entfernen (''' und '') — sie wuerden als wortwoertliche
                # Apostrophe gerendert, was im Verbatim haesslich aussieht
                zeile = zeilen[i].lstrip(" ")
                zeile = zeile.replace("'''", "").replace("''", "")
                # HTML-Tags innerhalb von Verbatim-Bloecken stehen sonst
                # woertlich da. <sup>7</sup> -> ^7 (akkordeontypisch),
                # <sub>1</sub> -> _1, andere Tags (u, b, i, em) ganz raus.
                zeile = re.sub(r"<sup>([^<]*)</sup>", r"^\1", zeile)
                zeile = re.sub(r"<sub>([^<]*)</sub>", r"_\1", zeile)
                zeile = re.sub(r"</?[a-zA-Z][^>]*>", "", zeile)
                # &nbsp; und andere HTML-Entities
                zeile = zeile.replace("&nbsp;", " ")
                zeile = zeile.replace("&amp;", "&")
                zeile = zeile.replace("&lt;", "<")
                zeile = zeile.replace("&gt;", ">")
                # Lange Zeilen vorab umbrechen
                for teilzeile in _bricht_lange_zeile(zeile):
                    block.append(teilzeile)
                i += 1
            if len(block) >= 2:
                ergebnis.append("")
                ergebnis.append("<pre>")
                ergebnis.extend(block)
                ergebnis.append("</pre>")
                ergebnis.append("")
            else:
                ergebnis.extend(block)
            continue
        ergebnis.append(z)
        i += 1
    text = "\n".join(ergebnis)

    # 12) Mehrfache Leerzeichen mitten im Text auf eines reduzieren
    #     (entstehen ebenfalls nach <ref>-Entfernung)
    text = re.sub(r"(?<=\S) {2,}(?=\S)", " ", text)

    # 12b) Verwaiste Medien-Verweise entfernen, die nach <ref>-Entfernung
    #      keine Funktion mehr haben. Pattern: kurze einzelne Zeilen, die
    #      mit "Video", "Audio", "Hoerbeispiel", "Video eines" anfangen und
    #      lediglich ankuendigen, was eigentlich verlinkt war (jetzt aber
    #      nicht mehr darstellbar im Buch).
    text = re.sub(
        r"(?m)^(?:Video(?: eines| beispiel| von)?|Audio(?: beispiel)?"
        r"|H(?:ö|oe)rbeispiel(?:e)?|Klangbeispiel(?:e)?)"
        r"\s+[^\n]{0,200}\.\s*$",
        "",
        text,
    )

    # 13) Leere Listenpunkte entfernen — entstehen, wenn der Inhalt eines
    #     `* {{Literatur|...}}` durch die Vorlagen-Entfernung leer wurde.
    #     Pattern: eine Zeile, die nur aus *, ** oder *** + Whitespace besteht
    text = re.sub(r"(?m)^\s*\*+\s*$", "", text)
    #     Auch Listenpunkte mit nur "()" oder ", " als Rest (Reste von Links)
    text = re.sub(r"(?m)^(\s*\*+\s*)[\(\)\.,;:\s\-–—]+$", "", text)

    # 14) Leere Sektionen entfernen — wenn nach Saeuberung unter einer
    #     "== Literatur ==" oder "== Weblinks ==" gar nichts mehr steht.
    #     Ein Header gefolgt von Leerzeilen und dem naechsten Header → weg.
    for _ in range(3):  # mehrfach, falls mehrere leere Sektionen aufeinander
        text = re.sub(
            r"(?m)^(==+)\s*([^=\n]+?)\s*\1\s*$\n+(?=^==+)",
            "",
            text,
        )
        # Auch leere Sektion am Dateiende
        text = re.sub(
            r"(?ms)^(==+)\s*([^=\n]+?)\s*\1\s*$\n+\Z",
            "",
            text,
        )

    # 15) Mehrfache Leerzeilen erneut reduzieren (nach den Loeschungen)
    text = re.sub(r"\n{3,}", "\n\n", text)

    return text.strip()


# ---------------------------------------------------------------------------
# Pandoc-Konvertierung
# ---------------------------------------------------------------------------

def wikitext_zu_latex(wikitext: str, vorhandene_bilder: List[str]) -> str:
    """Konvertiert Wikitext mit pandoc nach LaTeX."""
    proc = subprocess.run(
        ["pandoc", "-f", "mediawiki", "-t", "latex", "--wrap=preserve"],
        input=wikitext.encode("utf-8"),
        capture_output=True,
        check=False,
    )
    if proc.returncode != 0:
        print("      pandoc-Fehler:",
              proc.stderr.decode("utf-8", "ignore")[:300], file=sys.stderr)
        return ""
    return saeubere_latex(proc.stdout.decode("utf-8"), vorhandene_bilder)


def saeubere_latex(tex: str, vorhandene_bilder: List[str]) -> str:
    """LaTeX nach pandoc aufräumen und Bilder auf lokale Dateien biegen."""

    # 1) figure-Blöcke: Bild-Datei einbinden, falls lokal vorhanden
    def _figure_block(m: re.Match) -> str:
        block = m.group(0)
        gm = re.search(r"\\includegraphics(?:\[[^\]]*\])?\{([^}]+)\}", block)
        if not gm:
            return ""
        return _bild_block(gm.group(1), vorhandene_bilder)
    tex = re.sub(r"\\begin\{figure\}.*?\\end\{figure\}",
                 _figure_block, tex, flags=re.DOTALL)

    # 2) Reine \includegraphics{...} (ohne figure-Wrapper)
    tex = re.sub(r"\\includegraphics(\[[^\]]*\])?\{(?!bilder/)([^}]+)\}",
                 lambda m: _bild_block(m.group(2), vorhandene_bilder),
                 tex)

    # 3) Heading-Tiefen
    tex = re.sub(r"\\section\{", r"\\section*{", tex)
    tex = re.sub(r"\\subsection\{", r"\\subsection*{", tex)
    tex = re.sub(r"\\subsubsection\{", r"\\subsubsection*{", tex)

    # 4) Pandoc-\hypertarget vor Überschriften wegputzen
    tex = re.sub(r"\\hypertarget\{[^}]*\}\{%?\s*"
                 r"(\\(?:sub){0,2}section\*?\{[^}]*\})\\label\{[^}]*\}\}",
                 r"\1", tex)

    # 5) Reste von leeren externen Links bereinigen
    tex = re.sub(r"\bund\s*([.,;:])", r"\1", tex)
    tex = re.sub(r"\boder\s*([.,;:])", r"\1", tex)
    tex = re.sub(r"\(\s*\)", "", tex)
    tex = re.sub(r" {2,}", " ", tex)
    tex = re.sub(r" +([.,;:!?])", r"\1", tex)

    # 6) \includesvg verwerfen (pandoc baut das ein; LaTeX ohne svg-Paket
    #    versteht es nicht und der Dateiname rutscht in den Text)
    tex = re.sub(r"\\includesvg(?:\[[^\]]*\])?\{[^}]*\}", "", tex)

    # 7) "Findet sich der Satz:" oder "schreibt:" am Zeilenende entfernen,
    #    wenn das nachfolgende Zitat (eine Vorlage) bereits weggeputzt wurde.
    #    Wir erkennen das daran, dass nach dem Doppelpunkt ein Absatz folgt
    #    oder eine Leerzeile.
    tex = re.sub(
        r"([^\n]*?\bder Satz):\s*\n\s*\n",
        "",
        tex,
    )
    # Allgemeiner: Halbsätze, die mit "...Satz:" oder "...heißt es:" enden und
    # dann ohne Inhalt in eine Leerzeile/neuen Absatz übergehen
    tex = re.sub(
        r"(?m)^(.*?(?:der Satz|heißt es|schreibt|notiert|formuliert)):\s*$",
        "",
        tex,
    )

    # 8) Leere \item entfernen — entstehen, wenn der Inhalt eines
    #    Listenpunkts eine Vorlage war, die wir vorher entsorgt haben.
    #    Pandoc rendert leere Items als nur "\item\n" mit nichts dahinter.
    #    Achtung: Pandoc setzt den Inhalt eingerueckt UNTER das \item, also
    #    "\item\n  Inhalt\n". Daher ist \item nur leer, wenn die Zeile direkt
    #    danach selbst mit \item, \end{itemize}, \end{enumerate} oder einer
    #    Leerzeile + Sektionswechsel beginnt — NICHT wenn dort eingeruekter
    #    Inhalt steht.
    tex = re.sub(
        r"(?m)^\\item[ \t]*\n(?=\\item\b|\\end\{itemize\}|\\end\{enumerate\})",
        "",
        tex,
    )
    # \item gefolgt von einer Zeile, die nur Whitespace und Trennzeichen ist
    tex = re.sub(
        r"(?m)^\\item[ \t]*\n[ \t]+[\(\)\.,;:\-–—\s]+\n",
        "\n",
        tex,
    )
    # \item gefolgt von einer Zeile, die nur " -- Beschreibungstext" ist
    # (entsteht wenn der Anker/URL der Quelle vor dem " --" weg ist)
    tex = re.sub(
        r"(?m)^(\\item[ \t]*\n[ \t]+)--\s+",
        r"\1",
        tex,
    )

    # 9) Leere itemize-/enumerate-Bloecke entfernen
    #    \begin{itemize}\n\tightlist\n\end{itemize}  oder ohne tightlist
    for env in ("itemize", "enumerate"):
        tex = re.sub(
            r"\\begin\{" + env + r"\}\s*(?:\\tightlist\s*)?\\end\{" + env + r"\}",
            "",
            tex,
        )

    # 10) \subsection*{...}-Header, denen direkt ein anderer Header oder
    #     das Kapitelende folgt (also ohne Inhalt), entfernen.
    for _ in range(3):
        tex = re.sub(
            r"(?m)^\\subsection\*\{[^}]+\}\s*\n+(?=\\(?:sub)*section\*\{|\\medskip|\\end\{document\}|\\chapter\{)",
            "",
            tex,
        )
        tex = re.sub(
            r"(?m)^\\subsubsection\*\{[^}]+\}\s*\n+(?=\\(?:sub)*section\*\{|\\medskip|\\end\{document\}|\\chapter\{)",
            "",
            tex,
        )

    # 11) Mehrfache Leerzeilen
    tex = re.sub(r"\n{3,}", "\n\n", tex)

    # 12) Verbatim-Bloecke — kleinere Schrift + Auto-Umbruch wenn zu lang.
    #     \footnotesize ist gut lesbar. Bei \footnotesize passen ca. 70
    #     Zeichen pro Zeile. Laengere Zeilen werden automatisch via
    #     fancyvrb-breaklines umgebrochen.
    def _verbatim_klein(m: re.Match) -> str:
        return ("{\\footnotesize\n"
                "\\begin{Verbatim}[breaklines=true,breakanywhere=true,"
                "breaksymbol=]"
                + m.group(1) +
                "\\end{Verbatim}\n"
                "}")
    tex = re.sub(
        r"\\begin\{verbatim\}(.*?)\\end\{verbatim\}",
        _verbatim_klein,
        tex,
        flags=re.DOTALL,
    )

    return tex


def _bild_block(name: str, vorhandene_bilder: List[str]) -> str:
    """Erzeugt LaTeX-Block für ein Bild — nur wenn lokal vorhanden."""
    base = re.sub(r"^File:|^Datei:|^Bild:", "", name).strip().replace(" ", "_")
    sicher = slugify(Path(base).stem) + Path(base).suffix.lower()
    if sicher not in vorhandene_bilder:
        return ""  # Bild nicht lokal vorhanden → einfach weglassen
    pfad = f"bilder/{sicher}"
    # Bilder kompakt halten: width=65%, max height=40% Textheight.
    # Das laesst genug Platz fuer Text auf derselben Seite und reduziert
    # Float-Sprünge auf Folgeseiten mit Loechern.
    return (f"\\begin{{center}}\n"
            f"\\includegraphics[width=.65\\linewidth,"
            f"max height=0.40\\textheight,keepaspectratio]{{{pfad}}}\n"
            f"\\end{{center}}\n")


# ---------------------------------------------------------------------------
# LaTeX-Vorlagen
# ---------------------------------------------------------------------------

PREAMBLE = r"""% =============================================================
% preamble.tex — KDP-Format 6"x9" (Inter-Schrift)
% Buch: "Das Akkordeon und seine Geschichte"
% Kompilation: lualatex (zweimal für Inhaltsverzeichnis)
% =============================================================

\documentclass[11pt,twoside,openright]{book}

% --- Schriften (Inter / DejaVu Mono) --------------------------
\RequirePackage{fontspec}
\setmainfont{Inter}[
  Scale=1.02,
  UprightFont=*-Regular,
  BoldFont=*-Bold,
  ItalicFont=*-Italic,
  BoldItalicFont=*-BoldItalic,
  Ligatures=TeX
]
\setsansfont{Inter}[Scale=MatchLowercase,Ligatures=TeX]
\setmonofont{DejaVu Sans Mono}[Scale=0.88]

% --- Sprache --------------------------------------------------
\usepackage{babel}
\babelprovide[main,import]{ngerman}

% --- Mikrotypografie ------------------------------------------
\usepackage{microtype}
\frenchspacing
\emergencystretch=3em
\tolerance=2500
\hbadness=10000
\hfuzz=2pt
\clubpenalty=10000
\widowpenalty=10000

% --- Geometrie: KDP-Paperback 6" x 9" -------------------------
\usepackage[
  paperwidth=152.4mm,
  paperheight=228.6mm,
  inner=19mm,
  outer=14mm,
  top=18mm,
  bottom=18mm,
  headheight=14pt,
  headsep=6mm,
  footskip=10mm
]{geometry}

% --- Standardpakete -------------------------------------------
\usepackage{xcolor}
\usepackage{graphicx}
\usepackage[export]{adjustbox}
\graphicspath{{./}{bilder/}}
\usepackage{amsmath,amssymb}
\usepackage{enumitem}
\usepackage{longtable,booktabs,array,tabularx}
\usepackage{calc}
\usepackage{caption}
\captionsetup{font=small,labelfont=bf,labelsep=period}
\usepackage{parskip}

% --- Verbatim mit Auto-Umbruch ------------------------------------
% fancyvrb erlaubt \begin{Verbatim}[breaklines,...] mit Zeilenumbruch
\usepackage{fancyvrb}
\fvset{breaklines=true,breakanywhere=true}
\setlength{\parskip}{0.5ex plus 0.2ex minus 0.1ex}
\usepackage{setspace}
\setstretch{1.15}

% --- Farben (T0-Stil) -----------------------------------------
\definecolor{t0blue}{RGB}{0,70,127}
\definecolor{t0green}{RGB}{0,120,60}
\definecolor{t0red}{RGB}{180,0,0}
\definecolor{t0gray}{RGB}{80,80,80}

% --- Hyperlinks -----------------------------------------------
\usepackage[unicode,colorlinks=true,
            linkcolor=t0blue,
            citecolor=t0green,
            urlcolor=t0blue,
            bookmarksnumbered=true]{hyperref}
\usepackage{xurl}

% --- Kapitelüberschriften -------------------------------------
\usepackage{titlesec}
\titleformat{\chapter}[display]
  {\normalfont\sffamily\bfseries\color{t0blue}}
  {\filright \fontsize{32}{40}\selectfont \chaptertitlename\ \thechapter}
  {0.6ex}
  {\titlerule[0.6pt]\vspace{0.5ex}\filright\Large}
\titlespacing*{\chapter}{0pt}{-20pt}{18pt}

\titleformat{\part}[display]
  {\normalfont\sffamily\Huge\bfseries\filcenter\color{t0blue}}
  {\partname\ \thepart}{16pt}
  {\Huge}
\titlespacing*{\part}{0pt}{*4}{*4}

\titleformat{\section}{\normalfont\sffamily\large\bfseries\color{t0blue}}{\thesection}{0.7em}{}
\titleformat{\subsection}{\normalfont\sffamily\normalsize\bfseries\color{t0blue!85!black}}{\thesubsection}{0.6em}{}
\titleformat{\subsubsection}{\normalfont\sffamily\normalsize\bfseries}{\thesubsubsection}{0.5em}{}

% --- Kopf-/Fußzeilen ------------------------------------------
\usepackage{fancyhdr}
\pagestyle{fancy}
\fancyhf{}
\fancyhead[LE]{\small\textcolor{t0gray}{\thepage \quad \textit{Das Akkordeon und seine Geschichte}}}
\fancyhead[RO]{\small\textcolor{t0gray}{\nouppercase{\textit{\leftmark}} \quad \thepage}}
\renewcommand{\headrulewidth}{0.3pt}

\fancypagestyle{plain}{%
  \fancyhf{}%
  \fancyfoot[C]{\small\textcolor{t0gray}{\thepage}}%
  \renewcommand{\headrulewidth}{0pt}%
}

% --- Inhaltsverzeichnis ---------------------------------------
\setcounter{tocdepth}{1}
\setcounter{secnumdepth}{0}

% --- Pandoc-Hilfsbefehle --------------------------------------
\providecommand{\tightlist}{\setlength{\itemsep}{0pt}\setlength{\parskip}{0pt}}
\providecommand{\passthrough}[1]{\texttt{#1}}
\providecommand{\textsubscript}[1]{\ensuremath{_{#1}}}
\providecommand{\pandocbounded}[1]{#1}

\setkeys{Gin}{keepaspectratio}
\usepackage{etoolbox}
"""


TITELSEITE = r"""\thispagestyle{empty}
\begin{titlepage}
\centering
\vspace*{2.5cm}
{\sffamily\itshape\color{t0gray}\Large
Eine Sammlung von Wikipedia-Artikeln\par}
\vspace{1.6cm}
{\sffamily\bfseries\color{t0blue}\fontsize{28}{34}\selectfont
Das Akkordeon\\[0.3em]
und seine Geschichte\par}
\vspace{2.0cm}
{\centering\rule{6cm}{0.4pt}\par}
\vspace{1.2cm}
{\sffamily\large
zusammengestellt von\par}
\vspace{0.5em}
{\sffamily\bfseries\Large
Johann Pascher\par}
\vfill
{\sffamily\small\color{t0gray}
Inhaltliche Quelle: deutschsprachige Wikipedia\\[0.3em]
Lizenz: Creative Commons Attribution--ShareAlike 4.0 International\\[0.3em]
\textsc{cc-by-sa 4.0}\par}
\vspace{0.6cm}
{\sffamily\small\color{t0gray}
Zusammenstellung: \today\par}
\end{titlepage}

\clearpage
\thispagestyle{empty}
~\vfill
{\footnotesize\color{t0gray}
Diese Buchausgabe ist eine private Zusammenstellung von Artikeln aus der
deutschsprachigen Wikipedia. Die Artikel stehen unter der Creative-Commons-Lizenz
\textsc{cc-by-sa 4.0}; ihre Autorenschaft ist die Gemeinschaft der jeweiligen
Wikipedia-Bearbeiter. Die vollständige Versionsgeschichte ist über die
Quellenlinks am Ende jedes Kapitels zugänglich.\par}
\clearpage
"""


LIZENZSEITE = r"""\thispagestyle{empty}

\begin{center}
{\sffamily\Large\bfseries\color{t0blue} Über dieses Buch\par}
\end{center}
\vspace{0.5cm}

Dieses Buch versammelt eine Reihe von Artikeln aus der deutsch\-sprachigen
Wikipedia, die in der dortigen Buchfunktion unter dem Titel
\textit{Das Akkordeon und seine Geschichte} (Benutzer:Jpascher) als
Lese\-sammlung gespeichert sind. Da die Wikipedia-eigene PDF-Erzeugung
seit der Abschaltung des OCG-Renderers (2017) nur eingeschränkt
funk\-tioniert, ist diese Sammlung als eigenständiges, dauerhaft
verfügbares Buchwerk angelegt.

\vspace{0.6cm}
\begin{center}
{\sffamily\bfseries\color{t0blue} Lizenz und Nachnutzung\par}
\end{center}
\vspace{0.2cm}

Die Inhalte aller Artikel stehen unter der Lizenz \\
\textbf{Creative Commons Attribution -- ShareAlike 4.0 International} (\textbf{CC-BY-SA 4.0}).

Lizenzbedingungen: \\ \url{https://creativecommons.org/licenses/by-sa/4.0/deed.de}

Die Autorenschaft an einem Artikel ist die Gemeinschaft seiner
Wikipedia-Bearbeiter. Die vollständige Versionsgeschichte mit allen
Beitragenden ist über den jeweils am Kapitelende angegebenen Link
einsehbar. Bei Weitergabe oder Bearbeitung dieses Buches müssen die
Lizenzbedingungen erhalten bleiben.

\vspace{0.6cm}
\begin{center}
{\sffamily\bfseries\color{t0blue} Hinweise zur Darstellung\par}
\end{center}
\vspace{0.2cm}

Vorlagen, Infoboxen, Einzelnachweise und Tabellen, die in der
Wikipedia-Online-Darstellung enthalten sind, wurden bei der
Konvertierung gekürzt oder entfernt, weil sie sich nur eingeschränkt
in eine Buchform übertragen lassen. Der Fließtext der Artikel ist
unverändert übernommen. Bilder, deren Lizenz dies zuließ, wurden
heruntergeladen und im Buch eingebunden.

\clearpage
"""


# ---------------------------------------------------------------------------
# Build
# ---------------------------------------------------------------------------

def baue_kapitel(titel: str, quellen: Path, kapitel_dir: Path,
                 bilder_index: Dict[str, List[str]]) -> Optional[Path]:
    """Wandelt einen Wikitext-Artikel in ein LaTeX-Kapitel um."""
    slug = slugify(titel)
    wiki_pfad = quellen / "artikel" / f"{slug}.wiki"
    if not wiki_pfad.exists():
        print(f"      NICHT VORHANDEN: {wiki_pfad.name}", file=sys.stderr)
        return None

    print(f"   • {titel}", flush=True)

    wikitext = wiki_pfad.read_text(encoding="utf-8")
    aufbereitet = vorbereite_wikitext(wikitext)

    vorhandene_bilder = bilder_index.get(slug, [])
    latex_inhalt = wikitext_zu_latex(aufbereitet, vorhandene_bilder)
    if not latex_inhalt.strip():
        print("      Leere Konvertierung — übersprungen.", file=sys.stderr)
        return None

    quelle_url = "https://de.wikipedia.org/wiki/" + urllib.parse.quote(
        titel.replace(" ", "_"))
    pfad = kapitel_dir / f"{slug}.tex"
    pfad.write_text(
        f"\\chapter{{{escape_latex(titel)}}}\n"
        f"\\label{{chap:{slug}}}\n\n"
        f"{latex_inhalt}\n\n"
        f"\\medskip\n"
        f"\\begin{{flushright}}\\footnotesize\\itshape\n"
        f"Quelle: \\href{{{quelle_url}}}{{{escape_latex(titel)} -- Wikipedia}}, "
        f"Lizenz: CC-BY-SA 4.0.\n"
        f"\\end{{flushright}}\n",
        encoding="utf-8",
    )
    return pfad


def schreibe_buch_tex(buch_tex: Path, kapitel_pfade: List[Tuple[str, List[Path]]],
                      ausgabe: Path) -> None:
    """Erzeugt die Master-LaTeX-Datei buch.tex."""
    lines = [
        PREAMBLE,
        "",
        f"\\title{{{BUCH_TITEL}}}",
        f"\\author{{{BUCH_AUTOR}}}",
        "",
        "\\begin{document}",
        "",
        TITELSEITE,
        "",
        LIZENZSEITE,
        "",
        "\\tableofcontents",
        "\\cleardoublepage",
        "",
    ]
    for teil_titel, pfade in kapitel_pfade:
        if not pfade:
            continue
        lines.append(f"\\part{{{escape_latex(teil_titel)}}}")
        for p in pfade:
            rel = p.relative_to(ausgabe).as_posix()
            lines.append(f"\\input{{{rel}}}")
        lines.append("")
    lines.append("\\end{document}")
    buch_tex.write_text("\n".join(lines), encoding="utf-8")


def kompiliere(ausgabe: Path) -> bool:
    print("\n>>> LuaLaTeX-Lauf 1 von 2 ...")
    rc = subprocess.run(
        ["lualatex", "-interaction=nonstopmode", "buch.tex"],
        cwd=ausgabe,
    ).returncode
    if rc != 0:
        print(f"   (LuaLaTeX-Lauf 1 hatte Warnungen — siehe buch.log)",
              file=sys.stderr)
    print("\n>>> LuaLaTeX-Lauf 2 von 2 ...")
    rc = subprocess.run(
        ["lualatex", "-interaction=nonstopmode", "buch.tex"],
        cwd=ausgabe,
    ).returncode
    pdf = ausgabe / "buch.pdf"
    if pdf.exists():
        print(f"\n>>> Fertig: {pdf}")
        return True
    print("\nFEHLER: kein PDF erzeugt — siehe buch.log", file=sys.stderr)
    return False


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quellen", type=Path,
                        default=Path("akkordeon_quellen"),
                        help="Pfad zum Quellenverzeichnis (Default: ./akkordeon_quellen)")
    parser.add_argument("--ausgabe", type=Path,
                        default=Path("buch_ausgabe"),
                        help="Ausgabeverzeichnis (Default: ./buch_ausgabe)")
    parser.add_argument("--kompilieren", action="store_true",
                        help="Nach dem LaTeX-Erzeugen LuaLaTeX zweimal aufrufen.")
    args = parser.parse_args()

    quellen = args.quellen.resolve()
    ausgabe = args.ausgabe.resolve()

    if not quellen.exists():
        print(f"FEHLER: Quellenverzeichnis nicht gefunden: {quellen}",
              file=sys.stderr)
        print("Hinweis: Erst herunterladen.py laufen lassen, dann das Ergebnis "
              "(akkordeon_quellen/) hierher legen.", file=sys.stderr)
        return 1

    print("=" * 64)
    print(f" Quellen:  {quellen}")
    print(f" Ausgabe:  {ausgabe}")
    print("=" * 64)

    # Bilder-Index laden
    idx_pfad = quellen / "bilder_index.json"
    if not idx_pfad.exists():
        print(f"FEHLER: bilder_index.json fehlt in {quellen}", file=sys.stderr)
        return 1
    bilder_index = json.loads(idx_pfad.read_text(encoding="utf-8"))

    # Ausgabeverzeichnis vorbereiten
    kapitel_dir = ausgabe / "kapitel"
    bilder_dir = ausgabe / "bilder"
    kapitel_dir.mkdir(parents=True, exist_ok=True)
    bilder_dir.mkdir(parents=True, exist_ok=True)

    # Bilder kopieren
    quell_bilder = quellen / "bilder"
    if quell_bilder.exists():
        n = 0
        for src in quell_bilder.iterdir():
            if src.is_file():
                ziel = bilder_dir / src.name
                if not ziel.exists() or ziel.stat().st_size != src.stat().st_size:
                    shutil.copy2(src, ziel)
                n += 1
        print(f"Bilder bereitgestellt: {n} in {bilder_dir}")

    # Kapitel bauen
    kapitel_pfade: List[Tuple[str, List[Path]]] = []
    fehlend: List[str] = []

    for teil_titel, artikel in KAPITEL:
        print(f"\n=== {teil_titel} ===")
        pfade: List[Path] = []
        for art in artikel:
            p = baue_kapitel(art, quellen, kapitel_dir, bilder_index)
            if p:
                pfade.append(p)
            else:
                fehlend.append(art)
        kapitel_pfade.append((teil_titel, pfade))

    # Master-Datei
    buch_tex = ausgabe / "buch.tex"
    schreibe_buch_tex(buch_tex, kapitel_pfade, ausgabe)
    print(f"\nMaster-Datei: {buch_tex}")

    if fehlend:
        print(f"\n{len(fehlend)} Artikel fehlten und wurden übersprungen:")
        for f in fehlend:
            print(f"  - {f}")

    # Kompilieren
    if args.kompilieren:
        ok = kompiliere(ausgabe)
        return 0 if ok else 1

    print("\nZum Kompilieren:")
    print(f"  cd {ausgabe}")
    print("  lualatex buch.tex   # zweimal aufrufen für Inhaltsverzeichnis")
    print("oder gleich:")
    print(f"  python3 {Path(__file__).name} --kompilieren")
    return 0


if __name__ == "__main__":
    sys.exit(main())
