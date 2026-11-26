# Auswahl fuer "Harmonika erzeugen" einlesen und pruefen (Johann 26.11.2026)
# Quelle: Formular "Run workflow" (Umgebungsvariablen IN_<feld>) oder eine Modell-Anfrage als Issue (ISSUE_BODY).
# Nur gepruefte Werte aus festen Listen bzw. streng gepruefte Zahlen-/Kurzform-Texte werden weitergegeben.
# Ausgabe: Zeilen "name=wert" nach $GITHUB_OUTPUT; bei Fehlern ok=nein und grund=<Meldung>.
import os, re, sys, json, glob

HIER = os.path.dirname(os.path.abspath(__file__))
M3D = os.path.dirname(HIER)
LISTEN = {"teil": ["bass", "diskant"], "bass_variante": ["markus", "marie38", "eduard", "gerhard", "andreas", "v340"],
          "teilung": ["standard", "18", "19", "20"], "instrument": ["markus", "marie"],
          "griffstock": ["standard", "traditionell", "kugellager", "kolben", "kolben12", "winkel", "winkel_ital"],
          "knoepfe": ["50", "44", "46", "48", "52", "54"], "hebel": ["10", "9", "8"]}
STANDARD = {"teil": "bass", "bass_variante": "markus", "teilung": "standard", "laenge": "380", "instrument": "markus",
            "griffstock": "standard", "knoepfe": "50", "hebel": "10", "platten": "", "belegung": "", "vorschau": "ja"}
# Feldtitel im Issue-Formular (.github/ISSUE_TEMPLATE/modell_anfrage.yml) -> Feld
TITEL = {"Teil": "teil", "Bass: Bassteil": "bass_variante", "Bass: Knopfteilung": "teilung",
         "Bass: Gehäuselänge (nur Markus-Bass, mm)": "laenge", "Bass: Belegung (nur v340)": "belegung",
         "Diskant: Instrument": "instrument", "Diskant: Griffstock": "griffstock", "Diskant: Knopfzahl": "knoepfe",
         "Diskant: Hebelteilung (mm)": "hebel", "Diskant: Stimmplattenlängen": "platten", "Vorschau erzeugen": "vorschau"}

def aus_issue(text):
    werte, feld = {}, None
    for z in (text or "").splitlines():
        m = re.match(r"^###\s+(.+?)\s*$", z)
        if m: feld = TITEL.get(m.group(1).strip()); continue
        if feld and z.strip() and z.strip() != "_No response_":
            werte[feld] = (werte.get(feld, "") + " " + z.strip()).strip()
    return werte

def pruefe(w):
    a = dict(STANDARD); a.update({k: v.strip() for k, v in w.items() if v is not None and v.strip() != ""})
    for k, erlaubt in LISTEN.items():
        if a[k] not in erlaubt: raise ValueError(f"{k}: '{a[k]}' ist nicht erlaubt ({', '.join(erlaubt)})")
    a["vorschau"] = "true" if a["vorschau"].lower() in ("ja", "true", "1") else "false"
    try: l = float(a["laenge"].replace(",", "."))
    except ValueError: raise ValueError(f"Gehäuselänge '{a['laenge']}' ist keine Zahl")
    if not 300 <= l <= 500: raise ValueError(f"Gehäuselänge {l:g} mm außerhalb 300–500 mm")
    a["laenge"] = f"{l:g}"
    if not re.fullmatch(r"[Rr0-9:;,.\s]{0,600}", a["platten"]): raise ValueError("Stimmplattenlängen: nur R, Ziffern, Doppelpunkt, Strichpunkt, Komma, Punkt")
    if not re.fullmatch(r"[0-9.=A-Za-z,;\s]{0,800}", a["belegung"]): raise ValueError("Belegung: nur die Kurzform '1.1=V077, 1.2=BL152, ...'")
    a["platten"] = " ".join(a["platten"].split()); a["belegung"] = " ".join(a["belegung"].split())
    sys.path.insert(0, HIER); sys.path.insert(0, os.path.join(M3D, "bass"))
    if a["teil"] == "diskant":
        import tastenliste_bauen as tb
        tb.pruefe(a["instrument"], a["griffstock"], int(a["knoepfe"]), int(a["hebel"]))
        if a["platten"]: tb.lies_platten(a["platten"])
        a["variante"] = f'{a["instrument"]}:{a["griffstock"]}:{a["knoepfe"]}:{a["hebel"]}'
        a["titel"] = f'Diskant {a["instrument"]}, {a["griffstock"]}, {a["knoepfe"]} Knoepfe, Hebel {a["hebel"]} mm'
    else:
        a["variante"] = a["bass_variante"]
        if a["belegung"]:
            mj = glob.glob(os.path.join(M3D, "bass", a["bass_variante"], "belegung_modell_*.json"))
            if not mj: raise ValueError(f"Freie Belegung gibt es für {a['bass_variante']} noch nicht (bisher nur v340)")
            import belegung
            m, katalog, gassen, kx, reihen = belegung.lade_modell(mj[0])
            belegung.baue_tabelle(a["belegung"], katalog, gassen, kx, reihen, m.get("min_luecke", 12.0))
        a["titel"] = f'Bass {a["bass_variante"]}, Teilung {a["teilung"]}' + (", eigene Belegung" if a["belegung"] else "")
    return a

def main():
    out = open(os.environ.get("GITHUB_OUTPUT", "/dev/stdout"), "a", encoding="utf-8")
    if os.environ.get("QUELLE") == "issue": w = aus_issue(os.environ.get("ISSUE_BODY", ""))
    else: w = {k: os.environ.get("IN_" + k.upper(), "") for k in STANDARD}
    try:
        a = pruefe(w)
    except (ValueError, SystemExit) as e:
        grund = " ".join(str(e).split())
        print("Auswahl ungueltig:", grund)
        out.write("ok=nein\n"); out.write(f"grund={grund}\n"); return
    for k in ("teil", "variante", "teilung", "laenge", "vorschau", "platten", "belegung", "titel"):
        out.write(f"{k}={a[k]}\n")
    out.write("ok=ja\n")
    print("Auswahl:", json.dumps({k: a[k] for k in ("teil", "variante", "teilung", "laenge", "vorschau", "platten", "belegung")}, ensure_ascii=False))

if __name__ == "__main__":
    main()
