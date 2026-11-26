# Tastenliste aus der Auswahl auf der Webseite zusammensetzen (Johann 26.11.2026)
#   instrument : markus (3-choerig) | marie (2-choerig)
#   griffstock : standard | traditionell | kugellager | kolben | kolben12 | winkel | winkel_ital
#   knoepfe    : 44 46 48 50 52 54 (nur Markus; die Reihen kommen aus tastenliste_markus<N>_winkel.txt)
#   hebel      : 8 | 9 | 10 mm Hebelteilung (Knopfabstand 16 / 18 / 20 mm; nur Markus, Winkel nur 10 mm)
#   platten    : optional eigene Stimmplattenlaengen je Reihe, z. B. "R1: 32 33 34 ...; R3: 27 29 ..."
#                (genau so viele Werte wie die Reihe Klappen hat; nicht genannte Reihen bleiben wie gezeichnet)
# Unzulaessige Kombinationen werden sofort mit einer klaren Meldung abgelehnt (pruefe()).
# Grundlage sind die vorhandenen, geprueften Tastenlisten in diskant/disk; neu kombiniert werden nur Kopf (Griffstock,
# Teilung, Knopfversatz) und Reihen (Knopfzahl, Klappen, Platten). Die Knopfpyramide richtet instrument.py selbst mittig aus.
import os, re

DISK = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "diskant", "disk")
HEBEL_ERLAUBT = {"standard": (8, 9, 10), "traditionell": (8, 9, 10), "kugellager": (8, 9, 10), "kolben": (8, 9, 10),
                 "kolben12": (8, 9, 10), "winkel": (10,), "winkel_ital": (10,)}
KNOEPFE_ERLAUBT = {"winkel": (44, 46, 48, 50, 52, 54), "winkel_ital": (44, 46, 48, 50, 52, 54),
                   "standard": (44, 46, 48, 50, 52, 54)}   # sonst nur 50: Kugellager scheitert an den POM-Stuetzen (26.11.)
PLATTE_MIN, PLATTE_MAX = 15.0, 100.0
MARIE = {"standard": "marie", "traditionell": "marie_trad", "kugellager": "marie_kugel", "kolben": "marie_kolben",
         "kolben12": "marie_kolben12", "winkel": "marie_winkel"}

def lies(name):
    return open(os.path.join(DISK, f"tastenliste_{name}.txt"), encoding="utf-8").read().splitlines()

def teile(zeilen):
    """Kopfzeilen und Reihenbloecke {nr: [zeilen]}"""
    kopf, reihen, akt = [], {}, None
    for z in zeilen:
        m = re.match(r"\s*\[reihe\s+(\d+)\]", z)
        if m: akt = int(m.group(1)); reihen[akt] = [z]; continue
        (reihen[akt] if akt else kopf).append(z)
    return kopf, reihen

def setze(zeilen, schluessel, wert, oben=False):
    """Wert setzen (oder entfernen mit wert=None); neu eingefuegt nach 'choerig'/'knopfabstand' bzw. (oben) nach der Blockzeile"""
    rx = re.compile(rf"^\s*{schluessel}\s*=")
    neu = [z for z in zeilen if not rx.match(z)]
    if wert is None: return neu
    pos = 1 if oben else next((i + 1 for i, z in enumerate(neu) if re.match(r"^\s*(choerig|knopfabstand)\s*=", z)), len(neu))
    neu.insert(pos, f"{schluessel:<14} = {wert}")
    return neu

def wert(zeilen, schluessel):
    for z in zeilen:
        m = re.match(rf"^\s*{schluessel}\s*=\s*([^#]*)", z)
        if m: return m.group(1).strip()
    return None

def lies_platten(text):
    """'R1: 32 33; R3: 27 29' -> {1: [32.0, 33.0], 3: [27.0, 29.0]}"""
    erg = {}
    for teil in re.split(r"[;\n]", text or ""):
        teil = teil.strip()
        if not teil: continue
        m = re.match(r"^[Rr](?:eihe)?\s*(\d)\s*[:=]\s*(.+)$", teil)
        if not m: raise SystemExit(f"Plattenlaengen: '{teil}' nicht verstanden - Form 'R1: 32 33 34 ...; R2: ...'")
        try: werte = [float(v.replace(",", ".")) for v in m.group(2).split()]
        except ValueError: raise SystemExit(f"Plattenlaengen Reihe {m.group(1)}: nur Zahlen in mm, getrennt durch Leerzeichen")
        erg[int(m.group(1))] = werte
    return erg

def pruefe(instrument, griffstock, knoepfe, hebel):
    """unzulaessige Auswahl sofort ablehnen"""
    if instrument not in ("markus", "marie"): raise SystemExit("Instrument muss markus oder marie sein")
    if griffstock not in HEBEL_ERLAUBT: raise SystemExit(f"unbekannter Griffstock {griffstock}")
    if instrument == "marie":
        if griffstock not in MARIE: raise SystemExit(f"Marie gibt es nicht mit Griffstock {griffstock} (moeglich: {', '.join(MARIE)})")
        return
    if hebel not in HEBEL_ERLAUBT[griffstock]:
        raise SystemExit(f"Griffstock {griffstock}: Hebelteilung {hebel} mm geht nicht - moeglich: {', '.join(map(str, HEBEL_ERLAUBT[griffstock]))} mm "
                         f"(Knopfabstand {', '.join(str(2 * h) for h in HEBEL_ERLAUBT[griffstock])} mm)")
    erlaubt = KNOEPFE_ERLAUBT.get(griffstock, (50,))
    if knoepfe not in erlaubt:
        raise SystemExit(f"Griffstock {griffstock}: {knoepfe} Knoepfe gehen derzeit nicht - moeglich: {', '.join(map(str, erlaubt))}. "
                         "Andere Knopfzahlen gibt es bisher nur mit Winkel- und Standard-Griffstock.")

def baue(instrument, griffstock, knoepfe, hebel, ziel, platten=""):
    """schreibt die Tastenliste nach ziel und gibt (pfad, name, hinweise) zurueck"""
    hinweise = []
    knoepfe, hebel = int(knoepfe), int(hebel)
    pruefe(instrument, griffstock, knoepfe, hebel)
    eigene = lies_platten(platten)
    if instrument == "marie":
        if griffstock not in MARIE: raise SystemExit(f"Marie gibt es nicht mit Griffstock {griffstock}")
        if knoepfe != 50: hinweise.append("Knopfzahl gilt nur fuer Markus; Marie mit ihrer gezeichneten Knopfzahl gerechnet.")
        if hebel != 10: hinweise.append("Hebelteilung gilt nur fuer Markus; Marie mit 19 mm Knopfabstand gerechnet.")
        zeilen = lies(MARIE[griffstock]); name = f"Marie_{griffstock}"
        zeilen = setze(zeilen, "name", name)
    elif instrument == "markus":
        if griffstock.startswith("winkel"):
            quelle = "markus_winkel_ital" if griffstock == "winkel_ital" else "markus_winkel"
        else:
            basis = {"standard": "std", "traditionell": "trad", "kugellager": "kugel", "kolben": "std", "kolben12": "std"}[griffstock]
            quelle = f"markus_{basis}_{hebel}mm"
        kopf, reihen = teile(lies(quelle))
        if griffstock in ("kolben", "kolben12"):
            kopf = setze(kopf, "griffstock", griffstock)
            hinweise.append(f"Markus mit {griffstock}-Griffstock ist eine neue Kombination (bisher nur Marie gerechnet) - Protokoll pruefen.")
        if knoepfe != 50:
            _, reihen_n = teile(lies(f"markus{knoepfe}_winkel"))
            for nr, blk in reihen_n.items():                  # Knopfversatz der Griffstockvariante uebernehmen
                kv = wert(reihen.get(nr, []), "knopf_versatz")
                reihen_n[nr] = setze(blk, "knopf_versatz", kv, oben=True) if kv else setze(blk, "knopf_versatz", None)
            reihen = reihen_n
            if griffstock == "standard": hinweise.append(f"{knoepfe} Knoepfe mit Standard-Griffstock: rechnet (Probe 46 Knoepfe am 26.11.), aber noch nicht angesehen - Protokoll pruefen.")
            kopf = setze(kopf, "gehaeuse_x", "auto")
            if not eigene: hinweise.append("Klappen der Knopfvariante mit Platzhalter-Plattenlaengen - eigene Laengen im Feld 'Plattenlaengen' eintragen.")
        name = f"Markus{knoepfe}_{griffstock}_{hebel}mm"
        kopf = setze(kopf, "name", name)
        zeilen = kopf + [z for nr in sorted(reihen) for z in reihen[nr]]
    else:
        raise SystemExit("instrument muss markus oder marie sein")
    if eigene:                                                # eigene Plattenlaengen einsetzen
        kopf_, reihen_ = teile(zeilen)
        for nr, werte in eigene.items():
            if nr not in reihen_: raise SystemExit(f"Plattenlaengen: Reihe {nr} gibt es bei {name} nicht (Reihen {', '.join(map(str, sorted(reihen_)))})")
            alt = (wert(reihen_[nr], "platten") or "").split()
            if len(werte) != len(alt):
                raise SystemExit(f"Plattenlaengen Reihe {nr}: {len(werte)} Werte angegeben, die Reihe hat aber {len(alt)} Klappen "
                                 f"(gezeichnet: {' '.join(alt)})")
            schlecht = [v for v in werte if not PLATTE_MIN <= v <= PLATTE_MAX]
            if schlecht: raise SystemExit(f"Plattenlaengen Reihe {nr}: {schlecht} ausserhalb {PLATTE_MIN:g}-{PLATTE_MAX:g} mm")
            reihen_[nr] = [(f"platten   = {' '.join(f'{v:g}' for v in werte)}   # eigene Eingabe" if re.match(r"^\s*platten\s*=", z) else z)
                           for z in reihen_[nr]]
            hinweise.append(f"Reihe {nr}: eigene Plattenlaengen {' '.join(f'{v:g}' for v in werte)}")
        zeilen = kopf_ + [z for nr in sorted(reihen_) for z in reihen_[nr]]
    os.makedirs(os.path.dirname(os.path.abspath(ziel)), exist_ok=True)
    open(ziel, "w", encoding="utf-8").write("\n".join(zeilen) + "\n")
    return ziel, name, hinweise

if __name__ == "__main__":
    import sys
    p, n, h = baue(*sys.argv[1:5], sys.argv[5], sys.argv[6] if len(sys.argv) > 6 else "")
    print(p, n); [print("Hinweis:", x) for x in h]
