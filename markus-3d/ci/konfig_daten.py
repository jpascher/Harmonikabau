# Daten fuer die Webseite docs/konfigurator.html erzeugen (Johann 26.11.2026):
# Reihen, Knopfzahlen und gezeichnete Stimmplattenlaengen je Instrument/Knopfzahl aus den Tastenlisten,
# dazu die erlaubten Kombinationen aus tastenliste_bauen.py. Nach Aenderungen an Tastenlisten neu ausfuehren:
#   python3 markus-3d/ci/konfig_daten.py docs/konfig_daten.json
import sys, os, json
HIER = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HIER)
import tastenliste_bauen as tb

def reihen(name):
    _, r = tb.teile(tb.lies(name))
    return {nr: {"knoepfe": int(tb.wert(blk, "knoepfe") or 0), "platten": [float(x) for x in (tb.wert(blk, "platten") or "").split()]}
            for nr, blk in r.items()}

d = {"diskant": {"markus": {"50": reihen("markus_std_10mm")}, "marie": {}},
     "regeln": {"abstaende": list(tb.ABSTAENDE), "geprueft": sorted(f"{i}:{g}:{a}" for i, g, a in tb.GEPRUEFT), "knoepfe": tb.KNOEPFE_ERLAUBT, "marie_griffstock": list(tb.MARIE),
                "platte_min": tb.PLATTE_MIN, "platte_max": tb.PLATTE_MAX}}
for n in (44, 46, 48, 52, 54): d["diskant"]["markus"][str(n)] = reihen(f"markus{n}_winkel")
for g, f in tb.MARIE.items(): d["diskant"]["marie"][g] = reihen(f)
d["diskant"]["marie"]["standard"] = reihen("marie")
json.dump(d, open(sys.argv[1] if len(sys.argv) > 1 else "konfig_daten.json", "w", encoding="utf-8"), ensure_ascii=False)
print("geschrieben:", sys.argv[1] if len(sys.argv) > 1 else "konfig_daten.json")
