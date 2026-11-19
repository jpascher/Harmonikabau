# Schritt 1 Tastenliste: liest aus den Markus-Eingaben (boden.dxf, gs.dxf, var/BESASDES) je Taste Knopf-X, Klappe, Bodenloecher
# und Kanzellen aus und schreibt tastenliste_markus.json + eine Teilungs-Auswertung.
import ezdxf, json, math, numpy as np
from ezdxf import bbox
from collections import defaultdict
BXO = 380.0
def reihe_y(y): return 1 if y > 152 else 2 if y > 117 else 3 if y > 82 else 4
bm = ezdxf.readfile("boden.dxf").modelspace()
klappen = []
seen = set()
for ins in bm.query('INSERT[layer=="2Füllung-Öffnug"]'):
    p0 = (round(ins.dxf.insert.x, 2), round(ins.dxf.insert.y, 2))
    if p0 in seen: continue
    seen.add(p0)
    cs = [(v.dxf.center.x + BXO, v.dxf.center.y) for v in ins.virtual_entities() if v.dxftype() == "CIRCLE"]
    m = (sum(c[0] for c in cs) / len(cs), sum(c[1] for c in cs) / len(cs))
    klappen.append({"reihe": reihe_y(m[1]), "mitte": [round(m[0], 2), round(m[1], 2)], "loecher": [[round(a, 2), round(b, 2)] for a, b in sorted(cs)]})
gm = ezdxf.readfile("gs.dxf").modelspace()
lr = defaultdict(list)
for e in gm:
    if e.dxf.layer != "Lagerbock_aufsicht": continue
    b = bbox.extents([e])
    if b.has_data and abs(b.size.x - 5.0) < 0.2 and (abs(b.size.y - 14.2) < 0.3 or abs(b.size.y - 11.3) < 0.3):
        lr[round(b.center.y, 0)].append(round(b.center.x + BXO, 2))
ly = sorted(lr, reverse=True)
knoepfe = {i + 1: sorted(lr[y]) for i, y in enumerate(ly)}
knopf_y = {i + 1: y for i, y in enumerate(ly)}
out = {"knopfreihen": {r: {"y_aufsicht": knopf_y[r], "x": knoepfe[r]} for r in knoepfe},
       "klappen": sorted(klappen, key=lambda k: (k["reihe"], k["mitte"][0]))}
json.dump(out, open("tastenliste_markus_roh.json", "w"), indent=1)
print("Knopfreihen (Aufsicht-Y, Anzahl, Teilung):")
for r in sorted(knoepfe):
    xs = knoepfe[r]; print("  R%d y=%s n=%d dX=%s  X %.1f..%.1f" % (r, knopf_y[r], len(xs), sorted(set(np.round(np.diff(xs), 2))), xs[0], xs[-1]))
print("Klappen je Stimmstockreihe:")
for r in (1, 2, 3, 4):
    ks = sorted([k for k in klappen if k["reihe"] == r], key=lambda k: k["mitte"][0])
    xs = [k["mitte"][0] for k in ks]; ys = sorted(set(round(k["mitte"][1], 1) for k in ks))
    lo = ks[0]["loecher"]
    print("  R%d n=%d dX=%s  Y=%s  Loecher je Klappe=%d, Lochabstand %s" % (r, len(ks), sorted(set(np.round(np.diff(xs), 2))), ys, len(lo),
          [round(math.dist(lo[i], lo[i+1]), 2) for i in range(len(lo)-1)]))
print("Abstand Knopfreihen (Aufsicht):", [round(knopf_y[r] - knopf_y[r+1], 2) for r in (1, 2, 3)])
