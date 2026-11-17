import os, glob, math, re, cadquery as cq, ezdxf
from shapely.geometry import Polygon
from shapely.ops import unary_union
import schnitt_trafo as ST
B = "sec_brep/"
def lade(f): return cq.Shape.importBrep(B + f)
dateien = [os.path.basename(f) for f in glob.glob(B + "*.brep")]
statisch = [f for f in dateien if not f.startswith(("Hebel_", "Stimmplatte", "Verdeck"))]
hebel = [f for f in dateien if f.startswith("Hebel_")]
def zeich(Y, Z):                      # Modell (Y,Z) -> Zeichnung gschnitt.dxf
    u, s = -Z, ST.Y0 - Y
    return (ST.P_B[0] + s * ST.D[0] + u * ST.N[0], ST.P_B[1] + s * ST.D[1] + u * ST.N[1])
LAYER = {"Griffstock": 42, "Kunststoff": 8, "Dreieck": 34, "Diskantboden": 32, "Verstaerkung": 34, "Verdeck": 250,
         "Rahmen": 33, "Schraube": 9, "Stimmstock": 36, "Stimmzunge": 9}
ROT = {1: 1, 2: 5, 3: 3, 4: 30}
def layer_statisch(f):
    for k in LAYER:
        if f.startswith(k): return k
    return "Sonstiges"
CACHE = {}
def kanten(f, x):
    s = CACHE.setdefault(f, lade(f))
    bb = s.BoundingBox()
    if not (bb.xmin - 1e-6 <= x <= bb.xmax + 1e-6): return []
    try: sec = s.intersect(cq.Face.makePlane(800, 800, cq.Vector(x, 100, 0), cq.Vector(1, 0, 0)))
    except Exception: return []
    out = []
    for e in sec.Edges():
        if e.geomType() == "LINE": ts = [0, 1]
        else: ts = [i / max(2, int(e.Length() / 0.3)) for i in range(max(2, int(e.Length() / 0.3)) + 1)]
        out.append([e.positionAt(t) for t in ts])
    return out
def projektion(f):
    s = CACHE.setdefault(f, lade(f)); v, t = s.tessellate(0.03, 0.1)
    tris = [Polygon([(v[a].y, v[a].z), (v[b].y, v[b].z), (v[c].y, v[c].z)]) for a, b, c in t]
    tris = [p.buffer(0.001) for p in tris if p.area > 1e-9]
    return unary_union(tris).buffer(-0.001)
# je Reihe der mittlere Hebel
grp = {}
for f in hebel:
    g = f.split("__")[0]; grp.setdefault(int(g[7]), set()).add(g)
wahl = {}
for r, gs in grp.items():
    xs = sorted((lade(g + "__Knopf.brep").BoundingBox().center.x, g) for g in gs if os.path.exists(B + g + "__Knopf.brep"))
    wahl[r] = xs[len(xs) // 2]
print(wahl)
ges = ezdxf.new("R2010"); gm = ges.modelspace()
for r in sorted(wahl):
    x, g = wahl[r]
    doc = ezdxf.new("R2010"); ms = doc.modelspace()
    def lay(d, n, c):
        if n not in d.layers: d.layers.add(n, color=c)
    dx = 0.0
    def poly(n, c, pts, closed=False):
        n = f"Schnitt{r}_" + n
        for d, m, off in ((ges, gm, dx),):
            lay(d, n, c); m.add_lwpolyline([(p[0] + off, p[1]) for p in pts], close=closed, dxfattribs={"layer": n})
    # Schnitt aller festen Teile
    for f in statisch:
        n = layer_statisch(f)
        for k in kanten(f, x): poly(n + "_Schnitt", LAYER.get(n, 7), [zeich(p.y, p.z) for p in k])
    # Schnitt aller Hebel, die die Ebene kreuzen
    for f in hebel:
        rr = int(f[7]); teil = f.split("__")[1][:-5]
        if teil in ("Knopf", "Klappe", "Klappenpolster") or f.startswith(g + "__"): continue
        for k in kanten(f, x): poly(f"R{rr}_{teil}_Schnitt", ROT[rr], [zeich(p.y, p.z) for p in k])
    # Ansicht (Projektion) des gewaehlten Hebels: Knopf, Holzteil, Alu, Achse, Klappe, Polster
    for teil in ("Knopf", "Holzteil", "Hebel_Alu", "Achse", "Klappe", "Klappenpolster"):
        f = f"{g}__{teil}.brep"
        if not os.path.exists(B + f): continue
        P = projektion(f)
        for p in getattr(P, "geoms", [P]):
            for ring in [p.exterior, *p.interiors]:
                poly(f"Hebel_R{r}_Ansicht_{teil}", ROT[r] if teil in ("Hebel_Alu", "Achse") else (34 if teil == "Holzteil" else 7),
                     [zeich(a, b) for a, b in ring.coords], True)
    t0 = zeich(230, -90 - 6 * r)
    lay(ges, f"Schnitt{r}_Text", 7)
    gm.add_text(f"Schnitt {r}: X={x:.1f}, Hebel {g} als Ansicht", height=3.5,
                dxfattribs={"layer": f"Schnitt{r}_Text"}).set_placement(t0)
    print("Reihe", r, "fertig")
ges.saveas("Schnitt_Hebel_4_Layer.dxf")
