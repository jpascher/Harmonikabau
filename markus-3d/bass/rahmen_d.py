# Diskant-Rahmenteile Markus (wie Bass, Hoehe H, Verzierungen nur als Taschen, ohne Durchbrueche/Bohrungen)
# Rahmenteile Bass Markus aus Markus-Bass-Rahmen.dxf - liegend: X = Laenge, Y = Hoehe, Z = Dicke (Z = DICKE = Aussenseite)
# Aufruf: python3 teile.py Markus-Bass-Rahmen.dxf <ausgabe.step> front|rueck|seite
import sys, ezdxf, ezdxf.path, cadquery as cq
m = ezdxf.readfile(sys.argv[1]).modelspace()
TEIL = sys.argv[3]
DICKE, ZT, TIEFE = 5.5, 5.5, 2.0          # Dicke, Zinkentiefe, Tiefe der Verzierungen
H = float(sys.argv[4]) if len(sys.argv) > 4 else 70.0
Y0, Y1 = -4.24, -4.24 + H
if TEIL in ("front", "rueck"):
    X0, X1 = -382.0, 0.0
    NUTEN = [(0, 6), (12, 18), (24, 30), (36, 42), (48, 54), (60, 66)]
    ORN_LO = (-374.9, -28.85)               # Ornamente an den Enden (linke Kante in der Zeichnung)
else:
    X0, X1 = -207.0, 0.0
    NUTEN = [(-4.24, -0.25), (6.25, 11.75), (18.25, 23.75), (30.25, 35.75), (42.25, 47.75), (54.25, 59.75), (66.25, 71.75)]
    ORN_LO = (-199.9, -28.52)

def sx(x): return x - X0
def sy(y): return y - Y0
def kontur(e): return [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.05)]
f = cq.Workplane("XY").rect(X1-X0, Y1-Y0, centered=False).extrude(DICKE)
def durch(pts):
    global f; f = f.cut(cq.Workplane("XY").polyline([(sx(a), sy(b)) for a, b in pts]).close().extrude(DICKE+2).translate((0, 0, -1)))
def tasche(pts):
    global f; f = f.cut(cq.Workplane("XY").polyline([(sx(a), sy(b)) for a, b in pts]).close().extrude(TIEFE+1).translate((0, 0, DICKE-TIEFE)))
def kreis(x, y, r):
    global f; f = f.cut(cq.Workplane("XY").center(sx(x), sy(y)).circle(r).extrude(DICKE+2).translate((0, 0, -1)))

NUTEN = [(a, min(b, Y1)) for a, b in NUTEN if a < Y1 - 0.5]
for a, b in NUTEN:                                         # Zinken an beiden Enden
    for xa in (X0, X1-ZT):
        durch([(xa, a), (xa+ZT, a), (xa+ZT, b), (xa, b)])
n = {"Ornament-Tasche": 0, "Ranken-Tasche": 0, "Wellen-Tasche": 0, "Bohrung": 0}
# Ornamente an den Enden: 2 mm tief; je Ende senkrecht zentriert und bei niedrigerem Rahmen verkleinert (4 mm Rand)
_gr_orn = {}
for p in m.query('LWPOLYLINE[layer=="Layer2"]'):
    pts = kontur(p)[:-1]
    for lo in ORN_LO:
        if abs(min(a for a, _ in pts) - lo) < 0.1: _gr_orn.setdefault(lo, []).append(pts)
for lo, gl in _gr_orn.items():
    xa = min(a for q in gl for a, _ in q); xb = max(a for q in gl for a, _ in q)
    ya = min(b for q in gl for _, b in q); yb = max(b for q in gl for _, b in q)
    s = min(1.0, (H - 8.0) / (yb - ya))
    cx, cy = (xa + xb) / 2, (ya + yb) / 2
    xr = xa if lo == min(ORN_LO) else xb                  # Abstand zur Rahmenkante bleibt
    for pts in gl:
        tasche([(xr + (a - xr) * s, Y0 + H / 2 + (b - cy) * s) for a, b in pts]); n["Ornament-Tasche"] += 1
    print("Ornament", lo, "Massstab %.3f" % s, "Hoehe %.1f" % ((yb - ya) * s))
if False:                                                   # Ranken (nur Bassteil) - beim Diskant entfallen sie
    for e in m.query('SPLINE[layer=="LAYER_1"]'):
        pts = kontur(e)[:-1]
        if min(b for _, b in pts) > Y0: tasche(pts); n["Ranken-Tasche"] += 1
if False:                                                   # wellige Flaeche (nur Bassteil) - beim Diskant entfaellt sie
    for e in m.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]'):
        if not e.closed:
            pts = kontur(e); tasche(pts + [(pts[-1][0], Y0-1), (pts[0][0], Y0-1)]); n["Wellen-Tasche"] += 1
    for lay in ():
        for c in m.query(f'CIRCLE[layer=="{lay}"]'):
            x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
            if X0 < x < X1 and Y0 < y < Y1: kreis(x, y, r); n["Bohrung"] += 1
cq.exporters.export(f, sys.argv[2])
print(TEIL, n, "Groesse", round(X1-X0, 1), "x", round(Y1-Y0, 1), "x", DICKE)
