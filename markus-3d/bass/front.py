# Bass-Front Markus aus Markus-Bass-Rahmen.dxf - liegend, X = Laenge, Y = Hoehe, Z = Dicke
import sys, ezdxf, ezdxf.path, cadquery as cq
m = ezdxf.readfile(sys.argv[1]).modelspace()
TEIL = sys.argv[3] if len(sys.argv) > 3 else 'front'   # 'front' oder 'rueck'
X0, X1, Y0, Y1 = -382.0, 0.0, -4.24, 71.75     # Front in der Zeichnung
DICKE = 5.5                                      # = Zinkentiefe der Seitenteile
ZT = 5.5                                         # Zinkentiefe
NUTEN = [(0, 6), (12, 18), (24, 30), (36, 42), (48, 54), (60, 66)]   # Aussparungen an beiden Enden (y)

def sx(x): return x - X0
def sy(y): return y - Y0
f = cq.Workplane("XY").rect(X1-X0, Y1-Y0, centered=False).extrude(DICKE)
def cut2d(wp):
    global f; f = f.cut(wp.extrude(DICKE+2).translate((0, 0, -1)))
for a, b in NUTEN:                                  # Zinken links und rechts
    cut2d(cq.Workplane("XY").center(0, sy(a)).rect(ZT, b-a, centered=False))
    cut2d(cq.Workplane("XY").center(sx(X1)-ZT, sy(a)).rect(ZT, b-a, centered=False))
n_orn = n_loch = 0
for p in m.query('LWPOLYLINE[layer=="Layer2"]'):   # Ornament-Durchbrueche an den Enden
    pts = [(v.x, v.y) for v in ezdxf.path.make_path(p).flattening(0.05)][:-1]   # mit Boegen
    lo = min(a for a, _ in pts); hi = max(a for a, _ in pts)
    if lo >= X0 and hi <= X1 and min(lo-X0, X1-hi) < 40:
        cut2d(cq.Workplane("XY").polyline([(sx(a), sy(b)) for a, b in pts]).close()); n_orn += 1
for lay in (("Hörner-Bohrung", "bohrungen-bassdrücker") if TEIL == "front" else ()):
    for c in m.query(f'CIRCLE[layer=="{lay}"]'):
        x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
        if X0 < x < X1 and Y0 < y < Y1:
            cut2d(cq.Workplane("XY").center(sx(x), sy(y)).circle(r)); n_loch += 1
# Verzierungen: 2 mm tiefe Taschen von der Aussenseite (Z = DICKE), nicht durchbrochen
TIEFE = 2.0
def tasche(pts):
    global f
    f = f.cut(cq.Workplane("XY").polyline([(sx(a), sy(b)) for a, b in pts]).close()
              .extrude(TIEFE+1).translate((0, 0, DICKE-TIEFE)))
n_t = 0
for e in m.query('SPLINE[layer=="LAYER_1"]'):                       # Ranken = Durchbrueche im Rueckteil
    if TEIL != "rueck": continue
    pts = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.05)][:-1]
    if min(b for _, b in pts) > Y0:
        cut2d(cq.Workplane("XY").polyline([(sx(a), sy(b)) for a, b in pts]).close()); n_t += 1
for e in m.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]'):  # wellige Flaeche unten an der Front, 2 mm tief
    if TEIL == "front" and not e.closed:
        pts = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.05)]
        tasche(pts + [(pts[-1][0], Y0-1), (pts[0][0], Y0-1)]); n_t += 1
print(TEIL, "Ranken/Taschen", n_t)
cq.exporters.export(f, sys.argv[2])
print("Ornamente", n_orn, "Bohrungen", n_loch, "Volumen", round(f.val().Volume()))
