# Rahmen Bass Markus zusammenbauen: Front, Rueckteil, 2 Seitenteile (Zinken ineinander) + Boden
# X = Laenge (382), Y = Breite (207), Z = Hoehe (76); Z = 0 = Unterkante der Rahmenzeichnung
import sys, ezdxf, ezdxf.path, cadquery as cq
D = 5.5
L, B, H = 382.0, 207.0, 76.0
Z_BODEN, T_BODEN = 31.2, 6.0                 # Boden-Oberkante 31,2 mm ab Z=0 (aus Schnitt) - ANNAHME
OX, OY = 2924.3 - 1.0, 1521.7 - 1.0           # Draufsicht der Hauptzeichnung -> Rahmen (Aussenmass 380x205 mittig in 382x207)

imp = lambda f: cq.importers.importStep(f).val()
front, rueck, seite = imp("Markus_front_v4.step"), imp("Markus_rueck_v4.step"), imp("Markus_seite_v4.step")
X, Z = cq.Vector(1, 0, 0), cq.Vector(0, 0, 1); O = cq.Vector(0, 0, 0)
stehend = lambda s: s.rotate(O, X, 90)        # Platte aufstellen: Hoehe -> Z, Aussenseite -> -Y
rueck_s = stehend(rueck).translate(cq.Vector(0, D, 0))                        # Rueckteil bei Y = 0
front_s = stehend(front).rotate(O, Z, 180).translate(cq.Vector(L, B - D, 0))   # Front bei Y = B (auf Johanns Angabe getauscht)
seite_l = stehend(seite).rotate(O, Z, -90).translate(cq.Vector(D, B, 0))
seite_r = stehend(seite).rotate(O, Z, 90).translate(cq.Vector(L - D, 0, 0))

# Boden mit Oeffnungen aus der Hauptzeichnung
m = ezdxf.readfile(sys.argv[1]).modelspace()
boden = cq.Workplane("XY").box(L - 2*D, B - 2*D, T_BODEN, centered=False).translate((D, D, Z_BODEN - T_BODEN))
def loch(wp):
    global boden; boden = boden.cut(wp.extrude(T_BODEN + 2).translate((0, 0, Z_BODEN - T_BODEN - 1)))
n = 0
for lay in ("Füllung-DH-Öffnungen2", "Begleiterstimmstock-fuellung-öffnungen"):
    for c in m.query(f'CIRCLE[layer=="{lay}"]'):
        loch(cq.Workplane("XY").center(c.dxf.center.x - OX, c.dxf.center.y - OY).circle(c.dxf.radius)); n += 1
for lay in ("2Füllung-Öffnug", "Füllung-luftklappe"):
    for p in m.query(f'LWPOLYLINE[layer=="{lay}"]'):
        pts = [(v.x - OX, v.y - OY) for v in ezdxf.path.make_path(p).flattening(0.05)][:-1]
        if 0 < min(a for a, _ in pts) < L:
            loch(cq.Workplane("XY").polyline(pts).close()); n += 1
print("Bodenoeffnungen", n)

teile = [("Front", front_s, (0.80, 0.62, 0.40)), ("Rueckteil", rueck_s, (0.74, 0.56, 0.36)),
         ("Seitenteil_links", seite_l, (0.68, 0.50, 0.32)), ("Seitenteil_rechts", seite_r, (0.68, 0.50, 0.32)),
         ("Boden", boden.val(), (0.88, 0.80, 0.62))]
a = cq.Assembly(name="Rahmen_Bass_Markus")
for nme, s, col in teile:
    a.add(s, name=nme, color=cq.Color(*col))
a.save(sys.argv[2])
# Kontrolle: Ueberschneidungen der Rahmenteile (sollten ~0 sein, wenn die Zinken passen)
for i in range(4):
    for j in range(i + 1, 4):
        v = teile[i][1].intersect(teile[j][1]).Volume()
        if v > 1: print("Ueberschneidung", teile[i][0], teile[j][0], round(v, 1), "mm3")
for nme, s, _ in teile:
    bb = s.BoundingBox(); print(nme, [round(v, 1) for v in (bb.xmin, bb.xmax, bb.ymin, bb.ymax, bb.zmin, bb.zmax)])
    cq.exporters.export(cq.Workplane().add(s), f"stl_r_{nme}.stl")
