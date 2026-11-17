# Lagerboecke (Buchsen) einfach + doppelt aus lager-bohrung-beide.dxf / lager-2mm.cb
# Querschnitt (linke Ansicht) wird ueber die Laenge 8,6 extrudiert (Achsrichtung = X).
# Achsbohrung(en) Ø2,4 in Laengsrichtung; Querbohrung Ø1,5 durchgehend (Hoehe 5 mm) in Laengsmitte = Befestigung,
# Schraube M1,4 x 10 in den Boden. Grundflaeche = 8,6 x Breite, Hoehe 5.
import sys, ezdxf, cadquery as cq
m = ezdxf.readfile(sys.argv[1]).modelspace()
LANG = 8.6
def bock(layer_profil, n_achsen):
    p = [e for e in m.query(f'LWPOLYLINE[layer=="{layer_profil}"]') if e.closed][0]
    # Zeichnung: x = -13.6..-8.6 -> Hoehe z = x + 13.6 (0..5), y -> Breite (Y = -y)
    pts = [(-q[1], q[0] + 13.6) for q in p.get_points("xy")]
    s = cq.Workplane("YZ").polyline(pts).close().extrude(LANG).val()
    ys = [2.6, 5.8][:n_achsen]
    for yy in ys:
        s = s.cut(cq.Workplane("YZ").center(yy, -11.1 + 13.6).circle(1.2).extrude(LANG + 2).translate((-1, 0, 0)).val())
        s = s.cut(cq.Workplane("XY").center(LANG / 2, yy).circle(0.75).extrude(7).translate((0, 0, -1)).val())
    schrauben = []
    for yy in ys:
        sc = cq.Workplane("XY").center(LANG / 2, yy).circle(0.7).extrude(10).translate((0, 0, 5 - 10)).val() \
            .fuse(cq.Workplane("XY").center(LANG / 2, yy).circle(1.25).extrude(0.9).translate((0, 0, 5)).val())
        schrauben.append(sc)
    return s, schrauben
einf, s1 = bock("Layer1", 1)
dopp, s2 = bock("0", 2)
dopp = dopp.translate(cq.Vector(15, 0, 0)); s2 = [x.translate(cq.Vector(15, 0, 0)) for x in s2]
a = cq.Assembly(name="Lagerboecke")
a.add(einf, name="Lagerbock_einfach", color=cq.Color(0.85, 0.85, 0.80))
for i, x in enumerate(s1): a.add(x, name=f"Schraube_M1_4x10_e{i+1}", color=cq.Color(0.5, 0.5, 0.52))
a.add(dopp, name="Lagerbock_doppelt", color=cq.Color(0.85, 0.85, 0.80))
for i, x in enumerate(s2): a.add(x, name=f"Schraube_M1_4x10_d{i+1}", color=cq.Color(0.5, 0.5, 0.52))
a.save(sys.argv[2])
for n, s in [("einfach", einf), ("doppelt", dopp)]:
    b = s.BoundingBox(); print(n, round(b.xlen, 2), "x", round(b.ylen, 2), "x", round(b.zlen, 2), "Vol", round(s.Volume(), 1))
