# ---------- Hoerner (Horn-schnitt3.dxf, Schalltrichter) in den 6 Hornbohrungen der Front ----------
# Profil: Achse = Linie y=0, a = axial (0..11,5), r = radial (6,3..16,5); Schaft Ø18 x 4 steckt in der Bohrung,
# die Stufe bei a = 4 liegt an der Aussenseite der Front, der Trichter (Ø33) zeigt nach aussen.
HD = ezdxf.readfile("horn.dxf").modelspace()
segs = []
for e in HD.query('LINE ARC'):
    if e.dxftype() == "LINE" and abs(e.dxf.start.y) < 1e-6 and abs(e.dxf.end.y) < 1e-6: continue   # Achse
    p = [(round(q.x, 4), round(q.y, 4)) for q in ezdxf.path.make_path(e).flattening(0.02)]
    segs.append(p)
kette = segs.pop(0)
while segs:
    for i, s in enumerate(segs):
        if math.dist(s[0], kette[-1]) < 0.05: kette += s[1:]; segs.pop(i); break
        if math.dist(s[-1], kette[-1]) < 0.05: kette += s[::-1][1:]; segs.pop(i); break
    else:
        raise SystemExit("Hornprofil nicht geschlossen")
horn_prof = [(a_, r_) for a_, r_ in kette[:-1]] if math.dist(kette[0], kette[-1]) < 0.05 else kette
A_STUFE = 4.01
def horn(x, z, y_aussen):
    # Profil in der XY-Ebene (X = radial, Y = axial), um die Y-Achse drehen
    w = cq.Workplane("XY").polyline([(r_, a_) for a_, r_ in horn_prof]).close().revolve(360, (0, 0, 0), (0, 1, 0))
    return w.val().translate(cq.Vector(x, y_aussen - A_STUFE, z))
hoerner = []
for c in RD.query('CIRCLE[layer=="Hörner-Bohrung"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        X_, Z_ = -x, y + 4.24
        ring = cq.Workplane("XZ").center(X_, Z_).circle(r + 2.5).circle(r).extrude(-40).translate((0, B - 30, 0)).val()
        y_out = front.intersect(ring).BoundingBox().ymax
        hoerner.append((f"Horn_{len(hoerner)+1}", horn(X_, Z_, y_out)))
print("Hoerner", len(hoerner), "Profilpunkte", len(horn_prof), "Ø max", round(2*max(r_ for _, r_ in horn_prof), 1))
# ---------- Fuesse: Unterlage Ø15,5 x 4,5 mit Bohrung Ø3,5 (knopf-fuer-fuesse.cb), unter dem Bassdeckel ----------
# je Eckblock einer, mittig unter der Deckelschraube; eine Schraube Ø3 durch Fuss und Deckel in den Eckblock
D_FUSS, H_FUSS, L_SCHR = 15.5, 4.5, 20.0
fuesse = []
for k, (bx, by) in enumerate(ecken):
    cx, cy = bx + BLK/2, by + BLK/2
    fu = cq.Workplane("XY").center(cx, cy).circle(D_FUSS/2).circle(D_DS/2).extrude(H_FUSS).translate((0, 0, -H_FUSS)).val()
    sc = cq.Workplane("XY").center(cx, cy).circle(D_BL/2).extrude(L_SCHR).translate((0, 0, -H_FUSS)).val() \
        .fuse(cq.Workplane("XY").center(cx, cy).circle(3.0).extrude(1.5).translate((0, 0, -H_FUSS - 1.5)).val())
    fuesse += [(f"Fuss_{k+1}", fu), (f"Fuss_Schraube_{k+1}", sc)]
# Kollisionskontrolle der neuen Teile gegen alles
alle = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile
for n1, s1 in hoerner + fuesse:
    for n2, s2 in alle:
        vv = s1.intersect(s2).Volume()
        if vv > 0.5: print("Ueberschneidung", n1, n2, round(vv, 1))
