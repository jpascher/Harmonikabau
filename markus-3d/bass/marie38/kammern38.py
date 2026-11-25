# Bassteil Marie 38 cm (Zeichnung bassteil-alles.dxf / HB-trennwaende.dxf), aufgebaut auf dem Markus-Bassteil:
# Rahmen mit Zinken, Deckel mit Rosetten, hinterer Leiste, Leder, Eckbloecken und Fuessen, Deckelauflage, Gaze von Markus.
# Zeichnung 380 x 205 -> Rohmass 382 x 207 (Markus): Zeichnungskoordinaten + 1 mm (OX, OY).
# Achsen: X = Laenge, Y = Tiefe (0 = Rueckteil aussen, B = Front aussen), Z ab Deckel aussen Richtung Balg.
# Aufruf: python3 kammern38.py ausgabe.step
import sys, math, ezdxf, ezdxf.path, cadquery as cq
from ezdxf import bbox
OUT = sys.argv[1]
T_KNOPF = float(sys.argv[2]) if len(sys.argv) > 2 else 19.0   # Knopfteilung (Variable, Zeichnung 19)
D, L, B = 5.5, 382.0, 207.0
OX = OY = 1.0
T_DECKEL_BASS = 4.0
Z_BU, Z_BO = 41.0, 46.0                       # Boden Unterseite/Oberseite (Schnitt)
FI = B - D
def box(x0, x1, y0, y1, z0, z1):
    return cq.Solid.makeBox(x1 - x0, y1 - y0, z1 - z0, cq.Vector(x0, y0, z0))
def py(x0, x1, pts):                          # Profil (Y, Z) ueber X extrudiert
    return cq.Workplane("YZ").polyline(pts).close().extrude(x1 - x0).val().translate(cq.Vector(x0, 0, 0))
teile = []
# ---------- Rahmen mit Zinken (Markus v4); Front: Markus-Front ausgefuellt, dann Marie-Bohrungen ----------
imp = lambda f: cq.importers.importStep(f).val()
O, XA, ZA = cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), cq.Vector(0, 0, 1)
st = lambda s: s.rotate(O, XA, 90)
rueck = st(imp("Markus_rueck_v4.step")).translate(cq.Vector(0, D, 0))
front = st(imp("Markus_front_v4.step")).rotate(O, ZA, 180).translate(cq.Vector(L, B - D, 0))
seite = imp("Markus_seite_v4.step")
sl = st(seite).rotate(O, ZA, -90).translate(cq.Vector(D, B, 0))
sr = st(seite).rotate(O, ZA, 90).translate(cq.Vector(L - D, 0, 0))
H_R = front.BoundingBox().zmax
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v_ = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v_.Length, p0, v_.normalized())
A = ezdxf.readfile("bassteil-alles.dxf").modelspace()
X0, Y0P, ZF0 = -386.6, -624.0, -401.4         # Zeichnungs-Ursprung Draufsicht / Unterkante Frontansicht
# ---------- Front wie Markus: Knopfreihen-Abstufung (Reihe 1 vertieft, Aufdoppelung 2 + 5 mm), Hoerner, Blende, Knoepfe ----------
RD = ezdxf.readfile("rahmen.dxf").modelspace()
# Knopfreihen wie Standardmodell: Reihe 1 (unten, 9 Knoepfe, Bohrung Ø11 mit Filz) buendig an der Front, Reihe 2 (10 Knoepfe,
# Ø10) ueber eine 7-mm-Blende erhoeht; Hoehen und Teilung nach der Marie-Zeichnung (Teilung variabel), Mitte X 198.
X_KM = L / 2                                   # Knopfreihen mittig am Korpus
REIHEN_N = {1: 9, 2: 10}
Z_RH = {1: 401.4 - 387.5, 2: 401.4 - 373.5}                  # 13,9 / 27,9 ueber der Unterkante
R_RH = {1: 5.5, 2: 5.0}
kn = [(round(X_KM + (i_ - (n_ - 1) / 2) * T_KNOPF, 2), Z_RH[r_], R_RH[r_], r_) for r_, n_ in REIHEN_N.items() for i_ in range(n_)]
_wl = [e for e in RD.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]') if not e.closed][0]
_wp = [(q.x, q.y) for q in ezdxf.path.make_path(_wl).flattening(0.05)]
_wp = _wp + [(_wp[-1][0], -4.24 - 1), (_wp[0][0], -4.24 - 1)]
front = front.fuse(cq.Workplane("XZ").polyline([(-x, y + 4.24) for x, y in _wp]).close().extrude(-(D - 0.04)).translate((0, B - D + 0.02, 0)).val().intersect(box(0, L, B - D, B, 0, H_R)))   # Vertiefung Reihe 1 (Welle) zu
for c in RD.query('CIRCLE[layer=="bohrungen-bassdrücker"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        front = front.fuse(stab((-x, B - D + 0.02, y + 4.24), (-x, B - 0.02, y + 4.24), 2 * r + 0.2))
front = front.clean()
knopf_bohr = []
for X_, Z_, r_, rr in kn:
    cyl = stab((X_, FI - 15, Z_), (X_, B + 1, Z_), 2 * r_)
    front = front.cut(cyl); knopf_bohr.append(cyl)
wx0, wx1 = min(k[0] for k in kn) - 12.0, max(k[0] for k in kn) + 12.0
VERL_AUF = 0.0
HD = ezdxf.readfile("horn.dxf").modelspace()
segs = []
for e in HD.query('LINE ARC'):
    if e.dxftype() == "LINE" and abs(e.dxf.start.y) < 1e-6 and abs(e.dxf.end.y) < 1e-6: continue
    segs.append([(round(q.x, 4), round(q.y, 4)) for q in ezdxf.path.make_path(e).flattening(0.02)])
kette = segs.pop(0)
while segs:
    for i_, s_ in enumerate(segs):
        if math.dist(s_[0], kette[-1]) < 0.05: kette += s_[1:]; segs.pop(i_); break
        if math.dist(s_[-1], kette[-1]) < 0.05: kette += s_[::-1][1:]; segs.pop(i_); break
    else: raise SystemExit("Hornprofil nicht geschlossen")
horn_prof = [(a_, r_) for a_, r_ in kette[:-1]] if math.dist(kette[0], kette[-1]) < 0.05 else kette
aus_front = []
horn_pos = sorted({(round(-c.dxf.center.x, 2), round(c.dxf.center.y + 4.24, 2)) for c in RD.query('CIRCLE[layer=="Hörner-Bohrung"]')
                   if -382 < c.dxf.center.x < 0 and -4.24 < c.dxf.center.y < 71.75})     # wie Markus
for X_, Z_ in horn_pos:
        r = 9.0
        ring = cq.Workplane("XZ").center(X_, Z_).circle(r + 2.5).circle(r).extrude(-40).translate((0, B - 30, 0)).val()
        y_out = front.intersect(ring).BoundingBox().ymax
        w = cq.Workplane("XY").polyline([(r_, a_) for a_, r_ in horn_prof]).close().revolve(360, (0, 0, 0), (0, 1, 0))
        teile.append((f"Horn_{len(aus_front)+1}", w.val().translate(cq.Vector(X_, y_out - 4.01, Z_))))
        aus_front.append(cq.Workplane("XZ").center(X_, Z_).circle(r + 5.0).extrude(-(6.0 + 2)).translate((0, FI - 6.0 - 1, 0)).val())
# Blende ueber Reihe 2 (wie Standardmodell: Bogen um jede Bohrung, Einzug dazwischen, grosse Endboegen), 7 mm
D_ZYL, D_KOPF, H_KOPF, L_ZYL, SPALT, T_BL = 7.5, 15.5, 4.5, 20.0, 6.5, 7.0
def pts_of(e, tol=0.1):
    p_ = [(v_.x, v_.y) for v_ in ezdxf.path.make_path(e).flattening(tol)]
    q_ = [p_[0]] + [b_ for a_, b_ in zip(p_[:-1], p_[1:]) if math.dist(a_, b_) > 1e-6]
    if len(q_) > 2 and math.dist(q_[0], q_[-1]) < 1e-6: q_ = q_[:-1]
    return q_
def blende_kontur(xs_):
    xs_ = sorted(xs_, reverse=True); o_, u_ = 6.66, -4.42; pts = []
    pts.append((xs_[0] - 5.06, o_, -0.301))
    for h_ in xs_[1:-1]: pts += [(h_ + 5.06, o_, 0.388), (h_ - 5.06, o_, -0.301)]
    pts.append((xs_[-1] + 5.06, o_, 2.267)); pts.append((xs_[-1] + 5.06, u_, -0.301))
    for h_ in xs_[-2:0:-1]: pts += [(h_ - 5.06, u_, 0.388), (h_ + 5.06, u_, -0.301)]
    pts.append((xs_[0] - 5.06, u_, 2.267))
    e = ezdxf.new().modelspace().add_lwpolyline([(x_, y_, 0, 0, b_) for x_, y_, b_ in pts], format="xyseb", close=True)
    return pts_of(e, 0.05)
x2 = [k[0] for k in kn if k[3] == 2]
bl = cq.Workplane("XZ").polyline([(-x_, Z_RH[2] + y_) for x_, y_ in blende_kontur([-x_ for x_ in x2])]).close().extrude(-T_BL).val()
bl = bl.translate(cq.Vector(0, B - bl.BoundingBox().ymin, 0))
for X_, Z_, r_, rr in kn:
    if rr == 2: bl = bl.cut(stab((X_, B - 1, Z_), (X_, B + T_BL + 1, Z_), 2 * r_))
teile.append(("Frontblende", bl))
# Knoepfe in Ruhelage (wie Standardmodell), Filzeinlage in der aeussersten Lage (Front bzw. Blende)
for X_, Z_, r_, rr in kn:
    ya = B + (T_BL if rr == 2 else 0.0); dl_ = T_BL if rr == 2 else D; yk = ya + SPALT; lz = yk - (B + SPALT - L_ZYL)
    if rr == 1 and 207.5 + OX - D_ZYL / 2 - 1 < X_ < 285.5 + OX + D_ZYL / 2 + 1: lz -= 3.5   # Reihe 1 ueber der Lufttaste kuerzer (gedrueckt frei)
    teile.append((f"Knopf_R{rr}_{X_:.0f}", stab((X_, yk, Z_), (X_, yk + H_KOPF, Z_), D_KOPF)))
    teile.append((f"Knopfzylinder_R{rr}_{X_:.0f}", stab((X_, yk - lz, Z_), (X_, yk, Z_), D_ZYL)))
    teile.append((f"Filzeinlage_R{rr}_{X_:.0f}", stab((X_, ya - dl_, Z_), (X_, ya, Z_), 2 * r_).cut(stab((X_, ya - dl_ - 1, Z_), (X_, ya + 1, Z_), D_ZYL))))
print("Front wie Standardmodell (Reihe 2 mit 7-mm-Blende), Teilung", T_KNOPF, "Knoepfe", len(kn), "X", min(k[0] for k in kn), "-", max(k[0] for k in kn), "| Hoerner", len(aus_front), horn_pos)
# Aufdoppelung innen an der Front oberhalb (Auflage der vorderen Abdeckplatte), ueber den Markus-Aufdoppelungen
# zwischen den Abschlusswaenden der Vertiefung nur die erforderliche Leiste (Nut fuer das Bodenstueck mit der Luftklappe),
# an den Enden (Begleiter) vom Boden bis Z 82
afo = box(59.2, 324.2, FI - 5, FI, 75.0, 82).fuse(box(D, 59.2, FI - 5, FI, Z_BO, 82)).fuse(box(324.2, L - D, FI - 5, FI, Z_BO, 82))
# Im Bereich der Hoerner rechteckig ausgenommen (zwischen Boden und oberem Verleimstreifen), damit Gaze/Horn frei sind

for b_ in knopf_bohr: afo = afo.cut(b_)
teile.append(("Aufdoppelung_Front_oben", afo))
teile.append(("Aufdoppelung_Rueckteil", box(D, L - D, D, D + 5, Z_BO, 82)))
# ---------- Deckelauflage wie Markus (6 mm, vom Deckel bis zur Bodenunterseite), Aussparungen fuer Gaze am Rueckteil ----------
B_AUF, RAND_G, REST = 6.0, 5.0, 4.0
za0 = T_DECKEL_BASS
aus_rueck = []
for e in RD.query('SPLINE[layer=="LAYER_1"]'):
    pp = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.1)]
    if min(b for _, b in pp) > -4.24:
        xs = [a + 382.0 for a, _ in pp]
        aus_rueck.append(box(min(xs) - RAND_G, max(xs) + RAND_G, D - 1, D + B_AUF + 1, za0 + REST, Z_BU + 1))
dh = box(D, L - D, D, D + B_AUF, za0, Z_BU)
for a_ in aus_rueck: dh = dh.cut(a_)
za1 = za0 + 6.0
def mit_aus(s_, l_):
    for a_ in l_: s_ = s_.cut(a_)
    return s_
teile += [("Deckelauflage_hinten", dh), ("Deckelauflage_links", box(D, D + B_AUF, D + B_AUF, FI, za0, Z_BU)),
          ("Deckelauflage_rechts", box(L - D - B_AUF, L - D, D + B_AUF, FI, za0, Z_BU)),
          ("Deckelauflage_vorne_links", mit_aus(box(D + B_AUF, wx0 - VERL_AUF, FI - B_AUF, FI, za0, za1), aus_front)),
          ("Deckelauflage_vorne_rechts", mit_aus(box(wx1 + VERL_AUF, L - D - B_AUF, FI - B_AUF, FI, za0, za1), aus_front))]
# Im Knopfbereich zwei kleine Auflagekloetze (8 breit, 5 hoch) zwischen den Knopfzylindern der unteren Reihe (wie Standardmodell)
_r1 = sorted(k[0] for k in kn if k[3] == 1); _lu = [(a_ + b_) / 2 for a_, b_ in zip(_r1[:-1], _r1[1:])]
for j_, xk_ in enumerate([min(_lu, key=lambda x_: abs(x_ - (wx0 + (wx1 - wx0) * f_))) for f_ in (1 / 3, 2 / 3)]):
    teile.append((f"Deckelauflage_Klotz_vorne_{j_+1}", box(xk_ - 4.0, xk_ + 4.0, FI - B_AUF, FI, za0, za0 + 5.0)))
# ---------- Bassdeckel 4 mm zweigeteilt (hintere Leiste 15) + Eckbloecke + Rosetten (verdeck.dxf) ----------
B_DL, BLK, D_DS, D_BL = 15.0, 15.0, 3.5, 3.0
deckel_h = box(D, L - D, D, D + B_DL, 0.0, T_DECKEL_BASS)
deckel_r = box(D, L - D, D + B_DL, FI, 0.0, T_DECKEL_BASS)
xi0, xi1 = D + B_AUF, L - D - B_AUF
yi0, yi1 = D + B_AUF, FI - B_AUF
ecken = [(xi0, yi0), (xi1 - BLK, yi0), (xi0, yi1 - BLK), (xi1 - BLK, yi1 - BLK)]
for k, (bx, by) in enumerate(ecken):
    cx, cy = bx + BLK / 2, by + BLK / 2
    teile.append((f"Eckblock_{k+1}", box(bx, bx + BLK, by, by + BLK, T_DECKEL_BASS, Z_BU).cut(cq.Solid.makeCylinder(D_BL / 2, 30, cq.Vector(cx, cy, T_DECKEL_BASS - 1)))))
    loch_d = cq.Solid.makeCylinder(D_DS / 2, T_DECKEL_BASS + 2, cq.Vector(cx, cy, -1))
    if cy < D + B_DL: deckel_h = deckel_h.cut(loch_d)
    else: deckel_r = deckel_r.cut(loch_d)
VD = ezdxf.readfile("verdeck.dxf").modelspace()
_ros = sorted({(round(c.dxf.center.x + 381, 2), round(c.dxf.center.y + 206, 2), round(c.dxf.radius, 2)) for c in VD.query('CIRCLE[layer=="0"]')})
for (x_, y_, r_) in _ros: deckel_r = deckel_r.cut(cq.Solid.makeCylinder(r_, T_DECKEL_BASS + 2, cq.Vector(x_, y_, -1)))
# ---------- Boden 5 mm: Y 11..153 ganz, vorne nur an den Enden (Begleiter); Oeffnungen und Klappen aus der Zeichnung ----------
DY_KANAL = 0.2                                 # Kanalreihe der Helikonklappen an den Kanal (123..148)
boden = box(D, L - D, D, 153.0, Z_BU, Z_BO).fuse(box(D, 59.2, 153.0, FI, Z_BU, Z_BO)).fuse(box(324.2, L - D, 153.0, FI, Z_BU, Z_BO))   # an den Enden bis an Front und Rueckteil
oeff, kl_ = [], []
for e in A.query("LWPOLYLINE"):
    pts = [(p[0] - X0 + OX, p[1] - Y0P + OY) for p in e.get_points()]
    xs = [p[0] for p in pts]; ys = [p[1] for p in pts]
    if not (0 < min(xs) and max(xs) < L and 0 < min(ys) and max(ys) < B): continue
    dy_ = DY_KANAL if min(ys) > 100 else 0.0
    if e.dxf.layer.startswith("FULLUNG-RUND") and len(pts) == 4: oeff.append((min(xs), max(xs), min(ys) + dy_, max(ys) + dy_))
    if e.dxf.layer == "0KLAPPEN-RAND-AUSSEN": kl_.append((min(xs), max(xs), min(ys) + dy_, max(ys) + dy_))
# Klappen wie Markus: auf der Mechanikseite um jede Oeffnung eine Vertiefung (5 mm Rand, 2 mm tief), Klappe 0,5 mm
# kleiner als die Vertiefung (Leder 1, Filz 4, Epoxy 1,5), liegt im Grund der Vertiefung
RAND_V, TIEF_V, SPIEL = 5.0, 2.0, 0.5
klappen = []
oeff.sort(key=lambda k: (-k[2], k[0]))
for i, (x0, x1, y0, y1) in enumerate(oeff):
    boden = boden.cut(box(x0, x1, y0, y1, Z_BU - 1, Z_BO + 1)).cut(box(x0 - RAND_V, x1 + RAND_V, y0 - RAND_V, y1 + RAND_V, Z_BU - 1, Z_BU + TIEF_V))
    z_ = Z_BU + TIEF_V; o_ = RAND_V - SPIEL
    for nm, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
        klappen.append((f"Klappe_{i+1}_{nm}", box(x0 - o_, x1 + o_, y0 - o_, y1 + o_, z_ - dk, z_))); z_ -= dk
nb = 0
for c in sorted(A.query("CIRCLE"), key=lambda c: (c.dxf.center.x, c.dxf.center.y)):
    if c.dxf.layer.startswith("FULLUNG-RUND") and abs(c.dxf.radius - 9) < 0.1 and Y0P < c.dxf.center.y < Y0P + 205:
        cx_, cy_ = c.dxf.center.x - X0 + OX, c.dxf.center.y - Y0P + OY
        boden = boden.cut(cq.Solid.makeCylinder(9, 10, cq.Vector(cx_, cy_, Z_BU - 1))).cut(cq.Solid.makeCylinder(9 + RAND_V, TIEF_V + 1, cq.Vector(cx_, cy_, Z_BU - 1)))
        nb += 1; z_ = Z_BU + TIEF_V
        for nm, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
            klappen.append((f"Begleiterklappe_{nb}_{nm}", cq.Solid.makeCylinder(9 + RAND_V - SPIEL, dk, cq.Vector(cx_, cy_, z_ - dk)))); z_ -= dk
kl_ = oeff
print("Boden: Helikonoeffnungen", len(oeff), "Begleiter rund", nb, "| Klappen", len(kl_))
# ---------- Vordere Abdeckplatte 2,5 mm ueber dem Mechanikbereich (Vertiefung), Luftoeffnung gespiegelt (s. Luftklappe) ----------
Z_ABD = 79.5
abd = box(59.2 - 1.5, 324.2 + 1.5, 153.0 - 1.5, FI - 5 + 1.5, Z_ABD, Z_ABD + 2.5)   # ringsum 1,5 mm in Nuten eingelassen
lo_x0, lo_x1 = 75.5 + OX, 165.5 + OX
abd = abd.cut(box(lo_x0, lo_x1, 166.2 + OY, 186.2 + OY, Z_ABD - 1, Z_ABD + 4))
teile.append(("Abdeckplatte_vorne", abd))
# Abschlusswaende der Vertiefung an beiden Stirnseiten (5 mm, auf dem Boden bis Oberkante Abdeckplatte)
teile.append(("Abschlusswand_Vertiefung_links", box(59.2 - 5.0, 59.2, 153.0, FI - 5, Z_BO, Z_ABD + 2.5)))
teile.append(("Abschlusswand_Vertiefung_rechts", box(324.2, 324.2 + 5.0, 153.0, FI - 5, Z_BO, Z_ABD + 2.5)))
# ---------- Helikonkammern: 5 Bloecke (Teilung 54, Trennwaende 3,5 in die Laengsteile eingelassen, Nut 2 mm),
# untere Kammer + Verbindungskanal an der Rueckwand, obere Etage in der Tiefe zentriert. Innenraum je ca. 1 mm groesser
# als der Plattenkoerper (Stimmplatte 53 x 98 mit 1,5-mm-Randabstufung, liegt nur mit dem Rand auf).
XK0 = 53.5 + OX; XK1 = XK0 + 5 * 54 + 3.5
PL_L, PL_B, RST = 98.0, 53.0, 1.5
T_NUT = 2.0
Z_U1, Z_E0, Z_E1, Z_O1 = 85.0, 109.0, 112.0, 154.0
T_DICHT = 1.0
Z_KO = Z_E0 - T_DICHT                          # Oberkante Verbindungskanal (darauf die Dichtung)
YU_W = D + 5                                   # Rueckwand-Aufdoppelung innen
YU0 = YU_W + 6 - 1.0
YU1 = YU0 + PL_L
YKH0 = YU1 - 1.0                               # Stufe (vorne) innen
T_STUFE = 5.5
YK0 = YKH0 + T_STUFE + 5.0; YK1 = YK0 + 25.0     # Kanal innen; Kanal-Vorderwand endet an der Bodenkante 153
YO0 = B / 2 - PL_L / 2; YO1 = YO0 + PL_L       # obere Etage in der Tiefe mittig; steht vorne ca. 4,5 mm ueber den Kanal vor
wand_h = py(XK0, XK1, [(YU_W, Z_BO), (YU_W, Z_U1), (YU0 + 1, Z_U1), (YU0 + 1, Z_BO)])
wand_m = py(XK0, XK1, [(YKH0, Z_BO), (YKH0, Z_U1), (YKH0 + T_STUFE, Z_U1), (YKH0 + T_STUFE, Z_KO), (YK0, Z_KO), (YK0, Z_BO)])
wand_k = py(XK0, XK1, [(YK1, Z_BO), (YK1, Z_KO), (YK1 + 5, Z_KO), (YK1 + 5, Z_BO)])
oben_b = py(XK0, XK1, [(YO0 + 1 - 6, Z_E0), (YO0 + 1 - 6, Z_E1), (YO1 - 1 + 6, Z_E1), (YO1 - 1 + 6, Z_E0)])
oben_h = py(XK0, XK1, [(YO0 + 1 - 6, Z_E1), (YO0 + 1 - 6, Z_O1), (YO0 + 1, Z_O1), (YO0 + 1, Z_E1)])
oben_v = py(XK0, XK1, [(YO1 - 1, Z_E1), (YO1 - 1, Z_O1), (YO1 - 1 + 6, Z_O1), (YO1 - 1 + 6, Z_E1)])
for k in range(5):
    xa = XK0 + 3.5 + 54 * k
    oben_b = oben_b.cut(box(xa, xa + 50.5, YK0, YK1, Z_E0 - 1, Z_E1 + 1))
tw_l = []
for k in range(6):
    xa = XK0 + 54 * k
    tw_l.append((f"HB_Trennwand_unten_{k+1}", box(xa, xa + 3.5, YU0 + 1 - T_NUT, YKH0 + T_NUT, Z_BO, Z_U1)))
    tw_l.append((f"HB_Trennwand_Kanal_{k+1}", box(xa, xa + 3.5, YK0 - T_NUT, YK1 + T_NUT, Z_BO, Z_KO)))
    tw_l.append((f"HB_Trennwand_oben_{k+1}", box(xa, xa + 3.5, YO0 + 1 - T_NUT, YO1 - 1 + T_NUT, Z_E1, Z_O1)))
laengs = []
for n_, w_ in [("HB_Laengswand_hinten", wand_h), ("HB_Laengswand_Kanal_hinten", wand_m), ("HB_Laengswand_Kanal_vorne", wand_k),
               ("HB_Etagenboden", oben_b), ("HB_Etage_Wand_hinten", oben_h), ("HB_Etage_Wand_vorne", oben_v)]:
    for _, tw in tw_l: w_ = w_.cut(tw)                  # Nuten fuer die Trennwaende
    laengs.append((n_, w_))
teile += laengs
teile += tw_l
dicht = box(XK0, XK1, YK0 - 5, YK1 + 5, Z_KO, Z_E0)
for k in range(5):
    xa = XK0 + 3.5 + 54 * k
    dicht = dicht.cut(box(xa, xa + 50.5, YK0, YK1, Z_KO - 1, Z_E0 + 1))
teile.append(("Dichtung_Kanal", dicht))
# Formteile 8 mm nach HB-trennwaende.dxf (Konturen mit Rundungen R4 und Auge um die Bohrung Ø4), flach liegend:
# U am unteren Teil (Steg vor der Stirnseite, Schenkel umfassen den Kanal 35 mm, Oberkante = Kanaloberkante),
# L an der Etage (Unterkante = Etagenboden, kurzer Schenkel vor der Etagen-Vorderwand); je eine Schraube von oben durch L in U.
# Hinten in der Mitte: I unten an der Rueckseite der Kanalwand, I oben unter dem Etagenboden, Schraube von unten.
HBD = ezdxf.readfile("hb_trennwaende.dxf").modelspace()
def kontur(sel):
    for e in HBD.query("LWPOLYLINE"):
        if e.dxf.layer != "0": continue
        bb_ = bbox.extents([e])
        if sel(bb_): return [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.05)], bb_
K_I, _ = kontur(lambda b: abs(b.extmin.x - 61.73) < 0.2 and b.extmin.y > 70)
K_L, _ = kontur(lambda b: abs(b.extmin.x - 61.8) < 0.2 and 40 < b.extmin.y < 45)
K_U, _ = kontur(lambda b: abs(b.extmin.x - 65.8) < 0.2 and b.extmin.y < 15)
K_L = K_L + [K_L[0]] if math.dist(K_L[0], K_L[-1]) > 0.01 else K_L
K_U = [((74.3 if abs(x - 73.8) < 0.01 else 109.3 if abs(x - 109.8) < 0.01 and y < 24 else x), y) for x, y in K_U]   # Oeffnung 36 -> 35
def platte(pts, z0, z1):
    q_ = [pts[0]] + [b for a, b in zip(pts[:-1], pts[1:]) if math.dist(a, b) > 1e-4]
    if math.dist(q_[0], q_[-1]) < 1e-4: q_ = q_[:-1]
    return cq.Workplane("XY").polyline(q_).close().extrude(z1 - z0).translate((0, 0, z0)).val()
YC = (YK0 + YK1) / 2; T_FT = 8.0
UEB = (YO1 - 1 + 6) - (YK1 + 5)                 # Ueberstand der Etage vorne ueber den Kanal
YS = YC + UEB                                   # Schraubenlage: L-Schenkel liegt vor der Etagen-Vorderwand
ft_l = []
# U und L moeglichst gleich: gleicher Steg (10 breit vor der Stirnseite, Auge R11 um die Schraube, Ecken R4), von oben
# deckungsgleich im Schraubenbereich. L: ein Schenkel vorn vor der Etagen-Vorderwand. U: Schenkel hinten und vorn um den
# Kanal (35 mm); U zu den Begleiterstimmstoecken hin mit 5-mm-Fase an der unteren Aussenkante.
YA_F = YK0 - 5 - 8.0                         # Steg hinten (U-Schenkel hinter dem Kanal)
YB_F = YO1 - 1 + 6 + 8.0                     # Steg vorne (L-Schenkel vor der Etagen-Vorderwand)
def ft_form(sg, xe, z0, z1, schenkel):
    # u = Y, v = Abstand von der Stirnseite nach aussen (positiv), Schenkel nach innen (negativ)
    X_ = lambda v: xe + sg * v
    xa, xb = sorted((X_(0), X_(10)))
    t_ = box(xa, xb, YA_F, YB_F, z0, z1)
    for (u0, u1) in schenkel:
        xi0, xi1 = sorted((X_(-10), X_(0)))
        t_ = t_.fuse(box(xi0, xi1, u0, u1, z0, z1))
    auge = cq.Solid.makeCylinder(11.0, z1 - z0, cq.Vector(X_(5.0), YS, z0)).intersect(box(min(X_(0), X_(20)), max(X_(0), X_(20)), YS - 12, YS + 12, z0, z1))
    t_ = t_.fuse(auge).clean()
    ecken = [(X_(10), YA_F), (X_(10), YB_F)] + [(X_(-10), u_) for (u0, u1) in schenkel for u_ in (u0, u1) if abs(u_ - YA_F) < 0.01 or abs(u_ - YB_F) < 0.01]
    for (ex, ey) in ecken:
        try:
            e_ = [e for e in t_.Edges() if abs(e.Center().x - ex) < 0.05 and abs(e.Center().y - ey) < 0.05 and abs(e.Center().z - (z0 + z1) / 2) < 0.05]
            if e_: t_ = t_.fillet(3.99, e_)
        except Exception: pass
    return t_
for sg in (-1, 1):
    xe = XK0 if sg < 0 else XK1
    seite = "links" if sg < 0 else "rechts"
    u_ = ft_form(sg, xe, Z_KO - T_FT, Z_KO, [(YA_F, YK0 - 5), (YK1 + 5, YB_F)])
    l_ = ft_form(sg, xe, Z_E0, Z_E0 + T_FT, [(YO1 - 1 + 6, YB_F)])
    # Fase 5 mm (45 Grad) an der unteren Aussenkante des U (zu den Begleiterstimmstoecken): weg, was v - 5 > Z - Zu
    zu_ = Z_KO - T_FT
    pk_ = [(xe + sg * 5.0, zu_ - 0.01), (xe + sg * 30.0, zu_ - 0.01), (xe + sg * 30.0, zu_ + 25.0)]
    keil = cq.Workplane("XZ").polyline(pk_).close().extrude(-(YB_F - YA_F + 30)).translate((0, YA_F - 15, 0)).val()
    u_ = u_.cut(keil)
    xs = xe + sg * 5.0
    u_ = u_.cut(cq.Solid.makeCylinder(1.75, T_FT + 2, cq.Vector(xs, YS, Z_KO - T_FT - 1)))
    l_ = l_.cut(cq.Solid.makeCylinder(2.0, T_FT + 2, cq.Vector(xs, YS, Z_E0 - 1)))
    sc_ = cq.Solid.makeCylinder(3.5, 2.5, cq.Vector(xs, YS, Z_E0 + T_FT)).fuse(cq.Solid.makeCylinder(1.75, 14.0, cq.Vector(xs, YS, Z_E0 + T_FT - 14.0)))
    ft_l += [(f"Formteil_U_{seite}", u_), (f"Formteil_L_{seite}", l_), (f"Formteil_Schraube_{seite}", sc_)]
# I-Paar vorne in der Mitte an derselben Trennfuge (Dichtung) wie U/L: unteres I an der Kanal-Vorderwand (Oberkante =
# Kanaloberkante), um den Ueberstand der Etage breiter; oberes I an der Etagen-Vorderwand (Unterkante = Etagenboden);
# von oben in Deckung, Schraube von oben durch das obere in das untere I.
yk_v, ye_v = YK1 + 5.0, YO1 - 1 + 6
def mI_u(x, y):
    t = y - 72.28
    return (L / 2 + (x - 91.73), yk_v + (t * (10.0 + UEB) / 10.0 if t <= 10.0 else t + UEB))
mI_o = lambda x, y: (L / 2 + (x - 91.73), ye_v + (y - 72.28))
i_u = platte([mI_u(x, y) for x, y in K_I], Z_KO - T_FT, Z_KO)
i_o = platte([mI_o(x, y) for x, y in K_I], Z_E0, Z_E0 + T_FT)
yI = ye_v + 5.0
i_u = i_u.cut(cq.Solid.makeCylinder(1.75, T_FT + 2, cq.Vector(L / 2, yI, Z_KO - T_FT - 1)))
i_o = i_o.cut(cq.Solid.makeCylinder(2.0, T_FT + 2, cq.Vector(L / 2, yI, Z_E0 - 1)))
sc_i = cq.Solid.makeCylinder(3.5, 2.5, cq.Vector(L / 2, yI, Z_E0 + T_FT)).fuse(cq.Solid.makeCylinder(1.75, 14.0, cq.Vector(L / 2, yI, Z_E0 + T_FT - 14.0)))
ft_l += [("Formteil_I_unten", i_u), ("Formteil_I_oben", i_o), ("Formteil_Schraube_mitte", sc_i)]
print(f"Formteile: Etage steht vorne {UEB:.1f} mm ueber, unteres I um {UEB:.1f} breiter, Schrauben U/L bei Y {YS:.1f}, I bei Y {yI:.1f}")
teile += ft_l
print("Dichtung", T_DICHT, "mm auf dem Kanal")
T_FOLIE = 0.6
def folie(xa, xb, y_h, y_v, z0):
    lang = y_v - y_h; h, r = 16.5, 36.3; t_kn = lang - 23.2
    prof = [(23.2, h), (t_kn, h)] + [(t_kn + 23.2 * i_ / 12, -19.7 + math.sqrt(r * r - (23.2 * i_ / 12) ** 2)) for i_ in range(1, 13)]
    oben = [(y_h + t, z0 + z + T_FOLIE / 2) for t, z in prof]; unten = [(y_h + t, z0 + z - T_FOLIE / 2) for t, z in reversed(prof)]
    return py(xa, xb, oben + unten)
# Stimmplatten (Markus/Standard): Rand 2 mm dick, 1,5 mm Abstufung ringsum (liegt nur damit auf); Plattenkoerper darunter
# 3 mm kleiner, keilfoermig 5 -> 16 (dickes Ende vorne); Zungenschlitze 7,5/9,2/9,2/7,5 enden 5,5 vor dem dicken Ende
T_DUENN, T_DICK, T_RAND = 5.0, 16.0, 2.0
def stimmplatte(px0, px1, y0, y1, z1):
    pl = box(px0, px1, y0, y1, z1, z1 + T_RAND)
    a0, a1 = y0 + RST, y1 - RST
    pl = pl.fuse(py(px0 + RST, px1 - RST, [(a0, z1 + 0.01), (a1, z1 + 0.01), (a1, z1 + T_RAND - T_DICK), (a0, z1 + T_RAND - T_DUENN)]))
    for (o_, sb, sll) in [(5.1, 7.5, 51.0), (15.9, 9.2, 71.5), (27.9, 9.2, 71.5), (40.4, 7.5, 51.0)]:
        pl = pl.cut(box(px0 + o_, px0 + o_ + sb, y1 - 5.5 - sll, y1 - 5.5, z1 - T_DICK - 1, z1 + T_RAND + 1))
    return pl
# Geschwungene Trennfolie 0,6 durchgehend durch alle 5 Kammern (ein Stueck je Ebene), die Kammer-Trennwaende sind dafuer geschlitzt
fol_u = folie(XK0 + 1.75, XK1 - 1.75, YU0 + 1, YKH0, Z_BO); fol_o = folie(XK0 + 1.75, XK1 - 1.75, YO0 + 1, YO1 - 1, Z_E1)
teile.append(("HB_Folie_unten", fol_u)); teile.append(("HB_Folie_oben", fol_o))
teile = [(n_, s_.cut(fol_u).cut(fol_o) if n_.startswith(("HB_Trennwand_unten_", "HB_Trennwand_oben_")) else s_) for n_, s_ in teile]
for k in range(5):
    xa = XK0 + 3.5 + 54 * k; xb = xa + 50.5
    px0 = XK0 + 54 * k + 1.75 + 0.5; px1 = px0 + PL_B
    teile.append((f"HB_Stimmplatte_unten_{k+1}", stimmplatte(px0, px1, YU0, YU1, Z_U1)))
    teile.append((f"HB_Stimmplatte_oben_{k+1}", stimmplatte(px0, px1, YO0, YO1, Z_O1)))
# Befestigung wie Markus: je Platte vorne und hinten ein Schraubhaken (Draht 1,8, kurzer Schenkel 5 auf dem Plattenrand,
# langer Schenkel 15 ins Holz der Laengswand), mittig in der Plattenbreite
D_HAKEN, L_KURZ, L_LANG = 1.8, 5.0, 15.0
def zyl(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1)
    return cq.Solid.makeCylinder(d / 2, (p1 - p0).Length, p0, (p1 - p0).normalized())
def haken(px, py_, dy, z1):
    zh = z1 + T_RAND + D_HAKEN / 2
    return zyl((px, py_, zh - L_LANG + D_HAKEN / 2), (px, py_, zh), D_HAKEN).fuse(zyl((px, py_, zh), (px, py_ + dy * L_KURZ, zh), D_HAKEN)).fuse(cq.Solid.makeSphere(D_HAKEN / 2, cq.Vector(px, py_, zh)))
haken_l = []
off = 0.1 + D_HAKEN / 2
for k in range(5):
    xm = XK0 + 54 * k + 1.75 + 0.5 + PL_B / 2
    for (y0_, y1_, z1_, nm) in ((YU0, YU1, Z_U1, "unten"), (YO0, YO1, Z_O1, "oben")):
        haken_l.append((f"Haken_{nm}_{k+1}_hinten", haken(xm, y0_ - off, 1, z1_)))
        haken_l.append((f"Haken_{nm}_{k+1}_vorne", haken(xm, y1_ + off, -1, z1_)))
for k_ in range(len(teile)):
    n_, s_ = teile[k_]
    if n_.startswith("HB_Laengswand") or n_.startswith("HB_Etage_Wand"):
        for _, h_ in haken_l:
            if s_.intersect(h_).Volume() > 0.01: s_ = s_.cut(h_)
        teile[k_] = (n_, s_)
teile += haken_l
print(f"Kammern: unten Platte Y {YU0:.1f}-{YU1:.1f}, Kanal {YK0:.1f}-{YK1:.1f}, Etage Platte Y {YO0:.1f}-{YO1:.1f} (Mitte {B/2:.1f})")
# ---------- Begleiterstimmstoecke links/rechts im Detail (wie Markus, je 8 x 2 Kammern = 32 gesamt) ----------
# Sole 6 mm (Zeichnung 33 x 146) mit einer Oeffnung je Kammer; Block 20 breit, 50 hoch, Kammern 15 breit (Stege 2 wie
# Zeichnung), 45 hoch, 9,5 tief von beiden Seiten (Mittelwand 1); je Kammer eine Stimmplatte 2,5 x 16,8 x 47 aussen.
# Je zwei Kammern liegen ueber einer runden Bodenoeffnung Ø18 (Teilung 34).
T_SOLE, B_ST, H_ST, K_H, K_T = 6.0, 20.0, 50.0, 45.0, 9.5
T_SPL, B_SPL, H_SPL = 2.5, 16.8, 47.0
for nm, x0 in (("links", 14.6 + OX), ("rechts", 332.4 + OX)):
    sy0, sy1 = 27.1 + OY - 6.0, 173.1 + OY + 6.0           # Sole an beiden Enden 6 mm laenger: Schraube in den Boden
    xm = x0 + 16.5; bx0, bx1 = xm - B_ST / 2, xm + B_ST / 2
    sole = box(x0, x0 + 33, sy0, sy1, Z_BO, Z_BO + T_SOLE)
    blk = box(bx0, bx1, 31.1 + OY, 169.1 + OY, Z_BO + T_SOLE, Z_BO + T_SOLE + H_ST)
    plt_ = []
    for j in range(8):
        ky0 = 33.1 + OY + 17 * j; ky1 = ky0 + 15
        for (xa_, xb_, sg) in ((bx0 - 1, bx0 + K_T, -1), (bx1 - K_T, bx1 + 1, 1)):
            blk = blk.cut(box(xa_, xb_, ky0, ky1, Z_BO + T_SOLE - 1, Z_BO + T_SOLE + K_H))
            sole = sole.cut(box(max(xa_, bx0), min(xb_, bx1), ky0, ky1, Z_BO - 1, Z_BO + T_SOLE + 1))
            xp = bx0 - T_SPL if sg < 0 else bx1
            yc = (ky0 + ky1) / 2
            zb_ = Z_BO + T_SOLE + 0.5
            spl_ = box(xp, xp + T_SPL, yc - B_SPL / 2, yc + B_SPL / 2, zb_, zb_ + H_SPL)
            for dy_ in (-3.5, 3.5):                            # zwei Zungenschlitze 2,6 x 34 (ab 6 ueber Unterkante), Nietloch Ø1,6
                spl_ = spl_.cut(box(xp - 1, xp + T_SPL + 1, yc + dy_ - 1.3, yc + dy_ + 1.3, zb_ + 6.0, zb_ + 40.0))
                spl_ = spl_.cut(stab((xp - 1, yc + dy_, zb_ + 3.0), (xp + T_SPL + 1, yc + dy_, zb_ + 3.0), 1.6))
            plt_.append(spl_)
    for k_, ys_ in enumerate((sy0 + 5.0, sy1 - 5.0)):
        sole = sole.cut(cq.Solid.makeCylinder(1.75, T_SOLE + 2, cq.Vector(xm, ys_, Z_BO - 1)))
        sc_ = cq.Solid.makeCylinder(3.0, 2.0, cq.Vector(xm, ys_, Z_BO + T_SOLE)).fuse(cq.Solid.makeCylinder(1.5, 10.0, cq.Vector(xm, ys_, Z_BO + T_SOLE - 10.0)))
        boden = boden.cut(cq.Solid.makeCylinder(1.25, 4.5, cq.Vector(xm, ys_, Z_BO - 4.5)))
        teile.append((f"Begleiter_Schraube_{nm}_{k_+1}", sc_))
    teile.append((f"Begleiter_Sole_{nm}", sole)); teile.append((f"Begleiterstimmstock_{nm}", blk))
    for q_, p_ in enumerate(plt_): teile.append((f"Begleiter_Stimmplatte_{nm}_{q_+1}", p_))
lg = []
# --- Luftklappe wie Standardmodell, in der vorderen Vertiefung unter der Abdeckplatte (Brett waagrecht, Unterseite Z 79,5) ---
import math, ezdxf.path
def kugel(p, d): return cq.Solid.makeSphere(d / 2 + 0.01, cq.Vector(*p), angleDegrees1=-90, angleDegrees2=90)
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v_ = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v_.Length, p0, v_.normalized())
def pfad(pts, d):
    s_ = None
    for a_, b_ in zip(pts[:-1], pts[1:]):
        seg = stab(a_, b_, d).fuse(kugel(b_, d)); s_ = seg if s_ is None else s_.fuse(seg)
    return s_.fuse(kugel(pts[0], d))
def pts_of(e, tol=0.1):
    p_ = [(v_.x, v_.y) for v_ in ezdxf.path.make_path(e).flattening(tol)]
    q_ = [p_[0]] + [b_ for a_, b_ in zip(p_[:-1], p_[1:]) if math.dist(a_, b_) > 1e-6]
    if len(q_) > 2 and math.dist(q_[0], q_[-1]) < 1e-6: q_ = q_[:-1]
    return q_
T_LEDER_K, T_FILZ_K, T_KDECK = 0.5, 4.0, 0.5
T_KLK = T_FILZ_K + T_KDECK
DR_F, DI_F, N_F = 0.5, 2.8, 3
Z_B0 = Z_ABD
RAND_LK = 4.5
luft_kl = (75.5 + OX - RAND_LK, 165.5 + OX + RAND_LK, 166.2 + OY - RAND_LK, 186.2 + OY + RAND_LK)
x0_, x1_, y0_, y1_ = luft_kl
lg.append(("Luftklappe_Leder", box(x0_, x1_, y0_, y1_, Z_B0 - T_LEDER_K, Z_B0)))
lg.append(("Luftklappe_Filz", box(x0_ + 0.5, x1_ - 0.5, y0_ + 0.5, y1_ - 0.5, Z_B0 - T_LEDER_K - T_FILZ_K, Z_B0 - T_LEDER_K)))
lg.append(("Luftklappe", box(x0_ + 0.5, x1_ - 0.5, y0_ + 0.5, y1_ - 0.5, Z_B0 - T_LEDER_K - T_KLK, Z_B0 - T_LEDER_K - T_FILZ_K)))
# Lufttaste wie Markus/Standard: Langloch mit halbrunden Enden aus dem Def-Rechteck des Deckels, Filz 1, Hub 10
LT_X, LT_Y, LT_L, LT_B = (207.5 + 285.5) / 2 + OX, (174.0 + 188.0) / 2 + OY, 285.5 - 207.5, 188.0 - 174.0
T_FILZ_LT, HUB_LT = 1.0, 10.0
def langloch(l_, b_, z0_, h_):
    return cq.Workplane("XY").center(LT_X, LT_Y).slot2D(l_, b_).extrude(h_).translate((0, 0, z0_)).val()
lt_filz = langloch(LT_L, LT_B, 0, T_DECKEL_BASS).cut(langloch(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT, -1, T_DECKEL_BASS + 2))
taste = cq.Workplane("XY").center(LT_X, LT_Y).slot2D(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT).extrude(T_DECKEL_BASS + HUB_LT).faces("<Z").edges().fillet(1.0).val()
taste = taste.translate(cq.Vector(0, 0, -HUB_LT))
lg.append(("Lufttaste_Filz", lt_filz)); lg.append(("Lufttaste", taste))
LH = ezdxf.readfile("kb_lufthebel.dxf").modelspace()
LHX = LT_X + 12.65 + 5.5                                          # 5,5 nach rechts: Knick des Arms ausserhalb der Luftklappe
T_LH = 3.0
LY = (luft_kl[2] + luft_kl[3]) / 2
Z_LH = Z_B0 - (T_KLK + T_LEDER_K) - 0.05 + 5.35
def lh_xz(p_): return [(a_ + LHX, Z_LH - b_) for a_, b_ in p_]
def platte_xz(pts, ya, yb): return cq.Workplane("XZ").polyline(pts).close().extrude(-(yb - ya)).translate((0, ya, 0)).val()
hp_ = lh_xz(pts_of([e for e in LH.query('LWPOLYLINE[layer=="2-Gelenk-makierung"]')][0], 0.05))
kr = [(c.dxf.center.x + LHX, Z_LH - c.dxf.center.y, c.dxf.radius) for c in LH.query('CIRCLE[layer=="0"]')]
XP, ZP = [(k_[0], k_[1]) for k_ in kr if abs(k_[2] - 1.0) < 0.01 and abs(k_[0] - (-55.0 + LHX)) < 0.6][0]
XA, ZA_ = [(k_[0], k_[1]) for k_ in kr if abs(k_[2] - 0.54) < 0.01][0]
hebel = platte_xz(hp_, LY - T_LH / 2, LY + T_LH / 2).fuse(stab((XP, LY - T_LH / 2, ZP), (XP, LY + T_LH / 2, ZP), 5.0))
hebel = hebel.cut(stab((XP, LY - 3, ZP), (XP, LY + 3, ZP), 2.05)).cut(stab((XA, LY - 3, ZA_), (XA, LY + 3, ZA_), 1.08))
lg.append(("Lufthebel", hebel))
L_FL = N_F * (DR_F + 0.15) + DR_F; W_LB = 8.0
y_f = [(LY - T_LH / 2 - 0.15 - L_FL, LY - T_LH / 2 - 0.15), (LY + T_LH / 2 + 0.15, LY + T_LH / 2 + 0.15 + L_FL)]
y_w = [(y_f[0][0] - 0.15 - W_LB, y_f[0][0] - 0.15), (y_f[1][1] + 0.15, y_f[1][1] + 0.15 + W_LB)]
bk = [p_ for p_ in [pts_of(e, 0.05) for e in LH.query('LWPOLYLINE[layer=="2-Gelenke"]')] if min(a_ for a_, _ in p_) < -40][0]
xb0, xb1 = min(a_ for a_, _ in lh_xz(bk)), max(a_ for a_, _ in lh_xz(bk))
x_schr = -64.0 + LHX
schr = []
for k_, (ya_, yb_) in enumerate(y_w):
    w_ = platte_xz(lh_xz(bk), ya_, yb_).fuse(box(xb0, xb1, ya_, yb_, Z_LH - 0.5, Z_B0 + 3.0))
    w_ = w_.cut(box(xb0 - 1, xb1 + 1, ya_ - 1, yb_ + 1, Z_B0, Z_B0 + 10)).cut(stab((XP, LY - 30, ZP), (XP, LY + 30, ZP), 2.05))
    for dy_ in (-2.0, 2.0):
        yc_ = (ya_ + yb_) / 2 + dy_
        sc_ = stab((x_schr, yc_, Z_LH - 2.53 - 0.9), (x_schr, yc_, Z_LH - 2.53), 2.6).fuse(stab((x_schr, yc_, Z_LH - 2.53), (x_schr, yc_, Z_B0 + 2.0), 1.4))
        w_ = w_.cut(stab((x_schr, yc_, Z_LH - 4), (x_schr, yc_, Z_B0 + 1), 1.5)); schr.append(sc_)
    lg.append((f"Lufthebel_Stuetze_{k_+1}", w_))
for j_, sc_ in enumerate(schr): lg.append((f"Lufthebel_Stuetze_Schraube_{j_+1}", sc_))
lg.append(("Lufthebel_Bolzen", stab((XP, y_w[0][0] - 0.6, ZP), (XP, y_w[1][1] + 0.6, ZP), 2.0)))
def kante_oben(x_):
    zz_ = []
    for (xa_, za_), (xb_, zb_) in [(hp_[i], hp_[i + 1]) for i in range(len(hp_) - 1)]:
        if min(xa_, xb_) <= x_ <= max(xa_, xb_) and abs(xb_ - xa_) > 1e-6: zz_.append(za_ + (zb_ - za_) * (x_ - xa_) / (xb_ - xa_))
    return max(zz_)
rmf = (DI_F + DR_F) / 2
for k_, ((ya0, ya1), sg) in enumerate(((y_f[0], -1), (y_f[1], 1))):
    y_in = ya1 if sg < 0 else ya0; y_out = ya0 if sg < 0 else ya1
    h_ = cq.Wire.makeHelix((L_FL - DR_F) / N_F, L_FL - DR_F, rmf)
    pr = cq.Wire.makeCircle(DR_F / 2, cq.Vector(rmf, 0, 0), cq.Vector(0, 1, 0))
    sp = cq.Solid.sweep(pr, [], h_, isFrenet=True).rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -90 * sg)
    sp = sp.translate(cq.Vector(XP, y_in + sg * DR_F / 2, ZP))
    p0 = (XP + rmf, y_in + sg * DR_F / 2, ZP); p1 = (XP + rmf, y_out - sg * DR_F / 2, ZP)
    x_u = XP + 15.0 + 0.6 * sg
    z_u = max(kante_oben(XP + 14.4), kante_oben(XP + 15.6)) + DR_F / 2 + 0.05 + (DR_F + 0.15 if sg > 0 else 0.0)
    y_a, y_b = LY + sg * (T_LH / 2 + DR_F / 2 + 0.05), LY - sg * (T_LH / 2 + DR_F / 2 + 0.05)
    haken = pfad([p0, (XP + rmf + 3.0, y_a, ZP + 0.5), (XP + 13.4, y_a, z_u), (x_u, y_a, z_u), (x_u, y_b, z_u), (x_u, y_b, z_u - 2.0)], DR_F)
    zb_ = Z_B0 - DR_F / 2 - 0.05
    s1 = pfad([p1, (XP + rmf + 3.0, y_out - sg * DR_F / 2, zb_), (XP + rmf + 12.0, y_out - sg * DR_F / 2, zb_)], DR_F)
    lg.append((f"Lufthebel_Feder_{k_+1}", sp.fuse(haken).fuse(s1)))
G_L, D_GS = 10.0, 6.0
gg = stab((XA, LY, ZA_ + 3.0), (XA, LY, ZA_ - G_L), D_GS).cut(box(XA - 4, XA + 4, LY - T_LH / 2 - 0.4, LY + T_LH / 2 + 0.4, ZA_ - 3.0, ZA_ + 4))
gg = gg.cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6)).cut(stab((XA, LY, ZA_ - 3.5), (XA, LY, ZA_ - G_L - 1), 3.0))
E0 = cq.Vector(XA, LY, ZA_ - 3.6); Pb = cq.Vector(XA, LT_Y, T_DECKEL_BASS - HUB_LT - 0.5)   # Stange senkrecht zur Taste (Auge liegt ueber der Taste)
dw = math.degrees(math.atan2(LT_Y - LY, -(Pb.z - ZA_)))
gg = gg.rotate(cq.Vector(XA, LY, ZA_), cq.Vector(XA + 1, LY, ZA_), dw)
lg.append(("Lufttaste_Gabelgelenk", gg)); lg.append(("Lufttaste_Gelenkbolzen", stab((XA, LY - D_GS / 2, ZA_), (XA, LY + D_GS / 2, ZA_), 1.0)))
dvec = (Pb - cq.Vector(XA, LY, ZA_)).normalized(); E = cq.Vector(XA, LY, ZA_) + dvec * 3.6
lg.append(("Lufttaste_Schubstange", stab(E.toTuple(), Pb.toTuple(), 3.0).cut(taste)))
# Lufttaste gedrueckt (eben mit dem Deckel) halten, waehrend die Mechanik eingebaut wird: Hebel/Klappe um das Lager schwenken
LT_GEDRUECKT = True
if LT_GEDRUECKT:
    A0 = cq.Vector(XA, LY, ZA_); P0 = cq.Vector(XP, LY, ZP); Pb2 = Pb + cq.Vector(0, 0, HUB_LT); l_st = (Pb - A0).Length
    def _A(ph_):
        c_, s_ = math.cos(ph_), math.sin(ph_); dx_, dz_ = XA - XP, ZA_ - ZP
        return cq.Vector(XP + c_ * dx_ - s_ * dz_, LY, ZP + s_ * dx_ + c_ * dz_)
    best = min((abs((Pb2 - _A(math.radians(g_ / 20))).Length - l_st), g_ / 20) for g_ in range(-900, 901))
    PHI = best[1]
    rot = lambda s_: s_.rotate(cq.Vector(XP, LY, ZP), cq.Vector(XP, LY - 1, ZP), PHI)      # Drehung um Y durch das Lager (Vorzeichen wie _A)
    if abs(rot(cq.Vertex.makeVertex(XA, LY, ZA_)).Center().z - _A(math.radians(PHI)).z) > 0.01:
        rot = lambda s_: s_.rotate(cq.Vector(XP, LY, ZP), cq.Vector(XP, LY + 1, ZP), PHI)
    A2 = _A(math.radians(PHI))
    gg2 = stab((XA, LY, ZA_ + 3.0), (XA, LY, ZA_ - G_L), D_GS).cut(box(XA - 4, XA + 4, LY - T_LH / 2 - 0.4, LY + T_LH / 2 + 0.4, ZA_ - 3.0, ZA_ + 4))
    gg2 = gg2.cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6)).cut(stab((XA, LY, ZA_ - 3.5), (XA, LY, ZA_ - G_L - 1), 3.0))
    dw2 = math.degrees(math.atan2(LT_Y - LY, -(Pb2.z - A2.z)))
    gg2 = gg2.translate(A2 - A0).rotate(A2, A2 + cq.Vector(1, 0, 0), dw2)
    taste2 = taste.translate(cq.Vector(0, 0, HUB_LT))
    dv2 = (Pb2 - A2).normalized(); E2 = A2 + dv2 * 3.6
    neu = {"Lufttaste": taste2, "Lufttaste_Gabelgelenk": gg2, "Lufttaste_Schubstange": stab(E2.toTuple(), Pb2.toTuple(), 3.0).cut(taste2)}
    lg = [(n_, neu[n_] if n_ in neu else (rot(s_) if n_ in ("Lufthebel", "Lufttaste_Gelenkbolzen") or n_.startswith("Luftklappe") else s_)) for n_, s_ in lg]
    print(f"Lufttaste GEDRUECKT: Hebel {PHI:.2f} Grad, Stangenlaengenfehler {best[0]:.3f}, Klappe offen")
print(f"Luftklappe X {x0_:.1f}-{x1_:.1f} Y {y0_:.1f}-{y1_:.1f} | Lufttaste {LT_L:.1f} x {LT_B:.1f} | Lufthebel Lager X {XP:.1f} Z {ZP:.1f}, Auge X {XA:.1f} Z {ZA_:.1f}, Stange {dw:.1f} Grad")
# Luftgruppe (wie Standardmodell aufgebaut) an X = L/2 spiegeln: Einbau seitenverkehrt zur Zeichnung
spg = lambda s_: s_.mirror("YZ", basePointVector=(L / 2, 0, 0))
LT_XM = LT_X
def langloch_m(l_, b_, z0_, h_):
    return cq.Workplane("XY").center(LT_XM, LT_Y).slot2D(l_, b_).extrude(h_).translate((0, 0, z0_)).val()
deckel_r = deckel_r.cut(langloch_m(LT_L, LT_B, -1, T_DECKEL_BASS + 2))
schr_m = [s_ for n_, s_ in lg if n_.startswith("Lufthebel_Stuetze_Schraube")]
teile = [(n, (lambda t: [t := t.cut(sc) for sc in schr_m][-1])(t) if n == "Abdeckplatte_vorne" else t) for n, t in teile]
# Gedrueckte Taste: Hebelende und Gabelgelenk heben 10 mm und reichen in die Abdeckplatte -> Schlitz mit aufgesetzter Haube (dicht)
_abd = [t for n, t in teile if n == "Abdeckplatte_vorne"][0]
_bew = [s_ for n_, s_ in lg if n_ in ("Lufthebel", "Lufttaste_Gabelgelenk", "Lufttaste_Gelenkbolzen", "Lufttaste_Schubstange")]
_ib = [s_.intersect(box(0, L, 0, B, Z_ABD - 0.5, Z_ABD + 30)).BoundingBox() for s_ in _bew if s_.BoundingBox().zmax > Z_ABD - 0.5]
_ib = [b_ for b_ in _ib if b_.xlen > 0.01]
if _ib:
    hx0, hx1 = min(b_.xmin for b_ in _ib) - 1.0, max(b_.xmax for b_ in _ib) + 1.0
    hy0, hy1 = min(b_.ymin for b_ in _ib) - 1.0, max(b_.ymax for b_ in _ib) + 1.0
    hz1 = max(b_.zmax for b_ in _ib) + 1.0
    T_HB = 1.5
    teile = [(n, t.cut(box(hx0, hx1, hy0, hy1, Z_ABD - 1, Z_ABD + 5)) if n == "Abdeckplatte_vorne" else t) for n, t in teile]
    haube = box(hx0 - T_HB, hx1 + T_HB, hy0 - T_HB, hy1 + T_HB, Z_ABD + 2.5, max(hz1, Z_ABD + 2.5 + 1) + T_HB).cut(box(hx0, hx1, hy0, hy1, Z_ABD, max(hz1, Z_ABD + 2.5 + 1)))
    teile.append(("Abdeckplatte_Haube_Lufthebel", haube))
    print(f"Haube Lufthebel X {hx0:.1f}-{hx1:.1f} Y {hy0:.1f}-{hy1:.1f} bis Z {hz1 + T_HB:.1f}")
teile += lg
# ---------- Bassbodenleder 0,5 (wie Markus) zwischen den Rosetten, hinteres Drittel frei, in die Lufttaste hineingefuehrt ----------
T_LEDER_D, R_LEDER, A_SCHALL = 0.5, 10.0, 4.0
ly0, ly1 = D + (B - 2 * D) / 3.0, FI
lx0 = max(x_ + r_ for x_, _, r_ in _ros if x_ < L / 2) + A_SCHALL
lx1 = min(x_ - r_ for x_, _, r_ in _ros if x_ > L / 2) - A_SCHALL
feld = (cq.Workplane("XY").center((lx0 + lx1) / 2, (ly0 + ly1 + 2 * R_LEDER) / 2).rect(lx1 - lx0, ly1 - ly0 + 2 * R_LEDER)
        .extrude(T_DECKEL_BASS).edges("|Z").fillet(R_LEDER - 0.01).val()).intersect(box(lx0, lx1, ly0, ly1, 0, T_DECKEL_BASS))
leder_fl = feld.intersect(box(lx0, lx1, ly0, ly1, 0, T_LEDER_D))
leder_sl = langloch_m(LT_L + 2 * T_LEDER_D, LT_B + 2 * T_LEDER_D, 0, T_DECKEL_BASS).cut(langloch_m(LT_L, LT_B, -1, T_DECKEL_BASS + 2))
leder = leder_fl.cut(langloch_m(LT_L, LT_B, -1, T_DECKEL_BASS + 2)).fuse(leder_sl).clean()
deckel_r = deckel_r.cut(leder)
teile.append(("Bassbodenleder", leder))
# ---------- Rueckteil: Aussenkante (Deckelseite) zwischen den Schalloeffnungen gerundet (wie Markus) ----------
R_KANTE, L_AUSLAUF = 20.0, 30.0
def kantenprofil(x_, r_):
    r_ = max(r_, 0.05)
    w_ = cq.Workplane("YZ", origin=(x_, 0, 0)).moveTo(-0.01, -0.01).lineTo(r_, -0.01).threePointArc((r_ * (1 - 0.7071), r_ * (1 - 0.7071)), (-0.01, r_)).close()
    return w_.wires().val()
stationen = [(lx0 - L_AUSLAUF, 0.0), (lx0, R_KANTE), (lx1, R_KANTE), (lx1 + L_AUSLAUF, 0.0)]
kante = None
for (xa, ra), (xb, rb) in zip(stationen[:-1], stationen[1:]):
    seg_ = cq.Solid.makeLoft([kantenprofil(xa, ra), kantenprofil(xb, rb)], True)
    kante = seg_ if kante is None else kante.fuse(seg_)
rueck = rueck.cut(kante); deckel_h = deckel_h.cut(kante)
for k_ in range(len(teile)):
    if teile[k_][1].intersect(kante).Volume() > 0.05: teile[k_] = (teile[k_][0], teile[k_][1].cut(kante))
teile.append(("Bassdeckel_Leiste_hinten", deckel_h)); teile.append(("Bassdeckel", deckel_r))
_yfuge = D + B_DL
teile.append(("Bassdeckel_Auflage_hinten", box(xi0 + BLK + 0.5, xi1 - BLK - 0.5, _yfuge - 3.0, _yfuge + 3.0, T_DECKEL_BASS, T_DECKEL_BASS + 2.0)))
# ---------- Fuesse (fus2.dxf-Profil wie Markus) unter den Eckbloecken ----------
L_FUSS = 21.5
aussen = [(0, 0), (0, 7.0), (-0.5, 7.5), (-3.0, 7.5)]
_c, _r = (-8.889409, 11.283499), 7.0
for i in range(1, 30):
    t = math.radians(327.0 - (327.0 - 223.0) * i / 30)
    aussen.append((_c[0] + _r * math.cos(t), _c[1] + _r * math.sin(t)))
aussen += [(-14.0, 6.5), (-16.0, 6.5), (-16.0, 7.5), (-21.5, 7.5), (-21.5, 0)]
innen = [(0, 0), (0, 6.75), (-2.5, 6.75), (-2.5, 4.0), (-4.0, 4.0), (-4.0, 4.75), (-5.5, 4.75), (-5.5, 4.0), (-7.5, 2.0), (-21.5, 2.0), (-21.5, 0)]
def drehteil(profil):
    pts = [(r_, -(L_FUSS + x_)) for x_, r_ in profil]
    return cq.Workplane("XZ").polyline(pts).close().revolve(360, (0, 0, 0), (0, 1, 0)).val()
fuss_roh = drehteil(aussen).cut(drehteil(innen))
fuesse = []
for k, (bx, by) in enumerate(ecken):
    cx, cy = bx + BLK / 2, by + BLK / 2
    fuesse.append((f"Fuss_{k+1}", fuss_roh.translate(cq.Vector(cx, cy, 0))))
    z_k0, z_k1 = -(L_FUSS - 5.5), -(L_FUSS - 7.5)
    kopf = cq.Solid.makeCone(4.0, 1.5, z_k1 - z_k0, cq.Vector(cx, cy, z_k0), cq.Vector(0, 0, 1))
    sc = cq.Workplane("XY").center(cx, cy).circle(D_BL / 2).extrude(L_FUSS - 7.5 + T_DECKEL_BASS + 12).translate((0, 0, z_k1)).val().fuse(kopf)
    kuppe = [(0.4 + 0.6 * (1 - (r_ / 6.25) ** 2), r_) for r_ in [6.25 - 6.25 * i / 12 for i in range(13)]]
    gp = drehteil([(-5.5, 0), (-5.5, 4.75), (-4.0, 4.75), (-4.0, 4.0), (-2.5, 4.0), (-2.5, 6.75), (0, 6.75), (0, 6.25)] + kuppe)
    fuesse.append((f"Fuss_Gummipfropfen_{k+1}", gp.translate(cq.Vector(cx, cy, 0)))); fuesse.append((f"Fuss_Schraube_{k+1}", sc))
teile += fuesse
# ---------- Gaze innen ueber Rosetten (Deckel) und F-Loechern (Rueckteil), wie Markus ----------
T_GZ, R_GZ, RASTER, D_FADEN = 0.3, 3.0, 4.0, 0.4
from shapely.geometry import Polygon as _Poly, Point as _Pt
from gaze2d import gaze_platte
_huellen = []
for j_, (mx_, my_, _) in enumerate(sorted([r_ for r_ in _ros if r_[2] > 6.5])):
    rr_ = max(math.hypot(x_ - mx_, y_ - my_) + r_ for x_, y_, r_ in _ros if math.hypot(x_ - mx_, y_ - my_) < 30) + R_GZ
    _huellen.append((f"Gaze_Rosette_{j_+1}", cq.Solid.makeCylinder(rr_, T_GZ, cq.Vector(mx_, my_, T_DECKEL_BASS)), (_Pt(mx_, my_).buffer(rr_, 32), "XY", T_DECKEL_BASS)))
_jf = 0
for e in RD.query('SPLINE[layer=="LAYER_1"]'):
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.2)]
    if min(b for _, b in pp_) > -4.24:
        pg_ = _Poly([(a + 382.0, b + 4.24) for a, b in pp_]).buffer(R_GZ).simplify(0.05); _jf += 1
        _huellen.append((f"Gaze_FLoch_{_jf}", cq.Workplane("XZ").polyline(list(pg_.exterior.coords)[:-1]).close().extrude(-T_GZ).translate((0, D, 0)).val(), (pg_, "XZ", D)))
_jh = 0
for c in RD.query('CIRCLE[layer=="Hörner-Bohrung"]'):                # Gaze innen hinter den Hoernern (wie Markus)
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        _jh += 1
        _huellen.append((f"Gaze_Horn_{_jh}", cq.Solid.makeCylinder(r + R_GZ, T_GZ, cq.Vector(-x, FI - T_GZ, y + 4.24), cq.Vector(0, 1, 0)), (_Pt(-x, y + 4.24).buffer(r + R_GZ, 32), "XZ", FI - T_GZ)))
_aus = set()
for k_ in range(len(teile)):
    n_, s_ = teile[k_]
    if n_ in ("Bassdeckel",) or n_.startswith("Fuss") or n_.startswith("Horn_"): continue
    for gn_, gs_, _ in _huellen:
        if (s_.BoundingBox().zmin < 12 or gn_.startswith("Gaze_Horn")) and s_.intersect(gs_).Volume() > 1e-3: s_ = s_.cut(gs_); _aus.add(n_)
    teile[k_] = (n_, s_)
for gn_, gs_, _ in _huellen:
    if boden.intersect(gs_).Volume() > 1e-3: boden = boden.cut(gs_); _aus.add("Bassboden")
teile += [(gn_, gaze_platte(g2_[0], g2_[1], g2_[2], T_GZ, RASTER, D_FADEN)[0]) for gn_, gs_, g2_ in _huellen]
print("Gaze:", len(_huellen), "| Horn-Gaze", _jh, "| Aussparungen in", sorted(_aus))
for n_, s_ in teile:
    if n_.startswith("Horn_"): print(n_, "ymin", round(s_.BoundingBox().ymin, 2), "FI", FI)
# Auflagenuten fuer das Bodenstueck mit der Luftklappe (Abdeckplatte): Kanal-Laengsleiste, Front-Leiste, Abschlusswaende
_abd = [t for n, t in teile if n == "Abdeckplatte_vorne"][0]
_nut = box(59.2 - 1.5, 324.2 + 1.5, 153.0 - 1.5, FI - 5 + 1.5, Z_ABD, Z_ABD + 2.5)
for k_ in range(len(teile)):
    n_, t_ = teile[k_]
    if n_ in ("Aufdoppelung_Front_oben", "Abschlusswand_Vertiefung_links", "Abschlusswand_Vertiefung_rechts", "HB_Laengswand_Kanal_vorne"):
        teile[k_] = (n_, t_.cut(_nut))
# Innenleisten 2 mm an den schmalen Seitenteilen (neben den Begleiterstimmstoecken), vom Boden bis zur Hoehe der Front-Leiste
teile.append(("Innenleiste_Seite_links", box(D, D + 2.0, D + 5, FI - 5, Z_BO, 82)))
teile.append(("Innenleiste_Seite_rechts", box(L - D - 2.0, L - D, D + 5, FI - 5, Z_BO, 82)))
# ---------- Ausgabe mit Farben wie Markus/Standard ----------
teile = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + klappen
FARBEN = [("Gaze", (0.78, 0.80, 0.82)), ("Fuss_Gummi", (0.08, 0.08, 0.08)), ("Fuss_Schraube", (0.6, 0.6, 0.62)), ("Fuss", (0.28, 0.20, 0.14)),
          ("Frontblende", (0.70, 0.52, 0.30)), ("Front_Aufdoppelung", (0.78, 0.60, 0.38)), ("Front", (0.78, 0.60, 0.38)), ("Horn", (0.05, 0.05, 0.05)), ("Knopfzylinder", (0.95, 0.95, 0.92)), ("Knopf", (0.15, 0.15, 0.15)), ("Filzeinlage", (0.55, 0.55, 0.20)), ("Rueck", (0.74, 0.56, 0.36)), ("Seiten", (0.70, 0.52, 0.33)), ("Aufdoppelung", (0.78, 0.60, 0.38)),
          ("Bassbodenleder", (0.12, 0.08, 0.05)), ("Bassboden", (0.90, 0.82, 0.62)), ("Bassdeckel_Leiste", (0.30, 0.22, 0.15)), ("Bassdeckel_Auflage", (0.60, 0.45, 0.28)),
          ("Bassdeckel", (0.30, 0.22, 0.15)), ("Deckelauflage", (0.60, 0.45, 0.28)), ("Eckblock", (0.55, 0.40, 0.25)), ("Abdeckplatte", (0.55, 0.40, 0.25)), ("Abschlusswand", (0.55, 0.40, 0.25)), ("Innenleiste", (0.55, 0.40, 0.25)),
          ("Klappe_", None), ("HB_Folie", (0.85, 0.30, 0.30)), ("HB_Stimmplatte", (0.80, 0.82, 0.86)), ("HB_Trennwand", (0.50, 0.62, 0.80)), ("HB_", (0.65, 0.75, 0.88)),
          ("Begleiter_Stimmplatte", (0.80, 0.82, 0.86)), ("Begleiter_Sole", (0.55, 0.68, 0.45)), ("Begleiter_Schraube", (0.40, 0.40, 0.42)), ("Begleiter", (0.62, 0.72, 0.55)), ("Haken", (0.40, 0.40, 0.42)), ("Formteil_Schraube", (0.40, 0.40, 0.42)), ("Formteil", (0.60, 0.45, 0.28)), ("Dichtung", (0.12, 0.08, 0.05)), ("Lufttaste_Filz", (0.85, 0.15, 0.15)), ("Lufttaste_Gabelgelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)),
          ("Lufttaste_Schubstange", (0.10, 0.60, 0.60)), ("Lufttaste", (0.35, 0.25, 0.15)), ("Lufthebel_Feder", (0.75, 0.10, 0.75)), ("Lufthebel_Stuetze_Schraube", (0.40, 0.40, 0.42)),
          ("Lufthebel_Stuetze", (0.25, 0.75, 0.35)), ("Lufthebel_Bolzen", (0.95, 0.85, 0.20)), ("Lufthebel", (0.70, 0.72, 0.78)), ("Luftklappe_Leder", (0.45, 0.30, 0.20)),
          ("Luftklappe_Filz", (0.85, 0.85, 0.80)), ("Luftklappe", (0.20, 0.20, 0.22))]
def farbe(n):
    if n.startswith(("Klappe_", "Begleiterklappe_")): return (0.45, 0.30, 0.20) if "Leder" in n else (0.85, 0.85, 0.80) if "Filz" in n else (0.20, 0.20, 0.22)
    for p_, c_ in FARBEN:
        if n.startswith(p_) and c_: return c_
    print("FARBE ohne Zuordnung:", n); return (0.6, 0.6, 0.6)
exec(open("mech38.py").read())                # Mechanik (vor der Spiegelung)
a = cq.Assembly(name="Bassteil_Marie38")
SPIEGELN = True                                # wie Markus: links/rechts getauscht (X = L/2)
for n, s_ in teile:
    if SPIEGELN: s_ = s_.mirror("YZ", basePointVector=(L / 2, 0, 0))
    a.add(s_, name=n, color=cq.Color(*farbe(n)))
a.save(OUT)
print("gespeichert", OUT, len(teile), "Teile")
