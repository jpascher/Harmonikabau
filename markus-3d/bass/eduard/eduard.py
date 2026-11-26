# Bassteil Eduard Winkelbass (Zeichnungen Eduard-winkelbass-Mechanik.dxf, Eduard-WB_schnitt.dxf), Gehaeuse wie Markus:
# Rahmen mit Zinken, Deckel mit Rosetten, hinterer Leiste, Leder, Eckbloecken und Fuessen, Deckelauflage, Gaze von Markus.
# Zeichnung 380 x 205 -> Rohmass 382 x 207 (Markus): Zeichnungskoordinaten + 1 mm (OX, OY).
# Achsen: X = Laenge, Y = Tiefe (0 = Rueckteil aussen, B = Front aussen), Z ab Deckel aussen Richtung Balg.
# Aufruf: python3 eduard.py ausgabe.step [Knopfteilung]
import sys, math, ezdxf, ezdxf.path, cadquery as cq
from ezdxf import bbox
OUT = sys.argv[1]
T_KNOPF = float(sys.argv[2]) if len(sys.argv) > 2 else 20.0   # Knopfteilung (Variable, Zeichnung 20)
D, L, B = 5.5, 382.0, 207.0
OX = OY = 1.0
T_DECKEL_BASS = 4.0
Z_BU, Z_BO = 35.0, 38.5                       # Boden 3,5 (Schnitt: 1378,15..1381,65 ab Deckel aussen 1343,15)
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
XP0, YP0 = 721.8, 383.7                       # Draufsicht (Mechanik.dxf) -> Modell: X = x - XP0, Y = y - YP0 (Rohmass +1)
XS0, YS0 = 1343.15, 430.16                     # Schnitt -> Modell: Z = x - XS0, Y = y - YS0
# ---------- Front wie Markus: Knopfreihen-Abstufung (Reihe 1 vertieft, Aufdoppelung 2 + 5 mm), Hoerner, Blende, Knoepfe ----------
RD = ezdxf.readfile("rahmen.dxf").modelspace()
# Knopfreihen wie Standardmodell: Reihe 1 (unten, 9 Knoepfe, Bohrung Ø11 mit Filz) buendig an der Front, Reihe 2 (10 Knoepfe,
# Ø10) ueber eine 7-mm-Blende erhoeht; Hoehen und Teilung nach der Marie-Zeichnung (Teilung variabel), Mitte X 198.
X_KM = L / 2                                   # Knopfreihen mittig am Korpus
REIHEN_N = {1: 11, 2: 10}                       # Eduard: 21 Knoepfe (Schubstangen X 91..291 im 10er-Raster)
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
    teile.append((f"Knopf_R{rr}_{X_:.0f}", stab((X_, yk, Z_), (X_, yk + H_KOPF, Z_), D_KOPF)))
    teile.append((f"Knopfzylinder_R{rr}_{X_:.0f}", stab((X_, yk - lz, Z_), (X_, yk, Z_), D_ZYL)))
    teile.append((f"Filzeinlage_R{rr}_{X_:.0f}", stab((X_, ya - dl_, Z_), (X_, ya, Z_), 2 * r_).cut(stab((X_, ya - dl_ - 1, Z_), (X_, ya + 1, Z_), D_ZYL))))
print("Front wie Standardmodell (Reihe 2 mit 7-mm-Blende), Teilung", T_KNOPF, "Knoepfe", len(kn), "X", min(k[0] for k in kn), "-", max(k[0] for k in kn), "| Hoerner", len(aus_front), horn_pos)
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
# ======================= Eduard Winkelbass: Einbauten =======================
P = lambda x_, y_: (x_ - XP0, y_ - YP0)                 # Draufsicht -> Modell (vor der Spiegelung)
MD = ezdxf.readfile("eduard_mechanik.dxf").modelspace()
def rect(e):
    p_ = [P(q[0], q[1]) for q in e.get_points()]; xs_ = [q[0] for q in p_]; ys_ = [q[1] for q in p_]
    return (min(xs_), max(xs_), min(ys_), max(ys_))
# ---------- Boden 3,5 (Z 35..38,5) ueber die ganze Innenflaeche; Oeffnungen und Klappen aus der Draufsicht ----------
boden = box(D, L - D, D, FI, Z_BU, Z_BO)
oeff = sorted([rect(e) for e in MD.query('LWPOLYLINE[layer=="Füllung-DH-Öffnungen2"]')], key=lambda r_: (-r_[2], r_[0]))
LUFT_OEFF = [r_ for r_ in oeff if r_[1] - r_[0] > 60][0]
hb_oeff = [r_ for r_ in oeff if r_ is not LUFT_OEFF]
hb_oeff = [(x0, x1, y0 - 0.3, y1 - 0.3) if y1 < 110 else (x0, x1, y0, y1) for (x0, x1, y0, y1) in hb_oeff]   # hintere Reihe an die Trennwand (Zeichnung 0,2 ungenau)
# Klappen wie Markus: Vertiefung auf der Mechanikseite (Klappe + 0,5 ringsum, 2 tief); Klappe nach Zeichnung (Oeffnung + 3 ringsum)
RAND_K, SPIEL, TIEF_V = 3.0, 0.5, 2.0
klappen = []
def klappe_rect(nm, x0, x1, y0, y1, rand=RAND_K, durch=True):
    global boden                                   # durch=False: nur Vertiefung (Oeffnungen kommen aus dem Stimmstock)
    if durch: boden = boden.cut(box(x0, x1, y0, y1, Z_BU - 1, Z_BO + 1))
    boden = boden.cut(box(x0 - rand - SPIEL, x1 + rand + SPIEL, y0 - rand - SPIEL, y1 + rand + SPIEL, Z_BU - 1, Z_BU + TIEF_V))
    z_ = Z_BU + TIEF_V
    for t_, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
        klappen.append((f"{nm}_{t_}", box(x0 - rand, x1 + rand, y0 - rand, y1 + rand, z_ - dk, z_))); z_ -= dk
for i, (x0, x1, y0, y1) in enumerate(hb_oeff): klappe_rect(f"Klappe_{i+1}", x0, x1, y0, y1)
# ---------- Helikon-/Grundbasskammern (Winkelbass): 5 Stellungen x 2 Kammern, Teilung 54, Trennwaende 3,5 ----------
# Profil aus dem Schnitt (Y, Z): vordere Kammer mit stehender Platte an der Vorderseite, hintere Kammer mit liegender Platte
# oben, die genau im rechten Winkel zur stehenden steht (Zeichnung 2 Grad schief); Platten 94 x 53 wie Markus breit.
ZB = Z_BO
c_x = sorted({round((r_[0] + r_[1]) / 2, 2) for r_ in hb_oeff})
T_TW, TEIL = 3.5, 54.0
Y_PV = 158.05                                   # Innenseite der stehenden Platte (Schnitt 588,21)
Z_PO = 155.5                                    # Auflage der liegenden Platte (Schnitt 1498,65)
PL_L, PL_B, PL_T = 94.0, 53.0, 3.5
# Aufbau wie gezeichnet (Schnitt, Werte Z/Y): Langteile ueber alle Kammern, an den Knicken stumpf gestossen;
# Vierkantleisten unten/oben an der stehenden Platte und hinten unter der liegenden Platte; Trennwaende (vorne/hinten)
# in die Langteile und Vierkantleisten eingelassen (Nuten = Ueberdeckung mit der Trennwand); Trennfeder durchgehend.
y_po0 = (57.8 + 147.8) / 2 - PL_L / 2                           # liegende Platte mittig ueber der Oeffnung 57,8..147,8
def pzy(x0, x1, zy): return py(x0, x1, [(y_, z_) for z_, y_ in zy])          # Profil als (Z, Y) wie im Schnitt
XD = [c_x[0] - TEIL / 2 + TEIL * k for k in range(len(c_x) + 1)]            # Trennwandmitten
LX0, LX1 = XD[0] - T_TW / 2, XD[-1] + T_TW / 2
Z_VLU, Z_VLO = 49.0, 139.0                                                   # Plattenrand 2 auf den Vierkantleisten
lang = {
    "HB_Langteil_hinten_1": pzy(LX0, LX1, [(ZB, 68.0), (100.0, 68.0), (100.0, 72.0), (ZB, 72.0)]),
    "HB_Langteil_hinten_2": pzy(LX0, LX1, [(100.0, 68.0), (126.73, 51.03), (127.46, 54.96), (100.0, 72.0)]),
    "HB_Langteil_hinten_3": pzy(LX0, LX1, [(126.73, 51.0), (145.5, 51.0), (145.5, 55.0), (127.46, 55.0)]),
    "HB_Langteil_mitte_1": pzy(LX0, LX1, [(ZB, 103.43), (126.5, 103.43), (126.5, 107.43), (ZB, 107.43)]),
}
_d = (24.44 / 55.31, 49.62 / 55.31); _n = (-_d[1] * 4.0, _d[0] * 4.0)
_A, _Bp = (126.5, 103.43), (126.5 + (165.0 - 103.43) * 24.44 / 49.62, 165.0)
vl = {
    "HB_Vierkantleiste_unten": pzy(LX0, LX1, [(ZB, 148.05), (Z_VLU, 148.05), (Z_VLU, Y_PV), (ZB, Y_PV)]),
    "HB_Vierkantleiste_oben": pzy(LX0, LX1, [(Z_VLO, 147.8), (Z_PO, 147.8), (Z_PO, Y_PV), (Z_VLO, Y_PV)]),
    "HB_Vierkantleiste_hinten": pzy(LX0, LX1, [(145.5, 47.8), (Z_PO, 47.8), (Z_PO, 57.8), (145.5, 57.8)]),
}
m5 = pzy(LX0, LX1, [_A, _Bp, (_Bp[0] + _n[0], _Bp[1] + _n[1]), (_A[0] + _n[0], _A[1] + _n[1])])
m5 = m5.cut(lang["HB_Langteil_mitte_1"]).cut(vl["HB_Vierkantleiste_oben"]).cut(box(0, L, 0, B, Z_PO, 300)).cut(box(0, L, Y_PV, B, 0, 300))
lang["HB_Langteil_mitte_2"] = m5
TW_H = [(ZB, 70.0), (100.0, 70.0), (127.1, 53.0), (Z_PO, 53.0), (Z_PO, 153.05), (150.94, 153.05), (126.5, 103.43), (ZB, 103.43)]
TW_V = [(ZB, 105.53), (125.21, 105.53), (150.31, 154.84), (150.31, Y_PV), (43.5, Y_PV), (43.5, 148.05), (ZB, 148.05)]
feder = pzy(XD[0], XD[-1], [(99.93, 128.23), (58.5, 123.14), (47.23, 149.86), (45.33, 148.95), (57.15, 120.93), (100.18, 126.25)])
feder = feder.cut(vl["HB_Vierkantleiste_unten"])
tw = []
for k, xd in enumerate(XD):
    tw.append((f"HB_Trennwand_hinten_{k+1}", pzy(xd - T_TW / 2, xd + T_TW / 2, TW_H)))
    tw.append((f"HB_Trennwand_vorne_{k+1}", pzy(xd - T_TW / 2, xd + T_TW / 2, TW_V).cut(feder)))   # Feder durch alle Kammern
_tws = cq.Compound.makeCompound([s_ for _, s_ in tw])
for n_, s_ in list(lang.items()) + list(vl.items()): teile.append((n_, s_.cut(_tws)))         # Nuten fuer die Trennwaende
teile += tw + [("HB_Trennfeder", feder)]
# Helikon-Stimmplatten wie Markus (Bassboden2.dxf, Kammer 53 x 94): keilfoermig 5 / 16, Rand 2, Koerper je Seite 2 schmaeler,
# 4 Zungenschlitze (7,5 x 51 / 9,2 x 71,5 / 9,2 x 71,5 / 7,5 x 51), buendig am dicken Ende.
# Dickes Ende = Schraubenende (stehend unten, liegend hinten), duennes Ende mit Haken am Knick (oben).
T_DUENN, T_DICK, T_RAND, RUECK = 5.0, 16.0, 2.0, 2.0
SCHLITZE = [(5.1, 12.6, 37.5, 88.5), (15.9, 25.1, 17.0, 88.5), (27.9, 37.1, 17.0, 88.5), (40.4, 47.9, 37.5, 88.5)]   # quer, laengs ab duennem Ende
def platte_lokal():
    """lokal: u quer 0..53 (X), w laengs 0..94 (duennes Ende 0, dickes 94), t: 0 = Auflage, + aussen, - in die Kammer"""
    pl = box(0, PL_B, 0, PL_L, 0, T_RAND)
    prof = [(RUECK, 0.0), (PL_L - RUECK, 0.0), (PL_L - RUECK, T_RAND - T_DICK), (RUECK, T_RAND - T_DUENN)]
    keil = cq.Workplane("YZ").polyline(prof).close().extrude(PL_B - 2 * RUECK).val().translate(cq.Vector(RUECK, 0, 0))
    pl = pl.fuse(keil)
    for (u0, u1, w0, w1) in SCHLITZE: pl = pl.cut(box(u0, u1, w0, w1, -T_DICK - 1, T_RAND + 1))
    return pl
PL0 = platte_lokal()
def stehend(xc):      # Auflage Y = Y_PV, aussen +Y; duennes Ende oben (Z 141), dickes unten (Z 47)
    return PL0.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -90).translate(cq.Vector(xc - PL_B / 2, Y_PV, 47.0 + PL_L))
def liegend(xc):      # Auflage Z = Z_PO, aussen +Z; duennes Ende vorne (Y 149,8), dickes hinten (Y 55,8)
    return PL0.mirror("XZ").translate(cq.Vector(xc - PL_B / 2, y_po0 + PL_L, Z_PO))
D_HAKEN, D_KOPF, H_KOPF, D_SCHAFT = 1.8, 5.0, 1.5, 2.5
befest = []
for k, xc in enumerate(c_x):
    teile.append((f"HB_Stimmplatte_stehend_{k+1}", stehend(xc)))
    teile.append((f"HB_Stimmplatte_liegend_{k+1}", liegend(xc)))
    # stehend: Haken oben (Draht in die obere Vierkantleiste), Schraube unten (in die untere)
    yh, zh = Y_PV + T_RAND + D_HAKEN / 2, 47.0 + PL_L + 0.1 + D_HAKEN / 2
    befest.append((f"HB_Haken_stehend_{k+1}", stab((xc, yh, zh), (xc, Y_PV - 9.0, zh), D_HAKEN).fuse(stab((xc, yh, zh), (xc, yh, zh - 5.0), D_HAKEN)).fuse(cq.Solid.makeSphere(D_HAKEN / 2 + 0.01, cq.Vector(xc, yh, zh)))))
    zs = 47.0 - 1.5
    befest.append((f"HB_Schraube_stehend_{k+1}", stab((xc, Y_PV + T_RAND, zs), (xc, Y_PV + T_RAND + H_KOPF, zs), D_KOPF).fuse(stab((xc, Y_PV - 9.0, zs), (xc, Y_PV + T_RAND, zs), D_SCHAFT))))
    # liegend: Haken vorne am Knick (obere Vierkantleiste), Schraube hinten (hintere Vierkantleiste)
    yh, zh = y_po0 + PL_L + 0.1 + D_HAKEN / 2, Z_PO + T_RAND + D_HAKEN / 2
    befest.append((f"HB_Haken_liegend_{k+1}", stab((xc, yh, zh), (xc, yh, Z_PO - 12.0), D_HAKEN).fuse(stab((xc, yh, zh), (xc, yh - 5.0, zh), D_HAKEN)).fuse(cq.Solid.makeSphere(D_HAKEN / 2 + 0.01, cq.Vector(xc, yh, zh)))))
    ys = y_po0 - 1.5
    befest.append((f"HB_Schraube_liegend_{k+1}", stab((xc, ys, Z_PO + T_RAND), (xc, ys, Z_PO + T_RAND + H_KOPF), D_KOPF).fuse(stab((xc, ys, Z_PO - 9.0), (xc, ys, Z_PO + T_RAND), D_SCHAFT))))
_bf = cq.Compound.makeCompound([s_ for _, s_ in befest])
teile = [(n_, s_.cut(_bf) if n_.startswith("HB_Vierkantleiste") else s_) for n_, s_ in teile]  # Vorbohrungen
teile += befest
# ---------- Begleiterstimmstock: zweiseitig, 16 Kammern Teilung 17 (Trennwaende 2), Platten 42,5 x 16, Schieber ----------
BG_X = sorted(rect(e)[0] for e in MD.query('LWPOLYLINE[layer=="Def"]') if abs(rect(e)[1] - rect(e)[0] - 2) < 0.1 and abs(rect(e)[3] - rect(e)[2] - 24) < 0.1)
bx0, bx1 = BG_X[0], BG_X[-1] + 2.0
S_Z0, S_Z1, P_Z0, P_Z1, W_Z1, K_Z1 = ZB, ZB + 3.5, ZB + 3.5, ZB + 46.0, ZB + 44.0, ZB + 47.5   # Sole, Platten 42,5, Waende, Deckel
Y_SO0, Y_SO1, Y_PH0, Y_PH1, Y_PV0, Y_PV1, Y_MW0, Y_MW1 = 28.0, 63.0, 33.0, 36.5, 54.5, 58.0, 44.5, 46.5
sole = box(bx0, bx1, Y_SO0, Y_SO1, S_Z0, S_Z1)
bg = box(bx0, bx1, Y_MW0, Y_MW1, S_Z1, W_Z1).fuse(box(bx0, bx1, Y_PH1, Y_PV0, W_Z1, K_Z1))
for xd in BG_X: bg = bg.fuse(box(xd, xd + 2.0, Y_PH1, Y_PV0, S_Z1, W_Z1))
loecher = [P(c.dxf.center.x, c.dxf.center.y) + (c.dxf.radius,) for c in MD.query('CIRCLE[layer=="Füllung-DH-Öffnungen2"]')]
# Oeffnungen an die Kammern angepasst: Ø10 der Zeichnung reicht unter die Stimmplatten -> Rechteck zwischen Platte und
# Mittelwand (0,5 Abstand), 10 breit an der gezeichneten Stelle (innerhalb Kammer und Klappe)
for (x_, y_, r_) in loecher:
    ya_, yb_ = (Y_PH1 + 0.5, Y_MW0 - 0.5) if y_ < Y_MW0 else (Y_MW1 + 0.5, Y_PV0 - 0.5)
    sole = sole.cut(box(x_ - 5.0, x_ + 5.0, ya_, yb_, S_Z0 - 1, S_Z1 + 1))
# Boden: die vier Oeffnungen unter jeder Begleiterklappe zu einer vereint (Kreuzsteg entfernt)
_lx = sorted({round(x_, 2) for x_, _, _ in loecher})
for i_ in range(0, len(_lx), 2):
    boden = boden.cut(box(_lx[i_] - 5.0, _lx[i_ + 1] + 5.0, Y_PH1 + 0.5, Y_PV0 - 0.5, Z_BU - 1, Z_BO + 1))
# Schieber 1,5 in einer Nut an der Unterseite der Sole (Stellung "offen" wie gezeichnet)
sch = [rect(e) for e in MD.query('LWPOLYLINE[layer=="Schieber"]')]
sch_aussen = max(sch, key=lambda r_: r_[1] - r_[0])
x0, x1, y0, y1 = sch_aussen
sole = sole.cut(box(x0 - 0.3, x1 + 0.3, y0 - 0.3, y1 + 0.3, S_Z0 - 0.1, S_Z0 + 1.5 + 0.2))
schieber = box(x0, x1, y0, y1, S_Z0 + 0.05, S_Z0 + 1.55)
for r_ in sch:
    if r_ is sch_aussen or r_[1] - r_[0] < 5: continue
    schieber = schieber.cut(box(r_[0], r_[1], r_[2], r_[3], S_Z0 - 1, S_Z0 + 3))
teile += [("Begleiter_Sole", sole), ("Begleiter_Stimmstock", bg), ("Begleiter_Schieber", schieber)]
for j in range(len(BG_X) - 1):
    xc = (BG_X[j] + 2.0 + BG_X[j + 1]) / 2
    teile.append((f"Begleiter_Stimmplatte_h_{j+1}", box(xc - 8.0, xc + 8.0, Y_PH0, Y_PH1, P_Z0, P_Z1 - 1.0)))
    teile.append((f"Begleiter_Stimmplatte_v_{j+1}", box(xc - 8.0, xc + 8.0, Y_PV0, Y_PV1, P_Z0, P_Z1 - 1.0)))
# Begleiterklappen: 8 Quadrate 28 x 28 (je 4 Loecher), Vertiefung 0,5 groesser
for j, e in enumerate(sorted(MD.query('LWPOLYLINE[layer=="Klappen-rand-aussen"]'), key=lambda e: rect(e)[0])):
    x0, x1, y0, y1 = rect(e)
    klappe_rect(f"Begleiterklappe_{j+1}", x0, x1, y0, y1, rand=0.0, durch=False)   # Boden nur mit den 4 Oeffnungen des Stimmstocks
print("Begleiter: Kammern", len(BG_X) - 1, "X", round(bx0, 1), "-", round(bx1, 1), "| Loecher", len(loecher), "| Schieberfenster", len(sch) - 2)
# ---------- Lufttaste + Luftklappe + Lufthebel (zuerst, die Bassgestaenge werden spaeter daran vorbeigefuehrt) ----------
x0, x1, y0, y1 = LUFT_OEFF
klappe_rect("Luftklappe", x0, x1, y0, y1)
luft_kl = (x0 - RAND_K, x1 + RAND_K, y0 - RAND_K, y1 + RAND_K)
lh_r = [rect(e) for e in MD.query('LWPOLYLINE[layer=="FismenHebel"]')][0]
LT_X, LT_Y, LT_L, LT_B = (lh_r[0] + lh_r[1]) / 2, (lh_r[2] + lh_r[3]) / 2, lh_r[1] - lh_r[0], lh_r[3] - lh_r[2]
T_FILZ_LT, HUB_LT = 1.0, 10.0
def langloch_m(l_, b_, z0_, h_):
    return cq.Workplane("XY").center(LT_X, LT_Y).slot2D(l_, b_).extrude(h_).translate((0, 0, z0_)).val()
lt_filz = langloch_m(LT_L, LT_B, 0, T_DECKEL_BASS).cut(langloch_m(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT, -1, T_DECKEL_BASS + 2))
taste = cq.Workplane("XY").center(LT_X, LT_Y).slot2D(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT).extrude(T_DECKEL_BASS + HUB_LT).faces("<Z").edges().fillet(1.0).val()
taste = taste.translate(cq.Vector(0, 0, -HUB_LT))
deckel_r = deckel_r.cut(langloch_m(LT_L, LT_B, -1, T_DECKEL_BASS + 2))
teile += [("Lufttaste_Filz", lt_filz), ("Lufttaste", taste)]
# Lufthebel wie Marie/Standard (kb_lufthebel.dxf: Alu-Form, Gabelgelenk + Schubstange zur Taste, Lagerboecke mit Schrauben,
# zwei Schenkelfedern). An Eduard angepasst: Tastenarm und Klappenarm ausserhalb des Lagerbereichs (+-5) gestreckt,
# Auge ueber der Taste wie Marie, Drehpunkt so gelegt, dass die Klappe wie bei Marie oeffnet (Mitte ca. 13,7 bei 10 Hub)
# und das Klappenende gedrueckt frei vom Deckel bleibt. Lagerboecke am Boden (Unterseite Z_BU) angeschraubt.
import math as _m
def kugel(p, d): return cq.Solid.makeSphere(d / 2 + 0.01, cq.Vector(*p), angleDegrees1=-90, angleDegrees2=90)
def pfad(pts, d):
    s_ = None
    for a_, b_ in zip(pts[:-1], pts[1:]):
        seg = stab(a_, b_, d).fuse(kugel(b_, d)); s_ = seg if s_ is None else s_.fuse(seg)
    return s_.fuse(kugel(pts[0], d))
def pts_of(e, tol=0.1):
    p_ = [(v_.x, v_.y) for v_ in ezdxf.path.make_path(e).flattening(tol)]
    q_ = [p_[0]] + [b_ for a_, b_ in zip(p_[:-1], p_[1:]) if _m.dist(a_, b_) > 1e-6]
    if len(q_) > 2 and _m.dist(q_[0], q_[-1]) < 1e-6: q_ = q_[:-1]
    return q_
def platte_xz(pts, ya, yb): return cq.Workplane("XZ").polyline(pts).close().extrude(-(yb - ya)).translate((0, ya, 0)).val()
T_LEDER_K, T_KLK = 1.0, 5.5                                         # Luftklappe: Leder 1, Filz 4 + Epoxy 1,5 (klappe_rect)
DR_F, DI_F, N_F = 0.5, 2.8, 3
LH = ezdxf.readfile("kb_lufthebel.dxf").modelspace()
X_P0, X_A0 = -55.0, -12.65                                          # Lager und Auge in der Zeichnung
XA = LT_X + 5.5                                                     # Auge ueber der Taste wie Marie
KX = (luft_kl[0] + luft_kl[1]) / 2
XP = (1.37 * XA + KX) / 2.37                                 # Klappenmitte oeffnet 13,7 bei 10 Tastenhub (wie Marie)
K_A = (XA - XP - 5.0) / (X_A0 - X_P0 - 5.0)
K_K = (XP - KX - 5.0) / (X_P0 + 144.15 - 5.0)                       # Marie: Klappenmitte 89,15 links vom Lager
def mx(x_):
    d_ = x_ - X_P0
    if abs(d_) <= 5.0: return XP + d_
    return XP + (5.0 + (d_ - 5.0) * K_A) if d_ > 0 else XP - (5.0 + (-d_ - 5.0) * K_K)
z_kl_u = Z_BU + TIEF_V - T_LEDER_K - T_KLK                          # Unterseite Luftklappe 30,5
Z_LH = z_kl_u + 5.3                                                  # wie Marie: Hebeloberkante 1,2 unter der Klappe
T_LH = 3.0
LY = (luft_kl[2] + luft_kl[3]) / 2
DZ_A = 4.5                                                          # Tastenarm 4,5 tiefer gekroepft: gedrueckt bleibt das Auge unter dem Boden
def dz_a(x_): return DZ_A * min(1.0, max(0.0, (x_ - X_P0 - 5.0) / (X_A0 - X_P0 - 5.0)))
def lh_xz(p_): return [(mx(a_), Z_LH - b_ - dz_a(a_)) for a_, b_ in p_]
hp_ = lh_xz(pts_of([e for e in LH.query('LWPOLYLINE[layer=="2-Gelenk-makierung"]')][0], 0.05))
ZP = Z_LH - 5.03; ZA_ = Z_LH - 10.4 - DZ_A
lg = []
lg.append(("Luftklappe_Leder", box(luft_kl[0], luft_kl[1], luft_kl[2], luft_kl[3], z_kl_u + T_KLK, Z_BU + TIEF_V)))
teile = [(n_, s_) for n_, s_ in teile]
klappen = [(n_, s_) for n_, s_ in klappen if not n_.startswith("Luftklappe")]
lg.append(("Luftklappe_Filz", box(luft_kl[0], luft_kl[1], luft_kl[2], luft_kl[3], z_kl_u + 1.5, z_kl_u + T_KLK)))
lg.append(("Luftklappe", box(luft_kl[0], luft_kl[1], luft_kl[2], luft_kl[3], z_kl_u, z_kl_u + 1.5)))
hebel = platte_xz(hp_, LY - T_LH / 2, LY + T_LH / 2).fuse(stab((XP, LY - T_LH / 2, ZP), (XP, LY + T_LH / 2, ZP), 5.0))
hebel = hebel.cut(stab((XP, LY - 3, ZP), (XP, LY + 3, ZP), 2.05)).cut(stab((XA, LY - 3, ZA_), (XA, LY + 3, ZA_), 1.08))
lg.append(("Lufthebel", hebel))
L_FL = N_F * (DR_F + 0.15) + DR_F; W_LB = 8.0
y_f = [(LY - T_LH / 2 - 0.15 - L_FL, LY - T_LH / 2 - 0.15), (LY + T_LH / 2 + 0.15, LY + T_LH / 2 + 0.15 + L_FL)]
y_w = [(y_f[0][0] - 0.15 - W_LB, y_f[0][0] - 0.15), (y_f[1][1] + 0.15, y_f[1][1] + 0.15 + W_LB)]
bk = [p_ for p_ in [pts_of(e, 0.05) for e in LH.query('LWPOLYLINE[layer=="2-Gelenke"]')] if min(a_ for a_, _ in p_) < -40][0]
bk_xz = [(a_ - X_P0 + XP, Z_LH - b_) for a_, b_ in bk]               # Lagerbock unveraendert, nur verschoben
xb0, xb1 = min(a_ for a_, _ in bk_xz), max(a_ for a_, _ in bk_xz)
x_schr = -64.0 - X_P0 + XP
schr = []
for k_, (ya_, yb_) in enumerate(y_w):
    w_ = platte_xz(bk_xz, ya_, yb_).fuse(box(xb0, xb1, ya_, yb_, Z_LH - 0.5, Z_BU + 3.0))
    w_ = w_.cut(box(xb0 - 1, xb1 + 1, ya_ - 1, yb_ + 1, Z_BU, Z_BU + 10)).cut(stab((XP, LY - 30, ZP), (XP, LY + 30, ZP), 2.05))
    for dy_ in (-2.0, 2.0):
        yc_ = (ya_ + yb_) / 2 + dy_
        sc_ = stab((x_schr, yc_, Z_LH - 2.53 - 0.9), (x_schr, yc_, Z_LH - 2.53), 2.6).fuse(stab((x_schr, yc_, Z_LH - 2.53), (x_schr, yc_, Z_BU + 2.0), 1.4))
        w_ = w_.cut(stab((x_schr, yc_, Z_LH - 4), (x_schr, yc_, Z_BU + 1), 1.5)); schr.append(sc_)
    lg.append((f"Lufthebel_Lagerbock_{k_+1}", w_))
for j_, sc_ in enumerate(schr): lg.append((f"Lufthebel_Lagerbock_Schraube_{j_+1}", sc_))
lg.append(("Lufthebel_Bolzen", stab((XP, y_w[0][0] - 0.6, ZP), (XP, y_w[1][1] + 0.6, ZP), 2.0)))
def kante_oben(x_):
    zz_ = []
    for (xa_, za_), (xb_, zb_) in [(hp_[i], hp_[(i + 1) % len(hp_)]) for i in range(len(hp_))]:
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
    zb_ = Z_BU - DR_F / 2 - 0.05
    s1 = pfad([p1, (XP + rmf + 3.0, y_out - sg * DR_F / 2, zb_), (XP + rmf + 12.0, y_out - sg * DR_F / 2, zb_)], DR_F)
    lg.append((f"Lufthebel_Feder_{k_+1}", sp.fuse(haken).fuse(s1)))
G_L, D_GS = 10.0, 6.0
def gabel(A_, P_):                                                  # Gabelgelenk am Auge A_, Achse genau zur Stange (Richtung P_)
    g_ = stab((XA, LY, ZA_ + 3.0), (XA, LY, ZA_ - G_L), D_GS).cut(box(XA - 4, XA + 4, LY - T_LH / 2 - 0.4, LY + T_LH / 2 + 0.4, ZA_ - 3.0, ZA_ + 4))
    g_ = g_.cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6)).cut(stab((XA, LY, ZA_ - 3.5), (XA, LY, ZA_ - G_L - 1), 3.0))
    A0_ = cq.Vector(XA, LY, ZA_); d_ = (P_ - A_).normalized(); ziel = A_ + d_ * 10.0
    wx = _m.degrees(_m.atan2(abs(d_.y), abs(d_.z))); wy = _m.degrees(_m.atan2(abs(d_.x), _m.hypot(d_.y, d_.z)))
    best_ = None
    for sx in (-1, 1):
        for sy in (-1, 1):
            f_ = lambda s_: s_.translate(A_ - A0_).rotate(A_, A_ + cq.Vector(1, 0, 0), sx * wx).rotate(A_, A_ + cq.Vector(0, 1, 0), sy * wy)
            e_ = (f_(cq.Vertex.makeVertex(XA, LY, ZA_ - 10.0)).Center() - ziel).Length
            if best_ is None or e_ < best_[0]: best_ = (e_, f_)
    return best_[1](g_)
A0 = cq.Vector(XA, LY, ZA_)
Pb = cq.Vector(XA, LT_Y, T_DECKEL_BASS - HUB_LT - 0.5)              # Stange in die Taste eingesetzt
dw = _m.degrees(_m.atan2(LT_Y - LY, -(Pb.z - ZA_)))
lg.append(("Lufttaste_Gabelgelenk", gabel(A0, Pb))); lg.append(("Lufttaste_Gelenkbolzen", stab((XA, LY - D_GS / 2, ZA_), (XA, LY + D_GS / 2, ZA_), 1.0)))
E = A0 + (Pb - A0).normalized() * 3.6
lg.append(("Lufttaste_Schubstange", stab(E.toTuple(), Pb.toTuple(), 3.0).cut(taste)))
# gedrueckter Zustand (Taste eben mit dem Deckel): Hebel/Klappe um das Lager schwenken, Stangenlaenge gleich
Pb2 = Pb + cq.Vector(0, 0, HUB_LT); l_st = (Pb - A0).Length
def _A(ph_):
    c_, s_ = _m.cos(ph_), _m.sin(ph_); dx_, dz_ = XA - XP, ZA_ - ZP
    return cq.Vector(XP + c_ * dx_ - s_ * dz_, LY, ZP + s_ * dx_ + c_ * dz_)
best = min(((abs((Pb2 - _A(_m.radians(g_ / 50))).Length - l_st), g_ / 50) for g_ in range(-1500, 1501)), key=lambda t_: (t_[0] > 0.05, abs(t_[1]), t_[0]))
PHI = best[1]
rot = lambda s_: s_.rotate(cq.Vector(XP, LY, ZP), cq.Vector(XP, LY - 1, ZP), PHI)
if abs(rot(cq.Vertex.makeVertex(XA, LY, ZA_)).Center().z - _A(_m.radians(PHI)).z) > 0.01:
    rot = lambda s_: s_.rotate(cq.Vector(XP, LY, ZP), cq.Vector(XP, LY + 1, ZP), PHI)
A2 = _A(_m.radians(PHI))
taste2 = taste.translate(cq.Vector(0, 0, HUB_LT))
E2 = A2 + (Pb2 - A2).normalized() * 3.6
neu = {"Lufttaste": taste2, "Lufttaste_Gabelgelenk": gabel(A2, Pb2),
       "Lufttaste_Schubstange": stab(E2.toTuple(), Pb2.toTuple(), 3.0).cut(taste2)}
lg_g = [(n_ + "_gedrueckt", neu[n_] if n_ in neu else rot(s_)) for n_, s_ in lg + [("Lufttaste", taste)]
        if n_ in neu or n_ in ("Lufthebel", "Lufttaste_Gelenkbolzen") or n_.startswith("Luftklappe")]
_zmin = min(s_.BoundingBox().zmin for n_, s_ in lg_g if n_.startswith(("Lufthebel", "Luftklappe")))
_kl2 = [s_ for n_, s_ in lg_g if n_.startswith("Luftklappe")]
print(f"Lufthebel: Lager X {XP:.1f} Z {ZP:.1f}, Auge X {XA:.1f} Z {ZA_:.1f}, Klappenmitte X {KX:.1f} | gedrueckt {PHI:.2f} Grad, "
      f"Stangenfehler {best[0]:.3f}, tiefster Punkt Hebel/Klappe Z {_zmin:.1f} (Deckel innen {T_DECKEL_BASS:.1f})")
# Eckblock vorne links: Ausnehmung fuer die Luftklappe samt Oeffnungsweg (gedrueckt + 1)
_ek = box(luft_kl[0] - 1.0, luft_kl[1] + 1.0, luft_kl[2] - 1.0, luft_kl[3] + 1.0, min(b_.BoundingBox().zmin for b_ in _kl2) - 1.0, Z_BU + 0.01)
teile = [(n_, s_.cut(_ek) if n_.startswith("Eckblock") and s_.intersect(_ek).Volume() > 0.01 else s_) for n_, s_ in teile]
boden = boden.cut(cq.Compound.makeCompound(schr))
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
# ---------- Ecke vorne links (bei der Luftklappe): statt Eckblock ein L-Winkel (Alu 3) an der Seitenwand, unten frei ----------
_bx, _by = ecken[2]; _cx, _cy = _bx + BLK / 2, _by + BLK / 2; T_W, H_W = 3.0, 22.0
_w = box(_bx, _bx + BLK, _by, _by + BLK, T_DECKEL_BASS, T_DECKEL_BASS + T_W).fuse(box(_bx, _bx + T_W, _by, _by + BLK, T_DECKEL_BASS, T_DECKEL_BASS + H_W))
_w = _w.cut(cq.Solid.makeCylinder(D_BL / 2, 10, cq.Vector(_cx, _cy, T_DECKEL_BASS - 1)))        # Gewinde M3 fuer Fuss-/Deckelschraube
_ws = []
for _zz in (T_DECKEL_BASS + 9.0, T_DECKEL_BASS + 17.0):
    for _yy in (_cy - 3.5, _cy + 3.5) if False else (_cy,):
        _ws.append(stab((_bx + T_W + 1.0, _yy, _zz), (_bx + T_W, _yy, _zz), 5.0).fuse(stab((_bx + T_W, _yy, _zz), (_bx - 9.0, _yy, _zz), 2.5)))
        _w = _w.cut(stab((_bx + T_W + 1, _yy, _zz), (_bx - 1, _yy, _zz), 2.7))
teile = [(n_, s_) for n_, s_ in teile if n_ != "Eckblock_3"]
teile = [(n_, s_.intersect(box(0, L, 0, B, -100, T_DECKEL_BASS + T_W)) if n_ == "Fuss_Schraube_3" else s_) for n_, s_ in teile]
teile = [(n_, s_.cut(cq.Compound.makeCompound(_ws)) if n_ == "Deckelauflage_links" else s_) for n_, s_ in teile]
teile += [("Eckwinkel_3", _w)] + [(f"Eckwinkel_Schraube_{j_+1}", s_) for j_, s_ in enumerate(_ws)]

# ---------- Innenliegende Randleiste (4 Teile, 4 dick, 10 ueber den Rahmen) zur Aufnahme des Balgrahmens ----------
T_RL, U_RL = 4.0, 10.0
Z_RL0, Z_RL1 = Z_BO, H_R + U_RL
rl = [("Randleiste_hinten", box(D, L - D, D, D + T_RL, Z_RL0, Z_RL1)), ("Randleiste_vorne", box(D, L - D, FI - T_RL, FI, Z_RL0, Z_RL1)),
      ("Randleiste_links", box(D, D + T_RL, D + T_RL, FI - T_RL, Z_RL0, Z_RL1)), ("Randleiste_rechts", box(L - D - T_RL, L - D, D + T_RL, FI - T_RL, Z_RL0, Z_RL1))]
teile += rl
# Balgrahmen (wie Schnitt: auf dem Rahmen, innen an der Randleiste, 25 hoch) - in Fusion ausgeblendet
H_BALGR = 25.0
balgrahmen = box(0, L, 0, B, H_R, H_R + H_BALGR).cut(box(D, L - D, D, B - D, H_R - 1, H_R + H_BALGR + 1))
# ---------- Schieber des Begleiterstimmstocks: Hebelmechanik ----------
# Schieber liegt in der Vertiefung der Sole, Weg = eine Lochbreite (10) in Laengsrichtung. Oben ein langer Hebel neben den
# Kammern (Stift im Langloch 2 x 10 des Schiebers), Welle mit Lagerbuchse durch den Boden nach vorne versetzt (X 322,7 / Y 110),
# unten ein Hebel nach vorne zur Lufttasten-Seite; dessen Stift greift von innen in den Deckelschieber. Der Deckelschieber
# laeuft in einer geschlossenen Kammer im Deckel (Wandstaerken 0,5), aussen nur ein Schlitz mit Knopf.
SCH_HUB = 10.0                                      # eine Lochbreite
sl_ = [r_ for r_ in sch if r_[1] - r_[0] < 5][0]    # Langloch im Schieber
XS_, YS_ = (sl_[0] + sl_[1]) / 2, (sl_[2] + sl_[3]) / 2
XW_, YW_ = 322.7, 110.0                             # Welle
R_O = math.hypot(XS_ - XW_, YS_ - YW_)
_a0 = math.atan2(YS_ - YW_, XS_ - XW_)
_ALF = None
for _k in range(1, 4000):                            # Schwenkwinkel, bei dem der Stift oben um SCH_HUB in X wandert
    _t = math.radians(_k * 0.005)
    if XW_ + R_O * math.cos(_a0 + _t) - XS_ >= SCH_HUB: _ALF = _t; break
TH_ = math.degrees(_ALF)                            # gesamter Schwenk
R_U = 80.0                                          # unterer Hebel nach vorne (in der Mitte genau +Y)
_ar = math.radians(90.0) - _ALF / 2
XP_, YP_ = XW_ + R_U * math.cos(_ar), YW_ + R_U * math.sin(_ar)   # Deckelstift in Ruhe
X_HUB = 2 * R_U * math.sin(_ALF / 2)                # Weg des Deckelschiebers (X)
DY_ = R_U * (1 - math.cos(_ALF / 2))
Z_OH0, Z_OH1 = S_Z0 + 1.8, S_Z0 + 3.3               # oberer Hebel ueber dem Schieber
Z_UH0, Z_UH1 = 15.0, 16.5                           # unterer Hebel unter den Klappen (auch geoeffnet) und ueber dem Gestaenge
def hebel_flach(x0, y0, x1, y1, z0, z1, b_=6.0):
    a_ = math.atan2(y1 - y0, x1 - x0); l_ = math.hypot(x1 - x0, y1 - y0)
    h_ = box(0, l_, -b_ / 2, b_ / 2, z0, z1).fuse(cq.Solid.makeCylinder(b_ / 2, z1 - z0, cq.Vector(0, 0, z0))).fuse(cq.Solid.makeCylinder(b_ / 2, z1 - z0, cq.Vector(l_, 0, z0)))
    return h_.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), math.degrees(a_)).translate(cq.Vector(x0, y0, 0))
welle = stab((XW_, YW_, Z_UH0), (XW_, YW_, Z_OH1), 3.0)
buchse = stab((XW_, YW_, Z_BU - 0.5), (XW_, YW_, Z_BO), 5.0).cut(stab((XW_, YW_, Z_BU - 1), (XW_, YW_, Z_BO + 1), 3.1))
h_o = hebel_flach(XW_, YW_, XS_, YS_, Z_OH0, Z_OH1).cut(stab((XW_, YW_, 0), (XW_, YW_, 60), 3.0))
st_o = stab((XS_, YS_, S_Z0 + 0.2), (XS_, YS_, Z_OH0), 1.8)
h_u = hebel_flach(XW_, YW_, XP_, YP_, Z_UH0, Z_UH1).cut(stab((XW_, YW_, 0), (XW_, YW_, 60), 3.0))
# Deckel 0..4: Kammer Z 0,5..3,5 (aussen und innen 0,5 Haut), Deckelschieber 2,8 dick, Knopf durch einen Schlitz aussen
Z_K0, Z_K1 = 0.5, T_DECKEL_BASS - 0.5
YPM_ = YP_ + DY_ / 2
DS_L, DS_B = 16.0, 8.0
dschieber = box(XP_ - DS_L / 2, XP_ + DS_L / 2, YPM_ - DS_B / 2, YPM_ + DS_B / 2, Z_K0 + 0.1, Z_K1 - 0.1)
dschieber = dschieber.fuse(box(XP_ - 2.0, XP_ + 2.0, YPM_ - 1.5, YPM_ + 1.5, -2.5, Z_K0 + 0.2))           # Knopf
dschieber = dschieber.cut(cq.Workplane("XY").center(XP_, YPM_).slot2D(2.6 + DY_ + 0.4, 2.6, 90).extrude(5).val().translate(cq.Vector(0, 0, 1.5)))
st_u = stab((XP_, YP_, Z_UH0), (XP_, YP_, 1.6), 2.5)                                                  # greift von innen ein
kammer_d = box(XP_ - X_HUB - DS_L / 2 - 0.2, XP_ + DS_L / 2 + 0.2, YPM_ - DS_B / 2 - 0.2, YPM_ + DS_B / 2 + 0.2, Z_K0, Z_K1)
schlitz_a = cq.Workplane("XY").center(XP_ - X_HUB / 2, YPM_).slot2D(X_HUB + 7.4, 3.4).extrude(Z_K0 + 2).val().translate(cq.Vector(0, 0, -2))
schlitz_i = cq.Workplane("XY").center(XP_ - X_HUB / 2, YPM_).slot2D(X_HUB + 3.2, 3.2 + DY_).extrude(3).val().translate(cq.Vector(0, 0, Z_K1 - 0.1))
boden = boden.cut(stab((XW_, YW_, Z_BU - 1), (XW_, YW_, Z_BO + 1), 5.0))
aus_d = kammer_d.fuse(schlitz_a).fuse(schlitz_i)
sl_d = schlitz_i
teile = [(n_, (s_.cut(aus_d).cut(sl_d) if n_ == "Bassdeckel" else s_.cut(sl_d) if n_ == "Bassbodenleder"
               else s_.cut(box(sl_[0], sl_[1], sl_[2], sl_[3], S_Z0 - 1, S_Z0 + 3)) if n_ == "Begleiter_Schieber" else s_)) for n_, s_ in teile]
sh_teile = [("Schieberhebel_oben", h_o), ("Schieberhebel_Stift_oben", st_o), ("Schieberhebel_Welle", welle), ("Schieberhebel_Lagerbuchse", buchse),
            ("Schieberhebel_unten", h_u), ("Schieberhebel_Stift_Deckel", st_u), ("Deckelschieber", dschieber)]
teile += sh_teile
print(f"Schieberhebel: Welle X {XW_:.1f} Y {YW_:.1f}, Hebel oben {R_O:.1f} unten {R_U:.1f}, Schwenk {TH_:.2f} Grad, Schieberweg {SCH_HUB}, Deckelschieber {X_HUB:.1f}, Stift Deckel X {XP_:.1f} Y {YP_:.1f}")
# Pruefung Schieber geschlossen: Hebel, Stifte, Begleiterschieber, Deckelschieber in Endlage gegen alles andere
_rz = lambda s_: s_.rotate(cq.Vector(XW_, YW_, 0), cq.Vector(XW_, YW_, 1), TH_)
_end = {"Schieberhebel_oben": _rz(h_o), "Schieberhebel_Stift_oben": _rz(st_o), "Schieberhebel_unten": _rz(h_u), "Schieberhebel_Stift_Deckel": _rz(st_u),
        "Deckelschieber": dschieber.translate(cq.Vector(-X_HUB, 0, 0))}  # Stift unten wandert nach -X
_end["Begleiter_Schieber"] = [s_ for n_, s_ in teile if n_ == "Begleiter_Schieber"][0].translate(cq.Vector(SCH_HUB, 0, 0))
_fix = [(n_, s_) for n_, s_ in teile if n_ not in _end and n_ not in ("Schieberhebel_Welle", "Schieberhebel_Lagerbuchse")]
_ks = []
for n1, s1 in _end.items():
    b1 = s1.BoundingBox()
    for n2, s2 in _fix + [(m_, t_) for m_, t_ in _end.items() if m_ != n1]:
        b2 = s2.BoundingBox()
        if b1.xmax < b2.xmin or b2.xmax < b1.xmin or b1.ymax < b2.ymin or b2.ymax < b1.ymin or b1.zmax < b2.zmin or b2.zmax < b1.zmin: continue
        v_ = s1.intersect(s2).Volume()
        if v_ > 0.05: _ks.append((n1, n2, round(v_, 2)))
print("Schieber geschlossen, Kollisionen:", _ks)

# ---------- Knopfzylinder auf die minimal erforderliche Laenge: in Ruhe innen buendig mit dem Ende der Fuehrung (Front/Aufdoppelung) ----------
_fr = [front] + [s_ for n_, s_ in teile if n_.startswith(("Front", "Aufdoppelung", "Filzeinlage"))]
_kz = []
def _kuerzen(n_, s_):
    bb_ = s_.BoundingBox(); X_, Z_ = (bb_.xmin + bb_.xmax) / 2, (bb_.zmin + bb_.zmax) / 2
    ring_ = stab((X_, B - 40, Z_), (X_, B + 20, Z_), D_ZYL + 1.0)
    y_f = min([t_.intersect(ring_).BoundingBox().ymin for t_ in _fr if t_.intersect(ring_).Volume() > 0.01])
    _kz.append((round(bb_.ymax - y_f, 1), round(bb_.ylen, 1)))
    return stab((X_, y_f, Z_), (X_, bb_.ymax, Z_), D_ZYL)
teile = [(n_, _kuerzen(n_, s_) if n_.startswith("Knopfzylinder") else s_) for n_, s_ in teile]
print("Knopfzylinder gekuerzt (neu/alt):", sorted(set(_kz)))

# ---------- Ausgabe mit Farben wie Markus/Standard ----------
teile = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + klappen
FARBEN = [("Gaze", (0.78, 0.80, 0.82)), ("Fuss_Gummi", (0.08, 0.08, 0.08)), ("Fuss_Schraube", (0.6, 0.6, 0.62)), ("Fuss", (0.28, 0.20, 0.14)),
          ("Frontblende", (0.70, 0.52, 0.30)), ("Front_Aufdoppelung", (0.78, 0.60, 0.38)), ("Front", (0.78, 0.60, 0.38)), ("Horn", (0.05, 0.05, 0.05)), ("Knopfzylinder", (0.95, 0.95, 0.92)), ("Knopf", (0.15, 0.15, 0.15)), ("Filzeinlage", (0.55, 0.55, 0.20)), ("Rueck", (0.74, 0.56, 0.36)), ("Seiten", (0.70, 0.52, 0.33)), ("Aufdoppelung", (0.78, 0.60, 0.38)),
          ("Bassbodenleder", (0.12, 0.08, 0.05)), ("Bassboden", (0.90, 0.82, 0.62)), ("Bassdeckel_Leiste", (0.30, 0.22, 0.15)), ("Bassdeckel_Auflage", (0.60, 0.45, 0.28)),
          ("Bassdeckel", (0.30, 0.22, 0.15)), ("Randleiste", (0.62, 0.45, 0.28)), ("Schieberhebel_Lagerbuchse", (0.80, 0.65, 0.30)), ("Schieberhebel", (0.75, 0.76, 0.80)), ("Deckelschieber", (0.20, 0.15, 0.10)), ("HB_Langteil", (0.65, 0.75, 0.88)), ("HB_Vierkantleiste", (0.55, 0.40, 0.25)), ("HB_Trennfeder", (0.85, 0.30, 0.30)), ("HB_Haken", (0.40, 0.40, 0.42)), ("HB_Schraube", (0.40, 0.40, 0.42)), ("Balgrahmen", (0.55, 0.35, 0.20)), ("Deckelauflage", (0.60, 0.45, 0.28)), ("Eckblock", (0.55, 0.40, 0.25)), ("Eckwinkel_Schraube", (0.40, 0.40, 0.42)), ("Eckwinkel", (0.75, 0.76, 0.80)), ("Abdeckplatte", (0.55, 0.40, 0.25)), ("Abschlusswand", (0.55, 0.40, 0.25)), ("Innenleiste", (0.55, 0.40, 0.25)),
          ("Klappe_", None), ("HB_Folie", (0.85, 0.30, 0.30)), ("HB_Stimmplatte", (0.80, 0.82, 0.86)), ("HB_Trennwand", (0.50, 0.62, 0.80)), ("HB_", (0.65, 0.75, 0.88)),
          ("Begleiter_Stimmplatte", (0.80, 0.82, 0.86)), ("Begleiter_Sole", (0.55, 0.68, 0.45)), ("Begleiter_Schraube", (0.40, 0.40, 0.42)), ("Begleiter", (0.62, 0.72, 0.55)), ("Haken", (0.40, 0.40, 0.42)), ("Formteil_Schraube", (0.40, 0.40, 0.42)), ("Formteil", (0.60, 0.45, 0.28)), ("Dichtung", (0.12, 0.08, 0.05)), ("Lufttaste_Filz", (0.85, 0.15, 0.15)), ("Lufttaste_Gabelgelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)),
          ("Lufttaste_Schubstange", (0.10, 0.60, 0.60)), ("Lufttaste", (0.35, 0.25, 0.15)), ("Lufthebel_Feder", (0.75, 0.10, 0.75)), ("Lufthebel_Stuetze_Schraube", (0.40, 0.40, 0.42)),
          ("Lufthebel_Lagerbock_Schraube", (0.40, 0.40, 0.42)), ("Lufthebel_Lagerbock", (0.25, 0.75, 0.35)), ("Lufthebel_Achse", (0.55, 0.55, 0.60)), ("Begleiter_Schieber", (0.85, 0.30, 0.30)), ("Lufthebel_Stuetze", (0.25, 0.75, 0.35)), ("Lufthebel_Bolzen", (0.95, 0.85, 0.20)), ("Lufthebel", (0.70, 0.72, 0.78)), ("Luftklappe_Leder", (0.45, 0.30, 0.20)),
          ("Luftklappe_Filz", (0.85, 0.85, 0.80)), ("Luftklappe", (0.20, 0.20, 0.22))]
def farbe(n):
    if n.endswith("_gedrueckt"): return (0.95, 0.55, 0.55)
    if n.startswith(("Klappe_", "Begleiterklappe_")): return (0.45, 0.30, 0.20) if "Leder" in n else (0.85, 0.85, 0.80) if "Filz" in n else (0.20, 0.20, 0.22)
    for p_, c_ in FARBEN:
        if n.startswith(p_) and c_: return c_
    print("FARBE ohne Zuordnung:", n); return (0.6, 0.6, 0.6)
exec(open("mech_eduard.py").read())             # Mechanik (vor der Spiegelung)
teile.append(("Balgrahmen", balgrahmen))          # in Fusion ausblenden
NUR_LUFT = len(sys.argv) > 3 and sys.argv[3] == "luft"
if NUR_LUFT:                                   # nur Hebeleinbau zeigen: Deckel, Boden, Eckbloecke, Luftgruppe in Ruhe und gedrueckt
    teile = [(n, s_) for n, s_ in teile if n.startswith(("Bassboden", "Bassdeckel", "Eckblock", "Eckwinkel", "Fuss_Schraube_3", "Deckelauflage", "Luft", "Front", "Knopf")) and not n.startswith("Bassbodenleder")] + lg_g
a = cq.Assembly(name="Bassteil_Eduard_Winkelbass")
SPIEGELN = True                                # wie Markus: links/rechts getauscht (X = L/2)
for n, s_ in teile:
    if SPIEGELN: s_ = s_.mirror("YZ", basePointVector=(L / 2, 0, 0))
    a.add(s_, name=n, color=cq.Color(*farbe(n)))
a.save(OUT)
print("gespeichert", OUT, len(teile), "Teile")
