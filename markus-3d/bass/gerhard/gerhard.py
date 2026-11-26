# Bassteil Gerhard (Zeichnung BassteilV3-in_verwendung.dxf: Draufsicht + Schnitt), Einbau wie Eduard, Gehaeuse wie Markus:
# Rahmen mit Zinken, Deckel mit Rosetten, hinterer Leiste, Leder, Eckbloecken und Fuessen, Deckelauflage, Gaze von Markus.
# Zeichnung 380 x 205 -> Rohmass 382 x 207 (Markus): Zeichnungskoordinaten + 1 mm (OX, OY).
# Achsen: X = Laenge, Y = Tiefe (0 = Rueckteil aussen, B = Front aussen), Z ab Deckel aussen Richtung Balg.
# Aufruf: python3 gerhard.py ausgabe.step [Knopfteilung]
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
XP0, YP0 = 3403.8, 1536.6                     # Draufsicht (BassteilV3) -> Modell: X = x - XP0, Y = y - YP0 (Rohmass +1)
XS0, YS0 = 3861.0, 1536.6                      # Schnitt -> Modell: Z = x - XS0, Y = y - YS0
# ---------- Front wie Markus: Knopfreihen-Abstufung (Reihe 1 vertieft, Aufdoppelung 2 + 5 mm), Hoerner, Blende, Knoepfe ----------
RD = ezdxf.readfile("rahmen.dxf").modelspace()
# Knopfreihen wie Standardmodell: Reihe 1 (unten, 9 Knoepfe, Bohrung Ø11 mit Filz) buendig an der Front, Reihe 2 (10 Knoepfe,
# Ø10) ueber eine 7-mm-Blende erhoeht; Hoehen und Teilung nach der Marie-Zeichnung (Teilung variabel), Mitte X 198.
X_KM = L / 2                                   # Knopfreihen mittig am Korpus
REIHEN_N = {1: 8, 2: 11, 3: 10}                  # Gerhard: 29 Knoepfe in drei Reihen (Zeichnung)
Z_RH = {1: 14.0, 2: 28.0, 3: 42.0}               # Schnitt: Bassdruecker 3875 / 3889 / 3903 - 3861
R_RH = {1: 4.25, 2: 4.0, 3: 4.0}                  # Stift 6,5 (Zeichnung), Filzeinlage
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
_hp_alle = list(horn_pos)
horn_pos = [(X_, Z_) for X_, Z_ in horn_pos if all(math.hypot(X_ - k[0], Z_ - k[1]) > 9.0 + 2.5 + 15.5 / 2 + 1.0 for k in kn)]   # Gerhard: Hoerner an Reihe 3 entfallen
HORN_WEG = [p_ for p_ in _hp_alle if p_ not in horn_pos]
for X_, Z_ in HORN_WEG:                                 # Hornbohrung in der Front zu, Knopfbohrungen wieder auf
    front = front.fuse(stab((X_, B - D - 0.5, Z_), (X_, B + 0.5, Z_), 2 * 9.0 + 0.4)).intersect(box(0, L, B - D, B, 0, H_R))
for cyl in knopf_bohr: front = front.cut(cyl)
for X_, Z_ in horn_pos:
        r = 9.0
        ring = cq.Workplane("XZ").center(X_, Z_).circle(r + 2.5).circle(r).extrude(-40).translate((0, B - 30, 0)).val()
        y_out = front.intersect(ring).BoundingBox().ymax
        w = cq.Workplane("XY").polyline([(r_, a_) for a_, r_ in horn_prof]).close().revolve(360, (0, 0, 0), (0, 1, 0))
        teile.append((f"Horn_{len(aus_front)+1}", w.val().translate(cq.Vector(X_, y_out - 4.01, Z_))))
        aus_front.append(cq.Workplane("XZ").center(X_, Z_).circle(r + 5.0).extrude(-(6.0 + 2)).translate((0, FI - 6.0 - 1, 0)).val())
# Blende ueber Reihe 2 (wie Standardmodell: Bogen um jede Bohrung, Einzug dazwischen, grosse Endboegen), 7 mm
D_ZYL, D_KOPF, H_KOPF, L_ZYL, SPALT, T_BL = 6.5, 15.5, 4.5, 20.0, 6.5, 7.0
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
def T_BLR(r_): return T_BL * (r_ - 1)            # Reihe 2 auf 7-mm-Blende, Reihe 3 auf 14-mm-Blende
# Frontrahmen im Bereich der Reihen 2 und 3 ausgeschnitten: 2 mm Ueberlappung mit Stufenprofil (Falz halbe Wandstaerke im
# Rahmen, Lippe an der Blende) fuer die Verklebung
from shapely.geometry import Polygon as _SPoly
from shapely.ops import unary_union as _sunion
UEB, T_FALZ = 2.0, D / 2
_bk = {}
for rr_ in (2, 3):
    xr_ = [k[0] for k in kn if k[3] == rr_]
    _bk[rr_] = _SPoly([(-x_, Z_RH[rr_] + y_) for x_, y_ in blende_kontur([-x_ for x_ in xr_])]).buffer(0)
O_BL = _sunion(list(_bk.values())).buffer(2.0, join_style=2).buffer(-2.0, join_style=2)   # beide Blenden zu einem Ausschnitt
_bk[2] = O_BL.difference(_bk[3]).buffer(-0.001).buffer(0.001)          # Blende R2 deckt auch den Spalt bis zur Blende R3
if _bk[2].geom_type != "Polygon": _bk[2] = max(_bk[2].geoms, key=lambda g_: g_.area)
O_IN = O_BL.buffer(-UEB, join_style=2)
def xz_prisma(pg_, y0_, y1_):
    pg_ = pg_.simplify(0.02)
    return cq.Workplane("XZ").polyline(list(pg_.exterior.coords)[:-1]).close().extrude(-(y1_ - y0_)).val().translate(cq.Vector(0, y0_, 0))
front = front.cut(xz_prisma(O_IN, B - D - 1, B + 1)).cut(xz_prisma(O_BL, B - T_FALZ, B + 1))
for rr_ in (2, 3):
    bl = xz_prisma(_bk[rr_], B, B + T_BLR(rr_))
    lip_ = _bk[rr_].difference(O_IN)
    for g_ in (getattr(lip_, "geoms", None) or [lip_]):
        if g_.area > 0.5: bl = bl.fuse(xz_prisma(g_, B - T_FALZ, B + 0.01))
    for X_, Z_, r_, rr in kn:
        if rr == rr_: bl = bl.cut(stab((X_, B - D - 1, Z_), (X_, B + T_BLR(rr_) + 1, Z_), 2 * r_))
    teile.append((f"Frontblende_R{rr_}", bl.clean()))
FX0, FX1 = O_BL.bounds[0], O_BL.bounds[2]           # X-Bereich des Frontausschnitts
# Knoepfe in Ruhelage (wie Standardmodell), Filzeinlage in der aeussersten Lage (Front bzw. Blende)
for X_, Z_, r_, rr in kn:
    ya = B + T_BLR(rr); dl_ = T_BLR(rr) if rr > 1 else D; yk = ya + SPALT; lz = yk - (B + SPALT - L_ZYL)
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
B_AUF_H = 2.0                                   # Gerhard: Deckelauflage hinten 2 mm (Platz fuer die hinteren Achsreihen)
dh = box(D, L - D, D, D + B_AUF_H, za0, Z_BU)
for a_ in aus_rueck: dh = dh.cut(a_)
za1 = za0 + 6.0
def mit_aus(s_, l_):
    for a_ in l_: s_ = s_.cut(a_)
    return s_
teile += [("Deckelauflage_hinten", dh), ("Deckelauflage_links", box(D, D + B_AUF, D + B_AUF_H, FI, za0, Z_BU)),
          ("Deckelauflage_rechts", box(L - D - B_AUF, L - D, D + B_AUF_H, FI, za0, Z_BU)),
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
yi0, yi1 = D + B_AUF_H, FI - B_AUF
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
# ======================= Gerhard: Einbauten =======================
P = lambda x_, y_: (x_ - XP0, y_ - YP0)                 # Draufsicht -> Modell (vor der Spiegelung)
DOC = ezdxf.readfile("gerhard.dxf"); MD = DOC.modelspace()
def rect(e):
    p_ = [P(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.2)]; xs_ = [q[0] for q in p_]; ys_ = [q[1] for q in p_]
    return (min(xs_), max(xs_), min(ys_), max(ys_))
def im_plan(r_): return -5 < r_[0] and r_[1] < L + 5 and -5 < r_[2] and r_[3] < B + 5
VE = []                                                   # Bloecke (Begleiter) aufgeloest
for ins in MD.query("INSERT"):
    for v_ in ins.virtual_entities():
        lay_ = v_.dxf.layer if v_.dxf.layer != "0" else ins.dxf.layer
        VE.append((lay_, v_))
# ---------- Boden 3,5 (Z 35..38,5) ueber die ganze Innenflaeche; Oeffnungen und Klappen aus der Draufsicht ----------
boden = box(D, L - D, D, FI, Z_BU, Z_BO)
oeff = sorted([r_ for r_ in (rect(e) for e in MD.query('LWPOLYLINE[layer=="Füllung-Rund-groß-Begleiter"]')) if im_plan(r_)], key=lambda r_: (r_[2], r_[0]))
LUFT_OEFF = [r_ for r_ in oeff if r_[1] - r_[0] > 60][0]
hb_oeff = [r_ for r_ in oeff if r_ is not LUFT_OEFF]
RAND_K, SPIEL, TIEF_V = 3.0, 0.5, 2.0
klappen = []
def klappe_rect(nm, x0, x1, y0, y1, rand=RAND_K, durch=True):
    global boden                                   # durch=False: nur Vertiefung (Oeffnung separat)
    if durch: boden = boden.cut(box(x0, x1, y0, y1, Z_BU - 1, Z_BO + 1))
    boden = boden.cut(box(x0 - rand - SPIEL, x1 + rand + SPIEL, y0 - rand - SPIEL, y1 + rand + SPIEL, Z_BU - 1, Z_BU + TIEF_V))
    z_ = Z_BU + TIEF_V
    for t_, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
        klappen.append((f"{nm}_{t_}", box(x0 - rand, x1 + rand, y0 - rand, y1 + rand, z_ - dk, z_))); z_ -= dk
for i, (x0, x1, y0, y1) in enumerate(hb_oeff): klappe_rect(f"Klappe_{i+1}", x0, x1, y0, y1)
# ---------- Grundbass: 11 Winkelbass-Kammern versetzt (vorne 5 stehende, hinten 6 liegende Platten) + 1 liegende Kammer ----------
# Schnitt (Z = x - 3861, Y = y - YP0): vordere Kammer mit stehender Platte zwischen den Vierkantleisten, hintere Kammer mit
# liegender Platte oben; hinteres Langteil mit zwei Knicken (stumpf gestossen), gebogenes Mittelteil (gerades Stueck +
# Bogen), Trennfeder durchgehend. Vordere und hintere Trennwaende um eine halbe Teilung versetzt (Teilung 47,1).
# Platten wie die schmalen von Markus: 46 x 94, Keil 5/16, Rand 2, 4 Zungenschlitze.
ZB = Z_BO
def pzy(x0, x1, zy): return py(x0, x1, [(y_, z_) for z_, y_ in zy])
hinten = sorted([r_ for r_ in hb_oeff if r_[3] < 110], key=lambda r_: r_[0])
vorne_all = sorted([r_ for r_ in hb_oeff if r_[2] > 110], key=lambda r_: r_[0])
r_lieg = vorne_all[-1]; vorne = vorne_all[:-1]
cx_h = [(r_[0] + r_[1]) / 2 for r_ in hinten]; cx_v = [(r_[0] + r_[1]) / 2 for r_ in vorne]
TEIL = (cx_h[-1] - cx_h[0]) / (len(cx_h) - 1); T_TW = 3.5
XD_H = [cx_h[0] - TEIL / 2 + TEIL * k for k in range(len(cx_h) + 1)]
XD_V = [cx_v[0] - TEIL / 2 + TEIL * k for k in range(len(cx_v) + 1)]
HX0, HX1 = XD_H[0] - T_TW / 2, XD_H[-1] + T_TW / 2
VX0, VX1 = XD_V[0] - T_TW / 2, XD_V[-1] + T_TW / 2
Y_PV, Z_PO = 157.5, 153.5                                  # Auflage stehende Platte (Y) / liegende Platte (Z)
Z_VLU, Z_VLO = 53.5, 143.7
PL_L, PL_B = 94.0, 46.0
lang = {
    "HB_Langteil_hinten_1": pzy(HX0, HX1, [(ZB, 65.6), (94.6, 65.6), (96.1, 70.6), (ZB, 70.6)]),
    "HB_Langteil_hinten_2": pzy(HX0, HX1, [(94.6, 65.6), (121.5, 47.5), (123.0, 52.5), (96.1, 70.6)]),
    "HB_Langteil_hinten_3": pzy(HX0, HX1, [(121.5, 47.5), (143.5, 47.5), (143.5, 52.5), (123.0, 52.5)]),
    "HB_Langteil_mitte": pzy(HX0, HX1, [(ZB, 99.3), (109.6, 99.3), (109.6, 101.8), (ZB, 101.8)]),
}
BOG_A = [(109.6, 99.3), (120.1, 101.0), (129.3, 105.8), (136.5, 113.3), (141.1, 122.9), (148.5, 148.6)]
BOG_I = [(109.6, 101.8), (119.3, 103.4), (127.8, 107.8), (134.4, 114.8), (138.7, 123.6), (145.6, 147.5)]
vl = {
    "HB_Vierkantleiste_unten": pzy(VX0, VX1, [(ZB, 147.5), (Z_VLU, 147.5), (Z_VLU, Y_PV), (ZB, Y_PV)]),
    "HB_Vierkantleiste_oben": pzy(HX0, HX1, [(Z_VLO, 147.5), (Z_PO, 147.5), (Z_PO, Y_PV), (Z_VLO, Y_PV)]),
    "HB_Vierkantleiste_hinten": pzy(HX0, HX1, [(143.5, 46.5), (Z_PO, 46.5), (Z_PO, 56.5), (143.5, 56.5)]),
}
lang["HB_Langteil_gebogen"] = pzy(HX0, HX1, BOG_A + BOG_I[::-1]).cut(vl["HB_Vierkantleiste_oben"])
TW_H = [(ZB, 99.3), (ZB, 68.1), (95.3, 68.1), (122.3, 50.0), (Z_PO, 50.0), (Z_PO, 152.5), (148.5, 152.5), (148.5, 148.6)] + BOG_A[::-1][1:]
TW_V = [(ZB, 101.8)] + BOG_I + [(148.5, 147.5), (148.5, Y_PV), (48.5, Y_PV), (48.5, 147.5), (ZB, 147.5)]
feder = pzy(VX0 + T_TW / 2, VX1 - T_TW / 2, [(100.8, 126.8), (59.3, 121.7), (48.1, 148.4), (46.2, 147.5), (58.0, 119.5), (101.0, 124.8)]).cut(vl["HB_Vierkantleiste_unten"])
tw = [(f"HB_Trennwand_hinten_{k+1}", pzy(xd - T_TW / 2, xd + T_TW / 2, TW_H)) for k, xd in enumerate(XD_H)]
tw += [(f"HB_Trennwand_vorne_{k+1}", pzy(xd - T_TW / 2, xd + T_TW / 2, TW_V).cut(feder)) for k, xd in enumerate(XD_V)]
_tws = cq.Compound.makeCompound([s_ for _, s_ in tw])
for n_, s_ in list(lang.items()) + list(vl.items()): teile.append((n_, s_.cut(_tws)))        # Nuten fuer die Trennwaende
teile += tw + [("HB_Trennfeder", feder)]
# liegende Kammer (Querschnitt wie Standard): Kasten auf dem Boden, Platte oben, Hoehe 54 wie Markus
LX0, LX1 = (r_lieg[0] + r_lieg[1]) / 2 - 24.5, (r_lieg[0] + r_lieg[1]) / 2 + 24.5
LY0, LY1 = 102.0, FI - 4.0 - 0.1                     # hinter dem gebogenen Mittelteil .. vor der Randleiste
Z_LK = ZB + 54.0
lk = box(LX0, LX1, LY0, LY1, ZB, Z_LK).cut(box(LX0 + 3.0, LX1 - 3.0, LY0 + 2.0, LY1 - 2.0, ZB - 1, Z_LK + 1))
teile.append(("HB_Kammer_liegend", lk))
T_DUENN, T_DICK, T_RAND, RUECK = 5.0, 16.0, 2.0, 2.0
SCHLITZE = [(3.8, 11.3, 37.5, 88.5), (12.8, 22.0, 17.0, 88.5), (24.0, 33.2, 17.0, 88.5), (34.7, 42.2, 37.5, 88.5)]   # Markus 46er Platte
def platte_lokal():
    pl = box(0, PL_B, 0, PL_L, 0, T_RAND)
    prof = [(RUECK, 0.0), (PL_L - RUECK, 0.0), (PL_L - RUECK, T_RAND - T_DICK), (RUECK, T_RAND - T_DUENN)]
    keil = cq.Workplane("YZ").polyline(prof).close().extrude(PL_B - 2 * RUECK).val().translate(cq.Vector(RUECK, 0, 0))
    pl = pl.fuse(keil)
    for (u0, u1, w0, w1) in SCHLITZE: pl = pl.cut(box(u0, u1, w0, w1, -T_DICK - 1, T_RAND + 1))
    return pl
PL0 = platte_lokal()
Z_STP = (Z_VLU + Z_VLO) / 2 - PL_L / 2                # stehende Platte unten
Y_LP = (56.5 + 147.5) / 2 - PL_L / 2                  # liegende Platte hinten
def stehend(xc): return PL0.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -90).translate(cq.Vector(xc - PL_B / 2, Y_PV, Z_STP + PL_L))
def liegend(xc, y0=Y_LP, z0=Z_PO): return PL0.mirror("XZ").translate(cq.Vector(xc - PL_B / 2, y0 + PL_L, z0))
D_HAKEN, D_KOPF, H_KOPF, D_SCHAFT = 1.8, 5.0, 1.5, 2.5
def haken(p_, r_lang, r_kurz):
    p_ = cq.Vector(*p_)
    return stab(p_.toTuple(), (p_ + r_lang).toTuple(), D_HAKEN).fuse(stab(p_.toTuple(), (p_ + r_kurz).toTuple(), D_HAKEN)).fuse(cq.Solid.makeSphere(D_HAKEN / 2 + 0.01, p_))
def schraube(p_, r_):                                 # p_ = Kopfunterseite, r_ = Richtung in das Holz (Laenge = Schaft)
    p_ = cq.Vector(*p_); n_ = r_.normalized()
    return stab(p_.toTuple(), (p_ - n_ * H_KOPF).toTuple(), D_KOPF).fuse(stab(p_.toTuple(), (p_ + r_).toTuple(), D_SCHAFT))
befest = []
for k, xc in enumerate(cx_v):
    teile.append((f"HB_Stimmplatte_stehend_{k+1}", stehend(xc)))
    yh, zh = Y_PV + T_RAND + D_HAKEN / 2, Z_STP + PL_L + 0.1 + D_HAKEN / 2
    befest.append((f"HB_Haken_stehend_{k+1}", haken((xc, yh, zh), cq.Vector(0, -9.0 - T_RAND - D_HAKEN / 2, 0), cq.Vector(0, 0, -5.0))))
    befest.append((f"HB_Schraube_stehend_{k+1}", schraube((xc, Y_PV + T_RAND, Z_STP - 1.5), cq.Vector(0, -9.0 - T_RAND, 0))))
for k, xc in enumerate(cx_h):
    teile.append((f"HB_Stimmplatte_liegend_{k+1}", liegend(xc)))
    yh, zh = Y_LP + PL_L + 0.1 + D_HAKEN / 2, Z_PO + T_RAND + D_HAKEN / 2
    befest.append((f"HB_Haken_liegend_{k+1}", haken((xc, yh, zh), cq.Vector(0, 0, -9.0 - T_RAND - D_HAKEN / 2), cq.Vector(0, -5.0, 0))))
    befest.append((f"HB_Schraube_liegend_{k+1}", schraube((xc, Y_LP - 1.5, Z_PO + T_RAND), cq.Vector(0, 0, -9.0 - T_RAND))))
xl_ = (LX0 + LX1) / 2; yl0_ = (LY0 + LY1) / 2 - PL_L / 2
teile.append(("HB_Stimmplatte_liegend_Kammer", liegend(xl_, yl0_, Z_LK)))
befest.append(("HB_Haken_liegend_Kammer", haken((xl_, yl0_ + PL_L + 0.1 + D_HAKEN / 2, Z_LK + T_RAND + D_HAKEN / 2), cq.Vector(0, 0, -6.0 - T_RAND), cq.Vector(0, -5.0, 0))))
teile.append(("HB_Kammer_liegend_Aufdopplung", box(HX1 + 0.3, LX1, LY0 - 4.0, LY0 + 0.01, Z_LK - 12.0, Z_LK)))   # aussen an der Stirnwand, damit der Haken Holz hat
befest.append(("HB_Haken_liegend_Kammer_2", haken((xl_, yl0_ - 0.1 - D_HAKEN / 2, Z_LK + T_RAND + D_HAKEN / 2), cq.Vector(0, 0, -6.0 - T_RAND), cq.Vector(0, 5.0, 0))))
_bf = cq.Compound.makeCompound([s_ for _, s_ in befest])
teile = [(n_, s_.cut(_bf) if n_.startswith(("HB_Vierkantleiste", "HB_Kammer_liegend")) else s_) for n_, s_ in teile]
teile += befest
print("Helikon: Teilung", round(TEIL, 2), "hinten", [round(c, 1) for c in cx_h], "vorne", [round(c, 1) for c in cx_v], "liegende Kammer X", round(LX0, 1), "-", round(LX1, 1))
# ---------- Begleiterstimmstock (Block begleiter-165-4fach): 12 Begleiter, je 2 x 2 Kammern Teilung 16,5, Platten 16 x 42,5 x 3 ----------
bw_ = [rect(v_) for l_, v_ in VE if l_ == "Begleitstimmstock-Trennwände1" and v_.dxftype() == "LWPOLYLINE"]
bp_ = [rect(v_) for l_, v_ in VE if l_ == "Stimmplatten" and v_.dxftype() == "LWPOLYLINE"]
bl_ = [P(v_.dxf.center.x, v_.dxf.center.y) + (v_.dxf.radius,) for l_, v_ in VE if l_ == "1Füllung-öffnugen-rund-klein-begleiter" and v_.dxftype() == "CIRCLE"]
bo_ = [P(v_.dxf.center.x, v_.dxf.center.y) + (v_.dxf.radius,) for l_, v_ in VE if l_ == "Füllung-Rund-groß-Begleiter" and v_.dxftype() == "CIRCLE"]
bk_ = [rect(v_) for l_, v_ in VE if l_ == "Klappen-rand-aussen" and v_.dxftype() == "LWPOLYLINE"]
bw_ = [r_ for r_ in bw_ if im_plan(r_)]; bp_ = [r_ for r_ in bp_ if im_plan(r_)]
bl_ = [c_ for c_ in bl_ if 0 < c_[0] < L and 0 < c_[1] < B]; bo_ = [c_ for c_ in bo_ if 0 < c_[0] < L and 0 < c_[1] < B]
bk_ = [r_ for r_ in bk_ if im_plan(r_)]
S_Z0, S_Z1, P_Z1, W_Z1, K_Z1 = ZB, ZB + 3.5, ZB + 46.0, ZB + 44.0, ZB + 47.5
sole_all, bg_all, nb_ = [], [], 0
BG_SCHLITZE, BG_SZ0, BG_SZ1 = [(3.5, 6.5), (9.5, 12.5)], S_Z1 + 5.0, S_Z1 + 35.0   # Platte 16 breit: 2 Schlitze 3 x 30
_bgx = []
gruppen = {}
_yc = sorted({round((r_[2] + r_[3]) / 2, 1) for r_ in bp_})
_gk = {}
for y_ in _yc:                                        # je Stimmstock zwei Plattenlagen (hinten/vorne) < 40 auseinander
    k_ = next((k for k in _gk.values() if abs(y_ - k) < 40), y_); _gk[y_] = k_
for r_ in bp_:
    gruppen.setdefault(len({v for v in _gk.values() if v <= _gk[round((r_[2] + r_[3]) / 2, 1)]}), []).append(r_)
_bgs = []
for g_, pl_ in gruppen.items():
    y_h = min(r_[2] for r_ in pl_); y_v = max(r_[3] for r_ in pl_)
    x0_ = min(r_[0] for r_ in bw_ if y_h < r_[2] < y_v); x1_ = max(r_[1] for r_ in bw_ if y_h < r_[2] < y_v)
    w_ = [r_ for r_ in bw_ if y_h < r_[2] < y_v]
    y_mw0 = max(r_[3] for r_ in w_ if r_[3] < (y_h + y_v) / 2); y_mw1 = min(r_[2] for r_ in w_ if r_[2] > (y_h + y_v) / 2)
    y_pi0 = min(r_[3] for r_ in pl_ if r_[2] < y_mw0); y_pi1 = max(r_[2] for r_ in pl_ if r_[2] > y_mw1)
    sole = box(x0_, x1_, y_h, y_v, S_Z0, S_Z1)
    bg = box(x0_, x1_, y_mw0, y_mw1, S_Z1, W_Z1).fuse(box(x0_, x1_, y_pi0, y_pi1, W_Z1, K_Z1))
    for r_ in w_: bg = bg.fuse(box(r_[0], r_[1], r_[2], r_[3], S_Z1, W_Z1))
    for (x_, y_, rr_) in bo_:                          # Sole: ein Loch je Begleiter wie die Bodenoeffnung (nicht vier)
        if x0_ < x_ < x1_ and y_h < y_ < y_v: sole = sole.cut(cq.Solid.makeCylinder(rr_, 10, cq.Vector(x_, y_, S_Z0 - 1)))
    for r_ in pl_:
        pl_s = box(r_[0], r_[1], r_[2], r_[3], S_Z1, P_Z1 - 1.0)
        for (u0_, u1_) in BG_SCHLITZE: pl_s = pl_s.cut(box(r_[0] + u0_, r_[0] + u1_, r_[2] - 1, r_[3] + 1, BG_SZ0, BG_SZ1))   # 2 Zungenschlitze
        nb_ += 1; teile.append((f"Begleiter_Stimmplatte_{nb_}", pl_s))
    teile += [(f"Begleiter_Sole_{g_}", sole), (f"Begleiter_Stimmstock_{g_}", bg)]
    _bgs.append(box(x0_, x1_, y_h, y_v, S_Z0, K_Z1)); _bgx.append((g_, x0_, x1_, y_h, y_v, (y_mw0 + y_mw1) / 2, len(pl_)))
# --- Befestigung: langer Stimmstock je Ende ein Schraubhaken in der Randleiste (greift oben auf den Stimmstock);
#     kurzer Stimmstock: L-Plattfeder ueber den Stimmstock, an der Helikon-Kanzelwand verschraubt, U-Begrenzung am Boden ---
def schr_k(p_, r_, dk=4.0, hk=1.2, ds=2.0):            # kleine Schraube: p_ = Kopfunterseite, r_ = Schaft
    p_ = cq.Vector(*p_); n_ = r_.normalized()
    return stab(p_.toTuple(), (p_ - n_ * hk).toTuple(), dk).fuse(stab(p_.toTuple(), (p_ + r_).toTuple(), ds))
T_RIEG, B_RIEG, T_FED, B_FED, B_U, SP_U = 1.5, 6.0, 0.8, 8.0, 3.0, 0.2
_rl_schlitz, _bg_bohr, _tw_bohr = [], {}, []
for (g_, x0_, x1_, y_h, y_v, ymw_, n_) in _bgx:
    if n_ > 4:                                         # langer Stimmstock: je Ende ein Schraubhaken (wie Bassplatten), aussen in der Randleiste
        for sd_, xe_, sg_ in (("links", D + 4.0, 1), ("rechts", L - D - 4.0, -1)):
            hk_ = haken((xe_ + sg_ * (D_HAKEN / 2 + 0.1), ymw_, K_Z1 + D_HAKEN / 2 + 0.1), cq.Vector(-sg_ * 4.5, 0, 0), cq.Vector(sg_ * 5.0, 0, 0))
            teile.append((f"Begleiter_Haken_{sd_}", hk_)); _rl_schlitz.append(hk_)
    else:                                              # kurzer Stimmstock (Begleiter 12)
        yc_ = (y_h + y_v) / 2; xw_ = VX0                 # Aussenseite Trennwand vorne 1 (Helikon-Kanzel)
        fd = box(x0_ + 4.0, xw_ - T_FED, yc_ - B_FED / 2, yc_ + B_FED / 2, K_Z1, K_Z1 + T_FED)
        fd = fd.fuse(box(xw_ - T_FED, xw_, yc_ - B_FED / 2, yc_ + B_FED / 2, K_Z1 - 15.0, K_Z1 + T_FED))
        scs_ = [schr_k((xw_ - T_FED, yc_, z_), cq.Vector(3.0, 0, 0), 4.0, 1.0, 1.8) for z_ in (K_Z1 - 11.5, K_Z1 - 5.0)]
        teile.append(("Begleiter_Plattfeder", fd.cut(cq.Compound.makeCompound(scs_))))
        teile += [(f"Begleiter_Plattfeder_Schraube_{k+1}", s_) for k, s_ in enumerate(scs_)]
        _tw_bohr += scs_
        ua_, ub_, uc_, ud_ = x0_ - SP_U, x1_ + SP_U, y_h - SP_U, y_v + SP_U
        teile.append(("Begleiter_Begrenzung_Rahmen", box(ua_ - B_U, ub_ + B_U, uc_ - B_U, ud_ + B_U, S_Z0, S_Z1).cut(box(ua_, ub_, uc_, ud_, S_Z0 - 1, S_Z1 + 1))))   # rundum geschlossen
teile = [(n_, s_.cut(cq.Compound.makeCompound(_bg_bohr[int(n_.rsplit("_", 1)[1])]))) if n_.startswith("Begleiter_Stimmstock_") and int(n_.rsplit("_", 1)[1]) in _bg_bohr
         else (n_, s_.cut(cq.Compound.makeCompound(_tw_bohr))) if n_ == "HB_Trennwand_vorne_1" and _tw_bohr else (n_, s_) for n_, s_ in teile]
# Begleiterklappen rund Ø24 (Klappenrand rund), Vertiefung Ø25 x 2, Bodenoeffnung rund Ø18 wie gezeichnet
bkr_ = [P(v_.dxf.center.x, v_.dxf.center.y) + (v_.dxf.radius,) for l_, v_ in VE if l_ == "-Klappenrand-rund -begleiter" and v_.dxftype() == "CIRCLE"]
bkr_ = sorted([c_ for c_ in bkr_ if 0 < c_[0] < L and 0 < c_[1] < B], key=lambda c_: (c_[1], c_[0]))
def klappe_rund(nm, x_, y_, r_):
    global boden
    boden = boden.cut(cq.Solid.makeCylinder(r_ + SPIEL, TIEF_V + 1, cq.Vector(x_, y_, Z_BU - 1)))
    z_ = Z_BU + TIEF_V
    for t_, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
        klappen.append((f"{nm}_{t_}", cq.Solid.makeCylinder(r_, dk, cq.Vector(x_, y_, z_ - dk)))); z_ -= dk
for j, (x_, y_, r_) in enumerate(bkr_): klappe_rund(f"Begleiterklappe_{j+1}", x_, y_, r_)
for (x_, y_, rr_) in bo_: boden = boden.cut(cq.Solid.makeCylinder(rr_, 10, cq.Vector(x_, y_, Z_BU - 1)))
bk_ = bkr_
print("Begleiter: Stimmstoecke", len(gruppen), "Platten", nb_, "Loecher", len(bl_), "Klappen", len(bk_), "Bodenoeffnungen", len(bo_))
# ---------- Lufttaste + Luftklappe + Lufthebel (zuerst, die Bassgestaenge werden spaeter daran vorbeigefuehrt) ----------
x0, x1, y0, y1 = LUFT_OEFF
klappe_rect("Luftklappe", x0, x1, y0, y1)
luft_kl = (x0 - RAND_K, x1 + RAND_K, y0 - RAND_K, y1 + RAND_K)
lh_r = [rect(e) for e in MD.query('LWPOLYLINE[layer=="Ludttaste"]')][0]
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
    g_ = g_.cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6)).cut(stab((XA, LY, ZA_ - 3.5), (XA, LY, ZA_ - G_L - 1), 3.1))
    A0_ = cq.Vector(XA, LY, ZA_); d_ = (P_ - A_).normalized(); ziel = A_ + d_ * 10.0
    wx = _m.degrees(_m.atan2(abs(d_.y), abs(d_.z))); wy = _m.degrees(_m.atan2(abs(d_.x), _m.hypot(d_.y, d_.z)))
    best_ = None
    for sx in (-1, 1):
        for sy in (-1, 1):
            f_ = lambda s_, sx=sx, sy=sy: s_.translate(A_ - A0_).rotate(A_, A_ + cq.Vector(1, 0, 0), sx * wx).rotate(A_, A_ + cq.Vector(0, 1, 0), sy * wy)
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
teile.append(("Bassdeckel_Auflage_hinten", box(xi0 + BLK + 0.5, xi1 - BLK - 0.5, _yfuge - 3.0, _yfuge + 3.0, T_DECKEL_BASS, T_DECKEL_BASS + 1.2)))   # Gerhard: 1,2 hoch (Gelenke der hinteren Reihen)
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
    if -382 < x < 0 and -4.24 < y < 71.75 and (round(-x, 2), round(y + 4.24, 2)) not in HORN_WEG:
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
_rl_aus = [b_.translate(cq.Vector(0, 0, 0)) for b_ in _bgs] + _rl_schlitz + [stab((X_, FI - 10, Z_), (X_, B + 1, Z_), D_ZYL + 1.0) for X_, Z_, r_, rr in kn if Z_ > Z_RL0]
rl = [(n_, s_.cut(cq.Compound.makeCompound(_rl_aus))) for n_, s_ in rl]   # Begleiterstimmstock und Reihe-3-Knoepfe ausgespart
teile += rl
# Bereich hinter den Reihen 2/3: Boden ausgeschnitten (Platz fuer das Gestaenge), zum Balg hin mit einer Schraege und zwei
# Seitenwaenden abgedichtet (ersetzt den gezeichneten Winkel); um die Luftklappe bleibt eine Auflage im Boden.
T_SCH = 3.0
Y_S0, Y_S1, Z_S1 = 180.0, FI - T_RL, Z_RH[3] + 10.5       # Schraege: hinten auf dem Boden, vorne an der Randleiste
_sx0, _sx1 = FX0, FX1
_aus = box(_sx0, _sx1, Y_S0 + T_SCH + 0.5, FI, Z_BU - 1, Z_BO + 1)
_rest = box(luft_kl[0] - SPIEL - 3.0, luft_kl[1] + SPIEL + 3.0, luft_kl[2] - SPIEL - 3.0, luft_kl[3] + SPIEL + 3.0, Z_BU - 2, Z_BO + 2)
boden = boden.cut(_aus.cut(_rest))
_n = (Z_S1 - Z_BO, -(Y_S1 - Y_S0)); _nl = math.hypot(*_n); _n = (_n[0] / _nl * T_SCH, _n[1] / _nl * T_SCH)
schraege = py(_sx0, _sx1, [(Y_S0, Z_BO), (Y_S1, Z_S1), (Y_S1, Z_S1 + T_SCH * 1.3), (Y_S0 - _n[0], Z_BO - _n[1])]).cut(box(0, L, 0, B, 0, Z_BO))
teile.append(("Ausschnitt_Schraege", schraege))
for k_, xw_ in enumerate((_sx0 - T_SCH, _sx1)):
    teile.append((f"Ausschnitt_Seitenwand_{k_+1}", py(xw_, xw_ + T_SCH, [(Y_S0 - _n[0], Z_BO), (Y_S1, Z_BO), (Y_S1, Z_S1 + T_SCH * 1.3), (Y_S0 - _n[0], Z_BO - _n[1])])))
teile = [(n_, s_.cut(box(_sx0 - T_SCH, _sx1 + T_SCH, FI - T_RL - 0.1, FI + 0.1, Z_RL0 - 1, Z_S1)) if n_ == "Randleiste_vorne" else s_) for n_, s_ in teile]
# Balgrahmen (wie Schnitt: auf dem Rahmen, innen an der Randleiste, 25 hoch) - in Fusion ausgeblendet
H_BALGR = 25.0
balgrahmen = box(0, L, 0, B, H_R, H_R + H_BALGR).cut(box(D, L - D, D, B - D, H_R - 1, H_R + H_BALGR + 1))
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
          ("Bassdeckel", (0.30, 0.22, 0.15)), ("Ausschnitt_", (0.62, 0.45, 0.28)), ("HB_Kammer_liegend", (0.65, 0.75, 0.88)), ("Randleiste", (0.62, 0.45, 0.28)), ("Schieberhebel_Lagerbuchse", (0.80, 0.65, 0.30)), ("Schieberhebel", (0.75, 0.76, 0.80)), ("Deckelschieber", (0.20, 0.15, 0.10)), ("HB_Langteil", (0.65, 0.75, 0.88)), ("HB_Vierkantleiste", (0.55, 0.40, 0.25)), ("HB_Trennfeder", (0.85, 0.30, 0.30)), ("HB_Haken", (0.40, 0.40, 0.42)), ("HB_Schraube", (0.40, 0.40, 0.42)), ("Balgrahmen", (0.55, 0.35, 0.20)), ("Deckelauflage", (0.60, 0.45, 0.28)), ("Eckblock", (0.55, 0.40, 0.25)), ("Eckwinkel_Schraube", (0.40, 0.40, 0.42)), ("Eckwinkel", (0.75, 0.76, 0.80)), ("Abdeckplatte", (0.55, 0.40, 0.25)), ("Abschlusswand", (0.55, 0.40, 0.25)), ("Innenleiste", (0.55, 0.40, 0.25)),
          ("Klappe_", None), ("HB_Folie", (0.85, 0.30, 0.30)), ("HB_Stimmplatte", (0.80, 0.82, 0.86)), ("HB_Trennwand", (0.50, 0.62, 0.80)), ("HB_", (0.65, 0.75, 0.88)),
          ("Begleiter_Plattfeder_Schraube", (0.40, 0.40, 0.42)), ("Begleiter_Haken", (0.40, 0.40, 0.42)), ("Begleiter_Plattfeder", (0.35, 0.45, 0.70)), ("Begleiter_Begrenzung", (0.60, 0.45, 0.28)), ("Begleiter_Stimmplatte", (0.80, 0.82, 0.86)), ("Begleiter_Sole", (0.55, 0.68, 0.45)), ("Begleiter_Schraube", (0.40, 0.40, 0.42)), ("Begleiter", (0.62, 0.72, 0.55)), ("Haken", (0.40, 0.40, 0.42)), ("Formteil_Schraube", (0.40, 0.40, 0.42)), ("Formteil", (0.60, 0.45, 0.28)), ("Dichtung", (0.12, 0.08, 0.05)), ("Lufttaste_Filz", (0.85, 0.15, 0.15)), ("Lufttaste_Gabelgelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)),
          ("Lufttaste_Schubstange", (0.10, 0.60, 0.60)), ("Lufttaste", (0.35, 0.25, 0.15)), ("Lufthebel_Feder", (0.75, 0.10, 0.75)), ("Lufthebel_Stuetze_Schraube", (0.40, 0.40, 0.42)),
          ("Lufthebel_Lagerbock_Schraube", (0.40, 0.40, 0.42)), ("Lufthebel_Lagerbock", (0.25, 0.75, 0.35)), ("Lufthebel_Achse", (0.55, 0.55, 0.60)), ("Begleiter_Schieber", (0.85, 0.30, 0.30)), ("Lufthebel_Stuetze", (0.25, 0.75, 0.35)), ("Lufthebel_Bolzen", (0.95, 0.85, 0.20)), ("Lufthebel", (0.70, 0.72, 0.78)), ("Luftklappe_Leder", (0.45, 0.30, 0.20)),
          ("Luftklappe_Filz", (0.85, 0.85, 0.80)), ("Luftklappe", (0.20, 0.20, 0.22))]
def farbe(n):
    if n.endswith("_gedrueckt"): return (0.95, 0.55, 0.55)
    if n.startswith(("Klappe_", "Begleiterklappe_")): return (0.45, 0.30, 0.20) if "Leder" in n else (0.85, 0.85, 0.80) if "Filz" in n else (0.20, 0.20, 0.22)
    for p_, c_ in FARBEN:
        if n.startswith(p_) and c_: return c_
    print("FARBE ohne Zuordnung:", n); return (0.6, 0.6, 0.6)
exec(open("mech_gerhard.py").read())             # Mechanik (vor der Spiegelung)
teile.append(("Balgrahmen", balgrahmen))          # in Fusion ausblenden
NUR_LUFT = len(sys.argv) > 3 and sys.argv[3] == "luft"
if NUR_LUFT:                                   # nur Hebeleinbau zeigen: Deckel, Boden, Eckbloecke, Luftgruppe in Ruhe und gedrueckt
    teile = [(n, s_) for n, s_ in teile if n.startswith(("Bassboden", "Bassdeckel", "Eckblock", "Eckwinkel", "Fuss_Schraube_3", "Deckelauflage", "Luft", "Front", "Knopf")) and not n.startswith("Bassbodenleder")] + lg_g
a = cq.Assembly(name="Bassteil_Gerhard")
SPIEGELN = True                                # wie Markus: links/rechts getauscht (X = L/2)
for n, s_ in teile:
    if SPIEGELN: s_ = s_.mirror("YZ", basePointVector=(L / 2, 0, 0))
    a.add(s_, name=n, color=cq.Color(*farbe(n)))
a.save(OUT)
print("gespeichert", OUT, len(teile), "Teile")
