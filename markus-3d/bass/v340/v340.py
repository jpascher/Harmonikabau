# Bassteil "340 schmal leicht schraeg" (Zeichnung 340-schmal-leicht-schraeg.dxf): Gehaeuse 340 x 205 (Rohmass 342 x 207),
# Knopfreihen 9 / 9 (Teilung 19), 6 Hoerner, 9 schmale Helikonkammern (46 x 94) in zwei Reihen Ruecken an Ruecken mit leicht
# schraeg liegender Platte (aussen Z 89, innen Z 108,6), gerade Trennfeder; Begleiterstimmstoecke links/rechts wie Marie 38
# (je 4 Begleiter, rund O18); Luftklappe auf abgesenktem Bodenteil, oeffnet nach innen, in normalen Lagerboecken; Taste mit
# Achse und zwei Schenkeln (Fuehrung), direkte Schubstange Taste -> Klappe, alle Gelenke normale Gabelgelenke.
# Aufgebaut auf dem Markus-Bassteil (Rahmen gekuerzt). Achsen: X = Laenge, Y = Tiefe (0 = Rueckteil aussen), Z ab Deckel aussen.
# Knopfreihen (26.11.2026): Knopfzahl je Reihe frei (Kopf R1 = n / R2 = n der Belegung in BASS_BELEGUNG), jede Reihe mittig um
# die gezeichnete Reihenmitte mit der Knopfteilung verteilt; Knoepfe ohne Zuordnung bleiben blind (ohne Mechanik).
# Aufruf: python3 v340.py ausgabe.step [Knopfteilung]
import sys, os, math, ezdxf, ezdxf.path, cadquery as cq
from ezdxf import bbox
OUT = sys.argv[1]
T_KNOPF = float(sys.argv[2]) if len(sys.argv) > 2 else 19.0
D, L, B = 5.5, 342.0, 207.0
DL = L - 382.0                                 # Kuerzung gegenueber Markus (Rahmen an zwei Stellen gekuerzt)
OX = OY = 1.0
T_DECKEL_BASS = 4.0
Z_BU, Z_BO = 34.5, 38.5                       # Boden 4 mm (Schnitt 95,8..99,8 ab Deckel aussen 61,3)
FI = B - D
def box(x0, x1, y0, y1, z0, z1):
    return cq.Solid.makeBox(x1 - x0, y1 - y0, z1 - z0, cq.Vector(x0, y0, z0))
def py(x0, x1, pts):                          # Profil (Y, Z) ueber X extrudiert
    return cq.Workplane("YZ").polyline(pts).close().extrude(x1 - x0).val().translate(cq.Vector(x0, 0, 0))
teile = []
imp = lambda f: cq.importers.importStep(f).val()
O, XA, ZA = cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), cq.Vector(0, 0, 1)
st = lambda s: s.rotate(O, XA, 90)
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v_ = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v_.Length, p0, v_.normalized())
def kuerze(s_, schnitte):                       # an jeder Stelle (xc, d) das Stueck xc..xc+d herausnehmen (Querschnitt dort gleich)
    bb_ = s_.BoundingBox(); gr_ = lambda a_, b_: box(a_, b_, bb_.ymin - 1, bb_.ymax + 1, bb_.zmin - 1, bb_.zmax + 1)
    stuecke = []; x_alt = bb_.xmin - 1; off_ = 0.0
    for xc_, d_ in schnitte:
        stuecke.append(s_.intersect(gr_(x_alt, xc_)).translate(cq.Vector(-off_, 0, 0))); x_alt = xc_ + d_; off_ += d_
    stuecke.append(s_.intersect(gr_(x_alt, bb_.xmax + 1)).translate(cq.Vector(-off_, 0, 0)))
    r_ = stuecke[0]
    for t_ in stuecke[1:]: r_ = r_.fuse(t_)
    return r_.clean()
SCHN = [(110.0, -DL / 2), (252.0, -DL / 2)]     # im Markus-Modell (382), Querschnitt nach dem Fuellen der Front gleich
def m_k(X382):                                   # Markus-Modell-X -> gekuerztes X
    return X382 - sum(min(max(X382 - xc_, 0.0), d_) for xc_, d_ in SCHN)
m_r = m_k
RD = ezdxf.readfile("rahmen.dxf").modelspace()
L0 = 382.0
front0 = st(imp("Markus_front_v4.step")).rotate(O, ZA, 180).translate(cq.Vector(L0, B - D, 0))
H_R = front0.BoundingBox().zmax
_wl = [e for e in RD.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]') if not e.closed][0]
_wp = [(q.x, q.y) for q in ezdxf.path.make_path(_wl).flattening(0.05)]
_wp = _wp + [(_wp[-1][0], -4.24 - 1), (_wp[0][0], -4.24 - 1)]
front0 = front0.fuse(cq.Workplane("XZ").polyline([(-x, y + 4.24) for x, y in _wp]).close().extrude(-(D - 0.04)).translate((0, B - D + 0.02, 0)).val().intersect(box(0, L0, B - D, B, 0, H_R)))
for lay_ in ("bohrungen-bassdrücker", "Hörner-Bohrung"):                     # Markus-Bohrungen zu
    for c in RD.query(f'CIRCLE[layer=="{lay_}"]'):
        x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
        if -382 < x < 0 and -4.24 < y < 71.75:
            front0 = front0.fuse(stab((-x, B - D + 0.02, y + 4.24), (-x, B - 0.02, y + 4.24), 2 * r + 0.2))
front = kuerze(front0.clean(), SCHN)
rueck = kuerze(st(imp("Markus_rueck_v4.step")).translate(cq.Vector(0, D, 0)), SCHN)
seite = imp("Markus_seite_v4.step")
sl = st(seite).rotate(O, ZA, -90).translate(cq.Vector(D, B, 0))
sr = st(seite).rotate(O, ZA, 90).translate(cq.Vector(L - D, 0, 0))
print("Rahmen gekuerzt: Front X", round(front.BoundingBox().xmin, 2), "-", round(front.BoundingBox().xmax, 2), "Rueck", round(rueck.BoundingBox().xmax, 2))
# ---------- Zeichnung: Draufsicht X = x + 354,8, Y = y + 243,3 (Rohmass +1); Frontansicht Z = y + 16,1 ----------
AD = ezdxf.readfile("v340.dxf").modelspace()
PX = lambda x_: x_ + 354.8
PY = lambda y_: y_ + 243.3
# ---------- Front: Knopfreihen 9 / 9 wie gezeichnet (Reihe 1 buendig, Reihe 2 auf 7-mm-Blende), Hoerner nach Zeichnung ----------
_kc = [c for c in AD.query('CIRCLE[layer=="Front"]') if -16 < c.dxf.center.y < 60 and abs(c.dxf.radius - 8.5) > 0.1 and abs(c.dxf.radius - 16.5) > 0.1]
_R1Z = sorted({round(PX(c.dxf.center.x), 2) for c in _kc if abs(c.dxf.radius - 5.5) < 0.05})   # gezeichnet (9 / 9, Teilung 19)
_R2Z = sorted({round(PX(c.dxf.center.x), 2) for c in _kc if abs(c.dxf.radius - 5.0) < 0.05})
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath("v340.py"))))
import json as _json, belegung as _belegung
_BMJ = _json.load(open("belegung_modell_340.json", encoding="utf-8"))
REIHEN_N = {1: len(_R1Z), 2: len(_R2Z)}
if os.environ.get("BASS_BELEGUNG"): REIHEN_N = _belegung.reihen_aus(_belegung.lies(os.environ["BASS_BELEGUNG"])[0], REIHEN_N)
_KX = _belegung.knopf_lagen(_BMJ, REIHEN_N, T_KNOPF)          # mittig um die gezeichnete Reihenmitte, prueft den Platz
R1X, R2X = _KX[1], _KX[2]
Z_RH = {1: -2.3 + 16.1, 2: 11.7 + 16.1}
R_RH = {1: 5.5, 2: 5.0}
X_KM = L / 2
kn = [(X_, Z_RH[1], R_RH[1], 1) for X_ in R1X] + [(X_, Z_RH[2], R_RH[2], 2) for X_ in R2X]
def knx(r_, i_): return (R1X if r_ == 1 else R2X)[i_]
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
horn_pos = sorted({(round(PX(c.dxf.center.x), 2), round(c.dxf.center.y + 16.1, 2)) for c in AD.query('CIRCLE[layer=="Front"]')
                   if abs(c.dxf.radius - 8.5) < 0.05 and -16 < c.dxf.center.y < 60})
for X_, Z_ in horn_pos:
        r = 9.0
        front = front.cut(stab((X_, FI - 1, Z_), (X_, B + 1, Z_), 2 * r))
        ring = cq.Workplane("XZ").center(X_, Z_).circle(r + 2.5).circle(r).extrude(-40).translate((0, B - 30, 0)).val()
        y_out = front.intersect(ring).BoundingBox().ymax
        w = cq.Workplane("XY").polyline([(r_, a_) for a_, r_ in horn_prof]).close().revolve(360, (0, 0, 0), (0, 1, 0))
        teile.append((f"Horn_{len(aus_front)+1}", w.val().translate(cq.Vector(X_, y_out - 4.01, Z_))))
        aus_front.append(cq.Workplane("XZ").center(X_, Z_).circle(r + 5.0).extrude(-(6.0 + 2)).translate((0, FI - 6.0 - 1, 0)).val())
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
def T_BLR(r_): return T_BL * (r_ - 1)
bl = cq.Workplane("XZ").polyline([(-x_, Z_RH[2] + y_) for x_, y_ in blende_kontur([-x_ for x_ in R2X])]).close().extrude(-T_BL).val()
bl = bl.translate(cq.Vector(0, B - bl.BoundingBox().ymin, 0))
for X_, Z_, r_, rr in kn:
    if rr == 2: bl = bl.cut(stab((X_, B - 1, Z_), (X_, B + T_BL + 1, Z_), 2 * r_))
teile.append(("Frontblende", bl))
LT_ = [e for e in AD.query('LWPOLYLINE[layer=="Def"]')][0]
_lx = [PX(p[0]) for p in LT_.get_points()]; _ly = [PY(p[1]) for p in LT_.get_points()]
LT_X, LT_Y, LT_L, LT_B = (min(_lx) + max(_lx)) / 2, (min(_ly) + max(_ly)) / 2, max(_lx) - min(_lx), max(_ly) - min(_ly)   # Def-Rechteck der Lufttaste
for X_, Z_, r_, rr in kn:
    ya = B + T_BLR(rr); dl_ = T_BLR(rr) if rr > 1 else D; yk = ya + SPALT; lz = yk - (B + SPALT - L_ZYL)
    if rr == 1 and LT_X - LT_L / 2 - D_ZYL / 2 - 1 < X_ < LT_X + LT_L / 2 + D_ZYL / 2 + 1: lz -= 7.0   # ueber der Lufttaste kuerzer
    teile.append((f"Knopf_R{rr}_{X_:.0f}", stab((X_, yk, Z_), (X_, yk + H_KOPF, Z_), D_KOPF)))
    teile.append((f"Knopfzylinder_R{rr}_{X_:.0f}", stab((X_, yk - lz, Z_), (X_, yk, Z_), D_ZYL)))
    teile.append((f"Filzeinlage_R{rr}_{X_:.0f}", stab((X_, ya - dl_, Z_), (X_, ya, Z_), 2 * r_).cut(stab((X_, ya - dl_ - 1, Z_), (X_, ya + 1, Z_), D_ZYL))))
print("Front: Reihen", len(R1X), "/", len(R2X), "R1", R1X[0], "-", R1X[-1], "R2", R2X[0], "-", R2X[-1], "| Hoerner", len(aus_front), horn_pos, "| Lufttaste", round(LT_X, 1), round(LT_Y, 1), round(LT_L, 1), "x", round(LT_B, 1))
# Innenleisten 4 mm vorne/hinten vom Boden bis Z 89 (Schnitt), Anschlag der Kammer-Aussenwaende
Z_LST = 89.0
afo = box(D, L - D, FI - 4.0, FI, Z_BO, Z_LST)
for b_ in knopf_bohr: afo = afo.cut(b_)
teile.append(("Innenleiste_vorne", afo))
teile.append(("Innenleiste_hinten", box(D, L - D, D, D + 4.0, Z_BO, Z_LST)))

# ---------- Deckelauflage wie Markus (6 mm, vom Deckel bis zur Bodenunterseite), Aussparungen fuer Gaze am Rueckteil ----------
B_AUF, RAND_G, REST = 6.0, 5.0, 4.0
za0 = T_DECKEL_BASS
aus_rueck = []
for e in RD.query('SPLINE[layer=="LAYER_1"]'):
    pp = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.1)]
    if min(b for _, b in pp) > -4.24:
        xs = [m_r(a + 382.0) for a, _ in pp]
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
_ros = sorted({(round(c.dxf.center.x + 381 + DL / 2, 2), round(c.dxf.center.y + 206, 2), round(c.dxf.radius, 2)) for c in VD.query('CIRCLE[layer=="0"]')})   # mittig (+DL/2)
for (x_, y_, r_) in _ros: deckel_r = deckel_r.cut(cq.Solid.makeCylinder(r_, T_DECKEL_BASS + 2, cq.Vector(x_, y_, -1)))
# ---------- Boden 4 mm ueber die ganze Innenflaeche; Oeffnungen und Klappen aus der Zeichnung ----------
boden = box(D, L - D, D, FI, Z_BU, Z_BO)
oeff, luft_oeff = [], None
for e in AD.query('LWPOLYLINE[layer=="Füllung-Rund-groß-Begleiter"]'):
    xs = [PX(p[0]) for p in e.get_points()]; ys = [PY(p[1]) for p in e.get_points()]
    if not (0 < min(xs) and max(xs) < L and 0 < min(ys) and max(ys) < B): continue
    r_ = (round(min(xs), 2), round(max(xs), 2), round(min(ys), 2), round(max(ys), 2))
    if r_[3] - r_[2] > 50: luft_oeff = r_                  # Luftoeffnung 23 x 69
    else: oeff.append(r_)
oeff = sorted(set(oeff), key=lambda k: (-k[2], k[0]))
RAND_V, TIEF_V, SPIEL = 3.5, 2.0, 0.5                  # Klappenrand 3 (Zeichnung 44,4 x 30)
klappen = []
for i, (x0, x1, y0, y1) in enumerate(oeff):
    boden = boden.cut(box(x0, x1, y0, y1, Z_BU - 1, Z_BO + 1)).cut(box(x0 - RAND_V, x1 + RAND_V, y0 - RAND_V, y1 + RAND_V, Z_BU - 1, Z_BU + TIEF_V))
    z_ = Z_BU + TIEF_V; o_ = RAND_V - SPIEL
    for nm, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
        klappen.append((f"Klappe_{i+1}_{nm}", box(x0 - o_, x1 + o_, y0 - o_, y1 + o_, z_ - dk, z_))); z_ -= dk
# Begleiter: runde Oeffnungen O18, runde Klappen O24 (Klappenrand-rund), Vertiefung O25 x 2
nb = 0; bg_kr = []
for c in sorted(AD.query('CIRCLE[layer=="Füllung-Rund-groß-Begleiter"]'), key=lambda c: (c.dxf.center.x, c.dxf.center.y)):
    cx_, cy_ = PX(c.dxf.center.x), PY(c.dxf.center.y)
    if not (0 < cx_ < L and 0 < cy_ < B): continue
    rk_ = 12.0
    boden = boden.cut(cq.Solid.makeCylinder(c.dxf.radius, 10, cq.Vector(cx_, cy_, Z_BU - 1))).cut(cq.Solid.makeCylinder(rk_ + SPIEL, TIEF_V + 1, cq.Vector(cx_, cy_, Z_BU - 1)))
    nb += 1; z_ = Z_BU + TIEF_V; bg_kr.append((round(cx_, 2), round(cy_, 2)))
    for nm, dk in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
        klappen.append((f"Begleiterklappe_{nb}_{nm}", cq.Solid.makeCylinder(rk_, dk, cq.Vector(cx_, cy_, z_ - dk)))); z_ -= dk
kl_ = oeff
print("Boden: Helikonoeffnungen", len(oeff), "| Begleiter rund", nb, bg_kr, "| Luftoeffnung", luft_oeff)
# ---------- Helikonkammern: zwei Reihen Ruecken an Ruecken, Teilung 47 (Spalten ueber den Klappen), Platte 46 x 94 leicht schraeg
# (aussen Z 89, innen Z 108,6 wie Schnitt), Trennwaende 3,5, Stirnwaende 3, gerade Trennfeder 1 mm auf Z 58,5 von innen 71 lang ----------
PL_L, PL_B, RST = 94.0, 46.0, 1.5
T_RAND, T_DUENN, T_DICK = 2.0, 5.0, 16.0
TK, T_TW, T_STA, T_STI = 47.0, 3.5, 2.0, 2.5      # Stirnwand aussen (an der Innenleiste) / innen
XS_K = sorted({round((o[0] + o[1]) / 2, 2) for o in oeff})              # Spaltenmitten
reihe_h = sorted([o for o in oeff if o[2] < B / 2], key=lambda o: o[0]); reihe_v = sorted([o for o in oeff if o[2] > B / 2], key=lambda o: o[0])
Z_AUS, Z_INN = 89.0, 108.6
Y_H0, Y_H1 = D + 4.0, B / 2                    # hintere Reihe aussen .. innen (Ruecken an Ruecken)
Y_V0, Y_V1 = B / 2, FI - 4.0                   # vordere Reihe innen .. aussen
Z_FED, L_FED, T_FED = 58.5, 71.0, 1.0
def z_ob(y_, y_a, y_i):                        # Oberkante schraeg: aussen Z_AUS, innen Z_INN
    return Z_AUS + (Z_INN - Z_AUS) * (y_ - y_a) / (y_i - y_a)
from math import atan2, degrees
def stimmplatte0(px0, px1, l_):                # flach: Y 0..l_, Z 0 = Unterseite Rand; dickes Ende bei Y = l_
    pl = box(px0, px1, 0, l_, 0, T_RAND)
    a0, a1 = RST, l_ - RST
    pl = pl.fuse(py(px0 + RST, px1 - RST, [(a0, 0.01), (a1, 0.01), (a1, T_RAND - T_DICK), (a0, T_RAND - T_DUENN)]))
    for (o_, sb, sll) in [(3.8, 7.5, 51.0), (12.8, 9.2, 71.5), (24.0, 9.2, 71.5), (34.7, 7.5, 51.0)]:   # schmale Markus-Platte 46
        pl = pl.cut(box(px0 + o_, px0 + o_ + sb, l_ - 5.5 - sll, l_ - 5.5, -T_DICK - 1, T_RAND + 1))
    return pl
D_HAKEN = 1.8
def zyl(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1)
    return cq.Solid.makeCylinder(d / 2, (p1 - p0).Length, p0, (p1 - p0).normalized())
hb_teile, haken_l = [], []
for nm_r, y_a, y_i, spalten in (("hinten", Y_H0, Y_H1, reihe_h), ("vorne", Y_V1, Y_V0, reihe_v)):
    sg_ = 1 if y_i > y_a else -1                 # Richtung aussen -> innen
    xs_ = [(o[0] + o[1]) / 2 for o in spalten]
    xw_ = [x_ - TK / 2 for x_ in xs_] + [xs_[-1] + TK / 2]          # Trennwand-Mitten
    X0_, X1_ = xw_[0] - T_TW / 2, xw_[-1] + T_TW / 2
    ya_, yb_ = sorted((y_a, y_i))
    prof = [(y_a, Z_BO), (y_a, Z_AUS), (y_i, Z_INN), (y_i, Z_BO)]
    kasten = py(X0_, X1_, prof if sg_ > 0 else prof[::-1])
    for k_ in range(len(xs_)):                  # Kammern hohl (unten offen auf den Boden), Trennwaende bleiben
        kasten = kasten.cut(box(xw_[k_] + T_TW / 2, xw_[k_ + 1] - T_TW / 2, *sorted((y_a + sg_ * T_STA, y_i - sg_ * T_STI)), Z_BO - 1, Z_INN + 5))
    yf0_, yf1_ = sorted((y_i - sg_ * T_STI, y_i - sg_ * (T_STI + L_FED)))
    feder = box(X0_ + T_TW / 2, X1_ - T_TW / 2, yf0_, yf1_, Z_FED, Z_FED + T_FED)
    hb_teile += [(f"HB_Kammern_{nm_r}", kasten.cut(feder)), (f"HB_Trennfeder_{nm_r}", feder)]
    ang = degrees(atan2(Z_INN - Z_AUS, abs(y_i - y_a)))             # Neigung der Platte
    l_s = math.hypot(y_i - y_a, Z_INN - Z_AUS)
    for k_, x_ in enumerate(xs_):
        p0 = stimmplatte0(x_ - PL_B / 2, x_ + PL_B / 2, PL_L)     # flach, dickes Ende innen
        y_off = (l_s - PL_L) / 2
        p0 = p0.translate(cq.Vector(0, y_off, 0))
        if sg_ > 0: p0 = p0.rotate(O, cq.Vector(1, 0, 0), ang).translate(cq.Vector(0, y_a, Z_AUS))
        else: p0 = p0.mirror("XZ").rotate(O, cq.Vector(1, 0, 0), -ang).translate(cq.Vector(0, y_a, Z_AUS))
        hb_teile.append((f"HB_Stimmplatte_{nm_r}_{k_+1}", p0))
        # Haken: senkrecht in die Stirnwaende, kurzer Schenkel auf dem Plattenrand (schraeg)
        for (y_s, z_s, dr_) in ((y_a + sg_ * (y_off - 0.1 - D_HAKEN / 2) * math.cos(math.radians(ang)), Z_AUS, sg_), (y_i - sg_ * (y_off - 0.1 - D_HAKEN / 2) * math.cos(math.radians(ang)), Z_INN, -sg_)):
            zh_ = z_s + T_RAND + D_HAKEN / 2 + 0.2
            v_ = cq.Vector(0, dr_ * math.cos(math.radians(ang)) * 5.0, (1 if dr_ == sg_ else -1) * math.sin(math.radians(ang)) * 5.0)
            xh_ = x_ + (-8.0 if sg_ > 0 else 8.0)
            hk_ = zyl((xh_, y_s, zh_ - 15.0), (xh_, y_s, zh_), D_HAKEN).fuse(zyl((xh_, y_s, zh_), (xh_, y_s + v_.y, zh_ + v_.z), D_HAKEN)).fuse(cq.Solid.makeSphere(D_HAKEN / 2, cq.Vector(xh_, y_s, zh_)))
            haken_l.append((f"Haken_{nm_r}_{k_+1}_{'aussen' if z_s == Z_AUS else 'innen'}", hk_))
_hk = cq.Compound.makeCompound([s_ for _, s_ in haken_l])
hb_teile = [(n_, s_.cut(_hk) if n_.startswith("HB_Kammern") else s_) for n_, s_ in hb_teile]
teile = [(n_, s_.cut(_hk) if n_.startswith("Innenleiste") else s_) for n_, s_ in teile]
teile += hb_teile + haken_l
print("Kammern: Spalten", [round(x, 1) for x in XS_K], "| hinten", len(reihe_h), "vorne", len(reihe_v), f"| Neigung {ang:.1f} Grad")

# ---------- Begleiterstimmstoecke links/rechts wie Marie 38 (je 4 Begleiter, 8 x 2 Kammern, Sole 33 x 146 nach Zeichnung) ----------
T_SOLE, B_ST, H_ST, K_H, K_T = 6.0, 20.0, 50.0, 45.0, 9.5
T_SPL, B_SPL, H_SPL = 2.5, 16.8, 47.0
_sol = sorted([e for e in AD.query('LWPOLYLINE[layer=="0"]') if e.closed and 30 < bbox.extents([e]).size.x < 36 and bbox.extents([e]).size.y > 140], key=lambda e: bbox.extents([e]).extmin.x)
for nm, e_ in zip(("links", "rechts"), _sol):
    b_ = bbox.extents([e_]); x0 = PX(b_.extmin.x); ys0 = PY(b_.extmin.y)          # Sole 33 x 146 (Zeichnung)
    sy0, sy1 = ys0 - 6.0, ys0 + 146.0 + 6.0                                        # an beiden Enden 6 mm laenger: Schraube in den Boden
    xm = x0 + 16.5; bx0, bx1 = xm - B_ST / 2, xm + B_ST / 2
    sole = box(x0, x0 + 33, sy0, sy1, Z_BO, Z_BO + T_SOLE)
    blk = box(bx0, bx1, ys0 + 4.0, ys0 + 142.0, Z_BO + T_SOLE, Z_BO + T_SOLE + H_ST)
    plt_ = []
    for j in range(8):
        ky0 = ys0 + 6.0 + 17 * j; ky1 = ky0 + 15
        for (xa_, xb_, sg) in ((bx0 - 1, bx0 + K_T, -1), (bx1 - K_T, bx1 + 1, 1)):
            blk = blk.cut(box(xa_, xb_, ky0, ky1, Z_BO + T_SOLE - 1, Z_BO + T_SOLE + K_H))
            sole = sole.cut(box(max(xa_, bx0), min(xb_, bx1), ky0, ky1, Z_BO - 1, Z_BO + T_SOLE + 1))
            xp = bx0 - T_SPL if sg < 0 else bx1
            yc = (ky0 + ky1) / 2
            zb_ = Z_BO + T_SOLE + 0.5
            spl_ = box(xp, xp + T_SPL, yc - B_SPL / 2, yc + B_SPL / 2, zb_, zb_ + H_SPL)
            for dy_ in (-3.5, 3.5):                            # zwei Zungenschlitze 2,6 x 34, Nietloch O1,6
                spl_ = spl_.cut(box(xp - 1, xp + T_SPL + 1, yc + dy_ - 1.3, yc + dy_ + 1.3, zb_ + 6.0, zb_ + 40.0))
                spl_ = spl_.cut(stab((xp - 1, yc + dy_, zb_ + 3.0), (xp + T_SPL + 1, yc + dy_, zb_ + 3.0), 1.6))
            plt_.append(spl_)
    for k_, ys_ in enumerate((sy0 + 3.0, sy1 - 3.0)):
        sole = sole.cut(cq.Solid.makeCylinder(1.75, T_SOLE + 2, cq.Vector(xm, ys_, Z_BO - 1)))
        sc_ = cq.Solid.makeCylinder(3.0, 2.0, cq.Vector(xm, ys_, Z_BO + T_SOLE)).fuse(cq.Solid.makeCylinder(1.5, T_SOLE + 3.0, cq.Vector(xm, ys_, Z_BO - 3.0)))
        boden = boden.cut(cq.Solid.makeCylinder(1.25, 3.0, cq.Vector(xm, ys_, Z_BO - 3.0)))
        teile.append((f"Begleiter_Schraube_{nm}_{k_+1}", sc_))
    teile.append((f"Begleiter_Sole_{nm}", sole)); teile.append((f"Begleiterstimmstock_{nm}", blk))
    for q_, p_ in enumerate(plt_): teile.append((f"Begleiter_Stimmplatte_{nm}_{q_+1}", p_))
    print(f"Begleiterstimmstock {nm}: Sole X {x0:.1f}-{x0 + 33:.1f}, Y {ys0:.1f}-{ys0 + 146:.1f}")
# ---------- Luft: abgesenkter Bodenteil (4 mm, Z 69..73) vorne rechts, ringsum 2 mm in Nuten (lange Trennwand, neue Seitenwand,
# hintere Kammerreihe, vordere Innenleiste); im Hauptboden darunter ein Ausschnitt so gross wie die fehlende Basskammer.
# Luftklappe oben auf dem Bodenteil, an der Laengsseite gelagert (Achse in zwei Lagerboecken, zwei Schenkelfedern wie Markus),
# oeffnet nach innen in die Vertiefung. Lufttaste: Achse direkt unter dem Boden (Lagerboecke ohne Abstand), drei Schenkel:
# zwei zur Taste (vorne und hinten an der Taste, je Stange mit Gelenk oben und unten), einer schraeg hoch durch den Ausschnitt,
# Gelenk und Schubstange zur Klappe. ----------
lg = []
T_FILZ_LT, HUB_LT, RAND_LK = 1.0, 10.0, 3.0
def langloch_m(l_, b_, z0_, h_):
    return cq.Workplane("XY").center(LT_X, LT_Y).slot2D(l_, b_).extrude(h_).translate((0, 0, z0_)).val()
lt_filz = langloch_m(LT_L, LT_B, 0, T_DECKEL_BASS).cut(langloch_m(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT, -1, T_DECKEL_BASS + 2))
taste = cq.Workplane("XY").center(LT_X, LT_Y).slot2D(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT).extrude(T_DECKEL_BASS + HUB_LT).faces("<Z").edges().fillet(1.0).val()
taste = taste.translate(cq.Vector(0, 0, -HUB_LT))
Z_LB0, Z_LB1, T_NUT_L = 69.0, 73.0, 2.0
x0_, x1_, y0_, y1_ = luft_oeff
XW_L = reihe_v[-1][1] + (TK - (reihe_v[-1][1] - reihe_v[-1][0])) / 2
XL0 = XW_L + T_TW / 2
XS0 = x1_ + 13.0; XS1 = XS0 + 3.5
YL0, YL1 = B / 2, FI - 4.0
sw = box(XS0, XS1, YL0, YL1, Z_BO, Z_LB1)                  # Seitenwand nur bis zum Klappenboden
bt = box(XL0 - T_NUT_L, XS0 + T_NUT_L, YL0 - T_NUT_L, YL1 + T_NUT_L, Z_LB0, Z_LB1).cut(box(x0_, x1_, y0_, y1_, Z_LB0 - 1, Z_LB1 + 1))
teile = [(n_, s_.cut(bt) if n_ in ("HB_Kammern_vorne", "HB_Kammern_hinten", "Innenleiste_vorne") else s_) for n_, s_ in teile]
# Hauptboden: Ausschnitt wie die fehlende Basskammer (Innenmass der Spalte)
XA0, XA1 = XL0, XL0 + (TK - T_TW)
YA0, YA1 = YL0 + T_STI, YL1 - T_STA
boden = boden.cut(box(XA0, XA1, YA0, YA1, Z_BU - 1, Z_BO + 1))
# Klappe, Laengsseite (rechts, zur Seitenwand) gelagert
k0, k1, k2, k3 = x0_ - RAND_LK, x1_ + RAND_LK, y0_ - RAND_LK, y1_ + RAND_LK
lk_ = [("Luftklappe_Leder", box(k0, k1, k2, k3, Z_LB1, Z_LB1 + 1.0)), ("Luftklappe_Filz", box(k0, k1, k2, k3, Z_LB1 + 1.0, Z_LB1 + 5.0)),
       ("Luftklappe", box(k0, k1, k2, k3, Z_LB1 + 5.0, Z_LB1 + 6.5))]
X_H, Z_H, D_H = k1 + 3.0, Z_LB1 + 2.5, 2.4                # Klappenachse (Laengsseite), Hoehe wie normale Lagerboecke
for j_, (ya_, yb_) in enumerate(((k2 + 6.0, k2 + 14.0), (k3 - 14.0, k3 - 6.0))):   # Laschen an der Klappe um die Achse
    lk_.append((f"Luftklappe_Lasche_{j_+1}", box(k1, X_H + 1.8, ya_, yb_, Z_LB1 + 5.0, Z_LB1 + 6.5).fuse(stab((X_H, ya_, Z_H), (X_H, yb_, Z_H), 3.6)).cut(stab((X_H, ya_ - 1, Z_H), (X_H, yb_ + 1, Z_H), D_H + 0.1))))
lg += lk_
ach_k = stab((X_H, k2 - 4.6 + 1.0, Z_H), (X_H, k3 + 4.6 - 1.0, Z_H), D_H)   # Enden 1 mm vor den Schrauben
lg.append(("Luftklappe_Achse", ach_k))
_LBD = ezdxf.readfile("lagerbeide.dxf").modelspace()            # normale Lagerboecke (wie Mechanik), stehend auf dem Klappenboden
_pe = [e for e in _LBD.query('LWPOLYLINE[layer=="Layer1"]') if e.closed][0]
_pp = [(q.x, q.y) for q in ezdxf.path.make_path(_pe).flattening(0.02)]
_pp = _pp[:-1] if math.dist(_pp[0], _pp[-1]) < 1e-6 else _pp
for j_, yb_ in enumerate((k2 - 4.6, k3 + 4.6)):
    xz_ = [(X_H + (-q[1] - 2.6), Z_LB1 + (q[0] + 13.6)) for q in _pp]
    b_ = cq.Workplane("XZ").polyline(xz_).close().extrude(-8.6).val().translate(cq.Vector(0, yb_ - 4.3, 0))
    b_ = b_.cut(stab((X_H, yb_ - 5, Z_H), (X_H, yb_ + 5, Z_H), D_H + 0.1)).cut(stab((X_H, yb_, Z_LB1 - 1), (X_H, yb_, Z_LB1 + 6), 1.5))
    sc_ = stab((X_H, yb_, Z_LB1 + 5.0 + 0.9), (X_H, yb_, Z_LB1 + 5.0), 2.5).fuse(stab((X_H, yb_, Z_LB1 + 5.0), (X_H, yb_, Z_LB0 + 1.0), 1.4))
    lg += [(f"Luftklappe_Lagerbock_{j_+1}", b_), (f"Luftklappe_Lagerbock_Schraube_{j_+1}", sc_)]
# Schenkelfedern (3 Windungen) auf der Klappenachse neben den Lagerboecken, ein Schenkel auf der Klappe, einer auf dem Bodenteil
DR_K, RW_K, NW_K, PW_K = 0.5, D_H / 2 + 0.25 + 0.25, 3, 0.65
def draht_k(punkte):
    pts_ = [cq.Vector(*p) for p in punkte]
    prof_ = cq.Wire.makeCircle(DR_K / 2, pts_[0], (pts_[1] - pts_[0]).normalized())
    return cq.Solid.sweep(prof_, [], cq.Wire.makePolygon(pts_), transitionMode="round")
for j_, (yf_, sdir) in enumerate(((k2 + 0.2, 1.0), (k3 - 0.2, -1.0))):
    lsp = NW_K * PW_K
    hel = cq.Wire.makeHelix(PW_K, lsp, RW_K)
    sp = cq.Workplane("XZ").center(RW_K, 0).circle(DR_K / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -90)                 # Achse der Wendel in Y
    if sdir < 0: sp = sp.mirror("XZ")
    sp = sp.translate(cq.Vector(X_H, yf_ - (0 if sdir > 0 else 0), Z_H))
    ye_ = yf_ + sdir * lsp
    s1 = draht_k([(X_H - RW_K, yf_, Z_H), (X_H - RW_K, yf_, Z_LB1 + 6.5 + DR_K / 2 + 0.05), (X_H - RW_K - 10.0, yf_ + sdir * 3.0, Z_LB1 + 6.5 + DR_K / 2 + 0.05)])
    s2 = draht_k([(X_H + RW_K, ye_, Z_H), (X_H + RW_K + 1.5, ye_, Z_LB1 + DR_K / 2 + 0.05), (X_H + RW_K + 6.0, ye_, Z_LB1 + DR_K / 2 + 0.05)])
    lg.append((f"Luftklappe_Feder_{j_+1}", sp.fuse(s1).fuse(s2).clean()))
# Achse der Lufttaste direkt unter dem Boden (Lagerboecke ohne Abstand), etwas laenger als die Taste
Z_LA = Z_BU - 2.5; D_LA = 2.4
Y_T = LT_Y
ANG0 = -20.0                                                # Tastenschenkel in Ruhe 20 Grad nach unten, gedrueckt waagrecht
L_S1 = HUB_LT / math.sin(math.radians(-ANG0))
Y_LA = Y_T - L_S1 * math.cos(math.radians(ANG0))
Z_G1 = Z_LA + L_S1 * math.sin(math.radians(ANG0))
_kxa = sorted(k[0] for k in kn)                             # vordere Schenkel in der Luecke zwischen den Knopfstangen bei X 209
_xs1 = min(((a_ + b_) / 2 for a_, b_ in zip(_kxa[:-1], _kxa[1:])), key=lambda m_: abs(m_ - 208.95))
X_S1 = [208.95 if abs(_xs1 - 208.95) < 0.1 else _xs1, LT_X + LT_L / 2 - 6.0]   # wie gezeichnet bei 9 / 9, Teilung 19
X_S3 = k0 + 6.0                                             # Schubstange Taste -> Klappe (gegenueber der Klappenachse)
XB_ = [X_S1[0] - 9.0, XA1 + 5.5]
ach = stab((XB_[0] + 1.0, Y_LA, Z_LA), (XB_[1] - 1.0, Y_LA, Z_LA), D_LA)   # Enden 1 mm vor den Schrauben
for j_, xb_ in enumerate(XB_):
    b_ = box(xb_ - 4.3, xb_ + 4.3, Y_LA - 2.0, Y_LA + 2.0, Z_BU - 5.0, Z_BU).cut(stab((xb_ - 5, Y_LA, Z_LA), (xb_ + 5, Y_LA, Z_LA), D_LA + 0.1))
    b_ = b_.cut(stab((xb_, Y_LA, Z_BU - 6), (xb_, Y_LA, Z_BU + 1), 1.5))
    sc_ = stab((xb_, Y_LA, Z_BU - 5.9), (xb_, Y_LA, Z_BU + 3.0), 1.4).fuse(stab((xb_, Y_LA, Z_BU - 5.9), (xb_, Y_LA, Z_BU - 5.0), 2.5))
    lg += [(f"Lufthebel_Lagerbock_{j_+1}", b_), (f"Lufthebel_Lagerbock_Schraube_{j_+1}", sc_)]
    boden = boden.cut(stab((xb_, Y_LA, Z_BU - 1), (xb_, Y_LA, Z_BU + 3.1), 1.45))
lg.append(("Lufthebel_Achse", ach))
# normale Gabelgelenke (Gelenkkopf kb_gelenkkopf.stp + Gabel + Bolzen, wie Mechanik); Bolzen in X oder Y, Richtungen in der Ebene senkrecht dazu
_gkl = cq.importers.importStep("kb_gelenkkopf.stp").solids().vals()
GKL_H, GKL_L, GKL_R = 3.19, 13.0, 2.0
GKL0 = _gkl[0].fuse(_gkl[1]).clean().cut(stab((0, 0.5, 0), (0, -GKL_L - 1, 0), 2.4)).cut(stab((-3, GKL_H, 0), (3, GKL_H, 0), 1.5))
GKL0 = GKL0.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180).translate(cq.Vector(0, GKL_H, 0))   # Auge im Ursprung, Schaft +Y
def _kan(u_, pin):                                          # Weltrichtung -> Winkel um X im Grundsystem (Bolzen X, Schaft +Y)
    a_, b_ = (u_[1], u_[2]) if pin == "X" else (-u_[0], u_[2])
    return math.degrees(math.atan2(b_, a_))
def _welt(s_, P, pin):
    if pin == "Y": s_ = s_.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 90)
    return s_.translate(cq.Vector(*P))
def _rx(s_, w_): return s_.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), w_)
def _einh(a_, b_):
    v_ = cq.Vector(*b_) - cq.Vector(*a_); v_ = v_.normalized(); return (v_.x, v_.y, v_.z)
def g_kopf(P, u_, pin="X"):                                 # Gelenkkopf, Schaft in Richtung u_
    return _welt(_rx(GKL0, _kan(u_, pin)), P, pin)
def g_gabel(P, v_, k_, l_, pin="X", bohrung=True):          # Gabel, Schaft in Richtung v_ (Laenge l_), Schlitz zum Kopf in Richtung k_
    wv_, wk_ = _kan(v_, pin), _kan(k_, pin)
    w_ = (0, math.cos(math.radians(wv_)), math.sin(math.radians(wv_)))
    g_ = stab((-GKL_R, 0, 0), (GKL_R, 0, 0), 2 * GKL_R).fuse(stab((0, 0, 0), tuple(l_ * c for c in w_), 2 * GKL_R))
    g_ = g_.cut(_rx(box(-1.05, 1.05, -2.4, 6.0, -2.4, 2.4), wk_))
    if bohrung: g_ = g_.cut(stab(tuple(GKL_H * c for c in w_), tuple((l_ + 1) * c for c in w_), 2.4))
    return _welt(g_.cut(stab((-3, 0, 0), (3, 0, 0), 1.5)), P, pin)
def g_bolzen(P, pin="X"):
    d_ = (2.4, 0, 0) if pin == "X" else (0, 2.4, 0)
    return stab(tuple(p - q for p, q in zip(P, d_)), tuple(p + q for p, q in zip(P, d_)), 1.4)
def _pl(P, u_, t_): return tuple(p + t_ * c for p, c in zip(P, u_))
AX_T = lambda x_: (x_, Y_LA, Z_LA)
L_GAB_T = 6.0                                               # kurze Gabeln an der Tastenstange (Abstand der Gelenke nur 15)
def taste_zug(j_, xs_, PA, PB, suf=""):
    """Schenkel (Achse -> Gelenk A), Gelenk A (Kopf am Schenkel, Gabel an der Stange), Stange, Gelenk B (Gabel an der Stange, Kopf in der Taste)"""
    uA = _einh(PA, AX_T(xs_)); vAB = _einh(PA, PB); vBA = _einh(PB, PA); dn_ = (0, 0, -1.0)
    sch_ = stab(AX_T(xs_), _pl(PA, uA, GKL_H - 0.5), 2.4).cut(ach)
    st_ = stab(_pl(PA, vAB, GKL_H), _pl(PB, vBA, GKL_H), 2.4)
    return [(f"Lufthebel_Schenkel_Taste_{j_}{suf}", sch_), (f"Lufttaste_Gelenkkopf_{j_}a{suf}", g_kopf(PA, uA)),
            (f"Lufttaste_Gelenkgabel_{j_}a{suf}", g_gabel(PA, vAB, uA, L_GAB_T)), (f"Lufttaste_Gelenkbolzen_{j_}a{suf}", g_bolzen(PA)),
            (f"Lufttaste_Stange_{j_}{suf}", st_), (f"Lufttaste_Gelenkgabel_{j_}b{suf}", g_gabel(PB, vBA, dn_, L_GAB_T)),
            (f"Lufttaste_Gelenkbolzen_{j_}b{suf}", g_bolzen(PB)), (f"Lufttaste_Gelenkkopf_{j_}b{suf}", g_kopf(PB, dn_))]
Z_TB = T_DECKEL_BASS + 3.0                                  # Gelenke auf der Taste
for j_, xs_ in enumerate(X_S1):
    lg += taste_zug(j_ + 1, xs_, (xs_, Y_T, Z_G1), (xs_, Y_T, Z_TB))
# Schubstange Taste -> Klappe: unten Gelenk auf der Taste, oben Gelenk an der Klappenunterseite, beide Bolzen in Y
# (Taste geht gerade hoch, Klappe dreht um die Laengsachse -> Stange schwenkt in der XZ-Ebene)
P_KU, P_KO = (X_S3, Y_T, Z_TB), (X_S3, Y_T, Z_LB1 - 3.5)
def klappen_zug(PU, PO, gab_o, suf=""):
    vUO, vOU = _einh(PU, PO), _einh(PO, PU)
    return [("Luftklappe_Gelenkkopf_unten" + suf, g_kopf(PU, (0, 0, -1.0), "Y")), ("Luftklappe_Gelenkgabel_unten" + suf, g_gabel(PU, vUO, (0, 0, -1.0), 8.0, "Y")),
            ("Luftklappe_Gelenkbolzen_unten" + suf, g_bolzen(PU, "Y")), ("Luftklappe_Schubstange" + suf, stab(_pl(PU, vUO, GKL_H), _pl(PO, vOU, GKL_H - 0.5), 2.4)),
            ("Luftklappe_Gelenkkopf_oben" + suf, g_kopf(PO, vOU, "Y")), ("Luftklappe_Gelenkgabel_oben" + suf, gab_o), ("Luftklappe_Gelenkbolzen_oben" + suf, g_bolzen(PO, "Y"))]
gab_o = g_gabel(P_KO, (0, 0, 1.0), (0, 0, -1.0), 2.5, "Y", bohrung=False).fuse(box(X_S3 - 2.4, X_S3 + 2.4, Y_T - 2.4, Y_T + 2.4, Z_LB1 - 1.0, Z_LB1))
lg += klappen_zug(P_KU, P_KO, gab_o)
for n_, s_ in list(lg):                                     # Gelenkkoepfe unten sitzen eingeschraubt in der Taste
    if n_.startswith(("Lufttaste_Gelenkkopf_1b", "Lufttaste_Gelenkkopf_2b", "Luftklappe_Gelenkkopf_unten")): taste = taste.cut(s_)
lg.append(("Lufttaste_Filz", lt_filz)); lg.append(("Lufttaste", taste))
teile += [("Luft_Seitenwand", sw.cut(bt)), ("Luft_Bodenteil", bt.cut(cq.Compound.makeCompound([s_ for n_, s_ in lg if n_.startswith("Luftklappe_Lagerbock_Schraube")])))]
# gedrueckte Lage (nur Pruefung): Achse um -ANG0 gedreht, Taste 10 hoeher, Klappe um die Laengsachse geoeffnet
TH_L = -ANG0
_ra = lambda s_: s_.rotate(cq.Vector(0, Y_LA, Z_LA), cq.Vector(1, Y_LA, Z_LA), TH_L)
L_KZ = math.dist(P_KU, P_KO)
P_KUg = (P_KU[0], P_KU[1], P_KU[2] + HUB_LT)
def _rkw(w_, sg_):
    return lambda s_: s_.rotate(cq.Vector(X_H, 0, Z_H), cq.Vector(X_H, sg_, Z_H), w_)
_sg = 1.0 if _rkw(10.0, 1.0)(cq.Vertex.makeVertex(*P_KO)).Center().z > P_KO[2] else -1.0
lo_, hi_ = 0.0, 45.0                                         # Oeffnungswinkel so, dass die Stangenlaenge gleich bleibt
for _ in range(40):
    mi_ = (lo_ + hi_) / 2
    p_ = _rkw(mi_, _sg)(cq.Vertex.makeVertex(*P_KO)).Center()
    if math.dist(P_KUg, (p_.x, p_.y, p_.z)) < L_KZ: lo_ = mi_
    else: hi_ = mi_
PHI_LK = (lo_ + hi_) / 2
_rk = _rkw(PHI_LK, _sg)
_po = _rk(cq.Vertex.makeVertex(*P_KO)).Center(); P_KOg = (_po.x, _po.y, _po.z)
dz3 = P_KOg[2] - P_KO[2]
lg_g = [("Lufttaste_gedrueckt", taste.translate(cq.Vector(0, 0, HUB_LT)))]
for j_, xs_ in enumerate(X_S1):
    pA = _ra(cq.Vertex.makeVertex(xs_, Y_T, Z_G1)).Center()
    lg_g += taste_zug(j_ + 1, xs_, (pA.x, pA.y, pA.z), (xs_, Y_T, Z_TB + HUB_LT), "_gedrueckt")
lg_g += klappen_zug(P_KUg, P_KOg, _rk(gab_o), "_gedrueckt")
lg_g += [(n_ + "_gedrueckt", _rk(s_)) for n_, s_ in lg if n_ in ("Luftklappe", "Luftklappe_Filz", "Luftklappe_Leder") or n_.startswith("Luftklappe_Lasche")]
deckel_r = deckel_r.cut(langloch_m(LT_L, LT_B, -1, T_DECKEL_BASS + 2))
teile += lg
print(f"Luft: Bodenteil X {XL0 - T_NUT_L:.1f}-{XS0 + T_NUT_L:.1f}, Ausschnitt Hauptboden X {XA0:.1f}-{XA1:.1f} Y {YA0:.1f}-{YA1:.1f} | Tastenachse Y {Y_LA:.1f} Z {Z_LA}"
      f" X {XB_[0]:.1f}-{XB_[1]:.1f}, Schenkel Taste X {[round(float(v), 1) for v in X_S1]} (L {L_S1:.1f}), Schubstange Klappe X {X_S3:.1f} L {L_KZ:.1f} | gedrueckt: Klappe {PHI_LK:.1f} Grad, Hub {dz3:.1f}")

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
teile = [(n_, s_.intersect(box(0, L, 0, B, -100, T_DECKEL_BASS + 3.0)) if n_ in ("Fuss_Schraube_3", "Fuss_Schraube_4") else s_) for n_, s_ in teile]   # vorne nur in den Winkel
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
        pg_ = _Poly([(m_r(a + 382.0), b + 4.24) for a, b in pp_]).buffer(R_GZ).simplify(0.05); _jf += 1
        _huellen.append((f"Gaze_FLoch_{_jf}", cq.Workplane("XZ").polyline(list(pg_.exterior.coords)[:-1]).close().extrude(-T_GZ).translate((0, D, 0)).val(), (pg_, "XZ", D)))
_jh = 0
for (hx_, hz_) in horn_pos:                                            # Gaze innen hinter den Hoernern (wie Markus)
    r = 9.0; _jh += 1
    _huellen.append((f"Gaze_Horn_{_jh}", cq.Solid.makeCylinder(r + R_GZ, T_GZ, cq.Vector(hx_, FI - T_GZ, hz_), cq.Vector(0, 1, 0)), (_Pt(hx_, hz_).buffer(r + R_GZ, 32), "XZ", FI - T_GZ)))
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
# ---------- Ausgabe mit Farben wie Markus/Standard ----------
teile = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + klappen
FARBEN = [("Lufttaste_Gelenkkopf", (0.60, 0.20, 0.70)), ("Luftklappe_Gelenkkopf", (0.60, 0.20, 0.70)), ("Lufttaste_Gelenkgabel", (0.95, 0.45, 0.10)), ("Luftklappe_Gelenkgabel", (0.95, 0.45, 0.10)), ("Luftklappe_Gelenkbolzen", (0.95, 0.85, 0.20)), ("Lufttaste_Stange", (0.10, 0.60, 0.60)), ("Luftklappe_Schubstange", (0.10, 0.60, 0.60)), ("Gaze", (0.78, 0.80, 0.82)), ("Fuss_Gummi", (0.08, 0.08, 0.08)), ("Fuss_Schraube", (0.6, 0.6, 0.62)), ("Fuss", (0.28, 0.20, 0.14)),
          ("Frontblende", (0.70, 0.52, 0.30)), ("Front_Aufdoppelung", (0.78, 0.60, 0.38)), ("Front", (0.78, 0.60, 0.38)), ("Horn", (0.05, 0.05, 0.05)), ("Knopfzylinder", (0.95, 0.95, 0.92)), ("Knopf", (0.15, 0.15, 0.15)), ("Filzeinlage", (0.55, 0.55, 0.20)), ("Rueck", (0.74, 0.56, 0.36)), ("Seiten", (0.70, 0.52, 0.33)), ("Aufdoppelung", (0.78, 0.60, 0.38)),
          ("Bassbodenleder", (0.12, 0.08, 0.05)), ("Bassboden", (0.90, 0.82, 0.62)), ("Bassdeckel_Leiste", (0.30, 0.22, 0.15)), ("Bassdeckel_Auflage", (0.60, 0.45, 0.28)),
          ("Bassdeckel", (0.30, 0.22, 0.15)), ("Deckelauflage", (0.60, 0.45, 0.28)), ("Eckblock", (0.55, 0.40, 0.25)), ("Eckwinkel_Schraube", (0.40, 0.40, 0.42)), ("Eckwinkel", (0.75, 0.76, 0.80)), ("Abdeckplatte", (0.55, 0.40, 0.25)), ("Abschlusswand", (0.55, 0.40, 0.25)), ("Innenleiste", (0.55, 0.40, 0.25)),
          ("Klappe_", None), ("HB_Folie", (0.85, 0.30, 0.30)), ("HB_Trennfeder", (0.85, 0.30, 0.30)), ("HB_Kammern", (0.65, 0.75, 0.88)), ("Lufttaste_Stoessel", (0.10, 0.60, 0.60)), ("Luftklappe_Feder_Schraube", (0.40, 0.40, 0.42)), ("Luft_Seitenwand", (0.65, 0.75, 0.88)), ("Luft_Bodenteil", (0.90, 0.82, 0.62)), ("Lufthebel_Lagerbock_Schraube", (0.40, 0.40, 0.42)), ("Lufthebel_Lagerbock", (0.25, 0.75, 0.35)), ("Lufthebel_Achse", (0.55, 0.55, 0.60)), ("Lufthebel_Schenkel", (0.20, 0.40, 0.85)), ("Luftklappe_Lagerbock_Schraube", (0.40, 0.40, 0.42)), ("Luftklappe_Lagerbock", (0.25, 0.75, 0.35)), ("Luftklappe_Achse", (0.55, 0.55, 0.60)), ("Luftklappe_Gelenkbolzen", (0.95, 0.85, 0.20)), ("Luftklappe_Gelenk", (0.95, 0.45, 0.10)), ("Luftklappe_Schubstange", (0.10, 0.60, 0.60)), ("Luftklappe_Lasche", (0.20, 0.20, 0.22)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)), ("Lufttaste_Gelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Stange", (0.10, 0.60, 0.60)), ("Luftklappe_Stange", (0.10, 0.60, 0.60)), ("Luftklappe_Feder", (0.35, 0.45, 0.70)), ("HB_Stimmplatte", (0.80, 0.82, 0.86)), ("HB_Trennwand", (0.50, 0.62, 0.80)), ("HB_", (0.65, 0.75, 0.88)),
          ("Begleiter_Stimmplatte", (0.80, 0.82, 0.86)), ("Begleiter_Sole", (0.55, 0.68, 0.45)), ("Begleiter_Schraube", (0.40, 0.40, 0.42)), ("Begleiter_Haken", (0.40, 0.40, 0.42)), ("Begleiter", (0.62, 0.72, 0.55)), ("Haken", (0.40, 0.40, 0.42)), ("Formteil_Schraube", (0.40, 0.40, 0.42)), ("Formteil", (0.60, 0.45, 0.28)), ("Dichtung", (0.12, 0.08, 0.05)), ("Lufttaste_Filz", (0.85, 0.15, 0.15)), ("Lufttaste_Gabelgelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)),
          ("Lufttaste_Schubstange", (0.10, 0.60, 0.60)), ("Lufttaste", (0.35, 0.25, 0.15)), ("Lufthebel_Feder", (0.75, 0.10, 0.75)), ("Lufthebel_Stuetze_Schraube", (0.40, 0.40, 0.42)),
          ("Lufthebel_Stuetze", (0.25, 0.75, 0.35)), ("Lufthebel_Bolzen", (0.95, 0.85, 0.20)), ("Lufthebel", (0.70, 0.72, 0.78)), ("Luftklappe_Leder", (0.45, 0.30, 0.20)),
          ("Luftklappe_Filz", (0.85, 0.85, 0.80)), ("Luftklappe", (0.20, 0.20, 0.22))]
def farbe(n):
    if n.startswith(("Klappe_", "Begleiterklappe_")): return (0.45, 0.30, 0.20) if "Leder" in n else (0.85, 0.85, 0.80) if "Filz" in n else (0.20, 0.20, 0.22)
    for p_, c_ in FARBEN:
        if n.startswith(p_) and c_: return c_
    print("FARBE ohne Zuordnung:", n); return (0.6, 0.6, 0.6)
exec(open("mech_v340.py").read())             # Mechanik (vor der Spiegelung)
a = cq.Assembly(name="Bassteil_340_schmal_schraeg")
SPIEGELN = True                                # wie Markus: links/rechts getauscht (X = L/2)
for n, s_ in teile:
    if SPIEGELN: s_ = s_.mirror("YZ", basePointVector=(L / 2, 0, 0))
    a.add(s_, name=n, color=cq.Color(*farbe(n)))
a.save(OUT)
print("gespeichert", OUT, len(teile), "Teile")

