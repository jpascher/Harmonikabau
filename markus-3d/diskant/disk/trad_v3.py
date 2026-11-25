# Traditioneller Griffstock (Johann 20.11.2026, Zeichnung trad_schnitt.dxf, Layer "kammer 1 reihen"):
# zwei Teile mit je einer durchgehenden Achse Ø2; in jeder Kammer ein Holzhebel, der auf der Achse dreht.
# Teil 2 (vorne unten, Reihen 3/4) und Teil 1 (Reihen 1/2); jede Kammer hat einen 1,5-mm-Kammergrund (Kaemme gehen nicht
# durch), Kammern zur Knopfseite und nach unten offen. Kammerlaenge je Reihe: bis 1,5 mm unter die Oberkante der eigenen
# Stufe. Kammer-X und Kammerbreite wie im bestehenden Modell (Knopf-X inkl. Knopfversatz, Holz 4,8 in 5,4).
# Die Alu-Arme zu den Klappen bleiben (Fuehrung und Arm wie bisher), sie sitzen im Holzhebel und haengen an der Achse.
# Wird am Ende von hebel.py ausgefuehrt und ersetzt dort Kunststoffteil, Dreieckleiste, Federn, Holzteile, Buchsen,
# Einzelachsen und die beiden Griffstockteile.
import numpy as _np
from shapely.geometry import Polygon as _TP, Point as _TPt, LineString as _TL, box as _tbox
from shapely import affinity as _taff
from shapely.ops import unary_union as _tuu
import schnitt_trafo as _TST

_TDX = 58.17                                              # trad_schnitt.dxf liegt um +58,17 in x gegen Schnitt_Hebel_4_Layer.dxf
def _M(p):
    u_, Y_ = _TST.modell((p[0] - _TDX, p[1])); return (float(u_), float(Y_))
_src = _np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
_AFF = _np.linalg.solve(_np.c_[_src, _np.ones(3)], _np.array([_M(p_) for p_ in _src]))
def _Minv(q):
    return tuple(float(v_) for v_ in _np.linalg.solve(_AFF[:2, :].T, _np.array(q) - _AFF[2]))
def _Mpg(pg):                                             # Zeichnung -> (u, Y)
    pg = max(pg.geoms, key=lambda g_: g_.area) if hasattr(pg, "geoms") else pg
    return _TP([_M(p_) for p_ in list(pg.exterior.coords)[:-1]], [[_M(p_) for p_ in list(h_.coords)[:-1]] for h_ in pg.interiors]).buffer(0)
def _pp(pg):                                              # saubere Punktliste (u, Y) fuer prisma_uy
    pg = pg.simplify(0.02).buffer(0)
    pg = max(pg.geoms, key=lambda g_: g_.area) if hasattr(pg, "geoms") else pg
    out = []
    for a_, b_ in list(pg.exterior.coords)[:-1]:
        if not out or math.hypot(a_ - out[-1][0], b_ - out[-1][1]) > 0.02: out.append((float(a_), float(b_)))
    while len(out) > 3 and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) <= 0.02: out.pop()
    return out
def _tteile(pg):
    return [g_ for g_ in (pg.geoms if hasattr(pg, "geoms") else [pg]) if g_.geom_type == "Polygon" and g_.area > 0.05]

_tm = ezdxf.readfile("trad_schnitt.dxf").modelspace()
_tpl = [[(float(x_), float(y_)) for x_, y_, *_r in e.get_points()] for e in _tm.query('LWPOLYLINE[layer=="kammer 1 reihen"]')]
_tpoly = sorted([p_ for p_ in _tpl if len(p_) >= 6], key=lambda p_: -max(q_[1] for q_ in p_))
P12, P34 = _TP(_tpoly[0]).buffer(0), _TP(_tpoly[1]).buffer(0)          # Reihen 1/2 (hoch, hinten) / Reihen 3/4 (vorne unten)
_grund34 = [p_ for p_ in _tpl if len(p_) < 6][0]                         # Kammergrund Teil 2 (1,5 mm innen)
X_GR34 = max(q_[0] for q_ in _grund34); Y_GR34 = max(q_[1] for q_ in _grund34)
_ach = [(float(e.dxf.center.x), float(e.dxf.center.y)) for e in _tm.query('CIRCLE[layer=="2-BassgestängeQuer"]')]
A12 = [a_ for a_ in _ach if P12.contains(_TPt(a_))][0]; A34 = [a_ for a_ in _ach if P34.contains(_TPt(a_))][0]
WAND_T = 1.5
X_GR12 = P12.bounds[2] - WAND_T; Y_OB12 = P12.bounds[3] - WAND_T
# Stufen (Oberkanten) aus der Zeichnung: Reihe 4 bis 83,92, Reihe 3 bis Oberkante Teil 2, Reihe 2 bis zur Stufe 1, Reihe 1 bis Teil-1-Kammer
_ys34 = sorted(set(round(q_[1], 2) for q_ in _tpoly[1]))
_ys12 = sorted(set(round(q_[1], 2) for q_ in _tpoly[0]))
Y_ST4 = max(y_ for y_ in _ys34 if y_ < P34.bounds[3] - 1.0)               # Stufe 4/3 (83,92)
Y_ST2 = max(y_ for y_ in _ys12 if y_ < P12.bounds[3] - 1.0)               # Stufe 2/1 (116,45)
_KB = {}
_kb34 = P34.intersection(_tbox(-500, -500, X_GR34, Y_GR34))
_kb12 = P12.intersection(_tbox(-500, -500, X_GR12, Y_OB12))
# 20.11. (5): Kamm mit 1,5-mm-Boden hinten; Reihe 1 lang, Reihe 2 kuerzer (Teil 1), Reihe 3 lang, Reihe 4 kurz (Teil 2);
# jede Nut hat oben 1,5 mm Abschluss, zur Knopfseite und nach unten offen
_kb34 = P34.intersection(_tbox(-500, -500, X_GR34, Y_GR34))
_kb12 = P12.intersection(_tbox(-500, -500, X_GR12, Y_OB12))
_KB[1] = _kb12
_KB[2] = _kb12.intersection(_tbox(-500, -500, 500, Y_ST2 - WAND_T))
_KB[3] = _kb34
_KB[4] = _kb34.intersection(_tbox(-500, -500, 500, Y_ST4 - WAND_T))
for r_ in (1, 2, 3, 4):
    k_ = _KB[r_]
    _KB[r_] = _tuu([k_, _taff.translate(k_.intersection(_tbox(-500, -500, 500, k_.bounds[1] + 3.0)), 0, -6.0)]).buffer(0)
Y_OBEN_T1 = P12.bounds[3]                                 # Teil 1 endet oben mit deiner Kontur (132,1), darueber nichts
_ACHSE = {1: A12, 2: A12, 3: A34, 4: A34}
print(f"Traditionell: Achse Reihen 1/2 {tuple(round(v_, 2) for v_ in _M(A12))}, Reihen 3/4 {tuple(round(v_, 2) for v_ in _M(A34))} (u, Y); "
      f"Kammergrund Teil 1 x {X_GR12:.2f}, Teil 2 x {X_GR34:.2f}; Kammer-Oberkanten",
      {r_: round(_KB[r_].bounds[3], 2) for r_ in _KB})

# --- alte Mechanikteile entfernen ---
_weg_n = {"Griffstock_Koerper", "Griffstock_Teil2", "Griffstock_Abdeckung", "Kunststoffteil_Lagerblock", "Dreieckleiste"}
_weg_s = ("/Holzteil", "/Hebel_Alu", "/Messingbuchse", "/Achse", "/Feder")
_weg_p = ("Schraube_Kunststoff_", "Schraube_Boden_Dreieck_")
_vorher = len(teile)
teile = [t_ for t_ in teile if not (t_[0] in _weg_n or t_[0].endswith(_weg_s) or t_[0].startswith(_weg_p))]
print(f"Traditionell: {_vorher - len(teile)} Teile der Kunststoff-Lagerung entfernt")

# --- Holzhebel und Alu je Taste ---
_THUB = 3.6
AUS_ABST = 3.5                                            # Rechteckausschnitt Teil 2 bis 3,5 mm unter die Achse 3/4
X_KR0, X_KR1 = P34.bounds[2] + 0.3, P34.bounds[2] - 5.2  # Kroepfzone Alu 1/2 (Zeichnung x), im Ausschnitt von Teil 2
TRAD_HEBEL = []
_TRAD_2D = {}
_sweep_holz, _sweep_alu = {}, {}
def _rot(pg, wink, um):
    return _taff.rotate(pg, wink, origin=um, use_radians=False)
for n_, (r, x, k, kr) in enumerate(zuordnung):
    G = f"Hebel_R{r}_{n_+1:02d}"
    xk = x + KV[r]                                        # Holzhebel, Kammer und Knopf in Knopf-X (mit Knopfversatz)
    ax_ = _ACHSE[r]
    ku, kY = knopf_u[r - 1]
    if JP: ku, kY = ku + JP_KNOPF_D[r][0], kY + JP_KNOPF_D[r][1]
    bx, by = _Minv((ku, kY))                              # Knopfauflage in der Zeichnung (links der Stufe)
    face = bx + 3.5
    dist = math.hypot(bx - ax_[0], by - ax_[1])
    th = math.degrees(_THUB / dist)
    C = _KB[r].buffer(-0.4)
    # gedrueckt: oberes Ende nach +x -> Drehung im Uhrzeigersinn (-th); gedrueckt muss in der Kammer bleiben
    L0 = _tuu([_tbox(ax_[0] - 5.0, ax_[1] - 6.0, ax_[0] + 5.0, by + 4.0), _tbox(bx, by - 4.0, ax_[0] + 5.0, by + 4.0),
               _TPt(ax_).buffer(5.0, 48)])
    Lh = L0.intersection(C).intersection(_rot(C, th, ax_))
    Lh = _tuu([Lh, _tbox(bx, by - 4.0, face + 0.6, by + 4.0)]).buffer(0)
    Lh = max(Lh.geoms, key=lambda g_: g_.area) if hasattr(Lh, "geoms") else Lh
    if not Lh.contains(_TPt(ax_).buffer(2.5)): print(f"   ACHTUNG {G}: Holzhebel umschliesst die Achse nicht")
    Lh = Lh.difference(_TPt(ax_).buffer(1.0, 32))
    # Alu: Verbindungsstueck von der Achse (im Holz) bis zum Beginn des langen Arms bei ya, dann Arm wie bisher
    if kr != r:
        ya = 36.0; pfad = [(_M(ax_)[1], x), (ya, x), (46.0, k[0]), (k[1], k[0])]
    else:
        f = FUEHR[r]; v = k[0]; lane = v + f["lane"]
        if r == 1 and f["lane"] != 0.0 and not any(abs(k3[0] - (v + 11.5)) < 4.0 for k3 in kl_r[3]): lane = v + 5.0
        ya, yb = f["y_lane"]
        pfad = [(_M(ax_)[1], x), (ya, x), (yb, lane)]
        if f["y_fin"]:
            yf0, yf1 = f["y_fin"]
            if yf0 - yb > 0.5: pfad.append((yf0, lane))
            pfad.append((yf1, v))
        pfad.append((k[1], v))
    bd = arm_band(r, ya - 0.5) or arm_band(r, ya + 1)
    if r == 1 and ya < 40: bd = arm_band(1, 36.5)
    if os.environ.get("DEBUG_TRAD"):
        _DBG = globals().setdefault("_DBG_TRAD", {})
        if r not in _DBG:
            from shapely.geometry import Polygon as _dP
            _DBG[r] = dict(ya=ya, bd=bd, knopf=(ku, kY), bxby=(bx, by), Lh=Lh, th=th, hebel=_dP(HEBEL[r]).buffer(0), yk=k[1], U_KL=U_KL)
        if len(_DBG) == 4:
            import pickle as _pkT
            _pkT.dump({"rows": _DBG, "P12": P12, "P34": P34, "KB": _KB, "A12": A12, "A34": A34, "HUELLE_G": HUELLE_G, "AFF": _AFF}, open(os.environ["DEBUG_TRAD"], "wb"))
            raise SystemExit("DEBUG_TRAD geschrieben")
    _p0, _p1 = _Minv((bd[0], ya + 0.5)), _Minv((bd[1], ya + 0.5))
    if ax_ == A12:      # Reihen 1/2: unter dem Rechteckausschnitt von Teil 2 durch (Oberkante <= Ausschnitt - 1)
        _yo = A34[1] - AUS_ABST - 1.0
        _st = _TL([(X_KR0 + 1.0, _yo - 6.5), (X_KR0 + 1.0, _yo)]).buffer(0.01)
        Ca = _tuu([_tuu([_TPt(ax_).buffer(3.5, 48), _st]).convex_hull,
                   _tuu([_st, _TL([_p0, _p1]).buffer(0.01)]).convex_hull]).buffer(0).difference(_TPt(ax_).buffer(1.0, 32))
    else:
        Ca = _tuu([_TPt(ax_).buffer(3.5, 48), _TL([_p0, _p1]).buffer(0.01)]).convex_hull.difference(_TPt(ax_).buffer(1.0, 32))
    teile_arm = arm_durchgehend(r, ya, k[1], pfad)
    alu = prisma_uy(_pp(_Mpg(Ca)), xk, T_ALU)
    if abs(xk - x) > 1e-6:                                  # seitliche Kroepfung im Ausschnitt: xk -> x
        _Ya, _Yb = _M((X_KR0, 35.0))[1], _M((X_KR1, 35.0))[1]
        W_ = 1e4
        _p1_ = alu.intersect(box(-W_, W_, -W_, _Ya, -W_, W_))
        _p3_ = alu.intersect(box(-W_, W_, _Yb, W_, -W_, W_)).translate(cq.Vector(x - xk, 0, 0))
        _p2_ = alu.intersect(box(-W_, W_, _Ya, _Yb, -W_, W_))
        _k_ = (x - xk) / (_Yb - _Ya)
        _p2_ = _p2_.transformGeometry(cq.Matrix([[1.0, _k_, 0.0, -_k_ * _Ya], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]))
        alu = _p1_
        for t_ in (_p2_, _p3_):
            if t_.Volume() > 1e-6: alu = alu.fuse(t_)
    for s_ in teile_arm: alu = alu.fuse(s_)
    _uax, _Yax = _M(ax_)
    alu = alu.cut(cq.Solid.makeCylinder(1.0, T_ALU + 2, cq.Vector(xk - T_ALU / 2 - 1, _Yax, -_uax), cq.Vector(1, 0, 0)))
    holz = prisma_uy(_pp(_Mpg(Lh)), xk, T_HOLZ)
    holz = holz.cut(cq.Solid.makeCylinder(1.0, T_HOLZ + 2, cq.Vector(xk - T_HOLZ / 2 - 1, *(lambda q_: (q_[1], -q_[0]))(_M(ax_))), cq.Vector(1, 0, 0)))
    holz = holz.cut(prisma_uy(_pp(_Mpg(Ca.buffer(0.05))), xk, T_ALU + 0.1))     # Schlitz fuer das Alu
    teile.append((f"{G}/Holzhebel", holz, f"Holz_R{r}"))
    teile.append((f"{G}/Hebel_Alu", alu, f"Alu_R{r}"))
    DREHPUNKTE[G] = (xk, _M(ax_)[1], -_M(ax_)[0])
    # Schwenkbereich fuer die Kammern (Ruhe bis gedrueckt, 0,3 Luft)
    _sh = _tuu([_rot(Lh, -th * i_ / 6.0, ax_) for i_ in range(7)]).buffer(0.3)
    _Cs = Ca.intersection(_tbox(X_KR0, -500, 500, 500)) if abs(xk - x) > 1e-6 else Ca
    _sa = _tuu([_rot(_Cs, -th * i_ / 6.0, ax_) for i_ in range(7)]).buffer(0.3)
    _sweep_holz.setdefault(round(xk, 2), []).append(_sh); _sweep_alu.setdefault(round(xk, 2), []).append((r, _sa))
    if abs(xk - x) > 1e-6:                                  # Alu nach der Kroepfung (in der alten Hebel-X) ebenfalls frei schneiden
        _Cs2 = Ca.intersection(_tbox(-500, -500, X_KR1, 500))
        if not _Cs2.is_empty:
            _sweep_alu.setdefault(round(x, 2), []).append((r, _tuu([_rot(_Cs2, -th * i_ / 6.0, ax_) for i_ in range(7)]).buffer(0.3)))
    TRAD_HEBEL.append((G, r, xk, round(th, 2), round(dist, 1)))
    if kr == r and r not in _TRAD_2D: _TRAD_2D[r] = dict(L=Lh, Ca=Ca, th=th, ax=ax_, knopf=(bx, by), G=G)
print(f"Traditionell: {len(TRAD_HEBEL)} Holzhebel, Hub-Winkel {min(t_[3] for t_ in TRAD_HEBEL)}-{max(t_[3] for t_ in TRAD_HEBEL)} Grad")

# --- Griffstockteile ---
_T2uy = _Mpg(P34)
# 20.11. (3): kein hinterer Deckel - Teil 1 = genau Johanns Kontur (Kamm mit 1,5-mm-Boden), Teil 2 ebenso
_DAUMEN = _gr(HUELLE_G.intersection(_Mpg(_tbox(-500, P12.bounds[3] - 0.01, 500, 1000))))   # durchgehende Daumenleiste
_T1uy = _gr(_tuu([_Mpg(P12), _DAUMEN]).buffer(0.01).buffer(-0.01).difference(_T2uy))
_ENV = _gr(_tuu([_T1uy, _T2uy]).buffer(0.01).buffer(-0.01))
HUELLE_G = _ENV
teile = [t_ for t_ in teile if t_[0] != "Griffstock_Abdeckung"]
_LGS = GS_X1 - GS_X0
t1 = prisma_uy(_pp(_T1uy), (GS_X0 + GS_X1) / 2, _LGS)
t2 = prisma_uy(_pp(_gr(_T2uy)), (GS_X0 + GS_X1) / 2, _LGS)
t1 = schneide(t1, seitenschnitt); t2 = schneide(t2, seitenschnitt)
_kam = []
for n_, (r, x, k, kr) in enumerate(zuordnung):
    xk = x + KV[r]                                        # Holzhebel, Kammer und Knopf in Knopf-X (mit Knopfversatz)
    for g_ in _tteile(_Mpg(_KB[r])): _kam.append(prisma_uy(_pp(g_), xk, T_HOLZ + 0.6))
_extra = []
for xk_, lst_ in _sweep_holz.items():
    for g_ in _tteile(_tuu(lst_)):
        _extra.append(prisma_uy(_pp(_Mpg(g_)), xk_, T_HOLZ + 0.6))
for xk_, lst_ in _sweep_alu.items():
    for r_, g_ in lst_:
        for h_ in _tteile(g_): _extra.append(prisma_uy(_pp(_Mpg(h_)), xk_, T_ALU + 0.6))
_v1, _v2 = t1.Volume(), t2.Volume()
t1 = schneide(t1, _kam); t2 = schneide(t2, _kam)
_v1k, _v2k = t1.Volume(), t2.Volume()
t1 = schneide(t1, _extra); t2 = schneide(t2, _extra)
print(f"Traditionell: Kammern nehmen Teil 1 {_v1 - _v1k:.0f}, Teil 2 {_v2 - _v2k:.0f} mm3; Schwenkbereich Hebel/Alu zusaetzlich "
      f"Teil 1 {_v1k - t1.Volume():.0f}, Teil 2 {_v2k - t2.Volume():.0f} mm3 (u. a. Durchgang der Alu-Arme 1/2 durch Teil 2)")
# Ausschnitt (Rechteck) in Teil 2 von unten (Bodenseite) bis kurz vor die Achse 3/4, ueber den Hebelbereich; Enden massiv.
# Darunter laufen die Alu-Arme der Reihen 1/2 durch.
AUS_ABST = 3.5
_xs_h = [x_ for r_, x_, k_, kr_ in zuordnung]
_aus = P34.intersection(_tbox(-500, -500, 500, A34[1] - AUS_ABST))
if not _aus.is_empty:
    t2 = t2.cut(prisma_uy(_pp(_Mpg(_aus.buffer(0.5).intersection(_tbox(-500, -500, 500, A34[1] - AUS_ABST)))),
                          (min(_xs_h) + max(_xs_h)) / 2, max(_xs_h) - min(_xs_h) + 16.0))
    print(f"Traditionell: Rechteckausschnitt Teil 2 von unten bis {AUS_ABST} mm unter die Achse 3/4, X {min(_xs_h) - 8:.0f}-{max(_xs_h) + 8:.0f}")
# durchgehende Achsen Ø2 und ihre Bohrungen
for nm_, a_ in (("Reihen12", A12), ("Reihen34", A34)):
    u_, Y_ = _M(a_)
    _cyl = cq.Solid.makeCylinder(1.0, _LGS + 20.0, cq.Vector(GS_X0 - 10.0, Y_, -u_), cq.Vector(1, 0, 0))
    t1 = t1.cut(_cyl); t2 = t2.cut(_cyl)
    _rod = schneide(cq.Solid.makeCylinder(1.0, _LGS - 2.0, cq.Vector(GS_X0 + 1.0, Y_, -u_), cq.Vector(1, 0, 0)), seitenschnitt)
    teile.append((f"Achse_durchgehend_{nm_}", _rod, "Stahl"))
# Griffstock-Schrauben von unten wie bisher
if GS_LOECHER:
    _gsl = [w_ for x_, y_, n_, w_ in GS_LOECHER]
    t1 = schneide(t1, _gsl); t2 = schneide(t2, _gsl)
teile.append(("Griffstock_Teil2", t2, "Griffstock"))
teile.append(("Griffstock_Koerper", t1, "Griffstock"))
print(f"Traditionell: Griffstock Teil 1 {t1.Volume():.0f} mm3, Teil 2 {t2.Volume():.0f} mm3, Achsen 2 x {_LGS - 2.0:.0f} mm")

# --- Ergebnis-Schnitt als DXF/PNG im Koordinatensystem von trad_schnitt.dxf (je Reihe der erste Hebel) ---
try:
    _dout = ezdxf.new(); _dm = _dout.modelspace()
    _farb = {"Teil1": 1, "Teil2": 3, "Holzhebel": 30, "Holzhebel_gedrueckt": 8, "Alu": 5, "Achse": 7, "Knopf": 2, "Kammer": 4}
    def _dpl(pg, lay):
        if lay not in _dout.layers: _dout.layers.add(lay, color=_farb.get(lay.split("_R")[0], 7))
        for g_ in _tteile(pg):
            _dm.add_lwpolyline(list(g_.exterior.coords), close=True, dxfattribs={"layer": lay})
            for h_ in g_.interiors: _dm.add_lwpolyline(list(h_.coords), close=True, dxfattribs={"layer": lay})
    def _uy2d(pg):
        pg = max(pg.geoms, key=lambda g_: g_.area) if hasattr(pg, "geoms") else pg
        return _TP([_Minv(p_) for p_ in list(pg.exterior.coords)[:-1]]).buffer(0)
    _dpl(_uy2d(_T1uy), "Teil1"); _dpl(P34.difference(_aus), "Teil2"); 
    for r_, D_ in sorted(_TRAD_2D.items()):
        _dpl(_KB[r_].intersection(_tuu([P12, P34])), f"Kammer_R{r_}")
        _dpl(D_["L"], f"Holzhebel_R{r_}"); _dpl(_rot(D_["L"], -D_["th"], D_["ax"]), f"Holzhebel_gedrueckt_R{r_}")
        _dpl(D_["Ca"].intersection(_tbox(-200, -100, 100, 300)), f"Alu_R{r_}")
        bx_, by_ = D_["knopf"]; _dpl(_tbox(bx_ - 3.0, by_ - 7.5, bx_, by_ + 7.5), f"Knopf_R{r_}")
    for a_ in (A12, A34):
        if "Achse" not in _dout.layers: _dout.layers.add("Achse", color=7)
        _dm.add_circle(a_, 1.0, dxfattribs={"layer": "Achse"})
    _dn = OUT.replace(".step", "_trad_schnitt.dxf"); _dout.saveas(_dn)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as _plt
    _fg, _ax = _plt.subplots(1, 4, figsize=(20, 9))
    for i_, (r_, D_) in enumerate(sorted(_TRAD_2D.items())):
        a_ = _ax[i_]
        for pg_, c_, lw_ in ((_uy2d(_T1uy), "#b5651d", 1.0), (P34.difference(_aus), "#2e8b57", 1.0), (_KB[r_].intersection(_tuu([P12, P34])), "#4682b4", 0.6),
                             (D_["L"], "#d2691e", 1.6), (_rot(D_["L"], -D_["th"], D_["ax"]), "#999999", 0.8), (D_["Ca"], "#1f3fbf", 1.2)):
            for g_ in _tteile(pg_): a_.plot(*g_.exterior.xy, color=c_, lw=lw_)
        bx_, by_ = D_["knopf"]; a_.add_patch(_plt.Rectangle((bx_ - 3.0, by_ - 7.5), 3.0, 15.0, color="#e8d9a0"))
        for q_ in (A12, A34): a_.add_patch(_plt.Circle(q_, 1.0, color="k"))
        a_.set_aspect("equal"); a_.set_title(f"Reihe {r_} ({D_['G']}), Hub {D_['th']:.1f} Grad"); a_.set_xlim(-95, -35); a_.set_ylim(0, 150)
    _fg.tight_layout(); _fg.savefig(OUT.replace(".step", "_trad_schnitt.png"), dpi=90)
    print("Traditionell: Schnitt je Reihe ->", _dn)
except Exception as _ex:
    print("Traditionell: Schnitt-DXF nicht geschrieben:", _ex)
