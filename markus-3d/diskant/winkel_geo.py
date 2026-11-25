# 22.11.2026: Griffstock "winkel" - wie "kolben" (verdeckt, Zylinder unter den Knoepfen, Hebel nur Alu), aber statt
# Kunststoffteil/Dreieckleiste ein Alu-Winkel (Fuehrung 46753, Alu 2 mm) innen an der Rueckwand (Deckel):
#   Platte 79 Grad wie der Griffstock, unten 5 mm Umbiegung waagrecht auf dem Boden (verschraubt), oben 2 mm Einbiegung
#   zur Hebelseite (Randverzahnung fuer die 0,6-mm-Federdraehte), 50 Schlitze 2,4 breit / Teilung 4,8,
#   5 Laschen, Achse Ø3 durch die unteren Laschenbohrungen - alle Hebel nebeneinander auf dieser einen Achse.
# Hebel (Alu): Knopfarm wie Kolben-Variante oberhalb des Zylinders, darunter neu entlang der Rueckwand durch den Schlitz,
# Klappenarm als Steg 4,1..10,9 ueber der Achse (wie Lagerung Knick/V3), Ms-Lasche 0,5 mit Bohrung Ø3 um die Achse.
# Seitlich: im Winkel liegt jeder Hebel in seinem Schlitz, Z-Kroepfung zum Knopf (u 34..44) und zur Klappe (Y 36..ya-1).
import math as _wm
from shapely.geometry import Polygon as _WP, Point as _WPt, LineString as _WL, box as _Wbox
from shapely.ops import unary_union as _Wuu
W_T, W_L, W_FL, W_EIN, W_EIN_W, W_SCHLITZ_S = 2.0, 31.0, 10.0, 2.0, 45.0, 28.0     # 16:36 Johann: Umbiegung unten 10 mm
W_SCHLITZ_B, W_TEIL = 2.4, 4.8
# 22.11. 17:17 (Johann): je ein Hebel aussen neben den aeusseren Laschen, dazwischen 4 Bereiche (bei 50: 1+4x12+1);
# auf der Achse liegen die Ms-Laschen (4,8 breit) dicht an dicht, an einer Winkel-Lasche (2 dick)
# 22.11. 22:34: Knopfanzahl variabel (44/46/48/50/52): N-2 Hebel moeglichst gleich auf 4 Bereiche, groessere Bereiche innen
W_SCHLITZE = (sum(len(R_["klappen"]) for R_ in INSTR["reihen"].values())          # ein Hebel je Klappe
              if globals().get("INSTR") is not None else 50)
_wq, _wr = divmod(W_SCHLITZE - 2, 4)
W_BEREICHE = [_wq + (1 if k_ in ((1, 2, 0, 3)[:_wr]) else 0) for k_ in range(4)]
# 17:19 (Johann): je Bereich (Ms-Laschen 4,8 + Winkel-Lasche 2) 2,4 mm Spiel, gleichmaessig auf die n+1 Zwischenraeume
W_SPIELE = [2.4 / (n_ + 1.0) for n_ in W_BEREICHE]
W_SPIEL = W_SPIELE[0]
_wlz = [0.0]                                             # Mitten der Winkel-Laschen
for _k in range(4): _wlz.append(_wlz[-1] + 2.0 + W_BEREICHE[_k] * W_TEIL + 2.4)
W_SCHLITZ_X = [_wlz[0] - 1.0 - W_SPIELE[0] - W_TEIL / 2]
for _k in range(4):
    for _i in range(W_BEREICHE[_k]): W_SCHLITZ_X.append(_wlz[_k] + 1.0 + W_SPIELE[_k] + W_TEIL / 2 + _i * (W_TEIL + W_SPIELE[_k]))
W_SCHLITZ_X.append(_wlz[4] + 1.0 + W_SPIELE[3] + W_TEIL / 2)
W_LASCHE_NACH = tuple([0] + [sum(W_BEREICHE[:k_ + 1]) for k_ in range(4)])   # Winkel-Lasche zwischen Hebel i und i+1
_w0 = W_SCHLITZ_X[0]; W_SCHLITZ_X = [v_ - _w0 for v_ in W_SCHLITZ_X]
W_LASCHEN = tuple(v_ - _w0 - 1.0 for v_ in _wlz)        # Laschen (2 dick) ab x
print(f"Winkel: {W_SCHLITZE} Hebel = 1 + {'+'.join(map(str, W_BEREICHE))} + 1, Spiel je Zwischenraum {[round(v_, 3) for v_ in W_SPIELE]}")
W_LASCHE_PTS = [(4.0, 0.0), (21.683, 26.236), (21.832, 26.442), (21.994, 26.639), (22.168, 26.825), (22.354, 27.0), (22.55, 27.162), (22.756, 27.312), (22.972, 27.448), (23.196, 27.57), (23.426, 27.677), (23.664, 27.77), (23.907, 27.848), (24.154, 27.909), (24.404, 27.955), (24.657, 27.985), (24.912, 27.999), (25.167, 27.997), (25.421, 27.978), (25.673, 27.943), (25.923, 27.892), (26.169, 27.825), (26.41, 27.743), (26.645, 27.646), (26.874, 27.534), (27.095, 27.407), (27.308, 27.267), (27.511, 27.113), (27.705, 26.947), (27.887, 26.769), (28.057, 26.579), (28.215, 26.38), (28.36, 26.17), (28.492, 25.952), (28.609, 25.725), (28.711, 25.492), (28.799, 25.253), (28.871, 25.008), (28.927, 24.76), (28.968, 24.508), (28.992, 24.255), (29.0, 24.0), (29.0, 0.0)]
W_ACHSE_D = 3.0
_wro0 = min((q for q in GRIFFSTOCK if q[0] < 10), key=lambda q: q[1]); _wro1 = min(GRIFFSTOCK, key=lambda q: q[1])
_wdu, _wdy = _wro1[0] - _wro0[0], _wro1[1] - _wro0[1]; _wl = _wm.hypot(_wdu, _wdy)
W_ES = (_wdu / _wl, _wdy / _wl)                          # entlang der Rueckwand nach oben (u, Y)
W_EN = (-W_ES[1], W_ES[0])                               # nach vorne (+Y)
W_O = (0.0, _wro0[1] + (0.0 - _wro0[0]) * _wdy / _wdu + 3.0)   # Deckel-Innenseite am Boden
def W_P(s, n): return (W_O[0] + s * W_ES[0] + n * W_EN[0], W_O[1] + s * W_ES[1] + n * W_EN[1])
_wyf0 = _WL([W_P(-5, W_T), W_P(W_L, W_T)]).intersection(_WL([(0, -50), (0, 100)])).y
def _w_profil():
    pl = _WP([W_P(-5, 0), W_P(W_L, 0), W_P(W_L, W_T), W_P(-5, W_T)]).intersection(_Wbox(0, -50, 100, 100))
    fl = _WP([(0, W_O[1]), (0, _wyf0 + W_FL), (W_T, _wyf0 + W_FL), (W_T, W_O[1])])
    a = _wm.radians(W_EIN_W); d = (_wm.cos(a) * W_ES[0] + _wm.sin(a) * W_EN[0], _wm.cos(a) * W_ES[1] + _wm.sin(a) * W_EN[1])
    nd = (-d[1], d[0]) if (-d[1]) * W_EN[0] + d[0] * W_EN[1] > 0 else (d[1], -d[0])
    B, F = W_P(W_L, 0), W_P(W_L, W_T); B2 = (B[0] + W_EIN * d[0], B[1] + W_EIN * d[1]); F2 = (B2[0] + W_T * nd[0], B2[1] + W_T * nd[1])
    return pl.union(fl).union(_WP([B, B2, F2, F]).buffer(0)).buffer(0.001, join_style=2).buffer(-0.001, join_style=2)
W_PROFIL = _w_profil()
W_SCHLITZ = _WP([W_P(W_SCHLITZ_S, -1), W_P(W_L + 6, -1), W_P(W_L + 6, W_T + 6), W_P(W_SCHLITZ_S, W_T + 6)])
W_LOCH_Y = _wyf0 + W_FL / 2.0
# --- Hebel (u, Y) nach 4685 Lagerung Knick (Reihen 1-3) / 468401 V3 (Reihe 4): Zeichnungsrahmen x entlang der Rueckwand
# nach unten, y von der Wand weg; Bohrung (Achse) im Ursprung, Knopfarm y 4,1..10,9, Knick bei x 14,01 bzw. 21,21.
# 22.11. 14:57 (Johann): Lasche so nah wie moeglich an die Rueckwand, Hoehe nach den Klappen der 1. Reihe; der 101-Grad-
# Winkel wird je Reihe so angepasst, dass der Knopfarm unter dem Zylinder endet; Ende 2 mm schmal (Zapfen im Zylinder).
W_HEBEL_ZG = {"knick": [(-22, 4.1), (14.01, 4.1), (16.73, 19.57), (10.04, 20.75), (8.3, 10.9), (-22, 10.9)],
              "v3": [(-22, 4.1), (21.21, 4.1), (23.93, 19.57), (17.24, 20.75), (15.5, 10.9), (-22, 10.9)]}
W_B_HEBEL = 7.7                                          # Hebelbreite (Johann 22.11.)
W_N_ACHSE = W_T + 0.5 + 3.0                              # Achse: Ms-Lasche (R3) 0,5 vor der Platte
W_HEBEL_ALT = {r_: _WP(HEBEL[r_]).buffer(0) for r_ in (1, 2, 3, 4)}
def _w_rahmen(A): return lambda x, y: (A[0] - x * W_ES[0] + y * W_EN[0], A[1] - x * W_ES[1] + y * W_EN[1])
def _w_unterkante(A, Y):
    M = _w_rahmen(A); p = W_HEBEL_ZG["knick"]; a = M(*p[1]); b = M(*p[2]); t = (Y - a[1]) / (b[1] - a[1]); return a[0] + t * (b[0] - a[0])
_wy1 = FUEHR[1]["y_lane"][0] - 10.0
_wsoll = W_HEBEL_ALT[1].intersection(_Wbox(-500, _wy1, 500, _wy1 + 0.5)).bounds[0]     # Unterkante Klappenarm Reihe 1
_lo, _hi = 0.0, 40.0
for _ in range(50):
    _m = (_lo + _hi) / 2
    if _w_unterkante(W_P(_m, W_N_ACHSE), _wy1) < _wsoll: _lo = _m
    else: _hi = _m
W_S_ACHSE = _lo
# 22.11. 15:12 (Johann): Achsbohrung am Winkel 23 mm ueber dem Boden, 10 mm von der Rueckwand (senkrecht zur Wand)
W_N_ACHSE = 10.0
W_S_ACHSE = (23.0 - W_N_ACHSE * W_EN[0] - W_O[0]) / W_ES[0]
W_ACHSE = W_P(W_S_ACHSE, W_N_ACHSE)
_wA = W_ACHSE
W_LASCHE = _WL([W_P(4.0, W_T), W_P(29.0, W_T)]).union(_WPt(_wA).buffer(4.0, 64)).convex_hull.union(
    _WP([W_P(4.0, W_T - 0.01), W_P(29.0, W_T - 0.01), W_P(29.0, W_T + 0.5), W_P(4.0, W_T + 0.5)])).buffer(0)
W_LASCHE_MS, W_WINKEL_R, W_DD, W_NF, W_XKNICK = {}, {}, {}, {}, {}
# 22.11. 16:16-16:39 (Johann, Einzelhebel R1/R4 v4 ok): Knopfarm 7,7 x 2 gerade von der Lasche zum Armende der verdeckten
# Variante (Uebergang Zylinder unveraendert); Lasche (4685 Knick, Reihe 4: 468401 V3) ohne Luft an der inneren Knickecke,
# V3 9 mm weiter vom Knick; Klappenarm 6 breit, Unterkante gerade (kein Knick), zur Klappe hin auf 2 mm verjuengt;
# Zylinder bis zur Hebelkante, die unter dem ganzen Zylinderquerschnitt am weitesten hinten liegt
W_B_KLAPPE, W_B_ENDE, W_LUFT_KNICK, W_V3_ZUSATZ = 6.0, 2.0, 0.0, 9.0
def _w_schnitt(p, r, q, s):
    den = r[0] * s[1] - r[1] * s[0]; t = ((q[0] - p[0]) * s[1] - (q[1] - p[1]) * s[0]) / den; return (p[0] + t * r[0], p[1] + t * r[1])
for r_ in (1, 2, 3, 4):
    alt_ = W_HEBEL_ALT[r_]; KN_ = (knopf_u[r_ - 1][0] + (JP_KNOPF_D[r_][0] if JP else 0), knopf_u[r_ - 1][1] + (JP_KNOPF_D[r_][1] if JP else 0))
    uk_ = KOLBEN_UNTEN[r_][0] - 6.0
    end_ = alt_.intersection(_Wbox(uk_, -500, 500, 500)); end_ = max(end_.geoms, key=lambda q: q.area) if hasattr(end_, "geoms") else end_
    sl_ = alt_.intersection(_Wbox(uk_, -500, uk_ + 0.3, 500)).bounds
    Pk = (uk_ + 0.15, (sl_[1] + sl_[3]) / 2)
    yc_ = 4.1 + W_B_HEBEL / 2
    v_ = (Pk[0] - _wA[0], Pk[1] - _wA[1]); L_ = _wm.hypot(*v_); v_ = (v_[0] / L_, v_[1] / L_); a_ = _wm.asin(yc_ / L_)
    best_ = None
    for sg_ in (1, -1):
        c_, s_ = _wm.cos(sg_ * a_), _wm.sin(sg_ * a_)
        dd_ = (v_[0] * c_ - v_[1] * s_, v_[0] * s_ + v_[1] * c_)
        nf_ = (-dd_[1], dd_[0]) if (-dd_[1]) * W_EN[0] + dd_[0] * W_EN[1] > 0 else (dd_[1], -dd_[0])
        lot_ = (_wA[0] + yc_ * nf_[0], _wA[1] + yc_ * nf_[1])
        err_ = abs((Pk[0] - lot_[0]) * nf_[0] + (Pk[1] - lot_[1]) * nf_[1])
        if best_ is None or err_ < best_[0]: best_ = (err_, dd_, nf_)
    _, DD_, NF_ = best_
    W_DD[r_], W_NF[r_] = DD_, NF_
    M = lambda x, y, DD_=DD_, NF_=NF_: (_wA[0] - x * DD_[0] + y * NF_[0], _wA[1] - x * DD_[1] + y * NF_[1])
    b1_ = alt_.intersection(_Wbox(-500, 60, 500, 60.5)).bounds
    Y2_ = min(120.0, alt_.bounds[3] - 15.0); b2_ = alt_.intersection(_Wbox(-500, Y2_, 500, Y2_ + 0.5)).bounds
    kdir = (b2_[0] - b1_[0], Y2_ - 60.0); _l = _wm.hypot(*kdir); kdir = (kdir[0] / _l, kdir[1] / _l)
    kn_ = (-kdir[1], kdir[0]) if -kdir[1] > 0 else (kdir[1], -kdir[0])
    def _ecken(xk):
        ka_ = M(xk, 4.1)
        ki_ = _w_schnitt(M(0, 4.1 + W_B_HEBEL), (-DD_[0], -DD_[1]), (ka_[0] + W_B_KLAPPE * kn_[0], ka_[1] + W_B_KLAPPE * kn_[1]), kdir)
        return ka_, ki_, -((ki_[0] - _wA[0]) * DD_[0] + (ki_[1] - _wA[1]) * DD_[1])
    lo_, hi_ = 0.0, 40.0
    for _ in range(60):
        m_ = (lo_ + hi_) / 2
        if _ecken(m_)[2] < 6.75 + W_LUFT_KNICK: lo_ = m_
        else: hi_ = m_
    xk_ = lo_ + (W_V3_ZUSATZ if r_ == 4 else 0.0); W_XKNICK[r_] = xk_
    ka, ki, _xi = _ecken(xk_)
    ymax_ = alt_.bounds[3]; tend_ = (ymax_ - ka[1]) / kdir[1]
    klappe_ = _WP([ka, (ka[0] + tend_ * kdir[0], ka[1] + tend_ * kdir[1]),
                   (ka[0] + tend_ * kdir[0] + W_B_ENDE * kn_[0], ka[1] + tend_ * kdir[1] + W_B_ENDE * kn_[1]), ki])
    xo_ = -((Pk[0] - _wA[0]) * DD_[0] + (Pk[1] - _wA[1]) * DD_[1]) - 3.0
    stab_ = _WP([M(xo_, 4.1), ka, ki, M(xo_, 4.1 + W_B_HEBEL)])
    neu_ = _Wuu([stab_, klappe_, end_]).buffer(0.05).buffer(-0.05)
    neu_ = max(neu_.geoms, key=lambda q: q.area) if hasattr(neu_, "geoms") else neu_
    HEBEL[r_] = [(float(a_), float(b_)) for a_, b_ in list(neu_.simplify(0.02).exterior.coords)[:-1]]
    LAGER_R[r_] = (float(_wA[0]), float(_wA[1]))
    W_WINKEL_R[r_] = round(180.0 - _wm.degrees(_wm.acos(-(DD_[0] * kdir[0] + DD_[1] * kdir[1]))), 1)
    # Zylinder bis zur hintersten Hebelkante unter dem ganzen Zylinderquerschnitt
    _ax = _WL([KN_, (KN_[0] - 30 * KD[0], KN_[1] - 30 * KD[1])]).buffer(3.75, cap_style=2).intersection(neu_)
    _pp = [c_ for g_ in (getattr(_ax, "geoms", None) or [_ax]) if not g_.is_empty for c_ in g_.exterior.coords]
    _t = max((KN_[0] - c_[0]) * KD[0] + (KN_[1] - c_[1]) * KD[1] for c_ in _pp)
    KOLBEN_UNTEN[r_] = (KN_[0] - _t * KD[0], KN_[1] - _t * KD[1])
# Ms-Lasche nach 4685 (Ansicht LASCHE, Schnitt): Schenkel 0,5 beidseits, oben ueber die Hebelkante gefalzt (Kerbe R2),
# unter y 4,9 nach aussen gekroepft (lichte Weite 3,8), Bohrung Ø3 auf der Achse, Niete Ø3,2 bei (+-4 | 7,4)
import ezdxf as _wez
_wlas = [[(float(x_), float(y_)) for x_, y_, *_r in e_.get_points()] for e_ in _wez.readfile("lasche_4685_knick.dxf").modelspace().query('LWPOLYLINE[layer=="LASCHE"]')][0]
W_LAS_PG = _WP(_wlas).buffer(0)
def _w_platte(pg, z0, z1):
    erg = None
    for g in (pg.geoms if hasattr(pg, "geoms") else [pg]):
        if g.area < 0.01: continue
        t_ = cq.Workplane("XY").polyline([(float(a_), float(b_)) for a_, b_ in list(g.exterior.coords)[:-1]]).close().extrude(z1 - z0).translate((0, 0, z0)).val()
        erg = t_ if erg is None else erg.fuse(t_)
    return erg
def _w_lasche_lokal():
    oben, unten = W_LAS_PG.intersection(_Wbox(-20, 4.9, 20, 20)), W_LAS_PG.intersection(_Wbox(-20, -10, 20, 4.1))
    ms = None
    for sg in (1, -1):
        for pg, z0, z1 in ((oben, 1.0, 1.5), (unten, 1.9, 2.4)):
            p_ = _w_platte(pg, z0, z1) if sg > 0 else _w_platte(pg, -z1, -z0)
            ms = p_ if ms is None else ms.fuse(p_)
        kr = cq.Workplane("YZ").polyline([(4.1, 1.9 * sg), (4.9, 1.0 * sg), (4.9, 1.5 * sg), (4.1, 2.4 * sg)]).close().extrude(11.2).translate((-5.6, 0, 0)).val()
        ms = ms.fuse(kr.intersect(_w_platte(W_LAS_PG, -3, 3)))
    ms = ms.fuse(_w_platte(W_LAS_PG.intersection(_Wbox(-20, 10.9, 20, 20)), -1.5, 1.5))
    ms = ms.cut(cq.Solid.makeCylinder(1.5, 10, cq.Vector(0, 0, -5)))
    for xn in (-4.0, 4.0): ms = ms.cut(cq.Solid.makeCylinder(1.6, 10, cq.Vector(xn, 7.4, -5)))
    return ms.clean()
W_LASCHE_LOKAL = _w_lasche_lokal()
W_NIET_LOKAL = [cq.Solid.makeCylinder(1.6, 3.0, cq.Vector(xn, 7.4, -1.5)).fuse(cq.Solid.makeCylinder(2.6, 0.6, cq.Vector(xn, 7.4, 1.5))).fuse(
    cq.Solid.makeCylinder(2.6, 0.6, cq.Vector(xn, 7.4, -2.1))) for xn in (-4.0, 4.0)]
W_NIETLOCH_LOKAL = [cq.Solid.makeCylinder(1.6, 10, cq.Vector(xn, 7.4, -5)) for xn in (-4.0, 4.0)]
def W_LOKAL(sol, r, xs):     # Zeichnungsrahmen (x entlang Stab nach unten, y nach vorne, z = X) -> Modell
    DD_, NF_ = W_DD[r], W_NF[r]
    return sol.transformGeometry(cq.Matrix([[0.0, 0.0, 1.0, xs], [-DD_[1], NF_[1], 0.0, _wA[1]], [DD_[0], -NF_[0], 0.0, -_wA[0]]]))
print("Winkel-Hebel: Winkel Knopfarm/Klappenarm je Reihe", W_WINKEL_R, "Knick unter der Bohrung", {r_: round(v_, 2) for r_, v_ in W_XKNICK.items()}, "Zylinder unten", {r_: tuple(round(c_, 1) for c_ in v_) for r_, v_ in KOLBEN_UNTEN.items()})
print("Winkel-Griffstock: Achse Ø3 bei u %.2f Y %.2f (Platte s %.1f, %.1f vor der Platte), Platte ab Y %.2f (Deckel innen)" % (_wA[0], _wA[1], W_S_ACHSE, W_N_ACHSE - W_T, W_O[1]))
if os.environ.get("WINKEL_DEBUG"):
    import pickle as _wpk
    _wpk.dump({"neu": {r_: _WP(HEBEL[r_]) for r_ in (1, 2, 3, 4)}, "alt": W_HEBEL_ALT, "profil": W_PROFIL, "lasche": W_LASCHE, "achse": _wA,
               "schlitz": W_SCHLITZ, "kontakt": KOLBEN_KONTAKT, "unten": KOLBEN_UNTEN}, open("/tmp/winkel_dbg.pkl", "wb"))
    raise SystemExit("WINKEL_DEBUG")
# 22.11. 19:30 (Johann): jeder Hebel bekommt eine Spiralfeder, 3 Windungen, Draht 0,6, zwischen den Laschenschenkeln
# neben der Achse; ein Schenkel endet mittig am 2-mm-Hebel, abgewinkelt 2 mm tief in einer Bohrung Ø1, das andere Ende
# liegt in der Kerbe (Schlitz) der Randverzahnung des Winkels
W_D_DRAHT, W_FED_RM, W_FED_STEIG, W_FED_WDG = 0.6, 1.5 + 0.15 + 0.3, 0.65, 3.0
W_FED_X_LOCH = -9.0                                     # Bohrung Ø1 in der Hebelkante (zur Achse hin), oberhalb der Lasche
def _w_draht(p, q):
    v = q.sub(p); L = v.Length
    return cq.Solid.makeCylinder(W_D_DRAHT / 2, L, p, v.normalized()) if L > 1e-6 else None
def W_FEDER_LOKAL(r):
    """Feder und Bohrwerkzeug im Zeichnungsrahmen der Reihe r (x entlang Stab nach unten, y nach vorne, z = X)"""
    DD_, NF_ = W_DD[r], W_NF[r]
    Nu, NY = W_P(W_SCHLITZ_S + W_D_DRAHT / 2 + 0.05, W_T / 2)   # Kerbe: Draht liegt am Grund des Schlitzes der Randverzahnung, Plattenmitte
    nx = -((Nu - _wA[0]) * DD_[0] + (NY - _wA[1]) * DD_[1]); ny = (Nu - _wA[0]) * NF_[0] + (NY - _wA[1]) * NF_[1]
    p1 = (W_FED_X_LOCH, 4.1 - 0.9)                         # Schenkel 1 vor der Bohrung, dann 2 mm in die Bohrung
    _wx = -((W_P(W_SCHLITZ_S, W_T + 1.2)[0] - _wA[0]) * DD_[0] + (W_P(W_SCHLITZ_S, W_T + 1.2)[1] - _wA[1]) * DD_[1]); _wy = (W_P(W_SCHLITZ_S, W_T + 1.2)[0] - _wA[0]) * NF_[0] + (W_P(W_SCHLITZ_S, W_T + 1.2)[1] - _wA[1]) * NF_[1]
    th1 = math.atan2(p1[1], p1[0]); th2 = math.atan2(_wy, _wx)
    wdg = W_FED_WDG + ((th1 - th2) % (2 * math.pi)) / (2 * math.pi)
    h = W_FED_STEIG * wdg
    hel = cq.Wire.makeHelix(W_FED_STEIG, h, W_FED_RM, cq.Vector(0, 0, -h / 2), cq.Vector(0, 0, 1))
    hel = hel.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), math.degrees(th2))
    prof = cq.Wire.makeCircle(W_D_DRAHT / 2, hel.startPoint(), hel.tangentAt(0))
    fed = cq.Solid.sweep(prof, [], hel, True, True)
    S, E = hel.startPoint(), hel.endPoint()
    pk = [E, cq.Vector(p1[0], p1[1], 0.0), cq.Vector(p1[0], 4.1 + 2.0 - W_D_DRAHT / 2, 0.0)]
    Wu, WY = W_P(W_SCHLITZ_S + W_D_DRAHT / 2 + 0.05, W_T + 1.2)   # vor der Platte ueber der Kerbe, dann gerade in die Kerbe
    wx = -((Wu - _wA[0]) * DD_[0] + (WY - _wA[1]) * DD_[1]); wy = (Wu - _wA[0]) * NF_[0] + (WY - _wA[1]) * NF_[1]
    ps = [S, cq.Vector(wx, wy, 0.0), cq.Vector(nx, ny, 0.0)]
    for pts in (pk, ps):
        for a_, b_ in zip(pts, pts[1:]):
            d_ = _w_draht(a_, b_)
            if d_ is not None: fed = fed.fuse(d_)
        for q_ in pts[1:]: fed = fed.fuse(cq.Solid.makeSphere(W_D_DRAHT / 2 + 0.01, q_, angleDegrees1=-90, angleDegrees2=90))
    loch = cq.Solid.makeCylinder(0.5, 2.3, cq.Vector(W_FED_X_LOCH, 4.1 - 0.1, 0.0), cq.Vector(0, 1, 0))
    return fed.clean(), loch
W_FEDER = {r_: W_FEDER_LOKAL(r_) for r_ in (1, 2, 3, 4)}
print("Winkel-Federn: Spirale %.1f Windungen (Reihe 1), Bohrung Ø1 x 2 bei x %.1f" % (W_FED_WDG + 0, W_FED_X_LOCH))
