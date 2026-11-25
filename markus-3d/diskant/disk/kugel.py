# Kugellager-Griffstock (Johann 21.11.2026, Unterlagen C:\0-A-kugellager): Traeger aus POM-Stuetzen (4 mm) in zwei
# FR2-Leisten (2 mm), je Hebel ein Holzhebel mit zwei Kugellagern (5x2x2) nebeneinander auf einer kurzen 2-mm-Achse,
# die in den beiden Nachbarstuetzen sitzt; Alu-Schenkel seitlich in der 2-mm-Vertiefung des Holzhebels eingerastet
# (Reihen 1/2 und 3/4 auf gegenueberliegenden Seiten). Zwei Holzkaemme (Teil 1 hinten/oben fuer die Reihen 1/2, Teil 2
# vorne fuer die Reihen 3/4) aus dem alten Schnitt beide-Teile11a+schnitt.dxf (= kugel_schnitt.dxf), Hebel aus
# Hebel-kombination.dxf (= kugel_hebel.dxf).
# Lage: der alte Schnitt liegt um 180 Grad gedreht genau auf den Knoepfen des Modells (Abweichung < 0,15 mm).
# Hebelwinkel: Knopfrueckseite liegt auf der seitlichen Kopfflaeche des Holzhebels; die Kopfflaeche ist um den Hubwinkel
# geneigt, so dass der Knopf gedrueckt flach (parallel) auf seiner Stufe liegt. Wo ein Holzhebel das fremde Kammteil
# (Stufe) beruehrt, wird nur der duenne Schenkel verlaengert.
# Stuetzen und Leisten variabel: Stuetzen mittig zwischen den Hebel-X (halbe Knopfteilung), Leisten ueber alle Stuetzen.
# Wird am Ende von hebel.py ausgefuehrt (griffstock = kugellager) und ersetzt dort die Kunststoff-Lagerung.
import numpy as _np
from shapely.geometry import Polygon as _TP, Point as _TPt, LineString as _TL, box as _tbox
from shapely import affinity as _taff
from shapely.ops import unary_union as _tuu
from ezdxf import path as _ezp
import schnitt_trafo as _TST

_TDX = 58.17                                              # Schnittkoordinaten wie trad_schnitt.dxf (= Schnitt_Hebel_4_Layer + 58,17)
def _M(p):
    u_, Y_ = _TST.modell((p[0] - _TDX, p[1])); return (float(u_), float(Y_))
_src = _np.array([[0.0, 0.0], [100.0, 0.0], [0.0, 100.0]])
_AFF = _np.linalg.solve(_np.c_[_src, _np.ones(3)], _np.array([_M(p_) for p_ in _src]))
def _Minv(q):
    return tuple(float(v_) for v_ in _np.linalg.solve(_AFF[:2, :].T, _np.array(q) - _AFF[2]))
def _Mpg(pg):
    pg = max(pg.geoms, key=lambda g_: g_.area) if hasattr(pg, "geoms") else pg
    return _TP([_M(p_) for p_ in list(pg.exterior.coords)[:-1]], [[_M(p_) for p_ in list(h_.coords)[:-1]] for h_ in pg.interiors]).buffer(0)
def _pp(pg):
    pg = pg.simplify(0.02).buffer(0)
    pg = max(pg.geoms, key=lambda g_: g_.area) if hasattr(pg, "geoms") else pg
    out = []
    for a_, b_ in list(pg.exterior.coords)[:-1]:
        if not out or math.hypot(a_ - out[-1][0], b_ - out[-1][1]) > 0.02: out.append((float(a_), float(b_)))
    while len(out) > 3 and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) <= 0.02: out.pop()
    return out
def _tteile(pg):
    return [g_ for g_ in (pg.geoms if hasattr(pg, "geoms") else [pg]) if g_.geom_type == "Polygon" and g_.area > 0.05]
def _gross(pg): return max(pg.geoms, key=lambda g_: g_.area) if hasattr(pg, "geoms") else pg
def _rot(pg, wink, um):
    return _taff.rotate(pg, wink, origin=um, use_radians=False)
def _Ub(u0, u1=500.0):                                     # Bereich u0 <= u <= u1 (Hoehe ueber der Deckelflaeche) im Schnitt
    return _TP([_Minv((u0, -500.0)), _Minv((u1, -500.0)), _Minv((u1, 500.0)), _Minv((u0, 500.0))])

# --- Knoepfe des Modells im Schnitt ---
KN = {}
for r_ in (1, 2, 3, 4):
    ku_, kY_ = knopf_u[r_ - 1]
    if JP: ku_, kY_ = ku_ + JP_KNOPF_D[r_][0], kY_ + JP_KNOPF_D[r_][1]
    KN[r_] = _Minv((ku_, kY_))                            # Knopfrueckseite (Hebelseite), Knopf ragt nach -x

# --- alter Schnitt: Knopf-Einfuegepunkte (Reihe 1 = kleinstes y) -> starre Einpassung (Drehung + Verschiebung) ---
_ks = ezdxf.readfile("kugel_schnitt.dxf").modelspace()
_oldk = sorted([(float(e.dxf.insert.x), float(e.dxf.insert.y)) for e in _ks.query("INSERT")], key=lambda q_: q_[1])
_O = _np.array(_oldk); _N = _np.array([KN[r_] for r_ in (1, 2, 3, 4)])
_co, _cn = _O.mean(0), _N.mean(0)
_U, _S, _Vt = _np.linalg.svd((_O - _co).T @ (_N - _cn)); _RK = _Vt.T @ _U.T
if _np.linalg.det(_RK) < 0: raise SystemExit("Kugellager: alter Schnitt waere gespiegelt")
_tK = _cn - _RK @ _co
def _KT(p): return tuple(float(v_) for v_ in _RK @ _np.array(p, float) + _tK)
_res = max(float(_np.hypot(*(_np.array(_KT(o_)) - n_))) for o_, n_ in zip(_O, _N))
print(f"Kugellager: alter Schnitt auf die Knoepfe gelegt, Drehung {math.degrees(math.atan2(_RK[1, 0], _RK[0, 0])):.1f} Grad, "
      f"groesste Abweichung {_res:.2f} mm")
# Stuetze (Stuetze-fuer-Hebel.dxf, Teil B; im alten Schnitt um +96,21/+71,51 verschoben) mit den Achsbohrungen r 1,025
_DST = (96.21, 71.51)
def _B(p): return _KT((p[0] + _DST[0], p[1] + _DST[1]))
_st = [e for e in _ks.query('LWPOLYLINE[layer=="0"]') if len(list(e.get_points())) == 32][0]
STUETZE = _TP([_KT((x_, y_)) for x_, y_, *_r in _st.get_points()]).buffer(0)
ACHSE_K = {1: _B((-32.52, 33.70)), 2: _B((-32.52, 28.70)), 3: _B((-16.50, 15.56)), 4: _B((-16.50, 10.56))}
# 22.11. 07:54 (Johann): Bohrungen aus dem CamBam-Plan Stuetze-fuer-hebel.cb (Bohren1, in Teil-B-Koordinaten):
# neben jeder Achsbohrung eine zweite Ø2 um 1,436 versetzt - oeffnet die Achsbohrung, damit die Achse einschnappt;
# dazu Ø1,51-Freibohrungen in den Innenecken der Kontur
_ST_OEFFNUNG = [_TPt(_B((x_ + 1.436, y_))).buffer(1.0, 32) for x_, y_ in ((-32.52, 33.70), (-32.52, 28.70), (-16.50, 15.56), (-16.50, 10.56))]
_ST_ECKEN = [_TPt(_B(q_)).buffer(0.755, 32) for q_ in ((-25.543, 41.72), (-36.4, 39.599), (-26.715, 7.534), (-36.647, 4.586), (-42.544, 33.11), (-42.549, 8.054))]
STUETZE = _gross(STUETZE.difference(_tuu(_ST_OEFFNUNG + _ST_ECKEN)))
print(f"Kugellager: Stuetze nach CamBam-Plan - 4 Schnappoeffnungen Ø2 an den Achsbohrungen, 6 Eckfreibohrungen Ø1,51 (Flaeche {STUETZE.area:.1f} mm2)")
# Leiste2 (32 breit, Schlitze 24 lang) an der Stuetzen-Kante, Leiste (34 breit, Schlitze 10 lang bei 7..17, Loecher
# Ø3,5 bei 25,18) auf der schraegen Zunge
# 22.11. 03:04 (Johann): beide Epoxy-Leisten so duenn, dass die Nasen der Stuetzen 0,1 mm durch die Platte reichen
_NASE2 = 1.5                                              # seitliche Nase (x -43,51 .. -42,01 in Stuetze B)
T_L2 = _NASE2 - 0.1
LEISTE2 = _TP([_B(q_) for q_ in ((-42.01 - T_L2, 3.58), (-42.01, 3.58), (-42.01, 35.58), (-42.01 - T_L2, 35.58))])
_p0_, _p1_ = _np.array((-42.01, 37.97)), _np.array((-35.77, 39.19))
_ud = (_p1_ - _p0_) / _np.linalg.norm(_p1_ - _p0_); _nd = _np.array([-_ud[1], _ud[0]])
_NASE1 = min(float((_np.array(q_) - _p0_) @ _nd) for q_ in ((-36.09, 40.65), (-26.27, 42.56)))   # Hoehe der Zunge ueber der Kante
T_L1 = _NASE1 - 0.1
LEISTE1 = _TP([_B(tuple(_p1_ + _ud * a_ + _nd * b_)) for a_, b_ in ((-7, 0), (27, 0), (27, T_L1), (-7, T_L1))])
LEISTE1_LOCH = _B(tuple(_p1_ + _ud * 18.18 + _nd * 1.0))
# Holzkaemme aus dem alten Schnitt: Umriss, Teil 2 (vorne, Reihen 3/4), Teil 1 = Rest
_GES_O = [(50.76, -10.13), (71.96, -10.13), (71.96, 16.84), (75.96, 17.0), (75.96, 33.0), (79.96, 33.0), (79.96, 49.0),
          (83.96, 49.0), (83.96, 117.01), (50.76, 110.56)]
_T2_O = [(68.96, 33.0), (79.96, 33.0), (79.96, 49.0), (83.96, 49.0), (83.96, 117.01), (68.96, 114.1)]
KAMM_GES = _TP([_KT(q_) for q_ in _GES_O]).buffer(0)
KAMM2 = _TP([_KT(q_) for q_ in _T2_O]).buffer(0)
KAMM1 = _gross(KAMM_GES.difference(KAMM2))
# 22.11. 09:31 (Johann): hauchduenner Rest von Teil 1 auf der Oberseite von Teil 2 (Stufe 3 -> Stufe 2): die Oberkanten
# der alten Zeichnung sind nicht genau kollinear - Spaene unter 0,2 mm wegnehmen, Ecken bleiben scharf
KAMM1 = _gross(KAMM1.buffer(-0.1, join_style=2).buffer(0.1, join_style=2))
# 22.11. (Johann): die untere Epoxy-Leiste (34) bleibt, sie ist mit dem Fuss der Stuetzen in den Diskantboden eingelassen
# (im Modell liegt die Deckelflaeche 3 mm hoeher als beim alten Instrument); Kaemme enden an der Deckelflaeche
STUETZE_VOLL = STUETZE
KAMM_GES = _gross(KAMM_GES.intersection(_Ub(0.0))); KAMM1 = _gross(KAMM1.intersection(_Ub(0.0))); KAMM2 = _gross(KAMM2.intersection(_Ub(0.0)))
LEISTE1 = _gross(LEISTE1.intersection(_tbox(KAMM_GES.bounds[0] + 0.3, -500, 500, 500)))   # vorne buendig mit dem Kamm
# 22.11. 03:24 (Johann): die ganze Mechanik (ohne Knoepfe und Kammteile) parallel zur Rueckwand von Teil 1 anheben, bis
# die untere Leiste auf dem Diskantboden aufsteht
_v_rw = _RK @ _np.array([0.0, -1.0])                      # Rueckwand Teil 1 (alter Schnitt x = 50,76) nach oben
_du_ds = float(_AFF[0, 0] * _v_rw[0] + _AFF[1, 0] * _v_rw[1])
HUB_MECH = -min(_M(q_)[0] for q_ in LEISTE1.exterior.coords) / _du_ds
_dv = tuple(_v_rw * HUB_MECH)
STUETZE = _taff.translate(STUETZE, *_dv); STUETZE_VOLL = STUETZE
LEISTE1 = _taff.translate(LEISTE1, *_dv); LEISTE2 = _taff.translate(LEISTE2, *_dv)
ACHSE_K = {r_: (a_[0] + _dv[0], a_[1] + _dv[1]) for r_, a_ in ACHSE_K.items()}
print(f"Kugellager: Mechanik parallel zur Rueckwand von Teil 1 um {HUB_MECH:.2f} mm angehoben - untere Leiste steht auf dem Diskantboden")
_TRAEGER = _tuu([STUETZE, LEISTE1, LEISTE2]).buffer(0)
# Verdeck (Deckplatte u 28,5-32 und Verdeckleiste u 32-37) vor der Kammfront: dort darf kein Hebelteil hinaus
_ZV = _Ub(27.5, 38.5).intersection(_tbox(-500, -500, KAMM_GES.bounds[0] + 0.3, 500))
def _frei_machen(pg):                                     # Alu-Schenkel ueber dem Boden (das Verdeck wird hinten hoeher gesetzt)
    return _gross(pg.intersection(_Ub(0.4)))

# --- Hebel aus kugel_hebel.dxf: Holz T4/T3/T13/T10, Alu T8/T7/T14/T11, Vertiefungen T1/T2/T12/T9 ---
_hm = ezdxf.readfile("kugel_hebel.dxf").modelspace()
_hpl = []
for e in _hm.query("LWPOLYLINE"):
    _q = [(v_.x, v_.y) for v_ in _ezp.make_path(e).flattening(0.02)]
    _g = _TP(_q).buffer(0)
    if _g.area < 0.9 * abs(_TP(_q).area):                # Selbstueberschneidung durch doppelte Punkte (z. B. Vertiefung R2)
        from shapely.validation import make_valid as _mv
        _g = _gross(_tuu([x_ for x_ in getattr(_mv(_TP(_q)), "geoms", [_mv(_TP(_q))]) if x_.geom_type == "Polygon"]))
    _hpl.append((e.dxf.layer, len(list(e.get_points())), _g))
def _hk(lay, n, near=None):
    c_ = [p_ for l_, k_, p_ in _hpl if l_ == lay and k_ == n]
    if near is not None: c_.sort(key=lambda p_: p_.distance(_TPt(near)))
    return _gross(c_[0])
H_PIV = {1: (128.81, 31.82), 2: (167.65, -10.47), 3: (88.79, 98.84), 4: (46.11, 162.585)}
H_HOLZ = {1: _hk("HEBEL-ALU-1R", 16), 2: _hk("HEBEL-ALU-2R", 14), 3: _hk("HEBEL-ALU-3R", 11), 4: _hk("HEBEL-ALU-4R", 9)}
H_ALU = {1: _hk("AM_5", 31), 2: _hk("AM_5", 46), 3: _hk("HEBEL-ALU-3R", 16), 4: _hk("HEBEL-ALU-4R", 23)}
H_VERT = {1: _hk("LAYER9", 11, H_PIV[1]), 2: _hk("LAYER9", 11, H_PIV[2]), 3: _hk("0", 5, H_PIV[3]), 4: _hk("0", 9, H_PIV[4])}
# 22.11. 04:39 (Johann): die Vertiefung im Holz hat dieselbe Form wie ein Teil des Alu-Hebels -> Alu-Kontur starr auf die
# Vertiefung legen (ICP: Vertiefungsrand auf den naechsten Alu-Rand, 70 % der besten Paare)
from scipy.spatial import cKDTree as _KD
def _rand(pg, n):
    L_ = pg.exterior
    return _np.array([L_.interpolate(i_ / n, normalized=True).coords[0] for i_ in range(n)])
def _alu_einpassen(Al, Vt):
    P = _rand(Vt, 240); best = None
    for w0 in range(-30, 31, 5):
        A_ = _taff.rotate(Al, w0, origin=Vt.centroid)
        M_ = _np.eye(3); M_[:2, :2] = _np.array([[math.cos(math.radians(w0)), -math.sin(math.radians(w0))], [math.sin(math.radians(w0)), math.cos(math.radians(w0))]])
        c_ = _np.array(Vt.centroid.coords[0]); M_[:2, 2] = c_ - M_[:2, :2] @ c_
        for _it in range(80):
            Q = _rand(A_, 1500); d_, i_ = _KD(Q).query(P); keep = d_ <= _np.quantile(d_, 0.7)
            src_, dst_ = Q[i_[keep]], P[keep]          # Alu-Punkte -> Vertiefungspunkte
            cs_, cd_ = src_.mean(0), dst_.mean(0)
            U_, S_, Vt_ = _np.linalg.svd((src_ - cs_).T @ (dst_ - cd_)); R_ = Vt_.T @ U_.T
            if _np.linalg.det(R_) < 0: break
            t_ = cd_ - R_ @ cs_
            A_ = _taff.affine_transform(A_, [R_[0, 0], R_[0, 1], R_[1, 0], R_[1, 1], t_[0], t_[1]])
            T_ = _np.eye(3); T_[:2, :2] = R_; T_[:2, 2] = t_; M_ = T_ @ M_
            if _np.hypot(*t_) < 1e-4 and abs(math.atan2(R_[1, 0], R_[0, 0])) < 1e-6: break
        d_, _i = _KD(_rand(A_, 1500)).query(P); sc_ = float(_np.mean(_np.sort(d_)[:int(0.7 * len(d_))]))
        if A_.intersection(Vt).area < 0.85 * Vt.area: continue   # Alu muss die Vertiefung ausfuellen
        if best is None or sc_ < best[0]: best = (sc_, A_, M_)
    return best
# 22.11. 05:47 (Johann): Alu ganz nach Zeichnung - Nabe bis kurz nach der Kerbe (Nase) unveraendert, dann gerades Band
# (Winkel/Laenge so korrigiert, dass die Endform auf der Klappe liegt), Endform der Zeichnung unveraendert
H_KERBE = {1: (104.0, 50.8), 2: (142.0, 2.5), 3: (84.0, 89.5), 4: (38.0, 136.6)}     # Bandmitte kurz nach der Kerbe/dem Knick
H_ENDE = {1: (-27.0, 53.8), 2: (39.0, 1.5), 3: (28.0, 78.9), 4: (-12.0, 131.5)}      # Bandmitte am Beginn der Endform
for r_ in (1, 2, 3, 4):
    _sc, _A, _Mx = _alu_einpassen(H_ALU[r_], H_VERT[r_])
    print(f"Kugellager Reihe {r_}: Alu auf die Vertiefung eingepasst - Drehung {math.degrees(math.atan2(_Mx[1, 0], _Mx[0, 0])):+.2f} Grad, "
          f"Verschiebung {_np.hypot(_Mx[0, 2] - 0, 0):.0f}, mittlere Abweichung {_sc:.3f} mm")
    H_ALU[r_] = _A
H_FLAECHE = {1: ((135.88, 105.05), (141.69, 130.61)), 2: ((166.68, 59.07), (167.88, 68.94)),
             3: ((87.28, 134.57), (87.59, 149.13)), 4: ((41.59, 183.31), (42.79, 191.62))}   # Knopf-Anlageflaeche am Kopf
_THUB = 3.6
T_HZ = 4.0                                                # Holzhebel = zwei Lager 5x2x2 nebeneinander
D_LAGER, B_LAGER, D_ACHSE_K = 5.0, 2.0, 2.0
def _verlaengern(W, p, face, L, s_rel=0.5):
    """Holzhebel im duennen Schenkel um L verlaengern (Kopf wandert entlang der Schenkelachse)."""
    if abs(L) < 1e-9: return W, face
    f0 = _np.array(face[0]); v = f0 - _np.array(p); s = _np.linalg.norm(v); v = v / s
    c = _np.array(p) + v * s * s_rel; n = _np.array([-v[1], v[0]]); big = 500.0
    hi_ = _TP([c + n * big, c + n * big + v * big, c - n * big + v * big, c - n * big])
    lo = W.difference(hi_); hi = _taff.translate(W.intersection(hi_), *(v * L))
    strip = W.intersection(_TP([c + n * big - v * 0.01, c + n * big + v * 0.01, c - n * big + v * 0.01, c - n * big - v * 0.01]).buffer(0.2))
    band = _tuu([_taff.translate(strip, *(v * L * k_ / 20)) for k_ in range(21)]) if L > 0 else strip   # L < 0: Schenkel kuerzer
    return _gross(_tuu([lo, hi, band]).buffer(0.01).buffer(-0.01)), (tuple(f0 + v * L), tuple(_np.array(face[1]) + v * L))
def _drehung(face, p, A, K):
    a, b = map(_np.array, face)
    def dist(th):
        c, s = math.cos(th), math.sin(th); Mr = _np.array([[c, -s], [s, c]])
        aa = Mr @ (a - p) + A; bb = Mr @ (b - p) + A; d = bb - aa; n = _np.array([-d[1], d[0]]) / _np.linalg.norm(d)
        return float(n @ (K - aa))
    ths = _np.radians(_np.linspace(-15, 45, 601)); ds = [dist(t_) for t_ in ths]
    for i_ in range(len(ths) - 1):
        if ds[i_] * ds[i_ + 1] <= 0:
            lo, hi = ths[i_], ths[i_ + 1]
            for _k in range(50):
                m_ = (lo + hi) / 2
                if dist(lo) * dist(m_) <= 0: hi = m_
                else: lo = m_
            return lo
    raise SystemExit("Kugellager: keine Hebellage fuer den Knopf gefunden")
def _hebel_2d(r, L):
    p = _np.array(H_PIV[r]); A = _np.array(ACHSE_K[r]); K = _np.array(KN[r])
    W, face = _verlaengern(H_HOLZ[r], H_PIV[r], H_FLAECHE[r], L)
    rot = _drehung(face, p, A, K)
    g = lambda G_: _taff.translate(_taff.rotate(G_, rot, origin=tuple(p), use_radians=True), *(A - p))
    Wr = g(W); fa = _np.array(g(_TPt(face[0])).coords[0]); fb = _np.array(g(_TPt(face[1])).coords[0])
    dist = float(_np.linalg.norm(K - A)); th = _THUB / dist
    d = (fb - fa) / _np.linalg.norm(fb - fa); n = _np.array([-d[1], d[0]])
    tilt = math.degrees(math.atan2(-d[0], d[1]))          # Neigung der Kopfflaeche = Knopfneigung (Hebelform bleibt)
    if r == 3:                                            # 22.11. 06:58 (Johann): Stirnseite R3 kuerzen - stiess gedrueckt an das Kanalende
        _yl = _KT((76.46, 33.0))[1] - 1.5 - 0.3           # Kanalende Teil 2 (1,5-mm-Abschluss) minus 0,3 Luft
        Wr = _gross(Wr.intersection(_taff.rotate(_tbox(-500, -500, 500, _yl), math.degrees(th), origin=tuple(A))))
    return dict(W=Wr, Al=g(H_ALU[r]), KE=(g(_TPt(H_KERBE[r])).coords[0], g(_TPt(H_ENDE[r])).coords[0]), Vt=g(H_VERT[r]), A=tuple(A), K=tuple(K), th=math.degrees(th), dist=dist,
                rot=math.degrees(rot), tilt=tilt, L=L, d=d, n=n)
_FREMD = {r_: (KAMM2 if r_ in (1, 2) else KAMM1).difference(_TRAEGER.buffer(1.0)).intersection(_tbox(-500, max(a_[1] for a_ in ACHSE_K.values()) + 10.0, 500, 500))
          for r_ in (1, 2, 3, 4)}
H2D = {}
def _koll_fremd(h_, r_):
    return max(h_["W"].buffer(0.5).intersection(_FREMD[r_]).area, _rot(h_["W"], -h_["th"], h_["A"]).buffer(0.5).intersection(_FREMD[r_]).area)
for r_ in (1, 2, 3, 4):
    # 21.11. (Johann): Hebelform bleibt wie gezeichnet, nur der duenne Schenkel darf laenger/kuerzer werden
    _kand = []
    for L_ in (0.0,):                                     # 21.11. 18:46 (Johann): auch der Holzschenkel bleibt unveraendert
        h_ = _hebel_2d(r_, float(L_))
        if _koll_fremd(h_, r_) < 0.05: _kand.append((abs(h_["tilt"] - h_["th"]) + 0.15 * abs(L_), float(L_), h_))
    if not _kand: _kand = [(0, 0.0, _hebel_2d(r_, 0.0))]; print(f"   ACHTUNG Reihe {r_}: Kopf beruehrt die Stufe des anderen Teils")
    h_ = min(_kand, key=lambda q_: q_[0])[2]
    H2D[r_] = h_
    print(f"Kugellager Reihe {r_}: Achse {tuple(round(float(v_), 2) for v_ in h_['A'])}, Hebel um {h_['rot']:.1f} Grad gedreht, "
          f"duenner Schenkel {h_['L']:+.1f} mm; Knopf {h_['tilt']:+.2f} Grad geneigt, Hubwinkel {h_['th']:.2f} -> gedrueckt "
          f"{h_['tilt'] - h_['th']:+.2f} Grad gegen die Stufe")

# --- Spiralfeder unter jedem Hebel (wie traditionell, Johann 21.11.): Windung Ø3 aussen, Draht 0,6, 3 Windungen, an der
# Rueckseite des Schenkels (Seite, auf die der Hebel beim Druecken zugeht); ein Schenkel (14 mm) liegt auf der
# Hebelrueckseite, der andere in einer Nut im eigenen Kammteil. Platz dafuer wird im Kamm ausgespart.
F_D, D_DRAHT_T, F_SCHENKEL = 3.0, 0.6, 14.0
_F_RM = (F_D - D_DRAHT_T) / 2
_EIGEN = {r_: (KAMM1 if r_ in (1, 2) else KAMM2).difference(_TRAEGER.buffer(1.0)) for r_ in (1, 2, 3, 4)}
_PARTNER = {1: 3, 3: 1, 2: 4, 4: 2}                      # teilen sich die Luecke zwischen zwei Stuetzen
FED2D = {}
for r_ in (1, 2, 3, 4):
    H_ = H2D[r_]; A_ = _np.array(H_["A"]); K_ = _np.array(H_["K"]); Lh_ = H_["W"]
    t_ = (K_ - A_) / _np.linalg.norm(K_ - A_); nb_ = _np.array([t_[1], -t_[0]])
    if nb_[0] < 0: nb_ = -nb_                             # Rueckseite = +x (Druckrichtung)
    Hp_ = H2D[_PARTNER[r_]]
    _frei = _tuu([Hp_["W"], _rot(Hp_["W"], -Hp_["th"], Hp_["A"])]).buffer(0.5)
    def _kante(s_):
        P_ = A_ + t_ * s_
        g_ = Lh_.intersection(_TL([P_ - nb_ * 12.0, P_ + nb_ * 12.0]))
        if g_.is_empty: return None
        return max((_np.array(q_) - P_) @ nb_ for g2_ in getattr(g_, "geoms", [g_]) for q_ in g2_.coords)
    for s_ in _np.arange(12.0, 60.0, 0.5):
        _es = [_kante(s_ + d_) for d_ in _np.linspace(-F_D / 2, F_SCHENKEL, 12)]
        if any(e_ is None for e_ in _es): continue
        e1_ = e2_ = max(_es)                               # Feder und Schenkel ausserhalb der breitesten Stelle des Hebels
        C_ = A_ + t_ * s_ + nb_ * (e1_ + 0.4 + F_D / 2)
        Lz_ = _TPt(C_).buffer(F_D / 2 + 0.4).union(_TL([C_ + nb_ * _F_RM, C_ + nb_ * _F_RM + t_ * F_SCHENKEL]).buffer(0.6))
        if _EIGEN[r_].buffer(-0.8).contains(Lz_) and not Lz_.intersects(_frei):
            endA_ = A_ + t_ * (s_ + F_SCHENKEL) + nb_ * (e2_ + D_DRAHT_T / 2 + 0.05)
            FED2D[r_] = dict(C=tuple(C_), A0=tuple(C_ - nb_ * _F_RM), A1=tuple(endA_), B0=tuple(C_ + nb_ * _F_RM),
                             B1=tuple(C_ + nb_ * _F_RM + t_ * F_SCHENKEL),
                             cut=_tuu([_TPt(C_).buffer(F_D / 2 + 0.4), _TL([C_ - nb_ * _F_RM, endA_]).buffer(0.5),
                                       _TL([C_ + nb_ * _F_RM, C_ + nb_ * _F_RM + t_ * F_SCHENKEL]).buffer(0.5),
                                       _TL([A_ + t_ * s_ + nb_ * e1_, C_]).buffer(F_D / 2)]).buffer(0), s=s_)
            break
    print(f"Kugellager Reihe {r_}: Spiralfeder " + (f"{FED2D[r_]['s']:.1f} mm von der Achse an der Schenkelrueckseite" if r_ in FED2D else "- KEIN PLATZ gefunden"))
_feder_cut = []

# --- alte Mechanikteile entfernen ---
_weg_n = {"Griffstock_Koerper", "Griffstock_Teil2", "Griffstock_Abdeckung", "Kunststoffteil_Lagerblock", "Dreieckleiste"}
_weg_s = ("/Holzteil", "/Hebel_Alu", "/Messingbuchse", "/Achse", "/Feder")
_weg_p = ("Schraube_Kunststoff_", "Schraube_Boden_Dreieck_")
_vorher = len(teile)
_d_alt_pol = {n_.split("/")[0]: s_ for n_, s_, f_ in teile if n_.startswith("Hebel_R4_") and n_.endswith("/Klappenpolster")}
teile = [t_ for t_ in teile if not (t_[0] in _weg_n or t_[0].endswith(_weg_s) or t_[0].startswith(_weg_p))]
print(f"Kugellager: {_vorher - len(teile)} Teile der Kunststoff-Lagerung entfernt")

# --- Stuetzen-X: mittig zwischen den Hebelgruppen (halbe Knopfteilung, je Instrument aus den Knopf-X) ---
_hx = sorted(x_ + KV[r_] for r_, x_, k_, kr_ in zuordnung)
_gruppen = []
for x_ in _hx:
    if _gruppen and x_ - _gruppen[-1][-1] < 2.5: _gruppen[-1].append(x_)
    else: _gruppen.append([x_])
_gm = [sum(g_) / len(g_) for g_ in _gruppen]
STUETZEN_X = [_gm[0] - (_gm[1] - _gm[0]) / 2] + [(a_ + b_) / 2 for a_, b_ in zip(_gm, _gm[1:])] + [_gm[-1] + (_gm[-1] - _gm[-2]) / 2]
_teil_min = min(b_ - a_ for a_, b_ in zip(STUETZEN_X, STUETZEN_X[1:]))
T_STU = round(min(4.0, _teil_min - (2 * B_LAGER + 0.2) - 0.2), 2)   # 22.11.: bei 8 mm Teilung Stuetze 3,6, damit beide Lager (4 mm) ins Holz passen
_luecke = min(b_ - a_ for a_, b_ in zip(STUETZEN_X, STUETZEN_X[1:])) - T_STU
print(f"Kugellager: {len(STUETZEN_X)} POM-Stuetzen {T_STU} mm, X {STUETZEN_X[0]:.2f} ... {STUETZEN_X[-1]:.2f}, "
      f"Teilung {min(b_ - a_ for a_, b_ in zip(STUETZEN_X, STUETZEN_X[1:])):.2f}-{max(b_ - a_ for a_, b_ in zip(STUETZEN_X, STUETZEN_X[1:])):.2f}, "
      f"kleinste Luecke {_luecke:.2f} mm (Holzhebel {T_HZ})")
T_HZ = round(_luecke - 0.2, 2)                            # 22.11. 06:02 (Johann): Holzhebel so breit wie moeglich, 0,1 Spiel je Seite
print(f"Kugellager: Holzhebel {T_HZ} mm breit (Luecke {_luecke:.2f} mm - 2 x 0,1 Spiel), Lager 2 x 2 mm mittig, Alu {T_ALU} mm buendig in der Vertiefung")
X_TR0, X_TR1 = max(GS_X0 + 1.0, STUETZEN_X[0] - T_STU / 2 - 5.0), min(GS_X1 - 1.0, STUETZEN_X[-1] + T_STU / 2 + 5.0)

# --- Holzhebel, Lager, Achsen, Alu je Taste ---
TRAD_HEBEL = []
_TRAD_2D = {}
_ACHSE = ACHSE_K
_sweep_holz, _sweep_alu = {}, {}
_ALU_INFO = {}
X_KR0, X_KR1 = KAMM_GES.bounds[0] - 0.5, KAMM_GES.bounds[0] - 5.5   # Kroepfzone der Alus 1/2 (Schnitt-x) vor dem Kamm
def _p3(q, X_):
    u_, Y_ = _M(q); return cq.Vector(X_, Y_, -u_)
def _zyl_x(q, r_, x0, x1):
    return cq.Solid.makeCylinder(r_, x1 - x0, _p3(q, x0), cq.Vector(1, 0, 0))
for n_, (r, x, k, kr) in enumerate(zuordnung):
    G = f"Hebel_R{r}_{n_+1:02d}"
    xk = x + KV[r]
    xa = xk + (1.0 if r in (1, 2) else -1.0) * (T_HZ - T_ALU) / 2      # Alu buendig in der 2-mm-Vertiefung, 1/2 und 3/4 gegenueber
    xt = xa                                               # 21.11. 20:11 (Johann): Alu-Schenkel ohne Versatz, gerade durchgehend
                                                          # sonst laufen die Alus 1/3 bzw. 2/4 einer Luecke Flaeche an Flaeche
    H_ = H2D[r]; ax_ = H_["A"]; Lh = H_["W"]
    # Alu: in der Vertiefung, ausserhalb des Holzes nach Zeichnung, dann bis zum Armband bei ya (Schenkellaenge angepasst)
    if kr != r:
        ya = 36.0; pfad = [(_M(ax_)[1], xt), (ya, xt), (46.0, k[0]), (k[1], k[0])]
    else:
        f = FUEHR[r]; v = k[0]; lane = v + f["lane"]
        if r == 1 and f["lane"] != 0.0 and not any(abs(k3[0] - (v + 11.5)) < 4.0 for k3 in kl_r[3]): lane = v + 5.0
        ya, yb = f["y_lane"]
        pfad = [(_M(ax_)[1], xt), (ya, xt), (yb, lane)]
        if f["y_fin"]:
            yf0, yf1 = f["y_fin"]
            if yf0 - yb > 0.5: pfad.append((yf0, lane))
            pfad.append((yf1, v))
        pfad.append((k[1], v))
    bd = arm_band(r, ya - 0.5) or arm_band(r, ya + 1)
    if r == 1 and ya < 40: bd = arm_band(1, 36.5)
    _p0, _p1 = _Minv((bd[0], ya + 0.5)), _Minv((bd[1], ya + 0.5))
    _xb = min(_p0[0], _p1[0])
    # 22.11. 05:47 (Johann): Alu-Schenkel ganz nach Zeichnung in einem Stueck: Nabe bis kurz nach der Kerbe wie gezeichnet,
    # gerades Band bis zur Endform, Endform der Zeichnung auf der Klappe (Unterseite auf der Klappenoberseite, Spitze bei
    # Klappe + 4 mm); passt der Winkel nicht, wird er ab der Kerbe korrigiert - das Band laeuft gerade, ohne Absatz
    Al_ = H_["Al"]; Kp_, Ep_ = map(_np.array, H_["KE"])
    dKE = (Ep_ - Kp_) / _np.linalg.norm(Ep_ - Kp_); nKE = _np.array([-dKE[1], dKE[0]])
    _halb = lambda P_, sg_: _TP([P_ + nKE * 60, P_ + nKE * 60 + sg_ * dKE * 400, P_ - nKE * 60 + sg_ * dKE * 400, P_ - nKE * 60])
    nabe_ = _gross(Al_.intersection(_halb(Kp_, -1.0)))
    ende0_ = _gross(Al_.intersection(_halb(Ep_, +1.0)))
    _quer = lambda P_, G_: _gross(G_.buffer(0.01).intersection(_TL([P_ - nKE * 30, P_ + nKE * 30]).buffer(0.02)))
    cK_ = _quer(Kp_, Al_); cE0_ = _quer(Ep_, Al_)
    U_kl = U_KL4 if r == 4 else U_KL
    ytip = k[1] + 4.0
    _a0, _a1 = _np.array(_Minv((U_kl, ytip - 10.0))), _np.array(_Minv((U_kl, ytip)))
    dkl = (_a1 - _a0) / _np.linalg.norm(_a1 - _a0)       # Klappenoberseite im Schnitt, Richtung Klappenspitze
    # 22.11. 06:02 (Johann): Endform vom oberen Knickpunkt weg so drehen, dass ihre Unterkante flach auf der Klappe liegt
    _gu = _np.array([_AFF[0, 0], _AFF[1, 0]])              # Richtung zunehmender Hoehe u im Schnitt
    _co = list(ende0_.exterior.coords); _ccw = ende0_.exterior.is_ccw; _uk = None
    for a_, b_ in zip(_co, _co[1:]):
        a_, b_ = _np.array(a_), _np.array(b_); L_ = float(_np.linalg.norm(b_ - a_))
        if L_ < 2.0: continue
        d_ = (b_ - a_) / L_; n_ = _np.array([d_[1], -d_[0]]) if _ccw else _np.array([-d_[1], d_[0]])
        if n_ @ _gu < -0.7 * _np.linalg.norm(_gu) and (_uk is None or L_ > _uk[0]): _uk = (L_, d_)
    _wk = math.atan2(dkl[1], dkl[0])
    _w1 = [_wk - math.atan2(sg_ * _uk[1][1], sg_ * _uk[1][0]) for sg_ in (1.0, -1.0)]
    beta_ = math.degrees(min(((w_ + math.pi) % (2 * math.pi) - math.pi for w_ in _w1), key=abs))
    ende_ = _rot(ende0_, beta_, tuple(Ep_)); cE_ = _rot(cE0_, beta_, tuple(Ep_))
    _uy = [_M(q_) for q_ in ende_.exterior.coords]
    dU_, dY_ = U_kl - min(u_ for u_, Y_ in _uy), ytip - max(Y_ for u_, Y_ in _uy)
    _o = _np.array(_Minv((0.0, 0.0))); _t = _np.array(_Minv((dU_, dY_))) - _o
    ende_ = _taff.translate(ende_, *_t); cE_ = _taff.translate(cE_, *_t)
    band_ = _tuu([cK_, cE_]).convex_hull
    Afull = _gross(_tuu([nabe_, band_, ende_]).buffer(0.02).buffer(-0.02))
    Afull = _gross(Afull.intersection(_Ub(0.4)))
    # Profil (u, Y) seitlich entlang des Pfades geschert (wie bisher): bis ya auf xa, dann zu Spur und Klappe
    prof_ = _Mpg(Afull)
    _Ymin, _Ymax = prof_.bounds[1], prof_.bounds[3]
    pts_ = [(_Ymin - 1.0, xa)] + [p_ for p_ in pfad if p_[0] > _Ymin + 0.5]
    if pts_[-1][0] < _Ymax + 1.0: pts_.append((_Ymax + 1.0, pts_[-1][1]))
    for i_ in range(1, len(pts_)):
        if abs(pts_[i_][1] - pts_[i_ - 1][1]) < 0.05: pts_[i_] = (pts_[i_][0], pts_[i_ - 1][1])
    _m = [pts_[0]]
    for q_ in pts_[1:]:
        if q_[0] - _m[-1][0] < 0.05: continue
        if len(_m) >= 2:
            (y1_, x1_), (y2_, x2_) = _m[-2], _m[-1]
            if abs((x2_ - x1_) / (y2_ - y1_) - (q_[1] - x2_) / (q_[0] - y2_)) < 1e-6:
                _m[-1] = q_; continue
        _m.append(q_)
    alu = None
    for (ya_, xa_), (yb_, xb_) in zip(_m, _m[1:]):
        st = prof_.intersection(_tbox(-500, ya_ - 0.01, 500, yb_ + 0.01))
        for st1 in _tteile(st):
            sol = prisma_uy(_pp(st1), xa_, T_ALU)
            kk = (xb_ - xa_) / (yb_ - ya_)
            if abs(kk) > 1e-6:
                sol = sol.transformGeometry(cq.Matrix([[1.0, kk, 0.0, -kk * ya_], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]))
            alu = sol if alu is None else alu.fuse(sol)
    Ca = _gross(Afull.intersection(_tbox(KAMM_GES.bounds[0] - 40.0, -500, 500, 500)))   # Griffstockbereich (Schlitze, Verdeck)
    if r not in _ALU_INFO:
        _ALU_INFO[r] = (beta_, Afull)
        _wb = math.degrees(math.atan2(*(_np.array(band_.centroid.coords[0]) - Kp_)[::-1]) - math.atan2(dKE[1], dKE[0]))
        print(f"Kugellager Reihe {r}: Endform um {beta_:+.1f} Grad gedreht (Unterkante flach auf der Klappe, u {U_kl:.1f}, Spitze Y {ytip:.1f}), "
              f"Band ab der Kerbe um ca. {(_wb + 180) % 360 - 180:+.1f} Grad korrigiert")

    alu = alu.cut(_zyl_x(ax_, D_LAGER / 2, xa - T_ALU, xa + T_ALU))
    holz = prisma_uy(_pp(_Mpg(Lh)), xk, T_HZ)
    holz = holz.cut(prisma_uy(_pp(_Mpg(nabe_.buffer(0.05))), xa, T_ALU + 0.1))   # 2-mm-Vertiefung fuer das Alu
    holz = holz.cut(_zyl_x(ax_, D_LAGER / 2, xk - T_HZ, xk + T_HZ))           # Lagersitz Ø5
    teile.append((f"{G}/Holzhebel", holz, f"Holz_R{r}"))
    teile.append((f"{G}/Hebel_Alu", alu, f"Alu_R{r}"))
    for i_, xl_ in enumerate((xk - B_LAGER / 2, xk + B_LAGER / 2)):          # zwei Kugellager 5x2x2 nebeneinander
        lg_ = _zyl_x(ax_, D_LAGER / 2, xl_ - B_LAGER / 2, xl_ + B_LAGER / 2).cut(_zyl_x(ax_, D_ACHSE_K / 2, xl_ - B_LAGER, xl_ + B_LAGER))
        teile.append((f"{G}/Kugellager{i_+1}", lg_, "Stahl"))
    _sl = max([s_ for s_ in STUETZEN_X if s_ < xk], default=xk - 5.0); _sr = min([s_ for s_ in STUETZEN_X if s_ > xk], default=xk + 5.0)
    teile.append((f"{G}/Achse", _zyl_x(ax_, D_ACHSE_K / 2, _sl, _sr), "Stahl"))   # kurze Achse von Stuetzenmitte zu Stuetzenmitte
    DREHPUNKTE[G] = (xk, _M(ax_)[1], -_M(ax_)[0])
    if r in FED2D:                                         # Spiralfeder: 3 Windungen (Achse = X), zwei Schenkel
        F_ = FED2D[r]; _feder = None
        for j_ in range(3):
            _ring = cq.Solid.makeTorus(_F_RM, D_DRAHT_T / 2, _p3(F_["C"], xk - 0.7 + 0.7 * j_), cq.Vector(1, 0, 0))
            _feder = _ring if _feder is None else _feder.fuse(_ring)
        for (q0_, q1_), X_ in (((F_["A0"], F_["A1"]), xk + 0.7), ((F_["B0"], F_["B1"]), xk - 0.7)):
            a3_, b3_ = _p3(q0_, X_), _p3(q1_, X_); v_ = b3_.sub(a3_)
            _feder = _feder.fuse(cq.Solid.makeCylinder(D_DRAHT_T / 2, v_.Length, a3_, v_.normalized()))
        teile.append((f"{G}/Feder", _feder, "Stahl"))
        _feder_cut.append((r, prisma_uy(_pp(_Mpg(F_["cut"])), xk, T_HZ + 0.6)))
    th = H_["th"]
    _sh = _tuu([_rot(Lh, -th * i_ / 6.0, ax_) for i_ in range(7)]).buffer(0.3)
    _sa = _tuu([_rot(Ca, -th * i_ / 6.0, ax_) for i_ in range(7)]).buffer(0.3)
    _sweep_holz.setdefault(r, []).append((xk, _sh)); _sweep_alu.setdefault(r, []).append((xa, _sa))   # je Reihe gesammelt (innerhalb einer Reihe ueberlappen die Werkzeuge nicht)
    TRAD_HEBEL.append((G, r, xk, round(th, 2), round(H_["dist"], 1)))
    if kr == r and r not in _TRAD_2D: _TRAD_2D[r] = dict(L=Lh, Ca=Ca, th=th, ax=ax_, knopf=H_["K"], G=G)
try:                                                      # Vorschau Alu ganz (je Reihe der erste Hebel)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as _plt2
    _f2, _a2 = _plt2.subplots(figsize=(16, 9))
    for g_ in _tteile(KAMM1) + _tteile(KAMM2): _a2.fill(*g_.exterior.xy, color="#d9b48f", alpha=0.4)
    _cl2 = {1: "tab:red", 2: "tab:blue", 3: "tab:green", 4: "tab:orange"}
    for r_, (b_, A_) in sorted(_ALU_INFO.items()):
        _a2.plot(*A_.exterior.xy, color=_cl2[r_], lw=1.2, label=f"Alu R{r_} (Band {b_:+.1f} Grad)")
        _a2.plot(*H2D[r_]["W"].exterior.xy, color=_cl2[r_], lw=0.8, ls=":")
    _ub = [_Minv((0.0, Y_)) for Y_ in (-10, 200)]; _a2.plot([q_[0] for q_ in _ub], [q_[1] for q_ in _ub], "k", lw=0.8, label="Deckelflaeche")
    _a2.set_aspect("equal"); _a2.grid(alpha=0.3); _a2.legend(fontsize=8)
    _f2.savefig(OUT.replace(".step", "_kugel_alu.png"), dpi=70, bbox_inches="tight")
except Exception as _ex2:
    print("Alu-Vorschau:", _ex2)
# Knoepfe um den Hubwinkel neigen (Drehachse = X durch die Mitte der Knopfrueckseite)
def _neig_winkel(K, th):
    a_, b_ = _np.array(_M(K)), _np.array(_M((K[0], K[1] + 10.0)))
    c_ = _np.array(_M((K[0] - 10.0 * math.sin(math.radians(th)), K[1] + 10.0 * math.cos(math.radians(th)))))
    v_ = _np.array([b_[1] - a_[1], -(b_[0] - a_[0])]); w_ = _np.array([c_[1] - a_[1], -(c_[0] - a_[0])])
    return math.degrees(math.atan2(v_[0] * w_[1] - v_[1] * w_[0], v_ @ w_))
_neu = []
for n_t, s_t, f_t in teile:
    if n_t.startswith("Hebel_R") and n_t.endswith("/Knopf"):
        r_ = int(n_t[7]); K_ = H2D[r_]["K"]; X_ = DREHPUNKTE[n_t.split("/")[0]][0]
        P_ = _p3(K_, X_); s_t = s_t.rotate(P_, P_.add(cq.Vector(1, 0, 0)), _neig_winkel(K_, H2D[r_]["tilt"]))
    _neu.append((n_t, s_t, f_t))
teile = _neu
print("Kugellager: Knoepfe geneigt wie die Kopfflaeche (Grad) je Reihe:", {r_: round(H2D[r_]["tilt"], 2) for r_ in H2D})
print(f"Kugellager: {len(TRAD_HEBEL)} Holzhebel mit je 2 Kugellagern, Hub-Winkel {min(t_[3] for t_ in TRAD_HEBEL)}-{max(t_[3] for t_ in TRAD_HEBEL)} Grad")

# Klappen der Reihe 4 flacher (wie traditionell): Decke 0,6, Polster Filz 2 + Leder 0,6
_neu_t = []
for n_t, s_t, f_t in teile:
    if n_t.startswith("Hebel_R4_") and n_t.endswith(("/Klappe", "/Klappenpolster")):
        _G = n_t.split("/")[0]; _po = _d_alt_pol.get(_G)
        if _po is None: continue
        _sl = lambda z0_, z1_: cq.Solid.makeBox(1e4, 1e4, z1_ - z0_, cq.Vector(-5e3, -5e3, z0_))
        if n_t.endswith("/Klappenpolster"): s_t = _po.intersect(_sl(-T_POLSTER4, 0.0))
        else: s_t = _po.intersect(_sl(-T_KLAPPE4, 0.0)).translate(cq.Vector(0, 0, -T_POLSTER4))
    _neu_t.append((n_t, s_t, f_t))
teile = _neu_t

# --- Traeger: POM-Stuetzen und FR2-Leisten ---
_stu_pp = _pp(_Mpg(STUETZE))
_stuetzen = []
for i_, xs_ in enumerate(STUETZEN_X):
    s_ = prisma_uy(_stu_pp, xs_, T_STU)
    for a_ in ACHSE_K.values(): s_ = s_.cut(_zyl_x(a_, 1.025, xs_ - T_STU, xs_ + T_STU))
    _stuetzen.append(s_); teile.append((f"Stuetze_POM_{i_+1:02d}", s_, "POM"))
_schlitz = [prisma_uy(_stu_pp, xs_, T_STU + 0.02) for xs_ in STUETZEN_X]
# 22.11. 06:58 (Johann): untere Epoxy-Leiste mit dem Diskantboden verschraubt - Bohrungen Ø3,5 (Leiste, Stelle wie
# Zeichnung 18,18 vom Zungenende) und Ø2 Vorbohrung 5 mm tief im Boden, je in einer Stuetzenluecke
_lo = _np.array(LEISTE1_LOCH) + _np.array(_dv); _nl = _RK @ _nd
if _M(tuple(_lo + _nl))[0] > _M(tuple(_lo))[0]: _nl = -_nl      # _nl zeigt zum Boden
_lo_o = _lo - _nl * (T_L1 + 1.0); _lu = _lo + _nl * 1.0         # ueber der Leiste / Leistenunterseite (Boden)
_lu = _lo_o + _nl * (max(0.0, *( (_np.array(q_) - _lo_o) @ _nl for q_ in LEISTE1.exterior.coords)))
assert LEISTE1.contains(_TPt(tuple(_lo))), "Leistenbohrung liegt nicht in der Leiste"
_gaps = [(a_ + b_) / 2 for a_, b_ in zip(STUETZEN_X, STUETZEN_X[1:])]
LOCH_X = sorted(set([_gaps[0], _gaps[-1]] + _gaps[2:-2:4]))
def _zyl_r(q0, q1, r_, X_):
    a_, b_ = _p3(tuple(q0), X_), _p3(tuple(q1), X_); v_ = b_.sub(a_)
    return cq.Solid.makeCylinder(r_, v_.Length, a_, v_.normalized())
_loch_l = [_zyl_r(_lo_o, _lu + _nl * 0.2, 1.75, X_) for X_ in LOCH_X]
_loch_b = [_zyl_r(_lu - _nl * 0.2, _lu + _nl * 5.0, 1.0, X_) for X_ in LOCH_X]
print(f"Kugellager: untere Leiste {len(LOCH_X)} Schraubloecher Ø3,5 (Boden Ø2 x 5 vorgebohrt) bei X " + ", ".join(f"{x_:.1f}" for x_ in LOCH_X))
for nm_, pg_ in (("Leiste_FR2_34", LEISTE1), ("Leiste_FR2_32", LEISTE2)):
    l_ = schneide(prisma_uy(_pp(_Mpg(pg_)), (X_TR0 + X_TR1) / 2, X_TR1 - X_TR0), _schlitz)
    if nm_ == "Leiste_FR2_34": l_ = schneide(l_, _loch_l)
    teile.append((nm_, l_, "FR2"))
# Tasche im Diskantboden: untere Leiste ganz, Stuetzenfuesse an ihren X (0,05 Luft)
_tb_l = _gross(LEISTE1.buffer(0.05).intersection(_Ub(-50.0, 0.05)))
_tb_s = STUETZE.buffer(0.05).intersection(_Ub(-50.0, 0.05))
_boden_t = [prisma_uy(_pp(_Mpg(_tb_l)), (X_TR0 + X_TR1) / 2, X_TR1 - X_TR0 + 0.1)]
_boden_s = [prisma_uy(_pp(_Mpg(g_)), xs_, T_STU + 0.1) for xs_ in STUETZEN_X for g_ in _tteile(_tb_s)]
teile = [(n_, schneide(schneide(schneide(s_, _boden_t), _boden_s), _loch_b) if n_ == "Diskantboden" else s_, f_) for n_, s_, f_ in teile]
_ul = min(_M(q_)[0] for q_ in LEISTE1.exterior.coords)
print(f"Kugellager: Leisten FR2 X {X_TR0:.1f}-{X_TR1:.1f} ({X_TR1 - X_TR0:.0f} mm, {len(STUETZEN_X)} Schlitze), Dicke unten {T_L1:.2f} / "
      f"seitlich {T_L2:.2f} mm (Nasen 0,1 durch); untere Leiste im Diskantboden eingelassen bis u {_ul:.2f} (Boden {T_BODEN} mm)")

# --- Holzkaemme Teil 1 / Teil 2 ---
# 22.11. 07:50 (Johann): Daumenleiste wie bei den anderen Varianten - Mulde (R17, 1,5 tief) in der Stirnseite, Kanten R1,5
_st2 = sorted(KAMM1.exterior.coords[:-1], key=lambda q_: -q_[1])[:2]
_smk = _np.array([(_st2[0][0] + _st2[1][0]) / 2, (_st2[0][1] + _st2[1][1]) / 2])
_sdk = _np.array(_st2[0]) - _np.array(_st2[1]); _sdk = _sdk / _np.linalg.norm(_sdk)
_snk = _np.array([-_sdk[1], _sdk[0]]); _snk = _snk if _snk[1] > 0 else -_snk          # nach aussen (Daumen)
_mulde_k = _TPt(tuple(_smk + _snk * (17.0 - 1.5))).buffer(17.0, 64)
_y_rund = _smk[1] - 12.0
_k1o = KAMM1.intersection(_tbox(-500, _y_rund, 500, 500)).difference(_mulde_k).buffer(-1.5, join_style=1).buffer(1.5, join_style=1)
KAMM1 = _gross(_tuu([KAMM1.intersection(_tbox(-500, -500, 500, _y_rund + 3.0)).difference(_mulde_k), _k1o]).buffer(0.01).buffer(-0.01).simplify(0.02))
print(f"Kugellager: Daumenleiste mit Mulde R17 x 1,5 und Kanten R1,5 (Stirnseite {math.dist(*_st2):.1f} mm breit)")
_T1uy, _T2uy = _Mpg(KAMM1), _Mpg(KAMM2)
HUELLE_G = _gr(_tuu([_T1uy, _T2uy]).buffer(0.01).buffer(-0.01))
_LGS = GS_X1 - GS_X0
t1 = schneide(prisma_uy(_pp(_T1uy), (GS_X0 + GS_X1) / 2, _LGS), seitenschnitt)
t2 = schneide(prisma_uy(_pp(_T2uy), (GS_X0 + GS_X1) / 2, _LGS), seitenschnitt)
_v1, _v2 = t1.Volume(), t2.Volume()
# 21.11. 20:11 (Johann): im Bereich des Kunststoffkamms (Stuetzen) beide Holzteile ganz ausgeschnitten bis ueber die Stuetzen
Y_TR_OBEN = STUETZE.bounds[3] + 0.5
X_GRUND1 = _KT((50.76, 50.0))[0] - 1.5                    # Teil 1 hinten: 1,5-mm-Wand bis auf den Diskantboden (auch an der Leiste)
# 22.11. 09:02 (Johann): Rechteckausschnitt bleibt in voller Breite (keine duennen Kaemme bis zum Boden), nur die
# schmale Oberseite deckt sich mit den POM-Teilen: Teil 2 gerade bis zur Oberkante seines Stuetzenteils, Teil 1 folgt
# der Stuetzenstufe samt Laengsnut; die Stufe liegt genau an der Trennfuge Teil 1/Teil 2 (keine Reste am Uebergang)
from shapely.geometry import LineString as _LS
_TR_OBEN = _tuu([STUETZE, LEISTE2]).buffer(0.8, join_style=2).buffer(-0.8, join_style=2)   # Eckfreibohrungen zu
_xb = max(c_[0] for g_ in getattr(KAMM2.intersection(_LS([(-500, 57.0), (500, 57.0)])), "geoms",
                                  [KAMM2.intersection(_LS([(-500, 57.0), (500, 57.0)]))]) for c_ in g_.coords)   # Trennfuge
_tb = _TR_OBEN.bounds; _dach = []
for x_ in _np.arange(_tb[0] + 0.01, _tb[2] - 0.005, 0.05):
    s_ = _TR_OBEN.intersection(_LS([(x_, -500), (x_, 500)]))
    if not s_.is_empty: _dach.append((float(x_), s_.bounds[3]))
Y_AUS2 = max(y_ for x_, y_ in _dach if x_ < _xb - 0.5) + 0.1          # Teil 2: gerade Oberseite
_d1 = [(x_, y_) for x_, y_ in _dach if x_ >= _xb + 0.5]
Y_AUS1 = _d1[0][1]
_DECKE1 = _TP([(-500, -500), (-500, Y_AUS1)] + [(_xb, Y_AUS1)] + _d1 + [(500, _d1[-1][1]), (500, -500)]).buffer(0).buffer(0.1, join_style=2)
_rahmen = KAMM_GES.buffer(1.0).intersection(_tbox(-500, -500, X_GRUND1, 500))
_aus1 = _gross(_rahmen.intersection(_DECKE1)); _aus2 = _gross(_rahmen.intersection(_tbox(-500, -500, 500, Y_AUS2)))
print(f"Kugellager: Rechteckausschnitt volle Breite, Oberseite an den POM-Teilen: Teil 2 bis y {Y_AUS2:.1f}, Teil 1 ab der Fuge "
      f"(x {_xb:.2f}) {Y_AUS1 + 0.1:.1f} mit Laengsnut bis {max(y_ for x_, y_ in _d1) + 0.1:.1f} (vorher {Y_TR_OBEN:.1f})")
_L_ = X_TR1 - X_TR0 + 0.2
t1 = t1.cut(prisma_uy(_pp(_Mpg(_aus1)), (X_TR0 + X_TR1) / 2, _L_)); t2 = t2.cut(prisma_uy(_pp(_Mpg(_aus2)), (X_TR0 + X_TR1) / 2, _L_))
# Nuten wie traditionell: zur Knopfseite offen, nach unten in den Ausschnitt offen; Teil 1 mit 1,5-mm-Kammergrund hinten,
# zur Daumenleiste geschlossen (Reihe 1 lang, Reihe 2 bis 1,5 mm unter ihre Stufe); Teil 2 ohne Grund durch bis Teil 1,
# oben 1,5 mm Abschluss (Reihe 3 lang, Reihe 4 kurz)
_yo = lambda p_: _KT(p_)[1]
Y_DAUMEN = _yo((61.36, -1.16)); Y_ST2 = _yo((73.96, 17.0)); Y_T2 = _yo((76.46, 33.0)); Y_ST4 = _yo((81.96, 49.0))
X_GRUND2 = _KT((68.96, 60.0))[0] - 1.5                    # 21.11. 22:29: Teil 2 ist ein Kamm mit 1,5-mm-Nutboden zu Teil 1 hin
NUT_B, NUT_T = 0.8, 0.5                                   # mittige Vertiefung in jedem Nutboden (wie traditionell, 0,5 tief)
_KAMMER = {1: KAMM1.intersection(_tbox(-500, -500, X_GRUND1, Y_DAUMEN - 1.5)),
           2: KAMM1.intersection(_tbox(-500, -500, X_GRUND1, Y_ST2 - 1.5)),
           3: KAMM2.buffer(0.5).intersection(_tbox(-500, -500, X_GRUND2, Y_T2 - 1.5)),
           4: KAMM2.buffer(0.5).intersection(_tbox(-500, -500, X_GRUND2, Y_ST4 - 1.5))}
_GRUNDX = {1: X_GRUND1, 2: X_GRUND1, 3: X_GRUND2, 4: X_GRUND2}
_kam_w, _kam_v = {}, {}
for n_, (r, x, k, kr) in enumerate(zuordnung):
    for g_ in _tteile(_KAMMER[r]): _kam_w.setdefault(r, []).append(prisma_uy(_pp(_Mpg(g_)), x + KV[r], T_HZ + 0.6))
    _b_ = _KAMMER[r].bounds                                # mittige Vertiefung 0,8 x 0,5 im Nutboden
    _yv0 = (max(y_ for x_, y_ in _d1) + 0.1) if r in (1, 2) else Y_AUS2      # 22.11. 09:23: nicht im Epoxy-/Traegerausschnitt
    _v_ = _tbox(_GRUNDX[r] - 0.01, max(_b_[1], _yv0), _GRUNDX[r] + NUT_T, _b_[3])
    _kam_v.setdefault(r, []).append(prisma_uy(_pp(_Mpg(_v_)), x + KV[r], NUT_B))
for r_ in (1, 2, 3, 4):
    if _kam_w.get(r_): t1 = schneide(t1, _kam_w[r_]); t2 = schneide(t2, _kam_w[r_])
for r_ in (1, 2, 3, 4):
    if _kam_v.get(r_): t1 = schneide(t1, _kam_v[r_]); t2 = schneide(t2, _kam_v[r_])
print(f"Kugellager: Nuten bis y Reihe 1 {Y_DAUMEN - 1.5:.1f}, 2 {Y_ST2 - 1.5:.1f}, 3 {Y_T2 - 1.5:.1f}, 4 {Y_ST4 - 1.5:.1f}; "
      f"Nutboden 1,5 mm Teil 1 bei x {X_GRUND1:.1f}, Teil 2 bei x {X_GRUND2:.1f}, Vertiefung {NUT_B}x{NUT_T}; Ausschnitt Traegerbereich bis y {Y_TR_OBEN:.1f}")
_vt1, _vt2 = t1.Volume(), t2.Volume()
# Hebelschlitze und Federplatz: je Reihe und Art ein eigener Schnitt - ueberlappende Werkzeuge in einem Compound
# verfaelschen den Schnitt (Teil 2 war sonst fast ganz weg)
_gruppen_w = []
for r_ in (1, 2, 3, 4):
    _gruppen_w.append([prisma_uy(_pp(_Mpg(h_)), xk_, T_HZ + 0.6) for xk_, g_ in _sweep_holz.get(r_, []) for h_ in _tteile(g_)])
    _gruppen_w.append([prisma_uy(_pp(_Mpg(h_)), xk_, T_ALU + 0.6) for xk_, g_ in _sweep_alu.get(r_, []) for h_ in _tteile(g_)])
    _gruppen_w.append([w_ for rr_, w_ in _feder_cut if rr_ == r_])
for gw_ in _gruppen_w:
    if gw_: t1 = schneide(t1, gw_); t2 = schneide(t2, gw_)
# Hebelform bleibt: wo ein Hebel (Ruhe bis gedrueckt, 0,3 Luft) in den Diskantboden reicht, bekommt der Boden eine Tasche;
# das Verdeck wird in diskant.py ebenso freigeschnitten (KUGEL_FREI)
KUGEL_FREI = [gw_ for gw_ in _gruppen_w if gw_]
# Schwenkbereiche aller Hebel (Holz und Alu, Ruhe bis gedrueckt, 0,3 Luft) im Modellschnitt (u, Y) - diskant.py setzt
# damit das Verdeck hinten so weit hoeher, dass kein Hebel anstoesst
KUGEL_HEBEL_UY = [_Mpg(g_) for d_ in (_sweep_holz, _sweep_alu) for r_ in d_ for x_, g0_ in d_[r_] for g_ in _tteile(g0_)]
_vb = [s_.Volume() for n_, s_, f_ in teile if n_ == "Diskantboden"]
_nb = []
for n_, s_, f_ in teile:
    if n_ == "Diskantboden":
        for gw_ in KUGEL_FREI: s_ = schneide(s_, gw_)
    _nb.append((n_, s_, f_))
teile = _nb
print(f"Kugellager: Diskantboden fuer die Hebel freigeschnitten: {_vb[0] - [s_.Volume() for n_, s_, f_ in teile if n_ == 'Diskantboden'][0]:.1f} mm3" if _vb else "")
print(f"Kugellager: Traegertasche nimmt Teil 1 {_v1 - _vt1:.0f}, Teil 2 {_v2 - _vt2:.0f} mm3; Hebelschlitze Teil 1 {_vt1 - t1.Volume():.0f}, "
      f"Teil 2 {_vt2 - t2.Volume():.0f} mm3")
if GS_LOECHER:
    _gsl = [w_ for x_, y_, n_, w_ in GS_LOECHER]
    t1 = schneide(t1, _gsl); t2 = schneide(t2, _gsl)
# 22.11. 07:03 (Johann): Teil 1 und Teil 2 von hinten mit Senkkopfschrauben an vier Stellen verschraubt - die zwei
# mittleren durch die Stege (ueber einer Stuetze) mit 2,4-mm-Schrauben, die zwei aeusseren mit 3,5-mm-Schrauben
_nr = _np.array([-_v_rw[1], _v_rw[0]])
if _nr[0] > 0: _nr = -_nr                                   # von der Rueckwand (Teil 1) Richtung Teil 2
Y_SCHR = 71.0
_q0 = _np.array([KAMM_GES.bounds[2] + 5.0, Y_SCHR])
_sl = _TL([tuple(_q0), tuple(_q0 + _nr * 80.0)])
_e1 = _sl.intersection(KAMM1); _e2 = _sl.intersection(KAMM2)
_ein = _q0 + _nr * min(float((_np.array(c_) - _q0) @ _nr) for g_ in getattr(_e1, "geoms", [_e1]) for c_ in g_.coords)
_trenn = _q0 + _nr * min(float((_np.array(c_) - _q0) @ _nr) for g_ in getattr(_e2, "geoms", [_e2]) for c_ in g_.coords)
_t2end = _q0 + _nr * max(float((_np.array(c_) - _q0) @ _nr) for g_ in getattr(_e2, "geoms", [_e2]) for c_ in g_.coords)
_sx = sorted(STUETZEN_X, key=lambda x_: abs(x_ - (GS_X0 + GS_X1) / 3))[0], sorted(STUETZEN_X, key=lambda x_: abs(x_ - 2 * (GS_X0 + GS_X1) / 3))[0]
SCHRAUB_T12 = [((GS_X0 + X_TR0) / 2, 3.5), (_sx[0], 2.4), (_sx[1], 2.4), ((X_TR1 + GS_X1) / 2, 3.5)]
_schr = []
for X_, d_ in SCHRAUB_T12:
    _tief = min(12.0, float((_t2end - _trenn) @ _nr) - 2.0)
    _schr.append(_zyl_r(_ein - _nr * 1.0, _trenn + _nr * 0.1, (d_ + 0.2) / 2, X_))           # Durchgang Teil 1
    _schr.append(_zyl_r(_trenn - _nr * 0.1, _trenn + _nr * _tief, 0.35 * d_, X_))              # Vorbohrung Teil 2
    _a3 = _p3(tuple(_ein - _nr * 0.01), X_); _b3 = _p3(tuple(_ein + _nr * d_ / 2), X_)
    _schr.append(cq.Solid.makeCone(d_, (d_ + 0.2) / 2, (_b3.sub(_a3)).Length, _a3, (_b3.sub(_a3)).normalized()))  # Senkung 90 Grad
for w_ in _schr: t1 = schneide(t1, [w_]); t2 = schneide(t2, [w_])   # einzeln (ueberlappende Werkzeuge)
# 22.11. 09:24 (Johann): die Schrauben selbst zeichnen - Senkkopf 90 Grad buendig in der Rueckwand, 3,5 x 30 und 2,4 x 25
for i_, (X_, d_) in enumerate(SCHRAUB_T12):
    L_ = 30.0 if d_ >= 3.0 else 25.0
    _a3 = _p3(tuple(_ein), X_); _dv3 = (_p3(tuple(_ein + _nr), X_)).sub(_a3).normalized()
    _k3 = cq.Solid.makeCone(d_, d_ / 2, d_ / 2, _a3, _dv3)
    _s3 = cq.Solid.makeCylinder(d_ / 2, L_ - d_ / 2, _a3.add(_dv3.multiply(d_ / 2)), _dv3)
    teile.append((f"Schraube_Teil12_{i_+1}_{d_}x{L_:.0f}", _k3.fuse(_s3).clean(), "Schraube"))
print("Kugellager: Teil 1/Teil 2 verschraubt (Senkkopf von hinten, y %.0f): " % Y_SCHR
      + ", ".join(f"X {x_:.1f} Ø{d_}" for x_, d_ in SCHRAUB_T12) + f"; Teil 1 {float((_trenn - _ein) @ _nr):.1f} mm, Vorbohrung Teil 2 {_tief:.1f} mm")
teile.append(("Griffstock_Teil2", t2, "Griffstock"))
teile.append(("Griffstock_Koerper", t1, "Griffstock"))
print(f"Kugellager: Griffstock Teil 1 {t1.Volume():.0f} mm3, Teil 2 {t2.Volume():.0f} mm3")

# --- Schnitt als DXF/PNG (je Reihe der erste Hebel) ---
try:
    _dout = ezdxf.new(); _dm = _dout.modelspace()
    _farb = {"Teil1": 1, "Teil2": 3, "Holzhebel": 30, "Holzhebel_gedrueckt": 8, "Alu": 5, "Knopf": 2, "Stuetze": 6, "Leiste": 4}
    def _dpl(pg, lay):
        if lay not in _dout.layers: _dout.layers.add(lay, color=_farb.get(lay.split("_R")[0], 7))
        for g_ in _tteile(pg):
            _dm.add_lwpolyline(list(g_.exterior.coords), close=True, dxfattribs={"layer": lay})
            for h_ in g_.interiors: _dm.add_lwpolyline(list(h_.coords), close=True, dxfattribs={"layer": lay})
    def _knopf_pg(r_):
        H_ = H2D[r_]; K_ = _np.array(H_["K"]); d_, nn_ = H_["d"], H_["n"]
        return _TP([K_ - 7.75 * d_, K_ + 7.75 * d_, K_ + 7.75 * d_ + 3 * nn_, K_ - 7.75 * d_ + 3 * nn_])
    _dpl(KAMM1, "Teil1"); _dpl(KAMM2, "Teil2"); _dpl(STUETZE, "Stuetze"); _dpl(LEISTE1, "Leiste"); _dpl(LEISTE2, "Leiste")
    for r_, D_ in sorted(_TRAD_2D.items()):
        _dpl(D_["L"], f"Holzhebel_R{r_}"); _dpl(_rot(D_["L"], -D_["th"], D_["ax"]), f"Holzhebel_gedrueckt_R{r_}")
        _dpl(D_["Ca"].intersection(_tbox(-200, -100, 100, 300)), f"Alu_R{r_}"); _dpl(_knopf_pg(r_), f"Knopf_R{r_}")
        if "Achse" not in _dout.layers: _dout.layers.add("Achse", color=7)
        _dm.add_circle(D_["ax"], D_LAGER / 2, dxfattribs={"layer": "Achse"}); _dm.add_circle(D_["ax"], D_ACHSE_K / 2, dxfattribs={"layer": "Achse"})
    _dn = OUT.replace(".step", "_kugel_schnitt.dxf"); _dout.saveas(_dn)
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as _plt
    _fg, _ax = _plt.subplots(figsize=(9, 15))
    for pg_, c_, a_ in ((KAMM1, "#d9b48f", 0.5), (KAMM2, "#a7c7a0", 0.6), (STUETZE, "#cc88cc", 0.5), (LEISTE1, "#4682b4", 0.8), (LEISTE2, "#4682b4", 0.8)):
        for g_ in _tteile(pg_): _ax.fill(*g_.exterior.xy, color=c_, alpha=a_)
    _cl = {1: "tab:red", 2: "tab:blue", 3: "tab:green", 4: "tab:orange"}
    for r_, D_ in sorted(_TRAD_2D.items()):
        for g_ in _tteile(D_["L"]): _ax.plot(*g_.exterior.xy, color=_cl[r_], lw=1.6)
        for g_ in _tteile(_rot(D_["L"], -D_["th"], D_["ax"])): _ax.plot(*g_.exterior.xy, ":", color=_cl[r_], lw=1)
        for g_ in _tteile(D_["Ca"]): _ax.plot(*g_.exterior.xy, "--", color=_cl[r_], lw=1)
        _ax.fill(*_knopf_pg(r_).exterior.xy, color="k"); _ax.plot(*_rot(_knopf_pg(r_), -D_["th"], D_["ax"]).exterior.xy, color="gray", lw=0.8)
        _ax.plot(*_TPt(D_["ax"]).buffer(D_LAGER / 2).exterior.xy, "k", lw=0.8)
        if r_ in FED2D:
            F_ = FED2D[r_]; _ax.plot(*_TPt(F_["C"]).buffer(F_D / 2).exterior.xy, color="k", lw=1)
            for a_, b_ in ((F_["A0"], F_["A1"]), (F_["B0"], F_["B1"])): _ax.plot([a_[0], b_[0]], [a_[1], b_[1]], "k", lw=1.2)
            _dpl(_TPt(F_["C"]).buffer(F_D / 2), f"Feder_R{r_}")
        _ax.text(D_["knopf"][0] - 12, D_["knopf"][1], f"R{r_}", color=_cl[r_], fontsize=14, weight="bold")
    _ax.set_aspect("equal"); _ax.set_xlim(-110, -30); _ax.set_ylim(10, 150); _ax.grid(alpha=0.3)
    _ax.set_title("Kugellager-Griffstock: Teil 1 braun, Teil 2 gruen, Stuetze violett, Leisten blau, punktiert = gedrueckt")
    _fg.savefig(OUT.replace(".step", "_kugel_schnitt.png"), dpi=80, bbox_inches="tight")
    print("Kugellager: Schnitt ->", _dn)
except Exception as _ex:
    print("Kugellager: Schnitt-DXF nicht geschrieben:", _ex)

# --- Kollisionspruefung mit GEDRUECKTER Taste (3D) ---
if os.environ.get("KOLLISION"):
    import time as _tm
    _t0 = _tm.time()
    _d = {n_: s_ for n_, s_, f_ in teile}
    _fest = [n_ for n_ in _d if n_.startswith(("Griffstock_Koerper", "Griffstock_Teil2", "Stuetze_POM", "Leiste_FR2"))]
    _hebel = sorted(set(n_.split("/")[0] for n_ in _d if n_.startswith("Hebel_R") and "/" in n_))
    def _bb_ueber(a_, b_):
        A_, B_ = a_.BoundingBox(), b_.BoundingBox()
        return not (A_.xmax < B_.xmin or B_.xmax < A_.xmin or A_.ymax < B_.ymin or B_.ymax < A_.ymin or A_.zmax < B_.zmin or B_.zmax < A_.zmin)
    _kg = []
    for G in _hebel:
        r_ = int(G[7]); ax_ = ACHSE_K[r_]; th_ = H2D[r_]["th"]
        xk_, Y_, Z_ = DREHPUNKTE[G]
        bx_, by_ = _taff.rotate(_TPt(H2D[r_]["K"]), -th_, origin=ax_, use_radians=False).coords[0]
        u2_, Y2_ = _M((bx_, by_)); ku_, kY_ = _M(H2D[r_]["K"])
        best_ = None
        for sg_ in (1.0, -1.0):
            p_ = cq.Vector(xk_, kY_, -ku_).sub(cq.Vector(xk_, Y_, Z_)); a_ = math.radians(sg_ * th_)
            q_ = cq.Vector(xk_, Y_ + p_.y * math.cos(a_) - p_.z * math.sin(a_), Z_ + p_.y * math.sin(a_) + p_.z * math.cos(a_))
            e_ = math.hypot(q_.y - Y2_, q_.z + u2_)
            if best_ is None or e_ < best_[0]: best_ = (e_, sg_)
        sg_ = best_[1]
        beweg = [n_ for n_ in _d if n_.startswith(G + "/") and n_.split("/")[1] in ("Holzhebel", "Hebel_Alu", "Knopf")]
        gedr = {n_: _d[n_].rotate(cq.Vector(xk_, Y_, Z_), cq.Vector(xk_ + 1.0, Y_, Z_), sg_ * th_) for n_ in beweg}
        andere = _fest + [n_ for n_ in _d if n_.startswith("Hebel_R") and not n_.startswith(G + "/")]
        for n1_, s1_ in gedr.items():
            for n2_ in andere:
                s2_ = _d[n2_]
                if not _bb_ueber(s1_, s2_): continue
                try: _is = s1_.intersect(s2_); v_ = _is.Volume()
                except Exception: v_ = 0.0
                if v_ > 0.5:
                    _bb = _is.BoundingBox(); _c = (_bb.xmin + _bb.xmax) / 2, (_bb.ymin + _bb.ymax) / 2, (_bb.zmin + _bb.zmax) / 2
                    _sc = _Minv((-_c[2], _c[1]))
                    _kg.append((n1_, n2_, round(v_, 1), "X %.1f Schnitt (%.1f, %.1f)" % (_c[0], _sc[0], _sc[1])))
    print(f"Kugellager: Kollisionen bei GEDRUECKTER Taste (> 0,5 mm3): {len(_kg)}  ({_tm.time() - _t0:.0f} s)")
    for k_ in _kg: print("   KOLL_GEDRUECKT", k_)
