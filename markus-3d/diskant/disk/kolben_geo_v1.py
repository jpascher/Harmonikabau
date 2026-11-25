# Griffstock "kolben" (22.11.2026, Johann): verdeckter Griffstock mit Kunststoffzylindern (Kolben) unter den Knoepfen,
# ohne Nuten. Hebel nur Alu (keine Holzteile); Knopfarme aus griffstock-alternativ-50-tasten-dez-2018-45°-leiste-
# ueberarbeitet2 (kolben_schnitt.dxf). Die Zeichnung steht auf 76 Grad (Rueckwand zum Boden) - auf 79 Grad wie das
# bisherige Modell gedreht: Rueckwand deckungsgleich mit gschnitt. Lagerung, Kunststoffteil und Drehfedern wie beim
# Kunststoff-Griffstock; jeder Alu-Hebel wird nur auf das Lager des Kunststoffteils verschoben (ca. 3-4 mm), nicht gedreht.
# Eingebunden von hebel.py (exec) nach dem JP-Block, vor der Hebelschleife.
import ezdxf as _kz
from ezdxf import bbox as _kbb
from shapely.geometry import Polygon as _KP, LineString as _KL, Point as _KPt
from shapely import affinity as _kaff
from shapely.ops import unary_union as _kuu
import numpy as _knp
import schnitt_trafo as _KT_
_KM = _kz.readfile("kolben_schnitt.dxf").modelspace()
# Zeichnung -> gschnitt: Fusspunkt der Rueckwand (629,42 / -1,12) auf den Rueckwand-Fuss in gschnitt, +90 Grad
# (Rueckwand der Zeichnung waagrecht, in gschnitt senkrecht); danach entlang des Bodens so verschoben, dass die
# Rueckwand-Aussenseiten zusammenfallen (-38,08 mm in Modell-Y, ausgemessen an beiden Griffstock-Umrissen)
_KF = _knp.array((629.42, -1.12))
_KFg = _knp.array((-135.24, -10.2 + (280.56 - 135.24) * math.tan(math.radians(11.0))))
_KR = _knp.array([[0.0, -1.0], [1.0, 0.0]])
_KDY = -38.08
def _kmod(p):
    q_ = (_knp.array(p[:2], float) - _KF) @ _KR.T + _KFg
    u_, Y_ = _KT_.modell((float(q_[0]), float(q_[1]))); return (u_ - ABSENK, Y_ + _KDY)
def _kpg(e):
    pts_ = [_kmod((v_.x, v_.y)) for v_ in _kz.path.make_path(e).flattening(0.03)]
    return _KP(pts_).buffer(0)
_KHAUPT = lambda e: (lambda b_: b_.extmin.x > 590 and b_.extmin.y > 0)(_kbb.extents([e]))
KOLBEN_LAYER = {1: "Hebel-alu-1r", 2: "Hebel-alu-2r", 3: "Hebel-alu-3r", 4: "Hebel-4r"}
_KHEB = {r_: max((_kpg(e) for e in _KM.query(f'LWPOLYLINE[layer=="{l_}"]') if _KHAUPT(e)), key=lambda g_: g_.area)
         for r_, l_ in KOLBEN_LAYER.items()}
_KLAG = [_kmod((c_.dxf.center.x, c_.dxf.center.y)) for c_ in _KM.query('CIRCLE[layer=="Hebel-4r"]') if c_.dxf.center.x > 590]
D_KOLBEN, L_KOLBEN_MIN = 7.5, 6.0
U_KOLBEN_AUSTRITT = 30.0                                # Abdeckung erst ab u 30 (wie Zeichnung): darunter Austritt der Hebel
T_ABDECK = 4.0                                            # Abdeckung (Knopfseite, mit Stufen), fest mit den Seitenteilen
KOLBEN_ALU, KOLBEN_KONTAKT, KOLBEN_DREH, KOLBEN_UNTEN = {}, {}, {}, {}
HEBEL_ALT = {r_: list(HEBEL[r_]) for r_ in HEBEL}
for r_ in (1, 2, 3, 4):
    pg_ = _KHEB[r_]
    # eigenes Lager: der Kreis, um den die meisten Konturpunkte im Abstand < 7 liegen (Nabe)
    def _nah(c_): return -sum(1 for q_ in pg_.exterior.coords if math.dist(c_, q_) < 7.0)
    lk_ = min(_KLAG, key=_nah)
    up_, Yp_ = lager_von(r_)
    pg_ = _kaff.translate(pg_, up_ - lk_[0], Yp_ - lk_[1])          # Lager auf das Lager des Kunststoffteils
    ku_, kY_ = knopf_u[r_ - 1]
    if JP: ku_, kY_ = ku_ + JP_KNOPF_D[r_][0], kY_ + JP_KNOPF_D[r_][1]
    kd_ = _knp.array(KD)
    def _ende(g_):                                        # Armende = Konturpunkt mit dem groessten u (Knopfarm)
        return max(g_.exterior.coords, key=lambda q_: q_[0])
    def _abst(g_):                                        # Armende auf der Knopfachse, 0..30 mm unter dem Knopf
        e_ = _knp.array(_ende(g_)) - _knp.array((ku_, kY_))
        t_ = -(e_[0] * kd_[0] + e_[1] * kd_[1])
        return abs(e_[0] * kd_[1] - e_[1] * kd_[0]) + max(0.0, -t_) + max(0.0, t_ - 30.0)
    best_ = min(((_abst(_kaff.rotate(pg_, w_ / 10.0, origin=(up_, Yp_))), w_ / 10.0) for w_ in range(-150, 151)))
    best_ = (_abst(pg_), 0.0)      # 22.11.: Armenden liegen ohne Drehung unter den Knoepfen (Kolben 12-16 mm) - nicht drehen
    print(f"   Kolben Reihe {r_}: Lager Zeichnung {tuple(round(v_, 1) for v_ in lk_)} -> {round(up_, 1), round(Yp_, 1)}, Armende vorher {tuple(round(v_, 1) for v_ in _ende(pg_))}, Rest {best_[0]:.2f}")
    KOLBEN_DREH[r_] = best_[1]
    pg_ = _kaff.rotate(pg_, best_[1], origin=(up_, Yp_))
    # Knopfarm (Richtung +u vom Lager) aus der neuen Zeichnung, Klappenarm (Richtung +Y) aus dem bisherigen Hebel
    def _keil(w0_, w1_):
        pts_ = [(up_, Yp_)] + [(up_ + 400 * math.cos(math.radians(w_)), Yp_ + 400 * math.sin(math.radians(w_))) for w_ in range(w0_, w1_ + 1, 5)]
        return _KP(pts_).buffer(0)
    k_knopf = _keil(-120, 40); k_klappe = _keil(40, 240)
    alt_ = _KP(HEBEL[r_]).buffer(0)
    neu_ = _kuu([pg_.intersection(k_knopf), alt_.intersection(k_klappe), pg_.intersection(_KPt(up_, Yp_).buffer(6.0)),
                 alt_.intersection(_KPt(up_, Yp_).buffer(6.0))]).buffer(0.05).buffer(-0.05)
    neu_ = max(neu_.geoms, key=lambda g_: g_.area) if hasattr(neu_, "geoms") else neu_
    KOLBEN_ALU[r_] = neu_
    # Beruehrpunkt Kolben/Alu: von der Knopfunterseite entlang der Knopfachse nach innen bis aufs Alu
    s_ = _KL([(ku_, kY_), (ku_ - 80 * kd_[0], kY_ - 80 * kd_[1])]).intersection(neu_)
    pts_ = [c_ for g_ in (getattr(s_, "geoms", None) or [s_]) if not g_.is_empty for c_ in g_.coords]
    KOLBEN_KONTAKT[r_] = max(pts_, key=lambda c_: (c_[0] - ku_) * kd_[0] + (c_[1] - kY_) * kd_[1]) if pts_ else None
    # 22.11. 11:40 (Johann): der Zylinder reicht ueber den Hebel - der verjuengte Hebelteil sitzt in einer gleich grossen
    # Ausnehmung im Zylinder. Zylinder nach innen bis hinter das Alu im Zylinderquerschnitt (+1,5 mm Boden)
    _band = _KL([(ku_, kY_), (ku_ - 80 * kd_[0], kY_ - 80 * kd_[1])]).buffer(D_KOLBEN / 2, cap_style=2).intersection(neu_)
    _tt = [-((c_[0] - ku_) * kd_[0] + (c_[1] - kY_) * kd_[1]) for g_ in (getattr(_band, "geoms", None) or [_band]) if not g_.is_empty
           for c_ in g_.exterior.coords]
    _tb = (max(_tt) + 1.5) if _tt else None
    KOLBEN_UNTEN[r_] = (ku_ - _tb * kd_[0], kY_ - _tb * kd_[1]) if _tb else KOLBEN_KONTAKT[r_]
    HEBEL[r_] = [(float(a_), float(b_)) for a_, b_ in list(neu_.simplify(0.02).exterior.coords)[:-1]]
print("Kolben-Griffstock: Hebel um ihr Lager gedreht (Grad)", {r_: KOLBEN_DREH[r_] for r_ in KOLBEN_DREH},
      "Kolbenlaenge", {r_: (round(math.dist(KOLBEN_UNTEN[r_], (knopf_u[r_ - 1][0] + (JP_KNOPF_D[r_][0] if JP else 0),
                                                                 knopf_u[r_ - 1][1] + (JP_KNOPF_D[r_][1] if JP else 0))), 1)
                                     if KOLBEN_KONTAKT[r_] else None) for r_ in KOLBEN_KONTAKT})
if os.environ.get("KOLBEN_DEBUG"):
    import pickle as _kpk
    _kpk.dump({"neu": KOLBEN_ALU, "kontakt": KOLBEN_KONTAKT, "knopf": {r_: (knopf_u[r_ - 1][0] + (JP_KNOPF_D[r_][0] if JP else 0), knopf_u[r_ - 1][1] + (JP_KNOPF_D[r_][1] if JP else 0)) for r_ in (1, 2, 3, 4)},
               "lager": {r_: lager_von(r_) for r_ in (1, 2, 3, 4)}, "KD": KD, "HUELLE": globals().get("HUELLE"), "alt": {r_: _KP(HEBEL_ALT[r_]).buffer(0) for r_ in HEBEL_ALT}}, open("/tmp/kolben_dbg.pkl", "wb"))
    raise SystemExit("KOLBEN_DEBUG")
