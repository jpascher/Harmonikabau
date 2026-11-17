# Griffstock-Geometrie im Schnitt (u, Y): Kunststoffteil-Taschen und Huellform fuer die massiven Enden
if "KUNSTSTOFF" not in globals():
    exec(open("schnitt_teile.py").read().split('\nif __name__ == "__main__":')[0])
from shapely.geometry import Polygon as _P, LineString as _L
import math as _m
KPG = _P(KUNSTSTOFF).buffer(0)
_c = list(KPG.exterior.coords)[:-1]
_ed = [(_c[i], _c[(i+1) % 4]) for i in range(4)]
_lang = sorted(_ed, key=lambda e: -_m.dist(e[0], e[1]))[:2]
OFFEN = max(_lang, key=lambda e: (e[0][1] + e[1][1]) / 2)          # offene Seite = Seite zu den Knoepfen (+Y)
_a, _b = OFFEN
_t = ((_b[0]-_a[0]) / _m.dist(_a, _b), (_b[1]-_a[1]) / _m.dist(_a, _b))
_n = (-_t[1], _t[0])
_mid = KPG.centroid
if (_mid.x - _a[0]) * _n[0] + (_mid.y - _a[1]) * _n[1] < 0: _n = (-_n[0], -_n[1])   # n zeigt ins Teil
def kst_s(p): return (p[0]-_a[0]) * _t[0] + (p[1]-_a[1]) * _t[1]
def kst_rechteck(s0, s1, d0, d1):
    return [(_a[0] + s*_t[0] + d*_n[0], _a[1] + s*_t[1] + d*_n[1]) for s, d in ((s0, d0), (s1, d0), (s1, d1), (s0, d1))]
# Huellform des Griffstocks (massive Enden): Rueckwand + Stufen der Knopfreihen (Stufe = 3,5 mm unter der Knopfkappe)
_kd = (KNOPF_RICHT[0] / _m.hypot(*KNOPF_RICHT), KNOPF_RICHT[1] / _m.hypot(*KNOPF_RICHT))
_kt = (_kd[1], -_kd[0])                                           # entlang der Stufe (u waechst)
if _kt[0] < 0: _kt = (-_kt[0], -_kt[1])
_kn = sorted(KNOPF, key=lambda q: q[0])                           # innen (kleines u) -> aussen
_stufen = []
for i, k in enumerate(_kn):
    p = (k[0] - 3.5 * _kd[0], k[1] - 3.5 * _kd[1])
    ua = (_kn[i-1][0] + k[0]) / 2 if i > 0 else k[0] - 8.0
    ub = (_kn[i+1][0] + k[0]) / 2 if i + 1 < len(_kn) else max(q[0] for q in GRIFFSTOCK)
    sa, sb = (ua - p[0]) / _kt[0], (ub - p[0]) / _kt[0]
    _stufen.append(((p[0] + sa*_kt[0], p[1] + sa*_kt[1]), (p[0] + sb*_kt[0], p[1] + sb*_kt[1])))
_front = [q for st in _stufen for q in st]
_oben = max(KUNSTSTOFF + DREIECK, key=lambda q: q[1])
_g = GRIFFSTOCK
_hinten_aussen = max(_g, key=lambda q: q[0] - q[1])               # unten aussen
_rueck_oben = min(_g, key=lambda q: q[0] + 0 * q[1])
HUELLE = _P([(_rueck_oben[0], _rueck_oben[1]), (min(q[0] for q in DREIECK), _oben[1])] + _front + [_hinten_aussen]).convex_hull.buffer(0)
HUELLE = HUELLE.union(_P(GRIFFSTOCK).buffer(0)).union(KPG).union(_P(DREIECK).buffer(0))
# Stufen wieder einschneiden (die konvexe Huelle hat sie geglaettet)
for (p0, p1) in _stufen:
    aussen = [p0, p1, (p1[0] + 60*_kd[0], p1[1] + 60*_kd[1]), (p0[0] + 60*_kd[0], p0[1] + 60*_kd[1])]
    HUELLE = HUELLE.difference(_P(aussen).buffer(0))
if HUELLE.geom_type != "Polygon": HUELLE = max(HUELLE.geoms, key=lambda g: g.area)
# Unterkanten von Griffstock, Dreieckleiste und Kunststoffteil in eine Ebene mit der Polster-Unterkante (Boden, u = 0)
def _gr(g): return max(g.geoms, key=lambda q: q.area) if g.geom_type != "Polygon" else g
def _bx(u0, u1, y0, y1): return _P([(u0, y0), (u1, y0), (u1, y1), (u0, y1)])
# Rueckwand-Innenseite (2.8, 7.5) -> (72.5, -5.6)
_ri0 = min((q for q in GRIFFSTOCK if q[0] < 10), key=lambda q: -q[1])
_ri1 = min((q for q in GRIFFSTOCK if 60 < q[0] < 75), key=lambda q: q[1])
def _linie_u(p, q, u): t = (u - p[0]) / (q[0] - p[0]); return p[1] + t * (q[1] - p[1])
def _schnitt(p1, p2, p3, p4):
    d = (p2[0]-p1[0])*(p4[1]-p3[1]) - (p2[1]-p1[1])*(p4[0]-p3[0])
    t = ((p3[0]-p1[0])*(p4[1]-p3[1]) - (p3[1]-p1[1])*(p4[0]-p3[0])) / d
    return (p1[0] + t*(p2[0]-p1[0]), p1[1] + t*(p2[1]-p1[1]))
_A = (0.0, _linie_u(_ri0, _ri1, 0.0))
_gt = sorted(GRIFFSTOCK, key=lambda q: q[0])[:2]
_ro0 = min((q for q in GRIFFSTOCK if q[0] < 10), key=lambda q: q[1]); _ro1 = min(GRIFFSTOCK, key=lambda q: q[1])
GS_E = _gr(_P(GRIFFSTOCK).buffer(0).union(_P([_ri0, _ro0, (0.0, _linie_u(_ro0, _ro1, 0.0)), _A]).buffer(0)))
# Kunststoffteil (Plattendicke unveraendert) oben bis an die Rueckwand verlaengert; Holz-Dreieck bleibt ein reines Dreieck:
# Boden, Rueckwand-Innenseite, Unterseite des Kunststoffteils
_unten = [q for q in KUNSTSTOFF if q not in OFFEN]                # Unterseite (gegenueber der Lagerseite)
_P2 = min(_unten, key=lambda q: q[0]); _P3 = max(_unten, key=lambda q: q[0])
_P1 = min(OFFEN, key=lambda q: q[0])
_C = _schnitt(OFFEN[0], OFFEN[1], _ri0, _ri1)
_E = _schnitt(_P2, _P3, _ri0, _ri1)
KST_E = _P([_P1, _C, _E, _P2]).buffer(0)
DRE_E = _P([_A, (0.0, _P2[1]), _E]).buffer(0)
HUELLE = _gr(HUELLE.union(_bx(0.0, 3.2, min(q[1] for q in GS_E.exterior.coords if q[0] < 0.1), max(q[1] for q in KUNSTSTOFF))))
HUELLE = _gr(HUELLE.union(GS_E).union(KST_E).union(DRE_E).buffer(0.01).buffer(-0.01))
if __name__ == "__main__" and "teile" not in globals():
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(9, 12))
    def pl(q, c, lw=1.5):
        q = list(q) + [q[0]]; ax.plot([p[1] for p in q], [p[0] for p in q], color=c, lw=lw)
    pl(list(HUELLE.exterior.coords)[:-1], 'k', 2); pl(GRIFFSTOCK, 'brown'); pl(KUNSTSTOFF, 'blue'); pl(DREIECK, 'green')
    for r in (1, 2, 3, 4):
        pl(HEBEL[r], 'gray', 0.6)
        c = lager_von(r); s = kst_s(c); pl(kst_rechteck(s - 7.3, s + 7.3, 0, 11), 'red', 1)
    for k in KNOPF: ax.plot(k[1], k[0], 'm^')
    ax.set_aspect('equal'); ax.invert_yaxis(); ax.grid(True, lw=0.3); ax.set_xlim(-30, 70); ax.set_ylim(130, -5)
    fig.savefig("gs_huelle.png", dpi=70)
    print("offene Seite", OFFEN, "Lager s:", [round(kst_s(lager_von(r)), 1) for r in (1, 2, 3, 4)])
