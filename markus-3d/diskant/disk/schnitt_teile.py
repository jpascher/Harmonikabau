# Teile aus dem Griffstock-Schnitt in Modellkoordinaten (u, Y) bereitstellen
exec(open("schnitt_trafo.py").read().split('\nif __name__ == "__main__":')[0])
from ezdxf import bbox as _bb
def _pts(e, tol=0.05):
    q = [modell((p.x, p.y)) for p in ezdxf.path.make_path(e).flattening(tol)]
    out = []
    for p in q:
        if not out or math.hypot(p[0]-out[-1][0], p[1]-out[-1][1]) > 1e-3: out.append(p)
    if len(out) > 2 and math.hypot(out[0][0]-out[-1][0], out[0][1]-out[-1][1]) < 1e-3: out.pop()
    return out
_links = [e for e in SM if _bb.extents([e]).has_data and _bb.extents([e]).extmax.x < -80]
_heb = [e for e in _links if e.dxftype() == "LWPOLYLINE" and e.dxf.layer in ("Hebel-alu-1r", "Hebel-alu-2r", "Hebel-alu-3r", "Hebel-4r")]
LAGER = [modell((e.dxf.center.x, e.dxf.center.y)) for e in _links if e.dxftype() == "CIRCLE" and e.dxf.layer == "Hebel-4r"]
_groß = {}
for e in _heb:
    p = _pts(e); ys = [q[1] for q in p]; us = [q[0] for q in p]
    if max(us) > 50:                                   # Hebel mit kurzem Arm im Griffstock (je Layer der laengste in Y)
        if e.dxf.layer not in _groß or (max(ys)-min(ys)) > (max(q[1] for q in _groß[e.dxf.layer]) - min(q[1] for q in _groß[e.dxf.layer])):
            _groß[e.dxf.layer] = p
REIHE_LAYER = {1: "Hebel-alu-1r", 2: "Hebel-alu-2r", 3: "Hebel-alu-3r", 4: "Hebel-4r"}
HEBEL = {r: _groß[REIHE_LAYER[r]] for r in (1, 2, 3, 4)}
# Endstuecke (gekroepfter Teil zur Klappe): kleine Konturen u < 25, je Klappenreihe nach Y
_end = [_pts(e) for e in _heb if max(q[0] for q in _pts(e)) < 25]
def _reihe(y): return 1 if y > 152 else 2 if y > 117 else 3 if y > 82 else 4
ENDE = {}
for p in _end:
    ym = (min(q[1] for q in p) + max(q[1] for q in p)) / 2
    ENDE[_reihe(ym + 8)] = p
def _guete(c, p):
    d = sorted(math.hypot(c[0]-q[0], c[1]-q[1]) for q in p)
    nah = [x for x in d if x < 7.0]
    return -len(nah) if nah else 99
_frei = list(LAGER); LAGER_R = {}
for r in sorted((1, 2, 3, 4), key=lambda r: min(_guete(c, HEBEL[r]) for c in LAGER)):
    c = min(_frei, key=lambda c: _guete(c, HEBEL[r])); LAGER_R[r] = c; _frei.remove(c)
def lager_von(r): return LAGER_R[r]
# Griffstock-Koerper, Kunststoffteil, Dreieckleiste
_gk = [e for e in _links if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "Diskantgriffstockschnitt"]
GRIFFSTOCK = max((_pts(e) for e in _gk), key=len)
KUNSTSTOFF = [p for p in (_pts(e) for e in _gk) if len(p) == 4][0]
DREIECK = [_pts(e) for e in _links if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "Griffstock Schnitt"][0]
# Knoepfe: Einfuegepunkt + Achsrichtung (Schnitt -x)
_ax = modell((-1.0 + P_B[0], P_B[1])); _o = modell(P_B)
KNOPF_RICHT = (_ax[0] - _o[0], _ax[1] - _o[1])
KNOPF = sorted([modell((e.dxf.insert.x, e.dxf.insert.y)) for e in _links if e.dxftype() == "INSERT"], key=lambda q: q[0])
# Holzteile (Layer 0) direkt aus dem Schnitt
HOLZ0 = [_pts(e) for e in _links if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "0"]
# --- Mechanik absenken: Kunststoffteil (Dicke unveraendert) mit allen Hebeln, Lagern und Knoepfen um ABSENK zum Boden;
#     Hebel am Lager abgewinkelt (Klappenarm um THETA[r] gedreht), damit die Endstuecke wieder auf den Klappen liegen
ABSENK = min(q[0] for q in KUNSTSTOFF)
U_ENDE0 = {r: min(q[0] for q in ENDE[r]) for r in ENDE}
def _ab(pts): return [(q[0] - ABSENK, q[1]) for q in pts]
def _dicht2(pts, step=1.0):
    out = []
    for i in range(len(pts)):
        a, b = pts[i], pts[(i+1) % len(pts)]
        n = max(1, int(math.hypot(b[0]-a[0], b[1]-a[1]) / step))
        out += [(a[0]+(b[0]-a[0])*k/n, a[1]+(b[1]-a[1])*k/n) for k in range(n)]
    return out
HEBEL = {r: _ab(HEBEL[r]) for r in HEBEL}; ENDE = {r: _ab(ENDE[r]) for r in ENDE}
LAGER = _ab(LAGER); LAGER_R = {r: (c[0] - ABSENK, c[1]) for r, c in LAGER_R.items()}
KUNSTSTOFF = _ab(KUNSTSTOFF); KNOPF = _ab(KNOPF); HOLZ0 = [_ab(h) for h in HOLZ0]
def _rot(r, pts, th, gewicht=None):
    up, Yp = LAGER_R[r]; out = []
    for u, Y in pts:
        du, dY = u - up, Y - Yp
        w = 1.0 if gewicht is None else gewicht(du, dY)
        t = th * w; c, s = math.cos(t), math.sin(t)
        out.append((up + du*c + dY*s, Yp - du*s + dY*c))
    return out
def _w(du, dY):                       # 0 im Knopfarm (Richtung +u), 1 im Klappenarm (Richtung +Y); Uebergang am Lager
    phi = math.degrees(math.atan2(dY, du))
    return min(1.0, max(0.0, (phi - 25.0) / 40.0))
THETA = {}
for r in ENDE:
    lo, hi = 0.0, 0.3
    for _ in range(50):
        m = (lo + hi) / 2
        if min(q[0] for q in _rot(r, ENDE[r], m)) < U_ENDE0[r]: lo = m
        else: hi = m
    THETA[r] = (lo + hi) / 2
def biege(r, pts, voll=False): return _rot(r, pts, THETA[r], None if voll else _w)
HEBEL_GERADE = dict(HEBEL)
HEBEL = {r: biege(r, _dicht2(HEBEL[r], 1.0)) for r in HEBEL}
ENDE = {r: biege(r, ENDE[r], voll=True) for r in ENDE}
if __name__ == "__main__":
    for r in (1, 2, 3, 4):
        print(r, "Lager", [round(v, 1) for v in lager_von(r)], "Hebelpunkte", len(HEBEL[r]), "Ende", r in ENDE)
    print("Knoepfe", [[round(v, 1) for v in k] for k in KNOPF], "Richtung", [round(v, 3) for v in KNOPF_RICHT])
    print("Holz0", [(len(h), round(min(q[0] for q in h)), round(max(q[0] for q in h))) for h in HOLZ0])
