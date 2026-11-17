# DXF-Hebel (Holzrand-Zeichnungen) per Drehung/Spiegelung um den Lagerpunkt auf den Hebel im Schnitt legen
exec(open("schnitt_teile.py").read().split('\nif __name__ == "__main__":')[0])
exec(open("hebel2d.py").read().split('\nif __name__')[0])
def _nah(p, pts): return min(math.hypot(p[0]-q[0], p[1]-q[1]) for q in pts)
def _dicht(pts, step=1.0):
    out = []
    for i in range(len(pts)):
        a, b = pts[i], pts[(i+1) % len(pts)]
        n = max(1, int(math.hypot(b[0]-a[0], b[1]-a[1]) / step))
        out += [(a[0]+(b[0]-a[0])*k/n, a[1]+(b[1]-a[1])*k/n) for k in range(n)]
    return out
REG = {}
for r in (1, 2, 3, 4):
    alu, ring = hebel_profil(r)
    up, Yp = lager_von(r)
    ziel = _dicht(HEBEL_GERADE[r], 1.0)
    ziel = [q for q in ziel if q[0] > up + 8]                       # kurzer Arm im Griffstock
    quelle = [q for q in _dicht(alu, 2.0) if q[0] > 10]             # kurzer Arm im DXF (x > 10)
    best = None
    for sp in (1, -1):
        for k in range(0, 720):
            th = math.radians(k * 0.5)
            c_, s_ = math.cos(th), math.sin(th)
            f = [(up + x*c_ - sp*y*s_, Yp + x*s_ + sp*y*c_) for x, y in quelle]
            err = sum(_nah(p, ziel) for p in f[::3]) / len(f[::3])
            if best is None or err < best[0]: best = (err, sp, th)
    REG[r] = best
    if __name__ == "__main__": print(r, "Fehler %.2f mm  Spiegel %d  Drehung %.1f" % (best[0], best[1], math.degrees(best[2])))
def trafo(r, pts):
    err, sp, th = REG[r]; up, Yp = lager_von(r); c_, s_ = math.cos(th), math.sin(th)
    return biege(r, [(up + x*c_ - sp*y*s_, Yp + x*s_ + sp*y*c_) for x, y in pts])
