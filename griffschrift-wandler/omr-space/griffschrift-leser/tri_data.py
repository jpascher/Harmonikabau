# Arbeitsordner /tmp/philipp: pages/ (Seiten aus Philipps Heft, 150 dpi), norm/ (auf Linienabstand 18 gebracht), ph_pins.json (Griffe aus dem Wandler);
# die Seiten selbst sind urheberrechtlich geschützt und liegen nicht im Repository.
# Übungsdaten für die Dreieck-Köpfe (4. Reihe, Notation Elke Margetich) aus Philipps Heft
# Positiv: Stellen, an denen die bekannten Griffe einen Knopf der 4. Reihe verlangen und ein Dreieck steht.
# Negativ: richtig gelesene runde Köpfe und Kreuze, Stellen eine Lage daneben, zufällige Tintenstellen.
import sys, json, cv2, numpy as np, random, pickle
GS = '../gsleser'
sys.path.insert(0, GS); sys.path.insert(0, '/tmp/philipp')
import gsread as G
from train import patch
from tri_proto import tri_t
from philipp_pruefen import PAGES, PINS, pin_chords
random.seed(2)
OFF = {1: -1, 2: 0, 3: 1, 4: 2}
def page_ctx(p):
    g = cv2.imread(p if isinstance(p, str) else f'/tmp/philipp/norm/p-{p:02d}.png', 0); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    st = G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1]); ctx = []
    for ys in st:
        d = float(np.mean(np.diff(ys))); nol = G.stems_off(G.remove_lines(bw, ys, d), d)
        r = int(d * 0.9)
        ctx.append(dict(ys=ys, d=d, nol=nol, mid=ys[2], r=r, bwp=cv2.copyMakeBorder(bw, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0),
                        nolp=cv2.copyMakeBorder(nol, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0)))
    return bw, ctx
def _fp(pad, r, x, y):
    x0, y0 = int(round(x)) - r, int(round(y)) - r
    return cv2.resize(pad[y0 + 2 * r:y0 + 4 * r, x0 + 2 * r:x0 + 4 * r], (20, 20), interpolation=cv2.INTER_AREA).ravel()
def feat(bw, c, x, y): return np.concatenate([_fp(c['bwp'], c['r'], x, y), _fp(c['nolp'], c['r'], x, y)]).astype(np.float32)
def ypos(c, p): return int(round(c['mid'] + 0.5 - p * c['d'] / 2))
def best_tri(c, x, y):
    d = c['d']; best = (0, None, None)
    for f in (True, False):
        t = tri_t(d, f); h = t.shape[0] // 2
        y0, y1 = int(y - 0.75 * d) - h, int(y + 0.35 * d) + h + 1; x0, x1 = int(x - 1.6 * d) - h, int(x + 1.6 * d) + h
        win = c['nol'][max(0, y0):y1, max(0, x0):x1]
        if win.shape[0] < t.shape[0] or win.shape[1] < t.shape[1]: continue
        r = cv2.matchTemplate(win, t, cv2.TM_CCOEFF_NORMED); _, mv, _, ml = cv2.minMaxLoc(r)
        if mv > best[0]: best = (mv, max(0, x0) + ml[0] + h, f)
    return best
INK = 0.12
NEGIMG = None
def candidates(bw, c):
    d = c['d']; nol = c['nol']; ink = cv2.boxFilter((nol > 0).astype(np.float32), -1, (max(3, int(d * 0.9)), max(3, int(d * 0.6)))); out = []
    for p in range(-12, 13, 2):
        y = ypos(c, p)
        if not (0 <= y < bw.shape[0]): continue
        row = ink[y]
        for x in range(int(d * 2), bw.shape[1] - int(d), 2):
            if row[x] >= INK: out.append((x, y, p))
    return out
def build(pieces, R, prev=None):
    X, Y = [], []; npos = 0; POS = {}
    for k in pieces:
        res = R[k]['res']; B = pin_chords(PINS[k])
        cache = {}
        for i, j in R[k]['pairs']:
            ch = res[i]; p = ch['page']
            if p not in cache: cache[p] = page_ctx(p)
            bw, ctx = cache[p]
            if ch['st'] >= len(ctx): continue
            c = ctx[ch['st']]; d = c['d']; x = ch['x']
            want = {(r, n) for r, n, _ in B[j]}
            w4 = {2 * (n - 6) + 2 for r, n in want if r == 4}
            for (r, n) in want:
                pp = 2 * (n - 6) + OFF[r]; y = ypos(c, pp)
                if r == 4:
                    sc, bx, f = best_tri(c, x, y)
                    if sc < 0.3: continue
                    for dx, dy in ((0, 0), (-2, 0), (2, 0), (0, 1), (0, -1)):
                        X.append(feat(bw, c, bx + dx, y + dy)); Y.append(1 if f else 2)
                    npos += 1; POS.setdefault((p, ch['st']), []).append((bx, y))
                    for ay, dp in ((d / 2, -1), (-d / 2, 1), (d, -2), (-d, 2)):
                        if pp + dp not in w4: X.append(feat(bw, c, bx, y + ay)); Y.append(0)
            # richtig gelesene runde Köpfe / Kreuze als Negative; falsch gemeldete Dreiecke als harte Negative
            for (hp, ht, (r, n)) in ch['heads']:
                y = ypos(c, hp)
                if (r, n) in want and r < 4:
                    X.append(feat(bw, c, x, y)); Y.append(0)
                elif ht == 't' and (r, n) not in want:
                    for dx, dy in ((0, 0), (-2, 0), (2, 0), (0, 1), (0, -1)): X.append(feat(bw, c, x + dx, y + dy)); Y.append(0)
        # harte Negative: Stellen, an denen der bisherige Klassifikator ein Dreieck meldet, die Griffe aber keines verlangen
        for pv in (prev or []):
            for p, (bw, ctx) in cache.items():
                for si, c in enumerate(ctx):
                    cs = candidates(bw, c)
                    if not cs: continue
                    F = np.array([feat(bw, c, x, y) for x, y, _ in cs]); pt = pv.predict_proba(F); pos = POS.get((p, si), [])
                    for (x, y, pp), q, f in zip(cs, pt, F):
                        if q[1] + q[2] >= 0.3 and all(abs(x - bx) > 0.7 * c['d'] or abs(y - by) > 3 for bx, by in pos): X.append(f); Y.append(0)
        # zufällige Tintenstellen
        for p, (bw, ctx) in cache.items():
            for c in ctx:
                ink = np.argwhere(c['nol'] > 0)
                for _ in range(150):
                    yy, xx = ink[random.randrange(len(ink))]; X.append(feat(bw, c, xx, yy)); Y.append(0)
    # Bilder ohne Dreiecke (Wandler-Griffschrift in Fuchs-Schreibweise): jede Meldung ist falsch
    for f in (NEGIMG or []):
        bw, ctx = page_ctx(f)
        for c in ctx:
            cs = candidates(bw, c)
            if not cs: continue
            F = np.array([feat(bw, c, x, y) for x, y, _ in cs])
            hit = np.zeros(len(cs), bool)
            for pv in (prev or []):
                q = pv.predict_proba(F); hit |= (q[:, 1] + q[:, 2]) >= 0.2
            for f_ in F[hit]: X.append(f_); Y.append(0)
            for i in random.sample(range(len(cs)), min(60, len(cs))): X.append(F[i]); Y.append(0)
    print('Dreiecke', npos, 'Beispiele', len(Y), 'davon positiv', sum(1 for y in Y if y))
    return np.array(X), np.array(Y)
if __name__ == '__main__':
    R = json.load(open(sys.argv[1])); train = sys.argv[2].split(',')
    prev = [pickle.load(open(f, 'rb')) for f in sys.argv[4].split(',')] if len(sys.argv) > 4 else None
    import glob
    if len(sys.argv) > 5: NEGIMG = [f for pat in sys.argv[5].split(',') for f in sorted(glob.glob(pat))]
    X, Y = build(train, R, prev)
    from sklearn.ensemble import HistGradientBoostingClassifier
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True, random_state=1).fit(X, Y)
    pickle.dump(clf, open(sys.argv[3], 'wb')); print('gespeichert', sys.argv[3], clf.classes_)
