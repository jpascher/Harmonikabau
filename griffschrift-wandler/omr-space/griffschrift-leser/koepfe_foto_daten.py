# Lerndaten für ein zweites Kopfmodell: Wandler-Bilder (Köpfe aus dem SVG), zum Teil künstlich wie Handyfotos verschlechtert
# (blasse dünne Striche, Unschärfe, Rauschen); Merkmale genau wie detect.py (Ausschnitt mit und ohne Linien/Hälse, je 20×20)
import sys, json, glob, cv2, numpy as np, random
sys.path.insert(0, '/home/claude/gsw/m')
import gsread as G
random.seed(3); rng = np.random.default_rng(3)
CL = ['none', 'o_f', 'o_h', 'x_f']
def fp(pad, r, x, y):
    x0, y0 = int(round(x)) - r, int(round(y)) - r
    return cv2.resize(pad[y0 + 2 * r:y0 + 4 * r, x0 + 2 * r:x0 + 4 * r], (20, 20), interpolation=cv2.INTER_AREA).ravel()
def degrade(g):
    g = g.astype(np.float32)
    if random.random() < 0.8:
        g = 255 - (255 - g) * random.uniform(0.4, 0.9)   # blasser Druck
        g = cv2.GaussianBlur(g, (0, 0), random.uniform(0.4, 1.3))
        if random.random() < 0.5:   # dünnere Striche: leicht aufhellen an Kanten
            g = cv2.dilate(g, np.ones((2, 2), np.uint8)) * 0.5 + g * 0.5
        g = g + rng.normal(0, random.uniform(2, 9), g.shape)
    return np.clip(g, 0, 255).astype(np.uint8)
XS = []; YS = []; KS = []
files = sorted(glob.glob(sys.argv[1] + '/*.json'))
for fi, f in enumerate(files):
    gt = json.load(open(f)); o = cv2.imread(f[:-5] + '.png', 0)
    if o is None: continue
    st = G.staves(cv2.threshold(o, 240, 255, cv2.THRESH_BINARY_INV)[1])
    if not st or len(gt['staff']) < 5: continue
    ys0 = st[0]; s = float(np.mean(np.diff(ys0))) / float(np.mean(np.diff(gt['staff'][:5]))); oy = ys0[0] - s * gt['staff'][0]
    d0 = float(np.mean(np.diff(ys0))); k = 18.0 / d0
    g = cv2.resize(o, None, fx=k, fy=k, interpolation=cv2.INTER_CUBIC)
    for rep in range(2):
        gg = degrade(g); _, bw = cv2.threshold(gg, 200, 255, cv2.THRESH_BINARY_INV)
        ys = [y * k for y in ys0]; d = 18.0
        nol = G.stems_off(G.remove_lines(bw, ys, d), d); r = int(d * 0.9)
        bwp = cv2.copyMakeBorder(bw, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0); nolp = cv2.copyMakeBorder(nol, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0)
        F = lambda x, y: np.concatenate([fp(bwp, r, x, y), fp(nolp, r, x, y)])
        pos = []
        for h in gt['heads']:
            x = h['x'] * s * k; y = (h['y'] * s + oy) * k
            if not (2 * d < x < bw.shape[1] - d and 0 < y < bw.shape[0]): continue
            c = 'x_f' if h['t'] == 'x' else ('o_f' if h['fill'] else 'o_h'); pos.append((x, y))
            for dx, dy in ((0, 0), (-2, 0), (2, 1), (0, -1)):
                XS.append(F(x + dx, y + dy)); YS.append(CL.index(c)); KS.append(fi)
        mid = ys[2]
        for (x, y) in pos:
            for ax, ay in ((0, d / 2), (0, -d / 2), (d * 0.9, 0), (-d * 0.9, 0), (0, d), (0, -d)):
                xx, yy = x + ax, y + ay
                if any(abs(xx - a) < d * 0.5 and abs(yy - b) < d * 0.3 for a, b in pos): continue
                XS.append(F(xx, yy)); YS.append(0); KS.append(fi)
        for _ in range(len(pos) * 2):
            x = random.uniform(d * 3, bw.shape[1] - d); p = random.randint(-12, 12); y = mid + 0.5 - p * d / 2
            if any(abs(x - a) < d * 0.6 and abs(y - b) < d * 0.3 for a, b in pos): continue
            XS.append(F(x, y)); YS.append(0); KS.append(fi)
    if fi % 50 == 0: print(fi, len(YS), flush=True)
X = np.array(XS, np.uint8); Y = np.array(YS); np.savez_compressed(sys.argv[2], X=X, Y=Y, K=np.array(KS)); print(X.shape, np.bincount(Y))
