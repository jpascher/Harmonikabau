# Von Hand bestätigte Köpfe einer Fotoseite -> Merkmale wie detect.py (für das Nachlernen des Kopfmodells auf Fotos)
import sys, json, cv2, numpy as np, random, importlib
sys.path.insert(0, '/home/claude/gsw/m'); sys.path.insert(0, '/home/claude/lab')
import gsread as G
random.seed(5)
CL = ['none', 'o_f', 'o_h', 'x_f']
def run(key):
    lab = importlib.import_module('labels_' + key); g = cv2.imread(f'/home/claude/lab/{key}.png', 0); det = json.load(open(f'/home/claude/lab/{key}_det.json'))
    G.FOTO_LINES = True
    _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV); st = G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1])
    XS = []; YS = []; heads = []
    for si, ys in enumerate(st):
        d = float(np.mean(np.diff(ys))); mid = ys[2]; nol = G.stems_off(G.remove_lines(bw, ys, d), d); r = int(d * 0.9)
        bwp = cv2.copyMakeBorder(bw, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0); nolp = cv2.copyMakeBorder(nol, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0)
        def fp(pad, x, y):
            x0, y0 = int(round(x)) - r, int(round(y)) - r
            return cv2.resize(pad[y0 + 2 * r:y0 + 4 * r, x0 + 2 * r:x0 + 4 * r], (20, 20), interpolation=cv2.INTER_AREA).ravel()
        F = lambda x, y: np.concatenate([fp(bwp, x, y), fp(nolp, x, y)])
        ypos = lambda p: int(round(mid + 0.5 - p * d / 2))
        H = [dict(x=h['x'], p=h['p'], c=('x_f' if h['t'] == 'x' else ('o_f' if h['fill'] else 'o_h')) if h['t'] in 'ox' else 'o_f') for h in det if h['st'] == si]
        for h in H:
            if h['x'] > lab.XMAX or h['x'] < getattr(lab, 'XMIN', 0): h['c'] = 'none'
        for e in lab.L.get(si, []):
            x, p, c = e[:3]; fix = len(e) > 3   # vierter Eintrag: Lage genau so übernehmen
            cand = [h for h in H if h['p'] == p and abs(h['x'] - x) < 0.8 * d]
            if cand: min(cand, key=lambda h: abs(h['x'] - x))['c'] = c
            else:   # ergänzen: genaue Lage über die meiste Tinte in der Nähe
                y = ypos(p); best = x if fix else max(range(int(x - 0.6 * d), int(x + 0.6 * d)), key=lambda xx: int(nol[max(0, y - 4):y + 5, max(0, xx - 6):xx + 7].sum()))
                H.append(dict(x=float(best), p=p, c=c))
        pos = [(h['x'], ypos(h['p'])) for h in H if h['c'] != 'none']
        for h in H:
            y = ypos(h['p'])
            for dx, dy in ((0, 0), (-2, 0), (2, 1), (0, -1)) if h['c'] != 'none' else ((0, 0), (-2, 0), (2, 0)):
                XS.append(F(h['x'] + dx, y + dy)); YS.append(CL.index(h['c']))
            if h['c'] != 'none': heads.append((si, h['x'], h['p'], h['c']))
        for (x, y) in pos:
            for ax, ay in ((0, d / 2), (0, -d / 2), (d * 0.9, 0), (-d * 0.9, 0), (0, d), (0, -d)):
                xx, yy = x + ax, y + ay
                if any(abs(xx - a) < d * 0.5 and abs(yy - b) < d * 0.3 for a, b in pos): continue
                XS.append(F(xx, yy)); YS.append(0)
        for _ in range(len(pos) * 4):
            x = random.uniform(max(d * 3, getattr(lab, 'XMIN', 0)), min(lab.XMAX, bw.shape[1] - d)); p = random.randint(-12, 12); y = ypos(p)
            if any(abs(x - a) < d * 0.6 and abs(y - b) < d * 0.3 for a, b in pos): continue
            XS.append(F(x, y)); YS.append(0)
    return np.array(XS, np.uint8), np.array(YS), heads
if __name__ == '__main__':
    X, Y, heads = run(sys.argv[1]); np.savez_compressed(f'/home/claude/lab/{sys.argv[1]}_ds.npz', X=X, Y=Y)
    json.dump(heads, open(f'/home/claude/lab/{sys.argv[1]}_heads.json', 'w')); print(X.shape, np.bincount(Y, minlength=4), 'Köpfe', len(heads))
