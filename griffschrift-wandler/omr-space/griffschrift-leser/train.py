# Kopf-Klassifikator: Ausschnitte an den möglichen Lagen (halbe Linienabstände) -> keiner / rund voll / rund hohl / Kreuz voll / Kreuz hohl
import cv2, numpy as np, json, glob, random, sys, pickle
sys.path.insert(0, '/tmp/claude-0/gsr'); import gsread as G
random.seed(1); np.random.seed(1)
P = 20   # Ausschnittgröße nach Skalierung
CL = ['none', 'o_f', 'o_h', 'x_f', 'x_h']

def prep(path):
    g = cv2.imread(path, 0); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    st = G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1])  # Linien oft grau
    return bw, st

def patch(img, x, y, d):
    r = int(d * 0.9); x0, y0 = int(round(x)) - r, int(round(y)) - r
    pad = cv2.copyMakeBorder(img, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0)
    c = pad[y0 + 2 * r:y0 + 4 * r, x0 + 2 * r:x0 + 4 * r]
    return cv2.resize(c, (P, P), interpolation=cv2.INTER_AREA).ravel()

def feats(bw, ys, x, y, d):
    # zwei Ansichten: mit Linien (Original) und ohne Linien/Hälse
    return np.concatenate([patch(bw, x, y, d), patch(G.stems_off(G.remove_lines(bw, ys, d), d), x, y, d)])

if __name__ == '__main__':
    X, Y = [], []
    files = sorted(glob.glob('/tmp/claude-0/gs/z*.json') + glob.glob('/tmp/claude-0/gs2/f*.json')); random.shuffle(files)
    train = files[:int(len(files) * 0.8)]; json.dump(files[int(len(files) * 0.8):], open('/tmp/claude-0/gsr/test.json', 'w'))
    for f in train:
        gt = json.load(open(f)); bw, st = prep(gt['img']); sc = gt['scale']
        if not st: continue
        ys = st[0]; d = float(np.mean(np.diff(ys))); nol = G.stems_off(G.remove_lines(bw, ys, d), d); mid = ys[2]
        pos = []
        for h in gt['heads']:
            x, y = h['x'] * sc, h['y'] * sc; pos.append((x, y))
            c = h['t'] + ('_f' if h['t'] == 'x' or h.get('fill', True) else '_h')
            for dx, dy in ((-2, 0), (0, 0), (2, 1), (0, -1)):
                X.append(np.concatenate([patch(bw, x + dx, y + dy, d), patch(nol, x + dx, y + dy, d)])); Y.append(CL.index(c))
        # harte Negative: eine Lage darüber/darunter und knapp daneben, wenn dort kein Kopf ist
        for (x, y) in list(pos):
            for (ax, ay) in ((0, d / 2), (0, -d / 2), (d * 0.9, 0), (-d * 0.9, 0), (0, d), (0, -d), (1, d / 2 + 1), (-1, -d / 2 - 1), (2, d / 2 - 1), (-2, -d / 2 + 1)):
                xx, yy = x + ax, y + ay
                if any(abs(xx - a) < d * 0.5 and abs(yy - b) < d * 0.3 for a, b in pos): continue
                X.append(np.concatenate([patch(bw, xx, yy, d), patch(nol, xx, yy, d)])); Y.append(0)
        # Negative an Taktnummern, Klammerzahlen, Wiederholungspunkten, Pausen
        for q in gt.get('neg', []):
            for ax, ay in ((0, 0), (-3, 0), (3, 0), (0, -3), (0, 3)):
                X.append(np.concatenate([patch(bw, q['x'] * sc + ax, q['y'] * sc + ay, d), patch(nol, q['x'] * sc + ax, q['y'] * sc + ay, d)])); Y.append(0)
        # Negative an den Bassbuchstaben und Druckstrichen unter dem System
        for bsy in gt['bass']:
            for ay in (-6, -3, 0, 3):
                for ax in (-3, 0, 3):
                    X.append(np.concatenate([patch(bw, bsy['x'] * sc + ax, (bsy['y'] - 5) * sc + ay * sc, d), patch(nol, bsy['x'] * sc + ax, (bsy['y'] - 5) * sc + ay * sc, d)])); Y.append(0)
        # Negative: zufällige Stellen auf Lagen, nicht an Köpfen
        W = bw.shape[1]
        for _ in range(len(gt['heads']) * 3):
            x = random.uniform(d * 3, W - d); p = random.randint(-13, 12); y = mid + 0.5 - p * d / 2 + random.choice((-1, 0, 1))
            if any(abs(x - a) < d * 0.6 and abs(y - b) < d * 0.3 for a, b in pos): continue
            X.append(np.concatenate([patch(bw, x, y, d), patch(nol, x, y, d)])); Y.append(0)
    X = np.array(X, dtype=np.uint8); Y = np.array(Y)
    neg = np.where(Y == 0)[0]; posi = np.where(Y != 0)[0]; keepn = np.random.choice(neg, min(len(neg), int(1.4 * len(posi))), replace=False)
    sel = np.concatenate([posi, keepn]); X = X[sel]; Y = Y[sel]
    np.savez('/tmp/claude-0/gsr/synth_cache.npz', X=X, Y=Y)
    # echte, mit bekannten Griffen abgeglichene Fuchs-Seiten (realdata.py), dreifach gewichtet
    import os
    if os.path.exists('/tmp/claude-0/gsr/real_train.npz'):
        R = np.load('/tmp/claude-0/gsr/real_train.npz'); X = np.concatenate([X] + [R['X']] * 2); Y = np.concatenate([Y] + [R['Y']] * 2)
    Y[Y == 4] = 3   # hohle Kreuze sind zu selten; mit dieser Klasse lernt HistGradientBoosting gar nicht mehr
    print(X.shape, np.bincount(Y)); import gc; gc.collect()   # uint8 direkt übergeben spart Speicher
    if os.environ.get('NUR_CACHE'): sys.exit()
    from sklearn.ensemble import HistGradientBoostingClassifier
    clf = HistGradientBoostingClassifier(max_iter=300, learning_rate=0.1, early_stopping=True).fit(X, Y)
    pickle.dump(clf, open('/tmp/claude-0/gsr/heads.pkl', 'wb')); print('ok')
