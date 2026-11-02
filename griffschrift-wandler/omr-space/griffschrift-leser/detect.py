import cv2, numpy as np, json, sys, pickle
sys.path.insert(0, '/tmp/claude-0/gsr'); import gsread as G; from train import patch, CL
INK = 0.12
clf = pickle.load(open('/tmp/claude-0/gsr/heads.pkl', 'rb'))

def detect(path, scale=2, thr=0.5):
    g = cv2.imread(path, 0); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    out = []
    # Linien suchen mit hellerer Schwelle: gedruckte Notenlinien sind oft grau
    _, bwl = cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)
    for si, ys in enumerate(G.staves(bwl)):
        d = float(np.mean(np.diff(ys))); nol = G.stems_off(G.remove_lines(bw, ys, d), d); mid = ys[2]
        # Kandidaten: alle Lagen und Spalten, an denen nach Entfernen von Linien und Hälsen genug Tinte in Kopfgröße liegt
        ink = cv2.boxFilter((nol > 0).astype(np.float32), -1, (max(3, int(d * 0.9)), max(3, int(d * 0.6))))
        cand = []; F = []
        for p in range(-13, 13):
            y = int(round(mid + 0.5 - p * d / 2))
            if not (0 <= y < bw.shape[0]): continue
            row = ink[y]
            for x in range(int(d * 2), bw.shape[1] - int(d), 2):
                if row[x] >= INK:
                    F.append(np.concatenate([patch(bw, x, y, d), patch(nol, x, y, d)])); cand.append((x, y, p, float(row[x])))
        if not cand: continue
        pr = clf.predict_proba(np.array(F))
        hits = [(round(1 - pr[i][0], 3), cand[i][:3], int(np.argmax(pr[i][1:])) + 1) for i in range(len(cand)) if pr[i][0] < 1 - thr]
        # Formgüte: Übereinstimmung mit einer Kopfvorlage genau an der Stelle (trennt echte Köpfe von Treffern eine Halbstufe daneben)
        T = G.templates(d)
        def good(x, y, c):
            t = T['o_f' if CL[c].startswith('o') else 'x_f']; hh, ww = t.shape
            y0, x0 = int(round(y)) - hh // 2, int(round(x)) - ww // 2
            win = nol[max(0, y0):y0 + hh, max(0, x0):x0 + ww].astype(np.float32)
            if win.shape != t.shape: return 0.0
            tt = t.astype(np.float32); return float((win * tt).sum() / (np.sqrt((win * win).sum() * (tt * tt).sum()) + 1e-6))
        ref = sorted(((sc, good(x, y, c), (x, y, p), c) for sc, (x, y, p), c in hits), key=lambda r: (-r[0], -r[1]))
        ref = sorted(ref, key=lambda r: -(r[0] * 0.5 + r[1]))
        keep = []
        for sc, gq, (x, y, p), c in ref:
            if all(abs(x - k[2][0]) > d * (1.3 if k[2][2] == p else 0.8) or abs(k[2][2] - p) > 1 for k in keep): keep.append((sc, gq, (x, y, p), c))
        keep = [(sc, xyp, c) for sc, gq, xyp, c in keep]
        for sc, (x, y, p), c in keep:
            out.append({'x': x / scale, 'y': y / scale, 'p': p, 't': CL[clf.classes_[c]][0] if False else CL[c][0], 'fill': CL[c].endswith('f'), 'sc': round(float(sc), 2), 'st': si, 'd': d / scale, 'mid': mid / scale})
    return out

if __name__ == '__main__':
    files = json.load(open('/tmp/claude-0/gsr/test.json'))
    tot = dict(gt=0, det=0, hit=0, pos=0, typ=0)
    for f in files:
        gt = json.load(open(f[:-4] + '.json')); det = detect(f); mid = gt['staff'][2]; used = set()
        for h in gt['heads']:
            best = None
            for i, x in enumerate(det):
                if i in used: continue
                dd = abs(x['x'] - h['x']) + abs(x['y'] - h['y'])
                if dd < 6 and (best is None or dd < best[0]): best = (dd, i)
            if best:
                used.add(best[1]); x = det[best[1]]; tot['hit'] += 1
                tot['pos'] += x['p'] == round((mid - h['y']) / 4); tot['typ'] += x['t'] == h['t']
        tot['gt'] += len(gt['heads']); tot['det'] += len(det)
    print(tot, 'Trefferquote', round(tot['hit'] / tot['gt'], 3), 'Genauigkeit', round(tot['hit'] / max(1, tot['det']), 3))
