import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2, numpy as np, json, sys, pickle
sys.path.insert(0, HERE); import gsread as G; from train import patch, CL
INK = 0.12
MOD = HERE if _os.path.exists(HERE+'/heads.pkl') else _os.path.join(HERE, '..', 'gsleser')   # Forschungsordner: Modelle aus gsleser
clf = pickle.load(open(MOD+'/heads.pkl', 'rb'))
# Dreiecke (4. Reihe in der Notation von Elke Margetich, z. B. Philipp Ramoner): eigener kleiner Klassifikator auf denselben Ausschnitten
_TP = _os.environ.get('GS_TRI', MOD+'/tri.pkl')
TRI = pickle.load(open(_TP, 'rb')) if _os.path.exists(_TP) else None
TRI_THR, TRI_HEAD, TRI_SURE = 0.5, 0.3, 0.95   # Dreieck, wenn dessen Klassifikator ≥ 0,5 meldet und der Kopf-Klassifikator ≥ 0,3 (oder Dreieck ≥ 0,95)
CLX = CL + ['t_f', 't_h']
TRI_BOOST = float(_os.environ.get('GS_TRIB', '0.25'))   # Vorrang der Dreieck-Meldungen beim Ausdünnen

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
        # einmal je System rändern statt bei jedem Ausschnitt (gleiches Ergebnis wie train.patch, aber viel schneller)
        r = int(d * 0.9); bwp = cv2.copyMakeBorder(bw, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0); nolp = cv2.copyMakeBorder(nol, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0)
        def fpatch(pad, x, y):
            x0, y0 = int(round(x)) - r, int(round(y)) - r
            return cv2.resize(pad[y0 + 2 * r:y0 + 4 * r, x0 + 2 * r:x0 + 4 * r], (20, 20), interpolation=cv2.INTER_AREA).ravel()
        for p in range(-12, 13):   # Lage -13 (sechseinhalb Linien unter dem System) kommt nicht vor, dort stehen die Bassbuchstaben
            y = int(round(mid + 0.5 - p * d / 2))
            if not (0 <= y < bw.shape[0]): continue
            row = ink[y]
            for x in range(int(d * 2), bw.shape[1] - int(d), 2):
                if row[x] >= INK:
                    F.append(np.concatenate([fpatch(bwp, x, y), fpatch(nolp, x, y)])); cand.append((x, y, p, float(row[x])))
        if not cand: continue
        F = np.array(F); pr = clf.predict_proba(F)
        pt = TRI.predict_proba(F) if TRI is not None else None
        hits = []
        for i in range(len(cand)):
            tp = pt[i][1] + pt[i][2] if pt is not None else 0
            if cand[i][2] % 2 == 0 and tp >= TRI_THR and (1 - pr[i][0] >= TRI_HEAD or tp >= TRI_SURE):   # 4. Reihe liegt immer auf geraden Lagen
                hits.append((round(float(tp), 3), cand[i][:3], 5 if pt[i][1] >= pt[i][2] else 6))
            elif pr[i][0] < 1 - thr:
                hits.append((round(1 - pr[i][0], 3), cand[i][:3], int(np.argmax(pr[i][1:])) + 1))
        # Formgüte: Übereinstimmung mit einer Kopfvorlage genau an der Stelle (trennt echte Köpfe von Treffern eine Halbstufe daneben)
        T = G.templates(d)
        def good(x, y, c):
            t = T['o_f' if CLX[c].startswith('o') else ('t_f' if CLX[c].startswith('t') else 'x_f')]; hh, ww = t.shape
            y0, x0 = int(round(y)) - hh // 2, int(round(x)) - ww // 2
            win = nol[max(0, y0):y0 + hh, max(0, x0):x0 + ww].astype(np.float32)
            if win.shape != t.shape: return 0.0
            tt = t.astype(np.float32); return float((win * tt).sum() / (np.sqrt((win * win).sum() * (tt * tt).sum()) + 1e-6))
        ref = sorted(((sc, good(x, y, c), (x, y, p), c) for sc, (x, y, p), c in hits), key=lambda r: (-r[0], -r[1]))
        # Dreieck-Meldungen haben Vorrang vor runden Köpfen an derselben Stelle (der Kopf-Klassifikator hält Dreiecke oft für runde Köpfe)
        ref = sorted(ref, key=lambda r: -(r[0] * 0.5 + r[1] + (TRI_BOOST if CLX[r[3]].startswith('t') else 0)))
        keep = []
        for sc, gq, (x, y, p), c in ref:
            # Dreieck neben Kreuz oder rundem Kopf eine Lage tiefer/höher (Sekunde links und rechts vom Hals) darf näher stehen
            def apart(k):
                if abs(k[2][2] - p) > 1: return True
                lim = d * (1.3 if k[2][2] == p else (0.45 if (k[3] >= 5) != (c >= 5) else 0.8))
                return abs(x - k[2][0]) > lim
            if all(apart(k) for k in keep): keep.append((sc, gq, (x, y, p), c))
        keep = [(sc, xyp, c) for sc, gq, xyp, c in keep]
        for sc, (x, y, p), c in keep:
            out.append({'x': x / scale, 'y': y / scale, 'p': p, 't': CLX[c][0], 'fill': CLX[c].endswith('f'), 'sc': round(float(sc), 2), 'st': si, 'd': d / scale, 'mid': mid / scale})
    return out

if __name__ == '__main__':
    files = json.load(open(HERE+'/test.json'))
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
