# Dreieck-Vorlage und einfache Vorlagensuche (zum Beschriften der Übungsdaten in tri_data.py)
import sys, glob, cv2, numpy as np, json
sys.path.insert(0, '../gsleser')
import gsread as G
def tri_t(d, filled=True):
    w = int(d * 1.3); h = int(d * 1.05); s = int(d * 1.6) | 1; c = s // 2
    a = np.zeros((s, s), np.uint8)
    pts = np.array([[c - w // 2, c + h // 2], [c + w // 2, c + h // 2], [c, c - h // 2]], np.int32)
    if filled: cv2.fillPoly(a, [pts], 255)
    else: cv2.polylines(a, [pts], True, 255, max(1, int(d * 0.15)))
    return a
def scan(path, thr=0.6):
    g = cv2.imread(path, 0); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    st = G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1]); out = []
    for si, ys in enumerate(st):
        d = float(np.mean(np.diff(ys))); nol = G.stems_off(G.remove_lines(bw, ys, d), d)
        for f in (True, False):
            t = tri_t(d, f); r = cv2.matchTemplate(nol, t, cv2.TM_CCOEFF_NORMED)
            yy, xx = np.where(r >= thr)
            for y, x in zip(yy, xx):
                cy = y + t.shape[0] // 2; cx = x + t.shape[1] // 2
                if ys[0] - 7 * d <= cy <= ys[4] + 7 * d:
                    out.append((float(r[y, x]), cx, cy, f, si, (ys[2] - cy) / (d / 2)))
    out.sort(reverse=True); keep = []
    for o in out:
        if all(abs(o[1] - k[1]) > 10 or abs(o[2] - k[2]) > 6 for k in keep): keep.append(o)
    return keep
if __name__ == '__main__':
    tot = 0
    for p in sorted(glob.glob('/tmp/philipp/norm/p-*.png')):
        k = scan(p); tot += len(k)
        print(p[-9:], len(k), 'voll', sum(1 for o in k if o[3]), 'hohl', sum(1 for o in k if not o[3]), [round(o[0], 2) for o in k[:3]], [round(o[0], 2) for o in k[-3:]])
    print('gesamt', tot)
