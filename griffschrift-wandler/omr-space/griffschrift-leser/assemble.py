# Griffschrift-Zeile -> Griffe (Knopf, Reihe, Balgrichtung, Töne)
import cv2, numpy as np, json, sys
sys.path.insert(0, '/tmp/claude-0/gsr'); import gsread as G
from detect import detect
D446 = {'druck': {1: [98,98,0,4,7,12,16,19,24,28,31,36,40,43,48,98], 2: [98,98,98,5,9,12,17,21,24,29,33,36,41,45,48,98,98], 3: [98,98,98,10,14,17,22,26,29,34,38,41,46,50,98,98], 4: [98,98,98,15,19,22,27,31,34,39,43,46,51,98,98]},
        'zug': {1: [98,98,5,7,11,14,17,21,23,26,29,33,35,38,41,98], 2: [98,98,98,10,12,16,19,22,24,28,31,34,38,40,43,98,98], 3: [98,98,98,15,17,21,24,27,29,33,36,39,43,45,98,98], 4: [98,98,98,20,22,26,29,32,34,38,41,44,48,98,98]}}
REF = {1: 7, 2: 8, 3: 8, 4: 8}; OFF = {1: -1, 2: 0, 3: 1, 4: 2}; BASE = 43   # G-C-F-B

def button(p, t):
    rows = (1, 2) if t == 'o' else (3, 4)
    r = rows[0] if (p - OFF[rows[0]]) % 2 == 0 else rows[1]
    return r, (p - OFF[r]) // 2 + 6

def pitch(r, num, dr):
    i = num - 6 + REF[r]; arr = D446['druck' if dr == 'D' else 'zug'][r]
    return BASE + arr[i] if 0 <= i < len(arr) and arr[i] != 98 else None

_ocr = None
def ocr():
    global _ocr
    if _ocr is None:
        from rapidocr import RapidOCR
        _ocr = RapidOCR(params={"Global.log_level": "error"})
    return _ocr

def read_line(path, scale=2):
    g = cv2.imread(path, 0); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    st = G.staves(bw)
    if not st: return [], []
    ys = st[0]; d = float(np.mean(np.diff(ys)))
    hs = detect(path, scale)
    # Druckstriche: lange waagrechte dünne Striche unterhalb des Systems (unter den Bassbuchstaben)
    low = int(ys[-1] + 3 * d); band = bw[low:, :]
    hl = cv2.morphologyEx(band, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_RECT, (int(d * 2.5), 1)))
    n, lab, stats, _ = cv2.connectedComponentsWithStats(hl)
    under = [(stats[i, 0] / scale, (stats[i, 0] + stats[i, 2]) / scale, (low + stats[i, 1]) / scale) for i in range(1, n) if stats[i, 3] <= d * 0.4]
    # Akkorde: Köpfe nach x gruppieren (Sekunden sind seitlich versetzt)
    hs.sort(key=lambda h: h['x']); chords = []
    for h in hs:
        if chords and h['x'] - chords[-1]['x1'] < d * 1.5 / scale and abs(h['x'] - chords[-1]['x0']) < d * 1.6 / scale: chords[-1]['h'].append(h); chords[-1]['x1'] = h['x']
        else: chords.append({'x0': h['x'], 'x1': h['x'], 'h': [h]})
    out = []
    for c in chords:
        x = sum(h['x'] for h in c['h']) / len(c['h'])
        dr = 'D' if any(a - 4 <= x <= b + 8 for a, b, _ in under) else 'Z'
        btns = [button(h['p'], h['t']) for h in c['h']]
        out.append({'x': x, 'dir': dr, 'btns': btns, 'pitches': sorted(p for p in (pitch(r, n, dr) for r, n in btns) if p is not None)})
    return out, under

if __name__ == '__main__':
    files = json.load(open('/tmp/claude-0/gsr/test.json')) if len(sys.argv) < 2 else sys.argv[1:]
    tot = dict(gt=0, found=0, dir=0, btn=0, pitch=0)
    for f in files:
        gt = json.load(open(f[:-4] + '.json'))
        if 'grips' not in gt: continue
        out, _ = read_line(f)
        for g in gt['grips']:
            tot['gt'] += 1
            c = min(out, key=lambda o: abs(o['x'] - g['x'])) if out else None
            if not c or abs(c['x'] - g['x']) > 10: continue
            tot['found'] += 1; tot['dir'] += c['dir'] == g['dir']
            tot['btn'] += sorted(map(tuple, c['btns'])) == sorted(map(tuple, g['btns']))
            tot['pitch'] += c['pitches'] == sorted(g['pitches'])
    print(tot)
