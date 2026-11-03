# Lerndaten für die Bassbuchstaben: Philipps Heft (Buchstaben samt Lage aus dem PDF) und Wandler-Bilder (Buchstaben aus dem SVG)
import sys, json, glob, cv2, numpy as np, collections
sys.path.insert(0, '/home/claude/gsw/m'); sys.path.insert(0, '/home/claude/bass')
import gs2xml as X, gsread as G, bassread as BR
LET = set('ABCDEFGHabcdefgh')
def staves_of(g):
    return G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1])
def label_page(g, boxes, src):
    """boxes: Liste (x0,y0,x1,y1,text) in Koordinaten von g"""
    _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV); st = staves_of(g); XS = []; YS = []
    for si, ys in enumerate(st):
        d = (ys[4] - ys[0]) / 4; C = BR.candidates(bw, st, si, d)
        inside = collections.defaultdict(list)
        for ci, (box, f) in enumerate(C):
            cx, cy = box[0] + box[2] / 2, box[1] + box[3] / 2
            for bi, (x0, y0, x1, y1, t) in enumerate(boxes):
                if x0 - 2 <= cx <= x1 + 2 and y0 - 2 <= cy <= y1 + 2: inside[bi].append(ci); break
        lab = ['none'] * len(C)
        for bi, cis in inside.items():
            t = boxes[bi][4]
            if t == "'":
                for ci in cis: lab[ci] = 'mark'
                continue
            big = max(cis, key=lambda ci: C[ci][0][2] * C[ci][0][3])
            for ci in cis:
                if ci == big: lab[ci] = t[0] if t[0] in LET else 'none'
                elif len(t) > 1 and t[1] in "²'": lab[ci] = 'mark'
        for (box, f), l in zip(C, lab): XS.append(f); YS.append(l)
    return XS, YS
XS = []; YS = []; KS = []
# Philipp
PH = json.load(open('/home/claude/bass/ph_letters.json')); byp = collections.defaultdict(list)
for p, l, bb in PH: byp[p].append((l, bb))
for p in sorted(byp):
    g = cv2.imread(f'/home/claude/philipp/norm/p-{p:02d}.png', 0); s = g.shape[1] / 606.614
    boxes = [(bb[0] * s, bb[1] * s, bb[2] * s, bb[3] * s, l) for l, bb in byp[p] if l in LET or l == "'"]
    x, y = label_page(g, boxes, 'P'); XS += x; YS += y; KS += [f'P:{p}'] * len(y)
print('Philipp', len(YS), collections.Counter(YS).most_common(), flush=True)
# Wandler
for f in sorted(glob.glob('/home/claude/bass/wb/*.json')):
    d = json.load(open(f)); o = cv2.imread(f[:-5] + '.png', 0)
    g = X.load_gray(f[:-5] + '.png')
    if g is None: continue
    s = g.shape[1] / o.shape[1]
    boxes = [(b['x0'] * s, b['y0'] * s, b['x1'] * s, b['y1'] * s, b['t']) for b in d['L']]
    x, y = label_page(g, boxes, 'W'); XS += x; YS += y; KS += ['W:' + d['key']] * len(y)
print('gesamt', len(YS), collections.Counter(YS).most_common(), flush=True)
np.savez_compressed('/home/claude/bass/ds.npz', X=np.array(XS), Y=np.array(YS), K=np.array(KS))
