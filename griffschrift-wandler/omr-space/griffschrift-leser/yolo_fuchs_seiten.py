# Echte Lernseiten für den Kopf-Erkenner (YOLO) aus Fuchs-Blättern, deren Griffe aus den MIDI-Sätzen bekannt sind (pins.json):
# der bisherige Erkenner findet die Köpfe; Lernrahmen werden die Köpfe, deren Lage und Form zu den bekannten Griffen passt, und die
# übrigen, wenn der Erkenner sehr sicher ist (für die Form zählt die Tonhöhe nicht; nicht zugeordnet sind oft nur Wiederholungen,
# die der MIDI-Satz anders ausschreibt). Unsichere, nicht bestätigte Funde bleiben ohne Rahmen im Bild.
# Seiten, auf denen zu viele bekannte Griffe fehlen, werden ausgelassen. Nur jede zweite Seite (nach Namen) wird genommen,
# die übrigen bleiben für die Prüfung. Die Blätter selbst kommen nicht ins Repository.
# Aufruf: python3 yolo_fuchs_seiten.py FUCHSORDNER PINSORDNER AUSGABEORDNER [PRÜFLISTE]
#   FUCHSORDNER: <name>gs/GS*.gif bzw. u_<name>.png (hochskaliert), PINSORDNER: <name>.pins.json
import sys, os, glob, json
import numpy as np, cv2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'gsleser'))
import gs2xml as X, detect as D
from align import chords_det, chords_pin, align, upscale

# Rahmengröße je Form in Linienabständen (Breite, Höhe), wie bei den vom Wandler gesetzten Lernbildern
GR = {0: (1.45, 1.15), 1: (1.45, 1.15), 2: (1.2, 1.2), 3: (1.35, 1.35), 4: (1.35, 1.1), 5: (1.35, 1.1)}


def klasse(h):
    if h['t'] == 'o': return 0 if h['fill'] else 1
    if h['t'] == 'x': return 2 if h['fill'] else 3
    return 4 if h['fill'] else 5


def seite(img, pins, mindest=0.9, sicher=0.9):
    g = X.load_gray(img)
    if g is None: return None
    work = '/tmp/_yfs.png'; cv2.imwrite(work, g)
    det = D.detect(work, 1); C = chords_det(det); B = chords_pin(pins)
    A = [{(h['p'], h['t']) for h in c} for c in C]; pr = align(A, B)
    gut = set(); tp = 0
    for i, j in pr:
        for h in C[i]:
            if (h['p'], h['t']) in B[j]: gut.add(id(h)); tp += 1
    nB = sum(len(b) for b in B)
    if nB == 0 or tp / nB < mindest: return dict(skip=True, recall=round(tp / max(1, nB), 3))
    boxen = []; weg = 0
    for c in C:
        for h in c:
            d = h['d']; k = klasse(h); w, hh = GR[k][0] * d, GR[k][1] * d
            if id(h) in gut or h.get('sc', 0) >= sicher: boxen.append((k, h['x'] - w / 2, h['y'] - hh / 2, w, hh))
            else: weg += 1
    return dict(g=g, boxen=boxen, recall=round(tp / nB, 3), weg=weg)


if __name__ == '__main__':
    fu, pn, ziel = sys.argv[1:4]; os.makedirs(ziel, exist_ok=True)
    namen = sorted(os.path.basename(f.rstrip('/'))[:-2] for f in glob.glob(os.path.join(fu, '*gs/')))
    # mit Prüfliste (viertes Argument, z. B. yolo_pruefseiten.txt): alle Seiten außer denen der Liste; sonst jede zweite wie bisher
    pruef = {l.strip() for l in open(sys.argv[4]) if l.strip()} if len(sys.argv) > 4 else None
    for n, k in enumerate(namen):
        if (k in pruef) if pruef is not None else n % 2: continue   # bleibt für die Prüfung
        if not os.path.exists(os.path.join(pn, k + '.pins.json')): continue
        img = os.path.join(os.path.dirname(fu.rstrip('/')), f'u_{k}.png')
        if not os.path.exists(img): upscale(glob.glob(os.path.join(fu, k + 'gs', 'GS*'))[0], img)
        r = seite(img, json.load(open(os.path.join(pn, k + '.pins.json')))['n'])
        if r is None or r.get('skip'): print(k, 'ausgelassen', r and r['recall']); continue
        cv2.imwrite(os.path.join(ziel, f'fuchs_{k}.png'), r['g'])
        json.dump(dict(boxes=[[int(c), round(float(x), 1), round(float(y), 1), round(float(w), 1), round(float(h), 1)] for c, x, y, w, h in r['boxen']]), open(os.path.join(ziel, f'fuchs_{k}.json'), 'w'))
        print(k, 'Köpfe', len(r['boxen']), 'ohne Rahmen', r['weg'], 'gefunden', r['recall'])
