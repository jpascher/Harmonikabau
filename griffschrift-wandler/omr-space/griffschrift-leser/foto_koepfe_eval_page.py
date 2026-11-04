# Köpfe einer bestätigten Fotoseite gegen die Lesung (detect auf der gespeicherten, normierten Seite)
import sys, json, os, importlib
sys.path.insert(0, '/home/claude/gsw/m')
import gsread as G; G.FOTO_LINES = True
import detect as D
key = sys.argv[1]; T = json.load(open(f'/home/claude/lab/{key}_heads.json')); lab = importlib.import_module('labels_' + key) if False else None
sys.path.insert(0, '/home/claude/lab'); lab = importlib.import_module('labels_' + key)
det = [h for h in D.detect(f'/home/claude/lab/{key}.png', 1) if getattr(lab, 'XMIN', 0) <= h['x'] <= lab.XMAX]
used = set(); hit = typ = 0; hol = holhit = 0
for si, x, p, c in T:
    best = None
    for i, h in enumerate(det):
        if i in used or h['st'] != si or abs(h['p'] - p) > 0 or abs(h['x'] - x) > 0.8 * h['d']: continue
        if best is None or abs(h['x'] - x) < abs(det[best]['x'] - x): best = i
    if c == 'o_h': hol += 1
    if best is not None:
        used.add(best); hit += 1; h = det[best]; tc = 'x_f' if h['t'] == 'x' else ('o_f' if h['fill'] else 'o_h'); typ += tc == c
        if c == 'o_h' and tc == 'o_h': holhit += 1
print(key, 'Köpfe', len(T), 'gefunden', hit, 'Typ richtig', typ, '| hohl', hol, 'gefunden+hohl', holhit, '| erkannt gesamt', len(det), 'überzählig', len(det) - hit)
