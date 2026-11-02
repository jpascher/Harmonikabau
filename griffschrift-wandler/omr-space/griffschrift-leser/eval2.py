import json, sys
sys.path.insert(0, '/tmp/claude-0/gsr'); from detect import detect
tot = {}
for f in json.load(open('/tmp/claude-0/gsr/test.json')):
    gt = json.load(open(f)); kind = 'fuchsartig' if '/gs2/' in f else 'wandler'
    det = detect(gt['img'], gt['scale']); mid = gt['staff'][2]; used = set(); t = tot.setdefault(kind, dict(gt=0, det=0, hit=0, pos=0, typ=0))
    for h in gt['heads']:
        for i, x in enumerate(det):
            if i not in used and abs(x['x'] - h['x']) + abs(x['y'] - h['y']) < 6:
                used.add(i); t['hit'] += 1; t['pos'] += x['p'] == round((mid - h['y']) / 4); t['typ'] += x['t'] == h['t']; break
    t['gt'] += len(gt['heads']); t['det'] += len(det)
for k, t in tot.items(): print(k, t, 'Trefferquote', round(t['hit'] / t['gt'], 3), 'Genauigkeit', round(t['hit'] / max(1, t['det']), 3))
