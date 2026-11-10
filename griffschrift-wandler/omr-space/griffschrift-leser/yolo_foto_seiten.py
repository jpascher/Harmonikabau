# Echte Lernseiten für den Kopf-Erkenner (YOLO) aus Fotos und Scans ohne bekannte Lösung (dritter Lernlauf).
# Vorbeschriftet wird mit zwei Erkennern: Rahmen bekommen die Köpfe, bei denen der bisherige Erkenner (detect.py) sehr sicher ist,
# und die Köpfe, die nur das gelernte Modell (yolo_koepfe.py, falls vorhanden) mit hoher Sicherheit meldet. Wo einer der beiden
# etwas Unsicheres meldet, das der andere nicht bestätigt, wird die Stelle im Bild weiß abgedeckt: so lernt das Modell dort weder
# „Kopf“ noch „kein Kopf“. Vorgesetzte Doppelkreuze (Rosenzopf) kommen nur vom gelernten Modell; liest der bisherige Erkenner an
# derselben Stelle einen Kopf, gilt das Doppelkreuz. Für die Form zählt das Urteil des bisherigen Erkenners, weil das gelernte Modell
# auf Fotos volle Köpfe öfter als hohle liest.
# Ausgabe wie yolo_fuchs_seiten.py (PNG + JSON mit boxes), dazu "foto": true, damit yolo_daten.py die Seite nicht noch einmal
# wie ein Foto verschlechtert. Die Blätter selbst kommen nicht ins Repository.
# Aufruf: [GS_YOLOMODELL=…/griffschrift_koepfe.onnx] python3 yolo_foto_seiten.py AUSGABEORDNER BILD [BILD …]
import sys, os, json
import numpy as np, cv2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'gsleser'))
import gs2xml as X, detect as D, gsread as G

GR = {0: (1.45, 1.15), 1: (1.45, 1.15), 2: (1.2, 1.2), 3: (1.35, 1.35), 4: (1.35, 1.1), 5: (1.35, 1.1), 6: (1.0, 1.0)}
SICHER = float(os.environ.get('FS_SICHER', '0.9'))      # bisheriger Erkenner: ab hier Rahmen
YSICHER = float(os.environ.get('FS_YSICHER', '0.7'))    # gelerntes Modell: ab hier Rahmen, wo der bisherige nichts meldet
YUNTEN = float(os.environ.get('FS_YUNTEN', '0.3'))      # gelerntes Modell: darunter gilt eine Meldung als nichts


def klasse(h):
    if h['t'] == 'o': return 0 if h['fill'] else 1
    if h['t'] == 'x': return 2 if h['fill'] else 3
    return 4 if h['fill'] else 5


def seite(path):
    g = X.load_gray(path)
    if g is None: return None
    work = '/tmp/_yfo.png'; cv2.imwrite(work, g)
    det = [h for h in D.detect(work, 1) if abs(h['y'] - h['mid']) < 9 * h['d']]
    yk = []
    try:
        import yolo_koepfe as YK
        if YK.vorhanden(): yk = [k for k in YK.koepfe(g, conf=YUNTEN)]
    except Exception as e: print('  Modell nicht nutzbar:', str(e)[:80])
    d0 = float(np.median([h['d'] for h in det])) if det else 18.0
    nah = lambda x1, y1, x2, y2, f=0.6: abs(x1 - x2) < f * d0 and abs(y1 - y2) < f * d0
    vork = [(cx, cy) for cx, cy, w, h, k, sc in yk if k == 6 and sc >= 0.5]
    boxen = []; weiss = []
    for h in det:
        x, y, d = h['x'], h['y'], h['d']
        if any(nah(x, y, vx, vy) for vx, vy in vork): continue        # Doppelkreuz, kein Kopf
        bestaetigt = any(nah(x, y, cx, cy) and k != 6 for cx, cy, w, hh, k, sc in yk)
        if h.get('sc', 0) >= SICHER or (bestaetigt and h.get('sc', 0) >= 0.6):
            k = klasse(h); w, hh = GR[k][0] * d, GR[k][1] * d; boxen.append((k, x - w / 2, y - hh / 2, w, hh))
        else: weiss.append((x, y, d))
    for cx, cy, w, hh, k, sc in yk:
        if any(nah(cx, cy, h['x'], h['y']) for h in det) and k != 6: continue
        if sc >= YSICHER: boxen.append((k, cx - w / 2, cy - hh / 2, w, hh))
        else: weiss.append((cx, cy, d0))
    g2 = g.copy(); hell = int(np.percentile(g, 90))
    for x, y, d in weiss:   # unsichere Stellen mit Papierhelligkeit abdecken, aber nicht über einen beschrifteten Kopf
        x0, y0, x1, y1 = int(x - 0.8 * d), int(y - 0.65 * d), int(x + 0.8 * d), int(y + 0.65 * d)
        if any(bx < x1 and bx + bw > x0 and by < y1 and by + bh > y0 for _, bx, by, bw, bh in boxen): continue
        g2[max(0, y0):y1, max(0, x0):x1] = hell
    return dict(g=g2, boxen=boxen, weiss=len(weiss), alt=len(det), yolo=len(yk), foto=bool(G.FOTO_LINES))


if __name__ == '__main__':
    ziel = sys.argv[1]; os.makedirs(ziel, exist_ok=True)
    for p in sys.argv[2:]:
        k = os.path.splitext(os.path.basename(p))[0].replace(' ', '_')
        try: r = seite(p)
        except Exception as e: print(k, 'Fehler', str(e)[:80]); continue
        if r is None or len(r['boxen']) < 30: print(k, 'ausgelassen', r and len(r['boxen'])); continue
        cv2.imwrite(os.path.join(ziel, f'foto_{k}.png'), r['g'])
        json.dump(dict(foto=True, boxes=[[int(c), round(float(x), 1), round(float(y), 1), round(float(w), 1), round(float(h), 1)] for c, x, y, w, h in r['boxen']]),
                  open(os.path.join(ziel, f'foto_{k}.json'), 'w'))
        print(k, 'Köpfe', len(r['boxen']), 'abgedeckt', r['weiss'], '(bisher', r['alt'], 'gelernt', r['yolo'], ')', 'Foto' if r['foto'] else 'Scan', flush=True)
