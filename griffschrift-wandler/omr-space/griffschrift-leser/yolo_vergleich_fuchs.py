# Vergleich bisheriger Kopf-Erkenner (detect.py) gegen gelernten Erkenner (yolo_koepfe.py) auf Fuchs-Seiten mit bekannten Griffen
# (midi/vm/<stück>.pins.json). Geprüft wird nur an Seiten, die nicht gelernt wurden (yolo_pruefseiten.txt: 30 zufällig gewählte,
# keine der 19 Lernseiten aus yolo_fuchs_seiten.py). Kennzahlen: Akkorde ganz gleich, Köpfe gefunden, Genauigkeit.
# Die Fuchs-Blätter (GS*.gif von volksmusik.cc) liegen nicht im Repository: FUCHSORDNER/<stück>gs/GS*.gif.
# Aufruf aus diesem Ordner: GS_YOLOMODELL=…/griffschrift_koepfe.onnx python3 yolo_vergleich_fuchs.py FUCHSORDNER [yolo_pruefseiten.txt]
import sys, os, glob, json, cv2
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'gsleser'))
import gs2xml as X, detect as D, yolo_koepfe as YK
from align import chords_det, chords_pin, align, upscale
PINS = os.path.join(HERE, '..', '..', 'midi', 'vm')


def mess(det, B):
    C = chords_det(det); A = [{(h['p'], h['t']) for h in c} for c in C]; pr = align(A, B)
    ex = sum(A[i] == B[j] for i, j in pr); tp = sum(len(A[i] & B[j]) for i, j in pr)
    return dict(pins=len(B), exact=ex, tp=tp, fp=sum(len(a) for a in A) - tp, fn=sum(len(b) for b in B) - tp)


if __name__ == '__main__':
    fu = sys.argv[1]; liste = sys.argv[2] if len(sys.argv) > 2 else os.path.join(HERE, 'yolo_pruefseiten.txt')
    T = {m: dict(pins=0, exact=0, tp=0, fp=0, fn=0) for m in ('bisher', 'yolo')}
    for k in [l.strip() for l in open(liste) if l.strip()]:
        gs = glob.glob(os.path.join(fu, k + 'gs', 'GS*'))
        if not gs: print(k, 'Blatt fehlt'); continue
        u = os.path.join(fu, f'u_{k}.png')
        if not os.path.exists(u): upscale(gs[0], u)
        g = X.load_gray(u)
        if g is None: print(k, 'keine Notenlinien'); continue
        w = os.path.join(fu, f'w_{k}.png'); cv2.imwrite(w, g); st, _ = X.druck_bars(g)
        B = chords_pin(json.load(open(os.path.join(PINS, k + '.pins.json')))['n'])
        r = {'bisher': mess(D.detect(w, 1), B), 'yolo': mess(YK.detect(g, st), B)}
        for m in r:
            for x in r[m]: T[m][x] += r[m][x]
        print(k, '  '.join(f"{m}: ganz {r[m]['exact']}/{r[m]['pins']}" for m in r), flush=True)
    for m, t in T.items():
        print(m, 'Akkorde ganz %.3f' % (t['exact'] / t['pins']), 'Köpfe gefunden %.3f' % (t['tp'] / (t['tp'] + t['fn'])),
              'Genauigkeit %.3f' % (t['tp'] / (t['tp'] + t['fp'])))
