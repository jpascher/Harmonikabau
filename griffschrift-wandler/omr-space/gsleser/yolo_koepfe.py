# Gelernter Kopf-Erkenner (YOLO, in Colab mit yolo_colab.ipynb gelernt): findet in einer Griffschrift-Seite die Köpfe nach ihrer Form.
# Klassen: 0 rund voll, 1 rund hohl, 2 Kreuz, 3 Kreuz im Kreis, 4 Dreieck voll, 5 Dreieck hohl, 6 vorgesetztes Doppelkreuz.
# Liefert Einträge im Format von detect.detect (x, y, p, t, fill, sc, st, d, mid). Die Seite wird in Ausschnitte von 960 Bildpunkten
# mit Überlappung zerlegt, wie beim Lernen. Modell: griffschrift_koepfe.onnx neben dieser Datei (GS_YOLOMODELL: anderer Pfad).
# Liegen statt der ganzen Datei nur Teile griffschrift_koepfe.onnx.00, .01 … daneben (über die Weboberfläche von Hugging Face
# lassen sich nur kleinere Dateien hochladen), werden sie beim ersten Gebrauch im Temp-Ordner zusammengesetzt.
import os, glob, tempfile, numpy as np, cv2

HERE = os.path.dirname(os.path.abspath(__file__))
T = 960
_SESS = [None]


def _teile(fn):
    return sorted(glob.glob(fn + '.[0-9][0-9]'))


def _modell():
    fn = os.environ.get('GS_YOLOMODELL', os.path.join(HERE, 'griffschrift_koepfe.onnx'))
    if os.path.exists(fn) or not _teile(fn): return fn
    ziel = os.path.join(tempfile.gettempdir(), 'griffschrift_koepfe.onnx')
    if not os.path.exists(ziel):
        tmp = ziel + '.%d' % os.getpid()
        with open(tmp, 'wb') as o:
            for t in _teile(fn): o.write(open(t, 'rb').read())
        os.replace(tmp, ziel)
    return ziel


def _sess():
    if _SESS[0] is None:
        import onnxruntime as ort
        fn = _modell()
        if not os.path.exists(fn): return None
        _SESS[0] = ort.InferenceSession(fn, providers=['CPUExecutionProvider'])
    return _SESS[0]


def vorhanden():
    fn = os.environ.get('GS_YOLOMODELL', os.path.join(HERE, 'griffschrift_koepfe.onnx'))
    return os.path.exists(fn) or bool(_teile(fn))


def _nms(b, s, thr):
    idx = cv2.dnn.NMSBoxes([[float(x), float(y), float(w), float(h)] for x, y, w, h in b], [float(v) for v in s], 0.0, thr)
    return [int(i) for i in np.array(idx).reshape(-1)]


def koepfe(g, conf=None, ueberlapp=192):
    """g: Graubild der Seite. Liefert Liste (cx, cy, w, h, klasse, sicherheit) in Bildpunkten von g."""
    sess = _sess()
    if sess is None: return []
    conf = float(os.environ.get('GS_YOLOCONF', '0.35')) if conf is None else conf
    H, W = g.shape[:2]; boxes = []; scores = []; cls = []
    xs = list(range(0, max(1, W - T) + 1, T - ueberlapp)); ys = list(range(0, max(1, H - T) + 1, T - ueberlapp))
    if xs[-1] + T < W: xs.append(W - T)
    if ys[-1] + T < H: ys.append(H - T)
    name = sess.get_inputs()[0].name
    dt = np.float16 if 'float16' in sess.get_inputs()[0].type else np.float32   # Modell kann halbgenau (fp16) ausgegeben sein
    for y0 in ys:
        for x0 in xs:
            crop = np.full((T, T), 255, np.uint8); c_ = g[max(0, y0):y0 + T, max(0, x0):x0 + T]; crop[:c_.shape[0], :c_.shape[1]] = c_
            inp = np.repeat((crop[None, None].astype(np.float32) / 255.0).astype(dt), 3, axis=1)
            out = sess.run(None, {name: inp})[0][0].astype(np.float32)   # (4+K, N)
            pr = out[4:]; k = pr.argmax(0); sc = pr.max(0); m = sc >= conf
            for (cx, cy, w, h), kk, ss in zip(out[:4, m].T, k[m], sc[m]):
                # Köpfe am Rand eines Ausschnitts nur nehmen, wenn sie nicht im Überlappungsstreifen liegen, den der Nachbar ganz sieht
                if (x0 > 0 and cx < ueberlapp / 2) or (y0 > 0 and cy < ueberlapp / 2): continue
                if (x0 + T < W and cx > T - ueberlapp / 2) or (y0 + T < H and cy > T - ueberlapp / 2): continue
                boxes.append((x0 + cx - w / 2, y0 + cy - h / 2, w, h)); scores.append(ss); cls.append(int(kk))
    if not boxes: return []
    keep = _nms(boxes, scores, 0.5)
    return [(boxes[i][0] + boxes[i][2] / 2, boxes[i][1] + boxes[i][3] / 2, boxes[i][2], boxes[i][3], cls[i], float(scores[i])) for i in keep]


def detect(g, st):
    """Köpfe als Einträge wie detect.detect; vorgesetzte Doppelkreuze machen den runden Kopf rechts daneben zum Kreuzkopf."""
    hs = koepfe(g)
    if not hs or not st: return []
    mids = [(ys[2], (ys[4] - ys[0]) / 4) for ys in st]
    out = []; pre = []
    for cx, cy, w, h, k, sc in hs:
        si = int(np.argmin([abs(cy - m) for m, _ in mids])); mid, d = mids[si]
        if abs(cy - mid) > 9 * d: continue
        if k == 6: pre.append((cx, cy, si, d)); continue
        t = 'o' if k in (0, 1) else ('x' if k in (2, 3) else 't')
        fill = k in (0, 2, 4)
        out.append(dict(x=float(cx), y=float(cy), p=int(round((mid + 0.5 - cy) * 2 / d)), t=t, fill=fill, sc=round(sc, 2), st=si, d=d, mid=mid, yolo=k))
    for px, py, si, d in pre:
        c = [h for h in out if h['st'] == si and h['t'] == 'o' and 0.7 * d < h['x'] - px < 2.4 * d and abs(h['y'] - py) < 0.5 * d]
        if c:
            h = min(c, key=lambda h: h['x'] - px); h['t'] = 'x'; h['xmark'] = True
    return out
