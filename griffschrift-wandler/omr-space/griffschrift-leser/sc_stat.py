# Colab: mittlere Sicherheit (sc) der Kopf-Erkennung je Blatt – Merkmal für „Handschrift / unsauberer Satz“ (nur Messung)
import os, sys, glob, numpy as np
M = '/content/drive/MyDrive/griffschrift_modell'; S = '/content/Harmonikabau/griffschrift-wandler/omr-space/gsleser'
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx'); sys.path.insert(0, S); os.chdir(S)
import gs2xml as X, yolo_koepfe as YK, fitz, tempfile
from fuchsread import druck_bars
G = '/content/gs'
pdfs = sorted(f for f in glob.glob(G + '/*.pdf') + glob.glob(G + '/*.PDF') if 'Bassbelegung' not in f and 'Ziachorgl' not in f and 'Posch' not in f and 'Hornsteiner' not in f)
for f in pdfs:
    sc = []; xs = []
    for i, pg in enumerate(fitz.open(f)):
        if i >= 2: break
        fn = tempfile.mktemp(suffix='.png'); pg.get_pixmap(dpi=200).save(fn)
        g = X.load_gray(fn)
        if g is None: continue
        st, _ = druck_bars(g); det = YK.detect(g, st)
        sc += [h['sc'] for h in det]; xs += [h['sc'] for h in det if h['t'] == 'x']
    if sc: print(f"{os.path.basename(f)[:40]} | Köpfe {len(sc)} mittel {np.mean(sc):.3f} unter0.5 {np.mean(np.array(sc)<0.5):.3f} Kreuze {len(xs)} Kreuz-mittel {np.mean(xs) if xs else 0:.3f}", flush=True)
