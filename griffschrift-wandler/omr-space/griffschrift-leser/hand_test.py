# Colab: Bummelpetrus mit Schreibweise auto / handschrift (v169; nur Messung, Blatt nicht im Repository): Akkorde, Köpfe je Form, Druck/Zug, Takte (nur Messung)
import os, sys, glob, shutil, collections, fitz, tempfile
M = '/content/drive/MyDrive/griffschrift_modell'; S = '/content/Harmonikabau/griffschrift-wandler/omr-space/gsleser'
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx', GS_GCACHE='/content/gcache', GS_GERUEST='1')
sys.path.insert(0, S); os.chdir(S)
import gs2xml as X
f = glob.glob('/content/gs/Bummelpetrus*.pdf')[0]
fn = tempfile.mktemp(suffix='.png'); fitz.open(f)[0].get_pixmap(dpi=200).save(fn)
for stil in ('auto', 'handschrift'):
    X.set_style(stil); X.GERUEST = True
    m, L = X.read(fn, 0)
    ev = [e for t in m for e in t if not e['rest']]
    print(stil, 'Takte', len(m), 'Akkorde', len(ev), 'Druck', sum(e['dir'] == 'D' for e in ev), 'Zug', sum(e['dir'] == 'Z' for e in ev),
          'Reihen', dict(collections.Counter(h[0] for e in ev for h in e['heads'])), 'Erkenner', X.D.LAST.get('erkenner'), 'YOLO-Köpfe', X.D.LAST.get('yolo'), flush=True)
