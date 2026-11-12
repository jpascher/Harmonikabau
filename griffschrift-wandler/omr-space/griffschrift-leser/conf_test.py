# Colab: Bummelpetrus u. a. mit verschiedenen Sicherheitsschwellen des Kopf-Erkenners (GS_YOLOCONF): Akkorde, Takte, offene Takte
import os, sys, glob, importlib
M = '/content/drive/MyDrive/griffschrift_modell'; L = '/content/Harmonikabau/griffschrift-wandler/omr-space/griffschrift-leser'
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx', GS_GCACHE='/content/gcache', GS_GERUEST='1')
os.chdir(L); sys.path.insert(0, L)
import gerust_fehler as F, gs2xml as X
for name in ('Bummelpetrus', 'An Bohm Luki', 'Frei und Treu'):
    f = [x for x in glob.glob('/content/gs/*.pdf') if name in x][0]
    for c in ('0.35', '0.25', '0.15', '0.08'):
        os.environ['GS_YOLOCONF'] = c
        n = [0]; r0 = F._orig_read
        def rd(p, *a, **k):
            m, L_ = r0(p, *a, **k); n[0] += sum(1 for t in m for e in t if not e['rest']); return m, L_
        F._orig_read = rd
        r = F.lauf(f, True); F._orig_read = r0
        echt = [o for o in r.get('offen', []) if o['art'] not in ('auftakt', 'schlusstakt', 'teiltakt_wiederholung')]
        print(os.path.basename(f)[:24], 'conf', c, 'Akkorde', n[0], 'Takte', r.get('takte'), 'offen', len(r.get('offen', [])), 'echt', len(echt), flush=True)
