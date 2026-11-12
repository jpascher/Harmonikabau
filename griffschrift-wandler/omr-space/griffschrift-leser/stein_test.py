# Colab: Steinkellner-Blätter mit Schreibweise „steinkellner“ (Taktstriche über und unter dem System, v163) gegen „auto“; nur Messung
import os, sys, glob, shutil, functools
M = '/content/drive/MyDrive/griffschrift_modell'; L = '/content/Harmonikabau/griffschrift-wandler/omr-space/griffschrift-leser'
shutil.copy(M + '/gerust_fehler.py', L)
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx', GS_GCACHE='/content/gcache', GS_GERUEST='1')
os.chdir(L); sys.path.insert(0, L)
import gerust_fehler as F, gs2xml as X
orig = X.convert
for f in sorted(glob.glob('/content/gs/*Steinkellner*.PDF')):
    if 'Bassbelegung' in f: continue
    for stil in ('auto', 'steinkellner'):
        X.convert = functools.partial(orig, style=stil)
        r = F.lauf(f, True)
        echt = [o for o in r.get('offen', []) if o['art'] not in ('auftakt', 'schlusstakt', 'teiltakt_wiederholung')]
        print(os.path.basename(f)[:25], stil, 'Takte', r.get('takte'), 'offen', len(r.get('offen', [])), 'echt', len(echt), r.get('fehler', ''), [F.grob(o['art']) for o in echt][:12], flush=True)
