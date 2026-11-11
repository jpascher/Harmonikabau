# Colab: Philipps Heft „Hansi Hornsteiner – Bergwelten“ Seite für Seite durch den Leser (v161, nur Messung): Takte, offene Takte, Dreieck-Seite
import os, sys, glob, json, shutil, collections, fitz
M = '/content/drive/MyDrive/griffschrift_modell'; L = '/content/Harmonikabau/griffschrift-wandler/omr-space/griffschrift-leser'
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx', GS_GCACHE='/content/gcache', GS_GERUEST='1')
shutil.copy(M + '/gerust_fehler.py', L); os.chdir(L); sys.path.insert(0, L)
import gerust_fehler as F, gs2xml as X
W = '/content/hornsteiner'; os.makedirs(W, exist_ok=True)
d = fitz.open('/content/drive/MyDrive/griffschrift_daten/Hansi Hornsteiner Bergwelten.pdf')
S = collections.Counter()
for i, p in enumerate(d):
    fn = f'{W}/h{i+1:02d}.png'
    if not os.path.exists(fn): p.get_pixmap(dpi=200).save(fn)
    r = F.lauf(fn, True)
    if 'fehler' in r: print(i + 1, 'ohne Noten', flush=True); continue
    echt = [o for o in r['offen'] if o['art'] not in ('auftakt', 'schlusstakt', 'teiltakt_wiederholung')]
    S['seiten'] += 1; S['takte'] += r['takte']; S['offen'] += len(r['offen']); S['echt'] += len(echt)
    print(f"S{i+1}: Takte {r['takte']} offen {len(r['offen'])} echt {len(echt)} Dreieck {bool(X.D.LAST.get('tri_page'))} Quellen {r['quellen']} Arten {dict(collections.Counter(F.grob(o['art']) for o in echt))}", flush=True)
print('GESAMT', dict(S))
