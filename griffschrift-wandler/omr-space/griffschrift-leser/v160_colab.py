# Colab (CPU genügt): Messung v160 gegen den Messlauf vom 11.11.2026 (fehler_geruest.json, Gerüst aus homr, Stand v158).
# v160: strengere Dreieck-Einstufung (gs2xml.tri_seite), Dreiecke für homr rund (v159), überzählige homr-Pausen weg (gerust.py).
# Liest dieselben Blätter wie gerust_fehler_colab.py (PDFs, Einsendungen, Fotos) und Philipps Heft Seite für Seite (nur Messung).
# Braucht gerust_fehler.py im Drive-Ordner griffschrift_modell; homr-Lesungen aus /content/gcache werden wiederverwendet.
# Aufruf als eigener Prozess:  !cd /content && python3 -u .../griffschrift-leser/v160_colab.py
# Versuch Kreuze rund: !cd /content && GS_GKREUZRUND=1 V160_NAME=v160_kreuz V160_GRUPPEN=PDFs,Einsendungen python3 -u …
import os, sys, glob, json, shutil, subprocess, collections
M = '/content/drive/MyDrive/griffschrift_modell'; G = '/content/gs'
L = '/content/Harmonikabau/griffschrift-wandler/omr-space/griffschrift-leser'; S = L + '/../gsleser'
for n in ('gs2xml', 'gerust'):   # vor dem Push: neue Fassungen von Drive (v160k vor v160)
    for v in ('v160k', 'v160'):
        if os.path.exists(f'{M}/{n}_{v}.py'): shutil.copy(f'{M}/{n}_{v}.py', f'{S}/{n}.py'); break
shutil.copy(M + '/gerust_fehler.py', L)
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx', GS_GCACHE='/content/gcache', GS_GERUEST='1')
os.chdir(L); sys.path.insert(0, L)
import gerust_fehler as F, gs2xml as X
assert hasattr(X, 'tri_seite')

fotos = {}
for f in sorted(glob.glob(G + '/zip/attachments*/**/*.jpg', recursive=True)): fotos.setdefault(os.path.basename(f), f)
fotos = [fotos[k] for k in sorted(fotos)] + ([G + '/20261103_165454.jpg'] if os.path.exists(G + '/20261103_165454.jpg') else [])
pdfs = sorted(f for f in glob.glob(G + '/*.pdf') + glob.glob(G + '/*.PDF') if 'Bassbelegung' not in f and 'Ziachorgl' not in f)
eins = sorted(f for f in glob.glob(G + '/zip/*/**/*.pdf', recursive=True) if '/attachments' not in f)
ph = sorted(glob.glob('/content/philipp_seiten/philipp_s0*.png'))
liste = [('PDFs', f) for f in pdfs] + [('Einsendungen', f) for f in eins] + [('Fotos', f) for f in fotos] + [('Philipp', f) for f in ph]
_gr = os.environ.get('V160_GRUPPEN')   # z. B. PDFs,Einsendungen
if _gr: liste = [(g, f) for g, f in liste if g in _gr.split(',')]
alt = json.load(open(M + '/fehler_geruest.json')); altP = {}
if os.path.exists(M + '/tri_rund3.json'):
    for k, v in json.load(open(M + '/tri_rund3.json')).items(): altP[k] = v.get('rund', {})
print('Dateien', len(liste), flush=True)

_ST = []; _r0 = F._orig_read
def _rd(path, *a, **k):
    m, L_ = _r0(path, *a, **k); _ST.append(dict(X._DBG.get('gerust') or {}, tri=bool(X.D.LAST.get('tri_page')))); return m, L_
F._orig_read = _rd
def stat():
    c = collections.Counter()
    for d in _ST:
        for k, v in d.items():
            if isinstance(v, bool): c[k] += int(v)
            elif isinstance(v, (int, float)): c[k] += v
    _ST.clear(); return dict(c)
ECHT = lambda o: o['art'] not in ('auftakt', 'schlusstakt', 'teiltakt_wiederholung')

NAME = os.environ.get('V160_NAME', 'v160'); AUS = M + '/' + NAME + '.json'; erg = json.load(open(AUS)) if os.path.exists(AUS) else {}
for gr_, p in liste:
    key = gr_ + '|' + os.path.basename(p)
    if key in erg: continue
    _ST.clear(); r = F.lauf(p, True)
    if 'fehler' in r: print(key, 'Fehler', r['fehler'], flush=True); continue
    a = alt.get(key, {}).get('mit') if gr_ != 'Philipp' else None
    vorher = (len(a['offen']), len([o for o in a['offen'] if ECHT(o)])) if a and 'offen' in a else ((altP[key].get('offen'), altP[key].get('echt')) if key in altP else (None, None))
    erg[key] = dict(gruppe=gr_, takte=r['takte'], offen=len(r['offen']), echt=len([o for o in r['offen'] if ECHT(o)]), vorher=vorher,
                    quellen=r['quellen'], homr=stat(), arten=dict(collections.Counter(F.grob(o['art']) for o in r['offen'] if ECHT(o))))
    json.dump(erg, open(AUS, 'w'), ensure_ascii=False)
    e = erg[key]; print(f"{key}: Takte {e['takte']} offen/echt {e['offen']}/{e['echt']} vorher {vorher[0]}/{vorher[1]} homr {e['homr']}", flush=True)

zs = ['', '== Summe offen/echt: v160 gegen vorher (v158 bzw. bei Philipp v159) ==']
for g in ('PDFs', 'Einsendungen', 'Fotos', 'Philipp'):
    V = [v for v in erg.values() if v['gruppe'] == g and v['vorher'][0] is not None]
    zs.append(f"{g}: {len(V)} Dateien, Takte {sum(v['takte'] for v in V)}, jetzt {sum(v['offen'] for v in V)}/{sum(v['echt'] for v in V)}, vorher {sum(v['vorher'][0] for v in V)}/{sum(v['vorher'][1] for v in V)}, Pausen weg {sum(v['homr'].get('pausen_weg', 0) for v in V)}, Dreieck-Seiten {sum(v['homr'].get('tri', 0) for v in V)}")
print('\n'.join(zs), flush=True)
open(M + '/' + NAME + '.txt', 'w').write('\n'.join(f"{k}: {json.dumps(v, ensure_ascii=False)}" for k, v in erg.items()) + '\n' + '\n'.join(zs) + '\n')
print('fertig')
