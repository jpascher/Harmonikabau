# Colab (CPU genügt): Prüfung des Zwischenschritts „Dreiecke für homr in runde Köpfe wandeln“ (gs2xml.tri_rund, v159).
# Gelesen werden alle Blätter, die im Messlauf fehler_geruest.json als Dreieck-Seite eingestuft wurden, und Philipps Heft
# (griffschrift_daten/Ziachorgl_Sticklan_Philipp.pdf, alle Seiten; nur zum Messen, nicht ins Repo, nicht in den Space, nicht ins Lernen).
# Je Datei drei Lesungen: ohne Gerüst | Gerüst wie bisher (Dreieck-Seiten ohne homr, GS_GERUESTTRI=0) | Gerüst mit runden Köpfen.
# Braucht gerust_fehler.py im Drive-Ordner griffschrift_modell. Aufruf in einer Zelle, nachdem Drive verbunden ist, als eigener Prozess
# (im Notizbuch-Prozess scheitert homr, wenn dort schon eine andere Pillow-Fassung geladen ist):
#   !cd /content && python3 -u /content/Harmonikabau/griffschrift-wandler/omr-space/griffschrift-leser/tri_rund_colab.py
# Ergebnis: griffschrift_modell/tri_rund3.json (je Datei) und tri_rund3.txt (Zusammenfassung).
import os, sys, glob, json, shutil, subprocess, tempfile, collections
sh = lambda c: subprocess.run(c, shell=True)
M = '/content/drive/MyDrive/griffschrift_modell'; G = '/content/gs'
L = '/content/Harmonikabau/griffschrift-wandler/omr-space/griffschrift-leser'
os.chdir('/content')
if not os.path.isdir('/content/Harmonikabau'): sh('git clone -q --depth 1 https://github.com/jpascher/Harmonikabau.git')
try: import homr  # noqa
except Exception: sh('pip -q install homr onnxruntime scikit-learn==1.9.1 mido pymupdf pytesseract opencv-python-headless')
shutil.copy(M + '/gerust_fehler.py', L)
if os.path.exists(M + '/gs2xml_v159.py') and 'tri_rund' not in open(L + '/../gsleser/gs2xml.py').read():   # nur vor dem Push von v159 nötig
    shutil.copy(M + '/gs2xml_v159.py', L + '/../gsleser/gs2xml.py')
for _m in ('gs2xml', 'gerust', 'gerust_fehler'): sys.modules.pop(_m, None)   # frisch laden, falls früher schon importiert
if not os.path.isdir(G):
    os.makedirs(G + '/zip', exist_ok=True)
    sh(f'cp -rn /content/drive/MyDrive/griffschrift_daten/. {G}/')
    for z in sorted(glob.glob(G + '/*.zip')): subprocess.run(['unzip', '-q', '-o', z, '-d', G + '/zip/' + os.path.basename(z)[:-4]])
else:
    sh(f'cp -n /content/drive/MyDrive/griffschrift_daten/*.pdf {G}/')
os.environ.update(GS_YOLO='1', GS_YOLOMODELL=M + '/griffschrift_koepfe2.onnx', GS_GCACHE='/content/gcache', GS_GERUEST='1')
os.chdir(L); sys.path.insert(0, L)
import gerust_fehler as F
import gs2xml as X
assert hasattr(X, 'tri_rund'), 'gs2xml ohne tri_rund (gs2xml_v159.py fehlt auf Drive)'

# Dateien: Dreieck-Seiten aus dem Messlauf + Philipps Heft
alt = json.load(open(M + '/fehler_geruest.json'))
alle = {os.path.basename(f): f for f in glob.glob(G + '/**/*', recursive=True) if os.path.isfile(f)}
liste = []
for k, v in alt.items():
    q = (v.get('mit') or {}).get('quellen') or {}
    if q.get('leser_dreieckseite') and v['datei'] in alle: liste.append((v['gruppe'], alle[v['datei']]))
ph = G + '/Ziachorgl_Sticklan_Philipp.pdf'
if os.path.exists(ph):   # Philipps Heft Seite für Seite (Titel-/Textseiten ohne Notenlinien werden übersprungen)
    import fitz
    _pd = '/content/philipp_seiten'; os.makedirs(_pd, exist_ok=True)
    for i, pg in enumerate(fitz.open(ph)):
        fn = f'{_pd}/philipp_s{i:03d}.png'
        if not os.path.exists(fn): pg.get_pixmap(dpi=200).save(fn)
        liste.append(('Philipp', fn))
print('Dateien:', len(liste), flush=True)

_seiten4 = F.seiten
def seiten_alle(p):   # Philipps Heft: alle Seiten (sonst wie im Space nur die ersten vier)
    if 'Ziachorgl' not in p: return _seiten4(p)
    import fitz
    tmp = tempfile.mkdtemp(); out = []
    for i, pg in enumerate(fitz.open(p)):
        fn = os.path.join(tmp, f's{i:03d}.png'); pg.get_pixmap(dpi=200).save(fn); out.append(fn)
    return out
F.seiten = seiten_alle

_ST = []
_r0 = F._orig_read   # ursprüngliches gs2xml.read (gerust_fehler hat X.read schon umhüllt)
def _rd(path, *a, **k):
    m, L = _r0(path, *a, **k); _ST.append(dict(X._DBG.get('gerust') or {}, tri=bool(X.D.LAST.get('tri_page')))); return m, L
F._orig_read = _rd
def stat():
    c = collections.Counter()
    for d in _ST:
        for k, v in d.items():
            if isinstance(v, bool): c[k] += int(v)
            elif isinstance(v, (int, float)): c[k] += v
            elif k == 'fehler': c['fehler:' + str(v)[:40]] += 1
    _ST.clear(); return dict(c)
def kurz(r):
    if 'fehler' in r: return dict(fehler=r['fehler'])
    o = r['offen']; echt = [x for x in o if x['art'] not in ('auftakt', 'schlusstakt', 'teiltakt_wiederholung')]
    return dict(takte=r['takte'], offen=len(o), echt=len(echt), quellen=r['quellen'], arten=dict(collections.Counter(F.grob(x['art']) for x in o)))

AUS = M + '/tri_rund3.json'; erg = json.load(open(AUS)) if os.path.exists(AUS) else {}
for gr_, p in liste:
    key = gr_ + '|' + os.path.basename(p)
    if key in erg: continue
    r = {}
    os.environ['GS_GERUESTTRI'] = '0'; r['bisher'] = kurz(F.lauf(p, True))
    os.environ['GS_GERUESTTRI'] = '1'; _ST.clear(); r['rund'] = kurz(F.lauf(p, True)); r['rund']['homrstat'] = stat()
    r['ohne'] = kurz(F.lauf(p, False))
    if all('fehler' in r[k] for k in r): print(key, 'übersprungen:', r['rund']['fehler'], flush=True); continue
    erg[key] = r; json.dump(erg, open(AUS, 'w'), ensure_ascii=False)
    f = lambda d: f"{d.get('offen','F')}/{d.get('echt','F')}"
    print(f"{key}: Takte {r['rund'].get('takte')}  offen/echt  ohne {f(r['ohne'])}  bisher {f(r['bisher'])}  rund {f(r['rund'])}  Quellen rund {r['rund'].get('quellen')}  homr {r['rund'].get('homrstat')}", flush=True)

zs = ['', '== Summe (offen / davon echt offen) ==']
for g in sorted({k.split('|')[0] for k in erg}) + ['alle']:
    V = [v for k, v in erg.items() if g == 'alle' or k.split('|')[0] == g]
    s = lambda n, f: sum(v[n].get(f, 0) or 0 for v in V)
    zs.append(f"{g}: {len(V)} Dateien, Takte {s('rund','takte')}  ohne {s('ohne','offen')}/{s('ohne','echt')}  bisher {s('bisher','offen')}/{s('bisher','echt')}  rund {s('rund','offen')}/{s('rund','echt')}")
print('\n'.join(zs), flush=True)
open(M + '/tri_rund3.txt', 'w').write('\n'.join(f"{k}: {json.dumps(v, ensure_ascii=False)}" for k, v in erg.items()) + '\n' + '\n'.join(zs) + '\n')
print('fertig')
