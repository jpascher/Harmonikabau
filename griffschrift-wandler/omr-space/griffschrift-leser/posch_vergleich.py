# Philipps Heft „Franz Posch“ gibt es in Griffschrift und in Normalnotation (nur zur Messung, beide nicht im Repository/Space/Lernen).
# Die Normalnotation wird mit homr gelesen und dient als Vergleich („zweites Soll“, selbst nicht fehlerfrei): Für jedes Stück werden die
# vom Griffschrift-Leser gelesenen Töne (Diskant) und Abstände der Einsätze mit der homr-Lesung der Normalnoten verglichen.
# Stimmung und Oktavlage werden aus dem ersten Stück bestimmt (beste Übereinstimmung), danach für alle Stücke festgehalten.
# Aufruf: python3 posch_vergleich.py ../gsleser GRIFFSCHRIFT.pdf NORMAL.pdf  [Arbeitsordner /content/posch]
import sys, os, glob, subprocess, collections, json
GS = sys.argv[1]; sys.path.insert(0, GS)
import fitz, numpy as np, xml.etree.ElementTree as ET
import gs2xml as X
from align import align
PDF_G, PDF_N = sys.argv[2], sys.argv[3]; W = '/content/posch'; os.makedirs(W, exist_ok=True)
# Stück: (Name, Seiten Griffschrift, Seiten Normalnotation) – Griffschrift Seite 1 ist die Belegung
STUECKE = [('A bissl Boarischer', 2), ('Der Martina zur Ehr', 4), ('Ein Hoch auf unsere Helga', 6), ('Gemmals an', 8), ('Gnadenlos', 10),
           ('Herbstlaub', 12), ('Hopp Heinz', 14), ('Marschl', 16), ('Ohrwürmchen Walzer', 18), ('Robert seiner', 20), ('Walzerle', 22), ('Zitterpartie', 24)]
STUECKE = [(n, [g, g + 1], [g - 1, g]) for n, g in STUECKE]

def seite(pdf, i, tag, dpi):
    fn = f'{W}/{tag}_{i:02d}.png'
    if not os.path.exists(fn): fitz.open(pdf)[i - 1].get_pixmap(dpi=dpi).save(fn)
    return fn

STEP = dict(C=0, D=2, E=4, F=5, G=7, A=9, B=11)
def homr_noten(png):
    """Normalnoten mit homr lesen -> Liste (Einsatz in Vierteln, Menge MIDI) der obersten Notenzeile je System (Melodie/Diskant)."""
    xml = png[:-4] + '.musicxml'
    if not os.path.exists(xml): subprocess.run(['homr', png], capture_output=True)
    if not os.path.exists(xml): return []
    root = ET.parse(xml).getroot(); ev = collections.defaultdict(set)
    t0 = 0.0
    for part in root.iter('part'):
        div = 1; t = 0.0
        for meas in part.iter('measure'):
            for el in meas:
                if el.tag == 'attributes' and el.find('divisions') is not None: div = float(el.find('divisions').text)
                elif el.tag == 'backup': t -= float(el.find('duration').text) / div
                elif el.tag == 'forward': t += float(el.find('duration').text) / div
                elif el.tag == 'note':
                    d = float(el.find('duration').text) / div if el.find('duration') is not None else 0
                    ch = el.find('chord') is not None; st = el.find('staff'); st = int(st.text) if st is not None else 1
                    if ch: t -= last_d
                    if el.find('pitch') is not None and st == 1:
                        p = el.find('pitch'); m = 12 * (int(p.find('octave').text) + 1) + STEP[p.find('step').text] + int(float(p.find('alter').text)) if p.find('alter') is not None else 12 * (int(p.find('octave').text) + 1) + STEP[p.find('step').text]
                        ev[round(t0 + t, 3)].add(m)
                    t += d; last_d = d
        t0 += t + 100   # Teile (Systeme) hintereinander
    return sorted(ev.items())

def griffschrift(pngs, shift):
    out = []; t = 0
    for p in pngs:
        meas, L = X.read(p, shift)
        for m in meas:
            for e in m:
                if not e['rest']: out.append((t / 12.0, {h[2] for h in e['heads']}))
                t += e['u']
    return out

def vergleich(G, N, okt):
    A = [{p + okt for p in s} for _, s in G]; B = [s for _, s in N]
    pr = align(A, B); gleich = sum(A[i] == B[j] for i, j in pr); oben = sum(max(A[i]) == max(B[j]) for i, j in pr)
    ver = []
    for (i, j), (i2, j2) in zip(pr, pr[1:]):
        if i2 == i + 1 and j2 == j + 1:
            dl = G[i2][0] - G[i][0]; dp = N[j2][0] - N[j][0]
            if dl > 0 and 0 < dp < 50: ver.append(round(dl / dp, 2))
    rok = sum(v == 1.0 for v in ver)
    return dict(n=len(B), g=len(A), paar=len(pr), gleich=gleich, oben=oben, rpaare=len(ver), rok=rok)

N_all = {}
for name, gs_, ns in STUECKE:
    N_all[name] = [x for i in ns for x in homr_noten(seite(PDF_N, i, 'n', 300))]
# Stimmung/Oktave am ersten Stück bestimmen
name, gs_, ns = STUECKE[0]; G0 = [seite(PDF_G, i, 'g', 200) for i in gs_]
best = None
for sh in (0, -2, 2, -4, 4, -1, 1):
    try: Gs = griffschrift(G0, sh)
    except Exception: continue
    for okt in (-12, 0, 12):
        r = vergleich(Gs, N_all[name], okt)
        if best is None or r['oben'] > best[0]: best = (r['oben'], sh, okt)
_, SHIFT, OKT = best; print('Stimmung shift', SHIFT, 'Oktave', OKT, flush=True)
ARTEN = (('leser', dict(GS_GERUEST='0')), ('v160', dict(GS_GERUEST='1')))
T = {a: collections.Counter() for a, _ in ARTEN}
for name, gs_, ns in STUECKE:
    zeile = [name]
    for a, env in ARTEN:
        X.GERUEST = env['GS_GERUEST'] == '1'
        try: G = griffschrift([seite(PDF_G, i, 'g', 200) for i in gs_], SHIFT)
        except Exception as e: zeile.append(f'{a}: Fehler {str(e)[:40]}'); continue
        r = vergleich(G, N_all[name], OKT); T[a].update(r)
        zeile.append(f"{a}: Akkorde gleich {r['gleich']}/{r['n']} oben gleich {r['oben']}/{r['n']} Rhythmus {r['rok']}/{r['rpaare']}")
    print(' | '.join(zeile), flush=True)
for a, _ in ARTEN:
    s = T[a]
    print(f"{a}: Normalnoten-Akkorde {s['n']}, gelesen {s['g']}, ganz gleich {s['gleich']} ({s['gleich']/max(1,s['n']):.3f}), oberster Ton gleich {s['oben']} ({s['oben']/max(1,s['n']):.3f}), Rhythmus {s['rok']}/{s['rpaare']} ({s['rok']/max(1,s['rpaare']):.3f})")
