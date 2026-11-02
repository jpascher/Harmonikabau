# Arbeitsordner /tmp/philipp wie in tri_data.py beschrieben; Aufruf: python3 philipp_pruefen.py ../gsleser [stück,stück]
# Philipps Heft lesen und mit seinen bekannten Griffen (aus dem Wandler) vergleichen
import sys, json, os, collections, numpy as np, cv2
GS = sys.argv[1] if len(sys.argv) > 1 else '../gsleser'
sys.path.insert(0, GS)
import gs2xml as X, detect as D
from align import chords_det, align
from assemble import button
from fuchsread import druck_bars
PAGES = {'griffiglandler': [6, 7], 'furnvoter': [8, 9], 'polkaroyal': [10, 11], 'zommgschwoasste': [12, 13],
         'eiswind': [14, 15, 16], 'isabellile': [17], 'zehnfuenfer': [18, 19, 20], 'polarweis': [21],
         'boehmischerhauch': [22, 23, 24], 'bobblpolka': [25], 'sibirischegedanken': [26, 27, 28, 29]}
OFF = {1: -1, 2: 0, 3: 1, 4: 2}
PINS = json.load(open('/tmp/philipp/ph_pins.json'))
def pin_chords(n):
    m = {}
    for t, pi, r, num, dr in n: m.setdefault(round(t, 3), []).append((r, num, dr))
    return [m[t] for t in sorted(m)]
def read_piece(k, pagedir='/tmp/philipp/pages', cache=True):
    out = []
    for p in PAGES[k]:
        f = f'{pagedir}/p-{p:02d}.png'
        g = X.load_gray(f); w = f'/tmp/philipp/norm/p-{p:02d}.png'; os.makedirs(os.path.dirname(w), exist_ok=True); cv2.imwrite(w, g)
        st, bars = druck_bars(g); det = D.detect(w, 1); C = chords_det(det)
        for c in C:
            si = c[0]['st']; x = float(np.median([h['x'] for h in c]))
            dr = 'D' if any(a - 3 <= x <= b + 3 for a, b, _ in bars.get(si, [])) else 'Z'
            out.append(dict(page=p, st=si, x=x, dir=dr, heads=[(h['p'], h['t'], button(h['p'], h['t'])) for h in c]))
    return out
def compare(k, res):
    B = pin_chords(PINS[k])
    A = [{(r, n) for _, _, (r, n) in c['heads']} for c in res]
    Bs = [{(r, n) for r, n, _ in b} for b in B]
    pr = align([{(2 * (n - 6) + OFF[r], 'o' if r < 3 else 'x') for r, n in a} for a in A],
               [{(2 * (n - 6) + OFF[r], 'o' if r < 3 else 'x') for r, n in b} for b in Bs])
    st = collections.Counter()
    for i, j in pr:
        a, b = A[i], Bs[j]; st['paar'] += 1; st['gleich'] += a == b
        bd = collections.Counter(d for _, _, d in B[j]).most_common(1)[0][0]; st['balg_ok'] += res[i]['dir'] == bd
        for (r, n) in b:
            st[f'soll_r{r}'] += 1; st[f'ok_r{r}'] += (r, n) in a
        for (r, n) in a - b: st[f'falsch_r{r}'] += 1
    st['bild'] = len(A); st['pins'] = len(B)
    return st, pr
if __name__ == '__main__':
    T = collections.Counter(); allres = {}
    for k in (sys.argv[2].split(',') if len(sys.argv) > 2 else PAGES):
        res = read_piece(k); st, pr = compare(k, res); T.update(st); allres[k] = dict(res=res, pairs=pr)
        print(k, dict(bild=st['bild'], pins=st['pins'], paar=st['paar'], gleich=st['gleich'], balg=st['balg_ok']), flush=True)
    print('GESAMT Akkorde gepaart', T['paar'], 'ganz gleich', T['gleich'], f"{T['gleich']/max(1,T['paar']):.3f}", 'Balg', f"{T['balg_ok']/max(1,T['paar']):.3f}")
    for r in (1, 2, 3, 4): print(f'Reihe {r}: soll {T[f"soll_r{r}"]}, gefunden {T[f"ok_r{r}"]} ({T[f"ok_r{r}"]/max(1,T[f"soll_r{r}"]):.3f}), falsch gelesen {T[f"falsch_r{r}"]}')
    json.dump(allres, open('/tmp/philipp/phcheck_out.json', 'w'), default=lambda o: o.item() if hasattr(o, 'item') else str(o))
