# Aufruf: python3 philipp_pruefen_gs2xml.py ../gsleser [stück,stück]  (Arbeitsordner /tmp/philipp)
# Ende-zu-Ende-Prüfung mit gs2xml (wie im Space): Philipps Seiten lesen, Ereignisse (Akkorde) mit seinen Griffen vergleichen
import sys, json, collections, os
GS = sys.argv[1]; sys.path.insert(0, GS)
import gs2xml as X
from align import align
sys.path.insert(0, '/tmp/philipp')
from philipp_pruefen import PAGES, PINS, pin_chords
OFF = {1: -1, 2: 0, 3: 1, 4: 2}
D450 = {"druck": {1: [98,30,25,4,7,12,16,19,24,28,31,36,40,43,32,98], 2: [98,98,20,5,9,12,17,21,24,29,33,36,41,45,48,98,98], 3: [98,98,8,10,14,17,22,26,29,34,38,41,46,50,98,98], 4: [98,98,13,15,19,22,27,31,34,39,43,46,51,98,98]},
        "zug": {1: [98,25,30,7,11,14,17,21,23,26,29,33,35,38,41,98], 2: [98,98,18,10,12,16,19,22,24,28,31,34,38,40,43,98,98], 3: [98,98,9,15,17,21,24,27,29,33,36,39,43,45,98,98], 4: [98,98,14,20,22,26,29,32,34,38,41,44,48,98,98]}}
def key(r, n): return (2 * (n - 6) + OFF[r], 'o' if r < 3 else ('x' if r == 3 else 't'))
if __name__ == '__main__':
    pieces = sys.argv[2].split(',') if len(sys.argv) > 2 else list(PAGES)
    table = None if os.environ.get('T446') else D450
    T = collections.Counter()
    for k in pieces:
        ev = []
        for p in PAGES[k]:
            meas, L = X.read(f'/tmp/philipp/pages/p-{p:02d}.png', -2, table=table) if table else X.read(f'/tmp/philipp/pages/p-{p:02d}.png', -2)
            for m in meas:
                for e in m:
                    if not e['rest']: ev.append(e)
        A = [{(r, n) for r, n, _ in e['heads']} for e in ev]
        B = pin_chords(PINS[k]); Bs = [{(r, n) for r, n, _ in b} for b in B]
        pr = align([{key(r, n) for r, n in a} for a in A], [{key(r, n) for r, n in b} for b in Bs])
        st = collections.Counter(); st['pins'] = len(B); st['bild'] = len(A)
        # Tonhöhen: gelesene Töne gegen die bekannten (G-C-F-B)
        for i, j in pr:
            st['paar'] += 1; st['gleich'] += A[i] == Bs[j]
            bd = collections.Counter(d for _, _, d in B[j]).most_common(1)[0][0]; st['balg'] += ev[i]['dir'] == bd
            pa = sorted(p for _, _, p in ev[i]['heads']); pb = sorted(pi for _, pi, _, _, _ in [x for x in PINS[k] if False])
            for (r, n) in Bs[j]: st[f's{r}'] += 1; st[f'g{r}'] += (r, n) in A[i]
            for (r, n) in A[i] - Bs[j]: st[f'f{r}'] += 1
        T.update(st)
        print(k, dict(pins=st['pins'], bild=st['bild'], paar=st['paar'], gleich=st['gleich'], balg=st['balg']), flush=True)
    print('GESAMT gepaart', T['paar'], 'von', T['pins'], 'ganz gleich', T['gleich'], f"{T['gleich']/max(1,T['paar']):.3f}", 'Balg', f"{T['balg']/max(1,T['paar']):.3f}")
    for r in (1, 2, 3, 4): print(f'Reihe {r}: soll {T[f"s{r}"]}, gefunden {T[f"g{r}"]} ({T[f"g{r}"]/max(1,T[f"s{r}"]):.3f}), überzählig {T[f"f{r}"]}')
