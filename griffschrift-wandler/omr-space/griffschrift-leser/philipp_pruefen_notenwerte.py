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
    T = collections.Counter(); CONF = collections.Counter(); WR = []; EX = collections.Counter()
    for k in pieces:
        ev = []
        for p in PAGES[k]:
            meas, L = X.read(f'/tmp/philipp/pages/p-{p:02d}.png', -2, table=table) if table else X.read(f'/tmp/philipp/pages/p-{p:02d}.png', -2)
            for m in meas:
                for e in m:
                    if not e['rest']: ev.append(e)
        A = [{(r, n) for r, n, _ in e['heads']} for e in ev]
        B = pin_chords(PINS[k]); Bs = [{(r, n) for r, n, _ in b} for b in B]; TS = sorted({round(t_, 3) for t_, *_ in PINS[k]})
        pr = align([{key(r, n) for r, n in a} for a in A], [{key(r, n) for r, n in b} for b in Bs])
        st = collections.Counter(); st['pins'] = len(B); st['bild'] = len(A)
        # Tonhöhen: gelesene Töne gegen die bekannten (G-C-F-B)
        for i, j in pr:
            st['paar'] += 1
            if j+1 < len(TS): st['dsoll'] += 1; st['dok'] += abs(ev[i]['u']/12 - (TS[j+1]-TS[j])) < 0.02
            if j+1 < len(TS): CONF[(round(TS[j+1]-TS[j],3), ev[i]['u']/12)] += 1
            st['x']+=0; st['gleich'] += A[i] == Bs[j]
            bd = collections.Counter(d for _, _, d in B[j]).most_common(1)[0][0]; st['balg'] += ev[i]['dir'] == bd
            pa = sorted(p for _, _, p in ev[i]['heads']); pb = sorted(pi for _, pi, _, _, _ in [x for x in PINS[k] if False])
            pin_p = {pi for t_, pi, r_, n_, d_ in PINS[k] if round(t_, 3) == TS[j]}; got_p = {p for _, _, p in ev[i]['heads']}
            st['tsoll'] += len(pin_p); st['tok'] += len(pin_p & got_p); st['tfalsch'] += len(got_p - pin_p)
            if got_p - pin_p and len(WR) < 40: WR.append((k, sorted(got_p), sorted(pin_p), ev[i].get('dir'), sorted(A[i]), sorted(Bs[j])))
            for (r, n) in Bs[j]: st[f's{r}'] += 1; st[f'g{r}'] += (r, n) in A[i]
            for (r, n) in A[i] - Bs[j]: st[f'f{r}'] += 1; EX[(r, n, 'leer' if not (A[i] & Bs[j]) else 'mit')] += 1
        T.update(st)
        print(k, 'dauer',st['dok'],st['dsoll'], dict(pins=st['pins'], bild=st['bild'], paar=st['paar'], gleich=st['gleich'], balg=st['balg']), flush=True)
    print('GESAMT gepaart', T['paar'], 'von', T['pins'], 'ganz gleich', T['gleich'], f"{T['gleich']/max(1,T['paar']):.3f}", 'Balg', f"{T['balg']/max(1,T['paar']):.3f}")
    print(EX.most_common(25))
    print(sorted(CONF.items(), key=lambda x: -x[1])[:30]);print('Dauer richtig', T['dok'], 'von', T['dsoll']);print('Töne richtig', T['tok'], 'von', T['tsoll'], f"{T['tok']/max(1,T['tsoll']):.3f}", 'falsch', T['tfalsch'])
    for r in (1, 2, 3, 4): print(f'Reihe {r}: soll {T[f"s{r}"]}, gefunden {T[f"g{r}"]} ({T[f"g{r}"]/max(1,T[f"s{r}"]):.3f}), überzählig {T[f"f{r}"]}')
