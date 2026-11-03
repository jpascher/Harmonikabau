# Prüfung an Philipps Heft (Arbeitsordner /tmp/philipp wie in tri_data.py): Köpfe je Reihe UND Tonhöhen (450 Öllerer/Schaborak, G-C-F-B)
# Aufruf: python3 philipp_pruefen_toene.py ../gsleser stück,stück  [details]
import sys, json, collections, numpy as np
GS = sys.argv[1]; sys.path.insert(0, GS)
import importlib
from philipp_pruefen import PAGES, PINS, pin_chords, read_piece
from align import align
from assemble import pitch
OFF = {1: -1, 2: 0, 3: 1, 4: 2}
D450 = {"druck": {1: [98,30,25,4,7,12,16,19,24,28,31,36,40,43,32,98], 2: [98,98,20,5,9,12,17,21,24,29,33,36,41,45,48,98,98], 3: [98,98,8,10,14,17,22,26,29,34,38,41,46,50,98,98], 4: [98,98,13,15,19,22,27,31,34,39,43,46,51,98,98]},
        "zug": {1: [98,25,30,7,11,14,17,21,23,26,29,33,35,38,41,98], 2: [98,98,18,10,12,16,19,22,24,28,31,34,38,40,43,98,98], 3: [98,98,9,15,17,21,24,27,29,33,36,39,43,45,98,98], 4: [98,98,14,20,22,26,29,32,34,38,41,44,48,98,98]}}
def key(r, n): return (2 * (n - 6) + OFF[r], 'o' if r < 3 else ('x' if r == 3 else 't'))
def run(pieces, detail=False):
    T = collections.Counter(); miss = []; wrong = []
    for k in pieces:
        res = read_piece(k); B = pin_chords(PINS[k])
        A = [{(r, n) for _, _, (r, n) in c['heads']} for c in res]; Bs = [{(r, n) for r, n, _ in b} for b in B]
        pr = align([{key(r, n) for r, n in a} for a in A], [{key(r, n) for r, n in b} for b in Bs])
        # Tonhöhen der Pins (Ton aus Knopf + Balgrichtung der Pins), gelesen: Knopf + gelesene Balgrichtung
        pinp = {}
        for t, pi, r, num, dr in PINS[k]: pinp.setdefault(round(t, 3), set()).add(pi)
        PT = [pinp[t] for t in sorted(pinp)]
        for i, j in pr:
            a, b = A[i], Bs[j]; T['paar'] += 1; T['gleich'] += a == b
            T['ton_soll'] += len(PT[j])
            for r, n in b:
                T[f's{r}'] += 1; T[f'g{r}'] += (r, n) in a
                if (r, n) not in a: miss.append((k, res[i]['page'], res[i]['st'], round(res[i]['x']), r, n, sorted(a)))
            for r, n in a - b: T[f'f{r}'] += 1
        T['pins'] += len(B)
        # Tonhöhen-Vergleich (unabhängig von der Schreibweise): gleiche Töne im gepaarten Akkord
        for i, j in pr:
            pa = {p for p in (pitch(r, n, res[i]['dir'], D450) for _, _, (r, n) in res[i]['heads']) if p is not None}
            pb = PT[j]; T['ton_ok'] += len(pa & pb); T['ton_falsch'] += len(pa - pb); T['akk_ton_gleich'] += pa == pb
            bd = collections.Counter(d_ for _, _, d_ in B[j]).most_common(1)[0][0]; ok_d = res[i]['dir'] == bd; T['balg_ok'] += ok_d
            if not ok_d: T['ton_falsch_balg'] += len(pa - pb)
            if detail and pa != pb: wrong.append((k, res[i]['page'], res[i]['st'], round(res[i]['x']), res[i]['dir'], bd, sorted(pa), sorted(pb), sorted(A[i]), sorted(Bs[j])))
    print('Akkorde gepaart', T['paar'], 'von', T['pins'], '| Köpfe ganz gleich', T['gleich'], f"{T['gleich']/max(1,T['paar']):.3f}",
          '| Töne gleich (Akkord)', T['akk_ton_gleich'], f"{T['akk_ton_gleich']/max(1,T['paar']):.3f}",
          '| Töne richtig', T['ton_ok'], 'von', T['ton_soll'], f"{T['ton_ok']/max(1,T['ton_soll']):.3f}", 'falsch', T['ton_falsch'], 'davon wegen Balgrichtung', T['ton_falsch_balg'], '| Balg', f"{T['balg_ok']/max(1,T['paar']):.3f}")
    for r in (1, 2, 3, 4): print(f'  Reihe {r}: soll {T[f"s{r}"]}, gefunden {T[f"g{r}"]} ({T[f"g{r}"]/max(1,T[f"s{r}"]):.3f}), überzählig {T[f"f{r}"]}')
    if detail:
        json.dump(miss, open('/tmp/philipp/miss.json', 'w')); json.dump(wrong, open('/tmp/philipp/wrong.json', 'w'))
    return T
if __name__ == '__main__':
    run(sys.argv[2].split(','), len(sys.argv) > 3)
