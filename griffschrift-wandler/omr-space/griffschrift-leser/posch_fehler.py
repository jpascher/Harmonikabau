# Wo weichen die Töne bei „Franz Posch“ ab? (nur Messung) Für jedes zugeordnete Akkordpaar Griffschrift-Leser / homr-Normalnoten:
# gleich | Kopf fehlt | Kopf zu viel | Halbton daneben (Vorzeichen/Tonart?) | Oktave | sonst; dazu je Reihe und Kopfform die falschen
# Köpfe des Lesers und je Tonklasse die Halbton-Abweichungen (gleiche Tonklasse immer wieder = Vorzeichen der Tonart).
# Braucht die homr-Lesungen aus posch_vergleich.py (/content/posch/*.musicxml).
# Aufruf: python3 posch_fehler.py ../gsleser GRIFFSCHRIFT.pdf NORMAL.pdf
import sys, collections
import posch_vergleich as P
import gs2xml as X
from align import align
P.PDF_G, P.PDF_N = sys.argv[2], sys.argv[3]
X.GERUEST = True
NAMEN = 'C Cis D Dis E F Fis G Gis A B H'.split()
K = collections.Counter(); R = collections.Counter(); HT = collections.Counter(); TA = collections.Counter(); BSP = []
for name, gs_, ns in P.STUECKE:
    N = [x for i in ns for x in P.homr_noten(P.seite(P.PDF_N, i, 'n', 300))]
    G = []; t = 0
    for p in [P.seite(P.PDF_G, i, 'g', 200) for i in gs_]:
        meas, L = X.read(p, -2)
        for m in meas:
            for e in m:
                if not e['rest']: G.append((t, e['heads']))
                t += e['u']
    A = [{h[2] for h in hs} for _, hs in G]; B = [s for _, s in N]
    pr = align(A, B); K['ohne_partner_normal'] += len(B) - len(pr); K['ohne_partner_leser'] += len(A) - len(pr)
    for i, j in pr:
        a, b = A[i], B[j]
        if a == b: K['gleich'] += 1; continue
        if a < b: K['kopf_fehlt'] += 1
        elif a > b: K['kopf_zuviel'] += 1
        else:
            ra, rb = sorted(a - b), sorted(b - a)
            if len(ra) == len(rb) and all(abs(x - y) == 1 for x, y in zip(ra, rb)):
                K['halbton'] += 1
                for x, y in zip(ra, rb): HT[f'{NAMEN[y % 12]} gelesen als {NAMEN[x % 12]}'] += 1
            elif len(ra) == len(rb) and all(abs(x - y) == 12 for x, y in zip(ra, rb)): K['oktave'] += 1
            elif len(ra) == len(rb) and all(abs(x - y) in (2, 3, 4) for x, y in zip(ra, rb)): K['terz_sekunde'] += 1
            else: K['sonst'] += 1
        for (r, n, p) in G[i][1]:
            R[(r, 'falsch' if p not in b else 'richtig')] += 1
            if p not in b: TA[(r, n)] += 1
        if len(BSP) < 25: BSP.append((name, sorted(a), sorted(b), [(r, n) for r, n, p in G[i][1]]))
print('Paare:', dict(K))
print('Köpfe je Reihe (richtig/falsch):', {r: (R[(r, 'richtig')], R[(r, 'falsch')]) for r in (1, 2, 3, 4)})
print('Halbton-Abweichungen (Normalnote gelesen als):', HT.most_common(12))
print('häufigste falsche Knöpfe (Reihe, Knopf):', TA.most_common(12))
for b in BSP: print(b)
