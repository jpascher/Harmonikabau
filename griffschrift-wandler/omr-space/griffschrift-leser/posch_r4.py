# Colab: Posch – würden die Köpfe der 4. Reihe als 3. Reihe (eine Lage höher oder tiefer) zur Normalnotation passen?
import sys, collections
sys.argv = ['x', '../gsleser', '/content/drive/MyDrive/griffschrift_daten/Franz Posch Griffschrift.pdf', '/content/drive/MyDrive/griffschrift_daten/Franz Posch Normalnotation.pdf']
import posch_vergleich as P
import gs2xml as X
from align import align
from assemble import OFF, pitch
P.PDF_G, P.PDF_N = sys.argv[2], sys.argv[3]; X.GERUEST = True
C = collections.Counter()
for name, gs_, ns in P.STUECKE:
    N = [x for i in ns for x in P.homr_noten(P.seite(P.PDF_N, i, 'n', 300))]
    G = []
    for p in [P.seite(P.PDF_G, i, 'g', 200) for i in gs_]:
        meas, L = X.read(p, -2)
        for m in meas:
            for e in m:
                if not e['rest']: G.append(e)
    A = [{h[2] for h in e['heads']} for e in G]; B = [s for _, s in N]
    for i, j in align(A, B):
        e = G[i]
        for (r, n, pt) in e['heads']:
            if r != 4: continue
            C['r4'] += 1; C['r4_richtig'] += pt in B[j]
            p = (n - 6) * 2 + OFF[4]
            for d in (-1, 1):
                n3 = (p + d - OFF[3]) // 2 + 6; q = pitch(3, n3, e['dir'], None)
                if q is not None: q += 0
                C[f'r3_{d:+d}'] += (q is not None and (q - 2 * 0) in B[j])
print(dict(C))
