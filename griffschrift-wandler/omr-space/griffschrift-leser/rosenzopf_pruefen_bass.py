# Gelesene Bassbuchstaben (ganze Kette wie im Space) gegen die abgelesenen, je System per Abgleich
import sys, os, numpy as np, cv2
sys.path.insert(0, '/home/claude/gsw/m'); sys.path.insert(0, '/home/claude/rosen')
import gs2xml as X, gsread as G, bassread as BR
from truth import T
if len(sys.argv) > 1: BR._M = __import__('pickle').load(__import__('gzip').open(sys.argv[1], 'rb'))
keys = sys.argv[2].split(',') if len(sys.argv) > 2 else list(T)
def al(A, B):
    n, m = len(A), len(B); D = np.zeros((n + 1, m + 1)); D[:, 0] = np.arange(n + 1); D[0, :] = np.arange(m + 1)
    for i in range(1, n + 1):
        for j in range(1, m + 1): D[i, j] = min(D[i - 1, j - 1] + (A[i - 1] != B[j - 1]), D[i - 1, j] + 1, D[i, j - 1] + 1)
    i, j = n, m; ok = 0
    while i > 0 and j > 0:
        if D[i, j] == D[i - 1, j - 1] + (A[i - 1] != B[j - 1]): ok += A[i - 1] == B[j - 1]; i -= 1; j -= 1
        elif D[i, j] == D[i - 1, j] + 1: i -= 1
        else: j -= 1
    return ok
S = np.zeros(3)
for k in keys:
    g = X.load_gray(f'/home/claude/rosen/x/20261103_{k}.jpg'); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV); st = G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1])
    for si, tr in enumerate(T[k]):
        A = [c for _, c in BR.read_letters(bw, st, si, (st[si][4] - st[si][0]) / 4)]; ok = al(A, tr); S += [len(tr), len(A), ok]
        print(k, si, 'soll', len(tr), 'gelesen', len(A), 'richtig', ok, ' '.join(A))
print('GESAMT soll %d gelesen %d richtig %d (%.3f)' % (S[0], S[1], S[2], S[2] / S[0]))
