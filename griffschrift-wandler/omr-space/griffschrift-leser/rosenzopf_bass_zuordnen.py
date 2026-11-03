# Flecken der Buchstabenzeile (samt übereinander gesetzten) mit den von Hand abgelesenen Buchstaben paaren -> Lerndaten
import sys, numpy as np, cv2, collections
sys.path.insert(0, '/home/claude/gsw/m'); sys.path.insert(0, '/home/claude/rosen')
import gs2xml as X, gsread as G, bassread as BR
from truth import T
XS = []; YS = []; KS = []; M = BR.model()
for key, rows in T.items():
    g = X.load_gray(f'/home/claude/rosen/x/20261103_{key}.jpg'); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    st = G.staves(cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)[1]); assert len(st) == len(rows), (len(st), len(rows))
    for si, tr in enumerate(rows):
        d = (st[si][4] - st[si][0]) / 4; C = BR.candidates(bw, st, si, d)
        P = M.predict_proba(np.array([f for _, f in C])); cl = M.classes_[P.argmax(1)]
        L = [(b, f, c) for (b, f), c in zip(C, cl) if c not in ('none', 'mark')]
        yc = np.median([b[1] + b[3] / 2 for b, _, _ in L]) if L else 0
        cand = [(b, f, c) for (b, f), c in zip(C, cl) if abs(b[1] + b[3] / 2 - yc) <= 2.4 * d and 0.5 * d <= b[3] <= 2.4 * d and b[2] >= 0.3 * d and b[3] <= 3.5 * b[2]]
        cand.sort(key=lambda t: (round((t[0][0] + t[0][2] / 2) / (0.8 * d)), -t[0][1]))
        A = [str(c) for _, _, c in cand]; B = []
        for t in tr: B += t.split('+') if '+' in t else [t]
        n, m = len(A), len(B); D = np.zeros((n + 1, m + 1)); D[:, 0] = np.arange(n + 1) * 0.6; D[0, :] = np.arange(m + 1) * 1.0; bt = np.zeros((n + 1, m + 1), np.int8)
        for i in range(1, n + 1):
            for j in range(1, m + 1):
                sub = 0 if A[i - 1] == B[j - 1] else (0.3 if A[i - 1].lower() == B[j - 1].lower() else (0.8 if A[i - 1] == 'none' else 0.9))
                o = [D[i - 1, j - 1] + sub, D[i - 1, j] + 0.6, D[i, j - 1] + 1.0]; k = int(np.argmin(o)); D[i, j] = o[k]; bt[i, j] = k
        i, j = n, m; lab = ['none'] * n; miss = 0
        while i > 0 and j > 0:
            k = bt[i, j]
            if k == 0: lab[i - 1] = B[j - 1]; i -= 1; j -= 1
            elif k == 1: i -= 1
            else: j -= 1; miss += 1
        miss += j
        ok = sum(1 for a, l in zip(A, lab) if a == l and l != 'none'); nl = sum(1 for l in lab if l != 'none')
        print(key, si, 'Soll', m, 'gepaart', nl, 'davon schon richtig', ok, 'fehlend', miss)
        for (b, f, c), l in zip(cand, lab): XS.append(f); YS.append(l); KS.append(f'R:{key}:{si}')
np.savez_compressed('/home/claude/rosen/rosen.npz', X=np.array(XS), Y=np.array(YS), K=np.array(KS))
print(collections.Counter(YS))
