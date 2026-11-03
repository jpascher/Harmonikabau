# Bassbuchstaben auf Fuchs-Seiten gegen FUCHS_BASS: Buchstabenfolge (Abgleich) und Zeitpunkte (bei bester Verschiebung)
import sys, re, json, glob, numpy as np
sys.path.insert(0, sys.argv[1] if len(sys.argv) > 1 else '/home/claude/gsw/m')
import gs2xml as X
s = open('/home/claude/zipm/griffschrift-wandler/index.html', encoding='utf-8').read()
FB = json.loads(re.search(r'window\.FUCHS_BASS=(\{.*?\});', s).group(1))
def align(A, B):
    n, m = len(A), len(B); D = np.zeros((n + 1, m + 1)); D[:, 0] = np.arange(n + 1); D[0, :] = np.arange(m + 1)
    for i in range(1, n + 1):
        for j in range(1, m + 1): D[i, j] = min(D[i - 1, j - 1] + (A[i - 1] != B[j - 1]), D[i - 1, j] + 1, D[i, j - 1] + 1)
    i, j = n, m; ok = 0; ci = 0; P_ = []
    while i > 0 and j > 0:
        if D[i, j] == D[i - 1, j - 1] + (A[i - 1] != B[j - 1]): ok += A[i - 1] == B[j - 1]; P_.append((i - 1, j - 1)); ci += A[i - 1].lower() == B[j - 1].lower(); i -= 1; j -= 1
        elif D[i, j] == D[i - 1, j] + 1: i -= 1
        else: j -= 1
    return ok, ci, P_
T = np.zeros(5)
for key, tru in FB.items():
    f = glob.glob(f'/home/claude/fuchs/zip/x/fuchs/{key}gs/GS*') + glob.glob(f'/home/claude/fuchs/zip/x/fuchs/{key}_gs/GS*')
    xml, P = X.convert([sorted(f)[0]], 0); b = P.get('b', [])
    A = [c for _, c in b]; B = [l for _, l in tru]; ok, ci, pr = align(A, B); pr = pr[::-1]
    dif = [b[i][0] - tru[j][0] for i, j in pr]; loc = 0
    for k in range(len(dif)):
        w = dif[max(0, k - 4):k + 5]
        loc += abs(dif[k] - np.median(w)) < 0.05
    T[4] += 0
    tt = np.array([t for t, _ in tru]); best = 0
    for sh in np.arange(-16, 16.01, 0.25):
        best = max(best, sum(1 for t, c in b if np.min(abs(tt - (t + sh))) < 0.05 and tru[int(np.argmin(abs(tt - (t + sh))))][1] == c))
    T += [len(B), len(A), ok, ci, loc]; best = loc
    print(f'{key}: soll {len(B)} gelesen {len(A)} gleich {ok} (ohne Groß/klein {ci}) Zeit örtlich stimmig {best}', flush=True)
print('GESAMT soll %d gelesen %d Buchstabe richtig %d (%.3f) ohne Groß/klein %d, Zeit und Buchstabe richtig %d (%.3f)' % (T[0], T[1], T[2], T[2] / T[0], T[3], T[4], T[4] / T[0]))
