import numpy as np, pickle, collections, sys
from sklearn.ensemble import HistGradientBoostingClassifier as H
A = np.load('/home/claude/bass/ds.npz'); S = np.load('/home/claude/bass/synth.npz'); Fu = np.load('/home/claude/bass/fuchs.npz'); fm = Fu['Y'] != 'none'
X = np.concatenate([A['X'], S['X'], Fu['X'][fm]]); Y = np.concatenate([A['Y'], S['Y'], Fu['Y'][fm]]); K = np.concatenate([A['K'], S['K'], Fu['K'][fm]])
X = X.copy(); X[:, 402] = 0.0
w = np.where(np.char.startswith(K.astype(str), 'F:'), 2.0, 1.0)
if '--ohnefuchs' in sys.argv: keep=~np.char.startswith(K.astype(str), 'F:'); X,Y,K,w=X[keep],Y[keep],K[keep],w[keep]
HELD_P = {f'P:{p}' for p in (14, 15, 16, 17, 18, 19, 20, 25)}
wk = sorted({k for k in K if k.startswith('W:')}); HELD_W = set(wk[::4]) | {'F:sternpolka', 'F:bluemerl', 'F:warschaueraustragoess'}
te = np.isin(K, list(HELD_P | HELD_W)); final = 'final' in sys.argv
tr = np.ones(len(Y), bool) if final else ~te
c = H(max_iter=300, learning_rate=0.1, l2_regularization=1.0, random_state=0).fit(X[tr], Y[tr], sample_weight=w[tr])
if final: import gzip; pickle.dump(c, gzip.open('/home/claude/bass/bass.pkl.gz', 'wb')); print('gespeichert'); sys.exit()
for nm, m in [('Philipp', np.isin(K, list(HELD_P))), ('Wandler', np.isin(K, list(HELD_W)) & ~np.char.startswith(K.astype(str), 'F:')), ('Fuchs', np.char.startswith(K.astype(str), 'F:') & np.isin(K, list(HELD_W)))]:
    p = c.predict(X[m]); y = Y[m]; let = y != 'none'
    print(nm, 'n', int(m.sum()), 'alle', round((p == y).mean(), 3), 'Buchstaben richtig', round((p[let] == y[let]).mean(), 3), 'falsche Buchstaben aus none', int(((y == 'none') & (p != 'none')).sum()))
    print('  Fehler', collections.Counter((a, b) for a, b in zip(y, p) if a != b).most_common(10))
