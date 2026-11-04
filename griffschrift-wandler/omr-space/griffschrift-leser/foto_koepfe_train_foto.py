# Kopfmodell für Fotos: von Hand bestätigte Fotoseiten (stark gewichtet) + verschlechterte Wandler-Bilder
import numpy as np, pickle, gzip, sys
from sklearn.ensemble import HistGradientBoostingClassifier as H
keys = sys.argv[1].split(','); out = sys.argv[2]
Xs = []; Ys = []; Ws = []
for k in keys:
    z = np.load(f'/home/claude/lab/{k}_ds.npz'); Xs.append(z['X']); Ys.append(z['Y']); Ws.append(np.full(len(z['Y']), 6.0))
Z = np.load('/home/claude/hol/hol.npz'); X, Y = Z['X'], Z['Y']; rng = np.random.default_rng(1)
sel = np.concatenate([rng.choice(np.where(Y == k)[0], n, replace=False) for k, n in ((0, 20000), (1, 10000), (2, 6000), (3, 6000))])
Xs.append(X[sel]); Ys.append(Y[sel]); Ws.append(np.ones(len(sel))); del X, Z
X = np.concatenate(Xs).astype(np.float32); Y = np.concatenate(Ys); W = np.concatenate(Ws)
c = H(max_iter=300, learning_rate=0.1, l2_regularization=1.0, random_state=0).fit(X, Y, sample_weight=W)
pickle.dump(c, gzip.open(out, 'wb')); print('ok', X.shape, np.bincount(Y))
