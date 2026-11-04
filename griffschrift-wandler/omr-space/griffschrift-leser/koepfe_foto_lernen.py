import numpy as np, pickle, gzip, sys, collections
from sklearn.ensemble import HistGradientBoostingClassifier as H
Z = np.load('/home/claude/hol/hol.npz'); X, Y, K = Z['X'], Z['Y'], Z['K']
rng = np.random.default_rng(0); te = (K % 5) == 0
def pick(m, n): i = np.where(m)[0]; return rng.choice(i, min(n, len(i)), replace=False)
tr = np.concatenate([pick(~te & (Y == 0), 70000), pick(~te & (Y == 1), 35000), pick(~te & (Y == 2), 15000), pick(~te & (Y == 3), 25000)])
ts = np.concatenate([pick(te & (Y == k), 8000) for k in range(4)])
Xtr, Ytr, Xts, Yts = X[tr].astype(np.float32), Y[tr], X[ts], Y[ts]; del X, Z; import gc; gc.collect()
w = np.where(Ytr == 2, 3.0, 1.0)
c = H(max_iter=300, learning_rate=0.1, l2_regularization=1.0, random_state=0, early_stopping=True).fit(Xtr, Ytr, sample_weight=w); print('fertig', flush=True)
pickle.dump(c, gzip.open('/home/claude/hol/heads2.pkl.gz', 'wb'))
p = c.predict(Xts); y = Yts
old = pickle.load(open('/home/claude/gsw/m/heads.pkl', 'rb')); po = old.predict(Xts); po[po == 4] = 3
for k, nm in enumerate(['none', 'o_f', 'o_h', 'x']):
    m = y == k; print(nm, int(m.sum()), 'neu', round((p[m] == k).mean(), 3), 'alt', round((po[m] == k).mean(), 3))
