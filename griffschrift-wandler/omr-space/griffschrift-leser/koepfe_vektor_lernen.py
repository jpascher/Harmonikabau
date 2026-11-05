# Aufruf: python3 koepfe_vektor_lernen.py daten.npz ../gsleser/heads3.pkl.gz
# Drittes Kopfmodell (heads3.pkl.gz) aus den Vektor-Lerndaten; hohle Kreuze wie beim Hauptmodell zu vollen Kreuzen zusammengelegt
import sys, numpy as np, pickle, gzip
from sklearn.ensemble import HistGradientBoostingClassifier as H
z = np.load(sys.argv[1]); X, Y = z['X'].astype(np.float32), z['Y'].copy(); Y[Y == 4] = 3
c = H(max_iter=300, learning_rate=0.1, l2_regularization=1.0, random_state=0).fit(X, Y)
pickle.dump(c, gzip.open(sys.argv[2], 'wb')); print('ok', X.shape, np.bincount(Y))
