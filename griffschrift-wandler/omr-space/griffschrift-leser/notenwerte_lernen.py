import numpy as np,pickle
from sklearn.ensemble import HistGradientBoostingClassifier as H
Zs=[np.load(f) for f in ['rh/P.npz','rh/F.npz']];X=np.concatenate([z['X'] for z in Zs]);Y=np.concatenate([z['Y'] for z in Zs]);K=np.concatenate([z['K'] for z in Zs])
tr=~np.isin(K,['P:eiswind','P:isabellile','P:zehnfuenfer','P:bobblpolka'])
c=H(max_iter=400,learning_rate=0.08,l2_regularization=1.0,min_samples_leaf=10,random_state=0).fit(X[tr],Y[tr]);pickle.dump(c,open('rh/rhythm2.pkl','wb'));print('fertig',int(tr.sum()))
