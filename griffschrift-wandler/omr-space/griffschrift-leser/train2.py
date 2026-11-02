# Endgültiges Anlernen aus den Zwischenspeichern: synthetische Ausschnitte (train.py) + echte Fuchs-Seiten (realdata.py)
import numpy as np,pickle,sys
from sklearn.ensemble import HistGradientBoostingClassifier as H
S=np.load('/tmp/claude-0/gsr/synth_cache.npz');R=np.load('/tmp/claude-0/gsr/real_train.npz')
X=np.concatenate([S['X'],R['X'],R['X']]);Y=np.concatenate([S['Y'],R['Y'],R['Y']]);Y[Y==4]=3
_,inv=np.unique(X,axis=0,return_inverse=True);inv=inv.ravel()
first={};bad=set()
for k,y in zip(inv,Y):
    if first.setdefault(k,y)!=y: bad.add(k)
ok=np.array([k not in bad for k in inv]);X,Y=X[ok],Y[ok];print(X.shape,np.bincount(Y),'Widersprüche',len(bad),flush=True)
c=H(max_iter=300,learning_rate=0.1,early_stopping=True,n_iter_no_change=20,random_state=0,l2_regularization=1.0,max_leaf_nodes=31,min_samples_leaf=40).fit(X,Y)
print(c.n_iter_,c.validation_score_[-1]);pickle.dump(c,open('/tmp/claude-0/gsr/heads.pkl','wb'));print('ok')
