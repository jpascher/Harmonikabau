# Endgültiges Anlernen aus den Zwischenspeichern: synthetische Ausschnitte (train.py) + echte Fuchs-Seiten (realdata.py)
import numpy as np,pickle,sys
from sklearn.ensemble import HistGradientBoostingClassifier as H
S=np.load('/tmp/claude-0/gsr/synth_cache.npz');R=np.load('/tmp/claude-0/gsr/real_train.npz')
rng=np.random.default_rng(0);RY=R['Y'];keep=np.r_[np.where(RY!=0)[0],rng.choice(np.where(RY==0)[0],min((RY==0).sum(),35000),replace=False)]
SY=S['Y'];sk=np.r_[np.where(SY!=0)[0],rng.choice(np.where(SY==0)[0],int((SY==0).sum()*0.8),replace=False)]   # Speicher
X=np.concatenate([S['X'][sk],R['X'][keep]]);Y=np.concatenate([SY[sk],RY[keep]]);Y[Y==4]=3
_,inv=np.unique(X,axis=0,return_inverse=True);inv=inv.ravel()
first={};bad=set()
for k,y in zip(inv,Y):
    if first.setdefault(k,y)!=y: bad.add(k)
ok=np.array([k not in bad for k in inv]);X,Y=X[ok],Y[ok];print(X.shape,np.bincount(Y),'Widersprüche',len(bad),flush=True)
c=H(max_iter=300,learning_rate=0.1,early_stopping=True,n_iter_no_change=20,random_state=0,l2_regularization=1.0,max_leaf_nodes=31,min_samples_leaf=40).fit(X,Y)
print(c.n_iter_,c.validation_score_[-1]);pickle.dump(c,open('/tmp/claude-0/gsr/heads.pkl','wb'));print('ok')
