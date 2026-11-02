import numpy as np,pickle,collections
from sklearn.ensemble import HistGradientBoostingClassifier as H
Z=np.load('/tmp/claude-0/gsr/rhythm.npz');X,Y,K=Z['X'],Z['Y'],Z['K']
tr=np.isin(K,['S','R'])
w=np.where(K[tr]=='R',3.0,1.0)
c=H(max_iter=300,learning_rate=0.1,l2_regularization=1.0,min_samples_leaf=20,early_stopping=True,n_iter_no_change=20,random_state=0).fit(X[tr],Y[tr],sample_weight=w)
pickle.dump(c,open('/tmp/claude-0/gsr/rhythm.pkl','wb'))
DUR=['1/16','1/8T','1/8','1/8.','1/4','1/4.','1/2','1/2.','1/1']
for k in ['T','RT']:
    m=K==k;p=c.predict(X[m]);print(k,'Genauigkeit',round(float((p==Y[m]).mean()),3),'n',int(m.sum()))
    E=collections.Counter((DUR[a],DUR[b]) for a,b in zip(Y[m],p) if a!=b);print('  Fehler (soll, gelesen):',E.most_common(8))
