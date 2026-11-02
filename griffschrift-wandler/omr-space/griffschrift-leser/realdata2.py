# Übungsbeispiele aus allen Fuchs-Seiten, geprüft über die MIDI-Datei (Knopf + Balgrichtung -> Ton muss zur MIDI passen)
import cv2,sys,json,random,numpy as np
sys.path.insert(0,'/tmp/claude-0/gsr');import gsread as G;from train import patch,CL;from assemble import button,pitch
random.seed(3)
def build(src,keys,out):
    R=json.load(open(src));X=[];Y=[];n=dict(pos=0,neg_extra=0,fix=0,neg=0)
    for k in keys:
        r=R[k];sh=r['sh'];g=cv2.imread(r['img'],0);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV)
        st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]);V={}
        def views(si):
            if si not in V: ys=st[si];d=float(np.mean(np.diff(ys)));V[si]=(G.stems_off(G.remove_lines(bw,ys,d),d),d,ys[2])
            return V[si]
        def add(si,x,y,c,jit=True):
            nol,d,_=views(si)
            for dx,dy in (((-2,0),(0,0),(2,1),(0,-1)) if jit else ((0,0),)):
                X.append(np.concatenate([patch(bw,x+dx,y+dy,d),patch(nol,x+dx,y+dy,d)]));Y.append(c)
        sure=[]
        for i,j in r['pairs']:
            c=r['res'][i];b=set(r['B'][j]);hs=c['heads'];si=c['st'];nol,d,mid=views(si)
            P=[(h['pitch']+sh) if h['pitch'] else None for h in hs];a={p for p in P if p}
            if a==b and len(set(P))==len(P):
                for h in hs: add(si,h['x'],h['y'],CL.index(h['t']+('_f' if h['fill'] else '_h')));n['pos']+=1;sure.append((si,h['x'],h['y']))
            elif len(hs)>=2 and any(({p for q,p in enumerate(P) if q!=e and p}==b) for e in range(len(hs))):
                e=next(e for e in range(len(hs)) if {p for q,p in enumerate(P) if q!=e and p}==b)
                for q,h in enumerate(hs):
                    if q==e: add(si,h['x'],h['y'],0,False);n['neg_extra']+=1
                    else: add(si,h['x'],h['y'],CL.index(h['t']+('_f' if h['fill'] else '_h')));n['pos']+=1;sure.append((si,h['x'],h['y']))
            elif len(hs)==len(b) and len(a&b)==len(b)-1:
                # ein Kopf falsch: andere Lage oder Form suchen, die genau den fehlenden Ton gibt
                miss=next(iter(b-a));bad=[q for q,p in enumerate(P) if p not in b]
                if len(bad)!=1: continue
                h=hs[bad[0]];cand=[]
                for p in range(h['p']-3,h['p']+4):
                    for t in 'ox':
                        rr,nn=button(p,t);q=pitch(rr,nn,c['dir'])
                        if q and q+sh==miss: cand.append((abs(p-h['p'])+(t!=h['t']),p,t))
                if not cand: continue
                _,p,t=min(cand);y=mid+0.5-p*d/2
                if abs(p-h['p'])>1: continue
                add(si,h['x'],y,CL.index(t+('_f' if h['fill'] else '_h')));n['fix']+=1
                if p!=h['p']: add(si,h['x'],h['y'],0,False);n['neg']+=1
                for q,hh in enumerate(hs):
                    if q!=bad[0]: add(si,hh['x'],hh['y'],CL.index(hh['t']+('_f' if hh['fill'] else '_h')));n['pos']+=1
        allh=[(c['st'],h['x'],h['y']) for c in r['res'] for h in c['heads']]
        for si,x,y in sure:
            nol,d,mid=views(si)
            for ax,ay in ((0,d/2),(0,-d/2),(d*.9,0),(-d*.9,0),(0,d),(0,-d)):
                if any(s==si and abs(x+ax-a)<d*.5 and abs(y+ay-b)<d*.3 for s,a,b in allh): continue
                add(si,x+ax,y+ay,0,False);n['neg']+=1
        for si in range(len(st)):
            nol,d,mid=views(si);W=bw.shape[1]
            ink=cv2.boxFilter((nol>0).astype(np.float32),-1,(max(3,int(d*.9)),max(3,int(d*.6))));cand=[]
            for p in range(-13,13):
                y=int(round(mid+0.5-p*d/2))
                if not 0<=y<bw.shape[0]: continue
                for x in range(int(d*2),W-int(d),4):
                    if ink[y,x]>=0.12 and not any(s==si and abs(x-a)<d*1.1 and abs(y-b)<d*.8 for s,a,b in allh): cand.append((x,y))
            for x,y in random.sample(cand,min(len(cand),400)): add(si,x,y,0,False);n['neg']+=1
    X=np.array(X,np.uint8);Y=np.array(Y);Y[Y==4]=3
    _,inv=np.unique(X,axis=0,return_inverse=True);inv=inv.ravel();lab={}
    for k2,y in zip(inv,Y): lab.setdefault(k2,set()).add(int(y))
    ok=np.array([len(lab[k2])==1 for k2 in inv]);X,Y=X[ok],Y[ok]
    np.savez_compressed(out,X=X,Y=Y);print(out,n,np.bincount(Y))
if __name__=='__main__':
    R=json.load(open(sys.argv[1]));ks=sorted(R);TEST=ks[::5]
    json.dump(TEST,open('/tmp/claude-0/gsr/midi_test.json','w'));print('Prüfstücke',TEST)
    build(sys.argv[1],[k for k in ks if k not in TEST],'/tmp/claude-0/gsr/real_train.npz')
