# Echte Fuchs-Seiten -> Trainingsausschnitte aus den mit FUCHS_PINS abgeglichenen Akkorden
import cv2,sys,json,random,numpy as np
sys.path.insert(0,'/tmp/claude-0/gsr');import gsread as G;from train import patch,CL;from align import chords_det
random.seed(2)
def build(keys,out):
    R=json.load(open('/tmp/claude-0/gsr/align.json'));X=[];Y=[];n=dict(pos=0,fix=0,neg=0)
    for k in keys:
        r=R[k];g=cv2.imread(r['img'],0);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV)
        st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]);C=chords_det(r['det']);B=[{tuple(x) for x in b} for b in r['B']]
        NOL={}
        def views(si):
            if si not in NOL:
                ys=st[si];d=float(np.mean(np.diff(ys)));NOL[si]=(G.stems_off(G.remove_lines(bw,ys,d),d),d,ys[2])
            return NOL[si]
        def add(si,x,y,c,jit=True):
            nol,d,_=views(si)
            for dx,dy in (((-2,0),(0,0),(2,1),(0,-1)) if jit else ((0,0),)):
                X.append(np.concatenate([patch(bw,x+dx,y+dy,d),patch(nol,x+dx,y+dy,d)]));Y.append(c)
        sure=[]
        for i,j in r['pairs']:
            c=C[i];A={(h['p'],h['t']) for h in c};miss=B[j]-A;extra=A-B[j]
            if not (A==B[j] or (len(miss)==1 and len(extra)==1)): continue
            for h in c:
                si=h['st'];nol,d,mid=views(si)
                if (h['p'],h['t']) in B[j]:
                    add(si,h['x'],h['y'],CL.index(h['t']+('_f' if h['fill'] else '_h')));n['pos']+=1;sure.append((si,h['x'],h['y']))
                else:
                    p,t=next(iter(miss));y=mid+0.5-p*d/2
                    add(si,h['x'],y,CL.index(t+('_f' if h['fill'] else '_h')));n['fix']+=1;sure.append((si,h['x'],y))
                    add(si,h['x'],h['y'],0,False);n['neg']+=1
        # harte Negative neben gesicherten Köpfen, zufällige Negative fern aller Erkennungen
        allh=[(h['st'],h['x'],h['y']) for h in r['det']]+sure
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
            for x,y in random.sample(cand,min(len(cand),700)): add(si,x,y,0,False);n['neg']+=1
    X=np.array(X,np.uint8);Y=np.array(Y);Y[Y==4]=3   # hohle Kreuze zu selten -> Kreuz
    # gleiche Ausschnitte mit verschiedenen Klassen (Widersprüche) entfernen
    _,inv=np.unique(X,axis=0,return_inverse=True);inv=inv.ravel();lab={}
    for k,y in zip(inv,Y): lab.setdefault(k,set()).add(int(y))
    ok=np.array([len(lab[k])==1 for k in inv]);X,Y=X[ok],Y[ok]
    np.savez_compressed(out,X=X,Y=Y);print(out,n,np.bincount(Y))
if __name__=='__main__':
    TEST=['sternpolka','warschaueraustragoess']
    keys=[k for k in json.load(open('/tmp/claude-0/gsr/align.json')) if k not in TEST]
    build(keys,'/tmp/claude-0/gsr/real_train.npz')
