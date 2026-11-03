# Fuchs-Seiten lesen: Köpfe -> Knöpfe, Druckbalken unter den Bassbuchstaben -> Balgrichtung, Töne; Abgleich mit der MIDI-Datei
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2,sys,json,glob,numpy as np,mido
sys.path.insert(0,HERE);import detect as D, gsread as G
from align import chords_det,upscale
from assemble import button,pitch
def druck_bars(g):
    _,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV);st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]);out={}
    for si,ys in enumerate(st):
        d=float(np.mean(np.diff(ys)));y0=int(ys[4]+1.5*d);y1=int(min(ys[4]+14*d,(st[si+1][0]-2.5*d) if si+1<len(st) else g.shape[0]))
        if y1<=y0+2: out[si]=[];continue   # nächstes System zu nah
        band=bw[y0:y1];h=cv2.morphologyEx(band,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(int(d*3),1)))
        h=cv2.morphologyEx(h,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(1,max(2,int(d*.25)))))
        n,lab,s,_=cv2.connectedComponentsWithStats(h)
        L=[(s[i,0],s[i,0]+s[i,2],y0+s[i,1]) for i in range(1,n) if s[i,2]>2*d]
        # Nur Druckstriche zählen: deutlich unter dem System und direkt unter Bassbuchstaben (nicht Triolen- oder
        # Wiederholungsklammern, Bindebögen); Buchstaben = Teile in Schriftgröße knapp über dem Strich in seinem Bereich
        if L and _os.environ.get('GS_DRUCK','1')=='1':
            nb,_,sb,_=cv2.connectedComponentsWithStats(band)
            let=[(sb[i,0],sb[i,0]+sb[i,2],y0+sb[i,1]+sb[i,3]) for i in range(1,nb) if 0.6*d<=sb[i,3]<=3.2*d and 0.2*d<=sb[i,2]<=3.0*d]
            def ok(l):
                a,b,y=l
                if y<ys[4]+2.5*d: return False
                return any(x1>a-0.5*d and x0<b+0.5*d and y-1.6*d<=yb<=y+0.2*d for x0,x1,yb in let)
            L=[l for l in L if ok(l)]
        out[si]=L
    return st,out
def read_page(img):
    g=cv2.imread(img,0);st,bars=druck_bars(g);det=D.detect(img,1);C=chords_det(det);res=[]
    for c in C:
        si=c[0]['st'];x=float(np.median([h['x'] for h in c]));dr='D' if any(a-3<=x<=b+3 for a,b,_ in bars.get(si,[])) else 'Z'
        ks=[button(h['p'],h['t']) for h in c];ps=[pitch(r,n,dr) for r,n in ks]
        res.append(dict(st=si,x=x,dir=dr,keys=ks,pitches=sorted(p for p in ps if p),heads=[dict(x=h['x'],y=h['y'],p=h['p'],t=h['t'],fill=bool(h['fill']),pitch=q) for h,q in zip(c,ps)]))
    return res,bars
def midi_chords(path):
    m=mido.MidiFile(path);best=None
    for t in m.tracks:
        ns=[];tt=0
        for x in t:
            tt+=x.time
            if x.type=='note_on' and x.velocity>0: ns.append((tt,x.note))
        if ns and (best is None or np.mean([n for _,n in ns])>np.mean([n for _,n in best])): best=ns
    q=m.ticks_per_beat/8;g={}
    for tt,n in best: g.setdefault(round(tt/q),set()).add(n)
    return [g[k] for k in sorted(g)]
def dp(A,B):
    n,m=len(A),len(B);S=np.zeros((n+1,m+1));S[:,0]=-np.arange(n+1)*.3;S[0,:]=-np.arange(m+1)*.3;bt=np.zeros((n+1,m+1),np.int8)
    for i in range(1,n+1):
        for j in range(1,m+1):
            a,b=A[i-1],B[j-1];s=(2 if a==b else len(a&b)*.7-len(a^b)*.3) if a&b else -1.5
            c=(S[i-1,j-1]+s,S[i-1,j]-.3,S[i,j-1]-.3);k=int(np.argmax(c));S[i,j]=c[k];bt[i,j]=k
    i,j,pr=n,m,[]
    while i>0 and j>0:
        k=bt[i,j]
        if k==0: pr.append((i-1,j-1));i-=1;j-=1
        elif k==1: i-=1
        else: j-=1
    return S[n,m],pr[::-1]
def check(folder,quiet=False):
    gs=glob.glob(folder+'/GS*');ms=glob.glob(folder+'/*.mid')
    if not gs or not ms: return None
    k=folder.rstrip('/').split('/')[-1];img=upscale(gs[0],f'/tmp/claude-0/fuchs/{k}_u.png')
    res,_=read_page(img);B=midi_chords(ms[0]);best=None
    for sh in range(-14,15):
        A=[{p+sh for p in r['pitches']} for r in res];sc,pr=dp(A,B)
        if best is None or sc>best[0]: best=(sc,sh,pr,A)
    sc,sh,pr,A=best;ex=sum(A[i]==B[j] for i,j in pr)
    return dict(k=k,img=img,sh=sh,n=len(A),midi=len(B),paired=len(pr),exact=ex,pairs=pr,res=res,B=[sorted(b) for b in B])
if __name__=='__main__':
    out={};T=[0,0,0]
    for f in sorted(glob.glob('/tmp/claude-0/fuchs/fuchs/*/')):
        if 'iliana' in f: continue   # 4-Reiher mit Mollbässen, eigener Fall
        if not glob.glob(f+'GS*'): continue
        r=check(f)
        if not r: continue
        print(r['k'],'Akkorde',r['n'],'MIDI',r['midi'],'gepaart',r['paired'],'gleich',r['exact'],'Versatz',r['sh'],flush=True)
        T[0]+=r['n'];T[1]+=r['paired'];T[2]+=r['exact'];out[r['k']]=r
        json.dump(out,open((sys.argv[1] if len(sys.argv)>1 else HERE+'/midicheck.json')+'.part','w'),default=lambda o:o.item() if hasattr(o,'item') else str(o))
    print('gesamt Akkorde',T[0],'gepaart',T[1],'gleich',T[2],round(T[2]/max(1,T[1]),3))
    json.dump(out,open(sys.argv[1] if len(sys.argv)>1 else HERE+'/midicheck.json','w'),default=lambda o:o.item() if hasattr(o,'item') else str(o))
