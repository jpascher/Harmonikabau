# Fuchs-Bilder mit bekannten Griffen (FUCHS_PINS) abgleichen: erkannte Akkorde in Lesereihenfolge <-> Pin-Akkorde (DP)
import cv2,sys,json,glob,numpy as np
from PIL import Image
sys.path.insert(0,'/tmp/claude-0/gsr');import detect as D
OFF={1:-1,2:0,3:1,4:2}
def folder(k):
    for c in (k+'gs',k+'_gs'):
        fs=glob.glob(f'/tmp/claude-0/fuchs/fuchs/{c}/GS*')
        if fs: return fs[0]
def upscale(src,dst):
    g=np.array(Image.open(src).convert('L'));g=cv2.resize(g,None,fx=3,fy=3,interpolation=cv2.INTER_CUBIC);cv2.imwrite(dst,g);return dst
def chords_det(det):
    out=[]
    for st in sorted(set(h['st'] for h in det)):
        hs=sorted([h for h in det if h['st']==st],key=lambda h:h['x']);d=hs[0]['d']
        for h in hs:
            if out and out[-1][0]['st']==st and h['x']-np.median([q['x'] for q in out[-1]])<0.8*d: out[-1].append(h)
            else: out.append([h])
    return out
def chords_pin(n):
    m={}
    for t,pi,r,num,dr in n: m.setdefault(t,set()).add((2*(num-6)+OFF[r],'o' if r<3 else 'x'))
    return [m[t] for t in sorted(m)]
def align(A,B):
    # A: erkannte Akkorde (Mengen), B: Pin-Akkorde; Bewertung: Schnitt - halbe Differenz
    n,m=len(A),len(B);S=np.zeros((n+1,m+1));S[:,0]=-np.arange(n+1)*.5;S[0,:]=-np.arange(m+1)*.5;bt={}
    def sc(a,b): return len(a&b)*2-len(a^b)*.5 if a&b else -3
    for i in range(1,n+1):
        for j in range(1,m+1):
            c=[(S[i-1,j-1]+sc(A[i-1],B[j-1]),'m'),(S[i-1,j]-.5,'a'),(S[i,j-1]-.5,'b')];S[i,j],bt[i,j]=max(c)
    i,j,pairs=n,m,[]
    while i>0 and j>0:
        k=bt[i,j]
        if k=='m': pairs.append((i-1,j-1));i-=1;j-=1
        elif k=='a': i-=1
        else: j-=1
    return pairs[::-1]
if __name__=='__main__':
    P=json.load(open('/tmp/claude-0/gsr/pins.json'));res={}
    for k,v in P.items():
        img=upscale(folder(k),f'/tmp/claude-0/fuchs/{k}_u.png');det=D.detect(img,1)
        C=chords_det(det);A=[{(h['p'],h['t']) for h in c} for c in C];B=chords_pin(v['n'])
        pr=align(A,B);ex=sum(A[i]==B[j] for i,j in pr);good=sum(len(A[i]&B[j])>=max(1,len(B[j])/2) for i,j in pr)
        print(k,'Akkorde Bild',len(A),'Pins',len(B),'gepaart',len(pr),'ganz gleich',ex,'brauchbar',good)
        res[k]=dict(img=img,det=det,pairs=pr,B=[sorted(b) for b in B])
    json.dump(res,open('/tmp/claude-0/gsr/align.json','w'))
