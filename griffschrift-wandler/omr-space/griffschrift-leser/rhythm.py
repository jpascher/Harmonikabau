# Notenwerte: senkrechter Streifen um jeden Akkord (Kopf, Hals, Fähnchen/Balken, Punkt) -> Dauer in Vierteln
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2,sys,json,glob,random,numpy as np
sys.path.insert(0,HERE);import gsread as G
DUR=[0.25,1/3,0.5,0.75,1,1.5,2,3,4]
def dclass(q):
    return int(np.argmin([abs(np.log(q/d)) for d in DUR]))
SW,SH=24,144
def strip(nol_lines,cx,mid,d,fill):
    x0=int(cx-1.6*d);x1=int(cx+2.4*d);y0=int(mid-9*d);y1=int(mid+9*d)
    pad=cv2.copyMakeBorder(nol_lines,int(10*d),int(10*d),int(4*d),int(4*d),cv2.BORDER_CONSTANT,value=0)
    c=pad[y0+int(10*d):y1+int(10*d),x0+int(4*d):x1+int(4*d)]
    return np.concatenate([cv2.resize(c,(SW,SH),interpolation=cv2.INTER_AREA).ravel(),[fill*255]]).astype(np.uint8)
def page_views(g):
    _,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV);st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]);V={}
    def view(si):
        if si not in V: ys=st[si];d=float(np.mean(np.diff(ys)));V[si]=(G.remove_lines(bw,ys,d),d,ys[2])
        return V[si]
    return st,view

import os as _os
BEAMCONN=_os.environ.get('GS_BEAMCONN','1')=='1'   # nur am Hals hängende Balken/Fähnchen zählen (Fingersätze und Akzente über dem Hals, etwa bei Michlbauer, sonst als Fähnchen gelesen)
def beam_feats(nl,cx,ys_heads,d):
    """Hals suchen (längste senkrechte Linie nahe den Köpfen), an seinem Ende Balken/Fähnchen links und rechts zählen, Punkt rechts."""
    H,W=nl.shape;top=int(min(ys_heads));bot=int(max(ys_heads));best=(0,0,0,0)
    for x in range(int(cx-0.9*d),int(cx+0.9*d)+1):
        if not 0<=x<W: continue
        col=nl[:,x]>0
        for direction in (-1,1):
            y=(top if direction<0 else bot);start=y
            # vom Kopf weg laufen, solange Tinte
            gap=0;last=y
            while 0<=y+direction<H and abs(y-start)<8*d:
                y+=direction
                if col[y]: gap=0;last=y
                else:
                    gap+=1
                    if gap>max(2,0.3*d): break
            y=last;L=abs(y-start)   # Lücken bis 0,3 Linienabstand überspringen (dort lag eine Notenlinie)
            if L>best[0]: best=(L,x,y,direction)
    L,sx,ey,dr=best
    if L<1.8*d: return [L/d,0,0,0,0,0]
    out=[L/d,dr]
    for side in (-1,1):
        x=int(sx+side*0.45*d)   # Balken nur am Halsende zählen, nicht bis zum Kopf
        if not 0<=x<W: out.append(0);continue
        y0,y1=sorted((ey,int(ey-dr*max(0.6*d,min(1.7*d,L-1.3*d)))));seg=nl[max(0,y0-2):min(H,y1+2),x]>0
        if BEAMCONN and seg.any():
            # nur Striche zählen, die am Hals hängen (Balken, Fähnchen); Fingersätze und Akzente über dem Hals stehen frei
            ya=max(0,y0-2);xa,xb=sorted((sx,x));band=nl[ya:ya+len(seg),xa:xb+1]>0
            conn=band.all(1) if band.shape[1] else seg
            conn=conn|np.r_[False,conn[:-1]]|np.r_[conn[1:],False]
            seg=seg&conn
        runs=int(np.sum(np.diff(seg.astype(np.int8))==1)+(1 if len(seg) and seg[0] else 0));out.append(runs)
    # Punkt rechts vom Kopf (kleiner Fleck), auf Kopfhöhe oder knapp darüber
    hy=int(np.mean(ys_heads));win=nl[max(0,hy-int(.6*d)):hy+int(.3*d)+1,int(cx+0.9*d):int(cx+2*d)]
    out.append(float((win>0).mean()) if win.size else 0)
    out.append(len(ys_heads))
    return out

def dot_feat(nl,cx,ys_heads,d):
    """Verlängerungspunkt: kleiner runder Fleck rechts neben einem Kopf (bis 2,3 Linienabstände), auf Kopfhöhe oder knapp darüber."""
    H,W=nl.shape;x0=int(cx+0.6*d);x1=min(W,int(cx+2.3*d));best=0.0
    for hy in ys_heads:
        y0=max(0,int(hy-0.9*d));y1=min(H,int(hy+0.5*d))
        if x1<=x0 or y1<=y0: continue
        win=(nl[y0:y1,x0:x1]>0).astype(np.uint8);n,lab,s,_=cv2.connectedComponentsWithStats(win)
        for i in range(1,n):
            x,y,w,h,a=s[i]
            if x==0: continue   # hängt am Kopf
            if w<=0.75*d and h<=0.75*d and 0.08*d*d<=a<=0.45*d*d and a/(w*h)>0.5 and y>0 and y+h<win.shape[0]: best=max(best,1.0)
    return best
