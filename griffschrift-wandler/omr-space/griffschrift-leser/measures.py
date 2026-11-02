# Taktstriche und Pausen auf Fuchs-Seiten
import cv2,sys,numpy as np
sys.path.insert(0,'/tmp/claude-0/gsr');import gsread as G
def barlines(bw,ys,d,heads_x):
    """Senkrechte Striche genau von der obersten bis zur untersten Linie, ohne Kopf daneben."""
    top,bot=int(ys[0]),int(ys[4]);W=bw.shape[1];cols=[]
    inside=(bw[top:bot+1]>0).mean(0)
    above=(bw[max(0,top-int(.8*d)):top-2]>0).mean(0);below=(bw[bot+3:bot+int(.8*d)]>0).mean(0)
    for x in range(W):
        if inside[x]>0.92 and above[x]<0.3 and below[x]<0.3: cols.append(x)
    out=[]
    for x in cols:
        if out and x-out[-1][-1]<=2: out[-1].append(x)
        else: out.append([x])
    res=[]
    for g in out:
        x=float(np.mean(g))
        if len(g)>0.6*d: continue           # zu breit: Balken o. ä.
        if any(abs(x-h)<0.9*d for h in heads_x): continue   # Hals eines Kopfes
        if res and x-res[-1]<0.9*d: res[-1]=(res[-1]+x)/2;continue   # Doppelstrich zusammenfassen
        res.append(x)
    return res

UNIT=[3,4,6,9,12,18,24,36,48]   # Notenwerte aus rhythm.DUR in Zwölftel-Vierteln
def rests(bw,nol,ys,d,heads,bars):
    """Pausen: Tintenflecken im System, die kein Kopf, Hals, Balken oder Taktstrich sind; nach Größe eingeteilt."""
    top,bot=int(ys[0]-d),int(ys[4]+d);band=nol[top:bot].copy()
    band=cv2.morphologyEx(band,cv2.MORPH_CLOSE,cv2.getStructuringElement(cv2.MORPH_RECT,(1,max(3,int(d*.4)))))   # Lücken der entfernten Linien schließen
    n,lab,s,_=cv2.connectedComponentsWithStats((band>0).astype(np.uint8));out=[];comps=[]
    for i in range(1,n):
        x,y,w,h,a=s[i];comps.append((x+w/2,top+y+h/2,w/d,h/d,a/(d*d)))
    for cx,cy,w,h,a in comps:
        if any(abs(cx-hx)<1.0*d and abs(cy-hy)<1.2*d for hx,hy in heads): continue
        if any(abs(cx-b)<0.6*d for b in bars) or cx<4*d: continue
        c=(cy-ys[2])/d;fillr=a/(w*h)
        if 2.5<=h<=3.5 and 0.8<=w<=1.4 and 0.38<=fillr<=0.65 and abs(c)<1: u=12            # Viertelpause
        elif 1.5<=h<=2.2 and 0.9<=w<=1.4 and 0.35<=fillr<=0.6 and abs(c)<0.8: u=6        # Achtelpause
        elif 0.5<=h<=0.95 and 1.05<=w<=1.7 and fillr>0.6 and -1.1<c<0.1: u=24 if c>-0.5 else 48   # Halbe/Ganze Pause
        else: continue
        # Punkt rechts daneben verlängert um die Hälfte
        if any(0.6*d<px-cx<2.2*d and abs(py-cy)<0.8*d and 0.1<pa<0.3 and pw<0.7 for px,py,pw,ph,pa in comps): u=u*3//2
        out.append((cx,u))
    return out
def fit_measure(P,rest_units,L,maxloss=6.0):
    """P: Wahrscheinlichkeiten je Akkord über UNIT; wählt die Werte mit größter Wahrscheinlichkeit, deren Summe (mit Pausen) L ergibt."""
    T=L-sum(rest_units)
    if T<=0 or not len(P): return None
    lp=np.log(np.clip(P,1e-6,1));best={0:(0.0,[])}
    for row in lp:
        nb={}
        for t,(sc,ch) in best.items():
            for k,u in enumerate(UNIT):
                if t+u<=T and (t+u not in nb or nb[t+u][0]<sc+row[k]): nb[t+u]=(sc+row[k],ch+[k])
        best=nb
    if T not in best: return None
    sc,ch=best[T]
    if lp.max(1).sum()-sc>maxloss: return None
    return ch
