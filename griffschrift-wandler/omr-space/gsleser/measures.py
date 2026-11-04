# Taktstriche und Pausen auf Fuchs-Seiten
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2,sys,numpy as np
sys.path.insert(0,HERE);import gsread as G
def straight_bar(bw,top,bot,a,b,d):
    """Ist der Tintenstreifen a..b zwischen top und bot ein gerader Strich (gleich breit über die Höhe)?
    Ein dichter Akkord füllt die Höhe auch, ist aber an den Köpfen viel breiter als sein Hals."""
    x0=max(0,int(a-0.8*d));x1=min(bw.shape[1],int(b+0.8*d))
    w=(bw[top:bot+1,x0:x1]>0).sum(1)
    return float(np.percentile(w,80))<=(b-a)+0.35*d
def local_shift(bw,ys,d):
    """Je Spalte: um wie viele Pixel die fünf Linien dort höher (-) oder tiefer (+) liegen als ys (Fotos: gewölbte Seiten)"""
    W=bw.shape[1];r=int(d*0.6);top=max(0,int(ys[0])-r);bot=min(bw.shape[0],int(ys[4])+r+1)
    hl=cv2.morphologyEx(bw[top:bot],cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(max(5,int(1.5*d)),1)))>0
    step=max(4,int(d));xs=list(range(0,W,step));sv=[]
    for x0 in xs:
        blk=hl[:,max(0,x0-int(4*d)):x0+int(4*d)+1].mean(1)
        best=max(range(-r,r+1),key=lambda k:sum(blk[int(round(y))-top+k] for y in ys if 0<=int(round(y))-top+k<len(blk)))
        sv.append(best)
    sv=[int(np.median(sv[max(0,i-3):i+4])) for i in range(len(sv))]   # Ausreißer (Bögen, Balken) glätten
    return np.round(np.interp(np.arange(W),xs,sv)).astype(int)
def barlines(bw,ys,d,heads_x):
    """Senkrechte Striche genau von der obersten bis zur untersten Linie, ohne Kopf daneben."""
    top,bot=int(ys[0]),int(ys[4]);W=bw.shape[1];cols=[]
    inside=(bw[top:bot+1]>0).mean(0)
    above=(bw[max(0,top-int(.8*d)):top-2]>0).mean(0);below=(bw[bot+3:bot+int(.8*d)]>0).mean(0)
    sh=local_shift(bw,ys,d) if G.FOTO_LINES and not _os.environ.get('GS_NOBLSHIFT') else None
    for x in range(W):
        if inside[x]>0.92 and above[x]<0.3 and below[x]<0.3: cols.append(x);continue
        if sh is None: continue
        # Fotos: System an dieser Stelle etwas höher oder tiefer (gewölbte Seite); Linien dort neu ansetzen
        t,b_=top+sh[x],bot+sh[x]
        if t-int(.8*d)<0 or b_+int(.8*d)>=bw.shape[0]: continue
        i_=(bw[t:b_+1,x]>0).mean();a_=(bw[t-int(.8*d):t-2,x]>0).mean();u_=(bw[b_+3:b_+int(.8*d),x]>0).mean()
        # Wiederholungs- oder Klammerhaken (1., 2.) kommen von oben an den Taktstrich: dann genügt eine saubere Unterseite
        if i_>0.92 and (a_<0.3 and u_<0.3 or i_>=0.97 and min(a_,u_)<0.1): cols.append(x)
    out=[]
    for x in cols:
        if out and x-out[-1][-1]<=2: out[-1].append(x)
        else: out.append([x])
    if G.FOTO_LINES:   # dicke Striche (Wiederholung, Schluss) auch mit Zeichen darüber oder darunter, etwa Segno-Strichen bei Rosenzopf
        th=np.flatnonzero(np.diff(np.concatenate([[0],(inside>0.95).astype(np.int8),[0]])))
        for a,b in zip(th[::2],th[1::2]):
            if 0.3*d<=b-a<=0.9*d and not any(a-3<=q[0]<=b+3 for q in out) and straight_bar(bw,top,bot,a,b,d): out.append(list(range(a,b)))
        out.sort(key=lambda q:q[0])
    res=[]
    for g in out:
        x=float(np.mean(g))
        if len(g)>(0.9 if G.FOTO_LINES else 0.6)*d: continue           # zu breit: Balken o. ä. (auf Fotos sind dicke Wiederholungsstriche breiter)
        if len(g)<0.3*d and any(abs(x-h)<0.9*d for h in heads_x): continue   # Hals eines Kopfes (dicke Striche wie beim Wiederholungszeichen sind keine Hälse; ein Punkt daneben wird manchmal als Kopf gelesen)
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
        if (2.5<=h<=3.5 and 0.8<=w<=1.4 or G.FOTO_LINES and 2.0<=h<=3.5 and 0.6<=w<=1.4) and 0.38<=fillr<=0.65 and abs(c)<1: u=12            # Viertelpause (auf Fotos etwa bei Rosenzopf kleiner)
        elif 1.5<=h<=2.2 and (0.9<=w<=1.4 and 0.35<=fillr<=0.6 or G.FOTO_LINES and _os.environ.get('GS_R8','0')=='1' and 0.7<=w<=2.1 and 0.15<=fillr<0.35) and abs(c)<0.8: u=6        # Achtelpause (auf Fotos dünner gezeichnet, mit Linienresten breiter)
        elif 0.5<=h<=0.95 and 1.05<=w<=1.7 and fillr>0.6 and -1.1<c<0.1: u=24 if c>-0.5 else 48   # Halbe/Ganze Pause
        else: continue
        # Punkt rechts daneben verlängert um die Hälfte
        if any(0.6*d<px-cx<2.2*d and abs(py-cy)<0.8*d and 0.1<pa<0.3 and pw<0.7 for px,py,pw,ph,pa in comps): u=u*3//2
        out.append((cx,u))
    return out
def fit_measure(P,rest_units,L,maxloss=6.0,units=None):
    """P: Wahrscheinlichkeiten je Akkord über UNIT; wählt die Werte mit größter Wahrscheinlichkeit, deren Summe (mit Pausen) L ergibt."""
    T=L-sum(rest_units)
    if T<=0 or not len(P): return None
    lp=np.log(np.clip(P,1e-6,1));best={0:(0.0,[])}
    for row in lp:
        nb={}
        for t,(sc,ch) in best.items():
            for k,u in enumerate(UNIT if units is None else units[len(ch)]):
                if t+u<=T and (t+u not in nb or nb[t+u][0]<sc+row[k]): nb[t+u]=(sc+row[k],ch+[k])
        best=nb
    if T not in best: return None
    sc,ch=best[T]
    if lp.max(1).sum()-sc>maxloss: return None
    return ch
