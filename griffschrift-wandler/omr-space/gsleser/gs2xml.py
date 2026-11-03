# Griffschrift-Bild -> MusicXML (Normalnoten) + Griffe für den Wandler (Pins: Zeit, Ton, Reihe, Knopf, Balgrichtung)
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2,sys,json,pickle,collections,numpy as np,os as _os
from PIL import Image, ImageOps
sys.path.insert(0,HERE)
import gsread as G, detect as D
from align import chords_det
from assemble import button,pitch
from fuchsread import druck_bars
from rhythm import strip,beam_feats,dot_feat
from measures import barlines,rests,fit_measure,UNIT
import bassread as BR
_RD=HERE if _os.path.exists(HERE+'/rhythm.pkl') else _os.path.join(HERE,'..','gsleser')
RCLF=pickle.load(open(_RD+'/rhythm.pkl','rb'))
# zweites Notenwert-Modell, angelernt an echten Seiten (Philipps Übungsstücke, Fuchs-Seiten mit MIDI); beide Urteile werden gemischt
import gzip as _gz
_R2=_RD+'/rhythm2.pkl.gz'   # gepackt, damit die Datei klein bleibt
RCLF2=pickle.load(_gz.open(_R2,'rb')) if _os.path.exists(_R2) and not _os.environ.get('NORH2') else None
RW2=float(_os.environ.get('RW2','0.6'))
SHIFT={'GCFB':-2,'ADGC':0,'BEsAsDes':1,'HEAD':2,'CFBEs':3,'FBEsAsDes':-4,'FBEsAs':-4,'EADG':-5,'DGCF':-7}
def skew_angle(g):
    """Schräglage in Grad (Foto, schief eingelegter Scan): der Winkel, bei dem die Notenlinien am schärfsten waagrecht stehen."""
    sm=cv2.resize(g,None,fx=min(1.0,1200/max(g.shape)),fy=min(1.0,1200/max(g.shape)),interpolation=cv2.INTER_AREA)
    bw=(sm<160).astype(np.float32);h,w=bw.shape;c=(w/2,h/2)
    def score(a):
        r=cv2.warpAffine(bw,cv2.getRotationMatrix2D(c,a,1.0),(w,h),flags=cv2.INTER_NEAREST);p=r.sum(1);return float((p*p).sum())
    best=max(np.arange(-4,4.01,0.25),key=score)
    best=max(np.arange(best-0.25,best+0.26,0.05),key=score)
    return float(best)
# Triolen-Ziffer „3“ (16×20 Raster, aus Philipps Heft); Vergleich über Korrelation der Form
_T3=np.array([int(c) for c in bin(int('01fc03fe07fe0f1f1c0f1c0f001f003e03fc03f800fc003c001c001ce01ce01cf03cf87c7ff83fe0',16))[2:].zfill(320)],np.float32).reshape(20,16)
_T3=(_T3-_T3.mean())/_T3.std()
def triplet_marks(nol,ys,d,thr=0.7):
    """x-Lagen der Triolen-Ziffern 3 über oder unter dem System (nicht im System, dort steht die Taktangabe)."""
    a=max(0,int(ys[0]-6*d));b=min(nol.shape[0],int(ys[4]+6*d));sub=nol[a:b]
    n,lab,S,_=cv2.connectedComponentsWithStats((sub>0).astype(np.uint8));out=[]
    for i in range(1,n):
        x,y,w,h,_a=S[i]
        if not(0.9*d<=h<=1.6*d and 0.6*d<=w<=1.3*d): continue
        if abs(a+y+h/2-ys[2])<2.5*d: continue
        c=cv2.resize(((lab[y:y+h,x:x+w]==i)*255).astype(np.uint8),(16,20),interpolation=cv2.INTER_AREA).astype(np.float32)
        if c.std()==0: continue
        if float((((c-c.mean())/c.std())*_T3).mean())>thr: out.append(x+w/2)
    return out
_DBG={}
TRIP={3:2,4:4,6:4,9:6,12:8,18:12,24:16,36:24,48:32}   # Notenwert in einer Triole (zwei Drittel); 4 ist schon eine Triolenachtel
def dewarp(g,gs):
    """Gewölbte Seite (Foto eines Hefts): Notenlinien steigen und fallen über die Breite. Je waagrechtem Streifen der Seite wird die
    Verschiebung der Linien von Spalte zu Spalte gemessen und das Bild so verzogen, dass die Linien wieder waagrecht laufen."""
    h,w=g.shape;k=max(1.0,max(h,w)/1600);sm=cv2.resize(gs,None,fx=1/k,fy=1/k,interpolation=cv2.INTER_AREA);H,W=sm.shape
    bw=(sm<200).astype(np.uint8)*255;lin=cv2.morphologyEx(bw,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(max(15,W//40),1))).astype(np.float32)
    NS,NB=16,6;sw=W//NS;bh=H//NB;L=12;S=np.zeros((NB,NS),np.float32)
    for b in range(NB):
        y0=max(0,b*bh-L);y1=min(H,(b+1)*bh+L);prof=[lin[y0:y1,q*sw:(q+1)*sw].sum(1) for q in range(NS)]
        sh=[0.0]
        for q in range(NS-1):
            a_,b_=prof[q],prof[q+1]
            if a_.sum()<5*255 or b_.sum()<5*255: sh.append(0.0);continue
            best=max(range(-L,L+1),key=lambda t:float((a_[max(0,-t):len(a_)-max(0,t)]*b_[max(0,t):len(b_)-max(0,-t)]).sum()))
            sh.append(float(best))
        c=np.cumsum(sh);S[b]=c-c[NS//2]
    if np.abs(S).max()<1.5: return g
    S=cv2.GaussianBlur(S,(3,3),0.8)
    F=cv2.resize(S,(w,h),interpolation=cv2.INTER_LINEAR)*k
    mx,my=np.meshgrid(np.arange(w,dtype=np.float32),np.arange(h,dtype=np.float32))
    return cv2.remap(g,mx,my+F,cv2.INTER_LINEAR,borderValue=255)
def orient(g):
    """Hochkant fotografierte Seite (quer im Bild): mehr senkrechte als waagrechte lange Striche -> um 90° drehen"""
    sm=cv2.resize(g,None,fx=800/max(g.shape),fy=800/max(g.shape),interpolation=cv2.INTER_AREA)
    sm=cv2.divide(sm,cv2.GaussianBlur(sm,(0,0),10),scale=255);bw=(sm<190).astype(np.uint8)*255
    hz=cv2.morphologyEx(bw,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(40,1))).sum()
    vt=cv2.morphologyEx(bw,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(1,40))).sum()
    return cv2.rotate(g,cv2.ROTATE_90_CLOCKWISE) if vt>2*max(hz,1) else g
def upside_down(g,st):
    """Steht die Seite auf dem Kopf? Gezählt werden sicher erkannte Bassbuchstaben unter den Systemen, einmal so und einmal um 180° gedreht;
    auf dem Kopf stehende Buchstaben erkennt das Modell nicht."""
    def cnt(im):
        _,bw=cv2.threshold(im,200,255,cv2.THRESH_BINARY_INV);ss=G.staves(cv2.threshold(im,240,255,cv2.THRESH_BINARY_INV)[1]);n=0
        for si in range(len(ss)): n+=len(BR.read_letters(bw,ss,si,(ss[si][4]-ss[si][0])/4))
        return n
    a=cnt(g)
    if a>=3*len(st): return False
    r=cnt(cv2.rotate(g,cv2.ROTATE_180))
    return r>=max(8,1.5*a)
def load_gray(path,target_d=18.0,_turn=True):
    im=Image.open(path)
    try: im=ImageOps.exif_transpose(im)   # Handyfotos: Drehung aus den Bilddaten
    except Exception: pass
    im=im.convert('RGBA');bg=Image.new('RGBA',im.size,(255,255,255,255));bg.alpha_composite(im);g=np.array(bg.convert('L'))
    g=orient(g)
    a=skew_angle(g)
    if abs(a)>=0.2:   # schräge Vorlage gerade drehen, sonst findet die Liniensuche die Systeme nicht
        h,w=g.shape;g=cv2.warpAffine(g,cv2.getRotationMatrix2D((w/2,h/2),a,1.0),(w,h),flags=cv2.INTER_CUBIC,borderValue=255)
    hell=g[g>200].astype(np.float32)   # Foto erkennen: grauer Hintergrund oder Rauschen im hellen Bereich (robust gegen geglättete Kanten sauberer Bilder)
    foto=bool(np.median(g)<250 or (hell.size and float(np.median(np.abs(hell-np.median(hell))))>=2))
    # bei Fotos die Helligkeit schon vor der Liniensuche ausgleichen: auf grauem Hintergrund liegt sonst alles unter der Linienschwelle
    gs=cv2.divide(g,cv2.GaussianBlur(g,(0,0),max(8.0,min(g.shape)/60)),scale=255) if foto else g
    if foto and not _os.environ.get('NODEWARP'):
        g=dewarp(g,gs);gs=cv2.divide(g,cv2.GaussianBlur(g,(0,0),max(8.0,min(g.shape)/60)),scale=255)
    for f in (1,2,3):
        gg=cv2.resize(gs,None,fx=f,fy=f,interpolation=cv2.INTER_CUBIC) if f>1 else gs
        st=G.staves(cv2.threshold(gg,240,255,cv2.THRESH_BINARY_INV)[1])
        if st:
            d=float(np.median([np.mean(np.diff(s)) for s in st]))
            g=cv2.resize(g,None,fx=target_d/(d/f),fy=target_d/(d/f),interpolation=cv2.INTER_CUBIC)
            # Hintergrund ausgleichen, aber nur bei Fotos (grauer oder verrauschter Hintergrund); saubere Scans und Bilder bleiben unverändert
            if foto: g=cv2.divide(g,cv2.GaussianBlur(g,(0,0),20),scale=255)
            if foto and _os.environ.get('GAMMA'): g=(255*(g/255.0)**float(_os.environ['GAMMA'])).astype(np.uint8)   # blasse Striche kräftiger
            if _turn and foto:   # Foto auf dem Kopf: um 180° drehen und noch einmal (nur bei Fotos; saubere Bilder stehen richtig)
                st2=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1])
                if st2 and upside_down(g,st2): return _load_rot(path,target_d)
            return g
    return None
def _load_rot(path,target_d):
    import tempfile;im=Image.open(path)
    try: im=ImageOps.exif_transpose(im)
    except Exception: pass
    fn=tempfile.mktemp(suffix='.png');im.convert('RGB').rotate(180).save(fn)
    try: return load_gray(fn,target_d,_turn=False)
    finally:
        try: _os.remove(fn)
        except Exception: pass
def read(path,shift=0,work=None,table=None,bass=None):
    """shift: Stimmung wie im Wandler (A-D-G-C = 0, G-C-F-B = -2 …). bass: Liste, an die je Takt die gelesenen Bassbuchstaben (Lage im Takt in Zwölftel-Vierteln, Buchstabe) angehängt werden."""
    g=load_gray(path)
    if g is None: raise ValueError('Keine Notenlinien gefunden')
    _DBG['g']=g
    if work is None:
        import tempfile;work=_tmp=tempfile.mktemp(suffix='.png')
    cv2.imwrite(work,g);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV)
    st,bars_d=druck_bars(g);det=D.detect(work,1);C=chords_det(det)
    off=shift+2   # Tabelle in assemble.py gilt für G-C-F-B (shift -2)
    items=[]   # je System: Akkorde und Pausen
    feats=[]
    # Zeilenanfang je System (Taktangabe im ersten System nicht als Kopf lesen)
    x_start={};x_ext={};foto_like=bool(np.median(g)<253)
    for si,ys in enumerate(st):
        row=np.where(bw[int(round(ys[0]))]>0)[0];x_start[si]=float(row.min()) if len(row) else 0.0
        # Ausdehnung des Systems: Spalten, in denen (fast) alle fünf Linien Tinte haben; längstes Stück (Fotos: Seitenrand und Nachbarseite bleiben draußen)
        hl=cv2.morphologyEx(bw[max(0,int(ys[0])-8):int(ys[4])+9],cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(max(5,int(2*(ys[4]-ys[0])/4)),1)));o0=max(0,int(ys[0])-8)
        w_=max(1,int(0.25*(ys[4]-ys[0])/4));cols=sum((hl[max(0,int(round(y))-o0-w_):int(round(y))-o0+w_+1].max(0)>0).astype(np.int8) for y in ys)>=3
        cols=cv2.dilate(cols.astype(np.uint8)[None],np.ones((1,9),np.uint8))[0]>0
        runs=[];a_=None
        for xx,v in enumerate(list(cols)+[False]):
            if v and a_ is None: a_=xx
            if not v and a_ is not None: runs.append((a_,xx));a_=None
        dd=(ys[4]-ys[0])/4;mr=[]
        for r in runs:   # Lücken bis 8 Linienabstände überbrücken (blasse Linienstücke)
            if mr and r[0]-mr[-1][1]<8*dd: mr[-1]=(mr[-1][0],r[1])
            else: mr.append(r)
        runs=mr
        if runs:
            a_,b_=max(runs,key=lambda r:r[1]-r[0]);x_ext[si]=(float(a_),float(b_))
            if foto_like: x_start[si]=float(a_)
    for c in C:
        si=c[0]['st'];x=float(np.median([h['x'] for h in c]))
        if si==0 and x<x_start[0]+6*(st[0][4]-st[0][0])/4: continue
        if si in x_ext and not (x_ext[si][0]-0.5*(st[si][4]-st[si][0])/4<=x<=x_ext[si][1]+0.5*(st[si][4]-st[si][0])/4): continue   # außerhalb des Systems
        dr='D' if any(a-3<=x<=b+3 for a,b,_ in bars_d.get(si,[])) else 'Z'
        ks=[button(h['p'],h['t']) for h in c];ps=[pitch(r,n,dr,table) for r,n in ks]
        hs=[(h,k,p+off) for h,k,p in zip(c,ks,ps) if p]
        if not hs: continue
        ys=st[si];d=(ys[4]-ys[0])/4;nl=G.remove_lines(bw,ys,d)
        fill=float(np.mean([h['fill'] for h,_,_ in hs]))
        feats.append(np.concatenate([strip(nl,x,ys[2],d,fill),np.array(beam_feats(nl,x,[h['y'] for h,_,_ in hs],d)+[dot_feat(nl,x,[h['y'] for h,_,_ in hs],d)])*20]).astype(np.float32))
        items.append(dict(kind='n',st=si,x=x,dir=dr,heads=[(k[0],k[1],p) for _,k,p in hs],fi=len(feats)-1,hy=[(h['x'],h['y'],k[0]) for h,k,_ in hs]))
    _DBG['feats']=feats
    P=RCLF.predict_proba(np.array(feats)) if feats else np.zeros((0,len(RCLF.classes_)))
    P9=np.zeros((len(P),len(UNIT)));P9[:,RCLF.classes_]=P;lp=np.log(np.clip(P9,1e-6,1))
    if RCLF2 is not None and len(feats):
        Q9=np.zeros((len(P),len(UNIT)));Q9[:,RCLF2.classes_]=RCLF2.predict_proba(np.array(feats))
        lp=(1-RW2)*np.log(np.clip(P9,1e-4,1))+RW2*np.log(np.clip(Q9,1e-4,1))
    for it,row in zip(items,lp): it['lp']=row
    meas=[];medge=[];letters={}
    for si,ys in enumerate(st):
        d=(ys[4]-ys[0])/4;its=sorted([it for it in items if it['st']==si],key=lambda it:it['x'])
        hxy=[]
        for it in its:
            for r_,n_,p_ in it['heads']: hxy.append((it['x'],0))
        heads_xy=[(h['x'],h['y']) for c in C if c[0]['st']==si for h in c]
        bl=barlines(bw,ys,d,[x for x,_ in heads_xy]);nol=G.stems_off(G.remove_lines(bw,ys,d),d)
        rs=[dict(kind='r',st=si,x=x,u=u) for x,u in rests(bw,nol,ys,d,heads_xy,bl)]
        edges=[0]+bl+[1e9];allg=[]
        if bass is not None and not _os.environ.get('NOBASS'): letters[si]=BR.read_letters(bw,st,si,d)
        for a,b in zip(edges,edges[1:]):
            ev=sorted([it for it in its if a<it['x']<b]+[r for r in rs if a<r['x']<b],key=lambda e:e['x'])
            grp=[]
            for e in ev:
                if e['kind']=='n' and grp and grp[-1][0]['kind']=='n' and e['x']-grp[-1][-1]['x']<1.3*d: grp[-1].append(e)
                else: grp.append([e])
            if grp: meas.append((si,grp));allg+=grp;medge.append((a,b,d))
        # Triolen: die Gruppe unter/über der Ziffer 3 und ihre beiden Nachbarn
        if allg and not _os.environ.get('NOTRIP'):
            gx=[float(np.mean([e['x'] for e in g_])) for g_ in allg]
            for mx in triplet_marks(nol,ys,d):
                j=int(np.argmin([abs(x_-mx) for x_ in gx]))
                if abs(gx[j]-mx)>3*d: continue
                for g_ in allg[max(0,j-1):j+2]: g_[0]['trip']=True
    def gl(gr): return sum(e['lp'] for e in gr)
    def uu(gr,u): return TRIP.get(u,u) if gr[0].get('trip') else u
    def length(m,ch=None):
        t=0;k=0
        for q,gr in enumerate(m):
            if gr[0]['kind']=='r': t+=uu(gr,gr[0]['u'])
            else: t+=uu(gr,UNIT[ch[k] if ch else int(np.argmax(gl(gr)))]);k+=1
        return t
    cg=[c for c in collections.Counter(length(m) for _,m in meas).most_common() if c[0] in (18,24,36,48,72)];Lg=cg[0][0] if cg else 36
    # Taktlänge je System (Taktwechsel zwischen Teilen): häufigste Länge im System, wenn mindestens zwei Takte sie haben
    Ls={}
    for si in set(si for si,_ in meas):
        cs=[c for c in collections.Counter(length(m) for s_,m in meas if s_==si).most_common() if c[0] in (18,24,36,48,72)]
        Ls[si]=cs[0][0] if cs and cs[0][1]>=2 else Lg   # nur übliche Taktarten (3/8, 2/4, 3/4, 4/4, 6/4)
    out=[];lens=[]
    for mi,(si,m) in enumerate(meas):
        Lm=Ls[si]
        notes=[q for q,gr in enumerate(m) if gr[0]['kind']=='n'];ch=[int(np.argmax(gl(m[q]))) for q in notes]
        if length(m)!=Lm and not (mi==0 and length(m)<Lm):
            PG=np.array([np.exp(gl(m[q])-gl(m[q]).max()) for q in notes]);PG=PG/PG.sum(1,keepdims=True) if len(PG) else PG
            f=fit_measure(PG,[uu(gr,gr[0]['u']) for gr in m if gr[0]['kind']=='r'],Lm,units=[[uu(m[q],u) for u in UNIT] for q in notes])
            if f: ch=f
        evs=[];k=0
        for gr in m:
            if gr[0]['kind']=='r': evs.append(dict(rest=True,dbg=(gr[0]['st'],[],gr[0]['x']),u=uu(gr,gr[0]['u']),trip=bool(gr[0].get('trip'))))
            else:
                heads=sorted({h for e in gr for h in e['heads']},key=lambda h:-h[2]);evs.append(dict(rest=False,dbg=(gr[0]['st'],[q for e in gr for q in e['hy']],gr[0]['x']),fis=[e['fi'] for e in gr],mtrip=bool(gr[0].get('trip')),u=uu(gr,UNIT[ch[k]]),trip=bool(gr[0].get('trip')) or UNIT[ch[k]]==4,dir=gr[0]['dir'],heads=heads));k+=1
        out.append(evs);lens.append(Lm)
        if bass is not None:   # Bassbuchstaben des Takts: Lage über die Einsätze der Noten und Pausen (nächster Einsatz oder dazwischen auf Achtel gerundet)
            a,b,d=medge[mi];xs=[];pos=0
            for gr,e in zip(m,evs): xs.append((float(np.mean([q['x'] for q in gr])),pos));pos+=e['u']
            bl_=[];ins=[(x,c if isinstance(c,tuple) else str(c)) for x,c in letters.get(si,[]) if a<x<b];T=pos if 0<pos<Lm else Lm
            if ins and not _os.environ.get('BASS_ALL') and T%len(ins)==0 and T//len(ins) in (6,12,18,24,36):   # gleichmäßige Begleitung (häufigster Fall): Buchstaben auf gleich lange Schläge
                bass.append([(i*(T//len(ins)),c) for i,(x,c) in enumerate(ins)]);continue
            for x,c in ins:
                j=int(np.argmin([abs(x-gx) for gx,_ in xs]))
                if abs(x-xs[j][0])<0.9*d: q=xs[j][1]
                else:
                    L=[t for t in xs if t[0]<x];R=[t for t in xs if t[0]>x]
                    l=L[-1] if L else (a,0);r=R[0] if R else (min(b,xs[-1][0]+3*d),Lm)
                    q=l[1]+(x-l[0])/max(1e-6,r[0]-l[0])*(r[1]-l[1]);q=int(round(q/6)*6)
                q=max(0,min(Lm-1,q))
                if _os.environ.get('BASS_ALL') or not any(q==q2 for q2,_ in bl_): bl_.append((q,c))
            bass.append(bl_)
    try: _os.remove(work)
    except Exception: pass
    return out,lens
def parse_table(s):
    """Diskantbelegung vom Wandler als JSON-Text ({"druck":{"1":[...],...},"zug":{...}}); leer oder fehlerhaft -> None (446 Normal)."""
    if not s: return None
    try:
        t=json.loads(s) if isinstance(s,str) else s
        out={}
        for dr in ('druck','zug'):
            out[dr]={int(r):[int(v) for v in a][:40] for r,a in t[dr].items() if int(r) in (1,2,3,4)}
        return out if all(len(out[dr])==4 for dr in out) else None
    except Exception: return None
def convert(paths,shift=0,title='Gescannte Griffschrift',table=None):
    """Eine oder mehrere Seiten (Bilder) -> (MusicXML, Pins); Takte der Seiten hintereinander. table: Diskantbelegung des Wandlers (sonst 446 Normal)."""
    allm=[];Lm=None;bass=[]
    for p in paths:
        m,L=read(p,shift,table=table,bass=bass);allm+=m;Lm=(Lm or [])+L
    if not allm: raise ValueError('Keine Griffschrift gefunden')
    return to_xml(allm,Lm,title,shift,bass)
STEP=['C','C','D','D','E','F','F','G','G','A','A','B'];ALT=[0,1,0,1,0,0,1,0,1,0,1,0]
TRIPTYPE={2:('16th',0),4:('eighth',0),6:('eighth',1),8:('quarter',0),12:('quarter',1),16:('half',0),24:('half',1),32:('whole',0)}   # Notenart innerhalb einer Triole
TYPE={3:('16th',0),4:('eighth',0),6:('eighth',0),9:('eighth',1),12:('quarter',0),18:('quarter',1),24:('half',0),36:('half',1),48:('whole',0)}
def to_xml(meas,lens,title='Gescannte Griffschrift',shift=0,bass=None):
    if not isinstance(lens,list): lens=[lens]*len(meas)
    def ts(L): return (L//12,4) if L%12==0 else (L//6,8)
    x=['<?xml version="1.0" encoding="UTF-8"?>','<score-partwise version="3.1">',f'<work><work-title>{title}</work-title></work>',
       '<part-list><score-part id="P1"><part-name>Melodie</part-name></score-part></part-list>','<part id="P1">']
    pins=[];t0=0.0;bpins=[]
    for mi,evs in enumerate(meas):
        Lm=lens[mi];beats,btype=ts(Lm)
        imp=' implicit="yes"' if mi==0 and sum(e['u'] for e in evs)<Lm else ''
        if mi==0: first_imp=bool(imp)
        x.append(f'<measure number="{mi if first_imp else mi+1}"{imp}>')
        if mi==0: x.append(f'<attributes><divisions>12</divisions><key><fifths>0</fifths></key><time><beats>{beats}</beats><beat-type>{btype}</beat-type></time><clef><sign>G</sign><line>2</line></clef></attributes>')
        elif lens[mi]!=lens[mi-1]: x.append(f'<attributes><time><beats>{beats}</beats><beat-type>{btype}</beat-type></time></attributes>')
        pos=0
        for e in evs:
            tr=e.get('trip',e['u']==4);ty,dot=(TRIPTYPE if tr else TYPE).get(e['u'],TYPE.get(e['u'],('quarter',0)));tm='<time-modification><actual-notes>3</actual-notes><normal-notes>2</normal-notes></time-modification>' if tr else ''
            if e['rest']: x.append(f'<note><rest/><duration>{e["u"]}</duration><type>{ty}</type>{"<dot/>"*dot}{tm}</note>')
            else:
                for q,(r,n,p) in enumerate(e['heads']):
                    al=f'<alter>1</alter>' if ALT[p%12] else ''
                    x.append(f'<note>{"<chord/>" if q else ""}<pitch><step>{STEP[p%12]}</step>{al}<octave>{p//12-1}</octave></pitch><duration>{e["u"]}</duration><type>{ty}</type>{"<dot/>"*dot}{tm}</note>')
                    pins.append([round(t0+pos/12,4),p,r,n,e['dir']])
            pos+=e['u']
        if bass and mi<len(bass): bpins+=[[round(t0+q/12,4),c] for q,c in bass[mi]]
        x.append('</measure>');t0+=pos/12
    x+=['</part>','</score-partwise>']
    P={'sh':shift,'n':pins}
    if bpins: P['b']=bpins   # Bassbuchstaben wie in der Griffschrift (Zeit in Vierteln, Name)
    return '\n'.join(x),P
if __name__=='__main__':
    sh=SHIFT[sys.argv[2] if len(sys.argv)>2 else 'ADGC'];meas,Lm=read(sys.argv[1],sh)
    xml,pins=to_xml(meas,Lm,shift=sh)
    open(sys.argv[3] if len(sys.argv)>3 else HERE+'/out.musicxml','w').write(xml)
    json.dump(pins,open((sys.argv[3] if len(sys.argv)>3 else HERE+'/out.musicxml')+'.pins.json','w'))
    print('Takte',len(meas),'Taktlängen',sorted(set(l/12 for l in Lm)),'Viertel; Noten',sum(1 for m in meas for e in m if not e['rest']),'Pausen',sum(1 for m in meas for e in m if e['rest']))
