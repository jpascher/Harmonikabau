# Griffschrift-Bild -> MusicXML (Normalnoten) + Griffe für den Wandler (Pins: Zeit, Ton, Reihe, Knopf, Balgrichtung)
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2,sys,json,pickle,collections,re,numpy as np,os as _os
from PIL import Image, ImageOps
sys.path.insert(0,HERE)
import gsread as G, detect as D
from align import chords_det
from assemble import button,pitch
from fuchsread import druck_bars
from rhythm import strip,beam_feats,dot_feat
from measures import barlines,rests,fit_measure,UNIT,straight_bar
import bassread as BR
_RD=HERE if _os.path.exists(HERE+'/rhythm.pkl') else _os.path.join(HERE,'..','gsleser')
RCLF=pickle.load(open(_RD+'/rhythm.pkl','rb'))
# zweites Notenwert-Modell, angelernt an echten Seiten (Philipps Übungsstücke, Fuchs-Seiten mit MIDI); beide Urteile werden gemischt
import gzip as _gz
_R2=_os.environ.get('GS_RH2',_RD+'/rhythm2.pkl.gz')   # gepackt, damit die Datei klein bleibt
RCLF2=pickle.load(_gz.open(_R2,'rb')) if _os.path.exists(_R2) and not _os.environ.get('NORH2') else None
RW2=float(_os.environ.get('RW2','0.6'))
GERUEST=_os.environ.get('GS_GERUEST','0')=='1'   # Gerüst (Takte, Notenwerte, Pausen) aus homr; Versuch, voreingestellt aus
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
    NS,NB=int(_os.environ.get('DW_NS','16')),int(_os.environ.get('DW_NB','6'));sw=W//NS;bh=H//NB;L=12;S=np.zeros((NB,NS),np.float32)
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
    if vt>2*max(hz,1):
        # senkrechte Zierrahmen oder Randlinien: stehen schon so Notensysteme im Bild, nicht drehen
        if len(G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]))>=2: return g
        return cv2.rotate(g,cv2.ROTATE_90_COUNTERCLOCKWISE if quer_richtung(g)==270 else cv2.ROTATE_90_CLOCKWISE)
    return g
def quer_richtung(g):
    """Quer liegende Seite: in welche Richtung drehen? Beide Drehungen werden mit der Texterkennung (Tesseract) gelesen; die Richtung,
    in der mehr Wörter mit mindestens drei Buchstaben erkannt werden (Titel, Bearbeiter, Teilbezeichnungen), gewinnt. Gefunden am
    11.11.2026 bei den Steinkellner-PDFs, die quer und im Uhrzeigersinn gedreht im PDF liegen; bisher wurde immer im Uhrzeigersinn
    gedreht, die Seite stand dann auf dem Kopf. Ohne Tesseract oder bei Gleichstand: 90 (bisheriges Verhalten)."""
    try:
        import pytesseract,re
        sm=cv2.resize(g,None,fx=min(1.0,2400/max(g.shape)),fy=min(1.0,2400/max(g.shape)),interpolation=cv2.INTER_AREA)
        def n(im): return len(re.findall(r'[A-Za-zÄÖÜäöüß]{3,}',pytesseract.image_to_string(im,lang='deu',config='--psm 11')))
        a=n(cv2.rotate(sm,cv2.ROTATE_90_CLOCKWISE));b=n(cv2.rotate(sm,cv2.ROTATE_90_COUNTERCLOCKWISE))
        return 270 if b>a else 90
    except Exception: return 90
def upside_down(g,st):
    """Steht die Seite auf dem Kopf? Gezählt werden sicher erkannte Bassbuchstaben unter den Systemen, einmal so und einmal um 180° gedreht;
    auf dem Kopf stehende Buchstaben erkennt das Modell nicht."""
    def cnt(im):
        _,bw=cv2.threshold(im,200,255,cv2.THRESH_BINARY_INV);ss=G.staves(cv2.threshold(im,240,255,cv2.THRESH_BINARY_INV)[1]);n=0
        for si in range(len(ss)): n+=len(BR.read_letters(bw,ss,si,(ss[si][4]-ss[si][0])/4))
        return n
    a=cnt(g)
    if a>=8*len(st): return False   # viele sicher gelesene Buchstaben: steht richtig (bei weniger immer auch gedreht zählen; auf dem Kopf liest das Modell auch einige)
    r=cnt(cv2.rotate(g,cv2.ROTATE_180))
    if max(a,r)>=len(st): return r>=max(8,1.5*a)
    # kaum Bassbuchstaben unter den Systemen (etwa Steinkellner: Buchstaben über dem System): der Titel steht oben;
    # die größte Schrift außerhalb der Systeme nur unterhalb der Systeme heißt: Seite liegt auf dem Kopf
    if _os.environ.get('GS_TITLETURN','0')!='1' or not st: return r>=max(8,1.5*a)
    d=float(np.median([(s_[4]-s_[0])/4 for s_ in st]));_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV)
    def big(y0,y1):
        y0=max(0,int(y0));y1=min(bw.shape[0],int(y1))
        if y1-y0<2*d: return 0.0
        n,_,S,_=cv2.connectedComponentsWithStats((bw[y0:y1]>0).astype(np.uint8))
        hs=[S[i][3] for i in range(1,n) if 1.8*d<=S[i][3]<=6*d and S[i][2]<=6*d and S[i][4]>0.25*S[i][2]*S[i][3]]
        return float(np.percentile(hs,90)) if len(hs)>=3 else 0.0
    up=big(0,st[0][0]-1.5*d);dn=big(st[-1][4]+6*d,bw.shape[0])
    return dn>=2.2*d and dn>=1.3*up
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
    G.FOTO_LINES=foto
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
            if _turn and (foto or _os.environ.get('GS_TURNALL','0')=='1'):   # Foto auf dem Kopf (GS_TURNALL=1: auch saubere Scans prüfen, Versuch): um 180° drehen und noch einmal (nur bei Fotos; saubere Bilder stehen richtig)
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
# Schreibweisen: je eigene Regeln. „auto“ nimmt die allgemeinen Regeln; erkennt die Texterkennung einen Bearbeiter, wird dessen Schreibweise gewählt.
#   xmark: vorgesetztes Kreuz = nur Hinweis auf die Kreuzreihe (kein eigener Ton); druck: Druckstriche dick/dünn/automatisch
STYLE_RULES={
 'auto':        dict(xmark=True,  druck='auto'),
 'fuchs':       dict(xmark=False, druck='dick'),
 'michlbauer':  dict(xmark=True,  druck='auto'),
 'rosenzopf':   dict(xmark=True,  druck='duenn'),
 'margetich':   dict(xmark=False, druck='auto'),
 'bayerisch':   dict(xmark=False, druck='auto', kreuzrund=True),   # Kreuze für homr rund: Schaborak/Kiermeier 11.11.2026 besser
 'handschrift': dict(xmark=False, druck='auto', kreuzrund=True),   # Handschriften: ebenso (Die alte Lena 5 statt 15 offene Takte)
 'steinkellner':dict(xmark=False, druck='auto'),
 'pacher':      dict(xmark=False, druck='auto'),   # Günther Pacher (MC-Score): Kreuze vor runden Köpfen sind eigene Köpfe der 3. Reihe
 'klausner':    dict(xmark=False, druck='auto'),   # Hubert Klausner und Stefan Sirtl (Finale): ebenso
}
STYLE=['auto']
def set_style(name):
    name=(name or 'auto').strip().lower()
    STYLE[0]=name if name in STYLE_RULES else 'auto'
    import fuchsread as _F;_F.DRUCK_MODE[0]=STYLE_RULES[STYLE[0]]['druck']
    return STYLE[0]
# Texterkennung (Tesseract, deutsch): Titel, Bearbeiter, Teilbezeichnungen; nur Plausibles wird übernommen
TEIL_WOERTER=['Einleitung','Vorspiel','Tanz','Trio','Fine','D.S.','D.C.','Coda','Zwischenspiel','Schluss','Ende']
BEARB_STIL=[('rosenzopf','rosenzopf'),('michlbauer','michlbauer'),('wachtberg','michlbauer'),('fuchs','fuchs'),('steinkellner','steinkellner'),
            ('margetich','margetich'),('kiener','bayerisch'),('kiermeier','bayerisch'),('heimatpflege','bayerisch'),('knöpferl','bayerisch'),('schaborak','bayerisch'),
            ('pacher','pacher'),('klausner','klausner'),('sirtl','klausner')]
def _plaus(t):
    t=t.strip(" .,:;-–|_'\"")
    if len(t)<3: return ''
    let=sum(c.isalpha() for c in t);voc=sum(c.lower() in 'aeiouäöüy' for c in t)
    return t if let>=0.6*len(t.replace(' ','')) and voc>=1 and let>=3 else ''
def page_text(g,st):
    """Text oberhalb des ersten Systems (Titel, Bearbeiter) und Teilbezeichnungen auf der Seite. Ohne Tesseract: leer."""
    out=dict(titel='',bearbeiter='',teile=[],stil='')
    if _os.environ.get('GS_OCR','1')!='1': return out
    try: import pytesseract
    except Exception: return out
    try:
        d=(st[0][4]-st[0][0])/4 if st else 18.0;top=int(max(0,(st[0][0]-1.5*d) if st else g.shape[0]*0.2))
        lines={};words_all=[]
        if top>2*d:
            D_=pytesseract.image_to_data(g[:top],lang='deu',output_type=pytesseract.Output.DICT,config='--psm 3')
            for i,w in enumerate(D_['text']):
                w=w.strip()
                if not w: continue
                try: cf=float(D_['conf'][i])
                except Exception: cf=-1
                k=(D_['block_num'][i],D_['par_num'][i],D_['line_num'][i]);lines.setdefault(k,[]).append((w,cf,D_['height'][i],D_['left'][i]))
        L=[]
        for k,ws in lines.items():
            good=[w for w in ws if w[1]>=55]
            if not good: continue
            txt=' '.join(w[0] for w in sorted(good,key=lambda w:w[3]));h=float(np.median([w[2] for w in good]));L.append((txt,h,float(np.mean([w[1] for w in good]))))
        for txt,h,cf in L:
            low=txt.lower()
            if any(k in low for k in ('bearb','eingerichtet','satz:','arr.','gesetzt')) and not out['bearbeiter']:
                b_=re.sub(r'^.*?(bearbeitung|bearb\.?|eingerichtet von|eingerichtet für|satz|arr\.?)\s*:?\s*','',txt,flags=re.I).strip()
                out['bearbeiter']='' if re.search(r'eingerichtet\s+für',txt,flags=re.I) else (_plaus(b_) or '')   # „eingerichtet für die Steirische“ nennt keine Person
        cand=[(h,txt) for txt,h,cf in L if cf>=65 and not any(k in txt.lower() for k in ('bearb','eingerichtet','satz:','musik:','volksweise','komponist'))]
        if cand:
            h,txt=max(cand);txt=re.sub(r'^\s*(nr\.?\s*)?\d+\s*[.)]?\s*','',txt,flags=re.I).strip()
            if h>=1.2*d: out['titel']=_plaus(txt)
        # Teilbezeichnungen: ganze Seite, nur bekannte Wörter
        T=pytesseract.image_to_string(g,lang='deu',config='--psm 11')
        for w in TEIL_WOERTER:
            if re.search(r'(?<![A-Za-zÄÖÜäöü])'+re.escape(w)+r'(?![A-Za-zÄÖÜäöü])',T,flags=re.I) and w not in out['teile']: out['teile'].append(w)
        alltxt=(out['bearbeiter']+' '+' '.join(t for t,_,_ in L)).lower()
        for k,s_ in BEARB_STIL:
            if k in alltxt: out['stil']=s_;break
    except Exception: pass
    return out
XDROP=[]   # weggelassene vorgesetzte Kreuze (für die Pausensuche: ihre Tinte ist keine Pause)
def cross_marks(det):
    """Vorgesetzte Kreuze (Rosenzopf, Michlbauer): ein Kreuz knapp links vor einem runden Kopf in derselben Lage ist keine eigene Note,
    sondern zeigt nur die Reihe an. Das Kreuz fällt weg, der Kopf dahinter gilt als Kreuzkopf (3./4. Reihe). GS_XMARK=0: aus"""
    XDROP.clear()
    if _os.environ.get('GS_XMARK','1')!='1' or not STYLE_RULES[STYLE[0]]['xmark']: return det
    drop=set()
    for i,h in enumerate(det):
        if h['t']!='x': continue
        best=None
        for j,o in enumerate(det):
            if o['t']!='o' or o['st']!=h['st'] or o['p']!=h['p']: continue
            dx=o['x']-h['x']
            if 0.5*h['d']<dx<float(_os.environ.get('GS_XMARKD','2.4'))*h['d'] and (best is None or dx<best[0]): best=(dx,j)
        if best: drop.add(i);det[best[1]]=dict(det[best[1]],t='x',xmark=True);XDROP.append(h)
    # übrige Kreuze derselben Spalte ohne runden Kopf in ihrer Lage (Kopf übersehen): als Kreuzkopf in die Spalte des Akkords rücken
    W=float(_os.environ.get('GS_XMARKD','2.4'))
    for i,h in enumerate(det):
        if i in drop or h['t']!='x' or h.get('xmark'): continue
        cand=[o for o in det if o.get('xmark') and o['st']==h['st'] and 0.5*h['d']<o['x']-h['x']<W*h['d'] and abs(o['p']-h['p'])<=3]
        if cand: det[i]=dict(h,x=min(cand,key=lambda o:o['x']-h['x'])['x'],xmark=True)
    return [h for i,h in enumerate(det) if i not in drop]
def prefix_marks(det,bw,st):
    """Vorgesetzte Kreuze, die der Kopferkenner nicht als Kopf meldet (Rosenzopf: kräftiges Doppelkreuz knapp links vor dem runden Kopf).
    Ein runder Kopf gilt als Kreuzkopf, wenn 0,9–2,1 Linienabstände links in derselben Höhe ein Fleck in Kopfgröße liegt, der deutlich
    weniger gefüllt ist als der Kopf selbst und eine Kreuzform hat (beide Diagonalen tragen Tinte). Ein dort gemeldeter runder Kopf
    derselben Lage fällt weg. GS_XPREFIX=0: aus"""
    if _os.environ.get('GS_XPREFIX','1')!='1' or not STYLE_RULES[STYLE[0]]['xmark'] or not det: return det
    nl={};drop=set();out=[]
    def box(im,x,y,hw,hh):
        w_=im[max(0,int(y-hh)):int(y+hh)+1,max(0,int(x-hw)):int(x+hw)+1];return (w_>0).mean() if w_.size else 0
    def diag(im,x,y,r):
        c=0
        for t in np.linspace(-0.8,0.8,9):
            for sx,sy in ((1,1),(1,-1)):
                px,py=int(round(x+sx*t*r)),int(round(y+sy*t*r))
                if 0<=py<im.shape[0] and 0<=px<im.shape[1] and im[py,px]: c+=1
        return c/18
    nl0={}
    def sides(im,x,y,d):   # Seitenmitten des Flecks: beim Kreuz leer (Arme gehen über die Ecken), beim vollen Kopf Tinte
        c=0
        for px,py in ((x-0.42*d,y),(x+0.42*d,y),(x,y-0.36*d),(x,y+0.36*d)):
            w_=im[max(0,int(py)-1):int(py)+2,max(0,int(px)-1):int(px)+2]
            c+=bool(w_.size and (w_>0).mean()>0.5)
        return c
    for i,h in enumerate(det):
        if h['t']!='o' or h.get('xmark') or h['st']>=len(st): continue
        si=h['st'];ys=st[si];d=(ys[4]-ys[0])/4
        if si not in nl: nl0[si]=G.remove_lines(bw,ys,d);nl[si]=cv2.dilate(nl0[si],np.ones((3,3),np.uint8))
        im=nl[si];fh=box(im,h['x'],h['y'],0.45*d,0.32*d)
        best=None
        for dx in np.arange(0.9*d,2.15*d,max(1.0,d/9)):
            x=h['x']-dx;f=box(im,x,h['y'],0.45*d,0.32*d)
            kreuz=G.FOTO_LINES and sides(nl0[si],x,h['y'],d)<=1 and f<=min(0.88,fh-0.06) and dx>=1.0*d and not any(o is not h and o['st']==si and o['t']=='o' and abs(o['x']-x)<0.5*d and abs(o['y']-h['y'])<0.6*d for o in det)   # Kreuzform deutlich: auch dicht am Kopf, aber nicht dort, wo ein anderer Kopf gemeldet ist (dichte Akkordfolgen)
            if (0.3<=f<=min(0.68,fh-0.1) or kreuz) and diag(im,x,h['y'],0.45*d)>=0.75:
                # rechts davon bis zum Kopf eine Lücke (sonst ist es der Kopf selbst oder ein Balken)
                gap=box(im,(x+h['x'])/2,h['y'],0.08*d,0.25*d)
                lo_,hi_=max(0,int(x-0.8*d)),int(x+0.8*d)
                bar_=((bw[int(ys[0]):int(ys[4])+1,lo_:hi_]>0).mean(0)>0.9).any() if hi_>lo_ else False   # Taktstrich/Wiederholung daneben: Punkte, kein Kreuz
                lm_=box(im,x-0.7*d,h['y'],0.1*d,0.25*d)   # links vom Kreuz frei (dichte Akkordfolgen haben dort den nächsten Kopf)
                if (kreuz or gap<0.5 and lm_<0.5) and not bar_ and (best is None or f<best[1]): best=(x,f)
        if best:
            det[i]=dict(h,t='x',xmark=True)
            for j,o in enumerate(det):
                if j!=i and o['st']==si and o['p']==h['p'] and abs(o['x']-best[0])<0.6*d: drop.add(j)
    return [h for k,h in enumerate(det) if k not in drop]
def hollow_stack(det,bw,st):
    """Fotos: hohle Köpfe übereinander (Halbe-Akkorde in Terzen, die Köpfe berühren sich) meldet der Kopferkenner oft nur einmal.
    Über und unter einem gemeldeten hohlen Kopf wird in Terzabstand (zwei Lagen) mit der Vorlage des hohlen Kopfes nachgesehen;
    passt sie deutlich (und besser als ein voller Kopf), kommt dort ein hohler Kopf dazu. GS_HSTACK=0: aus"""
    if _os.environ.get('GS_HSTACK','1')!='1' or not G.FOTO_LINES or not det: return det
    nl={};add=[]
    for h in det:
        if h['t']!='o' or h['fill'] or h['st']>=len(st): continue
        si=h['st'];ys=st[si];d=(ys[4]-ys[0])/4;mid=ys[2]
        if si not in nl: nl[si]=(G.remove_lines(bw,ys,d),G.templates(d))
        im,T=nl[si]
        def good(x,y,key):
            t=T[key];hh,ww=t.shape;y0,x0=int(round(y))-hh//2,int(round(x))-ww//2;b_=0.0;tt=t.astype(np.float32)
            for dy in (-2,0,2):
                for dx in (-2,0,2):
                    win=im[max(0,y0+dy):y0+dy+hh,max(0,x0+dx):x0+dx+ww].astype(np.float32)
                    if win.shape!=t.shape: continue
                    b_=max(b_,float((win*tt).sum()/(np.sqrt((win*win).sum()*(tt*tt).sum())+1e-6)))
            return b_
        for step in (-2,2):
            p=h['p']
            for _ in range(3):
                p+=step;y=mid+0.5-p*d/2
                if any(o['st']==si and o['p']==p and abs(o['x']-h['x'])<0.8*d for o in det+add): continue
                gh,gf=good(h['x'],y,'o_h'),good(h['x'],y,'o_f')
                if gh<0.55 or gh<gf+0.03: break
                add.append(dict(h,y=float(y),p=p,sc=round(gh,2),stack=True))
    return det+add
def segno_filter(det,bw,st):
    """Segno-Zeichen (Rosenzopf: am Teilbeginn und am Schluss) wird als Köpfe gelesen. Erkennung: die Köpfe liegen in einem Tintenfleck
    von mindestens 2,6 Linienabständen Höhe ohne geraden Hals (längster senkrechter Strich unter 1,6 Linienabständen), direkt neben
    einem dicken Taktstrich (Wiederholung, Schluss). GS_SEGNO=0: aus"""
    if _os.environ.get('GS_SEGNO','1')!='1' or not det: return det   # nur neben einem dicken Taktstrich (sonst nicht sicher von dichten Akkorden zu trennen)
    out=[];cache={}
    for h in det:
        si=h['st']
        if si>=len(st): out.append(h);continue
        if si not in cache:
            ys=st[si];d=(ys[4]-ys[0])/4;nl=G.remove_lines(bw,ys,d);t=max(0,int(ys[0]-2*d));b=min(bw.shape[0],int(ys[4]+2*d))
            sub=cv2.morphologyEx(nl[t:b],cv2.MORPH_CLOSE,np.ones((max(3,int(0.4*d)),1),np.uint8))
            n_,lab,S_,_=cv2.connectedComponentsWithStats((sub>0).astype(np.uint8));bad=set()
            for i in range(1,n_):
                x,y,w,hh,a=S_[i]
                if hh<2.6*d or not (0.9*d<=w<=2.8*d): continue
                m=(lab[y:y+hh,x:x+w]==i);run=0;many=0
                for c in range(m.shape[1]):
                    col=np.concatenate([[0],m[:,c].astype(np.int8),[0]]);e=np.flatnonzero(np.diff(col))
                    if len(e): run=max(run,int((e[1::2]-e[::2]).max()));many+=len(e)//2>=4
                spiral=w>=1.4*d and many>=0.5*m.shape[1]   # Rosenzopf: gewundenes Zeichen, viele schräge Striche übereinander
                if run>=1.6*d and not spiral: continue
                # daneben (bis 2,5 Linienabstände) ein dicker Strich über die ganze Systemhöhe: Wiederholung oder Schluss
                def thick(lo,hi):
                    lo,hi=max(0,int(lo)),min(bw.shape[1],int(hi))
                    if hi<=lo: return False
                    cov=((bw[int(ys[0]):int(ys[4])+1,lo:hi]>0).mean(0)>0.9).astype(np.int8)
                    r_=np.flatnonzero(np.diff(np.concatenate([[0],cov,[0]])));return any(b2-a2>=0.3*d for a2,b2 in zip(r_[::2],r_[1::2]))
                if thick(x+w,x+w+2.5*d) or thick(x-2.5*d,x): bad.add(i)
            cache[si]=(t,lab,bad)
        t,lab,bad=cache[si];yy=int(h['y'])-t;xx=int(h['x'])
        if bad and 0<=yy<lab.shape[0] and 0<=xx<lab.shape[1] and lab[yy,xx] in bad: continue
        out.append(h)
    return out
def thick_bar_filter(det,bw):
    """Köpfe auf dicken Taktstrichen (Wiederholungszeichen, Schlussstrich) verwerfen: auf Fotos liest das Kopfmodell
    die dicke senkrechte Linie oft als Kopf oder Kreuz und verdeckt damit auch den Taktstrich (GS_THICK=0: aus)"""
    if _os.environ.get('GS_THICK','1')=='0' or not det or not G.FOTO_LINES: return det   # nur Fotos: auf sauberen Seiten kostet es echte Akkorde
    out=[];cache={}
    for h in det:
        k=h['st']
        if k not in cache:
            d=h['d'];top=int(round(h['mid']-2*d));bot=int(round(h['mid']+2*d))
            if top<0 or bot>=bw.shape[0]: cache[k]=[]
            else:
                ins=(bw[top:bot+1]>0).mean(0)>0.92;e=np.flatnonzero(np.diff(np.concatenate([[0],ins.astype(np.int8),[0]])))
                cache[k]=[(a+b_)/2 for a,b_ in zip(e[::2],e[1::2]) if 0.3*d<=b_-a<=0.9*d and straight_bar(bw,top,bot,a,b_,d)]
        if any(abs(h['x']-x)<0.8*h['d'] for x in cache[k]): continue
        out.append(h)
    return out
ERKENNER=''   # vom Space je Scan gesetzt: 'alt' oder 'yolo'; leer = Variable GS_YOLO entscheidet
def tri_seite(det):
    """Seite in Dreieck-Schreibweise (Margetich/Philipp)? Gemessen am 11.11.2026: Philipps Notenseiten haben mindestens 8 sichere
    Dreiecke (sc>=0.6) und mindestens 5 % Dreiecke unter allen Köpfen; die fälschlich eingestuften Seiten (Schaborak, Boarischer,
    Kiermeier) höchstens 12 sichere Dreiecke bei höchstens 3,8 %. Bisher genügten 3 Dreiecke. GS_TRIMIN/GS_TRIANTEIL zum Anpassen."""
    t=[h for h in det if h['t']=='t' and h.get('sc',1)>=0.6]
    return len(t)>=int(_os.environ.get('GS_TRIMIN_Y','6')) and len(t)>=float(_os.environ.get('GS_TRIANTEIL','0.045'))*max(1,len(det))
def _koepfe(work,g,st):
    """Köpfe der Seite: bisheriger Erkenner (detect.py) oder der gelernte Erkenner (yolo_koepfe.py), wenn das Modell vorhanden ist.
    Wahl: ERKENNER (aus dem Wandler) vor der Space-Variablen GS_YOLO=1. Fotos liest auch bei GS_YOLO=1 der bisherige Erkenner,
    weil der gelernte dort mehr Takte offen lässt (GS_YOLOFOTO=1 nimmt ihn auch für Fotos); ERKENNER='yolo' gilt immer."""
    D.LAST['erkenner']='alt'
    _foto_alt=bool(G.FOTO_LINES) and _os.environ.get('GS_YOLOFOTO','0')!='1'
    if ERKENNER=='yolo' or (ERKENNER!='alt' and _os.environ.get('GS_YOLO','0')=='1' and not _foto_alt):
        try:
            import yolo_koepfe as YK
            if YK.vorhanden():
                det=YK.detect(g,st);D.LAST['tri_page']=tri_seite(det);D.LAST['yolo']=len(det);D.LAST['erkenner']='yolo'
                return det
        except Exception as e: D.LAST['yolo']='Fehler: '+str(e)[:80]
    return D.detect(work,1)
def kreuz_rund():
    """Kreuze für homr rund zeichnen? GS_GKREUZRUND=1 erzwingt es, =0 schaltet es ab; sonst nach Schreibweise (bayerisch, Handschrift)."""
    v=_os.environ.get('GS_GKREUZRUND','')
    return v=='1' if v in ('0','1') else bool(STYLE_RULES[STYLE[0]].get('kreuzrund'))
def tri_rund(g,det,st):
    """Kopie des Graubilds für homr: jedes Dreieck weiß abdecken, die Notenlinien dort wieder ziehen und an derselben Stelle einen
    runden Kopf zeichnen (voll, wenn das Dreieck voll war, sonst hohl). Die übrigen Köpfe bleiben unverändert."""
    g2=g.copy();H,W=g2.shape[:2]
    arten=('t','x') if kreuz_rund() else ('t',)   # Kreuze nur je nach Schreibweise (kreuz_rund)
    for h in det:
        if h.get('t') not in arten: continue
        d=float(h.get('d') or 0)
        if d<=2: continue
        x,y=float(h['x']),float(h['y'])
        x0,x1=int(max(0,x-0.72*d)),int(min(W,x+0.72*d+1));y0,y1=int(max(0,y-0.58*d)),int(min(H,y+0.58*d+1))
        g2[y0:y1,x0:x1]=255
        lw=max(1,int(round(0.1*d)))
        si=h.get('st')
        if si is not None and 0<=si<len(st):
            for ly in st[si]:
                if y0-lw<=ly<=y1+lw: cv2.line(g2,(x0,int(round(ly))),(x1-1,int(round(ly))),0,lw)
        ax=(max(2,int(round(0.6*d))),max(2,int(round(0.42*d))))
        if h.get('fill') or h.get('t')=='x': cv2.ellipse(g2,(int(round(x)),int(round(y))),ax,-20,0,360,0,-1)
        else: cv2.ellipse(g2,(int(round(x)),int(round(y))),ax,-20,0,360,0,max(1,int(round(0.14*d))))
    return g2
def read(path,shift=0,work=None,table=None,bass=None):
    """shift: Stimmung wie im Wandler (A-D-G-C = 0, G-C-F-B = -2 …). bass: Liste, an die je Takt die gelesenen Bassbuchstaben (Lage im Takt in Zwölftel-Vierteln, Buchstabe) angehängt werden."""
    g=load_gray(path)
    if g is None: raise ValueError('Keine Notenlinien gefunden')
    _DBG['g']=g
    if work is None:
        import tempfile;work=_tmp=tempfile.mktemp(suffix='.png')
    cv2.imwrite(work,g);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV)
    st,bars_d=druck_bars(g)
    # Vektor-PDF (MC-Score, Finale): Köpfe aus der Textebene statt aus dem Bild (vektor.py, GS_VEKTOR=0 schaltet ab)
    _vk=None
    try:
        import vektor as VK;_vh=VK.load(path)
        if _vh:
            _W,_H=Image.open(path).size
            if abs(g.shape[0]/g.shape[1]-_H/_W)<0.01: _vk=VK.det_from_heads(_vh,g.shape[1]/_W,st)
    except Exception: _vk=None
    _DBG['vektor']=len(_vk) if _vk else 0
    if _vk:
        det=_vk;D.LAST['tri_page']=False;XDROP.clear();C=chords_det(det)
    else:
        det=hollow_stack(prefix_marks(cross_marks(segno_filter(thick_bar_filter(_koepfe(work,g,st),bw),bw,st)),bw,st),bw,st);C=chords_det(det)
    off=shift+2   # Tabelle in assemble.py gilt für G-C-F-B (shift -2)
    items=[]   # je System: Akkorde und Pausen
    fills=[];dots=[]
    _dc={}
    def dot_img(si,d,nl):   # Fotos: Punkte im Graubild ohne Linien suchen (remove_lines lässt dort Linienreste stehen, die am Punkt kleben)
        if not (G.FOTO_LINES and _os.environ.get('GS_DOTCLEAN')): return nl   # erprobt: im Graubild sind die Punkte zu klein, daher aus
        if si not in _dc: _dc[si]=G.clean_gray(g,d,int(_os.environ.get('GS_DOTT','60')))
        return _dc[si]
    feats=[]
    # Zeilenanfang je System (Taktangabe im ersten System nicht als Kopf lesen)
    x_start={};x_ext={};foto_like=bool(np.median(g)<253)
    for si,ys in enumerate(st):
        row=np.where(bw[int(round(ys[0]))]>0)[0];x_start[si]=float(row.min()) if len(row) else 0.0
        # Ausdehnung des Systems: Spalten, in denen (fast) alle fünf Linien Tinte haben; längstes Stück (Fotos: Seitenrand und Nachbarseite bleiben draußen)
        # Toleranz um jede Linie (GS_XEXTW, in Linienabständen): auf gewölbten Fotoseiten laufen die Linien bis zu einem halben Abstand von ihrer Höhe weg
        _xw=float(_os.environ.get('GS_XEXTW','0.45' if (foto_like or G.FOTO_LINES) else '0.25'));_pad=int(8+_xw*(ys[4]-ys[0])/4)
        hl=cv2.morphologyEx(bw[max(0,int(ys[0])-_pad):int(ys[4])+_pad+1],cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(max(5,int(2*(ys[4]-ys[0])/4)),1)));o0=max(0,int(ys[0])-_pad)
        w_=max(1,int(_xw*(ys[4]-ys[0])/4));cols=sum((hl[max(0,int(round(y))-o0-w_):int(round(y))-o0+w_+1].max(0)>0).astype(np.int8) for y in ys)>=3
        # volle Tintenblöcke breiter als 1,5 Linienabstände (Stempel, Aufkleber) sind kein System; Taktstriche und Wiederholungszeichen bleiben
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
    # Taktangabe im ersten System: statt pauschal 6 Linienabstände nach dem Zeilenanfang nur die Ziffern selbst (zwei übereinander stehende
    # Flecken von etwa 2 Linienabständen Höhe) und davor Schlüssel/TAB-Zeichen; ein Akkord gleich nach der Taktangabe bleibt (Rosenzopf, Philipp)
    _tsr=[None]
    def _ts_skip(x):
        d0=(st[0][4]-st[0][0])/4
        if _os.environ.get('GS_TSOLD')=='1': return x<x_start[0]+6*d0
        if _tsr[0] is None:
            ys0=st[0];x0=int(x_start[0]);x1=int(min(bw.shape[1],x_start[0]+7*d0));nl0=G.remove_lines(bw,ys0,d0)
            sub=cv2.morphologyEx(nl0[int(ys0[0]-0.3*d0):int(ys0[4]+0.3*d0),x0:x1],cv2.MORPH_CLOSE,np.ones((max(3,int(0.4*d0)),1),np.uint8));n_,_,S_,_=cv2.connectedComponentsWithStats((sub>0).astype(np.uint8))
            dig=[(x0+S_[i][0],x0+S_[i][0]+S_[i][2]) for i in range(1,n_) if 1.3*d0<=S_[i][3]<=4.6*d0 and 0.7*d0<=S_[i][2]<=2.7*d0]
            if dig:   # nur die linke Gruppe (Ziffern übereinander); ein Akkord gleich danach sieht ähnlich aus, steht aber rechts davon
                l0=min(dig);dig=[q for q in dig if q[0]<l0[1]]
            _tsr[0]=(min(a for a,_ in dig)-0.3*d0,max(b for _,b in dig)+0.3*d0) if dig else (x_start[0],x_start[0]+6*d0)
        a_,b_=_tsr[0];return x<x_start[0]+1.5*d0 or a_<=x<=b_
    for c in C:
        si=c[0]['st'];x=float(np.median([h['x'] for h in c]))
        xo=[h['x'] for h in c if h['t']!='x']   # Kreuze links vor dem Akkord: Lage des Akkords (und Hals) nach den runden Köpfen
        if xo and len(xo)<len(c) and all(h['x']<min(xo)-0.5*h['d'] for h in c if h['t']=='x'): x=float(np.median(xo))
        if si==0 and _ts_skip(x): continue
        if si in x_ext and not (x_ext[si][0]-0.5*(st[si][4]-st[si][0])/4<=x<=x_ext[si][1]+0.5*(st[si][4]-st[si][0])/4): continue   # außerhalb des Systems
        dr='D' if any(a-3<=x<=b+3 for a,b,_ in bars_d.get(si,[])) else 'Z'
        ks=[button(h['p'],h['t']) for h in c];ps=[pitch(r,n,dr,table) for r,n in ks]
        hs=[(h,k,p+off) for h,k,p in zip(c,ks,ps) if p]
        if not hs: continue
        ys=st[si];d=(ys[4]-ys[0])/4;nl=G.remove_lines(bw,ys,d)
        fill=float(np.mean([h['fill'] for h,_,_ in hs]));fills.append(fill)
        dots.append(dot_feat(dot_img(si,d,nl),x,[h['y'] for h,_,_ in hs],d));feats.append(np.concatenate([strip(nl,x,ys[2],d,fill),np.array(beam_feats(nl,x,[h['y'] for h,_,_ in hs],d)+[dots[-1]])*20]).astype(np.float32))
        items.append(dict(kind='n',st=si,x=x,dir=dr,heads=[(k[0],k[1],p) for _,k,p in hs],fi=len(feats)-1,hy=[(h['x'],h['y'],k[0]) for h,k,_ in hs],ps=[h['p'] for h,_,_ in hs]))
    _DBG['feats']=feats;_DBG['x_ext']=x_ext;_DBG['x_start']=x_start
    P=RCLF.predict_proba(np.array(feats)) if feats else np.zeros((0,len(RCLF.classes_)))
    P9=np.zeros((len(P),len(UNIT)));P9[:,RCLF.classes_]=P;lp=np.log(np.clip(P9,1e-6,1))
    if RCLF2 is not None and len(feats):
        Q9=np.zeros((len(P),len(UNIT)));Q9[:,RCLF2.classes_]=RCLF2.predict_proba(np.array(feats))
        lp=(1-RW2)*np.log(np.clip(P9,1e-4,1))+RW2*np.log(np.clip(Q9,1e-4,1))
    # hohle Köpfe sind Halbe oder länger, volle Köpfe kürzer (weiche Vorgabe, das Modell kann sie mit klaren Hinweisen überstimmen)
    FP=float(_os.environ.get('GS_FILLPEN','0'))
    DP=float(_os.environ.get('GS_DOTPEN','0'))   # erprobt, ohne Gewinn: aus   # erprobt: auf Rosenzopf-Fotos schlechter, daher aus
    for i_,(it,f_) in enumerate(zip(items,fills)):
        if G.FOTO_LINES and dots[i_]>0.5 and DP>0:   # Fotos: ein klar erkannter Punkt spricht für einen punktierten Wert
            for k_,u_ in enumerate(UNIT):
                if u_ not in (9,18,36): lp[i_][k_]-=DP
        if FP>0:
            for k_,u_ in enumerate(UNIT):
                if (f_<0.5 and u_<24) or (f_>=0.5 and u_>=24 and not _os.environ.get('FP_HOLONLY')): lp[i_][k_]-=FP
    for it,row in zip(items,lp): it['lp']=row
    meas=[];medge=[];letters={}
    for si,ys in enumerate(st):
        d=(ys[4]-ys[0])/4;its=sorted([it for it in items if it['st']==si],key=lambda it:it['x'])
        hxy=[]
        for it in its:
            for r_,n_,p_ in it['heads']: hxy.append((it['x'],0))
        heads_xy=[(h['x'],h['y']) for c in C if c[0]['st']==si for h in c]+[(h['x'],h['y']) for h in XDROP if h['st']==si]
        bl=barlines(bw,ys,d,[x for x,_ in heads_xy]);nol=G.stems_off(G.remove_lines(bw,ys,d),d)
        rs=[dict(kind='r',st=si,x=x,u=u) for x,u in rests(bw,nol,ys,d,heads_xy,bl)]
        edges=[0]+bl+[1e9];allg=[]
        if bass is not None and not _os.environ.get('NOBASS'):
            _vb=None
            if _vk:
                try:
                    import vektor as VK;_bb=VK.load_bass(path)
                    if _bb is not None: _vb=VK.letters_from_bass(_bb,g.shape[1]/Image.open(path).size[0],st,si)
                except Exception: _vb=None
            letters[si]=_vb if _vb else BR.read_letters(bw,st,si,d)   # Vektor-PDF: Bassbuchstaben aus der Textebene
        if _os.environ.get('GS_BARDOTS','1')=='1':   # „Köpfe“ direkt am Taktstrich sind Wiederholungspunkte (Fuchs: dicker Strich mit Punkten)
            its=[it for it in its if not any(abs(it['x']-x_)<0.8*d for x_ in bl)]
        if _os.environ.get('GS_LEDGER','1')=='1':   # Köpfe weit über oder unter dem System brauchen Hilfslinien; ohne sind es Zahlen oder Text (etwa „MM = 92“)
            def _hl(it,up):
                yy=int(round(ys[0]-d if up else ys[4]+d));x0,x1=int(it['x']-0.9*d),int(it['x']+0.9*d)
                return max(((bw[max(0,yy+k),max(0,x0):x1]>0).mean() for k in range(-2,3)),default=0)>0.6
            its2=[]
            for it in its:
                ps_=it.get('ps') or []
                if ps_ and len(ps_)==len(it.get('heads',[])):
                    up_=_hl(it,True) if max(ps_)>=8 else True;dn_=_hl(it,False) if min(ps_)<=-8 else True
                    hy_=it.get('hy') if isinstance(it.get('hy'),list) and len(it.get('hy'))==len(ps_) else None
                    nb_=lambda k:any(abs((hy_[k][0] if hy_ else it['x'])-b_)<1.5*d for b_ in bl)   # Zierstriche über/unter einem Wiederholungsstrich (Rosenzopf) sind keine Hilfslinien
                    keep_=[k for k,q in enumerate(ps_) if not ((q>=8 and not up_) or (q<=-8 and not dn_) or (G.FOTO_LINES and abs(q)>=6 and len(ps_)==1 and nb_(k)))]   # einzelne Köpfe ohne Hilfslinie weglassen (Text über dem Akkord)
                    if not keep_: continue
                    if len(keep_)<len(ps_):
                        it=dict(it)
                        for f_ in ('heads','ps','hy'):
                            if isinstance(it.get(f_),list) and len(it[f_])==len(ps_): it[f_]=[it[f_][k] for k in keep_]
                its2.append(it)
            its=its2
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
    _DBG['gerust']=None
    # Gerüst aus der Notenerkennung (GS_GERUEST=1): Takte, Notenwerte und Pausen von homr, Köpfe und Balgrichtung vom Leser (gerust.py)
    # Dreieck-Seiten (Philipp/Margetich): homr liest Dreiecke schlecht. Für homr wird daher eine Kopie gelesen, in der jedes Dreieck
    # durch einen runden Kopf ersetzt ist (voll/hohl wie das Dreieck; Hals, Fähnchen, Balken bleiben). Kopfform und Reihe kommen
    # weiter vom Original. GS_GERUESTTRI=0 schaltet homr auf Dreieck-Seiten ganz ab (Stand vor v159).
    _tri=bool(D.LAST.get('tri_page'))
    if GERUEST and (not G.FOTO_LINES or _os.environ.get('GS_GERUESTFOTO','1')=='1') and not (_tri and _os.environ.get('GS_GERUESTTRI','1')=='0'):   # Seiten mit Dreiecken (Philipp/Margetich): homr liest sie schlecht, daher Leser   # Fotos: mit homrs eigener Systemsuche (gerust.py, GS_GFOTOHOMR); GS_GERUESTFOTO=0 nimmt Fotos aus
        import gerust as GR
        GR.FOTO=bool(G.FOTO_LINES)
        _gh=tri_rund(g,det,st) if (_tri or kreuz_rund()) else g
        _DBG['gerust']={};meas,medge=GR.rebuild(_gh,st,x_ext,meas,medge,button,pitch,table,off,_DBG['gerust'])
        if _tri: _DBG['gerust']['tri_rund']=True
    def gl(gr): return sum(e['lp'] for e in gr) if 'lp' in gr[0] else np.zeros(len(UNIT))
    def uu(gr,u): return TRIP.get(u,u) if gr[0].get('trip') else u
    def letter_times(mi,si,m,Lm):
        """Fotos: Die Bassbuchstaben stehen auf den Schlägen und bilden ein Zeitraster. Einsätze der Noten und Pausen über ihre Lage
        zwischen den Buchstaben bestimmen (auf Achtel gerundet), Dauern als Abstand zum nächsten Einsatz. None, wenn es nicht aufgeht."""
        if any(gr[0].get('trip') for gr in m) or not letters.get(si): return None
        a,b,dd=medge[mi];L_=sorted(x for x,c in letters[si] if a<x<b);k=len(L_)
        if k<2 or Lm%k: return None
        step=Lm//k;sx=(L_[-1]-L_[0])/(k-1)
        if sx<1.5*dd: return None
        A=[(L_[0]-sx,-step)]+[(x,i*step) for i,x in enumerate(L_)]+[(L_[-1]+sx,Lm)]
        def tx(x):
            for (x1,t1),(x2,t2) in zip(A,A[1:]):
                if x<=x2: return t1+(x-x1)/max(1e-6,x2-x1)*(t2-t1)
            return Lm
        grid=6 if step%6==0 else 3
        ts=[max(0,int(round(tx(float(np.mean([e['x'] for e in gr])))/grid)*grid)) for gr in m]
        if any(t2<=t1 for t1,t2 in zip(ts,ts[1:])) or ts[0]!=0 or ts[-1]>=Lm: return None
        us=[t2-t1 for t1,t2 in zip(ts,ts[1:])]+[Lm-ts[-1]];out=[]
        for gr,u in zip(m,us):
            if gr[0]['kind']=='r':
                if u!=gr[0]['u']: return None
            else:
                if u not in UNIT: return None
                out.append(UNIT.index(u))
        return out
    def length(m,ch=None):
        t=0;k=0
        for q,gr in enumerate(m):
            if gr[0].get('hu') is not None: t+=gr[0]['hu'];k+=gr[0]['kind']!='r'
            elif gr[0]['kind']=='r': t+=uu(gr,gr[0]['u'])
            else: t+=uu(gr,UNIT[ch[k] if ch else int(np.argmax(gl(gr)))]);k+=1
        return t
    cg=[c for c in collections.Counter(length(m) for _,m in meas).most_common() if c[0] in (18,24,36,48,72)];Lg=cg[0][0] if cg else 36
    # Taktlänge je System (Taktwechsel zwischen Teilen): häufigste Länge im System, wenn mindestens zwei Takte sie haben
    Ls={}
    for si in set(si for si,_ in meas):
        cs=[c for c in collections.Counter(length(m) for s_,m in meas if s_==si).most_common() if c[0] in (18,24,36,48,72)]
        nm_=sum(1 for s_,_ in meas if s_==si);Ls[si]=cs[0][0] if cs and cs[0][1]>=2 and (not G.FOTO_LINES or cs[0][0]==Lg or cs[0][1]>=max(3,0.5*nm_)) else Lg   # Fotos: Abweichung vom Stück nur bei klarer Mehrheit im System   # nur übliche Taktarten (3/8, 2/4, 3/4, 4/4, 6/4)
    # Textzeichen über den Systemen (Einleitung, Trio, Fine, D.C. …) mit Lage: dem Takt zuordnen, über dem sie stehen (GS_WORTE=0: aus)
    _WORDS={};_VOLTA={};_VOLTA_END={}
    if _os.environ.get('GS_WORTE','1')=='1':
        try:
            import pytesseract
            D_=pytesseract.image_to_data(g,lang='deu',config='--psm 11',output_type=pytesseract.Output.DICT)
            want={w.lower().replace(' ','').rstrip('.'):w for w in TEIL_WOERTER}
            used_=set()
            for i_,txt in enumerate(D_['text']):
                t_=(txt or '').strip().replace(' ','')
                if not t_ or float(D_['conf'][i_])<55 or i_ in used_: continue
                key_=t_.lower().replace(',','.').rstrip(':,;').rstrip('.')
                if key_ in ('d','d.') and i_+1<len(D_['text']) and D_['line_num'][i_+1]==D_['line_num'][i_] and (D_['text'][i_+1] or '').strip().lower().rstrip('.') in ('s','c'):
                    used_.add(i_+1);key_='d.'+D_['text'][i_+1].strip().lower().rstrip('.');t_='D.'+key_[-1].upper()+'.'   # „D. S.“ in zwei Teilen gelesen
                if key_ in ('d.s','d.c','ds','dc'):   # Sprunganweisung mit Rest der Zeile: „D.S. al Fine dann Trio“
                    ph=[]   # ganze Zeile: „Trio D.S. al Fine“, „D.S. al Fine dann Trio“
                    for j_ in range(len(D_['text'])):
                        if D_['block_num'][j_]!=D_['block_num'][i_] or D_['line_num'][j_]!=D_['line_num'][i_]: continue
                        w_=(D_['text'][j_] or '').strip()
                        if j_==i_: ph.append('D.'+key_[-1].upper()+'.')
                        elif j_ in used_ and j_==i_+1: continue
                        elif w_ and float(D_['conf'][j_])>=40: ph.append(w_);used_.add(j_)
                    want[key_]=re.sub(r'\bTr[lI1]o\b','Trio',' '.join(ph))
                vm=re.match(r'^\[?([12])\.\]?$',t_)
                if vm and float(D_['conf'][i_])>=30:   # Volta-Klammer „1.“ / „2.“ knapp über dem System
                    cx=D_['left'][i_]+D_['width'][i_]/2;cy=D_['top'][i_]+D_['height'][i_]/2
                    cand=[(ys[0]-cy,si) for si,ys in enumerate(st) if 0<ys[0]-cy<((ys[4]-ys[0])/4)*4]
                    if cand:
                        si_=min(cand)[1];ms_i=[mi_ for mi_,(s_,_) in enumerate(meas) if s_==si_ and medge[mi_][0]-((st[si_][4]-st[si_][0])/4)<=cx<medge[mi_][1]]
                        if ms_i: _VOLTA[ms_i[0]]=int(vm.group(1))
                    continue
                if key_ not in want: continue
                cx=D_['left'][i_]+D_['width'][i_]/2;cy=D_['top'][i_]+D_['height'][i_]/2
                cand=[(ys[0]-cy,si) for si,ys in enumerate(st) if -((ys[4]-ys[0])/4)*1.0<ys[0]-cy<((ys[4]-ys[0])/4)*6]
                if not cand and key_ in ('fine','d.s','d.c','ds','dc'):   # Fine und Sprunganweisungen stehen oft unter dem System (unter den Bassbuchstaben)
                    cand=[(cy-ys[4],si) for si,ys in enumerate(st) if 0<cy-ys[4]<((ys[4]-ys[0])/4)*9]
                if not cand: continue
                si_=min(cand)[1]
                ms_i=[mi_ for mi_,(s_,_) in enumerate(meas) if s_==si_]
                if not ms_i: continue
                best=min(ms_i,key=lambda mi_:0 if medge[mi_][0]<=cx<medge[mi_][1] else min(abs(cx-medge[mi_][0]),abs(cx-medge[mi_][1])))
                if want[key_] not in _WORDS.setdefault(best,[]): _WORDS[best].append(want[key_])
        except Exception: _WORDS={}
    # Volta-Klammern geometrisch: waagrechter Strich 1–5 Linienabstände über dem System, links mit Haken nach unten; Ziffer darin einzeln lesen
    if _os.environ.get('GS_VOLTA','1')=='1':
        try:
            import pytesseract
            for si,ys in enumerate(st):
                d=(ys[4]-ys[0])/4;y0=max(0,int(ys[0]-4.5*d));y1=int(ys[0]-1.0*d)
                if y1-y0<d: continue
                band=bw[y0:y1];hz=cv2.morphologyEx(band,cv2.MORPH_OPEN,cv2.getStructuringElement(cv2.MORPH_RECT,(max(5,int(3*d)),1)))
                n_,lab_,S_,_=cv2.connectedComponentsWithStats((hz>0).astype(np.uint8))
                vols=[]
                for i in range(1,n_):
                    x,y,w,h,_a=S_[i]
                    if w<6*d or h>(0.45 if w>=10*d else 0.35)*d: continue   # Klammer: lang und dünn (Balken sind kürzer, Druckstriche dicker; lange schräge Klammern auf Fotos etwas dicker)
                    yy=y0+y;xl=x
                    hook=bw[yy:int(yy+1.2*d),max(0,xl-2):xl+int(0.3*d)]
                    if hook.size==0 or (hook.max(1)>0).mean()<(0.3 if w>=8*d else 0.7): continue   # kein Haken links: Balken o. Ä. (lange Klammern mit kurzem Haken, Michlbauer)
                    crop=g[int(yy+0.15*d):int(yy+1.7*d),int(xl+0.15*d):int(xl+2.2*d)]
                    num=None
                    if crop.size:
                        crop=cv2.resize(crop,None,fx=3,fy=3,interpolation=cv2.INTER_CUBIC)
                        t=pytesseract.image_to_string(crop,config='--psm 7 -c tessedit_char_whitelist=12.').strip()
                        mm_=re.search(r'[12]',t);num=int(mm_.group(0)) if mm_ else None
                    vols.append((xl,xl+w,num))
                vols=sorted(vols);vols=[(xa,xb,num if num else (2 if k_ and xa-vols[k_-1][1]<1.5*d else 1)) for k_,(xa,xb,num) in enumerate(vols)]   # Nummer: direkt anschließende Klammer ist die 2.
                for k_,(xa,xb,num) in enumerate(vols):
                    ms_i=[mi_ for mi_,(s_,_) in enumerate(meas) if s_==si]
                    if not ms_i: continue
                    ov=lambda mi_:min(xb,medge[mi_][1])-max(xa,medge[mi_][0])   # Überdeckung Klammer/Takt
                    start=min(ms_i,key=lambda mi_:0 if medge[mi_][0]<=xa+0.5*d<medge[mi_][1] else min(abs(xa-medge[mi_][0]),abs(xa-medge[mi_][1])))
                    if ov(start)<0.5*min(xb-xa,medge[start][1]-medge[start][0]):   # Anfang knapp vor dem Taktstrich: Takt mit der größten Überdeckung
                        start=max(ms_i,key=ov)
                    end=min(ms_i,key=lambda mi_:0 if medge[mi_][0]<xb-0.5*d<=medge[mi_][1] else min(abs(xb-medge[mi_][0]),abs(xb-medge[mi_][1])))
                    _VOLTA[start]=num;_VOLTA_END[start]=max(start,end)
        except Exception: pass
        ks_=sorted(_VOLTA)   # eine Klammer endet spätestens vor der nächsten
        for a_,b_ in zip(ks_,ks_[1:]):
            if _VOLTA_END.get(a_,a_)>=b_: _VOLTA_END[a_]=b_-1
    # Wiederholungszeichen über die Punkte: zwei Punkte in den beiden mittleren Zwischenräumen knapp rechts vom Taktanfang
    # (Wiederholung beginnt) oder knapp links vom Taktende (Wiederholung endet). homr übersieht sie vor allem auf Fotos (GS_REPDOTS=0: aus)
    _DBG['repdots']=0
    if _os.environ.get('GS_WIEDERHOLUNG','1')=='1' and _os.environ.get('GS_REPDOTS','1')=='1':
        try:
            _nl={}
            def _dots(si,ys,d,xa,xb,side):
                if si not in _nl: _nl[si]=G.remove_lines(bw,ys,d)
                y1,y2=(ys[1]+ys[2])/2,(ys[2]+ys[3])/2;xa,xb=max(0,int(xa)),min(bw.shape[1],int(xb))
                if xb-xa<3: return False
                win=_nl[si][int(ys[1]):int(ys[3])+1,xa:xb]
                n_,_,S_,C_=cv2.connectedComponentsWithStats((win>0).astype(np.uint8));hit=[]
                for i in range(1,n_):
                    w_,h_,a_=S_[i][2],S_[i][3],S_[i][4]
                    if 0.2*d<=w_<=0.75*d and 0.2*d<=h_<=0.75*d and a_>=0.5*w_*h_: hit.append((C_[i][0]+xa,C_[i][1]+int(ys[1])))
                up=[(x,y) for x,y in hit if abs(y-y1)<0.45*d];dn=[(x,y) for x,y in hit if abs(y-y2)<0.45*d]
                px=[(u+v)/2 for u,uy in up for v,vy in dn if abs(u-v)<0.4*d and 0.7*d<vy-uy<1.3*d]   # zwei Punkte genau übereinander, eine Linie Abstand
                if not px: return False
                # daneben (zur Taktgrenze hin) ein doppelter Strich über die ganze Systemhöhe: dünn und dick
                for x0 in px:
                    lo,hi=(int(x0-1.4*d),int(x0-0.2*d)) if side=='l' else (int(x0+0.2*d),int(x0+1.4*d))
                    lo,hi=max(0,lo),min(bw.shape[1],hi)
                    cov=(bw[int(ys[0]):int(ys[4])+1,lo:hi]>0).mean(0)>0.9
                    runs=np.flatnonzero(np.diff(np.concatenate([[0],cov.astype(np.int8),[0]])))
                    wid=[b2-a2 for a2,b2 in zip(runs[::2],runs[1::2])]
                    if len(wid)>=2 or any(w2>=0.3*d for w2 in wid): return True
                return False
            for mi,(si,m) in enumerate(meas):
                if not m: continue
                a_,b_,d=medge[mi];ys=st[si]
                xs_=[gr[0]['x'] for gr in m]
                if b_>bw.shape[1]:   # letzter Takt im System: Ende der Notenlinien
                    row=(bw[int(ys[0])-1:int(ys[0])+2].max(0)>0);nz=np.flatnonzero(row);b_=float(nz[-1]) if nz.size else max(xs_)+3*d
                fst=a_<=0 or min(xs_)-a_>3.4*d   # erster Takt im System (Anfang ist der Rand der Notenlinien)
                if fst: a_=min(xs_)-3.4*d   # erster Takt im System: Zeichen kurz vor der ersten Note
                mf=m[0][0].setdefault('mflags',{})
                if not mf.get('rep_start') and _dots(si,ys,d,a_-0.5*d,a_+1.8*d if not fst else min(xs_)-0.3*d,'l'): mf['rep_start']=True;_DBG['repdots']+=1
                if not mf.get('rep_end') and _dots(si,ys,d,b_-1.8*d,b_+0.5*d,'r'): mf['rep_end']=True;_DBG['repdots']+=1
                if not mf: del m[0][0]['mflags']
        except Exception as e_: _DBG['repdots']=str(e_)[:80]
        # Klammer direkt nach einem Wiederholungsende ist die 2.; fehlt davor die 1., gilt der Takt mit dem Wiederholungsende als 1. Klammer
        for k_ in sorted(_VOLTA):
            if _VOLTA[k_]==1 and k_>0 and k_-1 not in _VOLTA and meas[k_-1][0]==meas[k_][0] and meas[k_-1][1] and (meas[k_-1][1][0][0].get('mflags') or {}).get('rep_end'):
                _VOLTA[k_]=2;_VOLTA[k_-1]=1;_VOLTA_END[k_-1]=k_-1
    _DBG['worte']=dict(_WORDS);_DBG['volta']=dict(_VOLTA)
    out=[];lens=[];_VSTOP=set()
    for mi,(si,m) in enumerate(meas):
        Lm=Ls[si]
        notes=[q for q,gr in enumerate(m) if gr[0]['kind']=='n'];ch=[int(np.argmax(gl(m[q]))) for q in notes]
        hm=all(gr[0].get('hu') is not None for gr in m)   # Takt aus dem homr-Gerüst: Notenwerte stehen fest
        if hm and sum(gr[0]['hu'] for gr in m)!=Lm and not (mi==0 and sum(gr[0]['hu'] for gr in m)<Lm) and _os.environ.get('GS_GFIT','1')=='1':
            # geht ein homr-Takt nicht auf, wählt der Leser die Notenwerte selbst (Taktstriche bleiben von homr), homrs Wert als starker Vorschlag
            notes_=[gr for gr in m if gr[0]['kind']!='r']
            if notes_ and all('lp' in gr[0] or gr[0].get('nurhomr') for gr in notes_):
                PG=[]
                for gr in notes_:
                    v=np.array(gl(gr)) if 'lp' in gr[0] else np.full(len(UNIT),-8.0)
                    hu_=gr[0]['hu'];base=TRIP.get(hu_,hu_) if False else hu_
                    if gr[0].get('htrip'):
                        inv={vv:kk for kk,vv in TRIP.items()};base=inv.get(hu_,hu_)
                    if base in UNIT: v=v.copy();v[UNIT.index(base)]+=float(_os.environ.get('GS_GPRIOR','2.0'))
                    PG.append(np.exp(v-v.max()))
                PG=np.array(PG);PG=PG/PG.sum(1,keepdims=True)
                rest_u=[gr[0]['hu'] for gr in m if gr[0]['kind']=='r']
                f=fit_measure(maxloss=float(_os.environ.get('GS_FITLOSS','12' if G.FOTO_LINES else '6')),P=PG,rest_units=rest_u,L=Lm,units=[[(TRIP.get(u,u) if gr[0].get('htrip') else u) for u in UNIT] for gr in notes_])
                if f:
                    k_=0
                    for gr in m:
                        if gr[0]['kind']=='r': continue
                        u_=UNIT[f[k_]];gr[0]['hu']=TRIP.get(u_,u_) if gr[0].get('htrip') else u_;k_+=1
                    _DBG.setdefault('gerust',{})['nachgewaehlt']=_DBG.get('gerust',{}).get('nachgewaehlt',0)+1
        lch=letter_times(mi,si,m,Lm) if (not hm and (G.FOTO_LINES or _os.environ.get('GS_LETALL','1')=='1') and not _os.environ.get('NOLETTIME')) else None
        if hm: pass
        elif lch: ch=lch
        elif length(m)!=Lm and not (mi==0 and length(m)<Lm):
            PG=np.array([np.exp(gl(m[q])-gl(m[q]).max()) for q in notes]);PG=PG/PG.sum(1,keepdims=True) if len(PG) else PG
            f=fit_measure(maxloss=float(_os.environ.get('GS_FITLOSS','12' if G.FOTO_LINES else '6')),P=PG,rest_units=[uu(gr,gr[0]['u']) for gr in m if gr[0]['kind']=='r'],L=Lm,units=[[uu(m[q],u) for u in UNIT] for q in notes])
            if f: ch=f
        evs=[];k=0
        for gr in m:
            if gr[0].get('hu') is not None:
                if gr[0]['kind']=='r': evs.append(dict(rest=True,dbg=(gr[0]['st'],[],gr[0]['x']),u=gr[0]['hu'],trip=bool(gr[0].get('htrip'))))
                else:
                    heads=sorted({h for e in gr for h in e['heads']},key=lambda h:-h[2]);evs.append(dict(rest=False,dbg=(gr[0]['st'],[q for e in gr for q in e['hy']],gr[0]['x']),fis=[e['fi'] for e in gr],mtrip=bool(gr[0].get('htrip')),u=gr[0]['hu'],trip=bool(gr[0].get('htrip')),dir=gr[0]['dir'],heads=heads));k+=1
            elif gr[0]['kind']=='r': evs.append(dict(rest=True,dbg=(gr[0]['st'],[],gr[0]['x']),u=uu(gr,gr[0]['u']),trip=bool(gr[0].get('trip'))))
            else:
                heads=sorted({h for e in gr for h in e['heads']},key=lambda h:-h[2]);evs.append(dict(rest=False,dbg=(gr[0]['st'],[q for e in gr for q in e['hy']],gr[0]['x']),fis=[e['fi'] for e in gr],mtrip=bool(gr[0].get('trip')),u=uu(gr,UNIT[ch[k]]),trip=bool(gr[0].get('trip')) or UNIT[ch[k]]==4,dir=gr[0]['dir'],heads=heads));k+=1
        # Haltebögen und Taktzeichen (Wiederholung, Volta) aus dem homr-Gerüst an die Ereignisse hängen
        _ni=[gr for gr in m]
        for gr,e in zip(_ni,evs):
            if not e['rest']:
                if gr[0].get('tie_start'): e['tie_start']=True
                if gr[0].get('tie_stop'): e['tie_stop']=True
        if evs and m and m[0][0].get('mflags'): evs[0]['mflags']=dict(m[0][0]['mflags'])
        if evs and mi<len(medge) and medge[mi][1]<g.shape[1]: evs[0]['mx']=[float(medge[mi][0]),float(medge[mi][1])]   # Taktgrenzen im Bild (für Kontrollbild und Editor)
        if evs and mi in _VOLTA and not (evs[0].get('mflags') or {}).get('volta'):
            mf_=evs[0].setdefault('mflags',{});mf_['volta']=_VOLTA[mi]
            if _VOLTA_END.get(mi,mi)==mi: mf_['volta_stop']=True
            else: _VSTOP.add(_VOLTA_END[mi])
        if evs and mi in _VSTOP: evs[0].setdefault('mflags',{})['volta_stop']=True
        if evs and _WORDS.get(mi): evs[0]['words']=_WORDS[mi]
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
PB_COL={1:(220,90,0),2:(0,150,0),3:(0,0,220),4:(190,0,190)}   # Reihe 1 blau, 2 grün, 3 rot, 4 lila (BGR)
PB_UN={2:'16T',3:'16',4:'8T',6:'8',8:'4T',9:'8.',12:'4',16:'2T',18:'4.',24:'2',36:'2.',48:'1'}
def pruefbild_alt(g,meas,lens,num0,maxw=1600):
    """Kontrollbild einer Seite: erkannte Köpfe (Farbe = Reihe), Notenwert über dem System, Taktnummer darunter;
    Takte, deren Länge nicht aufgeht, rot hinterlegt. Ohne Wörter im Bild (Erklärung steht im Wandler). Rückgabe: JPEG-Bytes."""
    st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]);c=cv2.cvtColor(g,cv2.COLOR_GRAY2BGR);ov=c.copy()
    fs=lambda d:max(0.4,d/26);lab=[]
    for mi,(m,L) in enumerate(zip(meas,lens)):
        if not m: continue
        ok=sum(e['u'] for e in m)==L;by={}
        for e in m: by.setdefault(e['dbg'][0],[]).append(e)
        for si,es in by.items():
            if si>=len(st): continue
            ys=st[si];d=(ys[4]-ys[0])/4;xs=[e['dbg'][2] for e in es]
            if not ok: cv2.rectangle(ov,(int(min(xs)-0.9*d),int(ys[0]-1.2*d)),(int(max(xs)+0.9*d),int(ys[4]+1.2*d)),(150,150,255),-1)
            for e in es:
                for hx,hy,r in e['dbg'][1]: lab.append(('c',(int(hx),int(hy)),int(d*0.6),PB_COL.get(r,(0,0,0)),max(2,int(d/9))))
                t=('P' if e['rest'] else '')+PB_UN.get(e['u'],str(e['u']))
                lab.append(('t',t,(int(e['dbg'][2]-0.4*d),int(ys[0]-1.5*d)),fs(d),(0,0,200) if e['rest'] else (150,60,0)))
        si,_,x=m[-1]['dbg']
        if si<len(st):
            ys=st[si];d=(ys[4]-ys[0])/4
            lab.append(('t',str(num0+mi)+('' if ok else '!'),(int(x+0.6*d),int(ys[4]+2.3*d)),fs(d)*1.2,(0,130,0) if ok else (0,0,230)))
    c=cv2.addWeighted(ov,0.35,c,0.65,0)
    for k,*a in lab:
        if k=='c': cv2.circle(c,a[0],a[1],a[2],a[3])
        else: cv2.putText(c,a[0],a[1],cv2.FONT_HERSHEY_SIMPLEX,a[2],a[3],max(1,int(a[2]*2.2)))
    if c.shape[1]>maxw: c=cv2.resize(c,(maxw,int(c.shape[0]*maxw/c.shape[1])),interpolation=cv2.INTER_AREA)
    return cv2.imencode('.jpg',c,[cv2.IMWRITE_JPEG_QUALITY,72])[1].tobytes()
PB_DIR={'D':(200,90,0),'Z':(0,110,230)}   # Druck blau, Zug orange (BGR)
def _xq(tx,q,L,xe):
    """Lage eines Zeitpunkts q (Zwölftel-Viertel im Takt) im Bild: zwischen den Einsätzen, nach dem letzten bis zum Taktstrich verteilt."""
    if not tx: return xe
    for (p0,x0),(p1,x1) in zip(tx,tx[1:]):
        if p0<=q<=p1: return x0+(x1-x0)*(q-p0)/max(1,p1-p0)
    p0,x0=tx[-1]
    if q<=p0: return tx[0][1]
    return x0+(xe-x0)*(q-p0)/max(1,L-p0)
def pruefbild(g,meas,lens,num0,maxw=1600,bass=None):
    """Kontrollbild als Überlagerung: das entzerrte Original blass in Grau, darüber in Farbe, was im erzeugten Ergebnis steht
    (Köpfe an der Lage aus Reihe und Knopf, rund oder Kreuz, hohl ab der Halben; Hals, Fähnchen, Punkt; Pausen; Taktstriche grün;
    Bassbuchstaben). Farbe nach Balgrichtung: Druck blau, Zug orange. Takte, die nicht aufgehen, rot hinterlegt. GS_PB=alt: frühere Ansicht."""
    st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1])
    base=cv2.cvtColor(g,cv2.COLOR_GRAY2BGR);c=(255-(255-base.astype(np.float32))*0.35).astype(np.uint8)   # Original blass
    OFF_={1:-1,2:0,3:1,4:2};fs=lambda d:max(0.4,d/26)
    for mi,(m,L) in enumerate(zip(meas,lens)):
        if not m: continue
        ok=sum(e['u'] for e in m)==L;by={}
        for e in m: by.setdefault(e['dbg'][0],[]).append(e)
        for si,es in by.items():
            if si>=len(st): continue
            ys=st[si];d=(ys[4]-ys[0])/4;mid=ys[2];xs=[e['dbg'][2] for e in es]
            if not ok:
                ov=c.copy();cv2.rectangle(ov,(int(min(xs)-0.9*d),int(ys[0]-1.2*d)),(int(max(xs)+0.9*d),int(ys[4]+1.2*d)),(170,170,255),-1);c=cv2.addWeighted(ov,0.4,c,0.6,0)
            lw=max(2,int(d/8))
            pos=0;tx=[]
            for e in es:
                x=int(e['dbg'][2]);tx.append((pos,x));pos+=e['u']
                if e['rest']:
                    cv2.putText(c,'P'+PB_UN.get(e['u'],str(e['u'])),(int(x-0.5*d),int(mid+0.4*d)),cv2.FONT_HERSHEY_SIMPLEX,fs(d),(0,0,200),max(1,int(fs(d)*2.2)));continue
                col=PB_DIR.get(e.get('dir'),(80,80,80));u=e['u'];base_u=u if not e.get('trip') else {2:3,4:6,8:12,16:24,24:36,32:48}.get(u,u)
                hollow=base_u>=24;ysum=[]
                for r,n,pch in e['heads']:
                    p=2*(n-6)+OFF_.get(r,0);y=int(round(mid+0.5-p*d/2));ysum.append(y)
                    for lp in range(6,abs(p)//2*2+1,2):   # Hilfslinien
                        yl=int(round(mid+0.5-(lp if p>0 else -lp)*d/2));cv2.line(c,(int(x-0.9*d),yl),(int(x+0.9*d),yl),col,1)
                    if r==4 and D.LAST.get('tri_page'):   # Seiten mit Dreiecken (Margetich, Philipp): 4. Reihe als Dreieck, hohl ab der Halben
                        tp_=np.array([[x,y-int(0.58*d)],[x-int(0.66*d),y+int(0.5*d)],[x+int(0.66*d),y+int(0.5*d)]],np.int32)
                        if hollow: cv2.polylines(c,[tp_],True,col,lw)
                        else: cv2.fillPoly(c,[tp_],col)
                    elif r>=3:
                        a_=int(0.45*d);cv2.line(c,(x-a_,y-a_),(x+a_,y+a_),col,lw);cv2.line(c,(x-a_,y+a_),(x+a_,y-a_),col,lw)
                        if hollow: cv2.circle(c,(x,y),int(0.6*d),col,max(1,lw//2))
                    else:
                        cv2.ellipse(c,(x,y),(int(0.62*d),int(0.45*d)),-20,0,360,col,-1 if not hollow else lw)
                if base_u<48 and ysum:   # Hals nach oben, Fähnchen nach Notenwert
                    xs_=x+int(0.6*d);top=min(ysum)-int(3.3*d);cv2.line(c,(xs_,max(ysum)),(xs_,top),col,max(1,lw//2))
                    for k in range({3:2,6:1,9:1}.get(base_u,0)): cv2.line(c,(xs_,top+int(k*0.7*d)),(xs_+int(0.9*d),top+int(k*0.7*d)+int(0.7*d)),col,lw)
                    if base_u in (9,18,36): cv2.circle(c,(x+int(1.1*d),min(ysum)-int(0.2*d)),max(2,int(0.15*d)),col,-1)
                if e.get('trip'): cv2.putText(c,'3',(int(x-0.3*d),int(min(ysum)-3.8*d)),cv2.FONT_HERSHEY_SIMPLEX,fs(d)*0.9,col,1)
                if e.get('tie_start'): cv2.ellipse(c,(int(x+1.2*d),int(max(ysum)+0.8*d)),(int(1.2*d),int(0.5*d)),0,0,180,col,max(1,lw//2))
            # Taktstrich am Ende des Takts (Mitte zwischen letztem Ereignis und nächstem Takt) und Taktnummer
            mx_=es[0].get('mx');xe=int(mx_[1]-0.15*d) if mx_ and mx_[1]>max(xs)+0.5*d else int(max(xs)+1.1*d);cv2.line(c,(xe,int(ys[0])),(xe,int(ys[4])),(0,160,0),max(1,lw//2))   # Taktstrich an der Taktgrenze der Vorlage
            mf=es[0].get('mflags') or {}
            def _rep(xb_,sgn):   # Wiederholungszeichen: dicker grüner Strich, zwei Punkte auf der Seite des wiederholten Teils
                cv2.line(c,(int(xb_),int(ys[0])),(int(xb_),int(ys[4])),(0,160,0),max(3,int(0.3*d)))
                for yy_ in (mid-0.5*d,mid+0.5*d): cv2.circle(c,(int(xb_+sgn*0.6*d),int(yy_)),max(3,int(0.22*d)),(0,160,0),-1)
            if mf.get('rep_end'): _rep(xe+0.2*d,-1)
            if mf.get('rep_start'): _rep(min(xs)-1.4*d,1)
            if mf.get('volta'): cv2.putText(c,f"{mf['volta']}.",(int(min(xs)-0.5*d),int(ys[0]-2.6*d)),cv2.FONT_HERSHEY_SIMPLEX,fs(d),(0,160,0),max(1,int(fs(d)*2)))
            cv2.putText(c,str(num0+mi)+('' if ok else '!'),(int(xe-0.6*d),int(ys[4]+2.3*d)),cv2.FONT_HERSHEY_SIMPLEX,fs(d)*1.1,(0,130,0) if ok else (0,0,230),max(1,int(fs(d)*2.2)))
            # Bassbuchstaben aus dem Ergebnis: Lage über die Einsätze
            if bass and mi<len(bass) and bass[mi]:
                for q,cn in bass[mi]:
                    xq=_xq(tx,q,L,xe-0.9*d)
                    cv2.putText(c,str(cn),(int(xq-0.4*d),int(ys[4]+4.6*d)),cv2.FONT_HERSHEY_SIMPLEX,fs(d)*1.2,(150,0,150),max(1,int(fs(d)*2.4)))
    if c.shape[1]>maxw: c=cv2.resize(c,(maxw,int(c.shape[0]*maxw/c.shape[1])),interpolation=cv2.INTER_AREA)
    return cv2.imencode('.jpg',c,[cv2.IMWRITE_JPEG_QUALITY,75])[1].tobytes()
def edit_page(g,meas,lens,num0,a0,bass=None,maxw=1600):
    """Daten zum Nachbearbeiten einer Seite im Wandler: entzerrtes Bild (grau, JPEG), Systeme (Linienhöhen) und je Takt die Ereignisse
    mit Lage im Bild, Notenwert, Balgrichtung und Köpfen (Reihe, Knopf). Koordinaten in Bildpunkten des gelieferten Bildes."""
    import base64
    k=min(1.0,maxw/g.shape[1]);im=cv2.resize(g,None,fx=k,fy=k,interpolation=cv2.INTER_AREA) if k<1 else g
    st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1])
    ms=[]
    for mi,(m,L) in enumerate(zip(meas,lens)):
        evs=[]
        for e in m:
            ev=dict(x=round(float(e['dbg'][2])*k,1),s=int(e['dbg'][0]),u=int(e['u']),r=bool(e['rest']),t=bool(e.get('trip')))
            if not e['rest']: ev.update(d=e.get('dir','Z'),h=[[int(r),int(n)] for r,n,_ in e['heads']],ts=bool(e.get('tie_start')),te=bool(e.get('tie_stop')))
            evs.append(ev)
        mx_=(m[0].get('mx') if m else None);ms.append(dict(n=num0+mi,i=a0+mi,L=int(L),e=evs,mx=[round(v*k,1) for v in mx_] if mx_ else None,f=(m[0].get('mflags') if m else None) or {},w=(m[0].get('words') if m else None) or [],b=[[int(q),str(c)] for q,c in (bass[mi] if bass and mi<len(bass) else [])]))
    jpg=cv2.imencode('.jpg',im,[cv2.IMWRITE_JPEG_QUALITY,60])[1].tobytes()
    return dict(img='data:image/jpeg;base64,'+base64.b64encode(jpg).decode(),w=im.shape[1],h=im.shape[0],st=[[round(float(y)*k,1) for y in ys] for ys in st],m=ms,tri=bool(D.LAST.get('tri_page')))
def split_double(path):
    """Quer liegende Doppelseite (zwei Notenseiten nebeneinander, z. B. aufgeschlagenes Heft): findet die Liniensuche auf der ganzen
    Seite keine Systeme, wird an der hellsten Spalte nahe der Mitte geteilt und jede Hälfte als eigene Seite gelesen. GS_SPLIT=0: aus"""
    if _os.environ.get('GS_SPLIT','1')!='1': return [path]
    try:
        im=Image.open(path);im=ImageOps.exif_transpose(im);W,H=im.size
        if W<1.2*H: return [path]
        g=load_gray(path);n0=len(G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1])) if g is not None else 0
        if n0 and W<1.3*H: return [path]
        a=np.array(im.convert('L')).astype(np.float32);col=(a<128).mean(0);lo,hi=int(W*0.4),int(W*0.6)
        x=lo+int(np.argmin(np.convolve(col[lo:hi],np.ones(15)/15,'same')))
        import tempfile;out=[];n1=0;rel=[]
        for box in ((0,0,x,H),(x,0,W,H)):
            fn=tempfile.mktemp(suffix='.png');im.crop(box).save(fn)
            g2=load_gray(fn)
            if g2 is not None:
                ss=G.staves(cv2.threshold(g2,240,255,cv2.THRESH_BINARY_INV)[1])
                if ss: out.append(fn);n1+=len(ss);rel.append([s_[2]/g2.shape[0] for s_ in ss])
        if not n0: return out or [path]
        # Doppelseite: die Hälften haben zusammen deutlich mehr Systeme als das ganze Blatt, und die Systeme links und rechts stehen
        # auf verschiedener Höhe (bei einer breiten Einzelseite laufen sie durch und liegen in beiden Hälften gleich hoch)
        if len(rel)==2 and n1>n0*1.5:
            same=sum(1 for a_ in rel[0] if any(abs(a_-b_)<0.015 for b_ in rel[1]))
            if same<0.5*min(len(rel[0]),len(rel[1])) or n0<0.6*min(len(rel[0]),len(rel[1])): return out
        return [path]
    except Exception:
        return [path]
def convert(paths,shift=0,title='Gescannte Griffschrift',table=None,bilder=None,style='auto',text=None,seiten_out=None):
    """Eine oder mehrere Seiten (Bilder) -> (MusicXML, Pins); Takte der Seiten hintereinander. table: Diskantbelegung des Wandlers (sonst 446 Normal).
    bilder: Liste, an die je Seite ein Kontrollbild (JPEG-Bytes) angehängt wird; die Pins bekommen dann 'offen' (Nummern der Takte, die nicht aufgehen)."""
    allm=[];Lm=None;bass=[];seiten=[]
    paths=[q for p_ in paths for q in split_double(p_)]   # quer liegende Doppelseiten in zwei Seiten teilen
    st_=set_style(style)
    if text is not None:   # Texterkennung auf der ersten Seite; bei „automatisch“ wählt ein erkannter Bearbeiter die Schreibweise
        try:
            g0=load_gray(paths[0]);T=page_text(g0,G.staves(cv2.threshold(g0,240,255,cv2.THRESH_BINARY_INV)[1])) if g0 is not None else {}
        except Exception: T={}
        text.update(T or {})
        if st_=='auto' and text.get('stil'): set_style(text['stil'])
        text['schreibweise']=STYLE[0]
    for p in paths:
        m,L=read(p,shift,table=table,bass=bass)
        if bilder is not None: seiten.append((_DBG['g'],len(allm),m,L))
        allm+=m;Lm=(Lm or [])+L
    if not allm: raise ValueError('Keine Griffschrift gefunden')
    xml,P=to_xml(allm,Lm,title,shift,bass)
    if bilder is not None:
        n0=0 if sum(e['u'] for e in allm[0])<Lm[0] else 1   # Nummern wie im MusicXML (Auftakt = 0)
        for g,a,m,L in seiten:
            if seiten_out is not None:
                try: seiten_out.append(edit_page(g,m,L,n0+a,a,bass[a:a+len(m)] if bass else None))
                except Exception: pass
            try: bilder.append(pruefbild(g,m,L,n0+a,bass=bass[a:a+len(m)] if bass else None) if _os.environ.get('GS_PB','ueber')=='ueber' else pruefbild_alt(g,m,L,n0+a))
            except Exception: pass
        P['offen']=[n0+i for i,(m,L) in enumerate(zip(allm,Lm)) if sum(e['u'] for e in m)!=L and not (i==0 and n0==0)]
        P['takte']=len(allm)
    set_style('auto')
    return xml,P
STEP=['C','C','D','D','E','F','F','G','G','A','A','B'];ALT=[0,1,0,1,0,0,1,0,1,0,1,0]
TRIPTYPE={2:('16th',0),4:('eighth',0),6:('eighth',1),8:('quarter',0),12:('quarter',1),16:('half',0),24:('half',1),32:('whole',0)}   # Notenart innerhalb einer Triole
TYPE={3:('16th',0),4:('eighth',0),6:('eighth',0),9:('eighth',1),12:('quarter',0),18:('quarter',1),24:('half',0),36:('half',1),48:('whole',0)}
def to_xml(meas,lens,title='Gescannte Griffschrift',shift=0,bass=None):
    if not isinstance(lens,list): lens=[lens]*len(meas)
    def ts(L): return (L//12,4) if L%12==0 else (L//6,8)
    x=['<?xml version="1.0" encoding="UTF-8"?>','<score-partwise version="3.1">',f'<work><work-title>{title}</work-title></work>',
       '<part-list><score-part id="P1"><part-name>Melodie</part-name></score-part></part-list>','<part id="P1">']
    pins=[];t0=0.0;bpins=[];cur_volta=[None]
    REP=_os.environ.get('GS_WIEDERHOLUNG','1')=='1'   # Wiederholungen, Volten, Haltebögen, Textzeichen ins MusicXML
    for mi,evs in enumerate(meas):
        Lm=lens[mi];beats,btype=ts(Lm)
        imp=' implicit="yes"' if mi==0 and sum(e['u'] for e in evs)<Lm else ''
        if mi==0: first_imp=bool(imp)
        x.append(f'<measure number="{mi if first_imp else mi+1}"{imp}>')
        if mi==0: x.append(f'<attributes><divisions>12</divisions><key><fifths>0</fifths></key><time><beats>{beats}</beats><beat-type>{btype}</beat-type></time><clef><sign>G</sign><line>2</line></clef></attributes>')
        elif lens[mi]!=lens[mi-1]: x.append(f'<attributes><time><beats>{beats}</beats><beat-type>{btype}</beat-type></time></attributes>')
        mf=(evs[0].get('mflags') or {}) if (evs and REP) else {}
        if mf.get('volta')==1 and not mf.get('volta_stop'):   # 1. Klammer ohne Ende: bis zum nächsten Wiederholungsstrich (höchstens 4 Takte), sonst 1 Takt
            end_=next((j for j in range(mi,min(len(meas),mi+4)) if meas[j] and (meas[j][0].get('mflags') or {}).get('rep_end')),mi)
            if meas[end_]: meas[end_][0].setdefault('mflags',{})['volta_stop']=True
            mf=(evs[0].get('mflags') or {})
        if mf.get('rep_start') or mf.get('volta'):
            bl='<barline location="left">'
            if mf.get('volta'): cur_volta[0]=mf['volta'];bl+=f'<ending number="{mf["volta"]}" type="start"/>'
            if mf.get('rep_start'): bl+='<repeat direction="forward"/>'
            x.append(bl+'</barline>')
        if evs and REP:
            for w_ in evs[0].get('words',[]): x.append(f'<direction placement="above"><direction-type><words>{w_}</words></direction-type></direction>')
        pos=0
        for e in evs:
            tr=e.get('trip',e['u']==4);ty,dot=(TRIPTYPE if tr else TYPE).get(e['u'],TYPE.get(e['u'],('quarter',0)));tm='<time-modification><actual-notes>3</actual-notes><normal-notes>2</normal-notes></time-modification>' if tr else ''
            if e['rest']: x.append(f'<note><rest/><duration>{e["u"]}</duration><type>{ty}</type>{"<dot/>"*dot}{tm}</note>')
            else:
                for q,(r,n,p) in enumerate(e['heads']):
                    al=f'<alter>1</alter>' if ALT[p%12] else ''
                    ti=(('<tie type="stop"/>' if e.get('tie_stop') else '')+('<tie type="start"/>' if e.get('tie_start') else '')) if REP else ''
                    tn=(('<notations>'+('<tied type="stop"/>' if e.get('tie_stop') else '')+('<tied type="start"/>' if e.get('tie_start') else '')+'</notations>') if ti else '')
                    x.append(f'<note>{"<chord/>" if q else ""}<pitch><step>{STEP[p%12]}</step>{al}<octave>{p//12-1}</octave></pitch><duration>{e["u"]}</duration>{ti}<type>{ty}</type>{"<dot/>"*dot}{tm}{tn}</note>')
                    pins.append([round(t0+pos/12,4),p,r,n,e['dir']])
            pos+=e['u']
        if bass and mi<len(bass): bpins+=[[round(t0+q/12,4),c] for q,c in bass[mi]]
        if mf.get('rep_end') or mf.get('volta_stop'):
            bl='<barline location="right"><bar-style>light-heavy</bar-style>' if mf.get('rep_end') else '<barline location="right">'
            if mf.get('volta_stop') and cur_volta[0]: bl+=f'<ending number="{cur_volta[0]}" type="stop"/>';cur_volta[0]=None
            if mf.get('rep_end'): bl+='<repeat direction="backward"/>'
            x.append(bl+'</barline>')
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
