# Griffschrift-Bild -> MusicXML (Normalnoten) + Griffe für den Wandler (Pins: Zeit, Ton, Reihe, Knopf, Balgrichtung)
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2,sys,json,pickle,collections,numpy as np,os as _os
from PIL import Image
sys.path.insert(0,HERE)
import gsread as G, detect as D
from align import chords_det
from assemble import button,pitch
from fuchsread import druck_bars
from rhythm import strip,beam_feats,dot_feat
from measures import barlines,rests,fit_measure,UNIT
RCLF=pickle.load(open((HERE if _os.path.exists(HERE+'/rhythm.pkl') else _os.path.join(HERE,'..','gsleser'))+'/rhythm.pkl','rb'))
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
def load_gray(path,target_d=18.0):
    im=Image.open(path).convert('RGBA');bg=Image.new('RGBA',im.size,(255,255,255,255));bg.alpha_composite(im);g=np.array(bg.convert('L'))
    a=skew_angle(g)
    if abs(a)>=0.2:   # schräge Vorlage gerade drehen, sonst findet die Liniensuche die Systeme nicht
        h,w=g.shape;g=cv2.warpAffine(g,cv2.getRotationMatrix2D((w/2,h/2),a,1.0),(w,h),flags=cv2.INTER_CUBIC,borderValue=255)
    hell=g[g>200].astype(np.float32)   # Foto erkennen: grauer Hintergrund oder Rauschen im hellen Bereich (robust gegen geglättete Kanten sauberer Bilder)
    foto=bool(np.median(g)<250 or (hell.size and float(np.median(np.abs(hell-np.median(hell))))>=2))
    # bei Fotos die Helligkeit schon vor der Liniensuche ausgleichen: auf grauem Hintergrund liegt sonst alles unter der Linienschwelle
    gs=cv2.divide(g,cv2.GaussianBlur(g,(0,0),max(8.0,min(g.shape)/60)),scale=255) if foto else g
    for f in (1,2,3):
        gg=cv2.resize(gs,None,fx=f,fy=f,interpolation=cv2.INTER_CUBIC) if f>1 else gs
        st=G.staves(cv2.threshold(gg,240,255,cv2.THRESH_BINARY_INV)[1])
        if st:
            d=float(np.median([np.mean(np.diff(s)) for s in st]))
            g=cv2.resize(g,None,fx=target_d/(d/f),fy=target_d/(d/f),interpolation=cv2.INTER_CUBIC)
            # Hintergrund ausgleichen, aber nur bei Fotos (grauer oder verrauschter Hintergrund); saubere Scans und Bilder bleiben unverändert
            if foto: g=cv2.divide(g,cv2.GaussianBlur(g,(0,0),20),scale=255)
            return g
    return None
def read(path,shift=0,work=None,table=None):
    """shift: Stimmung wie im Wandler (A-D-G-C = 0, G-C-F-B = -2 …)."""
    g=load_gray(path)
    if g is None: raise ValueError('Keine Notenlinien gefunden')
    if work is None:
        import tempfile;work=_tmp=tempfile.mktemp(suffix='.png')
    cv2.imwrite(work,g);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV)
    st,bars_d=druck_bars(g);det=D.detect(work,1);C=chords_det(det)
    off=shift+2   # Tabelle in assemble.py gilt für G-C-F-B (shift -2)
    items=[]   # je System: Akkorde und Pausen
    feats=[]
    # Zeilenanfang je System (Taktangabe im ersten System nicht als Kopf lesen)
    x_start={}
    for si,ys in enumerate(st):
        row=np.where(bw[int(round(ys[0]))]>0)[0];x_start[si]=float(row.min()) if len(row) else 0.0
    for c in C:
        si=c[0]['st'];x=float(np.median([h['x'] for h in c]))
        if si==0 and x<x_start[0]+6*(st[0][4]-st[0][0])/4: continue
        dr='D' if any(a-3<=x<=b+3 for a,b,_ in bars_d.get(si,[])) else 'Z'
        ks=[button(h['p'],h['t']) for h in c];ps=[pitch(r,n,dr,table) for r,n in ks]
        hs=[(h,k,p+off) for h,k,p in zip(c,ks,ps) if p]
        if not hs: continue
        ys=st[si];d=(ys[4]-ys[0])/4;nl=G.remove_lines(bw,ys,d)
        fill=float(np.mean([h['fill'] for h,_,_ in hs]))
        feats.append(np.concatenate([strip(nl,x,ys[2],d,fill),np.array(beam_feats(nl,x,[h['y'] for h,_,_ in hs],d)+[dot_feat(nl,x,[h['y'] for h,_,_ in hs],d)])*20]).astype(np.float32))
        items.append(dict(kind='n',st=si,x=x,dir=dr,heads=[(k[0],k[1],p) for _,k,p in hs]))
    P=RCLF.predict_proba(np.array(feats)) if feats else np.zeros((0,len(RCLF.classes_)))
    P9=np.zeros((len(P),len(UNIT)));P9[:,RCLF.classes_]=P;lp=np.log(np.clip(P9,1e-6,1))
    for it,row in zip(items,lp): it['lp']=row
    meas=[]
    for si,ys in enumerate(st):
        d=(ys[4]-ys[0])/4;its=sorted([it for it in items if it['st']==si],key=lambda it:it['x'])
        hxy=[]
        for it in its:
            for r_,n_,p_ in it['heads']: hxy.append((it['x'],0))
        heads_xy=[(h['x'],h['y']) for c in C if c[0]['st']==si for h in c]
        bl=barlines(bw,ys,d,[x for x,_ in heads_xy]);nol=G.stems_off(G.remove_lines(bw,ys,d),d)
        rs=[dict(kind='r',st=si,x=x,u=u) for x,u in rests(bw,nol,ys,d,heads_xy,bl)]
        edges=[0]+bl+[1e9]
        for a,b in zip(edges,edges[1:]):
            ev=sorted([it for it in its if a<it['x']<b]+[r for r in rs if a<r['x']<b],key=lambda e:e['x'])
            grp=[]
            for e in ev:
                if e['kind']=='n' and grp and grp[-1][0]['kind']=='n' and e['x']-grp[-1][-1]['x']<1.3*d: grp[-1].append(e)
                else: grp.append([e])
            if grp: meas.append((si,grp))
    def gl(gr): return sum(e['lp'] for e in gr)
    def length(m,ch=None):
        t=0
        for q,gr in enumerate(m):
            t+=gr[0]['u'] if gr[0]['kind']=='r' else UNIT[ch[q] if ch else int(np.argmax(gl(gr)))]
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
            f=fit_measure(PG,[gr[0]['u'] for gr in m if gr[0]['kind']=='r'],Lm)
            if f: ch=f
        evs=[];k=0
        for gr in m:
            if gr[0]['kind']=='r': evs.append(dict(rest=True,u=gr[0]['u']))
            else:
                heads=sorted({h for e in gr for h in e['heads']},key=lambda h:-h[2]);evs.append(dict(rest=False,u=UNIT[ch[k]],dir=gr[0]['dir'],heads=heads));k+=1
        out.append(evs);lens.append(Lm)
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
    allm=[];Lm=None
    for p in paths:
        m,L=read(p,shift,table=table);allm+=m;Lm=(Lm or [])+L
    if not allm: raise ValueError('Keine Griffschrift gefunden')
    return to_xml(allm,Lm,title,shift)
STEP=['C','C','D','D','E','F','F','G','G','A','A','B'];ALT=[0,1,0,1,0,0,1,0,1,0,1,0]
TYPE={3:('16th',0),4:('eighth',0),6:('eighth',0),9:('eighth',1),12:('quarter',0),18:('quarter',1),24:('half',0),36:('half',1),48:('whole',0)}
def to_xml(meas,lens,title='Gescannte Griffschrift',shift=0):
    if not isinstance(lens,list): lens=[lens]*len(meas)
    def ts(L): return (L//12,4) if L%12==0 else (L//6,8)
    x=['<?xml version="1.0" encoding="UTF-8"?>','<score-partwise version="3.1">',f'<work><work-title>{title}</work-title></work>',
       '<part-list><score-part id="P1"><part-name>Melodie</part-name></score-part></part-list>','<part id="P1">']
    pins=[];t0=0.0
    for mi,evs in enumerate(meas):
        Lm=lens[mi];beats,btype=ts(Lm)
        imp=' implicit="yes"' if mi==0 and sum(e['u'] for e in evs)<Lm else ''
        if mi==0: first_imp=bool(imp)
        x.append(f'<measure number="{mi if first_imp else mi+1}"{imp}>')
        if mi==0: x.append(f'<attributes><divisions>12</divisions><key><fifths>0</fifths></key><time><beats>{beats}</beats><beat-type>{btype}</beat-type></time><clef><sign>G</sign><line>2</line></clef></attributes>')
        elif lens[mi]!=lens[mi-1]: x.append(f'<attributes><time><beats>{beats}</beats><beat-type>{btype}</beat-type></time></attributes>')
        pos=0
        for e in evs:
            ty,dot=TYPE.get(e['u'],('quarter',0));tm='<time-modification><actual-notes>3</actual-notes><normal-notes>2</normal-notes></time-modification>' if e['u']==4 else ''
            if e['rest']: x.append(f'<note><rest/><duration>{e["u"]}</duration><type>{ty}</type>{"<dot/>"*dot}{tm}</note>')
            else:
                for q,(r,n,p) in enumerate(e['heads']):
                    al=f'<alter>1</alter>' if ALT[p%12] else ''
                    x.append(f'<note>{"<chord/>" if q else ""}<pitch><step>{STEP[p%12]}</step>{al}<octave>{p//12-1}</octave></pitch><duration>{e["u"]}</duration><type>{ty}</type>{"<dot/>"*dot}{tm}</note>')
                    pins.append([round(t0+pos/12,4),p,r,n,e['dir']])
            pos+=e['u']
        x.append('</measure>');t0+=pos/12
    x+=['</part>','</score-partwise>']
    return '\n'.join(x),{'sh':shift,'n':pins}
if __name__=='__main__':
    sh=SHIFT[sys.argv[2] if len(sys.argv)>2 else 'ADGC'];meas,Lm=read(sys.argv[1],sh)
    xml,pins=to_xml(meas,Lm,shift=sh)
    open(sys.argv[3] if len(sys.argv)>3 else HERE+'/out.musicxml','w').write(xml)
    json.dump(pins,open((sys.argv[3] if len(sys.argv)>3 else HERE+'/out.musicxml')+'.pins.json','w'))
    print('Takte',len(meas),'Taktlängen',sorted(set(l/12 for l in Lm)),'Viertel; Noten',sum(1 for m in meas for e in m if not e['rest']),'Pausen',sum(1 for m in meas for e in m if e['rest']))
