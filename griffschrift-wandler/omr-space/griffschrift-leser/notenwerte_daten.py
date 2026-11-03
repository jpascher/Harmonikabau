# Notenwert-Lerndaten aus echten Seiten: Philipps Heft (Pins) und Fuchs-Seiten (MIDI); Merkmale wie gs2xml, Ziel = Abstand zum nächsten Einsatz
import sys,glob,json,collections,numpy as np,mido
GS=sys.argv[1];sys.path.insert(0,GS);src=sys.argv[2];out=sys.argv[3]
import gs2xml as X, fuchsread as F
from align import align
from rhythm import DUR
OFF={1:-1,2:0,3:1,4:2}
def key(r,n): return (2*(n-6)+OFF[r],'o' if r<3 else ('x' if r==3 else 't'))
def label(dur,trip):
    base=dur*1.5 if trip else dur
    k=int(np.argmin([abs(np.log(base/d)) for d in DUR]))
    return k if abs(np.log(base/DUR[k]))<0.08 else None
def page_events(path,shift,table=None):
    meas,_=X.read(path,shift,table=table);fe=X._DBG['feats'];seq=[]
    for m in meas:
        for e in m: seq.append(e)
    return seq,fe
XS=[];YS=[];KS=[];RS=[]
def emit(seq,fe,pairs,TS,piece):
    # pairs: Index in seq (nur Noten) -> Index der Zeitpunkte
    for q,e in enumerate(seq):
        if e['rest'] or q not in pairs: continue
        j=pairs[q]
        if q+1>=len(seq) or seq[q+1]['rest'] or (q+1) not in pairs or pairs[q+1]!=j+1: continue
        y=label(TS[j+1]-TS[j],e['mtrip'])
        if y is None: continue
        for fi in e['fis']: XS.append(fe[fi]);YS.append(y);KS.append(piece);RS.append(e['u'])
if src=='P':
    sys.path.insert(0,'/home/claude/philipp');from philipp_pruefen import PAGES,PINS,pin_chords
    D450={"druck":{1:[98,30,25,4,7,12,16,19,24,28,31,36,40,43,32,98],2:[98,98,20,5,9,12,17,21,24,29,33,36,41,45,48,98,98],3:[98,98,8,10,14,17,22,26,29,34,38,41,46,50,98,98],4:[98,98,13,15,19,22,27,31,34,39,43,46,51,98,98]},
          "zug":{1:[98,25,30,7,11,14,17,21,23,26,29,33,35,38,41,98],2:[98,98,18,10,12,16,19,22,24,28,31,34,38,40,43,98,98],3:[98,98,9,15,17,21,24,27,29,33,36,39,43,45,98,98],4:[98,98,14,20,22,26,29,32,34,38,41,44,48,98,98]}}
    for k in PAGES:
        seq=[];fe=[]
        for p in PAGES[k]:
            s_,f_=page_events(f'/home/claude/philipp/pages/p-{p:02d}.png',-2,D450)
            for e in s_:
                if not e['rest']: e['fis']=[i+len(fe) for i in e['fis']]
            seq+=s_;fe+=list(f_)
        notes=[q for q,e in enumerate(seq) if not e['rest']]
        A=[{key(r,n) for r,n,_ in seq[q]['heads']} for q in notes]
        B=pin_chords(PINS[k]);TS=sorted({round(t,3) for t,*_ in PINS[k]})
        pr=align(A,[{key(r,n) for r,n,_ in b} for b in B]);pairs={notes[i]:j for i,j in pr}
        emit(seq,fe,pairs,TS,'P:'+k);print(k,len(XS),flush=True)
else:
    for f in sorted(glob.glob('/home/claude/fuchs/zip/x/fuchs/*/')):
        if 'iliana' in f: continue
        gs=sorted(glob.glob(f+'GS*'));ms=glob.glob(f+'*.mid')
        if not gs or not ms: continue
        try: seq,fe=page_events(gs[0],0)
        except Exception as ex: print('FEHLER',f,ex);continue
        m=mido.MidiFile(ms[0]);best=None
        for t in m.tracks:
            ns=[];tt=0
            for x in t:
                tt+=x.time
                if x.type=='note_on' and x.velocity>0: ns.append((tt/m.ticks_per_beat,x.note))
            if ns and (best is None or np.mean([n for _,n in ns])>np.mean([n for _,n in best])): best=ns
        g={}
        for tt,n in best: g.setdefault(round(tt*24)/24,set()).add(n)
        TS=sorted(g);B=[g[t] for t in TS]
        notes=[q for q,e in enumerate(seq) if not e['rest']];A0=[{p for _,_,p in seq[q]['heads']} for q in notes];bb=None
        for sh in range(-14,15):
            A=[{p+sh for p in a} for a in A0];sc,pr=F.dp(A,B)
            if bb is None or sc>bb[0]: bb=(sc,pr)
        pairs={notes[i]:j for i,j in bb[1]}
        emit(seq,fe,pairs,TS,'F:'+f.rstrip('/').split('/')[-1]);print(f.rstrip('/').split('/')[-1],len(XS),flush=True)
np.savez_compressed(out,X=np.array(XS,np.float32),Y=np.array(YS),K=np.array(KS),R=np.array(RS))
