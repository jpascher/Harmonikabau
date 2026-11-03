# Fuchs-Seiten mit den MIDI-Dateien vergleichen (wie fuchsread.check): Akkorde genau gleich und Anteil der gelesenen Töne, die in der MIDI vorkommen
# Aufruf: python3 fuchs_pruefen_toene.py ../gsleser   (Seiten wie bisher in /tmp/claude-0/fuchs/fuchs/<stück>gs/GS*.gif|png samt .mid)
import sys,glob,os,json,numpy as np,cv2
GS=sys.argv[1];sys.path.insert(0,GS)
import fuchsread as F
from align import upscale
T=dict(n=0,paired=0,exact=0,tone_in=0,tone_tot=0);rows=[]
for f in sorted(glob.glob('/tmp/claude-0/fuchs/fuchs/*/')):
    if 'iliana' in f: continue
    gs=glob.glob(f+'GS*');ms=glob.glob(f+'*.mid')
    if not gs or not ms: continue
    k=f.rstrip('/').split('/')[-1];img=f'/tmp/claude-0/fuchs/u_{k}.png'
    if not os.path.exists(img): upscale(gs[0],img)
    res,_=F.read_page(img);B=F.midi_chords(ms[0]);best=None
    for sh in range(-14,15):
        A=[{p+sh for p in r['pitches']} for r in res];sc,pr=F.dp(A,B)
        if best is None or sc>best[0]: best=(sc,sh,pr,A)
    sc,sh,pr,A=best;ex=sum(A[i]==B[j] for i,j in pr)
    ti=sum(len(A[i]&B[j]) for i,j in pr);tt=sum(len(A[i]) for i,j in pr)
    T['n']+=len(A);T['paired']+=len(pr);T['exact']+=ex;T['tone_in']+=ti;T['tone_tot']+=tt
    rows.append((k,len(A),len(pr),ex,ti,tt))
print(os.path.basename(GS.rstrip('/')),T,'gleich',round(T['exact']/max(1,T['paired']),3),'Töne in MIDI',round(T['tone_in']/max(1,T['tone_tot']),3))
json.dump(rows,open('/tmp/claude-0/fuchs/rows_'+os.path.basename(GS.rstrip('/'))+'.json','w'))
