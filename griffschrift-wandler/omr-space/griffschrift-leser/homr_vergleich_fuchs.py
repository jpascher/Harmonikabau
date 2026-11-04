# Kopf-Vergleich gegen die bekannten Fuchs-Griffe (midi/vm/*.pins.json): ohne und mit homr-Zerlegungsnetz (GS_HOMR=1 in der Umgebung)
# Seiten wie gewohnt in /tmp/claude-0/fuchs/fuchs/<stück>gs/GS*.gif (von volksmusik.cc, nicht im Repo).
# Aufruf aus diesem Ordner: GS_HOMR=0 python3 homr_vergleich_fuchs.py ; GS_HOMR=1 GS_HOMRADD=0.8 GS_HOMRCACHE=/tmp/claude-0/hcache python3 homr_vergleich_fuchs.py
import sys,glob,os,json,time
import os as _o;GS=_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..','gsleser');sys.path.insert(0,GS)
import detect as D
from align import chords_det,chords_pin,align,upscale
PINS=_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..','..','midi','vm')
T=dict(pins=0,det=0,pairs=0,exact=0,tp=0,fp=0,fn=0);rows={};t0=time.time()
for f in sorted(glob.glob('/tmp/claude-0/fuchs/fuchs/*gs/')):
    k=f.rstrip('/').split('/')[-1][:-2];gs=glob.glob(f+'GS*')
    img=f'/tmp/claude-0/fuchs/u_{k}.png'
    if not os.path.exists(img): upscale(gs[0],img)
    det=D.detect(img,1);C=chords_det(det);A=[{(h['p'],h['t']) for h in c} for c in C]
    B=chords_pin(json.load(open(f'{PINS}/{k}.pins.json'))['n']);pr=align(A,B)
    ex=sum(A[i]==B[j] for i,j in pr);tp=sum(len(A[i]&B[j]) for i,j in pr)
    nA=sum(len(a) for a in A);nB=sum(len(b) for b in B)
    r=dict(pins=len(B),det=len(A),pairs=len(pr),exact=ex,tp=tp,fp=nA-tp,fn=nB-tp);rows[k]=r
    for x in r: T[x]+=r[x]
print('GS_HOMR',os.environ.get('GS_HOMR','0'),T,'Akkorde ganz gleich %.3f'%(T['exact']/T['pins']),'Köpfe gefunden %.3f'%(T['tp']/(T['tp']+T['fn'])),'Genauigkeit %.3f'%(T['tp']/(T['tp']+T['fp'])),'Zeit %.0fs'%(time.time()-t0))
json.dump(rows,open(f"/tmp/claude-0/fuchs/homr_rows_{os.environ.get('GS_HOMR','0')}{os.environ.get('TAG','')}.json",'w'))
