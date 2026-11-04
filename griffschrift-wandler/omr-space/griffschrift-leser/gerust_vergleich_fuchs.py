# Ende-zu-Ende gegen die bekannten Fuchs-Griffe (midi/vm/*.pins.json, mit Zeit und Ton): Akkorde ganz gleich, Töne gefunden,
# Töne richtig unter den gelesenen, und Rhythmus = Abstand zum nächsten Einsatz richtig (bei aufeinanderfolgenden zugeordneten Akkorden).
# Seiten wie gewohnt als /tmp/claude-0/fuchs/u_<stück>.png (hochgerechnete GS*.gif von volksmusik.cc, nicht im Repo).
# Aufruf: GS_GERUEST=0 python3 gerust_vergleich_fuchs.py leser.jsonl ; GS_GERUEST=1 GS_GCACHE=/tmp/claude-0/gcache python3 gerust_vergleich_fuchs.py geruest.jsonl
import sys,glob,os,json,collections,numpy as np
import os as _o;sys.path.insert(0,_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..','gsleser'))
import gs2xml as X
PINS=_o.path.join(_o.path.dirname(_o.path.abspath(__file__)),'..','..','midi','vm');out=sys.argv[1]
done=set()
if os.path.exists(out): done={json.loads(l)['k'] for l in open(out)}
def chords(n):
    m=collections.defaultdict(set)
    for t,p,*_ in n: m[round(float(t)*12)].add(int(p))
    ts=sorted(m);return [(t,m[t]) for t in ts]
def align(A,B):
    n,mm=len(A),len(B);S=np.zeros((n+1,mm+1));S[:,0]=-np.arange(n+1)*.5;S[0,:]=-np.arange(mm+1)*.5;bt={}
    sc=lambda a,b:(len(a&b)*2-len(a^b)*.5) if a&b else -3
    for i in range(1,n+1):
        for j in range(1,mm+1):
            c=[(S[i-1,j-1]+sc(A[i-1][1],B[j-1][1]),'m'),(S[i-1,j]-.5,'a'),(S[i,j-1]-.5,'b')];S[i,j],bt[i,j]=max(c)
    i,j,pr=n,mm,[]
    while i>0 and j>0:
        k=bt[i,j]
        if k=='m': pr.append((i-1,j-1));i-=1;j-=1
        elif k=='a': i-=1
        else: j-=1
    return pr[::-1]
ONLY=set(json.loads(l)['k'] for l in open(os.environ['ONLY'])) if os.environ.get('ONLY') else None
for f in sorted(glob.glob('/tmp/claude-0/fuchs/u_*.png')):
    k=os.path.basename(f)[2:-4]
    if ONLY is not None and k not in ONLY: continue
    if k in done: continue
    gt=json.load(open(f'{PINS}/{k}.pins.json'));B=chords(gt['n'])
    try: xml,P=X.convert([f],shift=gt.get('sh',0)); A0=chords(P['n'])
    except Exception as e: A0=[]
    best=None
    for sh in range(-3,4):
        A=[(t,{p+sh for p in s}) for t,s in A0];pr=align(A,B);v=sum(len(A[i][1]&B[j][1]) for i,j in pr)
        if best is None or v>best[0]: best=(v,pr,A)
    v,pr,A=best;ex=sum(A[i][1]==B[j][1] for i,j in pr)
    ioi=0;nio=0
    for (i1,j1),(i2,j2) in zip(pr,pr[1:]):
        if i2==i1+1 and j2==j1+1: nio+=1;ioi+=(A[i2][0]-A[i1][0])==(B[j2][0]-B[j1][0])
    r=dict(k=k,gt_ch=len(B),gt_n=sum(len(s) for _,s in B),out_n=sum(len(s) for _,s in A),paired=len(pr),exact=ex,tones=v,ioi_ok=ioi,ioi_n=nio,stats=X._DBG.get('gerust'))
    open(out,'a').write(json.dumps(r)+'\n')
