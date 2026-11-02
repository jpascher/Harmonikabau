# Vergleicht zwei Fassungen des Lesers auf Bildern ohne Dreiecke (Wandler-Griffschrift): Aufruf python3 vergleich_alt_neu.py 'wand/*.png' alter_gsleser neuer_gsleser
import sys,glob,collections,importlib
def run(gs):
    for m in ['gs2xml','detect','gsread','assemble','align','fuchsread','rhythm','measures','train']: sys.modules.pop(m,None)
    sys.path.insert(0,gs);import gs2xml as X;sys.path.pop(0)
    out={}
    for f in sorted(glob.glob(sys.argv[1])):
        try: meas,L=X.read(f,-2)
        except Exception as e: out[f]=[];continue
        out[f]=[sorted(h[:2] for h in e['heads']) for m in meas for e in m if not e['rest']]
    return out
a=run(sys.argv[2]);b=run(sys.argv[3])
same=diff=0;tri=0
for f in a:
    for x,y in zip(a[f],b[f]): same+=x==y;diff+=x!=y
    tri+=sum(1 for ev in b[f] for r,n in ev if r==4)-sum(1 for ev in a[f] for r,n in ev if r==4)
    if len(a[f])!=len(b[f]): print('Anzahl anders',f,len(a[f]),len(b[f]))
print('Ereignisse gleich',same,'verschieden',diff,'Reihe-4-Köpfe neu minus alt',tri)
