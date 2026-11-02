import json,sys,collections
R=json.load(open(sys.argv[1]));keys=json.load(open(sys.argv[2])) if len(sys.argv)>2 else list(R);T=collections.Counter()
for k in keys:
    r=R[k];sh=r['sh']
    for i,j in r['pairs']:
        a={p+sh for p in r['res'][i]['pitches']};b=set(r['B'][j]);T['gleich' if a==b else 'Teilmenge' if a and a<b else 'falsch']+=1
    T['Akkorde']+=r['n']
n=T['gleich']+T['Teilmenge']+T['falsch'];print(dict(T),'MIDI-stimmig',round((T['gleich']+T['Teilmenge'])/n,3))
