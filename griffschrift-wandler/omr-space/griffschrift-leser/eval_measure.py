import cv2,sys,json,glob,pickle,collections,numpy as np,mido
sys.path.insert(0,'/tmp/claude-0/gsr');from rhythm import *;from measures import *;import gsread as G
clf=pickle.load(open('/tmp/claude-0/gsr/rhythm.pkl','rb'))
R=json.load(open('/tmp/claude-0/gsr/midicheck_e.json'));keys=json.load(open('/tmp/claude-0/gsr/midi_test.json'))
def midi_len(k):
    m=mido.MidiFile(glob.glob(f'/tmp/claude-0/fuchs/fuchs/{k}/*.mid')[0]);best=None
    for t in m.tracks:
        on={};ns=[];tt=0
        for x in t:
            tt+=x.time
            if x.type=='note_on' and x.velocity>0: on[x.note]=tt
            elif x.type in('note_off','note_on') and x.note in on: s=on.pop(x.note);ns.append((s,x.note,tt-s))
        if ns and (best is None or np.mean([n for _,n,_ in ns])>np.mean([n for _,n,_ in best])): best=ns
    q=m.ticks_per_beat/8;grp={}
    for s,n,l in best: grp.setdefault(round(s/q),{})[n]=l/m.ticks_per_beat
    return [grp[k] for k in sorted(grp)]
T=collections.Counter()
for k in keys:
    r=R[k];sh=r['sh'];B=midi_len(k);g=cv2.imread(r['img'],0);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV);st,view=page_views(g)
    lab={}
    for i,j in r['pairs']:
        a={p+sh for p in r['res'][i]['pitches']};b=B[j]
        if a and a<=set(b):
            L=[b[p] for p in a]
            if max(L)-min(L)<=0.1: lab[i]=dclass(float(np.mean(L))/0.95)
    F=[];
    for c in r['res']:
        nl,d,mid=view(c['st']);fill=float(np.mean([h['fill'] for h in c['heads']]))
        F.append(np.concatenate([strip(nl,c['x'],mid,d,fill),np.array(beam_feats(nl,c['x'],[h['y'] for h in c['heads']],d)+[dot_feat(nl,c['x'],[h['y'] for h in c['heads']],d)])*20]).astype(np.float32))
    P=clf.predict_proba(np.array(F));pred=P.argmax(1).copy();P9=np.zeros((len(P),len(UNIT)));P9[:,clf.classes_]=P
    # Takte bilden
    meas=[]
    for si,ys in enumerate(st):
        d=(ys[4]-ys[0])/4;idx=[i for i,c in enumerate(r['res']) if c['st']==si]
        hs=[(h['x'],h['y']) for i in idx for h in r['res'][i]['heads']];bars=barlines(bw,ys,d,[x for x,_ in hs])
        nol=G.stems_off(G.remove_lines(bw,ys,d),d);rs=rests(bw,nol,ys,d,hs,bars);edges=[0]+bars+[1e9]
        for a,b in zip(edges,edges[1:]):
            ci=sorted([i for i in idx if a<r['res'][i]['x']<b],key=lambda i:r['res'][i]['x']);ru=[u for x,u in rs if a<x<b]
            grp=[]   # dicht nebeneinander liegende Akkorde (Sekunden) sind ein Ereignis
            for i in ci:
                if grp and r['res'][i]['x']-r['res'][grp[-1][-1]]['x']<1.3*d: grp[-1].append(i)
                else: grp.append([i])
            if ci or ru: meas.append((grp,ru))
    lp=np.log(np.clip(P9,1e-6,1))
    def gpred(gr): return int(np.argmax(sum(lp[i] for i in gr)))
    sums=[sum(UNIT[gpred(gr)] for gr in G_)+sum(ru) for G_,ru in meas];Lm=collections.Counter(sums).most_common(1)[0][0]
    new=pred.copy();fixed=0;ok=0
    for G_,ru in meas:
        for gr in G_:
            for i in gr: new[i]=gpred(gr)
        if sum(UNIT[gpred(gr)] for gr in G_)+sum(ru)==Lm: ok+=1;continue
        PG=np.array([np.exp(sum(lp[i] for i in gr)-max(sum(lp[i] for i in gr))) for gr in G_]) if G_ else np.zeros((0,9))
        PG=PG/PG.sum(1,keepdims=True) if len(PG) else PG
        ch=fit_measure(PG,ru,Lm)
        if ch:
            for gr,kk in zip(G_,ch):
                for i in gr: new[i]=kk
            fixed+=1;ok+=1
    if k==sys.argv[1] if len(sys.argv)>1 else False:
        v=cv2.cvtColor(g,cv2.COLOR_GRAY2BGR);NM=['16','8T','8','8.','4','4.','2','2.','1']
        for si,ys in enumerate(st):
            d=(ys[4]-ys[0])/4;idx=[i for i,c in enumerate(r['res']) if c['st']==si];hs=[(h['x'],h['y']) for i in idx for h in r['res'][i]['heads']]
            bars=barlines(bw,ys,d,[x for x,_ in hs]);nol=G.stems_off(G.remove_lines(bw,ys,d),d)
            for x in bars: cv2.line(v,(int(x),int(ys[0])-10),(int(x),int(ys[4])+10),(0,0,255),3)
            for x,u in rests(bw,nol,ys,d,hs,bars): cv2.putText(v,'P'+str(u),(int(x)-10,int(ys[0])-30),0,.7,(255,0,0),2)
            for i in idx:
                c=r['res'][i];cv2.putText(v,NM[new[i]],(int(c['x'])-10,int(ys[4])+50),0,.7,(0,0,255) if i in lab and lab[i]!=new[i] else (0,130,0),2)
                if i in lab and lab[i]!=new[i]: cv2.putText(v,NM[lab[i]],(int(c['x'])-10,int(ys[4])+75),0,.7,(200,0,200),2)
        cv2.imwrite('/tmp/claude-0/gsr/em.png',cv2.resize(v[:1600],None,fx=.5,fy=.5))
    a0=sum(pred[i]==y for i,y in lab.items());a1=sum(new[i]==y for i,y in lab.items())
    print(k,'Takt',Lm/12,'Takte',len(meas),'stimmig',ok,'angepasst',fixed,'vorher',round(a0/len(lab),3),'nachher',round(a1/len(lab),3),flush=True)
    T['n']+=len(lab);T['a0']+=a0;T['a1']+=a1
print('gesamt vorher',round(T['a0']/T['n'],3),'nachher',round(T['a1']/T['n'],3))
