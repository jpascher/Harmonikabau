import cv2,sys,json,glob,numpy as np,mido
sys.path.insert(0,'/tmp/claude-0/gsr');from rhythm import *;import detect as D;from align import chords_det
X=[];Y=[];K=[];NOTE={};PC=[]
# synthetisch: Griffe mit bekannter Dauer
test=set(json.load(open('/tmp/claude-0/gsr/test.json')))
for f in sorted(glob.glob('/tmp/claude-0/gs/z*.json')+glob.glob('/tmp/claude-0/gs2/f*.json')):
    gt=json.load(open(f));sc=gt['scale'];g=cv2.imread(gt['img'],0);st,view=page_views(g)
    if not st: continue
    nl,d,mid=view(0)
    for gr in gt['grips']:
        hs=[h for h in gt['heads'] if abs(h['x']-gr['x'])<3]
        fill=float(np.mean([h.get('fill',True) for h in hs])) if hs else 1.0
        X.append(np.concatenate([strip(nl,gr['x']*sc,mid,d,fill),np.array(beam_feats(nl,gr['x']*sc,[h['y']*sc for h in hs] or [mid],d)+[dot_feat(nl,gr['x']*sc,[h['y']*sc for h in hs] or [mid],d)])*20]).astype(np.float32));Y.append(dclass(gr['dur']));K.append('T' if f in test else 'S')
print('synthetisch',len(Y),flush=True)
# echt: Fuchs-Akkorde, die zur MIDI passen; Dauer = Länge der passenden MIDI-Töne
R=json.load(open('/tmp/claude-0/gsr/midicheck_e.json'));mt=set(json.load(open('/tmp/claude-0/gsr/midi_test.json')))
for k,r in R.items():
    mf=glob.glob(f'/tmp/claude-0/fuchs/fuchs/{k}/*.mid')[0];m=mido.MidiFile(mf);best=None
    for t in m.tracks:
        on={};ns=[];tt=0
        for x in t:
            tt+=x.time
            if x.type=='note_on' and x.velocity>0: on[x.note]=tt
            elif x.type in('note_off','note_on') and x.note in on: s=on.pop(x.note);ns.append((s,x.note,tt-s))
        if ns and (best is None or np.mean([n for _,n,_ in ns])>np.mean([n for _,n,_ in best])): best=ns
    q=m.ticks_per_beat/8;grp={}
    for s,n,l in best: grp.setdefault(round(s/q),{})[n]=l/m.ticks_per_beat
    keys=sorted(grp);sh=r['sh'];g=cv2.imread(r['img'],0);st,view=page_views(g)
    for i,j in r['pairs']:
        c=r['res'][i];a={p+sh for p in c['pitches']};b=grp[keys[j]]
        if not a or not a<=set(b): continue
        L=[b[p] for p in a]
        if max(L)-min(L)>0.1: continue
        dur=float(np.mean(L))/0.95   # MIDI-Töne klingen 95 % des Notenwerts
        # Abstand zum nächsten Akkord derselben Zeile (in der MIDI) ist verlässlicher als die Tonlänge
        nx=[(i2,j2) for i2,j2 in r['pairs'] if i2==i+1]
        if nx and r['res'][i+1]['st']==c['st']:
            ioi=(keys[nx[0][1]]-keys[j])*q/m.ticks_per_beat
            if 0<ioi<=4 and dclass(ioi)!=dclass(dur): NOTE[(DUR[dclass(dur)],DUR[dclass(ioi)],k)]=NOTE.get((DUR[dclass(dur)],DUR[dclass(ioi)],k),0)+1

        nl,d,mid=view(c['st']);fill=float(np.mean([h['fill'] for h in c['heads']]))
        X.append(np.concatenate([strip(nl,c['x'],mid,d,fill),np.array(beam_feats(nl,c['x'],[h['y'] for h in c['heads']],d)+[dot_feat(nl,c['x'],[h['y'] for h in c['heads']],d)])*20]).astype(np.float32));Y.append(dclass(dur));K.append('RT' if k in mt else 'R');PC.append(k)
print(sorted(NOTE.items(),key=lambda x:-x[1])[:15])
np.savez_compressed('/tmp/claude-0/gsr/rhythm.npz',X=np.array(X),Y=np.array(Y),K=np.array(K),PC=np.array(PC))
print('gesamt',len(Y),np.bincount(Y,minlength=len(DUR)),{k:K.count(k) for k in set(K)})
