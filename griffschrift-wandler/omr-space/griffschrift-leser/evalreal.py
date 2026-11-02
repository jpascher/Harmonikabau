import json,sys
sys.path.insert(0,'/tmp/claude-0/gsr');import detect as D;from align import chords_det,chords_pin,align
P=json.load(open('/tmp/claude-0/gsr/pins.json'));keys=sys.argv[1:] or list(P);T=[0,0,0]
for k in keys:
    img=f'/tmp/claude-0/fuchs/{k}_u.png';det=D.detect(img,1);C=chords_det(det);A=[{(h['p'],h['t']) for h in c} for c in C];B=chords_pin(P[k]['n'])
    pr=align(A,B);ex=sum(A[i]==B[j] for i,j in pr);print(k,'Bild',len(A),'gepaart',len(pr),'gleich',ex);T[0]+=len(pr);T[1]+=ex;T[2]+=len(det)
print('gesamt gepaart',T[0],'gleich',T[1],round(T[1]/T[0],3),'Köpfe',T[2])
