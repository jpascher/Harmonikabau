# Ansicht zum Prüfen der Köpfe: je System zwei Hälften, erkannte Köpfe mit Nummer, x-Raster alle 50 px, Lagen links
import sys,json,cv2,numpy as np
sys.path.insert(0,'/home/claude/gsw/m')
import gs2xml as X,gsread as G,detect as D
key=sys.argv[1];src=sys.argv[2]
g=X.load_gray(src);fn=f'/home/claude/lab/{key}.png';cv2.imwrite(fn,g)
det=D.detect(fn,1);json.dump(det,open(f'/home/claude/lab/{key}_det.json','w'))
st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1])
c=cv2.cvtColor(g,cv2.COLOR_GRAY2BGR);col={'o':(0,160,0),'x':(0,0,255),'t':(200,0,200)}
for i,h in enumerate(det):
    cc=col.get(h['t'],(0,0,0));cv2.circle(c,(int(h['x']),int(h['y'])),9,cc,1 if not h['fill'] else 2)
    cv2.putText(c,str(i),(int(h['x'])-6,int(h['y'])-11),cv2.FONT_HERSHEY_SIMPLEX,0.33,cc,1)
for si,ys in enumerate(st):
    d=(ys[4]-ys[0])/4;y0=int(ys[0]-4*d);y1=int(ys[4]+3*d);cr=c[y0:y1].copy()
    for x in range(0,cr.shape[1],50): cv2.line(cr,(x,0),(x,6),(255,0,0),1);cv2.putText(cr,str(x),(x+1,14),cv2.FONT_HERSHEY_SIMPLEX,0.3,(255,0,0),1)
    for p in range(-8,9,2):
        y=int(ys[2]+0.5-p*d/2)-y0;cv2.putText(cr,str(p),(2,y+3),cv2.FONT_HERSHEY_SIMPLEX,0.3,(255,0,0),1)
    W=cr.shape[1];h=W//2
    cv2.imwrite(f'/home/claude/lab/{key}_{si}a.png',cv2.resize(cr[:,:h+30],None,fx=1.6,fy=1.6));cv2.imwrite(f'/home/claude/lab/{key}_{si}b.png',cv2.resize(np.hstack([cr[:,h-30:]]),None,fx=1.6,fy=1.6))
print(len(st),len(det))
