# Aufruf: python3 koepfe_vektor_pruefen.py a.pdf …  (NOH3=1: ohne drittes Kopfmodell)
# Köpfe auf zurückgehaltenen Vektor-Seiten (Bild gelesen) gegen die Köpfe aus der Textebene
import sys,os,glob,numpy as np,cv2,fitz,collections,tempfile
sys.path.insert(0,os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','gsleser'))
import gs2xml as X, gsread as G, vektor as VK, detect as D
from PIL import Image
T=collections.Counter()
for pdf in sys.argv[1:]:
    doc=fitz.open(pdf)
    for i,pg in enumerate(doc):
        png=f'/tmp/claude-0/h3/{os.path.basename(pdf)[:8]}_{i}.png';pg.get_pixmap(dpi=200).save(png)
        hs=VK.heads_pdf_page(pg,200)
        g=X.load_gray(png);W,_=Image.open(png).size;fx=g.shape[1]/W
        st=G.staves(cv2.threshold(g,240,255,cv2.THRESH_BINARY_INV)[1]);gt=VK.det_from_heads(hs,fx,st)
        w=tempfile.mktemp(suffix='.png');cv2.imwrite(w,g);det=D.detect(w,1);os.remove(w)
        used=set()
        for h in gt:
            T['soll']+=1;T['soll_'+h['t']]+=1
            best=None
            for j,q in enumerate(det):
                if j in used or q['st']!=h['st'] or q['p']!=h['p']: continue
                dx=abs(q['x']-h['x'])
                if dx<0.7*h['d'] and (best is None or dx<best[0]): best=(dx,j)
            if best:
                used.add(best[1]);T['gef']+=1;T['gef_'+h['t']]+=1;T['typ']+=det[best[1]]['t']==h['t']
        T['det']+=len(det)
print(os.environ.get('NOH3','mitH3'),os.environ.get('GS_H3W',''),'gefunden %.3f | rund %.3f | Kreuz %.3f | Typ richtig %.3f | überzählig %d'%(T['gef']/T['soll'],T['gef_o']/max(1,T['soll_o']),T['gef_x']/max(1,T['soll_x']),T['typ']/max(1,T['gef']),T['det']-T['gef']),dict(T))
