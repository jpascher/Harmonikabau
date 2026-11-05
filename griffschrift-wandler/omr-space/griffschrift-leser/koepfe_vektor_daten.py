# Aufruf: python3 koepfe_vektor_daten.py daten.npz a.pdf b.pdf …  (PDFs aus MC-Score oder Finale; Seitenbilder in /tmp/claude-0/h3)
# Lerndaten für Köpfe aus Vektor-PDFs: Bild wie im Space (200 dpi, gs2xml.load_gray), Köpfe genau aus der Textebene,
# Merkmale wie detect.py (Ausschnitt mit und ohne Linien/Hälse an allen Kandidaten), dazu verschlechterte Fassungen (Scan-ähnlich)
import sys,os,json,glob,numpy as np,cv2,fitz
GS=os.path.join(os.path.dirname(os.path.abspath(__file__)),'..','gsleser');sys.path.insert(0,GS)
import gs2xml as X, gsread as G, vektor as VK, detect as D
CL=['none','o_f','o_h','x_f','x_h']
def degrade(g,k,rng):
    if k==0: return g
    g=g.astype(np.float32)
    if k>=1: g=cv2.GaussianBlur(g,(0,0),rng.uniform(0.6,1.4))
    if k>=2: g=g+rng.normal(0,rng.uniform(5,14),g.shape)
    if k>=3:
        ker=np.ones((2,2),np.uint8); g=cv2.erode(g,ker) if rng.random()<0.5 else cv2.dilate(g,ker)
    g=np.clip(g*rng.uniform(0.9,1.0)+rng.uniform(0,25),0,255).astype(np.uint8)
    if k>=2:
        ok,b=cv2.imencode('.jpg',g,[cv2.IMWRITE_JPEG_QUALITY,int(rng.uniform(35,80))]);g=cv2.imdecode(b,0)
    return g
def page_ds(png,heads,rng,variants=(0,1,2,3)):
    from PIL import Image
    g0=X.load_gray(png);W,H=Image.open(png).size;fx=g0.shape[1]/W
    Xs=[];Ys=[]
    for k in variants:
        g=degrade(g0,k,rng);_,bw=cv2.threshold(g,200,255,cv2.THRESH_BINARY_INV);_,bwl=cv2.threshold(g0,240,255,cv2.THRESH_BINARY_INV)
        st=G.staves(bwl)
        hs=VK.det_from_heads(heads,fx,st)
        for si,ys in enumerate(st):
            d=float(np.mean(np.diff(ys)));mid=ys[2];nol=G.stems_off(G.remove_lines(bw,ys,d),d)
            ink=cv2.boxFilter((nol>0).astype(np.float32),-1,(max(3,int(d*0.9)),max(3,int(d*0.6))))
            r=int(d*0.9);bwp=cv2.copyMakeBorder(bw,r*2,r*2,r*2,r*2,cv2.BORDER_CONSTANT,value=0);nolp=cv2.copyMakeBorder(nol,r*2,r*2,r*2,r*2,cv2.BORDER_CONSTANT,value=0)
            def fp(pad,x,y):
                x0,y0=int(round(x))-r,int(round(y))-r
                return cv2.resize(pad[y0+2*r:y0+4*r,x0+2*r:x0+4*r],(20,20),interpolation=cv2.INTER_AREA).ravel()
            my=[h for h in hs if h['st']==si]
            for p in range(-12,13):
                y=int(round(mid+0.5-p*d/2))
                if not(0<=y<bw.shape[0]): continue
                at=[h for h in my if h['p']==p]
                for x in range(int(d*2),bw.shape[1]-int(d),2):
                    if ink[y,x]<D.INK: continue
                    lab=0
                    for h in at:
                        if abs(h['x']-x)<=0.3*d: lab=CL.index(h['t']+('_f' if h['fill'] else '_h'));break
                        if abs(h['x']-x)<=0.6*d: lab=-1
                    if lab<0: continue
                    if lab==0 and rng.random()>0.08: continue   # Negative ausdünnen
                    Xs.append(np.concatenate([fp(bwp,x,y),fp(nolp,x,y)]));Ys.append(lab)
    return np.array(Xs,np.uint8),np.array(Ys)
if __name__=='__main__':
    out=sys.argv[1];pdfs=sys.argv[2:];rng=np.random.default_rng(1);AX=[];AY=[]
    for pdf in pdfs:
        doc=fitz.open(pdf)
        for i,pg in enumerate(doc):
            png=f'/tmp/claude-0/h3/{os.path.basename(pdf)[:8]}_{i}.png';pg.get_pixmap(dpi=200).save(png)
            hs=VK.heads_pdf_page(pg,200)
            if len(hs)<10: continue
            x,y=page_ds(png,hs,rng);AX.append(x);AY.append(y);print(png,x.shape,np.bincount(y,minlength=5),flush=True)
    np.savez_compressed(out,X=np.concatenate(AX),Y=np.concatenate(AY))
