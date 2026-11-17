import ezdxf, ezdxf.path, math
from ezdxf import bbox
def kontur_pts(e):
    pts=[(v.x,v.y) for v in ezdxf.path.make_path(e).flattening(0.05)]
    out=[]
    for q in pts:
        if not out or math.hypot(q[0]-out[-1][0],q[1]-out[-1][1])>1e-3: out.append(q)
    if len(out)>2 and math.hypot(out[0][0]-out[-1][0],out[0][1]-out[-1][1])<1e-3: out.pop()
    return out
HEBEL_DXF={1:("h1.dxf","Hebel-alu-1r"),2:("h2.dxf","Hebel-alu-2r"),3:("h3.dxf","Hebel-alu-3r"),4:("h4.dxf","Hebel-4r")}
def fl(q): return abs(sum(q[i][0]*q[(i+1)%len(q)][1]-q[(i+1)%len(q)][0]*q[i][1] for i in range(len(q))))/2
def hebel_profil(r):
    f,lay=HEBEL_DXF[r]; m_=ezdxf.readfile(f).modelspace()
    alu=max((e for e in m_ if e.dxftype()=="LWPOLYLINE" and e.dxf.layer==lay and len(e)>6),key=lambda e: bbox.extents([e]).size.y)
    offen=[kontur_pts(e) for e in m_ if e.dxftype()=="LWPOLYLINE" and e is not alu and not e.closed]
    offen=[o for o in offen if len(o)>=3]
    ring=None
    if len(offen)>=2:
        cand=[offen[0]+offen[1],offen[0]+list(reversed(offen[1]))]
        ring=max(cand,key=fl)
    return kontur_pts(alu),ring
def dreh(pts,th,c=(0,0)):
    c_,s_=math.cos(th),math.sin(th)
    return [((x-c[0])*c_-(y-c[1])*s_+c[0],(x-c[0])*s_+(y-c[1])*c_+c[1]) for x,y in pts]
# gestufter Griffstock-Schnitt (Mechanik-DXF, Schnitt Aussen), u = nach aussen, Y
GS=[(0.0,45.67),(74.01,32.62),(73.31,28.68),(89.07,25.91),(88.38,21.97),(104.13,19.19),(103.44,15.25),(131.01,10.39),(126.15,-17.19),(0.0,5.06)]
STUFE={4:((0.0,45.67),(74.01,32.62)),3:((73.31,28.68),(89.07,25.91)),2:((88.38,21.97),(104.13,19.19)),1:((103.44,15.25),(131.01,10.39))}
HUB={1:3.5}     # Annahme: aeusserste Reihe traegt die Taste -> 3,5 mm unter der Stufe
L0=(29.2,20.2); NEIG=math.radians(15.0); SD=(math.cos(NEIG),-math.sin(NEIG))
def oberkante(ring):
    # laengste nahezu waagrechte Kante des Holzteils mit groesstem y
    best=None
    for i in range(len(ring)):
        a,b=ring[i],ring[(i+1)%len(ring)]
        L=math.hypot(b[0]-a[0],b[1]-a[1])
        if L<8: continue
        ang=math.atan2(b[1]-a[1],b[0]-a[0])
        if abs(math.sin(ang))>0.3: continue
        sc=(a[1]+b[1])/2+L*0.01
        if best is None or sc>best[0]: best=(sc,a,b)
    return best[1],best[2]
def loese(r):
    alu,ring=hebel_profil(r)
    a,b=oberkante(ring)
    ka=math.atan2(b[1]-a[1],b[0]-a[0])
    if ka>math.pi/2: ka-=math.pi
    if ka<-math.pi/2: ka+=math.pi
    (pa,pb)=STUFE[r]; ks=math.atan2(pb[1]-pa[1],pb[0]-pa[0])
    th=ks-ka
    ring_r=dreh(ring,th); 
    # Hoehe ueber der Stufenlinie (senkrecht, positiv = ueber der Oberflaeche)
    n=(-math.sin(ks),math.cos(ks))
    def h(p,s):
        u=L0[0]+s*SD[0]+p[0]; Y=L0[1]+s*SD[1]+p[1]
        return (u-pa[0])*n[0]+(Y-pa[1])*n[1]
    # h linear in s: h(p,s)=h(p,0)+s*(SD.n)
    k=SD[0]*n[0]+SD[1]*n[1]
    hmax0=max(h(p,0) for p in ring_r)
    s=(-HUB.get(r,0.0)-hmax0)/k
    return dict(th=th,s=s,piv=(L0[0]+s*SD[0],L0[1]+s*SD[1]),alu=alu,ring=ring)
if __name__=="__main__":
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    fig,ax=plt.subplots(figsize=(11,9))
    ax.plot([p[0] for p in GS+[GS[0]]],[p[1] for p in GS+[GS[0]]],'k-',lw=1.5)
    cols={1:'r',2:'g',3:'b',4:'m'}
    for r in (1,2,3,4):
        d=loese(r); pv=d['piv']; th=d['th']
        print(r,"Drehung %.1f  s %.1f  Lager u %.1f Y %.1f"%(math.degrees(th),d['s'],pv[0],pv[1]))
        for pts,lw in ((d['alu'],1.0),(d['ring'],2.0)):
            q=[(pv[0]+x,pv[1]+y) for x,y in dreh(pts,th)]; q.append(q[0])
            ax.plot([p[0] for p in q],[p[1] for p in q],color=cols[r],lw=lw)
        ax.plot(pv[0],pv[1],'o',color=cols[r])
    ax.axhline(0,color='gray',lw=0.5); ax.set_aspect('equal'); ax.grid(True,lw=0.3); ax.set_xlabel('u (nach aussen ab Boden)'); ax.set_ylabel('Y')
    ax.set_xlim(-20,150); ax.set_ylim(-30,220)
    fig.savefig("hebel_schnitt.png",dpi=70)
