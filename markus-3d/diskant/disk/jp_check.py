import cadquery as cq, numpy as np, pyvista as pv, math, glob, os, re, json, warnings; warnings.filterwarnings("ignore")
B="jp_brep/"
def mesh(f, tol=0.05):
    s=cq.Shape.importBrep(B+f); v,t=s.tessellate(tol,0.2)
    v=np.array([(p.x,p.y,p.z) for p in v]); return v, np.array(t)
def poly(v,t): return pv.PolyData(v.astype(float), np.hstack([np.full((len(t),1),3),t]).astype(np.int64)).clean().triangulate()
koerper={n:poly(*mesh(n+".brep")) for n in ("Griffstock_Koerper","Griffstock_Teil2","Griffstock_Abdeckung","Kunststoffteil_Lagerblock","Dreieckleiste")}
DP=json.load(open("drehpunkte_tmp.json")) if os.path.exists("drehpunkte_tmp.json") else None
gruppen=sorted(set(f.split("__")[0] for f in os.listdir(B) if f.startswith("Hebel_")))
erg=[]
for g in gruppen:
    ax=cq.Shape.importBrep(B+g+"__Achse.brep").BoundingBox().center; Yp,Zp=ax.y,ax.z
    kn=mesh(g+"__Knopf.brep")[0].mean(0); th=3.5/math.hypot(kn[1]-Yp,kn[2]-Zp)
    for teil in ("Hebel_Alu","Holzteil","Knopf"):
        if not os.path.exists(B+f"{g}__{teil}.brep"): continue
        v,t=mesh(f"{g}__{teil}.brep",0.03)
        for f_ in (0.0,1.0):
            T_=th*f_; C,S=math.cos(T_),math.sin(T_); dY=v[:,1]-Yp; dZ=v[:,2]-Zp
            w=v.copy(); w[:,1]=Yp+C*dY+S*dZ; w[:,2]=Zp-S*dY+C*dZ
            m=poly(w,t).subdivide(1)
            for kn_,K in koerper.items():
                if teil=="Hebel_Alu" and kn_ in ("Kunststoffteil_Lagerblock",): continue   # Achse/Lager im Kunststoff
                b=K.bounds; P=m.points
                inb=(P[:,0]>b[0])&(P[:,0]<b[1])&(P[:,1]>b[2])&(P[:,1]<b[3])&(P[:,2]>b[4])&(P[:,2]<b[5])
                if not inb.any(): continue
                sub=pv.PolyData(P[inb]); sel=sub.select_enclosed_points(K,tolerance=0.0,check_surface=False)
                n=int(sel["SelectedPoints"].sum())
                if n: 
                    Q=sub.points[sel["SelectedPoints"].astype(bool)]
                    erg.append((g,teil,"gedrueckt" if f_ else "ruhe",kn_,n,np.round(Q.min(0),1).tolist(),np.round(Q.max(0),1).tolist()))
for e in erg: print(e)
print("Treffer:",len(erg))
