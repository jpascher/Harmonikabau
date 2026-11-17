import pickle, numpy as np, pyvista as pv, math, warnings; warnings.filterwarnings("ignore")
D=pickle.load(open("netze.pkl","rb")); T=D["teile"]; DP=D["drehpunkte"]
def poly(v,t): return pv.PolyData(v.astype(float), np.hstack([np.full((len(t),1),3),t]).astype(np.int64))
gk=poly(*T["Griffstock_Koerper"][1:]).clean().triangulate()
res={}
for g in sorted(DP):
    kn=T[g+"/Knopf"]; v=kn[1].astype(float)
    x,Yp,Zp=DP[g]; c=v.mean(0); th=3.5/math.hypot(c[1]-Yp,c[2]-Zp)
    worst=0
    for f in (0.0,0.5,1.0):
        t=th*f; C,S=math.cos(t),math.sin(t)
        dY=v[:,1]-Yp; dZ=v[:,2]-Zp
        w=v.copy(); w[:,1]=Yp+C*dY+S*dZ; w[:,2]=Zp-S*dY+C*dZ
        m=poly(w,kn[2]); 
        # dichter abtasten
        m=m.subdivide(2)
        sel=m.select_enclosed_points(gk,tolerance=0.0,check_surface=False)
        n=int(sel["SelectedPoints"].sum())
        if n: worst=max(worst,n); res.setdefault(g,[]).append((f,n))
    # minimaler Abstand gedrueckt
print(len(res),"Knoepfe mit Ueberschneidung"); 
for g in list(res)[:20]: print(g,res[g])
