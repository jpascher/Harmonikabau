import pickle, numpy as np, pyvista as pv, math, warnings; warnings.filterwarnings("ignore")
D=pickle.load(open("netze.pkl","rb")); T=D["teile"]; DP=D["drehpunkte"]
def poly(v,t): return pv.PolyData(v.astype(float), np.hstack([np.full((len(t),1),3),t]).astype(np.int64))
gk=poly(*T["Griffstock_Koerper"][1:]).clean().triangulate()
KD=(0.1908,0.9816); KE=(KD[1],-KD[0])
for g in ("Hebel_R1_05","Hebel_R2_20","Hebel_R3_30","Hebel_R4_45"):
    if g+"/Knopf" not in T: continue
    kn=T[g+"/Knopf"]; v=kn[1].astype(float)
    x,Yp,Zp=DP[g]; c=v.mean(0); th=3.5/math.hypot(c[1]-Yp,c[2]-Zp)
    for f in (0.0,1.0):
        t=th*f; C,S=math.cos(t),math.sin(t); dY=v[:,1]-Yp; dZ=v[:,2]-Zp
        w=v.copy(); w[:,1]=Yp+C*dY+S*dZ; w[:,2]=Zp-S*dY+C*dZ
        m=poly(w,kn[2]).subdivide(2)
        sel=m.select_enclosed_points(gk,tolerance=0.0,check_surface=False)
        P=m.points[sel["SelectedPoints"].astype(bool)]
        u=-P[:,2]; Y=P[:,1]; a=u*KD[0]+Y*KD[1]; b=u*KE[0]+Y*KE[1]
        uc=-w[:,2]; Yc=w[:,1]; ac=uc*KD[0]+Yc*KD[1]; bc=uc*KE[0]+Yc*KE[1]
        print(g,f,"Knopf a %.1f..%.1f b %.1f..%.1f"%(ac.min(),ac.max(),bc.min(),bc.max()), "Ueberschn.:", len(P), ("a %.1f..%.1f b %.1f..%.1f X %.1f..%.1f"%(a.min(),a.max(),b.min(),b.max(),P[:,0].min(),P[:,0].max())) if len(P) else "")
