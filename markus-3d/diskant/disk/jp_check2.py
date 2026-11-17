import cadquery as cq, math, os, sys, schnitt_trafo as ST
B="jp_brep/"; L=lambda f: cq.Shape.importBrep(B+f)
K={n:L(n+".brep") for n in ("Griffstock_Koerper","Griffstock_Teil2","Griffstock_Abdeckung","Dreieckleiste","Kunststoffteil_Lagerblock")}
def z(Y,Zv):
    u,s=-Zv,ST.Y0-Y; return (round(ST.P_B[0]+s*ST.D[0]+u*ST.N[0],1), round(ST.P_B[1]+s*ST.D[1]+u*ST.N[1],1))
gr=sys.argv[1:] or ["Hebel_R1_08","Hebel_R2_21","Hebel_R3_34","Hebel_R4_45","Hebel_R1_01"]
for g in gr:
    ax=L(g+"__Achse.brep").BoundingBox().center
    kn=L(g+"__Knopf.brep").BoundingBox().center; th=math.degrees(3.5/math.hypot(kn.y-ax.y,kn.z-ax.z))
    for teil in ("Hebel_Alu","Holzteil","Knopf"):
        if not os.path.exists(B+f"{g}__{teil}.brep"): continue
        s0=L(f"{g}__{teil}.brep")
        for f_ in (0,1):
            # gedrueckt: Y'=Yp+C dY+S dZ -> Drehung um X um -t (rechte Hand)
            s=s0.rotate(cq.Vector(ax.x,ax.y,ax.z),cq.Vector(ax.x+1,ax.y,ax.z),-th*f_) if f_ else s0
            for kn_,k in K.items():
                if teil=="Hebel_Alu" and kn_=="Kunststoffteil_Lagerblock": continue
                i=s.intersect(k); v=i.Volume()
                if v>0.01:
                    b=i.BoundingBox(); print(g,teil,"gedrueckt" if f_ else "ruhe",kn_,"V=%.2f"%v,"Zeichnung",z(b.ymin,b.zmin),z(b.ymax,b.zmax))
print("fertig")
