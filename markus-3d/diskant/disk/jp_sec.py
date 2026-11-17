import cadquery as cq, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, schnitt_trafo as ST, glob
def z(Y,Zv):
    u,s=-Zv,ST.Y0-Y; return (ST.P_B[0]+s*ST.D[0]+u*ST.N[0], ST.P_B[1]+s*ST.D[1]+u*ST.N[1])
xs=[190.0,195.0,200.0,205.0,197.5]
K=[("Griffstock_Koerper","peru"),("Griffstock_Teil2","green"),("Griffstock_Abdeckung","brown")]
fig,axs=plt.subplots(1,len(xs),figsize=(6*len(xs),9))
for ax,x in zip(axs,xs):
    for n,c in K:
        s=cq.Shape.importBrep(f"jp_brep/{n}.brep")
        sec=s.intersect(cq.Face.makePlane(800,800,cq.Vector(x,100,0),cq.Vector(1,0,0)))
        for e in sec.Edges():
            p=[z(q.y,q.z) for q in (e.positionAt(t/10) for t in range(11))]; ax.plot(*zip(*p),c)
    for f in glob.glob("jp_brep/Hebel_*__Holzteil.brep")+glob.glob("jp_brep/Hebel_*__Hebel_Alu.brep"):
        s=cq.Shape.importBrep(f); b=s.BoundingBox()
        if not (b.xmin<x<b.xmax): continue
        sec=s.intersect(cq.Face.makePlane(800,800,cq.Vector(x,100,0),cq.Vector(1,0,0)))
        for e in sec.Edges():
            p=[z(q.y,q.z) for q in (e.positionAt(t/10) for t in range(11))]; ax.plot(*zip(*p),'b' if 'Alu' in f else 'saddlebrown',lw=0.8)
    ax.set_xlim(-150,-95); ax.set_ylim(20,145); ax.set_aspect('equal'); ax.grid(alpha=.3); ax.set_title(f"X={x}")
fig.savefig("jp_sec.png",dpi=55)
