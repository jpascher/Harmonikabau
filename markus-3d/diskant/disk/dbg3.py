import cadquery as cq, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
L=lambda n: cq.Shape.importBrep(f"dbg_{n}.brep")
parts=[("Verdeck_Deckplatte","k"),("Verdeck_Seite_links","gray"),("Griffstock_Koerper","peru"),("Griffstock_Abdeckung","brown"),("Rahmen_vorne","g"),("Rahmen_hinten","b")]
fig,axs=plt.subplots(1,3,figsize=(18,7))
for ax,x in zip(axs,(12.5,30,200)):
  for n,c in parts:
    s=L(n)
    try: sec=s.intersect(cq.Face.makePlane(600,600,cq.Vector(x,0,0),cq.Vector(1,0,0)))
    except Exception: continue
    for e in sec.Edges():
      pts=[e.positionAt(t/20) for t in range(21)]
      ax.plot([p.y for p in pts],[p.z for p in pts],c)
  ax.set_xlim(-20,60); ax.set_ylim(-60,20); ax.set_aspect('equal'); ax.grid(); ax.set_title(f"x={x}")
plt.savefig("dbg3.png",dpi=80)
