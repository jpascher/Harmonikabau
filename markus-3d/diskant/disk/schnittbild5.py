import pyvista as pv, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, glob
farb={"Griffstock":"saddlebrown","Kunststoff":"blue","Alu":"gray","Hebelholz":"orange","Knopf":"m","Holz":"tan","Klappe":"k","Leder":"olive","Stahl":"r","Schraube":"gold"}
fig,axs=plt.subplots(1,3,figsize=(18,9))
for ax,x in zip(axs,(30,147,165)):
    for f in glob.glob("vor_*.stl"):
        n=f[4:-4]; m=pv.read(f); s=m.slice(normal="x",origin=(x,0,0))
        if s.n_points==0: continue
        for i in range(s.n_cells):
            p=s.get_cell(i).points; ax.plot(p[:,1],-p[:,2],color=farb.get(n,"c"),lw=0.8)
    ax.set_aspect("equal"); ax.set_title(f"X={x}"); ax.set_xlim(-5,45); ax.set_ylim(-3,40); ax.invert_yaxis(); ax.grid(lw=0.3)
fig.savefig("schnitt_v5.png",dpi=60)
