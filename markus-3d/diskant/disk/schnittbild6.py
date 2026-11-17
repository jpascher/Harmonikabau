import pyvista as pv, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt, glob
farb={"Griffstock":"saddlebrown","Kunststoff":"blue","Alu":"gray","Hebelholz":"orange","Knopf":"m","Holz":"tan","Klappe":"k","Leder":"olive","Stahl":"r","Schraube":"gold"}
fig,axs=plt.subplots(4,1,figsize=(16,22))
for ax,x in zip(axs,(165,155,160,150)):
    for f in glob.glob("vor_*.stl"):
        n=f[4:-4]; m=pv.read(f); s=m.slice(normal="x",origin=(x,0,0))
        if s.n_points==0: continue
        for i in range(s.n_cells):
            p=s.get_cell(i).points; ax.plot(p[:,1],-p[:,2],color=farb.get(n,"c"),lw=0.8)
    ax.set_aspect("equal"); ax.set_title(f"X={x}"); ax.set_xlim(-20,190); ax.set_ylim(-3,45); ax.invert_yaxis(); ax.grid(lw=0.3)
fig.savefig("schnitt_v6.png",dpi=60)
