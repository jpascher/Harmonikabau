import pyvista as pv, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
g=pv.read("vor_Griffstock.stl")
fig,axs=plt.subplots(2,3,figsize=(21,13)); axs=axs.ravel()
farb={"Griffstock":"saddlebrown","Alu":"gray","Hebelholz":"orange","Kunststoff":"blue","Knopf":"m"}
ms={n:pv.read(f"vor_{n}.stl") for n in farb}
for ax,x in zip(axs[:3],(65,70,75)):
    for n,m in ms.items():
        s=m.slice(normal="x",origin=(x,0,0))
        for i in range(s.n_cells):
            p=s.get_cell(i).points; ax.plot(p[:,1],-p[:,2],color=farb[n],lw=0.8)
    ax.set_aspect("equal"); ax.set_title(f"X={x}"); ax.set_xlim(-20,45); ax.set_ylim(-3,122); ax.invert_yaxis()
for ax,u in zip(axs[3:],(40,70,95)):
    for n,m in ms.items():
        s=m.slice(normal="z",origin=(0,0,-u))
        for i in range(s.n_cells):
            p=s.get_cell(i).points; ax.plot(p[:,0],p[:,1],color=farb[n],lw=0.6)
    ax.set_aspect("equal"); ax.set_title(f"u={u} (Aufsicht-Schnitt)"); ax.set_xlim(40,130); ax.set_ylim(-20,40)
fig.savefig("v10_schnitt.png",dpi=55)
