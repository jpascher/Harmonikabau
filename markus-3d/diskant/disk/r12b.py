import pyvista as pv, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
farb={"Griffstock":"saddlebrown","Alu_R1":"red","Alu_R2":"blue","Alu_R3":"green","Alu_R4":"orange","Holz_R1":"pink","Holz_R2":"lightblue","Holz_R3":"lightgreen","Holz_R4":"gold","Kunststoff":"gray","Knopf":"m","Klappe":"k"}
ms={n:pv.read(f"vor_{n}.stl") for n in farb}
fig,axs=plt.subplots(1,4,figsize=(20,14))
for ax,x in zip(axs,(125,135,140,130)):
    for n,m in ms.items():
        s=m.slice(normal="x",origin=(x,0,0))
        for i in range(s.n_cells):
            p=s.get_cell(i).points; ax.plot(p[:,1],-p[:,2],color=farb[n],lw=1)
    ax.set_aspect("equal"); ax.set_title(f"X={x}"); ax.set_xlim(-20,45); ax.set_ylim(-3,122); ax.invert_yaxis(); ax.grid(lw=.3)
fig.savefig("v12_schnitt.png",dpi=60)
