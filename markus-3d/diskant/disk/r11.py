import pyvista as pv
pv.global_theme.background="white"
p=pv.Plotter(off_screen=True,window_size=(1800,1000))
col={"Alu_R1":"red","Alu_R2":"blue","Alu_R3":"green","Alu_R4":"orange","Holz_R1":"pink","Holz_R2":"lightblue","Holz_R3":"lightgreen","Holz_R4":"khaki","Klappe":"saddlebrown","Leder":"peru","Kunststoff":"lightgray"}
for k,c in col.items():
    m=pv.read(f"vor_{k}.stl").clip("x",origin=(150,0,0),invert=True)
    p.add_mesh(m,color=c)
p.camera_position=[(60,250,-120),(110,80,-10),(0,0,-1)]
p.screenshot("v11_3d.png")
