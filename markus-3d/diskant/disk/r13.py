import pyvista as pv
pv.global_theme.background="white"
p=pv.Plotter(off_screen=True,window_size=(1800,900),shape=(2,1))
cols={"Griffstock":"burlywood","Knopf":"ivory","Holz_R1":"red","Holz_R2":"blue","Holz_R3":"green","Holz_R4":"orange"}
for i,clip in enumerate([None,(60,170)]):
    p.subplot(i,0)
    for k,c in cols.items():
        m=pv.read(f"vor_{k}.stl")
        if clip: m=m.clip_box([clip[0],clip[1],-100,200,-200,50],invert=False)
        p.add_mesh(m,color=c, smooth_shading=False)
    cx=190 if not clip else 115
    p.camera_position=[(cx,15+260*0.98,-80-260*0.19+(-60 if clip else 0)),(cx,15,-75),(0,-0.19,-0.98)]
    p.enable_parallel_projection() if not clip else None
p.screenshot("v13_stufen.png")
