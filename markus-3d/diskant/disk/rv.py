import pyvista as pv, glob
pv.global_theme.background="white"
cols={"Holz":"burlywood","Verdeck":"#222222","Griffstock":"saddlebrown","Knopf":"ivory","Klappe":"#333333","Leder":"peru","Kunststoff":"lightgray"}
p=pv.Plotter(off_screen=True,window_size=(1800,1000),shape=(1,2))
for i,cam in enumerate([[(-200,-260,-260),(190,100,10),(0,0,1)],[(-150,-200,250),(190,100,40),(0,0,1)]]):
    p.subplot(0,i)
    for f in glob.glob("vor_*.stl"):
        n=f[4:-4]; c=cols.get(n, cols.get(n.split("_")[0],"silver"))
        p.add_mesh(pv.read(f),color=c)
    p.camera_position=cam
p.screenshot("voll_v14.png")
