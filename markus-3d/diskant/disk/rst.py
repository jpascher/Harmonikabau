import pyvista as pv
pv.global_theme.background="white"
p=pv.Plotter(off_screen=True,window_size=(1800,900),shape=(2,1))
for i,(y0,y1) in enumerate([(150,190),(45,85)]):
    p.subplot(i,0)
    for n,c in (("Holz","burlywood"),("Alu","silver")):
        m=pv.read(f"vor_{n}.stl").clip_box([-5,385,y0,y1,5,75],invert=False)
        p.add_mesh(m,color=c)
    p.view_xz(); p.camera.zoom(1.6)
p.screenshot("stimmstock_seite.png")
