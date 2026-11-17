import pyvista as pv
pv.global_theme.background="white"
p=pv.Plotter(off_screen=True,window_size=(1600,1000),shape=(1,2))
m=pv.read("vor_Holz.stl")
p.subplot(0,0); p.add_mesh(m,color="burlywood"); p.camera_position=[(-150,-150,180),(60,40,60),(0,0,1)]
p.subplot(0,1); p.add_mesh(m.clip_box([-5,20,-5,20,50,90],invert=False),color="burlywood"); p.camera_position=[(60,60,110),(5,5,70),(0,0,1)]
p.screenshot("rahmen_leiste.png")
