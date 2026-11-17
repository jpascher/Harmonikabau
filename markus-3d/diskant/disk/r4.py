import pyvista as pv
pv.global_theme.background="white"
p=pv.Plotter(off_screen=True,window_size=(1600,1000),shape=(1,2))
p.subplot(0,0)
p.add_mesh(pv.read("vor_Kunststoff.stl").clip("x",origin=(200,0,0),invert=True),color="royalblue")
p.camera_position=[(100,90,30),(120,20,-18),(0,0,-1)]; p.add_text("Kunststoffteil",font_size=12,color="k")
p.subplot(0,1)
p.add_mesh(pv.read("vor_Griffstock.stl"),color="saddlebrown")
p.add_mesh(pv.read("vor_Kunststoff.stl"),color="royalblue")
p.add_mesh(pv.read("vor_Alu.stl"),color="silver",opacity=0.5)
p.camera_position=[(-60,160,-160),(190,10,-50),(0,0,-1)]; p.add_text("Griffstock",font_size=12,color="k")
p.screenshot("v4_3d.png")
