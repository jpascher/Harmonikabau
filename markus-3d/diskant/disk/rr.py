import cadquery as cq, pyvista as pv
pv.global_theme.background="white"
p=pv.Plotter(off_screen=True,window_size=(1800,1000),shape=(3,1))
for i,f in enumerate(("Dr_rueck.step","Dr_front.step","Dr_seite.step")):
    s=cq.importers.importStep(f).val(); cq.exporters.export(s,"_t.stl",tolerance=0.05)
    p.subplot(i,0); p.add_mesh(pv.read("_t.stl"),color="burlywood"); p.add_text(f,font_size=10,color='k')
    p.view_xy(); p.camera.zoom(1.3)
p.screenshot("rahmen_d_check.png")
