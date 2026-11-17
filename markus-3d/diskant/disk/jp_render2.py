import cadquery as cq, pyvista as pv, numpy as np
def mesh(f):
    s=cq.Shape.importBrep(f); v,t=s.tessellate(0.05,0.2); v=np.array([(p.x,p.y,p.z) for p in v]); return pv.PolyData(v,np.hstack([[3,*q] for q in t]))
for nm,f,c,vv in [("t2a","Griffstock_Teil2","#7a9a3a",(0.3,0.6,0.5)),("t2b","Griffstock_Teil2","#7a9a3a",(0.3,-0.4,0.8)),("gka","Griffstock_Koerper","#c8955a",(0.3,0.6,0.5))]:
    pl=pv.Plotter(off_screen=True,window_size=(1600,900)); pl.set_background("white")
    m=mesh(f"jp_brep/{f}.brep").clip_box((130,250,-100,300,-300,100),invert=False)
    pl.add_mesh(m,color=c,smooth_shading=True)
    pl.view_vector(vv,(0,0,-1)); pl.reset_camera(); pl.camera.zoom(1.3); pl.screenshot(f"jp_{nm}.png")
