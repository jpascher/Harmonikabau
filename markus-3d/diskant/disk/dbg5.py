import cadquery as cq, pyvista as pv, numpy as np, glob
def mesh(s):
    v,f=s.tessellate(0.05,0.2); v=np.array([(p.x,p.y,p.z) for p in v]); f=np.hstack([[3,*t] for t in f])
    return pv.PolyData(v,f)
col={"Verdeck":"#333333","Griffstock":"#c8955a","Rahmen":"#b07a40","Halbrund":"#555555"}
pl=pv.Plotter(off_screen=True,window_size=(1600,1000)); pl.set_background("white")
for fn in [f for f in glob.glob("dbg_*.brep") if "Griffstock" not in f]:
    n=fn[4:-5]; c=[v for k,v in col.items() if k in n]; c=c[-1] if c else "gray"
    pl.add_mesh(mesh(cq.Shape.importBrep(fn)),color=c,smooth_shading=True)

for nm,vv,zoom in [("d",(-0.4,-0.5,-0.75),1.2),("e",(0.6,-0.4,-0.7),4.0)]:
    pl.view_vector(vv,(0,0,-1)); pl.reset_camera(); pl.camera.zoom(zoom)
    if nm=="e": pl.camera.focal_point=(340,32,-32); pl.camera.position=(340+70,32-50,-32-90)
    pl.screenshot(f"dbg4{nm}.png")
