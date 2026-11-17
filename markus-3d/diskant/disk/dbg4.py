import cadquery as cq, pyvista as pv, numpy as np, glob
def mesh(s):
    v,f=s.tessellate(0.05,0.2); v=np.array([(p.x,p.y,p.z) for p in v]); f=np.hstack([[3,*t] for t in f])
    return pv.PolyData(v,f)
col={"Verdeck":"#333333","Griffstock":"#c8955a","Rahmen":"#b07a40","Halbrund":"#555555"}
pl=pv.Plotter(off_screen=True,window_size=(1600,1000)); pl.set_background("white")
for fn in glob.glob("dbg_*.brep"):
    n=fn[4:-5]; c=[v for k,v in col.items() if k in n]; c=c[-1] if c else "gray"
    pl.add_mesh(mesh(cq.Shape.importBrep(fn)),color=c,smooth_shading=True)

for nm,vv,zoom in [("a",(-0.5,-0.6,-0.6),1.0),("b",(0.0,-0.5,-0.85),1.0),("c",(-0.3,-0.7,-0.5),3.0)]:
    pl.view_vector(vv,(0,0,-1)); pl.reset_camera(); pl.camera.zoom(zoom)
    if nm=="c": pl.camera.focal_point=(20,30,-30); pl.camera.position=(20-60,30-140,-30-100)
    pl.screenshot(f"dbg4{nm}.png")
