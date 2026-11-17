import cadquery as cq, pyvista as pv, numpy as np, glob
def mesh(f):
    s=cq.Shape.importBrep(f); v,t=s.tessellate(0.05,0.2); v=np.array([(p.x,p.y,p.z) for p in v]); return pv.PolyData(v,np.hstack([[3,*q] for q in t]))
pl=pv.Plotter(off_screen=True,window_size=(1600,1000)); pl.set_background("white")
pl.add_mesh(mesh("jp_brep/Griffstock_Koerper.brep"),color="#c8955a",smooth_shading=True)
pl.add_mesh(mesh("jp_brep/Griffstock_Teil2.brep"),color="#7a9a3a",smooth_shading=True)
for f in glob.glob("jp_brep/Hebel_R*_2[0-3]__*.brep")+glob.glob("jp_brep/Hebel_R*_3[3-5]__*.brep"):
    if "Klappe" in f: continue
    c={'1':'red','2':'blue','3':'green','4':'orange'}[f.split("Hebel_R")[1][0]]
    if "Holz" in f: c="#8b5a2b"
    if "Knopf" in f: c="#eeeedd"
    pl.add_mesh(mesh(f),color=c)
pl.add_mesh(pv.Plane(center=(150,0,0)),opacity=0)
pl.camera_position=[(80,-160,-160),(160,20,-50),(0,0,-1)]
pl.screenshot("jp_r1.png")
# Schnitt bei X=190 (Hebel R4_45 x?) : clip
pl2=pv.Plotter(off_screen=True,window_size=(1400,1000)); pl2.set_background("white")
for f,c in [("jp_brep/Griffstock_Koerper.brep","#c8955a"),("jp_brep/Griffstock_Teil2.brep","#7a9a3a")]:
    pl2.add_mesh(mesh(f).clip(normal=(1,0,0),origin=(195,0,0)),color=c,smooth_shading=True)
for f in glob.glob("jp_brep/Hebel_R2_21__*.brep"):
    if "Klappe" in f: continue
    pl2.add_mesh(mesh(f),color="#8b5a2b" if "Holz" in f else "blue")
pl2.view_vector((1,0,0),(0,0,-1)); pl2.reset_camera(); pl2.camera.zoom(1.6)
pl2.screenshot("jp_r2.png")
