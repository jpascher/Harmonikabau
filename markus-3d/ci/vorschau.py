# Vorschau aus einem STEP-Modell: grobe GLB mit Farben (fuer den Browser / 3D-Viewer) und ein PNG-Standbild.
#   python3 vorschau.py modell.step vorschau.glb vorschau.png
import sys, os
os.environ.setdefault("PYVISTA_OFF_SCREEN", "true")
import cadquery as cq, trimesh, numpy as np, pyvista as pv

step, glb, png = sys.argv[1:4]
assy = cq.Assembly.importStep(step)
assy.export(glb, tolerance=0.8, angularTolerance=0.6)
print(f"GLB {os.path.getsize(glb) / 1e6:.1f} MB")
szene = trimesh.load(glb, force="scene")
gruppen = {}                                              # nach Farbe zusammenfassen: wenige grosse Netze rendern schnell
for m in szene.dump():
    if len(m.faces) == 0: continue
    farbe = (179, 179, 179)
    try:
        f_ = np.asarray(m.visual.material.baseColorFactor[:3], float)
        farbe = tuple(int(round(v)) for v in (f_ if f_.max() > 1 else f_ * 255))
    except Exception:
        try: farbe = tuple(int(v) for v in m.visual.main_color[:3])
        except Exception: pass
    gruppen.setdefault(farbe, []).append(trimesh.Trimesh(m.vertices, m.faces, process=False))
pl = pv.Plotter(off_screen=True, window_size=(1600, 1000))
pl.set_background("white")
for farbe, teile in gruppen.items():
    pl.add_mesh(pv.wrap(trimesh.util.concatenate(teile)), color=tuple(v / 255 for v in farbe))
pl.view_isometric(); pl.camera.zoom(1.3)
pl.screenshot(png)
print(f"PNG {png} ({len(gruppen)} Farben)")
