import pyvista as pv, numpy as np
g=pv.read("vor_Griffstock.stl"); w=pv.read("vor_Holz_R3.stl").triangulate()
# Punkte des Holzes, die im Griffstock liegen
sel=w.select_enclosed_points(g.triangulate().clean(),tolerance=0.0, check_surface=False)
ins=w.points[sel["SelectedPoints"].astype(bool)]
print(len(ins), ins[:10].round(1))
