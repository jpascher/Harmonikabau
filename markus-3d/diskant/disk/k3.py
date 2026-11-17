import pyvista as pv, itertools, warnings; warnings.filterwarnings("ignore")
n=["Alu_R1","Alu_R2","Alu_R3","Alu_R4","Holz_R1","Holz_R2","Holz_R3","Holz_R4","Klappe","Griffstock","Kunststoff"]
m={k:pv.read(f"vor_{k}.stl").triangulate().clean() for k in n}
for a in n[:8]:
    for b in n:
        if a==b or a[-2:]==b[-2:]: continue
        sel=m[a].select_enclosed_points(m[b],tolerance=0.0,check_surface=False)
        ins=m[a].points[sel["SelectedPoints"].astype(bool)]
        if len(ins): print(a,"in",b,len(ins),ins.min(0).round(1),ins.max(0).round(1))
