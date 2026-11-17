import pyvista as pv, itertools
n=["Alu_R1","Alu_R2","Alu_R3","Alu_R4","Holz_R1","Holz_R2","Holz_R3","Holz_R4","Klappe","Leder","Griffstock","Kunststoff","Holz"]
m={k:pv.read(f"vor_{k}.stl").triangulate() for k in n}
for a,b in itertools.combinations(n,2):
    if a[:4]==b[:4]=="Holz" and a!="Holz" and b!="Holz": pass
    if {a,b} <= {"Klappe","Leder","Griffstock","Kunststoff","Holz"}: continue
    if a[-2:]==b[-2:] and a[-3]=="R": continue  # gleiche Reihe Alu/Holz
    c,nc=m[a].collision(m[b],contact_mode=0)
    if nc: 
        pts=c.extract_cells(c.field_data["ContactCells"]).points if "ContactCells" in c.field_data else None
        bb=[round(v) for v in c.extract_cells(c.field_data["ContactCells"]).bounds] if pts is not None else ""
        print(a,b,nc,bb)
