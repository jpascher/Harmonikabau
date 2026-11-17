import cadquery as cq, glob, os, math
B="jp_brep/"; L=lambda f: cq.Shape.importBrep(f)
neu=[f for f in glob.glob(B+"*.brep") if any(k in f for k in ("Verdeckstuetze","Schraube_Stuetze","Raendelschraube"))]
rest=[f for f in glob.glob(B+"*.brep") if f not in neu]
print(len(neu), len(rest))
# gedrueckte Hebel dazu
ged=[]
for g in sorted(set(f.split("__")[0] for f in glob.glob(B+"Hebel_*__Achse.brep"))):
    ax=L(g+"__Achse.brep").BoundingBox().center; kn=L(g+"__Knopf.brep").BoundingBox().center
    th=math.degrees(3.5/math.hypot(kn.y-ax.y,kn.z-ax.z))
    for t in ("Hebel_Alu","Klappe","Klappenpolster"):
        ged.append((g+"__"+t+"(gedr)", L(f"{g}__{t}.brep").rotate(cq.Vector(ax.x,ax.y,ax.z),cq.Vector(ax.x+1,ax.y,ax.z),-th)))
andere=[(os.path.basename(f), L(f)) for f in rest]+ged
for f in neu:
    s=L(f); b=s.BoundingBox()
    for n,o in andere:
        ob=o.BoundingBox()
        if ob.xmin>b.xmax or ob.xmax<b.xmin or ob.ymin>b.ymax or ob.ymax<b.ymin or ob.zmin>b.zmax or ob.zmax<b.zmin: continue
        v=s.intersect(o).Volume()
        if v>0.01: print(os.path.basename(f), "x", n, "%.2f"%v)
    # Abstand zum naechsten Hebel/Klappe
print("fertig")
