import cadquery as cq, glob, numpy as np, math, pickle
from shapely.geometry import Polygon, box, Point
from shapely.ops import unary_union
exec(open("stuetz_frei.py").read().split("obst=[]")[0])
L=lambda f: cq.Shape.importBrep(f)
obst=[pickle.load(open("stuetz_obst.pkl","rb"))]
for g in sorted(set(f.split("__")[0] for f in glob.glob("jp_brep/Hebel_*__Achse.brep"))):
    ax=L(g+"__Achse.brep").BoundingBox().center; kn=L(g+"__Knopf.brep").BoundingBox().center
    th=math.degrees(3.5/math.hypot(kn.y-ax.y,kn.z-ax.z))
    for t in ("Hebel_Alu","Klappe","Klappenpolster"):
        s=L(f"{g}__{t}.brep").rotate(cq.Vector(ax.x,ax.y,ax.z),cq.Vector(ax.x+1,ax.y,ax.z),-th)
        obst.append(proj(s))
O=unary_union(obst); D=pickle.load(open("stuetz_deck.pkl","rb"))
pickle.dump(O,open("stuetz_obst2.pkl","wb"))
res=[]
for x in np.arange(12,368,0.25):
    for y0 in np.arange(36,150,0.5):
        if O.intersects(box(x-2.5,y0,x+2.5,y0+40)): continue
        m=Point(x,y0+20)
        if D.intersects(m.buffer(1.6)): continue      # Loch im Verdeck fuer die Raendelschraube
        res.append((round(float(x),2),round(float(y0),1)))
print(len(res)); 
import collections
by=collections.defaultdict(list)
for x,y in res: by[x].append(y)
for x in sorted(by): print(x, min(by[x]), max(by[x]), len(by[x]))
