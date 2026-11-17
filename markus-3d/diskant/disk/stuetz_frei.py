import cadquery as cq, glob, numpy as np, math
from shapely.geometry import Polygon, box
from shapely.ops import unary_union
def proj(s, zmin=-28.0, zmax=0.5):
    # Grundriss (X,Y) der Teile im Verdeckraum
    v,t=s.tessellate(0.1,0.3); V=np.array([(p.x,p.y,p.z) for p in v]); out=[]
    for a,b,c in t:
        P=V[[a,b,c]]
        if P[:,2].max()<zmin or P[:,2].min()>zmax: continue
        pg=Polygon(P[:,:2])
        if pg.area>1e-6: out.append(pg.buffer(0.01))
    return unary_union(out)
obst=[]
for f in glob.glob("jp_brep/Hebel_*__*.brep"):
    if "Knopf" in f or "Holzteil" in f: continue
    s=cq.Shape.importBrep(f); b=s.BoundingBox()
    if b.zmin>0.5: continue
    obst.append(proj(s))
O=unary_union(obst)
import pickle; pickle.dump(O,open("stuetz_obst.pkl","wb"))
deck=cq.Shape.importBrep("jp_brep/Verdeck_Deckplatte.brep")
sec=deck.intersect(cq.Solid.makeBox(400,300,0.4,cq.Vector(-10,-50,-31.0)))
# Material der Deckplatte im Grundriss
D=proj(sec,-40,10); pickle.dump(D,open("stuetz_deck.pkl","wb"))
print("Hindernisse", O.bounds, "Deck", D.bounds)
