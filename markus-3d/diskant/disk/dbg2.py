import cadquery as cq
L=lambda n: cq.Shape.importBrep(f"dbg_{n}.brep")
d=L("Verdeck_Deckplatte"); g=L("Griffstock_Koerper"); a=L("Griffstock_Abdeckung")
for n,s in [("deck",d),("gs",g),("abd",a),("sl",L("Verdeck_Seite_links")),("front",L("Verdeck_Front"))]:
    b=s.BoundingBox(); print(n, round(b.xmin,1),round(b.xmax,1),round(b.ymin,1),round(b.ymax,1),round(b.zmin,1),round(b.zmax,1))
# Schnitt bei x=200, x=30: Y-Z Profil von deck/gs/abd
for x in (30,200):
  pl=cq.Face.makePlane(400,400,cq.Vector(x,100,0),cq.Vector(1,0,0))
  for n,s in [("deck",d),("gs",g),("abd",a)]:
    try:
      sec=s.intersect(cq.Solid.makeBox(0.2,300,200,cq.Vector(x-0.1,-50,-100)))
      b=sec.BoundingBox(); print(x,n,round(b.ymin,2),round(b.ymax,2),round(b.zmin,2),round(b.zmax,2))
    except Exception as e: print(x,n,e)
