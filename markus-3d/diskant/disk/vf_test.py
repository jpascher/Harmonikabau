import cadquery as cq
_vf=cq.Shape.importBrep("dbg_Verdeck_Front.brep")
_ff = max([f for f in _vf.Faces() if abs(abs(f.normalAt().y) - 1) < 1e-6], key=lambda f: f.Area())
r=[]
for w_ in _ff.innerWires():
    fa_=cq.Face.makeFromWires(w_); r.append((round(fa_.Area(),1), round(fa_.Center().x,1), round(fa_.Center().z,1)))
import collections; print(collections.Counter(a for a,_,_ in r).most_common(12)); print(sorted([x for x in r if 13<x[0]<16], key=lambda q:q[1]))
