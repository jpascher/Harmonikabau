# Vereinfachtes Fliegengitter: duenne Platte mit quadratischen Durchbruechen (Raster/Faden), in 2D (shapely) aufgebaut
# und einmal extrudiert - schnell, ein Koerper.  ebene "XY": (u,v)=(X,Y), Dicke in +Z ab h0; "XZ": (u,v)=(X,Z), Dicke in +Y ab h0.
import cadquery as cq
from shapely.geometry import Polygon, box as sbox
def gaze_platte(umriss, ebene, h0, dicke=0.3, raster=2.0, faden=0.2):
    innen = umriss.buffer(-faden)
    x0, y0, x1, y1 = umriss.bounds; loecher = []
    a = raster - faden
    nx, ny = int((x1 - x0) / raster) + 2, int((y1 - y0) / raster) + 2
    for i in range(nx):
        for j in range(ny):
            u, v = x0 + faden / 2 + i * raster, y0 + faden / 2 + j * raster
            q = sbox(u, v, u + a, v + a)
            if innen.contains(q): loecher.append(q)
    def P(u, v): return cq.Vector(u, v, h0) if ebene == "XY" else cq.Vector(u, h0, v)
    def wire(ring): return cq.Wire.makePolygon([P(u, v) for u, v in list(ring.coords)[:-1]], close=True)
    aussen = wire(umriss.exterior)
    inn = [wire(q.exterior) for q in loecher]
    f = cq.Face.makeFromWires(aussen, inn)
    d = cq.Vector(0, 0, dicke) if ebene == "XY" else cq.Vector(0, dicke, 0)
    return cq.Solid.extrudeLinear(f, d), len(loecher)
