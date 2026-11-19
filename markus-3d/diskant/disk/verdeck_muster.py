"""Verdeckmuster aus einer DXF-Zeichnung (Deckel, Seitenteile, Front in der Abwicklung wie verdeck.dxf) lesen.
Liefert die Durchbrueche als shapely-Polygone (mit Inseln) in den Koordinaten von verdeck.dxf:
Deckel x -359,1..0 / y 0..168 (zentriert), Seitenteile links bis x -360,2 bzw. rechts ab x 1,1 (Abstand zur Kante bleibt),
Front y < -1,1 (Abstand zur Kante bleibt, in x zentriert wie der Deckel)."""
import ezdxf, math
from ezdxf import bbox
from shapely.geometry import Polygon, LineString, MultiPolygon
from shapely.ops import unary_union, polygonize
from shapely import affinity

M_DECKEL = (-359.1, 0.0, 0.0, 168.0)     # verdeck.dxf
M_XK_L, M_XK_R, M_YK_F = -360.2, 1.1, -1.1


def _linien(e):
    if e.dxftype() == "CIRCLE":
        c, r = e.dxf.center, e.dxf.radius
        return [LineString([(c.x + r * math.cos(2 * math.pi * k / 64), c.y + r * math.sin(2 * math.pi * k / 64)) for k in range(65)])]
    try:
        p = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.03)]
    except Exception:
        return []
    return [LineString(p)] if len(p) >= 2 else []


def lade(datei):
    m = ezdxf.readfile(datei).modelspace()
    ents = [e for e in m if e.dxftype() in ("LWPOLYLINE", "POLYLINE", "CIRCLE", "LINE", "ARC", "SPLINE", "ELLIPSE")]
    gr = sorted(((bbox.extents([e]), e) for e in ents), key=lambda q: -(q[0].size.x * q[0].size.y) if q[0].has_data else 0)
    deckel = gr[0][0]
    D = (deckel.extmin.x, deckel.extmax.x, deckel.extmin.y, deckel.extmax.y)
    umriss = {id(gr[0][1])}
    links = rechts = front = None
    for b, e in gr[1:8]:          # die naechstgroessten geschlossenen Umrisse: Seitenteile und Front
        if not b.has_data: continue
        if b.extmax.x <= D[0] + 0.5 and b.size.y > 80: links = b; umriss.add(id(e))
        elif b.extmin.x >= D[1] - 0.5 and b.size.y > 80: rechts = b; umriss.add(id(e))
        elif b.extmax.y <= D[2] + 0.5 and b.size.x > 0.6 * (D[1] - D[0]): front = b; umriss.add(id(e))
    lin = [l for e in ents if id(e) not in umriss for l in _linien(e)]
    flaechen = list(polygonize(unary_union(lin)))
    # gerade/ungerade Fuellung: Durchbruch = Flaechen in ungerader Verschachtelungstiefe
    loch = None
    for f in flaechen:
        loch = f if loch is None else loch.symmetric_difference(f)
    teile = [g for g in (loch.geoms if hasattr(loch, "geoms") else [loch]) if g.area > 0.5] if loch else []
    out = {"deckel": [], "links": [], "rechts": [], "front": []}
    dxm, dym = (M_DECKEL[0] + M_DECKEL[1]) / 2 - (D[0] + D[1]) / 2, (M_DECKEL[2] + M_DECKEL[3]) / 2 - (D[2] + D[3]) / 2
    for g in teile:
        c = g.centroid
        if D[0] < c.x < D[1] and D[2] < c.y < D[3]:
            out["deckel"].append(affinity.translate(g, dxm, dym))
        elif links is not None and c.x < D[0]:
            out["links"].append(affinity.translate(g, M_XK_L - links.extmax.x, dym))
        elif rechts is not None and c.x > D[1]:
            out["rechts"].append(affinity.translate(g, M_XK_R - rechts.extmin.x, dym))
        elif front is not None and c.y < D[2]:
            out["front"].append(affinity.translate(g, dxm, M_YK_F - front.extmax.y))
    return out, {"deckel": D, "links": links is not None, "rechts": rechts is not None, "front": front is not None}


if __name__ == "__main__":
    import sys
    o, info = lade(sys.argv[1])
    print(info, {k: len(v) for k, v in o.items()})
