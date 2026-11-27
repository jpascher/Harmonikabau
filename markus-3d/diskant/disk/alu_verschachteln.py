# Alle Alu-Hebel (Abwicklungen) gemeinsam auf eine Fraesplatte setzen (Johann 27.11.2026)
#   python3 alu_verschachteln.py <ordner_mit_*_Alu_abgewickelt.dxf> <ausgabe.dxf> [plattenbreite_mm] [fraeser_mm] [steg_mm]
# Regel: Schaftfraeser 2 mm, zwischen zwei Teilen bleiben mindestens 2 mm Material stehen -> Abstand der Konturen
# = 2 x Fraeser + Steg (6 mm); zum Plattenrand Fraeser + Steg. Die Teile werden gedreht (0/90/180/270 Grad und die Lage
# des kleinsten umschliessenden Rechtecks) und von unten links dicht gesetzt (Rasterverfahren 0,5 mm), damit die
# Flaeche moeglichst gut genutzt wird. Gespiegelt wird nie (Biegerichtungen bleiben richtig).
# Ausgabe fuer CamBam: Layer Kontur (aussen fraesen), Ausschnitt (innen), Bohrung_D6, Biegelinie (nur anreissen/
# gravieren), Text (Name je Hebel), Platte (Rohplatte, nur zur Orientierung).
import sys, os, glob, math
import numpy as np, ezdxf
from scipy.signal import fftconvolve
from shapely.geometry import Polygon, Point, LineString
from shapely import affinity
from shapely.ops import unary_union

RASTER = 0.5
KONTAKT = 2.0
WINKELSCHRITT = 15     # Grad; zusaetzlich die Lagen des kleinsten umschliessenden Rechtecks

def lies(pfad):
    d = ezdxf.readfile(pfad); m = d.modelspace()
    ringe = [Polygon([(p[0], p[1]) for p in e.get_points()]) for e in m.query('LWPOLYLINE[layer=="Kontur"]')]
    ringe = [r.buffer(0) for r in ringe if r.area > 1e-6]
    aussen = max(ringe, key=lambda r: r.area)
    innen = [r for r in ringe if r is not aussen and aussen.contains(r.representative_point())]
    teil = Polygon(aussen.exterior.coords, [r.exterior.coords for r in innen]).buffer(0)
    kreise = [(Point(c.dxf.center.x, c.dxf.center.y), c.dxf.radius) for c in m.query('CIRCLE')]
    linien = [LineString([(l.dxf.start.x, l.dxf.start.y), (l.dxf.end.x, l.dxf.end.y)]) for l in m.query('LINE[layer=="Biegelinie"]')]
    name = os.path.basename(pfad).replace("_Alu_abgewickelt.dxf", "")
    return {"name": name, "teil": teil, "kreise": kreise, "linien": linien}

def maske(poly, x0, y0, nx, ny):
    import shapely
    xs = x0 + (np.arange(nx) + 0.5) * RASTER; ys = y0 + (np.arange(ny) + 0.5) * RASTER
    X, Y = np.meshgrid(xs, ys)
    return shapely.contains_xy(poly, X, Y)

def verschachteln(teile, breite, abstand, rand, faktor=3.0, schritt=WINKELSCHRITT, regel="kontakt", ordnung=None):
    W = int(math.ceil((breite - 2 * rand) / RASTER))
    H = int(math.ceil(sum(t["teil"].minimum_rotated_rectangle.buffer(abstand).area for t in teile) * faktor / (breite - 2 * rand) / RASTER)) + 400
    belegt = np.zeros((H, W), dtype=np.float32)
    platz = []; oben = 0
    for t in (ordnung or sorted(teile, key=lambda t: -t["teil"].minimum_rotated_rectangle.area)):
        rr = t["teil"].minimum_rotated_rectangle
        c = list(rr.exterior.coords)
        kanten = [(c[i], c[i + 1]) for i in range(4)]
        (a, b) = max(kanten, key=lambda e: math.dist(*e))
        w0 = math.degrees(math.atan2(b[1] - a[1], b[0] - a[0]))
        winkel = sorted({round(w % 360, 3) for w in [-w0, -w0 + 90, -w0 + 180, -w0 + 270] + list(range(0, 360, int(schritt)))})
        best = None
        for w in winkel:
            p = affinity.rotate(t["teil"], w, origin=(0, 0))
            env = p.buffer(abstand / 2.0 + RASTER * 0.75, join_style=1)   # Rasterfehler abdecken
            bx0, by0, bx1, by1 = env.bounds
            nx = int(math.ceil((bx1 - bx0) / RASTER)) + 1; ny = int(math.ceil((by1 - by0) / RASTER)) + 1
            if nx > W: continue
            m = maske(env, bx0, by0, nx, ny).astype(np.float32)
            if regel == "kontakt":
                ring = maske(p.buffer(abstand / 2.0 + KONTAKT, join_style=1), bx0 - KONTAKT, by0 - KONTAKT, nx + 2 * int(KONTAKT / RASTER), ny + 2 * int(KONTAKT / RASTER)).astype(np.float32)
            if oben + ny + 2 > H: continue
            konf = fftconvolve(belegt[:oben + ny + 2], m[::-1, ::-1], mode="valid")          # Ueberlappung fuer jede Lage
            frei = np.argwhere(konf < 0.5)
            if not len(frei): continue
            if regel == "oben":
                iy, ix = min(frei.tolist(), key=lambda v: (v[0] + ny, v[1]))   # Oberkante moeglichst tief, dann links
                wert = ((iy + ny) * RASTER, ix * RASTER)
            elif regel == "kontakt":
                # an schon gesetzte Teile und den Plattenrand anschmiegen; wo die Platte nicht hoeher wird, zaehlt nur die Beruehrung
                k = int(KONTAKT / RASTER)
                umg = np.pad(np.minimum(belegt[:oben + ny + 2], 1.0), ((k, k), (k, k)), constant_values=1.0)
                umg[-k:, :] = 0.0
                ber = fftconvolve(umg, ring[::-1, ::-1], mode="valid")[:konf.shape[0], :konf.shape[1]]
                zuw = np.maximum(frei[:, 0] + ny - oben, 0) * RASTER
                bew = zuw * 3.0 - ber[frei[:, 0], frei[:, 1]] * RASTER * RASTER / KONTAKT + frei[:, 0] * RASTER * 0.05
                j = int(np.argmin(bew)); iy, ix = int(frei[j, 0]), int(frei[j, 1])
                wert = (float(bew[j]), 0.0)
            else:
                iy, ix = min(frei.tolist(), key=lambda v: (v[0] + ny / 2, v[1]))  # Schwerpunkt moeglichst tief
                wert = ((iy + ny / 2) * RASTER, ix * RASTER)
            if best is None or wert < best[0]: best = (wert, w, iy, ix, m, bx0, by0)
        if best is None:
            if faktor < 20 and max(min(affinity.rotate(t["teil"], w, origin=(0, 0)).bounds[2] - affinity.rotate(t["teil"], w, origin=(0, 0)).bounds[0] for w in winkel) for _ in [0]) + abstand <= breite - 2 * rand:
                return verschachteln(teile, breite, abstand, rand, faktor * 2, schritt, regel, ordnung)
            raise SystemExit(f"{t['name']} passt nicht auf {breite} mm Breite")
        _, w, iy, ix, m, bx0, by0 = best
        belegt[iy:iy + m.shape[0], ix:ix + m.shape[1]] += m
        oben = max(oben, iy + m.shape[0])
        dx = rand + ix * RASTER - bx0; dy = rand + iy * RASTER - by0
        platz.append((t, w, dx, dy))
    hoehe = max(affinity.translate(affinity.rotate(t["teil"], w, origin=(0, 0)), dx, dy).bounds[3] for t, w, dx, dy in platz) + rand
    # Kontrolle: echter Mindestabstand zwischen den Konturen und zum Rand
    pp = [affinity.translate(affinity.rotate(t["teil"], w, origin=(0, 0)), dx, dy) for t, w, dx, dy in platz]
    from shapely.strtree import STRtree
    baum = STRtree(pp); dmin = 1e9
    for i, q in enumerate(pp):
        for j in baum.query(q.buffer(abstand)):
            if j > i: dmin = min(dmin, q.distance(pp[j]))
        b = q.bounds; dmin_r = min(b[0], breite - b[2], b[1])
        if dmin_r < rand - 1e-6: raise SystemExit(f"{t['name']}: Randabstand {dmin_r:.2f} mm zu klein")
    if dmin < abstand - 1e-6: raise SystemExit(f"Teileabstand {dmin:.2f} mm < {abstand} mm")
    return platz, hoehe

def schreibe(platz, breite, hoehe, ziel, fraeser, steg):
    d = ezdxf.new(setup=True); d.units = ezdxf.units.MM; m = d.modelspace()
    for lay, col in (("Kontur", 7), ("Ausschnitt", 5), ("Bohrung_D6", 1), ("Biegelinie", 3), ("Text", 2), ("Platte", 8)): d.layers.add(lay, color=col)
    m.add_lwpolyline([(0, 0), (breite, 0), (breite, hoehe), (0, hoehe)], close=True, dxfattribs={"layer": "Platte"})
    tr = lambda g, w, dx, dy: affinity.translate(affinity.rotate(g, w, origin=(0, 0)), dx, dy)
    for t, w, dx, dy in platz:
        p = tr(t["teil"], w, dx, dy)
        m.add_lwpolyline(list(p.exterior.coords), close=True, dxfattribs={"layer": "Kontur"})
        for i in p.interiors: m.add_lwpolyline(list(i.coords), close=True, dxfattribs={"layer": "Ausschnitt"})
        for c, rad in t["kreise"]:
            q = tr(c, w, dx, dy); m.add_circle((q.x, q.y), rad, dxfattribs={"layer": "Bohrung_D6"})
        for l in t["linien"]:
            q = tr(l, w, dx, dy); m.add_line(q.coords[0], q.coords[-1], dxfattribs={"layer": "Biegelinie"})
        z = p.representative_point()
        m.add_text(t["name"], height=2.0, dxfattribs={"layer": "Text", "rotation": 0}).set_placement((z.x, z.y))
    m.add_text(f"Alu-Hebel {len(platz)} Stueck, Platte {breite:g} x {hoehe:.0f} mm, Fraeser {fraeser:g} mm, Steg {steg:g} mm, "
               "nicht gespiegelt", height=3.0, dxfattribs={"layer": "Text"}).set_placement((0, hoehe + 4))
    d.saveas(ziel)

def main(ordner, ziel, breite=300.0, fraeser=2.0, steg=2.0):
    teile = [lies(p) for p in sorted(glob.glob(os.path.join(ordner, "*_Alu_abgewickelt.dxf")))]
    if not teile: raise SystemExit("keine Alu-Abwicklungen gefunden")
    abstand = 2 * fraeser + steg; rand = fraeser + steg
    ordnungen = [sorted(teile, key=lambda t: -t["teil"].minimum_rotated_rectangle.area),
                 sorted(teile, key=lambda t: -t["teil"].convex_hull.area),
                 sorted(teile, key=lambda t: -t["teil"].area)]
    ergebnisse = [verschachteln(teile, float(breite), abstand, rand, ordnung=o) for o in ordnungen]
    platz, hoehe = min(ergebnisse, key=lambda e: e[1])          # kuerzeste Platte gewinnt
    schreibe(platz, float(breite), hoehe, ziel, fraeser, steg)
    flaeche = sum(t["teil"].area for t in teile)
    print(f"Alu-Fraesplatte: {len(teile)} Hebel auf {float(breite):g} x {hoehe:.0f} mm, Nutzung {flaeche / (float(breite) * hoehe) * 100:.0f} %")
    return hoehe

if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[0], a[1], *(float(x) for x in a[2:5]))
