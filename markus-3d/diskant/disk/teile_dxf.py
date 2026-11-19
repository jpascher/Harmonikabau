# Fertigungszeichnungen (Johann 19.11.2026), aus den Daten des Generators (TEILE_DUMP, KST_STEP):
#  1. je Alu-Hebel eine DXF der Abwicklung (gestreckt): Kontur, Bohrung Ø6 (Messingbuchse), Biegelinien mit Winkel
#     Biegezugabe: Halbrund 5 mm (Innenradius 2,5), neutrale Faser k = 0,4, Blech 2 mm
#  2. je Reihe eine DXF des Holzteils (alle Holzteile einer Reihe gleich) mit Stueckzahl, Ein-/Ausfahrbogen R3 an der
#     laengsten geraden Kante, darunter das Alu-Teil als Referenz mit Bohrung
#  3. Kunststoffteil: STEP als Fraesteil (Rohblock: Abschraegungen werden nicht gefraest, nur Nuten, Taschen,
#     Durchbrueche, Rastmarken) und Fertigteil in Fraeslage; DXF Draufsicht auf die Lagerseite mit Taschen je Tiefe
# Aufruf: python3 teile_dxf.py <teile.pkl> <kunststoff.step> <ausgabeordner>
import sys, os, math, pickle, collections
import ezdxf, cadquery as cq
from shapely.geometry import Polygon, LineString, Point, box as sbox
from shapely.ops import unary_union

K_FAKTOR, R_INNEN = 0.4, 2.5        # Johann 19.11.: gebogen ueber Halbrund 5 mm -> Innenradius 2,5

def _ausgleich(theta, t):
    """Laengenkorrektur je Biegung gegenueber der scharfkantigen Mittellinie des Modells"""
    th = abs(theta)
    return th * (R_INNEN + K_FAKTOR * t) - 2 * (R_INNEN + t / 2) * math.tan(th / 2)

def _poly_pts(p): return [(float(a), float(b)) for a, b in p.exterior.coords]

def abwicklung(h, t, kvu):
    pts = h["pts"]                                   # [(Y, X)] Pfad des Arms ab ya
    seg = []                                         # (Y0, Y1, Steigung kk)
    for (ya_, xa_), (yb_, xb_) in zip(pts, pts[1:]):
        if yb_ - ya_ > 1e-6: seg.append((ya_, yb_, (xb_ - xa_) / (yb_ - ya_)))
    # Biegungen in Laengsrichtung (um senkrechte Achsen): an jedem Richtungswechsel, beginnend mit der Platte (kk = 0)
    knicke = []; k_vor = 0.0
    for y0_, y1_, kk in seg:
        th = math.atan(kk) - math.atan(k_vor)
        if abs(th) > math.radians(0.5): knicke.append((y0_, th))
        k_vor = kk
    def s_von(Y):
        s = Y; korr = 0.0
        for y0_, y1_, kk in seg:
            if Y <= y0_: break
            dy = min(Y, y1_) - y0_
            s += dy * (math.sqrt(1 + kk * kk) - 1)
        for yk, th in knicke:
            if Y > yk: korr += _ausgleich(th, t)
        return s + korr
    # Kroepfung (Reihen mit Knopfversatz): Biegungen um Achsen in Laengsrichtung bei u = U0 und U1
    u0, u1 = kvu; kv = h["kv"]
    g = math.sqrt(1 + (kv / (u1 - u0)) ** 2) if kv else 1.0
    thk = math.atan(kv / (u1 - u0)) if kv else 0.0
    def u_von(u):
        if not kv or u <= u0: return u
        if u <= u1: return u0 + (u - u0) * g + _ausgleich(thk, t)
        return u + (u1 - u0) * (g - 1) + 2 * _ausgleich(thk, t)
    proj = h["proj"]
    if proj.geom_type != "Polygon": proj = max(proj.geoms, key=lambda q: q.area)
    ring = LineString(proj.exterior.coords).segmentize(0.25)
    flach = Polygon([(s_von(Y), u_von(u)) for Y, u in ring.coords]).buffer(0)
    if flach.geom_type != "Polygon": flach = max(flach.geoms, key=lambda q: q.area)     # Splitter (Flaeche ~0) weg
    loch = (s_von(h["Yp"]), u_von(h["up"]))
    linien = []
    for yk, th in knicke:
        sk = s_von(yk - 1e-6) + _ausgleich(th, t) / 2
        L_ = LineString([(sk, -500), (sk, 500)]).intersection(flach)
        if not L_.is_empty: linien.append((L_, math.degrees(th), "laengs"))
    if kv:
        for uu, sg in ((u0, 1), (u1, -1)):
            uk = u_von(uu - 1e-6) + _ausgleich(thk, t) / 2
            L_ = LineString([(-500, uk), (900, uk)]).intersection(flach)
            if not L_.is_empty: linien.append((L_, sg * math.degrees(thk), "kroepf"))
    return flach, loch, linien

def _linien(msp, geo, layer):
    for g in (geo.geoms if hasattr(geo, "geoms") else [geo]):
        if g.geom_type == "LineString": msp.add_line(g.coords[0], g.coords[-1], dxfattribs={"layer": layer})

def _neu():
    d = ezdxf.new(setup=True); d.units = ezdxf.units.MM
    return d, d.modelspace()

def alu_dxf(D, ordner):
    t = D["T_ALU"]; rb = D["D_BUCHSE"] / 2; n = 0
    for h in D["hebel"]:
        flach, loch, linien = abwicklung(h, t, D["KV_U"])
        doc, msp = _neu()
        for lay, c in (("Kontur", 7), ("Bohrung_D6_Messingbuchse", 1), ("Biegelinie", 3), ("Text", 2)): doc.layers.add(lay, color=c)
        msp.add_lwpolyline(_poly_pts(flach), close=True, dxfattribs={"layer": "Kontur"})
        for inner in flach.interiors: msp.add_lwpolyline(list(inner.coords), close=True, dxfattribs={"layer": "Kontur"})
        msp.add_circle(loch, rb, dxfattribs={"layer": "Bohrung_D6_Messingbuchse"})
        b = flach.bounds; ty = b[1] - 6
        msp.add_text(f"{h['G']} Alu {t:g} mm abgewickelt (k={K_FAKTOR:g}, Ri={R_INNEN:g})", height=2.5,
                     dxfattribs={"layer": "Text"}).set_placement((b[0], ty))
        for i, (L_, wink, art) in enumerate(linien):
            _linien(msp, L_, "Biegelinie")
            p = L_.centroid if not L_.is_empty else Point(0, 0)
            seite = ("zur Federseite (+X)" if wink > 0 else "von der Federseite weg (-X)") if art == "laengs" else ("Kroepfung " + ("hin" if wink > 0 else "zurueck"))
            msp.add_text(f"B{i+1}: {abs(wink):.1f} Grad {seite}", height=1.8, dxfattribs={"layer": "Text", "rotation": 90 if art == "laengs" else 0}).set_placement((p.x + 1.0, p.y))
        doc.saveas(os.path.join(ordner, f"{h['G']}_Alu_abgewickelt.dxf")); n += 1
    return n

def _ein_aus(msp, poly, R=3.0):
    co = list(poly.exterior.coords)
    i = max(range(len(co) - 1), key=lambda j: math.dist(co[j], co[j + 1]))
    (x0, y0), (x1, y1) = co[i], co[i + 1]
    L = math.dist(co[i], co[i + 1]); dx, dy = (x1 - x0) / L, (y1 - y0) / L
    M = ((x0 + x1) / 2, (y0 + y1) / 2)
    n_ = (dy, -dx)                                                     # eine Normale; nach aussen pruefen
    if poly.contains(Point(M[0] + n_[0] * 0.5, M[1] + n_[1] * 0.5)): n_ = (-n_[0], -n_[1])
    C = (M[0] + R * n_[0], M[1] + R * n_[1])
    def ang(v): return math.degrees(math.atan2(v[1], v[0]))
    a_m, a_ein, a_aus = ang((-n_[0], -n_[1])), ang((-dx, -dy)), ang((dx, dy))
    for a0, a1, lay in ((a_ein, a_m, "Einfahren_R3"), (a_m, a_aus, "Ausfahren_R3")):
        # kuerzeren Bogen zwischen a0 und a1 (90 Grad) zeichnen
        if (a1 - a0) % 360 > 180: a0, a1 = a1, a0
        msp.add_arc(C, R, a0, a1, dxfattribs={"layer": lay})
    msp.add_point(M, dxfattribs={"layer": "Einfahren_R3"})

def holz_dxf(D, ordner):
    anz = collections.Counter(h["r"] for h in D["hebel"]); n = 0
    for r in sorted(anz):
        hp = D["holz"].get(r)
        if hp is None: continue
        hp = Polygon([(Y, u) for u, Y in hp.exterior.coords]).buffer(0)                   # (u, Y) -> Zeichnung (Y, u)
        h0 = [h for h in D["hebel"] if h["r"] == r][0]
        doc, msp = _neu()
        for lay, c in (("Kontur_Holz", 7), ("Einfahren_R3", 3), ("Ausfahren_R3", 1), ("Alu_Referenz", 8), ("Bohrung_Referenz_D6", 1), ("Text", 2)): doc.layers.add(lay, color=c)
        ap = h0["proj"] if h0["proj"].geom_type == "Polygon" else max(h0["proj"].geoms, key=lambda q: q.area)
        msp.add_lwpolyline(_poly_pts(ap), close=True, dxfattribs={"layer": "Alu_Referenz"})
        msp.add_circle((h0["Yp"], h0["up"]), D["D_BUCHSE"] / 2, dxfattribs={"layer": "Bohrung_Referenz_D6"})
        msp.add_lwpolyline(_poly_pts(hp), close=True, dxfattribs={"layer": "Kontur_Holz"})
        _ein_aus(msp, hp)
        b = hp.bounds
        msp.add_text(f"Holzteil Reihe {r}: {anz[r]} Stueck, Dicke {D['T_HOLZ']:g} mm (Alu darunter nur Referenz)", height=2.5,
                     dxfattribs={"layer": "Text"}).set_placement((b[0], b[1] - 6))
        doc.saveas(os.path.join(ordner, f"Holzteil_Reihe{r}_{anz[r]}Stueck.dxf")); n += 1
    return n

def _schnitt(teil, z):
    sec = teil.intersect(cq.Solid.makeBox(2000, 2000, 0.002, cq.Vector(-1000, -1000, z - 0.001)))
    pg = []
    for f in sec.Faces():
        if abs(f.normalAt().z) < 0.9: continue
        def ring(w):
            try: return [(p.x, p.y) for p in (w.positionAt(k / 300) for k in range(301))]
            except Exception:
                return [(v.X, v.Y) for v in w.Vertices()]           # entartete Kante: nur die Eckpunkte
        try:
            inn = [r for r in (ring(w) for w in f.innerWires()) if len(r) > 2]
            pg.append(Polygon(ring(f.outerWire()), inn).buffer(0))
        except Exception: pass
    return unary_union(pg) if pg else Polygon()

def kugel_bahnen(doc, msp, teil, aus=3.0):
    """Laengsbohrungen Ø2 unter den Lagernuten: Kugelkopffraeser Ø2 mit duennem Hals (Lollipop), Hals laeuft in der
    1,5-mm-Nut, Kugel in Achsmitte ueber die ganze Laenge; je Bahn eine Linie mit Ein-/Ausfahren an den Stirnseiten."""
    bb = teil.BoundingBox(); achsen = []
    for f in teil.Faces():
        if f.geomType() != "CYLINDER": continue
        b = f.BoundingBox()
        if b.xlen < 5 or max(b.ylen, b.zlen) > 2.2: continue
        y, z = (b.ymin + b.ymax) / 2, (b.zmin + b.zmax) / 2
        if b.zmax > -0.01: z = b.zmax - 1.0                      # oben offen (Nut): Mitte = Oberkante Rundung - Radius
        if not any(abs(y - a[0]) < 0.3 and abs(z - a[1]) < 0.3 for a in achsen): achsen.append((y, z))
    for y, z in sorted(achsen):
        lay = f"Kugelfraeser_D2_Laengsbahn_T{-z:.2f}"
        if lay not in doc.layers: doc.layers.add(lay, color=4)
        msp.add_line((bb.xmin - aus, y), (bb.xmax + aus, y), dxfattribs={"layer": lay})
    return achsen

def kunststoff(kst_step, ordner, name):
    """kst_step: Fraesteil in Fraeslage (vom Generator, Lagerseite oben Z = 0, Abschraegungen als Rohblock aufgefuellt);
    daneben <...>_fertig.step = Fertigteil in gleicher Lage."""
    import shutil
    teil = cq.importers.importStep(kst_step).val(); bb = teil.BoundingBox()
    shutil.copy(kst_step, os.path.join(ordner, f"Kunststoffteil_{name}_fraesen.step"))
    shutil.copy(kst_step.replace(".step", "_fertig.step"), os.path.join(ordner, f"Kunststoffteil_{name}_fertig_in_Fraeslage.step"))
    roh = sbox(bb.xmin, bb.ymin, bb.xmax, bb.ymax)
    stufen = []; vor = None; z = -0.1
    while z > bb.zmin + 0.05:
        q = _schnitt(teil, z)
        F = roh.difference(q).buffer(-0.01).buffer(0.01)
        F = unary_union([g for g in (getattr(F, "geoms", None) or [F]) if g.geom_type == "Polygon" and g.area > 0.2])
        neu_ = vor is None or len(getattr(F, "geoms", [F])) != len(getattr(vor[1], "geoms", [vor[1]])) or abs(F.area - vor[1].area) > 0.1 * max(vor[1].area, 1.0)
        if neu_: stufen.append([-z, F]); vor = stufen[-1]
        else: vor[0] = -z; vor[1] = vor[1].union(F)
        z -= 0.25
    doc, msp = _neu()
    doc.layers.add("Rohblock_nicht_fraesen", color=8)
    msp.add_lwpolyline([(bb.xmin, bb.ymin), (bb.xmax, bb.ymin), (bb.xmax, bb.ymax), (bb.xmin, bb.ymax)], close=True, dxfattribs={"layer": "Rohblock_nicht_fraesen"})
    info = []
    for k, (tiefe, F) in enumerate(stufen):
        if F.is_empty: continue
        lay = f"Tasche_{k+1}_bis_T{tiefe:.2f}" + ("_durch" if tiefe > -bb.zmin - 0.5 else "")
        doc.layers.add(lay, color=(k % 6) + 1)
        for g in (getattr(F, "geoms", None) or [F]):
            if g.geom_type == "Polygon": msp.add_lwpolyline(list(g.simplify(0.02).exterior.coords), close=True, dxfattribs={"layer": lay})
        info.append((round(tiefe, 2), len(getattr(F, "geoms", [F]))))
    kugel_bahnen(doc, msp, teil)
    doc.saveas(os.path.join(ordner, f"Kunststoffteil_{name}_fraesen.dxf"))
    return info, (round(bb.xlen, 1), round(bb.ylen, 1), round(-bb.zmin, 1))

def main(dump, kst_step, ordner):
    os.makedirs(ordner, exist_ok=True)
    D = pickle.load(open(dump, "rb"))
    ad = os.path.join(ordner, "Alu_Hebel_abgewickelt"); os.makedirs(ad, exist_ok=True)
    hd = os.path.join(ordner, "Holzteile"); os.makedirs(hd, exist_ok=True)
    print("Alu-Abwicklungen:", alu_dxf(D, ad))
    print("Holzteile:", holz_dxf(D, hd))
    print("Kunststoffteil:", kunststoff(kst_step, ordner, D["instrument"]))

if __name__ == "__main__":
    main(*sys.argv[1:4])
