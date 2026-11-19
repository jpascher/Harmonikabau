# Fraesdateien fuer den Griffstock (Teil 1 = Koerper, Teil 2 = Kamm mit Boden), Johann 19.11.2026:
# - Teil in Fraeslage: Deckelflaeche waagrecht, Z nach oben = vom Deckel weg; Oberseite bei Z = 0.
# - Aussenkontur wird nicht gefraest; gefraest wird von oben (vom Deckel weg) und von der Deckelseite.
# - Durchgang 1 (6-mm-Schaftfraeser): Flaechen/Absaetze; Innenraum und Deckel-Ausnehmung als Rechteck.
# - Durchgang 2 (2,5-mm-Schaftfraeser): jeder Kanal einzeln als Rechteck, zur offenen Seite (Freiraum)
#   um den Fraeserdurchmesser verlaengert, damit an der Muendung keine Rundung stehen bleibt.
# Aufruf: python3 gs_fraes.py <griffstock.step> <ausgabeordner> <name>
#   (griffstock.step = Diskantteil mit NUR="Griffstock" erzeugt: Koerper, Teil2, Abdeckung, Schrauben)
import sys, os, math, json
import cadquery as cq, ezdxf
from shapely.geometry import Polygon, box as sbox
from shapely.ops import unary_union

D_KANAL, D_FLAECHE = 2.5, 6.0
D_KANAL_TEIL = {"Teil1": 3.175, "Teil2": 2.5}   # Johann 19.11.: Teil 1 (Kanaele 18,9 tief) mit 3,175-mm-Fraeser
DZ = 0.25                     # Schnittabstand in Z fuer die Kanalerkennung

def lage(pfad):
    a = cq.importers.importStep(pfad)
    S = sorted(a.solids().vals(), key=lambda s: -s.Volume())
    koerper, teil2, abd = S[0], S[1], S[2]
    nv = sorted(abd.Faces(), key=lambda f: -f.Area())[0].normalAt()
    if nv.y > 0: nv = -nv
    rot = 180 - math.degrees(math.atan2(nv.y, nv.z))          # Deckelnormale (nach aussen) -> -Z
    tr = lambda s: s.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -rot)
    return {"Teil1": tr(koerper), "Teil2": tr(teil2)}, rot

def schnitt(teil, z):
    sec = teil.intersect(cq.Solid.makeBox(1000, 1000, 0.002, cq.Vector(-300, -500, z - 0.001)))
    pg = []
    for f in sec.Faces():
        if abs(f.normalAt().z) < 0.9: continue
        ring = lambda w: [(p.x, p.y) for p in (w.positionAt(k / 600) for k in range(601))]
        pg.append(Polygon(ring(f.outerWire()), [ring(w) for w in f.innerWires()]).buffer(0))
    return unary_union(pg) if pg else Polygon()

import statistics as _st
from shapely.geometry import Point

def analyse(secs, d_kanal=D_KANAL):
    """Aus den Z-Schnitten: Innenraum-Rechteck, Deckel-Ausnehmung, Kanal-Rechtecke (je Kanal eines, offen zur Innenseite)."""
    zs = sorted(secs); alle = unary_union([secs[z] for z in zs if not secs[z].is_empty])
    X0, Y0, X1, Y1 = alle.bounds; xm = (X0 + X1) / 2
    yb = alle.intersection(sbox(xm - 60, Y0 - 1, xm + 60, Y1 + 1)).bounds[1]          # Balken-Innenkante
    xl = alle.intersection(sbox(X0 - 1, Y0 - 1, xm - 60, yb - 0.5)).bounds[2]          # Schenkel-Innenkanten
    xr = alle.intersection(sbox(xm + 60, Y0 - 1, X1 + 1, yb - 0.5)).bounds[0]
    innen = (xl, Y0 - D_FLAECHE, xr, yb)
    roh = []
    for z in zs:
        s = secs[z]
        if s.is_empty: continue
        zu = s.buffer(3.4, join_style=2).buffer(-3.4, join_style=2); n = zu.difference(s)
        for g in (n.geoms if hasattr(n, "geoms") else [n]):
            if g.geom_type != "Polygon" or g.area < 1.5: continue
            a, b, c, d = g.bounds
            if not (3.0 < c - a < 6.5) or d - b < 8.0 or d < yb + 3.0: continue
            oo = not zu.contains(Point((a + c) / 2, d + 0.6)) and not s.contains(Point((a + c) / 2, d + 0.6))
            roh.append((z, a, b, c, d, oo))
    roh.sort(key=lambda r: (r[1] + r[3]) / 2)
    gruppen = []
    for r in roh:
        xm_ = (r[1] + r[3]) / 2
        if gruppen and abs(xm_ - gruppen[-1][-1][0]) < 2.0: gruppen[-1].append((xm_, r))
        else: gruppen.append([(xm_, r)])
    kan = []
    for G in gruppen:
        L = [r for _, r in G]
        if len(L) < 3: continue
        x0 = _st.median([r[1] for r in L]); x1 = _st.median([r[3] for r in L])
        zu_ = [r for r in L if not r[5]]
        y1 = max(r[4] for r in (zu_ or L))
        y0 = min(min(r[2] for r in L), yb) - d_kanal            # zur offenen Seite: ab der Balken-Innenkante (frei in jeder Hoehe) + Fraeserdurchmesser
        kan.append({"r": (x0, y0, x1, y1), "z0": min(r[0] for r in L) - DZ / 2, "z1": max(r[0] for r in L) + DZ / 2})
    s0 = secs[zs[0]]; da = None
    m0 = s0.intersection(sbox(xm - 60, Y0 - 1, xm + 60, Y1 + 1))
    if not m0.is_empty and m0.bounds[1] > yb + 5:                 # unterste Schicht im Mittelbereich frei: Deckel-Ausnehmung
        yd = m0.bounds[1]
        xdl = s0.intersection(sbox(X0 - 1, Y0 - 1, xm - 60, yd - 0.5)).bounds[2]
        xdr = s0.intersection(sbox(xm + 60, Y0 - 1, X1 + 1, yd - 0.5)).bounds[0]
        zt = next(z for z in zs if not secs[z].intersection(sbox(xm - 60, yb + 1, xm + 60, yd - 1)).is_empty)
        da = {"r": (xdl, Y0 - D_FLAECHE, xdr, yd), "tiefe": zt - DZ / 2 - (zs[0] - DZ / 2)}
    return {"innen": innen, "kanaele": kan, "deckel": da, "umriss": alle}

SCHLICHT_T = 0.5     # Johann 19.11.: Stufen der Vorderseite zuletzt mit 0,5 mm Tiefe schlichten

def stufen(teil):
    """Knopfreihen-Stufen der Vorderseite (oben): waagrechte Flaechen, die ueber die Mitte laufen und mindestens 10 mm breit sind.
    Schlicht-Rechteck: in X (seitlich) um den Fraeserdurchmesser ueber die Stufe hinaus, in Y genau die Stufe (Stufen bleiben)."""
    bb = teil.BoundingBox(); xm = (bb.xmin + bb.xmax) / 2
    gr = {}
    for f in teil.Faces():
        if f.geomType() != "PLANE" or f.normalAt().z < 0.99 or f.Area() < 5: continue
        gr.setdefault(round(f.Center().z, 2), []).append(f.BoundingBox())
    out = []
    for z, B in sorted(gr.items(), reverse=True):
        x0, x1 = min(b.xmin for b in B), max(b.xmax for b in B)
        y0, y1 = min(b.ymin for b in B), max(b.ymax for b in B)
        if y1 - y0 < 10.0 or z < bb.zmax - 10.0: continue
        if not any(b.xmin < xm < b.xmax for b in B) and not (x0 < xm < x1 and len(B) > 3): continue
        out.append({"z": z, "r": (x0 - D_FLAECHE, y0, x1 + D_FLAECHE, y1)})
    return out[:2]          # je Teil zwei Knopfreihen-Stufen (Teil 1: Reihe 1/2, Teil 2: Reihe 3/4), 4 mm Hoehenunterschied

def flaechen(teil):
    """Alle waagrechten Flaechen, die der 6-mm-Fraeser erreicht (mind. 6,5 mm breit), getrennt nach Seite:
    oben (Normale +Z) und unten/Deckelseite (Normale -Z). Je Hoehe Rechtecke (benachbarte Flaechen gleicher Y-Lage
    zusammengefasst, Kanaele dazwischen ueberbrueckt). In X (seitlich) um den Fraeserdurchmesser verlaengert, aber nur
    dort, wo neben der Flaeche nichts hoeher steht (Stufen bleiben erhalten). Johann 19.11.: Schlichtlauf 0,5 mm."""
    bb = teil.BoundingBox(); out = {"oben": [], "unten": []}
    for seite, sg in (("oben", 1), ("unten", -1)):
        gr = {}
        for f in teil.Faces():
            if f.geomType() != "PLANE" or f.normalAt().z * sg < 0.99: continue
            b = f.BoundingBox()
            if b.ylen < 2.0 or b.xlen < 1.0: continue           # schmale Zaehne bleiben drin: zusammengefasst bilden sie eine Auflage
            z_ = f.Center().z
            k_ = next((k for k in gr if abs(k - z_) < 0.06), None)                # gleiche Ebene (z. B. -7,79 / -7,81)
            gr.setdefault(k_ if k_ is not None else round(z_, 2), []).append([b.xmin, b.ymin, b.xmax, b.ymax])
        for z in sorted(gr, reverse=(sg > 0)):
            R = sorted(gr[z], key=lambda q: (round(q[1]), round(q[3]), q[0]))
            zus = []
            for q in R:
                if zus and abs(zus[-1][1] - q[1]) < 1.0 and abs(zus[-1][3] - q[3]) < 1.0 and q[0] - zus[-1][2] < D_FLAECHE + 0.5:
                    zus[-1] = [zus[-1][0], min(zus[-1][1], q[1]), max(zus[-1][2], q[2]), max(zus[-1][3], q[3])]
                else: zus.append(list(q))
            zus = [q for q in zus if q[2] - q[0] >= D_FLAECHE + 0.2 and q[3] - q[1] >= D_FLAECHE + 0.2]   # Fraeser muss hineinpassen
            if not zus: continue
            rechtecke = []
            for x0, y0, x1, y1 in zus:
                def frei(xa, xb):
                    za, zb = (z + 0.05, bb.zmax + 1) if sg > 0 else (bb.zmin - 1, z - 0.05)
                    try: return teil.intersect(cq.Solid.makeBox(xb - xa, y1 - y0 - 0.2, zb - za, cq.Vector(xa, y0 + 0.1, za))).Volume() < 0.5
                    except Exception: return False
                xa = x0 - D_FLAECHE if frei(x0 - D_FLAECHE, x0 - 0.05) else x0
                xb = x1 + D_FLAECHE if frei(x1 + 0.05, x1 + D_FLAECHE) else x1
                rechtecke.append((xa, y0, xb, y1))
            out[seite].append({"z": z, "r": rechtecke})
    return out

def schnitte(teil):
    bb = teil.BoundingBox()
    return {bb.zmin + DZ / 2 + i * DZ: schnitt(teil, bb.zmin + DZ / 2 + i * DZ) for i in range(int((bb.zmax - bb.zmin) / DZ))}

def main(pfad, ordner, name):
    os.makedirs(ordner, exist_ok=True)
    teile, rot = lage(pfad)
    info = {"name": name, "drehung_um_X_grad": round(rot, 3), "fraeser": {"flaeche": D_FLAECHE, "kanal": D_KANAL}, "teile": {}}
    for tn, t in teile.items():
        bb = t.BoundingBox()
        t = t.translate(cq.Vector(-bb.xmin, 0, -bb.zmax))          # Oberseite Z = 0, Anfang X = 0
        bb = t.BoundingBox()
        dk = D_KANAL_TEIL.get(tn, D_KANAL)
        A = analyse(schnitte(t), dk)
        if tn != "Teil1": A["deckel"] = None          # nur Teil 1 liegt am Deckel (Teil 2 beginnt an der Trennebene)
        doc = ezdxf.new(); msp = doc.modelspace()
        def rect(r, layer):
            x0, y0, x1, y1 = r
            msp.add_lwpolyline([(x0, y0), (x1, y0), (x1, y1), (x0, y1)], close=True, dxfattribs={"layer": layer})
        um = A["umriss"]
        for g in (um.geoms if hasattr(um, "geoms") else [um]):
            msp.add_lwpolyline(list(g.exterior.coords), close=True, dxfattribs={"layer": "Umriss_nicht_fraesen"})
        dicke = bb.zmax - bb.zmin
        rect(A["innen"], f"Innenraum_Rechteck_durch_T{dicke:.1f}_D{D_FLAECHE:g}")
        if A["deckel"]:     # Johann 19.11.: Deckel-Ausnehmung vollflaechig ausraeumen, ausser dem Innenraum (wird ohnehin ausgeschnitten)
            _ta = sbox(*A["deckel"]["r"]).difference(sbox(*A["innen"]))
            for g_ in (_ta.geoms if hasattr(_ta, "geoms") else [_ta]):
                msp.add_lwpolyline(list(g_.exterior.coords), close=True,
                                   dxfattribs={"layer": f"Deckel_Ausnehmung_vollflaechig_T{A['deckel']['tiefe']:.1f}_Deckelseite_D{D_FLAECHE:g}"})
        for k_ in A["kanaele"]:
            rect(k_["r"], f"Kanal_D{dk:g}_von_oben_T{-k_['z0']:.1f}")
        ST = stufen(t)
        for st_ in ST:
            rect(st_["r"], f"Stufe_schlichten_Z{st_['z']:.2f}_T{SCHLICHT_T:g}_D{D_FLAECHE:g}")
        doc.saveas(os.path.join(ordner, f"Griffstock_{tn}_{name}_fraesen.dxf"))
        t.exportStep(os.path.join(ordner, f"Griffstock_{tn}_{name}_fraesen.step"))
        r2 = lambda q: [round(v, 2) for v in q]
        FL = flaechen(t)
        for seite_, L_ in FL.items():
            for g_ in L_:
                for r_ in g_["r"]: rect(r_, f"Flaeche_{seite_}_Z{g_['z']:.2f}_schlichten_T{SCHLICHT_T:g}_D{D_FLAECHE:g}")
        info["teile"][tn] = {"flaechen": {k_: [{"z": g_["z"], "r": [r2(r_) for r_ in g_["r"]], "tiefe": SCHLICHT_T} for g_ in v_] for k_, v_ in FL.items()},
                             "stufen_schlichten": [{"z": st_["z"], "r": r2(st_["r"]), "tiefe": SCHLICHT_T} for st_ in ST], "fraeser_kanal": dk, "dicke": round(dicke, 2), "laenge": round(bb.xlen, 2), "breite": round(bb.ylen, 2),
                             "innenraum": r2(A["innen"]),
                             "deckel_ausnehmung": ({"r": r2(A["deckel"]["r"]), "tiefe": round(A["deckel"]["tiefe"], 2)} if A["deckel"] else None),
                             "kanaele": [{"r": r2(k_["r"]), "tiefe_von_oben": round(-k_["z0"], 2), "oben_bis": round(-k_["z1"], 2)} for k_ in A["kanaele"]]}
        print(f"{tn}: Dicke {dicke:.1f}, Kanaele {len(A['kanaele'])}, Innenraum {r2(A['innen'])}, Deckel {info['teile'][tn]['deckel_ausnehmung']}")
    json.dump(info, open(os.path.join(ordner, f"Griffstock_{name}_fraesen.json"), "w"), indent=1)

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
