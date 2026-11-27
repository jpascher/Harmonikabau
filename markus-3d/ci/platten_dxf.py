# Alle 2D-Fraesteile (Platten, Leisten, Sohlen, Begleiterstimmstoecke ...) eines Modells als DXF fuer CamBam
# (Johann 27.11.2026: "alles frästeile als dxf für cambam ... in fusion werden nur diskantstimmstöcke und griffstöcke
#  im cam bearbeitet, alle 2d im cambam wie auch die stimmstock solen und begleiterstimmstöcke")
#
#   python3 platten_dxf.py modell.step zielordner [--teil bass|diskant]
#
# Je Teil: Dickenachse = Normale der groessten ebenen Flaeche. Ein Teil gilt als 2D-Fraesteil (2,5D), wenn alle
# Flaechen entweder waagrecht (Normale = Dickenachse) oder senkrecht dazu stehen (Waende, Bohrungen). Dann wird es in
# Scheiben zwischen den waagrechten Ebenen zerlegt und daraus abgeleitet:
#   Kontur_aussen_durch      Aussenkontur, ganz durch (Dicke steht im Layernamen der Datei/Text)
#   Durchbruch               innere Ausschnitte ganz durch
#   Bohrung_D<d>_durch       runde Durchbrueche als Kreis (CamBam: Bohren oder Kreistasche)
#   Tasche_oben_T<t>         Taschen/Absaetze von oben, Tiefe t ab Oberseite
#   Tasche_unten_T<t>        von der Unterseite (eigene Datei ..._Unterseite.dxf, gespiegelt = Ansicht von unten)
# Gleiche Teile werden zusammengefasst (Stueckzahl im Dateinamen), Dateien nach Dicke in Ordnern Platten_<d>mm.
# Nicht erfasst (eigene Wege): Diskant-Stimmstoecke und Griffstoecke (Fusion-CAM), Alu-Hebel und Holz-Hebelteile
# (teile_dxf.py), Schrauben, Achsen, Federn, Leder, Filz, Stimmplatten usw. Teile mit Schraegen/Rundungen in
# Dickenrichtung stehen in nicht_2D.txt.
import sys, os, re, math, argparse, collections
import cadquery as cq, ezdxf
from shapely.geometry import Polygon, MultiPolygon
from shapely.ops import unary_union
from shapely import affinity

AUSSCHLUSS = re.compile(r"schraube|mutter|scheibe_m|leder|filz|epoxy|achse|feder|knopf|schubstange|klappenhebel|stimmplatte|"
                        r"folie|gaze|stift|niete|gummi|pfropfen|ventil|zunge|haken|gelenk|bolzen|druckstange|kugel|magnet|"
                        r"polster|buchse|draht|stange|_arm$|luftklappe_arm|hebel_.*alu|alu_|kunststoff", re.I)
NUR_FUSION = re.compile(r"griffstock|^stimmstock|diskant_stimmstock", re.I)       # Fusion-CAM (Diskant)
EIGENE_DXF = re.compile(r"^holzteil", re.I)                                      # Diskant: teile_dxf.py (je Reihe)
SOHLE = re.compile(r"sohle|sole", re.I)

def blaetter(a, loc=None, pfad=""):
    loc = (loc * a.loc) if loc is not None else a.loc
    p = (pfad + "/" + a.name) if pfad else a.name
    out = []
    if a.obj is not None:
        objs = a.obj.vals() if isinstance(a.obj, cq.Workplane) else [a.obj]
        for o in objs:
            if isinstance(o, cq.Shape):
                for s in o.moved(loc).Solids(): out.append((p, s))
    for c in a.children: out += blaetter(c, loc, p)
    return out

def kurzname(pfad):
    n = pfad.split("/")[-1]
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", n)

def ausgeschlossen(name, teil):
    n = name.split("/")[-1]
    if SOHLE.search(n): return None
    if teil == "diskant" and NUR_FUSION.search(n): return "Fusion-CAM"
    if teil == "diskant" and EIGENE_DXF.search(n): return "Holzteile (eigene DXF)"
    if AUSSCHLUSS.search(n): return "kein Fraesteil"
    return None

def richtungen(s):
    """Kandidaten fuer die Dickenachse: Normalen der ebenen Flaechen, nach Flaeche sortiert (je Richtung einmal)."""
    fl = sorted([f for f in s.Faces() if f.geomType() == "PLANE"], key=lambda f: -f.Area())
    out = []
    for f in fl:
        n = f.normalAt()
        if not any(abs(abs(n.dot(m.normalAt())) - 1) < 1e-4 for m in out): out.append(f)
        if len(out) == 3: break
    return out

def lage(s, f=None):
    """Lokales System: Z = Dickenachse (Normale der Flaeche f, sonst der groessten ebenen), X = laengste gerade Kante darin."""
    if f is None:
        fl = richtungen(s)
        if not fl: return None
        f = fl[0]
    n = f.normalAt()
    kanten = [e for e in f.Edges() if e.geomType() == "LINE"]
    if kanten:
        e = max(kanten, key=lambda e: e.Length()); xd = e.endPoint() - e.startPoint()
    else:
        xd = cq.Vector(1, 0, 0) if abs(n.x) < 0.9 else cq.Vector(0, 1, 0)
    xd = xd - n * xd.dot(n)
    if xd.Length < 1e-6: xd = cq.Vector(1, 0, 0) - n * n.x
    pl = cq.Plane(origin=s.Center(), xDir=xd.normalized(), normal=n)
    return pl.toLocalCoords(s)

def ist_25d(s):
    """Waagrechte Ebenen (z) und Abweichungen von reinem 2,5D (Schraegen, quer liegende Rundungen/Bohrungen)."""
    ebenen = []; abw = []
    for f in s.Faces():
        g = f.geomType()
        if g == "PLANE":
            n = f.normalAt()
            if abs(abs(n.z) - 1) < 1e-4: ebenen.append(round(f.Center().z, 3)); continue
            if abs(n.z) < 1e-4: continue
            abw.append("Schraege/Fase"); continue
        if g == "CYLINDER":
            try:
                ax = f._geomAdaptor().Cylinder().Axis().Direction()
                if abs(abs(ax.Z()) - 1) < 1e-4: continue
            except Exception: pass
            abw.append("Rundung/Bohrung quer zur Dicke"); continue
        try:                                   # z. B. BSPLINE, die in Wahrheit eben oder senkrecht ist
            pts = [f.positionAt((i + .5) / 5, (j + .5) / 5) for i in range(5) for j in range(5)]
            nz = [abs(f.normalAt(p).z) for p in pts]
            if all(x > 1 - 1e-4 for x in nz): ebenen.append(round(f.Center().z, 3)); continue
            if all(x < 1e-3 for x in nz): continue
        except Exception: pass
        abw.append(f"Flaeche {g}")
    return sorted(set(ebenen)), sorted(set(abw))

def schnitt(s, zz):
    """Querschnitt in Hoehe zz: Schnittkanten mit der Ebene, zu Flaechen zusammengesetzt (gerade/ungerade = Material/Loch)."""
    import shapely
    from OCP.BRepAlgoAPI import BRepAlgoAPI_Section
    from OCP.gp import gp_Pln, gp_Pnt, gp_Dir
    from shapely.geometry import LineString
    from shapely.ops import polygonize
    sec = BRepAlgoAPI_Section(s.wrapped, gp_Pln(gp_Pnt(0, 0, zz), gp_Dir(0, 0, 1)), True)
    linien = []
    for e in cq.Shape.cast(sec.Shape()).Edges():
        n = 1 if e.geomType() == "LINE" else max(8, int(e.Length() / 0.3))
        linien.append(LineString([(p.x, p.y) for p in (e.positionAt(k / n) for k in range(n + 1))]))
    if not linien: return Polygon()
    r = Polygon()
    ringe = [Polygon(q.exterior) for q in polygonize(shapely.set_precision(unary_union(linien), 0.001))]
    for q in sorted(ringe, key=lambda q: -q.area):          # gerade/ungerade: jede Schnittkurve wechselt Material/Luft
        r = r.symmetric_difference(q)
    return r.buffer(0.005).buffer(-0.005)

def polys(g):
    if g.is_empty: return []
    return [p for p in (g.geoms if hasattr(g, "geoms") else [g]) if p.geom_type == "Polygon" and p.area > 0.05]

def kreis(ring):
    p = Polygon(ring); c = p.centroid
    r = [math.hypot(x - c.x, y - c.y) for x, y in ring.coords]
    rm = sum(r) / len(r)
    if rm > 0 and (max(r) - min(r)) < max(0.02, rm * 0.02) and len(r) > 12 and abs(p.area - math.pi * rm * rm) < 0.02 * math.pi * rm * rm:
        return (c.x, c.y, rm)
    return None

def zerlege(s):
    lv, abw = ist_25d(s)
    bb = s.BoundingBox(); z1, z0 = bb.zmax, bb.zmin
    if abw:                                    # nicht rein 2,5D: zusaetzlich fein in Scheiben (0,5 mm)
        n = min(8, max(2, int((z1 - z0) / 0.5))); lv = lv + [z0 + (z1 - z0) * k / n for k in range(1, n)]
    lv = sorted(set([round(z0, 3), round(z1, 3)] + [round(v, 3) for v in lv if z0 - 1e-3 <= v <= z1 + 1e-3]), reverse=True)
    lv2 = [lv[0]]
    for v in lv[1:]:
        if lv2[-1] - v > 0.02: lv2.append(v)
    scheiben = [(a, b, schnitt(s, (a + b) / 2)) for a, b in zip(lv2, lv2[1:])]      # oben -> unten
    O = unary_union([R for _, _, R in scheiben])
    if O.is_empty: return None, "kein Querschnitt"
    O = rein(unary_union([Polygon(p.exterior) for p in polys(rein(O))]))   # Aussenumriss, Loecher gefuellt
    weg = [rein(O.difference(rein(R))) for _, _, R in scheiben]
    durch = weg[0]
    for w in weg[1:]: durch = rein(durch.intersection(w))
    oben = []; acc = None
    for k, w in enumerate(weg):
        acc = w if acc is None else rein(acc.intersection(w)); oben.append(acc)
    unten = [None] * len(weg); acc = None
    for k in range(len(weg) - 1, -1, -1):
        acc = weg[k] if acc is None else rein(acc.intersection(weg[k])); unten[k] = acc
    hinw = list(abw)
    if not abw:
        for k, w in enumerate(weg):
            if rein(w.difference(rein(oben[k].union(unten[k])))).area > 0.5: hinw.append("Hinterschnitt"); break
    if hinw:                                   # nur Umriss + Durchbrueche, Rest von Hand
        return {"T": round(z1 - z0, 2), "O": O, "durch": durch, "oben": [], "unten": [], "hinweis": ", ".join(hinw)}, None
    T = z1 - z0
    taschen_o = []; taschen_u = []
    for k in range(len(weg)):
        tief = z1 - scheiben[k][1]
        nur = rein(rein(oben[k].difference(oben[k + 1] if k + 1 < len(weg) else durch)).difference(durch))
        if nur.area > 0.05 and tief < T - 0.02: taschen_o.append((round(tief, 2), nur))
    for k in range(len(weg) - 1, -1, -1):
        tief = scheiben[k][0] - z0
        nur = rein(rein(rein(unten[k].difference(unten[k - 1] if k > 0 else durch)).difference(durch)).difference(oben[k]))
        if nur.area > 0.05 and tief < T - 0.02: taschen_u.append((round(tief, 2), nur))
    return {"T": round(T, 2), "O": O, "durch": durch, "oben": taschen_o, "unten": taschen_u, "hinweis": ""}, None

def ausrichten(d):
    """Flach legen: X liegt schon entlang der laengsten geraden Kante (lage); laengere Seite waagrecht, links unten 0/0."""
    b = d["O"].bounds
    w = 90.0 if (b[3] - b[1]) > (b[2] - b[0]) + 0.01 else 0.0
    f = lambda g: affinity.rotate(g, w, origin=(0, 0))
    O = f(d["O"]); x0, y0 = O.bounds[0], O.bounds[1]
    g = lambda q: rein(affinity.translate(f(q), -x0, -y0))
    return {"T": d["T"], "O": g(d["O"]), "durch": g(d["durch"]), "oben": [(t, g(q)) for t, q in d["oben"]],
            "unten": [(t, g(q)) for t, q in d["unten"]], "hinweis": d["hinweis"]}

def gleich(a, b):
    """Gleiches Teil? Dicke gleich, Umriss deckt sich (auch um 180 Grad gedreht; ohne Taschen auch gewendet)."""
    if abs(a["T"] - b["T"]) > 0.05 or a["hinweis"] != b["hinweis"]: return False
    if abs(a["O"].area - b["O"].area) > 1e-3 * a["O"].area + 0.5: return False
    if abs(a["durch"].area - b["durch"].area) > 0.3: return False
    if len(a["oben"]) != len(b["oben"]) or len(a["unten"]) != len(b["unten"]): return False
    ba = a["O"].bounds; bb = b["O"].bounds
    if abs((ba[2] - ba[0]) - (bb[2] - bb[0])) > 0.2 or abs((ba[3] - ba[1]) - (bb[3] - bb[1])) > 0.2: return False
    W, H = bb[2], bb[3]
    lagen = [lambda g: g, lambda g: affinity.rotate(g, 180, origin=(W / 2, H / 2))]
    if not (a["oben"] or a["unten"]):
        lagen += [lambda g: affinity.scale(g, -1, 1, origin=(W / 2, H / 2)), lambda g: affinity.scale(g, 1, -1, origin=(W / 2, H / 2))]
    def deckt(x, y): return x.symmetric_difference(y).area < 1e-4 * max(x.area, 1.0) + 0.3
    for L in lagen:
        if not deckt(a["O"], L(b["O"])) or not deckt(a["durch"], L(b["durch"])): continue
        if all(abs(t1 - t2) < 0.05 and deckt(q1, L(q2)) for (t1, q1), (t2, q2) in zip(sorted(a["oben"], key=lambda e: e[0]), sorted(b["oben"], key=lambda e: e[0]))):
            return True
    return False

def rein(g):
    """Nur Flaechen, auf 0,001 mm gerundet (macht die Booleschen Operationen robust)."""
    import shapely
    if g is None or g.is_empty: return Polygon()
    g = shapely.set_precision(g, 0.001)
    p = [q for q in (g.geoms if hasattr(g, "geoms") else [g]) if q.geom_type in ("Polygon", "MultiPolygon")]
    return unary_union(p).buffer(0) if p else Polygon()

def zeichne(msp, g, lay, kreise=True):
    for p in polys(g):
        msp.add_lwpolyline(list(p.exterior.coords), close=True, dxfattribs={"layer": lay})
        for i in p.interiors:
            msp.add_lwpolyline(list(i.coords), close=True, dxfattribs={"layer": lay})

def schreibe(d, name, anzahl, ziel, unterseite=False):
    doc = ezdxf.new(setup=True); doc.units = ezdxf.units.MM; msp = doc.modelspace()
    spiegel = (lambda g: affinity.scale(g, -1, 1, origin=(0, 0))) if unterseite else (lambda g: g)
    O = spiegel(d["O"]); bx = O.bounds
    doc.layers.add("Kontur_aussen_durch", color=7); doc.layers.add("Durchbruch", color=5); doc.layers.add("Text", color=2)
    for p in polys(O): msp.add_lwpolyline(list(p.exterior.coords), close=True, dxfattribs={"layer": "Kontur_aussen_durch"})
    for p in polys(spiegel(d["durch"])):
        k = kreis(p.exterior)
        if k:
            lay = f"Bohrung_D{2 * k[2]:.1f}_durch".replace(".", ",")
            if lay not in doc.layers: doc.layers.add(lay, color=1)
            msp.add_circle((k[0], k[1]), k[2], dxfattribs={"layer": lay})
        else:
            msp.add_lwpolyline(list(p.exterior.coords), close=True, dxfattribs={"layer": "Durchbruch"})
    seite = "unten" if unterseite else "oben"
    for i, (t, q) in enumerate(sorted(d[seite])):
        lay = f"Tasche_{seite}_T{t:.2f}".replace(".", ",")
        if lay not in doc.layers: doc.layers.add(lay, color=3 + i % 4)
        zeichne(msp, spiegel(q), lay)
    txt = f"{name}: {anzahl} Stueck, Dicke {d['T']:g} mm, {bx[2] - bx[0]:.1f} x {bx[3] - bx[1]:.1f} mm"
    if d["hinweis"]: txt += f" - NUR UMRISS, von Hand: {d['hinweis']}"
    if unterseite: txt += " - UNTERSEITE (gespiegelt, Teil umgedreht aufspannen)"
    elif d["unten"]: txt += " - Oberseite, Unterseite siehe eigene Datei"
    msp.add_text(txt, height=2.5, dxfattribs={"layer": "Text"}).set_placement((bx[0], bx[1] - 6))
    doc.saveas(ziel)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("step"); ap.add_argument("ordner"); ap.add_argument("--teil", default="bass")
    a = ap.parse_args()
    os.makedirs(a.ordner, exist_ok=True)
    asm = cq.Assembly.importStep(a.step)
    gruppen = collections.OrderedDict(); nicht = []; weg = collections.Counter()
    for pfad, s in blaetter(asm):
        n = kurzname(pfad)
        if os.environ.get("PLATTEN_LOG"): print("  ", n, flush=True)
        grund = ausgeschlossen(pfad, a.teil)
        if grund: weg[grund] += 1; continue
        d, warum = None, "keine ebene Flaeche"
        rr = richtungen(s); bb = s.BoundingBox(); dm = sorted([bb.xlen, bb.ylen, bb.zlen])
        if dm[0] < 0.3 * dm[1]: rr = rr[:1]            # eindeutig eine Platte: nur die Hauptrichtung
        for i, f in enumerate(rr):                      # sonst beste Aufspannung: rein 2D, einseitig, wenig Stufen
            try:
                e, w = zerlege(lage(s, f))
            except Exception as x:
                e, w = None, f"Fehler {type(x).__name__}"
            if e is None:
                if d is None: warum = w
                continue
            wert = (bool(e["hinweis"]), bool(e["unten"]), e["T"], len(e["oben"]) + len(e["unten"]))
            if d is None or wert < d[0]: d = (wert, e)
            if not e["hinweis"] and not e["unten"]: break
        d = d[1] if d else None
        if d is None: nicht.append((n, warum)); continue
        d = ausrichten(d)
        for g in gruppen.values():
            if gleich(g["d"], d): g["namen"].append(n); break
        else: gruppen[len(gruppen)] = {"d": d, "namen": [n]}
    zeilen = ["Datei;Stueck;Dicke_mm;Laenge_mm;Breite_mm;Seiten;Hinweis;Teile"]
    for g in gruppen.values():
        d = g["d"]; n = len(g["namen"]); name = re.sub(r"(_\d+)+$", "", g["namen"][0])
        ordn = os.path.join(a.ordner, f"Platten_{d['T']:g}mm".replace(".", ","))
        os.makedirs(ordn, exist_ok=True)
        datei = f"{name}_{n}Stk.dxf"; i = 2
        while os.path.exists(os.path.join(ordn, datei)): datei = f"{name}_{i}_{n}Stk.dxf"; i += 1
        schreibe(d, name, n, os.path.join(ordn, datei))
        if d["unten"]: schreibe(d, name, n, os.path.join(ordn, datei.replace(".dxf", "_Unterseite.dxf")), True)
        b = d["O"].bounds
        zeilen.append(f"{os.path.basename(ordn)}/{datei};{n};{d['T']:g};{b[2] - b[0]:.1f};{b[3] - b[1]:.1f};"
                      f"{'beidseitig' if d['unten'] else 'oben'};{d['hinweis']};{' '.join(g['namen'])}")
    open(os.path.join(a.ordner, "Uebersicht_2D_Teile.csv"), "w").write("\n".join(zeilen) + "\n")
    open(os.path.join(a.ordner, "nicht_2D.txt"), "w").write(
        "Teile, die nicht als 2D (2,5D) gefraest werden koennen - eigene Bearbeitung noetig:\n" +
        "".join(f"{n}: {w}\n" for n, w in nicht))
    print(f"2D-Fraesteile: {len(gruppen)} Dateien ({sum(len(g['namen']) for g in gruppen.values())} Teile), "
          f"nicht 2D: {len(nicht)}, ausgelassen: {dict(weg)}")

if __name__ == "__main__":
    main()
