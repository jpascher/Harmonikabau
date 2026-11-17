# Diskantmechanik Markus nach griffstock-alternativ1-schnitt.dxf (eingebunden von diskant.py per exec)
# Schnitt -> Modell: u = Abstand nach aussen von der Bodenaussenseite (Z = -u), Y wie Modell; X aus der Griffstock-Aufsicht.
# Hebel: Kontur aus dem Schnitt (Alu 2 mm), langer Arm bis zur Kroepfung, dort seitlich zur Klappenmitte versetzt,
# Endstueck mit Stufe auf die Klappe (Kontur aus dem Schnitt). Holzteile aus den Holzrand-Zeichnungen, auf den Hebel
# im Schnitt eingepasst (Abweichung < 0,3 mm); kein Alu ueber das Holz hinaus. Feste Achse Ø5 x 11 im Lagerpunkt.
exec(open("hebel_reg.py").read().split('\nif __name__')[0])

_tmpdoc = ezdxf.new()
def _flach(vs):
    pl = _tmpdoc.modelspace().add_lwpolyline([(x_, y_, 0, 0, b_) for x_, y_, b_ in vs], format="xyseb")
    q = [(v.x, v.y) for v in ezdxf.path.make_path(pl).flattening(0.05)]
    out = []
    for p_ in q:
        if not out or math.hypot(p_[0]-out[-1][0], p_[1]-out[-1][1]) > 1e-3: out.append(p_)
    return out
def _umkehr(v):
    # Polylinie umdrehen: Bulge gehoert zum Segment ab dem Punkt -> verschieben und Vorzeichen wechseln
    r_ = list(reversed(v)); out = []
    for i_ in range(len(r_)):
        b_ = -r_[i_ + 1][2] if i_ + 1 < len(r_) else 0.0
        out.append([r_[i_][0], r_[i_][1], b_])
    return out
def holz_ring(r):
    """Holzteil-Umriss ohne die Fraes-Einlaeufe: Bogen-Enden der beiden offenen Linien weglassen, gerade schliessen."""
    f_, lay_ = HEBEL_DXF[r]
    m_ = ezdxf.readfile(f_).modelspace()
    alu_e = max((e for e in m_ if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == lay_ and len(e) > 6), key=lambda e: bbox.extents([e]).size.y)
    linien = []
    for e in m_:
        if e.dxftype() != "LWPOLYLINE" or e is alu_e or e.closed: continue
        vs = []
        for x_, y_, s_, e_, b_ in e.get_points():
            if not vs or math.hypot(x_ - vs[-1][0], y_ - vs[-1][1]) > 1e-3: vs.append([x_, y_, b_])
            else: vs[-1][2] = b_
        if len(vs) < 3: continue
        linien.append(vs)
    if len(linien) < 2: return None
    from shapely.geometry import LineString as _LS, Point as _P
    from shapely.ops import polygonize as _pz, unary_union as _uu
    A, B = linien[0], linien[1]
    def _ende(v, anfang): return (v[0][0], v[0][1]) if anfang else (v[-1][0], v[-1][1])
    # Paarung der Enden: A-Ende mit naechstem B-Ende
    if math.dist(_ende(A, False), _ende(B, True)) > math.dist(_ende(A, False), _ende(B, False)):
        B = [[x_, y_, 0.0] for x_, y_, b_ in reversed(B)] if all(abs(v[2]) < 0.01 for v in B) else _umkehr(B)
    la, lb = _LS(_flach(A)), _LS(_flach(B))
    ip = la.intersection(lb)
    ipts = [ip] if ip.geom_type == "Point" else list(getattr(ip, "geoms", []))
    def _kreuzt(pa, pb):
        return any(min(math.dist((q.x, q.y), pa), math.dist((q.x, q.y), pb)) < 8.0 for q in ipts if q.geom_type == "Point")
    def _ohne_einlauf_anfang(v):
        for i_ in range(min(2, len(v) - 2)):
            if abs(v[i_][2]) > 0.1: return v[i_ + 1:]
        return v
    def _ohne_einlauf_ende(v):
        n_ = len(v)
        for i_ in range(n_ - 2, max(n_ - 4, 0), -1):
            if abs(v[i_][2]) > 0.1:
                w = v[:i_ + 1]; w[-1] = [w[-1][0], w[-1][1], 0.0]; return w
        return v
    zu = []
    # Ende A-Ende / B-Anfang
    if not _kreuzt(_ende(A, False), _ende(B, True)):
        A = _ohne_einlauf_ende(A); B = _ohne_einlauf_anfang(B); zu.append((_ende(A, False), _ende(B, True)))
    if not _kreuzt(_ende(A, True), _ende(B, False)):
        A = _ohne_einlauf_anfang(A); B = _ohne_einlauf_ende(B); zu.append((_ende(A, True), _ende(B, False)))
    netz = _uu([_LS(_flach(A)), _LS(_flach(B))] + [_LS(list(z)) for z in zu])
    polys = list(_pz(netz))
    if not polys: return None
    pg = max(polys, key=lambda g: g.area)
    # Holzteil umschliesst das Alu in seinem Bereich vollstaendig (Alu mittig im Kanal, nicht sichtbar)
    from shapely.geometry import Polygon as _PG, box as _bx
    alu_pg = _PG(hebel_profil(r)[0]).buffer(0)
    x0_, y0_, x1_, y1_ = pg.bounds
    if r == 3: pg = pg.union(alu_pg.intersection(_bx(x0_, y0_, x1_, y1_)).buffer(0.3).buffer(-0.3))   # nur Reihe 3: Holz zu knapp gezeichnet
    if pg.geom_type != "Polygon": pg = max(pg.geoms, key=lambda g: g.area)
    pg = pg.simplify(0.02)
    out = []
    for x_, y_ in list(pg.exterior.coords)[:-1]:
        if not out or math.hypot(x_ - out[-1][0], y_ - out[-1][1]) > 0.05: out.append((float(x_), float(y_)))
    if math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) < 0.05: out.pop()
    return out

def prisma_uy(pts_uy, x0, dicke):
    return (cq.Workplane("YZ").polyline([(Y, -u) for u, Y in pts_uy]).close()
            .extrude(dicke).translate((x0 - dicke / 2, 0, 0)).val())
# --- Klappenmitten aus dem Boden ---
klappen_m, _ges = [], []
for ins in bm.query('INSERT[layer=="2Füllung-Öffnug"]'):
    p0 = (round(ins.dxf.insert.x, 2), round(ins.dxf.insert.y, 2))
    if p0 in _ges: continue
    _ges.append(p0)
    cs = [(v.dxf.center.x + BXO, v.dxf.center.y) for v in ins.virtual_entities() if v.dxftype() == "CIRCLE"]
    klappen_m.append((sum(c[0] for c in cs) / len(cs), sum(c[1] for c in cs) / len(cs)))
# --- Hebelplaetze aus der Griffstock-Aufsicht ---
gm = ezdxf.readfile("gs.dxf").modelspace()
from collections import defaultdict as _dd
_lr = _dd(list)
for e in gm:
    if e.dxf.layer != "Lagerbock_aufsicht": continue
    b = bbox.extents([e])
    if b.has_data and abs(b.size.x - 5.0) < 0.2 and (abs(b.size.y - 14.2) < 0.3 or abs(b.size.y - 11.3) < 0.3):
        _lr[round(b.center.y, 0)].append(b.center.x + BXO)
_ly = sorted(_lr.keys(), reverse=True)                 # 14er-Reihe = aeussere Knopfreihe = Reihe 1
HEBEL_X = {i + 1: sorted(_lr[y]) for i, y in enumerate(_ly)}
# --- feste Zuordnung Knopfreihe r -> Stimmstockreihe r, eine Des-Klappe an der aeusseren Reihe ---
kl_r = {r: sorted([k for k in klappen_m if _reihe(k[1]) == r], key=lambda k: k[0]) for r in (1, 2, 3, 4)}
frei_r = {r: list(HEBEL_X[r]) for r in (1, 2, 3, 4)}
zuordnung = []
for r in (1, 2, 3):
    for k in kl_r[r]:
        x = min(frei_r[r], key=lambda h: abs(h - k[0])); frei_r[r].remove(x); zuordnung.append((r, x, k, r))
if frei_r[1]:
    xf = frei_r[1][0]; kx = min(kl_r[4], key=lambda k: abs(k[0] - xf)); kl_r[4].remove(kx); frei_r[1].remove(xf)
    zuordnung.append((1, xf, kx, 4)); print(f"Ausnahme: Des-Klappe X {kx[0]:.1f} an aeusserer Knopfreihe (Hebel X {xf:.1f})")
for k in kl_r[4]:
    x = min(frei_r[4], key=lambda h: abs(h - k[0])); frei_r[4].remove(x); zuordnung.append((4, x, k, 4))

T_ALU, T_HOLZ = 2.0, 4.8          # Holzteil fuellt die 5-mm-Nut mit minimalem Spiel, Alu mittig im Kanal
U_KL = T_POLSTER + T_KLAPPE                            # Klappenoberseite (u)
def schnitt_y(pts, y0, y1):
    q = [p for p in pts if y0 <= p[1] <= y1]
    return (min(p[0] for p in q), max(p[0] for p in q)) if q else None
knopf_u = sorted(KNOPF, key=lambda q: -q[0])          # aeusserster Knopf = Reihe 1
DREHPUNKTE = {}                                       # Gruppe -> (X, Y, Z) Lagerpunkt, Drehachse = X
for n_, (r, x, k, kr) in enumerate(zuordnung):
    up, Yp = lager_von(r)
    G = f"Hebel_R{r}_{n_+1:02d}"
    DREHPUNKTE[G] = (x, Yp, -up)
    alu_s = HEBEL[r]
    end = ENDE[kr]
    y_kr = min(p[1] for p in end)                       # Beginn des gekroepften Endstuecks
    # Arm in der Hebelebene bis zur Kroepfung
    hb = prisma_uy(alu_s, x, T_ALU).cut(box(x - 5, x + 5, y_kr + 3.5, 400, -300, 50))
    # kein Alu ueber das Holzteil hinaus
    ring = holz_ring(r)
    if ring:
        holz_m = trafo(r, ring)
        xw0 = min(q[0] for q in ring)
        zone_m = trafo(r, [(xw0, -200), (400, -200), (400, 200), (xw0, 200)])
        zone = prisma_uy(zone_m, x, T_ALU + 2); holzp = prisma_uy(holz_m, x, T_ALU + 2)
        hb = hb.cut(zone).fuse(hb.intersect(zone).intersect(holzp))
        teile.append((f"{G}/Holzteil", prisma_uy(holz_m, x, T_HOLZ), "Hebelholz"))
    # Hoehe des Arms an der Kroepfung
    ua = schnitt_y(alu_s, y_kr, y_kr + 3.5) or schnitt_y(end, y_kr, y_kr + 3.5)
    _c = (ua[0] + ua[1]) / 2; _h = max((ua[1] - ua[0]) / 2, 1.0); ua = (_c - _h, _c + _h)
    if kr == r:
        ende = prisma_uy(end, k[0], T_ALU)
    else:   # Ausnahme: Arm der aeusseren Reihe zur Des-Klappe -> gerade Stufe hinunter auf die Klappe
        ende = box(k[0] - T_ALU/2, k[0] + T_ALU/2, y_kr, k[1], -ua[1], -ua[0])
        ende = ende.fuse(box(k[0] - T_ALU/2, k[0] + T_ALU/2, k[1] - 3, k[1], -ua[1], -U_KL))
        ende = ende.fuse(box(k[0] - 4, k[0] + 4, k[1] - 6, k[1] + 2, -U_KL - T_ALU, -U_KL))
    xa, xb = sorted((x, k[0]))
    bruecke = box(xa - T_ALU/2, xb + T_ALU/2, y_kr, y_kr + 3.5, -ua[1], -ua[0])     # seitliche Kroepfung
    teile.append((f"{G}/Hebel_Alu", hb.fuse(bruecke).fuse(ende), "Alu"))
    km = min(KL_TEILE, key=lambda t_: math.hypot(t_[0][0] - k[0], t_[0][1] - k[1]))
    teile.append((f"{G}/Klappe", km[1], "Klappe")); teile.append((f"{G}/Klappenpolster", km[2], "Leder"))
    teile.append((f"{G}/Achse", cq.Workplane("YZ").center(Yp, -up).circle(2.5).extrude(11.0).translate((x - 5.5, 0, 0)).val(), "Stahl"))
    # Knopf (Ø15,5, Kappe 3 mm) auf der Knopfreihe r, Achse wie im Schnitt
    ku, kY = knopf_u[r - 1]
    d_ = KNOPF_RICHT; L_ = math.hypot(*d_)
    knopf = cq.Solid.makeCylinder(15.5 / 2, 3.0, cq.Vector(x, kY, -ku), cq.Vector(0, d_[1] / L_, -d_[0] / L_))
    teile.append((f"{G}/Knopf", knopf, "Knopf"))
print("Hebel:", len(zuordnung))
# --- Griffstock aus dem Schnitt: Koerper, Kunststoffteil (schraeg, Lagerpunkte), Dreieckleiste ---
GS_X0, GS_X1 = 0.0, 380.0
gk = prisma_uy(GRIFFSTOCK, (GS_X0 + GS_X1) / 2, GS_X1 - GS_X0)
kst = prisma_uy(KUNSTSTOFF, 190.0, 290.0)
dre = prisma_uy(DREIECK, 190.0, 360.0)
schl_g, schl_k = [], []
for r, x, k, kr in zuordnung:
    up, Yp = lager_von(r)
    schl_k.append(prisma_uy(HEBEL[r], x, 3.0))
    schl_k.append(cq.Workplane("YZ").center(Yp, -up).circle(2.6).extrude(12.0).translate((x - 6.0, 0, 0)).val())
    schl_g.append(prisma_uy(HEBEL[r], x, 3.0))
    ring = holz_ring(r)
    if ring: schl_g.append(prisma_uy(trafo(r, ring), x, 5.0))
# Seitenform links/rechts nach der Griffstock-Aufsicht (Diskantteil-top): Stirnseiten zum Knopfende hin eingezogen
_m2 = ezdxf.readfile("gs2.dxf").modelspace()
_u_max = max(q[0] for q in GRIFFSTOCK)
def _seite(e):
    q = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.1)]
    y_hi = max(p[1] for p in q); y_lo = min(p[1] for p in q)
    return [(p[0], (y_hi - p[1]) / (y_hi - y_lo) * (_u_max + 2)) for p in q]
_kont = [e for e in _m2.query("LWPOLYLINE") if e.dxf.layer == "Diskantteil-top" and len(e) == 5]
_r = max(_kont, key=lambda e: max(v[0] for v in e.get_points()))      # rechte Aussenkontur (x bis 0)
_l = min(_kont, key=lambda e: min(v[0] for v in e.get_points()))      # linke Aussenkontur (x bis -380)
def _schnitt_kontur(e, rechts):
    q = _seite(e)
    pts = [((380.0 + x_) if rechts else (x_ + 380.0), -u_) for x_, u_ in q]   # (X, Z)
    aussen = 420.0 if rechts else -40.0
    poly = pts + [(aussen, pts[-1][1]), (aussen, pts[0][1])]
    return cq.Workplane("XZ").polyline(poly).close().extrude(-400).translate((0, -100, 0)).val()
seitenschnitt = [_schnitt_kontur(_r, True), _schnitt_kontur(_l, False)]
gk = schneide(gk, seitenschnitt); dre = schneide(dre, seitenschnitt)
teile.append(("Griffstock_Koerper", schneide(gk, schl_g), "Griffstock"))
teile.append(("Kunststoffteil_Lagerblock", schneide(kst, schl_k), "Kunststoff"))
teile.append(("Dreieckleiste", dre, "Griffstock"))
# Verschraubung von unten (innen): durch Verstaerkungsleiste und Boden in Dreieckleiste und Kunststoffteil
_du = sum(q[0] for q in DREIECK) / 3; _dY = sum(q[1] for q in DREIECK) / 3
_ku_max = max(q[0] for q in KUNSTSTOFF)
schr_l = []
for i_, xs_ in enumerate((45.0, 95.0, 145.0, 235.0, 285.0, 335.0)):
    sch = cq.Solid.makeCylinder(1.5, T_BODEN + 5.0 + _du + 10.0, cq.Vector(xs_, _dY, T_BODEN + 5.0), cq.Vector(0, 0, -1))
    kopf = cq.Solid.makeCone(3.0, 1.5, 1.5, cq.Vector(xs_, _dY, T_BODEN + 5.0), cq.Vector(0, 0, -1))
    teile.append((f"Schraube_Griffstock_{i_+1}", sch.fuse(kopf), "Schraube"))
    schr_l.append(cq.Solid.makeCylinder(1.6, 200, cq.Vector(xs_, _dY, 50), cq.Vector(0, 0, -1)))
verst = schneide(verst, schr_l)
teile.append(("Verstaerkungsleiste_innen", verst, "Holz"))
import json as _js
open("drehpunkte.json", "w").write(_js.dumps(DREHPUNKTE))
