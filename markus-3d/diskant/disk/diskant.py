# Diskantteil Markus (ohne Diskantmechanik) aus den CamBam-DXF in 0-Markus/Diskant
# Modell: X = x_boden + 380 (0..380), Y = y_boden (0..205, Y=0 Griffstockseite), Z = 0 Bodenaussenseite,
#         Z > 0 nach innen (Balgseite), Verdeck bei Z < 0.
# Aufruf: python3 diskant.py <variante> <ausgabe.step>     variante = Ordnername unter var/ (z.B. BESASDES)
import sys, os, glob, math, ezdxf, ezdxf.path, cadquery as cq
from ezdxf import bbox

VAR = sys.argv[1] if len(sys.argv) > 1 else "BESASDES"
OUT = sys.argv[2] if len(sys.argv) > 2 else f"Diskantteil_Markus_{VAR}.step"
BXO = 380.0
L, B = 380.0, 205.0
T_BODEN = 4.0             # Diskantboden 4 mm (Angabe Johann)
T_LEDER = 1.0             # Bodenleder auf der Innenseite, gelocht
Z_ST = T_BODEN + T_LEDER  # Auflage der Stimmstoecke
H_RAHMEN, D_RAHMEN = 51.0, 5.5   # (Johann 16.11.: niedriger, Randleiste oben 61 >= Stimmstock 60,5)  # Rahmen wie Bass (Zinken, Verzierungen als Taschen), 382 x 207 aussen, Boden 1 mm eingelassen
RI = D_RAHMEN - 1.0             # Rahmen-Innenkante ab Bodenkante
T_ZARGE = RI
T_SOLE, T_WAND, H_WAND, NEIG = 3.0, 10.0, 60.0, 4.0   # Sohle 3 mm (CamBam -6 nur zum sicheren Durchfraesen)
T_PLATTE, TASCHE = 2.0, 2.0
R_GEH, R_VD, VD_DY = 2.0, 1.5, -3.0   # Kantenradien Gehaeuse/Verdeck, Verdeck 3 mm Richtung Griffstock
T_VD, H_VD = 3.5, 32.0    # Verdeck-Wandstaerke 3,5 mm (Johann 16.11.), Seitenstreifen 31,5

def box(x0, x1, y0, y1, z0, z1):
    return cq.Workplane("XY").box(x1-x0, y1-y0, z1-z0, centered=False).translate((x0, y0, z0)).val()
def vereinigt(liste):
    c = liste[0]
    for s in liste[1:]: c = c.fuse(s)
    return c
def bohrungen(punkte, d, z0, z1):
    return [cq.Workplane("XY").center(x, y).circle(d/2).extrude(z1-z0).translate((0, 0, z0)).val() for x, y in punkte]
def schneide(s, werkzeuge):
    if not werkzeuge: return s
    return s.cut(cq.Compound.makeCompound(werkzeuge))

def kreise_der_inserts(m, layer):
    out = []
    for ins in m.query(f'INSERT[layer=="{layer}"]'):
        for e in ins.virtual_entities():
            if e.dxftype() == "CIRCLE":
                out.append((e.dxf.center.x, e.dxf.center.y, e.dxf.radius))
    return out
def dedupe(pts, tol=0.3):
    r = []
    for p in pts:
        if all(abs(p[0]-q[0]) > tol or abs(p[1]-q[1]) > tol for q in r): r.append(p)
    return r

teile = []   # (name, solid, farbe)
# ---------------- Diskantboden ----------------
bm = ezdxf.readfile("boden.dxf").modelspace()
loecher = dedupe(kreise_der_inserts(bm, "2Füllung-Öffnug"))
boden = box(RI, L - RI, RI, B - RI, 0, T_BODEN)          # Boden liegt innerhalb des Rahmens (zwischen den Rahmenwaenden)
werk = bohrungen([(x + BXO, y) for x, y, r in loecher], 10.5, -1, T_BODEN + 1)
for c in [c for c in bm.query('CIRCLE') if c.dxf.layer in ("GH-füllung-innenabmessungen", "0Griffstockleiste")]:
    werk += bohrungen([(c.dxf.center.x + BXO, c.dxf.center.y)], 2*c.dxf.radius, -1, T_BODEN + 1)
boden = schneide(boden, werk)
teile.append(("Diskantboden", boden, "Holz"))
leder = box(T_ZARGE, L - T_ZARGE, T_ZARGE, B - T_ZARGE, T_BODEN, Z_ST)
leder = schneide(leder, bohrungen([(x + BXO, y) for x, y, r in loecher], 10.5, T_BODEN - 1, Z_ST + 1))
# Verstaerkung innen im Griffstockbereich: 5-mm-Leiste auf dem Boden von der Frontwand bis zu den Stimmstocksohlen
Y_VERST = 47.0
leder = leder.cut(box(0, L, 0, Y_VERST, T_BODEN - 1, Z_ST + 1))
teile.append(("Bodenleder", leder, "Leder"))
verst = box(T_ZARGE, L - T_ZARGE, T_ZARGE, Y_VERST, T_BODEN, T_BODEN + 5.0)
print("Boden: Ventilloecher", len(loecher), "Bohrungen", len(werk) - len(loecher))

# ---------------- Diskantklappen (aussen auf dem Boden, unter dem Verdeck) ----------------
# Umriss je Klappe aus dem Bodenblock: "Klappen-rand-aussen" = Klappe, "Klappen-rand-3mm" = Filz/Leder-Polster.
T_POLSTER, T_KLAPPE = 3.4, 4.0     # Polster (Filz+Leder) gleich gross wie die Klappe, liegt am Boden auf; Klappe 4 mm (Schnitt)
n_kl = 0
gesehen = []
KL_TEILE = []        # (Mitte, Klappe, Polster) - mit Mechanik kommen sie in die Hebelgruppe
MIT_MECH = os.environ.get("MECHANIK", "1") == "1"
for ins in bm.query('INSERT[layer=="2Füllung-Öffnug"]'):
    p0 = (round(ins.dxf.insert.x, 2), round(ins.dxf.insert.y, 2))
    if p0 in gesehen: continue                   # doppelte Inserts in der Zeichnung
    gesehen.append(p0)
    for v in ins.virtual_entities():
        if v.dxftype() != "LWPOLYLINE": continue
        pts = [(q.x + BXO, q.y) for q in ezdxf.path.make_path(v).flattening(0.05)]
        pts = [q for i_, q in enumerate(pts) if i_ == 0 or math.hypot(q[0]-pts[i_-1][0], q[1]-pts[i_-1][1]) > 1e-4]
        if math.hypot(pts[0][0]-pts[-1][0], pts[0][1]-pts[-1][1]) < 1e-4: pts = pts[:-1]
        if v.dxf.layer == "Klappen-rand-aussen":
            kl_ = cq.Workplane("XY").polyline(pts).close().extrude(T_KLAPPE).translate((0, 0, -T_POLSTER - T_KLAPPE)).val()
            po_ = cq.Workplane("XY").polyline(pts).close().extrude(T_POLSTER).translate((0, 0, -T_POLSTER)).val()
            mi_ = (sum(q[0] for q in pts) / len(pts), sum(q[1] for q in pts) / len(pts))
            KL_TEILE.append((mi_, kl_, po_))
            if not MIT_MECH:
                teile.append((f"Klappe_{n_kl+1}", kl_, "Klappe")); teile.append((f"Klappenpolster_{n_kl+1}", po_, "Leder"))
    n_kl += 1
print("Klappen:", n_kl)

# ---------------- Rahmen (aus Dr_front/rueck/seite.step, erzeugt mit rahmen_d.py wie die Bass-Rahmenteile) ----------------
imp = lambda f: cq.importers.importStep(f).val()
O_, XA, ZA = cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), cq.Vector(0, 0, 1)
st = lambda s_: s_.rotate(O_, XA, 90)
LR, BR = L + 2.0, B + 2.0
rah = [("Rahmen_hinten", st(imp("Dr_rueck_51.step")).translate(cq.Vector(0, D_RAHMEN, 0))),
       ("Rahmen_vorne", st(imp("Dr_front_51.step")).rotate(O_, ZA, 180).translate(cq.Vector(LR, BR - D_RAHMEN, 0))),
       ("Rahmen_links", st(imp("Dr_seite_51.step")).rotate(O_, ZA, -90).translate(cq.Vector(D_RAHMEN, BR, 0))),
       ("Rahmen_rechts", st(imp("Dr_seite_51.step")).rotate(O_, ZA, 90).translate(cq.Vector(LR - D_RAHMEN, 0, 0)))]
falz = box(1.0, LR - 1.0, 1.0, BR - 1.0, -1, T_BODEN)          # Falz fuer den eingelassenen Boden
for n_, s_ in rah:
    s_ = s_.rotate(cq.Vector(LR/2, BR/2, 0), cq.Vector(LR/2, BR/2, 1), 180)   # Front zur Griffstockseite (Y = 0)
    s_ = s_.translate(cq.Vector(-1, -1, 0))
    # Aussenkanten R 0,5 wie das Verdeck: senkrechte Gehaeuseecken und untere Aussenkanten (Bodenseite),
    # nicht zum Balg hin und nicht die Kante, die mit dem Griffstock buendig ist (Y = -1)
    _x0, _x1, _y0, _y1 = -1.0, LR - 1.0, -1.0, BR - 1.0
    def _rund(e):
        p0, p1 = e.startPoint(), e.endPoint()
        ecke = lambda p: (abs(p.x - _x0) < 0.05 or abs(p.x - _x1) < 0.05) and (abs(p.y - _y0) < 0.05 or abs(p.y - _y1) < 0.05)
        if ecke(p0) and ecke(p1) and abs(p0.x - p1.x) < 0.05 and abs(p0.y - p1.y) < 0.05: return True
        if abs(p0.z) < 0.05 and abs(p1.z) < 0.05:
            return (abs(p0.x - p1.x) < 0.05 and abs(p0.x - _x0) < 0.05 or abs(p0.x - p1.x) < 0.05 and abs(p0.x - _x1) < 0.05
                    or abs(p0.y - p1.y) < 0.05 and abs(p0.y - _y1) < 0.05)
        return False
    _ed = [e for e in s_.Edges() if e.geomType() == "LINE" and _rund(e)]
    for _r in (R_GEH, 1.5, 1.0, 0.5):
        try:
            if _ed: s_ = s_.fillet(_r, _ed)
            print(n_, f"Aussenkanten R{_r}:", len(_ed)); break
        except Exception as ex_:
            print(n_, f"Abrundung R{_r} nicht moeglich:", ex_)
    s_ = s_.mirror("YZ", cq.Vector(L / 2, 0, 0))          # seitenverkehrt wie der ganze Diskantteil
    teile.append((n_.replace("hinten", "X").replace("vorne", "hinten").replace("X", "vorne"), s_, "Holz"))
# Umlaufende Leiste innen zum Balg: 3 mm dick, 5 mm im Rahmen, 10 mm vorstehend, Ecken auf Gehrung,
# untere Innenkante abgeschraegt (2 mm), uebrige Kanten leicht gebrochen (0,3 mm)
T_LE, Z_LE0, Z_LE1 = 3.0, H_RAHMEN - 5.0, H_RAHMEN + 10.0
_xi0, _xi1, _yi0, _yi1 = RI, L - RI, RI, B - RI
_leisten = {"Rahmenleiste_vorne": [(_xi0, _yi0), (_xi1, _yi0), (_xi1 - T_LE, _yi0 + T_LE), (_xi0 + T_LE, _yi0 + T_LE)],
            "Rahmenleiste_hinten": [(_xi0, _yi1), (_xi0 + T_LE, _yi1 - T_LE), (_xi1 - T_LE, _yi1 - T_LE), (_xi1, _yi1)],
            "Rahmenleiste_links": [(_xi0, _yi0), (_xi0 + T_LE, _yi0 + T_LE), (_xi0 + T_LE, _yi1 - T_LE), (_xi0, _yi1)],
            "Rahmenleiste_rechts": [(_xi1, _yi0), (_xi1, _yi1), (_xi1 - T_LE, _yi1 - T_LE), (_xi1 - T_LE, _yi0 + T_LE)]}
for n_, pts_ in _leisten.items():
    le = cq.Workplane("XY").polyline(pts_).close().extrude(Z_LE1 - Z_LE0).translate((0, 0, Z_LE0)).val()
    try:
        cx, cy = (_xi0 + _xi1) / 2, (_yi0 + _yi1) / 2
        innen = [e for e in le.Edges() if abs(e.Center().z - Z_LE0) < 0.01 and
                 math.hypot(e.Center().x - cx, e.Center().y - cy) < max(math.hypot(q[0] - cx, q[1] - cy) for q in pts_) - 1.0
                 and (abs(e.startPoint().x - e.endPoint().x) < 0.01 or abs(e.startPoint().y - e.endPoint().y) < 0.01)]
        innen = [e for e in innen if not (abs(e.Center().x - _xi0) < 0.01 or abs(e.Center().x - _xi1) < 0.01 or abs(e.Center().y - _yi0) < 0.01 or abs(e.Center().y - _yi1) < 0.01)]
        le = le.chamfer(2.0, None, innen)
        oben = [e for e in le.Edges() if abs(e.Center().z - Z_LE1) < 0.01]
        le = le.chamfer(0.3, None, oben)
    except Exception as ex_:
        print(n_, "Fase nicht moeglich:", ex_)
    teile.append((n_, le, "Holz"))

# ---------------- Stimmstoecke ----------------
bl = [(x + BXO, y) for x, y, r in loecher]
def kammern(m):
    """Plattenfelder nach Layer: Vorder* -> vorne, R*ck* -> hinten; Rechtecke ~15 breit, Hoehe = Plattenlaenge."""
    vorn, hint = {}, {}
    for p in m.query("LWPOLYLINE"):
        lay = p.dxf.layer.lower()
        ziel = vorn if lay.startswith("vorder") else hint if (lay.startswith("rück") or lay.startswith("rueck")) else None
        if ziel is None: continue
        pts = [(a, b) for a, b, *_ in p.get_points()]
        if len(pts) < 4: continue
        xs = [a for a, _ in pts]; ys = [b for _, b in pts]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if not (13 < w < 17 and 15 < h < 60): continue
        xm = round((max(xs) + min(xs)) / 2, 1)
        ziel[xm] = max(ziel.get(xm, 0), h)
    return vorn, hint
def sohle(m):
    best = None
    for p in m.query("LWPOLYLINE"):
        if p.dxf.layer != "Sole": continue
        xs = [a for a, b, *_ in p.get_points()]; ys = [b for a, b, *_ in p.get_points()]
        r = (min(xs), max(xs), min(ys), max(ys))
        if (r[1]-r[0]) < 100 or (r[3]-r[2]) < 20: continue
        if best is None or (r[1]-r[0])*(r[3]-r[2]) > (best[1]-best[0])*(best[3]-best[2]): best = r
    return best or (-372.0, 0.0, 0.0, 35.0)     # Sohle nur als offene Linien gezeichnet -> Standard 372 x 35

def passung(sl, tr):
    n = 0
    for x, y in sl:
        X, Y = tr(x, y)
        if any(abs(X - a) < 0.6 and abs(Y - b) < 0.6 for a, b in bl): n += 1
    return n

# Neigung je Stimmstock laut Seitenansicht Diskant in der Mechanik-DXF (Schnitt Aussen): Waende parallel,
# Gesamtbreite 24 (Wand 10 + Spalt 4 + Wand 10), Grad positiv = Oberkante Richtung hinten (+Y).
# Reihenfolge von hinten (1. Reihe) nach vorne (4. Reihe, Griffstockseite).
NEIGUNG = {"1.reihe": math.degrees(math.atan(5/55)), "2.reihe": 0.0,
           "3.reihe": -math.degrees(math.atan(6/55)), "4.reihe": -math.degrees(math.atan(12/55))}
XW0, XW1 = RI, L - RI          # Rahmen-Innenwaende in X
B_LEISTE, T_TW, H_NASE = 10.0, 2.5, 20.0
B_RASTE_U, B_RASTE_O = 6.0, 4.0   # Rastleiste fuer die Nasen an der Wand der tiefen Seite: unten breiter als oben   # Schraubleiste an der hohen Seite, Trennwand, Nasenhoehe ueber Sohle
SEITE_TIEF, Z_ENDE_HOCH, STOCK_INFO = [], [], []
B_STOCK, B_OBEN = 24.0, 19.0     # Stimmstock unten 24, oben 19 breit (Angabe Johann)
dateien = sorted(glob.glob(f"var/{VAR}/*.reihe/*.dxf"))
assert dateien, "Variante nicht gefunden"
holz2 = []
for f in dateien:
    reihe = f.split("/")[-2]; ton = f.split("/")[-1][:-4]
    m = ezdxf.readfile(f).modelspace()
    sl = [(x, y) for x, y, r in dedupe(kreise_der_inserts(m, "2Füllung-Öffnugen-ori"))]
    # beste Lage auf dem Boden suchen: Spiegelung in y, Verschiebung aus Lochdifferenzen
    best = None
    Y_BAND = {'1.reihe': (152.8, 187.8), '2.reihe': (117.8, 152.8), '3.reihe': (82.8, 117.8), '4.reihe': (47.2, 82.2)}
    SX = -1      # Stimmstock-Zeichnung ist fuer das Fraesen aufgeklappt: beim Einbau gewendet -> in X gespiegelt
    for sy in (1, -1):
        for a, b in bl:
            x0, y0 = sl[0]
            dx = a - SX*x0; dy = b - sy*y0
            if not (Y_BAND[reihe][0] - 3 <= sy*y0 + dy <= Y_BAND[reihe][1] + 3): continue
            tr = lambda x, y, dx=dx, dy=dy, sy=sy: (SX*x + dx, sy*y + dy)
            n = passung(sl, tr)
            if best is None or n > best[0]: best = (n, sy, dx, dy)
    n, sy, dx, dy = best
    tr = lambda x, y: (SX*x + dx, sy*y + dy)
    vorn, hint = kammern(m)
    sx0, sx1, sy0, sy1 = sohle(m)
    xa, ya = tr(sx0, sy0); xb, yb = tr(sx1, sy1)
    X0, X1, Y0, Y1 = max(min(xa, xb), 4.0), min(max(xa, xb), 376.0), min(ya, yb), max(ya, yb)
    print(f"{reihe} {ton}: Loecher {n}/{len(sl)}  Sohle X {X0:.1f}-{X1:.1f} Y {Y0:.1f}-{Y1:.1f}  Kammern vorn {len(vorn)} hinten {len(hint)}")
    zs = Z_ST
    # Laengsaufteilung: tiefe Seite (lange Platten) -> Stock endet nach der letzten Trennwand, Nase 20 hoch bis zur Wand;
    # hohe Seite -> alle Stoecke gleich lang, davor eine Schraubleiste an der Wand.
    felder_x = [(tr(xm, 0)[0], h) for xm, h in list(vorn.items()) + list(hint.items())]
    xmid = (X0 + X1) / 2
    lo = [h for x, h in felder_x if x < xmid]; hi = [h for x, h in felder_x if x >= xmid]
    tief_links = sum(lo) / max(len(lo), 1) > sum(hi) / max(len(hi), 1)
    if tief_links:
        X0, X1 = XW0 + B_RASTE_O, XW1 - B_LEISTE
        x_tw = min(x for x, h in felder_x) - 7.5 - T_TW          # Ende nach der letzten Trennwand
    else:
        X0, X1 = XW0 + B_LEISTE, XW1 - B_RASTE_O
        x_tw = max(x for x, h in felder_x) + 7.5 + T_TW
    SEITE_TIEF.append(tief_links)
    sole = box(XW0 + B_RASTE_U if tief_links else X0, X1 if tief_links else XW1 - B_RASTE_U, Y0, Y1, zs, zs + T_SOLE)
    sole = schneide(sole, bohrungen([tr(x, y) for x, y in sl], 10.5, zs - 1, zs + T_SOLE + 1))
    teile.append((f"Stimmstock_{reihe}_{ton}_Sohle", sole, "Holz"))
    # Vorderwand an der Sohlenkante, unter der die Loecher der vorderen Kammern liegen
    yv = [min(sl, key=lambda q: abs(q[0] - xm))[1] for xm in vorn] or [sy1]
    ym = sum(yv) / len(yv)
    kv = sy1 if abs(ym - sy1) < abs(ym - sy0) else sy0
    kh = sy0 if kv == sy1 else sy1
    A = NEIGUNG.get(reihe, 0.0)
    yc = sum(tr(x, y)[1] for x, y in sl) / len(sl)            # Stockmitte = Mitte der Lochreihen
    kv_y = tr(0, kv)[1]
    # Stimmstock = ein keilfoermiger Vollkoerper: unten 24, oben 19 breit, Hoehe H_WAND;
    # Vorder- und Rueckseite je leicht gegeneinander geneigt, Kammern 15 x h, 9,5 mm tief (CamBam), Stimmplatten liegen aussen auf.
    phi = math.degrees(math.atan((B_STOCK - B_OBEN) / 2 / H_WAND))
    def gesamt(s_):
        return s_.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -A).translate(cq.Vector(0, yc, zs + T_SOLE))
    _t = math.tan(math.radians(phi)); _u = 20.0                # nach unten verlaengert, spaeter waagrecht auf der Sohle abgeschnitten
    koerper = (cq.Workplane("YZ").polyline([(-B_STOCK/2 - _u*_t, -_u), (B_STOCK/2 + _u*_t, -_u), (B_OBEN/2, H_WAND), (-B_OBEN/2, H_WAND)])
               .close().extrude(X1 - X0).translate((X0, 0, 0)).val())
    cuts_all, platten_all = [], []
    for seite, felder in (("vorn", vorn), ("hinten", hint)):
        vorne_hoch = kv_y > yc
        auf_plus = vorne_hoch if seite == "vorn" else not vorne_hoch
        aussen_y = B_STOCK/2 if auf_plus else -B_STOCK/2
        innen = -1 if auf_plus else 1
        def stelle(s_, aussen_y=aussen_y, innen=innen):
            s_ = s_.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), 90)            # v -> Z, w -> -Y
            if innen == 1: s_ = s_.mirror("XZ")                                 # w -> +Y
            s_ = s_.rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -innen * phi)  # Seitenflaeche geneigt
            return s_.translate(cq.Vector(0, aussen_y, 0))
        for xm, h in felder.items():
            Xm = tr(xm, 0)[0]
            cuts_all.append(stelle(box(Xm - 7.5, Xm + 7.5, -25, 1.0 + h, -1, 9.5)))              # Kammer 15 x h, 9,5 tief (CamBam), unten zur Sohle offen
            pl = box(Xm - 8.0, Xm + 8.0, -25.0, 2.0 + h, -T_PLATTE, 0.0)                         # Stimmplatte liegt aussen auf
            for s2 in (-2.2, 2.2):
                pl = pl.cut(box(Xm + s2 - 1.0, Xm + s2 + 1.0, 4.0, h - 3.0, -5, 5))   # Stimmschlitze durch
            platten_all.append((seite, stelle(pl)))
    koerper = schneide(koerper, cuts_all)
    koerper = gesamt(koerper)
    unten = box(-50, 450, yc - 60, yc + 60, -50, zs + T_SOLE)     # alles unter der Sohlenoberkante weg: Stock, Kammern, Platten buendig
    koerper = koerper.cut(unten)
    # Oberkante: parallel zum Boden, gleichmaessig verlaufend = obere Huelle ueber alle Plattenoberkanten (+2 mm)
    zb = zs + T_SOLE
    kw = math.cos(math.radians(max(0.0, abs(A) - phi)))
    pts = sorted(set((round(x + d, 2), zb + (h + 3.0) * kw) for x, h in felder_x for d in (-7.5, 7.5)))
    huelle = []
    for q in pts:                                                  # obere konvexe Huelle (monotone chain)
        while len(huelle) >= 2 and ((huelle[-1][0]-huelle[-2][0])*(q[1]-huelle[-2][1]) - (huelle[-1][1]-huelle[-2][1])*(q[0]-huelle[-2][0])) >= 0:
            huelle.pop()
        huelle.append(q)
    prof = [(X0 - 1, huelle[0][1])] + huelle + [(X1 + 1, huelle[-1][1])]
    prof = [q for i, q in enumerate(prof) if i == 0 or abs(q[0] - prof[i-1][0]) > 1e-6]
    # am hohen Ende auf Schieberlaenge (40 mm) eben, damit der Riegel plan aufliegt
    def z_at(x):
        for (xa, za), (xb, zb2) in zip(prof, prof[1:]):
            if xa <= x <= xb: return za + (zb2 - za) * (x - xa) / (xb - xa)
        return prof[-1][1] if x > prof[-1][0] else prof[0][1]
    if tief_links:
        xf = X1 - 40.0; zf = max(z_at(xf), z_at(X1))
        prof = [q for q in prof if q[0] < xf] + [(xf, zf), (X1 + 1, zf)]
    else:
        xf = X0 + 40.0; zf = max(z_at(xf), z_at(X0))
        prof = [(X0 - 1, zf), (xf, zf)] + [q for q in prof if q[0] > xf]
    ztop = 200.0
    deckel = (cq.Workplane("XZ").polyline(prof + [(X1 + 1, ztop), (X0 - 1, ztop)]).close()
              .extrude(-80).translate((0, yc - 40, 0)).val())
    koerper = koerper.cut(deckel)
    # Nase an der tiefen Seite: ab der letzten Trennwand bis zur Wand nur 20 mm hoch
    nx0, nx1 = (X0 - 1, x_tw) if tief_links else (x_tw, X1 + 1)
    koerper = koerper.cut(box(nx0, nx1, yc - 40, yc + 40, zb + H_NASE - T_SOLE, ztop))
    # Nasenstirn passend zur schraegen Rastleiste (unten 6, oben 4 breit) abschraegen
    xr_u, xr_o = (XW0 + B_RASTE_U, XW0 + B_RASTE_O) if tief_links else (XW1 - B_RASTE_U, XW1 - B_RASTE_O)
    zr0, zr1 = Z_ST, Z_ST + T_SOLE + H_NASE
    keil = cq.Workplane("XZ").polyline([(xr_u, zr0 - 1), (xr_o - (xr_u - xr_o) * 0.05, zr1 + 1),
                                        ((XW0 - 5) if tief_links else (XW1 + 5), zr1 + 1),
                                        ((XW0 - 5) if tief_links else (XW1 + 5), zr0 - 1)]).close().extrude(-80).translate((0, yc - 40, 0)).val()
    koerper = koerper.cut(keil)
    zz = [q[1] for q in prof]
    Z_ENDE_HOCH.append(prof[-1][1] if tief_links else prof[0][1])
    STOCK_INFO.append((reihe, yc, A, zb, Z_ENDE_HOCH[-1], X1 if tief_links else X0, tief_links))
    print(f"   Oberkante Z {min(zz):.1f} .. {max(zz):.1f}")
    teile.append((f"Stimmstock_{reihe}_{ton}", koerper, "Holz"))
    for i, (seite, pl) in enumerate(platten_all):
        teile.append((f"Stimmplatte_{reihe}_{ton}_{seite}_{i+1}", gesamt(pl).cut(unten), "Alu"))
    print(f"   Neigung {A:+.1f} Grad, Mitte Y {yc:.1f}, Seiten je {phi:.1f} Grad gegeneinander")

# Rastleiste fuer die Nasen an der Wand der tiefen Seite (Trapez: unten 6, oben 4 breit, Hoehe = Sohle + Nase)
_x0 = XW0 if all(SEITE_TIEF) else XW1
_d = 1 if all(SEITE_TIEF) else -1
_h = T_SOLE + H_NASE
raste = (cq.Workplane("XZ").polyline([(_x0, Z_ST), (_x0 + _d * B_RASTE_U, Z_ST), (_x0 + _d * B_RASTE_O, Z_ST + _h), (_x0, Z_ST + _h)])
         .close().extrude(-(B - 2 * RI)).translate((0, RI, 0)).val())
teile.append(("Rastleiste_Nasen", raste, "Holz"))
# Schraubleiste an der Wand der hohen Seite, etwas niedriger als die Stimmstock-Enden
zl = min(Z_ENDE_HOCH) - 5.0
xl0, xl1 = (XW1 - B_LEISTE, XW1) if all(SEITE_TIEF) else (XW0, XW0 + B_LEISTE)
teile.append(("Schraubleiste_Stimmstoecke", box(xl0, xl1, RI, B - RI, Z_ST, zl), "Holz"))
print(f"Schraubleiste X {xl0:.1f}-{xl1:.1f}, Oberkante Z {zl:.1f} (Stockenden min {min(Z_ENDE_HOCH):.1f}); tiefe Seite links: {SEITE_TIEF}")

# Verriegelungsschieber nach Rigel.dxf (0-Eduard/Diskant/Stimmstoecke/Halter, 2 mm Blech, 45 x 9):
# offener Schlitz an einem Ende (rastet unter der langen Schraube in der Schraubleiste), zwei Langloecher 10 x 3
# fuer die beiden Schrauben im Stock; hinteres Ende an der Einschnuerung (u = -7,5) 90 Grad nach oben gebogen.
T_SB, U_BIEG, U_SCHLITZ = 2.0, -7.5, -38.5
rg = ezdxf.readfile("rigel_j.dxf").modelspace()
def rg_kontur(e):
    pts = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.05)]
    out = []
    for q in pts:
        if not out or math.hypot(q[0] - out[-1][0], q[1] - out[-1][1]) > 1e-4: out.append(q)
    if math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) < 1e-4: out.pop()
    return out
rg_pl = sorted(rg.query("LWPOLYLINE"), key=lambda e: -bbox.extents([e]).size.x)
rg_aussen, rg_loecher = rg_kontur(rg_pl[0]), [rg_kontur(e) for e in rg_pl[1:]]
blech = cq.Workplane("XY").polyline(rg_aussen).close().extrude(T_SB).val()
blech = schneide(blech, [cq.Workplane("XY").polyline(l).close().extrude(T_SB + 2).translate((0, 0, -1)).val() for l in rg_loecher])
flach = blech.cut(box(U_BIEG, 5, -10, 10, -1, 5))
schenkel = blech.cut(box(-50, U_BIEG, -10, 10, -1, 5)).rotate(cq.Vector(U_BIEG, 0, 0), cq.Vector(U_BIEG, 1, 0), -90)
def zyl_z(x, y, z0, z1, d):
    return cq.Workplane("XY").center(x, y).circle(d / 2).extrude(z1 - z0).translate((0, 0, z0)).val()
def schraube(x, y, z_auf, laenge, d_k, d_s=3.0, h_k=1.8):
    """Linsenkopf mit Kreuzschlitz, steht senkrecht auf z_auf (lokal), Schaft nach unten."""
    k = cq.Workplane("XY").center(x, y).circle(d_k / 2).extrude(h_k).translate((0, 0, z_auf)).edges(">Z").fillet(h_k * 0.6).val()
    for dx, dy in ((d_k, 0.7), (0.7, d_k)):
        k = k.cut(box(x - dx / 2, x + dx / 2, y - dy / 2, y + dy / 2, z_auf + h_k - 1.0, z_auf + h_k + 1))
    return k.fuse(zyl_z(x, y, z_auf - laenge, z_auf, d_s))
for reihe, yc_, A_, zb_, zt, xe, tl in STOCK_INFO:
    d = 1 if tl else -1                                   # Richtung zur Wand der hohen Seite
    yt = yc_ + (zt - zb_) * math.tan(math.radians(A_))    # Mitte der Stockoberkante
    xs_ = XW1 - 7.5 if tl else XW0 + 7.5                  # lange Rastschraube in der Schraubleiste
    welt = lambda u: xs_ - d * (u - U_SCHLITZ)            # Rigel-Koordinate u -> Modell-X
    def platz(s_):
        s_ = s_.mirror("YZ") if d == 1 else s_
        return s_.translate(cq.Vector(welt(0), yt, zt))
    # Neigung wie der Stimmstock: Schenkel und Schrauben um die X-Achse durch die Auflage geneigt
    neig = lambda s_: s_.rotate(cq.Vector(0, yt, zt), cq.Vector(1, yt, zt), -A_)
    teile.append((f"Schieber_{reihe}", platz(flach).fuse(neig(platz(schenkel))), "Stahl"))
    for i, u in enumerate((-30.0 + 3.5, -15.0 + 3.5)):    # verriegelt: Schrauben am hinteren Ende der Langloecher
        teile.append((f"Schraube_Schieber_{reihe}_{i+1}", neig(schraube(welt(u), yt, zt + T_SB, 12.0, 5.5)), "Schraube"))
    teile.append((f"Rastschraube_{reihe}", neig(schraube(xs_, yt, zt + T_SB + 0.2, zt + T_SB + 0.2 - zl + 12, 6.0)), "Schraube"))
print("Schieber:", len(STOCK_INFO))

# ---------------- Verdeck (aussen, Z < 0) ----------------
vm = ezdxf.readfile("verdeck.dxf").modelspace()
VX, VY = 369.55, 28.0          # Deckel x -359.1..0 -> X 10.45..369.55 ; y 0..168 -> Y 28..196 (hinten ca. 1 cm zurueck, Lage bestimmt spaeter die Mechanik)
def kontur(e):
    if e.dxftype() == "CIRCLE":
        c = e.dxf.center; r = e.dxf.radius
        return [(c.x + r*math.cos(2*math.pi*k/48), c.y + r*math.sin(2*math.pi*k/48)) for k in range(48)]
    return [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.05)][:-1]
ornamente = [kontur(e) for e in vm.query('LWPOLYLINE CIRCLE') if e.dxf.layer == "LAYER2"]
ornamente = [o for o in ornamente if len(o) >= 3]
def prisma(pts2d, z0, z1):
    return cq.Workplane("XY").polyline(pts2d).close().extrude(z1 - z0).translate((0, 0, z0)).val()
# Deckplatte
dk = box(10.45, 369.55, VY, VY + 168.0, -H_VD, -H_VD + T_VD)
werk = [prisma([(x + VX, y + VY) for x, y in o], -H_VD - 1, -H_VD + T_VD + 1)
        for o in ornamente if -359.1 < min(a for a, _ in o) and max(a for a, _ in o) < 0 and min(b for _, b in o) > 0]
teile.append(("Verdeck_Deckplatte", schneide(dk, werk), "Verdeck"))
nd = len(werk)
# Seitenwaende (Abwicklung: Abstand d von der Deckelkante -> Z = -H_VD + d)
def wand_x(xk, richt, X):     # xk Kante in der Zeichnung, richt +1 rechts/-1 links
    sel = [o for o in ornamente if (min(a for a, _ in o) > xk if richt > 0 else max(a for a, _ in o) < xk)
           and abs(min(a for a, _ in o) - xk) < 40]
    s = box(X, X + T_VD, VY + T_VD, VY + 168.0, -H_VD + T_VD, 0)   # bis zur Deckelkante (Stirnseite ausgefuellt)
    ws = []
    for o in sel:
        pts = [(VY + y, -H_VD + abs(x - xk)) for x, y in o]      # (Y, Z)
        ws.append(cq.Workplane("YZ").polyline(pts).close().extrude(T_VD + 2).translate((X - 1, 0, 0)).val())
    return schneide(s, ws), len(ws)
sl_, nl = wand_x(-360.2, -1, 10.45)
sr_, nr = wand_x(1.1, 1, 369.55 - T_VD)
teile += [("Verdeck_Seite_links", sl_, "Verdeck"), ("Verdeck_Seite_rechts", sr_, "Verdeck")]
fs = box(10.45, 369.55, VY, VY + T_VD, -H_VD + T_VD, 0)
sel = [o for o in ornamente if max(b for _, b in o) < -1.1]
ws = [cq.Workplane("XZ").polyline([(x + VX, -H_VD + (-1.1 - y)) for x, y in o]).close().extrude(-(T_VD + 2)).translate((0, VY - 1, 0)).val()
      for o in sel]
teile.append(("Verdeck_Front", schneide(fs, ws), "Verdeck"))
# Verdeck 180 Grad drehen (vorne/hinten vertauscht), gleicher Platz Y 28..196
_m = (10.45 + 369.55) / 2, VY + 84.0
teile = [(n_, s_.rotate(cq.Vector(_m[0], _m[1], 0), cq.Vector(_m[0], _m[1], 1), 180) if f_ == "Verdeck" else s_, f_) for n_, s_, f_ in teile]
# ganzer Diskantteil seitenverkehrt (links <-> rechts): Verdeck in X gespiegelt
teile = [(n_, s_.mirror("YZ", cq.Vector(L / 2, 0, 0)).translate(cq.Vector(0, VD_DY, 0)) if f_ == "Verdeck" else s_, f_) for n_, s_, f_ in teile]
print("Verdeck: Durchbrueche Deckel", nd, "Seiten", nl, nr, "Front", len(ws))

# ---------------- Griffstock + Diskantmechanik (hebel.py, nach griffstock-alternativ1-schnitt.dxf) ----------------
if os.environ.get("MECHANIK", "1") == "1":
    exec(open("hebel.py").read())
    # Verdeck an Griffstock, Abdeckung, Kunststoffteil und Keil anpassen (keine Ueberschneidung),
    # umlaufende Aussenkante der Deckplatte an den freien Seiten R 0,5
    _gsf = [s_ for n_, s_, f_ in teile if n_ in ("Griffstock_Koerper", "Griffstock_Teil2", "Griffstock_Abdeckung", "Kunststoffteil_Lagerblock", "Dreieckleiste")]
    _gsf.append(prisma_uy(_pts(_gr(HUELLE_G.buffer(1.0))), 190.0, 400.0))   # ganzer Griffstockbereich inkl. Holzhebel + 1 mm Luft
    _neu = []
    for n_, s_, f_ in teile:
        if f_ == "Verdeck":
            v0 = s_.Volume()
            for g_ in _gsf:
                if s_.intersect(g_).Volume() > 1e-3: s_ = s_.cut(g_)
            if n_ == "Verdeck_Deckplatte":
                _bb = s_.BoundingBox()
                _ed = [e for e in s_.Edges() if abs(e.Center().z - _bb.zmin) < 0.05 and
                       (abs(e.Center().x - _bb.xmin) < 0.3 or abs(e.Center().x - _bb.xmax) < 0.3 or abs(e.Center().y - _bb.ymax) < 0.3)]
                for _r in (R_VD, 1.0, 0.5):
                    try:
                        s_ = s_.fillet(_r, _ed); print(f"Verdeck: Kante R{_r} an", len(_ed), "Kanten"); break
                    except Exception as ex_: print(f"Verdeck: Abrundung R{_r} nicht moeglich", ex_)
                _sec = s_.intersect(cq.Solid.makeBox(0.2, 300, 1.0, cq.Vector((KST_X0 + KST_X1) / 2 - 0.1, -50, -H_VD)))
                Y_VD_KANTE = _sec.BoundingBox().ymin
            if abs(s_.Volume() - v0) > 1e-3: print(f"Verdeck {n_}: angepasst, {v0 - s_.Volume():.0f} mm3 weggenommen")
        _neu.append((n_, s_, f_))
    teile = _neu
    # Leiste in der Ecke Verdeck/Griffstock (X wie die Ausnehmung KST_X0..KST_X1): flache Seite auf dem Verdeck,
    # flache Seite buendig am Griffstock, Vorderseite und beide Enden gerundet
    W_L, H_L, R_L = 6.0, 5.0, 3.0
    from shapely.geometry import Polygon as _LP, Point as _LPt
    _u0, _yk = H_VD, Y_VD_KANTE
    _pr = _LP([(_u0, _yk - 15), (_u0 + H_L, _yk - 15), (_u0 + H_L, _yk + W_L - R_L), (_u0, _yk + W_L - R_L)]).union(
          _LP([(_u0, _yk - 15), (_u0 + H_L - R_L, _yk - 15), (_u0 + H_L - R_L, _yk + W_L), (_u0, _yk + W_L)])).union(
          _LPt(_u0 + H_L - R_L, _yk + W_L - R_L).buffer(R_L, 64)).difference(HUELLE_G)
    _pr = [g for g in (_pr.geoms if hasattr(_pr, "geoms") else [_pr]) if g.contains(_LPt(_u0 + 0.5, _yk + 1.0))][0]
    def _y_hull(u_):                 # Griffstockflaeche (kleinstes Y der Leiste) in Hoehe u_
        s_ = _pr.intersection(_LP([(u_ - 0.01, -100), (u_ + 0.01, -100), (u_ + 0.01, 200), (u_ - 0.01, 200)]))
        return s_.bounds[1]
    _ya, _yb = _y_hull(_u0 + 0.02), _y_hull(_u0 + H_L - 0.02)
    _wp = (cq.Workplane("YZ").moveTo(_ya, -_u0).lineTo(_yk + W_L, -_u0).lineTo(_yk + W_L, -(_u0 + H_L - R_L))
           .threePointArc((_yk + W_L - R_L + R_L * 0.7071, -(_u0 + H_L - R_L + R_L * 0.7071)), (_yk + W_L - R_L, -(_u0 + H_L)))
           .lineTo(_yb, -(_u0 + H_L)).close())
    hl = _wp.extrude(KST_X1 - KST_X0).translate((KST_X0, 0, 0)).val()
    def _ist_end(e): return all(min(abs(p_.x - xe) for xe in (KST_X0, KST_X1)) < 1e-3 for p_ in (e.startPoint(), e.endPoint()))
    def _auf_deck(e): return all(abs(p_.z + _u0) < 1e-3 for p_ in (e.startPoint(), e.endPoint(), e.positionAt(0.5)))
    def _am_gs(e):
        m = e.positionAt(0.5); u_ = -m.z
        return abs(m.y - (_ya + (_yb - _ya) * (u_ - _u0) / H_L)) < 1e-3
    _ek = [e for e in hl.Edges() if _ist_end(e) and not _auf_deck(e) and not _am_gs(e)]
    for _r in (2.0, 1.5, 1.0):
        try:
            hl = hl.fillet(_r, _ek); print(f"Verdeckleiste: Enden R{_r} ({len(_ek)} Kanten)"); break
        except Exception as ex_: print(f"Verdeckleiste: Enden R{_r} nicht moeglich ({len(_ek)} Kanten)", ex_)
    for n_g, g_ in [(n_, s_) for n_, s_, f_ in teile if n_ in ("Griffstock_Koerper", "Griffstock_Teil2", "Griffstock_Abdeckung", "Kunststoffteil_Lagerblock", "Dreieckleiste")]:
        v_ = hl.intersect(g_).Volume()
        if v_ > 1e-3: hl = hl.cut(g_); print(f"Verdeckleiste: an {n_g} {v_:.2f} mm3 beschnitten")
    teile.append(("Verdeck_Leiste", hl, "Verdeck"))
    print(f"Verdeckleiste {W_L}x{H_L} R{R_L}: X {KST_X0:.0f}-{KST_X1:.0f}, Kante Y {Y_VD_KANTE:.1f}, Volumen {hl.Volume():.0f}")

# ---- Verdeckstuetzen (Johann 16.11.): 3-mm-Platte nach Verdeck-Stuetzen.dxf, 25 hoch, innen an den Seitenteilen des Verdecks,
# in der Mitte der Seitenteile. Je Fuss eine M2,4x8 Sechskantschraube von aussen durch den Fuss in den Boden;
# das Verdeck wird von der Seite verschraubt: Raendelschraube M2,4 waagrecht durch die Seitenwand in die Stuetze.
_st_dxf = [(0, 10), (0, 0), (-3, 0), (-3, 5), (-25, 10), (-25, 30), (-3, 35), (-3, 40), (0, 40), (0, 30), (-8, 26), (-8, 14)]
def _raendel_x(x_aussen, y, z, richtung, dk, hk, d, l):     # Kopf aussen an der Wand, Schaft in Richtung +-X
    x0 = x_aussen - richtung * hk
    kopf = cq.Solid.makeCylinder(dk / 2, hk, cq.Vector(x0, y, z), cq.Vector(richtung, 0, 0))
    for i in range(24):
        a = 2 * math.pi * i / 24
        kopf = kopf.cut(cq.Solid.makeCylinder(0.25, hk + 0.2, cq.Vector(x0 - richtung * 0.1, y + dk / 2 * math.cos(a), z + dk / 2 * math.sin(a)), cq.Vector(richtung, 0, 0)))
    return kopf.fuse(cq.Solid.makeCylinder(d / 2, l, cq.Vector(x_aussen, y, z), cq.Vector(richtung, 0, 0)))
_vd = {n_: s_ for n_, s_, f_ in teile if n_.startswith("Verdeck_Seite")}
_wb = [(n_, s_.BoundingBox()) for n_, s_ in _vd.items()]
_boden_l = []; _wand_l = {}; STUETZEN = []
_Z_SCHR = -(H_VD - T_VD) / 2                                   # Mitte der Seitenwand (Hoehe)
for n_w, b_ in _wb:
    links = b_.xmin < L / 2
    s_n = "links" if links else "rechts"
    xi = b_.xmax if links else b_.xmin                          # Innenseite der Wand
    xa = b_.xmin if links else b_.xmax                          # Aussenseite
    ri = 1.0 if links else -1.0                                 # nach innen
    xs_ = xi + ri * 2.6; ym_ = (b_.ymin + b_.ymax) / 2          # 1,1 mm Abstand zur Wand (Sechskantkoepfe der Fussschrauben)
    STUETZEN.append((round(xs_, 2), round(ym_, 2)))
    pts_ = [(ym_ - 20 + yy, float(xx)) for xx, yy in _st_dxf]   # (Y, Z): Fuesse am Boden (Z=0), oben Z=-25
    stz = cq.Workplane("YZ").polyline(pts_).close().extrude(3.0).translate((xs_ - 1.5, 0, 0)).val()
    for yf in (2.5, 37.5):
        _z_auf = -3.0
        _kopf = cq.Workplane("XY").polygon(6, 5.0 / math.cos(math.pi / 6)).extrude(1.7).val().rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 30).translate(cq.Vector(xs_, ym_ - 20 + yf, _z_auf - 1.7))
        teile.append((f"Schraube_Stuetze_{s_n}_{int(yf)}", _kopf.fuse(cq.Solid.makeCylinder(1.2, 8.0, cq.Vector(xs_, ym_ - 20 + yf, _z_auf), cq.Vector(0, 0, 1))), "Schraube"))
        _boden_l.append(cq.Solid.makeCylinder(1.2, 9.0, cq.Vector(xs_, ym_ - 20 + yf, _z_auf - 0.5), cq.Vector(0, 0, 1)))
        stz = stz.cut(cq.Solid.makeCylinder(1.3, 4.0, cq.Vector(xs_, ym_ - 20 + yf, _z_auf - 0.5), cq.Vector(0, 0, 1)))
    _l_schr = abs(xi - xa) + 1.1 + 3.0                          # durch Wand, Luft und Stuetze
    stz = stz.cut(cq.Solid.makeCylinder(1.2, 6.0, cq.Vector(xi, ym_, _Z_SCHR), cq.Vector(ri, 0, 0)))     # Gewinde M2,4
    _wand_l[n_w] = cq.Solid.makeCylinder(1.3, abs(xi - xa) + 1.0, cq.Vector(xa - ri * 0.5, ym_, _Z_SCHR), cq.Vector(ri, 0, 0))
    teile.append((f"Verdeckstuetze_{s_n}", stz, "Holz"))
    teile.append((f"Verdeck_Raendelschraube_{s_n}", _raendel_x(xa, ym_, _Z_SCHR, ri, 8.0, 3.0, 2.4, _l_schr), "Schraube"))
teile = [(n_, schneide(s_, _boden_l) if n_ == "Diskantboden" else (s_.cut(_wand_l[n_]) if n_ in _wand_l else s_), f_) for n_, s_, f_ in teile]
print("Verdeckstuetzen (X, Y Mitte):", STUETZEN, "Schraube seitlich bei Z", _Z_SCHR)

# ---- Verdeck-Front: in der Mitte fehlte eine Blume -> Rosette aus Blume1.dxf in die Luecke (durchgebrochen wie die anderen)
_vf = [s_ for n_, s_, f_ in teile if n_ == "Verdeck_Front"][0]
_bv = _vf.BoundingBox(); _ymv = (_bv.ymin + _bv.ymax) / 2
_ff = max([f for f in _vf.Faces() if abs(abs(f.normalAt().y) - 1) < 1e-6], key=lambda f: f.Area())
_ros = []
for w_ in _ff.innerWires():
    fa_ = cq.Face.makeFromWires(w_); a_ = fa_.Area()
    if 13.0 < a_ < 16.0: c_ = fa_.Center(); _ros.append((round(c_.x, 2), round(c_.z, 2)))
_ros = sorted(set(_ros))
print("Rosetten-Mitten Verdeck-Front:", _ros)
_gp = max(zip(_ros, _ros[1:]), key=lambda q: q[1][0] - q[0][0])
_xc, _zc = (_gp[0][0] + _gp[1][0]) / 2, _gp[0][1]
_b1 = ezdxf.readfile("blume1.dxf").modelspace(); _rx, _ry = -765.9, 3.3
_rw = []
for e in _b1:
    if e.dxf.layer != "LAYER2": continue
    if e.dxftype() == "CIRCLE" and math.hypot(e.dxf.center.x - _rx, e.dxf.center.y - _ry) < 0.2:
        _rw.append(cq.Solid.makeCylinder(e.dxf.radius, 20, cq.Vector(_xc, _bv.ymin - 5, _zc), cq.Vector(0, 1, 0)))
    elif e.dxftype() == "LWPOLYLINE" and e.closed:
        p_ = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.03)][:-1]
        if max(math.hypot(x_ - _rx, y_ - _ry) for x_, y_ in p_) < 9.5:
            _rw.append(cq.Workplane("XZ").polyline([(_xc + (x_ - _rx), _zc + (y_ - _ry)) for x_, y_ in p_]).close().extrude(-20).translate((0, _bv.ymin - 5, 0)).val())
teile = [(n_, schneide(s_, _rw) if n_ == "Verdeck_Front" else s_, f_) for n_, s_, f_ in teile]
print(f"Verdeck-Front: Rosette in der Mitte X {_xc:.1f} Z {_zc:.1f} ({len(_rw)} Teile)")

# ---- Rahmen vorne: Einlage mittig = Venus-Blume (Front-einlage, Fron-venus2.dxf): Ring (R16,5-17,5) und Blumenmuster
# schwarz eingelegt, 1,5 mm tief, buendig, auf der Front zentriert
T_EINL, R_RING_I = 1.5, 16.5
_vm = ezdxf.readfile("venus.dxf").modelspace()
_kr = list(_vm.query('CIRCLE[layer=="Layer9"]'))[0]
_cx, _cy, _r_out = _kr.dxf.center.x, _kr.dxf.center.y, _kr.dxf.radius
_rv = [s_ for n_, s_, f_ in teile if n_ == "Rahmen_vorne"][0]; _br = _rv.BoundingBox()
_yv, _xm, _zm = _br.ymax, (_br.xmin + _br.xmax) / 2, H_RAHMEN / 2
_ring = cq.Solid.makeCylinder(_r_out, T_EINL, cq.Vector(_xm, _yv - T_EINL, _zm), cq.Vector(0, 1, 0)).cut(
        cq.Solid.makeCylinder(R_RING_I, T_EINL + 1, cq.Vector(_xm, _yv - T_EINL - 0.5, _zm), cq.Vector(0, 1, 0)))
_muster = []
for e in _vm.query('LWPOLYLINE[layer=="Layer10"]'):
    if not e.closed: continue
    p_ = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.03)][:-1]
    _muster.append(cq.Workplane("XZ").polyline([(_xm - (x_ - _cx), _zm + (y_ - _cy)) for x_, y_ in p_]).close().extrude(-T_EINL).translate((0, _yv - T_EINL, 0)).val())
_rv = schneide(_rv, [_ring] + _muster)
teile = [(n_, _rv if n_ == "Rahmen_vorne" else s_, f_) for n_, s_, f_ in teile]
teile.append(("Rahmeneinlage/Ring", _ring, "Einlage"))
for i_, s_ in enumerate(_muster): teile.append((f"Rahmeneinlage/Blume_{i_ + 1:02d}", s_, "Einlage"))
print(f"Rahmen vorne: Venus-Blume Ring R{R_RING_I}-{_r_out} + {len(_muster)} Musterteile schwarz, Mitte X {_xm:.1f} Z {_zm:.1f}")

# ---------------- Ausgabe ----------------
farben = {"Einlage": cq.Color(0.12, 0.07, 0.04), "Holz": cq.Color(0.80, 0.62, 0.40), "Alu": cq.Color(0.78, 0.80, 0.83), "Verdeck": cq.Color(0.15, 0.15, 0.15),
          "Leder": cq.Color(0.55, 0.30, 0.15), "Stahl": cq.Color(0.55, 0.57, 0.60), "Schraube": cq.Color(0.80, 0.65, 0.20), "Klappe": cq.Color(0.30, 0.20, 0.12), "Hebelholz": cq.Color(0.70, 0.45, 0.25), "Kunststoff": cq.Color(0.92, 0.92, 0.88), "Knopf": cq.Color(0.95, 0.93, 0.85), "Griffstock": cq.Color(0.45, 0.28, 0.15),
          # Hebel je Knopfreihe eigene Farbe (Alu kraeftig, Holzteil heller)
          "Alu_R1": cq.Color(0.85, 0.15, 0.15), "Holz_R1": cq.Color(0.95, 0.55, 0.55),
          "Alu_R2": cq.Color(0.15, 0.35, 0.90), "Holz_R2": cq.Color(0.55, 0.70, 0.98),
          "Alu_R3": cq.Color(0.10, 0.65, 0.20), "Holz_R3": cq.Color(0.55, 0.88, 0.55),
          "Alu_R4": cq.Color(0.95, 0.60, 0.05), "Holz_R4": cq.Color(0.99, 0.85, 0.45)}
import os
asm = cq.Assembly(name=f"Diskantteil_Markus_{VAR}")
import re as _re
_filt = os.environ.get("NUR")          # z.B. NUR="Diskantboden|Klappe|Griffstock|Rahmen_vorne" fuer leichte Ansicht
if _filt: teile = [t_ for t_ in teile if _re.search(_filt, t_[0])]
_grp = {}
for name, s, f in teile:
    if name.startswith("Verdeck_"): name = "Verdeck/" + name[len("Verdeck_"):]     # Verdeck als eine Einheit (Gruppe)
    if "/" in name:
        g_, n2 = name.split("/", 1)
        if g_ not in _grp:
            _grp[g_] = cq.Assembly(name=g_)
        _grp[g_].add(cq.Workplane().add(s), name=n2, color=farben[f])
    else:
        asm.add(cq.Workplane().add(s), name=name, color=farben[f])
for g_, ga_ in _grp.items():            # Untergruppen erst fuellen, dann einhaengen (add kopiert die Gruppe)
    asm.add(ga_, name=g_)
asm.save(OUT)
print("Teile", len(teile), "->", OUT)
# Netze je Teil fuer Animationen (ANIMDIR=Ordner)
if os.environ.get("ANIMDIR"):
    import numpy as _np, pickle as _pk
    _out = {}
    for n, s, f in teile:
        try:
            v, t = s.tessellate(0.15, 0.4)
            _out[n] = (f, _np.array([(p.x, p.y, p.z) for p in v], dtype=_np.float32), _np.array(t, dtype=_np.int32))
        except Exception as ex: print("Netz fehlt", n, ex)
    _pk.dump({"teile": _out, "farben": {k: tuple(c.toTuple()[:3]) for k, c in farben.items()}, "drehpunkte": DREHPUNKTE if "DREHPUNKTE" in globals() else {}},
             open(os.path.join(os.environ["ANIMDIR"], "netze.pkl"), "wb"))
    print("Netze:", len(_out))
# Vorschau-STL je Farbe
import os
if os.environ.get("STL"):
    for f in farben:
        cs = [s for n, s, ff in teile if ff == f]
        cq.exporters.export(cq.Compound.makeCompound(cs), f"vor_{f}.stl", tolerance=0.1, angularTolerance=0.3)
