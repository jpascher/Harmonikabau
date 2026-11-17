# Diskantmechanik Markus (wird von diskant.py per exec eingebunden; nutzt box, teile, bm, BXO, cq, ezdxf, math)
# Hebel: Alu 2 mm nach 0-Markus/Diskant/hebel/<n>.reihe (Holzrand-Zeichnungen), Holzteil 5 mm am kurzen Arm,
# kurze feste Achse Ø5 x 11 im Lagerpunkt, Lagerpunkte = tiefe Nuten im schraeg eingesetzten Kunststoffteil.
# Langer Arm in der Hebelebene bis vor die Klappe, dann seitlich gekroepft zur Klappenmitte und mit einer Stufe
# auf die Klappe (damit sich offene Klappen und Nachbarhebel nicht beruehren).
# Koordinaten: u = Abstand nach aussen von der Bodenaussenseite (Z = -u), Y wie Modell, X wie Modell.

# --- Klappenmitten aus dem Boden (Mittel der drei Loecher je Block) ---
klappen_m = []
_ges = []
for ins in bm.query('INSERT[layer=="2Füllung-Öffnug"]'):
    p0 = (round(ins.dxf.insert.x, 2), round(ins.dxf.insert.y, 2))
    if p0 in _ges: continue
    _ges.append(p0)
    cs = [(v.dxf.center.x + BXO, v.dxf.center.y) for v in ins.virtual_entities() if v.dxftype() == "CIRCLE"]
    klappen_m.append((sum(c[0] for c in cs) / len(cs), sum(c[1] for c in cs) / len(cs)))
# Klappenreihen nach Y (hinten = weit weg vom Griffstock)
def kl_reihe(y):
    return 1 if y > 152 else 2 if y > 117 else 3 if y > 82 else 4
# --- Hebelpositionen aus der Griffstock-Aufsicht (Layer Lagerbock_aufsicht, Schlitze 5 x 14,2) ---
gm = ezdxf.readfile("gs.dxf").modelspace()
from collections import defaultdict as _dd
_lr = _dd(list)
for e in gm:
    if e.dxf.layer != "Lagerbock_aufsicht": continue
    b = bbox.extents([e])
    if b.has_data and abs(b.size.x - 5.0) < 0.2 and abs(b.size.y - 14.2) < 0.3 or (b.has_data and abs(b.size.x - 5.0) < 0.2 and abs(b.size.y - 11.3) < 0.3):
        _lr[round(b.center.y, 0)].append(b.center.x + BXO)
LAGER_Y = sorted(_lr.keys(), reverse=True)            # 77, 69, 58, 50  -> Hebelreihe 1..4 (Annahme)
HEBEL_X = {i + 1: sorted(_lr[y]) for i, y in enumerate(LAGER_Y)}
# --- Lagerpunkte je Reihe: Kreis der Seitenansicht (u 29,2 / Y 20,2) fuer Reihe 1, die anderen entlang der Schraege
NEIG_GS = math.radians(15.0)
def lager(r):
    d = LAGER_Y[0] - LAGER_Y[r - 1]
    return 29.2 + d * math.cos(NEIG_GS), 20.2 - d * math.sin(NEIG_GS)
HEBEL_DXF = {1: ("h1.dxf", "Hebel-alu-1r"), 2: ("h2.dxf", "Hebel-alu-2r"), 3: ("h3.dxf", "Hebel-alu-3r"), 4: ("h4.dxf", "Hebel-4r")}
def kontur_pts(e):
    pts = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.05)]
    out = []
    for q in pts:
        if not out or math.hypot(q[0] - out[-1][0], q[1] - out[-1][1]) > 1e-3: out.append(q)
    if len(out) > 2 and math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) < 1e-3: out.pop()
    return out
def hebel_profil(r):
    f, lay = HEBEL_DXF[r]
    m_ = ezdxf.readfile(f).modelspace()
    alu = max((e for e in m_ if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == lay and len(e) > 6),
              key=lambda e: bbox.extents([e]).size.y)
    holz = [kontur_pts(e) for e in m_ if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "0" and not e.closed]
    ring = None
    if len(holz) >= 2:                                  # zwei offene Linien -> geschlossener Holzumriss (groessere Flaeche)
        def fl(q): return abs(sum(q[i][0]*q[(i+1) % len(q)][1] - q[(i+1) % len(q)][0]*q[i][1] for i in range(len(q)))) / 2
        cand = [holz[0] + holz[1], holz[0] + list(reversed(holz[1]))]
        ring = max(cand, key=fl)
    return kontur_pts(alu), ring
# Oberseite des Griffstocks (Abstufungen vorerst als Gerade): Y = 45,67 - 0,2693 u
def oben_Y(u): return 45.67 - (45.67 - 10.39) / 131.01 * u
def drehwinkel(alu2d, up, Yp, abstand=7.0):
    """Hebel so um den Lagerpunkt drehen, dass die Spitze des kurzen Arms 'abstand' unter der Griffstock-Oberseite liegt."""
    sp = max(alu2d, key=lambda q: q[0])
    best = None
    for k in range(-400, 401):
        th = math.radians(k * 0.1)
        xr, yr = dreh([sp], th)[0]
        err = abs((Yp + yr) - (oben_Y(up + xr) - abstand))
        if best is None or err < best[0]: best = (err, th)
    return best[1]
def dreh(pts2d, th):
    c_, s_ = math.cos(th), math.sin(th)
    return [(x_ * c_ - y_ * s_, x_ * s_ + y_ * c_) for x_, y_ in pts2d]
def in_ebene(pts2d, x0, dicke, up, Yp, th=0.0):
    """Profil (xl, yl) um th gedreht in die Hebelebene: u = up + xl (nach aussen, Z = -u), Y = Yp + yl; Dicke in X."""
    return (cq.Workplane("YZ").polyline([(Yp + yl, -(up + xl)) for xl, yl in dreh(pts2d, th)]).close()
            .extrude(dicke).translate((x0 - dicke / 2, 0, 0)).val())

T_ALU, T_HOLZ, Z_KLAPPE_OBEN = 2.0, 5.0, -(T_POLSTER + T_KLAPPE)
zuordnung = []
# Feste Zuordnung Knopfreihe r -> Stimmstockreihe r (1 = aussen = B ... 4 = innen = Des);
# Ausnahme: eine Des-Klappe (4. Reihe) haengt an der aeusseren Knopfreihe (dort 14 Hebelplaetze).
kl_r = {r: sorted([k for k in klappen_m if kl_reihe(k[1]) == r], key=lambda k: k[0]) for r in (1, 2, 3, 4)}
frei_r = {r: list(HEBEL_X[r]) for r in (1, 2, 3, 4)}
for r in (1, 2, 3, 4):
    for k in kl_r[r][:] if r != 4 else []:
        x = min(frei_r[r], key=lambda h: abs(h - k[0])); frei_r[r].remove(x); zuordnung.append((r, x, k))
# Ausnahme-Klappe: Des-Klappe am naechsten zum freien Platz der aeusseren Reihe
if frei_r[1]:
    xf = frei_r[1][0]
    kx = min(kl_r[4], key=lambda k: abs(k[0] - xf)); kl_r[4].remove(kx); frei_r[1].remove(xf)
    zuordnung.append((1, xf, kx))
    print(f"Ausnahme: Des-Klappe X {kx[0]:.1f} an aeusserer Knopfreihe, Hebel X {xf:.1f}")
for k in kl_r[4]:
    x = min(frei_r[4], key=lambda h: abs(h - k[0])); frei_r[4].remove(x); zuordnung.append((4, x, k))
for r in (1, 2, 3, 4):
    up, Yp = lager(r)
    alu2d, holz2d = hebel_profil(r)
    TH = drehwinkel(alu2d, up, Yp)
    print(f"Hebelreihe {r}: Drehung {math.degrees(TH):+.1f} Grad")
    print(f"Hebelreihe {r}: Lager u {up:.1f} Y {Yp:.1f}, Hebel {sum(1 for z in zuordnung if z[0]==r)}, freie Plaetze {[round(f_,1) for f_ in frei_r[r]]}")
    # Hebelprofil: langer Arm wird vor der Klappe abgeschnitten (Y_klappe - 14), dann Kroepfung
    for i, (rr, x, k) in enumerate([z for z in zuordnung if z[0] == r]):
        hb = in_ebene(alu2d, x, T_ALU, up, Yp, TH)
        if holz2d:   # am kurzen Arm steht kein Alu ueber das Holzteil (nachtraeglich mit Kontur gefraest)
            xw0 = min(q[0] for q in holz2d)
            zone = in_ebene([(xw0, -200), (400, -200), (400, 200), (xw0, 200)], x, T_ALU + 2, up, Yp, TH)
            holzform = in_ebene(holz2d, x, T_ALU + 2, up, Yp, TH)
            hb = hb.cut(zone).fuse(hb.intersect(zone).intersect(holzform))
        y_end = k[1] - 14.0
        hb = hb.cut(box(x - 5, x + 5, y_end, 400, -300, 50))
        # Ende des Arms in der Ebene: u-Lage dort, wo der Arm Y = y_end schneidet (lineare Naeherung ueber die Armspitze)
        tip = max(dreh(alu2d, TH), key=lambda q: q[1])
        t_ = (y_end - Yp) / max(tip[1], 1e-6)
        u_end = up + tip[0] * t_
        z_arm = -u_end
        b_ = 3.5                                           # Armbreite an der Kroepfung
        xa, xb = sorted((x, k[0]))
        stuecke = [box(xa - T_ALU / 2, xb + T_ALU / 2, y_end - b_, y_end, z_arm - T_ALU / 2, z_arm + T_ALU / 2),       # seitliche Kroepfung
                   box(k[0] - T_ALU / 2, k[0] + T_ALU / 2, y_end - b_, k[1], z_arm - T_ALU / 2, z_arm + T_ALU / 2)]       # weiter zur Klappenmitte
        # Stufe hinunter auf die Klappe
        z0, z1 = sorted((z_arm, Z_KLAPPE_OBEN))
        stuecke.append(box(k[0] - T_ALU / 2, k[0] + T_ALU / 2, k[1] - 3.0, k[1], z0 - T_ALU / 2, z1))
        stuecke.append(box(k[0] - 4.0, k[0] + 4.0, k[1] - 6.0, k[1] + 2.0, Z_KLAPPE_OBEN - T_ALU, Z_KLAPPE_OBEN))  # Fuss auf der Klappe
        for s_ in stuecke: hb = hb.fuse(s_)
        teile.append((f"Hebel_R{r}_{i+1}", hb, "Alu"))
        if holz2d:
            try:
                teile.append((f"Hebelholz_R{r}_{i+1}", in_ebene(holz2d, x, T_HOLZ, up, Yp, TH), "Hebelholz"))
            except Exception as ex_:
                print("Holzteil", r, i, ex_)
        # kurze feste Achse
        teile.append((f"Hebelachse_R{r}_{i+1}", cq.Workplane("YZ").center(Yp, -up).circle(2.5).extrude(11.0).translate((x - 5.5, 0, 0)).val(), "Stahl"))
print("Hebel:", len(zuordnung))

# --- Kunststoffteil (Lagerblock) schraeg im Griffstock: 290 lang, 44,6 breit, 11 dick; Lagerpunkte in den tiefen Nuten
u1, Y1_ = lager(1); u4, Y4_ = lager(4)
s_dir = (math.cos(NEIG_GS), -math.sin(NEIG_GS)); n_dir = (math.sin(NEIG_GS), math.cos(NEIG_GS))   # (u, Y)
a0 = -8.0; a1 = a0 + 44.6
ecken = []
for a, n in ((a0, 0.0), (a1, 0.0), (a1, -11.0), (a0, -11.0)):
    ecken.append((u1 + a * s_dir[0] + (n + 2.5) * n_dir[0], Y1_ + a * s_dir[1] + (n + 2.5) * n_dir[1]))
kst = cq.Workplane("YZ").polyline([(Y, -u) for u, Y in ecken]).close().extrude(290.0).translate((190 - 145, 0, 0)).val()
schl = []
for r, x, k in zuordnung:
    up, Yp = lager(r)
    schl.append(cq.Workplane("YZ").center(Yp, -up).circle(2.6).extrude(12.0).translate((x - 6.0, 0, 0)).val())   # Achsnut
    schl.append(box(x - 1.5, x + 1.5, -40, 80, -150, 0))                                                         # Hebelschlitz
kst = schneide(kst, schl)
teile.append(("Kunststoffteil_Lagerblock", kst, "Kunststoff"))
