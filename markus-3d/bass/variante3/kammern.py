# Bassteil Variante 3 "traditionell, flacher Einbau": Markus-Gehaeuse, flacher Boden, Helikon-Kammern 4/5, Knoepfe 9 + 8
# X = Laenge (382), Y = Breite (207), Z = Hoehe; Boden-Oberkante Z_BODEN, Kammern stehen darauf (Z steigt Richtung Balg)
# Aufruf: python3 kammern.py Bassboden2.dxf ausgabe.step
import sys, math, ezdxf, ezdxf.path, cadquery as cq
from ezdxf import bbox
D, L, B = 5.5, 382.0, 207.0
Z_BODEN, T_BODEN = 31.2, 3.5                # Boden 3,5 mm (Oberkante wie Markus)
H_KAM   = 54.0      # Kammerhoehe = Leistenhoehe (Streifen 54 mm in der Trennwandzeichnung)
T_LEIST = 6.0       # gemeinsame Laengsleiste
T_WAND  = 3.6       # Trennwand
T_NUT   = 2.0       # Nuttiefe in der Leiste
T_FOLIE = 0.6
BX, BY  = 377.0, 5.0            # Bassboden2.dxf -> Modell: X = x + BX, Y = y + BY
v = lambda x, y: (x + BX, y + BY)

def box(x0, x1, y0, y1, z0, z1):
    return cq.Workplane("XY").box(x1-x0, y1-y0, z1-z0, centered=False).translate((x0, y0, z0)).val()

# ---------- Rahmen aus den Einzelteilen (wie Rahmen_Bass_Markus_v2) ----------
imp = lambda f: cq.importers.importStep(f).val()
O, XA, ZA = cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), cq.Vector(0, 0, 1)
st = lambda s: s.rotate(O, XA, 90)
rueck = st(imp("Markus_rueck_v4.step")).translate(cq.Vector(0, D, 0))
front = st(imp("Markus_front_v4.step")).rotate(O, ZA, 180).translate(cq.Vector(L, B - D, 0))
seite = imp("Markus_seite_v4.step")
sl = st(seite).rotate(O, ZA, -90).translate(cq.Vector(D, B, 0))
sr = st(seite).rotate(O, ZA, 90).translate(cq.Vector(L - D, 0, 0))

# ---------- Boden aus Bassboden2.dxf (inkl. Klappenoeffnungen in den Bloecken) ----------
m = ezdxf.readfile(sys.argv[1]).modelspace()
boden = cq.Workplane("XY").box(L - 2*D, B - 2*D, T_BODEN, centered=False).translate((D, D, Z_BODEN - T_BODEN))
def loch(wp):
    global boden; boden = boden.cut(wp.extrude(T_BODEN + 2).translate((0, 0, Z_BODEN - T_BODEN - 1)))
ents = list(m)
for ins in m.query("INSERT"): ents += list(ins.virtual_entities())
n_rect = n_rund = 0
# Auf der Mechanikseite (Unterseite) um jede Klappenoeffnung eine Vertiefung: 5 mm Rand, 2 mm tief,
# damit die Klappen etwas tiefer liegen. Klappenaufbau: 1 mm Leder, 4 mm Filz, 1,5 mm Epoxy.
RAND_V, TIEF_V, SPIEL = 5.0, 2.0, 0.5
T_LEDER, T_FILZ, T_EPOXY = 1.0, 4.0, 1.5
zu = Z_BODEN - T_BODEN                      # Unterseite Boden (Mechanikseite)
klappen = []
def vertiefung_und_klappe(form2d, name):
    global boden
    v_ = form2d(RAND_V)
    boden = boden.cut(v_.extrude(TIEF_V + 1).translate((0, 0, zu - 1)))
    z0_ = zu + TIEF_V
    lagen = []
    for nm, dk in (("Leder", T_LEDER), ("Filz", T_FILZ), ("Epoxy", T_EPOXY)):
        lagen.append((f"{name}_{nm}", form2d(RAND_V - SPIEL).extrude(dk).translate((0, 0, z0_ - dk)).val())); z0_ -= dk
    klappen.extend(lagen)
# Variante 3 (traditionell, flacher Einbau): Bassklappen-Oeffnungen nach der Helikon-Zeichnung (X wie dort, um 1 mm
# auf die Markus-Mitte 191 verschoben), Hinterkante im selben Abstand zur Rueckseite wie bei Markus (Y 68.7), damit die
# Mechanik zwischen Bass- und Begleiterklappen Platz hat. Begleiter-Rundoeffnungen wie Markus, um DYB nach vorne.
DYB = 8.0                                                         # Platz fuer drei Achsreihen zwischen Bass- und Begleiterklappen
DXB = -23.5                                                       # Begleiter (Zwischenteil inkl. 10-mm-Verlaengerung, Stimmstock, Klappen) mittig im Gehaeuse
HM = ezdxf.readfile("kb_haupt.dxf").modelspace()
def _rechtecke(layer):
    out = set()
    for e in HM.query(f'LWPOLYLINE[layer=="{layer}"]'):
        if not e.closed: continue
        p = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.1)]; xs = [a for a, _ in p]; ys = [b for _, b in p]
        if -5 < min(xs) and max(xs) < 390 and 250 < min(ys) and max(ys) < 470:
            out.add((round(min(xs), 2), round(max(xs), 2), round(min(ys), 2), round(max(ys), 2)))
    return sorted(out)
Y_OEFF0 = 68.7
bo_oeff = [o for o in _rechtecke("2Füllung-Öffnug") if 340 < o[2] < 380 and o[1] - o[0] < 60]
_y0min = min(o[2] for o in bo_oeff)
for o in bo_oeff:
    x0_, x1_ = o[0] + 6.0, o[1] + 6.0
    y0_ = Y_OEFF0 + (o[2] - _y0min); y1_ = y0_ + (o[3] - o[2])
    loch(cq.Workplane("XY").center((x0_ + x1_) / 2, (y0_ + y1_) / 2).rect(x1_ - x0_, y1_ - y0_)); n_rect += 1
    vertiefung_und_klappe(lambda o_, a=x0_, b=x1_, c=y0_, d=y1_: cq.Workplane("XY").center((a+b)/2, (c+d)/2).rect(b-a+2*o_, d-c+2*o_), f"Bassklappe_{n_rect}")
for e in ents:
    if e.dxftype() == "CIRCLE" and e.dxf.layer == "Begleiterstimmstock-fuellung-öffnungen":
        c_ = v(e.dxf.center.x, e.dxf.center.y); c_ = (c_[0] + DXB, c_[1] + DYB); r_ = e.dxf.radius
        loch(cq.Workplane("XY").center(*c_).circle(r_)); n_rund += 1
        vertiefung_und_klappe(lambda o, c=c_, r=r_: cq.Workplane("XY").center(*c).circle(r + o), f"Begleiterklappe_{n_rund}")
# Zusatzkammer vorne rechts (quer): Oeffnung 27 x 19 wie die schmalen Kammern, 1 mm hinter der Endleiste am dicken Ende
XZ0_, XZ1_ = 272.5, L - D
FY0_, FY1_ = B - D - T_LEIST, B - D
_ycz = ((FY0_ + 1.5 - 33.0 - 0.5 + 1.75) + FY0_) / 2
OEFF_Z = (XZ0_ + 7.0, XZ0_ + 34.0, _ycz - 9.5, _ycz + 9.5)
loch(cq.Workplane("XY").center((OEFF_Z[0] + OEFF_Z[1]) / 2, (OEFF_Z[2] + OEFF_Z[3]) / 2).rect(27.0, 19.0)); n_rect += 1
vertiefung_und_klappe(lambda o_, a=OEFF_Z: cq.Workplane("XY").center((a[0]+a[1])/2, (a[2]+a[3])/2).rect(27.0 + 2*o_, 19.0 + 2*o_), "Zusatzklappe")

print("Boden: Klappenoeffnungen", n_rect, "runde Oeffnungen", n_rund, "| Vertiefungen + Klappen", len(klappen)//3)
boden = boden.val()

# ---------- Kammern ----------
Z0, Z1 = Z_BODEN, Z_BODEN + H_KAM
def leiste(x0, x1, y0, y1, nuten_x=None, nuten_y=None, innen=None):
    """Laengsleiste mit 2-mm-Nuten (T_WAND breit) auf der Innenseite"""
    s = box(x0, x1, y0, y1, Z0, Z1)
    for nx in (nuten_x or []):                       # Leiste laeuft in X, Nut quer, Innenseite +Y oder -Y
        ya, yb = (y1 - T_NUT, y1 + 1) if innen == "+" else (y0 - 1, y0 + T_NUT)
        s = s.cut(box(nx - T_WAND/2, nx + T_WAND/2, ya, yb, Z0 - 1, Z1 + 1))
    for ny in (nuten_y or []):
        xa, xb = (x1 - T_NUT, x1 + 1) if innen == "+" else (x0 - 1, x0 + T_NUT)
        s = s.cut(box(xa, xb, ny - T_WAND/2, ny + T_WAND/2, Z0 - 1, Z1 + 1))
    return s

def folie_y(xa, xb, y_vorne, y_hinten, umgedreht=False, achse="y"):
    """Gekruemmte Folie (0,6 mm) nach dem Schnitt: flach 16,5 mm ueber dem Boden, Bogen R 36,3 zur hinteren Leiste.
       Profil in der Y-Z-Ebene, quer ueber die Kammerbreite (X) gezogen."""
    lang = y_hinten - y_vorne                       # lichte Laenge zwischen den Leisten
    s = lambda t: (y_hinten - t) if umgedreht else (y_vorne + t)   # Abstand von der vorderen Leiste -> Y
    h, r = 16.5, 36.3
    t_kn = lang - 23.2                              # Uebergang flach -> Bogen (Bogenmitte 23,2 vor der hinteren Leiste)
    pts = [(23.2, h)]                               # Beginn flacher Teil (aus dem Schnitt)
    pts += [(t_kn + r*math.sin(a), h + 19.7 + r*math.cos(a) - r - 0) for a in []]
    prof = [(23.2, h), (t_kn, h)]
    for i in range(1, 13):                          # Bogen bis an die hintere Leiste
        dy = 23.2 * i / 12
        prof.append((t_kn + dy, -19.7 + math.sqrt(r*r - dy*dy)))
    oben = [(s(t), Z0 + z + T_FOLIE/2) for t, z in prof]
    unten = [(s(t), Z0 + z - T_FOLIE/2) for t, z in reversed(prof)]
    if achse == "x":   # Profil laeuft in X, quer ueber die Breite (Y) gezogen
        w = cq.Workplane("XZ").polyline(oben + unten).close().extrude(-(xb - xa)).translate((0, xa, 0))
    else:
        w = cq.Workplane("YZ").polyline(oben + unten).close().extrude(xb - xa).translate((xa, 0, 0))
    return w.val()

teile = []
# ---------- Knopfreihen der Front: Abstufung je 7 mm ----------
# Reihe 1 (unten, in der welligen Flaeche): 7 mm tief vertieft; Reihe 2: Frontflaeche; Reihe 3: Blende 7 mm aussen.
# Dafuer innen aufgedoppelt: ueber den ganzen Knopfbereich 2 mm, hinter der welligen Flaeche zusaetzlich 5 mm.
# Die Aufdoppelung beginnt 4 mm ueber der Rahmenunterkante = Auflage fuer den 4-mm-Bassdeckel (innen).
RD = ezdxf.readfile("rahmen.dxf").modelspace()
T_DECKEL_BASS, VERT_R1 = 4.0, 7.0
welle = [e for e in RD.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]') if not e.closed][0]
wp_ = [(q.x, q.y) for q in ezdxf.path.make_path(welle).flattening(0.05)]
wp_ = wp_ + [(wp_[-1][0], -4.24 - 1), (wp_[0][0], -4.24 - 1)]
welle_xz = [(-x, y + 4.24) for x, y in wp_]                 # Front ist um 180 Grad gedreht eingebaut
wx0, wx1 = min(a for a, _ in welle_xz), max(a for a, _ in welle_xz)
FI = B - D                                                    # Innenseite der Front
def prisma_xz(pts, y0, y1):
    return cq.Workplane("XZ").polyline(pts).close().extrude(-(y1 - y0)).translate((0, y0, 0)).val()
auf2 = box(wx0, wx1, FI - 2.0, FI, T_DECKEL_BASS, 30.24)
auf5 = prisma_xz(welle_xz, FI - 7.0, FI - 2.0).intersect(box(wx0, wx1, FI - 7.0, FI - 2.0, T_DECKEL_BASS, 60))
for c in RD.query('CIRCLE[layer=="bohrungen-bassdrücker"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        cyl = cq.Workplane("XZ").center(-x, y + 4.24).circle(r).extrude(-20).translate((0, FI - 12, 0)).val()
        auf2 = auf2.cut(cyl); auf5 = auf5.cut(cyl)
# Aufdoppelungen innen an beiden Enden um 5 mm verlaengert (Leimflaeche an der Front)
VERL_AUF = 5.0
for xa_, xb_, xs_ in ((wx0 - VERL_AUF, wx0, wx0 + 0.5), (wx1, wx1 + VERL_AUF, wx1 - 0.5)):
    def _hoehe(sol_, y0_, y1_, x_, std_):
        try: return sol_.intersect(box(x_ - 1.0, x_ + 1.0, y0_, y1_, 0, 80)).BoundingBox().zmax
        except Exception: return std_
    h2_ = 30.24
    h5_ = _hoehe(auf5, FI - 7.0, FI - 2.0, xs_ + (1.0 if xs_ < L / 2 else -1.0), 20.0)
    auf2 = auf2.fuse(box(xa_, xb_, FI - 2.0, FI, T_DECKEL_BASS, h2_))
    auf5 = auf5.fuse(box(xa_, xb_, FI - 7.0, FI - 2.0, T_DECKEL_BASS, h5_))
teile.append(("Front_Aufdoppelung_2mm", auf2))
print("Knopfreihen: Vertiefung Reihe 1", VERT_R1, "mm, Aufdoppelung 2 + 5 mm, X", round(wx0, 1), "-", round(wx1, 1))
# ---------- Auflageleiste fuer den Bassdeckel rundum (4 mm ueber der Unterkante) ----------
# 6 x 6 mm, an Seitenteilen, Rueckteil und Front (an der Front ausserhalb der Aufdoppelung).
# Um die Oeffnungen (Ranken im Rueckteil, Hornbohrungen in der Front) grosszuegig ausgespart (+5 mm),
# damit innen Gaze ueber die Oeffnungen geklebt werden kann.
B_AUF, H_AUF, RAND_G = 6.0, 6.0, 5.0
za0, za1 = T_DECKEL_BASS, T_DECKEL_BASS + H_AUF
aus_rueck, aus_front = [], []
zb_unten = Z_BODEN - T_BODEN                 # Unterseite Boden
REST = 4.0                                   # Rest an der Deckelseite unter den Aussparungen
for e in RD.query('SPLINE[layer=="LAYER_1"]'):
    pp = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.1)]
    if min(b for _, b in pp) > -4.24:
        xs = [a + 382.0 for a, _ in pp]
        aus_rueck.append(box(min(xs) - RAND_G, max(xs) + RAND_G, D - 1, D + B_AUF + 1, za0 + REST, zb_unten + 1))
for c in RD.query('CIRCLE[layer=="Hörner-Bohrung"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        aus_front.append(cq.Workplane("XZ").center(-x, y + 4.24).circle(r + RAND_G).extrude(-(B_AUF + 2)).translate((0, FI - B_AUF - 1, 0)).val())
def mit_aussparung(s, liste):
    for a_ in liste:
        s = s.cut(a_)
    return s
# Seitenteile und Rueckteil: Einsaetze vom Deckel bis zum Boden; Rueckteil ueber die ganze Laenge,
# Aussparungen fuer die Gaze zum Boden hin, an der Deckelseite bleibt ein 4-mm-Rest als Auflage.
auflage = [
    ("Deckelauflage_hinten", mit_aussparung(box(D, L - D, D, D + B_AUF, za0, zb_unten), aus_rueck)),
    ("Deckelauflage_links",  box(D, D + B_AUF, D + B_AUF, FI, za0, zb_unten)),
    ("Deckelauflage_rechts", box(L - D - B_AUF, L - D, D + B_AUF, FI, za0, zb_unten)),
    ("Deckelauflage_vorne_links",  mit_aussparung(box(D + B_AUF, wx0 - VERL_AUF, FI - B_AUF, FI, za0, za1), aus_front)),
    ("Deckelauflage_vorne_rechts", mit_aussparung(box(wx1 + VERL_AUF, L - D - B_AUF, FI - B_AUF, FI, za0, za1), aus_front)),
]
teile.extend(auflage)
# ---------- Bassdeckel (4 mm, innen) zweigeteilt + Eckbloecke ----------
# Hintere Leiste ca. 15 mm breit, Rest des Deckels; in jeder Ecke ein Block 15 x 15 mm vom Deckel bis zum Boden
# (zum Verschrauben des Deckels und fuer die Fuesse). Je Block eine Schraube von unten durch den Deckel.
B_DL, BLK, D_DS, D_BL = 15.0, 15.0, 3.5, 3.0
deckel_h = box(D, L - D, D, D + B_DL, 0.0, T_DECKEL_BASS)
deckel_r = box(D, L - D, D + B_DL, B - D, 0.0, T_DECKEL_BASS)
xi0, xi1 = D + B_AUF, L - D - B_AUF                  # innerhalb der Einsaetze
yi0, yi1 = D + B_AUF, FI - B_AUF
ecken = [(xi0, yi0), (xi1 - BLK, yi0), (xi0, yi1 - BLK), (xi1 - BLK, yi1 - BLK)]
for k, (bx, by) in enumerate(ecken):
    cx, cy = bx + BLK/2, by + BLK/2
    blk_e = box(bx, bx + BLK, by, by + BLK, T_DECKEL_BASS, zb_unten)
    blk_e = blk_e.cut(cq.Workplane("XY").center(cx, cy).circle(D_BL/2).extrude(20).translate((0, 0, T_DECKEL_BASS - 1)).val())
    teile.append((f"Eckblock_{k+1}", blk_e))
    loch_d = cq.Workplane("XY").center(cx, cy).circle(D_DS/2).extrude(T_DECKEL_BASS + 2).translate((0, 0, -1)).val()
    if cy < D + B_DL: deckel_h = deckel_h.cut(loch_d)
    else: deckel_r = deckel_r.cut(loch_d)
# Im Bereich der Knopfreihen endet der Deckel an der Innenseite der Aufdoppelung (kein Ueberstand)
# Deckel reicht bis zur Vorderkante Reihe 1 (Taschenboden); ausgespart nur, wo die Front-Tasche wirklich bis zur Unterkante geht
# Lufttaste und Luftoeffnungen aus Markus-Bass-Verdeck.dxf (Draufsicht wie Bassboden: X = x + 381, Y = y + 206)
VD = ezdxf.readfile("verdeck.dxf").modelspace()
n_lo = 0
for c in VD.query('CIRCLE[layer=="0"]'):
    deckel_r = deckel_r.cut(cq.Workplane("XY").center(c.dxf.center.x + 381, c.dxf.center.y + 206).circle(c.dxf.radius).extrude(T_DECKEL_BASS + 2).translate((0, 0, -1)).val()); n_lo += 1
lt = [e for e in VD.query('LWPOLYLINE[layer=="0"]') if e.closed and 80 < bbox.extents([e], fast=True).size.x < 100][0]
lt_pts = [(q.x + 381, q.y + 206) for q in ezdxf.path.make_path(lt).flattening(0.05)][:-1]
# Lufttasten-Ausnehmung mit halbrunden Enden (Langloch 92 x 14,2), Filzeinlage 1 mm, Taste 90 x 12,2
_ltx = [a_ for a_, _ in lt_pts]; _lty = [b_ for _, b_ in lt_pts]
LT_CX, LT_CY = (min(_ltx) + max(_ltx)) / 2, (min(_lty) + max(_lty)) / 2
LT_DY = -3.0                                       # Lufttaste 3 mm nach hinten versetzt (Platz fuer Knopfreihe 1)
LT_CY += LT_DY
LT_L, LT_B = max(_ltx) - min(_ltx), max(_lty) - min(_lty)
T_FILZ_LT, T_LEDER_LT = 1.0, 0.5
def langloch(l_, b_, z0_, h_):
    return cq.Workplane("XY").center(LT_CX, LT_CY).slot2D(l_, b_).extrude(h_).translate((0, 0, z0_)).val()
deckel_r = deckel_r.cut(langloch(LT_L, LT_B, -1, T_DECKEL_BASS + 2))
lt_pts_alt = lt_pts
print("Bassdeckel: Luftoeffnungen", n_lo, "+ Lufttaste", [round(min(a for a, _ in lt_pts), 1), round(max(a for a, _ in lt_pts), 1)], "| Knopfbereich ausgespart")
# ---------- Bassbodenleder (0,5 mm) auf der Aussenseite des Bassdeckels ----------
# Feld zwischen den beiden Schalloeffnungs-Rosetten (X wie Knopfbereich), vorne bis an die Kante und um die
# Stirnseite herum, hinteres Drittel frei; hintere Ecken gerundet; Leder in die Lufttasten-Oeffnung hineingefuehrt.
T_LEDER_D, R_LEDER = 0.5, 10.0
ly0, ly1 = D + (B - 2*D) / 3.0, FI                          # hinteres Drittel frei .. Deckel-Vorderkante
A_SCHALL = 4.0                                                    # Mindestabstand Leder -> Schalloeffnungen
_lo0 = [(c.dxf.center.x + 381, c.dxf.radius) for c in VD.query('CIRCLE[layer=="0"]')]
lx0 = max(x_ + r_ for x_, r_ in _lo0 if x_ < L / 2) + A_SCHALL
lx1 = min(x_ - r_ for x_, r_ in _lo0 if x_ > L / 2) - A_SCHALL
feld = (cq.Workplane("XY").center((lx0 + lx1)/2, (ly0 + ly1 + 2*R_LEDER)/2).rect(lx1 - lx0, ly1 - ly0 + 2*R_LEDER)
        .extrude(T_DECKEL_BASS).edges("|Z").fillet(R_LEDER - 0.01).val()).intersect(box(lx0, lx1, ly0, ly1, 0, T_DECKEL_BASS))
_lo_ = [c for c in VD.query('CIRCLE[layer=="0"]')]
abst = min(min(abs(c.dxf.center.x + 381 - lx0), abs(c.dxf.center.x + 381 - lx1)) - c.dxf.radius for c in _lo_)
leder_flaeche = feld.intersect(box(lx0, lx1, ly0, ly1, 0, T_LEDER_D))                       # Oberflaeche
leder_schlitz = langloch(LT_L + 2*T_LEDER_D, LT_B + 2*T_LEDER_D, 0, T_DECKEL_BASS).cut(langloch(LT_L, LT_B, -1, T_DECKEL_BASS + 2)).intersect(box(lx0, lx1, ly0, ly1, 0, T_DECKEL_BASS))
leder = leder_flaeche.cut(langloch(LT_L, LT_B, -1, T_DECKEL_BASS + 2)).fuse(leder_schlitz).clean()
deckel_r = deckel_r.cut(leder)
teile.append(("Bassbodenleder", leder))
# Filzeinlage in der Ausnehmung (1 mm, ueber die Deckelstaerke) und Lufttaste (gedrueckt: aussen buendig, 10 mm nach innen)
filz_lt = langloch(LT_L, LT_B, 0, T_DECKEL_BASS).cut(langloch(LT_L - 2*T_FILZ_LT, LT_B - 2*T_FILZ_LT, -1, T_DECKEL_BASS + 2))
H_TASTE, R_TASTE = T_DECKEL_BASS + 10.0, 1.0
taste = cq.Workplane("XY").center(LT_CX, LT_CY).slot2D(LT_L - 2*T_FILZ_LT, LT_B - 2*T_FILZ_LT).extrude(H_TASTE).faces("<Z").edges().fillet(R_TASTE).val().translate(cq.Vector(0, 0, -10.0))   # Ruhelage: 10 mm aussen
teile.append(("Lufttaste_Filz", filz_lt)); teile.append(("Lufttaste", taste))
print("Lufttaste: Langloch", round(LT_L, 1), "x", round(LT_B, 1), "halbrund, Filz", T_FILZ_LT, "| Taste", round(LT_L - 2*T_FILZ_LT, 1), "x", round(LT_B - 2*T_FILZ_LT, 1), "x", H_TASTE, "R", R_TASTE)
print("Bassbodenleder: Feld X", round(lx0, 1), "-", round(lx1, 1), "Y", round(ly0, 1), "-", round(ly1, 1), "| Abstand zu Schalloeffnungen", round(abst, 1), "mm | Ecken R", R_LEDER)
# ---------- Rueckteil: Aussenkante (Deckelseite) zwischen den Schalloeffnungen gerundet, Uebergaenge auslaufend ----------
R_KANTE, L_AUSLAUF = 20.0, 30.0
def kantenprofil(x_, r_):
    r_ = max(r_, 0.05)
    w_ = cq.Workplane("YZ", origin=(x_, 0, 0)).moveTo(-0.01, -0.01).lineTo(r_, -0.01).threePointArc((r_ * (1 - 0.7071), r_ * (1 - 0.7071)), (-0.01, r_)).close()
    return w_.wires().val()
stationen = [(lx0 - L_AUSLAUF, 0.0), (lx0, R_KANTE), (lx1, R_KANTE), (lx1 + L_AUSLAUF, 0.0)]
kante = None
for (xa, ra), (xb, rb) in zip(stationen[:-1], stationen[1:]):
    seg_ = cq.Solid.makeLoft([kantenprofil(xa, ra), kantenprofil(xb, rb)], True)
    kante = seg_ if kante is None else kante.fuse(seg_)
rueck = rueck.cut(kante)
deckel_h = deckel_h.cut(kante)
_gerundet = ["Rueckteil", "Bassdeckel_Leiste_hinten"]
for k_ in range(len(teile)):
    if teile[k_][1].intersect(kante).Volume() > 0.05:
        teile[k_] = (teile[k_][0], teile[k_][1].cut(kante)); _gerundet.append(teile[k_][0])
print("Rundung trifft:", _gerundet)
print("Rueckteil: Kantenrundung R", R_KANTE, "X", round(lx0, 1), "-", round(lx1, 1), "Auslauf je", L_AUSLAUF, "mm")
teile.append(("Bassdeckel_Leiste_hinten", deckel_h))
# Auflage fuer den Deckel hinten (wie Standard-Bassteil 4/5): Leiste 2 dick x 6 breit innen unter der Fuge zwischen
# schmaler hinterer Deckelleiste und Deckel eingeleimt, halb unter der Leiste, halb unter dem Deckel; zwischen den Eckbloecken.
_yfuge = D + B_DL
teile.append(("Bassdeckel_Auflage_hinten", box(xi0 + BLK + 0.5, xi1 - BLK - 0.5, _yfuge - 3.0, _yfuge + 3.0, T_DECKEL_BASS, T_DECKEL_BASS + 2.0)))
print("Deckelauflage hinten 2 x 6 bei Y", round(_yfuge - 3.0, 1), "-", round(_yfuge + 3.0, 1), "X", round(xi0 + BLK + 0.5, 1), "-", round(xi1 - BLK - 0.5, 1))
teile.append(("Bassdeckel", deckel_r))
print("Bassdeckel: hintere Leiste", B_DL, "mm + Rest; Eckbloecke", len(ecken), "je", BLK, "x", BLK, "bis Boden")
print("Deckelauflage rundum: Aussparungen Rueckteil", len(aus_rueck), "Front (Hornbohrungen)", len(aus_front))
# ---------- Basskammern wie Standard-Helikon (4/5-Anordnung): 9 Trennwaende, 8 Kammern 98 lang, Laengsleisten 6 mm ----------
T_TW, B_LL = 3.5, 6.0
TW_X = [x_ + 6.0 for x_ in (7.75, 61.75, 95.75, 149.75, 183.75, 237.75, 271.75, 305.75, 359.75)]
Y_LH0, Y_LH1 = D, D + B_LL
PY0_ = Y_LH1 - 3.0; PY1_ = PY0_ + 98.0
Y_LV0, Y_LV1 = PY1_ - 3.0, PY1_ - 3.0 + B_LL
lh = box(D, L - D, Y_LH0, Y_LH1, Z0, Z1); lv = box(D, L - D, Y_LV0, Y_LV1, Z0, Z1)
for i, x_ in enumerate(TW_X):
    tw_ = box(x_, x_ + T_TW, PY0_, PY1_, Z0, Z1)
    teile.append((f"Basskammer_Trennwand_{i+1}", tw_)); lh = lh.cut(tw_); lv = lv.cut(tw_)
teile += [("Basskammer_Leiste_hinten", lh), ("Basskammer_Leiste_vorne", lv)]
teile.append(("Fuellleiste_Kammer_links", box(D, TW_X[0], Y_LH1, Y_LV0, Z0, Z1)))
teile.append(("Fuellleiste_Kammer_rechts", box(TW_X[-1] + T_TW, L - D, Y_LH1, Y_LV0, Z0, Z1)))
# Gekruemmte Teiler (HB-trennwaende.dxf), 1.5 mm, beidseits 1 mm in die Trennwaende genutet; Fuss vorne (dickes Ende)
TD = ezdxf.readfile("kb_hb_trennwand_detail.dxf").modelspace()
def _pz(e): return [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.05)]
kurve = sorted([p_ for p_ in [_pz(e) for e in TD.query('LWPOLYLINE') if not e.closed] if min(b for _, b in p_) < 10][0])
T_TEIL, T_NUTT = 1.5, 1.0
k_yz = [(PY1_ - (a + 94.88), Z0 + (b - 0.63)) for a, b in kurve]; k_yz[0] = (k_yz[0][0], Z0 + 0.01)
prof_t = k_yz + [(y_, z_ + T_TEIL) for y_, z_ in reversed(k_yz)]
def teiler_xy(x0_, x1_): return cq.Workplane("YZ").polyline(prof_t).close().extrude(x1_ - x0_).translate((x0_, 0, 0)).val()
for i in range(len(TW_X) - 1):
    teile.append((f"Basskammer_Teiler_{i+1}", teiler_xy(TW_X[i] + T_TW - T_NUTT, TW_X[i + 1] + T_NUTT)))
for k_ in [k for k, (n_, _) in enumerate(teile) if n_.startswith("Basskammer_Trennwand_")]:
    s_ = teile[k_][1]
    for i in range(len(TW_X) - 1): s_ = s_.cut(teiler_xy(TW_X[i] + T_TW - T_NUTT - 0.05, TW_X[i + 1] + T_NUTT + 0.05))
    teile[k_] = (teile[k_][0], s_)
# Helikon-Stimmplatten wie Markus (Keil 5 -> 16, Rand 2), doppelt 53 / einfach 33 x 98; Haken hinten und vorne
T_DUENN, T_DICK, T_RAND, RUECK = 5.0, 16.0, 2.0, 2.0
D_HAKEN, L_KURZ = 1.8, 5.0
ZH = Z1 + T_RAND + D_HAKEN / 2
def stab3(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); w_ = p1 - p0
    return cq.Solid.makeCylinder(d / 2, w_.Length, p0, w_.normalized())
def _haken(xc, py_, dy_, l_lang, k_lei, ax="y"):
    if ax == "y":
        hk = stab3((xc, py_, ZH - l_lang), (xc, py_, ZH + 0.3), D_HAKEN).fuse(stab3((xc, py_ - dy_ * 0.3, ZH), (xc, py_ + dy_ * L_KURZ, ZH), D_HAKEN)).clean()
        lo = stab3((xc, py_, ZH - l_lang - 0.5), (xc, py_, Z1 + 0.1), D_HAKEN)
    else:
        hk = stab3((xc, py_, ZH - l_lang), (xc, py_, ZH + 0.3), D_HAKEN).fuse(stab3((xc - dy_ * 0.3, py_, ZH), (xc + dy_ * L_KURZ, py_, ZH), D_HAKEN)).clean()
        lo = stab3((xc, py_, ZH - l_lang - 0.5), (xc, py_, Z1 + 0.1), D_HAKEN)
    teile[k_lei] = (teile[k_lei][0], teile[k_lei][1].cut(lo))
    return hk
_ilh = [k for k, (n_, _) in enumerate(teile) if n_ == "Basskammer_Leiste_hinten"][0]
_ilv = [k for k, (n_, _) in enumerate(teile) if n_ == "Basskammer_Leiste_vorne"][0]
n_pl = 0
for i in range(len(TW_X) - 1):
    PX0, PX1 = TW_X[i] + T_TW / 2 + 0.5, TW_X[i + 1] + T_TW / 2 - 0.5; breit = PX1 - PX0
    pl = box(PX0, PX1, PY0_, PY1_, Z1, Z1 + T_RAND)
    a0, a1, b0, b1 = Y_LH1 + 0.5, Y_LV0 - 0.5, PX0 + RUECK, PX1 - RUECK
    prof = [(a0, Z1), (a1, Z1), (a1, Z1 + T_RAND - T_DICK), (a0, Z1 + T_RAND - T_DUENN)]
    pl = pl.fuse(cq.Workplane("YZ").polyline(prof).close().extrude(b1 - b0).translate((b0, 0, 0)).val())
    if breit > 45: muster = [(5.1, 7.5, 51.0), (15.9, 9.2, 71.5), (27.9, 9.2, 71.5), (40.4, 7.5, 51.0)]
    else: muster = [(breit / 2 - 1.4 - 9.2, 9.2, 71.5), (breit / 2 + 1.4, 9.2, 71.5)]
    for (o_, sb, sll) in muster:
        pl = pl.cut(box(PX0 + o_, PX0 + o_ + sb, PY1_ - 5.5 - sll, PY1_ - 5.5, Z1 - T_DICK - 1, Z1 + T_RAND + 1))
    n_pl += 1; teile.append((f"Stimmplatte_{n_pl}", pl))
    xc = (PX0 + PX1) / 2; off = 0.1 + D_HAKEN / 2
    teile.append((f"Haken_P{n_pl}_hinten", _haken(xc, PY0_ - off, 1, 15.0, _ilh)))
    teile.append((f"Haken_P{n_pl}_vorne", _haken(xc, PY1_ + off, -1, 8.5, _ilv)))
# ---------- Zusatzkammer vorne rechts quer (4/5): an Front (Frontleiste = Vorderwand) und Seitenteil ----------
FY1 = B - D; FY0 = FY1 - T_LEIST
PZX0, PZX1 = XZ0_ + 3.0, XZ1_ - 3.0
PZY1 = FY0 + 1.5; PZY0 = PZY1 - 33.0                 # Keil endet 0.5 vor der Frontleiste
WZ1 = PZY0 - 0.5 + T_TW / 2; WZ0 = WZ1 - T_TW
z_wand = box(XZ0_, XZ1_, WZ0, WZ1, Z0, Z1)
z_la = box(XZ0_, XZ0_ + B_LL, WZ1, FY0, Z0, Z1); z_lb = box(XZ1_ - B_LL, XZ1_, WZ1, FY0, Z0, Z1)
pz = box(PZX0, PZX1, PZY0, PZY1, Z1, Z1 + T_RAND)
a0z, a1z = XZ0_ + B_LL + 0.5, XZ1_ - B_LL - 0.5
profz = [(a0z, Z1), (a1z, Z1), (a1z, Z1 + T_RAND - T_DUENN), (a0z, Z1 + T_RAND - T_DICK)]   # dickes Ende bei X klein (ueber der Klappe)
pz = pz.fuse(cq.Workplane("XZ").polyline(profz).close().extrude(-(PZY1 - PZY0 - 2 * RUECK)).translate((0, PZY0 + RUECK, 0)).val())
for (o_, sb, sll) in [(33.0 / 2 - 1.4 - 9.2, 9.2, 71.5), (33.0 / 2 + 1.4, 9.2, 71.5)]:
    pz = pz.cut(box(PZX0 + 5.5, PZX0 + 5.5 + sll, PZY0 + o_, PZY0 + o_ + sb, Z1 - T_DICK - 1, Z1 + T_RAND + 1))
prof_z = [(PZX0 + (a + 94.88), Z0 + (b - 0.63)) for a, b in kurve]; prof_z[0] = (prof_z[0][0], Z0 + 0.01)
prof_z = prof_z + [(x_, z_ + T_TEIL) for x_, z_ in reversed(prof_z)]
def _tz(y0_, y1_): return cq.Workplane("XZ").polyline(prof_z).close().extrude(-(y1_ - y0_)).translate((0, y0_, 0)).val()
teiler_z = _tz(WZ1 - T_NUTT, FY0 + T_NUTT)
z_wand = z_wand.cut(_tz(WZ1 - T_NUTT - 0.05, FY0 + T_NUTT + 0.05))
teile += [("Zusatzkammer_Wand", z_wand), ("Zusatzkammer_Leiste_a", z_la), ("Zusatzkammer_Leiste_b", z_lb), ("Zusatzkammer_Teiler", teiler_z), ("Stimmplatte_Zusatz", pz)]
_iza = len(teile) - 4; _izb = len(teile) - 3
_ycz2 = (PZY0 + PZY1) / 2
teile.append(("Haken_Zusatz_a", _haken(PZX0 - 0.1 - D_HAKEN / 2, _ycz2, 1, 8.5, _iza, "x")))
teile.append(("Haken_Zusatz_b", _haken(PZX1 + 0.1 + D_HAKEN / 2, _ycz2, -1, 15.0, _izb, "x")))
print("Basskammern Helikon 4/5: Trennwaende", len(TW_X), "Platten", n_pl, "+ Zusatzkammer, Leisten 6 mm, Kammern Y", PY0_, "-", PY1_, "| Zusatzkammer X", XZ0_, "-", XZ1_, "Y", round(WZ0, 1), "-", FY1)
# ---------- Vorderer Bereich (Luftklappe): Laengswand hinter dem Begleiterbereich, Zwischendeckel bis zur Zusatzkammer ----------
T_DECKEL, T_UNTEN = 3.0, 4.0
gy0 = v(0, 136.0)[1] + DYB                      # Vorderkante Begleiterbereich
gx0, gx1 = D, XZ0_
H_FRONT, UEBER = 76.0, 5.0
fl = box(D, L - D, FY0, FY1, Z0, Z1).cut(box(D - 1, XZ0_, FY0 - 1, FY1 + 1, Z0 - 1, H_FRONT - UEBER))   # bei der Zusatzkammer volle Hoehe (Vorderwand)
fl = fl.cut(_tz(WZ1 - T_NUTT - 0.05, FY0 + T_NUTT + 0.05))
teile.append(("Frontleiste", fl))
zd0, zd1 = Z1 - T_DECKEL, Z1
deckel = box(gx0, gx1 + 1.0, gy0, FY0 + T_NUT, zd0, zd1).fuse(box(gx1 - 1.0, L - D, gy0, WZ0 + 1.0, zd0, zd1)).clean()   # auch ueber den Streifen zwischen Begleiter und Zusatzkammer
for k_ in range(len(teile)):
    n_, s_ = teile[k_]
    if n_.startswith(("Stimmplatte", "Haken")): continue
    if s_.intersect(deckel).Volume() > 0.01: teile[k_] = (n_, s_.cut(deckel))
# Luftklappenoeffnung (Standardloesung): Lage aus der Lufthebel-Zeichnung - Hebelauge ueber der Tastenmitte
_LH = ezdxf.readfile("kb_lufthebel.dxf").modelspace()
_ob = sorted([[(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.05)] for e in _LH.query('LWPOLYLINE[layer=="2Füllung-Öffnug"]')], key=lambda p_: max(a for a, _ in p_) - min(a for a, _ in p_))[0]
_lhx = LT_CX + 12.65
DX_LK = -31.1                                                     # Oeffnung wie beim Standard-Bassteil relativ zum Hebeldrehpunkt (weiter vom Drehpunkt weg)
LK_X0, LK_X1 = min(a for a, _ in _ob) + _lhx + DX_LK, max(a for a, _ in _ob) + _lhx + DX_LK
LK_CY = (gy0 + T_UNTEN + 3.0 + 12.5 + 1.0)                       # Klappe 25 breit, 1 mm hinter der Laengswand + Stuetzenplatz
_ohh = max(b for _, b in _ob) - min(b for _, b in _ob)
deckel = deckel.cut(box(LK_X0, LK_X1, LK_CY - _ohh / 2, LK_CY + _ohh / 2, zd0 - 1, zd1 + 1))
print("Luftklappenoeffnung im Deckel (Standard): X", round(LK_X0, 1), "-", round(LK_X1, 1), "Y", round(LK_CY - _ohh / 2, 1), "-", round(LK_CY + _ohh / 2, 1))
teile.append(("Zwischen_Deckel", deckel))
teile.append(("Zwischen_Laengswand", box(gx0, L - D, gy0, gy0 + T_UNTEN, Z0, zd0)))
# Stirnseiten schliessen (Rahmen nur 76 hoch, Deckel bei 82): links ueber die ganze Tiefe des vorderen Bereichs,
# rechts im Streifen zwischen Begleiterbereich und Zusatzkammer
teile.append(("Fuellleiste_Vorne_links", box(D, D + 6.0, gy0 + T_UNTEN, FY0, Z0, zd0)))
teile.append(("Fuellleiste_Vorne_rechts", box(L - D - 6.0, L - D, gy0 + T_UNTEN, WZ0, Z0, zd0)))
boden = boden.cut(box(gx0, gx1, gy0 + T_UNTEN, FY1, Z_BODEN - T_BODEN - 1, Z_BODEN + 1))
boden = boden.cut([s_ for n_, s_ in teile if n_ == "Front_Aufdoppelung_2mm"][0])        # Boden endet an der Aufdoppelung
gy1 = FY0

# ---------- Begleiter: Zwischenteil + Stimmstock ----------
# Boden: je Gruppe eine grosse runde Oeffnung (Ø17). Zwischenteil (Markus.Begleiterzwischenteil.dxf) 304 x 34 mm
# mit 8 Fenstern 30 x 22 (je vier Stimmplatten zusammengefasst), Hoehe 31,2 mm laut Schnitt.
# Stimmstock (Begleiter.dxf): Sole 3,6 mm mit 32 Oeffnungen 14 x 9,5 (eine je Stimmplatte),
# zwei Aussenwaende 3,6 mm und Mittelsteg 2,1 mm, 42,5 mm hoch (Schnitt), Deckleiste oben, Abschluesse an den Enden.
H_ZW, T_SOLE, H_WAND, T_AW, T_MS = 31.2, 3.6, 42.5, 3.6, 2.1
bz = ezdxf.readfile("bzw.dxf").modelspace(); bg = ezdxf.readfile("begl.dxf").modelspace()
zx0, zx1, zy0, zy1 = -309.5, -5.5, 102.0, 136.0                # Zwischenteil/Sole in Bassboden-Koordinaten
(ZX0, ZY0), (ZX1, ZY1) = v(zx0 + DXB, zy0 + DYB), v(zx1 + DXB, zy1 + DYB)
ZY1 = min(ZY1, gy0)                                             # nicht in die Wand der schmalen Kammern
zb0, zb1 = Z_BODEN, Z_BODEN + H_ZW
ZYW = Y_LV1                                                     # Zwischenteil liegt an der vorderen Kammerleiste an                                            # Zwischenteil 1,5 mm breiter: liegt an der Reihenleiste an
zw = cq.Workplane("XY").box(ZX1 - ZX0 + 10.0, ZY1 - ZYW, H_ZW, centered=False).translate((ZX0 - 10.0, ZYW, zb0))   # links 10 mm laenger (zweite Stimmstock-Schraube)
n_f = 0
for e in bz.query('LWPOLYLINE[layer=="LAYER_1"]'):
    pts = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.05)][:-1]
    w_ = max(a for a, _ in pts) - min(a for a, _ in pts); h_ = max(b for _, b in pts) - min(b for _, b in pts)
    if abs(w_ - 30) < 0.5 and abs(h_ - 22) < 0.5:               # Fenster (separat gezeichnet: um -5 / +135,5 verschoben)
        pp = [v(a - 5.0 + DXB, b + 135.5 + DYB) for a, b in pts]
        zw = zw.cut(cq.Workplane("XY").polyline(pp).close().extrude(H_ZW + 2).translate((0, 0, zb0 - 1))); n_f += 1
teile.append(("Begleiter_Zwischenteil", zw.val()))
# offener Bereich am Schraubenende: zwischen Zwischenteil-Ende und Seitenteil bis Kammeroberkante fuellen
teile.append(("Fuellleiste_Begleiter_rechts", box(L - D - 5.0, L - D, ZYW, ZY1, Z0, Z1)))   # beidseits nur 5-mm-Leiste (Begleiter mittig)
teile.append(("Fuellleiste_Begleiter_links", box(D, D + 5.0, ZYW, ZY1, Z0, Z1)))   # nur 5-mm-Leiste: Aussenring geschlossen
# Sole und Stimmstock nach Begleiter.cb (CamBam): Sole 288,4 x 34 x 6 mit 32 Oeffnungen 14 x 9,5;
# Stimmstock aus einem Holzblock 288,4 x 25 x 50, je Seite 16 gleiche Kammern 16 x 45 mm, 9,5 mm tief (Taschen), unten offen zur Sole.
T_SOLE, L_ST, B_ST, H_ST = 6.0, 288.4, 20.0, 50.0       # Block nur 20 mm: Kammern 9,5 tief, Mittelwand 1 mm
L_BLK = 278.4                                              # Block: letzte Trennwand 5 mm; Sole 10 mm laenger (Schraube)
K_B, K_H, K_T, K_PITCH, K_X0 = 16.0, 45.0, 9.5, 17.0, 2.4
s0, s1 = zb1, zb1 + T_SOLE
SX0b = v(-309.5 + DXB, 0)[0]; SX1b = SX0b + L_ST
sole = cq.Workplane("XY").box(L_ST + 10.0, ZY1 - ZY0, T_SOLE, centered=False).translate((SX0b - 10.0, ZY0, s0))   # links 10 mm fuer die zweite Schraube
n_s = 0
# Sole-Oeffnungen so gross wie moeglich = Grundflaeche jeder Kammer (16 x 9,5); werden nach dem Verleimen
# zusammen mit den Kammern gefraest (3 Arbeitsgaenge: Kammern vorne, Kammern hinten, Sole-Oeffnungen)
_by0 = v(0, 102.0 + 17.0)[1] + DYB - B_ST/2
for i in range(16):
    kx = SX0b + K_X0 + i * K_PITCH
    for (ya, yb) in ((_by0, _by0 + K_T), (_by0 + B_ST - K_T, _by0 + B_ST)):
        sole = sole.cut(cq.Workplane("XY").box(K_B, yb - ya, T_SOLE + 2, centered=False).translate((kx, ya, s0 - 1))); n_s += 1
teile.append(("Begleiter_Sole", sole.val()))
by0 = v(0, 102.0 + 17.0)[1] + DYB - B_ST/2                              # Block mittig ueber der Trennung der zwei Oeffnungsreihen
blk = cq.Workplane("XY").box(L_BLK, B_ST, H_ST, centered=False).translate((SX0b, by0, s1))
for i in range(16):
    kx = SX0b + K_X0 + i * K_PITCH
    for (ya, yb) in ((by0 - 1, by0 + K_T), (by0 + B_ST - K_T, by0 + B_ST + 1)):   # Vorder- und Rueckseite
        blk = blk.cut(cq.Workplane("XY").box(K_B, yb - ya, K_H + 1, centered=False).translate((kx, ya, s1 - 1)))
# --- Befestigung des Stimmstocks ---
# Freies Ende (rechts): Nase in Hoehe der Sole verlaengert, mit einer Schraube in das Zwischenteil/den Boden geschraubt.
# Anderes Ende (links): kleine Nase, die von oben in einen Schlitz der Einzelkammer-Wand eingeschoben wird;
# (keine Aufdoppelung: innen ist ueber der Folie kein Platz).
# Freies Ende: keine breitere Nase mehr (kollidiert mit dem Balg). Die Sole steht 10 mm ueber den Block hinaus,
# dort eine Schraube durch die Sole ins Zwischenteil.
D_LOCH, D_SCHR, D_SKOPF, H_SKOPF, L_SCHR = 3.5, 3.0, 6.0, 2.5, 20.0
ymid = by0 + B_ST/2
sx_, sy_ = SX0b + L_BLK + (L_ST - L_BLK)/2, ymid
nz1 = s1
i_so = [n for n, _ in teile].index("Begleiter_Sole")
teile[i_so] = ("Begleiter_Sole", teile[i_so][1].cut(cq.Workplane("XY").center(sx_, sy_).circle(D_LOCH/2).extrude(T_SOLE + 2).translate((0, 0, s0 - 1)).val()))
schr = cq.Workplane("XY").center(sx_, sy_).circle(D_SKOPF/2).extrude(H_SKOPF).translate((0, 0, nz1)).val().fuse(
       cq.Workplane("XY").center(sx_, sy_).circle(D_SCHR/2).extrude(L_SCHR).translate((0, 0, nz1 - L_SCHR)).val())
i_zw = [n for n, _ in teile].index("Begleiter_Zwischenteil")
teile[i_zw] = ("Begleiter_Zwischenteil", teile[i_zw][1].cut(cq.Workplane("XY").center(sx_, sy_).circle(D_SCHR/2).extrude(L_SCHR).translate((0, 0, nz1 - L_SCHR)).val()))
# (Variante 3: keine Einzelkammer -> keine Nase; linkes Ende liegt an der Fuellleiste an)
assert len(blk.solids().vals()) == 1, "Stimmstock nicht aus einem Stueck"
teile.append(("Begleiter_Stimmstock", blk.val()))
teile.append(("Begleiter_Schraube", schr))
# zweite Schraube am linken Ende durch die Sole ins verlaengerte Zwischenteil
sxl_ = SX0b - 5.0
i_so = [n for n, _ in teile].index("Begleiter_Sole")
teile[i_so] = ("Begleiter_Sole", teile[i_so][1].cut(cq.Workplane("XY").center(sxl_, sy_).circle(D_LOCH/2).extrude(T_SOLE + 2).translate((0, 0, s0 - 1)).val()))
schr2 = cq.Workplane("XY").center(sxl_, sy_).circle(D_SKOPF/2).extrude(H_SKOPF).translate((0, 0, s1)).val().fuse(
        cq.Workplane("XY").center(sxl_, sy_).circle(D_SCHR/2).extrude(L_SCHR).translate((0, 0, s1 - L_SCHR)).val())
i_zw2 = [n for n, _ in teile].index("Begleiter_Zwischenteil")
teile[i_zw2] = ("Begleiter_Zwischenteil", teile[i_zw2][1].cut(cq.Workplane("XY").center(sxl_, sy_).circle(D_SCHR/2).extrude(L_SCHR).translate((0, 0, s1 - L_SCHR)).val()))
teile.append(("Begleiter_Schraube_2", schr2))
w1 = s1 + H_ST
print("Begleiter: Fenster im Zwischenteil", n_f, "Sole-Oeffnungen", n_s, "| Oberkante Stimmstock Z", round(w1, 1))

befest = []
# Kollisionskontrolle Kammern untereinander (Nut/Feder sollte 0 ergeben)
k = [s for n, s in teile]
kol = 0
for i in range(len(k)):
    for j in range(i+1, len(k)):
        vol = k[i].intersect(k[j]).Volume()
        if vol > 0.5: kol += 1; print("Ueberschneidung", teile[i][0], teile[j][0], round(vol, 1))
print("Kammerteile", len(teile), "Ueberschneidungen", kol)

# ---------- Hoerner (Horn-schnitt3.dxf, Schalltrichter) in den 6 Hornbohrungen der Front ----------
# Profil: Achse = Linie y=0, a = axial (0..11,5), r = radial (6,3..16,5); Schaft Ø18 x 4 steckt in der Bohrung,
# die Stufe bei a = 4 liegt an der Aussenseite der Front, der Trichter (Ø33) zeigt nach aussen.
HD = ezdxf.readfile("horn.dxf").modelspace()
segs = []
for e in HD.query('LINE ARC'):
    if e.dxftype() == "LINE" and abs(e.dxf.start.y) < 1e-6 and abs(e.dxf.end.y) < 1e-6: continue   # Achse
    p = [(round(q.x, 4), round(q.y, 4)) for q in ezdxf.path.make_path(e).flattening(0.02)]
    segs.append(p)
kette = segs.pop(0)
while segs:
    for i, s in enumerate(segs):
        if math.dist(s[0], kette[-1]) < 0.05: kette += s[1:]; segs.pop(i); break
        if math.dist(s[-1], kette[-1]) < 0.05: kette += s[::-1][1:]; segs.pop(i); break
    else:
        raise SystemExit("Hornprofil nicht geschlossen")
horn_prof = [(a_, r_) for a_, r_ in kette[:-1]] if math.dist(kette[0], kette[-1]) < 0.05 else kette
A_STUFE = 4.01
def horn(x, z, y_aussen):
    # Profil in der XY-Ebene (X = radial, Y = axial), um die Y-Achse drehen
    w = cq.Workplane("XY").polyline([(r_, a_) for a_, r_ in horn_prof]).close().revolve(360, (0, 0, 0), (0, 1, 0))
    return w.val().translate(cq.Vector(x, y_aussen - A_STUFE, z))
hoerner = []
for c in RD.query('CIRCLE[layer=="Hörner-Bohrung"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        X_, Z_ = -x, y + 4.24
        ring = cq.Workplane("XZ").center(X_, Z_).circle(r + 2.5).circle(r).extrude(-40).translate((0, B - 30, 0)).val()
        y_out = front.intersect(ring).BoundingBox().ymax
        hoerner.append((f"Horn_{len(hoerner)+1}", horn(X_, Z_, y_out)))
print("Hoerner", len(hoerner), "Profilpunkte", len(horn_prof), "Ø max", round(2*max(r_ for _, r_ in horn_prof), 1))
# ---------- Fuesse aus fus2.dxf (Drehteil, Halbschnitt): Achse = Linie y = 0, Laenge 21,5, Ø15 ----------
# Profil x = -21,5 liegt am Bassdeckel (Z = 0), x = 0 ist die Standflaeche (Z = -21,5).
# Innen: Senkung Ø13,5 x 2,5 (Einlage), Ø8 fuer den Schraubenkopf, Bohrung Ø4 durchgehend.
# Je Eckblock ein Fuss, mittig unter der Deckelschraube; Schraube Ø3 durch Fuss und Deckel in den Eckblock.
L_FUSS = 21.5
aussen = [(0, 0), (0, 7.0), (-0.5, 7.5), (-3.0, 7.5)]
import math as _m
_c, _r = (-8.889409, 11.283499), 7.0
for i in range(1, 30):                                   # Hohlkehle (Bogen) von -3 nach -14
    t = _m.radians(327.0 - (327.0 - 223.0) * i / 30)
    aussen.append((_c[0] + _r * _m.cos(t), _c[1] + _r * _m.sin(t)))
aussen += [(-14.0, 6.5), (-16.0, 6.5), (-16.0, 7.5), (-21.5, 7.5), (-21.5, 0)]
innen = [(0, 0), (0, 6.75), (-2.5, 6.75), (-2.5, 4.0), (-4.0, 4.0), (-4.0, 4.75), (-5.5, 4.75), (-5.5, 4.0), (-7.5, 2.0), (-21.5, 2.0), (-21.5, 0)]
def drehteil(profil):
    pts = [(r_, -(L_FUSS + x_)) for x_, r_ in profil]          # Ebene XZ: X = radial, Z = Hoehe
    return cq.Workplane("XZ").polyline(pts).close().revolve(360, (0, 0, 0), (0, 1, 0)).val()
fuss_roh = drehteil(aussen).cut(drehteil(innen))
fuesse = []
for k, (bx, by) in enumerate(ecken):
    cx, cy = bx + BLK/2, by + BLK/2
    fuesse.append((f"Fuss_{k+1}", fuss_roh.translate(cq.Vector(cx, cy, 0))))
    # Senkkopfschraube: Kopf im Konus (x = -5,5 .. -7,5, Ø8 -> Ø4), Schaft Ø3 durch Fuss und Deckel in den Eckblock
    z_k0, z_k1 = -(L_FUSS - 5.5), -(L_FUSS - 7.5)
    kopf = cq.Solid.makeCone(4.0, 1.5, z_k1 - z_k0, cq.Vector(cx, cy, z_k0), cq.Vector(0, 0, 1))
    sc = cq.Workplane("XY").center(cx, cy).circle(D_BL/2).extrude(L_FUSS - 7.5 + T_DECKEL_BASS + 12).translate((0, 0, z_k1)).val().fuse(kopf)
    # Gummipfropfen ueber der Schraube: fuellt Senkung Ø13,5 x 2,5, Ø8 und Haltenut Ø9,5 bis x = -5,5
    # steht ca. 1 mm ueber die Standflaeche, dort etwas kleiner (Ø12,5) und nach aussen gewoelbt (Rand 0,4 / Mitte 1,0)
    kuppe = [(0.4 + 0.6 * (1 - (r_ / 6.25) ** 2), r_) for r_ in [6.25 - 6.25 * i / 12 for i in range(13)]]   # vom Rand zur Mitte
    gp = drehteil([(-5.5, 0), (-5.5, 4.75), (-4.0, 4.75), (-4.0, 4.0), (-2.5, 4.0), (-2.5, 6.75), (0, 6.75), (0, 6.25)] + kuppe)
    fuesse.append((f"Fuss_Gummipfropfen_{k+1}", gp.translate(cq.Vector(cx, cy, 0))))
    fuesse.append((f"Fuss_Schraube_{k+1}", sc))
print("Fuesse", len(ecken), "Ø15 x", L_FUSS, "aus fus2.dxf | Volumen", round(fuss_roh.Volume()))
# Kollisionskontrolle der neuen Teile gegen alles
alle = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile
for n1, s1 in hoerner + fuesse:
    for n2, s2 in alle:
        vv = s1.intersect(s2).Volume()
        if vv > 0.5: print("Ueberschneidung", n1, n2, round(vv, 1))

# ---------- Fliegengitter (Gaze) innen ueber Rosetten (Deckel), F-Loechern (Rueckteil) und Hornbohrungen (Front) ----------
# Gewebe: Faeden 0,2 x 0,2 im Raster 2 mm, Kett-/Schussfaeden versetzt, 0,3 dick; je 3 mm groesser als die Oeffnung (Kleberand).
# Angrenzende feste Teile bekommen eine 0,3-mm-Aussparung.
# vereinfacht: Platte 0,3 mit quadratischen Durchbruechen, Raster 4 mm / Steg 0,4 (echt: Raster 2 mm / Faden 0,2)
T_GZ, R_GZ, RASTER, D_FADEN = 0.3, 3.0, 4.0, 0.4
from shapely.geometry import Polygon as _Poly, Point as _Pt
from gaze2d import gaze_platte
_huellen = []
_ros = sorted({(round(c.dxf.center.x + 381, 2), round(c.dxf.center.y + 206, 2), round(c.dxf.radius, 2)) for c in VD.query('CIRCLE[layer=="0"]')})
for j_, (mx_, my_, _) in enumerate(sorted([r_ for r_ in _ros if r_[2] > 6.5])):
    rr_ = max(math.hypot(x_ - mx_, y_ - my_) + r_ for x_, y_, r_ in _ros if math.hypot(x_ - mx_, y_ - my_) < 30) + R_GZ
    _huellen.append((f"Gaze_Rosette_{j_+1}", cq.Solid.makeCylinder(rr_, T_GZ, cq.Vector(mx_, my_, T_DECKEL_BASS)), (_Pt(mx_, my_).buffer(rr_, 32), "XY", T_DECKEL_BASS)))
_jf = 0
for e in RD.query('SPLINE[layer=="LAYER_1"]'):
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.2)]
    if min(b for _, b in pp_) > -4.24:
        pg_ = _Poly([(a + 382.0, b + 4.24) for a, b in pp_]).buffer(R_GZ).simplify(0.05); _jf += 1
        _huellen.append((f"Gaze_FLoch_{_jf}", cq.Workplane("XZ").polyline(list(pg_.exterior.coords)[:-1]).close().extrude(-T_GZ).translate((0, D, 0)).val(), (pg_, "XZ", D)))
_jh = 0
for c in RD.query('CIRCLE[layer=="Hörner-Bohrung"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if -382 < x < 0 and -4.24 < y < 71.75:
        _jh += 1
        _huellen.append((f"Gaze_Horn_{_jh}", cq.Solid.makeCylinder(r + R_GZ, T_GZ, cq.Vector(-x, FI - T_GZ, y + 4.24), cq.Vector(0, 1, 0)), (_Pt(-x, y + 4.24).buffer(r + R_GZ, 32), "XZ", FI - T_GZ)))
_aus = set()
for k_ in range(len(teile)):
    n_, s_ = teile[k_]
    if n_ == "Bassdeckel": continue
    for gn_, gs_, _ in _huellen:
        if s_.intersect(gs_).Volume() > 1e-3: s_ = s_.cut(gs_); _aus.add(n_)
    teile[k_] = (n_, s_)
for gn_, gs_, _ in _huellen:
    if boden.intersect(gs_).Volume() > 1e-3: boden = boden.cut(gs_); _aus.add("Bassboden")
teile += [(gn_, gaze_platte(g2_[0], g2_[1], g2_[2], T_GZ, RASTER, D_FADEN)[0]) for gn_, gs_, g2_ in _huellen]
print("Gaze:", len(_huellen), "Stueck (Rosetten, F-Loecher, Hoerner), Aussparungen in", sorted(_aus))

exec(open("mechanik.py", encoding="utf-8").read())
a = cq.Assembly(name="Bassteil_Markus")
# Links/rechts getauscht (auf Johanns Wunsch): alles an der Ebene X = L/2 spiegeln
SPIEGELN = True
_add_orig = a.add
def _add_gespiegelt(obj, *args, **kw):
    if SPIEGELN and isinstance(obj, cq.Shape):
        obj = obj.mirror("YZ", basePointVector=(L / 2, 0, 0))
    # Farbgebung wie Standard-Bassteil 4/5 (kugel_bass2.py, FARBEN/farbe), ergaenzt um die Teile dieser Variante
    if isinstance(obj, cq.Shape) and kw.get("name"):      # nur beim eigentlichen Teil (cadquery ruft add intern nochmals ohne Namen auf)
        kw["color"] = cq.Color(*farbe_std(kw["name"]))
    return _add_orig(obj, *args, **kw)
FARBEN_STD = [   # (Praefix, Farbe) - wie Standardmodell, erste passende Zeile gilt
    ("Gaze", (0.78, 0.80, 0.82)), ("Frontblende", (0.45, 0.25, 0.12)), ("Blende", (0.45, 0.25, 0.12)), ("Frontleiste", (0.60, 0.45, 0.28)), ("Front", (0.78, 0.60, 0.38)), ("Rueck", (0.74, 0.56, 0.36)), ("Seiten", (0.70, 0.52, 0.33)),
    ("Horn", (0.70, 0.55, 0.30)), ("Bassbodenleder", (0.12, 0.08, 0.05)),
    ("Bassboden", (0.90, 0.82, 0.62)), ("Deckelauflage", (0.60, 0.45, 0.28)), ("Eckblock", (0.55, 0.40, 0.25)), ("Lagerblock", (0.55, 0.40, 0.25)), ("Fuss", (0.28, 0.20, 0.14)),
    ("Fuellleiste", (0.60, 0.45, 0.28)), ("Zwischen_Laengswand", (0.60, 0.45, 0.28)), ("Zwischen_Deckel", (0.84, 0.86, 0.60)), ("Begleiter_Zwischenteil", (0.84, 0.86, 0.60)),
    ("Bassdeckel_Leder", (0.12, 0.08, 0.05)), ("Bassdeckel_Auflage", (0.60, 0.45, 0.28)), ("Bassdeckel", (0.30, 0.22, 0.15)),
    ("Haken", (0.40, 0.40, 0.42)), ("Zusatzkammer_Teiler", (0.95, 0.80, 0.45)), ("Zusatzkammer", (0.50, 0.62, 0.80)), ("Stimmplatte_", (0.80, 0.82, 0.86)),
    ("Basskammer_Teiler", (0.95, 0.80, 0.45)), ("Basskammer_Trennwand", (0.65, 0.75, 0.88)), ("Basskammer", (0.50, 0.62, 0.80)),
    ("Begleiter_Stimmplatte", (0.80, 0.82, 0.86)), ("Begleiter_Sole", (0.55, 0.68, 0.45)), ("Begleiter_Stimmstock", (0.62, 0.72, 0.55)), ("Begleiterstimmstock", (0.62, 0.72, 0.55)),
    ("Lufttaste_Filz", (0.85, 0.15, 0.15)), ("Lufttaste_Gabelgelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)),
    ("Lufttaste_Schubstange", (0.10, 0.60, 0.60)), ("Lufttaste", (0.35, 0.25, 0.15)),
    ("Lufthebel_Feder", (0.75, 0.10, 0.75)), ("Lufthebel_Stuetze_Schraube", (0.40, 0.40, 0.42)), ("Lufthebel_Stuetze", (0.25, 0.75, 0.35)), ("Lufthebel_Bolzen", (0.95, 0.85, 0.20)), ("Lufthebel", (0.70, 0.72, 0.78)),
    ("Knopfzylinder", (0.95, 0.95, 0.92)), ("Knopf", (0.12, 0.12, 0.12)), ("Filzeinlage", (0.85, 0.15, 0.15)), ("Filzring", (0.85, 0.15, 0.15)),
    ("Achse", (0.55, 0.55, 0.60)), ("Schenkel_Knopf", (0.20, 0.40, 0.85)), ("Schenkel_Klappe", (0.10, 0.65, 0.30)), ("Klappenhebel", (0.10, 0.65, 0.30)),
    ("Druckstange", (0.10, 0.60, 0.60)), ("Schubstange", (0.10, 0.60, 0.60)), ("Gelenkgabel", (0.95, 0.45, 0.10)), ("Gelenkkopf", (0.60, 0.20, 0.70)),
    ("Gelenkbolzen", (0.95, 0.85, 0.20)), ("Feder", (0.75, 0.10, 0.75)),
    ("Lagerschelle", (0.25, 0.75, 0.35)), ("Lagerbuchse", (0.25, 0.75, 0.35)),
]
def farbe_std(n):
    if n.startswith("Stimmplatte") and ("_Haken" in n or n.endswith("_Schraube")): return (0.40, 0.40, 0.42)
    if n.startswith("Stimmplatte"): return (0.80, 0.82, 0.86)
    if n.endswith("_Filz"): return (0.85, 0.15, 0.15)              # Klappenfilz rot
    if n.endswith("_Schraube") or "_Schraube_" in n: return (0.40, 0.40, 0.42)
    if n.endswith("_Gummipfropfen"): return (0.08, 0.08, 0.08)
    if "klappe" in n.lower(): return (0.22, 0.22, 0.25)
    for p_, c_ in FARBEN_STD:
        if n.startswith(p_): return c_
    print("FARBE ohne Zuordnung:", n)
    return (0.60, 0.60, 0.63)
a.add = _add_gespiegelt
farben = {"Luftklappe": (0.22, 0.22, 0.25), "Haken": (0.40, 0.40, 0.42), "Teiler": (0.95, 0.80, 0.45), "Zusatzkammer": (0.50, 0.62, 0.80), "Gaze": (0.78, 0.80, 0.82), "Rahmen": (0.78, 0.60, 0.38), "Boden": (0.88, 0.80, 0.62), "Leiste": (0.55, 0.40, 0.25),
          "Trennwand": (0.65, 0.75, 0.88), "Folie": (0.85, 0.30, 0.30), "Stimmplatte": (0.80, 0.82, 0.86), "Zwischen": (0.70, 0.55, 0.35), "Frontleiste": (0.55, 0.40, 0.25), "Begleiter": (0.62, 0.72, 0.55), "Fuellleiste": (0.55, 0.40, 0.25), "Front_Aufdoppelung": (0.78, 0.60, 0.38), "Deckelauflage": (0.60, 0.45, 0.28), "Eckblock": (0.55, 0.40, 0.25), "Bassdeckel": (0.30, 0.22, 0.15), "Bassbodenleder": (0.12, 0.08, 0.05), "Lufttaste_Filz": (0.85, 0.85, 0.80), "Lufttaste": (0.35, 0.25, 0.15), "Aufdoppelung": (0.65, 0.75, 0.88)}
for n, s in [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr)]:
    a.add(s, name=n, color=cq.Color(*farben["Rahmen"]))
a.add(boden, name="Bassboden", color=cq.Color(*farben["Boden"]))
for n, s in befest:
    a.add(s, name=n, color=cq.Color(0.25, 0.25, 0.28))
for n, s in teile:
    f = next((c for key, c in farben.items() if key in n), (0.6, 0.6, 0.6))
    a.add(s, name=n, color=cq.Color(*f))
# Frontblende (Markus-20-mm-Blende.dxf): Aufsatz aussen auf der Front ueber der kurzen Bassdruecker-Reihe,
# 4 Bohrungen Ø10 deckungsgleich mit der Front, 7 mm hoch
# Frontblende ueber Reihe 2 (wie Markus 3. Reihe, sinngemaess auf 8 Bohrungen erweitert, geschwungene Kontur), 7 mm
T_BL, Z_REIHE = 7.0, 24.0 + 4.24
def _blende_kontur(xs_):
    xs_ = sorted(xs_, reverse=True); o_, u_ = 6.66, -4.42; pts_ = []
    pts_.append((xs_[0] - 5.06, o_, -0.301))
    for h_ in xs_[1:-1]: pts_ += [(h_ + 5.06, o_, 0.388), (h_ - 5.06, o_, -0.301)]
    pts_.append((xs_[-1] + 5.06, o_, 2.267)); pts_.append((xs_[-1] + 5.06, u_, -0.301))
    for h_ in xs_[-2:0:-1]: pts_ += [(h_ - 5.06, u_, 0.388), (h_ + 5.06, u_, -0.301)]
    pts_.append((xs_[0] - 5.06, u_, 2.267))
    e_ = ezdxf.new().modelspace().add_lwpolyline([(x_, y_, 0, 0, b_) for x_, y_, b_ in pts_], format="xyseb", close=True)
    return [(q.x, q.y) for q in ezdxf.path.make_path(e_).flattening(0.05)][:-1]
_x2 = sorted(-c.dxf.center.x for c in RD.query('CIRCLE[layer=="bohrungen-bassdrücker"]') if abs(c.dxf.center.y - 24.0) < 0.1 and -382 < c.dxf.center.x < 0)
_bk = [(-x_, Z_REIHE + y_) for x_, y_ in _blende_kontur([-x_ for x_ in _x2])]
bl = cq.Workplane("XZ").polyline(_bk).close().extrude(-T_BL).translate((0, B, 0))
for x_ in _x2:
    bl = bl.cut(cq.Workplane("XZ").center(x_, Z_REIHE).circle(5.0).extrude(-T_BL - 2).translate((0, B - 1, 0)))
for n_, s_ in mech:
    if n_.startswith(("Knopfzylinder", "Schubstange")) and s_.intersect(bl.val()).Volume() > 0.05: print("Blende kollidiert mit", n_)
# Kontrolle Knoepfe/Filz zur Abstufung: Material entlang der Knopfachse (Ring r 5..8 um die Bohrung)
_stack = [("Front", front), ("Blende", bl.val())] + [(n_, s_) for n_, s_ in teile if n_.startswith("Front_Aufdoppelung")]
for n_, s_ in mech:
    if not n_.startswith("Knopf_"): continue
    bbk = s_.BoundingBox(); Xk_, Zk_ = (bbk.xmin + bbk.xmax) / 2, (bbk.zmin + bbk.zmax) / 2
    ring = cq.Workplane("XZ").center(Xk_, Zk_).circle(8.0).circle(5.2).extrude(-60).translate((0, 160, 0)).val()
    lagen_ = []
    for m_, t_ in _stack:
        iv = t_.intersect(ring)
        if iv.Volume() > 0.5: b_ = iv.BoundingBox(); lagen_.append((m_, round(b_.ymin, 1), round(b_.ymax, 1)))
    aussen_ = max(l_[2] for l_ in lagen_); innen_ = min(l_[1] for l_ in lagen_)
    fz_ = [t2 for n2, t2 in mech if n2 == n_.replace("Knopf_", "Filzring_")][0].BoundingBox()
    fl_ok = all(t_.intersect(fz_ and [t2 for n2, t2 in mech if n2 == n_.replace("Knopf_", "Filzring_")][0]).Volume() < 0.05 for m_, t_ in _stack)
    print("FILZ", n_, "Y", round(fz_.ymin, 1), "-", round(fz_.ymax, 1), "frei" if fl_ok else "KOLLISION")
    print("KNOPF", n_, "X", round(Xk_, 1), "Z", round(Zk_, 1), "| Material Y", innen_, "-", aussen_, lagen_, "| Kopf ab Y", round(bbk.ymin, 1), "Abstand", round(bbk.ymin - aussen_, 1))
a.add(bl.val(), name="Frontblende", color=cq.Color(0.70, 0.52, 0.30))
for n, s in mech:
    a.add(s, name=n, color=cq.Color(0.15, 0.15, 0.15) if n.startswith("Knopf_") else cq.Color(0.95, 0.95, 0.92) if n.startswith("Knopfzyl") else cq.Color(0.55, 0.55, 0.20) if n.startswith("Filz") else cq.Color(0.20, 0.45, 0.75) if n.startswith("Gabel") else cq.Color(0.80, 0.20, 0.20) if n.startswith("Gelenkbolzen") else cq.Color(0.85, 0.85, 0.80) if n.startswith("Lagerbock_") and "Schraube" not in n else cq.Color(0.60, 0.60, 0.63))
for n, s in hoerner:
    a.add(s, name=n, color=cq.Color(0.05, 0.05, 0.05))
for n, s in fuesse:
    a.add(s, name=n, color=cq.Color(0.08, 0.08, 0.08) if "Gummi" in n else cq.Color(0.45, 0.30, 0.18) if "Schraube" not in n else cq.Color(0.6, 0.6, 0.62))
for n, s in klappen:
    col = (0.45, 0.30, 0.20) if "Leder" in n else (0.85, 0.85, 0.80) if "Filz" in n else (0.20, 0.20, 0.22)
    a.add(s, name=n, color=cq.Color(*col))
_bb = bl.val().BoundingBox(); print("Frontblende X", round(_bb.xmin, 1), "-", round(_bb.xmax, 1), "Z", round(_bb.zmin, 1), "-", round(_bb.zmax, 1), "Y", round(_bb.ymin, 1), "-", round(_bb.ymax, 1))
# Kontrolle Klappen untereinander
kk = [s for _, s in klappen if True]
_n = 0
for i in range(0, len(klappen), 3):
    for j in range(i + 3, len(klappen), 3):
        if klappen[i+2][1].intersect(klappen[j+2][1]).Volume() > 0.1: _n += 1
print("Klappen ueberschneiden sich:", _n)
a.save(sys.argv[2])
# STL fuer Vorschau
import os
os.makedirs("stl_k", exist_ok=True)
for n, s in [("Rahmen", front.fuse(rueck).fuse(sl).fuse(sr)), ("Boden", boden)] + teile:
    cq.exporters.export(cq.Workplane().add(s), f"stl_k/{n}.stl")
