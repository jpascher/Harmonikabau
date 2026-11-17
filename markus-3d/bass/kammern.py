# Bassteil Markus: Rahmen (v2) + Boden aus Bassboden2.dxf + Basskammern aus Einzelteilen
# X = Laenge (382), Y = Breite (207), Z = Hoehe; Boden-Oberkante Z_BODEN, Kammern stehen darauf (Z steigt Richtung Balg)
# Aufruf: python3 kammern.py Bassboden2.dxf ausgabe.step
import sys, math, ezdxf, ezdxf.path, cadquery as cq
from ezdxf import bbox
D, L, B = 5.5, 382.0, 207.0
Z_BODEN, T_BODEN = 31.2, 6.0
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
for e in ents:
    if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "2Füllung-Öffnug" and e.closed:
        pts = [v(p.x, p.y) for p in ezdxf.path.make_path(e).flattening(0.05)][:-1]
        loch(cq.Workplane("XY").polyline(pts).close()); n_rect += 1
        xs = [a for a, _ in pts]; ys = [b for _, b in pts]
        x0_, x1_, y0_, y1_ = min(xs), max(xs), min(ys), max(ys)
        vertiefung_und_klappe(lambda o, a=x0_, b=x1_, c=y0_, d=y1_: cq.Workplane("XY").center((a+b)/2, (c+d)/2).rect(b-a+2*o, d-c+2*o), f"Klappe_{n_rect}")
    elif e.dxftype() == "CIRCLE" and e.dxf.layer == "Begleiterstimmstock-fuellung-öffnungen":   # Begleiter: je Gruppe eine grosse Oeffnung
        c_ = v(e.dxf.center.x, e.dxf.center.y); r_ = e.dxf.radius
        loch(cq.Workplane("XY").center(*c_).circle(r_)); n_rund += 1
        vertiefung_und_klappe(lambda o, c=c_, r=r_: cq.Workplane("XY").center(*c).circle(r + o), f"Begleiterklappe_{n_rund}")
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
tasche7 = prisma_xz(welle_xz, B - VERT_R1, B + 1)
front = front.cut(tasche7)
auf2 = auf2.cut(tasche7)
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
teile.append(("Front_Aufdoppelung_5mm", auf5))
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
kerbe = tasche7.intersect(box(wx0, wx1, B - VERT_R1, B + 1, -1, T_DECKEL_BASS + 1))
deckel_r = deckel_r.cut(kerbe)
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
leder_stirn = kerbe.translate(cq.Vector(0, -T_LEDER_D, 0)).cut(kerbe).intersect(box(lx0, lx1, 0, B, 0, T_DECKEL_BASS))   # um die Vorderkante (Taschenrand)
leder_schlitz = langloch(LT_L + 2*T_LEDER_D, LT_B + 2*T_LEDER_D, 0, T_DECKEL_BASS).cut(langloch(LT_L, LT_B, -1, T_DECKEL_BASS + 2)).intersect(box(lx0, lx1, ly0, ly1, 0, T_DECKEL_BASS))
leder = leder_flaeche.cut(kerbe).cut(langloch(LT_L, LT_B, -1, T_DECKEL_BASS + 2)).fuse(leder_stirn).fuse(leder_schlitz).clean()
deckel_r = deckel_r.cut(leder)
teile.append(("Bassbodenleder", leder))
# Filzeinlage in der Ausnehmung (1 mm, ueber die Deckelstaerke) und Lufttaste (gedrueckt: aussen buendig, 10 mm nach innen)
filz_lt = langloch(LT_L, LT_B, 0, T_DECKEL_BASS).cut(langloch(LT_L - 2*T_FILZ_LT, LT_B - 2*T_FILZ_LT, -1, T_DECKEL_BASS + 2))
H_TASTE, R_TASTE = T_DECKEL_BASS + 10.0, 1.0
taste = cq.Workplane("XY").center(LT_CX, LT_CY).slot2D(LT_L - 2*T_FILZ_LT, LT_B - 2*T_FILZ_LT).extrude(H_TASTE).faces("<Z").edges().fillet(R_TASTE).val()
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
teile.append(("Bassdeckel", deckel_r))
print("Bassdeckel: hintere Leiste", B_DL, "mm + Rest; Eckbloecke", len(ecken), "je", BLK, "x", BLK, "bis Boden")
print("Deckelauflage rundum: Aussparungen Rueckteil", len(aus_rueck), "Front (Hornbohrungen)", len(aus_front))
# Kammer-Innenmasse laut Bassboden2.dxf (HB-Rahmen innen): Reihe 49 x 90, Einzelkammer, schmale Kammern 91 x 43.
# Die Leisten liegen ausserhalb dieser Innenmasse; Trennwaende stecken 2 mm in den Nuten.
def nuten_dazu(s, xs, y0, y1):
    for nx in xs:
        s = s.cut(box(nx - T_WAND/2, nx + T_WAND/2, y0, y1, Z0 - 1, Z1 + 1))
    return s
# Reihe: 8 Trennwaende, 7 Kammern 53/53/53/46/53/46/53
wx = [v(x, 0)[0] for x in (-364.5, -311.5, -258.5, -205.5, -159.5, -106.5, -60.5, -7.5)]
ri0, ri1 = v(0, 6.5)[1], v(0, 96.5)[1]               # Innenraum in Y
lx0, lx1 = D, L - D                                   # Reihenleisten bis an die Seitenteile
teile.append(("Reihe_Leiste_vorne", leiste(lx0, lx1, ri0 - T_LEIST, ri0, nuten_x=wx, innen="+")))
ex = [v(-364.5, 0)[0], v(-311.5, 0)[0]]
T_MITTE = v(0, 100.5)[1] - ri1        # gemeinsame Leiste Reihe/Einzelkammer: 4 mm (Innenmasse 96,5 / 100,5 laut Zeichnung)
T_NUT1 = 1.0                           # flache Nut, wo die Wand an einer duennen Gegenwand anstoesst
hint = leiste(lx0, lx1, ri1, ri1 + T_MITTE, nuten_x=wx, innen="-")
hint = nuten_dazu(hint, ex, ri1 + T_MITTE - T_NUT1, ri1 + T_MITTE + 1)   # 1-mm-Nuten fuer die Einzelkammer-Waende
teile.append(("Reihe_Leiste_hinten", hint))
for i, x in enumerate(wx):
    teile.append((f"Reihe_Trennwand_{i+1}", box(x - T_WAND/2, x + T_WAND/2, ri0 - T_NUT, ri1 + T_NUT, Z0, Z1)))
for i in range(7):
    teile.append((f"Reihe_Folie_{i+1}", folie_y(wx[i] + T_WAND/2, wx[i+1] - T_WAND/2, ri0, ri1)))
# Ausgleichsbereiche an den Enden (zwischen aeusserer Trennwand und Seitenteil) mit Fuellleisten schliessen,
# damit oben ein gleichmaessig umlaufender Rand entsteht
teile.append(("Fuellleiste_Reihe_links", box(D, wx[0] - T_WAND/2, ri0, ri1, Z0, Z1)))
teile.append(("Fuellleiste_Reihe_rechts", box(wx[-1] + T_WAND/2, L - D, ri0, ri1, Z0, Z1)))
# Einzelkammer: vorne die hintere Reihenleiste, hinten eigene Leiste bis an die Front
ei0, ei1 = ri1 + T_MITTE, v(0, 190.5)[1]
teile.append(("Einzel_Leiste_hinten", leiste(ex[0] - 6, ex[1] + T_WAND/2, ei1, ei1 + T_LEIST, nuten_x=ex, innen="-")))
for i, x in enumerate(ex):
    teile.append((f"Einzel_Trennwand_{i+1}", box(x - T_WAND/2, x + T_WAND/2, ei0 - T_NUT1, ei1 + T_NUT, Z0, Z1)))
teile.append(("Fuellleiste_Einzel_links", box(D, ex[0] - T_WAND/2, ei0, ei1, Z0, Z1)))
_i = [n for n, _ in teile].index("Einzel_Trennwand_2")
_w = teile[_i][1]
for _y in (v(0, 139.5)[1] - T_WAND/2, v(0, 182.5)[1] + T_WAND/2):
    _w = _w.cut(box(ex[1] + T_WAND/2 - T_NUT1, ex[1] + T_WAND/2 + 1, _y - T_WAND/2, _y + T_WAND/2, Z0 - 1, Z1 + 1))
teile[_i] = ("Einzel_Trennwand_2", _w)
teile.append(("Einzel_Folie", folie_y(ex[0] + T_WAND/2, ex[1] - T_WAND/2, ei0, ei1, umgedreht=True)))
# Zwei schmale Kammern (innen 91 x 43, quer): Waende in X-Richtung, Leisten an den X-Enden
si0, si1 = v(0, 139.5)[1], v(0, 182.5)[1]
wy = [si0 - T_WAND/2, si1 + T_WAND/2]
for k, (ix0, ix1) in enumerate([(-309.9, -218.9), (-97.0, -6.0)]):
    I0, I1 = v(ix0, 0)[0], v(ix1, 0)[0]
    ly0, ly1 = wy[0] - T_WAND/2, wy[1] + T_WAND/2   # Leisten buendig mit den Waenden (Ueberstaende beseitigt)
    if k == 0:
        I0 = ex[1] + T_WAND/2            # links schliesst die Wand der Einzelkammer ab (mit 1-mm-Nuten)
        w0 = I0 - T_NUT1
    else:
        teile.append((f"Schmal{k+1}_Leiste_links", leiste(I0 - T_LEIST, I0, ly0, ly1, nuten_y=wy, innen="+")))
        w0 = I0 - T_NUT
    r1 = min(I1 + T_LEIST, L - D)        # nicht in das Seitenteil
    teile.append((f"Schmal{k+1}_Leiste_rechts", leiste(I1, r1, ly0, ly1, nuten_y=wy, innen="-")))
    for i, y in enumerate(wy):
        teile.append((f"Schmal{k+1}_Trennwand_{i+1}", box(w0, I1 + T_NUT, y - T_WAND/2, y + T_WAND/2, Z0, Z1)))
    teile.append((f"Schmal{k+1}_Folie", folie_y(si0, si1, I0, I1, umgedreht=True, achse="x")))

# ---------- Raum zwischen und oberhalb der beiden schmalen Kammern schliessen ----------
# Frontleiste 6 mm ueber die ganze Laenge an der Front (zugleich hintere Leiste der Einzelkammer).
# Eingeschlossen wird der Raum zwischen den beiden schmalen Kammern UND der Streifen oberhalb beider Kammern
# bis zur Frontleiste. Unten zwischen den Kammern 4-mm-Laengswand (buendig mit den Kammerwaenden).
# Deckel 3 mm, oben buendig, steckt in Nuten aller angrenzenden Teile (2 mm in Leisten, 1 mm in duennen Waenden).
# Unter dem abgedeckelten Bereich ist der Boden ausgenommen.
T_DECKEL, T_UNTEN = 3.0, 4.0
names = [n for n, _ in teile]
i1 = names.index("Schmal1_Leiste_rechts"); i2 = names.index("Schmal2_Leiste_links")
b1, b2 = teile[i1][1].BoundingBox(), teile[i2][1].BoundingBox()
gx0, gx1 = b1.xmax, b2.xmin
gy0, gy1 = wy[0] - T_WAND/2, wy[1] + T_WAND/2       # Aussenseiten der Kammerwaende unten / oben
FY1 = B - D; FY0 = FY1 - T_LEIST
SX0 = ex[1] + T_WAND/2                             # links: Einzelkammer-Wand
SX1 = L - D                                        # rechts: Seitenteil
# Einzelkammer-Hinterleiste durch durchgehende Frontleiste ersetzen
k_e = names.index("Einzel_Leiste_hinten"); teile.pop(k_e); names.pop(k_e)
i1 = names.index("Schmal1_Leiste_rechts"); i2 = names.index("Schmal2_Leiste_links")
fl = leiste(D, L - D, FY0, FY1, nuten_x=ex, innen="-")
# Frontleiste reicht nicht bis zum Boden: nur 5 mm Ueberstand nach unten zum Verleimen mit der Front.
# Ausnahme: im Bereich der Einzelkammer bleibt sie voll hoch (dort ist sie die Rueckwand der Kammer).
H_FRONT, UEBER = 76.0, 5.0
fl = fl.cut(box(SX0, L, FY0 - 1, FY1 + 1, Z0 - 1, H_FRONT - UEBER))
# Nuten fuer die Laengswand unten in den beiden zugewandten Kammerleisten
for k_, (xa, xb) in ((i1, (gx0 - T_NUT, gx0 + 1)), (i2, (gx1 - 1, gx1 + T_NUT))):
    teile[k_] = (names[k_], teile[k_][1].cut(box(xa, xb, gy0 - 1, gy0 + T_UNTEN, Z0 - 1, Z1 + 1)))
teile.append(("Frontleiste", fl))
# Deckel: Zwischenraum + Streifen oberhalb der Kammern, mit Einstand in die Nachbarteile
zd0, zd1 = Z1 - T_DECKEL, Z1
deckel = box(gx0 - T_NUT, gx1 + T_NUT, gy0, gy1, zd0, zd1).fuse(
         box(SX0 - 1.0, SX1, gy1 - 1.0, FY0 + T_NUT, zd0, zd1))
# Nuten = Deckelvolumen aus allen beruehrten Teilen ausschneiden
for k_ in range(len(teile)):
    n_, s_ = teile[k_]
    if n_.startswith(("Stimmplatte", "Haken", "Schraube")): continue
    if s_.intersect(deckel).Volume() > 0.01:
        teile[k_] = (n_, s_.cut(deckel))
# Luftklappenoeffnung im mittleren Teil des Deckels (Kontur "Füllung-luftklappe" aus der Hauptzeichnung)
mm = ezdxf.readfile("master.dxf").modelspace()
for e in mm.query('LWPOLYLINE[layer=="Füllung-luftklappe"]'):
    pts = [(p_.x - 2923.3, p_.y - 1520.7) for p_ in ezdxf.path.make_path(e).flattening(0.05)][:-1]
    if gx0 < min(a for a, _ in pts) and max(a for a, _ in pts) < gx1:
        deckel = deckel.cut(cq.Workplane("XY").polyline(pts).close().extrude(T_DECKEL + 2).translate((0, 0, zd0 - 1)).val())
        bb_ = cq.Workplane("XY").polyline(pts).close().extrude(1).val().BoundingBox()
        LK_BB = (bb_.xmin, bb_.xmax, bb_.ymin, bb_.ymax)
        print("Luftklappenoeffnung im Deckel: X", round(bb_.xmin, 1), "-", round(bb_.xmax, 1), "Y", round(bb_.ymin, 1), "-", round(bb_.ymax, 1))
print("Deckel aus einem Stueck:", len(cq.Workplane().add(deckel).solids().vals()) == 1)
teile.append(("Zwischen_Deckel", deckel))
teile.append(("Zwischen_Laengswand", box(gx0 - T_NUT, gx1 + T_NUT, gy0, gy0 + T_UNTEN, Z0, zd0)))
# Boden unter dem abgedeckelten Bereich ausnehmen
boden = boden.cut(box(gx0, gx1, gy0 + T_UNTEN, gy1, Z_BODEN - T_BODEN - 1, Z_BODEN + 1))
boden = boden.cut(box(SX0, SX1, gy1, FY1, Z_BODEN - T_BODEN - 1, Z_BODEN + 1))   # bis an die Front (Leiste endet 5 mm unter Oberkante Front)
# Stirnseite rechts: der Streifen zwischen Kammerwand und Frontleiste ist ueber der Rahmen-Oberkante offen -> Fuellstueck
_sl2 = [s_ for n_, s_ in teile if n_ == "Schmal2_Leiste_rechts"][0].BoundingBox()
fs_x0 = _sl2.xmin
teile.append(("Fuellleiste_Streifen_rechts", box(fs_x0, L - D, gy1, FY0, Z0, zd0)))
boden = boden.cut(box(fs_x0, L - D, gy1, FY0, Z0 - 0.01, Z0 + 0.01)) if False else boden
print("Fuellstueck Stirnseite rechts: X", round(fs_x0, 1), "-", L - D, "Y", round(gy1, 1), "-", round(FY0, 1), "Z", Z0, "-", round(zd0, 1))
print("Eingeschlossen: Zwischenraum X", round(gx0, 1), "-", round(gx1, 1), "und Streifen X", round(SX0, 1), "-", round(SX1, 1), "Y", round(gy1, 1), "-", round(FY0, 1))

# ---------- Begleiter: Zwischenteil + Stimmstock ----------
# Boden: je Gruppe eine grosse runde Oeffnung (Ø17). Zwischenteil (Markus.Begleiterzwischenteil.dxf) 304 x 34 mm
# mit 8 Fenstern 30 x 22 (je vier Stimmplatten zusammengefasst), Hoehe 31,2 mm laut Schnitt.
# Stimmstock (Begleiter.dxf): Sole 3,6 mm mit 32 Oeffnungen 14 x 9,5 (eine je Stimmplatte),
# zwei Aussenwaende 3,6 mm und Mittelsteg 2,1 mm, 42,5 mm hoch (Schnitt), Deckleiste oben, Abschluesse an den Enden.
H_ZW, T_SOLE, H_WAND, T_AW, T_MS = 31.2, 3.6, 42.5, 3.6, 2.1
bz = ezdxf.readfile("bzw.dxf").modelspace(); bg = ezdxf.readfile("begl.dxf").modelspace()
zx0, zx1, zy0, zy1 = -309.5, -5.5, 102.0, 136.0                # Zwischenteil/Sole in Bassboden-Koordinaten
(ZX0, ZY0), (ZX1, ZY1) = v(zx0, zy0), v(zx1, zy1)
ZY1 = min(ZY1, gy0)                                             # nicht in die Wand der schmalen Kammern
zb0, zb1 = Z_BODEN, Z_BODEN + H_ZW
ZYW = ri1 + T_MITTE                                            # Zwischenteil 1,5 mm breiter: liegt an der Reihenleiste an
zw = cq.Workplane("XY").box(ZX1 - ZX0, ZY1 - ZYW, H_ZW, centered=False).translate((ZX0, ZYW, zb0))
n_f = 0
for e in bz.query('LWPOLYLINE[layer=="LAYER_1"]'):
    pts = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.05)][:-1]
    w_ = max(a for a, _ in pts) - min(a for a, _ in pts); h_ = max(b for _, b in pts) - min(b for _, b in pts)
    if abs(w_ - 30) < 0.5 and abs(h_ - 22) < 0.5:               # Fenster (separat gezeichnet: um -5 / +135,5 verschoben)
        pp = [v(a - 5.0, b + 135.5) for a, b in pts]
        zw = zw.cut(cq.Workplane("XY").polyline(pp).close().extrude(H_ZW + 2).translate((0, 0, zb0 - 1))); n_f += 1
teile.append(("Begleiter_Zwischenteil", zw.val()))
# offener Bereich am Schraubenende: zwischen Zwischenteil-Ende und Seitenteil bis Kammeroberkante fuellen
teile.append(("Fuellleiste_Begleiter_rechts", box(ZX1, L - D, ZYW, ZY1, Z0, Z1)))
# Sole und Stimmstock nach Begleiter.cb (CamBam): Sole 288,4 x 34 x 6 mit 32 Oeffnungen 14 x 9,5;
# Stimmstock aus einem Holzblock 288,4 x 25 x 50, je Seite 16 gleiche Kammern 16 x 45 mm, 9,5 mm tief (Taschen), unten offen zur Sole.
T_SOLE, L_ST, B_ST, H_ST = 6.0, 288.4, 20.0, 50.0       # Block nur 20 mm: Kammern 9,5 tief, Mittelwand 1 mm
L_BLK = 278.4                                              # Block: letzte Trennwand 5 mm; Sole 10 mm laenger (Schraube)
K_B, K_H, K_T, K_PITCH, K_X0 = 16.0, 45.0, 9.5, 17.0, 2.4
s0, s1 = zb1, zb1 + T_SOLE
SX0b = v(-309.5, 0)[0]; SX1b = SX0b + L_ST
sole = cq.Workplane("XY").box(L_ST, ZY1 - ZY0, T_SOLE, centered=False).translate((SX0b, ZY0, s0))
n_s = 0
# Sole-Oeffnungen so gross wie moeglich = Grundflaeche jeder Kammer (16 x 9,5); werden nach dem Verleimen
# zusammen mit den Kammern gefraest (3 Arbeitsgaenge: Kammern vorne, Kammern hinten, Sole-Oeffnungen)
_by0 = v(0, 102.0 + 17.0)[1] - B_ST/2
for i in range(16):
    kx = SX0b + K_X0 + i * K_PITCH
    for (ya, yb) in ((_by0, _by0 + K_T), (_by0 + B_ST - K_T, _by0 + B_ST)):
        sole = sole.cut(cq.Workplane("XY").box(K_B, yb - ya, T_SOLE + 2, centered=False).translate((kx, ya, s0 - 1))); n_s += 1
teile.append(("Begleiter_Sole", sole.val()))
by0 = v(0, 102.0 + 17.0)[1] - B_ST/2                              # Block mittig ueber der Trennung der zwei Oeffnungsreihen
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
# linke Nase + Schlitz in der Einzelkammer-Wand + Aufdoppelung
SCHL_T, SCHL_B = 1.5, 8.0                     # minimale Nase; innen aufdoppeln geht nicht (Folie/Stimmplatte)
wand_x = ex[1] + T_WAND/2                     # Aussenseite der Wand (zum Stimmstock)
nase_l = box(wand_x - SCHL_T, SX0b, ymid - SCHL_B/2, ymid + SCHL_B/2, s1, Z1)
blk = blk.union(cq.Workplane().add(nase_l))
i_w = [n for n, _ in teile].index("Einzel_Trennwand_2")
schlitz = box(wand_x - SCHL_T - 0.1, wand_x + 1, ymid - SCHL_B/2 - 0.1, ymid + SCHL_B/2 + 0.1, s1, Z1 + 1)
teile[i_w] = ("Einzel_Trennwand_2", teile[i_w][1].cut(schlitz))
assert len(blk.solids().vals()) == 1, "Stimmstock nicht aus einem Stueck"
teile.append(("Begleiter_Stimmstock", blk.val()))
teile.append(("Begleiter_Schraube", schr))
w1 = s1 + H_ST
print("Begleiter: Fenster im Zwischenteil", n_f, "Sole-Oeffnungen", n_s, "| Oberkante Stimmstock Z", round(w1, 1))

# ---------- Helikon-Stimmplatten ----------
# Keilfoermig: 5 mm am Schraubenende, 16 mm am anderen Ende (Ende, an dem die Schlitze buendig enden);
# umlaufender Rand 2 mm dick, Plattenkoerper darueber 2 mm schmaeler (je Seite). Liegt auf der Kammer (Z1).
T_DUENN, T_DICK, T_RAND, RUECK = 5.0, 16.0, 2.0, 2.0
schlitze, schrauben = [], []
for e in ents:
    if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == "HB-Rahmen" and e.closed:
        xs = [p[0] for p in e.get_points("xy")]; ys = [p[1] for p in e.get_points("xy")]
        w, h = max(xs) - min(xs), max(ys) - min(ys)
        if sorted((round(w, 1), round(h, 1))) in ([9.2, 71.5], [7.5, 51.0]):
            schlitze.append((min(xs), max(xs), min(ys), max(ys)))
    elif e.dxftype() == "CIRCLE" and e.dxf.layer == "HB-Schrauben" and e.dxf.radius < 1:
        schrauben.append((e.dxf.center.x, e.dxf.center.y, e.dxf.radius))
schlitze = sorted(set(tuple(round(c, 2) for c in s) for s in schlitze))
kammern_rect = [(-364.5, -311.5, 4.5, 98.5), (-311.5, -258.5, 4.5, 98.5), (-258.5, -205.5, 4.5, 98.5), (-205.5, -159.5, 4.5, 98.5),
                (-159.5, -106.5, 4.5, 98.5), (-106.5, -60.5, 4.5, 98.5), (-60.5, -7.5, 4.5, 98.5),
                (-364.5, -311.5, 98.5, 192.5), (-311.4, -217.4, 138.0, 184.0), (-98.5, -4.5, 138.0, 184.0)]
n_pl = 0
platten_info = []
for (kx0, kx1, ky0, ky1) in kammern_rect:
    ss = [s for s in schlitze if kx0 < (s[0]+s[1])/2 < kx1 and ky0 < (s[2]+s[3])/2 < ky1]
    if not ss: continue
    x0, x1 = min(s[0] for s in ss), max(s[1] for s in ss); y0, y1 = min(s[2] for s in ss), max(s[3] for s in ss)
    laengs_y = (y1 - y0) > (x1 - x0)
    # Aussenmass der Platte = Aussenmass der Kammer (53 bzw. 46 x 94 mm), 2-mm-Rand liegt auf der Kammer auf
    px0, px1, py0, py1 = kx0, kx1, ky0, ky1
    if laengs_y:   # buendiges Schlitzende = dickes Ende
        dick_hoch = sum(abs(s[3] - y1) < 0.2 for s in ss) >= sum(abs(s[2] - y0) < 0.2 for s in ss)
    else:
        dick_hoch = sum(abs(s[1] - x1) < 0.2 for s in ss) >= sum(abs(s[0] - x0) < 0.2 for s in ss)
    (PX0, PY0), (PX1, PY1) = v(px0, py0), v(px1, py1)
    # Rand
    pl = cq.Workplane("XY").box(PX1 - PX0, PY1 - PY0, T_RAND, centered=False).translate((PX0, PY0, Z1))   # Rand liegt auf der Kammer
    # Keilkoerper: Profil entlang der Laengsachse, ueber die Breite gezogen
    if laengs_y:
        a0, a1, b0, b1 = PY0 + RUECK, PY1 - RUECK, PX0 + RUECK, PX1 - RUECK
        t0, t1 = (T_DUENN, T_DICK) if dick_hoch else (T_DICK, T_DUENN)
        prof = [(a0, Z1), (a1, Z1), (a1, Z1 + T_RAND - t1), (a0, Z1 + T_RAND - t0)]   # Keil nach unten in die Kammer
        keil = cq.Workplane("YZ").polyline(prof).close().extrude(b1 - b0).translate((b0, 0, 0))
    else:
        a0, a1, b0, b1 = PX0 + RUECK, PX1 - RUECK, PY0 + RUECK, PY1 - RUECK
        t0, t1 = (T_DUENN, T_DICK) if dick_hoch else (T_DICK, T_DUENN)
        prof = [(a0, Z1), (a1, Z1), (a1, Z1 + T_RAND - t1), (a0, Z1 + T_RAND - t0)]   # Keil nach unten in die Kammer
        keil = cq.Workplane("XZ").polyline(prof).close().extrude(-(b1 - b0)).translate((0, b0, 0))
    pl = pl.union(keil)
    for s in ss:   # Zungenschlitze durch
        (sx0, sy0), (sx1, sy1) = v(s[0], s[2]), v(s[1], s[3])
        pl = pl.cut(cq.Workplane("XY").box(sx1 - sx0, sy1 - sy0, T_DICK + 4, centered=False).translate((sx0, sy0, Z1 - T_DICK)))
    for (cx, cy, r) in schrauben:
        if px0 < cx < px1 and py0 < cy < py1:
            pl = pl.cut(cq.Workplane("XY").center(*v(cx, cy)).circle(r).extrude(T_DICK + 4).translate((0, 0, Z1 - T_DICK)))
    n_pl += 1
    teile.append((f"Stimmplatte_{n_pl}", pl.val()))
    platten_info.append((n_pl, PX0, PX1, PY0, PY1, laengs_y))
print("Stimmplatten", n_pl, "Schlitze gesamt", len(schlitze))

# ---------- Befestigung der Stimmplatten ----------
# Je Platte oben und unten ein Schraubhaken (Draht 1,8 mm, Winkel 5 / 15 mm, Gewindeende ins Holz),
# der kurze Schenkel liegt auf dem Rand der Platte. Wo zwei Platten aneinanderstossen: eine Schraube fuer beide.
D_HAKEN, L_KURZ, L_LANG = 1.8, 5.0, 15.0
D_KOPF, H_KOPF, D_SCHAFT = 5.0, 1.5, 2.5       # Schraube an den Stossstellen - ANNAHME
ZH = Z1 + T_RAND + D_HAKEN/2                   # Hoehe der Drahtachse ueber dem Rand
def zyl(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1)
    return cq.Solid.makeCylinder(d/2, (p1 - p0).Length, p0, (p1 - p0).normalized())
def haken(px, py, dx, dy):
    """Senkrechter Schenkel bei (px,py), kurzer Schenkel 5 mm in Richtung (dx,dy) ueber den Plattenrand"""
    senk = zyl((px, py, ZH - L_LANG + D_HAKEN/2), (px, py, ZH), D_HAKEN)
    quer = zyl((px, py, ZH), (px + dx*L_KURZ, py + dy*L_KURZ, ZH), D_HAKEN)
    ecke = cq.Solid.makeSphere(D_HAKEN/2, cq.Vector(px, py, ZH))
    return senk.fuse(quer).fuse(ecke)
def schraube(px, py):
    kopf = zyl((px, py, Z1 + T_RAND), (px, py, Z1 + T_RAND + H_KOPF), D_KOPF)
    schaft = zyl((px, py, Z1 + T_RAND - L_LANG), (px, py, Z1 + T_RAND), D_SCHAFT)
    return kopf.fuse(schaft)
# Enden jeder Platte (Laengsachse), Mittelpunkt in der Breite
enden = []   # (platte, x, y, Richtung nach innen dx, dy)
for (nr, x0, x1, y0, y1, ly) in platten_info:
    if ly:
        enden += [(nr, (x0+x1)/2, y0, 0, 1), (nr, (x0+x1)/2, y1, 0, -1)]
    else:
        enden += [(nr, x0, (y0+y1)/2, 1, 0), (nr, x1, (y0+y1)/2, -1, 0)]
befest, geteilt = [], set()
for i, (nr, ex_, ey_, dx, dy) in enumerate(enden):
    # Nachbarplatte direkt vor diesem Ende (Spalt <= 6 mm)?
    for (nr2, a0, a1, b0, b1, _) in platten_info:
        if nr2 == nr: continue
        if dx == 0:   # Ende quer zu Y
            kante = b1 if dy == 1 else b0
            if 0 <= (ey_ - kante) * dy <= 6 and a0 - 1 < ex_ < a1 + 1:
                key = tuple(sorted((nr, nr2)))
                if key not in geteilt:
                    geteilt.add(key); befest.append((f"Schraube_P{key[0]}_P{key[1]}", schraube(ex_, (ey_ + kante)/2)))
                break
            # Seitenkante der Nachbarplatte (quer liegende Platte)
        else:
            kante = a1 if dx == 1 else a0
            if 0 <= (ex_ - kante) * dx <= 6 and b0 - 1 < ey_ < b1 + 1:
                key = tuple(sorted((nr, nr2)))
                if key not in geteilt:
                    geteilt.add(key); befest.append((f"Schraube_P{key[0]}_P{key[1]}", schraube((ex_ + kante)/2, ey_)))
                break
    else:
        off = 0.1 + D_HAKEN/2               # Draht knapp neben der Plattenkante
        befest.append((f"Haken_P{nr}_{'a' if i % 2 == 0 else 'b'}", haken(ex_ - dx*off, ey_ - dy*off, dx, dy)))
for n_b, s_b in befest:
    if n_b.startswith("Schraube"):
        bbk = s_b.BoundingBox(); cx, cy = (bbk.xmin + bbk.xmax)/2, (bbk.ymin + bbk.ymax)/2
        loch_b = cq.Workplane("XY").center(cx, cy).circle(D_SCHAFT/2 + 0.25).extrude(T_DICK + 10).translate((0, 0, Z1 - T_DICK - 2)).val()
        for k_ in range(len(teile)):
            if teile[k_][0].startswith("Stimmplatte"):
                teile[k_] = (teile[k_][0], teile[k_][1].cut(loch_b))
print("Haken", sum(n.startswith("Haken") for n, _ in befest), "Schrauben", sum(n.startswith("Schraube") for n, _ in befest), sorted(geteilt))
# Kontrolle: Befestigung darf keine Stimmplatte schneiden
for n, s in befest:
    for m_, p_ in teile:
        if m_.startswith("Stimmplatte") and s.intersect(p_).Volume() > 0.05:
            print("Befestigung schneidet", n, m_, round(s.intersect(p_).Volume(), 2))

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

exec(open("mechanik.py", encoding="utf-8").read())
a = cq.Assembly(name="Bassteil_Markus")
# Links/rechts getauscht (auf Johanns Wunsch): alles an der Ebene X = L/2 spiegeln
SPIEGELN = True
_add_orig = a.add
def _add_gespiegelt(obj, *args, **kw):
    if SPIEGELN and isinstance(obj, cq.Shape):
        obj = obj.mirror("YZ", basePointVector=(L / 2, 0, 0))
    return _add_orig(obj, *args, **kw)
a.add = _add_gespiegelt
farben = {"Rahmen": (0.78, 0.60, 0.38), "Boden": (0.88, 0.80, 0.62), "Leiste": (0.55, 0.40, 0.25),
          "Trennwand": (0.65, 0.75, 0.88), "Folie": (0.85, 0.30, 0.30), "Stimmplatte": (0.80, 0.82, 0.86), "Zwischen": (0.70, 0.55, 0.35), "Frontleiste": (0.55, 0.40, 0.25), "Begleiter": (0.62, 0.72, 0.55), "Fuellleiste": (0.55, 0.40, 0.25), "Front_Aufdoppelung": (0.78, 0.60, 0.38), "Deckelauflage": (0.60, 0.45, 0.28), "Eckblock": (0.55, 0.40, 0.25), "Bassdeckel": (0.30, 0.22, 0.15), "Bassbodenleder": (0.12, 0.08, 0.05), "Lufttaste_Filz": (0.85, 0.85, 0.80), "Lufttaste": (0.35, 0.25, 0.15), "Aufdoppelung": (0.65, 0.75, 0.88)}
for n, s in [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr)]:
    a.add(s, name=n, color=cq.Color(*farben["Rahmen"]))
a.add(boden, name="Bassboden", color=cq.Color(*farben["Boden"]))
for n, s in befest:
    a.add(s, name=n, color=cq.Color(0.25, 0.25, 0.28))
for n, s in teile:
    f = next(c for key, c in farben.items() if key in n)
    a.add(s, name=n, color=cq.Color(*f))
# Frontblende (Markus-20-mm-Blende.dxf): Aufsatz aussen auf der Front ueber der kurzen Bassdruecker-Reihe,
# 4 Bohrungen Ø10 deckungsgleich mit der Front, 7 mm hoch
fb = ezdxf.readfile("blende.dxf").modelspace()
T_BL, Z_REIHE, X_REF = 7.0, 38.0 + 4.24, 161.0      # deckungsgleich mit den Frontbohrungen der kurzen Reihe (Rahmen y = 38, X 161..221)
umr = [e for e in fb.query("LWPOLYLINE") if e.closed][0]
pts = [(X_REF - q.x, Z_REIHE + q.y) for q in ezdxf.path.make_path(umr).flattening(0.05)][:-1]
bl = cq.Workplane("XZ").polyline(pts).close().extrude(-T_BL).translate((0, B, 0))
for c in fb.query("CIRCLE"):
    bl = bl.cut(cq.Workplane("XZ").center(X_REF - c.dxf.center.x, Z_REIHE + c.dxf.center.y).circle(c.dxf.radius).extrude(-T_BL - 2).translate((0, B - 1, 0)))
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
