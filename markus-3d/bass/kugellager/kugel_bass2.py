# Bassteil Standardmechanik mit 4/4 Helikon (2-reihig, 9 + 8 Knoepfe; vormals "Kugellager-Bassteil v2") nach der Hauptzeichnung Bassteil8-gederht-ausgebssert.dxf
# Aufruf: python3 kugel_bass2.py <teilung 18|19|20> <laenge_oder_diskant> <ausgabe.step>
# Quellen: kb_haupt.dxf (Grundriss: Klappen, Oeffnungen, Trennwaende, Achslinien, Knoepfe; Schnitt: Rahmen, Boden,
# Basskammern, schraeges Begleiterbrett, Begleiter-Stimmstock), rahmen.dxf + horn.dxf (Markus: Zinken, Verzierungen,
# Hoerner), kb_front.dxf (Knopf-/Hornhoehen), kb_deckel.dxf (Bassdeckel, Lufttastenschlitz).
# Achsen: X = Laenge, Y = Tiefe (0 = Rueckteil aussen, B = Front aussen), Z = Hoehe ab Deckelseite (Z 0) Richtung Balg.
# Grundriss: X = x + D, Y = y - 272.1 + D (D = 5).  Schnitt: Z = x - 401.84, Y = y - 267.06 (+0.46 auf den Grundriss).
# Mechanik: 6 Achslinien parallel zur Front (Grundriss y 379.1 ... 394.1), Zuordnung Knopf -> Klappe nach dem Grundriss.
# Bassachsen: Knopfarm nach oben (Schubstange Z 20.5), Hebel nach hinten unter die Bassklappe (oeffnet nach unten).
# Begleiterachsen: Knopfarm nach unten (Schubstange Z 7) = Drehung andersherum; Hebel steigt nach oben unter die
# Begleiterklappe am schraegen Brett und oeffnet sie nach unten ("von der anderen Seite bedient").
import sys, os, math, json, subprocess, ezdxf, ezdxf.path, cadquery as cq
H = os.path.dirname(os.path.abspath(__file__))
T = float(sys.argv[1]); LD = float(sys.argv[2]); OUT = sys.argv[3]
L_END = max(380.0, math.ceil(LD / 10.0 - 1e-9) * 10.0)
D, L, B = 5.0, 380.0, 205.0                                   # Gehaeuse aussen 380 x 205, Wandstaerke 5 (Fertigmass)
H_RAHMEN, H_ZW = 66.0, 28.0                                   # Rahmen 66, Zwischenrahmen 66..94 (Schnitt)
T_DECKEL = 4.5
VAR45 = os.environ.get("KB_VAR", "") == "4-5"                 # Variante 4/5: zusaetzliche halbe Kammer vorne quer
NAME_V = "4-5" if VAR45 else "4-4"
Z_UB, Z_OB = 25.5, 29.0                                       # Bassboden 3.5 (Unterseite wie Schnitt, Klappen/Mechanik darunter)
Z_KT = Z_OB + 47.31                                           # Basskammer-Oberkante (Trennwandhoehe 47.31)
T_LEDER, T_FILZ, T_KDECK = 0.5, 4.0, 0.5                     # Klappenaufbau: Leder, Filz (rot), Deckel
T_KL = T_FILZ + T_KDECK
Z_KL_U = Z_UB - T_KL - T_LEDER                                # Bassklappen-Unterkante 20.5
Z_A = 15.0                                                    # Achsmitte
Z_RB, Z_RG = 20.5, 7.0                                        # Schubstangen Bass / Begleiter
Y_BR0, Z_BR0 = 385.1 - 267.1, 64.84                           # schraeges Brett: Innenflaeche beginnt hier (= Vorderleiste)
ALPHA = math.degrees(math.atan2(10.0, 82.0)); TA = math.tan(math.radians(ALPHA)); CA = math.cos(math.radians(ALPHA))
T_BRETT = 3.5                                                 # Begleiterboden 3.5
def X(px): return px + D                                     # Grundriss: Innenkante bei x 0 / y 272.1
def Y(py): return py - 272.1 + D
def XK(px): return L / 2 + (px - 185.5) * T / 19.0            # Knopf-X bei Teilung T (Zeichnung 19, Mitte 185.5)
def box(x0, x1, y0, y1, z0, z1):
    return cq.Solid.makeBox(x1 - x0, y1 - y0, z1 - z0, cq.Vector(x0, y0, z0))
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v.Length, p0, v.normalized())
def kugel(p, d): return cq.Solid.makeSphere(d / 2 + 0.01, cq.Vector(*p), angleDegrees1=-90, angleDegrees2=90)
def pfad(pts, d):
    s = None
    for a, b in zip(pts[:-1], pts[1:]):
        seg = stab(a, b, d).fuse(kugel(b, d)); s = seg if s is None else s.fuse(seg)
    return s.fuse(kugel(pts[0], d))
def pts_of(e, tol=0.1):
    p = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(tol)]
    q = [p[0]] + [b_ for a_, b_ in zip(p[:-1], p[1:]) if math.dist(a_, b_) > 1e-6]
    if len(q) > 2 and math.dist(q[0], q[-1]) < 1e-6: q = q[:-1]
    return q
def prisma_xy(pts, z0, z1): return cq.Workplane("XY").polyline(pts).close().extrude(z1 - z0).translate((0, 0, z0)).val()
def platte_yz(pts, x0, x1):                                  # Polygon in (Y, Z), in X extrudiert
    return cq.Workplane("YZ").polyline(pts).close().extrude(x1 - x0).translate((x0, 0, 0)).val()
def mod(n): return os.path.join(H, n)
def schraeg(s): return s.rotate(cq.Vector(0, Y_BR0, Z_BR0), cq.Vector(1, Y_BR0, Z_BR0), ALPHA)
def lbox(x0, x1, y0, y1, z0, z1):                             # Box im Brett-System (y laengs der Schraege ab Y_BR0, z normal)
    return schraeg(box(x0, x1, Y_BR0 + y0, Y_BR0 + y1, Z_BR0 + z0, Z_BR0 + z1))
def yl(Yp): return (Yp - Y_BR0) / CA                          # Grundriss-Y -> Laenge auf dem Brett
def z_innen(Yp): return Z_BR0 + (Yp - Y_BR0) * TA             # Brett-Innenflaeche ueber Grundriss-Y
def z_kl_u(Yp): return z_innen(Yp) - (T_KL + T_LEDER) / CA    # Unterkante Begleiter-/Luftklappe
teile, mech = [], []
INNEN = box(D, L - D, D, B - D, -20, 300)
HM = ezdxf.readfile(mod("kb_haupt.dxf")).modelspace()
def rechtecke(layer, xr=(-5, 390), yr=(250, 470)):
    out = set()
    for e in HM.query(f'LWPOLYLINE[layer=="{layer}"]'):
        if not e.closed: continue
        p = pts_of(e); xs = [a for a, _ in p]; ys = [b for _, b in p]
        if xr[0] < min(xs) and max(xs) < xr[1] and yr[0] < min(ys) and max(ys) < yr[1]:
            out.add((round(min(xs), 2), round(max(xs), 2), round(min(ys), 2), round(max(ys), 2)))
    return sorted(out)

# ---------------- Knoepfe (Hoehen aus kb_front.dxf) ----------------
FR = ezdxf.readfile(mod("kb_front.dxf")).modelspace()
Z_R = {}; hoerner_xz = []
for c in FR.query('CIRCLE[layer=="Front"]'):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if r in (5.0, 5.5): Z_R[1 if r == 5.5 else 2] = (-y, r)
    elif r == 9.0: hoerner_xz.append((L / 2 + x, -y))
KN_PX = {1: [109.5 + 19 * i for i in range(9)], 2: [119.0 + 19 * i for i in range(8)]}
kn = [(XK(px), Z_R[rr][0], Z_R[rr][1], rr, px) for rr in (1, 2) for px in KN_PX[rr]]
def knopf(rr, px): return [k for k in kn if k[3] == rr and abs(k[4] - px) < 0.1][0]

# ---------------- Rahmen (Markus: Zinken, Verzierungen) + Zwischenrahmen ----------------
bohr = [[-(L - X_), Z_ - 4.24, r_] for X_, Z_, r_, rr, px in kn] + [[-(L - X_), Z_ - 4.24, 9.0] for X_, Z_ in hoerner_xz]   # gespiegelt (s. Zusammenbau)
json.dump(bohr, open(mod("kb_bohrungen.json"), "w"))
rt = {}
for t_ in ("front", "rueck", "seite"):
    fn = mod(f"_kb_{t_}.step")
    r_ = subprocess.run([sys.executable, mod("kb_rahmen.py"), mod("rahmen.dxf"), fn, t_], cwd=H, capture_output=True, text=True,
                        env=dict(os.environ, KB_L=str(L), KB_B=str(B), KB_D=str(D)))
    print(r_.stdout.strip());
    if r_.returncode: raise SystemExit(r_.stderr)
    rt[t_] = cq.importers.importStep(fn).val()
O, XA, ZA = cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), cq.Vector(0, 0, 1)
st = lambda s: s.rotate(O, XA, 90)
rueck = st(rt["rueck"]).translate(cq.Vector(0, D, 0))
front = st(rt["front"]).rotate(O, ZA, 180).translate(cq.Vector(L, B - D, 0))
sl = st(rt["seite"]).rotate(O, ZA, -90).translate(cq.Vector(D, B, 0))
sr = st(rt["seite"]).rotate(O, ZA, 90).translate(cq.Vector(L - D, 0, 0))
teile += [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr)]
RAHMEN_ = ("Front", "Rueckteil", "Seitenteil_links", "Seitenteil_rechts")
z0, z1 = H_RAHMEN, H_RAHMEN + H_ZW
# Balgrahmen (66..94) gezinkt wie der Rahmen: Zinken 6 hoch; an den Ecken abwechselnd Laengs- und Seitenteil.
_zA = [(0, 6), (12, 18), (24, H_ZW)]; _zB = [(6, 12), (18, 24)]
def _ecken_weg(s_, xr_, yr_, zl_):
    for (xa_, xb_) in xr_:
        for (ya_, yb_) in yr_:
            for (za_, zb_) in zl_: s_ = s_.cut(box(xa_, xb_, ya_, yb_, z0 + za_, z0 + zb_))
    return s_
_ex = [(-1, D), (L - D, L + 1)]; _ey = [(-1, D), (B - D, B + 1)]
bg_v = _ecken_weg(box(0, L, B - D, B, z0, z1), _ex, [(B - D - 1, B + 1)], _zA)
bg_h = _ecken_weg(box(0, L, 0, D, z0, z1), _ex, [(-1, D + 1)], _zA)
bg_l = _ecken_weg(box(0, D, 0, B, z0, z1), [(-1, D + 1)], _ey, _zB)
bg_r = _ecken_weg(box(L - D, L, 0, B, z0, z1), [(L - D - 1, L + 1)], _ey, _zB)
# Oben (zum Balg) vier Einsatzleisten 10 breit x 6 hoch, schmale Seite (6) am Rahmen; jede um ihre Breite kuerzer
# (Windmuehlen-Anordnung): jede Leiste stoesst an genau einer Ecke an den Rahmen.
B_EL, H_EL = 10.0, 6.0
_xi0, _xi1, _yi0, _yi1 = D, L - D, D, B - D
teile += [("Balgrahmen_vorne", bg_v), ("Balgrahmen_hinten", bg_h), ("Balgrahmen_links", bg_l), ("Balgrahmen_rechts", bg_r),
          ("Balgrahmen_Einsatzleiste_hinten", box(_xi0, _xi1 - B_EL, _yi0, _yi0 + B_EL, z1 - H_EL, z1)),
          ("Balgrahmen_Einsatzleiste_rechts", box(_xi1 - B_EL, _xi1, _yi0, _yi1 - B_EL, z1 - H_EL, z1)),
          ("Balgrahmen_Einsatzleiste_vorne", box(_xi0 + B_EL, _xi1, _yi1 - B_EL, _yi1, z1 - H_EL, z1)),
          ("Balgrahmen_Einsatzleiste_links", box(_xi0, _xi0 + B_EL, _yi0 + B_EL, _yi1, z1 - H_EL, z1))]
# Hoerner (Markus horn.dxf, Drehprofil, Stufe an der Frontaussenseite)
HD = ezdxf.readfile(mod("horn.dxf")).modelspace()
segs = []
for e in HD.query('LINE ARC'):
    if e.dxftype() == "LINE" and abs(e.dxf.start.y) < 1e-6 and abs(e.dxf.end.y) < 1e-6: continue
    segs.append([(round(q.x, 4), round(q.y, 4)) for q in ezdxf.path.make_path(e).flattening(0.08)])
kette = segs.pop(0)
while segs:
    for i, s in enumerate(segs):
        if math.dist(s[0], kette[-1]) < 0.05: kette += s[1:]; segs.pop(i); break
        if math.dist(s[-1], kette[-1]) < 0.05: kette += s[::-1][1:]; segs.pop(i); break
    else: raise SystemExit("Hornprofil nicht geschlossen")
horn_prof = kette[:-1] if math.dist(kette[0], kette[-1]) < 0.05 else kette
for i, (X_, Z_) in enumerate(hoerner_xz):
    w = cq.Workplane("XY").polyline([(r_, a_) for a_, r_ in horn_prof]).close().revolve(360, (0, 0, 0), (0, 1, 0))
    teile.append((f"Horn_{i+1}", w.val().translate(cq.Vector(X_, B - 4.01, Z_))))
print("Rahmen Markus, Hoerner", len(hoerner_xz))

# ---------------- Bassboden mit Oeffnungen, Bassklappen ----------------
off = rechtecke("2Füllung-Öffnug")
bo_oeff = [o for o in off if 340 < o[2] < 380 and o[1] - o[0] < 60]
be_oeff = [o for o in off if 390 < o[2] < 420 and o[1] - o[0] < 60]
lu_oeff = [o for o in off if o[2] > 425 and o[1] - o[0] > 60][0]
B_SL = 6.0                                                    # seitliche Leisten (Seitenteile), Kammereinbau beidseits um B_SL kuerzer
boden = box(D + B_SL, L - D - B_SL, D, Y(385.1), Z_UB, Z_OB)
for o in bo_oeff: boden = boden.cut(box(X(o[0]), X(o[1]), Y(o[2]), Y(o[3]), Z_UB - 1, Z_OB + 1))
teile.append(("Bassboden", boden))
kl = rechtecke("HB-Klappe")
bass_kl = sorted([(X(k[0]), X(k[1]), Y(k[2]), Y(k[3])) for k in kl if k[2] < 390], key=lambda k: k[0])
begl_kl = sorted([(X(k[0]), X(k[1]), Y(k[2]), Y(k[3])) for k in kl if k[2] >= 390], key=lambda k: k[0])
for i, (x0, x1, y0, y1) in enumerate(bass_kl):
    teile.append((f"Bassklappe_{i+1}_Leder", box(x0, x1, y0, y1, Z_UB - T_LEDER, Z_UB)))
    teile.append((f"Bassklappe_{i+1}_Filz", box(x0 + 0.5, x1 - 0.5, y0 + 0.5, y1 - 0.5, Z_UB - T_LEDER - T_FILZ, Z_UB - T_LEDER)))
    teile.append((f"Bassklappe_{i+1}", box(x0 + 0.5, x1 - 0.5, y0 + 0.5, y1 - 0.5, Z_KL_U, Z_UB - T_LEDER - T_FILZ)))
print("Bassboden: Oeffnungen", len(bo_oeff), "| Klappen Bass", len(bass_kl), "Begleiter", len(begl_kl))

# ---------------- Basskammern (Detailzeichnungen HB-trennwaende.dxf / tremnnwand2.dxf) ----------------
# 9 Trennwaende 98 lang x 47.3 hoch x 3.5, durch die Leisten (Nuten) gesteckt; Vorderleiste 6 = Stirnflaeche der
# Mechanik (Y 118.5), Hinterleiste 12 (hinten 6 mm ueber das Plattenende hinaus fuer die Haken).
L_TW, T_TW = 98.0, 3.5
Y_LV1 = Y(385.1); Y_LV0 = max(Y(o[3]) for o in bo_oeff) + 1.0     # Vorderleiste so breit wie moeglich: bis 1 mm vor die Klappenoeffnungen
Y_TW1 = Y_LV1 - 3.0; Y_TW0 = Y_TW1 - L_TW                          # Trennwaende 17.5 .. 115.5, 3 mm Nut in der Vorderleiste
Y_LH0, Y_LH1 = D, Y_TW0 + 3.0                                      # Hinterleiste bis an den Rahmen (kein Zwischenraum)
TW_X = [x_ + D for x_ in (7.75, 61.75, 95.75, 149.75, 183.75, 237.75, 271.75, 305.75, 359.75)]
lh = box(D + B_SL, L - D - B_SL, Y_LH0, Y_LH1, Z_OB, Z_KT); lv = box(D + B_SL, L - D - B_SL, Y_LV0, Y_LV1, Z_OB, Z_KT)
teile += [("Basskammer_Leiste_links", box(D, D + B_SL, Y_LH0, Y_LV1, Z_UB, Z_KT)), ("Basskammer_Leiste_rechts", box(L - D - B_SL, L - D, Y_LH0, Y_LV1, Z_UB, Z_KT))]
for i, x_ in enumerate(TW_X):
    tw_ = box(x_, x_ + T_TW, Y_TW0, Y_TW1, Z_OB, Z_KT)
    teile.append((f"Basskammer_Trennwand_{i+1}", tw_)); lh = lh.cut(tw_); lv = lv.cut(tw_)
teile += [("Basskammer_Leiste_hinten", lh), ("Basskammer_Leiste_vorne", lv)]
# Gekruemmte Teiler (HB-trennwaende.dxf: Kurve in der Trennwand-Seitenansicht) je Kammer, 1.5 mm stark, beidseits
# 1 mm tief in genutete Trennwaende eingesetzt. Zeichnung x -94.88 = vorderes Wandende (Y_TW1), y 0.63 = Wand-Unterkante.
TD = ezdxf.readfile(mod("kb_hb_trennwand_detail.dxf")).modelspace()
kurve = [p_ for p_ in [pts_of(e, 0.05) for e in TD.query('LWPOLYLINE') if not e.closed] if min(b for _, b in p_) < 10][0]
kurve = sorted(kurve)                                               # vom Fuss (unten) zum freien Ende (oben)
T_TEIL, T_NUT = 1.5, 1.0
k_yz = [(Y_TW1 - (a + 94.88), Z_OB + (b - 0.63)) for a, b in kurve]   # von der anderen Seite: Fuss vorne, steigt nach hinten
k_yz[0] = (k_yz[0][0], Z_OB + 0.01)                                 # Fuss steht auf dem Boden
prof_t = k_yz + [(y_, z_ + T_TEIL) for y_, z_ in reversed(k_yz)]
def teiler_xy(x0_, x1_):
    return cq.Workplane("YZ").polyline(prof_t).close().extrude(x1_ - x0_).translate((x0_, 0, 0)).val()
for i in range(len(TW_X) - 1):
    xa_, xb_ = TW_X[i] + T_TW - T_NUT, TW_X[i + 1] + T_NUT
    teile.append((f"Basskammer_Teiler_{i+1}", teiler_xy(xa_, xb_)))
nut = teiler_xy(-1, 1).BoundingBox()
for k_ in [k for k, (n_, _) in enumerate(teile) if n_.startswith("Basskammer_Trennwand_")]:
    s_ = teile[k_][1]
    for i in range(len(TW_X) - 1):
        s_ = s_.cut(teiler_xy(TW_X[i] + T_TW - T_NUT - 0.05, TW_X[i + 1] + T_NUT + 0.05))
    teile[k_] = (teile[k_][0], s_)
print("Teiler", len(TW_X) - 1, "Fuss Y", round(k_yz[0][0], 1), "freies Ende Y", round(k_yz[-1][0], 1), "Z bis", round(k_yz[-1][1], 1))
# Helikon-Stimmplatten wie Markus: Keil 5 -> 16 mm (dickes Ende vorne, buendiges Schlitzende), Rand 2 mm liegt auf
# Trennwaenden und Leisten; doppelt 53 x 98, einfach 33 x 98; Platte = Wandmitte bis Wandmitte minus 1 mm Spalt.
T_DUENN, T_DICK, T_RAND, RUECK = 5.0, 16.0, 2.0, 2.0
Z1 = Z_KT
pl_info = []
for i in range(len(TW_X) - 1):
    xm0, xm1 = TW_X[i] + T_TW / 2, TW_X[i + 1] + T_TW / 2
    PX0, PX1 = xm0 + 0.5, xm1 - 0.5
    PY0, PY1 = Y_TW0, Y_TW1
    breit = PX1 - PX0
    pl = box(PX0, PX1, PY0, PY1, Z1, Z1 + T_RAND)
    a0, a1, b0, b1 = Y_LH1 + 0.5, Y_LV0 - 0.5, PX0 + RUECK, PX1 - RUECK      # Keil zwischen den Leisten
    prof = [(a0, Z1), (a1, Z1), (a1, Z1 + T_RAND - T_DICK), (a0, Z1 + T_RAND - T_DUENN)]
    pl = pl.fuse(cq.Workplane("YZ").polyline(prof).close().extrude(b1 - b0).translate((b0, 0, 0)).val())
    # Zungenschlitze nach der Markus-Zeichnung (boden2.dxf, Platte 53 x 94): Lage ab Plattenkante x 5.1 / 15.9 / 27.9 / 40.4,
    # Breiten 7.5 / 9.2 / 9.2 / 7.5, Laengen 51 / 71.5 / 71.5 / 51, alle enden 5.5 vor der Kante am dicken Ende.
    # Halbe Platte (33): zwei lange Zungen 9.2 x 71.5, Abstand wie bei Markus (2.8), mittig.
    if breit > 45: muster = [(5.1, 7.5, 51.0), (15.9, 9.2, 71.5), (27.9, 9.2, 71.5), (40.4, 7.5, 51.0)]
    else: muster = [(breit / 2 - 1.4 - 9.2, 9.2, 71.5), (breit / 2 + 1.4, 9.2, 71.5)]
    for (o_, sb, sll) in muster:
        pl = pl.cut(box(PX0 + o_, PX0 + o_ + sb, PY1 - 5.5 - sll, PY1 - 5.5, Z1 - T_DICK - 1, Z1 + T_RAND + 1))
    teile.append((f"Stimmplatte_{i+1}", pl)); pl_info.append((i + 1, PX0, PX1, PY0, PY1, breit))
# Befestigung wie Markus, nur Schraubhaken (Draht 1.8, kurzer Schenkel 5 auf dem Plattenrand), hinten und vorne
# je einer in der Leiste; vorne langer Schenkel 8.5 (unter dem Haken liegt die Nut fuer das Begleiterbrett)
D_HAKEN, L_KURZ = 1.8, 5.0
ZH = Z1 + T_RAND + D_HAKEN / 2
def haken_pl(xc, py_, dy_, l_lang, leiste):
    hk = stab((xc, py_, ZH - l_lang), (xc, py_, ZH + 0.3), D_HAKEN).fuse(stab((xc, py_ - dy_ * 0.3, ZH), (xc, py_ + dy_ * L_KURZ, ZH), D_HAKEN)).clean()
    _i = [k for k, (n_, _) in enumerate(teile) if n_ == leiste][0]
    teile[_i] = (teile[_i][0], teile[_i][1].cut(stab((xc, py_, ZH - l_lang - 0.5), (xc, py_, Z1 + 0.1), D_HAKEN)))
    return hk
for (nr, x0_, x1_, y0_, y1_, br_) in pl_info:
    xc = (x0_ + x1_) / 2; off = 0.1 + D_HAKEN / 2
    teile.append((f"Stimmplatte_{nr}_Haken_hinten", haken_pl(xc, y0_ - off, 1, 15.0, "Basskammer_Leiste_hinten")))
    teile.append((f"Stimmplatte_{nr}_Haken_vorne", haken_pl(xc, y1_ + off, -1, 8.5, "Basskammer_Leiste_vorne")))
print("Basskammern: Trennwaende", len(TW_X), "Kammern", len(pl_info), "| Stimmplatten", [round(p[5], 1) for p in pl_info],
      "x", L_TW, "| Kammer innen Y", round(Y_LV0 - Y_LH1, 1), "Platte Y", Y_TW0, "-", Y_TW1)

# ---------------- schraeges Begleiterbrett mit Begleiter- und Luftklappen ----------------
brett = lbox(D, L - D, -3.0, yl(B - D) + 5, 0, T_BRETT).intersect(INNEN)
for o in be_oeff + [lu_oeff]: brett = brett.cut(box(X(o[0]), X(o[1]), Y(o[2]), Y(o[3]), 50, 100))
brett = brett.cut(box(0, L, 0, Y(385.1) - 2.0, 0, 200))          # 2 mm Nut in der Kammer-Laengswand (Vorderleiste)
_iv2 = [k for k, (n_, _) in enumerate(teile) if n_ == "Basskammer_Leiste_vorne"][0]
teile[_iv2] = (teile[_iv2][0], teile[_iv2][1].cut(brett))
for _n in ("Basskammer_Leiste_links", "Basskammer_Leiste_rechts"):
    _k = [k for k, (n_, _) in enumerate(teile) if n_ == _n][0]
    teile[_k] = (_n, teile[_k][1].cut(brett))
teile.append(("Begleiterbrett", brett))
# Seitliche Leisten ueber die ganze Gehaeusebreite (Y): ueber den Basskammern bis zur Kammeroberkante, davor in der
# Seitenansicht schraeg bis unter das Begleiterbrett (stuetzt es). Front: Leiste unter dem Brett als Abschluss zum Rahmen,
# ueber der oberen Knopfreihe; ausserhalb des Knopfbereichs bis auf den Bassboden-Horizont herunter.
_ueber = lbox(-10, L + 10, 0, yl(B - D) + 20, 0, 60).intersect(box(-10, L + 10, Y_LV1, B + 10, 0, 200))   # alles ueber der Brettunterseite vor den Kammern
for _n in ("Basskammer_Leiste_links", "Basskammer_Leiste_rechts"):
    _k = [k for k, (n_, _) in enumerate(teile) if n_ == _n][0]
    _x0 = D if _n.endswith("links") else L - D - B_SL
    _sl = box(_x0, _x0 + B_SL, Y_LH0, B - D, Z_UB, Z_KT).cut(_ueber).cut(brett)
    teile[_k] = (_n, _sl)
_zfr = max(k[1] for k in kn) + 3.75 + 1.5                            # ueber den Knopfzylindern der oberen Reihe
_fx0, _fx1 = min(k[0] for k in kn) - 9.0, max(k[0] for k in kn) + 9.0
_fl = box(D + B_SL, L - D - B_SL, B - D - B_SL, B - D, _zfr, 120)
_fl = _fl.fuse(box(D + B_SL, _fx0, B - D - B_SL, B - D, Z_UB, 120)).fuse(box(_fx1, L - D - B_SL, B - D - B_SL, B - D, Z_UB, 120))
teile.append(("Begleiterbrett_Leiste_vorne", _fl.cut(_ueber).cut(brett).clean()))
print("Seitenleisten bis zur Front, Frontleiste ab Z", round(_zfr, 1), "(ausserhalb X", round(_fx0, 1), "-", round(_fx1, 1), "ab Z", Z_UB, ")")
def schraeg_klappe(nm, x0, x1, y0, y1):
    teile.append((nm + "_Leder", lbox(x0, x1, yl(y0), yl(y1), -T_LEDER, 0)))
    teile.append((nm + "_Filz", lbox(x0 + 0.5, x1 - 0.5, yl(y0) + 0.5, yl(y1) - 0.5, -T_LEDER - T_FILZ, -T_LEDER)))
    teile.append((nm, lbox(x0 + 0.5, x1 - 0.5, yl(y0) + 0.5, yl(y1) - 0.5, -T_LEDER - T_KL, -T_LEDER - T_FILZ)))
for i, k in enumerate(begl_kl): schraeg_klappe(f"Begleiterklappe_{i+1}", *k)
k2 = [k for k in rechtecke("2Klappen")][0]
luft_kl = (X(k2[0]), X(k2[1]), Y(k2[2]), Y(k2[3]))
schraeg_klappe("Luftklappe", *luft_kl)
# Begleiter-Sole + Stimmstock nach Begleiter.dxf (Kugellager): je Klappenoeffnung 4 Kammern (2 je Seite, 15 breit, 44 hoch,
# Kammertiefe 9.5 wie Markus), Block 20 breit / 50 hoch, Sole 6 mit einer Oeffnung je Kammer, Stimmplatte je Kammer.
BG = ezdxf.readfile(mod("kb_begleiter.dxf")).modelspace(); BGX = 349.0 - (5.5 - D)          # Zeichnung x -> X (Oeffnungen deckungsgleich)
yc = yl((Y(be_oeff[0][2]) + Y(be_oeff[0][3])) / 2)
T_SOLE, B_ST, H_ST, K_T, K_H = 3.0, 20.0, 50.0, 9.5, 44.0
T_SPL, B_SPL, H_SPL = 2.5, 16.8, 47.0
A_SL, B_SL_Z, L_SL_Z, Z_SL0, Z_NIET, D_NIET = 7.0, 2.6, 34.0, 6.0, 3.0, 1.6   # Schlitzabstand, -breite, -laenge, Beginn ueber Plattenunterkante; Nietloch
sd = pts_of([e for e in BG.query('LWPOLYLINE[layer=="DEFPOINTS"]')][0])
so0, so1 = min(a for a, _ in sd) + BGX, max(a for a, _ in sd) + BGX
B_SO = max(b for _, b in sd) - min(b for _, b in sd)
bo = [e for e in BG.query('LWPOLYLINE[layer=="0"]') if (lambda p_: max(a for a, _ in p_) - min(a for a, _ in p_) > 200)(pts_of(e))][0]
bx0, bx1 = min(a for a, _ in pts_of(bo)) + BGX, max(a for a, _ in pts_of(bo)) + BGX
L_VERL = 10.0                                                  # Schraubansatz je Blockende; Sole laeuft darunter durch
so0, so1 = min(so0, bx0 - L_VERL), max(so1, bx1 + L_VERL)
sole = lbox(so0, so1, yc - B_SO / 2, yc + B_SO / 2, T_BRETT, T_BRETT + T_SOLE)
blk = lbox(bx0, bx1, yc - B_ST / 2, yc + B_ST / 2, T_BRETT + T_SOLE, T_BRETT + T_SOLE + H_ST)
kz = sorted({(round(min(a for a, _ in p_), 2), round(max(a for a, _ in p_), 2)) for p_ in [pts_of(e) for e in BG.query('LWPOLYLINE[layer=="KANZELLE"]')]})
n_spl = 0
for (k0, k1) in kz:
    x0_, x1_ = k0 + BGX, k1 + BGX; xm_ = (x0_ + x1_) / 2
    for sg in (-1, 1):
        ya, yb = (yc - B_ST / 2, yc - B_ST / 2 + K_T) if sg < 0 else (yc + B_ST / 2 - K_T, yc + B_ST / 2)
        blk = blk.cut(lbox(x0_, x1_, ya - (1 if sg < 0 else 0), yb + (1 if sg > 0 else 0), T_BRETT + T_SOLE - 1, T_BRETT + T_SOLE + K_H))
        sole = sole.cut(lbox(x0_ + 0.5, x1_ - 0.5, ya, yb, T_BRETT - 1, T_BRETT + T_SOLE + 1))
        yp0, yp1 = (yc - B_ST / 2 - T_SPL, yc - B_ST / 2) if sg < 0 else (yc + B_ST / 2, yc + B_ST / 2 + T_SPL)
        n_spl += 1
        zb_ = T_BRETT + T_SOLE + 0.5
        spl_ = lbox(xm_ - B_SPL / 2, xm_ + B_SPL / 2, yp0, yp1, zb_, zb_ + H_SPL)
        for dx_ in (-A_SL / 2, A_SL / 2):                       # zwei Zungenschlitze nebeneinander, je ein Nietloch an der Wurzel
            spl_ = spl_.cut(lbox(xm_ + dx_ - B_SL_Z / 2, xm_ + dx_ + B_SL_Z / 2, yp0 - 1, yp1 + 1, zb_ + Z_SL0, zb_ + Z_SL0 + L_SL_Z))
            spl_ = spl_.cut(schraeg(stab((xm_ + dx_, Y_BR0 + yp0 - 1, Z_BR0 + zb_ + Z_NIET), (xm_ + dx_, Y_BR0 + yp1 + 1, Z_BR0 + zb_ + Z_NIET), D_NIET)))
        teile.append((f"Begleiter_Stimmplatte_{n_spl}", spl_))
print("Begleiterstimmstock: Kammern", 2 * len(kz), "Stimmplatten", n_spl, "Block X", round(bx0, 1), "-", round(bx1, 1))
# Seitliche Verlaengerung an beiden Blockenden nur fuer die Verschraubung (auf dem Brett, neben der Sole): 10 lang / 10 hoch;
# je eine Schraube Ø3 (Kopf Ø5.5) ins Brett (endet 0.5 vor der Unterseite).
L_ST1, H_ST1, L_ST2, H_ST2 = 0.0, 10.0, L_VERL, 10.0
sole_huelle = lbox(so0, so1, yc - B_SO / 2, yc + B_SO / 2, T_BRETT, T_BRETT + T_SOLE)
schr_st = []
for sg, xe in ((-1, bx0), (1, bx1)):
    xa1, xb1 = sorted((xe, xe + sg * L_ST1)); xa2, xb2 = sorted((xe + sg * L_ST1, xe + sg * (L_ST1 + L_ST2)))
    ext = lbox(xa2, xb2, yc - B_ST / 2, yc + B_ST / 2, T_BRETT + T_SOLE, T_BRETT + H_ST2)   # auf der Sole
    xs_ = (xa2 + xb2) / 2
    sc_ = schraeg(stab((xs_, Y_BR0 + yc, Z_BR0 + 0.5), (xs_, Y_BR0 + yc, Z_BR0 + T_BRETT + H_ST2), 3.0).fuse(
        stab((xs_, Y_BR0 + yc, Z_BR0 + T_BRETT + H_ST2), (xs_, Y_BR0 + yc, Z_BR0 + T_BRETT + H_ST2 + 2.0), 5.5)))
    blk = blk.fuse(ext.cut(sole_huelle).cut(sc_)).clean()
    schr_st.append(sc_)
_ib3 = [k for k, (n_, _) in enumerate(teile) if n_ == "Begleiterbrett"][0]
_b3 = teile[_ib3][1]
for sc_ in schr_st: _b3 = _b3.cut(sc_)
teile[_ib3] = ("Begleiterbrett", _b3)
sole = sole.cut(schr_st[0]).cut(schr_st[1])
teile += [("Begleiter_Sole", sole), ("Begleiterstimmstock", blk)]
teile += [(f"Begleiterstimmstock_{j_+1}_Schraube", sc_) for j_, sc_ in enumerate(schr_st)]
print("Stimmstock-Verlaengerung je Seite", L_ST2, "x", H_ST2, "mm, Sole", T_SOLE, "mm, Schrauben bei X", [round(v_, 1) for v_ in (bx0 - L_ST1 - L_ST2 / 2, bx1 + L_ST1 + L_ST2 / 2)])

# ---------------- Bassdeckel ----------------
DM = ezdxf.readfile(mod("kb_deckel.dxf")).modelspace()
dx = lambda x: x + 370.75 + D; dy = lambda y: y - 0.62 + D
rp = pts_of([e for e in DM.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]')][0])
dk = box(dx(min(a for a, _ in rp)), dx(max(a for a, _ in rp)), dy(min(b for _, b in rp)), dy(max(b for _, b in rp)), 0, T_DECKEL)
for e in DM.query('LWPOLYLINE[layer=="0"]'):
    p = pts_of(e); xs = [a for a, _ in p]; ys = [b for _, b in p]
    if max(ys) - min(ys) < 20 and max(xs) - min(xs) > 300:
        teile.append(("Bassdeckel_Leiste_hinten", box(dx(min(xs)), dx(max(xs)), dy(min(ys)), dy(max(ys)), 0, T_DECKEL).intersect(INNEN)))
        dk = dk.fuse(box(dx(min(a for a, _ in rp)), dx(max(a for a, _ in rp)), dy(max(ys)), dy(min(b for _, b in rp)) + 0.01, 0, T_DECKEL)).clean()   # kein Spalt zur hinteren Leiste
    # (die beiden Zierausschnitte aus kb_deckel.dxf entfallen - Deckeloeffnungen wie bei Markus, siehe unten)
slp = [(dx(a), dy(b)) for a, b in pts_of([e for e in DM.query('LWPOLYLINE[layer=="Def"]')][0])]
dk = dk.cut(prisma_xy(slp, -1, T_DECKEL + 1))
# Deckeloeffnungen wie Markus (Markus-Bass-Verdeck.dxf: 4 Rosetten zu je 9 Loechern, Ø12 / Mitte Ø14), auf das
# Gehaeuse 380 x 205 umgelegt (Markus 382 x 207: je 1 mm zur Mitte). Lufttasten-Schlitz bleibt der eigene (kb_deckel.dxf).
VD_ = ezdxf.readfile(mod("kb_verdeck_markus.dxf")).modelspace()
_ros = sorted({(round(c.dxf.center.x + 381 - 1.0, 2), round(c.dxf.center.y + 206 - 1.0, 2), round(c.dxf.radius, 2)) for c in VD_.query('CIRCLE[layer=="0"]')})
for cx_, cy_, r_ in _ros:
    dk = dk.cut(cq.Solid.makeCylinder(r_, T_DECKEL + 2, cq.Vector(cx_, cy_, -1)))
print("Deckeloeffnungen wie Markus:", len(_ros), "Loecher, X", min(a for a, _, _ in _ros), "-", max(a for a, _, _ in _ros), "Y", min(b for _, b, _ in _ros), "-", max(b for _, b, _ in _ros))
# Deckelauflage rundum (wie Markus 6 x 6): Leisten an Rueckteil, beiden Seitenteilen und - ausserhalb des Knopfbereichs -
# an der Front, vom Deckel (Z T_DECKEL) bis unter den Bassboden (Z_UB): gleichmaessiger Abschluss rundherum.
B_AUF = 6.0
FI_ = B - D
_kx0 = min(k[0] for k in kn) - 9.0; _kx1 = max(k[0] for k in kn) + 9.0
# F-Loecher im Rueckteil (rahmen.dxf LAYER_1, mittig): Auflageleiste dahinter wie Markus grosszuegig ausgespart (+5 mm),
# fuer Gaze innen; an der Deckelseite bleibt ein 4-mm-Rest als Auflage.
RAND_G, REST_F = 5.0, 4.0
_dl_hinten = box(D, L - D, D, D + B_AUF, T_DECKEL, Z_UB)
_RD = ezdxf.readfile(mod("rahmen.dxf")).modelspace(); _nf = 0
for e in _RD.query('SPLINE[layer=="LAYER_1"]'):
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.1)]
    if min(b for _, b in pp_) > -4.24 and max(b for _, b in pp_) < 61.76 - 2:
        xs_ = [L - (a + 382.0 - (382.0 - L) / 2) for a, _ in pp_]        # Rueckteil wird nicht gespiegelt, die Leiste schon
        _dl_hinten = _dl_hinten.cut(box(min(xs_) - RAND_G, max(xs_) + RAND_G, D - 1, D + B_AUF + 1, T_DECKEL + REST_F, Z_UB + 1)); _nf += 1
print("F-Loecher im Rueckteil:", _nf, "- Auflageleiste hinten ausgespart")
teile += [
    ("Deckelauflage_hinten", _dl_hinten),
    ("Deckelauflage_links", box(D, D + B_AUF, D + B_AUF, FI_, T_DECKEL, Z_UB)),
    ("Deckelauflage_rechts", box(L - D - B_AUF, L - D, D + B_AUF, FI_, T_DECKEL, Z_UB)),
    ("Deckelauflage_vorne_links", box(D + B_AUF, _kx0, FI_ - B_AUF, FI_, T_DECKEL, Z_UB)),
    ("Deckelauflage_vorne_rechts", box(_kx1, L - D - B_AUF, FI_ - B_AUF, FI_, T_DECKEL, Z_UB)),
]
# Im Knopfbereich zwei kleine Auflagekloetze (8 breit, 5 hoch) zwischen den Knopfzylindern der unteren Reihe
_r1 = sorted(k[0] for k in kn if k[3] == 1); _lu = [(a_ + b_) / 2 for a_, b_ in zip(_r1[:-1], _r1[1:])]
_kl_x = [min(_lu, key=lambda x_: abs(x_ - (_kx0 + (_kx1 - _kx0) * f_))) for f_ in (1 / 3, 2 / 3)]
for j_, xk_ in enumerate(_kl_x):
    teile.append((f"Deckelauflage_Klotz_vorne_{j_+1}", box(xk_ - 4.0, xk_ + 4.0, FI_ - B_AUF, FI_, T_DECKEL, T_DECKEL + 5.0)))
print("Auflagekloetze vorne bei X", [round(x_, 1) for x_ in _kl_x])
# Fuesse mit Kloetzen (Eckbloecke) wie Markus: in jeder Ecke innerhalb der Leisten ein Block 15 x 15 vom Deckel bis unter
# den Boden; darunter aussen der Fuss (fus2.dxf, Drehteil Ø15 x 21,5) mit Senkkopfschraube Ø3 durch Fuss und Deckel in den
# Block, Gummipfropfen in der Senkung. Fussmitte je 18,5 vom Seitenteil und 15,5 von Rueck-/Frontkante (Deckelleiste hinten).
BLK, D_DS, D_BL, L_FUSS = 15.0, 3.5, 3.0, 21.5
_xi0, _xi1, _yi0, _yi1 = D + B_AUF, L - D - B_AUF, D + B_AUF, FI_ - B_AUF
_ecken = [(_xi0, _yi0, 15.5), (_xi1 - BLK, _yi0, 15.5), (_xi0, _yi1 - BLK, B - 15.5), (_xi1 - BLK, _yi1 - BLK, B - 15.5)]
_aussen = [(0, 0), (0, 7.0), (-0.5, 7.5), (-3.0, 7.5)]
_c, _r = (-8.889409, 11.283499), 7.0
for i in range(1, 30):
    t_ = math.radians(327.0 - (327.0 - 223.0) * i / 30)
    _aussen.append((_c[0] + _r * math.cos(t_), _c[1] + _r * math.sin(t_)))
_aussen += [(-14.0, 6.5), (-16.0, 6.5), (-16.0, 7.5), (-21.5, 7.5), (-21.5, 0)]
_innen = [(0, 0), (0, 6.75), (-2.5, 6.75), (-2.5, 4.0), (-4.0, 4.0), (-4.0, 4.75), (-5.5, 4.75), (-5.5, 4.0), (-7.5, 2.0), (-21.5, 2.0), (-21.5, 0)]
def _drehteil(profil):
    pts = [(r_, -(L_FUSS + x_)) for x_, r_ in profil]
    return cq.Workplane("XZ").polyline(pts).close().revolve(360, (0, 0, 0), (0, 1, 0)).val()
_fuss_roh = _drehteil(_aussen).cut(_drehteil(_innen))
_kuppe = [(0.4 + 0.6 * (1 - (r_ / 6.25) ** 2), r_) for r_ in [6.25 - 6.25 * i / 12 for i in range(13)]]
_gp_roh = _drehteil([(-5.5, 0), (-5.5, 4.75), (-4.0, 4.75), (-4.0, 4.0), (-2.5, 4.0), (-2.5, 6.75), (0, 6.75), (0, 6.25)] + _kuppe)
_ilh = [i for i, (n_, _) in enumerate(teile) if n_ == "Bassdeckel_Leiste_hinten"][0]
for k, (bx, by, cy) in enumerate(_ecken):
    cx = bx + BLK / 2
    bohr = cq.Solid.makeCylinder(D_BL / 2, 20, cq.Vector(cx, cy, T_DECKEL - 1))
    teile.append((f"Eckblock_{k+1}", box(bx, bx + BLK, by, by + BLK, T_DECKEL, Z_UB).cut(bohr)))
    loch = cq.Solid.makeCylinder(D_DS / 2, T_DECKEL + 2, cq.Vector(cx, cy, -1))
    dk = dk.cut(loch)
    teile[_ilh] = (teile[_ilh][0], teile[_ilh][1].cut(loch))
    teile.append((f"Fuss_{k+1}", _fuss_roh.translate(cq.Vector(cx, cy, 0))))
    z_k0, z_k1 = -(L_FUSS - 5.5), -(L_FUSS - 7.5)
    kopf = cq.Solid.makeCone(4.0, 1.5, z_k1 - z_k0, cq.Vector(cx, cy, z_k0), cq.Vector(0, 0, 1))
    sc = cq.Solid.makeCylinder(D_BL / 2, L_FUSS - 7.5 + T_DECKEL + 12, cq.Vector(cx, cy, z_k1)).fuse(kopf).clean()
    teile.append((f"Fuss_{k+1}_Gummipfropfen", _gp_roh.translate(cq.Vector(cx, cy, 0))))
    teile.append((f"Fuss_{k+1}_Schraube", sc))
# Auflage fuer den Deckel hinten: Leiste 2 dick x 6 breit innen unter der schmalen hinteren Deckelleiste eingeleimt,
# halb unter der Leiste, halb unter der Deckelkante; zwischen den hinteren Eckbloecken (nicht ueber die ganze Laenge).
_yfuge = [bb_ for n_, bb_ in [(n_, s_.BoundingBox()) for n_, s_ in teile if n_ == "Bassdeckel_Leiste_hinten"]][0].ymax
teile.append(("Bassdeckel_Auflage_hinten", box(_xi0 + BLK + 0.5, _xi1 - BLK - 0.5, _yfuge - 3.0, _yfuge + 3.0, T_DECKEL, T_DECKEL + 2.0)))
print("Deckelauflage hinten 2 x 6 bei Y", round(_yfuge - 3, 1), "-", round(_yfuge + 3, 1), "X", _xi0 + BLK + 0.5, "-", _xi1 - BLK - 0.5)
print("Deckelauflage 6 x 6 rundum, Front ausserhalb X", round(_kx0, 1), "-", round(_kx1, 1), "| Eckbloecke + Fuesse", len(_ecken), "Fussmitten", [(round(b_ + BLK / 2, 1), round(c_, 1)) for b_, _, c_ in _ecken])
teile.append(("Bassdeckel", dk.intersect(INNEN)))
# Hintere Gehaeusekante (Deckelseite, aussen) wie Markus gerundet: R 20, an den Enden 30 mm auslaufend
R_KANTE, L_AUSLAUF = 20.0, 30.0
lx0, lx1 = 60.0, L - 60.0
def kantenprofil(x_, r_):
    r_ = max(r_, 0.05)
    w_ = cq.Workplane("YZ", origin=(x_, 0, 0)).moveTo(-0.01, -0.01).lineTo(r_, -0.01).threePointArc((r_ * (1 - 0.7071), r_ * (1 - 0.7071)), (-0.01, r_)).close()
    return w_.wires().val()
stationen = [(lx0 - L_AUSLAUF, 0.0), (lx0, R_KANTE), (lx1, R_KANTE), (lx1 + L_AUSLAUF, 0.0)]
kante = None
for (xa, ra), (xb, rb) in zip(stationen[:-1], stationen[1:]):
    seg_ = cq.Solid.makeLoft([kantenprofil(xa, ra), kantenprofil(xb, rb)], True)
    kante = seg_ if kante is None else kante.fuse(seg_)
_gerundet = []
for k_ in range(len(teile)):
    if teile[k_][1].intersect(kante).Volume() > 0.05:
        teile[k_] = (teile[k_][0], teile[k_][1].cut(kante)); _gerundet.append(teile[k_][0])
print("Hintere Kante R", R_KANTE, "X", lx0, "-", lx1, "trifft:", _gerundet)
LT_X = (min(a for a, _ in slp) + max(a for a, _ in slp)) / 2; LT_Y = (min(b for _, b in slp) + max(b for _, b in slp)) / 2

# ---------------- Knoepfe ----------------
D_ZYL, D_KOPF, H_KOPF, L_ZYL, SPALT = 7.5, 15.5, 4.5, 20.0, 6.5
# Blende ueber der inneren Reihe (R2) wie Markus 3. Reihe (blende.dxf, 4 Bohrungen) sinngemaess auf 8 Bohrungen erweitert:
# geschwungene Kontur - Bogen um jede Bohrung, Einzug dazwischen, grosse Endboegen. 7 mm stark, Knoepfe 7 mm erhoeht.
T_BL = 7.0
def blende_kontur(xs_):
    from ezdxf.entities import LWPolyline
    xs_ = sorted(xs_, reverse=True); o_, u_ = 6.66, -4.42; pts = []
    pts.append((xs_[0] - 5.06, o_, -0.301))
    for i, h_ in enumerate(xs_[1:-1]):
        pts += [(h_ + 5.06, o_, 0.388), (h_ - 5.06, o_, -0.301)]
    pts.append((xs_[-1] + 5.06, o_, 2.267)); pts.append((xs_[-1] + 5.06, u_, -0.301))
    for h_ in xs_[-2:0:-1]:
        pts += [(h_ - 5.06, u_, 0.388), (h_ + 5.06, u_, -0.301)]
    pts.append((xs_[0] - 5.06, u_, 2.267))
    e = ezdxf.new().modelspace().add_lwpolyline([(x_, y_, 0, 0, b_) for x_, y_, b_ in pts], format="xyseb", close=True)
    return pts_of(e, 0.05)
x2 = [k[0] for k in kn if k[3] == 2]; zb2 = Z_R[2][0]
bk = [(-x_, zb2 + y_) for x_, y_ in blende_kontur([-x_ for x_ in x2])]   # Ansicht von aussen: x = -X
bl = cq.Workplane("XZ").polyline(bk).close().extrude(-T_BL).val()
bb_ = bl.BoundingBox(); bl = bl.translate(cq.Vector(0, B - bb_.ymin, 0))
for X_, Z_, r_, rr, px in kn:
    if rr == 2: bl = bl.cut(stab((X_, B - 1, Z_), (X_, B + T_BL + 1, Z_), 2 * r_))
teile.append(("Frontblende", bl))
# Knoepfe; Filzeinlage in der Durchfuehrung der aeussersten Lage (wie Markus: Huelse Bohrung/Zylinder)
for X_, Z_, r_, rr, px in kn:
    ya = B + (T_BL if rr == 2 else 0.0); dl_ = T_BL if rr == 2 else D; yk = ya + SPALT; lz = yk - (B + SPALT - L_ZYL)
    teile.append((f"Knopf_R{rr}_{X_:.0f}", stab((X_, yk, Z_), (X_, yk + H_KOPF, Z_), D_KOPF)))
    teile.append((f"Knopfzylinder_R{rr}_{X_:.0f}", stab((X_, yk - lz, Z_), (X_, yk, Z_), D_ZYL)))
    teile.append((f"Filzeinlage_R{rr}_{X_:.0f}", stab((X_, ya - dl_, Z_), (X_, ya, Z_), 2 * r_).cut(stab((X_, ya - dl_ - 1, Z_), (X_, ya + 1, Z_), D_ZYL))))
Y_KN_END = B + SPALT - L_ZYL

# ---------------- Mechanik ----------------
# (Knoepfe [(reihe, x)], Klappenmitte x im Grundriss, Achslinie y im Grundriss)
BASS = [([(1, 109.5)], 36.5, 382.1), ([(2, 119)], 80.5, 379.1), ([(1, 147.5)], 124.5, 382.1), ([(2, 157)], 168.5, 379.1),
        ([(1, 185.5)], 212.5, 382.1), ([(2, 195)] + ([] if VAR45 else [(1, 261.5)]), 256.5, 379.1), ([(2, 233)], 290.5, 385.1), ([(1, 223.5)], 334.5, 388.1)]
BEGL = [([(1, 128.5)], 59.5, 385.1), ([(2, 138)], 95.5, 388.1), ([(1, 166.5)], 131.5, 391.1), ([(1, 204.5)], 167.5, 388.1),
        ([(2, 176)], 203.5, 394.1), ([(2, 214)], 239.5, 391.1), ([(1, 242.5)], 275.5, 394.1), ([(2, 252)], 311.5, 391.1)]
KX = [k[0] for k in kn]
# U-Achsen an der Stirnflaeche der Basskammer-Vorderleiste (Y 118.5): Achse in Einzelschellen (lager-schelle.dxf:
# 8.4 x 5.29, Lagerpunkt 2.5 ueber der Auflage, eine Befestigungsschraube), Schenkel an den Enden aufgebogen.
# Ebenen: Bass unten (kurze Wege zur Druckstange oben), Begleiter oben. Einheitlicher Klappenhub: Uebersetzung
# Klappenweg/Knopfweg = R_UE fuer alle Achsen - Bass ueber die Lage des Klappen-Angriffspunkts, Begleiter ueber die Armlaenge.
Y_FL = Y(385.1); Y_AX = Y_FL + 2.5; Y_SK = Y_AX + 4.0          # Schenkel verlassen das Achspaket schraeg (ohne rechten Winkel)
Z_BASS = [28.5, 35.5, 40.0]; Z_BEGL = [50.2, 55.0]
Z_RB, R_UE, HUB = 51.0, 1.25, 6.0                            # Bass-Druckstangen unter den offenen Begleiterklappen
ach = []
for art, tab, kls in (("Bass", BASS, bass_kl), ("Begleiter", BEGL, begl_kl)):
    for knl, pxk, pyl in tab:
        ki = min(range(len(kls)), key=lambda j: abs((kls[j][0] + kls[j][1]) / 2 - X(pxk)))
        ach.append(dict(art=art, kn=[knopf(rr, px) for rr, px in knl], kl=kls[ki], nm=f"{art}klappe_{ki+1}"))
hx_list = []
KX_BASS = [XK(px) for knl, _, _ in BASS for rr, px in knl]
def frei(x, abst, kx_=None): return all(abs(x - k) >= abst for k in (KX if kx_ is None else kx_)) and all(abs(x - h) >= 5.0 for h in hx_list)
for a in ach:
    k = a["kl"]; xm = (k[0] + k[1]) / 2
    # Begleiter-Klappenschenkel steigen nach oben, Begleiter-Knopfarme gehen nach unten -> nur Bassknoepfe meiden
    kx_ = None
    for d_ in [0] + [s_ * j_ * 0.5 for j_ in range(1, 60) for s_ in (1, -1)]:
        if k[0] + 4 <= xm + d_ <= k[1] - 4 and frei(xm + d_, 5.0, kx_): a["xh"] = xm + d_; break
    else: a["xh"] = xm; print("WARNUNG Hebel ohne freie Lage", a["nm"])
    kx_ = [kk[0] for kk in a["kn"]]
    if max(kx_ + [a["xh"]]) - min(kx_ + [a["xh"]]) < 24:
        kand = [x_ for x_ in [k[0] + 4 + 0.5 * j_ for j_ in range(int((k[1] - k[0] - 8) / 0.5) + 1)] if frei(x_, 3.5)]
        if kand: a["xh"] = max(kand, key=lambda x_: min(abs(x_ - q_) for q_ in kx_))
    hx_list.append(a["xh"])
    xs_ = kx_ + [a["xh"]]; a["x0"], a["x1"] = min(xs_), max(xs_)
for art, lv in (("Bass", Z_BASS), ("Begleiter", Z_BEGL)):     # Ebenen zuteilen: tiefste freie Ebene ohne Ueberlappung
    for a in sorted([a for a in ach if a["art"] == art], key=lambda a: a["x0"]):
        for z_ in lv:
            if all(b_.get("z") != z_ or a["x1"] < b_["x0"] - 3 or b_["x1"] < a["x0"] - 3 for b_ in ach if b_ is not a):
                a["z"] = z_; break
        else: a["z"] = lv[-1]; print("WARNUNG keine freie Ebene", a["nm"])
belegt = []                                                  # (x, z0, z1): Schenkelanfang im Bereich der Schellen (Y < Y_SK)
for a in ach:
    s0_ = 1 if a["art"] == "Bass" else -1
    for kk in a["kn"]: belegt.append((kk[0], min(a["z"], a["z"] + 3 * s0_) - 1.2, max(a["z"], a["z"] + 3 * s0_) + 1.2))
    belegt.append((a["xh"], min(a["z"], a["z"] - 3 * s0_) - 1.2, max(a["z"], a["z"] - 3 * s0_) + 1.2))
# Federn: 3 Windungen, Draht 0.5, an der Aufbiegung; Haken umschliesst den Schenkel, Gegenschenkel liegt an der Stirnflaeche an
DR_F, DI_F, N_F = 0.5, 2.8, 3
D_W = 2.4                                                    # Achsen, Schenkel, Druckstangen 2.4
L_FED = N_F * (DR_F + 0.15) + DR_F
federn = []
for a in ach:
    for seite, x_e in ((1, a["x0"]), (-1, a["x1"])):
        if not any(abs(kk[0] - x_e) < 0.01 for kk in a["kn"]): continue      # Feder nur am Knopfschenkel (Haken ums Gelenk)
        federn.append((x_e + seite * 1.3, x_e + seite * (1.3 + L_FED), a["z"], 1 if a["art"] == "Bass" else -1)); a.setdefault("fed", []).append((federn[-1], x_e, seite))
# Schellen: Block 8.4 (in Z) x 5.29 (in Y) x 3 (in X), Lagerpunkt 2.6 vom Blockende; je nach Platz nach unten oder oben
sch = []
def sch_ok(xs_, z0_, z1_, a):
    if any(abs(xs_ - b_[0]) < 2.95 and b_[1] < z1_ and b_[2] > z0_ for b_ in belegt): return False
    if any(min(f_[0], f_[1]) - 1.7 < xs_ < max(f_[0], f_[1]) + 1.7 and z0_ - 0.8 < max(f_[2] - 2.2, f_[2] + f_[3] * 6.6) and min(f_[2] - 2.2, f_[2] + f_[3] * 6.6, f_[2] + 2.2) < z1_ + 0.8 for f_ in federn): return False
    if any(b_ is not a and b_["x0"] - 1.6 < xs_ < b_["x1"] + 1.6 and z0_ - 1.05 < b_["z"] < z1_ + 1.05 for b_ in ach): return False
    if z0_ < Z_UB or z1_ > z_innen(Y_FL) - 0.5: return False          # Stirnflaeche: Bodenkante + Vorderleiste
    return all(abs(xs_ - s_[0]) >= 3.2 or s_[2] < z0_ or s_[1] > z1_ for s_ in sch)
for a in ach:
    for seite, x_e in ((1, a["x0"]), (-1, a["x1"])):
        ff = [f_ for f_, xe2, s2 in a.get("fed", []) if s2 == seite]
        x_s = (max(ff[0][0], ff[0][1]) if seite > 0 else min(ff[0][0], ff[0][1])) if ff else x_e + seite * 1.3
        ok = False
        for j_ in range(0, 120):
            xs_ = x_s + seite * (1.7 + 0.5 * j_)
            if not (a["x0"] + 2 < xs_ < a["x1"] - 2): break
            for richt in (-1, 1):
                z0_, z1_ = (a["z"] - 5.8, a["z"] + 2.6) if richt < 0 else (a["z"] - 2.6, a["z"] + 5.8)
                if sch_ok(xs_, z0_, z1_, a): sch.append((xs_, z0_, z1_, richt, a["z"])); a.setdefault("sch", []).append(xs_); ok = True; break
            if ok: break
        if not ok: print("WARNUNG keine Schellenlage", a["nm"])
for i, (xs_, z0_, z1_, richt, za) in enumerate(sch):
    s_ = box(xs_ - 1.5, xs_ + 1.5, Y_FL, Y_FL + 5.29, z0_, z1_)
    s_ = s_.cut(stab((xs_ - 3, Y_AX, za), (xs_ + 3, Y_AX, za), D_W + 0.1))
    s_ = s_.cut(box(xs_ - 3, xs_ + 3, Y_FL - 1, Y_AX, za - 0.38, za + 0.38))                 # Klemmschlitz zur Auflage
    zs_ = za + 3.9 * (-1 if richt < 0 else 1)
    s_ = s_.cut(stab((xs_, Y_FL - 1, zs_), (xs_, Y_FL + 6, zs_), 1.5))
    mech.append((f"Lagerschelle_{i+1}", s_))
    mech.append((f"Lagerschelle_{i+1}_Schraube", stab((xs_, Y_FL - 2.8, zs_), (xs_, Y_FL + 5.29, zs_), 1.4).fuse(
        stab((xs_, Y_FL + 5.29, zs_), (xs_, Y_FL + 5.9, zs_), 2.6))))
print("Lagerschellen", len(sch))
_iv = [i for i, (n_, _) in enumerate(teile) if n_ == "Basskammer_Leiste_vorne"][0]
_lv = teile[_iv][1]
for n_, s2 in mech:
    if n_.endswith("_Schraube"): _lv = _lv.cut(s2)
teile[_iv] = ("Basskammer_Leiste_vorne", _lv)
# Gelenk aus zwei Gelenkkoepfen (Fusion "gelenkkopf": Schaft Oe4 x 13 mit Bohrung, Kopf 2 mm stark, Querbohrung Oe1.5
# 3.19 ueber dem Schaft): Gelenkkopf an der Druckstange (Bohrung auf 2.4 fuer den Draht), Gegenstueck = Gabel am
# Schenkeldraht mit je einer Wange (Zylinder) beidseits des Kopfes; Bolzen durch alle drei, beidseits vorstehend.
_gk = cq.importers.importStep(mod("kb_gelenkkopf.step")).solids().vals()
GK = _gk[0].fuse(_gk[1]).clean()
GK_H, GK_L, GK_R = 3.19, 13.0, 2.0                           # Querbohrung ueber Schaftanfang, Schaftlaenge, Schaftradius
GK = GK.cut(stab((0, 0.5, 0), (0, -GK_L - 1, 0), D_W)).cut(stab((-3, GK_H, 0), (3, GK_H, 0), 1.5))
L_GABEL = 10.0
def gelenk(P, s_, draht_s, draht_st, nm):
    x, y, z = P
    kopf = GK.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180).translate(cq.Vector(x, y + GK_H, z)).cut(draht_st)
    mech.append((f"Gelenkkopf_{nm}", kopf))
    g = stab((x - GK_R, y, z), (x + GK_R, y, z), 2 * GK_R).fuse(stab((x, y, z), (x, y, z - s_ * L_GABEL), 2 * GK_R))
    g = g.cut(box(x - 1.05, x + 1.05, y - 2.4, y + 6, *sorted((z - s_ * 2.4, z + s_ * 2.4))))
    g = g.cut(stab((x, y, z - s_ * GK_H), (x, y, z - s_ * (L_GABEL + 1)), D_W)).cut(stab((x - 3, y, z), (x + 3, y, z), 1.5)).cut(draht_s)
    mech.append((f"Gelenkgabel_{nm}", g))
    mech.append((f"Gelenkbolzen_{nm}", stab((x - 2.4, y, z), (x + 2.4, y, z), 1.4)))
def feder(f_, x_e, seite, z, P, s_):
    # Schenkelfeder 3 Windungen: Haken umfasst den Gelenkzylinder (Gabelschaft) des Knopfschenkels,
    # Gegenschenkel liegt an der Stirnflaeche an.
    xa_, xb_ = min(f_[0], f_[1]), max(f_[0], f_[1]); rm = (DI_F + DR_F) / 2
    h_ = cq.Wire.makeHelix((xb_ - xa_ - DR_F) / N_F, xb_ - xa_ - DR_F, rm)
    pr = cq.Wire.makeCircle(DR_F / 2, cq.Vector(rm, 0, 0), cq.Vector(0, 1, 0))
    sp = cq.Solid.sweep(pr, [], h_, isFrenet=True).rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90).translate(cq.Vector(xa_ + DR_F / 2, Y_AX, z))
    x_nah = xa_ + DR_F / 2 if seite > 0 else xb_ - DR_F / 2
    x_fern = xb_ - DR_F / 2 if seite > 0 else xa_ + DR_F / 2
    zr_ = P[2] - s_ * 6.0                                          # U-Haken um den Gabelzylinder des Knopfschenkels,
    r_ring = GK_R + 0.1 + DR_F / 2                                 # auf der Rueckseite (-Y): drueckt den Schenkel nach vorne zurueck
    u_ = [(P[0] + seite * r_ring * math.cos(math.radians(w_)), P[1] - r_ring * math.sin(math.radians(w_)), zr_) for w_ in range(0, 181, 15)]
    ring = pfad(u_, DR_F)
    # Federschenkel = Windungsende naeher beim Arm; Gegenschenkel (fernes Ende) liegt an der Stirnflaeche an
    haken = pfad([(x_nah, Y_AX, z - rm), (P[0] + seite * r_ring, Y_AX + 2.5, z - rm), u_[0]], DR_F)
    gs = pfad([(x_fern, Y_AX, z - rm), (x_fern, Y_FL + DR_F / 2 + 0.05, z - rm + s_ * 0.8), (x_fern, Y_FL + DR_F / 2 + 0.05, z + s_ * 6.0)], DR_F)   # in Hakenrichtung an der Stirnflaeche
    return sp.fuse(ring).fuse(haken).fuse(gs)
for i, a in enumerate(ach):
    z = a["z"]; nm = f"{i+1}"; s_ = 1 if a["art"] == "Bass" else -1          # Richtung des Knopfschenkels (+Z / -Z)
    achse = stab((a["x0"], Y_AX, z), (a["x1"], Y_AX, z), D_W)
    mech.append((f"Achse_{nm}", achse))
    for kk in a["kn"]:
        X_, Zk = kk[0], kk[1]
        if a["art"] == "Bass": P = (X_, Y_SK + 3.0, Z_RB); arm_l = Z_RB - z
        else: arm_l = ((a["kl"][2] + a["kl"][3]) / 2 + 3.0 - Y_AX) / R_UE; P = (X_, Y_SK + 1.5, z - arm_l)
        Q = (X_, P[1], P[2] - s_ * 4.5); A = (X_, P[1], P[2] - s_ * (L_GABEL + 1.5))
        arm = pfad([(X_, Y_AX, z), (X_, Y_SK, z + s_ * 3.0), A, Q], D_W).cut(achse)
        mech.append((f"Schenkel_Knopf_{nm}_R{kk[3]}_{X_:.0f}", arm)); a.setdefault("arm_l", arm_l)
        y_e = Y_KN_END - 1.25; y_r = y_e - abs(Zk - P[2])
        stange = pfad([(X_, P[1] + 4.6, P[2]), (X_, y_r, P[2]), (X_, y_e, Zk)], D_W)   # steckt im Gelenkkopf-Schaft
        mech.append((f"Druckstange_R{kk[3]}_{X_:.0f}", stange))
        gelenk(P, s_, arm, stange, f"{nm}_R{kk[3]}_{X_:.0f}")
    k = a["kl"]; xh = a["xh"]
    # Schenkelende verlaengert, mit kurzem flachem Ende (6 x 4 x 1) unter der Klappe: Druckpunkt verteilt
    L_FUSS, B_FUSS, T_FUSS = 6.0, 4.0, 1.0
    if a["art"] == "Bass":
        zt = Z_KL_U - 0.05 - T_FUSS - 0.9; yt = (k[2] + k[3]) / 2 - 3.0              # Fuss etwas ueber die Klappenmitte
        hp = [(xh, Y_AX, z), (xh, Y_SK, z - 3.0)] + ([(xh, Y_AX + 2.7, Z_UB)] if z - 3.0 > Z_UB + 0.1 else []) + [(xh, Y_FL + 3.0, Z_UB - 4.5), (xh, k[3] + 1.0, zt), (xh, yt - L_FUSS / 2 + 1.2, zt)]
        fuss = box(xh - B_FUSS / 2, xh + B_FUSS / 2, yt - L_FUSS / 2, yt + L_FUSS / 2, Z_KL_U - 0.05 - T_FUSS, Z_KL_U - 0.05)
    else:
        yt = (k[2] + k[3]) / 2 + 3.0
        dzf = (0.05 + D_W / 2) / CA                              # Ende zusaetzlich geknickt: liegt 12 mm flach unter der Klappe
        hp = [(xh, Y_AX, z), (xh, Y_SK, z + 2.0), (xh, yt - 9.0, z_kl_u(yt - 9.0) - dzf), (xh, yt + 3.0, z_kl_u(yt + 3.0) - dzf)]
        fuss = None
    a["yt"] = yt
    _sk = pfad(hp, D_W)
    mech.append((f"Schenkel_Klappe_{nm}", (_sk.fuse(fuss) if fuss is not None else _sk).cut(achse)))
    for (f_, x_e, seite) in a.get("fed", []):                  # Haken ums Gelenk dieses Knopfschenkels
        P_ = (x_e, Y_SK + 3.0, Z_RB) if a["art"] == "Bass" else (x_e, Y_SK + 1.5, z - a["arm_l"])
        mech.append((f"Feder_{nm}_{'ab'[0 if seite > 0 else 1]}", feder(f_, x_e, seite, z, P_, s_)))
    print(f"   Achse {nm}: {' + '.join(f'R{kk[3]} X {kk[0]:.1f}' for kk in a['kn'])} -> {a['nm']} (Schenkel X {xh:.1f}), "
          f"X {a['x0']:.1f}-{a['x1']:.1f}, Z {z:.1f}, Arm {a['arm_l']:.1f}, Angriff Y {yt:.1f}, Klappenhub {abs(yt - Y_AX) / a['arm_l'] * HUB:.1f}, Schellen {[round(x_, 1) for x_ in a.get('sch', [])]}")
# Lufttaste wie Markus (Langloch mit halbrunden Enden, Filzeinlage 1 mm, Taste Kante R1, Hub 10), Masse aus der
# Deckelausnehmung. Lufthebel = Alu-Profil nach Lufthebel-lagerbock-kurz146mm.dxf (Seitenansicht): Lagerbock am
# Begleiterbrett (Brettunterseite), Hebel haengt darunter, sein langer Arm liegt unter der Luftklappe; zwei
# Schenkelfedern spiegelgleich auf dem Lagerbolzen druecken den Arm zur Klappe (zu). Taste -> Schubstange Oe3 mit
# Gabelgelenk ans Hebelauge: Taste hoch -> Auge hoch -> langer Arm runter -> Klappe auf.
LT_L, LT_B = max(a for a, _ in slp) - min(a for a, _ in slp), max(b for _, b in slp) - min(b for _, b in slp)
T_FILZ_LT, HUB_LT = 1.0, 10.0
def langloch(l_, b_, z0_, h_):
    return cq.Workplane("XY").center(LT_X, LT_Y).slot2D(l_, b_).extrude(h_).translate((0, 0, z0_)).val()
mech_lt_filz = langloch(LT_L, LT_B, 0, T_DECKEL).cut(langloch(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT, -1, T_DECKEL + 2))
H_TASTE = T_DECKEL + HUB_LT
taste = cq.Workplane("XY").center(LT_X, LT_Y).slot2D(LT_L - 2 * T_FILZ_LT, LT_B - 2 * T_FILZ_LT).extrude(H_TASTE).faces("<Z").edges().fillet(1.0).val()
taste = taste.translate(cq.Vector(0, 0, -HUB_LT))                    # Ruhelage: steht 10 mm aus dem Deckel
teile.append(("Lufttaste_Filz", mech_lt_filz)); mech.append(("Lufttaste", taste))
_idk = [i for i, (n_, _) in enumerate(teile) if n_ == "Bassdeckel"][0]
teile[_idk] = ("Bassdeckel", teile[_idk][1].cut(langloch(LT_L, LT_B, -1, T_DECKEL + 2)))   # Ausnehmung = Langloch der Filzeinlage
# Ledereinsatz 0.5 wie Markus: Feld zwischen den Schalloeffnungs-Rosetten (4 mm Abstand), hinteres Drittel frei, vorne bis
# an die Front, hintere Ecken R10; Leder in die Lufttasten-Ausnehmung hineingefuehrt; im Deckel eingelassen.
T_LEDER_D, R_LEDER, A_SCHALL = 0.5, 10.0, 4.0
ly0, ly1 = D + (B - 2 * D) / 3.0, B - D
lx0 = max(x_ + r_ for x_, _, r_ in _ros if x_ < L / 2) + A_SCHALL
lx1 = min(x_ - r_ for x_, _, r_ in _ros if x_ > L / 2) - A_SCHALL
feld = (cq.Workplane("XY").center((lx0 + lx1) / 2, (ly0 + ly1 + 2 * R_LEDER) / 2).rect(lx1 - lx0, ly1 - ly0 + 2 * R_LEDER)
        .extrude(T_DECKEL).edges("|Z").fillet(R_LEDER - 0.01).val()).intersect(box(lx0, lx1, ly0, ly1, 0, T_DECKEL))
leder_fl = feld.intersect(box(lx0, lx1, ly0, ly1, 0, T_LEDER_D))
leder_sl = langloch(LT_L + 2 * T_LEDER_D, LT_B + 2 * T_LEDER_D, 0, T_DECKEL).cut(langloch(LT_L, LT_B, -1, T_DECKEL + 2))
leder_d = leder_fl.cut(langloch(LT_L, LT_B, -1, T_DECKEL + 2)).fuse(leder_sl).clean()
teile[_idk] = ("Bassdeckel", teile[_idk][1].cut(leder_d))
teile.append(("Bassdeckel_Leder", leder_d))
print("Deckel-Ledereinsatz X", round(lx0, 1), "-", round(lx1, 1), "Y", round(ly0, 1), "-", round(ly1, 1), "R", R_LEDER)
LH = ezdxf.readfile(mod("kb_lufthebel.dxf")).modelspace()
# Lufthebel-Gruppe (Hebel, Stuetzen, Bolzen, Federn, Gabelgelenk) wird flach aufgebaut (Brettunterseite waagrecht bei
# Z_B0) und dann um die Brettneigung gekippt: im Schnitt senkrecht zum Begleiterbrett.
LHX = LT_X + 12.65                                                 # Zeichnung x -> X: Hebelauge (x -12.65) ueber der Tastenmitte
T_LH = 3.0
def luft_gruppe(LY):
    grp = []
    Z_B0 = z_innen(LY)
    Z_LH = Z_B0 - (T_KL + T_LEDER) - 0.05 + 5.35                      # Zeichnung y -> Z = Z_LH - y (langer Arm an der Klappe)
    def lh_xz(p_): return [(a + LHX, Z_LH - b) for a, b in p_]
    def platte_xz(pts, y0, y1): return cq.Workplane("XZ").polyline(pts).close().extrude(-(y1 - y0)).translate((0, y0, 0)).val()
    hp_ = lh_xz(pts_of([e for e in LH.query('LWPOLYLINE[layer=="2-Gelenk-makierung"]')][0], 0.05))
    kr = [(c.dxf.center.x + LHX, Z_LH - c.dxf.center.y, c.dxf.radius) for c in LH.query('CIRCLE[layer=="0"]')]
    XP, ZP = [(k[0], k[1]) for k in kr if abs(k[2] - 1.0) < 0.01 and abs(k[0] - (-55.0 + LHX)) < 0.6][0]
    XA, ZA_ = [(k[0], k[1]) for k in kr if abs(k[2] - 0.54) < 0.01][0]
    hebel = platte_xz(hp_, LY - T_LH / 2, LY + T_LH / 2).fuse(stab((XP, LY - T_LH / 2, ZP), (XP, LY + T_LH / 2, ZP), 5.0))
    hebel = hebel.cut(stab((XP, LY - 3, ZP), (XP, LY + 3, ZP), 2.05)).cut(stab((XA, LY - 3, ZA_), (XA, LY + 3, ZA_), 1.08))
    grp.append(("Lufthebel", hebel))
    # Federn direkt neben dem Alu-Hebel auf dem Bolzen, Stuetzen (Profil "2-Gelenke", je 8 breit, 2 Schrauben) aussen daneben
    L_FL = N_F * (DR_F + 0.15) + DR_F
    W_LB = 8.0
    y_f = [(LY - T_LH / 2 - 0.15 - L_FL, LY - T_LH / 2 - 0.15), (LY + T_LH / 2 + 0.15, LY + T_LH / 2 + 0.15 + L_FL)]
    y_w = [(y_f[0][0] - 0.15 - W_LB, y_f[0][0] - 0.15), (y_f[1][1] + 0.15, y_f[1][1] + 0.15 + W_LB)]
    bk = [p_ for p_ in [pts_of(e, 0.05) for e in LH.query('LWPOLYLINE[layer=="2-Gelenke"]')] if min(a for a, _ in p_) < -40][0]
    xb0, xb1 = min(a for a, _ in lh_xz(bk)), max(a for a, _ in lh_xz(bk))
    x_schr = -64.0 + LHX
    schr = []
    for k_, (y0_, y1_) in enumerate(y_w):
        w_ = platte_xz(lh_xz(bk), y0_, y1_).fuse(box(xb0, xb1, y0_, y1_, Z_LH - 0.5, Z_B0 + 3.0))
        w_ = w_.cut(box(xb0 - 1, xb1 + 1, y0_ - 1, y1_ + 1, Z_B0, Z_B0 + 10)).cut(stab((XP, LY - 30, ZP), (XP, LY + 30, ZP), 2.05))
        for dy_ in (-2.0, 2.0):
            yc_ = (y0_ + y1_) / 2 + dy_
            sc_ = stab((x_schr, yc_, Z_LH - 2.53 - 0.9), (x_schr, yc_, Z_LH - 2.53), 2.6).fuse(stab((x_schr, yc_, Z_LH - 2.53), (x_schr, yc_, Z_B0 + 3.5), 1.4))
            w_ = w_.cut(stab((x_schr, yc_, Z_LH - 4), (x_schr, yc_, Z_B0 + 1), 1.5)); schr.append(sc_)
        grp.append((f"Lufthebel_Stuetze_{k_+1}", w_))
    for j_, sc_ in enumerate(schr): grp.append((f"Lufthebel_Stuetze_Schraube_{j_+1}", sc_))
    grp.append(("Lufthebel_Bolzen", stab((XP, y_w[0][0] - 0.6, ZP), (XP, y_w[1][1] + 0.6, ZP), 2.0)))
    def kante_oben(x_):
        zz_ = []
        for (xa_, za_), (xb_, zb_) in [(hp_[i], hp_[i + 1]) for i in range(len(hp_) - 1)]:
            if min(xa_, xb_) <= x_ <= max(xa_, xb_) and abs(xb_ - xa_) > 1e-6: zz_.append(za_ + (zb_ - za_) * (x_ - xa_) / (xb_ - xa_))
        return max(zz_)
    rmf = (DI_F + DR_F) / 2
    for k_, ((ya0, ya1), sg) in enumerate(((y_f[0], -1), (y_f[1], 1))):
        y_in = ya1 if sg < 0 else ya0; y_out = ya0 if sg < 0 else ya1
        h_ = cq.Wire.makeHelix((L_FL - DR_F) / N_F, L_FL - DR_F, rmf)
        pr = cq.Wire.makeCircle(DR_F / 2, cq.Vector(rmf, 0, 0), cq.Vector(0, 1, 0))
        sp = cq.Solid.sweep(pr, [], h_, isFrenet=True).rotate(cq.Vector(0, 0, 0), cq.Vector(1, 0, 0), -90 * sg)
        sp = sp.translate(cq.Vector(XP, y_in + sg * DR_F / 2, ZP))
        p0 = (XP + rmf, y_in + sg * DR_F / 2, ZP); p1 = (XP + rmf, y_out - sg * DR_F / 2, ZP)
        x_u = XP + 15.0 + 0.6 * sg                                   # U-Haken ueber den kurzen Schenkel (zur Taste)
        z_u = max(kante_oben(XP + 14.4), kante_oben(XP + 15.6)) + DR_F / 2 + 0.05 + (DR_F + 0.15 if sg > 0 else 0.0)   # aeusserer Haken ueber dem inneren
        y_a, y_b = LY + sg * (T_LH / 2 + DR_F / 2 + 0.05), LY - sg * (T_LH / 2 + DR_F / 2 + 0.05)
        haken = pfad([p0, (XP + rmf + 3.0, y_a, ZP + 0.5), (XP + 13.4, y_a, z_u), (x_u, y_a, z_u), (x_u, y_b, z_u), (x_u, y_b, z_u - 2.0)], DR_F)
        zb_ = Z_B0 - DR_F / 2 - 0.05                                 # Gegenschenkel: liegt an der Brettunterseite an, verlaengert
        s1 = pfad([p1, (XP + rmf + 3.0, y_out - sg * DR_F / 2, zb_), (XP + rmf + 12.0, y_out - sg * DR_F / 2, zb_)], DR_F)   # gleiche Richtung wie der Haken
        grp.append((f"Lufthebel_Feder_{k_+1}", sp.fuse(haken).fuse(s1)))
    G_L, D_GS = 10.0, 6.0
    gg = stab((XA, LY, ZA_ + 3.0), (XA, LY, ZA_ - G_L), D_GS).cut(box(XA - 4, XA + 4, LY - T_LH / 2 - 0.4, LY + T_LH / 2 + 0.4, ZA_ - 3.0, ZA_ + 4))
    gg = gg.cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6)).cut(stab((XA, LY, ZA_ - 3.5), (XA, LY, ZA_ - G_L - 1), 3.0))   # Bolzenloch mit Spiel, Bohrung fuer die Schubstange
    grp.append(("Lufttaste_Gabelgelenk", gg))
    grp.append(("Lufttaste_Gelenkbolzen", stab((XA, LY - D_GS / 2, ZA_), (XA, LY + D_GS / 2, ZA_), 1.0)))
    ax0 = cq.Vector(0, LY, Z_B0); ax1 = cq.Vector(1, LY, Z_B0)
    kipp = lambda s_: s_.rotate(ax0, ax1, ALPHA)
    grp = [(n_, kipp(s_)) for n_, s_ in grp]
    E = cq.Vector(XA, LY, ZA_ - 3.6); E = E.rotateX(0) if False else E
    def kp(v_):                                                      # Punkt um die Achse kippen
        return cq.Vector(*cq.Solid.makeSphere(0.01, v_).rotate(ax0, ax1, ALPHA).Center().toTuple())
    E = kp(cq.Vector(XA, LY, ZA_ - 3.6)); F_ = kp(cq.Vector(XA, LY, ZA_ - 50.0)); d_ = (F_ - E).normalized()
    z_end = T_DECKEL - HUB_LT - 0.5
    t_end = (E.z - z_end) / (-d_.z)
    Pb = E + d_ * t_end
    return grp, E, Pb, (XP, ZP, XA, ZA_)
LY = (luft_kl[2] + luft_kl[3]) / 2                                 # Hebel trifft die Luftklappe mittig
grp, E, Pb, _info = luft_gruppe(LY)
# Schubstange zur Tastenmitte geneigt; Gabelgelenk um das Hebelauge mitgedreht
Pb = cq.Vector(Pb.x, LT_Y, Pb.z)
_XP, _ZP, _XA, _ZA = _info
_E0 = cq.Solid.makeSphere(0.01, cq.Vector(_XA, LY, _ZA)).rotate(cq.Vector(0, LY, z_innen(LY)), cq.Vector(1, LY, z_innen(LY)), ALPHA).Center()
_d_alt = (E - _E0).normalized(); _d_neu = (Pb - _E0).normalized()
_dw = math.degrees(math.atan2(_d_neu.y, -_d_neu.z) - math.atan2(_d_alt.y, -_d_alt.z))
grp = [(n_, s_.rotate(_E0, _E0 + cq.Vector(1, 0, 0), _dw) if n_ == "Lufttaste_Gabelgelenk" else s_) for n_, s_ in grp]
E = _E0 + _d_neu * (E - _E0).Length
_gb = [s_ for n_, s_ in grp if n_ == "Lufttaste_Gabelgelenk"][0]
print("Gabel-Mitte", [round(v, 2) for v in _gb.Center().toTuple()], "Stangenrichtung", [round(v, 3) for v in _d_neu.toTuple()])
mech.extend(grp)
schr_lh = [s_ for n_, s_ in grp if "Schraube" in n_]
_ib = [i for i, (n_, _) in enumerate(teile) if n_ == "Begleiterbrett"][0]
_br = teile[_ib][1]
for sc_ in schr_lh: _br = _br.cut(sc_)
teile[_ib] = ("Begleiterbrett", _br)
mech.append(("Lufttaste_Schubstange", stab(E.toTuple(), Pb.toTuple(), 3.0).cut(taste)))
print(f"Lufttaste {LT_L:.1f} x {LT_B:.1f}, Hebelebene Y {LY:.1f} = Klappenmitte (gekippt {ALPHA:.1f} Grad), Schubstange zur Tastenmitte Y {Pb.y:.1f} um {_dw:.1f} Grad geneigt")

# ---------------- Variante 4/5: zusaetzliche halbe Kammer vorne quer in der Ecke (Front + Seitenteil) ----------------
# Begleiterbrett, seitliche Leiste und Frontleiste dort ausgeschnitten; Kammer liegt an Front und Seitenteil an.
# Einfache Helikonplatte 33 x 98 (zwei lange Zungen), Keil, Teiler, Haken wie bei den Basskammern; eigener Boden mit
# Klappenoeffnung und Klappe. Der bisher an Achse 6 gekoppelte Knopf (Reihe 1, Zeichnung x 261.5) bekommt eine eigene
# Mechanik an der Rueckseite der Kammer: Achse in zwei Schellen, Knopfschenkel nach unten (Drehung wie Begleiter),
# Klappenschenkel unter die Klappe, Feder mit Haken um den Gelenkzylinder.
if VAR45:
    XK0_5, XK1_5, YK0_5, YK1_5 = 271.0, L - D, 162.5, B - D
    T_W5, T_E5, ZB5, T_B5 = 3.5, 6.0, 30.0, 3.5
    ZO5 = ZB5 + T_B5; ZT5 = ZO5 + 47.3
    # Klappe wie bei den schmalen Kammern (25 x 33, Oeffnung 19 x 27, Rand 3), quer: 33 laengs der Kammer;
    # Oeffnung beginnt 1 mm hinter der Endleiste am dicken Plattenende (wie bei den Basskammern)
    _yc5 = ((YK0_5 + 3.5) + (YK1_5 - 3.5)) / 2
    OEFF5 = (XK0_5 + 6.0 + 1.0, XK0_5 + 6.0 + 1.0 + 27.0, _yc5 - 9.5, _yc5 + 9.5)
    KL5 = (OEFF5[0] - 3.0, OEFF5[1] + 3.0, OEFF5[2] - 3.0, OEFF5[3] + 3.0)
    for _n in ("Begleiterbrett", "Basskammer_Leiste_rechts", "Begleiterbrett_Leiste_vorne"):
        _k = [k for k, (n_, _) in enumerate(teile) if n_ == _n][0]
        teile[_k] = (_n, teile[_k][1].cut(box(XK0_5, L + 1, YK0_5, B + 1, -1, 200)))
    boden5 = box(XK0_5, XK1_5, YK0_5, YK1_5, ZB5, ZO5).cut(box(OEFF5[0], OEFF5[1], OEFF5[2], OEFF5[3], ZB5 - 1, ZO5 + 1))
    wand_h5 = box(XK0_5, XK1_5, YK0_5, YK0_5 + T_W5, ZO5, ZT5)
    wand_v5 = box(XK0_5, XK1_5, YK1_5 - T_W5, YK1_5, ZO5, ZT5)
    lei_a5 = box(XK0_5, XK0_5 + T_E5, YK0_5 + T_W5, YK1_5 - T_W5, ZO5, ZT5)
    lei_b5 = box(XK1_5 - T_E5, XK1_5, YK0_5 + T_W5, YK1_5 - T_W5, ZO5, ZT5)
    # Platte 33 x 98: X 274..372 (je 3 mm auf den Endleisten), Y zwischen den Wandmitten minus 0.5; gedreht: dickes Ende
    # und Teilerfuss bei X klein, ueber der Klappe (wie bei den Basskammern vorne)
    PX0_5, PX1_5 = XK0_5 + 3.0, XK1_5 - 3.0
    PY0_5, PY1_5 = YK0_5 + T_W5 / 2 + 0.5, YK1_5 - T_W5 / 2 - 0.5
    pl5 = box(PX0_5, PX1_5, PY0_5, PY1_5, ZT5, ZT5 + T_RAND)
    a0_, a1_ = XK0_5 + T_E5 + 0.5, XK1_5 - T_E5 - 0.5
    prof5 = [(a0_, ZT5), (a1_, ZT5), (a1_, ZT5 + T_RAND - T_DUENN), (a0_, ZT5 + T_RAND - T_DICK)]
    pl5 = pl5.fuse(cq.Workplane("XZ").polyline(prof5).close().extrude(-(PY1_5 - PY0_5 - 2 * RUECK)).translate((0, PY0_5 + RUECK, 0)).val())
    br5 = PY1_5 - PY0_5
    for (o_, sb, sll) in [(br5 / 2 - 1.4 - 9.2, 9.2, 71.5), (br5 / 2 + 1.4, 9.2, 71.5)]:
        pl5 = pl5.cut(box(PX0_5 + 5.5, PX0_5 + 5.5 + sll, PY0_5 + o_, PY0_5 + o_ + sb, ZT5 - T_DICK - 1, ZT5 + T_RAND + 1))
    # Teiler (gleiche Kurve wie in den Basskammern) laengs der Kammer, in die Laengswaende genutet
    prof_t5 = [(PX0_5 + (a + 94.88), ZO5 + (b - 0.63)) for a, b in kurve]
    prof_t5[0] = (prof_t5[0][0], ZO5 + 0.01)
    prof_t5 = prof_t5 + [(x_, z_ + T_TEIL) for x_, z_ in reversed(prof_t5)]
    def _teiler5(y0_, y1_): return cq.Workplane("XZ").polyline(prof_t5).close().extrude(-(y1_ - y0_)).translate((0, y0_, 0)).val()
    teiler5 = _teiler5(YK0_5 + T_W5 - T_NUT, YK1_5 - T_W5 + T_NUT)
    _nut5 = _teiler5(YK0_5 + T_W5 - T_NUT - 0.05, YK1_5 - T_W5 + T_NUT + 0.05)
    wand_h5 = wand_h5.cut(_nut5); wand_v5 = wand_v5.cut(_nut5)
    # Haken an beiden Endleisten (langer Schenkel in die Leiste, kurzer Schenkel auf dem Plattenrand)
    haken5 = []
    for xh_, dx_, l_lang, nm_ in ((PX0_5 - 0.1 - D_HAKEN / 2, 1, 8.5, "a"), (PX1_5 + 0.1 + D_HAKEN / 2, -1, 15.0, "b")):
        yc_ = (PY0_5 + PY1_5) / 2
        ZH5 = ZT5 + T_RAND + D_HAKEN / 2
        hk = stab((xh_, yc_, ZH5 - l_lang), (xh_, yc_, ZH5 + 0.3), D_HAKEN).fuse(stab((xh_ - dx_ * 0.3, yc_, ZH5), (xh_ + dx_ * L_KURZ, yc_, ZH5), D_HAKEN)).clean()
        loch_ = stab((xh_, yc_, ZH5 - l_lang - 0.5), (xh_, yc_, ZT5 + 0.1), D_HAKEN)
        if nm_ == "a": lei_a5 = lei_a5.cut(loch_)
        else: lei_b5 = lei_b5.cut(loch_)
        haken5.append((f"Kammer5_Stimmplatte_Haken_{nm_}", hk))
    # Klappe unter dem Kammerboden
    Z_KL_U5 = ZB5 - T_KL - T_LEDER
    kl5 = [("Kammer5_Klappe_Leder", box(KL5[0], KL5[1], KL5[2], KL5[3], ZB5 - T_LEDER, ZB5)),
           ("Kammer5_Klappe_Filz", box(KL5[0] + 0.5, KL5[1] - 0.5, KL5[2] + 0.5, KL5[3] - 0.5, ZB5 - T_LEDER - T_FILZ, ZB5 - T_LEDER)),
           ("Kammer5_Klappe", box(KL5[0] + 0.5, KL5[1] - 0.5, KL5[2] + 0.5, KL5[3] - 0.5, Z_KL_U5, ZB5 - T_LEDER - T_FILZ))]
    # eigene Mechanik an der Kammer-Rueckseite (Wandflaeche Y YK0_5)
    Y_W5 = YK0_5; Y_AX5 = Y_W5 - 2.5; Z_AX5 = ZO5 + 6.5
    kk5 = knopf(1, 261.5); X_k5, Zk5 = kk5[0], kk5[1]
    P5 = (X_k5, Y_AX5, Z_KL_U5 - 1.0); arm5 = Z_AX5 - P5[2]
    xh5 = (KL5[0] + KL5[1]) / 2; yt5 = (KL5[2] + KL5[3]) / 2
    achse5 = stab((X_k5, Y_AX5, Z_AX5), (xh5, Y_AX5, Z_AX5), D_W)
    mech.append(("Achse_K5", achse5))
    A5 = (X_k5, Y_AX5, P5[2] + L_GABEL + 1.5); Q5 = (X_k5, Y_AX5, P5[2] + 4.5)
    arm5_s = pfad([(X_k5, Y_AX5, Z_AX5), A5, Q5], D_W).cut(achse5)
    mech.append((f"Schenkel_Knopf_K5_R1_{X_k5:.0f}", arm5_s))
    y_e5 = Y_KN_END - 1.25; y_r5 = y_e5 - abs(Zk5 - P5[2])
    stange5 = pfad([(X_k5, P5[1] + 4.6, P5[2]), (X_k5, y_r5, P5[2]), (X_k5, y_e5, Zk5)], D_W)
    mech.append((f"Druckstange_R1_{X_k5:.0f}", stange5))
    gelenk(P5, -1, arm5_s, stange5, f"K5_R1_{X_k5:.0f}")
    L_FUSS5, B_FUSS5, T_FUSS5 = 6.0, 4.0, 1.0
    zt5 = Z_KL_U5 - 0.05 - T_FUSS5 - 0.9
    hp5 = [(xh5, Y_AX5, Z_AX5), (xh5, Y_AX5, zt5), (xh5, KL5[2] + 1.0, zt5), (xh5, yt5 + L_FUSS5 / 2 - 1.2, zt5)]
    fuss5 = box(xh5 - B_FUSS5 / 2, xh5 + B_FUSS5 / 2, yt5 - L_FUSS5 / 2, yt5 + L_FUSS5 / 2, Z_KL_U5 - 0.05 - T_FUSS5, Z_KL_U5 - 0.05)
    mech.append(("Schenkel_Klappe_K5", pfad(hp5, D_W).fuse(fuss5).cut(achse5)))
    # Feder am Knopfschenkel: 3 Windungen, Haken um den Gelenkzylinder (Rueckseite), Gegenschenkel an der Kammerwand
    xa_, xb_ = X_k5 + 1.3, X_k5 + 1.3 + L_FED; rm5 = (DI_F + DR_F) / 2
    h_ = cq.Wire.makeHelix((xb_ - xa_ - DR_F) / N_F, xb_ - xa_ - DR_F, rm5)
    pr_ = cq.Wire.makeCircle(DR_F / 2, cq.Vector(rm5, 0, 0), cq.Vector(0, 1, 0))
    sp5 = cq.Solid.sweep(pr_, [], h_, isFrenet=True).rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90).translate(cq.Vector(xa_ + DR_F / 2, Y_AX5, Z_AX5))
    r_ring5 = GK_R + 0.1 + DR_F / 2; zr5 = P5[2] + 6.0
    u5 = [(X_k5 + r_ring5 * math.cos(math.radians(w_)), P5[1] - r_ring5 * math.sin(math.radians(w_)), zr5) for w_ in range(0, 181, 15)]
    haken_f5 = pfad([(xa_ + DR_F / 2, Y_AX5, Z_AX5 - rm5), (X_k5 + r_ring5, Y_AX5, Z_AX5 - rm5 - 1.5), u5[0]], DR_F)
    gs5 = pfad([(xb_ - DR_F / 2, Y_AX5, Z_AX5 - rm5), (xb_ - DR_F / 2, Y_W5 - DR_F / 2 - 0.05, Z_AX5 - rm5 - 0.8), (xb_ - DR_F / 2, Y_W5 - DR_F / 2 - 0.05, Z_AX5 - 6.0)], DR_F)
    mech.append(("Feder_K5", sp5.fuse(pfad(u5, DR_F)).fuse(haken_f5).fuse(gs5)))
    # zwei Schellen an der Kammer-Rueckwand (je eine Schraube in die Wand)
    for j_, xs_ in enumerate((xb_ + 1.7 + 1.5, xh5 - 3.0)):
        s_ = box(xs_ - 1.5, xs_ + 1.5, Y_W5 - 5.29, Y_W5, Z_AX5 - 5.8, Z_AX5 + 2.6)
        s_ = s_.cut(stab((xs_ - 3, Y_AX5, Z_AX5), (xs_ + 3, Y_AX5, Z_AX5), D_W + 0.1))
        s_ = s_.cut(box(xs_ - 3, xs_ + 3, Y_AX5, Y_W5 + 1, Z_AX5 - 0.38, Z_AX5 + 0.38))
        zs_ = Z_AX5 - 3.9
        s_ = s_.cut(stab((xs_, Y_W5 - 6, zs_), (xs_, Y_W5 + 1, zs_), 1.5))
        sc_ = stab((xs_, Y_W5 + 2.3, zs_), (xs_, Y_W5 - 5.29, zs_), 1.4).fuse(stab((xs_, Y_W5 - 5.29, zs_), (xs_, Y_W5 - 5.9, zs_), 2.6))
        wand_h5 = wand_h5.cut(sc_)
        mech.append((f"Lagerschelle_K5_{j_+1}", s_)); mech.append((f"Lagerschelle_K5_{j_+1}_Schraube", sc_))
    teile += [("Kammer5_Boden", boden5), ("Kammer5_Wand_hinten", wand_h5), ("Kammer5_Wand_vorne", wand_v5),
              ("Kammer5_Leiste_a", lei_a5), ("Kammer5_Leiste_b", lei_b5), ("Kammer5_Teiler", teiler5),
              ("Kammer5_Stimmplatte", pl5)] + haken5 + kl5
    print(f"Variante 4/5: Zusatzkammer X {XK0_5}-{XK1_5} Y {YK0_5}-{YK1_5} Z {ZB5}-{ZT5 + T_RAND:.1f}, Platte {PX1_5 - PX0_5:.0f} x {br5:.0f}, "
          f"Knopf R1 X {X_k5:.1f} eigene Achse Y {Y_AX5} Z {Z_AX5}, Arm {arm5:.1f}, Klappenhub {(yt5 - Y_AX5) / arm5 * HUB:.1f}")

# ---------------- Fliegengitter (Gaze) innen ueber Rosetten, F-Loechern und Hornbohrungen ----------------
# 0.3 dick, je 3 mm groesser als die Oeffnung (Kleberand); angrenzende feste Teile bekommen eine 0.3-mm-Aussparung.
T_GZ, R_GZ = 0.3, 3.0
from shapely.geometry import Polygon as _Poly, Point as _Pt
sys.path.insert(0, H); from gaze2d import gaze_platte
gaze = []; _g2d = {}
_mitten = [(x_, y_) for x_, y_, r_ in _ros if r_ > 6.5]
for j_, (mx_, my_) in enumerate(sorted(_mitten)):
    rr_ = max(math.hypot(x_ - mx_, y_ - my_) + r_ for x_, y_, r_ in _ros if math.hypot(x_ - mx_, y_ - my_) < 30) + R_GZ
    gaze.append((f"Gaze_Rosette_{j_+1}", cq.Solid.makeCylinder(rr_, T_GZ, cq.Vector(mx_, my_, T_DECKEL)), "Bassdeckel"))
    _g2d[gaze[-1][0]] = (_Pt(mx_, my_).buffer(rr_, 32), "XY", T_DECKEL)
_RD2 = ezdxf.readfile(mod("rahmen.dxf")).modelspace(); _jf = 0
for e in _RD2.query('SPLINE[layer=="LAYER_1"]'):
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(0.2)]
    if min(b for _, b in pp_) > -4.24 and max(b for _, b in pp_) < 61.76 - 2:
        pg_ = _Poly([(L - (a + 382.0 - (382.0 - L) / 2), b + 4.24) for a, b in pp_]).buffer(R_GZ).simplify(0.05)
        _jf += 1
        gaze.append((f"Gaze_FLoch_{_jf}", cq.Workplane("XZ").polyline(list(pg_.exterior.coords)[:-1]).close().extrude(-T_GZ).translate((0, D, 0)).val(), "Rueckteil"))
        _g2d[gaze[-1][0]] = (pg_, "XZ", D)
for j_, (X_, Z_) in enumerate(hoerner_xz):
    if VAR45 and XK0_5 - 9 < X_ < L and ZO5 - 9 < Z_ < ZT5: continue      # von der Zusatzkammer verschlossen: kein Gitter
    gaze.append((f"Gaze_Horn_{j_+1}", cq.Solid.makeCylinder(9.0 + R_GZ, T_GZ, cq.Vector(X_, B - D - T_GZ, Z_), cq.Vector(0, 1, 0)), "Front"))
    _g2d[gaze[-1][0]] = (_Pt(X_, Z_).buffer(9.0 + R_GZ, 32), "XZ", B - D - T_GZ)
_aus = []
for k_ in range(len(teile)):
    n_, s_ = teile[k_]
    if n_ in RAHMEN_ or n_ == "Bassdeckel": continue
    for gn_, gs_, _ in gaze:
        if s_.intersect(gs_).Volume() > 1e-3:
            s_ = s_.cut(gs_); _aus.append((n_, gn_))
    teile[k_] = (n_, s_)
# Gaze vereinfacht: Platte 0.3 mit quadratischen Durchbruechen, Raster 4 mm / Steg 0.4 (echt: 2 mm / Faden 0.2)
RASTER, D_FADEN = 4.0, 0.4
teile += [(gn_, gaze_platte(_g2d[gn_][0], _g2d[gn_][1], _g2d[gn_][2], T_GZ, RASTER, D_FADEN)[0]) for gn_, gs_, tr_ in gaze]
print("Gaze:", len(gaze), "Stueck, Aussparungen in", sorted({a_ for a_, _ in _aus}))

# ---------------- Laenge und Zusammenbau ----------------
BASS_DL = L_END - 380.0
exec(open(os.path.join(H, "streck.py"), encoding="utf-8").read().replace('FRONT_EBENE = float(_os.environ.get("BASS_FRONT_EBENE", "50.0"))', 'FRONT_EBENE = 23.0'))
a_ = cq.Assembly(name=f"Bassteil_Standardmechanik_{NAME_V}_Helikon")
FARBEN = [   # (Praefix, Farbe) - jede Teileart eigene Farbe, erste passende Zeile gilt
    ("Gaze", (0.78, 0.80, 0.82)), ("Frontblende", (0.45, 0.25, 0.12)), ("Front", (0.78, 0.60, 0.38)), ("Rueck", (0.74, 0.56, 0.36)), ("Seiten", (0.70, 0.52, 0.33)),
    ("Balgrahmen_Einsatzleiste", (0.60, 0.45, 0.28)), ("Balgrahmen", (0.86, 0.70, 0.48)), ("Horn", (0.70, 0.55, 0.30)),
    ("Bassboden", (0.90, 0.82, 0.62)), ("Deckelauflage", (0.60, 0.45, 0.28)), ("Eckblock", (0.55, 0.40, 0.25)), ("Fuss", (0.28, 0.20, 0.14)), ("Begleiterbrett_Leiste", (0.60, 0.45, 0.28)), ("Begleiterbrett", (0.84, 0.86, 0.60)), ("Bassdeckel_Leder", (0.12, 0.08, 0.05)), ("Bassdeckel_Auflage", (0.60, 0.45, 0.28)), ("Bassdeckel", (0.30, 0.22, 0.15)),
    ("Kammer5_Stimmplatte_Haken", (0.40, 0.40, 0.42)), ("Kammer5_Stimmplatte", (0.80, 0.82, 0.86)), ("Kammer5_Teiler", (0.95, 0.80, 0.45)), ("Kammer5", (0.50, 0.62, 0.80)), ("Stimmplatte_", (0.80, 0.82, 0.86)), ("Basskammer_Teiler", (0.95, 0.80, 0.45)), ("Basskammer_Trennwand", (0.65, 0.75, 0.88)), ("Basskammer", (0.50, 0.62, 0.80)),
    ("Begleiter_Stimmplatte", (0.80, 0.82, 0.86)), ("Begleiter_Sole", (0.55, 0.68, 0.45)), ("Begleiterstimmstock", (0.62, 0.72, 0.55)),
    ("Lufttaste_Filz", (0.85, 0.15, 0.15)), ("Lufttaste_Gabelgelenk", (0.95, 0.45, 0.10)), ("Lufttaste_Gelenkbolzen", (0.95, 0.85, 0.20)),
    ("Lufttaste_Schubstange", (0.10, 0.60, 0.60)), ("Lufttaste", (0.35, 0.25, 0.15)),
    ("Lufthebel_Feder", (0.75, 0.10, 0.75)), ("Lufthebel_Stuetze_Schraube", (0.40, 0.40, 0.42)), ("Lufthebel_Stuetze", (0.25, 0.75, 0.35)), ("Lufthebel_Bolzen", (0.95, 0.85, 0.20)), ("Lufthebel", (0.70, 0.72, 0.78)),
    ("Knopfzylinder", (0.95, 0.95, 0.92)), ("Knopf", (0.12, 0.12, 0.12)), ("Filzeinlage", (0.85, 0.15, 0.15)),
    ("Achse", (0.55, 0.55, 0.60)), ("Schenkel_Knopf", (0.20, 0.40, 0.85)), ("Schenkel_Klappe", (0.10, 0.65, 0.30)),
    ("Druckstange", (0.10, 0.60, 0.60)), ("Gelenkgabel", (0.95, 0.45, 0.10)), ("Gelenkkopf", (0.60, 0.20, 0.70)),
    ("Gelenkbolzen", (0.95, 0.85, 0.20)), ("Feder", (0.75, 0.10, 0.75)),
    ("Lagerschelle", (0.25, 0.75, 0.35)), ("Lagerschelle_", (0.25, 0.75, 0.35)),
]
def farbe(n):
    if n.startswith("Stimmplatte") and ("_Haken" in n or n.endswith("_Schraube")): return (0.40, 0.40, 0.42)
    if n.startswith("Stimmplatte"): return (0.80, 0.82, 0.86)
    if n.endswith("_Filz"): return (0.85, 0.15, 0.15)              # Klappenfilz rot
    if n.endswith("_Schraube"): return (0.40, 0.40, 0.42)
    if n.endswith("_Gummipfropfen"): return (0.08, 0.08, 0.08)
    if "klappe" in n.lower(): return (0.22, 0.22, 0.25)
    for p_, c_ in FARBEN:
        if n.startswith(p_): return c_
    return (0.60, 0.60, 0.63)
RAHMEN = ("Front", "Rueckteil", "Seitenteil_links", "Seitenteil_rechts")
def spiegel(s, n):                                            # Grundriss/Deckel/Front-Zeichnung sind von aussen gezeichnet -> in X spiegeln
    return s if n in RAHMEN else s.mirror("YZ", (L / 2, 0, 0))
for n, s in teile + mech:
    s = spiegel(s, n)
    a_.add(strecke(s, n), name=n, color=cq.Color(*farbe(n)))
if BASS_DL > 0: print("Strecken um", BASS_DL, "mm - Querschnitt nicht gleichmaessig bei:", _streck_warn)
print(f"Bassteil Standardmechanik 4/4 Helikon: Teilung {T:g}, Gehaeuse {L_END:g} mm (Rohmass {L_END + 2:g}), Teile {len(teile) + len(mech)}")
a_.save(OUT)
ANSICHT = ("Rueckteil", "Seitenteil", "Stimmplatte", "Basskammer", "Front", "Frontblende", "Filzeinlage", "Knopf", "Bassboden", "Basskammer_Leiste_vorne", "Begleiterbrett", "Kammer5", "Begleiterstimmstock", "Begleiter_Sole", "Gaze", "Bassdeckel", "Deckelauflage", "Eckblock", "Fuss",
           "Bassklappe", "Begleiterklappe", "Luftklappe", "Lufttaste", "Lufthebel")
a2 = cq.Assembly(name=f"Bassteil_Standardmechanik_{NAME_V}_Helikon_Mechanik")
for n, s in teile + mech:
    if (n, s) in mech or n.startswith(ANSICHT):
        a2.add(strecke(spiegel(s, n), n), name=n, color=cq.Color(*farbe(n)))
a2.save(OUT.replace(".step", "_Mechanik.step"))
print("Ansicht Mechanik:", OUT.replace(".step", "_Mechanik.step"))
