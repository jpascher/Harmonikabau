# Kugellager-Bassteil (2-reihig, 9 + 8 Knoepfe) - erstes Modell, 23.11.2026
# Aufruf: python3 kugel_bass.py <teilung 18|19|20> <laenge_oder_diskant> <ausgabe.step>
# Quellen (aus C:\ACAD-DWG\0-A-kugellager\Bassteil): kb_front.dxf (Frontbohrungen), kb_boden.dxf (Bassboden mit
# Klappen), kb_deckel.dxf (Bassdeckel), kb_begleiter.dxf (Begleiterstimmstock), kb_trennwand.dxf / kb_hb_trennwand.dxf
# (Basskammern), kb_schelle2.dxf (Lagerschelle), kb_lufthebel.dxf (Luftklappenhebel).
# Achsen: X = Laenge, Y = Tiefe (0 = Rueckteil aussen, B = Front aussen), Z = Hoehe ab Deckelseite (Z 0) Richtung Balg.
# Mechanik (eigener Entwurf nach dem Markus-Prinzip, die DWG-Hauptzeichnung lag nicht vor): je Knopf eine Achse Ø2
# parallel zur Front, hinter den Bassklappen unter dem Boden in Lagerschellen; vom Knopf eine Schubstange nach hinten
# zum Arm der Achse, am anderen Achsende ein Klappenhebel unter die Klappe.
import sys, os, math, ezdxf, ezdxf.path, cadquery as cq
H = os.path.dirname(os.path.abspath(__file__))
T = float(sys.argv[1]); LD = float(sys.argv[2]); OUT = sys.argv[3]
L_END = max(380.0, math.ceil(LD / 10.0 - 1e-9) * 10.0)
D, L, B = 5.5, 382.0, 207.0              # Rohmass wie Markus (2 mm Uebermass auf 380 x 205)
H_RAHMEN = 66.0                          # Front laut Frontbohrungen.dxf
X0B, Y0B = D, D                          # Boden/Deckel-Zeichnung 0..371 x 0..196 -> innen
Z_UB, T_BODEN = 36.0, 6.0                # Boden-Unterseite (Deckelseite) / Dicke
Z_OB = Z_UB + T_BODEN
T_DECKEL = 4.0
T_KL, T_LEDER = 5.5, 1.0                 # Klappe Holz + Leder (liegt unter dem Boden, Richtung Deckel)
Z_KL_U = Z_UB - T_KL - T_LEDER           # Klappen-Unterkante 29.5
Z_ACHSE = Z_UB - 2.0                     # Achsmitte (Lagerschelle: Bohrung 2 mm unter der Auflage)
D_ACHSE, D_STANGE, D_HEBEL = 2.0, 2.0, 2.0
Z_STANGE = 10.0                          # Schubstangen-Hoehe (unter allen Klappen, auch geoeffnet)
H_KAM = 47.3                             # Basskammer-Hoehe (HB-Trennwand)
teile = []
def box(x0, x1, y0, y1, z0, z1):
    return cq.Solid.makeBox(x1 - x0, y1 - y0, z1 - z0, cq.Vector(x0, y0, z0))
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v.Length, p0, v.normalized())
def kugel(p, d):
    return cq.Solid.makeSphere(d / 2 + 0.01, cq.Vector(*p), angleDegrees1=-90, angleDegrees2=90)
def pfad(pts, d):
    s = None
    for a, b in zip(pts[:-1], pts[1:]):
        seg = stab(a, b, d).fuse(kugel(b, d))
        s = seg if s is None else s.fuse(seg)
    return s.fuse(kugel(pts[0], d))
def pts_of(e, tol=0.1):
    p = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(tol)]
    q = [p[0]] + [b_ for a_, b_ in zip(p[:-1], p[1:]) if math.dist(a_, b_) > 1e-6]
    if len(q) > 2 and math.dist(q[0], q[-1]) < 1e-6: q = q[:-1]
    return q
def prisma_xy(pts, z0, z1):
    return cq.Workplane("XY").polyline(pts).close().extrude(z1 - z0).translate((0, 0, z0)).val()
def mod(n): return os.path.join(H, n)

# ---------------- Rahmen ----------------
FR = ezdxf.readfile(mod("kb_front.dxf")).modelspace()
XM = L / 2                                                    # Frontmitte (Zeichnung x 0)
fr = box(0, L, B - D, B, 0, H_RAHMEN)
ru = box(0, L, 0, D, 0, H_RAHMEN)
sl = box(0, D, D, B - D, 0, H_RAHMEN); sr = box(L - D, L, D, B - D, 0, H_RAHMEN)
knoepfe = []                                                  # (X, Z, r_bohrung, reihe)
hoerner = []
for c in FR.query("CIRCLE"):
    x, y, r = c.dxf.center.x, c.dxf.center.y, c.dxf.radius
    if c.dxf.layer == "Front" and r in (5.0, 5.5): knoepfe.append((x, -y, r, 1 if r == 5.5 else 2))
    elif c.dxf.layer == "Front" and r == 9.0: hoerner.append((XM + x, -y))
# Knopfteilung: Reihen mittig, Teilung T (Zeichnung 19 mm)
kn = []
for rr in (1, 2):
    lst = sorted([k for k in knoepfe if k[3] == rr]); n = len(lst)
    for i, k in enumerate(lst): kn.append((XM + (i - (n - 1) / 2) * T, k[1], k[2], rr))
for X_, Z_, r_, rr in kn:
    fr = fr.cut(stab((X_, B + 1, Z_), (X_, B - D - 1, Z_), 2 * r_))
for X_, Z_ in hoerner:
    fr = fr.cut(stab((X_, B + 1, Z_), (X_, B - D - 1, Z_), 18.0))
for e in FR.query('LWPOLYLINE[layer=="EINLAGE-KONTUR_ANGEPASST"]'):     # Einlagen an den Enden, 1 mm tief
    p = [(XM + a, -b) for a, b in pts_of(e)]
    xs_, zs_ = [a for a, _ in p], [b for _, b in p]
    fr = fr.cut(box(min(xs_), max(xs_), B - 1.0, B + 1, max(min(zs_), 0.5), min(max(zs_), H_RAHMEN - 0.5)))
teile += [("Front", fr), ("Rueckteil", ru), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr)]
print("Front: Knoepfe", len(kn), "(Reihe 1:", sum(1 for k in kn if k[3] == 1), "Reihe 2:", sum(1 for k in kn if k[3] == 2), ") Teilung", T, "| Hoerner", len(hoerner))

# ---------------- Boden mit Oeffnungen, Klappen ----------------
BM = ezdxf.readfile(mod("kb_boden.dxf")).modelspace()
def bod(p): return (p[0] + X0B, p[1] - 276.06 + Y0B)
platten, oeffn = [], []
for e in BM.query('LWPOLYLINE[layer=="2Füllung-Öffnug"]'):
    p = [bod(q) for q in pts_of(e)]
    xs = [a for a, _ in p]
    (platten if max(xs) - min(xs) > 300 else oeffn).append(p)
boden = None
for p in platten:
    s = prisma_xy(p, Z_UB, Z_OB)
    boden = s if boden is None else boden.fuse(s)
for p in oeffn: boden = boden.cut(prisma_xy(p, Z_UB - 1, Z_OB + 1))
INNEN = box(D, L - D, D, B - D, -20, 200)
boden = boden.intersect(INNEN)
teile.append(("Bassboden", boden))
klappen = []                                                  # (name, X0, X1, Y0, Y1, art)
for e in BM.query('LWPOLYLINE[layer=="HB-Klappe"]'):
    if not e.closed: continue
    p = [bod(q) for q in pts_of(e)]; xs = [a for a, _ in p]; ys = [b for _, b in p]
    klappen.append([min(xs), max(xs), min(ys), max(ys)])
klappen.sort(key=lambda k: (k[2], k[0]))
bass_kl = sorted([k for k in klappen if k[2] < 115], key=lambda k: k[0])
begl_kl = sorted([k for k in klappen if k[2] >= 115], key=lambda k: k[0])
lo = [p for p in oeffn if min(b for _, b in p) > 160][0]                 # Lufttasten-Oeffnung
lx0, lx1 = min(a for a, _ in lo), max(a for a, _ in lo); ly0, ly1 = min(b for _, b in lo), max(b for _, b in lo)
luft_kl = [lx0 - 3, lx1 + 3, ly0 - 3, ly1 + 3]
def klappe(nm, k):
    x0, x1, y0, y1 = k
    teile.append((nm + "_Leder", box(x0, x1, y0, y1, Z_UB - T_LEDER, Z_UB)))
    teile.append((nm, box(x0 + 0.5, x1 - 0.5, y0 + 0.5, y1 - 0.5, Z_KL_U, Z_UB - T_LEDER)))
for i, k in enumerate(bass_kl): klappe(f"Bassklappe_{i+1}", k)
for i, k in enumerate(begl_kl): klappe(f"Begleiterklappe_{i+1}", k)
klappe("Luftklappe", luft_kl)
print("Klappen: Bass", len(bass_kl), "Begleiter", len(begl_kl), "+ Luftklappe")

# ---------------- Basskammern (Trennwaende + Laengsleisten) ----------------
TW = ezdxf.readfile(mod("kb_trennwand.dxf")).modelspace()
Y_KAM0, Y_KAM1, T_LEIST = 8.0, 114.0, 6.0                    # Rueck-/Frontleiste, dazwischen Trennwaende (98 tief)
lei_h = box(X0B, X0B + 371, Y_KAM0, Y_KAM0 + T_LEIST, Z_OB, Z_OB + H_KAM)
lei_v = box(X0B, X0B + 371, Y_KAM1 - T_LEIST, Y_KAM1, Z_OB, Z_OB + H_KAM)
teile += [("Basskammer_Leiste_hinten", lei_h), ("Basskammer_Leiste_vorne", lei_v)]
for i, e in enumerate(TW.query('LWPOLYLINE[layer=="Layer1"]')):
    xs = [a for a, _ in pts_of(e)]
    tw_ = box(min(xs) + X0B, max(xs) + X0B, Y_KAM0 + T_LEIST - 2, Y_KAM1 - T_LEIST + 2, Z_OB, Z_OB + H_KAM)
    teile.append((f"Basskammer_Trennwand_{i+1}", tw_))
    for j_, (n_, s_) in enumerate(teile):                     # 2-mm-Nuten in den Laengsleisten
        if n_.startswith("Basskammer_Leiste"): teile[j_] = (n_, s_.cut(tw_))

# ---------------- Begleiterstimmstock (8 Toene, Kanzellen beidseitig) ----------------
BG = ezdxf.readfile(mod("kb_begleiter.dxf")).modelspace()
YB = (begl_kl[0][2] + begl_kl[0][3]) / 2                     # Mitte der Begleiteroeffnungen
sole = [e for e in BG.query('LWPOLYLINE[layer=="DEFPOINTS"]')][0]
sp = pts_of(sole); sx0, sx1 = min(a for a, _ in sp), max(a for a, _ in sp); sy0, sy1 = min(b for _, b in sp), max(b for _, b in sp)
# Zeichnung x -> Boden x: Oeffnung -45.48 <-> 298 (Boden), also +343.5
bx = lambda x: x + 343.48 + X0B
H_BEGL = 44.0
stock = box(bx(sx0), bx(sx1), YB + sy0, YB + sy1, Z_OB, Z_OB + H_BEGL)
for e in BG.query('LWPOLYLINE[layer=="KANZELLE"]'):
    p = pts_of(e); xs = [a for a, _ in p]; ys = [b for _, b in p]
    seite = 1 if max(ys) > 1 else -1                            # Kanzellen auf beiden Laengsseiten, 3 mm tief
    y_ = YB + (sy1 if seite > 0 else sy0)
    stock = stock.cut(box(bx(min(xs)), bx(max(xs)), y_ - 3.0 if seite > 0 else y_, y_ if seite > 0 else y_ + 3.0, Z_OB + 2, Z_OB + H_BEGL - 2))
teile.append(("Begleiterstimmstock", stock))

# ---------------- Deckel ----------------
DM = ezdxf.readfile(mod("kb_deckel.dxf")).modelspace()
dx = lambda x: x + 371.0 - 0.75 + X0B - 371.0 + 371.0      # Zeichnung -370.75..-0.75 -> X 5.5..375.5
dx = lambda x: x + 376.25
dy = lambda y: y - 0.62 + Y0B
rect = [e for e in DM.query('LWPOLYLINE[layer=="GHInnenabmessungRechteck"]')][0]
rp = pts_of(rect); dk = box(dx(min(a for a, _ in rp)), dx(max(a for a, _ in rp)), dy(min(b for _, b in rp)), dy(max(b for _, b in rp)), 0, T_DECKEL)
for e in DM.query('LWPOLYLINE[layer=="0"]'):
    p = pts_of(e); xs = [a for a, _ in p]; ys = [b for _, b in p]
    if max(ys) - min(ys) < 20 and max(xs) - min(xs) > 300:
        teile.append(("Bassdeckel_Leiste_hinten", box(dx(min(xs)), dx(max(xs)), dy(min(ys)), dy(max(ys)), 0, T_DECKEL).intersect(INNEN)))
    elif len(p) > 100:                                         # Schalloeffnungen (Ornamente)
        dk = dk.cut(prisma_xy([(dx(a), dy(b)) for a, b in p], -1, T_DECKEL + 1))
slot = [e for e in DM.query('LWPOLYLINE[layer=="Def"]')][0]
slp = [(dx(a), dy(b)) for a, b in pts_of(slot)]
dk = dk.cut(prisma_xy(slp, -1, T_DECKEL + 1))
teile.append(("Bassdeckel", dk.intersect(INNEN)))
LT_X = (min(a for a, _ in slp) + max(a for a, _ in slp)) / 2; LT_Y = (min(b for _, b in slp) + max(b for _, b in slp)) / 2
print("Deckel: Lufttaste X", round(LT_X, 1), "Y", round(LT_Y, 1), "| Luftklappe X", round((luft_kl[0] + luft_kl[1]) / 2, 1), "Y", round((luft_kl[2] + luft_kl[3]) / 2, 1))

# ---------------- Knoepfe ----------------
D_ZYL, D_KOPF, H_KOPF, L_ZYL, SPALT = 7.5, 15.5, 4.5, 20.0, 6.5
for X_, Z_, r_, rr in kn:
    yk = B + SPALT
    teile.append((f"Knopf_R{rr}_{X_:.0f}", stab((X_, yk, Z_), (X_, yk + H_KOPF, Z_), D_KOPF)))
    teile.append((f"Knopfzylinder_R{rr}_{X_:.0f}", stab((X_, yk - L_ZYL, Z_), (X_, yk, Z_), D_ZYL)))
    teile.append((f"Filzring_R{rr}_{X_:.0f}", stab((X_, B, Z_), (X_, B + 2.0, Z_), 2 * r_ + 3).cut(stab((X_, B - 1, Z_), (X_, B + 3, Z_), D_ZYL + 0.4))))
for i, (X_, Z_) in enumerate(hoerner):                         # Hoerner: Rohr Ø18/14 mit Kragen
    h_ = stab((X_, B - D, Z_), (X_, B, Z_), 18.0).cut(stab((X_, B - D - 1, Z_), (X_, B + 1, Z_), 14.0))   # Hornrohr nur in der Front
    h_ = h_.fuse(stab((X_, B, Z_), (X_, B + 1.5, Z_), 33.0).cut(stab((X_, B - 1, Z_), (X_, B + 3, Z_), 14.0)))
    teile.append((f"Horn_{i+1}", h_))

# ---------------- Mechanik ----------------
Y_KN_END = B + SPALT - L_ZYL                                 # Ende Knopfzylinder innen (Y 193.5)
zuord = []                                                   # (Knopf X, Knopf Z, Reihe, Klappe-Name, Klappe)
r1 = sorted([k for k in kn if k[3] == 1]); r2 = sorted([k for k in kn if k[3] == 2])
for k, kl in zip(r2, begl_kl): zuord.append((k, f"Begleiterklappe_{begl_kl.index(kl)+1}", kl))
for k, kl in zip(r1[:len(bass_kl)], bass_kl): zuord.append((k, f"Bassklappe_{bass_kl.index(kl)+1}", kl))
frei = r1[len(bass_kl):]
# Klappenhebel-Ansatz: Bassklappen direkt in Klappenmitte; Begleiterklappen durch die naechste Luecke zwischen den Bassklappen
luecken = [((a[1] + b[0]) / 2, b[0] - a[1]) for a, b in zip(bass_kl[:-1], bass_kl[1:])] + [(bass_kl[0][0] - 6, 10), (bass_kl[-1][1] + 6, 10)]
KX = [k[0] for k in kn]
def frei_von_knoepfen(x): return all(abs(x - k) >= 4.0 for k in KX)
def hebel_x(nm, kl):
    xm = (kl[0] + kl[1]) / 2
    if nm.startswith("Bass") or not any(b[0] - 2 < xm < b[1] + 2 for b in bass_kl):
        lo_, hi_ = kl[0] + 4, kl[1] - 4                       # Hebel irgendwo unter der Klappe, nicht auf Hoehe eines Arms
        if nm.startswith("Begleiter"):
            lo_ = max(lo_, max([b[1] + 2 for b in bass_kl if b[1] < xm] + [lo_])); hi_ = min(hi_, min([b[0] - 2 for b in bass_kl if b[0] > xm] + [hi_]))
    else:
        g_ = min(luecken, key=lambda l: abs(l[0] - xm)); xm = g_[0]; lo_, hi_ = xm - g_[1] / 2 + 1.5, xm + g_[1] / 2 - 1.5
    for d_ in [0] + [s_ * j_ * 0.5 for j_ in range(1, 40) for s_ in (1, -1)]:
        if lo_ <= xm + d_ <= hi_ and frei_von_knoepfen(xm + d_): return xm + d_
    return xm
ach = []
for (k, nm, kl) in zuord:
    xa = hebel_x(nm, kl); x0, x1 = min(k[0], xa) - 8.0, max(k[0], xa) + 8.0
    ach.append(dict(k=k, nm=nm, kl=kl, xh=xa, x0=x0, x1=x1))
# Achsen auf Linien verteilen (gleiche Linie nur ohne Ueberlappung, 6 mm Abstand), hinter den Bassklappen
Y_LIN0, DY_LIN = min(b[2] for b in bass_kl) - 8.0, 6.0
linien = []
for a in sorted(ach, key=lambda a: (a["x1"] - a["x0"]), reverse=True):
    for li, lst in enumerate(linien):
        if all(a["x0"] > b["x1"] + 6 or a["x1"] < b["x0"] - 6 for b in lst): lst.append(a); a["li"] = li; break
    else: linien.append([a]); a["li"] = len(linien) - 1
for a in ach: a["y"] = Y_LIN0 - a["li"] * DY_LIN
print("Mechanik: Achsen", len(ach), "auf", len(linien), "Linien, Y", round(Y_LIN0, 1), "...", round(Y_LIN0 - (len(linien) - 1) * DY_LIN, 1))
mech = []
SCH_B, SCH_T, SCH_H = 5.0, 4.0, 5.3                           # Lagerschelle je Achsende (Y-Breite, X-Dicke, Hoehe)
for i, a in enumerate(ach):
    y, (X_, Z_, r_, rr) = a["y"], a["k"]
    mech.append((f"Achse_{i+1}", stab((a["x0"], y, Z_ACHSE), (a["x1"], y, Z_ACHSE), D_ACHSE)))
    for xs_ in (a["x0"] + SCH_T / 2 + 0.5, a["x1"] - SCH_T / 2 - 0.5):
        s_ = box(xs_ - SCH_T / 2, xs_ + SCH_T / 2, y - SCH_B / 2, y + SCH_B / 2, Z_UB - SCH_H, Z_UB).cut(stab((xs_ - 3, y, Z_ACHSE), (xs_ + 3, y, Z_ACHSE), D_ACHSE + 0.1))
        mech.append((f"Lagerschelle_{i+1}_{'a' if xs_ < a['x1'] - 3 else 'b'}", s_))
    # Arm an der Knopfstelle: von der Achse senkrecht nach unten auf Stangenhoehe, Gabelgelenk
    mech.append((f"Arm_{i+1}", box(X_ - 0.75, X_ + 0.75, y - 1.5, y + 1.5, Z_STANGE - 1.75, Z_ACHSE + 1.0)))
    gk = box(X_ - 1.85, X_ + 1.85, y + 1.5, y + 5.0, Z_STANGE - 2.0, Z_STANGE + 2.0).fuse(stab((X_, y + 5.0, Z_STANGE), (X_, y + 11.0, Z_STANGE), 4.0))
    mech.append((f"Gabelkopf_{i+1}", gk))
    mech.append((f"Gelenkbolzen_{i+1}", stab((X_ - 2.3, y, Z_STANGE), (X_ + 2.3, y, Z_STANGE), 1.4)))
    # Schubstange: vom Gabelkopf nach vorne auf Stangenhoehe, zuletzt 45 Grad auf Knopfhoehe
    y_r = Y_KN_END - abs(Z_ - Z_STANGE)
    pp = [(X_, y + 11.0, Z_STANGE), (X_, y_r, Z_STANGE), (X_, Y_KN_END, Z_)]
    mech.append((f"Schubstange_R{rr}_{X_:.0f}", pfad([q for j_, q in enumerate(pp) if j_ == 0 or math.dist(q, pp[j_ - 1]) > 0.05], D_STANGE)))
    # Klappenhebel: von der Achse nach unten knapp unter die Klappe, dann nach vorne bis Klappenmitte (Begleiter mit X-Versatz)
    kl, xh = a["kl"], a["xh"]; xm, ym = (kl[0] + kl[1]) / 2, (kl[2] + kl[3]) / 2
    z_kh = Z_KL_U - 1.2
    if abs(xh - xm) < 0.1:
        hp = [(xh, y, Z_ACHSE), (xh, y, z_kh), (xm, ym, z_kh)]
    else:
        y_q = (bass_kl[0][3] + begl_kl[0][2]) / 2                # zwischen Bass- und Begleiterreihe in X versetzen
        hp = [(xh, y, Z_ACHSE), (xh, y, z_kh), (xh, y_q, z_kh), (xm, y_q, z_kh), (xm, ym, z_kh)]
    mech.append((f"Klappenhebel_{i+1}", pfad(hp, D_HEBEL)))
    print(f"   Achse {i+1}: Knopf R{rr} X {X_:.1f} -> {a['nm']} (Hebel X {xh:.1f}), Achse X {a['x0']:.1f}-{a['x1']:.1f}, Y {y:.1f}")
for (X_, Z_, r_, rr) in frei: print(f"   Knopf R{rr} X {X_:.1f}: noch ohne Zuordnung (9. Knopf der Reihe 1)")
# Lufttaste: Taste im Deckelschlitz, Hebel (Ø2.5) von der Taste zur Luftklappe, Lagerbock dazwischen unter dem Boden
lkx, lky = (luft_kl[0] + luft_kl[1]) / 2, (luft_kl[2] + luft_kl[3]) / 2
xp = lkx + (LT_X - lkx) * (114.0 - 55.0) / (114.0 - 12.65)     # Drehpunkt wie im Lufthebel-Plan (Verhaeltnis)
# Taste: Kappe im Deckelschlitz, schmaler Stoessel (zwischen zwei Schubstangen) bis ueber die Stangenebene;
# Hebel liegt darueber unter dem Boden, Drehpunkt zwischen Taste und Luftklappe
xs_frei = sorted(KX); xst = min([(xs_frei[j_] + xs_frei[j_ + 1]) / 2 for j_ in range(len(xs_frei) - 1)], key=lambda x_: abs(x_ - LT_X))
mech.append(("Lufttaste", box(LT_X - 30, LT_X + 30, LT_Y - 4.5, LT_Y + 4.5, -6.0, 3.0).fuse(box(xst - 2.5, xst + 2.5, LT_Y - 2, LT_Y + 2, 2.9, 27.0))))
xp = lkx + (xst - lkx) * (114.0 - 55.0) / (114.0 - 12.65)
mech.append(("Lufthebel", pfad([(xst, LT_Y, 28.3), (xp, LT_Y, Z_ACHSE), (lkx, LT_Y, Z_KL_U - 1.25)], 2.5).cut(stab((xp, LT_Y - 5, Z_ACHSE), (xp, LT_Y + 5, Z_ACHSE), 1.6))))
_lb = box(xp - 4, xp + 4, LT_Y - 4.0, LT_Y - 1.6, Z_UB - 5.0, Z_UB).fuse(box(xp - 4, xp + 4, LT_Y + 1.6, LT_Y + 4.0, Z_UB - 5.0, Z_UB))
mech.append(("Lufthebel_Lagerbock", _lb.cut(stab((xp, LT_Y - 5, Z_ACHSE), (xp, LT_Y + 5, Z_ACHSE), 1.6))))
mech.append(("Lufthebel_Bolzen", stab((xp, LT_Y - 4.3, Z_ACHSE), (xp, LT_Y + 4.3, Z_ACHSE), 1.5)))

# ---------------- Laenge (Strecken wie beim Markus-Bass) und Zusammenbau ----------------
BASS_DL = L_END - 380.0
exec(open(os.path.join(H, "streck.py"), encoding="utf-8").read().replace('FRONT_EBENE = float(_os.environ.get("BASS_FRONT_EBENE", "50.0"))', 'FRONT_EBENE = 23.0'))
a_ = cq.Assembly(name="Kugellager_Bassteil")
farbe = lambda n: (0.78, 0.60, 0.38) if n.startswith(("Front", "Rueck", "Seiten")) else (0.88, 0.80, 0.62) if n.startswith("Bassboden") else \
    (0.30, 0.22, 0.15) if n.startswith("Bassdeckel") else (0.45, 0.30, 0.20) if n.endswith("Leder") else (0.20, 0.20, 0.22) if "klappe" in n.lower() else \
    (0.65, 0.75, 0.88) if n.startswith("Basskammer") else (0.62, 0.72, 0.55) if n.startswith("Begleiter") else (0.15, 0.15, 0.15) if n.startswith("Knopf_") else \
    (0.95, 0.95, 0.92) if n.startswith("Knopfzyl") else (0.55, 0.55, 0.20) if n.startswith("Filz") else (0.70, 0.55, 0.30) if n.startswith("Horn") else \
    (0.85, 0.85, 0.80) if n.startswith(("Lagerschelle", "Lufthebel_Lager")) else (0.20, 0.45, 0.75) if n.startswith("Gabel") else (0.60, 0.60, 0.63)
for n, s in teile + mech:
    a_.add(strecke(s, n), name=n, color=cq.Color(*farbe(n)))
if BASS_DL > 0: print("Strecken um", BASS_DL, "mm - Querschnitt nicht gleichmaessig bei:", _streck_warn)
print(f"Kugellager-Bassteil: Teilung {T:g}, Gehaeuse {L_END:g} mm (Rohmass {L_END + 2:g}), Teile {len(teile) + len(mech)}")
a_.save(OUT)
