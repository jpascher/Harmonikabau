# ---------- Bassmechanik Marie 38 (wird von kammern38.py per exec eingebunden, vor der Spiegelung) ----------
# Bauteile wie Markus/V3: Achsen Ø2,4 unter dem Boden in Lagerboecken (lagerbeide.dxf), Schenkel Ø2,4 von der Achse
# zum Gelenk (Fusion-Gelenkkopf + Gabel + Bolzen), Schubstange Ø2,4 vom Knopfzylinder zum Gelenkkopf, Klappenhebel
# Ø2,4 von der Achse unter die Klappenmitte, je Achse eine Schenkelfeder (3 Windungen). Anordnung nach bassteil-alles.dxf
# (Linien 0HB_FEDER / 0BG_QUER, Rest nach Plan). Koordinaten vor der Spiegelung, Zeichnung +1.
Z_UB = Z_BU
H_LB, Z_ACHSE = 5.0, Z_BU - 2.5
D_ACHSE = D_STANGE = 2.4
GK_H, GK_L, GK_R, L_GABEL = 3.19, 13.0, 2.0, 8.0
Z_GG = T_DECKEL_BASS + 0.3 + GK_R + 1.5
Z_KL_U = min(s.BoundingBox().zmin for n_, s in teile if n_.startswith(("Klappe_", "Begleiterklappe_")))
Z_KH = Z_KL_U - D_STANGE / 2 - 0.05
HUB_K = 6.0
mech = []
R_Z = {"A": Z_RH[1], "B": Z_RH[2]}
X_R = {"A": [round(X_KM + (i_ - 4) * T_KNOPF, 2) for i_ in range(9)], "B": [round(X_KM + (i_ - 4.5) * T_KNOPF, 2) for i_ in range(10)]}
def kn(r_, x_): return min(X_R[r_], key=lambda q_: abs(q_ - x_))
KMY = 135.5; KUY = 76.3                                   # Klappenmitten Kanalreihe / untere Reihe
# (Name, Achs-Y, Arme [(Reihe, Knopf-X)], Hebel [(X am Achsansatz, Klappenmitte)])
ACHSEN = [
    ("H1", 52.0, [("B", 124.5)], [(83.0, (83.0, KUY))]),
    ("H2", 52.0, [("B", 162.5)], [(139.5, (137.0, KUY))]),
    ("H3", 52.0, [("B", 200.5)], [(191.0, (191.0, KUY))]),
    ("H4", 49.3, [("B", 238.5)], [(245.0, (245.0, KUY))]),
    ("H5", 49.3, [("B", 276.5)], [(299.0, (299.0, KUY))]),
    ("K1", 111.0, [("A", 115.0)], [(83.0, (83.0, KMY))]),
    ("K2", 111.0, [("A", 153.0)], [(137.0, (137.0, KMY))]),
    ("K3", 112.3, [("A", 191.0)], [(196.0, (196.0, KMY))]),
    ("K4", 100.3, [("B", 105.5), ("A", 229.0)], [(245.0, (245.0, KMY))]),
    ("K5", 108.3, [("A", 267.0)], [(299.0, (299.0, KMY))]),
    ("BL1", 37.3, [("B", 181.5)], [(56.0, (32.6, 50.1))]),
    ("BL2", 43.3, [("A", 172.0)], [(60.0, (32.6, 84.1))]),
    ("BL3", 96.3, [("B", 143.5)], [(58.0, (32.6, 118.1))]),
    ("BL4", 115.3, [("A", 134.0)], [(58.0, (32.6, 152.1))]),
    ("BR1", 37.3, [("B", 219.5)], [(326.0, (350.4, 50.1))]),
    ("BR2", 43.3, [("A", 210.0)], [(322.0, (350.4, 84.1))]),
    ("BR3", 96.3, [("B", 257.5)], [(324.0, (350.4, 118.1))]),
    ("BR4", 115.3, [("A", 248.0)], [(324.0, (350.4, 152.1))]),
]
ACHSEN = [(n_, y_, [(r_, kn(r_, x_)) for r_, x_ in a_], h_) for n_, y_, a_, h_ in ACHSEN]
# Hindernisse fuer Lagerboecke (Grundriss unter dem Boden)
rechtecke = [(x0 - RAND_V + 0.2, x1 + RAND_V - 0.2, y0 - RAND_V + 0.2, y1 + RAND_V - 0.2) for (x0, x1, y0, y1) in oeff]
kreise = [((s_.BoundingBox().xmin + s_.BoundingBox().xmax) / 2, (s_.BoundingBox().ymin + s_.BoundingBox().ymax) / 2, s_.BoundingBox().xlen / 2 + 0.3)
          for n_, s_ in teile if n_.startswith("Begleiterklappe_") and n_.endswith("Leder")]
rechtecke += [(b_.xmin - 0.3, b_.xmax + 0.3, b_.ymin - 0.3, b_.ymax + 0.3) for b_ in [s_.BoundingBox() for n_, s_ in teile if n_.startswith(("Eckblock_", "Deckelauflage_hinten"))]]
HL, HB = 4.3, 2.6                                         # halbe Bocklaenge / -breite
def _seg_d(p, a, b):
    ax, ay = a; bx, by = b; dx, dy = bx - ax, by - ay; l2 = dx * dx + dy * dy
    t = max(0, min(1, ((p[0] - ax) * dx + (p[1] - ay) * dy) / l2)) if l2 > 0 else 0
    return math.dist(p, (ax + t * dx, ay + t * dy))
def _rect_seg(r, a, b):                                   # Abstand Rechteck (x0,x1,y0,y1) zu Strecke a-b (grob ueber Stuetzpunkte)
    x0, x1, y0, y1 = r; best = 99
    for k in range(21):
        p = (a[0] + (b[0] - a[0]) * k / 20, a[1] + (b[1] - a[1]) * k / 20)
        dx = max(x0 - p[0], 0, p[0] - x1); dy = max(y0 - p[1], 0, p[1] - y1); best = min(best, math.hypot(dx, dy))
    return best
# Arme und Hebel als feste Hindernisse
arme_alle = [(n_, x_, y_) for n_, y_, a_, h_ in ACHSEN for r_, x_ in a_]
hebel_pa = {n_: [list(h) for h in h_] for n_, y_, a_, h_ in ACHSEN}
boecke = []                                               # (Achse, x, y)
def bock_frei(nm, x, y, eigene):
    r = (x - HL, x + HL, y - HB, y + HB)
    if y - HB < 11.8 or y + HB > 153.0 - 0.3: return False
    for q in rechtecke:
        if not (r[1] < q[0] or q[1] < r[0] or r[3] < q[2] or q[3] < r[2]): return False
    for cx, cy, rr in kreise:
        if math.hypot(max(r[0] - cx, 0, cx - r[1]), max(r[2] - cy, 0, cy - r[3])) < rr: return False
    for (n2, bx, by) in boecke:
        if abs(bx - x) < 2 * HL + 0.3 and abs(by - y) < 2 * HB + 0.3: return False
    for n2, y2, a2, h2 in ACHSEN:                          # fremde Achskoerper
        if n2 == nm: continue
        xs2 = achs_x.get(n2)
        if xs2 and abs(y2 - y) < HB + D_ACHSE / 2 + 0.15 and r[0] < xs2[1] and xs2[0] < r[1]: return False
    for (n2, ax, ay) in arme_alle:
        if math.hypot(max(r[0] - ax, 0, ax - r[1]), max(r[2] - ay, 0, ay - r[3])) < D_STANGE / 2 + (1.8 if n2 != nm else 0.5):
            return False
    for n2, y2, a2, h2 in ACHSEN:                          # Klappenhebel (Grundriss), eigene tauchen unter dem Bock durch
        if n2 == nm: continue
        for (pax, pk) in hebel_pa[n2]:
            if _rect_seg(r, (pax, y2), pk) < D_STANGE / 2 + 0.3: return False
    for t_ in eigene:
        if abs(t_ - x) < HL + D_STANGE / 2 + 0.5: return False
    return True
# Lagerboecke nach Regeln: Ende hinter dem aeussersten Arm/Hebel (6 mm), bei schraeg nach aussen laufendem Begleiterhebel
# innen neben dem Hebelansatz; lange Achsen bekommen einen Mittelbock; stossen zwei Achsen in einer Linie zusammen,
# teilen sie sich einen Bock (Achsenden je 1 mm vor der Schraube). Feder neben dem Arm auf der freien Seite.
A_FREI = HL + D_STANGE / 2 + 0.6                          # Bock neben Hebel
A_ARM = A_FREI + 3 * 0.6 + D_STANGE / 2 + 0.6               # Bock neben Arm (Platz fuer die Feder)
FEAT = {}
for n_, y_, a_, h_ in ACHSEN:
    FEAT[n_] = sorted([(x_, "A") for r_, x_ in a_] + [(h[0], "H" if abs(h[1][0] - h[0]) < 10 else "S") for h in hebel_pa[n_]])
def _abst(kind): return A_ARM if kind == "A" else A_FREI
def _erlaubt(n_, x_, art):
    f_ = FEAT[n_]
    if any(abs(x_ - fx) < _abst(fk) for fx, fk in f_): return False
    if art == "lo": return x_ < f_[0][0]
    if art == "hi": return x_ > f_[-1][0]
    return f_[0][0] < x_ < f_[-1][0]
bock_soll = {}
for n_, y_, a_, h_ in ACHSEN:
    lo_, hi_ = FEAT[n_][0], FEAT[n_][-1]
    bl_ = lo_[0] - _abst(lo_[1])                          # Endboecke immer aussen; Schraeghebel taucht unter dem Bock durch
    bh_ = hi_[0] + _abst(hi_[1])
    bock_soll[n_] = [[bl_, "lo"], [bh_, "hi"]]
    if False:                                             # kein Mittelbock (Schraube wuerde die Achse kreuzen)
        frei_ = [x_ for x_ in [bl_ + 0.5 * k for k in range(int((bh_ - bl_) * 2))] if _erlaubt(n_, x_, "mid")]
        if frei_: bock_soll[n_].append([min(frei_, key=lambda q: abs(q - (bl_ + bh_) / 2)), "mid"])
def _achs_x_neu():
    for n_, y_, a_, h_ in ACHSEN:
        xs_ = [f[0] for f in FEAT[n_]] + [b[0] for b in bock_soll[n_]]
        achs_x[n_] = (min(xs_), max(xs_))
achs_x = {}; _achs_x_neu()
# gemeinsame Boecke kollinearer Achsen
geteilt = set()
namen = [q[0] for q in ACHSEN]; ay = {q[0]: q[1] for q in ACHSEN}
for i_, n1 in enumerate(namen):
    for n2 in namen[i_ + 1:]:
        if abs(ay[n1] - ay[n2]) > 0.05: continue
        for e1 in bock_soll[n1]:
            for e2 in bock_soll[n2]:
                if abs(e1[0] - e2[0]) < 2 * HL + 0.3:
                    m_ = (e1[0] + e2[0]) / 2; e1[0] = m_; e2[0] = m_; e1[1] = e2[1] = "geteilt"; geteilt.add((n2, round(m_, 2)))
                    print("Gemeinsamer Bock", n1, n2, "bei X", round(m_, 1))
_achs_x_neu()
# Konflikte loesen: Bock entlang der Achse verschieben (Enden nach aussen bzw. beim Schraeghebel nach innen)
boecke = []
for n_, y_, a_, h_ in sorted(ACHSEN, key=lambda q: (q[0][0] != "H", q[0][0] != "K")):
    for e in bock_soll[n_]:
        if (n_, round(e[0], 2)) in geteilt: continue
        x0_ = e[0]; gef = None
        for k in range(0, 61):
            for sg_ in ((1, -1) if e[1] == "mid" else (-1, 1)):
                x_ = x0_ + sg_ * 0.5 * k
                if e[1] != "geteilt" and not _erlaubt(n_, x_, e[1]): continue
                if e[1] == "geteilt" and k > 0: continue
                if bock_frei(n_, x_, y_, []): gef = x_; break
            if gef is not None: break
        if gef is None: print("HINWEIS Bock nicht frei", n_, e[1], round(x0_, 1), round(y_, 1)); gef = x0_
        elif abs(gef - x0_) > 0.01: print("Bock verschoben", n_, e[1], round(x0_, 1), "->", round(gef, 1))
        e[0] = gef; boecke.append((n_, gef, y_)); _achs_x_neu()
print("Mechanik: Achsen", len(ACHSEN), "Lagerboecke", len(boecke))
# ---- Lagerbock-Koerper (lagerbeide.dxf) ----
_LBD = ezdxf.readfile("lagerbeide.dxf").modelspace()
def _prof(lay):
    e_ = [e for e in _LBD.query(f'LWPOLYLINE[layer=="{lay}"]') if e.closed][0]
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e_).flattening(0.02)]
    return pp_[:-1] if math.dist(pp_[0], pp_[-1]) < 1e-6 else pp_
PROF_E = _prof("Layer1")
def bock_solid(g_=1):
    yz_ = [(-g_ * (-q[1] - 2.6), Z_UB - (q[0] + 13.6)) for q in PROF_E]
    lb_ = cq.Workplane("YZ").polyline(yz_).close().extrude(8.6).val().translate(cq.Vector(-4.3, 0, 0))
    lb_ = lb_.cut(stab((-5, 0, Z_ACHSE), (5, 0, Z_ACHSE), D_ACHSE + 0.1)).cut(stab((0, 0, Z_UB - H_LB - 1), (0, 0, Z_UB + 1), 1.5))
    sc_ = stab((0, 0, Z_UB - H_LB - 0.9), (0, 0, Z_UB - H_LB + 10), 1.4).fuse(stab((0, 0, Z_UB - H_LB - 0.9), (0, 0, Z_UB - H_LB), 2.5))
    return lb_, sc_
for k_, (n_, bx_, by_) in enumerate(boecke):
    nb_ = [y2 for n2, y2, a2, h2 in ACHSEN if n2 != n_ and achs_x[n2][0] - 4 < bx_ < achs_x[n2][1] + 4 and abs(y2 - by_) < 8]
    g_ = 1 if (nb_ and min(nb_, key=lambda q: abs(q - by_)) > by_) else -1
    lb, sc = bock_solid(g_)
    mech.append((f"Lagerbock_{n_}_{k_+1}", lb.translate(cq.Vector(bx_, by_, 0))))
    mech.append((f"Lagerbock_Schraube_{n_}_{k_+1}", sc.translate(cq.Vector(bx_, by_, 0))))
    boden = boden.cut(stab((bx_, by_, Z_UB - 1), (bx_, by_, Z_UB + 4.0), 1.1))      # Vorbohrung im Boden
for n_, y_, a_, h_ in ACHSEN:
    x0_, x1_ = achs_x[n_]
    mech.append((f"Achse_{n_}", stab((x0_ + 1.0, y_, Z_ACHSE), (x1_ - 1.0, y_, Z_ACHSE), D_ACHSE)))   # Enden 1 mm vor der Schraube
# ---- Gelenk (wie Markus/V3) ----
_gk = cq.importers.importStep("kb_gelenkkopf.stp").solids().vals()
GK = _gk[0].fuse(_gk[1]).clean()
GK = GK.cut(stab((0, 0.5, 0), (0, -GK_L - 1, 0), D_STANGE)).cut(stab((-3, GK_H, 0), (3, GK_H, 0), 1.5))
def gelenk(P, s_, draht_s, draht_st, nm):
    x, y, z = P
    mech.append((f"Gelenkkopf_{nm}", GK.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180).translate(cq.Vector(x, y + GK_H, z)).cut(draht_st)))
    g = stab((x - GK_R, y, z), (x + GK_R, y, z), 2 * GK_R).fuse(stab((x, y, z), (x, y, z - s_ * L_GABEL), 2 * GK_R))
    g = g.cut(box(x - 1.05, x + 1.05, y - 2.4, y + 6, *sorted((z - s_ * 2.4, z + s_ * 2.4))))
    g = g.cut(stab((x, y, z - s_ * GK_H), (x, y, z - s_ * (L_GABEL + 1)), D_STANGE)).cut(stab((x - 3, y, z), (x + 3, y, z), 1.5)).cut(draht_s)
    mech.append((f"Gelenkgabel_{nm}", g))
    mech.append((f"Gelenkbolzen_{nm}", stab((x - 2.4, y, z), (x + 2.4, y, z), 1.4)))
_lt = [s_ for n_, s_ in teile if n_ == "Lufttaste"][0].BoundingBox()     # gedrueckt (eben mit dem Deckel)
def pfad_stab(pfad):
    st = None
    for p_a, p_b in zip(pfad[:-1], pfad[1:]):
        seg = stab(p_a, p_b, D_STANGE).fuse(cq.Solid.makeSphere(D_STANGE / 2 + 0.01, cq.Vector(*p_b)))   # Kugel etwas groesser: keine tangentialen Flaechen
        st = seg if st is None else st.fuse(seg)
    st = st.clean()
    if len(st.Solids()) > 1: print("WARNUNG mehrteiliger Stab", [[round(v, 1) for v in p] for p in pfad])
    return st
arm_info = []
for n_, y_, a_, h_ in ACHSEN:
    for r_, Xk in a_:
        Zr = R_Z[r_]; zyl_ = [s_ for nn_, s_ in teile if nn_.startswith(f"Knopfzylinder_R{1 if r_ == 'A' else 2}_{Xk:.0f}")][0]
        y_end = zyl_.BoundingBox().ymin + 0.6
        gx, gy = Xk, y_
        sk_ = stab((gx, gy, Z_ACHSE), (gx, gy, Z_GG + 4.5), D_STANGE).cut(stab((gx - 3, gy, Z_ACHSE), (gx + 3, gy, Z_ACHSE), D_ACHSE))
        mech.append((f"Schenkel_Knopf_{n_}_{r_}{Xk:.0f}", sk_))
        ueber = _lt.xmin - 2.5 < Xk < _lt.xmax + 2.5
        z_mid = max(Zr, _lt.zmax + D_STANGE / 2 + 1.0) if ueber else Zr      # ueber die gedrueckte Lufttaste
        z_mid = min(z_mid, Zr + D_ZYL / 2 - D_STANGE / 2 - 0.3) if ueber else z_mid
        y_a = (min(_lt.ymin - 1.5, y_end) if ueber else y_end) - (z_mid - Z_GG)
        pfad = [(gx, gy + 4.6, Z_GG), (gx, y_a, Z_GG), (gx, y_a + (z_mid - Z_GG), z_mid), (Xk, y_end, z_mid)]
        pfad = [p for k_, p in enumerate(pfad) if k_ == 0 or math.dist(p, pfad[k_ - 1]) > 0.05]
        st = pfad_stab(pfad)
        nm_ = f"{n_}_{r_}{Xk:.0f}"
        mech.append((f"Schubstange_{nm_}", st))
        gelenk((gx, gy, Z_GG), -1, sk_, st, nm_)
        arm_info.append((n_, nm_, gx, gy, r_, Xk))
# ---- Klappenhebel: von der Achse senkrecht hinunter, waagrecht unter die Klappenmitte; unter fremden Boecken tiefer ----
for n_, y_, a_, h_ in ACHSEN:
    for j_, (pax, pk) in enumerate(hebel_pa[n_]):
        pa = (pax, y_); L_ = math.dist(pa, pk); ux_, uy_ = (pk[0] - pa[0]) / L_, (pk[1] - pa[1]) / L_
        unter_ = []
        for (n2, bx_, by_) in boecke:
            t_ = (bx_ - pa[0]) * ux_ + (by_ - pa[1]) * uy_
            if _rect_seg((bx_ - HL, bx_ + HL, by_ - HB, by_ + HB), pa, pk) < D_STANGE / 2 + 0.3: unter_.append(t_)
        if not unter_:
            pts_ = [(pa[0], pa[1], Z_ACHSE), (pa[0], pa[1], Z_KH), (pk[0], pk[1], Z_KH)]
        else:
            z_lo = Z_UB - H_LB - 0.9 - D_STANGE / 2 - 0.15; t_r = min(L_ - 1.0, max(unter_) + HL + 6)
            pts_ = [(pa[0], pa[1], Z_ACHSE), (pa[0], pa[1], z_lo), (pa[0] + ux_ * t_r, pa[1] + uy_ * t_r, z_lo),
                    (pa[0] + ux_ * min(L_, t_r + (Z_KH - z_lo)), pa[1] + uy_ * min(L_, t_r + (Z_KH - z_lo)), Z_KH), (pk[0], pk[1], Z_KH)]
        pts_ = [q_ for k_, q_ in enumerate(pts_) if k_ == 0 or math.dist(q_, pts_[k_ - 1]) > 0.05]
        mech.append((f"Klappenhebel_{n_}", pfad_stab(pts_)))
# ---- Rueckdruckfedern (wie Markus): 3 Windungen neben dem Arm zur naechsten Buchse ohne Hebel dazwischen ----
D_DRAHT, R_WIND, N_WIND, P_WIND = 0.5, D_ACHSE / 2 + 0.2 + 0.25, 3.5, 0.6   # 3,5 Windungen: Anfang unten, Ende oben
def draht(punkte):                                        # Draht als Sweep mit runden Ecken (robuster als Zylinder+Kugeln)
    pts_ = [cq.Vector(*p) for p in punkte]
    prof_ = cq.Wire.makeCircle(D_DRAHT / 2, pts_[0], (pts_[1] - pts_[0]).normalized())
    return cq.Solid.sweep(prof_, [], cq.Wire.makePolygon(pts_), transitionMode="round")
feder_ersatz = {}
for n_, y_, a_, h_ in ACHSEN:
    # Feder am offenen Ende: zwischen einem Arm und dem naechsten Lagerbock, ohne Hebel oder anderen Arm dazwischen
    bx_ = [b[1] for b in boecke if abs(b[2] - y_) < 0.05 and achs_x[n_][0] - 0.5 <= b[1] <= achs_x[n_][1] + 0.5]
    hx_ = [h[0] for h in hebel_pa[n_]]; ax_ = [x_ for r_, x_ in a_]
    kand = []
    for gx_ in ax_:
        for sd in (1.0, -1.0):
            nb_ = [b for b in bx_ if (b - gx_) * sd > 0]
            if not nb_: continue
            bn_ = min(nb_, key=lambda b: abs(b - gx_))
            if any(min(gx_, bn_) < h < max(gx_, bn_) for h in hx_ + [x for x in ax_ if x != gx_]): continue
            kand.append((-abs(bn_ - gx_), gx_, sd))
    if kand:
        _, gx, s_dir = max(kand)
    else:
        gx, s_dir = ax_[0], 1.0; print("WARNUNG Feder ohne offenes Ende", n_)
    t0 = gx + s_dir * (D_STANGE / 2 + 0.4); L_SP = N_WIND * P_WIND
    hel = cq.Wire.makeHelix(P_WIND, L_SP, R_WIND)
    sp = cq.Workplane("XZ").center(R_WIND, 0).circle(D_DRAHT / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90)
    if s_dir < 0: sp = sp.mirror("YZ")                  # Gegenseite: Feder spiegeln (Links- statt Rechtswindung)
    sp = sp.translate(cq.Vector(t0, y_, Z_ACHSE))
    nx_, ny_ = 0.0, -1.0
    z_h = Z_GG + 4.5 + 2.5; ab = GK_R + 0.1 + D_DRAHT / 2
    s1 = draht([(t0, y_, Z_ACHSE - R_WIND), (t0, y_ - ab, Z_ACHSE - R_WIND - 1.0), (t0, y_ - ab, z_h),
                (gx - s_dir * ab, y_ - ab, z_h), (gx - s_dir * ab, y_ + 0.8, z_h)])
    tE = t0 + s_dir * L_SP
    s2 = draht([(tE, y_, Z_ACHSE + R_WIND), (tE, y_ + 0.8, Z_UB - 0.3)])
    gel = [cq.Solid.makeSphere(D_DRAHT / 2 + 0.06, cq.Vector(t0, y_, Z_ACHSE - R_WIND)), cq.Solid.makeSphere(D_DRAHT / 2 + 0.06, cq.Vector(tE, y_, Z_ACHSE + R_WIND))]
    fd_ = sp.fuse(gel[0], gel[1], s1, s2).clean()
    if len(fd_.Solids()) > 1: print("WARNUNG Feder mehrteilig", n_)
    mech.append((f"Feder_Achse_{n_}", fd_))
    # Ersatzkoerper nur fuer die Pruefung (Boolesche Operationen mit der Wendel haengen gelegentlich)
    xa_, xb_ = sorted((t0, tE))
    hk_ = stab((xa_ - 0.25, y_, Z_ACHSE), (xb_ + 0.25, y_, Z_ACHSE), 2 * (R_WIND + D_DRAHT / 2)).cut(stab((xa_ - 1, y_, Z_ACHSE), (xb_ + 1, y_, Z_ACHSE), 2 * (R_WIND - D_DRAHT / 2)))
    pp1_ = [(t0, y_, Z_ACHSE - R_WIND), (t0, y_ - ab, Z_ACHSE - R_WIND - 1.0), (t0, y_ - ab, z_h), (gx - s_dir * ab, y_ - ab, z_h), (gx - s_dir * ab, y_ + 0.8, z_h)]
    for q0_, q1_ in list(zip(pp1_[:-1], pp1_[1:])) + [((tE, y_, Z_ACHSE + R_WIND), (tE, y_ + 0.8, Z_UB - 0.3))]: hk_ = hk_.fuse(stab(q0_, q1_, D_DRAHT))
    feder_ersatz[f"Feder_Achse_{n_}"] = hk_
print("Mechanik Teile", len(mech), "| Achsen", {n_: [round(v, 1) for v in achs_x[n_]] for n_, *_ in ACHSEN})
print("Hebelansaetze", {n_: [round(h[0], 1) for h in hebel_pa[n_]] for n_, *_ in ACHSEN})
# ---- Pruefung: Mechanik gegen Gehaeuse in Ruhe, dann jede Taste gedrueckt (Lufttaste gedrueckt) ----
def _bbt(b1, b2): return not (b1.xmax < b2.xmin or b2.xmax < b1.xmin or b1.ymax < b2.ymin or b2.ymax < b1.ymin or b1.zmax < b2.zmin or b2.zmax < b1.zmin)
gehaeuse = [("Bassboden_neu", boden)] + [(n_, s_) for n_, s_ in teile if n_ != "Bassboden"]
mech_pr = [(n_, feder_ersatz.get(n_, s_)) for n_, s_ in mech]
kol = {}
for n1, s1 in mech_pr:
    if "Schraube" in n1: continue
    b1 = s1.BoundingBox()
    for n2, s2 in gehaeuse + [m_ for m_ in mech_pr if m_[0] != n1]:
        if not _bbt(b1, s2.BoundingBox()): continue
        if (n1.startswith("Achse_") and n2.startswith("Lagerbock")) or (n2.startswith("Achse_") and n1.startswith("Lagerbock")): continue
        if {n1.split("_")[0], n2.split("_")[0]} == {"Achse", "Klappenhebel"} and n1.split("_")[-1] == n2.split("_")[-1]: continue
        if n1.startswith("Schubstange_") and n2.startswith("Knopfzylinder_") and n2.split("_")[-1] == n1.split("_")[-1][1:]: continue
        v = s1.intersect(s2).Volume()
        if v > 0.1: kol.setdefault(n1, []).append((n2, round(v, 2)))
for k, v in sorted(kol.items()): print("KOLLISION", k, v[:6])
mech_d = dict(mech_pr); ged = {}
ARMLEN = Z_ACHSE - Z_GG; th_ = math.degrees(HUB_K / ARMLEN)
for n_, y_, a_, h_ in ACHSEN:
    x0_, x1_ = achs_x[n_]; ax_p, ax_q = cq.Vector(x0_, y_, Z_ACHSE), cq.Vector(x1_, y_, Z_ACHSE)
    pt_ = cq.Vertex.makeVertex(a_[0][1], y_, Z_GG).rotate(ax_p, ax_q, 1.0)
    sg_ = -1.0 if pt_.Center().y > y_ else 1.0
    dreh = lambda s_: s_.rotate(ax_p, ax_q, sg_ * th_)
    beweg = {f"Achse_{n_}": dreh(mech_d[f"Achse_{n_}"]), f"Feder_Achse_{n_}": dreh(mech_d[f"Feder_Achse_{n_}"]), f"Klappenhebel_{n_}": dreh(mech_d[f"Klappenhebel_{n_}"])}
    for (n2, nm_, gx_, gy_, r_, Xk_) in [q for q in arm_info if q[0] == n_]:
        for p_ in ("Schenkel_Knopf_", "Gelenkgabel_", "Gelenkbolzen_"): beweg[p_ + nm_] = dreh(mech_d[p_ + nm_])
        for p_ in ("Gelenkkopf_", "Schubstange_"): beweg[p_ + nm_] = mech_d[p_ + nm_].translate(cq.Vector(0, -HUB_K, 0))
        for nn_, s_ in teile:
            if nn_.startswith((f"Knopf_R{1 if r_ == 'A' else 2}_{Xk_:.0f}", f"Knopfzylinder_R{1 if r_ == 'A' else 2}_{Xk_:.0f}")): beweg[nn_] = s_.translate(cq.Vector(0, -HUB_K, 0))
    for (pax, pk) in hebel_pa[n_]:
        for nn_, s_ in teile:
            if nn_.startswith(("Klappe_", "Begleiterklappe_")):
                bb_ = s_.BoundingBox()
                if bb_.xmin < pk[0] < bb_.xmax and bb_.ymin < pk[1] < bb_.ymax: beweg[nn_] = dreh(s_)
    fest_ = [(nn_, s_) for nn_, s_ in gehaeuse + mech_pr if nn_ not in beweg]
    for n1, s1 in beweg.items():
        b1 = s1.BoundingBox()
        for n2, s2 in fest_:
            if "Schraube" in n2 or (n1.startswith("Achse_") and n2.startswith("Lagerbock")) or n2.startswith(("Filzeinlage", "Luftklappe")): continue
            if n1.startswith(("Knopf_", "Knopfzylinder_")) and n2.startswith(("Front", "Frontblende", "Front_Aufdoppelung", "Filz")): continue
            if not _bbt(b1, s2.BoundingBox()): continue
            v_ = s1.intersect(s2).Volume()
            if v_ > 0.1: ged.setdefault(n_, []).append((n1, n2, round(v_, 2)))
print("Gedrueckt-Pruefung: Winkel", round(th_, 1), "Grad, Armlaenge", round(ARMLEN, 1))
for k_, v_ in ged.items(): print("GEDRUECKT", k_, v_[:8])
teile = [(n_, boden if n_ == "Bassboden" else s_) for n_, s_ in teile] + mech
FARBEN[:0] = [("Lagerbock_Schraube", (0.40, 0.40, 0.42)), ("Lagerbock", (0.25, 0.75, 0.35)), ("Achse", (0.55, 0.55, 0.60)),
              ("Schenkel_Knopf", (0.20, 0.40, 0.85)), ("Klappenhebel", (0.10, 0.65, 0.30)), ("Schubstange", (0.10, 0.60, 0.60)),
              ("Gelenkgabel", (0.95, 0.45, 0.10)), ("Gelenkkopf", (0.60, 0.20, 0.70)), ("Gelenkbolzen", (0.95, 0.85, 0.20)), ("Feder_Achse", (0.75, 0.10, 0.75))]
