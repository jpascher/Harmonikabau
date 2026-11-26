# ---------- Bassmechanik Eduard Winkelbass (nach mech38 Marie; von eduard.py per exec eingebunden, vor der Spiegelung) ----------
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
def _bbt(b1, b2): return not (b1.xmax < b2.xmin or b2.xmax < b1.xmin or b1.ymax < b2.ymin or b2.ymax < b1.ymin or b1.zmax < b2.zmin or b2.zmax < b1.zmin)
R_Z = Z_RH                                                # Reihe 1 (Z 14), Reihe 2 (Z 28), Reihe 3 (Z 42)
_kl_m = [((s_.BoundingBox().xmin + s_.BoundingBox().xmax) / 2, (s_.BoundingBox().ymin + s_.BoundingBox().ymax) / 2)
         for n_, s_ in klappen if n_.endswith("_Leder") and n_.startswith(("Klappe_", "Begleiterklappe_"))]
def km(x_, y_): return min(_kl_m, key=lambda p_: math.dist(p_, (x_, y_)))   # Klappenmitte zum Zeichnungspunkt
# Gerhard (BassteilV3: 0-GB-/2-Begleiter-Schubstangen, 2-BassgestaengeQuer, Begleitergestaenge-quer). Reihe 1 und 3 teilen die
# X-Positionen 121..261: die Reihe-1-Stange laeuft um 5 versetzt (wie gezeichnet, 181 -> 186 statt 176 wegen der Lagerboecke).
# Achsreihen in den Streifen zwischen den Klappen, schmale Lagerboecke (4 mm): Begleiterbereich 5 Reihen, hinten 4, vorne 4.
YB = [15.6, 19.6, 23.6, 27.6]                             # hinter den Begleiterklappen (Deckelauflage hinten 2 mm)
YH = [57.6, 61.1, 64.6, 68.1]                             # zwischen Begleiter- und hinteren Grundbassklappen
YV = [104.9, 108.2, 111.5, 115.5]                         # zwischen hinteren und vorderen Grundbassklappen
# (Name, Y, [(Reihe, Knopf-X, Stangen-X)], [(Hebelansatz X, Klappe nach Zeichnungspunkt)], Wegpunkte Hebel)
# Zuordnung getauscht: Begleiter 12 (links vorne) von Knopf 121 (Reihe 3), Begleiterklappe 26 von Knopf 161 (Reihe 3).
_A = [
    ("B161k", YB[3], [(3, 161, 161)], [(38.0, (26, 42.6))], None),
    ("B201b", YB[3], [(1, 201, 206)], [(191.0, (191, 42.6))], None),
    ("B241", YB[3], [(3, 241, 241)], [(257.0, (257, 42.6))], None),
    ("B291", YB[3], [(2, 291, 291)], [(356.0, (356, 42.6))], None),
    ("B131", YB[2], [(2, 131, 131)], [(59.0, (59, 42.6))], None),
    ("B201", YB[2], [(3, 201, 201)], [(158.0, (158, 42.6))], None),
    ("B251", YB[2], [(2, 251, 251)], [(296.0, (290, 42.6))], None),      # Hebel neben dem Arm 291 vorbei
    ("B161", YB[1], [(1, 161, 166)], [(92.0, (92, 42.6))], None),
    ("B211", YB[1], [(2, 211, 211)], [(224.0, (224, 42.6))], None),
    ("B261", YB[1], [(1, 261, 266)], [(323.0, (323, 42.6))], None),
    ("B171", YB[0], [(2, 171, 171)], [(125.0, (125, 42.6))], None),
    # hinter den hinteren Grundbassklappen (Hg/Hf mit Torsionsfeder, gemeinsamer Bock zwischen den Armen 176/191)
    ("Hg", YH[0], [(1, 181, 176)], [(165.3, (167.5, 86.6))], None),
    ("Hf", YH[0], [(2, 191, 191)], [(214.6, (214.6, 86.6))], None),
    ("Hc", YH[0], [(2, 231, 231)], [(308.8, (308.8, 86.6))], None),
    ("HA1", YH[1], [(2, 111, 111), (3, 261, 261)], [(115.8, (120.4, 86.6))], None),
    ("HA2", YH[2], [(1, 121, 126), (2, 271, 271)], [(256.0, (261.7, 86.6))], None),
    ("Hj", YH[3], [(3, 101, 104), (1, 241, 246)], [(73.4, (73.3, 86.6))], None),   # Stange 101 -> 104 (Auflage Luftklappe)
    # hinter den vorderen Grundbassklappen
    ("B12", YV[0], [(3, 121, 121)], [(43.0, (43, 131.5))], None),          # Begleiter 12
    ("Hl", YV[0], [(1, 141, 146)], [(337.4, (337.4, 133.3))], None),        # liegende Kammer
    ("He", YV[1], [(2, 91, 91), (3, 221, 221)], [(191.0, (191.0, 133.3))], None),
    ("Hi", YV[2], [(3, 141, 141)], [(96.8, (96.8, 133.3))], None),
    ("Hd", YV[2], [(3, 181, 181)], [(210.0, (238.1, 133.3))], None),
    ("Hk", YV[2], [(1, 221, 226), (3, 281, 281)], [(285.3, (285.2, 133.3))], None),
    ("Hh", YV[3], [(2, 151, 151)], [(156.0, (143.9, 133.3))], None),
]
ACHSEN, KNOPF_X, HEBEL_WEG = [], {}, {}
for n_, y_, ks_, hs_, weg_ in _A:
    ACHSEN.append((n_, y_, [(r_, xs_) for r_, xk_, xs_ in ks_], [(xa_, km(*pz_)) for xa_, pz_ in hs_]))
    for r_, xk_, xs_ in ks_: KNOPF_X[(n_, r_, xs_)] = xk_
    if weg_: HEBEL_WEG[n_] = weg_
for q_ in ACHSEN: print("Achse", q_[0], q_[1], q_[2], [(round(h[0], 1), tuple(round(v, 1) for v in h[1])) for h in q_[3]])
# Hindernisse fuer Lagerboecke (Grundriss unter dem Boden)
_rund = lambda n_: n_.startswith("Begleiterklappe_")
rechtecke = [(b_.xmin - 0.5 - 0.3, b_.xmax + 0.5 + 0.3, b_.ymin - 0.5 - 0.3, b_.ymax + 0.5 + 0.3) for b_ in
             [s_.BoundingBox() for n_, s_ in klappen + lg if n_.endswith("_Leder") and not _rund(n_)]]
kreise = [((b_.xmin + b_.xmax) / 2, (b_.ymin + b_.ymax) / 2, (b_.xmax - b_.xmin) / 2 + 0.5 + 0.3) for b_ in
          [s_.BoundingBox() for n_, s_ in klappen if n_.endswith("_Leder") and _rund(n_)]]
rechtecke += [(b_.xmin - 0.3, b_.xmax + 0.3, b_.ymin - 0.3, b_.ymax + 0.3) for b_ in [s_.BoundingBox() for n_, s_ in teile if n_.startswith(("Eckblock_", "Eckwinkel", "Deckelauflage", "Lufthebel_Lagerbock"))]]
HL, HB = 4.3, 2.0                                         # halbe Bocklaenge / -breite (Gerhard: schmale Boecke 4 mm)
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
def hebel_seg(n2, pax, y2, pk):                           # Klappenhebel im Grundriss als Strecken (mit Wegpunkten)
    pts_ = [(pax, y2)] + HEBEL_WEG.get(n2, []) + [pk]
    return list(zip(pts_[:-1], pts_[1:]))
# Arme und Hebel als feste Hindernisse
arme_alle = [(n_, x_, y_) for n_, y_, a_, h_ in ACHSEN for r_, x_ in a_]
hebel_pa = {n_: [list(h) for h in h_] for n_, y_, a_, h_ in ACHSEN}
boecke = []                                               # (Achse, x, y)
def bock_frei(nm, x, y, eigene):
    r = (x - HL, x + HL, y - HB, y + HB)
    if y - HB < 11.8 or y + HB > FI - 0.3: return False
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
        if math.hypot(max(r[0] - ax, 0, ax - r[1]), max(r[2] - ay, 0, ay - r[3])) < D_STANGE / 2 + (0.8 if n2 != nm else 0.5):
            return False
    for n2, y2, a2, h2 in ACHSEN:                          # Klappenhebel (Grundriss), eigene tauchen unter dem Bock durch
        if n2 == nm: continue
        for (pax, pk) in hebel_pa[n2]:
            if any(_rect_seg(r, a_, b_) < D_STANGE / 2 + 0.3 for a_, b_ in hebel_seg(n2, pax, y2, pk)): return False
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
TORSION = {"Hg", "Hf"}                                   # Federn dort als Torsionsfeder neben dem Arm (Seite zum Hebel)
def _abst(kind, n_=None): return A_FREI if (kind != "A" or n_ in TORSION) else A_ARM
def _erlaubt(n_, x_, art):
    f_ = FEAT[n_]
    if any(abs(x_ - fx) < _abst(fk, n_) for fx, fk in f_): return False
    if art == "lo": return x_ < f_[0][0]
    if art == "hi": return x_ > f_[-1][0]
    return f_[0][0] < x_ < f_[-1][0]
bock_soll = {}
for n_, y_, a_, h_ in ACHSEN:
    lo_, hi_ = FEAT[n_][0], FEAT[n_][-1]
    bl_ = lo_[0] - _abst(lo_[1], n_)                          # Endboecke immer aussen; Schraeghebel taucht unter dem Bock durch
    bh_ = hi_[0] + _abst(hi_[1], n_)
    bock_soll[n_] = [[bl_, "lo"], [bh_, "hi"]]
    if False:                                             # kein Mittelbock (Schraube wuerde die Achse kreuzen)
        frei_ = [x_ for x_ in [bl_ + 0.5 * k for k in range(int((bh_ - bl_) * 2))] if _erlaubt(n_, x_, "mid")]
        if frei_: bock_soll[n_].append([min(frei_, key=lambda q: abs(q - (bl_ + bh_) / 2)), "mid"])
# hintere Reihe (Begleiter, Y 21,1 und 26,1): jede Reihe in einer Linie auf gemeinsamen Lagerboecken, keine Doppelboecke mit
# der Nachbarreihe (5 mm Reihenabstand) -> Boecke der beiden Reihen gegeneinander versetzt
FEST_BOCK = {"B201b": {"hi": 222.3}, "B241": {"lo": 231.3}, "B201": {"hi": 210.7}, "Hi": {"hi": 150.7}, "Hd": {"lo": 171.3},
             "Hh": {"lo": 141.3, "hi": 162.1}, "He": {"hi": 233.0}}
FEDER_ARM = {"Hk": 281}                                   # Feder an diesem Arm (Platz zur Nachbarreihe)     # Boecke benachbarter Reihen gegeneinander versetzt
for n_, d_ in FEST_BOCK.items():
    for e in bock_soll[n_]:
        if e[1] in d_: e[0] = d_[e[1]]; e[1] = "fest"
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
                if e[1] not in ("geteilt", "fest") and not _erlaubt(n_, x_, e[1]): continue
                if e[1] in ("geteilt", "fest") and k > 0: continue
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
    yz_ = [(-g_ * (-q[1] - 2.6) * HB / 2.6, Z_UB - (q[0] + 13.6)) for q in PROF_E]
    lb_ = cq.Workplane("YZ").polyline(yz_).close().extrude(8.6).val().translate(cq.Vector(-4.3, 0, 0))
    lb_ = lb_.cut(stab((-5, 0, Z_ACHSE), (5, 0, Z_ACHSE), D_ACHSE + 0.1)).cut(stab((0, 0, Z_UB - H_LB - 1), (0, 0, Z_UB + 1), 1.5))
    sc_ = stab((0, 0, Z_UB - H_LB - 0.9), (0, 0, Z_UB + 3.0), 1.4).fuse(stab((0, 0, Z_UB - H_LB - 0.9), (0, 0, Z_UB - H_LB), 2.5))
    return lb_, sc_
# Doppelboecke: zwei Boecke benachbarter paralleler Achsen (5 mm Abstand), die sich ueberschneiden -> ein gemeinsamer Bock
doppel = []
for i_ in range(len(boecke)):
    for j_ in range(i_ + 1, len(boecke)):
        (n1_, x1_, y1_), (n2_, x2_, y2_) = boecke[i_], boecke[j_]
        if False and n1_ != n2_ and abs(x1_ - x2_) < 2 * HL + 0.3 and abs(y1_ - y2_) < 2 * HB + 0.3 and abs(y1_ - y2_) > 0.05:
            doppel.append((i_, j_))
for i_, j_ in doppel:
    xm_ = (boecke[i_][1] + boecke[j_][1]) / 2
    for k_ in (i_, j_):
        n_, xo_, yo_ = boecke[k_]
        for e in bock_soll[n_]:
            if abs(e[0] - xo_) < 0.01: e[0] = xm_
        boecke[k_] = (n_, xm_, yo_)
    print("Doppelbock", boecke[i_][0], boecke[j_][0], "bei X", round(xm_, 1))
_achs_x_neu()
_dop = {k_ for p_ in doppel for k_ in p_}
for i_, j_ in doppel:
    xm_, ya_, yb_ = boecke[i_][1], min(boecke[i_][2], boecke[j_][2]), max(boecke[i_][2], boecke[j_][2])
    db_ = box(xm_ - HL, xm_ + HL, ya_ - HB, yb_ + HB, Z_UB - H_LB, Z_UB)
    for y_ in (ya_, yb_): db_ = db_.cut(stab((xm_ - 5, y_, Z_ACHSE), (xm_ + 5, y_, Z_ACHSE), D_ACHSE + 0.1))
    ym_ = (ya_ + yb_) / 2
    db_ = db_.cut(stab((xm_, ym_, Z_UB - H_LB - 1), (xm_, ym_, Z_UB + 1), 1.5))
    sc_ = stab((xm_, ym_, Z_UB - H_LB - 0.9), (xm_, ym_, Z_UB + 3.0), 1.4).fuse(stab((xm_, ym_, Z_UB - H_LB - 0.9), (xm_, ym_, Z_UB - H_LB), 2.5))
    mech.append((f"Lagerbock_Doppel_{boecke[i_][0]}_{boecke[j_][0]}", db_)); mech.append((f"Lagerbock_Schraube_Doppel_{boecke[i_][0]}_{boecke[j_][0]}", sc_))
    boden = boden.cut(stab((xm_, ym_, Z_UB - 1), (xm_, ym_, Z_UB + 3.0), 1.1))
for k_, (n_, bx_, by_) in enumerate(boecke):
    if k_ in _dop: continue
    nb_ = [y2 for n2, y2, a2, h2 in ACHSEN if n2 != n_ and achs_x[n2][0] - 4 < bx_ < achs_x[n2][1] + 4 and abs(y2 - by_) < 8]
    g_ = 1 if (nb_ and min(nb_, key=lambda q: abs(q - by_)) > by_) else -1
    lb, sc = bock_solid(g_)
    mech.append((f"Lagerbock_{n_}_{k_+1}", lb.translate(cq.Vector(bx_, by_, 0))))
    mech.append((f"Lagerbock_Schraube_{n_}_{k_+1}", sc.translate(cq.Vector(bx_, by_, 0))))
    boden = boden.cut(stab((bx_, by_, Z_UB - 1), (bx_, by_, Z_UB + 3.0), 1.1))      # Vorbohrung im Boden
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
# Luftgruppe in beiden Stellungen als Hindernis: Gestaenge laeuft zwischen gedrueckter und ungedrueckter Lage durch
LUFT_HI = [s_ for n_, s_ in lg + lg_g + [("Lufttaste_Filz", lt_filz)] if not n_.startswith("Lufthebel_Feder")]
LUFT_HI += [box(*(lambda b_: (b_.xmin, b_.xmax, b_.ymin, b_.ymax, b_.zmin, b_.zmax))(s_.BoundingBox())) for n_, s_ in lg + lg_g if n_.startswith("Lufthebel_Feder")]
Y_BD0, Y_BD1 = 152.0, 188.0                               # Band der Luftgruppe (Klappe, Hebel, Boecke, Taste)
Z_UNTEN = T_DECKEL_BASS + 0.5                              # Deckel + Bassbodenleder
R_FREI = D_STANGE / 2 + 0.5
def z_band(x_):
    """niedrigste freie Hoehe fuer eine Stange bei X im Band (None: kein Durchgang)"""
    sl_ = box(x_ - R_FREI, x_ + R_FREI, Y_BD0, Y_BD1, -20, Z_BU)
    iv_ = []
    for h_ in LUFT_HI:
        if not _bbt(sl_.BoundingBox(), h_.BoundingBox()): continue
        c_ = sl_.intersect(h_)
        if c_.Volume() > 1e-4:
            for so_ in c_.Solids(): b_ = so_.BoundingBox(); iv_.append((b_.zmin, b_.zmax))
    z_ = max(Z_GG, Z_UNTEN + R_FREI)
    while z_ < Z_BU - 3:
        if all(z_ + R_FREI <= a_ or z_ - R_FREI >= b_ for a_, b_ in iv_): return round(z_, 2)
        z_ += 0.25
    return None
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
    for r_, Xs in a_:
        Xk = KNOPF_X[(n_, r_, Xs)]
        Zr = R_Z[r_]; zyl_ = [s_ for nn_, s_ in teile if nn_ == f"Knopfzylinder_R{r_}_{Xk:.0f}"][0]
        y_end = zyl_.BoundingBox().ymin + 0.6
        gx, gy = Xs, y_
        sk_ = stab((gx, gy, Z_ACHSE), (gx, gy, Z_GG + 4.5), D_STANGE).cut(stab((gx - 3, gy, Z_ACHSE), (gx + 3, gy, Z_ACHSE), D_ACHSE))
        mech.append((f"Schenkel_Knopf_{n_}_{r_}{Xk:.0f}", sk_))
        Z_BMAX = 20.5                                     # Band zwischen gedrueckter Taste und ungedruecktem Hebel
        jx_, zb_ = 0.0, z_band(Xs)
        if zb_ is not None and zb_ > Z_BMAX: zb_ = None
        if zb_ is None:                                   # Gabelgelenk/Stange der Lufttaste: seitlich ausweichen
            for d_ in [0.5 * k for k in range(1, 17)]:
                for sg_ in (-1, 1):
                    zz_ = z_band(Xs + sg_ * d_)
                    if zz_ is not None and zz_ <= Z_BMAX: jx_, zb_ = sg_ * d_, zz_; break
                if zb_ is not None: break
        if abs(zb_ - Z_GG) < 0.01 and jx_ == 0:          # frei: wie Standard (auf Gelenkhoehe, vor dem Knopf hinauf)
            y_s = min(y_end - 6.0, FI - 1.0)              # erst innerhalb des Kastens absenken (Reihe 2/3: Frontausschnitt)
            y_a = y_s - (Zr - Z_GG) * (0.364 if Zr > 35 else 0.577 if Zr > 20 else 1.0)   # Reihe 3 steil (Luftgruppe)
            pfad = [(gx, gy + 4.6, Z_GG), (gx, y_a, Z_GG), (gx, y_s, Zr), (Xk, y_end, Zr)]   # Versatz im waagrechten Stueck
        else:                                             # im Band zwischen gedrueckter und ungedrueckter Luftgruppe
            xb_ = Xs + jx_; d1_ = max(abs(zb_ - Z_GG), abs(jx_)); d2_ = max(abs(Zr - zb_), abs(jx_ + Xs - Xk))
            y_k = min(Y_BD1, y_end - 2.0 - d2_)            # Reihe 3: steigt schon im Band zur Knopfhoehe
            pfad = [(gx, gy + 4.6, Z_GG), (gx, Y_BD0 - d1_, Z_GG), (xb_, Y_BD0, zb_), (xb_, y_k, zb_), (Xk, y_k + d2_, Zr), (Xk, y_end, Zr)]
            print(f"Stange {n_} Knopf {Xk:.0f}: Band Z {zb_:.2f}" + (f", seitlich {jx_:+.1f}" if jx_ else ""))
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
        if n_ in HEBEL_WEG:                                # Hebel mit Wegpunkten (zwischen den Begleiterklappen durch)
            pts_ = [(pa[0], pa[1], Z_ACHSE), (pa[0], pa[1], Z_KH)] + [(wx_, wy_, Z_KH) for wx_, wy_ in HEBEL_WEG[n_]] + [(pk[0], pk[1], Z_KH)]
        elif not unter_:
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
            nf_ = [h for h in hx_ + [x for x in ax_ if x != gx_] if min(gx_, bn_) < h < max(gx_, bn_)]
            if not nf_ and abs(bn_ - gx_) >= A_ARM - 0.1:
                kand.append((1, -abs(bn_ - gx_), gx_, sd))          # offenes Ende zum Lagerbock
            elif nf_ and min(abs(h - gx_) for h in nf_) >= 6.5:
                kand.append((0, -min(abs(h - gx_) for h in nf_), gx_, sd))   # Torsionsfeder zur Hebelseite
    if n_ in FEDER_ARM: kand = [k for k in kand if abs(k[2] - FEDER_ARM[n_]) < 0.5] or kand
    if kand:
        _, _, gx, s_dir = max(kand)
    else:
        gx, s_dir = ax_[0], 1.0; print("WARNUNG Feder ohne offenes Ende", n_)
    t0 = gx + s_dir * (D_STANGE / 2 + 0.4); L_SP = N_WIND * P_WIND
    hel = cq.Wire.makeHelix(P_WIND, L_SP, R_WIND)
    sp = cq.Workplane("XZ").center(R_WIND, 0).circle(D_DRAHT / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90)
    if s_dir < 0: sp = sp.mirror("YZ")                  # Gegenseite: Feder spiegeln (Links- statt Rechtswindung)
    sp = sp.translate(cq.Vector(t0, y_, Z_ACHSE))
    nx_, ny_ = 0.0, -1.0
    z_h = Z_GG + L_GABEL + 1.0; ab = D_STANGE / 2 + 0.1 + D_DRAHT / 2   # Schenkel ueber der Gabel am Arm: knapp neben dem Arm
    s1 = draht([(t0, y_, Z_ACHSE - R_WIND), (t0, y_ - 0.6, Z_ACHSE - R_WIND - 3.0), (t0, y_ - ab, Z_ACHSE - R_WIND - 5.0), (t0, y_ - ab, z_h),
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
    pp1_ = [(t0, y_, Z_ACHSE - R_WIND), (t0, y_ - 0.6, Z_ACHSE - R_WIND - 3.0), (t0, y_ - ab, Z_ACHSE - R_WIND - 5.0), (t0, y_ - ab, z_h), (gx - s_dir * ab, y_ - ab, z_h), (gx - s_dir * ab, y_ + 0.8, z_h)]
    for q0_, q1_ in list(zip(pp1_[:-1], pp1_[1:])) + [((tE, y_, Z_ACHSE + R_WIND), (tE, y_ + 0.8, Z_UB - 0.3))]: hk_ = hk_.fuse(stab(q0_, q1_, D_DRAHT))
    feder_ersatz[f"Feder_Achse_{n_}"] = hk_
print("Mechanik Teile", len(mech), "| Achsen", {n_: [round(v, 1) for v in achs_x[n_]] for n_, *_ in ACHSEN})
print("Hebelansaetze", {n_: [round(h[0], 1) for h in hebel_pa[n_]] for n_, *_ in ACHSEN})
# ---- Pruefung: Mechanik gegen Gehaeuse in Ruhe, dann jede Taste gedrueckt (Lufttaste gedrueckt) ----
# Deckel-Auflagestreifen unter den Gelenken der hinteren Reihe ausgenommen (1,3 bleibt stehen)
_ausn = [box(gx_ - 4.5, gx_ + 4.5, gy_ - HUB_K - 4.0, gy_ + 17.5, T_DECKEL_BASS - 0.1, Z_GG + 5) for (n2, nm_, gx_, gy_, r_, Xk_) in arm_info if gy_ < 40]
teile = [(n_, s_.cut(cq.Compound.makeCompound(_ausn)) if n_ == "Bassdeckel_Auflage_hinten" else s_) for n_, s_ in teile]
gehaeuse = [("Bassboden_neu", boden)] + [(n_, s_) for n_, s_ in teile + lg_g if n_ != "Bassboden"]
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
mech_d = dict(mech_pr); ged = {}; alle_ged = {}
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
            if nn_.startswith((f"Knopf_R{r_}_{Xk_:.0f}", f"Knopfzylinder_R{r_}_{Xk_:.0f}")): beweg[nn_] = s_.translate(cq.Vector(0, -HUB_K, 0))
    for (pax, pk) in hebel_pa[n_]:
        for nn_, s_ in klappen:
            if nn_.startswith(("Klappe_", "Begleiterklappe_")):
                bb_ = s_.BoundingBox()
                if bb_.xmin < pk[0] < bb_.xmax and bb_.ymin < pk[1] < bb_.ymax: beweg[nn_] = dreh(s_)
    alle_ged[n_] = dict(beweg)
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
# ---- alle Tasten gleichzeitig gedrueckt (dazu Lufttaste gedrueckt): bewegte Teile untereinander und gegen das Feste ----
_bw = {nn_: (n_, s_) for n_, d_ in alle_ged.items() for nn_, s_ in d_.items()}
_fest = [(nn_, s_) for nn_, s_ in gehaeuse + mech_pr if nn_ not in _bw and not nn_.startswith(("Lufttaste", "Lufthebel", "Luftklappe")) or nn_.endswith("_gedrueckt")]
_alle = []
_bl = list(_bw.items())
for i_, (n1, (a1, s1)) in enumerate(_bl):
    b1 = s1.BoundingBox()
    geg_ = [(n2, s2) for n2, s2 in _fest] + [(n2, s2) for n2, (a2, s2) in _bl[i_ + 1:] if a2 != a1]
    for n2, s2 in geg_:
        if "Schraube" in n2 or (n1.startswith("Achse_") and n2.startswith("Lagerbock")) or n2.startswith("Filzeinlage"): continue
        if n1.startswith(("Knopf_", "Knopfzylinder_")) and n2.startswith(("Front", "Frontblende", "Front_Aufdoppelung", "Filz")): continue
        if not _bbt(b1, s2.BoundingBox()): continue
        v_ = s1.intersect(s2).Volume()
        if v_ > 0.1: _alle.append((n1, n2, round(v_, 2)))
print("ALLE GEDRUECKT (inkl. Lufttaste): bewegte Teile", len(_bl), "| Kollisionen", len(_alle))
for q_ in _alle: print("ALLE GEDRUECKT", q_)
teile = [(n_, boden if n_ == "Bassboden" else s_) for n_, s_ in teile] + mech
FARBEN[:0] = [("Lagerbock_Schraube", (0.40, 0.40, 0.42)), ("Lagerbock", (0.25, 0.75, 0.35)), ("Achse", (0.55, 0.55, 0.60)),
              ("Schenkel_Knopf", (0.20, 0.40, 0.85)), ("Klappenhebel", (0.10, 0.65, 0.30)), ("Schubstange", (0.10, 0.60, 0.60)),
              ("Gelenkgabel", (0.95, 0.45, 0.10)), ("Gelenkkopf", (0.60, 0.20, 0.70)), ("Gelenkbolzen", (0.95, 0.85, 0.20)), ("Feder_Achse", (0.75, 0.10, 0.75))]
