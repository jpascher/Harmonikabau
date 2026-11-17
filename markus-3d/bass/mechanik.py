# ---------- Bassmechanik aus Markus.Basteil-Mechanik3-fehler behoben.dxf (Aufsicht) ----------
# Wird von kammern.py per exec() eingebunden (nutzt box, cq, ezdxf, math, B, D, L, teile, klappen, boden, front).
# Plan -> Modell: X = x - 2923.3, Y = y - 1520.7 (wie Hauptzeichnung).
# Hoehen (Z ab Deckelseite): Boden-Unterseite 25.2; Lagerboecke 5 mm hoch darunter angeschraubt, Achsmitte 22.7;
# Klappenhebel von der Achse zur Klappenmitte knapp unter der Klappe; zweiter Arm nach unten Richtung Deckel,
# dort Gabelgelenk; Schubstange Ø2 vom Gabelgelenk zum Knopfzylinder; Knopf: Zylinder Ø7.5 x 20 + Kopf Ø15.5.
MX, MY = 2923.3, 1520.7
MD = ezdxf.readfile("mech.dxf").modelspace()
# Hintere Achsgruppe (Plan Y < 70) Richtung Rueckteil versetzen, damit Lagerboecke/Achsen nicht an die hinteren Klappen stossen
DY_HINTEN = -3.5
vy = lambda y: y + DY_HINTEN if y < 70 else y
Z_UB = Z_BODEN - T_BODEN                 # Boden-Unterseite 25.2
H_LB, Z_ACHSE = 5.0, Z_BODEN - T_BODEN - 2.5
Z_GG = 10.0                              # Hoehe Gabelgelenk (Armende Richtung Deckel)
D_ACHSE, D_STANGE = 2.0, 2.0
Z_KL_U = min(s.BoundingBox().zmin for _, s in klappen)   # Unterkante Klappen
# Abstufung je 7 mm: Oberflaeche Reihe 1 = Taschenboden Y 200, Reihe 2 = Front Y 207, Reihe 3 = Frontblende Y 214.
# Knopf-Ruhelage wie im Knopfschnitt der Mechanikzeichnung: 6 mm Hub -> Kopf 6.5 mm vor der jeweiligen Oberflaeche.
# Filzring sitzt in der aeussersten Lage der jeweiligen Stufe.
HUB, SPALT = 6.0, 6.5
STUFE = {"A": (200.0, 5.5), "B": (207.0, 5.5), "C": (214.0, 7.0)}            # Oberflaeche Y, Dicke der aeussersten Lage
REIHEN = {r_: (z_, STUFE[r_][0] + SPALT) for r_, z_ in (("A", 10.0 + 4.24), ("B", 24.0 + 4.24), ("C", 38.0 + 4.24))}
mech, mech_teile = [], []
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v.Length, p0, v.normalized())
# --- Achsen ---
achsen = []
for e in MD.query('LINE[layer=="Gestänge quer"]'):
    a0 = (e.dxf.start.x - MX, vy(e.dxf.start.y - MY)); a1 = (e.dxf.end.x - MX, vy(e.dxf.end.y - MY))
    achsen.append((a0, a1))
def fuss(p, a):
    (x0, y0), (x1, y1) = a; dx, dy = x1 - x0, y1 - y0; t = max(0, min(1, ((p[0]-x0)*dx + (p[1]-y0)*dy) / (dx*dx + dy*dy)))
    f = (x0 + t*dx, y0 + t*dy); return f, math.dist(p, f)
# --- Lagerboecke (lg2) -> Profil aus lager-bohrung-beide.dxf, einfach 8.6 x 5.2 x 5 ---
lbs = []
for e in MD.query('INSERT[name=="lg2"]'):
    a = math.radians(e.dxf.rotation); px, py = 301.84, 138.06
    lbs.append((e.dxf.insert.x + px*math.cos(a) - py*math.sin(a) - MX, vy(e.dxf.insert.y + px*math.sin(a) + py*math.cos(a) - MY), e.dxf.rotation))
n_lb_neu = 0
for (a0, a1) in achsen:                       # jedes Achsende braucht einen Bock auf der Achslinie (innerhalb 15 mm)
    dx_, dy_ = a1[0] - a0[0], a1[1] - a0[1]; Ln_ = math.hypot(dx_, dy_)
    ts_ = [((b_[0]-a0[0])*dx_ + (b_[1]-a0[1])*dy_) / Ln_ for b_ in lbs if abs((b_[0]-a0[0])*dy_ - (b_[1]-a0[1])*dx_) / Ln_ < 2.5]
    for t_end, p_ in ((0.0, a0), (Ln_, a1)):
        if not any(abs(t_ - t_end) < 15.0 for t_ in ts_):
            rot_ = math.degrees(math.atan2(dy_, dx_))
            rot_ = rot_ - 180 if rot_ > 90 else rot_ + 180 if rot_ < -90 else rot_
            lbs.append((p_[0], p_[1], rot_)); n_lb_neu += 1
            print("Lagerbock ergaenzt bei X", round(L - p_[0], 1), "(Modellansicht) Y", round(p_[1], 1))
# --- Doppelte Lagerboecke = minimaler Achsabstand 3.2 mm (lager-bohrung-beide.dxf): Paare zusammenfassen, Achsen nachziehen ---
A_DOPPEL = 3.2
paare, in_paar = [], set()
for i in range(len(lbs)):
    for j in range(i + 1, len(lbs)):
        if i in in_paar or j in in_paar: continue
        if math.dist(lbs[i][:2], lbs[j][:2]) < 6.0:
            paare.append((i, j)); in_paar |= {i, j}
ziel = {}                                       # Bock-Index -> Bohrungsmitte im Doppelbock
for i, j in paare:
    (x1, y1, _), (x2, y2, _) = lbs[i], lbs[j]
    mx_, my_ = (x1 + x2) / 2, (y1 + y2) / 2; d_ = math.dist((x1, y1), (x2, y2)); ux, uy = (x2 - x1) / d_, (y2 - y1) / d_
    ziel[i] = (mx_ - ux * A_DOPPEL / 2, my_ - uy * A_DOPPEL / 2); ziel[j] = (mx_ + ux * A_DOPPEL / 2, my_ + uy * A_DOPPEL / 2)
n_versch = 0
for ai, (a0, a1) in enumerate(achsen):
    dx_, dy_ = a1[0] - a0[0], a1[1] - a0[1]; Ln_ = math.hypot(dx_, dy_); nx_, ny_ = -dy_ / Ln_, dx_ / Ln_
    neu_ = [list(a0), list(a1)]
    for k_, (t_end, p_) in enumerate(((0.0, a0), (Ln_, a1))):
        for bi, b_ in enumerate(lbs):
            if bi not in ziel: continue
            t_ = ((b_[0]-a0[0])*dx_ + (b_[1]-a0[1])*dy_) / Ln_; dn_ = (b_[0]-a0[0])*nx_ + (b_[1]-a0[1])*ny_
            if abs(dn_) < 2.5 and abs(t_ - t_end) < 15.0:
                sh_ = (ziel[bi][0] - b_[0]) * nx_ + (ziel[bi][1] - b_[1]) * ny_
                neu_[k_][0] += sh_ * nx_; neu_[k_][1] += sh_ * ny_; n_versch += 1
                break
    achsen[ai] = (tuple(neu_[0]), tuple(neu_[1]))
for i, (a0, a1) in enumerate(achsen):
    mech.append((f"Achse_{i+1}", stab((a0[0], a0[1], Z_ACHSE), (a1[0], a1[1], Z_ACHSE), D_ACHSE)))
def bock_solid(breite, bohrungen):
    lb_ = box(-4.3, 4.3, -breite / 2, breite / 2, Z_UB - H_LB, Z_UB); sc_ = None
    for yb in bohrungen:
        lb_ = lb_.cut(stab((-5, yb, Z_ACHSE), (5, yb, Z_ACHSE), 2.4)).cut(stab((0, yb, Z_UB - H_LB - 1), (0, yb, Z_UB + 1), 1.5))
        s_ = stab((0, yb, Z_UB - H_LB - 0.9), (0, yb, Z_UB - H_LB + 10), 1.4).fuse(stab((0, yb, Z_UB - H_LB - 0.9), (0, yb, Z_UB - H_LB), 2.5))
        sc_ = s_ if sc_ is None else sc_.fuse(s_)
    return lb_, sc_
nb = 0
for i, (cx, cy, rot) in enumerate(lbs):
    if i in in_paar: continue
    lb, sc = bock_solid(5.2, [0.0])
    tr = lambda s: s.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), rot).translate(cq.Vector(cx, cy, 0))
    nb += 1; mech.append((f"Lagerbock_{nb}", tr(lb))); mech.append((f"Lagerbock_Schraube_{nb}", tr(sc)))
for i, j in paare:
    (x1, y1, _), (x2, y2, _) = lbs[i], lbs[j]
    mx_, my_ = (x1 + x2) / 2, (y1 + y2) / 2
    rot = math.degrees(math.atan2(y2 - y1, x2 - x1)) - 90.0          # Achsrichtung senkrecht zur Paar-Verbindung
    lb, sc = bock_solid(8.4, [-A_DOPPEL / 2, A_DOPPEL / 2])
    tr = lambda s: s.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), rot).translate(cq.Vector(mx_, my_, 0))
    nb += 1; mech.append((f"Lagerbock_doppelt_{nb}", tr(lb))); mech.append((f"Lagerbock_Schraube_{nb}", tr(sc)))
print("Lagerboecke: einfach", len(lbs) - 2 * len(paare), "doppelt", len(paare), "| Achsenden auf 3.2 mm nachgezogen:", n_versch)
# --- Knoepfe (Knopfmitte X wie Frontbohrungen) ---
knoepfe = []
for c in RD.query('CIRCLE[layer=="bohrungen-bassdrücker"]'):
    x, y = c.dxf.center.x, c.dxf.center.y
    if -382 < x < 0 and -5 < y < 72:
        r = {10.0: "A", 24.0: "B", 38.0: "C"}[round(y, 1)]
        if (round(-x, 1), r) not in [(k[0], k[1]) for k in knoepfe]: knoepfe.append((round(-x, 1), r))
D_ZYL, D_KOPF, H_KOPF = 7.5, 15.5, 4.5
L_ZYL_R = {"A": 14.0, "B": 20.0, "C": 20.0}     # Reihe 1 gekuerzt (Lufttaste)
for X_, r in knoepfe:
    Zr, Yk = REIHEN[r]
    zyl = stab((X_, Yk - L_ZYL_R[r], Zr), (X_, Yk, Zr), D_ZYL)
    kopf = stab((X_, Yk, Zr), (X_, Yk + H_KOPF, Zr), D_KOPF).fuse(
        cq.Solid.makeSphere(50.27, cq.Vector(X_, Yk + H_KOPF - 50.27 + 0.6, Zr)).intersect(stab((X_, Yk + H_KOPF - 0.01, Zr), (X_, Yk + H_KOPF + 1, Zr), D_KOPF)))
    mech.append((f"Knopf_{r}_{X_:.0f}", kopf)); mech.append((f"Knopfzylinder_{r}_{X_:.0f}", zyl))
    yo_, dl_ = STUFE[r]
    filz = stab((X_, yo_ - dl_, Zr), (X_, yo_, Zr), 10.0).cut(stab((X_, yo_ - dl_ - 1, Zr), (X_, yo_ + 1, Zr), D_ZYL))
    mech.append((f"Filzring_{r}_{X_:.0f}", filz))
# --- Schubstangen + Arm + Gabelgelenk ---
stangen = []
for e in MD.query('*[layer=="GB-Schubstangen" | layer=="Begleiter-Schubstangen"]'):
    pts = [(q.x - MX, q.y - MY) for q in ezdxf.path.make_path(e).flattening(0.1)]
    hinten, vorne = min(pts, key=lambda p: p[1]), max(pts, key=lambda p: p[1])
    hinten = (hinten[0], vy(hinten[1]))
    r = "A" if vorne[1] < 215 else "B" if vorne[1] < 220.5 else "C"
    kn = min([k for k in knoepfe if k[1] == r], key=lambda k: abs(k[0] - vorne[0]))
    stangen.append((hinten, r, kn[0]))
n_ok = 0
# Klappen-Grundrisse (fuer "nicht mit Klappe kollidieren") und Lufttasten-Langloch im Deckel
kl_bb = [s_.BoundingBox() for n_, s_ in klappen if n_.endswith("Leder")]
LT = (lt_x0, lt_x1, lt_y0, lt_y1) = (LT_CX - LT_L / 2 - 2, LT_CX + LT_L / 2 + 2, LT_CY - LT_B / 2 - 1.5, LT_CY + LT_B / 2 + 1.5)
Z_UNTER = Z_UB - 1.4                    # Stangenhoehe knapp unter dem Boden
stange_achse = {}
stange_pfad = {}
arme = {}
for i, (P, r, Xk) in enumerate(stangen):
    ai = min(range(len(achsen)), key=lambda j: fuss(P, achsen[j])[1]); F, dist = fuss(P, achsen[ai])
    Zr, Yk = REIHEN[r]
    gx, gy = P[0], F[1]
    stange_achse.setdefault(ai, []).append((r, Xk))
    arme.setdefault(ai, []).append((gx, gy))
    # Gabelgelenk zweiteilig: Gabelkopf (an der Schubstange) + Gelenkbolzen; Drehpunkt = Bolzenachse am Armende
    bolzen_loch = stab((gx - 3.0, gy, Z_GG), (gx + 3.0, gy, Z_GG), 1.5)
    arm_ = box(gx - 0.75, gx + 0.75, gy - 1.5, gy + 1.5, Z_GG - 1.75, Z_ACHSE + 1.0).cut(bolzen_loch)
    mech.append((f"Arm_{i+1}", arm_))
    wange = lambda xa, xb: box(xa, xb, gy - 2.5, gy + 3.0, Z_GG - 2.5, Z_GG + 2.5)
    gk = wange(gx - 1.85, gx - 0.85).fuse(wange(gx + 0.85, gx + 1.85)).fuse(box(gx - 1.85, gx + 1.85, gy + 3.0, gy + 5.0, Z_GG - 2.5, Z_GG + 2.5))
    gk = gk.fuse(stab((gx, gy + 5.0, Z_GG), (gx, gy + 11.0, Z_GG), 4.0)).cut(bolzen_loch)
    mech.append((f"Gabelkopf_{i+1}", gk))
    bz = stab((gx - 2.3, gy, Z_GG), (gx + 2.1, gy, Z_GG), 1.4).fuse(stab((gx - 2.3, gy, Z_GG), (gx - 1.9, gy, Z_GG), 2.4))
    mech.append((f"Gelenkbolzen_{i+1}", bz))
    # Verlauf: unten (Z_GG) bis hinter die letzte Klappe im Weg, 45 Grad hoch bis knapp unter den Boden,
    # dort zum Knopf (X-Versatz), zuletzt 45 Grad hoch auf Knopfhoehe; Lufttasten-Langloch wird unten nicht ueberquert.
    y0 = gy + 11.0; y_end = Yk - L_ZYL_R[r]
    xs = (min(gx, Xk) - 2.5, max(gx, Xk) + 2.5)
    y_frei = max([b_.ymax + 2.0 for b_ in kl_bb if b_.xmax > xs[0] and b_.xmin < xs[1] and b_.ymax > y0] + [y0])
    z_mid = min(Zr, Z_UNTER)
    if xs[1] > lt_x0 and xs[0] < lt_x1 and y0 < lt_y1:
        z_mid = max(z_mid, H_TASTE + 2.0)                                # Bogen ueber die gedrueckte Lufttaste
        if y_frei > lt_y0 - (z_mid - Z_GG) - 1.0 and z_mid < Z_KL_U - 1.5:
            y_frei = lt_y0 - (z_mid - Z_GG) - 1.0                        # schon unter der Klappe hochgehen (unter Klappen-Unterkante)
    y1 = max(y_frei, y0)
    y2 = y1 + (z_mid - Z_GG)
    y3 = y_end - (Zr - z_mid)
    if y3 < y2: y3 = y2
    pfad = [(gx, y0, Z_GG), (gx, y1, Z_GG), (gx, y2, z_mid), (Xk, y3, z_mid), (Xk, y_end, Zr)]
    pfad = [p for k_, p in enumerate(pfad) if k_ == 0 or math.dist(p, pfad[k_ - 1]) > 0.05]
    st = None
    for p_a, p_b in zip(pfad[:-1], pfad[1:]):
        seg = stab(p_a, p_b, D_STANGE).fuse(cq.Solid.makeSphere(D_STANGE / 2, cq.Vector(*p_b)))
        st = seg if st is None else st.fuse(seg)
    mech.append((f"Schubstange_{r}_{Xk:.0f}", st)); stange_pfad[f"Schubstange_{r}_{Xk:.0f}"] = pfad
    if dist < 6: n_ok += 1
print("Mechanik: Achsen", len(achsen), "Lagerboecke", len(lbs), "Knoepfe", len(knoepfe), "Schubstangen", len(stangen), "davon an Achse (<6 mm)", n_ok)
# --- Klappenhebel: Plan-Linien Achse -> Klappenmitte; von der Achse senkrecht nach unten unter die Klappe, dann waagrecht ---
gesehen = set(); n_kh = 0; kh_liste = []; kh_namen = []
Z_KH = Z_KL_U - 1.2
for e in MD.query('*[layer=="Klappen-hebel"]'):
    pts = [(q.x - MX, q.y - MY) for q in ezdxf.path.make_path(e).flattening(0.1)]
    if len(pts) < 2: continue
    p0, p1 = pts[0], pts[-1]
    if math.dist(p0, p1) < 5 or abs(p0[1] - p1[1]) < 3: continue
    key = (round(p0[0]), round(p0[1]), round(p1[0]), round(p1[1]))
    if key in gesehen or key[2:] + key[:2] in gesehen: continue
    gesehen.add(key)
    da = min(fuss(p0, a_)[1] for a_ in achsen); db = min(fuss(p1, a_)[1] for a_ in achsen)
    if da < db: pk, pa = p1, pts[-2]                     # Hebel = letztes Segment zur Klappe
    else: pk, pa = p0, pts[1]
    if math.dist(pa, pk) < 3: pa = p0 if da < db else p1
    pa = (pa[0], vy(pa[1]))
    pa = fuss(pa, min(achsen, key=lambda a_: fuss(pa, a_)[1]))[0]
    n_kh += 1
    kh_liste.append((pa, pk)); kh_namen.append(f"Klappenhebel_{n_kh}")
    kh = stab((pa[0], pa[1], Z_ACHSE + 1.0), (pa[0], pa[1], Z_KH), 2.0).fuse(stab((pa[0], pa[1], Z_KH), (pk[0], pk[1], Z_KH), 2.0))
    mech.append((f"Klappenhebel_{n_kh}", kh))
print("Klappenhebel", n_kh, "| Klappen-Unterkante Z", round(Z_KL_U, 2))
# Zuordnung Knopf -> Achse -> Klappe(n); Klappennamen aus den Texten der Mechanikzeichnung
texte = [(t.dxf.insert.x - MX, t.dxf.insert.y - MY, t.plain_text().strip()) for t in MD.query("MTEXT")]
kh_namen_used = []
kl_c = [(n_.rsplit("_", 1)[0], ((s_.BoundingBox().xmin + s_.BoundingBox().xmax) / 2, (s_.BoundingBox().ymin + s_.BoundingBox().ymax) / 2)) for n_, s_ in klappen if n_.endswith("Leder")]
zuordnung = {}
for (pa_, pk_) in kh_liste:
    ai = min(range(len(achsen)), key=lambda j: fuss(pa_, achsen[j])[1])
    kn_, kc_ = min(kl_c, key=lambda k_: math.dist(k_[1], pk_))
    tx_ = min(texte, key=lambda t_: math.dist((t_[0], t_[1]), kc_))
    zuordnung.setdefault(ai, []).append((kn_, tx_[2] if math.dist((tx_[0], tx_[1]), kc_) < 30 else "?", round(math.dist(kc_, pk_), 1), kh_namen[len(kh_namen_used)]))
    kh_namen_used.append(1)
# Klappe an zwei Achsen (alter, in der Zeichnung nicht geloeschter Hebel): dort entfernen, wo die Achse noch eine andere Klappe hat
weg = []
for ai, lst_ in zuordnung.items():
    if len(lst_) < 2: continue
    for ent in list(lst_):
        if any(ent[0] in [x_[0] for x_ in l2] for a2, l2 in zuordnung.items() if a2 != ai and len(l2) == 1):
            lst_.remove(ent); weg.append(ent[3])
mech[:] = [m_ for m_ in mech if m_[0] not in weg]
print("Doppelte Klappenhebel entfernt:", weg)
# Achsen mit Knopf, aber ohne Klappenhebel in der Zeichnung: naechste Klappe ohne Hebel ergaenzen
belegt = {k_[0] for v_ in zuordnung.values() for k_ in v_}
for ai in sorted(stange_achse):
    if ai in zuordnung: continue
    frei = [k_ for k_ in kl_c if k_[0] not in belegt]
    if not frei: continue
    kn_, kc_ = min(frei, key=lambda k_: fuss(k_[1], achsen[ai])[1])
    F_, d_ = fuss(kc_, achsen[ai])
    if d_ > 40: continue
    n_kh += 1
    mech.append((f"Klappenhebel_{n_kh}_ergaenzt", stab((kc_[0], F_[1], Z_ACHSE + 1.0), (kc_[0], F_[1], Z_KH), 2.0).fuse(stab((kc_[0], F_[1], Z_KH), (kc_[0], kc_[1], Z_KH), 2.0))))
    tx_ = min(texte, key=lambda t_: math.dist((t_[0], t_[1]), kc_))
    zuordnung[ai] = [(kn_, tx_[2], "ERGAENZT")]; belegt.add(kn_)
for ai in sorted(set(stange_achse) | set(zuordnung)):
    print("ZUORDNUNG Achse", ai + 1, "Knopf", [(r_, round(L - x_)) for r_, x_ in stange_achse.get(ai, [])], "-> Klappe", [(k_[0], k_[1], k_[2]) for k_ in zuordnung.get(ai, [])] or "KEINE")
# --- Rueckdruckfedern: je Achse eine Schenkelfeder, 3 Windungen um die Achse neben dem Arm,
#     ein Schenkel mit Haken um den Gabelgelenk-Arm (drueckt ihn nach vorne zurueck), der andere stuetzt sich am Boden ab ---
D_DRAHT, R_WIND, N_WIND, P_WIND = 0.5, 1.0 + 0.2 + 0.25, 3, 0.6
def draht(punkte):
    w_ = None
    for p_a, p_b in zip(punkte[:-1], punkte[1:]):
        s_ = stab(p_a, p_b, D_DRAHT).fuse(cq.Solid.makeSphere(D_DRAHT / 2, cq.Vector(*p_b)))
        w_ = s_ if w_ is None else w_.fuse(s_)
    return w_
n_feder = 0
for ai, al in sorted(arme.items()):
    (a0, a1) = achsen[ai]; gx, gy = al[0]
    Ln_ = math.dist(a0, a1); dx_, dy_ = (a1[0] - a0[0]) / Ln_, (a1[1] - a0[1]) / Ln_
    nx_, ny_ = -dy_, dx_
    if ny_ > 0: nx_, ny_ = -nx_, -ny_                             # n zeigt nach hinten (hinter den Arm)
    F_, _ = fuss((gx, gy), achsen[ai]); t_arm = (F_[0] - a0[0]) * dx_ + (F_[1] - a0[1]) * dy_
    L_SP = N_WIND * P_WIND
    s_dir = 1.0 if t_arm + 1.2 + L_SP + 5.0 < Ln_ else -1.0          # Seite mit Platz zum naechsten Bock
    t0 = t_arm + s_dir * 1.15                                      # Federanfang neben dem Arm
    hel = cq.Wire.makeHelix(P_WIND, L_SP, R_WIND)
    sp = cq.Workplane("XZ").center(R_WIND, 0).circle(D_DRAHT / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90)     # Wendelachse Z -> X
    if s_dir < 0: sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180)
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), math.degrees(math.atan2(dy_, dx_))).translate(cq.Vector(a0[0] + dx_ * t0, a0[1] + dy_ * t0, Z_ACHSE))
    pA = (a0[0] + dx_ * t0, a0[1] + dy_ * t0)                      # Federende am Arm
    pE = (a0[0] + dx_ * (t0 + s_dir * L_SP), a0[1] + dy_ * (t0 + s_dir * L_SP))
    z_h = Z_ACHSE - 7.0; ab = 1.5 + 0.35                          # Arm ist 3 mm tief -> 1.5 + Luft
    s1 = draht([(pA[0], pA[1], Z_ACHSE - R_WIND),
                (pA[0] + nx_ * ab, pA[1] + ny_ * ab, Z_ACHSE - R_WIND - 1.0),
                (pA[0] + nx_ * ab, pA[1] + ny_ * ab, z_h),
                (F_[0] - dx_ * s_dir * 1.1 + nx_ * ab, F_[1] - dy_ * s_dir * 1.1 + ny_ * ab, z_h),
                (F_[0] - dx_ * s_dir * 1.1 - nx_ * 0.6, F_[1] - dy_ * s_dir * 1.1 - ny_ * 0.6, z_h)])
    s2 = draht([(pE[0], pE[1], Z_ACHSE + R_WIND),
                (pE[0] + nx_ * 4.0, pE[1] + ny_ * 4.0, Z_UB - 0.3)])
    mech.append((f"Feder_Achse_{ai+1}", sp.fuse(s1).fuse(s2))); n_feder += 1
print("Rueckdruckfedern", n_feder)
# --- Luftklappe innen auf dem Zwischendeckel (Balgseite): Achse mit zwei Armen zur Klappe, zwei Lagerboecke
#     (so weit aussen wie moeglich), zwei Schenkelfedern druecken zu. Die Lufttaste (Hub 10 mm) drueckt ueber eine
#     Druckstange die Klappe auf. Dargestellt: Taste gedrueckt -> Klappe an der Druckstange um 10 mm angehoben.
lx0_, lx1_, ly0_, ly1_ = LK_BB
RAND_LK, HUB_LT = 4.5, 10.0
kl_x0, kl_x1, kl_y0, kl_y1 = lx0_ - RAND_LK, lx1_ + RAND_LK, ly0_ - RAND_LK, ly1_ + RAND_LK
T_KL = 6.5
Z_AX_LK = Z1 + 2.5                                                   # Lagerboecke stehen auf dem Zwischendeckel
Y_AX_LK = kl_y0 - 4.0                                                # Achse hinter der Klappe (Scharnierseite)
x_lb0, x_lb1 = gx0 + 4.3 + 0.5, gx1 - 4.3 - 0.5                      # Boecke ganz aussen am Zwischendeckel
x_ds = (max(lx0_, LT_CX - LT_L / 2 + 6) + min(lx1_, LT_CX + LT_L / 2 - 6)) / 2
y_ds_oben = (ly0_ + ly1_) / 2
r_kontakt = y_ds_oben - Y_AX_LK
_dz0 = Z1 - Z_AX_LK                                                   # Klappen-Unterseite relativ zur Achse
_rr = math.hypot(r_kontakt, _dz0); _ph = math.atan2(_dz0, r_kontakt)
winkel = math.degrees(math.asin(min(0.99, (_dz0 + HUB_LT) / _rr)) - _ph)   # Unterseite an der Druckstange um HUB hoch
def um_achse(s_):
    return s_.rotate(cq.Vector(0, Y_AX_LK, Z_AX_LK), cq.Vector(1, Y_AX_LK, Z_AX_LK), winkel)
z_ = Z1
for nm_, dk_ in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
    mech.append((f"Luftklappe_{nm_}", um_achse(box(kl_x0, kl_x1, kl_y0, kl_y1, z_, z_ + dk_)))); z_ += dk_
z_kl_oben = Z1 + T_KL
mech.append(("Luftklappe_Achse", stab((x_lb0 + 1.0, Y_AX_LK, Z_AX_LK), (x_lb1 - 1.0, Y_AX_LK, Z_AX_LK), D_ACHSE)))   # Achse steckt bis zur Schraube in den Boecken
for k_, xa in enumerate((kl_x0 + 8.0, kl_x1 - 8.0)):                 # zwei Arme Achse -> Klappenoberseite
    arm_ = box(xa - 0.75, xa + 0.75, Y_AX_LK - 1.5, kl_y0 + 15.0, z_kl_oben, z_kl_oben + 1.5).fuse(
           box(xa - 0.75, xa + 0.75, Y_AX_LK - 1.5, Y_AX_LK + 1.5, Z_AX_LK - 1.0, z_kl_oben + 1.5))
    mech.append((f"Luftklappe_Arm_{k_+1}", um_achse(arm_)))
for k_, xb in enumerate((x_lb0, x_lb1)):
    lb_, sc_ = bock_solid(5.2, [0.0])
    lb_ = lb_.mirror("XY", (0, 0, (Z_UB - H_LB + Z_UB) / 2)).translate(cq.Vector(0, 0, Z1 - (Z_UB - H_LB)))   # stehend auf dem Deckel
    sc_ = sc_.mirror("XY", (0, 0, (Z_UB - H_LB + Z_UB) / 2)).translate(cq.Vector(0, 0, Z1 - (Z_UB - H_LB)))
    lb_ = lb_.translate(cq.Vector(xb, Y_AX_LK, 0)); sc_ = sc_.translate(cq.Vector(xb, Y_AX_LK, 0))
    mech.append((f"Luftklappe_Lagerbock_{k_+1}", lb_)); mech.append((f"Luftklappe_Lagerbock_Schraube_{k_+1}", sc_))
    # Schenkelfeder neben dem Bock: 3 Windungen um die Achse, ein Schenkel auf dem Deckel, einer auf der Klappe
    s_dir = 1.0 if k_ == 0 else -1.0
    t0 = xb + s_dir * (4.3 + 0.4)
    hel = cq.Wire.makeHelix(P_WIND, N_WIND * P_WIND, R_WIND)
    sp = cq.Workplane("XZ").center(R_WIND, 0).circle(D_DRAHT / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90)
    if s_dir < 0: sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180)
    sp = sp.translate(cq.Vector(t0, Y_AX_LK, Z_AX_LK))
    xe = t0 + s_dir * N_WIND * P_WIND
    s1 = draht([(t0, Y_AX_LK, Z_AX_LK - R_WIND), (t0, Y_AX_LK - 6.0, Z1 + D_DRAHT / 2)])
    s2 = draht([(xe, Y_AX_LK, Z_AX_LK + R_WIND), (xe, Y_AX_LK + 1.0, z_kl_oben + 1.0), (xe, kl_y0 + 12.0, z_kl_oben + D_DRAHT / 2 + 0.05)])
    mech.append((f"Luftklappe_Feder_{k_+1}", sp.fuse(s1).fuse(um_achse(s2)) if False else sp.fuse(s1).fuse(um_achse(s2))))
# Druckstange Ø3 von der Tastenoberseite (gedrueckt) bis an die angehobene Klappe
mech.append(("Luftklappe_Druckstange", stab((x_ds, LT_CY - 2.0, H_TASTE), (x_ds, y_ds_oben, Z1 + HUB_LT), 3.0)))
print("Luftklappe: Achse Y", round(Y_AX_LK, 1), "Z", round(Z_AX_LK, 1), "| Boecke X", round(x_lb0, 1), round(x_lb1, 1),
      "| Taste Hub", HUB_LT, "-> Klappe", round(winkel, 1), "Grad offen")
_lk = [(n_, s_) for n_, s_ in mech if n_.startswith("Luftklappe")]
for i_ in range(len(_lk)):
    for j_ in range(i_ + 1, len(_lk)):
        v_ = _lk[i_][1].intersect(_lk[j_][1]).Volume()
        if v_ > 0.3: print("LK intern", _lk[i_][0], _lk[j_][0], round(v_, 1))
# --- Freischnitte: wo eine Schubstange Boden oder Kammerteile schneidet, dort 1 mm um die Stange ausnehmen ---
frei_n = []
for n1, s1 in mech:
    if not n1.startswith("Schubstange"): continue
    pf_ = stange_pfad[n1]; dick = None
    for p_a, p_b in zip(pf_[:-1], pf_[1:]):
        sg_ = stab(p_a, p_b, D_STANGE + 2.0).fuse(cq.Solid.makeSphere(D_STANGE / 2 + 1.0, cq.Vector(*p_b)))
        dick = sg_ if dick is None else dick.fuse(sg_)
    if boden.intersect(s1).Volume() > 0.05:
        boden = boden.cut(dick); frei_n.append((n1, "Bassboden"))
    for k_ in range(len(teile)):
        n2, s2 = teile[k_]
        if n2.startswith("Lufttaste"): continue
        if s2.intersect(s1).Volume() > 0.05:
            teile[k_] = (n2, s2.cut(dick)); frei_n.append((n1, n2))
print("Freischnitte:", frei_n)
# --- Kollisionen Mechanik gegen Gehaeuse/Kammern ---
gehaeuse = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + hoerner + fuesse + klappen
kol = {}
for n1, s1 in mech:
    bb1 = s1.BoundingBox()
    for n2, s2 in gehaeuse:
        bb2 = s2.BoundingBox()
        if bb1.xmax < bb2.xmin or bb2.xmax < bb1.xmin or bb1.ymax < bb2.ymin or bb2.ymax < bb1.ymin or bb1.zmax < bb2.zmin or bb2.zmax < bb1.zmin: continue
        v = s1.intersect(s2).Volume()
        if v > 0.3 and not n1.startswith(("Lagerbock_Schraube", "Luftklappe_Lagerbock_Schraube")): kol.setdefault(n1, []).append((n2, round(v, 1)))
for k, v in sorted(kol.items()): print("Kollision", k, v[:5])
