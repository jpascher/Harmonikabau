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
D_ACHSE, D_STANGE = 2.4, 2.4             # wie Standard / Variante 3: Achsen, Schenkel, Druckstangen, Hebel Oe2.4
GK_H, GK_L, GK_R, L_GABEL = 3.19, 13.0, 2.0, 8.0   # Fusion-Gelenkkopf (Querbohrung, Schaft, Radius), Gabel gekuerzt
Z_GG = T_DECKEL_BASS + 0.3 + GK_R + 1.5   # Gelenk so tief wie moeglich (1.5 ueber Gaze/Deckel): langer Arm, kleiner Klappenschwenk
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
# --- Buchsen nicht hinter einem Gelenkarm und nicht im Federbereich neben dem Arm: entlang ihrer Achse verschieben ---
_arm_pl = []
for e in MD.query('*[layer=="GB-Schubstangen" | layer=="Begleiter-Schubstangen"]'):
    pts_ = [(q.x - MX, q.y - MY) for q in ezdxf.path.make_path(e).flattening(0.1)]
    h_ = min(pts_, key=lambda p: p[1]); h_ = (h_[0], vy(h_[1]))
    aj_ = min(range(len(achsen)), key=lambda j: fuss(h_, achsen[j])[1])
    _arm_pl.append((h_[0], fuss(h_, achsen[aj_])[0][1], aj_))
def _achse_von(b_):
    # Achse, auf deren Linie die Buchse liegt (Buchse am Achsende darf bis 15 mm ueber das Ende hinaus liegen)
    def _m(j):
        (p0_, p1_) = achsen[j]; l_ = math.dist(p0_, p1_); ux_, uy_ = (p1_[0] - p0_[0]) / l_, (p1_[1] - p0_[1]) / l_
        t_ = (b_[0] - p0_[0]) * ux_ + (b_[1] - p0_[1]) * uy_; q_ = abs((b_[0] - p0_[0]) * uy_ - (b_[1] - p0_[1]) * ux_)
        return q_ + (0 if -15 <= t_ <= l_ + 15 else 100)
    return min(range(len(achsen)), key=_m)
def _stoert(bx_, by_, aj_):
    for (gx_, gy_, ai_) in _arm_pl:
        if ai_ == aj_:
            if abs(bx_ - gx_) < 4.3 + D_STANGE / 2 + 0.4 + 3 * 0.6 + 0.3: return True       # Feder zwischen Arm und eigener Buchse
        elif 0.0 < gy_ - by_ < 7.0 and abs(bx_ - gx_) < 4.3 + GK_R + 1.5: return True     # hinter dem Arm
        elif abs(gy_ - by_) <= 4.5 and abs(bx_ - gx_) < 4.3 + D_STANGE / 2 + 0.4 + 3 * 0.6 + 0.5: return True   # Feder der Nachbarachse
    return False
n_bv = 0
for bi_, b_ in enumerate(lbs):
    aj_ = _achse_von(b_)
    if not _stoert(b_[0], b_[1], aj_): continue
    (a0_, a1_) = achsen[aj_]; Ln_ = math.dist(a0_, a1_); ux_, uy_ = (a1_[0] - a0_[0]) / Ln_, (a1_[1] - a0_[1]) / Ln_
    t_ = (b_[0] - a0_[0]) * ux_ + (b_[1] - a0_[1]) * uy_
    arm_t = [ (g_[0] - a0_[0]) * ux_ + (g_[1] - a0_[1]) * uy_ for g_ in _arm_pl if g_[2] == aj_]
    rich_ = [-1] if t_ < 15 else [1] if t_ > Ln_ - 15 else [-1, 1]
    # gemeinsame Buchse zweier Achsen in einer Linie (Achsen stecken von beiden Seiten): Partnerende mitziehen
    partner_ = [(j_, k_) for j_, (p0_, p1_) in enumerate(achsen) if j_ != aj_ for k_, pe_ in enumerate((p0_, p1_)) if math.dist(pe_, b_[:2]) < 2.5]
    if partner_: rich_ = [-1, 1]
    best_ = None
    for st_ in [0.5 * k for k in range(1, 31)]:
        for sg_ in rich_:
            tn_ = t_ + sg_ * st_
            if tn_ > 15 and tn_ < Ln_ - 15 and any(abs(tn_ - ta_) < 9 for ta_ in arm_t): continue
            if (t_ < 15 and any(tn_ > ta_ - 8.5 for ta_ in arm_t)) or (t_ > Ln_ - 15 and any(tn_ < ta_ + 8.5 for ta_ in arm_t)): continue
            nx_, ny_ = a0_[0] + ux_ * tn_, a0_[1] + uy_ * tn_
            if _stoert(nx_, ny_, aj_): continue
            ok_p = True
            for (j_, k_) in partner_:
                arm_j = [g_[0] for g_ in _arm_pl if g_[2] == j_]
                if any(abs(nx_ - gx_) < 8.5 for gx_ in arm_j) or _stoert(nx_, ny_, j_): ok_p = False
                pe_ = achsen[j_][1 - k_]
                if (nx_ - pe_[0]) * (b_[0] - pe_[0]) <= 0: ok_p = False                    # Partner darf nicht umklappen
            if not ok_p: continue
            # nicht in die Spanne einer anderen Achse auf derselben Linie hinein
            if any(j_ != aj_ and j_ not in [q_[0] for q_ in partner_] and fuss((nx_, ny_), achsen[j_])[1] < 2.0 and min(achsen[j_][0][0], achsen[j_][1][0]) - 4.3 < nx_ < max(achsen[j_][0][0], achsen[j_][1][0]) + 4.3 for j_ in range(len(achsen))): continue
            if any(k_ != bi_ and math.dist((nx_, ny_), q_[:2]) < 6.0 and abs(q_[1] - ny_) > 1.0 for k_, q_ in enumerate(lbs)): continue
            best_ = (tn_, nx_, ny_); break
        if best_: break
    if best_ is None: print("WARNUNG Buchse nicht frei", round(b_[0], 1), round(b_[1], 1)); continue
    tn_, nx_, ny_ = best_
    lbs[bi_] = (nx_, ny_, b_[2]); n_bv += 1
    for (j_, k_) in partner_:
        achsen[j_] = ((nx_, ny_), achsen[j_][1]) if k_ == 0 else (achsen[j_][0], (nx_, ny_))
    if t_ < 15: achsen[aj_] = ((nx_, ny_), a1_)                        # Buchse am Achsende: Ende mitnehmen (verlaengern oder kuerzen)
    elif t_ > Ln_ - 15: achsen[aj_] = (a0_, (nx_, ny_))
    print("Buchse verschoben X", round(b_[0], 1), "->", round(nx_, 1), "Y", round(b_[1], 1), "(Achse", aj_ + 1, ")")
print("Buchsen verschoben:", n_bv)
# Achsen in einer Linie duerfen sich nicht ueberlappen: Enden auf die gemeinsame Buchse setzen
for i_ in range(len(achsen)):
    for j_ in range(i_ + 1, len(achsen)):
        (p0_, p1_), (q0_, q1_) = achsen[i_], achsen[j_]
        _lp = lambda q_: abs((q_[0] - p0_[0]) * (p1_[1] - p0_[1]) - (q_[1] - p0_[1]) * (p1_[0] - p0_[0])) / math.dist(p0_, p1_)
        if max(_lp(q0_), _lp(q1_)) > 1.5: continue
        ilo_, ihi_ = sorted((p0_[0], p1_[0])); jlo_, jhi_ = sorted((q0_[0], q1_[0]))
        if ihi_ <= jlo_ + 0.01 or jhi_ <= ilo_ + 0.01: continue
        kand_ = [b_[0] for b_ in lbs if max(ilo_, jlo_) - 0.5 <= b_[0] <= min(ihi_, jhi_) + 0.5 and abs(b_[1] - (p0_[1] + q0_[1]) / 2) < 3.0]
        m_ = kand_[0] if kand_ else (max(ilo_, jlo_) + min(ihi_, jhi_)) / 2
        def _setz(ax_, lo_, hi_, m_, links_):
            a0_, a1_ = ax_; ym_ = fuss((m_, a0_[1]), ax_)[0][1]
            if links_: return ((m_, ym_), a1_) if a0_[0] < a1_[0] else (a0_, (m_, ym_))
            return (a0_, (m_, ym_)) if a0_[0] < a1_[0] else ((m_, ym_), a1_)
        if ilo_ < jlo_: achsen[i_] = _setz(achsen[i_], ilo_, ihi_, m_, False); achsen[j_] = _setz(achsen[j_], jlo_, jhi_, m_, True)
        else: achsen[j_] = _setz(achsen[j_], jlo_, jhi_, m_, False); achsen[i_] = _setz(achsen[i_], ilo_, ihi_, m_, True)
        print("Achsen", i_ + 1, "und", j_ + 1, "treffen sich in der Buchse bei X", round(m_, 1))
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
# --- Achsenden in die Buchsenbohrung: Achse bis zur Bohrungsmitte verlaengern/kuerzen; einfache Buchse seitlich auf die
#     Achslinie ruecken (bei Doppelbuchsen das Achsende auf die Bohrung setzen); Ende ohne Buchse bekommt eine eigene ---
def _bohrungen():
    bl_ = []
    for i_, (cx_, cy_, rot_) in enumerate(lbs):
        if i_ not in in_paar: bl_.append(("E", i_, cx_, cy_, rot_))
    for i_, j_ in paare:
        (x1_, y1_, _), (x2_, y2_, _) = lbs[i_], lbs[j_]
        mx_, my_ = (x1_ + x2_) / 2, (y1_ + y2_) / 2; rot_ = math.degrees(math.atan2(y2_ - y1_, x2_ - x1_)) - 90.0
        nx_, ny_ = -math.sin(math.radians(rot_)), math.cos(math.radians(rot_))
        for sg_ in (-1, 1): bl_.append(("D", i_, mx_ + sg_ * nx_ * A_DOPPEL / 2, my_ + sg_ * ny_ * A_DOPPEL / 2, rot_))
    return bl_
# Doppelbuchse: beide Bohrungen auf die Linien ihrer zwei Achsen legen (auch bei schraegen Achsen), Abstand 3.2
def _y_auf(j_, x_):
    (p0_, p1_) = achsen[j_]
    return p0_[1] + (p1_[1] - p0_[1]) * (x_ - p0_[0]) / (p1_[0] - p0_[0]) if abs(p1_[0] - p0_[0]) > 1e-6 else p0_[1]
for i_, j_ in paare:
    ai_a, ai_b = _achse_von(lbs[i_]), _achse_von(lbs[j_])
    if ai_a == ai_b: continue
    x_ = (lbs[i_][0] + lbs[j_][0]) / 2; ya_, yb_ = _y_auf(ai_a, x_), _y_auf(ai_b, x_); m_ = (ya_ + yb_) / 2
    sg_ = 1 if ya_ > yb_ else -1
    lbs[i_] = (x_, m_ + sg_ * A_DOPPEL / 2, lbs[i_][2]); lbs[j_] = (x_, m_ - sg_ * A_DOPPEL / 2, lbs[j_][2])
n_fix = [0, 0, 0]; belegt_b = {}
for ai_ in range(len(achsen)):
    for k_ in (0, 1):
        a0_, a1_ = achsen[ai_]; l_ = math.dist(a0_, a1_); ux_, uy_ = (a1_[0] - a0_[0]) / l_, (a1_[1] - a0_[1]) / l_
        e_ = a0_ if k_ == 0 else a1_; aus_ = -1 if k_ == 0 else 1                   # Richtung nach aussen
        best_ = None
        for (art_, bi_, bx_, by_, br_) in _bohrungen():
            t_ = (bx_ - e_[0]) * ux_ + (by_ - e_[1]) * uy_; q_ = (bx_ - e_[0]) * (-uy_) + (by_ - e_[1]) * ux_
            if art_ == "D" and any(ux_ * aus_ * vx_ + uy_ * aus_ * vy_ > 0 for vx_, vy_ in belegt_b.get((round(bx_, 2), round(by_, 2)), [])): continue   # je Bohrung ein Achsende von jeder Seite
            if abs(q_) < (6.5 if art_ == "D" else 3.0) and -15.0 < t_ * aus_ < 15.0 and abs(math.sin(math.radians(br_) - math.atan2(uy_, ux_))) < 0.3:
                if best_ is None or abs(t_) + 5 * abs(q_) < abs(best_[2]) + 5 * abs(best_[3]): best_ = (art_, bi_, t_, q_, bx_, by_)   # quer zaehlt staerker: in der Doppelbuchse die richtige Bohrung
        if best_ is None:
            lbs.append((e_[0], e_[1], math.degrees(math.atan2(uy_, ux_)))); n_fix[2] += 1
            print("Buchse ergaenzt am Ende von Achse", ai_ + 1, [round(v_, 1) for v_ in e_]); continue
        art_, bi_, t_, q_, bx_, by_ = best_
        if art_ == "D": belegt_b.setdefault((round(bx_, 2), round(by_, 2)), []).append((ux_ * aus_, uy_ * aus_))
        if abs(q_) > 0.1:
            geteilt_ = any(j_ != ai_ and math.dist(pe_, (bx_, by_)) < 4.0 for j_, (p0_, p1_) in enumerate(achsen) for pe_ in (p0_, p1_))
            if art_ == "E" and geteilt_: art_ = "D"                        # gemeinsame Buchse zweier Achsen: Achsende auf die Bohrung
            if art_ == "E":                                              # einfache Buchse seitlich auf die Achslinie
                cx_, cy_, rot_ = lbs[bi_]; lbs[bi_] = (cx_ + uy_ * q_, cy_ - ux_ * q_, rot_); n_fix[1] += 1
                bx_, by_ = bx_ + uy_ * q_, by_ - ux_ * q_
            else: e_n = (bx_, by_)
        ne_ = (e_[0] + ux_ * t_, e_[1] + uy_ * t_) if (art_ == "E" or abs(q_) <= 0.1) else (bx_, by_)
        if math.dist(ne_, e_) > 0.05:
            achsen[ai_] = (ne_, a1_) if k_ == 0 else (a0_, ne_); n_fix[0] += 1
for _ai in (4, 17): print("DBG Achse", _ai + 1, [[round(v_, 2) for v_ in p_] for p_ in achsen[_ai]])
print("DBG Bohrungen", [(b_[0], b_[1], round(b_[2], 1), round(b_[3], 1)) for b_ in _bohrungen() if b_[0] == "D" and (b_[2] < 150 or b_[2] > 350)])
print("Achsenden in die Buchse gesetzt:", n_fix[0], "| Buchsen seitlich gerueckt:", n_fix[1], "| Buchsen ergaenzt:", n_fix[2])
for i, (a0, a1) in enumerate(achsen):
    _l = math.dist(a0, a1); _ux, _uy = (a1[0] - a0[0]) / _l, (a1[1] - a0[1]) / _l
    mech.append((f"Achse_{i+1}", stab((a0[0] + _ux, a0[1] + _uy, Z_ACHSE), (a1[0] - _ux, a1[1] - _uy, Z_ACHSE), D_ACHSE)))   # Enden 1 mm vor der Schraube
_LBD = ezdxf.readfile("lagerbeide.dxf").modelspace()
def _prof(lay):
    e_ = [e for e in _LBD.query(f'LWPOLYLINE[layer=="{lay}"]') if e.closed][0]
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e_).flattening(0.02)]
    return pp_[:-1] if math.dist(pp_[0], pp_[-1]) < 1e-6 else pp_
PROF_E, PROF_D = _prof("Layer1"), _prof("0")
def bock_solid(breite, bohrungen, g_=1):
    # Lagerbuchse nach lagerbeide.dxf: einfach 5.2 breit mit Hohlkehle R1.2 auf der Seite g_ (Nachbarachse), doppelt 8.4 breit
    pp_, mi_ = (PROF_E, 2.6) if breite < 6 else (PROF_D, 4.2)
    yz_ = [(-g_ * (-q[1] - mi_), Z_UB - (q[0] + 13.6)) for q in pp_]
    lb_ = cq.Workplane("YZ").polyline(yz_).close().extrude(8.6).val().translate(cq.Vector(-4.3, 0, 0)); sc_ = None
    for yb in bohrungen:
        lb_ = lb_.cut(stab((-5, yb, Z_ACHSE), (5, yb, Z_ACHSE), D_ACHSE + 0.1)).cut(stab((0, yb, Z_UB - H_LB - 1), (0, yb, Z_UB + 1), 1.5))
        s_ = stab((0, yb, Z_UB - H_LB - 0.9), (0, yb, Z_UB - H_LB + 10), 1.4).fuse(stab((0, yb, Z_UB - H_LB - 0.9), (0, yb, Z_UB - H_LB), 2.5))
        sc_ = s_ if sc_ is None else sc_.fuse(s_)
    return lb_, sc_
nb = 0
def kehle(cx, cy, rot):
    nx_, ny_ = -math.sin(math.radians(rot)), math.cos(math.radians(rot)); best_ = (99, 1)
    for (a0, a1) in achsen:
        f_, d_ = fuss((cx, cy), (a0, a1))
        if 0.5 < d_ < best_[0]: best_ = (d_, 1 if (f_[0] - cx) * nx_ + (f_[1] - cy) * ny_ > 0 else -1)
    return best_[1]
for i, (cx, cy, rot) in enumerate(lbs):
    if i in in_paar: continue
    lb, sc = bock_solid(5.2, [0.0], kehle(cx, cy, rot))
    tr = lambda s: s.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), rot).translate(cq.Vector(cx, cy, 0))
    nb += 1; mech.append((f"Lagerbock_{nb}", tr(lb))); mech.append((f"Lagerbock_Schraube_{nb}", tr(sc)))
for i, j in paare:
    (x1, y1, _), (x2, y2, _) = lbs[i], lbs[j]
    mx_, my_ = (x1 + x2) / 2, (y1 + y2) / 2
    rot = math.degrees(math.atan2(y2 - y1, x2 - x1)) - 90.0          # Achsrichtung senkrecht zur Paar-Verbindung
    lb, sc = bock_solid(8.4, [-A_DOPPEL / 2, A_DOPPEL / 2])
    tr = lambda s: s.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), rot).translate(cq.Vector(mx_, my_, 0))
    nb += 1; mech.append((f"Lagerbock_doppelt_{nb}", tr(lb)))
    for k_, sv_ in enumerate(tr(sc).Solids()): mech.append((f"Lagerbock_Schraube_{nb}_{k_+1}", sv_))
bohr_ = []
for i_, (cx_, cy_, rot_) in enumerate(lbs):
    if i_ in in_paar: continue
    bohr_.append(("E%d" % i_, cx_, cy_, rot_))
for i_, j_ in paare:
    (x1_, y1_, _), (x2_, y2_, _) = lbs[i_], lbs[j_]
    mx_, my_ = (x1_ + x2_) / 2, (y1_ + y2_) / 2; rot_ = math.degrees(math.atan2(y2_ - y1_, x2_ - x1_)) - 90.0
    nx_, ny_ = -math.sin(math.radians(rot_)), math.cos(math.radians(rot_))
    for sg_ in (-1, 1): bohr_.append(("D%d" % i_, mx_ + sg_ * nx_ * A_DOPPEL / 2, my_ + sg_ * ny_ * A_DOPPEL / 2, rot_))
for ai_, (a0_, a1_) in enumerate(achsen):
    l_ = math.dist(a0_, a1_); ux_, uy_ = (a1_[0] - a0_[0]) / l_, (a1_[1] - a0_[1]) / l_
    for k_, e_ in enumerate(((a0_[0] + ux_, a0_[1] + uy_), (a1_[0] - ux_, a1_[1] - uy_))):
        best_ = None
        for (nm_, bx_, by_, br_) in bohr_:
            vx_, vy_ = math.cos(math.radians(br_)), math.sin(math.radians(br_))
            along_ = (e_[0] - bx_) * vx_ + (e_[1] - by_) * vy_; quer_ = abs((e_[0] - bx_) * vy_ - (e_[1] - by_) * vx_)
            d_ = math.dist(e_, (bx_, by_))
            if best_ is None or d_ < best_[0]: best_ = (d_, nm_, along_, quer_)
        ok_ = best_[3] < 0.3 and abs(best_[2]) <= 4.3 - 0.5
        print("ACHSENDE", ai_ + 1, "Ende", k_, "OK" if ok_ else "FEHLT", "Buchse", best_[1], "laengs", round(best_[2], 2), "quer", round(best_[3], 2))
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
arm_info = []
# Gelenk wie Standard-Bassteil: Fusion "gelenkkopf" an der Druckstange, Gabel (8 mm) am Schenkeldraht, Bolzen Oe1.4
import os as _os
_gkp = "kb_gelenkkopf.stp" if _os.path.exists("kb_gelenkkopf.stp") else "/tmp/bass/v0/kb_gelenkkopf.stp"
_gk = cq.importers.importStep(_gkp).solids().vals()
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
arme = {}
for i, (P, r, Xk) in enumerate(stangen):
    ai = min(range(len(achsen)), key=lambda j: fuss(P, achsen[j])[1]); F, dist = fuss(P, achsen[ai])
    Zr, Yk = REIHEN[r]
    gx, gy = P[0], F[1]
    stange_achse.setdefault(ai, []).append((r, Xk))
    arme.setdefault(ai, []).append((gx, gy))
    # Gabelgelenk zweiteilig: Gabelkopf (an der Schubstange) + Gelenkbolzen; Drehpunkt = Bolzenachse am Armende
    gy = fuss((gx, P[1]), achsen[ai])[0][1]                          # Schenkel genau unter der Achse
    arm_info.append((ai, i + 1, gx, gy, r, Xk))
    ax0_, ax1_ = achsen[ai]
    sk_ = stab((gx, gy, Z_ACHSE), (gx, gy, Z_GG + 4.5), D_STANGE).cut(stab((ax0_[0], ax0_[1], Z_ACHSE), (ax1_[0], ax1_[1], Z_ACHSE), D_ACHSE))
    mech.append((f"Schenkel_Knopf_{i+1}", sk_))
    # Verlauf: unten (Z_GG) bis hinter die letzte Klappe im Weg, 45 Grad hoch bis knapp unter den Boden,
    # dort zum Knopf (X-Versatz), zuletzt 45 Grad hoch auf Knopfhoehe; Lufttasten-Langloch wird unten nicht ueberquert.
    y0 = gy + GK_H + GK_L + 0.5; y_end = Yk - L_ZYL_R[r]
    xs = (min(gx, Xk) - 2.5, max(gx, Xk) + 2.5)
    y_frei = max([b_.ymax + 2.0 for b_ in kl_bb if b_.xmax > xs[0] and b_.xmin < xs[1] and b_.ymax > y0] + [y0])
    z_mid = min(Zr, Z_UNTER)
    if xs[1] > lt_x0 and xs[0] < lt_x1 and y0 < lt_y1:
        z_mid = max(z_mid, H_TASTE + 2.0)                                # Bogen ueber die gedrueckte Lufttaste
        if y_frei > lt_y0 - (z_mid - Z_GG) - 1.0 and z_mid < Z_KL_U - 1.5:
            y_frei = lt_y0 - (z_mid - Z_GG) - 1.0                        # schon unter der Klappe hochgehen (unter Klappen-Unterkante)
    y1 = max(y_frei, y0)
    hind_ = []                                                     # Achsen und Buchsen ueber der Stange
    for (q0_, q1_) in achsen:
        if min(q0_[0], q1_[0]) - 1 < gx < max(q0_[0], q1_[0]) + 1 and abs(q1_[0] - q0_[0]) > 1e-6:
            hind_.append((q0_[1] + (q1_[1] - q0_[1]) * (gx - q0_[0]) / (q1_[0] - q0_[0]), D_ACHSE / 2, Z_ACHSE - D_ACHSE / 2))
    for (bx_, by_, _) in lbs:
        if abs(bx_ - gx) < 4.3 + D_STANGE / 2 + 0.5: hind_.append((by_, 4.3, Z_UB - H_LB - 0.9))
    for (yo_, hb_, zb_) in hind_:
        if yo_ > gy + 1.0 and zb_ - D_STANGE / 2 - 0.3 < z_mid:
            y1 = max(y1, yo_ + hb_ + D_STANGE / 2 + 0.5 + HUB - max(0.0, zb_ - D_STANGE / 2 - 0.3 - Z_GG))
    y2 = y1 + (z_mid - Z_GG)
    y3 = y_end - (Zr - z_mid)
    if y3 < y2: y3 = y2
    pfad = [(gx, gy + 4.6, Z_GG), (gx, y1, Z_GG), (gx, y2, z_mid), (Xk, y3, z_mid), (Xk, y_end, Zr)]
    pfad = [p for k_, p in enumerate(pfad) if k_ == 0 or math.dist(p, pfad[k_ - 1]) > 0.05]
    st = None
    for p_a, p_b in zip(pfad[:-1], pfad[1:]):
        seg = stab(p_a, p_b, D_STANGE).fuse(cq.Solid.makeSphere(D_STANGE / 2, cq.Vector(*p_b)))
        st = seg if st is None else st.fuse(seg)
    mech.append((f"Schubstange_{r}_{Xk:.0f}", st)); stange_pfad[f"Schubstange_{r}_{Xk:.0f}"] = pfad
    gelenk((gx, gy, Z_GG), -1, sk_, st, f"{i+1}")
    if dist < 6: n_ok += 1
print("Mechanik: Achsen", len(achsen), "Lagerboecke", len(lbs), "Knoepfe", len(knoepfe), "Schubstangen", len(stangen), "davon an Achse (<6 mm)", n_ok)
# --- Klappenhebel: Plan-Linien Achse -> Klappenmitte; von der Achse senkrecht nach unten unter die Klappe, dann waagrecht ---
gesehen = set(); n_kh = 0; kh_liste = []; kh_namen = []
def hebel_pfad(pa, pk):
    # von der Achse senkrecht hinunter, waagrecht zur Klappenmitte; liegt eine Buchse ueber dem Weg, tiefer darunter durch
    # und erst mit Abstand (auch beim Zurueckschwenken) wieder auf Klappenhoehe
    aj_ = min(range(len(achsen)), key=lambda j: fuss(pa, achsen[j])[1]); (q0_, q1_) = achsen[aj_]
    lq_ = math.dist(q0_, q1_); vx_, vy_ = (q1_[0] - q0_[0]) / lq_, (q1_[1] - q0_[1]) / lq_
    _frei = lambda p_: all(math.dist(p_, b_[:2]) >= 4.3 + D_STANGE / 2 + 1.2 for b_ in lbs)
    if not _frei(pa):
        for st_ in [0.25 * k for k in range(1, 41)]:
            c_ = [(pa[0] + sg_ * st_ * vx_, pa[1] + sg_ * st_ * vy_) for sg_ in (1, -1)]
            c_ = [p_ for p_ in c_ if _frei(p_)]
            if c_: pa = min(c_, key=lambda p_: math.dist(p_, pk)); break
        else: print("WARNUNG Klappenhebel nahe Buchse", [round(v_, 1) for v_ in pa])
    L_ = math.dist(pa, pk); ux_, uy_ = (pk[0] - pa[0]) / L_, (pk[1] - pa[1]) / L_
    unter_ = []
    for (bx_, by_, _) in lbs:
        t_ = (bx_ - pa[0]) * ux_ + (by_ - pa[1]) * uy_; d_ = abs((bx_ - pa[0]) * uy_ - (by_ - pa[1]) * ux_)
        if 1.5 < t_ < L_ and d_ < 4.3 + D_STANGE / 2 + 0.8: unter_.append(t_)
    if not unter_:
        pts_ = [(pa[0], pa[1], Z_ACHSE), (pa[0], pa[1], Z_KH), (pk[0], pk[1], Z_KH)]
    else:
        z_lo = Z_UB - H_LB - 0.9 - D_STANGE / 2 - 0.15; t_r = min(L_ - 1.0, max(unter_) + 4.3 + 1.2 + 4.5)
        pts_ = [(pa[0], pa[1], Z_ACHSE), (pa[0], pa[1], z_lo), (pa[0] + ux_ * t_r, pa[1] + uy_ * t_r, z_lo),
                (pa[0] + ux_ * min(L_, t_r + (Z_KH - z_lo)), pa[1] + uy_ * min(L_, t_r + (Z_KH - z_lo)), Z_KH), (pk[0], pk[1], Z_KH)]
    pts_ = [q_ for k_, q_ in enumerate(pts_) if k_ == 0 or math.dist(q_, pts_[k_ - 1]) > 0.05]
    w_ = None
    for q0_, q1_ in zip(pts_[:-1], pts_[1:]):
        sg_ = stab(q0_, q1_, D_STANGE).fuse(cq.Solid.makeSphere(D_STANGE / 2 + 0.01, cq.Vector(*q1_)))
        w_ = sg_ if w_ is None else w_.fuse(sg_)
    return w_
Z_KH = Z_KL_U - D_STANGE / 2 - 0.05
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
    kh = hebel_pfad(pa, pk)
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
    mech.append((f"Klappenhebel_{n_kh}_ergaenzt", hebel_pfad((kc_[0], F_[1]), kc_)))
    tx_ = min(texte, key=lambda t_: math.dist((t_[0], t_[1]), kc_))
    zuordnung[ai] = [(kn_, tx_[2], "ERGAENZT", f"Klappenhebel_{n_kh}_ergaenzt")]; belegt.add(kn_)
for ai in sorted(set(stange_achse) | set(zuordnung)):
    print("ZUORDNUNG Achse", ai + 1, "Knopf", [(r_, round(L - x_)) for r_, x_ in stange_achse.get(ai, [])], "-> Klappe", [(k_[0], k_[1], k_[2]) for k_ in zuordnung.get(ai, [])] or "KEINE")
if __import__("os").environ.get("DUMP_M") and False:
    import json as _json
    _nb = []; _k = 0
    for _i, (_cx, _cy, _rot) in enumerate(lbs):
        if _i in in_paar: continue
        _k += 1; _nb.append(("Lagerbock_%d" % _k, _cx, _cy, _rot))
    for _i, _j in paare:
        _k += 1; _nb.append(("Lagerbock_doppelt_%d" % _k, (lbs[_i][0] + lbs[_j][0]) / 2, (lbs[_i][1] + lbs[_j][1]) / 2, 0))
    _kl = {n_: [round(v_, 1) for v_ in (s_.BoundingBox().xmin, s_.BoundingBox().xmax, s_.BoundingBox().ymin, s_.BoundingBox().ymax)] for n_, s_ in klappen if n_.endswith("Leder")}
    _json.dump(dict(achsen=achsen, lbs=lbs, boecke=_nb, arme=arm_info, zuord={str(k_): v_ for k_, v_ in zuordnung.items()}, kh=kh_liste, klappen=_kl, taste=[s_.BoundingBox().__dict__ for n_, s_ in teile if n_ == "Lufttaste"] and str([s_.BoundingBox().ymin for n_, s_ in teile if n_ == "Lufttaste"])), open("/tmp/bass/mdump.json", "w"), default=str)
    raise SystemExit("DUMP")
# --- Rueckdruckfedern: je Achse eine Schenkelfeder, 3 Windungen um die Achse neben dem Arm,
#     ein Schenkel mit Haken um den Gabelgelenk-Arm (drueckt ihn nach vorne zurueck), der andere stuetzt sich am Boden ab ---
D_DRAHT, R_WIND, N_WIND, P_WIND = 0.5, D_ACHSE / 2 + 0.2 + 0.25, 3, 0.6
def draht(punkte):
    w_ = None
    for p_a, p_b in zip(punkte[:-1], punkte[1:]):
        s_ = stab(p_a, p_b, D_DRAHT).fuse(cq.Solid.makeSphere(D_DRAHT / 2 + 0.01, cq.Vector(*p_b))).fuse(cq.Solid.makeSphere(D_DRAHT / 2 + 0.01, cq.Vector(*p_a)))
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
    tb_ = [((b_[0] - a0[0]) * dx_ + (b_[1] - a0[1]) * dy_) for b_ in lbs if abs((b_[0] - a0[0]) * dy_ - (b_[1] - a0[1]) * dx_) < 3.0] or [0.0, Ln_]
    kh_t_ = [((q_[0][0] - a0[0]) * dx_ + (q_[0][1] - a0[1]) * dy_) for q_ in kh_liste if fuss(q_[0], achsen[ai])[1] < 3.0]
    wahl_ = []                                                     # Feder ans freie Ende: zwischen Arm und Buchse kein Klappenhebel
    for (gx2_, gy2_) in al:
        F2_, _ = fuss((gx2_, gy2_), achsen[ai]); t2_ = (F2_[0] - a0[0]) * dx_ + (F2_[1] - a0[1]) * dy_
        for sd_ in (1.0, -1.0):
            tb2_ = [t_ for t_ in tb_ if (t_ - t2_) * sd_ > 0]
            if not tb2_: continue
            tn_ = min(tb2_, key=lambda t_: abs(t_ - t2_))
            if abs(tn_ - t2_) < 8.0: continue
            if any(min(t2_, tn_) < th_ < max(t2_, tn_) for th_ in kh_t_): continue
            if any(min(t2_, tn_) < ta_ < max(t2_, tn_) for ta_ in [((fuss(g_, achsen[ai])[0][0] - a0[0]) * dx_ + (fuss(g_, achsen[ai])[0][1] - a0[1]) * dy_) for g_ in al if g_ != (gx2_, gy2_)]): continue
            wahl_.append((abs(tn_ - t2_), sd_, (gx2_, gy2_), F2_, t2_))
    if wahl_:
        _, s_dir, (gx, gy), F_, t_arm = min(wahl_)
    else:
        t_nah = min(tb_, key=lambda t_: abs(t_ - t_arm))
        s_dir = 1.0 if t_nah > t_arm else -1.0
        print("WARNUNG Feder Achse", ai + 1, "ohne freies Ende")
    t0 = t_arm + s_dir * (D_STANGE / 2 + 0.4)
    hel = cq.Wire.makeHelix(P_WIND, L_SP, R_WIND)
    sp = cq.Workplane("XZ").center(R_WIND, 0).circle(D_DRAHT / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90)     # Wendelachse Z -> X
    if s_dir < 0: sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180)
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), math.degrees(math.atan2(dy_, dx_))).translate(cq.Vector(a0[0] + dx_ * t0, a0[1] + dy_ * t0, Z_ACHSE))
    pA = (a0[0] + dx_ * t0, a0[1] + dy_ * t0)                      # Federende am Arm
    pE = (a0[0] + dx_ * (t0 + s_dir * L_SP), a0[1] + dy_ * (t0 + s_dir * L_SP))
    z_h = Z_GG + 4.5 + 2.5; ab = GK_R + 0.1 + D_DRAHT / 2            # Haken um den Gabelschaft
    s1 = draht([(pA[0], pA[1], Z_ACHSE - R_WIND),
                (pA[0] + nx_ * ab, pA[1] + ny_ * ab, Z_ACHSE - R_WIND - 1.0),
                (pA[0] + nx_ * ab, pA[1] + ny_ * ab, z_h),
                (F_[0] - dx_ * s_dir * ab + nx_ * ab, F_[1] - dy_ * s_dir * ab + ny_ * ab, z_h),
                (F_[0] - dx_ * s_dir * ab - nx_ * 0.8, F_[1] - dy_ * s_dir * ab - ny_ * 0.8, z_h)])
    s2 = draht([(pE[0], pE[1], Z_ACHSE + R_WIND),
                (pE[0] - nx_ * 0.8, pE[1] - ny_ * 0.8, Z_UB - 0.3)])
    mech.append((f"Feder_Achse_{ai+1}", sp.fuse(s1).fuse(s2))); n_feder += 1
    _kh_t = [((q_[0][0] - a0[0]) * dx_ + (q_[0][1] - a0[1]) * dy_) for q_ in kh_liste if fuss(q_[0], achsen[ai])[1] < 3.0]
    print("FEDER Achse", ai + 1, "Klappe", [z_[0] for z_ in zuordnung.get(ai, [])], "Arm t", round(t_arm, 1), "Feder", "+" if s_dir > 0 else "-", "Hebel t", [round(t_, 1) for t_ in _kh_t], "Buchsen t", [round(t_, 1) for t_ in tb_], "Laenge", round(Ln_, 1))
print("Rueckdruckfedern", n_feder)
if __import__("os").environ.get("STOP_DIAG"): raise SystemExit("DIAG")
# --- Luftklappe innen auf dem Zwischendeckel (Balgseite): Achse mit zwei Armen zur Klappe, zwei Lagerboecke
#     (so weit aussen wie moeglich), zwei Schenkelfedern druecken zu. Die Lufttaste (Hub 10 mm) drueckt ueber eine
#     Druckstange die Klappe auf. Dargestellt: Taste gedrueckt -> Klappe an der Druckstange um 10 mm angehoben.
lx0_, lx1_, ly0_, ly1_ = LK_BB
RAND_LK, HUB_LT = 4.5, 10.0
kl_x0, kl_x1, kl_y0, kl_y1 = lx0_ - RAND_LK, lx1_ + RAND_LK, ly0_ - RAND_LK, ly1_ + RAND_LK
T_KL = 6.5
Z_AX_LK = Z1 + 2.5                                                   # Lagerboecke stehen auf dem Zwischendeckel
# Gelenk nicht neben dem Begleiterstimmstock, sondern so weit aussen wie moeglich: an der Frontseite am Rand des
# Zwischendeckels (langer Hebelweg Achse -> Druckpunkt). SG_LK = Richtung Achse -> Klappe (-1: Klappe liegt in -Y).
_zd = [s_ for n_, s_ in teile if n_ == "Zwischen_Deckel"][0].BoundingBox()
Y_AX_LK = _zd.ymax - 2.6 - 0.4                                       # Bock 5.2 breit, 0.4 vom Deckelrand
SG_LK = -1.0
x_lb0, x_lb1 = gx0 + 4.3 + 0.5, gx1 - 4.3 - 0.5                      # Boecke ganz aussen am Zwischendeckel
x_ds = (max(lx0_, LT_CX - LT_L / 2 + 6) + min(lx1_, LT_CX + LT_L / 2 - 6)) / 2
y_ds_oben = (ly0_ + ly1_) / 2
r_kontakt = abs(y_ds_oben - Y_AX_LK)
_dz0 = Z1 - Z_AX_LK                                                   # Klappen-Unterseite relativ zur Achse
_rr = math.hypot(r_kontakt, _dz0); _ph = math.atan2(_dz0, r_kontakt)
winkel = (1 if SG_LK > 0 else -1) * math.degrees(math.asin(min(0.99, (_dz0 + HUB_LT) / _rr)) - _ph)   # Unterseite an der Druckstange um HUB hoch
def um_achse(s_):
    return s_.rotate(cq.Vector(0, Y_AX_LK, Z_AX_LK), cq.Vector(1, Y_AX_LK, Z_AX_LK), winkel)
z_ = Z1
for nm_, dk_ in (("Leder", 1.0), ("Filz", 4.0), ("Epoxy", 1.5)):
    mech.append((f"Luftklappe_{nm_}", um_achse(box(kl_x0, kl_x1, kl_y0, kl_y1, z_, z_ + dk_)))); z_ += dk_
z_kl_oben = Z1 + T_KL
mech.append(("Luftklappe_Achse", stab((x_lb0 + 1.0, Y_AX_LK, Z_AX_LK), (x_lb1 - 1.0, Y_AX_LK, Z_AX_LK), D_ACHSE)))   # Achse steckt bis zur Schraube in den Boecken
for k_, xa in enumerate((kl_x0 + 8.0, kl_x1 - 8.0)):                 # zwei Arme Achse -> Klappenoberseite
    _ya, _yb = sorted((Y_AX_LK - SG_LK * 1.5, (kl_y0 + 15.0) if SG_LK > 0 else (kl_y1 - 15.0)))
    arm_ = box(xa - 0.75, xa + 0.75, _ya, _yb, z_kl_oben, z_kl_oben + 1.5).fuse(
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
    s1 = draht([(t0, Y_AX_LK, Z_AX_LK - R_WIND), (t0, Y_AX_LK - SG_LK * 5.0, Z1 + D_DRAHT / 2)])
    s2 = draht([(xe, Y_AX_LK, Z_AX_LK + R_WIND), (xe, Y_AX_LK + SG_LK * 1.0, z_kl_oben + 1.0), (xe, (kl_y0 + 12.0) if SG_LK > 0 else (kl_y1 - 12.0), z_kl_oben + D_DRAHT / 2 + 0.05)])
    mech.append((f"Luftklappe_Feder_{k_+1}", sp.fuse(s1).fuse(um_achse(s2)) if False else sp.fuse(s1).fuse(um_achse(s2))))
# Druckstange Ø3 von der Tastenoberseite (gedrueckt) bis an die angehobene Klappe
mech.append(("Luftklappe_Druckstange", stab((x_ds, LT_CY - 2.0, H_TASTE), (x_ds, y_ds_oben, Z1 + HUB_LT - 0.4), 3.0)))
print("Luftklappe: Achse Y", round(Y_AX_LK, 1), "Z", round(Z_AX_LK, 1), "| Boecke X", round(x_lb0, 1), round(x_lb1, 1),
      "| Taste Hub", HUB_LT, "-> Klappe", round(winkel, 1), "Grad offen")
_lk = [(n_, s_) for n_, s_ in mech if n_.startswith("Luftklappe")]
for i_ in range(len(_lk)):
    for j_ in range(i_ + 1, len(_lk)):
        v_ = _lk[i_][1].intersect(_lk[j_][1]).Volume()
        if v_ > 0.3: print("LK intern", _lk[i_][0], _lk[j_][0], round(v_, 1))
def _bb_ueber(b1, b2):
    return not (b1.xmax < b2.xmin or b2.xmax < b1.xmin or b1.ymax < b2.ymin or b2.ymax < b1.ymin or b1.zmax < b2.zmin or b2.zmax < b1.zmin)
# --- Lufttaste-Federstange: Drehstab Ø2.4 links (hohes X) an der Seitenwandverstaerkung knapp unter dem Deckel, in zwei
#     Lagerboecken an der Wand; hinteres Ende abgewinkelt in die Wand (verdrehfest). Vorne 90 Grad Richtung Taste, hinter dem
#     Eckblock vorbei bis zum Tastenbeginn, dort nach oben (Richtung Deckel) fest in die Taste eingelassen. Taste gedrueckt dargestellt.
D_FS, Z_FS = 2.4, T_DECKEL_BASS + 3.0                                # Stabmitte 1.8 mm unter dem Deckel (Stab-Oberkante)
_wand = [s_ for n_, s_ in teile if n_ == "Deckelauflage_rechts"][0].BoundingBox()
_tbb = [s_ for n_, s_ in teile if n_ == "Lufttaste"][0].BoundingBox()
X_FS = _wand.xmin - 2.2                                              # Stabachse 1 mm vor der Wandverstaerkung
Y_FS_ARM = min(LT_CY, yi1 - BLK - D_FS / 2 - 1.3)                    # Arm hinter dem vorderen Eckblock vorbei, unter der Taste
X_FS_T = _tbb.xmax - 5.0                                             # Tastenbeginn (Seite zur Wand)
Y_FS_H = Y_FS_ARM - 82.0                                             # wirksame Drehstablaenge ~82 mm
r_fs = X_FS - X_FS_T
w_fs = math.degrees(math.asin(HUB_LT / r_fs))                        # Armdrehung fuer 10 mm Tastenhub
def _fs_dreh(s_):
    return s_.rotate(cq.Vector(X_FS, 0, Z_FS), cq.Vector(X_FS, 1, Z_FS), w_fs)
def _zug(pkte):
    t_ = None
    for a_, b_ in zip(pkte[:-1], pkte[1:]):
        g_ = stab(a_, b_, D_FS)
        t_ = g_ if t_ is None else t_.fuse(g_)
    for q_ in pkte[1:-1]: t_ = t_.fuse(cq.Solid.makeSphere(D_FS / 2, cq.Vector(*q_)))
    return t_
_ecke = (X_FS, Y_FS_ARM, Z_FS)
_stab_ = _zug([(X_FS + 3.0, Y_FS_H, Z_FS), (X_FS, Y_FS_H, Z_FS), _ecke])            # Drehstab + Wandende
# Bassklappen unter dem Arm (oeffnen Richtung Deckel): Arm nach dem Eckblock nach vorne versetzen, vor die Klappen
X_FS_J = xi1 - BLK - 3.5
_kly = [s_.BoundingBox().ymax for n_, s_ in klappen if _bb_ueber(s_.BoundingBox(), box(X_FS_T - 2, X_FS_J + 2, _tbb.ymin, _tbb.ymax, 0, 40).BoundingBox())]
_hw = math.sqrt(max(0.0, (LT_B / 2 - T_FILZ_LT) ** 2 - max(0.0, X_FS_T - (_tbb.xmax - (LT_B / 2 - T_FILZ_LT))) ** 2))
Y_FS_V = max([Y_FS_ARM] + [y_ + D_FS / 2 + 1.0 for y_ in _kly])
Y_FS_V = min(Y_FS_V, LT_CY + _hw - D_FS / 2 - 0.6)
print("Lufttaste-Federstange: Klappen unter dem Arm bis Y", [round(y_, 1) for y_ in _kly], "-> Arm vorne bei Y", round(Y_FS_V, 1))
_arm_ = _zug([_ecke, (X_FS_J, Y_FS_ARM, Z_FS), (X_FS_J, Y_FS_V, Z_FS), (X_FS_T, Y_FS_V, Z_FS), (X_FS_T, Y_FS_V, Z_FS - 8.0)])   # Arm + Schenkel in die Taste
_arm_ = _fs_dreh(_arm_)
_fs = _stab_.fuse(_arm_).clean()
mech.append(("Lufttaste_Federstange", _fs))
for k_ in range(len(teile)):                                         # Taste: Bohrung fuer den Schenkel (fest eingeklebt)
    n_, s_ = teile[k_]
    if n_ == "Lufttaste": teile[k_] = (n_, s_.cut(_arm_))
    if n_ == "Deckelauflage_rechts": teile[k_] = (n_, s_.cut(stab((X_FS, Y_FS_H, Z_FS), (X_FS + 4.5, Y_FS_H, Z_FS), D_FS + 0.1)))
for k_, yb_ in enumerate((Y_FS_ARM - 7.0, Y_FS_H + 6.0)):
    lb_ = box(X_FS - 2.3, _wand.xmin, yb_ - 3.0, yb_ + 3.0, T_DECKEL_BASS + 0.3, Z_FS + 5.5)
    lb_ = lb_.cut(stab((X_FS, yb_ - 4, Z_FS), (X_FS, yb_ + 4, Z_FS), D_FS + 0.1))
    z_sc = Z_FS + 3.5
    lb_ = lb_.cut(stab((X_FS - 3, yb_, z_sc), (_wand.xmin + 1, yb_, z_sc), 1.5))
    sc_ = stab((X_FS - 3.2, yb_, z_sc), (_wand.xmin + 5.0, yb_, z_sc), 1.4).fuse(stab((X_FS - 3.2, yb_, z_sc), (X_FS - 2.3, yb_, z_sc), 2.5))
    mech.append((f"Lufttaste_Lagerbock_{k_+1}", lb_)); mech.append((f"Lufttaste_Lagerbock_Schraube_{k_+1}", sc_))
print("Lufttaste-Federstange: X", round(X_FS, 1), "Z", round(Z_FS, 1), "| Drehstab Y", round(Y_FS_H, 1), "-", round(Y_FS_ARM, 1),
      "| Arm bis X", round(X_FS_T, 1), "r", round(r_fs, 1), "| gedrueckt", round(w_fs, 1), "Grad")
# --- Freischnitte: wo eine Schubstange Boden oder Kammerteile schneidet, dort 1 mm um die Stange ausnehmen ---
KAMMER_PRAEFIX = ("Reihe_", "Schmal", "Einzel_")
fase_kam = {}
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
for n1, s1 in mech:
    if not n1.startswith(("Knopfzylinder", "Schubstange")): continue
    kopien_ = [s1.translate(cq.Vector(0, -sh_, 0)) for sh_ in (HUB / 3, 2 * HUB / 3, HUB)]
    if n1.startswith("Knopfzylinder"):
        r_ = n1.split("_")[1]; xk_ = min([k_[0] for k_ in knoepfe if k_[1] == r_], key=lambda q: abs(q - float(n1.split("_")[2])))
        bb_ = s1.BoundingBox(); Zr_ = REIHEN[r_][0]
        frei_ = [stab((xk_, bb_.ymin - HUB - 0.5, Zr_), (xk_, bb_.ymax + 0.5, Zr_), D_ZYL + 1.0)]
    else:
        pf_ = stange_pfad[n1]; frei_ = []
        for sh_ in (0.0, HUB / 3, 2 * HUB / 3, HUB):
            for p_a, p_b in zip(pf_[:-1], pf_[1:]):
                pa_, pb_ = (p_a[0], p_a[1] - sh_, p_a[2]), (p_b[0], p_b[1] - sh_, p_b[2])
                frei_.append(stab(pa_, pb_, D_STANGE + 1.0)); frei_.append(cq.Solid.makeSphere(D_STANGE / 2 + 0.5, cq.Vector(*pb_)))
    def _trifft(t_):
        return any(_bb_ueber(t_.BoundingBox(), k_.BoundingBox()) and t_.intersect(k_).Volume() > 0.05 for k_ in kopien_)
    def _ausnehmen(t_):
        for f_ in frei_:
            if _bb_ueber(t_.BoundingBox(), f_.BoundingBox()): t_ = t_.cut(f_)
        return t_
    fase_ok = n1.startswith("Knopfzylinder") and _trifft(boden) and any(n2.startswith(KAMMER_PRAEFIX) and _trifft(s2) for n2, s2 in teile)   # nur Kante Boden/Kammer
    if fase_ok:
        # Kammer mit Boden: statt Schlitz eine 45-Grad-Fase knapp ueber den Zylinderweg; innen schliesst ein Fuelldreieck die Oeffnungen
        R_f = (D_ZYL + 1.0) / 2 + 0.5
        Y_c, Z_c = bb_.ymin - HUB - 1.0, Zr_ + R_f
        w_f = 2 * R_f + 0.5
        fase_ = cq.Workplane("YZ").polyline([(Y_c - 40, Z_c - 40), (Y_c + 40, Z_c + 40), (Y_c + 40, Z_c - 40)]).close().extrude(w_f).val().translate(cq.Vector(xk_ - w_f / 2, 0, 0))
        if _bb_ueber(boden.BoundingBox(), fase_.BoundingBox()) and boden.intersect(fase_).Volume() > 0.01:
            boden = boden.cut(fase_); frei_n.append((n1 + " (Fase)", "Bassboden"))
        kam_g = set()
        for k_ in range(len(teile)):
            n2, s2 = teile[k_]
            if n2.startswith(KAMMER_PRAEFIX) and _bb_ueber(s2.BoundingBox(), fase_.BoundingBox()) and s2.intersect(fase_).Volume() > 0.01:
                teile[k_] = (n2, s2.cut(fase_)); frei_n.append((n1 + " (Fase)", n2)); kam_g.add(n2.split("_")[0])
        for kam_ in kam_g: fase_kam.setdefault(kam_, []).append((Y_c - Z_c, fase_))
    if _trifft(boden) and not fase_ok:
        boden = _ausnehmen(boden); frei_n.append((n1 + " (Weg)", "Bassboden"))
    for k_ in range(len(teile)):
        n2, s2 = teile[k_]
        if n2.startswith("Lufttaste") or (fase_ok and n2.startswith(KAMMER_PRAEFIX)): continue
        if _trifft(s2):
            teile[k_] = (n2, _ausnehmen(s2)); frei_n.append((n1 + " (Weg)", n2))
# Abdeckung innen je Gruppe benachbarter Fasen (aehnliche Tiefe): 2 mm Platte unter 45 Grad, aussen buendig mit der tiefsten
# Fase der Gruppe, nur ueber den Bereich dieser Fasen; liegt ganz in der Kammer (nichts ragt ueber die Kammer hinaus).
for kam_, fl_ in fase_kam.items():
    tl_ = [(n_, s_) for n_, s_ in teile if n_.startswith(kam_ + "_")]
    wand_ = max([s_ for n_, s_ in tl_ if "_Trennwand" in n_], key=lambda s_: s_.BoundingBox().ymax).BoundingBox()
    fol_ = [s_.BoundingBox() for n_, s_ in tl_ if n_.endswith("_Folie")]
    z_u = min(f_.zmin for f_ in fol_) - 0.3 if fol_ else Z_BODEN + 15
    Y_w = wand_.ymin
    fl_ = sorted(fl_, key=lambda q_: q_[1].BoundingBox().xmin)
    gruppen_ = []
    for q_ in fl_:
        if gruppen_ and q_[1].BoundingBox().xmin - gruppen_[-1][-1][1].BoundingBox().xmax < 12.0 and abs(q_[0] - gruppen_[-1][-1][0]) < 4.0: gruppen_[-1].append(q_)
        else: gruppen_.append([q_])
    for g_i, g_ in enumerate(gruppen_):
        k_ = min(q_[0] for q_ in g_); t_s = 2.0 * math.sqrt(2)
        x0_ = min(q_[1].BoundingBox().xmin for q_ in g_) - 1.5; x1_ = max(q_[1].BoundingBox().xmax for q_ in g_) + 1.5
        ab_ = cq.Workplane("YZ").polyline([(k_ + Z_BODEN - 5, Z_BODEN - 5), (Y_w + 5, Y_w + 5 - k_), (Y_w + 5 - t_s, Y_w + 5 - k_), (k_ + Z_BODEN - 5 - t_s, Z_BODEN - 5)]).close().extrude(x1_ - x0_).val().translate(cq.Vector(x0_, 0, 0))
        ab_ = ab_.intersect(box(x0_, x1_, Y_w - 80, Y_w, Z_BODEN, z_u)).cut(boden)
        for n_, s_ in tl_: ab_ = ab_.cut(s_)
        for _, f_ in fl_: ab_ = ab_.cut(f_)
        for ns_, pf_ in stange_pfad.items():                           # Schubstangen-Weg freihalten
            for sh_ in (0.0, HUB / 2, HUB):
                for p_a, p_b in zip(pf_[:-1], pf_[1:]):
                    t_ = stab((p_a[0], p_a[1] - sh_, p_a[2]), (p_b[0], p_b[1] - sh_, p_b[2]), D_STANGE + 1.0)
                    if _bb_ueber(ab_.BoundingBox(), t_.BoundingBox()) and ab_.intersect(t_).Volume() > 0.01: ab_ = ab_.cut(t_)
        so_ = sorted(ab_.Solids(), key=lambda q_: -q_.Volume())
        if so_:
            teile.append((f"{kam_}_Abdeckung_{g_i+1}", so_[0]))
            b_ = so_[0].BoundingBox(); zw_ = Y_w - k_
            print(f"Kammer {kam_}: Abdeckung {g_i+1} 2 mm, X {b_.xmin:.1f}-{b_.xmax:.1f} ({len(g_)} Fasen), Y {b_.ymin:.1f}-{b_.ymax:.1f} Z {b_.zmin:.1f}-{b_.zmax:.1f}",
                  "| an Vorderwand Z", round(zw_, 1), "Folie ab", round(z_u + 0.3, 1), "OK" if zw_ <= z_u else "STOESST AN FOLIE", "| Teile", len(so_))
print("Freischnitte:", frei_n)
for k_ in range(len(teile)):                                       # Splitter nach Freischnitten entfernen
    n2, s2 = teile[k_]
    so_ = s2.Solids() if hasattr(s2, "Solids") else [s2]
    if len(so_) > 1:
        vol_ = sorted([round(x_.Volume(), 2) for x_ in so_], reverse=True)
        gross_ = [x_ for x_ in so_ if x_.Volume() > 5.0]
        teile[k_] = (n2, gross_[0] if len(gross_) == 1 else cq.Compound.makeCompound(gross_))
        print("Splitter entfernt", n2, vol_, "" if len(gross_) == 1 else "WARNUNG mehrere grosse Teile")
# --- Pruefung gedrueckt: Lufttaste (10 mm) und jede Basstaste einzeln (6 mm); Schenkel, Gabel, Bolzen, Achse, Klappenhebel,
#     Feder und Klappe drehen um die Achse (Winkel = Hub / Armlaenge), Knopf, Zylinder, Stange und Gelenkkopf 6 mm hinein ---
gehaeuse_alle = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + hoerner + fuesse + klappen + mech
def _bbt(b1, b2):
    return not (b1.xmax < b2.xmin or b2.xmax < b1.xmin or b1.ymax < b2.ymin or b2.ymax < b1.ymin or b1.zmax < b2.zmin or b2.zmax < b1.zmin)
_tg = [s_ for n_, s_ in teile if n_ == "Lufttaste"][0]                    # im Markus-Modell schon gedrueckt dargestellt
tg_kol = [(n1, round(s1.intersect(_tg).Volume(), 2)) for n1, s1 in gehaeuse_alle if not n1.startswith(("Lufttaste", "Luftklappe")) and _bbt(s1.BoundingBox(), _tg.BoundingBox()) and s1.intersect(_tg).Volume() > 0.01]
print("Lufttaste gedrueckt Kollisionen:", tg_kol)
mech_d = dict(mech); ged_kol = {}
ARMLEN = Z_ACHSE - Z_GG
for ai in sorted(set(stange_achse)):
    a0, a1 = achsen[ai]; th_ = math.degrees(HUB / ARMLEN)
    ax_p, ax_q = cq.Vector(a0[0], a0[1], Z_ACHSE), cq.Vector(a1[0], a1[1], Z_ACHSE)
    arms_ = [q for q in arm_info if q[0] == ai]
    pt_ = cq.Vertex.makeVertex(arms_[0][2], arms_[0][3], Z_GG).rotate(ax_p, ax_q, 1.0)
    sg_ = -1.0 if pt_.Center().y > arms_[0][3] else 1.0                # Drehsinn: Gelenk bewegt sich nach hinten (-Y)
    dreh = lambda s_: s_.rotate(ax_p, ax_q, sg_ * th_)
    beweg = {}
    for nm_ in [f"Achse_{ai+1}", f"Feder_Achse_{ai+1}"] + [z_[3] for z_ in zuordnung.get(ai, []) if len(z_) > 3]:
        if nm_ in mech_d: beweg[nm_] = dreh(mech_d[nm_])
    for (_, idx, gx_, gy_, r_, Xk_) in arms_:
        for nm_ in (f"Schenkel_Knopf_{idx}", f"Gelenkgabel_{idx}", f"Gelenkbolzen_{idx}"): beweg[nm_] = dreh(mech_d[nm_])
        for nm_ in (f"Gelenkkopf_{idx}", f"Schubstange_{r_}_{Xk_:.0f}", f"Knopf_{r_}_{Xk_:.0f}", f"Knopfzylinder_{r_}_{Xk_:.0f}"):
            if nm_ in mech_d: beweg[nm_] = mech_d[nm_].translate(cq.Vector(0, -HUB, 0))
    for kn_ in [z_[0] for z_ in zuordnung.get(ai, [])]:
        for n_, s_ in klappen:
            if n_.startswith(kn_ + "_"): beweg[n_] = dreh(s_)
    fest_ = [(n_, s_) for n_, s_ in gehaeuse_alle if n_ not in beweg] + [("Lufttaste_gedrueckt", _tg)]
    for n1, s1 in beweg.items():
        b1 = s1.BoundingBox()
        for n2, s2 in fest_:
            if "Schraube" in n2 or (n1.startswith("Achse_") and n2.startswith("Lagerbock")) or n2.startswith("Luftklappe"): continue
            if not _bbt(b1, s2.BoundingBox()): continue
            v_ = s1.intersect(s2).Volume()
            if v_ > 0.05:
                bi_ = s1.intersect(s2).BoundingBox()
                ged_kol.setdefault(f"Achse {ai+1} Knopf " + ",".join(f"{q[4]}{L - q[5]:.0f}" for q in arms_), []).append((n1, n2, round(v_, 2), [round(bi_.ymin, 1), round(bi_.ymax, 1), round(bi_.zmin, 1), round(bi_.zmax, 1)]))
print("Gedrueckt-Pruefung Achsen", len(set(stange_achse)), "Armlaenge", round(ARMLEN, 1), "Winkel", round(math.degrees(HUB / ARMLEN), 1))
for k_, v_ in ged_kol.items(): print("GEDRUECKT", k_, v_)
for n_, s_ in mech:
    so_ = s_.Solids() if hasattr(s_, "Solids") else [s_]
    if len(so_) > 1: print("MEHRTEILIG", n_, [[round(v_, 1) for v_ in (x_.BoundingBox().xmin, x_.BoundingBox().xmax, x_.BoundingBox().ymin, x_.BoundingBox().ymax, x_.BoundingBox().zmin, x_.BoundingBox().zmax, x_.Volume())] for x_ in so_])
# --- Kollisionen Mechanik gegen Gehaeuse/Kammern ---
gehaeuse = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + hoerner + fuesse + klappen
kol = {}
for n1, s1 in mech:
    bb1 = s1.BoundingBox()
    for n2, s2 in gehaeuse:
        bb2 = s2.BoundingBox()
        if bb1.xmax < bb2.xmin or bb2.xmax < bb1.xmin or bb1.ymax < bb2.ymin or bb2.ymax < bb1.ymin or bb1.zmax < bb2.zmin or bb2.zmax < bb1.zmin: continue
        v = s1.intersect(s2).Volume()
        if v > 0.3 and not n1.startswith(("Lagerbock_Schraube", "Luftklappe_Lagerbock_Schraube", "Lufttaste_Lagerbock_Schraube")): kol.setdefault(n1, []).append((n2, round(v, 1)))
for k, v in sorted(kol.items()): print("Kollision", k, v[:5])
