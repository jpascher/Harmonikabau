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
D_ACHSE, D_STANGE = 2.4, 2.4                 # wie Standard: Achsen, Schenkel, Druckstangen Oe2.4
Z_KL_U = min(s.BoundingBox().zmin for _, s in klappen)   # Unterkante Klappen
# Abstufung je 7 mm: Oberflaeche Reihe 1 = Taschenboden Y 200, Reihe 2 = Front Y 207, Reihe 3 = Frontblende Y 214.
# Knopf-Ruhelage wie im Knopfschnitt der Mechanikzeichnung: 6 mm Hub -> Kopf 6.5 mm vor der jeweiligen Oberflaeche.
# Filzring sitzt in der aeussersten Lage der jeweiligen Stufe.
HUB, SPALT = 6.0, 6.5
STUFE = {"A": (207.0, 5.5), "B": (214.0, 7.0), "C": (221.0, 7.0)}   # Variante 3: Reihe 1 in der Frontflaeche, Reihe 2 mit Blende            # Oberflaeche Y, Dicke der aeussersten Lage
REIHEN = {r_: (z_, STUFE[r_][0] + SPALT) for r_, z_ in (("A", 10.0 + 4.24), ("B", 24.0 + 4.24), ("C", 38.0 + 4.24))}
mech, mech_teile = [], []
def stab(p0, p1, d):
    p0, p1 = cq.Vector(*p0), cq.Vector(*p1); v = p1 - p0
    return cq.Solid.makeCylinder(d / 2, v.Length, p0, v.normalized())
# Variante 3: Bassmechanik wird neu nach Markus-Muster angelegt (Schritt 2); hier nur Knoepfe und Luftklappe.
def bock_solid(breite, bohrungen):
    lb_ = box(-4.3, 4.3, -breite / 2, breite / 2, Z_UB - H_LB, Z_UB); sc_ = None
    for yb in bohrungen:
        lb_ = lb_.cut(stab((-5, yb, Z_ACHSE), (5, yb, Z_ACHSE), D_ACHSE + 0.1)).cut(stab((0, yb, Z_UB - H_LB - 1), (0, yb, Z_UB + 1), 1.5))
        s_ = stab((0, yb, Z_UB - H_LB - 0.9), (0, yb, Z_UB - H_LB + 10), 1.4).fuse(stab((0, yb, Z_UB - H_LB - 0.9), (0, yb, Z_UB - H_LB), 2.5))
        sc_ = s_ if sc_ is None else sc_.fuse(s_)
    return lb_, sc_
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
stange_pfad = {}
# --- Rueckdruckfedern: je Achse eine Schenkelfeder, 3 Windungen um die Achse neben dem Arm,
#     ein Schenkel mit Haken um den Gabelgelenk-Arm (drueckt ihn nach vorne zurueck), der andere stuetzt sich am Boden ab ---
D_DRAHT, R_WIND, N_WIND, P_WIND = 0.5, 1.0 + 0.2 + 0.25, 3, 0.6
def draht(punkte):
    w_ = None
    for p_a, p_b in zip(punkte[:-1], punkte[1:]):
        s_ = stab(p_a, p_b, D_DRAHT).fuse(cq.Solid.makeSphere(D_DRAHT / 2 + 0.01, cq.Vector(*p_b)))
        w_ = s_ if w_ is None else w_.fuse(s_)
    return w_

# --- Bassmechanik Variante 3 nach Markus-Muster (unter dem flachen Boden) ---
# Je Klappe eine Achse quer (X) hinter der Klappe in Lagerboecken an der Bodenunterseite (Achsmitte 2.5 unter dem Boden);
# am Knopf-X ein Arm nach unten zum Gabelgelenk (zweiteiliger Gabelkopf + Bolzen), Schubstange Oe2 gerade nach vorne
# zum Knopfzylinder; Klappenhebel von der Achse senkrecht unter die Klappe, dann waagrecht zur Klappenmitte.
# Knopf gedrueckt -> Armende nach hinten -> Achse dreht -> Klappenhebel nach unten -> Klappe offen.
# Achsen, die sich in X ueberschneiden, liegen in versetzten Reihen (5.8 mm); in einer Reihe stossen zwei Achsen
# an einem gemeinsamen Lagerbock aneinander, wenn der Abstand klein ist.
def fuss(p, a):
    (x0, y0), (x1, y1) = a; dx, dy = x1 - x0, y1 - y0; t = max(0, min(1, ((p[0]-x0)*dx + (p[1]-y0)*dy) / (dx*dx + dy*dy)))
    f = (x0 + t*dx, y0 + t*dy); return f, math.dist(p, f)
def pfad3(pts, d):
    s_ = None
    for a_, b_ in zip(pts[:-1], pts[1:]):
        seg = stab(a_, b_, d).fuse(cq.Solid.makeSphere(d / 2, cq.Vector(*b_)))
        s_ = seg if s_ is None else s_.fuse(seg)
    return s_
XKn = lambda px: L / 2 + (px - 185.5) * 20.0 / 19.0
BASS_T = [([(1, 109.5)], 36.5), ([(2, 119)], 80.5), ([(1, 147.5)], 124.5), ([(2, 157)], 168.5),
          ([(1, 185.5)], 212.5), ([(2, 195)], 256.5), ([(2, 233)], 290.5), ([(1, 223.5)], 334.5)]
BEGL_T = [(1, 128.5), (2, 138), (1, 166.5), (1, 204.5), (2, 176), (2, 214), (1, 242.5), (2, 252)]
kl_c = {}
for n_, s_ in klappen:
    if n_.endswith("_Leder"):
        b_ = s_.BoundingBox(); kl_c[n_[:-6]] = ((b_.xmin + b_.xmax) / 2, (b_.ymin + b_.ymax) / 2, b_.xmax - b_.xmin, b_.ymax - b_.ymin, b_.ymin, b_.ymax)
bass_kl = sorted([k for k in kl_c if k.startswith("Bassklappe")], key=lambda k: kl_c[k][0])
begl_kl = sorted([k for k in kl_c if k.startswith("Begleiterklappe")], key=lambda k: kl_c[k][0])
RB = {1: "A", 2: "B"}
ach = []
for knl, pxk in BASS_T:
    k_ = min(bass_kl, key=lambda k: abs(kl_c[k][0] - (pxk + 6.0)))
    ach.append(dict(art="Bass", kn=[(RB[r_], XKn(px)) for r_, px in knl], kl=k_))
for i, (r_, px) in enumerate(BEGL_T):
    ach.append(dict(art="Begleiter", kn=[(RB[r_], XKn(px))], kl=begl_kl[i]))
ach.append(dict(art="Zusatz", kn=[("A", XKn(261.5))], kl="Zusatzklappe"))
taste_bb = [s_ for n_, s_ in teile if n_ == "Lufttaste"][0].BoundingBox()
TX0, TX1, TY0, TY1 = taste_bb.xmin - 2.5, taste_bb.xmax + 2.5, taste_bb.ymin - 2.5, taste_bb.ymax + 2.5
# Armlage: am Knopf-X; Ausnahmen: Luft-Schubstange (Knopf bei X des Hebelauges) und Zusatzklappe (Boden erst ab X 272.5)
X_AUGE = LT_CX
Y_VERT = {k: (kl_c[k][4] - 0.5) for k in kl_c}                     # Hinterkante der Bodenvertiefung (Klappe + 0.5)
for a in ach:
    a["arm"] = []
    for r_, x_ in a["kn"]:
        xa_ = x_
        if a["art"] == "Zusatz": xa_ = x_ + 4.5                             # Schenkel seitlich versetzt: Stange geht im U neben der Gabel vorbei
        elif abs(x_ - X_AUGE) < 4.0: xa_ = x_ - 6.0
        a["arm"].append((r_, x_, xa_))
    a["xh"] = kl_c[a["kl"]][0]
    if a["art"] == "Zusatz": a["xh"] = max(a["xh"], taste_bb.xmax + D_STANGE / 2 + 0.6)   # Zusatzhebel neben der Lufttaste
import os as _os, json as _json
if _os.environ.get("DUMP_ACH"):
    _json.dump(dict(ach=[{k: v for k, v in a_.items()} for a_ in ach], kl_c=kl_c, Y_VERT=Y_VERT, TY0=TY0, Z_UB=Z_UB, Z_ACHSE=Z_ACHSE, Z_KL_U=Z_KL_U), open("/tmp/bass/ach.json", "w")); raise SystemExit("DUMP")
# Reihen: Achsabstand 3.2 wie die doppelte Lagerbuchse (Bohrungen 2.6 / 5.8); die einfache Buchse hat seitlich eine
# Hohlkehle R1.2, an der die Achse der Nachbarreihe vorbeilaeuft. Zwei Achsen einer Reihe stecken von beiden Seiten
# in derselben Buchse (Schraube M1.4 quer durch die Bohrung = Anschlag).
DY_R = 3.2
def reihen(lst, y_start, dy_):
    rows = []
    for a in sorted(lst, key=lambda a: min([a["xh"]] + [q[2] for q in a["arm"]])):
        a["x0"] = min([a["xh"]] + [q[2] for q in a["arm"]]); a["x1"] = max([a["xh"]] + [q[2] for q in a["arm"]])
        for ri, row in enumerate(rows):
            if ri < a.get("rmin", 0): continue
            if not row or a["x0"] - row[-1]["x1"] >= 15.0: row.append(a); a["ri"] = ri; break
        else:
            while len(rows) < a.get("rmin", 0): rows.append([])
            rows.append([a]); a["ri"] = len(rows) - 1
    rows = [r_ for r_ in rows]
    for a in lst: a["y"] = y_start + dy_ * a["ri"]
    return [r_ for r_ in rows if r_]
GK_H, GK_L, GK_R = 3.19, 13.0, 2.0                          # Fusion-Gelenkkopf: Querbohrung ueber Schaftanfang, Schaftlaenge, Radius
GK_L_KURZ = 5.0                                              # Zusatztaste: Gelenkkopf-Schaft gekuerzt (Stange geht gleich seitlich weg)
Z_TASTE_GEDR = T_DECKEL_BASS + 10.0 + D_STANGE / 2 + 0.8          # gedrueckte Taste (Oberkante 14) + Stangenradius + Luft
Z_MIN_GG = T_DECKEL_BASS + 0.3 + GK_R + 1.5                    # Gelenk 1.5 mm ueber Gaze/Deckel
Z_ST_TIEF = T_DECKEL_BASS + 0.3 + D_STANGE / 2 + 0.3          # Stangen unter den Klappen ganz tief
FREI_KL = Z_KL_U - (Z_ST_TIEF + D_STANGE / 2 + 0.5)            # so weit darf eine Klappenkante beim Oeffnen hinunter
Y_AX_Z = Y_VERT["Zusatzklappe"] - 2.6 - 2.3                      # Zusatzachse knapp hinter der Zusatzklappe (Buchse 2.3 vor der Vertiefung)
y_bass0 = min(Y_VERT[a["kl"]] for a in ach if a["art"] == "Bass") - 2.6 - 1.0     # Platz fuer die oeffnende Klappe
y_begl0 = min(Y_VERT[a["kl"]] for a in ach if a["art"] == "Begleiter") - 2.6 - 0.4
def layout():
    global rows_b, rows_g, rows_z
    rows_b = reihen([a for a in ach if a["art"] == "Bass"], y_bass0, -DY_R)
    rows_g = reihen([a for a in ach if a["art"] == "Begleiter"], y_begl0, -DY_R)
    rows_z = reihen([a for a in ach if a["art"] == "Zusatz"], Y_AX_Z, -DY_R)
# Buchsenlage: am Achsende neben dem aeussersten Element (Arm: Platz fuer die Feder zwischen Arm und Buchse,
# Feder immer auf der Buchsenseite -> linke/rechte Feder; Klappenhebel: nur Drahtradius + Luft)
L_FED_ = 3 * 0.6 + 0.5
def rand(a, x_):
    return (D_STANGE / 2 + 0.4 + L_FED_ + 0.3 + 4.3) if any(abs(q[2] - x_) < 0.05 for q in a["arm"]) else (D_STANGE / 2 + 1.5 + 4.3)
B_HALB = 4.3
def verteilen():
    global gruppen, bl, kl_
    gruppen = [rows_b, rows_g, rows_z]
    def lage(a): return a["y"]
    for rows in gruppen:
        for row in rows:
            for a in row:
                a["e0"] = a["x0"] - rand(a, a["x0"]); a["e1"] = a["x1"] + rand(a, a["x1"]); a["gem1"] = None
            for j in range(len(row) - 1):
                if row[j + 1]["e0"] - row[j]["e1"] < 12.0:                       # gemeinsame Buchse, beide Achsen von den Seiten
                    m_ = (row[j]["e1"] + row[j + 1]["e0"]) / 2
                    if row[j]["x1"] + rand(row[j], row[j]["x1"]) - 0.01 <= m_ + 0.01 or True:
                        row[j]["e1"] = m_; row[j + 1]["e0"] = m_; row[j]["gem1"] = m_
    def alle_buchsen():
        bl = []
        for gi, rows in enumerate(gruppen):
            for row in rows:
                for j, a in enumerate(row):
                    if j == 0 or row[j - 1]["gem1"] is None: bl.append(dict(x=a["e0"], y=a["y"], g=gi, ach=[(id(a), 0)]))
                    b_ = dict(x=a["e1"], y=a["y"], g=gi, ach=[(id(a), 1)])
                    if a["gem1"] is not None: b_["ach"].append((id(row[j + 1]), 0))
                    bl.append(b_)
        return bl
    def achse_bei(y_, x_, m_=0.0):                                   # laeuft in Reihe y_ eine Achse durch X x_ (+-m_)?
        return any(abs(a["y"] - y_) < 0.01 and a["e0"] - m_ < x_ < a["e1"] + m_ for a in ach)
    def konflikte(bl):
        # je Buchse: Nachbarreihen vorne/hinten; Arm-Regel: hinter einem Gelenkarm keine Buchse
        erg = []
        for b in bl:
            k_ = []
            for sg in (1, -1):
                yn_ = b["y"] + sg * DY_R
                nb_ = [c for c in bl if c is not b and abs(c["y"] - yn_) < 0.01 and abs(c["x"] - b["x"]) < 2 * B_HALB + 0.5]
                if nb_:
                    if all(abs(c["x"] - b["x"]) < 0.01 for c in nb_): k_.append(("doppel", sg))
                    else: k_.append(("ueberlapp", sg))
                elif achse_bei(yn_, b["x"], B_HALB): k_.append(("vorbei", sg))
                # uebernaechste Reihe darf nicht an einer Doppelbuchse vorbei (keine Hohlkehle aussen)
            for a in ach:
                for q in a["arm"]:
                    if a["y"] - 2 * DY_R - 0.1 < b["y"] < a["y"] - 0.01 and abs(q[2] - b["x"]) < B_HALB + GK_R0 + 1.0: k_.append(("hinterArm", 0))
                # Nachbarreihe: an Arm, Feder oder Klappenhebel der vorbeilaufenden Achse darf keine Buchse stehen
                if abs(abs(a["y"] - b["y"]) - DY_R) < 0.01:
                    for xe_ in [q[2] for q in a["arm"]] + [a["xh"]]:
                        if abs(xe_ - b["x"]) < B_HALB + D_STANGE / 2 + 0.4: k_.append(("anElement", 0))
                    for q in a["arm"]:
                        sd_ = -1 if (q[2] - a["e0"]) < (a["e1"] - q[2]) else 1
                        f0_, f1_ = sorted((q[2] - sd_ * (D_STANGE / 2 + 0.4), q[2] - sd_ * (D_STANGE / 2 + 0.4 + L_FED_ + 0.6)))
                        if f0_ < b["x"] + B_HALB and b["x"] - B_HALB < f1_: k_.append(("anFeder", 0))
            erg.append(k_)
        return erg
    GK_R0 = 2.0
    def bewertung(bl):
        n_ = 0
        for b, k_ in zip(bl, konflikte(bl)):
            arten = [t for t, _ in k_]
            if "ueberlapp" in arten or "hinterArm" in arten or "anElement" in arten or "anFeder" in arten: n_ += 10
            if arten.count("vorbei") + arten.count("doppel") > 1: n_ += 5
        return n_
    # Konflikte loesen: Buchse (= Achsende) in 0.5-Schritten nach aussen schieben, gemeinsame Buchsen in beide Richtungen
    by_id = {id(a): a for a in ach}
    for runde in range(6):
        bl = alle_buchsen(); kl_ = konflikte(bl); bes = False
        for b, k_ in zip(bl, kl_):
            arten = [t for t, _ in k_]
            if not ("ueberlapp" in arten or "hinterArm" in arten or "anElement" in arten or "anFeder" in arten or arten.count("vorbei") + arten.count("doppel") > 1): continue
            best = None
            for st_ in [0.5 * k for k in range(1, 61)]:
                for sg in (-1, 1):
                    ok_ = True; alt = []
                    for (ida, en) in b["ach"]:
                        a = by_id[ida]
                        if len(b["ach"]) == 1 and ((en == 0 and sg > 0) or (en == 1 and sg < 0)): ok_ = False
                        alt.append((a, en, a["e0"] if en == 0 else a["e1"]))
                    if not ok_: continue
                    for a, en, v_ in alt:
                        if en == 0: a["e0"] = v_ + sg * st_
                        else: a["e1"] = v_ + sg * st_
                        if a["gem1"] is not None and en == 1: a["gem1"] = a["e1"]
                    gueltig = all(a["e0"] <= a["x0"] - rand(a, a["x0"]) + 0.01 and a["e1"] >= a["x1"] + rand(a, a["x1"]) - 0.01 for a, _, _ in alt)
                    # gleiche Reihe: keine Ueberschneidung mit der naechsten Achse
                    for a, en, _ in alt:
                        for c in ach:
                            if c is not a and abs(c["y"] - a["y"]) < 0.01 and id(c) not in [x_ for x_, _ in b["ach"]] and c["e0"] - 0.01 < a["e1"] and a["e0"] < c["e1"] + 0.01 and not (abs(c["e0"] - a["e1"]) < 0.01 or abs(c["e1"] - a["e0"]) < 0.01): gueltig = False
                    sc_ = bewertung(alle_buchsen()) if gueltig else 999
                    for a, en, v_ in alt:
                        if en == 0: a["e0"] = v_
                        else: a["e1"] = v_
                        if a["gem1"] is not None and en == 1: a["gem1"] = a["e1"]
                    if best is None or sc_ < best[0]: best = (sc_, sg, st_)
                if best is not None and best[0] < bewertung(alle_buchsen()) and best[0] == 0: break
            if best is not None and best[0] < bewertung(alle_buchsen()):
                sc_, sg, st_ = best
                for (ida, en) in b["ach"]:
                    a = by_id[ida]
                    if en == 0: a["e0"] += sg * st_
                    else:
                        a["e1"] += sg * st_
                        if a["gem1"] is not None: a["gem1"] = a["e1"]
                bes = True
        if not bes: break
    bl = alle_buchsen(); kl_ = konflikte(bl)

    bl = alle_buchsen(); kl_ = konflikte(bl)
def schlecht(k_):
    arten = [t for t, _ in k_]
    return "ueberlapp" in arten or "hinterArm" in arten or "anElement" in arten or "anFeder" in arten or arten.count("vorbei") + arten.count("doppel") > 1
by_id = {id(a): a for a in ach}
Y_MIN_R = {"Bass": 44.0, "Begleiter": 104.0, "Zusatz": 140.0}      # hinterste zulaessige Achslage je Gruppe
# Reihenzuordnung (systematisch gesucht: alle Kombinationen ohne Ueberschneidung in einer Reihe, keine Buchsenkonflikte)
REIHE_FEST = {'Bassklappe_1': 0, 'Bassklappe_2': 2, 'Bassklappe_3': 0, 'Bassklappe_4': 1, 'Bassklappe_5': 0, 'Bassklappe_6': 2, 'Bassklappe_8': 3, 'Bassklappe_7': 0, 'Begleiterklappe_1': 0, 'Begleiterklappe_2': 1, 'Begleiterklappe_3': 2, 'Begleiterklappe_4': 3, 'Begleiterklappe_5': 0, 'Begleiterklappe_6': 0, 'Begleiterklappe_7': 2, 'Begleiterklappe_8': 0, 'Zusatzklappe': 0}
layout()
y_st_ = {"Bass": y_bass0, "Begleiter": y_begl0, "Zusatz": Y_AX_Z}
SCHRAEG = sorted([a for a in ach if a["art"] == "Bass"], key=lambda a: a["xh"])[5:]      # Bassklappen 6, 7, 8 (Seite Zusatzkammer)
for a in ach: a["ri"] = REIHE_FEST[a["kl"]]; a["y"] = y_st_[a["art"]] - DY_R * a["ri"]
for a in SCHRAEG: a["y"] = -1000.0
def _rows(lst):
    rr_ = {}
    for a in lst: rr_.setdefault(a["ri"], []).append(a)
    return [sorted(rr_[k], key=lambda a: a["x0"]) for k in sorted(rr_)]
rows_b = _rows([a for a in ach if a["art"] == "Bass" and a not in SCHRAEG]); rows_g = _rows([a for a in ach if a["art"] == "Begleiter"]); rows_z = _rows([a for a in ach if a["art"] == "Zusatz"])
verteilen()
boecke = []
for b, k_ in zip(bl, kl_):
    arten = [t for t, _ in k_]
    if "ueberlapp" in arten or "hinterArm" in arten or "anElement" in arten or "anFeder" in arten or arten.count("vorbei") + arten.count("doppel") > 1:
        print("WARNUNG Lagerbuchse", round(b["x"], 1), round(b["y"], 1), k_)
    dop = [sg for t, sg in k_ if t == "doppel"]; vor = [sg for t, sg in k_ if t == "vorbei"]
    if dop and dop[0] < 0: continue                                    # Doppelbuchse wird von der vorderen Reihe angelegt
    boecke.append(dict(x=b["x"], y=b["y"], art="doppelt" if dop else "einfach", seite=(dop[0] if dop else (vor[0] if vor else -1))))
for a in ach:
    if a not in SCHRAEG: a["p0"], a["p1"] = (a["e0"], a["y"]), (a["e1"], a["y"])
for b in boecke: b.setdefault("ang", 0.0)
# --- Bass: links und rechts je drei Achsen leicht schraeg (Klappen 1-3 und 6-8): Lagerbuchse am Klappenende knapp
#     hinter der Klappe (vorderste Reihe), zum Knopf hin aufgefaechert nach hinten, ohne Kreuzung; Mitte (4, 5) gerade ---
def ya(a, x_):
    (xa_, ya_), (xb_, yb_) = a["p0"], a["p1"]
    return ya_ if abs(xb_ - xa_) < 1e-9 else ya_ + (yb_ - ya_) * (x_ - xa_) / (xb_ - xa_)
y_f = y_bass0
gesetzt = []
for a in ach:
    if a["art"] == "Bass" and a not in SCHRAEG:
        gesetzt.append((a, [a["e0"], a["e1"]] + ([a["gem1"]] if a.get("gem1") else [])))
SEP_A, SEP_B, SEP_ARM = 3.3, 4.1, 6.5
for trio in (SCHRAEG,):
    for a in trio:
        if a["xh"] < a["arm"][0][2]: xk_, xn_ = a["x0"] - rand(a, a["x0"]), a["x1"] + rand(a, a["x1"])
        else: xk_, xn_ = a["x1"] + rand(a, a["x1"]), a["x0"] - rand(a, a["x0"])
        yn_ = y_f; lo_, hi_ = min(xk_, xn_), max(xk_, xn_)
        xs_ = [lo_ + 0.25 * k for k in range(int((hi_ - lo_) / 0.25) + 1)] + [xn_]
        for (c, bx_l) in gesetzt:
            clo_, chi_ = min(c["p0"][0], c["p1"][0]), max(c["p0"][0], c["p1"][0])
            for x_ in xs_:
                if not (clo_ - 6.0 <= x_ <= chi_ + 6.0): continue
                sep_ = SEP_A if clo_ <= x_ <= chi_ else 0.0
                if any(abs(x_ - bx_) <= B_HALB + 1.5 for bx_ in bx_l): sep_ = SEP_B
                if abs(x_ - xn_) <= B_HALB + 1.5 and clo_ - 1 <= x_ <= chi_ + 1: sep_ = max(sep_, SEP_B)
                if abs(x_ - xn_) < 0.01 and any(abs(q[2] - xn_) < B_HALB + GK_R + 1.0 for q in c["arm"]): sep_ = max(sep_, SEP_ARM)
                if sep_ == 0.0: continue
                ymax_ = ya(c, min(max(x_, clo_), chi_)) - sep_
                t_ = (x_ - xk_) / (xn_ - xk_)
                if t_ < 0.02:
                    if ymax_ < y_f - 0.01: print("WARNUNG Schraegachse Klappenende", a["kl"], round(x_, 1))
                    continue
                yn_ = min(yn_, y_f + (ymax_ - y_f) / t_)
        a["p0"], a["p1"] = (xk_, y_f), (xn_, yn_); a["e0"], a["e1"] = lo_, hi_; a["y"] = y_f
        boecke.append(dict(x=xk_, y=y_f, art="einfach", seite=-1, ang=math.degrees(math.atan2(yn_ - y_f, xn_ - xk_))))
        boecke.append(dict(x=xn_, y=yn_, art="einfach", seite=1, ang=math.degrees(math.atan2(yn_ - y_f, xn_ - xk_))))
        gesetzt.append((a, [xk_, xn_]))
        print(f"Schraegachse {a['kl']}: Klappenende X {xk_:.1f} Y {y_f:.1f} -> Knopfende X {xn_:.1f} Y {yn_:.1f} ({math.degrees(math.atan2(abs(yn_ - y_f), abs(xn_ - xk_))):.1f} Grad)")
print("Lagerbuchsen", len(boecke), "einfach", sum(b["art"] == "einfach" for b in boecke), "doppelt", sum(b["art"] == "doppelt" for b in boecke))
# Hub ausgleichen: Armlaenge = Abstand Achse -> Klappenmitte * 6 / 9 (Ziel-Klappenhub 9 mm), 10 .. 16.7 mm
for a in ach:
    yk_ = ya(a, a["xh"])
    d_ = kl_c[a["kl"]][1] - yk_
    e_ = kl_c[a["kl"]][5] - yk_                                  # Achse -> Klappenvorderkante
    a["armlen"] = Z_ACHSE - Z_MIN_GG                                 # alle Arme so lang wie moeglich (Gelenk 1.5 ueber Gaze): kleiner Klappenschwenk, Gelenkkoepfe tief
    a["zgg"] = Z_ACHSE - a["armlen"]; a["hub"] = d_ / a["armlen"] * HUB
    a["ygg"] = a["y"]; a["gkl"] = GK_L_KURZ if a["art"] == "Zusatz" else GK_L
Z_KH = Z_KL_U - D_STANGE / 2 - 0.05
achsen, arme = [], {}
for i, a in enumerate(ach):
    (xa_, ya_), (xb_, yb_) = a["p0"], a["p1"]; ll_ = math.dist(a["p0"], a["p1"]); ux_, uy_ = (xb_ - xa_) / ll_, (yb_ - ya_) / ll_
    q0_, q1_ = (xa_ + ux_, ya_ + uy_), (xb_ - ux_, yb_ - uy_)
    achsen.append((q0_, q1_))
    a["sol"] = stab((q0_[0], q0_[1], Z_ACHSE), (q1_[0], q1_[1], Z_ACHSE), D_ACHSE)
    a["sol_l"] = stab((xa_ - 2 * ux_, ya_ - 2 * uy_, Z_ACHSE), (xb_ + 2 * ux_, yb_ + 2 * uy_, Z_ACHSE), D_ACHSE)
    mech.append((f"Achse_{i+1}", a["sol"]))
# Lagerbuchsen nach lager-bohrung-beide.dxf: Querschnitt (YZ) ueber 8.6 in X extrudiert; einfach 5.2 breit mit
# Hohlkehle R1.2 auf einer Seite (Nachbarachse), doppelt 8.4 breit mit zwei Bohrungen 3.2 auseinander; je Bohrung
# eine Schraube M1.4 x 10 in Laengsmitte von unten in den Boden (Anschlag fuer die Achsen).
_LBD = ezdxf.readfile("lagerbeide.dxf").modelspace()
def _prof(lay):
    e_ = [e for e in _LBD.query(f'LWPOLYLINE[layer=="{lay}"]') if e.closed][0]
    pp_ = [(q.x, q.y) for q in ezdxf.path.make_path(e_).flattening(0.02)]
    return pp_[:-1] if math.dist(pp_[0], pp_[-1]) < 1e-6 else pp_
PROF_E, PROF_D = _prof("Layer1"), _prof("0")
def buchse(b):
    pp_ = PROF_E if b["art"] == "einfach" else PROF_D
    g_ = b["seite"] if b["art"] == "doppelt" else -b["seite"]           # einfach: Hohlkehle zur Nachbarachse; doppelt: zweite Bohrung zur Nachbarreihe
    yz_ = [(g_ * (-q[1] - 2.6), Z_UB - (q[0] + 13.6)) for q in pp_]
    sol_ = cq.Workplane("YZ").polyline(yz_).close().extrude(8.6).val().translate(cq.Vector(-4.3, 0, 0))
    ys_ = [0.0] + ([g_ * DY_R] if b["art"] == "doppelt" else [])
    sc_ = None
    for yb in ys_:
        sol_ = sol_.cut(stab((-5, yb, Z_ACHSE), (5, yb, Z_ACHSE), D_ACHSE + 0.1)).cut(stab((0, yb, Z_UB - 6), (0, yb, Z_UB + 1), 1.5))
        s_ = stab((0, yb, Z_UB - 5.0 - 0.9), (0, yb, Z_UB + 5.0), 1.4).fuse(stab((0, yb, Z_UB - 5.0 - 0.9), (0, yb, Z_UB - 5.0), 2.5))
        sc_ = s_ if sc_ is None else sc_.fuse(s_)
    w_ = b.get("ang", 0.0)
    sol_ = sol_.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), w_).translate(cq.Vector(b["x"], b["y"], 0))
    sc_ = sc_.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), w_).translate(cq.Vector(b["x"], b["y"], 0))
    return sol_, sc_
for j, b in enumerate(boecke):
    print(f"BUCHSE {j+1} X {b['x']:.1f} Y {b['y']:.1f} {b['art']} Winkel {b.get('ang', 0):.1f}")
    lb_, sc_ = buchse(b)
    mech.append((f"Lagerbuchse_{b['art']}_{j+1}", lb_))
    for k_, sv_ in enumerate(cq.Shape.cast(sc_.wrapped).Solids() if hasattr(sc_, "Solids") else [sc_]):
        mech.append((f"Lagerbuchse_Schraube_{j+1}" + (f"_{k_+1}" if b["art"] == "doppelt" else ""), sv_))
bpos = []
for b in boecke:
    bpos.append((b["x"], b["y"]))
    if b["art"] == "doppelt": bpos.append((b["x"], b["y"] + b["seite"] * DY_R))
# Boden links der Zusatzkammer (vorderer Bereich offen) unter den Lagerbuchsen der Zusatzachse verlaengern
for bx_, by_ in bpos:
    if abs(by_ - Y_AX_Z) < 0.01 and bx_ < XZ0_ + 4.3:
        # Lagerblock fuer die Buchse im offenen Bereich: an die Bodenkante und die Seitenleiste der Zusatzkammer geleimt
        _bv = box(bx_ - 4.3 - 1.5, XZ0_, by_ - 2.6 - 1.5, by_ + 2.6 + 1.5, Z_UB, Z_BODEN + 8.0)
        for n2_, s2_ in teile:
            if s2_.BoundingBox().xmin < XZ0_ + 0.5 and _bv.intersect(s2_).Volume() > 0.01: _bv = _bv.cut(s2_)
        teile.append(("Lagerblock_Zusatz", _bv)); print("Lagerblock Zusatz: X", round(bx_ - 5.8, 1), "-", XZ0_, "Y", round(by_ - 4.1, 1), "-", round(by_ + 4.1, 1), "Z", round(Z_UB, 1), "-", round(Z_BODEN + 8.0, 1))
# Klappenhebel
for i, a in enumerate(ach):
    c_, w_ = kl_c[a["kl"]][0], kl_c[a["kl"]][2]; ykl_ = kl_c[a["kl"]][1]; yh0_ = ya(a, c_)
    stoer = [bx_ for bx_, by_ in bpos if yh0_ + 1.0 < by_ < ykl_] + [q[2] for b in ach if b is not a for q in b["arm"] if yh0_ + 1.0 < ya(b, q[2]) < ykl_]
    kand = sorted([c_ + 0.25 * k * sg for k in range(0, int(w_ / 2 - 3.0) * 4 + 1) for sg in (1, -1)], key=lambda x_: abs(x_ - c_))
    nah_ = [(bx_, 6.0 if abs(by_ - yh0_) < 0.01 else 6.5) for bx_, by_ in bpos if abs(by_ - yh0_) < 5.5]   # Buchsen gleiche Reihe / Nachbarreihe
    arm_x = [q[2] for b in ach if b is not a for q in b["arm"] if yh0_ + 1.0 < ya(b, q[2]) < ykl_]
    bock_v = [bx_ for bx_, by_ in bpos if yh0_ + 1.0 < by_ < ykl_]
    x_min_ = (taste_bb.xmax + D_STANGE / 2 + 0.6) if a["art"] == "Zusatz" else -1e9    # Zusatzhebel nicht ueber der Lufttaste
    ok_ = lambda x_: x_ >= x_min_ and all(abs(x_ - sx_) > 5.6 for sx_ in bock_v) and all(abs(x_ - sx_) >= 4.8 for sx_ in arm_x) and all(abs(x_ - bx_) >= d_ for bx_, d_ in nah_)
    a["xh"] = next((x_ for x_ in kand if ok_(x_)), c_)
    if not ok_(a["xh"]): print("WARNUNG Klappenhebel ohne freie Lage", i + 1)
    xh_ = a["xh"]
    yh_ = ya(a, xh_)
    unter_ = [by_ for bx_, by_ in bpos if yh_ + 1.0 < by_ < ykl_ and abs(bx_ - xh_) < 4.3 + 1.3]
    if unter_:                                                         # unter einem Lagerbock durch: tiefer, danach wieder hoch
        z_lo = Z_UB - H_LB - 0.9 - D_STANGE / 2 - 0.15; y_r = max(unter_) + 2.6 + 1.2 + 4.5            # Luft fuer das Zurueckschwenken beim Druecken
        kh = pfad3([(xh_, yh_, Z_ACHSE), (xh_, yh_, z_lo), (xh_, y_r, z_lo), (xh_, y_r + (Z_KH - z_lo), Z_KH), (xh_, ykl_, Z_KH)], D_STANGE)
    else:
        kh = pfad3([(xh_, yh_, Z_ACHSE), (xh_, yh_, Z_KH), (xh_, ykl_, Z_KH)], D_STANGE)
    mech.append((f"Klappenhebel_{i+1}", kh.cut(a["sol_l"]))); print(f"HEBEL {i+1} X {xh_:.1f} Y {yh_:.1f} unter {['%.1f' % u_ for u_ in unter_]}")
# Knoepfe sind oben angelegt; Arm, Gabelkopf, Bolzen, Schubstange je Knopf
Y_BODEN_OFFEN = gy0 + T_UNTEN + 0.5                                 # ab hier ist der Boden ausgenommen (vorderer Bereich)
# Gelenk wie Standard-Bassteil: Fusion "gelenkkopf" (Schaft Oe4 x 13, Kopf, Querbohrung Oe1.5 3.19 ueber dem Schaft)
# an der Druckstange, Gegenstueck = Gabel am Schenkeldraht Oe2.4 (Wangen beidseits des Kopfes), Bolzen Oe1.4.
_gk = cq.importers.importStep("kb_gelenkkopf.step").solids().vals()
GK = _gk[0].fuse(_gk[1]).clean()
GK = GK.cut(stab((0, 0.5, 0), (0, -GK_L - 1, 0), D_STANGE)).cut(stab((-3, GK_H, 0), (3, GK_H, 0), 1.5))
L_GABEL = 8.0                                                # Gabel gekuerzt (Fusion-Gabel 10)
def gelenk(P, s_, draht_s, draht_st, nm, gkl=GK_L, rueck=False):
    x, y, z = P
    gk_ = GK if gkl >= GK_L else GK.cut(box(-5, 5, -GK_L - 2, -gkl, -5, 5))
    if rueck: kopf = gk_.translate(cq.Vector(x, y - GK_H, z)).cut(draht_st)                 # Schaft nach hinten (vom Knopf weg)
    else: kopf = gk_.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180).translate(cq.Vector(x, y + GK_H, z)).cut(draht_st)
    mech.append((f"Gelenkkopf_{nm}", kopf))
    g = stab((x - GK_R, y, z), (x + GK_R, y, z), 2 * GK_R).fuse(stab((x, y, z), (x, y, z - s_ * L_GABEL), 2 * GK_R))
    g = g.cut(box(x - 1.05, x + 1.05, (y - 6) if rueck else (y - 2.4), (y + 2.4) if rueck else (y + 6), *sorted((z - s_ * 2.4, z + s_ * 2.4))))
    g = g.cut(stab((x, y, z - s_ * GK_H), (x, y, z - s_ * (L_GABEL + 1)), D_STANGE)).cut(stab((x - 3, y, z), (x + 3, y, z), 1.5)).cut(draht_s)
    mech.append((f"Gelenkgabel_{nm}", g))
    mech.append((f"Gelenkbolzen_{nm}", stab((x - 2.4, y, z), (x + 2.4, y, z), 1.4)))
schenkel = {}
n_st = 0
for i, a in enumerate(ach):
    for (r_, Xk, gx) in a["arm"]:
        n_st += 1; gy = ya(a, gx) if a["art"] != "Zusatz" else a["ygg"]; zg = a["zgg"]
        Zr, Yk = REIHEN[r_]
        arme.setdefault(i, []).append((gx, ya(a, gx)))
        # Schenkel: Draht Oe2.4 von der Achse nach unten in den Gabelschaft (bei der Zusatzachse schraeg nach hinten)
        sp_ = [(gx, gy, Z_ACHSE), (gx, gy, zg + 4.5)]
        sk_ = pfad3(sp_, D_STANGE).cut(a["sol_l"])
        mech.append((f"Schenkel_Knopf_{n_st}_{r_}_{Xk:.0f}", sk_)); schenkel[i] = sp_
        y0 = gy + GK_H + a["gkl"] + 0.5; y_end = Yk - L_ZYL_R[r_]
        ueber_taste = (TX0 < max(gx, Xk) and min(gx, Xk) < TX1)
        z_mid = max(Zr, Z_TASTE_GEDR) if ueber_taste else Zr
        if z_mid > Z_UB - 1.4: y_hoch = max(Y_BODEN_OFFEN, y0)            # ueber die Boden-Unterkante erst im offenen Bereich
        else: y_hoch = (TY0 - (z_mid - zg) - 1.0) if ueber_taste else y_end - (Zr - zg) - 8.0
        y_hoch = max(y_hoch, y0)
        xs_ = (min(gx, Xk) - 2.5, max(gx, Xk) + 2.5)
        y_frei = max([kl_c[k_][5] + 2.0 for k_ in kl_c if kl_c[k_][0] + kl_c[k_][2] / 2 > xs_[0] and kl_c[k_][0] - kl_c[k_][2] / 2 < xs_[1] and kl_c[k_][5] > y0 and k_ != "Zusatzklappe"] + [y0])
        if y_frei > y0 + 0.01:                                          # unter Klappen: ganz tief
            pf = [(gx, gy + 4.6, zg), (gx, y0, zg), (gx, y0 + (zg - Z_ST_TIEF), Z_ST_TIEF)]; z_t = Z_ST_TIEF
        else: pf = [(gx, gy + 4.6, zg), (gx, y0, zg)]; z_t = zg
        y_h = max(y_frei, pf[-1][1])
        if z_mid > Z_UB - D_STANGE: y_h = max(y_h, Y_BODEN_OFFEN)                 # ueber Boden-Unterkante erst im offenen Bereich
        for bx_, by_ in bpos:                                           # unter einer Buchse durch: erst dahinter steigen (auch gedrueckt)
            if abs(bx_ - gx) < B_HALB + D_STANGE / 2 + 0.5 and by_ > pf[-1][1] - 3.0:
                y_h = max(y_h, by_ + 2.6 + HUB + 1.0 - ((Z_UB - H_LB - 0.9) - D_STANGE / 2 - 0.5 - z_t))
        if a["art"] != "Zusatz" and ueber_taste and y_h + (min(z_mid, Z_TASTE_GEDR) - z_t) > TY0 - 1.0: print("WARNUNG Stange steigt erst ueber der Lufttaste", Xk)
        pf.append((gx, y_h, z_t)); pf.append((gx, y_h + (z_mid - z_t), z_mid))
        y_j = max(pf[-1][1], LT_CY + 11.0) if abs(Xk - X_AUGE) < 4.0 else pf[-1][1]   # hinter der Luft-Schubstange auch bei gedruecktem Knopf
        if abs(gx - Xk) > 0.05:                                         # seitlich zum Knopf: gerader Querabschnitt
            pf += [(gx, y_j, z_mid), (Xk, y_j, z_mid)]
        y_last = y_end - abs(Zr - z_mid)
        if y_last < pf[-1][1]: y_last = pf[-1][1]
        pf += [(Xk, y_last, z_mid), (Xk, y_end, Zr)]
        if a["art"] == "Zusatz":                                       # Gelenkkopf nach hinten, Stange im senkrechten U nach oben, ueber die Lufttaste zum Knopf
            yb_ = gy - GK_H - a["gkl"] - 0.5; zo_ = Z_TASTE_GEDR
            pf = [(gx, gy - 4.6, zg), (gx, yb_, zg), (Xk, yb_ - 1.5, zg + 1.5), (Xk, yb_ - 1.5, zo_ - 1.5), (Xk, yb_, zo_),
                  (Xk, y_end - (zo_ - Zr), zo_), (Xk, y_end, Zr)]
            print(f"Zusatzstange: Gelenk X {gx:.1f} Y {gy:.1f} Z {zg:.1f}, U hinten bis Y {yb_ - 1.5:.1f}, oben Z {zo_:.1f}")
        pf = [p for k_, p in enumerate(pf) if k_ == 0 or math.dist(p, pf[k_ - 1]) > 0.05]
        st = None
        for p_a, p_b in zip(pf[:-1], pf[1:]):
            seg = stab(p_a, p_b, D_STANGE).fuse(cq.Solid.makeSphere(D_STANGE / 2 + 0.01, cq.Vector(*p_b)))
            st = seg if st is None else st.fuse(seg)
        nm_ = f"Schubstange_{r_}_{Xk:.0f}"
        mech.append((nm_, st)); stange_pfad[nm_] = pf
        gelenk((gx, gy, zg), -1, sk_, st, f"{n_st}_{r_}_{Xk:.0f}", a["gkl"], rueck=(a["art"] == "Zusatz"))
for i, a in enumerate(ach):
    print(f"MECH Achse {i+1} {a['art']:9s} Reihe {a['ri']+1} Y {a['y']:.1f} X {a['e0']:.1f}-{a['e1']:.1f} | Knoepfe {[(r_, round(x_)) for r_, x_, _ in a['arm']]} -> {a['kl']} | Arm {a['armlen']:.1f} Hub {a['hub']:.1f}")
print("Lagerbuchsen", len(boecke), "| Reihen Bass", len(rows_b), "Begleiter", len(rows_g))
# Rueckdruckfedern (wie Markus): je Achse eine Schenkelfeder neben dem ersten Arm, Haken hinter den Arm, Schenkel an den Boden
n_feder = 0
for ai, al in sorted(arme.items()):
    (a0, a1) = achsen[ai]; gx, gy = al[0]
    Ln_ = math.dist(a0, a1); dx_, dy_ = (a1[0] - a0[0]) / Ln_, (a1[1] - a0[1]) / Ln_
    nx_, ny_ = -dy_, dx_
    if ny_ > 0: nx_, ny_ = -nx_, -ny_                             # n zeigt nach hinten (hinter den Arm)
    F_, _ = fuss((gx, gy), achsen[ai]); t_arm = (F_[0] - a0[0]) * dx_ + (F_[1] - a0[1]) * dy_
    L_SP = N_WIND * P_WIND
    s_dir = 1.0 if (Ln_ - t_arm) < t_arm else -1.0                  # Feder immer zwischen Arm und naechster Buchse (linke/rechte Feder)
    t0 = t_arm + s_dir * (D_STANGE / 2 + 0.4)                       # Federanfang neben dem Schenkeldraht
    hel = cq.Wire.makeHelix(P_WIND, L_SP, R_WIND)
    sp = cq.Workplane("XZ").center(R_WIND, 0).circle(D_DRAHT / 2).sweep(cq.Workplane().add(hel), isFrenet=True).val()
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 1, 0), 90)     # Wendelachse Z -> X
    if s_dir < 0: sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), 180)
    sp = sp.rotate(cq.Vector(0, 0, 0), cq.Vector(0, 0, 1), math.degrees(math.atan2(dy_, dx_))).translate(cq.Vector(a0[0] + dx_ * t0, a0[1] + dy_ * t0, Z_ACHSE))
    pA = (a0[0] + dx_ * t0, a0[1] + dy_ * t0)                      # Federende am Arm
    pE = (a0[0] + dx_ * (t0 + s_dir * L_SP), a0[1] + dy_ * (t0 + s_dir * L_SP))
    _sp = schenkel[ai]; _yh = _sp[-1][1]                          # Gabelschaft (Y) unter der Achse bzw. hinten versetzt
    z_h = _sp[-1][2] + 2.5                                         # Haken um den Gabelschaft (wie Standard)
    ab = (gy - _yh) + GK_R + 0.1 + D_DRAHT / 2; sd_ = GK_R + 0.1 + D_DRAHT / 2
    s1 = draht([(pA[0], pA[1], Z_ACHSE - R_WIND),
                (pA[0] + nx_ * ab, pA[1] + ny_ * ab, Z_ACHSE - R_WIND - 1.0),
                (pA[0] + nx_ * ab, pA[1] + ny_ * ab, z_h),
                (F_[0] - dx_ * s_dir * sd_ + nx_ * ab, F_[1] - dy_ * s_dir * sd_ + ny_ * ab, z_h),
                (F_[0] - dx_ * s_dir * sd_ + nx_ * (gy - _yh - 0.8), F_[1] - dy_ * s_dir * sd_ + ny_ * (gy - _yh - 0.8), z_h)])
    s2 = draht([(pE[0], pE[1], Z_ACHSE + R_WIND),
                (pE[0] - nx_ * 0.8, pE[1] - ny_ * 0.8, Z_UB - 0.3)])      # Gegenschenkel steil an den Boden (zwischen den Reihen)
    mech.append((f"Feder_Achse_{ai+1}", sp.fuse(s1).fuse(s2))); n_feder += 1

print("Rueckdruckfedern", n_feder)

# --- Luftklappe wie beim Standard-Bassteil: Klappe unter dem Zwischendeckel (oeffnet nach unten in den vorderen
#     Bereich), Lufthebel = Alu-Profil nach Lufthebel-lagerbock-kurz146mm.dxf, Lagerbock (zwei Stuetzen 8 breit, je 2
#     Schrauben) am Zwischendeckel, zwei spiegelgleiche Schenkelfedern; Taste -> Schubstange Ø3 mit Gabelgelenk ans
#     Hebelauge. Hebelauge ueber der Tastenmitte, Klappe (Lage aus der Hebelzeichnung) ueber dem langen Arm.
#     Dargestellt in Ruhelage (Taste aussen, Klappe zu).
LH = ezdxf.readfile("kb_lufthebel.dxf").modelspace()
def _pts(e, tol=0.05):
    p = [(q.x, q.y) for q in ezdxf.path.make_path(e).flattening(tol)]
    q = [p[0]] + [b_ for a_, b_ in zip(p[:-1], p[1:]) if math.dist(a_, b_) > 1e-6]
    if len(q) > 2 and math.dist(q[0], q[-1]) < 1e-6: q = q[:-1]
    return q
def pfad(pts, d):
    s_ = None
    for a_, b_ in zip(pts[:-1], pts[1:]):
        seg = stab(a_, b_, d).fuse(cq.Solid.makeSphere(d / 2 + 0.01, cq.Vector(*b_)))
        s_ = seg if s_ is None else s_.fuse(seg)
    return s_.fuse(cq.Solid.makeSphere(d / 2 + 0.01, cq.Vector(*pts[0])))
LT_LEDER, LT_FILZ, LT_KDECK = 0.5, 4.0, 0.5
DR_F, DI_F, N_F = 0.5, 2.8, 3
HUB_LT = 10.0
LHX = LT_CX + 12.65                                               # Hebelauge (Zeichnung x -12.65) ueber der Tastenmitte
_kl = [e for e in LH.query('LWPOLYLINE[layer=="2Füllung-Öffnug"]')]
_kb = sorted([(min(a for a, _ in p_), max(a for a, _ in p_), min(b for _, b in p_), max(b for _, b in p_)) for p_ in [_pts(e) for e in _kl]], key=lambda r_: r_[1] - r_[0])
_oe, _kk = _kb[0], _kb[-1]                                        # Oeffnung (kleiner), Klappe (groesser)
LY = LK_CY                                                        # Hebelebene = Klappenmitte (Oeffnung aus kammern.py)
Z_B0 = zd0                                                        # Unterseite Zwischendeckel
klx0, klx1 = _kk[0] + LHX + DX_LK, _kk[1] + LHX + DX_LK; kly0, kly1 = LY - (_kk[3] - _kk[2]) / 2, LY + (_kk[3] - _kk[2]) / 2
z_ = Z_B0
for nm_, dk_ in (("Leder", LT_LEDER), ("Filz", LT_FILZ), ("Deckel", LT_KDECK)):
    mech.append((f"Luftklappe_{nm_}", box(klx0 + (0 if nm_ == "Leder" else 0.5), klx1 - (0 if nm_ == "Leder" else 0.5), kly0 + (0 if nm_ == "Leder" else 0.5), kly1 - (0 if nm_ == "Leder" else 0.5), z_ - dk_, z_))); z_ -= dk_
T_KLL = LT_FILZ + LT_KDECK
T_LH = 3.0
Z_LH = Z_B0 - (T_KLL + LT_LEDER) - 0.05 + 3.55                     # Hebeloberkante beruehrt die Klappe knapp nicht                     # Zeichnung y -> Z = Z_LH - y (langer Arm unter der Klappe)
def lh_xz(p_): return [(a + LHX, Z_LH - b) for a, b in p_]
def platte_xz(pts, y0, y1): return cq.Workplane("XZ").polyline(pts).close().extrude(-(y1 - y0)).translate((0, y0, 0)).val()
hp_ = lh_xz(_pts([e for e in LH.query('LWPOLYLINE[layer=="2-Gelenk-makierung"]')][0], 0.05))
kr = [(c.dxf.center.x + LHX, Z_LH - c.dxf.center.y, c.dxf.radius) for c in LH.query('CIRCLE[layer=="0"]')]
XP, ZP = [(k[0], k[1]) for k in kr if abs(k[2] - 1.0) < 0.01 and abs(k[0] - (-55.0 + LHX)) < 0.6][0]
XA, ZA_ = [(k[0], k[1]) for k in kr if abs(k[2] - 0.54) < 0.01][0]
hebel = platte_xz(hp_, LY - T_LH / 2, LY + T_LH / 2).fuse(stab((XP, LY - T_LH / 2, ZP), (XP, LY + T_LH / 2, ZP), 5.0))
hebel = hebel.cut(stab((XP, LY - 3, ZP), (XP, LY + 3, ZP), 2.05)).cut(stab((XA, LY - 3, ZA_), (XA, LY + 3, ZA_), 1.08))
mech.append(("Lufthebel", hebel))
L_FL = N_F * (DR_F + 0.15) + DR_F
W_LB = 8.0
y_f = [(LY - T_LH / 2 - 0.15 - L_FL, LY - T_LH / 2 - 0.15), (LY + T_LH / 2 + 0.15, LY + T_LH / 2 + 0.15 + L_FL)]
y_w = [(y_f[0][0] - 0.15 - W_LB, y_f[0][0] - 0.15), (y_f[1][1] + 0.15, y_f[1][1] + 0.15 + W_LB)]
bk = [p_ for p_ in [_pts(e, 0.05) for e in LH.query('LWPOLYLINE[layer=="2-Gelenke"]')] if min(a for a, _ in p_) < -40][0]
xb0, xb1 = min(a for a, _ in lh_xz(bk)), max(a for a, _ in lh_xz(bk))
x_schr = -64.0 + LHX
schr = []
for k_, (y0_, y1_) in enumerate(y_w):
    w_ = platte_xz(lh_xz(bk), y0_, y1_).fuse(box(xb0, xb1, y0_, y1_, Z_LH - 0.5, Z_B0 + 3.0))
    w_ = w_.cut(box(xb0 - 1, xb1 + 1, y0_ - 1, y1_ + 1, Z_B0, Z_B0 + 10)).cut(stab((XP, LY - 30, ZP), (XP, LY + 30, ZP), 2.05))
    for dy_ in (-2.0, 2.0):
        yc_ = (y0_ + y1_) / 2 + dy_
        sc_ = stab((x_schr, yc_, Z_LH - 2.53 - 0.9), (x_schr, yc_, Z_LH - 2.53), 2.6).fuse(stab((x_schr, yc_, Z_LH - 2.53), (x_schr, yc_, Z_B0 + T_DECKEL - 0.5), 1.4))
        w_ = w_.cut(stab((x_schr, yc_, Z_LH - 4), (x_schr, yc_, Z_B0 + 1), 1.5)); schr.append(sc_)
    mech.append((f"Lufthebel_Stuetze_{k_+1}", w_))
for j_, sc_ in enumerate(schr): mech.append((f"Lufthebel_Stuetze_Schraube_{j_+1}", sc_))
mech.append(("Lufthebel_Bolzen", stab((XP, y_w[0][0] - 0.6, ZP), (XP, y_w[1][1] + 0.6, ZP), 2.0)))
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
    x_u = XP + 15.0 + 0.6 * sg
    z_u = max(kante_oben(XP + 14.4), kante_oben(XP + 15.6)) + DR_F / 2 + 0.05 + (DR_F + 0.15 if sg > 0 else 0.0)
    y_a, y_b = LY + sg * (T_LH / 2 + DR_F / 2 + 0.05), LY - sg * (T_LH / 2 + DR_F / 2 + 0.05)
    haken = pfad([p0, (XP + rmf + 3.0, y_a, ZP + 0.5), (XP + 13.4, y_a, z_u), (x_u, y_a, z_u), (x_u, y_b, z_u), (x_u, y_b, z_u - 2.0)], DR_F)
    zb_ = Z_B0 - DR_F / 2 - 0.05
    s1 = pfad([p1, (XP + rmf + 3.0, y_out - sg * DR_F / 2, zb_), (XP + rmf + 12.0, y_out - sg * DR_F / 2, zb_)], DR_F)
    mech.append((f"Lufthebel_Feder_{k_+1}", sp.fuse(haken).fuse(s1)))
G_L, D_GS = 10.0, 6.0
gg = stab((XA, LY, ZA_ + 3.0), (XA, LY, ZA_ - G_L), D_GS).cut(box(XA - 4, XA + 4, LY - T_LH / 2 - 0.4, LY + T_LH / 2 + 0.4, ZA_ - 3.0, ZA_ + 4))
gg = gg.cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6)).cut(stab((XA, LY, ZA_ - 3.5), (XA, LY, ZA_ - G_L - 1), 3.0))
E0 = cq.Vector(XA, LY, ZA_); Pb = cq.Vector(XA, LT_CY, T_DECKEL_BASS - HUB_LT - 0.5)
_d = (Pb - E0).normalized(); _dw = math.degrees(math.atan2(_d.y, -_d.z))
gg = gg.rotate(E0, E0 + cq.Vector(1, 0, 0), _dw).cut(stab((XA, LY - 4, ZA_), (XA, LY + 4, ZA_), 1.6))
mech.append(("Lufttaste_Gabelgelenk", gg))
mech.append(("Lufttaste_Gelenkbolzen", stab((XA, LY - D_GS / 2, ZA_), (XA, LY + D_GS / 2, ZA_), 1.0)))
E = E0 + _d * 3.6
_it = [k_ for k_, (n_, _) in enumerate(teile) if n_ == "Lufttaste"][0]
mech.append(("Lufttaste_Schubstange", stab(E.toTuple(), Pb.toTuple(), 3.0).cut(teile[_it][1])))
_iz = [k_ for k_, (n_, _) in enumerate(teile) if n_ == "Zwischen_Deckel"][0]
_zd = teile[_iz][1]
for sc_ in schr: _zd = _zd.cut(sc_)
teile[_iz] = ("Zwischen_Deckel", _zd)
print(f"Luftklappe (Standard): Klappe X {klx0:.1f}-{klx1:.1f} Y {kly0:.1f}-{kly1:.1f}, Hebelachse X {XP:.1f}, Auge X {XA:.1f} ueber Taste X {LT_CX:.1f}, Schubstange um {_dw:.1f} Grad geneigt")

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
# --- Lufttaste gedrueckt (10 mm nach innen): darf das Gestaenge nicht beruehren ---
_tg = [s_ for n_, s_ in teile if n_ == "Lufttaste"][0].translate(cq.Vector(0, 0, 10.0))
_bt = _tg.BoundingBox(); tg_kol = []
for n1, s1 in mech + teile:
    if n1.startswith(("Lufttaste", "Lufthebel", "Luftklappe")): continue
    b1 = s1.BoundingBox()
    if b1.xmax < _bt.xmin - 0.5 or _bt.xmax + 0.5 < b1.xmin or b1.ymax < _bt.ymin - 0.5 or _bt.ymax + 0.5 < b1.ymin or b1.zmax < _bt.zmin or _bt.zmax + 0.5 < b1.zmin: continue
    v_ = s1.intersect(_tg).Volume()
    if v_ > 0.001: tg_kol.append((n1, round(v_, 2)))
print("Lufttaste gedrueckt (Z", round(_bt.zmin, 1), "-", round(_bt.zmax, 1), ") Kollisionen:", tg_kol)
# --- Basstasten gedrueckt (je eine, 6 mm), zusammen mit gedrueckter Lufttaste ---
# Knopf, Zylinder, Stange, Gelenkkopf: 6 mm nach innen; Schenkel, Gabel, Bolzen, Achse, Klappenhebel, Feder: um die Achse
# gedreht (Winkel = Hub / Armlaenge); Klappe offen (vereinfacht: um den Klappenhub parallel abgesenkt).
mech_d = dict(mech); n_st = 0; ged_kol = {}
fest_alle = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + hoerner + fuesse + klappen + mech
def _bbt(b1, b2, m_=0.0):
    return not (b1.xmax < b2.xmin - m_ or b2.xmax < b1.xmin - m_ or b1.ymax < b2.ymin - m_ or b2.ymax < b1.ymin - m_ or b1.zmax < b2.zmin - m_ or b2.zmax < b1.zmin - m_)
for i, a in enumerate(ach):
    th_ = math.degrees(HUB / a["armlen"])
    beweg = {}
    for nm_ in (f"Achse_{i+1}", f"Klappenhebel_{i+1}", f"Feder_Achse_{i+1}"):
        if nm_ in mech_d: beweg[nm_] = mech_d[nm_].rotate(cq.Vector(a["p0"][0], a["p0"][1], Z_ACHSE), cq.Vector(a["p1"][0], a["p1"][1], Z_ACHSE), -th_ * (1 if a["p1"][0] > a["p0"][0] else -1))
    for (r_, Xk, gx) in a["arm"]:
        n_st += 1; tag_ = f"{n_st}_{r_}_{Xk:.0f}"
        for nm_ in (f"Schenkel_Knopf_{tag_}", f"Gelenkgabel_{tag_}", f"Gelenkbolzen_{tag_}"):
            beweg[nm_] = mech_d[nm_].rotate(cq.Vector(a["p0"][0], a["p0"][1], Z_ACHSE), cq.Vector(a["p1"][0], a["p1"][1], Z_ACHSE), -th_ * (1 if a["p1"][0] > a["p0"][0] else -1))
        for nm_ in (f"Gelenkkopf_{tag_}", f"Schubstange_{r_}_{Xk:.0f}", f"Knopf_{r_}_{Xk:.0f}", f"Knopfzylinder_{r_}_{Xk:.0f}"):
            beweg[nm_] = mech_d[nm_].translate(cq.Vector(0, -HUB, 0))
    for n_, s_ in klappen:
        if n_ == a["kl"] or n_.startswith(a["kl"] + "_"): beweg[n_] = s_.rotate(cq.Vector(a["p0"][0], a["p0"][1], Z_ACHSE), cq.Vector(a["p1"][0], a["p1"][1], Z_ACHSE), -th_ * (1 if a["p1"][0] > a["p0"][0] else -1))   # Klappe sitzt am Klappenhebel
    fest_ = [(n_, s_) for n_, s_ in fest_alle if n_ not in beweg] + [("Lufttaste_gedrueckt", _tg)]
    for n1, s1 in beweg.items():
        b1 = s1.BoundingBox()
        for n2, s2 in fest_:
            if n2.startswith(("Lagerbuchse_Schraube",)) or (n1.startswith("Achse_") and n2.startswith("Lagerbuchse")): continue
            if not _bbt(b1, s2.BoundingBox()): continue
            v_ = s1.intersect(s2).Volume()
            if v_ > 0.05:
                bi_ = s1.intersect(s2).BoundingBox()
                ged_kol.setdefault(f"Taste {a['arm'][0][0]}{a['arm'][0][1]:.0f}", []).append((n1, n2, round(v_, 2), [round(bi_.xmin, 1), round(bi_.xmax, 1), round(bi_.ymin, 1), round(bi_.ymax, 1), round(bi_.zmin, 1), round(bi_.zmax, 1)]))
print("Gedrueckt-Pruefung Achsen", len(ach))
for k_, v_ in ged_kol.items(): print("GEDRUECKT", k_, v_)
# --- Kollisionen Mechanik gegen Gehaeuse/Kammern ---
gehaeuse = [("Front", front), ("Rueckteil", rueck), ("Seitenteil_links", sl), ("Seitenteil_rechts", sr), ("Bassboden", boden)] + teile + hoerner + fuesse + klappen
kol = {}
for n1, s1 in mech:
    bb1 = s1.BoundingBox()
    for n2, s2 in gehaeuse:
        bb2 = s2.BoundingBox()
        if bb1.xmax < bb2.xmin or bb2.xmax < bb1.xmin or bb1.ymax < bb2.ymin or bb2.ymax < bb1.ymin or bb1.zmax < bb2.zmin or bb2.zmax < bb1.zmin: continue
        v = s1.intersect(s2).Volume()
        if v > 0.3 and not n1.startswith(("Lagerbuchse_Schraube", "Luftklappe_Lagerbock_Schraube")): kol.setdefault(n1, []).append((n2, round(v, 1)))
for k, v in sorted(kol.items()): print("Kollision", k, v[:5])
