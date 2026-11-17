# Diskantmechanik Markus nach griffstock-alternativ1-schnitt.dxf (eingebunden von diskant.py per exec)
# Schnitt -> Modell: u = Abstand nach aussen von der Bodenaussenseite (Z = -u), Y wie Modell; X aus der Griffstock-Aufsicht.
# Hebel: Kontur aus dem Schnitt (Alu 2 mm), langer Arm bis zur Kroepfung, dort seitlich zur Klappenmitte versetzt,
# Endstueck mit Stufe auf die Klappe (Kontur aus dem Schnitt). Holzteile aus den Holzrand-Zeichnungen, auf den Hebel
# im Schnitt eingepasst (Abweichung < 0,3 mm); kein Alu ueber das Holz hinaus. Feste Achse Ø2 x 11 im Lagerpunkt.
exec(open("hebel_reg.py").read().split('\nif __name__')[0])

_tmpdoc = ezdxf.new()
def _flach(vs):
    pl = _tmpdoc.modelspace().add_lwpolyline([(x_, y_, 0, 0, b_) for x_, y_, b_ in vs], format="xyseb")
    q = [(v.x, v.y) for v in ezdxf.path.make_path(pl).flattening(0.05)]
    out = []
    for p_ in q:
        if not out or math.hypot(p_[0]-out[-1][0], p_[1]-out[-1][1]) > 1e-3: out.append(p_)
    return out
def _umkehr(v):
    # Polylinie umdrehen: Bulge gehoert zum Segment ab dem Punkt -> verschieben und Vorzeichen wechseln
    r_ = list(reversed(v)); out = []
    for i_ in range(len(r_)):
        b_ = -r_[i_ + 1][2] if i_ + 1 < len(r_) else 0.0
        out.append([r_[i_][0], r_[i_][1], b_])
    return out
def holz_ring(r):
    """Holzteil-Umriss ohne die Fraes-Einlaeufe: Bogen-Enden der beiden offenen Linien weglassen, gerade schliessen."""
    f_, lay_ = HEBEL_DXF[r]
    m_ = ezdxf.readfile(f_).modelspace()
    alu_e = max((e for e in m_ if e.dxftype() == "LWPOLYLINE" and e.dxf.layer == lay_ and len(e) > 6), key=lambda e: bbox.extents([e]).size.y)
    linien = []
    for e in m_:
        if e.dxftype() != "LWPOLYLINE" or e is alu_e or e.closed: continue
        vs = []
        for x_, y_, s_, e_, b_ in e.get_points():
            if not vs or math.hypot(x_ - vs[-1][0], y_ - vs[-1][1]) > 1e-3: vs.append([x_, y_, b_])
            else: vs[-1][2] = b_
        if len(vs) < 3: continue
        linien.append(vs)
    if len(linien) < 2: return None
    from shapely.geometry import LineString as _LS, Point as _P
    from shapely.ops import polygonize as _pz, unary_union as _uu
    A, B = linien[0], linien[1]
    def _ende(v, anfang): return (v[0][0], v[0][1]) if anfang else (v[-1][0], v[-1][1])
    # Paarung der Enden: A-Ende mit naechstem B-Ende
    if math.dist(_ende(A, False), _ende(B, True)) > math.dist(_ende(A, False), _ende(B, False)):
        B = [[x_, y_, 0.0] for x_, y_, b_ in reversed(B)] if all(abs(v[2]) < 0.01 for v in B) else _umkehr(B)
    la, lb = _LS(_flach(A)), _LS(_flach(B))
    ip = la.intersection(lb)
    ipts = [ip] if ip.geom_type == "Point" else list(getattr(ip, "geoms", []))
    def _kreuzt(pa, pb):
        return any(min(math.dist((q.x, q.y), pa), math.dist((q.x, q.y), pb)) < 8.0 for q in ipts if q.geom_type == "Point")
    def _ohne_einlauf_anfang(v):
        for i_ in range(min(2, len(v) - 2)):
            if abs(v[i_][2]) > 0.1: return v[i_ + 1:]
        return v
    def _ohne_einlauf_ende(v):
        n_ = len(v)
        for i_ in range(n_ - 2, max(n_ - 4, 0), -1):
            if abs(v[i_][2]) > 0.1:
                w = v[:i_ + 1]; w[-1] = [w[-1][0], w[-1][1], 0.0]; return w
        return v
    zu = []
    # Ende A-Ende / B-Anfang
    if not _kreuzt(_ende(A, False), _ende(B, True)):
        A = _ohne_einlauf_ende(A); B = _ohne_einlauf_anfang(B); zu.append((_ende(A, False), _ende(B, True)))
    if not _kreuzt(_ende(A, True), _ende(B, False)):
        A = _ohne_einlauf_anfang(A); B = _ohne_einlauf_ende(B); zu.append((_ende(A, True), _ende(B, False)))
    netz = _uu([_LS(_flach(A)), _LS(_flach(B))] + [_LS(list(z)) for z in zu])
    polys = list(_pz(netz))
    if not polys: return None
    pg = max(polys, key=lambda g: g.area)
    # Holzteil umschliesst das Alu in seinem Bereich vollstaendig (Alu mittig im Kanal, nicht sichtbar)
    from shapely.geometry import Polygon as _PG, box as _bx
    alu_pg = _PG(hebel_profil(r)[0]).buffer(0)
    x0_, y0_, x1_, y1_ = pg.bounds
    if r == 3: pg = pg.union(alu_pg.intersection(_bx(x0_, y0_, x1_, y1_)).buffer(0.3).buffer(-0.3))   # nur Reihe 3: Holz zu knapp gezeichnet
    if pg.geom_type != "Polygon": pg = max(pg.geoms, key=lambda g: g.area)
    pg = pg.simplify(0.02)
    out = []
    for x_, y_ in list(pg.exterior.coords)[:-1]:
        if not out or math.hypot(x_ - out[-1][0], y_ - out[-1][1]) > 0.05: out.append((float(x_), float(y_)))
    if math.hypot(out[0][0] - out[-1][0], out[0][1] - out[-1][1]) < 0.05: out.pop()
    return out

def prisma_uy(pts_uy, x0, dicke):
    return (cq.Workplane("YZ").polyline([(Y, -u) for u, Y in pts_uy]).close()
            .extrude(dicke).translate((x0 - dicke / 2, 0, 0)).val())
# --- Klappenmitten aus dem Boden ---
klappen_m, _ges = [], []
for ins in bm.query('INSERT[layer=="2Füllung-Öffnug"]'):
    p0 = (round(ins.dxf.insert.x, 2), round(ins.dxf.insert.y, 2))
    if p0 in _ges: continue
    _ges.append(p0)
    cs = [(v.dxf.center.x + BXO, v.dxf.center.y) for v in ins.virtual_entities() if v.dxftype() == "CIRCLE"]
    klappen_m.append((sum(c[0] for c in cs) / len(cs), sum(c[1] for c in cs) / len(cs)))
# --- Hebelplaetze aus der Griffstock-Aufsicht ---
gm = ezdxf.readfile("gs.dxf").modelspace()
from collections import defaultdict as _dd
_lr = _dd(list)
for e in gm:
    if e.dxf.layer != "Lagerbock_aufsicht": continue
    b = bbox.extents([e])
    if b.has_data and abs(b.size.x - 5.0) < 0.2 and (abs(b.size.y - 14.2) < 0.3 or abs(b.size.y - 11.3) < 0.3):
        _lr[round(b.center.y, 0)].append(b.center.x + BXO)
_ly = sorted(_lr.keys(), reverse=True)                 # 14er-Reihe = aeussere Knopfreihe = Reihe 1
HEBEL_X = {i + 1: sorted(_lr[y]) for i, y in enumerate(_ly)}
# --- feste Zuordnung Knopfreihe r -> Stimmstockreihe r, reihenweise in X-Reihenfolge (keine Ueberkreuzung);
#     die zusaetzliche (14.) Taste der aeusseren Reihe nimmt die letzte Des-Klappe ---
kl_r = {r: sorted([k for k in klappen_m if _reihe(k[1]) == r], key=lambda k: k[0]) for r in (1, 2, 3, 4)}
zuordnung = []
# Ausnahme: der zusaetzliche Knopf der aeusseren Reihe sitzt am Ende mit dem kleinsten X und geht auf die Klappe der
# 4. Reihe an diesem Ende (von der Spielerseite gesehen: Knopf 14 -> Klappe 12)
des_x = kl_r[4].pop(0)
h1 = HEBEL_X[1]
zuordnung.append((1, h1[0], des_x, 4))
for x, k in zip(h1[1:], kl_r[1]): zuordnung.append((1, x, k, 1))
print(f"Ausnahme: Klappe X {des_x[0]:.1f} (4. Reihe) an aeusserer Knopfreihe (Hebel X {h1[0]:.1f})")
for r in (2, 3):
    for x, k in zip(HEBEL_X[r], kl_r[r]): zuordnung.append((r, x, k, r))
# Reihe 4: 11 Knoepfe ohne Luecke; Fenster mit der kleinsten Summe der seitlichen Versaetze zu den 11 Klappen
_h4, _k4 = HEBEL_X[4], kl_r[4]
_i0 = 1                                  # Knoepfe X 90-290 (Plaetze 2..12 der Aufsicht)
for x, k in zip(_h4[_i0:_i0 + len(_k4)], _k4): zuordnung.append((4, x, k, 4))
print("Reihe 4: Knoepfe X", round(_h4[_i0]), "-", round(_h4[_i0 + len(_k4) - 1]), "Versatz zur Klappe", [round(k[0] - x, 1) for x, k in zip(_h4[_i0:], _k4)])

T_ALU, T_HOLZ = 2.0, 4.8          # Holzteil fuellt die 5-mm-Nut mit minimalem Spiel, Alu mittig im Kanal
U_KL = T_POLSTER + T_KLAPPE                            # Klappenoberseite (u)
def schnitt_y(pts, y0, y1):
    q = [p for p in pts if y0 <= p[1] <= y1]
    return (min(p[0] for p in q), max(p[0] for p in q)) if q else None
knopf_u = sorted(KNOPF, key=lambda q: -q[0])          # aeusserster Knopf = Reihe 1
D_KNOPF = 15.0                                        # Stufe 16 mm tief: gedrueckt 0,5 mm Luft zur Stufenkante
# --- Griffstock-Stufen: je Knopfreihe eine ebene Stufe senkrecht zur Knopfachse, von Reihe zu Reihe genau 4 mm;
#     Knopf 3,5 mm ueber seiner Stufe (Hub). a = entlang Knopfachse, b = entlang der Stufe (waechst mit u)
_Lk = math.hypot(*KNOPF_RICHT); KD = (KNOPF_RICHT[0] / _Lk, KNOPF_RICHT[1] / _Lk); KE = (KD[1], -KD[0])
def _ab(p): return (p[0] * KD[0] + p[1] * KD[1], p[0] * KE[0] + p[1] * KE[1])
def _uy(a, b): return (a * KD[0] + b * KE[0], a * KD[1] + b * KE[1])
_A0 = _ab(knopf_u[0])[0] - 3.5                        # Stufe der aeusseren Reihe (Reihe 1)
A_STUFE = {r: _A0 + 4.0 * (r - 1) for r in (1, 2, 3, 4)}
_kb = {r: _ab(knopf_u[r - 1])[1] for r in (1, 2, 3, 4)}
knopf_u = [_uy(A_STUFE[r] + 3.5, _kb[r]) for r in (1, 2, 3, 4)]   # Knoepfe genau 3,5 ueber ihrer Stufe
U_VD_IN = H_VD - T_VD                                 # Verdeck-Innenkante (u)
def _b_bei_u(a, u): return (u - a * KD[0]) / KE[0]
def terrasse(versenkt=None, tiefe=0.0):
    """Stufenlinie von innen (u = Verdeck-Innenkante) nach aussen; Stufen ausser 'versenkt' um 'tiefe' abgesenkt."""
    pts = []
    # innerste Stufe wieder in voller Hoehe bis zur Verdeck-Innenkante (keine abgesenkte Flaeche)
    _a4 = A_STUFE[4] - (tiefe if (versenkt is not None and versenkt != 4) else 0.0)
    b4_0 = _b_bei_u(_a4, U_VD_IN)
    for r in (4, 3, 2, 1):
        a = A_STUFE[r] - (tiefe if (versenkt is not None and r != versenkt) else 0.0)
        b0 = b4_0 if r == 4 else (_kb[r + 1] + _kb[r]) / 2
        b1 = _kb[r] + 9.0 if r == 1 else (_kb[r] + _kb[r - 1]) / 2
        pts += [_uy(a, b0), _uy(a, b1)]
    return pts
TERRASSE = terrasse()
from shapely.geometry import Polygon as _SPGc
def _clip(pts):
    a4 = _ab(pts[0])[0]; a1 = _ab(pts[-1])[0]
    return _SPGc([_uy(a4, _b_bei_u(a4, -50))] + pts + [_uy(a1, _kb[1] + 200), _uy(-300, _kb[1] + 200), _uy(-300, _b_bei_u(a4, -50))]).buffer(0)
FRONT_CLIP = _clip(TERRASSE)
from shapely import affinity as _aff0
from shapely.ops import unary_union as _uu0
# Vorderkanten der Stufen: an jeder Stufenkante bleibt ein 1,5 mm breiter Steg durchgehend geschlossen
def _nasen():
    out = []
    for r_hi in (4, 3, 2):                                   # hoehere (innere) Stufe, Kante zur naechst tieferen
        b_k = (_kb[r_hi] + _kb[r_hi - 1]) / 2; a_k = A_STUFE[r_hi]
        out.append(_SPGc([_uy(a_k - 5.5, b_k - 1.5), _uy(a_k + 1.0, b_k - 1.5), _uy(a_k + 1.0, b_k), _uy(a_k - 5.5, b_k)]))
    return out
NASEN = None
def holz_geclippt(r, ring):
    """Holzteil reicht als Steg durch den Griffstock bis zum Knopf: Oberkante genau 3,5 mm ueber der eigenen Stufe
    (Knopfunterseite); ein zu kurzer Steg wird verlaengert, ein zu langer gekuerzt. Der Griffstock wird ausgenommen."""
    g = _SPGc(trafo(r, ring)).buffer(0)
    a_top = A_STUFE[r] + 3.5
    a_max = max(_ab(p)[0] for p in g.exterior.coords)
    if a_max < a_top:
        steg = g.intersection(_clip_a(A_STUFE[r] - 1.0))
        g = _uu0([g] + [_aff0.translate(steg, KD[0] * t, KD[1] * t) for t in [0.25 * i for i in range(1, int((a_top - a_max) / 0.25) + 2)]])
    g = g.difference(_clip_a(a_top))
    # wo das Holz bis 2 mm unter einer Stufe liegt, reicht es bis knapp (0,5 mm) ueber die Stufe -> in den Durchbruechen sichtbar
    _fc_up = _clip(terrasse(99, -0.5))
    _band = _fc_up.difference(FRONT_CLIP.buffer(-2.0))
    _hoch = _uu0([_aff0.translate(g, KD[0] * t, KD[1] * t) for t in [0.25 * i for i in range(1, 9)]])
    g = g.union(_hoch.intersection(_band)).difference(_clip_a(a_top)).buffer(0.01).buffer(-0.01)
    # nirgends mehr als 0,5 mm ueber die Griffstockoberflaeche, ausser dem Steg zum Knopf
    _steg = _SPGc([_uy(-300, _kb[r] - 5.0), _uy(300, _kb[r] - 5.0), _uy(300, _kb[r] + 5.0), _uy(-300, _kb[r] + 5.0)])
    g = g.intersection(_fc_up.union(_steg)).buffer(0)
    g = max(g.geoms, key=lambda q: q.area) if g.geom_type != "Polygon" else g
    g = g.simplify(0.02).buffer(0)
    g = g.difference(_uu0(_nasen()).buffer(0.3)).buffer(0)      # Hebel unter der Stufenkante etwas gekuerzt
    return max(g.geoms, key=lambda q: q.area) if g.geom_type != "Polygon" else g
def _clip_a(a0):                                      # Halbebene a >= a0 (zur Knopfseite)
    return _SPGc([_uy(a0, -500), _uy(a0, 500), _uy(a0 + 300, 500), _uy(a0 + 300, -500)])
def steg_zone(r, ring):                               # b-Bereich, in dem das Holz durch die eigene Stufe nach oben geht
    b0 = -500 if r == 4 else (_kb[r + 1] + _kb[r]) / 2       # nur innerhalb der eigenen Stufe
    b1 = 500 if r == 1 else (_kb[r] + _kb[r - 1]) / 2
    g = holz_geclippt(r, ring).intersection(_clip_a(A_STUFE[r] - 0.5)).intersection(
        _SPGc([_uy(-300, b0), _uy(300, b0), _uy(300, b1), _uy(-300, b1)]))
    import shapely as _shp
    bs = [_ab(p)[1] for p in _shp.get_coordinates(g)] or [_kb[r] - 2.4, _kb[r] + 2.4]
    return _SPGc([_uy(-300, min(bs) - 0.6), _uy(300, min(bs) - 0.6), _uy(300, max(bs) + 0.6), _uy(-300, max(bs) + 0.6)])
def alu_eff(r, ring):                                 # Alu-Kontur nach dem Konturfraesen (nicht ueber das Holz hinaus)
    alu = _SPGc(HEBEL[r]).buffer(0)
    if not ring: return alu
    xw0 = min(q[0] for q in ring)
    zone = _SPGc(trafo(r, [(xw0, -200), (400, -200), (400, 200), (xw0, 200)])).buffer(0)
    g = alu.difference(zone).union(alu.intersection(zone).intersection(holz_geclippt(r, ring)))
    return max(g.geoms, key=lambda q: q.area) if g.geom_type != "Polygon" else g
print("Stufen a:", {r: round(A_STUFE[r], 2) for r in A_STUFE}, "Knoepfe:", [tuple(round(v, 1) for v in k) for k in knopf_u])
DREHPUNKTE = {}                                       # Gruppe -> (X, Y, Z) Lagerpunkt, Drehachse = X
from shapely.geometry import Polygon as _SPG, LineString as _SLS
def arm_band(r, y):
    """u-Bereich des langen Arms der Reihe r bei Y = y (Schnittkontur)."""
    g = _SPG(HEBEL[r]).buffer(0).intersection(_SLS([(-50, y), (40, y)]))
    segs = [g] if g.geom_type == "LineString" else list(getattr(g, "geoms", []))
    segs = [sg for sg in segs if sg.length > 0.5]
    if not segs: return None
    sg = min(segs, key=lambda q: q.bounds[0])
    return sg.bounds[0], sg.bounds[2]
def arm_seg(xa, ya, xb, yb, band):
    u0, u1 = band
    pl = [(xa - T_ALU/2, ya), (xa + T_ALU/2, ya), (xb + T_ALU/2, yb), (xb - T_ALU/2, yb)]
    return cq.Workplane("XY").polyline(pl).close().extrude(u1 - u0).translate((0, 0, -u1)).val()
def arm_durchgehend(r, ya, yk, pfad):
    """Langer Arm als ein Stueck: Hoehe (u) verlaeuft stetig vom Hebelquerschnitt bei ya bis zur Klappe;
    Unterseite liegt flach auf der Klappenoberseite (u = U_KL), Oberseite laeuft keilfoermig zur Spitze zu.
    In der Aufsicht wird das Profil abschnittsweise entlang des Pfades geschert (schraege Abwinkelungen)."""
    bd = arm_band(r, ya - 0.5) or arm_band(r, ya + 1)
    if r == 1 and ya < 40: bd = arm_band(1, 36.5)
    y0 = ya - 0.3
    ytip = yk + 4.0
    if r == 4:          # innere Reihe: ohne Stufe, stetig keilfoermig zur Klappe
        yc0 = max(yk - 12.0, ya + 6.0)
        prof = _SPG([(bd[0], y0), (U_KL, yc0), (U_KL, ytip), (U_KL + 2.5, ytip), (bd[1], y0)])
    else:               # Reihen 1-3: Arm hoch (Freiraum bei gedrueckter Taste), Z-foermige Stufe hinunter auf die Klappe
        yz0 = max(yk - 18.0, ya + 2.0); yz1 = yz0 + 7.0
        prof = _SPG([(bd[0], y0), (bd[0], yz0), (U_KL, yz1), (U_KL, ytip), (U_KL + 2.0, ytip),
                     (U_KL + 3.5, yz1), (bd[1], yz0 + 1.5), (bd[1], y0)])
    pts = [(y0, pfad[[p[0] for p in pfad].index(ya)][1])] + [p for p in pfad if p[0] > ya] 
    if pts[-1][0] < ytip: pts.append((ytip, pts[-1][1]))
    out = []
    for (ya_, xa_), (yb_, xb_) in zip(pts, pts[1:]):
        if yb_ - ya_ < 0.05: continue
        st = prof.intersection(_SPG([(-50, ya_ - 0.01), (100, ya_ - 0.01), (100, yb_ + 0.01), (-50, yb_ + 0.01)]))
        if st.is_empty or st.area < 0.05: continue
        st = max(st.geoms, key=lambda g: g.area) if st.geom_type != "Polygon" else st
        sol = prisma_uy([(float(a), float(b)) for a, b in list(st.exterior.coords)[:-1]], xa_, T_ALU)
        kk = (xb_ - xa_) / (yb_ - ya_)
        if abs(kk) > 1e-6:
            sol = sol.transformGeometry(cq.Matrix([[1.0, kk, 0.0, -kk * ya_], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0], [0.0, 0.0, 0.0, 1.0]]))
        out.append(sol)
    return out
# Fuehrung in der Aufsicht: Reihen 1-3 nebeneinander mit schraegen Abwinkelungen (keine 90 Grad),
# gemeinsames Versatzfenster Y 50..80; Reihe 4 darf unter den anderen durch.
FUEHR = {1: dict(lane=+17.0, y_lane=(58.0, 88.0), y_fin=(140.0, 152.0)),
         2: dict(lane=0.0, y_lane=(58.0, 88.0), y_fin=None),
         3: dict(lane=0.0, y_lane=(58.0, 88.0), y_fin=None),
         4: dict(lane=0.0, y_lane=(43.0, 49.0), y_fin=None)}   # Reihe 4 unter den anderen (u < 13) seitlich verkroepft, vor der Klappe
PFADE = []                                            # fuer die Kollisionspruefung: (Reihe, [(Y, X)...])
for n_, (r, x, k, kr) in enumerate(zuordnung):
    up, Yp = lager_von(r)
    G = f"Hebel_R{r}_{n_+1:02d}"
    DREHPUNKTE[G] = (x, Yp, -up)
    alu_s = HEBEL[r]
    v = k[0]
    # Pfad in der Aufsicht (Y, X): gerade bis ya, dann schraeg abgewinkelt zur Klappenmitte
    if kr != r:            # Ausnahme: aeussere Reihe zur Des-Klappe der 4. Reihe
        ya = 36.0
        pfad = [(Yp, x), (ya, x), (46.0, v), (k[1], v)]
    else:
        f = FUEHR[r]
        lane = v + f["lane"]
        if r == 1 and not any(abs(k3[0] - (v + 11.5)) < 4.0 for k3 in kl_r[3]):
            lane = v + 5.0      # kein gruener Hebel dazwischen: enger am blauen (Platz fuer die Verdeckstuetze)
        ya, yb = f["y_lane"]
        pfad = [(Yp, x), (ya, x), (yb, lane)]
        if f["y_fin"]:
            yf0, yf1 = f["y_fin"]
            if yf0 - yb > 0.5: pfad.append((yf0, lane))
            pfad.append((yf1, v))
        pfad.append((k[1], v))
    PFADE.append((r if kr == r else 1, pfad))
    hb = prisma_uy(alu_s, x, T_ALU).cut(box(x - 5, x + 5, ya, 400, -300, 50))
    teile_arm = arm_durchgehend(r, ya, k[1], pfad)
    ende = None
    # kein Alu ueber das Holzteil hinaus
    ring = holz_ring(r)
    if ring:
        _hc = holz_geclippt(r, ring)                                     # buendig mit der eigenen Stufe, sonst darunter
        holz_m = [(float(a_), float(b_)) for a_, b_ in list(_hc.exterior.coords)[:-1]]
        xw0 = min(q[0] for q in ring)
        zone_m = trafo(r, [(xw0, -200), (400, -200), (400, 200), (xw0, 200)])
        zone = prisma_uy(zone_m, x, T_ALU + 2); holzp = prisma_uy(holz_m, x, T_ALU + 2)
        hb = hb.cut(zone).fuse(hb.intersect(zone).intersect(holzp))
        teile.append((f"{G}/Holzteil", prisma_uy(holz_m, x, T_HOLZ), f"Holz_R{r}"))
    for s_ in teile_arm: hb = hb.fuse(s_)
    teile.append((f"{G}/Hebel_Alu", hb, f"Alu_R{r}"))
    km = min(KL_TEILE, key=lambda t_: math.hypot(t_[0][0] - k[0], t_[0][1] - k[1]))
    teile.append((f"{G}/Klappe", km[1], "Klappe")); teile.append((f"{G}/Klappenpolster", km[2], "Leder"))
    teile.append((f"{G}/Achse", cq.Workplane("YZ").center(Yp, -up).circle(1.0).extrude(11.0).translate((x - 5.5, 0, 0)).val(), "Stahl"))
    ku, kY = knopf_u[r - 1]
    d_ = KNOPF_RICHT; L_ = math.hypot(*d_)
    knopf = cq.Solid.makeCylinder(D_KNOPF / 2, 3.0, cq.Vector(x, kY, -ku), cq.Vector(0, d_[1] / L_, -d_[0] / L_))
    teile.append((f"{G}/Knopf", knopf, "Knopf"))
print("Hebel:", len(zuordnung))
# Kollisionspruefung in der Aufsicht (Reihen 1-3 untereinander), Mindestabstand Mitte-Mitte 3,5 mm
def _x_bei(pf, y):
    for (ya_, xa_), (yb_, xb_) in zip(pf, pf[1:]):
        if ya_ <= y <= yb_: return xa_ + (xb_ - xa_) * (y - ya_) / max(yb_ - ya_, 1e-9)
    return None
_kn = 0; _min = 99
for yy in [30 + 0.5 * i for i in range(300)]:
    xs_ = sorted(x_ for r_, pf in PFADE if r_ != 4 for x_ in [_x_bei(pf, yy)] if x_ is not None)
    for a_, b_ in zip(xs_, xs_[1:]):
        _min = min(_min, b_ - a_)
        if b_ - a_ < 3.5: _kn += 1
print(f"Aufsicht Reihen 1-3: kleinster Abstand Mitte-Mitte {_min:.1f} mm, zu knappe Stellen {_kn}")
_m4 = 99.0
for yy in [49.5 + 0.5 * i for i in range(57)]:
    x4 = [x_ for r_, pf in PFADE if r_ == 4 for x_ in [_x_bei(pf, yy)] if x_ is not None]
    xo = [x_ for r_, pf in PFADE if r_ != 4 for x_ in [_x_bei(pf, yy)] if x_ is not None]
    for a_ in x4:
        for b_ in xo: _m4 = min(_m4, abs(a_ - b_))
print(f"Reihe 4 gegen andere (Y 49,5-78): kleinster Abstand Mitte-Mitte {_m4:.1f} mm")
# --- Griffstock aus dem Schnitt: massiv links/rechts bis zu den Hebeln, dazwischen ausgenommen (Hebel frei),
#     hinten 3-mm-Abdeckung; Kunststoffteil (schraeg, Lagerpunkte), Dreieckleiste so lang wie das Kunststoffteil
GS_X0, GS_X1 = 0.0, 380.0
exec(open("gs_geom.py").read())            # HUELLE (Griffstock-Querschnitt mit Stufen), KPG/kst_rechteck/kst_s (Kunststoffteil)
# Huellform neu: Rueckseite, Boden, saubere Stufen (TERRASSE), Stirnseite
_st2 = sorted(GRIFFSTOCK, key=lambda q: -q[0])[:2]
_st_o = max(_st2, key=lambda q: q[1]); _st_u = min(_st2, key=lambda q: q[1])
_ro_top = (0.0, _linie_u(_ro0, _ro1, 0.0))
HUELLE = _P([_ro_top, (0.0, max(q[1] for q in KUNSTSTOFF))] + TERRASSE + [_st_o, _st_u]).buffer(0)
HUELLE = _gr(HUELLE.union(KST_E).union(DRE_E).union(GS_E.intersection(FRONT_CLIP)).buffer(0.01).buffer(-0.01).simplify(0.02))   # alte Schnittkontur nie ueber die Stufen
from shapely.geometry import Polygon as _PG2
from shapely import affinity as _aff
import json as _js0
def _pts(pg): return [(float(a_), float(b_)) for a_, b_ in list(pg.exterior.coords)[:-1]]
# Kunststoffteil nach CamBam "Kunstoffunterteil-5-fraesen-neuer4": Platte 290 lang, X = CB-x + 335;
# CB-y (0 .. -43,26 quer) -> s entlang der Lagerseite: s = -y - 1,72 (Nutmitten auf die Lagerpunkte des Schnitts)
KST_X0, KST_X1 = 55.0, 333.0          # von hinten gesehen rechts (kleines X) 10 mm, links 2 mm gekuerzt
_cb = {eval(k): v for k, v in _js0.load(open("kst_cb.json")).items()}
def _s(y): return -y - 1.72
kst = prisma_uy(_pts(KST_E), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0)
dre = prisma_uy(_pts(DRE_E), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0)       # Dreieckleiste so lang wie das Kunststoffteil
schl_k = []
for (op, y), xs in _cb.items():
    if op == "Ausnehmung-H":      # Hebeldurchbruch 5,41 x 14,64, durchgehend
        q, b_ = kst_rechteck(_s(y) - 7.32, _s(y) + 7.32, -1.0, 20.0), 5.41
    elif y < -20.0 and y > -30.0 or y > -10.0:   # Zusatztaschen 3,41 x 4,65, 3 tief
        q, b_ = kst_rechteck(max(_s(y) - 2.325, 1.2), _s(y) + 2.325, -1.0, 3.0), 3.41   # Vorderkante bleibt geschlossen
    else: continue
    for x in sorted(set(xs)):
        schl_k.append(prisma_uy(q, x + 335.0, b_))
for r_ in (1, 2, 3, 4):           # Lagernut 1,5 breit bis zur Laengsbohrung Ø2 (Achse), daneben je eine Entlastungsnut 1,5 tief
    up, Yp = lager_von(r_); s_ = kst_s(lager_von(r_))
    d_ = (up - _a[0]) * _n[0] + (Yp - _a[1]) * _n[1]
    schl_k.append(prisma_uy(kst_rechteck(s_ - 0.75, s_ + 0.75, -1.0, d_), 190.0, 300.0))
    schl_k.append(cq.Workplane("YZ").center(Yp, -up).circle(1.0).extrude(300.0).translate((40.0, 0, 0)).val())

# Entlastungsnuten 1,5 breit, 1,5 tief wie im CamBam: je Achspaar aussen 2,6 neben der Achse und eine gemeinsame in der Mitte
for r1_, r2_ in ((1, 2), (3, 4)):
    s1_, s2_ = sorted((kst_s(lager_von(r1_)), kst_s(lager_von(r2_))))
    for sm_ in (s1_ - 2.6, (s1_ + s2_) / 2, s2_ + 2.6):
        schl_k.append(prisma_uy(kst_rechteck(sm_ - 0.75, sm_ + 0.75, -1.0, 1.5), 190.0, 300.0))
# Griffstock: ein Teil. Links/rechts bis zu den Hebeln massiv (voller Querschnitt mit Stufen),
# dazwischen ausgenommen bis auf die 3-mm-Abdeckung hinten -> Hebel frei
_hx = sorted(x_ for r_ in HEBEL_X for x_ in HEBEL_X[r_]) + [x_ + 335.0 for (op, y), xs in _cb.items() if op == "Ausnehmung-H" for x_ in xs]
XM0, XM1 = KST_X0, KST_X1              # Griffstock-Enden massiv bis an Kunststoffteil/Keil (keine Vertiefungen mehr)
# Griffstock massiv ueber die ganze Laenge; im Hebelbereich je Hebel eine Nute (Hebel samt Hub geschwenkt), die nach hinten
# durchgeht; hinten 3-mm-Abdeckung, die links, rechts und an der Bodenseite in einer Nute (1,5 breit, 2 tief) sitzt
from shapely.ops import unary_union as _uu
_ro_a, _ro_b = _ro0, _ro1                                    # Rueckwand-Aussenseite
def _band(d0, d1):                                           # Streifen parallel zur Rueckseite, d0..d1 mm nach innen (+Y)
    p, q = (-5.0, _linie_u(_ro_a, _ro_b, -5.0)), (130.0, _linie_u(_ro_a, _ro_b, 130.0))
    return _PG2([(p[0], p[1] + d0), (q[0], q[1] + d0), (q[0], q[1] + d1), (p[0], p[1] + d1)])
def _u_zw(u0, u1): return _PG2([(u0, -80), (u1, -80), (u1, 80), (u0, 80)])
# Stirnseite (Daumenleiste): leichte Mulde fuer den Daumen, Kanten gerundet
_st = sorted(GRIFFSTOCK, key=lambda q: -q[0])[:2]                  # (119.9, 1.8), (116.2, -17.4)
_sm = ((_st[0][0] + _st[1][0]) / 2, (_st[0][1] + _st[1][1]) / 2)
_sl = math.dist(_st[0], _st[1]); _sd = ((_st[0][0] - _st[1][0]) / _sl, (_st[0][1] - _st[1][1]) / _sl)
_sn = (_sd[1], -_sd[0]) if _sd[1] > 0 else (-_sd[1], _sd[0])       # nach aussen (+u)
if _sn[0] < 0: _sn = (-_sn[0], -_sn[1])
_T_MULDE, _R_MULDE = 1.5, 17.0
from shapely.geometry import Point as _PT
_mulde = _PT(_sm[0] + (_R_MULDE - _T_MULDE) * _sn[0], _sm[1] + (_R_MULDE - _T_MULDE) * _sn[1]).buffer(_R_MULDE, 64)
_U_RUND = max(q[0] for q in GRIFFSTOCK) - 12.0
_aussen = HUELLE.intersection(_u_zw(_U_RUND, 200)).difference(_mulde).buffer(-1.5, join_style=1).buffer(1.5, join_style=1)
HUELLE_G = _gr(HUELLE.intersection(_u_zw(-10, _U_RUND + 3.0)).difference(_mulde).union(_aussen).buffer(0.01).buffer(-0.01).simplify(0.02))
# Abdeckung hinten: vom Boden bis 2 mm vor die Stirnseite; innen 1 mm Luft, ringsum 3 mm Auflage
_U_ABD = min(q[0] for q in _st) - 5.0
_abd = _band(-1.0, 3.0).intersection(HUELLE_G)
_abd_platte = _abd.intersection(_u_zw(-1.0, _U_ABD))
_feder = _band(0.75, 2.25).intersection(HUELLE_G).intersection(_u_zw(-1.0, _U_ABD + 2.0))
# Deckel links, rechts und an der oberen langen Seite um 3 mm groesser: dort liegt er auf (3-mm-Auflage, ohne Luft)
_abd_gross = _gr(_abd.intersection(_u_zw(-1.0, _U_ABD + 3.0)))
abdeckung = prisma_uy(_pts(_abd_gross), (XM0 + XM1) / 2, XM1 - XM0 + 6.0)
gk = prisma_uy(_pts(HUELLE_G), (GS_X0 + GS_X1) / 2, GS_X1 - GS_X0)
gk = gk.cut(kst).cut(dre)
gk = gk.cut(prisma_uy(_pts(_gr(_abd_gross.buffer(0.05))), (XM0 + XM1) / 2, XM1 - XM0 + 6.1))
_abd_luft = _band(-1.0, 3.0 + 1.0 / KD[1]).intersection(HUELLE_G).intersection(_u_zw(-1.0, _U_ABD))   # ca. 1 mm Luft zwischen Deckel und Griffstock
gk = gk.cut(prisma_uy(_pts(_gr(_abd_luft.buffer(0.05))), (XM0 + XM1) / 2, XM1 - XM0 + 0.1))
# im Bereich von Dreieckleiste und Kunststoffteil ganz offen bis zur Hoehe der Verdeck-Innenkante
U_VD_INNEN = H_VD - T_VD
# Schnitt im Hebelbereich nach griffstock-alternativ2-schnitt: oben (aussen) voller Querschnitt bis zur Abdeckung,
# unten nur ein keilfoermiger Vorderteil unter der innersten Stufe; dahinter hohl (Hebel frei), hinten die Abdeckung
_a_back = _ab(_ro0)[0]
A_ABD_INNEN = _a_back + 3.0 * KD[1]                  # Innenseite der 3-mm-Abdeckung
_b_abd_ende = _ab((_U_ABD, _linie_u(_ro0, _ro1, _U_ABD)))[1]
B_HOHL_OBEN = _kb[3] + 1.9                            # bis hierher hohl (Schnitt: ca. 2 mm ueber der 3. Knopfreihe)
B_KEIL_UNTEN, KEIL_T, B_KEIL_KNICK = _kb[4] - 25.2, 15.8, _kb[4] - 8.0
A4_ = A_STUFE[4]
_Q = _SPGc([_uy(*p) for p in [(A4_, B_KEIL_UNTEN), (A4_ - KEIL_T, B_KEIL_KNICK), (A4_ - KEIL_T, B_HOHL_OBEN),
                                 (A_ABD_INNEN, B_HOHL_OBEN), (A_ABD_INNEN, _b_abd_ende), (-50, _b_abd_ende), (-50, 400),
                                 (200, 400), (200, B_KEIL_UNTEN)]]).buffer(0)
KORPUS_M = _gr(HUELLE_G.intersection(_Q))
_hohl = HUELLE_G.difference(_Q)
gk = gk.cut(prisma_uy(_pts(_gr(_hohl.buffer(0.02))), (XM0 + XM1) / 2, XM1 - XM0))
if _hohl.geom_type != "Polygon":
    for g_ in list(_hohl.geoms)[1:]: gk = gk.cut(prisma_uy(_pts(g_.buffer(0.02)), (XM0 + XM1) / 2, XM1 - XM0))
print(f"Griffstock-Schnitt: Keil unten b {B_KEIL_UNTEN:.1f}, hohl bis b {B_HOHL_OBEN:.1f}, Abdeckung bis b {_b_abd_ende:.1f}")
def _schwenk(r, poly):
    up, Yp = lager_von(r); ku, kY = knopf_u[r - 1]
    th = 3.6 / math.hypot(ku - up, kY - Yp)                 # 3,5 mm Hub am Knopf
    return _uu([_PG2(_rot(r, list(poly.exterior.coords), th * i / 6.0)).buffer(0) for i in range(7)])
def _nach_hinten(poly):                                     # Nute geht bis zur Abdeckung durch
    return _uu([_aff.translate(poly, 0, -d) for d in range(0, 70, 1)]).intersection(HUELLE_G).difference(_band(-1.0, 3.0)).simplify(0.3)
NUT = {}
for r in (1, 2, 3, 4):
    ring = holz_ring(r)
    alu_pg = alu_eff(r, ring).buffer(0.3)
    holz_pg = holz_geclippt(r, ring).buffer(0.3) if ring else None
    # Oberflaeche sauber: nur im Bereich des eigenen Holzstegs geht die Nute durch die Stufe, sonst bleibt >= 1 mm stehen
    _innen = FRONT_CLIP.buffer(-0.2)                  # Oberflaeche bleibt geschlossen (Haut >= 0,2 mm), Hebel unveraendert
    _zone = steg_zone(r, ring) if ring else None
    def _sauber(g):
        if g is None: return None
        return g.intersection(_innen).union(g.intersection(_zone)) if _zone is not None else g.intersection(_innen)
    na = _nach_hinten(_schwenk(r, alu_pg)).intersection(_innen)      # Alu-Nute nie durch die Oberflaeche
    if holz_pg is not None:                            # Durchbruch ueber dem Holzhebel bis zur Stufenflaeche (bis 6 mm Tiefe):
        _tmax = 16.0 if r == 3 else 6.0              # 3. Reihe liegt tiefer: Nute bis oben offen
        holz_pg = _uu([_aff.translate(holz_pg, KD[0] * t, KD[1] * t) for t in [0.5 * i for i in range(int(_tmax * 2) + 1)]])   # Stufe bleibt hoch, Hebel sichtbar
    nh = _nach_hinten(_schwenk(r, holz_pg)).difference(_uu(_nasen())) if holz_pg is not None else None
    na = na.difference(_uu(_nasen()))
    NUT[r] = (na, nh)
schl_g = []
for r, x, k, kr in zuordnung:
    na, nh = NUT[r]
    for g_, b_ in ((na, T_ALU + 0.6), (nh, T_HOLZ + 0.6)):
        if g_ is None or g_.is_empty: continue
        for pg_ in (g_.geoms if g_.geom_type != "Polygon" else [g_]):
            if pg_.area > 0.5: schl_g.append(prisma_uy(_pts(pg_.simplify(0.3)), x, b_))
for w_ in schl_g: gk = gk.cut(w_)
# Durchbrueche 10 x 12 (Schnitt alternativ2, "Ausnehmungen eckig") in jeder Stufe um den Holzsteg, durch bis in den Hebelkanal
for r, x, k, kr in zuordnung:
    a_ = A_STUFE[r]; b_ = _kb[r]
    rr = _SPGc([_uy(A_ABD_INNEN + 0.5, b_ - 6.0), _uy(a_ + 2.0, b_ - 6.0), _uy(a_ + 2.0, b_ + 6.0), _uy(A_ABD_INNEN + 0.5, b_ + 6.0)])
    gk = gk.cut(prisma_uy(_pts(_gr(rr.intersection(HUELLE_G.buffer(1.0)))), x, 10.0))
teile.append(("Griffstock_Abdeckung", abdeckung, "Griffstock"))
print(f"Griffstock massiv X 0-{XM0:.0f} und {XM1:.0f}-380, Kunststoff Durchbrueche {sum(len(set(v)) for (o, y), v in _cb.items() if o == 'Ausnehmung-H')}")
# Seitenform links/rechts nach der Griffstock-Aufsicht (Diskantteil-top): Stirnseiten zum Knopfende hin eingezogen
_m2 = ezdxf.readfile("gs2.dxf").modelspace()
_u_max = max(q[0] for q in GRIFFSTOCK)
def _seite(e):
    q = [(v.x, v.y) for v in ezdxf.path.make_path(e).flattening(0.1)]
    y_hi = max(p[1] for p in q); y_lo = min(p[1] for p in q)
    return [(p[0], (y_hi - p[1]) / (y_hi - y_lo) * (_u_max + 2)) for p in q]
_kont = [e for e in _m2.query("LWPOLYLINE") if e.dxf.layer == "Diskantteil-top" and len(e) == 5]
_r = max(_kont, key=lambda e: max(v[0] for v in e.get_points()))      # rechte Aussenkontur (x bis 0)
_l = min(_kont, key=lambda e: min(v[0] for v in e.get_points()))      # linke Aussenkontur (x bis -380)
def _schnitt_kontur(e, rechts):
    q = _seite(e)
    pts = [((380.0 + x_) if rechts else (x_ + 380.0), -u_) for x_, u_ in q]   # (X, Z)
    aussen = 420.0 if rechts else -40.0
    poly = pts + [(aussen, pts[-1][1]), (aussen, pts[0][1])]
    return cq.Workplane("XZ").polyline(poly).close().extrude(-400).translate((0, -100, 0)).val()
seitenschnitt = [_schnitt_kontur(_r, True), _schnitt_kontur(_l, False)]
gk = schneide(gk, seitenschnitt); dre = schneide(dre, seitenschnitt)
teile.append(("Griffstock_Koerper", gk, "Griffstock"))
# --- Verschraubung: Blechschrauben 2,6 mm, Halbrundkopf mit Kreuzschlitz (Kopf Ø5, 2 hoch) ---
D_BS, DK_BS, HK_BS = 2.6, 5.0, 2.0
def _bs(p0, dv, laenge, name):
    """Kopf liegt bei p0 auf, Schaft in Richtung dv; liefert Bohrung (Ø2,6) als Werkzeug."""
    dv = dv.normalized()
    xd = cq.Vector(1, 0, 0) if abs(dv.x) < 0.9 else cq.Vector(0, 1, 0)
    loc = cq.Location(cq.Plane(origin=p0, xDir=xd, normal=-dv))
    teile.append((name, schraube(0, 0, 0, laenge, DK_BS, D_BS, HK_BS).moved(loc), "Schraube"))
    return cq.Solid.makeCylinder(D_BS / 2, laenge, p0, dv)
def _v(u, Y, x=0.0): return cq.Vector(x, Y, -u)
_nv = cq.Vector(0, _n[1], -_n[0])                 # senkrecht zur Kunststoffoberflaeche, ins Teil
_cbA = sorted(set(x_ + 335.0 for (op, y), xs in _cb.items() if op == "Ausnehmung-H" and abs(y + 11.32) < 0.1 for x_ in xs))
_cbD = sorted(set(x_ + 335.0 for (op, y), xs in _cb.items() if op == "Ausnehmung-H" and abs(y + 35.94) < 0.1 for x_ in xs))
# 1) Kunststoffteil -> Holzdreieck: Mitte des Kunststoffteils, senkrecht zur Oberflaeche
_sm = kst_s(_P1) + 21.65
_pm = (_a[0] + _sm * _t[0], _a[1] + _sm * _t[1])
loch_k = [_bs(_v(*_pm, x=x_), _nv, 8.5 + 9.0, f"Schraube_Kunststoff_Dreieck_{x_:.0f}") for x_ in _cbA[0::2]]
# 2) Holzdreieck von unten durch Verstaerkungsleiste und Boden
_yb = 11.5
loch_b = [_bs(cq.Vector(xs_, _yb, T_BODEN + 5.0), cq.Vector(0, 0, -1), T_BODEN + 5.0 + 9.0, f"Schraube_Boden_Dreieck_{i_+1}")
          for i_, xs_ in enumerate((62.0, 95.0, 135.0, 175.0, 215.0, 255.0, 295.0, 326.0))]
# 3) kurze Schrauben am unteren Ende des Kunststoffteils in den Boden, gleicher Winkel wie 1)
_su = kst_s(_P1) + 4.0
_pu = (_a[0] + _su * _t[0], _a[1] + _su * _t[1])
_lu = (_pu[0] + 3.5) / (-_n[0])                   # bis 3,5 mm in den Boden (4 mm)
loch_u = [_bs(_v(*_pu, x=x_), _nv, _lu, f"Schraube_Kunststoff_Boden_{x_:.0f}") for x_ in _cbD[0::2]]
for w_ in schl_k + loch_k + loch_u: kst = kst.cut(w_)
for w_ in loch_k + loch_b: dre = dre.cut(w_)
print(f"Kunststoffteil Volumen {kst.Volume():.0f} mm3, Schrauben {len(loch_k)} + {len(loch_b)} + {len(loch_u)}, Holzdreieck Tiefe an Schraube 1: {_pm}")
teile.append(("Kunststoffteil_Lagerblock", kst, "Kunststoff"))
teile.append(("Dreieckleiste", dre, "Griffstock"))
verst = schneide(verst, loch_b)
teile.append(("Verstaerkungsleiste_innen", verst, "Holz"))
import json as _js
open("drehpunkte.json", "w").write(_js.dumps(DREHPUNKTE))
