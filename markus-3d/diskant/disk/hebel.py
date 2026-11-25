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
if globals().get("INSTR") is not None:                 # Instrument-Modus: Klappen und Knoepfe aus der Tastenliste
    _dy2 = 6.5 if INSTR["choerig"] == 2 else 0.0          # 2-choerig: Klappenmitte zwischen den Loechern (0 / 13)
    klappen_m = [(k_["x"], k_["y"] + _dy2) for R_ in INSTR["reihen"].values() for k_ in R_["klappen"]]
    HEBEL_X = {r_: list(R_["knoepfe"]) for r_, R_ in INSTR["reihen"].items()}
# --- feste Zuordnung Knopfreihe r -> Stimmstockreihe r, reihenweise in X-Reihenfolge (keine Ueberkreuzung);
#     die zusaetzliche (14.) Taste der aeusseren Reihe nimmt die letzte Des-Klappe ---
kl_r = {r: sorted([k for k in klappen_m if _reihe(k[1]) == r], key=lambda k: k[0]) for r in (1, 2, 3, 4)}
zuordnung = []
if globals().get("INSTR") is not None:
    # allgemein: je Knopfreihe ab Platz knopf_ab; "extra = z" -> erster Knopf treibt die erste Klappe der Reihe z
    _frei = {r_: list(kl_r[r_]) for r_ in kl_r}
    _plan = {}
    for r_, R_ in sorted(INSTR["reihen"].items()):
        _kn = HEBEL_X[r_][R_["knopf_ab"] - 1:]
        if R_["extra"]:
            k_ = _frei[R_["extra"]].pop(0); zuordnung.append((r_, _kn[0], k_, R_["extra"]))
            print(f"Ausnahme: Klappe X {k_[0]:.1f} ({R_['extra']}. Reihe) an Knopfreihe {r_} (Hebel X {_kn[0]:.1f})")
            _kn = _kn[1:]
        _plan[r_] = _kn
    for r_ in sorted(_plan):
        if len(_plan[r_]) < len(_frei[r_]): print(f"ACHTUNG Reihe {r_}: {len(_frei[r_])} Klappen, nur {len(_plan[r_])} Knoepfe")
        for x, k in zip(_plan[r_], _frei[r_]): zuordnung.append((r_, x, k, r_))
    kl_r = {r_: v_ for r_, v_ in _frei.items()}
else:
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
T_KLAPPE4, T_POLSTER4 = 0.6, 2.0 + 0.6                 # 21.11. traditionell: Reihe 4 Decke 0,6, Filz 2 + Leder 0,6
U_KL4 = T_POLSTER4 + T_KLAPPE4
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
    zone = _SPGc(trafo(r, [(xw0, -200), (L + 20, -200), (L + 20, 200), (xw0, 200)])).buffer(0)
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
    # 19.11. (Johann): flacher Endteil auf der Klappe nur 2 mm hoch, der lange Arm laeuft bis dorthin auf halbe Hoehe zu
    h_ = bd[1] - bd[0]
    if r == 4 and globals().get("TRAD"):
        # 21.11. (Johann, traditionell): Klappe Reihe 4 flach (0,6 Decke + 2 Filz + 0,6 Leder = 3,2); der Schenkel zur Klappe
        # wird in der Materialhoehe verringert und laeuft spitz aus (Spitze 0,5 mm) - liegt so tiefer unter den Armen 1/2
        _u4 = U_KL4
        yc0 = max(yk - 12.0, ya + 6.0)
        # Oberkante gleich nach dem Griffstock auf 3 mm ueber der Unterkante (gedrueckt frei unter den Armen 1/2)
        prof = _SPG([(bd[0], y0), (_u4, yc0), (_u4, ytip), (_u4 + 0.5, ytip), (bd[0] + 3.0, y0 + 1.0), (bd[1], y0)])
    elif r == 4:        # innere Reihe: ohne Stufe, keilfoermig zur Klappe
        yc0 = max(yk - 12.0, ya + 6.0)
        prof = _SPG([(bd[0], y0), (U_KL, yc0), (U_KL, ytip), (U_KL + 2.0, ytip), (U_KL + 2.0, yc0),
                     (max(U_KL + 2.0, bd[0] + h_ / 2), yc0 - 1.0), (bd[1], y0)])
    else:               # Reihen 1-3: Arm hoch (Freiraum bei gedrueckter Taste), Z-foermige Stufe hinunter auf die Klappe
        yz0 = max(yk - 18.0, ya + 2.0); yz1 = yz0 + 7.0
        oben_ = max(bd[0] + h_ / 2, bd[0] + 2.0)
        prof = _SPG([(bd[0], y0), (bd[0], yz0), (U_KL, yz1), (U_KL, ytip), (U_KL + 2.0, ytip),
                     (U_KL + 2.0, yz1), (oben_, yz0 + 1.5), (bd[1], y0)])
    pts = [(y0, pfad[[p[0] for p in pfad].index(ya)][1])] + [p for p in pfad if p[0] > ya] 
    if pts[-1][0] < ytip: pts.append((ytip, pts[-1][1]))
    # aufeinanderfolgende Abschnitte gleicher Richtung zusammenfassen (sonst stossen gleich laufende Prismen nur
    # flaechig aneinander, und die Vereinigung verliert Stuecke - bei Marie fehlte so der lange Schenkel der Reihe 1)
    # fast gleiche X-Lagen (< 0,05 mm) angleichen: beinahe deckungsgleiche Flaechen bringen die Vereinigung zum Scheitern
    for i_ in range(1, len(pts)):
        if abs(pts[i_][1] - pts[i_ - 1][1]) < 0.05: pts[i_] = (pts[i_][0], pts[i_ - 1][1])
    _m = [pts[0]]
    for q_ in pts[1:]:
        if len(_m) >= 2:
            (y1_, x1_), (y2_, x2_) = _m[-2], _m[-1]
            if y2_ - y1_ > 1e-9 and q_[0] - y2_ > 1e-9 and abs((x2_ - x1_) / (y2_ - y1_) - (q_[1] - x2_) / (q_[0] - y2_)) < 1e-6:
                _m[-1] = q_; continue
        _m.append(q_)
    pts = _m
    global ARM_PTS_LETZTE
    ARM_PTS_LETZTE = list(pts)
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
if globals().get("INSTR") is not None and INSTR["choerig"] == 2:
    # 2-choerig (Marie): Klappenteilung = Knopfteilung, Hebel laufen fast gerade -> keine Ausweichspuren wie beim Markus
    FUEHR = {r_: dict(lane=0.0, y_lane=(58.0, 88.0) if r_ != 4 else (43.0, 49.0), y_fin=None) for r_ in (1, 2, 3, 4)}
# ---- Johanns Schnitt (jp_schnitt.dxf, 16.11.): Holzteile/Alu/Knopf je Reihe, Kammern, Teilung Teil 2 ----
import schnitt_trafo as _ST
from shapely.geometry import Polygon as _JPoly
from shapely.ops import unary_union as _juu
JP = os.path.exists("jp_schnitt.dxf")
if JP:
    _jm = ezdxf.readfile("jp_schnitt.dxf").modelspace(); _om = ezdxf.readfile("Schnitt_Hebel_4_Layer.dxf").modelspace()
    def _jpg(m, lay, alle=False):
        pg = [_JPoly([_ST.modell((x_, y_)) for x_, y_, *_r in e.get_points()]).buffer(0) for e in m.query(f'LWPOLYLINE[layer=="{lay}"]')]
        return _juu(pg) if alle else (pg[0] if pg else None)
    JP_HOLZ = {r: _jpg(_jm, f"Schnitt{r}_Hebel_R{r}_Ansicht_Holzteil") for r in (1, 2, 3, 4)}
    JP_ALU = {r: _jpg(_jm, f"Schnitt{r}_Hebel_R{r}_Ansicht_Hebel_Alu") for r in (1, 2, 3, 4)}
    _ALU_ALT = {r: _jpg(_om, f"Schnitt{r}_Hebel_R{r}_Ansicht_Hebel_Alu") for r in (1, 2, 3, 4)}
    JP_ALU_ZONE = {}
    for r in (1, 2, 3, 4):
        d_ = JP_ALU[r].symmetric_difference(_ALU_ALT[r])
        d_ = [g for g in (d_.geoms if hasattr(d_, "geoms") else [d_]) if g.area > 0.3]
        JP_ALU_ZONE[r] = _juu(d_).buffer(1.0) if d_ else None
    JP_KNOPF_D = {}
    for r in (1, 2, 3, 4):
        a_, b_ = _jpg(_om, f"Schnitt{r}_Hebel_R{r}_Ansicht_Knopf").centroid, _jpg(_jm, f"Schnitt{r}_Hebel_R{r}_Ansicht_Knopf").centroid
        JP_KNOPF_D[r] = (b_.x - a_.x, b_.y - a_.y)          # (du, dY)
    JP_KAMMER = {1: _jpg(_jm, "kammer 1 reihen", True), 2: _jpg(_jm, "kammer 2 reihe"), 3: _jpg(_jm, "kammer 3reihe"), 4: _jpg(_jm, "kammer 4reihe")}
    def _jpts(pg):
        pg = pg.buffer(0).simplify(0.02)
        pg = max(pg.geoms, key=lambda q: q.area) if hasattr(pg, "geoms") else pg
        return [(float(a_), float(b_)) for a_, b_ in list(pg.exterior.coords)[:-1]]
    def _dm(pts): return _JPoly([_ST.modell(p_) for p_ in pts]).buffer(0)
    # Teil 2 (Reihen 3/4): bis x -118 (2 mm Boden hinter den Kammern 3/4 bei x -120), oben bis y 100,2 (deckt Kammer 2 teilweise ab)
    JP_TEIL2 = _dm([(-134.2, 46.9), (-118.0, 67.9), (-118.0, 100.2), (-170.0, 100.2), (-170.0, 46.9)])
    JP_OBEN = _dm([(-170.0, 100.2), (-118.0, 100.2), (-118.0, 93.9), (-80.0, 93.9), (-80.0, 220.0), (-170.0, 220.0)])
    # Kammern: Reihe 1 lange, Reihe 2 kurze Schlitze im oberen Koerper, beide unten offen (Kamm);
    # Reihe 3 lange, Reihe 4 kurze Schlitze in Teil 2 mit Boden; Reihen 1/2 beruehren Teil 2 nicht
    JP_KAMMER[1] = _dm([(-105.8, 132.1), (-121.8, 132.1), (-121.8, 116.5), (-125.9, 116.5), (-125.8, 100.2), (-118.0, 100.2), (-118.0, 93.9), (-105.8, 93.9)])
    JP_KAMMER[2] = _dm([(-125.8, 114.5), (-105.8, 114.5), (-105.8, 93.9), (-125.8, 93.9)])
    for r in (3, 4):      # Kamm unten offen: Schlitze 3/4 durch die schraege Unterseite von Teil 2
        JP_KAMMER[r] = JP_KAMMER[r].union(_dm([(-134.6, 46.4), (-120.0, 65.7), (-112.0, 56.0), (-134.6, 36.0)])).buffer(0)
    JP_BEREICH = {1: JP_OBEN, 2: JP_OBEN, 3: JP_TEIL2, 4: JP_TEIL2}
    for r in (1, 2, 3, 4): JP_KAMMER[r] = JP_KAMMER[r].intersection(JP_BEREICH[r].buffer(0.5)).buffer(0)
    for r in (3, 4):      # 19.11.: Schlitze 3/4 nach oben (weg vom Deckel) sicher offen - 3 mm ueber die Aussenflaeche hinaus
        JP_KAMMER[r] = _juu([JP_KAMMER[r]] + [_aff0.translate(JP_KAMMER[r], 0, t_ * 0.5) for t_ in range(1, 7)]).buffer(0)
    print("JP-Schnitt: Knopf-Verschiebung", {r: tuple(round(v, 2) for v in d) for r, d in JP_KNOPF_D.items()},
          "Alu-Zonen", {r: round(z.area, 1) if z else 0 for r, z in JP_ALU_ZONE.items()})
if os.environ.get("DEBUG_LAGER"):
    for r_ in (1, 2, 3, 4):
        u_, Y_ = lager_von(r_); pg_ = _SPG(HEBEL[r_]).buffer(0); bb_ = pg_.bounds
        print("LAGER", r_, "u", round(u_, 2), "Yp", round(Y_, 2), "Hebel u", round(bb_[0], 1), round(bb_[2], 1), "Y", round(bb_[1], 1), round(bb_[3], 1),
              "Kammer", [round(v_, 1) for v_ in JP_KAMMER[r_].bounds] if JP else "")
        for Yq in range(int(bb_[1]), int(Y_) + 12, 3):
            g_ = pg_.intersection(_SLS([(-80, Yq), (120, Yq)]))
            st_ = [q_.bounds for q_ in (getattr(g_, "geoms", None) or [g_]) if not q_.is_empty]
            print("   Y", Yq, [(round(q_[0], 1), round(q_[2], 1)) for q_ in st_])
    raise SystemExit("DEBUG_LAGER")
# 22.11.2026: Griffstock "kolben" - verdeckt, Kunststoffzylinder unter den Knoepfen, Hebel nur Alu (kolben_geo.py)
_GSART0 = ((INSTR.get("griffstock", "") if globals().get("INSTR") is not None else "") or os.environ.get("GRIFFSTOCK", ""))
_WINKEL = _GSART0 == "winkel"       # 22.11.: wie kolben, aber Alu-Winkel 46753 mit einer Achse fuer alle Hebel (winkel_geo.py)
_KOLBEN = _GSART0 in ("kolben", "kolben12", "winkel")
# 22.11. 11:45 (Johann): "kolben12" = Reihen 1/2 verdeckt mit Kolben, Reihen 3/4 offen mit Kanaelen und Holzhebeln (Teil 2 wie bisher)
KOLBEN_R = (1, 2) if _GSART0 == "kolben12" else (1, 2, 3, 4)
_KOLBEN_ALLE = _KOLBEN and 3 in KOLBEN_R
if _KOLBEN: exec(open("kolben_geo.py").read())
if _WINKEL: exec(open("winkel_geo.py").read())
# Knopfversatz (Johann 18.11.): Reihen 1/2 Knopf + Knopfschenkel um eine halbe Teilung nach +X, der Alu-Hebel bekommt
# zwischen u 45 und 62 (Knopfschenkel unterhalb der Kammern, oberhalb Teil 2 / Reihen 3/4) eine seitliche Z-Kroepfung
KV = {r_: (INSTR["reihen"][r_].get("knopf_versatz", 0.0) if globals().get("INSTR") is not None and r_ in INSTR["reihen"] else 0.0) for r_ in (1, 2, 3, 4)}
KV_U0, KV_U1 = 45.0, 62.0
def _kroepfen(s_, dx):
    """Koerper s_ (Alu) ueber u (= -z): unten unveraendert, oben um dx in X, dazwischen geschert."""
    if not dx: return s_
    W = 1e4
    unten = s_.intersect(box(-W, W, -W, W, -KV_U0, W))                 # u < KV_U0  (z > -u0)
    oben = s_.intersect(box(-W, W, -W, W, -W, -KV_U1)).translate(cq.Vector(dx, 0, 0))
    mitte = s_.intersect(box(-W, W, -W, W, -KV_U1, -KV_U0))
    k_ = dx / (KV_U1 - KV_U0)                                             # X' = X + k*(u - u0), u = -z
    mat = cq.Matrix([[1.0, 0.0, -k_, -k_ * KV_U0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]])
    erg = unten
    for t_ in (mitte.transformGeometry(mat) if mitte.Volume() > 1e-6 else None, oben if oben.Volume() > 1e-6 else None):
        if t_ is not None: erg = erg.fuse(t_)
    return erg.clean()
if any(KV.values()): print("Knopfversatz:", {r_: v_ for r_, v_ in KV.items() if v_}, "Kroepfung u", KV_U0, "-", KV_U1)
D_BUCHSE = 6.0                                        # 19.11.: Messingbuchse aussen 6, innen 2 (Achse Ø2), Laenge = Alu 2 mm
WAND_BUCHSE = []
HEBEL_DATEN = []                                      # fuer die Abwicklungen (teile_dxf.py)
ARM_PTS_LETZTE = []
def min_wand_buchse(r_, hb_, x_, Yp_, up_):
    """kleinste Alu-Wandstaerke rund um die 6-mm-Bohrung (Schnitt in der Hebelebene)"""
    try:
        sec_ = hb_.intersect(cq.Solid.makeBox(0.2, 60, 60, cq.Vector(x_ - 0.1, Yp_ - 30, -up_ - 30)))
        from shapely.geometry import Point as _Pt, Polygon as _PG2
        from shapely.ops import unary_union as _uu2
        pg_ = []
        for f_ in sec_.Faces():
            if abs(f_.normalAt().x) < 0.9: continue
            pg_.append(_PG2([(v_.y, v_.z) for v_ in (f_.outerWire().positionAt(k_ / 200) for k_ in range(201))]).buffer(0))
        a_ = _uu2(pg_)
        c_ = _Pt(Yp_, -up_)
        lo_, hi_ = 0.0, 10.0
        for _ in range(25):
            m_ = (lo_ + hi_) / 2
            ring_ = c_.buffer(D_BUCHSE / 2 + m_, 64).difference(c_.buffer(D_BUCHSE / 2 + 0.01, 64))
            if ring_.difference(a_).area < 1e-3: lo_ = m_
            else: hi_ = m_
        return lo_
    except Exception as ex_:
        print("Wandpruefung Buchse:", ex_); return -1.0
PFADE = []                                            # fuer die Kollisionspruefung: (Reihe, [(Y, X)...])
if _WINKEL:
    # 22.11. 16:14 (Johann): von hinten gesehen weichen die Hebel den Zylindern aus - der Freiraum zwischen den Zylindern
    # wird gleichmaessig fuer die Durchfuehrung der Hebel aufgeteilt. Je Zylinderreihe b eine Zone (u-Bereich, in dem
    # Hebel hoeherer Reihen a < b durch die Zylinder b laufen); dort liegen die durchlaufenden Hebel gleichmaessig
    # verteilt in den Luecken; unterhalb im Schlitz des Winkels, oberhalb beim eigenen Knopf.
    from shapely.geometry import LineString as _WLS2, Polygon as _WPG2
    _wH = {r_: _WPG2(HEBEL[r_]).buffer(0) for r_ in (1, 2, 3, 4)}
    _wKN = {r_: (knopf_u[r_ - 1][0] + (JP_KNOPF_D[r_][0] if JP else 0), knopf_u[r_ - 1][1] + (JP_KNOPF_D[r_][1] if JP else 0)) for r_ in (1, 2, 3, 4)}
    _wC = {r_: _WLS2([KOLBEN_UNTEN[r_], _wKN[r_]]).buffer(3.75, cap_style=2) for r_ in (1, 2, 3, 4)}
    W_ZONE = {}
    for b_ in (1, 2, 3, 4):
        z_ = [_wH[a_].intersection(_wC[b_].buffer(0.5)) for a_ in (1, 2, 3, 4) if a_ < b_]
        z_ = [g_ for g_ in z_ if not g_.is_empty]
        if z_: W_ZONE[b_] = (min(g_.bounds[0] for g_ in z_) - 0.3, max(g_.bounds[2] for g_ in z_) + 0.3)
    W_OWN_U = {r_: _wC[r_].intersection(_wH[r_]).bounds[0] - 0.3 for r_ in (1, 2, 3, 4)}
    print("Winkel: Zylinderzonen (u)", {b_: tuple(round(v_, 1) for v_ in z_) for b_, z_ in W_ZONE.items()}, "eigener Zylinder ab u", {r_: round(v_, 1) for r_, v_ in W_OWN_U.items()})
    _wlev = [(r_, x_, x_ + KV[r_]) for r_, x_, k_, kr_ in zuordnung]          # (Reihe, Hebelebene, Knopf-X)
    # 22.11. 17:48 (Johann): Fuehrung nach ital.df.dxf (Ansicht von hinten, 20 mm Knopfabstand): Reihen 3/4 gerade auf dem
    # eigenen X; Reihen 1/2 oben (Achse) +5, an der naechsttieferen Reihe +7, dann auf das eigene X (Reihe 2 zwischen den
    # Zylindern 4 und 3, Reihe 1 zwischen 3 und 2) -> oben gleichmaessige Teilung 5 mm
    # 18:11 (Johann): in jeder Engstelle die durchlaufenden Hebel gleichmaessig zwischen die beiden Zylinder verteilen
    # (welche Luecke und Reihenfolge wie in der Zeichnung); von der Lagerung nur eine Kroepfung bis zur ersten Engstelle
    _R_ZYL = 3.75
    _ZG = {4: {1: 5.0, 2: 7.0, 3: 0.0}, 3: {1: 7.0, 2: 0.0}, 2: {1: 0.0}}      # Lage nach Zeichnung (relativ zum eigenen X)
    W_POS = {}
    for b_ in (4, 3, 2):
        cyl_ = sorted(kx_ for r_, x_, kx_ in _wlev if r_ == b_)
        pas_ = [(kx_ + _ZG[b_][r_], r_, x_) for r_, x_, kx_ in _wlev if r_ < b_]
        cyx_ = [cyl_[0] - 20.0 * k_ for k_ in (4, 3, 2, 1)] + cyl_ + [cyl_[-1] + 20.0 * k_ for k_ in (1, 2, 3, 4)]   # ausserhalb gedachte Zylinder
        for i_ in range(len(cyx_) - 1):
            e0_ = cyx_[i_] + _R_ZYL
            e1_ = cyx_[i_ + 1] - _R_ZYL
            l_ = sorted(q_ for q_ in pas_ if e0_ < q_[0] < e1_)
            if not l_: continue
            if e0_ < -1e8 or e1_ > 1e8: xs_ = [q_[0] for q_ in l_]
            else:
                n_ = len(l_); g_ = (e1_ - e0_ - 2.0 * n_) / (n_ + 1)
                if g_ < 0.4: print(f"   ACHTUNG Winkel: Engstelle {e0_:.1f}-{e1_:.1f} (Zylinder Reihe {b_}) mit {n_} Hebeln - nur {g_:.2f} mm Luft")
                xs_ = [e0_ + g_ * (k_ + 1) + 2.0 * k_ + 1.0 for k_ in range(n_)]
            for (p_, r_, x_), v_ in zip(l_, xs_): W_POS.setdefault((r_, round(x_, 3)), {})[b_] = v_
        print(f"Winkel: Engstelle Zylinder Reihe {b_}: Luft Hebel-Zylinder/Hebel-Hebel min "
              f"{min((abs(v_[b_] - c_) - _R_ZYL - 1.0) for v_ in W_POS.values() if b_ in v_ for c_ in cyl_):.2f} mm")
    W_OBEN = {(r_, round(x_, 3)): (W_POS[(r_, round(x_, 3))][4] if r_ < 4 else kx_) for r_, x_, kx_ in _wlev}
    _wx_unten = sorted((W_OBEN[(r_, round(x_, 3))], r_, x_) for r_, x_, kx_ in _wlev)
    _wsm = [W_SCHLITZ_X[i_] for i_ in range(len(_wx_unten))]           # 1 + 4 x 12 + 1 Hebel zwischen den Winkel-Laschen
    W_XOFF = sum(a_[0] - b_ for a_, b_ in zip(_wx_unten, _wsm)) / len(_wx_unten)
    W_SLOT = {(r_, round(x_, 3)): W_XOFF + m_ for (v_, r_, x_), m_ in zip(_wx_unten, _wsm)}
    print("Winkel: Schlitz - Lage oben min %.2f max %.2f mm, Platten-Nullpunkt X %.2f" % (
        min(W_SLOT[(r_, round(x_, 3))] - v_ for v_, r_, x_ in _wx_unten), max(W_SLOT[(r_, round(x_, 3))] - v_ for v_, r_, x_ in _wx_unten), W_XOFF))
    W_U_SCHLITZ = 36.0
    def _w_knoten(r_, x_, kx_):
        sl_ = W_SLOT[(r_, round(x_, 3))]; ps_ = W_POS.get((r_, round(x_, 3)), {})
        kn_ = [(W_U_SCHLITZ, sl_)]
        for b_ in (4, 3, 2):
            if b_ <= r_: break
            z0_, z1_ = W_ZONE[b_]
            kn_ += [(z0_ - (1.0 if b_ == 4 else 0.3), ps_[b_]), (z1_, ps_[b_])]
        kn_.append((W_OWN_U[r_] - 0.3 if r_ < 4 else W_ZONE[4][0] - 1.0, kx_))
        kn_ = [k_ for i_, k_ in enumerate(kn_) if i_ == 0 or abs(k_[1] - kn_[i_ - 1][1]) > 1e-6 or i_ == len(kn_) - 1 or True]
        for (u0_, _a), (u1_, _b) in zip(kn_, kn_[1:]):
            if u1_ - u0_ < 2.5 and abs(_a - _b) > 0.05: print(f"   ACHTUNG Winkel R{r_} x {x_:.1f}: Uebergang u {u0_:.1f}-{u1_:.1f} nur {u1_ - u0_:.1f} mm")
        return kn_
    def _w_versatz(s_, x_, ya_, knoten_):
        """u < 36: im Schlitz, Klappenarm Y 36..ya-1 zurueck auf die Hebelebene x; u >= 36: X nach den Knoten (stueckweise geschert)"""
        W = 1e4; U0 = W_U_SCHLITZ
        dx = knoten_[0][1] - x_
        lo_ = s_.intersect(box(-W, W, -W, W, -U0, W)).translate(cq.Vector(dx, 0, 0)); hi_ = s_.intersect(box(-W, W, -W, W, -W, -U0))
        st_ = [lo_]
        gr_ = [k_ for k_ in knoten_] + [(1e3, knoten_[-1][1])]
        for (u0_, a_), (u1_, b_) in zip(gr_, gr_[1:]):
            sl_ = hi_.intersect(box(-W, W, -W, W, -u1_, -u0_))
            if sl_.Volume() < 1e-6: continue
            a_, b_ = a_ - x_, b_ - x_
            if abs(b_ - a_) < 1e-6: st_.append(sl_.translate(cq.Vector(a_, 0, 0)))
            else:
                k_ = (b_ - a_) / (u1_ - u0_)       # X' = X + a + k (u - u0), u = -z
                st_.append(sl_.transformGeometry(cq.Matrix([[1.0, 0.0, -k_, a_ - k_ * u0_], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]])))
        erg_ = None
        for t_ in st_:
            if t_ is None or t_.Volume() < 1e-6: continue
            erg_ = t_ if erg_ is None else erg_.fuse(t_)
        return erg_.clean()
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
        if r == 1 and f["lane"] != 0.0 and not any(abs(k3[0] - (v + 11.5)) < 4.0 for k3 in kl_r[3]):
            lane = v + 5.0      # kein gruener Hebel dazwischen: enger am blauen (Platz fuer die Verdeckstuetze)
        ya, yb = f["y_lane"]
        pfad = [(Yp, x), (ya, x), (yb, lane)]
        if f["y_fin"]:
            yf0, yf1 = f["y_fin"]
            if yf0 - yb > 0.5: pfad.append((yf0, lane))
            pfad.append((yf1, v))
        pfad.append((k[1], v))
    if _WINKEL:          # Klappenarm bleibt bis ya im Schlitz-X des Winkels, biegt dann zur Klappe
        _wxs0 = W_SLOT[(r, round(x, 3))]
        pfad = [(pfad[0][0], _wxs0), (pfad[1][0], _wxs0)] + pfad[2:]
    PFADE.append((r if kr == r else 1, pfad))
    hb = prisma_uy(alu_s, x, T_ALU).cut(box(x - 5, x + 5, ya, L + 20, -300, 50))
    if _WINKEL:
        _wkn = _w_knoten(r, x, x + KV[r]); _xs = _wkn[0][1]
        for _nl in W_NIETLOCH_LOKAL: hb = hb.cut(W_LOKAL(_nl, r, x))           # Nietloecher (vor dem Versatz, in der Ebene x)
        hb = hb.cut(W_LOKAL(W_FEDER[r][1], r, x))                               # Bohrung Ø1 fuer den Federschenkel
        _v0 = hb.Volume(); hb = _w_versatz(hb, x, ya, _wkn)
    teile_arm = arm_durchgehend(r, ya, k[1], pfad)
    if os.environ.get("DEBUG_ARM") and n_ in (0, 1, 15, 28, 40):
        print("DEBUG_ARM", G, "ya", ya, "yk", round(k[1], 1), "pfad", [(round(a_, 1), round(b_, 1)) for a_, b_ in pfad],
              "Stuecke", [tuple(round(v_, 1) for v_ in (q_.BoundingBox().ymin, q_.BoundingBox().ymax, q_.BoundingBox().zmin)) for q_ in teile_arm])
    ende = None
    # kein Alu ueber das Holzteil hinaus
    ring = None if (_KOLBEN and r in KOLBEN_R) else holz_ring(r)
    if ring:
        _hc = holz_geclippt(r, ring)                                     # buendig mit der eigenen Stufe, sonst darunter
        holz_m = [(float(a_), float(b_)) for a_, b_ in list(_hc.exterior.coords)[:-1]]
        xw0 = min(q[0] for q in ring)
        zone_m = trafo(r, [(xw0, -200), (L + 20, -200), (L + 20, 200), (xw0, 200)])
        zone = prisma_uy(zone_m, x, T_ALU + 2); holzp = prisma_uy(holz_m, x, T_ALU + 2)
        hb = hb.cut(zone).fuse(hb.intersect(zone).intersect(holzp))
        if JP: holz_m = _jpts(JP_HOLZ[r])
        teile.append((f"{G}/Holzteil", prisma_uy(holz_m, x + KV[r], T_HOLZ), f"Holz_R{r}"))
    _vsoll = hb.Volume() + sum(s_.Volume() for s_ in teile_arm)
    for s_ in teile_arm: hb = hb.fuse(s_)
    if len(hb.Solids()) > 1 or hb.Volume() < 0.9 * _vsoll:
        print(f"   ACHTUNG {G}: Hebel nicht zusammenhaengend ({len(hb.Solids())} Koerper, Volumen {hb.Volume():.0f} statt <= {_vsoll:.0f})")
    _zz = JP_ALU_ZONE[r].intersection(_JPoly([(-300, -300), (300, -300), (300, ya - 1.0), (-300, ya - 1.0)])) if (JP and not (_KOLBEN and r in KOLBEN_R) and JP_ALU_ZONE[r] is not None) else None
    if _zz is not None and not _zz.is_empty:      # nur im geraden Stueck (Griffstock), nicht am verkroepften Arm
        for z_ in [g for g in (_zz.geoms if hasattr(_zz, "geoms") else [_zz]) if g.geom_type == "Polygon" and g.area > 0.05]:
            hb = hb.cut(prisma_uy(_jpts(z_), x, 30.0))
            neu_ = JP_ALU[r].intersection(z_)
            for g_ in (neu_.geoms if hasattr(neu_, "geoms") else [neu_]):
                if g_.geom_type == "Polygon" and g_.area > 0.05: hb = hb.fuse(prisma_uy(_jpts(g_), x, T_ALU))
    if _WINKEL:
        if len(hb.Solids()) > 1: print(f"   ACHTUNG {G}: Winkel-Hebel {len(hb.Solids())} Koerper")
    if not _WINKEL: hb = _kroepfen(hb, KV[r])          # Winkel: Knopfversatz schon in der Kroepfung u 34-44
    # 19.11.: Messingbuchse zwischen Achse und Alu: aussen 6, Bohrung 2 (Achse), eingepresst, so breit wie das Alu (2 mm) -
    # Aussenmasse unveraendert; Alu bekommt die 6-mm-Bohrung
    _bz = lambda r_, l_: cq.Solid.makeCylinder(r_, l_, cq.Vector(x - l_ / 2, Yp, -up), cq.Vector(1, 0, 0))
    if not _WINKEL:
        hb = hb.cut(_bz(D_BUCHSE / 2, T_ALU + 2.0))
        _wand = min_wand_buchse(r, hb, x, Yp, up)
        if _wand < 0.8: print(f"   ACHTUNG {G}: um die Messingbuchse nur {_wand:.2f} mm Alu")
        WAND_BUCHSE.append(_wand)
    teile.append((f"{G}/Hebel_Alu", hb, f"Alu_R{r}"))
    HEBEL_DATEN.append({"G": G, "r": r, "x": x, "Yp": Yp, "up": up, "kv": KV[r], "ya": ya, "pts": ARM_PTS_LETZTE, "hb": hb})
    if _WINKEL:      # Ms-Lasche 0,5 nach 4685 / V3, Bohrung Ø3 auf der gemeinsamen Achse, 2 Niete
        teile.append((f"{G}/Ms_Lasche", W_LOKAL(W_LASCHE_LOKAL, r, _xs), "Messing"))
        for _i, _n in enumerate(W_NIET_LOKAL): teile.append((f"{G}/Niete_{_i + 1}", W_LOKAL(_n, r, _xs), "Messing"))
        teile.append((f"{G}/Feder", W_LOKAL(W_FEDER[r][0], r, _xs), "Stahl"))
    else: teile.append((f"{G}/Messingbuchse", _bz(D_BUCHSE / 2, T_ALU).cut(_bz(1.0, T_ALU + 1.0)), "Messing"))
    km = min(KL_TEILE, key=lambda t_: math.hypot(t_[0][0] - k[0], t_[0][1] - k[1]))
    if _WINKEL and kr in (1, 2, 3):
        # 22.11. 19:37 (Johann): Gummi 031404 V3 auf der Klappe (Reihen 1-3), Hebelende steckt im Loch 2 x 2,5;
        # Klappendeckel 1 mm niedriger, weil das Gummi mit seinem 1-mm-Fuss darunter liegt (Hebelhoehe bleibt)
        _kl = km[1].cut(box(-1e4, 1e4, -1e4, 1e4, -(U_KL + 60.0), -(U_KL - 1.0)))
        teile.append((f"{G}/Klappe", _kl, "Klappe")); teile.append((f"{G}/Klappenpolster", km[2], "Leder"))
        _gx, _gy = pfad[-1][1], k[1]
        _gp = [(-4.0, 0.0), (4.0, 0.0), (2.0, 1.0), (2.0, 4.3)] + [(2.0 * math.cos(math.pi * i_ / 24), 4.3 + 2.0 * math.sin(math.pi * i_ / 24)) for i_ in range(1, 24)] + [(-2.0, 4.3), (-2.0, 1.0)]   # Profil nach 031404_Gummi_V3.dxf
        _gu = cq.Workplane("XZ").polyline([(_gx + a_, -(U_KL - 1.0 + b_)) for a_, b_ in _gp]).close().extrude(-7.0).translate((0, _gy - 3.5, 0)).val()
        _gu = _gu.cut(box(_gx - 1.0, _gx + 1.0, _gy - 10, _gy + 10, -(U_KL - 1.0 + 3.5), -(U_KL - 1.0 + 1.0)))
        teile.append((f"{G}/Gummi_031404", _gu, "Leder"))
    else:
        teile.append((f"{G}/Klappe", km[1], "Klappe")); teile.append((f"{G}/Klappenpolster", km[2], "Leder"))
    if not _WINKEL: teile.append((f"{G}/Achse", cq.Workplane("YZ").center(Yp, -up).circle(1.0).extrude(11.0).translate((x - 5.5, 0, 0)).val(), "Stahl"))
    ku, kY = knopf_u[r - 1]
    d_ = KNOPF_RICHT; L_ = math.hypot(*d_)
    if JP: ku, kY = ku + JP_KNOPF_D[r][0], kY + JP_KNOPF_D[r][1]
    knopf = cq.Solid.makeCylinder(D_KNOPF / 2, 3.0, cq.Vector(x + KV[r], kY, -ku), cq.Vector(0, d_[1] / L_, -d_[0] / L_))
    teile.append((f"{G}/Knopf", knopf, "Knopf"))
    if _KOLBEN and r in KOLBEN_R and KOLBEN_UNTEN.get(r):                   # Kunststoffzylinder ueber das Armende bis unter den Knopf; der verjuengte
        cu_, cY_ = KOLBEN_UNTEN[r]; _lk = math.dist((cu_, cY_), (ku, kY))     # Hebelteil sitzt in einer gleich grossen Ausnehmung
        _kol = cq.Solid.makeCylinder(D_KOLBEN / 2, _lk, cq.Vector(x + KV[r], cY_, -cu_), cq.Vector(0, d_[1] / L_, -d_[0] / L_)).cut(hb)
        _thk = math.degrees(3.6 / math.hypot(ku - up, kY - Yp))      # Ausnehmung ueber den ganzen Hub (Hebel dreht im Zylinder)
        for _wk in (() if _WINKEL else (_thk / 2, _thk, -_thk / 2, -_thk)):     # Winkel: Zylinder starr mit dem Hebel (Zapfen)
            _kol = _kol.cut(hb.rotate(cq.Vector(x - 5, Yp, -up), cq.Vector(x + 5, Yp, -up), _wk))
        teile.append((f"{G}/Kolben", _kol, "Kunststoff"))
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
if WAND_BUCHSE: print(f"Messingbuchsen Ø{D_BUCHSE:g}/2: {len(WAND_BUCHSE)}, Alu-Wand rundum min {min(WAND_BUCHSE):.2f} mm, max {max(WAND_BUCHSE):.2f} mm")
# --- Griffstock aus dem Schnitt: massiv links/rechts bis zu den Hebeln, dazwischen ausgenommen (Hebel frei),
#     hinten 3-mm-Abdeckung; Kunststoffteil (schraeg, Lagerpunkte), Dreieckleiste so lang wie das Kunststoffteil
GS_X0, GS_X1 = 0.0, L
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
KST_X0, KST_X1 = 55.0, 333.0 + DL
if globals().get("INSTR") is not None:   # Kunststoffteil immer nach den Knoepfen der Tastenliste (Markus: 55/333 wie bisher)
    _kx = [x_ for r_ in HEBEL_X for x_ in HEBEL_X[r_]]
    KST_X0, KST_X1 = min(_kx) - 9.0, max(_kx) + 6.0  # 19.11.: hinter dem letzten Hebel 3 mm kuerzer (Johann)    # wie beim Markus: 10 mm vor dem ersten, 8 mm nach dem letzten Knopf         # von hinten gesehen rechts (kleines X) 10 mm, links 2 mm gekuerzt
_cb = {eval(k): v for k, v in _js0.load(open("kst_cb.json")).items()}
if globals().get("INSTR") is not None:
    # Durchbrueche/Taschen je Knopfreihe auf die Knoepfe der Tastenliste umlegen (gleiche Anzahl Ueberstand an den Enden)
    _ta = INSTR.get("knopfabstand", 20.0)
    _MK = {1: (64.99, 14), 2: (74.99, 13), 3: (80.02, 12), 4: (70.02, 13)}     # Markus-Knopfreihen (Start, Anzahl), Teilung 20
    def _umlegen(xs):
        xs = sorted(set(round(x_ + 335.0, 2) for x_ in xs))
        r_ = min(_MK, key=lambda r_: abs(((xs[0] - _MK[r_][0] + 10.0) % 20.0) - 10.0))   # Reihe ueber die Lage im 20er-Raster
        kn = sorted(HEBEL_X[r_])
        if len(xs) == 1:                                   # Einzeltasche: gleicher Platz in der Reihe (vom Ende gezaehlt)
            j_ = round((xs[0] - _MK[r_][0]) / 20.0) - (_MK[r_][1] - 1)
            return [kn[max(0, min(len(kn) - 1, len(kn) - 1 + j_))] - 335.0]
        vor = max(0, round((_MK[r_][0] - xs[0]) / 20.0))
        nach = max(0, round((xs[-1] - (_MK[r_][0] + 20.0 * (_MK[r_][1] - 1))) / 20.0))
        a_ = [kn[0] - _ta * j for j in range(vor, 0, -1)] + kn + [kn[-1] + _ta * j for j in range(1, nach + 1)]
        return [x_ - 335.0 for x_ in a_]
    _cb = {k_: _umlegen(v_) for k_, v_ in _cb.items()}
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
    schl_k.append(prisma_uy(kst_rechteck(s_ - 0.75, s_ + 0.75, -1.0, d_), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0 + 30.0))
    schl_k.append(cq.Workplane("YZ").center(Yp, -up).circle(1.0).extrude(KST_X1 - KST_X0 + 30.0).translate((KST_X0 - 15.0, 0, 0)).val())

# Entlastungsnuten 1,5 breit, 1,5 tief wie im CamBam: je Achspaar aussen 2,6 neben der Achse und eine gemeinsame in der Mitte
for r1_, r2_ in ((1, 2), (3, 4)):
    s1_, s2_ = sorted((kst_s(lager_von(r1_)), kst_s(lager_von(r2_))))
    for sm_ in (s1_ - 2.6, (s1_ + s2_) / 2, s2_ + 2.6):
        schl_k.append(prisma_uy(kst_rechteck(sm_ - 0.75, sm_ + 0.75, -1.0, 1.5), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0 + 30.0))
# Griffstock: ein Teil. Links/rechts bis zu den Hebeln massiv (voller Querschnitt mit Stufen),
# dazwischen ausgenommen bis auf die 3-mm-Abdeckung hinten -> Hebel frei
# Drehfedern (Johann 18.11.): Feder links vom Hebel (von den Knoepfen Richtung Klappen gesehen, Aussenseite oben = +X);
# dafuer Kunststoffteil, Dreieckleiste, Durchbrueche und Schrauben um KST_DX verschieben, Hebel liegen dann nicht mehr mittig
FEDER = (INSTR.get("feder", "") if globals().get("INSTR") is not None else os.environ.get("FEDER", "")) or ""
# 20.11.: traditioneller Griffstock (trad.py) - Holzhebel auf zwei durchgehenden Achsen statt Kunststoff-Lagerung
_GSART = ((INSTR.get("griffstock", "") if globals().get("INSTR") is not None else "") or os.environ.get("GRIFFSTOCK", ""))
TRAD = _GSART in ("traditionell", "kugellager")         # 21.11.: kugellager (kugel.py) nutzt dieselben Vorgaben wie traditionell
KUGEL = _GSART == "kugellager"
if TRAD: FEDER = ""
KST_DX = 1.5 if FEDER else 0.0      # Durchbruch 5,41: Hebel 2 mm mit 0,2 mm Luft an einer Seite, 3,2 mm frei fuer die Feder
if KST_DX:
    _cb = {k_: [x_ + KST_DX for x_ in v_] for k_, v_ in _cb.items()}
    KST_X0, KST_X1 = KST_X0 + KST_DX, KST_X1 + KST_DX
    kst = prisma_uy(_pts(KST_E), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0)
    dre = prisma_uy(_pts(DRE_E), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0)
    schl_k = []
    for (op, y), xs in _cb.items():
        if op == "Ausnehmung-H":
            q, b_ = kst_rechteck(_s(y) - 7.32, _s(y) + 7.32, -1.0, 20.0), 5.41
        elif y < -20.0 and y > -30.0 or y > -10.0:
            q, b_ = kst_rechteck(max(_s(y) - 2.325, 1.2), _s(y) + 2.325, -1.0, 3.0), 3.41
        else: continue
        for x in sorted(set(xs)):
            schl_k.append(prisma_uy(q, x + 335.0, b_))
    for r_ in (1, 2, 3, 4):
        up, Yp = lager_von(r_); s_ = kst_s(lager_von(r_))
        d_ = (up - _a[0]) * _n[0] + (Yp - _a[1]) * _n[1]
        schl_k.append(prisma_uy(kst_rechteck(s_ - 0.75, s_ + 0.75, -1.0, d_), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0 + 30.0))
        schl_k.append(cq.Workplane("YZ").center(Yp, -up).circle(1.0).extrude(KST_X1 - KST_X0 + 30.0).translate((KST_X0 - 15.0, 0, 0)).val())
    for r1_, r2_ in ((1, 2), (3, 4)):
        s1_, s2_ = sorted((kst_s(lager_von(r1_)), kst_s(lager_von(r2_))))
        for sm_ in (s1_ - 2.6, (s1_ + s2_) / 2, s2_ + 2.6):
            schl_k.append(prisma_uy(kst_rechteck(sm_ - 0.75, sm_ + 0.75, -1.0, 1.5), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0 + 30.0))
    print(f"Feder: Kunststoffteil/Dreieckleiste/Schrauben um {KST_DX:+.1f} mm in X verschoben")
_hx = sorted(x_ for r_ in HEBEL_X for x_ in HEBEL_X[r_]) + [x_ + 335.0 for (op, y), xs in _cb.items() if op == "Ausnehmung-H" for x_ in xs]
XM0, XM1 = KST_X0, KST_X1              # Griffstock-Enden massiv bis an Kunststoffteil/Keil (keine Vertiefungen mehr)
if any(KV.values()):                    # versetzte Knoepfe (Reihe 1/2) brauchen ihre Kammern auch ausserhalb des Kunststoffteils
    _kxv = [x_ + KV[r_] for r_, x_, k_, kr_ in zuordnung]
    XM0, XM1 = min(XM0, min(_kxv) - 9.0), max(XM1, max(_kxv) + 9.0)
if _WINKEL:                              # 22.11.: seitlich gekroepfte Winkel-Hebel brauchen innen Platz bis zur aeussersten Lage (+3 mm)
    _wall = [v_ for d_ in W_POS.values() for v_ in d_.values()] + list(W_SLOT.values())
    XM0, XM1 = min(XM0, min(_wall) - 3.0), max(XM1, max(_wall) + 3.0)
if globals().get("INSTR") is not None:   # 19.11.: beide Griffstock-Enden gleich breit (die schmaelere Seite gilt)
    _wgs = min(XM0 - GS_X0, GS_X1 - XM1)
    XM0, XM1 = GS_X0 + _wgs, GS_X1 - _wgs
    print(f"Griffstock-Enden gleich breit: je {_wgs:.1f} mm (Kunststoffteil X {KST_X0:.1f}-{KST_X1:.1f})")
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
if JP and not _KOLBEN:     # Deckel nur so weit wie die Schlitze + 3 mm Auflage (nicht bis an die Daumenleiste)
    _U_ABD = max(JP_KAMMER[r].intersection(_band(-2.0, 4.0)).bounds[2] for r in (1, 2, 3, 4) if not JP_KAMMER[r].intersection(_band(-2.0, 4.0)).is_empty) + 0.5
    print(f"Abdeckung: Schlitze bis u {_U_ABD - 0.5:.1f}, Deckel mit Auflage bis u {_U_ABD + 3.0:.1f}")
_abd = _band(-1.0, 3.0).intersection(HUELLE_G)
_abd_platte = _abd.intersection(_u_zw(-1.0, _U_ABD))
_feder = _band(0.75, 2.25).intersection(HUELLE_G).intersection(_u_zw(-1.0, _U_ABD + 2.0))
# Deckel links, rechts und an der oberen langen Seite um 3 mm groesser: dort liegt er auf (3-mm-Auflage, ohne Luft)
_abd_gross = _gr(_abd.intersection(_u_zw(-1.0, _U_ABD + 3.0)))
abdeckung = prisma_uy(_pts(_abd_gross), (XM0 + XM1) / 2, XM1 - XM0 + 6.0)
dre = dre.cut(abdeckung)              # 19.11.: Dreieckleiste lag ca. 0,05 mm im Deckel - jetzt nur Beruehrung
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
if JP: _hohl = HUELLE_G.difference(JP_TEIL2.union(JP_OBEN))      # Johann: Koerper = Teil 2 + oberer Block, dazwischen leer
if _KOLBEN:            # verdeckt: innen ganz hohl, Knopfseite 4-mm-Abdeckung mit den Stufen, oben die Daumenleiste
    _kdeck = HUELLE_G.difference(_aff.translate(HUELLE_G, -T_ABDECK * KD[0], -T_ABDECK * KD[1]))
    if _KOLBEN_ALLE:   # 22.11. 19:04 (Johann, 19:27 auch Kolben): Innenkontur der Abdeckung ohne Stufen wie im Schnitt (griffstock-alternativ DEZ 2018,
        # "4mm abdeckung"): an der untersten und obersten Stufe 4 mm, dazwischen eine gerade Schraege
        _fl = []
        _hc = list(HUELLE_G.exterior.coords)
        for _pa, _pb in zip(_hc, _hc[1:]):
            _dx, _dy = _pb[0] - _pa[0], _pb[1] - _pa[1]; _ll = math.hypot(_dx, _dy)
            if _ll < 5.0 or abs((_dx * KD[0] + _dy * KD[1]) / _ll) > 0.05: continue
            _nn = _pa[0] * KD[0] + _pa[1] * KD[1]
            if _nn < 20.0 or max(_pa[0] * KD[1] - _pa[1] * KD[0], _pb[0] * KD[1] - _pb[1] * KD[0]) < 30.0: continue    # nur die Knopfstufen
            _tt = sorted((_pa[0] * KD[1] - _pa[1] * KD[0], _pb[0] * KD[1] - _pb[1] * KD[0]))
            _fl.append((_nn, _tt[0], _tt[1]))
        _nf = sorted(set(round(f_[0], 1) for f_ in _fl), reverse=True)            # Stufenflaechen von unten (groesstes n) nach oben
        _f0 = max((f_ for f_ in _fl if round(f_[0], 1) == _nf[0]), key=lambda f_: f_[2])
        _f3 = min((f_ for f_ in _fl if round(f_[0], 1) == _nf[-1]), key=lambda f_: f_[1])
        _t1, _t2 = _f0[2] - 7.3, _f3[1] + 8.3                                    # Zeichnung: 18,5 - 11,18 bzw. 58,81 - 50,5
        _n1, _n2 = _nf[0] - T_ABDECK, _nf[-1] - T_ABDECK
        def _tn(t_, n_): return (n_ * KD[0] + t_ * KD[1], n_ * KD[1] - t_ * KD[0])
        _innen = _PG2([_tn(-300, -300), _tn(-300, _n1), _tn(_t1, _n1), _tn(_t2, _n2), _tn(400, _n2), _tn(400, -300)])
        _kdeck = HUELLE_G.difference(_innen)
        print(f"Verdeckter Griffstock: Abdeckung innen glatt - {T_ABDECK} mm an Stufe unten (n {_nf[0]:.1f}) bis t {_t1:.1f}, Schraege bis t {_t2:.1f}, dann {T_ABDECK} mm an Stufe oben (n {_nf[-1]:.1f})")
    _kdeck = _kdeck.intersection(_u_zw(U_KOLBEN_AUSTRITT, 400.0))   # unten offen: dort gehen die Alu-Hebel zu den Klappen hinaus
    U_KTOP = max(max(q_[0] for q_ in knopf_u) + D_KOLBEN / 2 + 3.0, max(g_.bounds[2] for g_ in KOLBEN_ALU.values()) + 4.0)   # 4 mm ueber dem hoechsten Armende (gedrueckt)
    _ktop = HUELLE_G.intersection(_u_zw(U_KTOP, 400.0))
    _hohl = HUELLE_G.difference(_kdeck).difference(_ktop)
    if not _KOLBEN_ALLE:   # nur Reihen 1/2 verdeckt: Teil 2 (Kanaele 3/4) wie bisher, oberer Block innen hohl bis auf Abdeckung/Daumenleiste
        _hohl = HUELLE_G.difference(JP_TEIL2.union(JP_OBEN.intersection(_kdeck.union(_ktop))))
    _hohl = max(_hohl.geoms, key=lambda g_: g_.area) if hasattr(_hohl, "geoms") else _hohl
    print(f"Kolben-Griffstock: Abdeckung {T_ABDECK} mm, Daumenleiste ab u {U_KTOP:.1f}, innen hohl {_hohl.area:.0f} mm2")
gk = gk.cut(prisma_uy(_pts(_gr(_hohl.buffer(0.02))), (XM0 + XM1) / 2, XM1 - XM0))
if _hohl.geom_type != "Polygon":
    for g_ in list(_hohl.geoms)[1:]: gk = gk.cut(prisma_uy(_pts(g_.buffer(0.02)), (XM0 + XM1) / 2, XM1 - XM0))
if os.environ.get("DEBUG_TEIL2"):
    import pickle as _pk
    _pk.dump({"HUELLE_G": HUELLE_G, "HUELLE": HUELLE, "T2": JP_TEIL2, "OBEN": JP_OBEN, "K": JP_KAMMER, "HOLZ": JP_HOLZ, "ALU": JP_ALU, "KST": KST_E, "DRE": DRE_E, "abd": _abd_gross, "GS": (GS_X0, GS_X1, XM0, XM1)}, open("/tmp/teil2.pkl", "wb"))
    print("DEBUG_TEIL2 geschrieben"); raise SystemExit(0)
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
if JP:                                                 # Kammern statt der Bewegungs-Nuten: je Hebel ein Schlitz (Holz + 0,6)
    schl_g = []
    JP_KAM_EXT = {}
    for r in (1, 2, 3, 4):     # eigener Hebel (Holz + gerades Alu), Ruhe bis gedrueckt, 0,3 Luft: Kammer notfalls ergaenzen
        _ya = FUEHR[r]["y_lane"][0] - 1.0
        _eig = JP_HOLZ[r].union(JP_ALU[r].intersection(_JPoly([(-300, -300), (300, -300), (300, _ya), (-300, _ya)]))).buffer(0.3)
        _sw = _schwenk(r, _gr(_eig)).buffer(0)
        _neu = _sw.difference(JP_KAMMER[r].buffer(0.01)).intersection(JP_BEREICH[r]).intersection(HUELLE_G)
        _neu = _juu([g for g in (_neu.geoms if hasattr(_neu, "geoms") else [_neu]) if g.geom_type == "Polygon" and g.area > 0.2])
        JP_KAM_EXT[r] = _neu
        print(f"Kammer {r}: ergaenzt um {_neu.area:.1f} mm2", [tuple(round(v, 1) for v in _neu.bounds)] if not _neu.is_empty else "")
        JP_KAMMER[r] = JP_KAMMER[r].union(_neu.buffer(0.1)).buffer(0) if not _neu.is_empty else JP_KAMMER[r]
    for r, x, k, kr in zuordnung:
        _kk = JP_KAMMER[r].intersection(HUELLE_G.buffer(2.0)).buffer(0)
        for g_ in (_kk.geoms if hasattr(_kk, "geoms") else [_kk]):
            if g_.geom_type == "Polygon" and g_.area > 0.5: schl_g.append(prisma_uy(_jpts(g_), x + KV[r], T_HOLZ + 0.6))
if _KOLBEN:
    schl_g = []
    for r, x, k, kr in zuordnung:                     # offene Kanaele nur fuer die Reihen ohne Kolben (Holzhebel)
        if r in KOLBEN_R: continue
        _kk = JP_KAMMER[r].intersection(HUELLE_G.buffer(2.0)).buffer(0)
        for g_ in (_kk.geoms if hasattr(_kk, "geoms") else [_kk]):
            if g_.geom_type == "Polygon" and g_.area > 0.5: schl_g.append(prisma_uy(_jpts(g_), x + KV[r], T_HOLZ + 0.6))
for w_ in schl_g: gk = gk.cut(w_)
if _KOLBEN:            # Durchfuehrung je Kolben durch die Abdeckung
    for n_, (r, x, k, kr) in enumerate(zuordnung):
        if r not in KOLBEN_R: continue
        ku, kY = knopf_u[r - 1]
        if JP: ku, kY = ku + JP_KNOPF_D[r][0], kY + JP_KNOPF_D[r][1]
        # 22.11. 12:23 (Johann): Durchfuehrung als Rechteck 10 x 12 mit runden Ecken (Draufsicht "Ausnehmungen eckig"),
        # 10 quer (X), 12 entlang der Stufe
        _pl = cq.Plane(origin=(x + KV[r], kY + 5 * KD[1], -(ku + 5 * KD[0])), xDir=(1, 0, 0), normal=(0, -KD[1], KD[0]))
        _ak = cq.Workplane(_pl).sketch().rect(10.0, 12.0).vertices().fillet(2.0).finalize().extrude(30.0).val()
        gk = gk.cut(_ak)
    if os.environ.get("WINKEL_GS_DUMP"):
        import pickle as _gpk
        _gpk.dump({"HUELLE_G": HUELLE_G, "hohl": _hohl, "abd": _abd_gross, "HEBEL": {r_: HEBEL[r_] for r_ in (1, 2, 3, 4)}, "LAGER": {r_: lager_von(r_) for r_ in (1, 2, 3, 4)},
                   "KNOPF": KOLBEN_KNOPF, "KD": KD}, open("/tmp/winkel_gs.pkl", "wb"))
        raise SystemExit("WINKEL_GS_DUMP")
    # Freischnitt: Knopfarm Ruhe bis gedrueckt darf Abdeckung und Deckel nicht beruehren (wo der Platz knapp ist, z. B. Reihe 1)
    _nfrei = 0
    if _KOLBEN_ALLE:  # 22.11. 19:00 (Johann, 19:27 auch Kolben): keine Freischnitte - Hebel halten 2 mm Abstand zur Griffstock-Innenseite (Ruhe und gedrueckt)
        from shapely.geometry import LineString as _WLS3
        from shapely.affinity import rotate as _wrot
        _wwand = HUELLE_G.difference(_hohl).buffer(2.0)
        _nk = 0
        for _ti, (_tn, _ts, _tf) in enumerate(teile):
            if not (_tn.startswith("Hebel_R") and _tn.endswith("/Hebel_Alu")): continue
            _r = int(_tn[7]); _L = lager_von(_r); _K = KOLBEN_KNOPF[_r]
            _tg = math.degrees(3.6 / math.hypot(_K[0] - _L[0], _K[1] - _L[1]))
            _band = _WLS3([(KOLBEN_UNTEN[_r][0] - KD[0], KOLBEN_UNTEN[_r][1] - KD[1]), _K]).buffer(D_KOLBEN / 2, cap_style=2)
            _vb = _wwand.union(_wrot(_wwand, _tg, origin=_L)).difference(_band).intersection(_u_zw(_L[0] + 12.0, 400.0))
            _vb = [g_ for g_ in (_vb.geoms if hasattr(_vb, "geoms") else [_vb]) if g_.geom_type == "Polygon" and g_.area > 0.05]
            _ns = _ts
            for g_ in _vb:
                _w = prisma_uy(_pts(g_), (GS_X0 + GS_X1) / 2, GS_X1 - GS_X0 + 20)
                if _ns.intersect(_w).Volume() > 1e-3: _ns = _ns.cut(_w); _nk += 1
            if _ns is not _ts:
                teile[_ti] = (_tn, _ns, _tf)
                for _hd in HEBEL_DATEN:
                    if _hd["G"] == _tn.split("/")[0]: _hd["hb"] = _ns
        print(f"Verdeckter Griffstock: Hebel auf 2 mm Abstand zur Innenseite gekuerzt: {_nk} Stellen")
    for _hd in ([] if _KOLBEN_ALLE else HEBEL_DATEN):
        if _hd["r"] not in KOLBEN_R: continue
        _hx, _hu, _hY = _hd["x"], _hd["up"], _hd["Yp"]
        _kkh = KOLBEN_KNOPF[_hd["r"]]; _tg = math.degrees(3.6 / math.hypot(_kkh[0] - _hu, _kkh[1] - _hY))
        _arm3 = _hd["hb"].intersect(box(-1e4, 1e4, -1e4, 1e4, -400.0, -(_hu + 8.0)))
        for _wk in (-_tg, -_tg / 2, _tg / 2, _tg):
            _a3 = _arm3.rotate(cq.Vector(_hx - 5, _hY, -_hu), cq.Vector(_hx + 5, _hY, -_hu), _wk)
            if gk.intersect(_a3).Volume() > 1e-3: gk = gk.cut(_a3); _nfrei += 1
            if abdeckung.intersect(_a3).Volume() > 1e-3: abdeckung = abdeckung.cut(_a3); _nfrei += 1
    print(f"Kolben-Griffstock: Freischnitte fuer den Hebelweg: {_nfrei}")
# Durchbrueche 10 x 12 (Schnitt alternativ2, "Ausnehmungen eckig") in jeder Stufe um den Holzsteg, durch bis in den Hebelkanal
for r, x, k, kr in ([] if JP else zuordnung):
    a_ = A_STUFE[r]; b_ = _kb[r]
    rr = _SPGc([_uy(A_ABD_INNEN + 0.5, b_ - 6.0), _uy(a_ + 2.0, b_ - 6.0), _uy(a_ + 2.0, b_ + 6.0), _uy(A_ABD_INNEN + 0.5, b_ + 6.0)])
    gk = gk.cut(prisma_uy(_pts(_gr(rr.intersection(HUELLE_G.buffer(1.0)))), x, 10.0))
teile.append(("Griffstock_Abdeckung", abdeckung, "Griffstock"))
print(f"Griffstock massiv X 0-{XM0:.0f} und {XM1:.0f}-{L:.0f}, Kunststoff Durchbrueche {sum(len(set(v)) for (o, y), v in _cb.items() if o == 'Ausnehmung-H')}")
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
    pts = [((L + x_) if rechts else (x_ + L_GRUND), -u_) for x_, u_ in q]   # (X, Z); rechte Stirnseite mit dem Gehaeuse verschoben
    aussen = L + 40.0 if rechts else -40.0
    poly = pts + [(aussen, pts[-1][1]), (aussen, pts[0][1])]
    return cq.Workplane("XZ").polyline(poly).close().extrude(-400).translate((0, -100, 0)).val()
seitenschnitt = [_schnitt_kontur(_r, True), _schnitt_kontur(_l, False)]
gk = schneide(gk, seitenschnitt); dre = schneide(dre, seitenschnitt)
if JP and not _KOLBEN_ALLE:                            # Teilung: Teil 2 (unten, Reihen 3/4) ueber die ganze Laenge, oberer Block = Koerper
    # 19.11. (2): Teil 2 = nur der vom Deckel abgewandte Teil unter dem oberen Block; Trennung eben, parallel zur Deckelflaeche,
    # etwa in halber Hoehe (beim Fraesen liegen beide Teile auf der Deckelseite). Vorderer Block mit Stufe 1/2 bleibt bei Teil 1.
    _R2 = HUELLE_G.difference(JP_OBEN)
    _lo, _hi = 0.0, 200.0
    for _ in range(40):
        _m = (_lo + _hi) / 2
        if _band(_m, 400.0).intersection(_R2).area > 1e-3: _lo = _m
        else: _hi = _m
    # 19.11. (3): Trennebene genau am Grund der Kammern 3/4 (+0,05), damit die Schlitze in Teil 2 durchgehen (Kamm mit Seitenteilen)
    def _dmin_(P):
        lo_, hi_ = -50.0, 200.0
        for _ in range(40):
            m_ = (lo_ + hi_) / 2
            if _band(-50.0, m_).intersection(P).area > 1e-4: hi_ = m_
            else: lo_ = m_
        return hi_
    _dk34 = min(_dmin_(JP_KAMMER[r_].intersection(HUELLE_G)) for r_ in (3, 4))
    _vd = []                                               # Kammergrund ist schraeg: hoechster Punkt des Grundes zaehlt
    for r_ in (3, 4):
        k_ = JP_KAMMER[r_].intersection(HUELLE_G)
        for g_ in (k_.geoms if hasattr(k_, "geoms") else [k_]):
            if g_.geom_type != "Polygon": continue
            _vd += [y_ - _linie_u(_ro_a, _ro_b, x_) for x_, y_ in g_.exterior.coords]
    _dg = max(v_ for v_ in _vd if v_ <= _dk34 + 3.0)
    # 19.11. (4): Teil 2 = urspruenglicher Kamm MIT Boden (JP_TEIL2), ueber die ganze Laenge; Trennebene genau
    # unter dessen Boden, damit im Hebelbereich kein Streifen von Teil 1 unter Reihe 3/4 stehen bleibt
    _dt2 = _dmin_(JP_TEIL2.intersection(HUELLE_G))
    D_TEIL2 = round(_dt2 - 0.01, 2)
    print(f"Griffstock Teilung: Boden Teil 2 (Kamm 3/4) bei {_dt2:.2f} mm, Kammergrund {_dk34:.2f}-{_dg:.2f} mm ueber der Deckelflaeche")
    _t2q = _R2.intersection(_band(D_TEIL2, 400.0))
    _t2q = max(_t2q.geoms, key=lambda q_: q_.area) if hasattr(_t2q, "geoms") else _t2q
    print(f"Griffstock Teilung: Hoehe ueber Deckelflaeche {_lo:.1f} mm, Trennebene bei {D_TEIL2} mm")
    if os.environ.get("DEBUG_EBENE"):
        def _dmin(P):
            lo_, hi_ = -50.0, 200.0
            for _ in range(40):
                m_ = (lo_ + hi_) / 2
                if _band(-50.0, m_).intersection(P).area > 1e-4: hi_ = m_
                else: lo_ = m_
            return hi_
        def _dmax(P):
            lo_, hi_ = -50.0, 200.0
            for _ in range(40):
                m_ = (lo_ + hi_) / 2
                if _band(m_, 400.0).intersection(P).area > 1e-4: lo_ = m_
                else: hi_ = m_
            return lo_
        for r_ in (1, 2, 3, 4):
            k_ = JP_KAMMER[r_].intersection(HUELLE_G)
            print(f"DEBUG Kammer {r_}: d {_dmin(k_):.2f} .. {_dmax(k_):.2f}")
        print("DEBUG T2 d", round(_dmin(_R2), 2), round(_dmax(_R2), 2), "OBEN d", round(_dmin(JP_OBEN.intersection(HUELLE_G)), 2), round(_dmax(JP_OBEN.intersection(HUELLE_G)), 2))
        import pickle as _pk2; _pk2.dump({"R2": _R2, "K": JP_KAMMER, "H": HUELLE_G, "t2q": _t2q, "abd": _abd_gross, "dre": DRE_E, "b0": _band(D_TEIL2 - 1.0, D_TEIL2), "band0": _band(-1.0, 0.0)}, open("/tmp/ebene.pkl", "wb"))
        raise SystemExit(0)
    _t2 = prisma_uy(_pts(_gr(_t2q)), (GS_X0 + GS_X1) / 2, GS_X1 - GS_X0 + 2.0)
    teile.append(("Griffstock_Teil2", gk.intersect(_t2), "Griffstock"))
    gk = gk.cut(_t2)
    print("Griffstock Teil 2 Volumen %.0f, Rest %.0f" % (teile[-1][1].Volume(), gk.Volume()))
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
          for i_, xs_ in enumerate((lambda q_: [v_ + KST_DX for v_ in q_])((62.0, 95.0, 135.0, 175.0, 215.0, 255.0, 295.0, 323.5)) if abs(DL) < 1e-6 else
                                    [KST_X0 + 7.0 + (KST_X1 - KST_X0 - 14.0) * j / (7 + round(DL / 40.0)) for j in range(8 + round(DL / 40.0))])]
# 3) kurze Schrauben am unteren Ende des Kunststoffteils in den Boden, gleicher Winkel wie 1)
_su = kst_s(_P1) + 4.0
_pu = (_a[0] + _su * _t[0], _a[1] + _su * _t[1])
_lu = (_pu[0] + 3.5) / (-_n[0])                   # bis 3,5 mm in den Boden (4 mm)
loch_u = [_bs(_v(*_pu, x=x_), _nv, _lu, f"Schraube_Kunststoff_Boden_{x_:.0f}") for x_ in _cbD[0::2]]
for w_ in schl_k + loch_k + loch_u: kst = kst.cut(w_)
for w_ in loch_k + loch_b: dre = dre.cut(w_)
print(f"Kunststoffteil Volumen {kst.Volume():.0f} mm3, Schrauben {len(loch_k)} + {len(loch_b)} + {len(loch_u)}, Holzdreieck Tiefe an Schraube 1: {_pm}")
teile.append(("Kunststoffteil_Lagerblock", kst, "Kunststoff"))
# 4) 19.11.: Griffstock von unten durch Verstaerkungsleiste und Boden angeschraubt, je Ende 2 Schrauben
#    (eine in Teil 1 an der Deckelseite, eine in Teil 2); Y je Teil so gewaehlt, dass der Schaft moeglichst viel Fleisch hat
GS_LOECHER = []
if JP and globals().get("INSTR") is not None:
    from shapely.geometry import LineString as _LS
    _T_GS = 12.0                                          # Einschraubtiefe in den Griffstock
    _P1 = HUELLE_G if _KOLBEN_ALLE else HUELLE_G.difference(_t2q.buffer(0.01)); _P1 = max(_P1.geoms, key=lambda q_: q_.area) if hasattr(_P1, "geoms") else _P1
    def _y_best(P):
        Pe = P.union(_aff.translate(P, -10.0, 0)).buffer(0)
        best = None
        for i_ in range(-200, 400):
            y_ = i_ * 0.25
            sg = _LS([(0.0, y_), (_T_GS, y_)])
            if not Pe.contains(sg.buffer(D_BS / 2)): continue
            d_ = Pe.boundary.distance(sg)
            if best is None or d_ > best[1]: best = (y_, d_)
        return best
    _yt = {n_: _y_best(P_) for n_, P_ in ((("Teil1", _P1),) if _KOLBEN_ALLE else (("Teil1", _P1), ("Teil2", _t2q)))}
    _xs = [GS_X0 + _wgs / 2, GS_X1 - _wgs / 2]
    _sl = [s_.BoundingBox() for n_, s_, f_ in teile if n_ == "Schraubleiste_Stimmstoecke"]
    for j_ in range(2):                                   # Schraubenkopf (Oe5) nicht unter die Schraubleiste der Stimmstoecke
        for bb_ in _sl:
            if bb_.xmin - DK_BS / 2 - 0.5 < _xs[j_] < bb_.xmax + DK_BS / 2 + 0.5:
                _xs[j_] = bb_.xmin - DK_BS / 2 - 0.5 if j_ == 1 else bb_.xmax + DK_BS / 2 + 0.5
                print(f"Griffstock-Schraube {'rechts' if j_ else 'links'} wegen Schraubleiste nach X {_xs[j_]:.1f}")
    for n_, (y_, d_) in _yt.items():
        for j_, x_ in enumerate(_xs):
            GS_LOECHER.append((x_, y_, n_, _bs(cq.Vector(x_, y_, T_BODEN + 5.0), cq.Vector(0, 0, -1), T_BODEN + 5.0 + _T_GS,
                                                f"Schraube_Griffstock_{n_}_{'links' if j_ == 0 else 'rechts'}")))
    print("Griffstock-Schrauben von unten:", [(round(x_, 1), round(y_, 1), n_) for x_, y_, n_, w_ in GS_LOECHER],
          {n_: round(d_, 1) for n_, (y_, d_) in _yt.items()}, "mm Fleisch neben dem Schaft")
    _gsl = [w_ for x_, y_, n_, w_ in GS_LOECHER]
    teile = [(n_, (schneide(s_, _gsl) if n_ in ("Griffstock_Teil2", "Griffstock_Koerper") else
                   schneide(s_, [cq.Solid.makeCylinder(1.5, T_BODEN + 2.0, cq.Vector(x_, y_, -1.0), cq.Vector(0, 0, 1)) for x_, y_, nn_, w_ in GS_LOECHER])
                   if n_ == "Diskantboden" else s_), f_) for n_, s_, f_ in teile]
# ---- Drehfeder je Hebel: Durchmesser 6 aussen, Draht 0,6, 3 Windungen, auf der Achse neben dem Hebel (+X);
#      Schenkel 1 liegt auf dem Kunststoffteil, Schenkel 2 mit U-Haken ueber die Kante des 2-mm-Alu-Hebels
D_FED, D_DRAHT, N_WDG = 6.0, 0.6, 3
def _draht(p0, p1):
    v_ = p1.sub(p0)
    return cq.Solid.makeCylinder(D_DRAHT / 2, v_.Length, p0, v_.normalized())
def _kugel(p): return cq.Solid.makeSphere(D_DRAHT / 2, p)
def _hebel_u_bei(r_, Y):
    """u-Bereich des Alu-Hebels der Reihe r_ bei Y (Querschnitt)."""
    g_ = _SPG(HEBEL[r_]).buffer(0).intersection(_SLS([(-60, Y), (60, Y)]))
    return (g_.bounds[0], g_.bounds[2]) if not g_.is_empty else None
FEDERN = []
RASTKERBE_DATEN = []
for G_, (xh_, Yp_, zp_) in sorted(DREHPUNKTE.items()):
    if not FEDER or (FEDER != "alle" and G_ != FEDER): continue
    r_ = int(G_.split("_R")[1][0]); up_ = -zp_
    rm_ = (D_FED - D_DRAHT) / 2; steig_ = D_DRAHT + 0.1
    x0_ = xh_ + T_ALU / 2 + D_DRAHT / 2 + 0.1                         # Draht 0,1 mm neben dem Hebel
    import math as _m_
    # --- Hakenpunkt am Knopfschenkel (U an der Unterkante, 8 mm vor der Achse, Oeffnung zum Blech) ---
    _pg = _SPG(HEBEL[r_]).buffer(0)
    def _u_unten(Y):                     # Unterkante (kleinstes u) des Knopfschenkels bei Y (Stueck mit dem groessten u)
        g_ = _pg.intersection(_SLS([(-80, Y), (120, Y)]))
        st_ = [q_ for q_ in (getattr(g_, "geoms", None) or [g_]) if not q_.is_empty]
        if not st_: return None
        q_ = max(st_, key=lambda q_: q_.bounds[2])
        return q_.bounds[0]
    from shapely.geometry import Point as _Pt
    _holz = JP_HOLZ.get(r_) if JP else None
    def _ok(Y):                          # Knopfschenkel bei Y vorhanden, Unterkante deutlich ueber der Achse, nicht im Holzteil
        if not all(_u_unten(Y + e_) is not None for e_ in (-0.5, 0.0, 0.5)) or _u_unten(Y) <= up_ + 1.5: return False
        return _holz is None or not _Pt(_u_unten(Y), Y).buffer(2.0).intersects(_holz)
    # Knopfschenkel: Reihen 1/2 zeigen zu kleinerem Y, Reihen 3/4 zu groesserem Y -> beide Richtungen probieren
    Yh_ = next((Y_ for d_ in (8.0, 7.0, 6.0, 9.0, 10.0, 5.0, 4.0, 12.0) for Y_ in (Yp_ - d_, Yp_ + d_) if _ok(Y_)), None)
    if Yh_ is None:
        print(f"   Feder {G_}: kein Knopfschenkel fuer den Haken gefunden - uebersprungen"); continue
    ue_ = _u_unten(Yh_)
    _tu = (_u_unten(Yh_ + 0.5) - _u_unten(Yh_ - 0.5)); _l = (_tu ** 2 + 1.0) ** 0.5
    _du, _dY = -1.0 / _l, _tu / _l                                      # Kantennormale nach aussen (kleineres u)
    Qu, QY = ue_ + _du * (D_DRAHT / 2 + 0.05), Yh_ + _dY * (D_DRAHT / 2 + 0.05)
    d3_ = cq.Vector(0, _dY, -_du)
    xa_, xi_ = xh_ + T_ALU / 2 + D_DRAHT / 2 + 0.05, xh_ - T_ALU / 2 - D_DRAHT / 2 - 0.05
    q1_ = cq.Vector(xa_, QY, -Qu); q2_ = cq.Vector(xi_, QY, -Qu)
    q0_ = q1_.sub(d3_.multiply(3.0)); q3_ = q2_.sub(d3_.multiply(3.0))  # 3 mm an beiden Blechseiten
    # --- Windungen: Schenkel gehen tangential (ohne Knick) aus der Spirale, nie durch die Achse ---
    def _wink(y, z): return _m_.atan2(z - zp_, y - Yp_)
    hel0 = cq.Wire.makeHelix(steig_, steig_, rm_, cq.Vector(x0_, Yp_, zp_), cq.Vector(1, 0, 0))
    _sgn = 1.0 if ((_wink(hel0.positionAt(0.02).y, hel0.positionAt(0.02).z) - _wink(hel0.startPoint().y, hel0.startPoint().z)
                    + 3 * _m_.pi) % (2 * _m_.pi) - _m_.pi) > 0 else -1.0                  # Drehsinn der Wendel
    # Windungsanfang (Hebelseite): Tangente von q0 an den Windungskreis
    _Ey, _Ez = q0_.y - Yp_, q0_.z - zp_; _D = (_Ey ** 2 + _Ez ** 2) ** 0.5; _al = _m_.atan2(_Ez, _Ey)
    _kand = [_al + _m_.acos(min(1.0, rm_ / _D)), _al - _m_.acos(min(1.0, rm_ / _D))]
    def _rueck(th): return (_sgn * _m_.sin(th), -_sgn * _m_.cos(th))   # Richtung vom Windungsanfang weg (rueckwaerts)
    th2 = max(_kand, key=lambda th: (_Ey - rm_ * _m_.cos(th)) * _rueck(th)[0] + (_Ez - rm_ * _m_.sin(th)) * _rueck(th)[1])
    # Windungsende (Kunststoffseite): dort, wo die Wendel zur Kunststoffoberflaeche zeigt, Schenkel tangential weiter
    n3_ = cq.Vector(0, _n[1], -_n[0])
    th1 = _m_.atan2(-n3_.z, -n3_.y)                                     # Richtung Oberflaeche (wie bisher -n)
    frac = ((_sgn * (th1 - th2)) % (2 * _m_.pi)) / (2 * _m_.pi)
    wdg_ = (N_WDG - 1 + frac if frac > 0.5 else N_WDG + frac) - 0.5    # 19.11.: halbe Windung weniger -> gerader Schenkel liegt vorne (2,0 ... 3,0)
    hel_ = cq.Wire.makeHelix(steig_, steig_ * wdg_, rm_, cq.Vector(x0_, Yp_, zp_), cq.Vector(1, 0, 0))
    _drh = _m_.degrees(th2 - _wink(hel_.startPoint().y, hel_.startPoint().z))
    hel_ = hel_.rotate(cq.Vector(x0_, Yp_, zp_), cq.Vector(x0_ + 1, Yp_, zp_), _drh)
    prof_ = cq.Wire.makeCircle(D_DRAHT / 2, hel_.startPoint(), hel_.tangentAt(0))
    spirale = cq.Solid.sweep(prof_, [], hel_, True, True)
    pe_, pa_ = hel_.startPoint(), hel_.endPoint()                      # Anfang am Hebel, Ende auf der Kunststoffseite
    # Schenkel 1: tangential aus dem Windungsende, 9 mm, liegt auf dem Kunststoffteil
    t_end = hel_.tangentAt(1.0); t1_ = cq.Vector(0, t_end.y, t_end.z).normalized()
    for l1_ in (9.0, 8.0, 7.0, 6.0, 5.0, 4.0):                         # 19.11.: Schenkel nicht in die Griffstock-Abdeckung
        sch1 = [_draht(pa_, pa_.add(t1_.multiply(l1_)))]
        try:
            if sch1[0].intersect(abdeckung).Volume() < 1e-4: break
        except Exception: break
    # Schenkel 2: tangential aus dem Windungsanfang gerade bis zum U-Haken
    sch2 = [_draht(pe_, q0_), _kugel(q0_), _draht(q0_, q1_), _kugel(q1_), _draht(q1_, q2_), _kugel(q2_), _draht(q2_, q3_), _kugel(q3_)]
    fed_ = spirale
    for t_ in sch1 + sch2:
        try: fed_ = fed_.fuse(t_)
        except Exception: pass
    teile.append((f"{G_}/Feder", fed_, "Stahl"))
    RASTKERBE_DATEN.append((pa_, t1_, n3_, xh_))
    FEDERN.append(G_)
if FEDER: print("Federn:", len(FEDERN), FEDERN[:3])
# 19.11. (Johann): Rastmarken im Kunststoffteil in Flucht mit jedem Hebel, nur an der hinteren Kante der Flaeche, auf der
# Federschenkel 1 liegt: ca. 2 mm lang, 0,5 mm tief, 0,8 breit (Draht 0,6). Schenkel werden in der Praxis bis dorthin
# verlaengert (im Modell nicht).
RASTKERBEN = []
RAST_VERSCHOBEN = []
if FEDERN:
    from shapely.geometry import LineString as _LSr
    _ik = [i_ for i_, t_ in enumerate(teile) if t_[0] == "Kunststoffteil_Lagerblock"][0]
    _kk = teile[_ik][1]
    if os.environ.get("DEBUG_RAST"):
        import pickle as _pk
        cq.exporters.export(_kk, os.environ["DEBUG_RAST"] + ".brep") if False else _kk.exportBrep(os.environ["DEBUG_RAST"] + ".brep")
        _pk.dump(([((p_.x, p_.y, p_.z), (t_.x, t_.y, t_.z), (n_.x, n_.y, n_.z), x_) for p_, t_, n_, x_ in RASTKERBE_DATEN], KST_E, D_DRAHT),
                 open(os.environ["DEBUG_RAST"] + ".pkl", "wb"))
        print("DEBUG_RAST gespeichert"); sys.exit(0)
    for (pa1_, t1k_, n3k_, xk_) in RASTKERBE_DATEN:
        ein_ = cq.Vector(n3k_.x, n3k_.y, n3k_.z).normalized()               # in den Kunststoff hinein (n zeigt ins Teil)
        S_ = pa1_.add(ein_.multiply(D_DRAHT / 2))                           # Punkt auf der Flaeche unter dem Schenkel
        r_dir = t1k_ if t1k_.y < 0 else t1k_.multiply(-1)                   # nach hinten (zum Deckel, kleineres Y)
        a2 = (-S_.z, S_.y); d2 = (-r_dir.z, r_dir.y)                        # (u, Y)
        ray_ = _LSr([(a2[0] - d2[0] * 0.5, a2[1] - d2[1] * 0.5), (a2[0] + d2[0] * 80, a2[1] + d2[1] * 80)])
        g_ = KST_E.buffer(0.02).intersection(ray_)
        st_ = [q_ for q_ in (getattr(g_, "geoms", None) or [g_]) if not q_.is_empty and q_.geom_type == "LineString"]
        if not st_: continue
        q_ = min(st_, key=lambda q_: q_.distance(_SPG([a2, (a2[0] + 0.01, a2[1]), (a2[0], a2[1] + 0.01)]).centroid))
        e_ = max(q_.coords, key=lambda c_: (c_[0] - a2[0]) * d2[0] + (c_[1] - a2[1]) * d2[1])
        E_ = cq.Vector(xk_, e_[1], -e_[0])                                   # hintere Kante, in Flucht mit dem Hebel
        quer_ = cq.Vector(1, 0, 0)
        def _kerbe(dx_):
            def _pkt(a_, b_, c_): return E_.add(r_dir.multiply(a_)).add(quer_.multiply(b_ + dx_)).add(ein_.multiply(c_))
            _ecken = [_pkt(a_, b_, -1.0) for a_, b_ in ((-2.0, -0.4), (1.0, -0.4), (1.0, 0.4), (-2.0, 0.4))]
            _w = cq.Wire.makePolygon(_ecken + [_ecken[0]])
            return cq.Solid.extrudeLinear(cq.Face.makeFromWires(_w), ein_.multiply(1.5))
        # 20.11. (Johann): mindestens ein Teil der 2 mm muss an der hinteren oberen Kante Material treffen. Liegt in Flucht
        # mit dem Hebel der Hebelschlitz (Reihe 4: Schlitz laeuft durch die Kante), Kerbe seitlich zum naechsten Material.
        _k_ = None
        for dx_ in [0.0] + [s_ * 0.25 * j_ for j_ in range(1, 33) for s_ in (1, -1)]:
            _kt = _kerbe(dx_)
            try: _vt = _kt.intersect(_kk).Volume()
            except Exception: _vt = 0.0
            if _vt >= 0.8:                                                  # Sollwert voll im Material ca. 1,6
                _k_ = _kt
                if dx_ != 0.0: RAST_VERSCHOBEN.append((round(xk_, 1), dx_))
                break
        if _k_ is None: _k_ = _kerbe(0.0); RAST_VERSCHOBEN.append((round(xk_, 1), None))
        RASTKERBEN.append(_k_)
    _kk_alt = _kk
    _kk = schneide(_kk, RASTKERBEN)
    teile[_ik] = (teile[_ik][0], _kk, teile[_ik][2])
    _vk = teile[_ik][1].Volume()
    print(f"   Rastmarken: {sum(1 for k_ in RASTKERBEN if k_.intersect(_kk_alt).Volume() > 0.05)} schneiden ins Kunststoffteil")
    if RAST_VERSCHOBEN: print("   Rastmarken seitlich verschoben (Hebelschlitz an der Kante), x: dx =", RAST_VERSCHOBEN)
    print(f"Kunststoffteil: {len(RASTKERBEN)} Rastmarken an der hinteren Kante (2 lang, 0,5 tief, 0,8 breit, in Flucht mit den Hebeln)")
teile.append(("Dreieckleiste", dre, "Griffstock"))
verst = schneide(verst, loch_b + [cq.Solid.makeCylinder(1.5, 7.0, cq.Vector(x_, y_, T_BODEN - 1.0), cq.Vector(0, 0, 1)) for x_, y_, n_, w_ in GS_LOECHER])
teile.append(("Verstaerkungsleiste_innen", verst, "Holz"))
if _WINKEL: exec(open("winkel_teile.py").read())
if TRAD: exec(open("kugel.py" if KUGEL else "trad.py").read())
import json as _js
open("drehpunkte.json", "w").write(_js.dumps(DREHPUNKTE))
if not TRAD and os.environ.get("KST_STEP"):          # 19.11.: Kunststoffteil in Fraeslage (Lagerseite oben, Z = 0): Fraesteil + Fertigteil
    _kf = [t_[1] for t_ in teile if t_[0] == "Kunststoffteil_Lagerblock"][0]
    _sv = [kst_s(p_) for p_ in KST_E.exterior.coords]
    _dv = [(p_[0] - _a[0]) * _n[0] + (p_[1] - _a[1]) * _n[1] for p_ in KST_E.exterior.coords]
    _smin, _smax, _dmax = min(_sv), max(_sv), max(_dv)
    # Rohblock in der Lagerseiten-Lage: was ausserhalb des Querschnitts liegt (Abschraegungen), bleibt stehen
    _rb = prisma_uy(kst_rechteck(_smin, _smax, 0.0, _dmax), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0)
    _fr = _kf.fuse(_rb.cut(prisma_uy(_pts(KST_E), (KST_X0 + KST_X1) / 2, KST_X1 - KST_X0 + 0.02))).clean()
    _wk = math.degrees(math.atan2(-_n[1], _n[0]))                      # Aussennormale der Lagerseite (0, -n1, n0) -> +Z
    from OCP.gp import gp_Trsf as _gT, gp_Ax1 as _gA1, gp_Pnt as _gP, gp_Dir as _gD, gp_Vec as _gV
    from OCP.BRepBuilderAPI import BRepBuilderAPI_Transform as _BT
    def _kst_lage(so_, dv_=None):
        tr_ = _gT(); tr_.SetRotation(_gA1(_gP(0, 0, 0), _gD(1, 0, 0)), math.radians(_wk))
        if dv_ is not None:
            t2_ = _gT(); t2_.SetTranslation(_gV(dv_.x, dv_.y, dv_.z)); tr_ = t2_.Multiplied(tr_)
        return cq.Shape.cast(_BT(so_.wrapped, tr_, True).Shape())
    _zl = _kst_lage(cq.Vertex.makeVertex(0.0, _a[1], -_a[0])).Center().z
    _bbr = _kst_lage(_fr).BoundingBox()
    _verschieb = cq.Vector(-KST_X0, -_bbr.ymin, -_zl)
    for _nm, _so in (("", _fr), ("_fertig", _kf)):
        cq.exporters.export(_kst_lage(_so, _verschieb), os.environ["KST_STEP"].replace(".step", _nm + ".step"))
    print(f"Kunststoffteil in Fraeslage: {KST_X1 - KST_X0:.1f} x {_smax - _smin:.1f} x {_dmax:.1f} mm -> {os.environ['KST_STEP']}")
if not TRAD and os.environ.get("TEILE_DUMP"):          # 19.11.: Daten fuer teile_dxf.py (Alu-Abwicklungen, Holzteile, Kunststoffteil)
    import pickle as _pkd
    _hd = []
    for d_ in HEBEL_DATEN:
        _pr = []
        for f_ in d_["hb"].Faces():
            try:
                if abs(f_.normalAt().x) < 0.2: continue
                _pr.append(_PG2([(v_.y, -v_.z) for v_ in (f_.outerWire().positionAt(k_ / 300) for k_ in range(301))]).buffer(0))
            except Exception: pass
        _hd.append({k_: v_ for k_, v_ in d_.items() if k_ != "hb"} | {"proj": _juu(_pr).buffer(0.001).buffer(-0.001)})
    _holz = {r_: (JP_HOLZ[r_] if JP else None) for r_ in (1, 2, 3, 4)}
    _pkd.dump({"hebel": _hd, "holz": _holz, "kst_name": "Kunststoffteil_Lagerblock", "T_ALU": T_ALU, "T_HOLZ": T_HOLZ,
               "D_BUCHSE": D_BUCHSE, "KV_U": (KV_U0, KV_U1), "instrument": (INSTR["name"] if globals().get("INSTR") else "Markus")},
              open(os.environ["TEILE_DUMP"], "wb"))
    print("Teile-Daten ->", os.environ["TEILE_DUMP"])
