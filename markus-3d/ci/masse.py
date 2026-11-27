# Gehaeusemasse eines fertigen Modells aendern (Johann 27.11.2026): Breite (Tiefe, Y) und Laenge (X) als Aussenmass.
#   python3 masse.py ein.step aus.step --teil bass|diskant [--breite 210] [--laenge 400] [--bericht bericht.json]
# Regeln:
#  - Verbreitern: Bass - der ganze Einbau bekommt vorne (Front, hohes Y) mehr Raum; Diskant - vorne (Seite mit der
#    Rahmeneinlage, hohes Y) wird der Freiraum und der Boden breiter, der Griffstock bleibt hinten an der Aussenkante.
#  - Schmaeler (bis 200): nur in freiem Raum, an Laengsleisten (bleiben mindestens 2 mm) und Gehaeusewaenden (bleiben
#    mindestens 4 mm); Stimmplatten, Stimmstoecke, Kammern, Klappen und die Mechanik werden nie beschnitten.
#  - Verlaengern: Leerraum zwischen Aussenwand und Einbau an beiden Enden (je die Haelfte); Deckel, Boden, Balg-/
#    Rahmenleisten, Fuesse und Waende bleiben an den Aussenmassen ausgerichtet.
# Umsetzung wie streck.py (Markus-Bass): jedes Teil, das die Schnittebene kreuzt, muss dort prismatisch sein (gleicher
# Querschnitt) und wird aufgeschnitten und mit seinem Querschnitt gefuellt bzw. um den Streifen gekuerzt; kleine Teile
# wandern als Ganzes; Randbretter (Front, Rueckwand, Seitenwaende, Rahmen) mit Ornamenten, Hornbohrungen oder Zinken an
# der Ebene werden an ihrer eigenen naechsten glatten Stelle geschnitten, was in ihnen sitzt (Hoerner, Einlagen) wandert
# mit. Danach werden alle neu entstandenen Beruehrungen gesucht; gibt es welche, bricht das Programm mit Meldung ab.
import sys, os, re, math, json, argparse, time
import cadquery as cq
from OCP.gp import gp_GTrsf, gp_Mat, gp_XYZ
from OCP.BRepBuilderAPI import BRepBuilderAPI_GTransform

VERBOTEN = re.compile(r"Stimmplatte|Stimmstock|Kammer|Klappe|Trennfeder|Sole|Hebel_R|Achse|Lagerbock|Begleiterstimmstock|Gelenk|Schenkel|Feder")
GANZ = re.compile(r"(^|/)(Horn|Knopf_|Gaze_|Rahmeneinlage|Filzeinlage|Fuss|Haken|Schraube|Raendel)")
# Deckel/Boden mit Rosetten oder Durchbruechen an den Enden: wie Randbretter an ihrer eigenen glatten Stelle (nach innen)
# schneiden; was am Deckel haengt (Gaze, Eckbloecke, Fuesse) wandert dann mit dem Rand - Schachtel an den Aussenmassen
STANGE = re.compile(r"(^|/)(Schubstange|Luftklappe_Schubstange|Lufttaste_Stange)")   # Draehte duerfen auch schraeg gestreckt werden
DECKEL = re.compile(r"^[^/]*/(Bassdeckel|Bassdeckel_Leiste_hinten|Bassdeckel_Leiste_vorne|Balgrahmen)$")
ECKE = re.compile(r"(^|/)(Eckblock|Eckwinkel|Fuss)")       # Ecken und Fuesse bleiben an den Aussenmassen
_WAND = re.compile(r"^[^/]*/(Front|Rueckteil|Seitenteil_links|Seitenteil_rechts|Rahmen_hinten|Rahmen_vorne|Rahmen_links|Rahmen_rechts)$")
class _W:                                                   # Gehaeusewand = direktes Teil des Modells mit Wandnamen
    @staticmethod
    def search(n): return _WAND.match(n)
WAND = _W
KLEIN = 12.0                 # Teile kuerzer als das (in Aenderungsrichtung) wandern als Ganzes
MIN_WAND, MIN_LEISTE = 4.0, 2.0
TOL_V, TOL_B = 2e-3, 0.03

class MassFehler(SystemExit):
    pass

def blaetter(a, loc=None, pfad=""):
    loc = (loc * a.loc) if loc is not None else a.loc
    p = (pfad + "/" + a.name) if pfad else a.name
    out = []
    if a.obj is not None:
        objs = a.obj.vals() if isinstance(a.obj, cq.Workplane) else [a.obj]
        for o in objs:
            if isinstance(o, cq.Shape): out.append({"name": p, "s": o.moved(loc), "knoten": a, "loc": loc})
    for c in a.children: out += blaetter(c, loc, p)
    return out

def bbox(s):
    b = s.BoundingBox(); return (b.xmin, b.ymin, b.zmin, b.xmax, b.ymax, b.zmax)

def kasten(k, a, b, bb, rand=5.0):
    lo = [bb[0] - rand, bb[1] - rand, bb[2] - rand]; hi = [bb[3] + rand, bb[4] + rand, bb[5] + rand]
    lo[k], hi[k] = a, b
    return cq.Solid.makeBox(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2], cq.Vector(*lo))

def scheibe(s, k, u, bb, d=0.1):
    try:
        c = s.intersect(kasten(k, u - d / 2, u + d / 2, bb))
        v = c.Volume()
        if v < 1e-9: return (0.0, None)
        cb = bbox(c); return (v, tuple(round(x, 2) for i, x in enumerate(cb) if i % 3 != k))
    except Exception:
        return (None, None)

_CACHE = {}
def prismatisch(t, k, a, b):
    key = (id(t["s"]), k, round(a, 2), round(b, 2))
    if key not in _CACHE: _CACHE[key] = _prismatisch(t, k, a, b)
    return _CACHE[key]

def _prismatisch(t, k, a, b):
    """gleicher Querschnitt bei a, Mitte, b (und dazwischen alle 2 mm), nur soweit das Teil reicht"""
    a, b = max(a, t["bb"][k] + 0.06), min(b, t["bb"][k + 3] - 0.06)
    if b <= a: return True
    n = max(2, int(math.ceil((b - a) / 2.0)))
    ref = None
    for i in range(n + 1):
        u = a + (b - a) * i / n
        v, q = scheibe(t["s"], k, u, t["bb"])
        if v is None: return False
        if ref is None: ref = (v, q); continue
        if abs(v - ref[0]) > TOL_V * max(ref[0], 1e-6) + 1e-6: return False
        if (q is None) != (ref[1] is None): return False
        if q is not None and any(abs(x - y) > TOL_B for x, y in zip(q, ref[1])): return False
    return True

def ausdehnung(t, k): return t["bb"][k + 3] - t["bb"][k]

def pruefe_ebene(teile, k, p, r, schrumpf):
    """Ebene p (r = 0, Einfuegen) bzw. Streifen [p, p+r] (Entfernen). Gibt (ok, {index: art}, grund) zurueck."""
    arten = {}
    a_, b_ = (p - 0.3, p + 0.3) if r == 0 else (p - 0.2, p + r + 0.2)
    kreuz = [i for i, t in enumerate(teile) if not (t["bb"][k + 3] <= a_ + 0.2 or t["bb"][k] >= b_ - 0.2)]
    for i in kreuz:                                          # zuerst die billigen Pruefungen (nur Namen und Lage)
        t = teiles = teile[i]; n = t["name"]; lo, hi = t["bb"][k], t["bb"][k + 3]; e = hi - lo
        if r > 0 and VERBOTEN.search(n): return False, None, f"{n} wuerde beschnitten"
        if r == 0 and not (GANZ.search(n) or ECKE.search(n) or (e < KLEIN and not WAND.search(n))):
            if VERBOTEN.search(n): return False, None, f"{n} wuerde geschnitten"
            if e < 8.0: return False, None, f"{n} (duenn, {e:.1f} mm) wuerde dicker"
    kreuz.sort(key=lambda i: len(teile[i]["s"].Faces()))      # einfache Teile zuerst pruefen
    for i in kreuz:
        t = teile[i]
        lo, hi = t["bb"][k], t["bb"][k + 3]
        n = t["name"]; e = hi - lo; mitte = (lo + hi) / 2
        if r > 0:                                            # Streifen entfernen
            if VERBOTEN.search(n): return False, None, f"{n} wuerde beschnitten"
            if GANZ.search(n):
                if lo > p - 0.05 and hi < p + r + 0.05: return False, None, f"{n} liegt im Streifen"
                # durch die Wand gesteckte Teile (Hoerner, Filzeinlagen): wo sie gleichmaessig sind, mit kuerzen
                if hi > p + r + 0.3 and prismatisch(t, k, a_, b_): arten[i] = "schnitt"; continue
                if "Horn" in n and hi > p - 0.05 and lo < p + r + 0.05: arten[i] = "horn"; continue   # Johann 27.11.: Horn versetzen/kuerzen
                arten[i] = "ganz"; continue
            if STANGE.search(n) or prismatisch(t, k, a_, b_):
                mind = MIN_WAND if WAND.search(n) else (MIN_LEISTE if e < 30 else 10.0)
                if e - r < mind - 1e-6: return False, None, f"{n} wuerde {e - r:.1f} mm (mindestens {mind:g})"
                arten[i] = "schnitt"; continue
            if WAND.search(n) and e > 60: arten[i] = "rand"; continue
            return False, None, f"{n} ist bei {p:.1f} nicht gleichmaessig"
        if GANZ.search(n) or ECKE.search(n) or (e < KLEIN and not WAND.search(n)):
            arten[i] = "ganz"; continue
        if VERBOTEN.search(n): return False, None, f"{n} wuerde geschnitten"
        if e < 8.0: return False, None, f"{n} (duenn, {e:.1f} mm) wuerde dicker"
        if STANGE.search(n) or prismatisch(t, k, a_, b_): arten[i] = "schnitt"; continue
        if (WAND.search(n) or DECKEL.match(n)) and e > 60: arten[i] = "rand"; continue
        return False, None, f"{n} ist bei {p:.1f} nicht gleichmaessig"
    return True, arten, ""

def eigene_ebene(t, k, p, r, seite):
    """naechste glatte Stelle eines Randbretts bei p (seite +1: Richtung hoeheres u)"""
    lo, hi = t["bb"][k] + 6.0, t["bb"][k + 3] - 6.0
    for j in range(0, 240):
        for sg in ((seite, -seite) if j else (1,)):
            q = p + sg * j * 0.5
            if not lo <= q <= hi - r: continue
            if prismatisch(t, k, q - 0.3, q + r + 0.3): return q
    return None

def gstreck(s, k, u0, f):
    """affin in Richtung k um u0 mit Faktor f strecken (Scheibe -> laengere/kuerzere Scheibe, Enden passen)"""
    m = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]; m[k][k] = f
    g = gp_GTrsf(); g.SetVectorialPart(gp_Mat(*[v for z in m for v in z]))
    tr = [0.0, 0.0, 0.0]; tr[k] = u0 * (1.0 - f); g.SetTranslationPart(gp_XYZ(*tr))
    if s.Volume() < 1e-9: return s
    try: return cq.Shape.cast(BRepBuilderAPI_GTransform(s.wrapped, g, True).Shape())
    except Exception:                                        # einzelne Koerper getrennt versuchen
        teile_ = []
        for so in s.Solids():
            try: teile_.append(cq.Shape.cast(BRepBuilderAPI_GTransform(so.wrapped, g, True).Shape()))
            except Exception: raise MassFehler("Strecken eines Querschnitts nicht moeglich (OCC GTransform)")
        return cq.Compound.makeCompound(teile_)

def verschiebe(s, k, d):
    v = [0.0, 0.0, 0.0]; v[k] = d; return s.translate(cq.Vector(*v))

def schneide(t, k, q, r, d):
    """r = 0: bei q um d strecken; r > 0: Streifen [q, q+r] entfernen"""
    s, bb = t["s"], t["bb"]
    def stueck(a, b):
        if b - a < 1e-6: return None
        c = s.intersect(kasten(k, a, b, bb)); return c if c.Volume() > 1e-9 else None
    if r == 0:                                               # Scheibe [u0, q] auf (q - u0) + d strecken (auch schraege Stangen bleiben stetig)
        u0 = max(q - 0.3, bb[k] + 1e-3)                     # nur der gepruefte Bereich wird gestreckt
        li = stueck(bb[k] - 1, u0)
        sl = stueck(u0, q)
        mi = gstreck(sl, k, u0, (q - u0 + d) / (q - u0)) if sl is not None else None
        re_ = stueck(q, bb[k + 3] + 1); re_ = verschiebe(re_, k, d) if re_ is not None else None
        teile_ = [x for x in (li, mi, re_) if x is not None]
    elif not STANGE.search(t["name"]):                       # gleichmaessiges Teil: Streifen herausnehmen, Rest zusammenschieben
        li = stueck(bb[k] - 1, q)
        re_ = stueck(q + r, bb[k + 3] + 1); re_ = verschiebe(re_, k, -r) if re_ is not None else None
        teile_ = [x for x in (li, re_) if x is not None]
    else:                                                    # Draht: Scheibe [q, q+r+1] auf 1 mm stauchen (bleibt stetig)
        u1 = min(q + r + 1.0, bb[k + 3] - 1e-3)
        li = stueck(bb[k] - 1, q)
        sl = stueck(q, u1)
        mi = gstreck(sl, k, q, max(u1 - q - r, 1e-3) / (u1 - q)) if sl is not None and u1 - q > r else None
        re_ = stueck(u1, bb[k + 3] + 1); re_ = verschiebe(re_, k, -r) if re_ is not None else None
        teile_ = [x for x in (li, mi, re_) if x is not None]
    teile_ = [x for x in teile_ if x.Volume() > 1e-9]
    neu = teile_[0]
    for x in teile_[1:]:
        neu = neu.fuse(x)
    try: neu = neu.clean()
    except Exception: pass
    sol = neu.Solids()
    if len(sol) == 1: neu = sol[0]
    return neu

def teiles_name(teile, i): return teile[i]["name"]

def aendere(teile, k, p, r, d, arten, bericht, was):
    """wendet eine Ebene an; d > 0 Einfuegen bei p, r > 0 Streifen [p, p+r] entfernen.
    Teile mit hoeherer Koordinate als die Ebene wandern (um d bzw. -r), die anderen bleiben."""
    _CACHE.clear()
    eigene = {}
    hoch = p + (r if r > 0 else 0.0)
    beruehrt = lambda b, o, tol=1.0: all(b[m] < o[m + 3] + tol and b[m + 3] > o[m] - tol for m in range(3))
    # 1) Randgruppen: was an einer Querwand haengt (Leisten, Kloetze, Eckbloecke, Fuesse samt Schrauben), bleibt bei ihr
    quer = lambda t: WAND.search(t["name"]) and ausdehnung(t, k) < 8.0
    seiten = {}
    for seite, wand in (("unten", [i for i, t in enumerate(teile) if quer(t) and t["bb"][k + 3] <= p + 1e-6]),
                        ("oben", [i for i, t in enumerate(teile) if quer(t) and t["bb"][k] >= hoch - 1e-6])):
        if not wand: continue
        w0 = min(teile[i]["bb"][k] for i in wand); w1 = max(teile[i]["bb"][k + 3] for i in wand)
        nahe = lambda t: (t["bb"][k] < w1 + 45.0) if seite == "unten" else (t["bb"][k + 3] > w0 - 45.0)
        kann = lambda i, t: (i not in arten or arten[i] == "ganz") and not VERBOTEN.search(t["name"]) and nahe(t) and "Gaze_" not in t["name"] and \
                            (GANZ.search(t["name"]) or ECKE.search(t["name"]) or ausdehnung(t, k) < KLEIN)
        grp = set(wand) | {i for i, t in enumerate(teile) if ECKE.search(t["name"]) and nahe(t) and not VERBOTEN.search(t["name"])}
        neu_ = list(grp)
        while neu_:
            j = neu_.pop()
            for i, t in enumerate(teile):
                if i in grp or not kann(i, t): continue
                if beruehrt(t["bb"], teile[j]["bb"], 0.3): grp.add(i); neu_.append(i)
        g_ = grp - set(wand)
        if r > 0:                                            # Streifen entfernen: nur was ganz ausserhalb des Streifens liegt, wandert mit
            g_ = {i for i in g_ if (teile[i]["bb"][k] >= hoch - 0.01 if seite == "oben" else teile[i]["bb"][k + 3] <= p + 0.01)}
        seiten[seite] = g_
    mit = {i: s_ for s_, g in seiten.items() for i in g}
    deck = [i for i, a in arten.items() if a in ("schnitt", "rand") and DECKEL.match(teiles_name(teile, i))] if r == 0 else []
    blatt = lambda t: "leder" in t["name"].lower() or min(t["bb"][m + 3] - t["bb"][m] for m in range(3)) < 1.5
    if deck:                                                 # Leder/Blaetter am Deckel werden mit ihm geschnitten
        deck += [i for i, a in arten.items() if a == "schnitt" and i not in deck and blatt(teile[i]) and not GANZ.search(teile[i]["name"])
                 and any(beruehrt(teile[i]["bb"], teile[j]["bb"], 0.3) for j in deck)]
    # Teile des Einbaus, die durch den Deckel gehen (Lufttaste, Schieber): zwischen Ebene und Deckelschnitt darf keines liegen
    durch = [i for i, t in enumerate(teile) if i not in mit and i not in deck and not GANZ.search(t["name"]) and not WAND.search(t["name"])
             and arten.get(i) not in ("schnitt", "rand") and not blatt(t)   # Leder/Gaze-Blaetter nicht
             and any(beruehrt(t["bb"], teile[j]["bb"], -0.2) for j in deck)]
    def deckel_ebene(innen_mit, haengt):
        innen = 1 if (any(mit[i] == "unten" for i in innen_mit) or (not innen_mit and p < (teile[deck[0]]["bb"][k] + teile[deck[0]]["bb"][k + 3]) / 2)) else -1
        sperre = [(teile[j]["bb"][k] - 0.5, teile[j]["bb"][k + 3] + 0.5) for j in haengt + innen_mit]
        for jj in range(0, 500):
            u = p + innen * jj * 0.5
            if any(a_ < u < b_ for a_, b_ in sperre): continue
            if innen_mit and ((innen > 0 and u < max(teile[i]["bb"][k + 3] for i in innen_mit) + 1.0) or
                              (innen < 0 and u > min(teile[i]["bb"][k] for i in innen_mit) - 1.0)): continue
            lo_, hi_ = sorted((p, u))
            if any(teile[j]["bb"][k + 3] > lo_ + 0.1 and teile[j]["bb"][k] < hi_ - 0.1 for j in durch): return None
            if all(teile[i]["bb"][k] + 1.0 < u < teile[i]["bb"][k + 3] - 1.0 and prismatisch(teile[i], k, u - 0.3, u + r + 0.3) for i in deck):
                return u
        return None
    jj_ = 1 - k                                              # Querrichtung (Baender des gestuften Deckelschnitts)
    def baender(i, innen_mit, haengt):
        """gestufter Schnitt eines Deckelteils (Johann 27.11.): in jedem Querbereich eine eigene Schnittstelle -
        hinter Fuessen/Eckbloecken (die wandern mit der Wand), vor Lufttaste/Schieber (die bleiben beim Einbau)"""
        t = teile[i]; tb = t["bb"]
        innen = 1 if p < (tb[k] + tb[k + 3]) / 2 else -1
        M = [m for m in innen_mit if beruehrt(teile[m]["bb"], tb, 1.0)]
        U = [u for u in durch if beruehrt(teile[u]["bb"], tb, -0.2)]
        H = [h for h in haengt if beruehrt(teile[h]["bb"], tb, 1.0)]
        kanten = sorted({tb[jj_], tb[jj_ + 3]} | {teiles_bb(x)[jj_] for x in M + U + H} | {teiles_bb(x)[jj_ + 3] for x in M + U + H})
        kanten = [x for x in kanten if tb[jj_] - 1e-6 <= x <= tb[jj_ + 3] + 1e-6]
        def stueck_t(a, b):
            lo = list(tb[:3]); hi = list(tb[3:]); lo[jj_], hi[jj_] = a, b
            box_ = cq.Solid.makeBox(hi[0] - lo[0], hi[1] - lo[1], hi[2] - lo[2] + 10, cq.Vector(lo[0], lo[1], lo[2] - 5))
            c = t["s"].intersect(box_); return {"s": c, "bb": bbox(c), "name": t["name"]} if c.Volume() > 1e-6 else None
        def aktiv(liste, a, b): return [x for x in liste if teiles_bb(x)[jj_] < b - 1e-6 and teiles_bb(x)[jj_ + 3] > a + 1e-6]
        def passt(u, a, b):
            for m in aktiv(M, a, b):
                if (innen < 0 and u > teiles_bb(m)[k] - 1.0) or (innen > 0 and u < teiles_bb(m)[k + 3] + 1.0): return False
            for x in aktiv(U, a, b):
                if (innen < 0 and u < teiles_bb(x)[k + 3] + 0.5) or (innen > 0 and u > teiles_bb(x)[k] - 0.5): return False
            ha = aktiv(H, a, b)                              # Gaze/Rosetten eines Bereichs nicht trennen: alle auf einer Seite
            if ha and not (u < min(teiles_bb(h)[k] for h in ha) - 0.5 or u > max(teiles_bb(h)[k + 3] for h in ha) + 0.5): return False
            return True
        roh = []
        for a, b in zip(kanten[:-1], kanten[1:]):
            if b - a < 1e-3: continue
            st = stueck_t(a, b)
            if st is None: continue
            q = None
            for n_ in range(0, 500):
                u = p + innen * n_ * 0.5
                if not (tb[k] + 1.0 < u < tb[k + 3] - 1.0): continue
                if passt(u, a, b) and prismatisch(st, k, u - 0.3, u + 0.3): q = u; break
            if q is None:
                nm = lambda x: teile[x]["name"].split("/")[-1]
                print(f"  Stufenschnitt {t['name'].split('/')[-1]} {('XYZ')[jj_]} {a:.1f}-{b:.1f}: keine Stelle; Wand:",
                      [nm(x) for x in aktiv(M, a, b)], "Einbau:", [nm(x) for x in aktiv(U, a, b)], "Gaze:", [nm(x) for x in aktiv(H, a, b)], flush=True)
                return None
            roh.append([a, b, q])
        ba = []                                              # gleiche Schnittstellen zusammenfassen
        for e in roh:
            if ba and abs(ba[-1][2] - e[2]) < 1e-6 and abs(ba[-1][1] - e[0]) < 1e-3: ba[-1][1] = e[1]
            else: ba.append(e)
        for n_ in range(len(ba) - 1):                       # Stufe an eine glatte Stelle legen (keine Bohrung zerreissen)
            a0, b0, q0 = ba[n_]; a1, b1, q1 = ba[n_ + 1]
            ylo, yhi = min(q0, q1) - 0.5, max(q0, q1) + abs(d) + 0.5
            def glatt(x):
                lo = list(tb[:3]); hi = list(tb[3:]); lo[k], hi[k] = ylo, yhi
                box_ = cq.Solid.makeBox(hi[0] - lo[0] + (0 if jj_ else 0), hi[1] - lo[1], hi[2] - lo[2] + 10, cq.Vector(lo[0], lo[1], lo[2] - 5))
                c = t["s"].intersect(box_)
                return c.Volume() > 1e-6 and prismatisch({"s": c, "bb": bbox(c)}, jj_, x - 0.3, x + 0.3)
            ok = False
            for sch in [0.5 * z for z in range(0, 80)]:
                for x in ((b0 - sch, b0 + sch) if sch else (b0,)):
                    if not (a0 + 0.5 < x < b1 - 0.5): continue
                    if x < b0 and not passt(q1, x, b0): continue
                    if x > b0 and not passt(q0, b0, x): continue
                    if glatt(x): ba[n_][1] = x; ba[n_ + 1][0] = x; ok = True; break
                if ok: break
            if not ok:
                print(f"  Stufenschnitt {t['name'].split('/')[-1]}: Stufe bei {('XYZ')[jj_]} {b0:.1f} findet keine glatte Stelle", flush=True)
                return None
        return ba
    def teiles_bb(x): return teile[x]["bb"]
    band = {}
    # 2) Deckel: gemeinsam innerhalb der Randgruppen schneiden (Fuss-/Eckblockbohrungen, Rosetten bleiben am Rand)
    innen_mit = [] if r > 0 else [i for i in mit if ((mit[i] == "unten" and teile[i]["bb"][k + 3] > p) or (mit[i] == "oben" and teile[i]["bb"][k] < hoch))
                 and (ECKE.search(teile[i]["name"]) or GANZ.search(teile[i]["name"]))]   # nur was im Deckel verschraubt ist
    if innen_mit:                                            # Leisten, durch die Eckwinkel-/Fussschrauben gehen, wie den Deckel behandeln
        for i, a_ in list(arten.items()):
            if a_ == "schnitt" and i not in deck and not WAND.search(teile[i]["name"]) and not VERBOTEN.search(teile[i]["name"]) \
               and not STANGE.search(teile[i]["name"]) and not ECKE.search(teile[i]["name"]) and not GANZ.search(teile[i]["name"]) and any(beruehrt(teile[i]["bb"], teile[m]["bb"], -0.2) for m in innen_mit):
                deck.append(i)
    if deck and (any(arten[i] == "rand" for i in deck) or innen_mit):
        haengt = [j for j, t in enumerate(teile) if GANZ.search(t["name"]) and j not in mit and any(beruehrt(t["bb"], teile[i]["bb"]) for i in deck)]
        q = deckel_ebene(innen_mit, haengt)
        if q is None and innen_mit:                          # gestufter Schnitt: Ecken am Rand, Lufttaste/Schieber beim Einbau
            bb_ = {i: baender(i, innen_mit, haengt) for i in deck}
            if all(v is not None for v in bb_.values()):
                for i, v in bb_.items():
                    band[i] = v; arten[i] = "band"
                    bericht.append(f"{teile[i]['name']}: gestuft geschnitten " + ", ".join(f"{('XYZ')[jj_]} {a:.0f}-{b:.0f} bei {q_:.1f}" for a, b, q_ in v))
                q = "band"
        if q is None and innen_mit:                          # geht nicht: Ecken/Fuesse bleiben beim Einbau (Deckelbohrungen passen)
            for i in [i for i in innen_mit if ECKE.search(teile[i]["name"]) or GANZ.search(teile[i]["name"])]: mit.pop(i, None)
            innen_mit = [i for i in innen_mit if i in mit]
            bericht.append(f"{was}: Eckbloecke/Fuesse bleiben beim Einbau (Deckel hat dort Durchbrueche fuer den Einbau)")
            haengt = [j for j, t in enumerate(teile) if GANZ.search(t["name"]) and j not in mit and any(beruehrt(t["bb"], teile[i]["bb"]) for i in deck)]
            q = deckel_ebene(innen_mit, haengt)
        if q is None:
            nm = lambda i: teile[i]["name"].split("/")[-1] + f"[{teile[i]['bb'][k]:.1f}..{teile[i]['bb'][k + 3]:.1f}]"
            print("  Deckel:", [nm(i) for i in deck], "\n  innen bei der Wand:", [nm(i) for i in innen_mit],
                  "\n  durch den Deckel:", [nm(i) for i in durch], "\n  haengt am Deckel:", [nm(i) for i in haengt], flush=True)
            raise MassFehler(f"{was}: Deckel/Balgrahmen haben keine passende glatte Stelle fuer den Schnitt")
        if q != "band":
            for i in deck: eigene[i] = q; arten[i] = "rand"
    for i, a in arten.items():
        if a == "rand" and i not in eigene:
            t = teile[i]; mitte = (t["bb"][k] + t["bb"][k + 3]) / 2
            q = eigene_ebene(t, k, p, r, 1 if p < mitte else -1)
            if q is None: raise MassFehler(f"{was}: {t['name']} hat keine glatte Stelle fuer den Schnitt")
            eigene[i] = q
    folgt = dict(eigene)                                     # Gaze, Hoerner, Einlagen folgen ihrem eigens geschnittenen Teil
    neu_ = list(eigene)
    while neu_:
        j = neu_.pop(); q = folgt[j]; lo_, hi_ = sorted((p, q))
        for i, t in enumerate(teile):
            if i in folgt or i in mit or arten.get(i) not in ("ganz", None) or not GANZ.search(t["name"]): continue
            mi = (t["bb"][k] + t["bb"][k + 3]) / 2
            if lo_ - 0.5 < mi < hi_ + r + 0.5 and beruehrt(t["bb"], teile[j]["bb"]):
                folgt[i] = q; neu_.append(i)
    for i_b, v in band.items():                              # Gaze usw. folgt dem Band des Deckels, in dem sie liegt
        for i, t in enumerate(teile):
            if i in folgt or i in mit or arten.get(i) not in ("ganz", None) or not GANZ.search(t["name"]): continue
            if not beruehrt(t["bb"], teile[i_b]["bb"]): continue
            mj = (t["bb"][jj_] + t["bb"][jj_ + 3]) / 2
            for a_, b_, q_ in v:
                if a_ <= mj <= b_:
                    lo_, hi_ = sorted((p, q_)); mi = (t["bb"][k] + t["bb"][k + 3]) / 2
                    if lo_ - 0.5 < mi < hi_ + 0.5: folgt[i] = q_
                    break
    dd = d if r == 0 else -r
    for i, t in enumerate(teile):
        q = folgt.get(i, p); a = arten.get(i)
        lo, hi = t["bb"][k], t["bb"][k + 3]; mi = (lo + hi) / 2
        grenze = q + (r if r > 0 else 0.0)
        if i in mit:
            if mit[i] == "oben":
                t["s"] = verschiebe(t["s"], k, dd); t["art"] = t.get("art", "") + f"M{k}"
            if (mit[i] == "oben") != (mi > grenze) or a == "ganz": bericht.append(f"{t['name']}: " + ("wandert mit der Wand" if mit[i] == "oben" else "bleibt bei der Wand"))
        elif a == "band":
            tb = t["bb"]; stuecke = []
            for a_, b_, q_ in band[i]:
                lo_ = list(tb[:3]); hi_ = list(tb[3:]); lo_[jj_], hi_[jj_] = a_, b_
                c = t["s"].intersect(cq.Solid.makeBox(hi_[0] - lo_[0], hi_[1] - lo_[1], hi_[2] - lo_[2] + 10, cq.Vector(lo_[0], lo_[1], lo_[2] - 5)))
                if c.Volume() > 1e-6: stuecke.append(schneide({"s": c, "bb": bbox(c), "name": t["name"]}, k, q_, 0.0, d))
            neu = stuecke[0]
            for x in stuecke[1:]: neu = neu.fuse(x)
            try: neu = neu.clean()
            except Exception: pass
            sol = neu.Solids(); t["s"] = sol[0] if len(sol) == 1 else neu
            t["art"] = t.get("art", "") + f"B{k}"
        elif a in ("schnitt", "rand"):
            t["s"] = schneide(t, k, q, r, d); t["art"] = t.get("art", "") + f"S{k}{q:.1f}"
            bericht.append(f"{t['name']}: geschnitten bei {('XYZ')[k]} {q:.1f}" + (f" (eigene Stelle statt {p:.1f})" if abs(q - p) > 0.01 else ""))
        elif a == "horn":                                    # mit der Aussenseite versetzen, innen auf das alte Ende kuerzen
            neu = verschiebe(t["s"], k, dd)
            if dd < 0:
                kb = kasten(k, lo, t["bb"][k + 3] + 1, t["bb"]); neu = neu.intersect(kb)
            t["s"] = neu; t["art"] = t.get("art", "") + f"S{k}{q:.1f}"
            bericht.append(f"{t['name']}: versetzt und innen gekuerzt")
        elif (a == "ganz" and mi > grenze) or (a is None and lo >= q - 1e-6):
            t["s"] = verschiebe(t["s"], k, dd); t["art"] = t.get("art", "") + f"M{k}"
        t["bb"] = bbox(t["s"])

def suche_einfuegen(teile, k, von, bis, schritt):
    u = von
    gruende = []
    while (u >= bis if schritt < 0 else u <= bis):
        ok, arten, grund = pruefe_ebene(teile, k, u, 0.0, False)
        if ok: return u, arten
        gruende.append(grund); u += schritt
    raise MassFehler(f"keine Stelle zum Einfuegen zwischen {von:.0f} und {bis:.0f} gefunden (z. B. {gruende[0] if gruende else '-'})")

def suche_streifen(teile, k, von, bis, bedarf):
    """freie, gleichmaessige Bereiche in [von, bis] -> Liste (p, r) mit Summe r = bedarf (groesste zuerst)"""
    gut = []
    u = von
    while u < bis:
        ok, _, _ = pruefe_ebene(teile, k, u, 0.5, True); gut.append((u, ok)); u += 0.5
    laeufe, start = [], None
    for u, ok in gut + [(bis, False)]:
        if ok and start is None: start = u
        if not ok and start is not None: laeufe.append((start, u - start + 0.0)); start = None
    kandidaten = []
    for p0, l0 in sorted(laeufe, key=lambda e: -e[1]):
        r = l0
        while r > 0.4:                                       # groessten pruefbaren Streifen im Lauf (Mindestmasse)
            ok, arten, grund = pruefe_ebene(teile, k, p0, r, True)
            if ok: kandidaten.append((p0, r)); break
            r -= 0.5
    gewaehlt, summe, weg = [], 0.0, {}                       # weg: schon entfernte mm je Teil (Mindestmasse ueber alle Streifen)
    moeglich = 0.0
    for p0, r in sorted(kandidaten, key=lambda e: -e[1]):
        ok, arten, _ = pruefe_ebene(teile, k, p0, r, True)
        if not ok: continue
        grenze = r
        for i, a in arten.items():
            if a != "schnitt": continue
            t = teile[i]; e = ausdehnung(t, k); mind = MIN_WAND if WAND.search(t["name"]) else (MIN_LEISTE if e < 30 else 10.0)
            grenze = min(grenze, e - weg.get(i, 0.0) - mind)
        grenze = math.floor(max(grenze, 0.0) * 2) / 2
        if grenze < 0.5: continue
        moeglich += grenze
        if summe >= bedarf - 1e-6: continue
        nimm = min(grenze, bedarf - summe); gewaehlt.append((p0, round(nimm, 2))); summe += nimm
        for i, a in arten.items():
            if a == "schnitt": weg[i] = weg.get(i, 0.0) + nimm
    return gewaehlt, summe, moeglich

def kollisionen(teile, alt, nur):
    """Beruehrungen (> 0,5 mm3) der veraenderten Teile mit allen anderen - vorher (alt) und nachher"""
    def paare(liste):
        bbs = [bbox(s) for s in liste]; erg = {}
        for i in nur:
            for j in range(len(liste)):
                if j == i or (j in nur and j < i): continue
                a, b = bbs[i], bbs[j]
                if any(a[m] > b[m + 3] or b[m] > a[m + 3] for m in range(3)): continue
                if teile[i].get("art") == teile[j].get("art") and teile[i].get("art"): continue   # gleich bewegt
                try: v = liste[i].intersect(liste[j]).Volume()
                except Exception: v = 0.0
                if v > 0.5 and not ("Gaze_" in teile[i]["name"] + teile[j]["name"] and v < 20.0):   # Gaze-Tuch: Modellueberlappung
                    erg[tuple(sorted((teile[i]["name"], teile[j]["name"])))] = round(v, 2)
        return erg
    vor = paare(alt); nach = paare([t["s"] for t in teile])
    return {p: v for p, v in nach.items() if p not in vor or v > 2 * vor[p] + 0.5}

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("ein"); ap.add_argument("aus")
    ap.add_argument("--teil", required=True); ap.add_argument("--breite", type=float); ap.add_argument("--laenge", type=float)
    ap.add_argument("--bericht")
    ar = ap.parse_args()
    t0 = time.time()
    assy = cq.Assembly.importStep(ar.ein)
    teile = blaetter(assy)
    for t in teile: t["bb"] = bbox(t["s"])
    alt = [t["s"] for t in teile]
    wand = [t for t in teile if WAND.search(t["name"])]
    if not wand: raise MassFehler("keine Gehaeusewaende im Modell gefunden")
    X0, X1 = min(t["bb"][0] for t in wand), max(t["bb"][3] for t in wand)
    Y0, Y1 = min(t["bb"][1] for t in wand), max(t["bb"][4] for t in wand)
    L0, B0 = round(X1 - X0 - 2.0, 1), round(Y1 - Y0 - 2.0, 1)       # Fertigmass (Rohteile 2 mm Uebermass)
    bericht = {"vorher": {"laenge": L0, "breite": B0}, "schritte": []}
    print(f"Modell {os.path.basename(ar.ein)}: {len(teile)} Teile, Gehaeuse {L0:g} x {B0:g} mm (fertig)", flush=True)
    schritte = []
    if ar.breite is not None and abs(ar.breite - B0) > 0.05:
        if not 200.0 - 1e-6 <= ar.breite <= 215.0 + 1e-6: raise MassFehler(f"Breite {ar.breite:g} mm ausserhalb 200-215 mm")
        dB = round(ar.breite - B0, 2)
        if dB > 0:
            p, arten = suche_einfuegen(teile, 1, Y1 - 6.0, Y1 - 120.0, -0.5)
            aendere(teile, 1, p, 0.0, dB, arten, schritte, "Breite")
            bericht["schritte"].append({"breite": f"+{dB:g} bei Y {p:.1f}"})
            print(f"Breite +{dB:g} mm: Ebene Y {p:.1f} (vorne)", flush=True)
        else:
            # gesucht wird nahe den Laengswaenden (Leisten, Waende, Freiraum vor/hinter den Stimmplatten); Diskant nur vorne
            zonen = [(Y1 - 25.0, Y1 - 0.5)] + ([(Y0 + 0.5, Y0 + 25.0)] if ar.teil == "bass" else [])
            streifen, summe, moeglich = [], 0.0, 0.0
            for von, bis in zonen:
                st_, su_, mo_ = suche_streifen(teile, 1, von, bis, -dB - summe)
                streifen += st_; summe += su_; moeglich += mo_
            if summe < -dB - 1e-6:
                raise MassFehler(f"Breite {ar.breite:g} mm geht bei diesem Modell nicht: hoechstens {moeglich:.1f} mm schmaeler "
                                 f"(= {B0 - moeglich:.1f} mm), ohne Stimmplatten, Kammern oder Mechanik zu beschneiden")
            for p, r in sorted(streifen, key=lambda e: -e[0]):              # von vorne nach hinten (Koordinaten bleiben gueltig)
                ok, arten, grund = pruefe_ebene(teile, 1, p, r, True)
                if not ok: raise MassFehler(f"Streifen Y {p:.1f}+{r:g}: {grund}")
                aendere(teile, 1, p, r, 0.0, arten, schritte, "Breite")
                bericht["schritte"].append({"breite": f"-{r:g} bei Y {p:.1f}"})
                print(f"Breite -{r:g} mm: Streifen Y {p:.1f}-{p + r:.1f}", flush=True)
    if ar.laenge is not None and abs(ar.laenge - L0) > 0.05:
        dL = round(ar.laenge - L0, 2)
        if dL < 0: raise MassFehler(f"Laenge {ar.laenge:g} mm ist kuerzer als gezeichnet ({L0:g} mm) - kuerzer geht nicht")
        if dL > 150: raise MassFehler(f"Laenge {ar.laenge:g} mm: hoechstens 150 mm laenger als gezeichnet")
        h = dL / 2
        X1n = max(t["bb"][3] for t in teile if WAND.search(t["name"]))
        p, arten = suche_einfuegen(teile, 0, X1n - 6.0, X1n - 120.0, -0.5)          # rechts zuerst (links bleibt gueltig)
        aendere(teile, 0, p, 0.0, h, arten, schritte, "Laenge rechts")
        X0n = min(t["bb"][0] for t in teile if WAND.search(t["name"]))
        p2, arten = suche_einfuegen(teile, 0, X0n + 6.0, X0n + 120.0, 0.5)
        aendere(teile, 0, p2, 0.0, h, arten, schritte, "Laenge links")
        for t in teile: t["s"] = verschiebe(t["s"], 0, -h); t["bb"] = bbox(t["s"]); t["art"] = t.get("art", "") + "V"
        bericht["schritte"].append({"laenge": f"+{h:g} bei X {p2:.1f} und {p:.1f}"})
        print(f"Laenge +{dL:g} mm: je {h:g} bei X {p2:.1f} und X {p:.1f}", flush=True)
    for z in schritte: print("  ", z)
    if not schritte:
        print("keine Aenderung"); bericht["kollisionen_neu"] = {}
    else:
        nur = [i for i, t in enumerate(teile) if "S" in t.get("art", "")] + \
              [i for i, t in enumerate(teile) if "S" not in t.get("art", "") and t.get("art")]
        nur = sorted(set(nur))
        neu = kollisionen(teile, alt, nur)
        bericht["kollisionen_neu"] = {" | ".join(k): v for k, v in neu.items()}
        if neu:
            for k_, v in neu.items(): print("NEUE BERUEHRUNG", k_[0], "|", k_[1], v)
        else: print("neue Beruehrungen: 0", flush=True)
    wand = [t for t in teile if WAND.search(t["name"])]
    bericht["nachher"] = {"laenge": round(max(t["bb"][3] for t in wand) - min(t["bb"][0] for t in wand) - 2.0, 1),
                          "breite": round(max(t["bb"][4] for t in wand) - min(t["bb"][1] for t in wand) - 2.0, 1)}
    print(f"Gehaeuse nachher {bericht['nachher']['laenge']:g} x {bericht['nachher']['breite']:g} mm (fertig)")
    if ar.bericht: json.dump(bericht, open(ar.bericht, "w", encoding="utf-8"), ensure_ascii=False, indent=1)
    if bericht.get("kollisionen_neu"): raise MassFehler("Masse geaendert, aber es gibt neue Beruehrungen (siehe oben) - nicht gespeichert")
    je = {}
    for t in teile: je.setdefault(id(t["knoten"]), (t["knoten"], []))[1].append(t["s"].moved(t["loc"].inverse))
    for kn, objs in je.values():
        kn.obj = objs[0] if len(objs) == 1 else cq.Workplane("XY").newObject(objs)
        for att in ("_subshape_names", "_subshape_colors", "_subshape_layers"):
            if hasattr(kn, att) and getattr(kn, att): setattr(kn, att, {})
    assy.save(ar.aus)
    print(f"gespeichert {ar.aus} ({time.time() - t0:.0f} s)")

if __name__ == "__main__":
    main()
