#!/usr/bin/env python3
"""Instrument-Generator: aus einer Tastenliste (Textdatei) werden je Reihe
Knopf-X, Klappen (Mitte, Lage), Bodenloecher und Stimmstock-Kammern erzeugt.

Aufruf:  python3 instrument.py <tastenliste.txt> [ausgabe.json]

Format der Tastenliste (Kommentare mit #):
    name           = Markus
    choerig        = 3          # 3: Hebel auf Trippel (3 Loecher), 2: zwei Loecher, Kammern gegenueber
    knopfabstand   = 20         # 18 / 19 / 20 ...
    klappenteilung = 26.25      # bei choerig=2 = knopfabstand
    gehaeuse_x     = 380        # Innenlaenge fix, oder "auto"
    rand           = 17.5       # Abstand Klappenrand -> Gehaeuse bei auto
    rand_links / rand_rechts    # optional getrennt (sonst = rand)
    verdeck_muster = keines     # oder verdeck.dxf (Markus) / geschwungen.dxf (Marie, im Ordner muster); Standard: glatt
    boden_dxf      = muster/marie_boden.dxf   # 2-choerig: Klappenumriss und Griffstockfuesse aus dieser Bodenzeichnung
    kammer_breite  = 15 / kammer_tiefe = 9.5  # 2-choerig
    [reihe n]: ton = G, neigung = 0           # optional
    streck_x       = 190        # optional: wo Rahmen/Verdeck verlaengert werden (Koordinate im 380er-Grundmodell)
    [reihe 1]
    knopf_x0   = 64.99          # erster Knopf
    knoepfe    = 14
    klappe_x0  = 32.36          # erste Klappe (Mitte)
    klappe_y   = 169.66         # Reihenmitte der Klappen (Zickzack +-2)
    erste      = unten          # erste Klappe unten/oben (3-choerig: Bodenlage *U4 = unten)
    kammer     = h              # Seite der ersten Kammer (3-choerig)
    platten    = 32 33 34 ...   # eine Plattenlaenge je Klappe (von links) -> Anzahl Klappen
    knopf_versatz = 10          # optional: Knopf samt Knopfschenkel seitlich (+X) versetzt, Hebel Z-gekroepft
    extra      = 4              # optional: erster Knopf dieser Reihe treibt die erste Klappe von Reihe 4
    knopf_ab   = 2              # optional: erster benutzter Knopfplatz (1 = erster)
Vorne/hinten anfuegen: Wert in "platten" vorne/hinten ergaenzen, bei vorne
klappe_x0 um eine Klappenteilung verkleinern; mit gehaeuse_x = auto wird
das Gehaeuse automatisch angepasst.
"""
import json, re, sys

ZICK = 2.0                 # Klappen-Y Zickzack +-2 mm (3-choerig)
KAMMER_3 = 8.75            # Kammerteilung im Trippel
LOCH_3 = {"unten": [(0, 8), (-6, -4), (6, -4)],    # *U4
          "oben":  [(0, -8), (-6, 4), (6, 4)]}     # *U3
LOCH_2 = [(0, 0), (0, 13)]                         # 2chor-klappe
KLAPPE_HALB = {3: 14.77, 2: 9.0}


def _klappen_umrisse():
    """Klappenumriss je Lage relativ zur Klappenmitte (aus boden.dxf: *U4 = unten, *U3 = oben); 2-choerig Rechteck."""
    import os, ezdxf
    from shapely.geometry import Polygon, box as sbox
    u = {"-": sbox(-9.0, -9.0, 9.0, 22.0)}
    f = os.path.join(os.path.dirname(os.path.abspath(__file__)), "boden.dxf")
    if not os.path.exists(f): return u
    doc = ezdxf.readfile(f)
    for name, lage, dy in (("*U4", "unten", 8.0), ("*U3", "oben", -8.0)):
        if name not in doc.blocks: continue
        for e in doc.blocks[name]:
            if e.dxf.layer == "Klappen-rand-aussen":
                p = [(q.x, q.y + dy) for q in ezdxf.path.make_path(e).flattening(0.05)]
                u[lage] = Polygon(p).buffer(0)
    return u


def lies(pfad):
    inst, reihe = {"reihen": {}}, None
    for z in open(pfad, encoding="utf-8"):
        z = z.split("#")[0].strip()
        if not z:
            continue
        m = re.match(r"\[\s*reihe\s+(\d+)\s*\]", z, re.I)
        if m:
            reihe = int(m.group(1)); inst["reihen"][reihe] = {}; continue
        k, v = [s.strip() for s in z.split("=", 1)]
        ziel = inst if reihe is None else inst["reihen"][reihe]
        if k == "platten":
            ziel[k] = [float(x) for x in v.replace(",", " ").split()]
        else:
            try: ziel[k] = float(v)
            except ValueError: ziel[k] = v
    inst["choerig"] = int(inst.get("choerig", 3))
    if "klappenteilung" not in inst:
        inst["klappenteilung"] = inst["knopfabstand"] if inst["choerig"] == 2 else 26.25
    return inst


def erzeuge(inst):
    ch, T = inst["choerig"], inst["klappenteilung"]
    out = {"name": inst.get("name", ""), "choerig": ch, "reihen": {}, "knopfabstand": float(inst["knopfabstand"])}
    if "streck_x" in inst: out["streck_x"] = float(inst["streck_x"])
    out["verdeck_muster"] = str(inst.get("verdeck_muster", "keines"))
    for k_ in ("boden_dxf",):
        if k_ in inst: out[k_] = str(inst[k_])
    out["griffstock"] = str(inst.get("griffstock", ""))   # 20.11.: "traditionell" = Holzhebel auf zwei durchgehenden Achsen
    out["feder"] = str(inst.get("feder", ""))        # Drehfedern an den Achsen: "" = keine, "alle", oder Hebelname z. B. Hebel_R1_07
    out["kammer_breite"] = float(inst.get("kammer_breite", 15.0))     # nur 2-choerig (3-choerig: Regeln nach Jakob)
    out["kammer_tiefe"] = float(inst.get("kammer_tiefe", 9.5))
    for r, R in sorted(inst["reihen"].items()):
        n = len(R["platten"])
        unten = str(R.get("erste", "unten")) == "unten"
        seite = str(R.get("kammer", "h"))
        kl, lo, ka = [], [], []
        for i in range(n):
            x = R["klappe_x0"] + i * T
            L = R["platten"][i]
            if ch == 3:
                lage = "unten" if (i % 2 == 0) == unten else "oben"
                y = R["klappe_y"] + (-ZICK if lage == "unten" else ZICK)
                kl.append({"x": x, "y": y, "lage": lage, "L": L})
                lo += [(x + a, y + b) for a, b in LOCH_3[lage]]
                for j, d in enumerate((-KAMMER_3, 0.0, KAMMER_3)):
                    s = seite if (3 * i + j) % 2 == 0 else ("v" if seite == "h" else "h")
                    ka.append({"x": x + d, "seite": s, "L": L, "klappe": i})
            else:
                y = R["klappe_y"]
                kl.append({"x": x, "y": y, "lage": "-", "L": L})
                lo += [(x + a, y + b) for a, b in LOCH_2]
                ka += [{"x": x, "seite": "v", "L": L, "klappe": i},
                       {"x": x, "seite": "h", "L": L, "klappe": i}]
        kn = [R["knopf_x0"] + j * inst["knopfabstand"] for j in range(int(R.get("knoepfe", n)))]
        out["reihen"][r] = {"klappen": kl, "loecher": lo, "kammern": ka, "knoepfe": kn,
                            "extra": int(R.get("extra", 0)), "knopf_ab": int(R.get("knopf_ab", 1)),
                            "ton": str(R.get("ton", "R%d" % r)), "neigung": float(R.get("neigung", 0.0)),
                            "knopf_versatz": float(R.get("knopf_versatz", 0.0))}
    # Gehaeuse
    h = KLAPPE_HALB[ch]
    xmin = min(k["x"] for R in out["reihen"].values() for k in R["klappen"]) - h
    xmax = max(k["x"] for R in out["reihen"].values() for k in R["klappen"]) + h
    # Verdeckstuetzen innen an den Seitenwaenden (Mitte 16,55 vom Rand, Platte 3 mm, Fussplatten bei Y 90,8-101,8 und
    # 119,8-130,8 bis Z -8 = Hoehe der Klappen, Sechskantkoepfe der Fussschrauben bei Y 93,3 / 128,3): Klappenumriss
    # (aus boden.dxf) gegen diese Hindernisse pruefen und den Rand so weit vergroessern, dass 0,5 mm Luft bleibt.
    ST_X, LUFT = 16.55, 0.5
    st_l = st_r = None
    try:
        from shapely.geometry import Polygon, box as sbox, Point
        from shapely.ops import unary_union
        from shapely import affinity
        hind = unary_union([sbox(-1.5, 90.8, 1.5, 101.8), sbox(-1.5, 119.8, 1.5, 130.8),
                            Point(0, 93.31).buffer(2.9), Point(0, 128.31).buffer(2.9)]).buffer(LUFT)
        umr = _klappen_umrisse()
        kp = [affinity.translate(umr.get(k["lage"], umr.get("-")), k["x"], k["y"])
              for R in out["reihen"].values() for k in R["klappen"]]
        kp = unary_union([q for q in kp if q is not None])
        def frei(t): return not kp.intersects(affinity.translate(hind, t, 0))
        def suche(a, b):                  # kleinstes/groesstes t zwischen a (frei) und b (belegt)
            for _ in range(40):
                m = (a + b) / 2
                if frei(m): a = m
                else: b = m
            return a
        xs = [k["x"] for R in out["reihen"].values() for k in R["klappen"]]
        xr, xl = max(xs) + 60, min(xs) - 60
        if not frei(xr) or not frei(xl): raise ValueError("Suchbereich")
        st_r = suche(xr, max(xs) - 60) if not frei(max(xs) - 60) else None    # Stuetzenmitte rechts mindestens hier
        st_l = suche(xl, min(xs) + 60) if not frei(min(xs) + 60) else None    # Stuetzenmitte links hoechstens hier
    except Exception as ex_:
        print("Hinweis: Stuetzenpruefung nicht moeglich:", ex_)
    rand = float(inst.get("rand", 17.5))
    rl = float(inst.get("rand_links", rand)); rr = float(inst.get("rand_rechts", rand))
    if str(inst.get("gehaeuse_x", "auto")) == "auto":
        verschub = rl - xmin
        if st_l is not None: verschub = max(verschub, ST_X - st_l)            # linke Stuetze liegt bei X 16,55
        if st_r is not None: rr = max(rr, st_r + ST_X - xmax)                 # rechte Stuetze bei L - 16,55
        rl = verschub + xmin
        for R in out["reihen"].values():
            for k in R["klappen"]: k["x"] += verschub
            for k in R["kammern"]: k["x"] += verschub
            R["loecher"] = [(x + verschub, y) for x, y in R["loecher"]]
            R["knoepfe"] = [x + verschub for x in R["knoepfe"]]
        out["gehaeuse_x"] = xmax - xmin + rl + rr
        out["verschub"] = verschub
    else:
        out["gehaeuse_x"] = float(inst["gehaeuse_x"])
        out["verschub"] = 0.0
        out["rand_links"], out["rand_rechts"] = xmin, out["gehaeuse_x"] - xmax
        if xmin < 0 or xmax > out["gehaeuse_x"]:
            print("WARNUNG: Klappen ragen ueber das Gehaeuse", xmin, xmax)
        if st_l is not None and st_l < ST_X: print("WARNUNG: Klappe links an der Verdeckstuetze")
        if st_r is not None and st_r > out["gehaeuse_x"] - ST_X: print("WARNUNG: Klappe rechts an der Verdeckstuetze")
    # 22.11.2026 (Johann): Griffstock/Knopfpyramide immer mittig am Korpus, abhaengig von der konfigurierten Knopfanzahl
    # (wirksame Knopf-X inkl. Knopfversatz; Klappen bleiben, wo sie sind)
    _kx = [x + R["knopf_versatz"] for R in out["reihen"].values() for x in R["knoepfe"]]
    _dm = out["gehaeuse_x"] / 2 - (min(_kx) + max(_kx)) / 2
    if abs(_dm) > 0.01:
        for R in out["reihen"].values(): R["knoepfe"] = [x + _dm for x in R["knoepfe"]]
        print("   Knopfpyramide mittig am Korpus: um %+.2f mm verschoben" % _dm)
    return out


if __name__ == "__main__":
    inst = lies(sys.argv[1])
    g = erzeuge(inst)
    print(g["name"], "%d-choerig" % g["choerig"], "Gehaeuse X %.2f" % g["gehaeuse_x"])
    for r, R in g["reihen"].items():
        print(" Reihe %d: %d Knoepfe, %d Klappen, %d Loecher, %d Kammern"
              % (r, len(R["knoepfe"]), len(R["klappen"]), len(R["loecher"]), len(R["kammern"])))
    if len(sys.argv) > 2:
        json.dump(g, open(sys.argv[2], "w"), indent=1)
