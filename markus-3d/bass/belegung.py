# Freie Bassbelegung (Johann 26.11.2026, Regeln in REGELN_Bassbelegung.md)
# Liest eine Belegungsdatei (welcher Knopf welche Klappe oeffnet; gleiche Klappe bei zwei Knoepfen = Kopplung),
# prueft sie und verteilt die Achsen automatisch auf die Achsgassen des Modells. Ergebnis ist die Achstabelle _A,
# die die vorhandene Mechanik-Berechnung (mech_*.py) wie bisher baut und auf Kollisionen prueft.
#
# Format der Belegungsdatei:
#   modell = v340
#   1.1 = V077        # Reihe 1, Knopf 1 (von links) oeffnet Klappe V077
#   2.1 = H218        # zwei Knoepfe auf derselben Klappe = Kopplung
#   R1 = 10           # optional: Knoepfe in Reihe 1 (ohne Angabe wie gezeichnet); die Reihe wird mittig neu verteilt
#   R2 = 9
# Kurzform (Webseite / GitHub-Formular): "R1=10, R2=9, 1.1=V077, 1.2=BL152, 2.1=H218, ..."
# Knoepfe ohne Zuordnung sind blind (Knopf ohne Mechanik) und zaehlen nicht: verdrahtete Knoepfe = Klappen + Kopplungen.
# Die Klappen, Gassen und Knopflagen eines Modells stehen in <modell>/belegung_modell_<name>.json.
#
# Aufruf zur Pruefung ohne CAD:  python3 belegung.py <modell.json> <belegungsdatei oder Kurzform>
import re, os, json, sys

class BelegungsFehler(SystemExit):
    pass

def lies(pfad):
    kopf, zuo = {}, {}
    text = open(pfad, encoding="utf-8").read() if os.path.exists(pfad) else pfad
    zeilen = [t for z in text.splitlines() for t in (z.split("#", 1)[0].replace(";", ",").split(",") if "=" in z.split("#", 1)[0] and z.count("=") > 1 else [z])]
    for nr, z in enumerate(zeilen, 1):
        z = z.split("#", 1)[0].strip()
        if not z: continue
        m = re.match(r"^(\d+)\s*\.\s*(\d+)\s*=\s*(\S+)$", z)
        if m:
            k = (int(m.group(1)), int(m.group(2)))
            if k in zuo: raise BelegungsFehler(f"Belegung Zeile {nr}: Knopf {k[0]}.{k[1]} doppelt")
            zuo[k] = m.group(3).upper(); continue
        m = re.match(r"^(\w+)\s*=\s*(.+)$", z)
        if m: kopf[m.group(1).lower()] = m.group(2).strip(); continue
        raise BelegungsFehler(f"Belegung Zeile {nr}: '{z}' nicht verstanden (Form: 1.3 = H124)")
    return kopf, zuo

def reihen_aus(kopf, reihen):
    """Knopfzahl je Reihe: wie gezeichnet, ueberschrieben durch R1 = n / R2 = n im Kopf der Belegung"""
    neu = dict(reihen)
    for k, v in kopf.items():
        m = re.fullmatch(r"r(?:eihe)?(\d+)", k)
        if not m: continue
        r = int(m.group(1))
        if r not in reihen: raise BelegungsFehler(f"Reihe {r} gibt es nicht (Reihen {', '.join(map(str, sorted(reihen)))})")
        if not re.fullmatch(r"\d+", v) or not 1 <= int(v) <= 20: raise BelegungsFehler(f"R{r} = {v}: Knopfzahl 1-20 erwartet")
        neu[r] = int(v)
    return neu

def knopf_lagen(m, reihen, teilung=None):
    """Knopf-X je Reihe. Die Reihen bleiben gegeneinander versetzt wie gezeichnet (Raster je Reihe, Versatz waechst mit der
    Teilung); jede Reihe sitzt so mittig wie auf ihrem Raster moeglich, dann wird das ganze Knopffeld auf die gezeichnete
    Mitte geschoben. Bei 9/9 ergibt das genau die Zeichnung, bei 10/9 sind beide Reihen exakt mittig."""
    t0 = float(m["teilung"]); t = float(teilung) if teilung not in (None, "", "standard") else t0
    mi = {int(r): x for r, x in m["reihen_mitte"].items()}
    c = (min(mi.values()) + max(mi.values())) / 2
    rs = sorted(reihen); o = {r: (mi[r] - c) * t / t0 for r in rs}     # Rasterversatz jeder Reihe gegen die Mitte
    folge = sorted(rs, key=lambda r: mi[r])                            # gezeichnete Reihenfolge von links
    def wertung(wahl):                                                 # schmal, jede Reihe mittig, Reihenfolge wie gezeichnet
        lagen = {r: [o[r] + (k0 + i) * t for i in range(reihen[r])] for r, k0 in zip(rs, wahl)}
        alle = [x for xs in lagen.values() for x in xs]; mitte = (min(alle) + max(alle)) / 2
        asym = sum(abs((xs[0] + xs[-1]) / 2 - mitte) for xs in lagen.values())
        anders = sorted(rs, key=lambda r: (lagen[r][0] + lagen[r][-1]) / 2) != folge
        return (round(max(alle) - min(alle), 6), round(asym, 6), anders), lagen, mitte
    import itertools
    _, kx, mitte = min((wertung(w) for w in itertools.product(*[range(-reihen[r] - 1, 2) for r in rs])), key=lambda e: e[0])
    sh = -mitte
    kx = {r: [round(c + x + sh, 2) for x in xs] for r, xs in kx.items()}
    x0, x1 = m.get("knopf_bereich", [-1e9, 1e9])
    for r, xs in kx.items():
        if xs and (xs[0] < x0 or xs[-1] > x1):
            raise BelegungsFehler(f"Reihe {r}: {len(xs)} Knoepfe mit {t:g} mm Teilung reichen von X {xs[0]:.0f} bis {xs[-1]:.0f} - "
                                  f"an der Front ist nur Platz von {x0:g} bis {x1:g} (Hoerner, Seitenwaende); weniger Knoepfe oder kleinere Teilung")
    return kx

def pruefe(zuo, katalog, reihen):
    """Zahlen pruefen: jeder verdrahtete Knopf genau eine Klappe, jede Klappe mindestens einmal,
    verdrahtete Knoepfe = Klappen + Kopplungen; Knoepfe ohne Zuordnung sind blind und zaehlen nicht"""
    fehler = []
    alle = {(r, i) for r, n in reihen.items() for i in range(1, n + 1)}
    for k in sorted(set(zuo) - alle): fehler.append(f"Knopf {k[0]}.{k[1]} gibt es nicht (Reihen: {', '.join(f'{r}: {n} Knoepfe' for r, n in reihen.items())})")
    for k, kl in zuo.items():
        if kl not in katalog: fehler.append(f"Knopf {k[0]}.{k[1]}: Klappe {kl} gibt es nicht (vorhanden: {' '.join(sorted(katalog))})")
    benutzt = set(zuo.values())
    for kl in sorted(set(katalog) - benutzt): fehler.append(f"Klappe {kl} wird von keinem Knopf geoeffnet")
    if fehler: raise BelegungsFehler("Belegung nicht moeglich:\n  " + "\n  ".join(fehler))
    n_kopp = len(zuo) - len(katalog)
    return n_kopp, sorted(alle - set(zuo))

def verteile(zuo, katalog, gassen, knopf_x, min_luecke=12.0):
    """Achsen auf Gassen verteilen. katalog[kl] = {"ziel": (x, y), "hebel_x": x, "gassen": [namen in Vorzugsreihenfolge]}"""
    achsen = {}
    for (r, i), kl in zuo.items(): achsen.setdefault(kl, []).append((r, i))
    def bereich(kl):
        xs = [knopf_x(r, i - 1) for r, i in achsen[kl]] + [katalog[kl]["hebel_x"]]
        return min(xs), max(xs)
    reihenfolge = sorted(achsen, key=lambda kl: (len(katalog[kl]["gassen"]), bereich(kl)[1] - bereich(kl)[0]))
    # wenig Auswahl zuerst, dann kurze Achsen (sie bekommen die Gasse nahe der Klappe, lange/gekoppelte weichen aus);
    # geht es nicht auf, wird zurueckgegangen und eine andere Verteilung probiert (Suche mit Ruecksprung)
    def frei(belegt, g, x0, x1):
        return [n for a, b, n in belegt.get(g, []) if not (x1 + min_luecke <= a or b + min_luecke <= x0)]
    tiefste = [0, None]
    def suche(k, belegt, wahl):
        if k == len(reihenfolge): return wahl
        kl = reihenfolge[k]; x0, x1 = bereich(kl)
        for g in katalog[kl]["gassen"]:
            if frei(belegt, g, x0, x1): continue
            belegt.setdefault(g, []).append((x0, x1, kl)); wahl[kl] = g
            erg = suche(k + 1, belegt, wahl)
            if erg: return erg
            belegt[g].pop(); del wahl[kl]
        if k >= tiefste[0]: tiefste[0], tiefste[1] = k, kl
        return None
    wahl = suche(0, {}, {})
    if not wahl:
        kl = tiefste[1]; x0, x1 = bereich(kl)
        raise BelegungsFehler(f"Belegung nicht baubar: Achse {kl} (Knoepfe {' '.join(f'{r}.{i}' for r, i in sorted(achsen[kl]))}, "
                              f"X {x0:.0f}-{x1:.0f}) findet in keiner Verteilung eine freie Gasse ({', '.join(katalog[kl]['gassen'])}); "
                              "andere Knoepfe fuer diese Klappe waehlen")
    bericht, _A = [], []
    for kl in reihenfolge:
        x0, x1 = bereich(kl); g = wahl[kl]; kn = sorted(achsen[kl])
        _A.append((kl, gassen[g], [(r, knopf_x(r, i - 1), knopf_x(r, i - 1)) for r, i in kn],
                   [(katalog[kl]["hebel_x"], katalog[kl]["ziel"])], None))
        bericht.append(f"Achse {kl:6s} Gasse {g:4s} X {x0:6.1f}-{x1:6.1f}  Knoepfe {' '.join(f'{r}.{i}' for r, i in kn)}"
                       + ("  (gekoppelt)" if len(kn) > 1 else ""))
    return _A, bericht

def baue_tabelle(pfad, katalog, gassen, knopf_x, reihen, min_luecke=12.0):
    """reihen = Knopfzahl je Reihe wie gebaut (Kopfangaben R1/R2 der Belegung muessen schon darin stecken)"""
    kopf, zuo = lies(pfad)
    n_kopp, blind = pruefe(zuo, katalog, reihen)
    _A, bericht = verteile(zuo, katalog, gassen, knopf_x, min_luecke)
    print(f"Belegung: {len(zuo)} verdrahtete Knoepfe = {len(katalog)} Klappen + {n_kopp} Kopplungen"
          + (f"; blind (ohne Mechanik): {' '.join(f'{r}.{i}' for r, i in blind)}" if blind else ""))
    for z in bericht: print(z)
    return _A

def lade_modell(pfad_json, belegung=None, teilung=None):
    """Modell laden; mit belegung (Datei/Kurzform) gelten deren Knopfzahlen R1/R2, mit teilung die Knopfteilung"""
    m = json.load(open(pfad_json, encoding="utf-8"))
    katalog = {k: {"ziel": tuple(v["ziel"]), "hebel_x": v["hebel_x"], "gassen": v["gassen"]} for k, v in m["klappen"].items()}
    reihen = {int(r): n for r, n in m["reihen"].items()}
    if belegung: reihen = reihen_aus(lies(belegung)[0], reihen)
    kx = knopf_lagen(m, reihen, teilung)
    return m, katalog, m["gassen"], (lambda r, i: kx[r][i]), reihen

def pruefe_alles(pfad_json, belegung, teilung=None):
    m, katalog, gassen, knopf_x, reihen = lade_modell(pfad_json, belegung, teilung)
    return baue_tabelle(belegung, katalog, gassen, knopf_x, reihen, m.get("min_luecke", 12.0))

if __name__ == "__main__":
    # python3 belegung.py <modell.json> <belegung> [teilung]
    pruefe_alles(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else None)
    print("Belegung baubar.")
