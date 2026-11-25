# Bassteil-Variante vorbereiten (23.11.2026): Knopfteilung je Reihe, gleicher Reihenabstand.
# Aufruf: python3 bass_variante.py <quellordner> <zielordner> <teilung> [reihenabstand]
# Schreibt in den Zielordner eine Kopie der Bass-Quellen, in der rahmen.dxf (Bohrungen Bassdruecker) und
# blende.dxf (4 Bohrungen der kurzen Reihe) auf die neue Teilung umgesetzt sind; Knopfreihen bleiben mittig
# (Mitte wie im Original, X 191). mechanik.py ordnet die Schubstangen je Reihe nach Rang statt nach Abstand zu.
import sys, os, shutil, ezdxf
Q, Z, T = sys.argv[1], sys.argv[2], float(sys.argv[3])
RA = float(sys.argv[4]) if len(sys.argv) > 4 else 14.0
if os.path.exists(Z): shutil.rmtree(Z)
shutil.copytree(Q, Z, ignore=shutil.ignore_patterns("*.step", "__pycache__"))
for f in ("Markus_front_v4.step", "Markus_rueck_v4.step", "Markus_seite_v4.step"):
    pass
# --- rahmen.dxf: Bohrungen Bassdruecker neu setzen (Zeichnung: x = -X, y = Reihenhoehe) ---
d = ezdxf.readfile(os.path.join(Q, "rahmen.dxf")); m = d.modelspace()
alt = [c for c in m.query('CIRCLE[layer=="bohrungen-bassdrücker"]')]
reihen = {}
for c in alt:
    x, y = round(-c.dxf.center.x, 1), round(c.dxf.center.y, 1)
    if 0 < x < 382 and -4.24 < y < 71.75: reihen.setdefault(y, set()).add((x, c.dxf.radius))
ys = sorted(reihen)                                     # 10, 24, 38
y0 = ys[0]
neu = {}
for i, y in enumerate(ys):
    xs = sorted(reihen[y]); n = len(xs); mitte = (xs[0][0] + xs[-1][0]) / 2; r = xs[0][1]
    neu[y] = (y0 + i * RA, [mitte + (k - (n - 1) / 2) * T for k in range(n)], r)
for c in alt:
    if 0 < -c.dxf.center.x < 382 and -4.24 < c.dxf.center.y < 71.75: m.delete_entity(c)
for y, (yn, xs, r) in neu.items():
    for x in xs: m.add_circle((-x, yn), r, dxfattribs={"layer": "bohrungen-bassdrücker"})
d.saveas(os.path.join(Z, "rahmen.dxf"))
# --- blende.dxf: Bohrungen der kurzen Reihe (relativ zur ersten, x = -k*T) ---
b = ezdxf.readfile(os.path.join(Q, "blende.dxf")); bm = b.modelspace()
bc = sorted(bm.query("CIRCLE"), key=lambda c: -c.dxf.center.x)
for k, c in enumerate(bc): c.dxf.center = (-k * T, c.dxf.center.y)
b.saveas(os.path.join(Z, "blende.dxf"))
yk = max(neu)                                          # kurze Reihe
xk0 = neu[yk][1][0]
ka = open(os.path.join(Z, "kammern.py"), encoding="utf-8").read()
alt_s = 'T_BL, Z_REIHE, X_REF = 7.0, 38.0 + 4.24, 161.0'
assert ka.count(alt_s) == 1
ka = ka.replace(alt_s, f'T_BL, Z_REIHE, X_REF = 7.0, {neu[yk][0]:.1f} + 4.24, {xk0:.2f}   # Variante Teilung {T:g}')
a_ = """def _add_gespiegelt(obj, *args, **kw):
    if SPIEGELN and isinstance(obj, cq.Shape):
        obj = obj.mirror("YZ", basePointVector=(L / 2, 0, 0))
"""
assert ka.count(a_) == 1
ka = ka.replace(a_, """exec(open("streck.py", encoding="utf-8").read())      # Variante: Laenge (BASS_DL)
def _add_gespiegelt(obj, *args, **kw):
    if SPIEGELN and isinstance(obj, cq.Shape):
        obj = obj.mirror("YZ", basePointVector=(L / 2, 0, 0))
    if isinstance(obj, cq.Shape): obj = strecke(obj, kw.get("name", ""))
""")
ka = ka.replace("a.save(sys.argv[2])", """if BASS_DL > 0: print("Strecken um", BASS_DL, "mm - Querschnitt nicht gleichmaessig bei:", _streck_warn)
a.save(sys.argv[2])""", 1)
open(os.path.join(Z, "kammern.py"), "w", encoding="utf-8").write(ka)
# --- mechanik.py: Reihenhoehen aus dem Reihenabstand, Schubstange -> Knopf nach Rang in der Reihe ---
me = open(os.path.join(Z, "mechanik.py"), encoding="utf-8").read()
a1 = 'r = {10.0: "A", 24.0: "B", 38.0: "C"}[round(y, 1)]'
assert me.count(a1) == 1
me = me.replace(a1, 'r = {%.1f: "A", %.1f: "B", %.1f: "C"}[round(y, 1)]' % tuple(neu[y][0] for y in ys))
a2 = 'REIHEN = {r_: (z_, STUFE[r_][0] + SPALT) for r_, z_ in (("A", 10.0 + 4.24), ("B", 24.0 + 4.24), ("C", 38.0 + 4.24))}'
assert me.count(a2) == 1
me = me.replace(a2, 'REIHEN = {r_: (z_, STUFE[r_][0] + SPALT) for r_, z_ in (("A", %.1f + 4.24), ("B", %.1f + 4.24), ("C", %.1f + 4.24))}' % tuple(neu[y][0] for y in ys))
a3 = '''    kn = min([k for k in knoepfe if k[1] == r], key=lambda k: abs(k[0] - vorne[0]))
    stangen.append((hinten, r, kn[0]))'''
assert me.count(a3) == 1
me = me.replace(a3, '''    stangen.append((hinten, r, vorne[0]))
# Variante: Zuordnung je Reihe nach Rang (Reihenfolge in X) - funktioniert auch bei anderer Teilung
for r_ in ("A", "B", "C"):
    _ix = sorted([i_ for i_, s_ in enumerate(stangen) if s_[1] == r_], key=lambda i_: stangen[i_][2])
    _kx = sorted(k[0] for k in knoepfe if k[1] == r_)
    assert len(_ix) == len(_kx), (r_, len(_ix), len(_kx))
    for i_, x_ in zip(_ix, _kx): stangen[i_] = (stangen[i_][0], r_, x_)''')
open(os.path.join(Z, "mechanik.py"), "w", encoding="utf-8").write(me)
print("Variante Teilung", T, "Reihenabstand", RA, {chr(65 + i): (neu[y][0], [round(v, 1) for v in neu[y][1]]) for i, y in enumerate(ys)})
