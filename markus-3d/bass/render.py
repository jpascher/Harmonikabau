import pickle, math, sys, numpy as np, pyvista as pv
pv.OFF_SCREEN = True
alle, d = pickle.load(open("anim/daten.pkl", "rb"))
L = d["L"]; ZA = d["Z_ACHSE"]; alpha = d["alpha"]; HUB = d["HUB"]
TEST = len(sys.argv) > 1 and sys.argv[1] == "test"
def faces(t): return np.hstack([np.full((len(t), 1), 3), t]).ravel()
def farbe(n):
    if n.startswith(("Knopf_",)): return (0.08, 0.08, 0.08)
    if n.startswith("Knopfzyl"): return (0.95, 0.95, 0.9)
    if n.startswith("Filz") or n.endswith("_Filz"): return (0.9, 0.88, 0.8)
    if n.endswith("_Leder"): return (0.45, 0.28, 0.15)
    if n.endswith("_Epoxy"): return (0.15, 0.15, 0.18)
    if n.startswith("Gabelkopf"): return (0.2, 0.45, 0.85)
    if n.startswith("Gelenkbolzen"): return (0.85, 0.2, 0.2)
    if n.startswith(("Lagerbock", "Luftklappe_Lagerbock")) and "Schraube" not in n: return (0.92, 0.9, 0.82)
    if n in ("Bassboden", "Front", "Front_Aufdoppelung_2mm", "Front_Aufdoppelung_5mm", "Zwischen_Deckel"): return (0.86, 0.72, 0.5)
    if n.startswith("Lufttaste"): return (0.35, 0.25, 0.15)
    return (0.7, 0.7, 0.74)
pl = pv.Plotter(off_screen=True, window_size=(1280, 720))
pl.set_background((0.97, 0.97, 0.97))
basis, akt = {}, {}
for n, (v, t) in alle.items():
    if n.startswith(("Achse_", "Luftklappe_Achse")) is False and n == "Zwischen_Deckel": pass
    v = v.copy(); v[:, 0] = L - v[:, 0]                      # Modell ist gespiegelt
    m = pv.PolyData(v, faces(t))
    basis[n] = v; akt[n] = m
    op = 0.35 if n in ("Bassboden", "Zwischen_Deckel") else 1.0
    pl.add_mesh(m, color=farbe(n), opacity=op, smooth_shading=False)
def rot_mat(axis, ang):
    axis = axis / np.linalg.norm(axis); a = math.cos(ang / 2); b, c, e = -axis * math.sin(ang / 2)
    return np.array([[a*a+b*b-c*c-e*e, 2*(b*c+a*e), 2*(b*e-a*c)], [2*(b*c-a*e), a*a+c*c-b*b-e*e, 2*(c*e+a*b)], [2*(b*e+a*c), 2*(c*e-a*b), a*a+e*e-b*b-c*c]])
def drehe(pts, p0, axis, ang):
    R = rot_mat(axis, ang); return (pts - p0) @ R.T + p0
# Achsen gespiegelt
gr = []
for ai, drehen, schieben, knoepfe in d["gruppen"]:
    a0, a1 = d["achsen"][ai]
    p0 = np.array([L - a0[0], a0[1], ZA]); ax = np.array([-(a1[0] - a0[0]), a1[1] - a0[1], 0.0])
    test = drehe(np.array([[p0[0], p0[1], ZA - 12.7]]), p0, ax, 0.2)[0]
    sg = 1.0 if test[1] < p0[1] else -1.0
    gr.append((p0, ax, sg, [n for n in drehen if n in basis], [n for n in schieben if n in basis], knoepfe))
lk = d["lk"]; lk_p0 = np.array([0.0, lk["Y_AX_LK"], lk["Z_AX_LK"]]); lk_ax = np.array([-1.0, 0, 0])   # gespiegelt: X umgedreht
lk_dreh = [n for n in basis if n.startswith("Luftklappe_") and not n.startswith(("Luftklappe_Lagerbock", "Luftklappe_Druckstange"))]
lk_schieb = [n for n in basis if n in ("Lufttaste", "Luftklappe_Druckstange")]
def lk_test():  # Vorzeichen: gedrueckt (f=1) = Modell; geschlossen = um -winkel zurueck
    return 1.0
def zustand(gruppe_f, lk_f):
    for n in akt: akt[n].points = basis[n]
    for (p0, ax, sg, dr, sc, _), f in zip(gr, gruppe_f):
        if f <= 0: continue
        for n in dr: akt[n].points = drehe(basis[n], p0, ax, sg * alpha * f)
        for n in sc: akt[n].points = basis[n] + np.array([0, -HUB * f, 0])
    w = math.radians(lk["winkel"]) * (1 - lk_f)
    # Klappe zu bei lk_f=0: Rueckdrehung um die Achse; Richtung so, dass die Klappe sinkt
    for n in lk_dreh:
        b = basis[n]; c1 = drehe(b, lk_p0, lk_ax, w); c2 = drehe(b, lk_p0, lk_ax, -w)
        akt[n].points = c1 if c1[:, 2].mean() < c2[:, 2].mean() else c2
    for n in lk_schieb: akt[n].points = basis[n] + np.array([0, 0, -lk["HUB_LT"] * (1 - lk_f)])
pl.camera_position = [(L / 2, 420, -260), (L / 2, 100, 25), (0, -1, 0)]
pl.reset_camera(); pl.camera.zoom(1.35)
pl.add_text("Bassmechanik Markus", font_size=12, color="black")
lab = pl.add_text("", position="lower_left", font_size=11, color="black")
# Reihenfolge: Reihe 1, 2, 3 von links (Modellansicht)
folge = []
for r in "ABC":
    ks = sorted({(L - X, gi) for gi, g in enumerate(gr) for (rr, X) in g[5] if rr == r})
    folge += [(r, x, gi) for x, gi in ks]
kurve = [i / 6 for i in range(1, 7)] + [1.0] * 3 + [1 - i / 6 for i in range(1, 7)]
frames = 0
if TEST:
    gf = [0] * len(gr); gf[folge[3][2]] = 1.0; zustand(gf, 1.0); pl.screenshot("anim/test.png"); sys.exit(0)
pl.open_movie("anim/Bassmechanik_Markus.mp4", framerate=20, quality=8)
zustand([0] * len(gr), 0.0)
for _ in range(15): pl.write_frame()
for r, x, gi in folge:
    lab.SetText(0, f"Reihe {'ABC'.index(r) + 1}, Knopf bei {x:.0f} mm")
    for f in kurve:
        gf = [0] * len(gr); gf[gi] = f; zustand(gf, 0.0); pl.write_frame(); frames += 1
lab.SetText(0, "Lufttaste -> Luftklappe")
for f in [i / 10 for i in range(1, 11)] + [1.0] * 8 + [1 - i / 10 for i in range(1, 11)]:
    zustand([0] * len(gr), f); pl.write_frame(); frames += 1
pl.close(); print("Frames", frames)
