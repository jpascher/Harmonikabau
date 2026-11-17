import pickle, math, numpy as np, pyvista as pv
pv.OFF_SCREEN = True
alle, L = pickle.load(open("anim/explo.pkl", "rb"))
def off(n):
    s = lambda *p: n.startswith(p)
    if s("Knopf_"): return (0, 120, 0)
    if s("Knopfzylinder", "Filzring"): return (0, 95, 0)
    if s("Frontblende"): return (0, 80, 0)
    if s("Front_Aufdoppelung"): return (0, 35, 0)
    if s("Front", "Horn_") and not s("Frontleiste"): return (0, 60, 0)
    if s("Rueckteil"): return (0, -60, 0)
    if s("Deckelauflage_hinten"): return (0, -30, -15)
    if s("Seitenteil_links"): return (-60, 0, 0)
    if s("Seitenteil_rechts"): return (60, 0, 0)
    if s("Deckelauflage"): return (0, 0, -15)
    if s("Fuss"): return (0, 0, -150)
    if s("Bassbodenleder"): return (0, 0, -125)
    if s("Bassdeckel", "Lufttaste"): return (0, 0, -110)
    if s("Eckblock"): return (0, 0, -90)
    if s("Klappe", "Begleiterklappe"): return (0, 0, -30)
    if s("Achse_", "Lagerbock", "Arm_", "Gabel", "Gelenkbolzen", "Schubstange", "Klappenhebel", "Feder_Achse"): return (0, 0, -60)
    if s("Luftklappe"): return (0, 0, 110)
    if s("Stimmplatte", "Haken", "Schraube_P"): return (0, 0, 80)
    if s("Begleiter_Stimmstock", "Begleiter_Sole", "Begleiter_Schraube", "Begleiter_Abschluss", "Begleiter_Aussenwand", "Begleiter_Mittelsteg", "Begleiter_Deckleiste"): return (0, 0, 90)
    if s("Zwischen_Deckel"): return (0, 0, 95)
    if s("Bassboden"): return (0, 0, 0)
    return (0, 0, 40)            # Kammerteile, Leisten, Folien, Begleiter-Zwischenteil
def farbe(n):
    s = lambda *p: n.startswith(p)
    if s("Knopf_", "Fuss_Gummi", "Horn_", "Bassbodenleder") or n.endswith("_Epoxy"): return (0.1, 0.1, 0.1)
    if n.endswith("_Leder"): return (0.45, 0.28, 0.15)
    if n.endswith("_Filz") or s("Filz", "Lufttaste_Filz"): return (0.9, 0.88, 0.8)
    if s("Folie"): return (0.85, 0.3, 0.3)
    if "Folie" in n: return (0.85, 0.3, 0.3)
    if s("Stimmplatte"): return (0.8, 0.82, 0.86)
    if s("Gabel"): return (0.2, 0.45, 0.85)
    if s("Gelenkbolzen"): return (0.85, 0.2, 0.2)
    if s("Achse_", "Arm_", "Schubstange", "Klappenhebel", "Feder", "Haken", "Schraube", "Lagerbock_Schraube", "Fuss_Schraube") or "Schraube" in n: return (0.6, 0.6, 0.64)
    if s("Knopfzyl", "Lagerbock"): return (0.95, 0.94, 0.88)
    if "Trennwand" in n: return (0.7, 0.78, 0.9)
    if s("Bassdeckel", "Lufttaste"): return (0.35, 0.25, 0.15)
    if s("Begleiter"): return (0.62, 0.72, 0.55)
    return (0.84, 0.68, 0.45)
pl = pv.Plotter(off_screen=True, window_size=(1280, 900)); pl.set_background((0.97, 0.97, 0.97))
basis, akt, verschub = {}, {}, {}
for n, (v, t) in alle.items():
    v = v.copy(); v[:, 0] = L - v[:, 0]
    o = np.array(off(n), float); o[0] = -o[0]
    m = pv.PolyData(v, np.hstack([np.full((len(t), 1), 3), t]).ravel())
    basis[n], akt[n], verschub[n] = v, m, o
    pl.add_mesh(m, color=farbe(n))
def setze(f):
    e = 0.5 - 0.5 * math.cos(math.pi * f)
    for n in akt: akt[n].points = basis[n] + verschub[n] * e
setze(1.0)
pl.camera_position = [(L / 2 + 420, -520, 380), (L / 2, 110, 10), (0, 0, 1)]
pl.reset_camera(); pl.camera.zoom(1.25)
pl.add_text("Bassteil Markus - Explosionsdarstellung", font_size=12, color="black")
pl.screenshot("anim/Explosion_Bassteil_Markus.png")
pl.camera_position = [(L / 2 + 420, 650, -380), (L / 2, 110, 10), (0, 0, 1)]
pl.reset_camera(); pl.camera.zoom(1.25); pl.screenshot("anim/Explosion_Bassteil_Markus_unten.png")
# Film: zusammen -> auseinander (mit Kamerafahrt) -> halten -> zusammen
pl.open_movie("anim/Explosion_Bassteil_Markus.mp4", framerate=20, quality=8)
N = 60
for i in range(N * 3 + 40):
    t = i / (N * 3 + 40)
    if i < N: f = i / N
    elif i < N + 40: f = 1.0
    elif i < 2 * N + 40: f = 1.0
    elif i < 3 * N + 40 - 0: f = 1 - (i - 2 * N - 40) / N
    ang = 2 * math.pi * 0.75 * t
    r = 680
    pl.camera_position = [(L / 2 + r * math.sin(ang + 0.6), 110 - r * math.cos(ang + 0.6), 380 * math.cos(ang * 0.5)), (L / 2, 110, 10), (0, 0, 1)]
    setze(f); pl.write_frame()
pl.close(); print("ok")
