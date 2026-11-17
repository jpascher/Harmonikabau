# Animationen Diskantmechanik Markus: (1) mit Griffstock, (2) ohne Griffstock, (3) Explosion
import pickle, sys, os, math, numpy as np, pyvista as pv, subprocess, shutil
pv.global_theme.background = "white"
D = pickle.load(open("netze.pkl", "rb"))
T, F, DP = D["teile"], D["farben"], D["drehpunkte"]
MODUS = sys.argv[1]
W, H = 1280, 720
def poly(v, t):
    return pv.PolyData(v, np.hstack([np.full((len(t), 1), 3), t]).astype(np.int64))
pl = pv.Plotter(off_screen=True, window_size=(W, H))
pl.enable_anti_aliasing("ssaa")
akt = {}
for n, (f, v, t) in T.items():
    if MODUS == "ohne" and n.startswith(("Griffstock_Koerper", "Griffstock_Teil2", "Griffstock_Abdeckung")): continue
    if os.environ.get("OHNE_VD") and n.startswith("Verdeck"): continue
    op = 0.35 if (MODUS == "ohne" and n.startswith(("Kunststoff", "Dreieck"))) else 1.0
    akt[n] = pl.add_mesh(poly(v, t), color=F[f], opacity=op, smooth_shading=False)
# Hebelgruppen
grp = {}
for n in akt:
    if "/" in n: grp.setdefault(n.split("/")[0], []).append(n)
def mitte(n):
    v = T[n][1]; return v.mean(0)
WINKEL = {}
for g, ns in grp.items():
    x, Yp, Zp = DP[g]
    kn = [n for n in ns if n.endswith("/Knopf")]
    k = mitte(kn[0]) if kn else np.array([x, Yp, Zp - 80])
    # Druck so weit, dass der Knopf 0,3 mm ueber seiner Stufe stehen bleibt (Knopf in Ruhe 3,5 ueber der Stufe)
    KD = (0.1908, 0.9816)
    vk = T[kn[0]][1].astype(float) if kn else None
    def amin(t):
        c_, s_ = math.cos(t), math.sin(t); dY = vk[:, 1] - Yp; dZ = vk[:, 2] - Zp
        Y = Yp + c_ * dY + s_ * dZ; Z = Zp - s_ * dY + c_ * dZ
        return (-Z * KD[0] + Y * KD[1]).min()
    if vk is not None:
        ziel = amin(0.0) - 3.2; lo, hi = 0.0, 0.2
        for _ in range(40):
            m = (lo + hi) / 2
            if amin(m) > ziel: lo = m
            else: hi = m
        WINKEL[g] = lo
    else:
        WINKEL[g] = 3.5 / math.hypot(k[1] - Yp, k[2] - Zp)
def drehmatrix(g, t):
    x, Yp, Zp = DP[g]; c, s = math.cos(t), math.sin(t)
    M = np.eye(4)
    M[1, 1], M[1, 2], M[2, 1], M[2, 2] = c, s, -s, c
    M[1, 3] = Yp - (c * Yp + s * Zp); M[2, 3] = Zp - (-s * Yp + c * Zp)
    return M
reihe = lambda g: int(g.split("_")[1][1])
reihenfolge = sorted(grp, key=lambda g: (reihe(g), DP[g][0]))
def kamera(modus):
    if modus == "mit":
        pl.camera_position = [(-120, 330, -260), (190, 60, -30), (0, 0, -1)]
    elif modus == "ohne":
        pl.camera_position = [(-90, -260, -230), (190, 70, -30), (0, 0, -1)]
    else:
        pl.camera_position = [(-420, 620, -560), (190, 40, -60), (0, 0, -1)]
kamera(MODUS)
os.makedirs(f"fr_{MODUS}{os.environ.get('OHNE_VD','')}", exist_ok=True)
fr = 0
def bild():
    global fr
    pl.render(); pl.screenshot(f"fr_{MODUS}{os.environ.get('OHNE_VD','')}/{fr:05d}.png"); fr += 1
if MODUS in ("mit", "ohne"):
    DAUER, VERSATZ = 30, 10
    N = VERSATZ * len(reihenfolge) + DAUER + 12
    for i in range(N):
        for j, g in enumerate(reihenfolge):
            p = (i - j * VERSATZ) / DAUER
            h = math.sin(math.pi * p) if 0 <= p <= 1 else 0.0
            M = drehmatrix(g, WINKEL[g] * h)
            for n in grp[g]: akt[n].user_matrix = M
        bild()
    # zum Schluss alle Tasten gleichzeitig
    for i in range(72):
        h = math.sin(math.pi * i / 71)
        for g in reihenfolge:
            M = drehmatrix(g, WINKEL[g] * h)
            for n in grp[g]: akt[n].user_matrix = M
        bild()
else:
    def ziel(n):
        v = np.zeros(3)
        if n.startswith("Diskantboden"): v = [0, 0, 60]
        elif n.startswith("Verstaerkung"): v = [0, 0, 95]
        elif n.startswith("Schraube_Boden"): v = [0, 0, 135]
        elif n.startswith("Griffstock_Koerper"): v = [0, -40, -150]
        elif n.startswith("Griffstock_Abdeckung"): v = [0, -110, -60]
        elif n.startswith("Kunststoff"): v = [0, -60, -35]
        elif n.startswith("Schraube_Kunststoff_Dreieck"): v = [0, -45, -60]
        elif n.startswith("Schraube_Kunststoff_Boden"): v = [0, -60, -65]
        elif n.startswith("Dreieck"): v = [0, -85, 10]
        elif "/" in n:
            r = reihe(n.split("/")[0]); t = n.split("/")[1]
            basis = np.array([0, 0, -40 - 22 * (5 - r)])
            if t in ("Klappe",): v = [0, 0, -12]
            elif t == "Klappenpolster": v = [0, 0, -4]
            elif t == "Knopf": v = basis + [0, 0, -35]
            elif t == "Achse": v = basis + [0, -25, 0]
            elif t == "Holzteil": v = basis + [0, 0, -15]
            else: v = basis
        return np.array(v, dtype=float)
    Z = {n: ziel(n) for n in akt}
    seq = [0.0] * 24 + [0.5 - 0.5 * math.cos(math.pi * i / 143) for i in range(144)] + [1.0] * 60
    seq += list(reversed(seq[24:168])) + [0.0] * 24
    for k, s in enumerate(seq):
        for n, a in akt.items():
            M = np.eye(4); M[:3, 3] = Z[n] * s; a.user_matrix = M
        pl.camera.azimuth = 0
        bild()
pl.close()
ff = shutil.which("ffmpeg")
subprocess.run([ff, "-y", "-loglevel", "error", "-framerate", "24", "-i", f"fr_{MODUS}{os.environ.get('OHNE_VD','')}/%05d.png", "-c:v", "libx264",
                "-preset", "veryfast", "-crf", "23", "-pix_fmt", "yuv420p", sys.argv[2]], check=True)
print("fertig", sys.argv[2], fr, "Bilder")
