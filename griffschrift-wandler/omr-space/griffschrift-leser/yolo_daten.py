# Trainingsdaten für den Kopf-Erkenner (YOLO) aus den Systembildern von yolo_erzeuge_wandlerbilder.js:
# mehrere Systeme werden zu einer Seite gestapelt, die Seite wird zufällig wie ein Scan, eine Kopie oder ein Handyfoto verschlechtert
# (Unschärfe, Rauschen, ungleiche Helligkeit, leichte Drehung und Verzerrung, dünnere oder dickere Striche, JPEG) und in Ausschnitte
# von 960 × 960 Bildpunkten geschnitten. Ausgabe im YOLO-Format (images/, labels/, data.yaml).
# Jedes System wird auf einen Linienabstand um 18 Bildpunkte gebracht, wie ihn der Leser (load_gray) für jede Seite herstellt.
# Ein Teil der Seiten wird wie die Fuchs-Blätter aus einer niedrig aufgelösten Vorlage hochgerechnet. Liegen echte Seiten mit
# Rahmen vor (yolo_fuchs_seiten.py, yolo_foto_seiten.py), wird jede vierte Seite (ECHT_JEDE) eine davon, jedesmal anders verschlechtert;
# echte Fotos ("foto": true) nur leicht, weil sie schon wie Fotos aussehen.
# Aufruf: python3 yolo_daten.py ROHORDNER AUSGABEORDNER [SEITEN] [ECHTORDNER]
import sys, os, json, glob, random
import numpy as np, cv2

KLASSEN = ['rund_voll', 'rund_hohl', 'kreuz', 'kreuz_kreis', 'dreieck_voll', 'dreieck_hohl', 'vorkreuz']
T = 960   # Ausschnittgröße


def lade(raw):
    out = []
    for f in sorted(glob.glob(os.path.join(raw, '*.json'))):
        j = json.load(open(f)); im = cv2.imread(f[:-5] + '.png', cv2.IMREAD_GRAYSCALE)
        if im is None or not j.get('boxes'): continue
        s = im.shape[1] / j['w']; ys = j.get('staff') or []
        d = (ys[4] - ys[0]) / 4 * s if len(ys) >= 5 else 8 * s
        out.append((im, [(c, x * s, y * s, w * s, h * s) for c, x, y, w, h in j['boxes']], d))
    return out


def lade_echt(ordner):
    out = []
    for f in sorted(glob.glob(os.path.join(ordner or '', '*.json'))):
        im = cv2.imread(f[:-5] + '.png', cv2.IMREAD_GRAYSCALE); j = json.load(open(f))
        if im is not None and j.get('boxes'): out.append((im, [tuple(b) for b in j['boxes']], bool(j.get('foto'))))
    return out


def echte_seite(echt, rnd):
    im, bx, foto = rnd.choice(echt); k = rnd.uniform(0.88, 1.12)
    im = cv2.resize(im, None, fx=k, fy=k, interpolation=cv2.INTER_AREA if k < 1 else cv2.INTER_LINEAR)
    return im, [(c, x * k, y * k, w * k, h * k) for c, x, y, w, h in bx], foto


def leicht(g, rnd):
    """echte Fotos sind schon unscharf und ungleich hell: nur etwas Rauschen, Helligkeit und JPEG"""
    g = g.astype(np.float32) * rnd.uniform(0.85, 1.05) + rnd.uniform(-15, 10) + np.random.normal(0, rnd.uniform(0, 4), g.shape)
    g = np.clip(g, 0, 255).astype(np.uint8)
    return cv2.imdecode(cv2.imencode('.jpg', g, [cv2.IMWRITE_JPEG_QUALITY, rnd.randint(50, 92)])[1], cv2.IMREAD_GRAYSCALE)


def grob(pg, rnd):
    """wie eine niedrig aufgelöste Vorlage (GIF), die hochgerechnet wurde: Punkte werden Kästchen, Linien grau und unscharf"""
    f = rnd.uniform(2.3, 3.2); H, W = pg.shape
    kl = cv2.resize(pg, (max(1, int(W / f)), max(1, int(H / f))), interpolation=cv2.INTER_AREA)
    if rnd.random() < 0.6: kl = np.where(kl < rnd.randint(110, 170), 0, 255).astype(np.uint8)
    return cv2.resize(kl, (W, H), interpolation=rnd.choice([cv2.INTER_LINEAR, cv2.INTER_CUBIC, cv2.INTER_NEAREST]))


def seite(sys_, rnd):
    """4–8 Systeme untereinander auf weißem Blatt, auf eine gemeinsame Breite gebracht."""
    dz = rnd.uniform(15.5, 21.0); teile = []; boxen = []; y = rnd.randint(60, 250)
    wahl = rnd.sample(sys_, k=min(len(sys_), rnd.randint(4, 8)))
    W = max(int(im.shape[1] * dz / d) for im, _, d in wahl) + 40
    for im, bx, d in wahl:
        k = dz / d * rnd.uniform(0.97, 1.03); im2 = cv2.resize(im, None, fx=k, fy=k, interpolation=cv2.INTER_AREA if k < 1 else cv2.INTER_LINEAR)
        x0 = rnd.randint(20, max(21, W - im2.shape[1] + 20))
        teile.append((im2, x0, y)); boxen += [(c, x0 + x * k, y + yy * k, w * k, h * k) for c, x, yy, w, h in bx]
        y += im2.shape[0] + rnd.randint(10, 120)
    H = y + rnd.randint(60, 250); pg = np.full((H, W + 60), 255, np.uint8)
    for im2, x0, y0 in teile:
        h, w = im2.shape; pg[y0:y0 + h, x0:x0 + w] = np.minimum(pg[y0:y0 + h, x0:x0 + w], im2[:, :min(w, pg.shape[1] - x0)])
    return pg, boxen


def verschlechtern(pg, boxen, rnd):
    art = rnd.choice(['sauber', 'scan', 'kopie', 'foto', 'foto'])
    g = pg.astype(np.float32)
    if art in ('kopie', 'foto'):   # Striche dünner oder dicker (auf Fotos auch deutlich dicker)
        k = np.ones((3, 3) if art == 'foto' and rnd.random() < 0.4 else (2, 2), np.uint8)
        g = cv2.erode(g, k) if rnd.random() < 0.55 else cv2.dilate(g, k)
    if art != 'sauber':
        g = cv2.GaussianBlur(g, (0, 0), rnd.uniform(0.3, 2.2 if art == 'foto' else 0.9))
    H, W = g.shape
    M = np.eye(3, dtype=np.float32)
    if art in ('scan', 'kopie', 'foto'):   # leichte Drehung, bei Fotos auch Verzerrung
        a = rnd.uniform(-2.5, 2.5) if art == 'foto' else rnd.uniform(-1, 1)
        R = cv2.getRotationMatrix2D((W / 2, H / 2), a, 1.0); M = np.vstack([R, [0, 0, 1]]).astype(np.float32)
        if art == 'foto':
            d = 0.03 * min(W, H); src = np.float32([[0, 0], [W, 0], [W, H], [0, H]])
            dst = src + np.float32([[rnd.uniform(-d, d), rnd.uniform(-d, d)] for _ in range(4)])
            M = cv2.getPerspectiveTransform(src, dst) @ M
        g = cv2.warpPerspective(g, M, (W, H), borderValue=255)
    if art == 'foto':   # ungleiche Helligkeit, grauer Hintergrund, Rauschen
        yy, xx = np.mgrid[0:H, 0:W].astype(np.float32)
        grad = 0.55 + 0.4 * (np.cos(xx / W * rnd.uniform(1, 4) + rnd.uniform(0, 6)) * 0.5 + 0.5) * (np.cos(yy / H * rnd.uniform(1, 4)) * 0.5 + 0.5)
        g = g * grad + rnd.uniform(-20, 10)
        g = g + np.random.normal(0, rnd.uniform(3, 12), g.shape)
        g = (g - g.min()) / max(1e-3, g.max() - g.min()) * rnd.uniform(110, 255) + rnd.uniform(0, 70)   # auch flau: graue Köpfe auf grauem Papier
    elif art in ('scan', 'kopie'):
        g = g + np.random.normal(0, rnd.uniform(1, 6), g.shape)
    g = np.clip(g, 0, 255).astype(np.uint8)
    if art != 'sauber':
        q = rnd.randint(35, 90); g = cv2.imdecode(cv2.imencode('.jpg', g, [cv2.IMWRITE_JPEG_QUALITY, q])[1], cv2.IMREAD_GRAYSCALE)
    nb = []
    for c, x, y, w, h in boxen:   # Rahmen mitbewegen (Ecken umrechnen, umschließendes Rechteck)
        pts = np.float32([[x, y], [x + w, y], [x + w, y + h], [x, y + h]]).reshape(-1, 1, 2)
        p2 = cv2.perspectiveTransform(pts, M).reshape(-1, 2)
        x0, y0 = p2.min(0); x1, y1 = p2.max(0); nb.append((c, x0, y0, x1 - x0, y1 - y0))
    return g, nb


def ausschnitte(g, boxen, rnd, n):
    H, W = g.shape; out = []
    for _ in range(n):
        x0 = rnd.randint(0, max(0, W - T)); y0 = rnd.randint(0, max(0, H - T))
        crop = np.full((T, T), 255, np.uint8); c_ = g[y0:y0 + T, x0:x0 + T]; crop[:c_.shape[0], :c_.shape[1]] = c_
        lab = []
        for c, x, y, w, h in boxen:
            cx, cy = x + w / 2 - x0, y + h / 2 - y0
            if 0 <= cx < T and 0 <= cy < T and x - x0 > -0.3 * w and x + w - x0 < T + 0.3 * w:
                lab.append(f"{c} {cx / T:.6f} {cy / T:.6f} {min(w, T) / T:.6f} {min(h, T) / T:.6f}")
        if lab: out.append((crop, lab))
    return out


if __name__ == '__main__':
    raw, ziel = sys.argv[1], sys.argv[2]; n_seiten = int(sys.argv[3]) if len(sys.argv) > 3 else 600
    rnd = random.Random(7); np.random.seed(7)
    sys_ = lade(raw); echt = lade_echt(sys.argv[4] if len(sys.argv) > 4 else None); print('Systeme', len(sys_), 'echte Seiten', len(echt))
    for teil in ('train', 'val'):
        os.makedirs(os.path.join(ziel, 'images', teil), exist_ok=True); os.makedirs(os.path.join(ziel, 'labels', teil), exist_ok=True)
    k = 0; jede = int(os.environ.get('ECHT_JEDE', '4'))   # jede n-te Seite ist eine echte (dritter Lernlauf: 2)
    for i in range(n_seiten):
        foto = False
        if echt and i % jede == jede - 1: pg, bx, foto = echte_seite(echt, rnd)
        else:
            pg, bx = seite(sys_, rnd)
            if rnd.random() < 0.2: pg = grob(pg, rnd)
        if foto: g = leicht(pg, rnd)
        else: g, bx = verschlechtern(pg, bx, rnd)
        teil = 'val' if i % 10 == 0 else 'train'
        for crop, lab in ausschnitte(g, bx, rnd, 3):
            fn = f"s{i:05d}_{k:06d}"; cv2.imwrite(os.path.join(ziel, 'images', teil, fn + '.jpg'), crop, [cv2.IMWRITE_JPEG_QUALITY, 92])
            open(os.path.join(ziel, 'labels', teil, fn + '.txt'), 'w').write('\n'.join(lab) + '\n'); k += 1
    with open(os.path.join(ziel, 'data.yaml'), 'w') as f:
        f.write(f"path: {os.path.abspath(ziel)}\ntrain: images/train\nval: images/val\nnames:\n" + ''.join(f"  {i}: {n}\n" for i, n in enumerate(KLASSEN)))
    print('Ausschnitte', k)
