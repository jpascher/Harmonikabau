# Erster Versuch eines Griffschrift-Lesers: Notensystem, Notenköpfe (rund/Kreuz, voll/hohl), Lage, Bässe mit Druckstrich
import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2, numpy as np, json, sys

def staves(bw):
    """5-Linien-Systeme über die waagrechte Projektion finden; liefert Liste von 5 y-Werten."""
    h, w = bw.shape
    # nur lange waagrechte Striche zählen (Balken und Köpfe fallen weg, sonst verschmelzen sie mit einer Linie)
    lin = cv2.morphologyEx(bw, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_RECT, (max(20, int(w * 0.08)), 1)))
    prof = lin.sum(axis=1) / 255
    rows = [y for y in range(h) if prof[y] > 0.3 * prof.max() and prof[y] > 0.15 * w]   # 0,3: auch kurze Systeme (z. B. vor „Fine“); fünf gleich weite Linien filtern den Rest
    # benachbarte Zeilen zusammenfassen (Liniendicke)
    lines = []
    for y in rows:
        if lines and y - lines[-1][-1] <= 1: lines[-1].append(y)
        else: lines.append([y])
    ys = [sum(l) / len(l) for l in lines]
    # alle Fünfergruppen bewerten (Schwankung der Abstände), die gleichmäßigsten zuerst nehmen, ohne Überlappung
    cand = []
    for i in range(len(ys) - 4):
        g = ys[i:i + 5]; d = np.diff(g)
        if d.mean() > 4 and d.max() - d.min() <= max(2, 0.25 * d.mean()): cand.append(((d.max() - d.min()) / d.mean(), i))
    used = set(); out = []
    for _, i in sorted(cand):
        if any(j in used for j in range(i, i + 5)): continue
        used.update(range(i, i + 5)); g = ys[i:i + 5]
        # gleichmäßig ausgleichen (eine Linie kann durch einen anliegenden Balken leicht verschoben sein)
        k = np.arange(5); a1, a0 = np.polyfit(k, g, 1); out.append([float(a0 + a1 * q) for q in k])
    out.sort(key=lambda g: g[0])
    return out

def remove_lines(bw, ys, d):
    """Notenlinien (und Hilfslinien) entfernen: waagrechte dünne Strukturen löschen, die nicht von Köpfen bedeckt sind."""
    k = cv2.getStructuringElement(cv2.MORPH_RECT, (int(d * 2.5), 1))
    horiz = cv2.morphologyEx(bw, cv2.MORPH_OPEN, k)
    thin = cv2.subtract(horiz, cv2.morphologyEx(horiz, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_RECT, (1, max(3, int(d * 0.35))))))
    return cv2.subtract(bw, thin)

def stems_off(img, d):
    """Senkrechte Hälse entfernen (lange dünne senkrechte Striche)."""
    k = cv2.getStructuringElement(cv2.MORPH_RECT, (1, int(d * 2.2)))
    vert = cv2.morphologyEx(img, cv2.MORPH_OPEN, k)
    return cv2.subtract(img, vert)

def templates(d):
    """Kopfvorlagen in der Größe des Linienabstands d (Wandler-Stil): rund voll/hohl, Kreuz, Kreuz im Kreis."""
    s = int(d * 1.6) | 1; c = s // 2; T = {}
    def blank(): return np.zeros((s, s), np.uint8)
    a = blank(); cv2.ellipse(a, (c, c), (int(d * 0.68), int(d * 0.49)), -20, 0, 360, 255, -1); T['o_f'] = a
    a = blank(); cv2.ellipse(a, (c, c), (int(d * 0.68), int(d * 0.49)), -20, 0, 360, 255, max(1, int(d * 0.19))); T['o_h'] = a
    a = blank(); r = int(d * 0.56); t = max(2, int(d * 0.27))
    cv2.line(a, (c - r, c - r), (c + r, c + r), 255, t); cv2.line(a, (c - r, c + r), (c + r, c - r), 255, t); T['x_f'] = a
    a = blank(); t2 = max(1, int(d * 0.16)); cv2.line(a, (c - r, c - r), (c + r, c + r), 255, t2); cv2.line(a, (c - r, c + r), (c + r, c - r), 255, t2)
    cv2.circle(a, (c, c), int(d * 0.7), 255, max(1, int(d * 0.12))); T['x_h'] = a
    # Dreieck mit der Spitze nach oben (4. Reihe bei Elke Margetich)
    w, h = int(d * 1.3), int(d * 1.05); pts = np.array([[c - w // 2, c + h // 2], [c + w // 2, c + h // 2], [c, c - h // 2]], np.int32)
    a = blank(); cv2.fillPoly(a, [pts], 255); T['t_f'] = a
    a = blank(); cv2.polylines(a, [pts], True, 255, max(1, int(d * 0.15))); T['t_h'] = a
    return T

def heads(img, d, top, bottom):
    T = templates(d); found = []
    for name, t in T.items():
        r = cv2.matchTemplate(img, t, cv2.TM_CCOEFF_NORMED)
        if name.startswith('t'): continue   # Dreiecke nur über detect.py
        thr = {'o_f': 0.5, 'o_h': 0.42, 'x_f': 0.3, 'x_h': 0.33}[name]
        ys, xs = np.where(r >= thr)
        for y, x in zip(ys, xs):
            cy = y + t.shape[0] // 2
            if top - 6 * d <= cy <= bottom + 6 * d:
                found.append((float(r[y, x]), x + t.shape[1] / 2, cy, name))
    found.sort(reverse=True); keep = []
    for sc, x, y, n in found:   # Nicht-Maximum-Unterdrückung
        if all(abs(x - k[1]) > d * 0.8 or abs(y - k[2]) > d * 0.4 for k in keep): keep.append((sc, x, y, n))
    return keep

def read(path, scale=2):
    g = cv2.imread(path, cv2.IMREAD_GRAYSCALE)
    _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    st = staves(bw)
    res = []
    for ys in st:
        d = float(np.mean(np.diff(ys)))
        nol = stems_off(remove_lines(bw, ys, d), d)
        hs = heads(nol, d, ys[0], ys[-1])
        mid = ys[2]
        for sc, x, y, n in hs:
            p = round((mid - y) / (d / 2))   # Halbstufen über der Mittellinie
            res.append({'x': x / scale, 'y': y / scale, 'p': p, 't': n[0], 'fill': n.endswith('f'), 'sc': round(sc, 2)})
    return st, res

if __name__ == '__main__':
    st, res = read(sys.argv[1]); print(len(st), 'Systeme', len(res), 'Köpfe')
