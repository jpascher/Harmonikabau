# Erster Versuch eines Griffschrift-Lesers: Notensystem, Notenköpfe (rund/Kreuz, voll/hohl), Lage, Bässe mit Druckstrich
import cv2, numpy as np, json, sys

def staves(bw):
    """5-Linien-Systeme über die waagrechte Projektion finden; liefert Liste von 5 y-Werten."""
    h, w = bw.shape
    prof = bw.sum(axis=1) / 255
    rows = [y for y in range(h) if prof[y] > 0.6 * prof.max() and prof[y] > 0.15 * w]
    # benachbarte Zeilen zusammenfassen (Liniendicke)
    lines = []
    for y in rows:
        if lines and y - lines[-1][-1] <= 1: lines[-1].append(y)
        else: lines.append([y])
    ys = [sum(l) / len(l) for l in lines]
    out, i = [], 0
    while i + 4 < len(ys):
        g = ys[i:i + 5]; d = np.diff(g)
        if d.max() - d.min() <= 2 and d.mean() > 4: out.append(g); i += 5
        else: i += 1
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
    return T

def heads(img, d, top, bottom):
    T = templates(d); found = []
    for name, t in T.items():
        r = cv2.matchTemplate(img, t, cv2.TM_CCOEFF_NORMED)
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
