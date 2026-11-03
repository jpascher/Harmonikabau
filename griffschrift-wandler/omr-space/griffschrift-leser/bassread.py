# Bassbuchstaben unter den Notenlinien lesen: Kandidaten (zusammenhängende Flecken in einem Band unter dem System), Merkmale, Zuordnung
import os as _os, cv2, numpy as np
HERE = _os.path.dirname(_os.path.abspath(__file__))

import pickle as _pk, gzip as _gz
_M = None
def model():
    global _M
    if _M is None:
        f = _os.path.join(HERE, 'bass.pkl.gz')
        _M = _pk.load(_gz.open(f, 'rb')) if _os.path.exists(f) else False
    return _M

def read_letters(bw, st, si, d):
    """Bassbuchstaben unter System si: Liste (x, Buchstabe) von links nach rechts; Zusatzzeichen (' oder ²) wird als ' angehängt"""
    M = model()
    if not M: return []
    C = candidates(bw, st, si, d)
    if not C: return []
    P = M.predict_proba(np.array([f for _, f in C])); cl = M.classes_[P.argmax(1)]; pr = P.max(1)
    let = [(b, c) for (b, _), c, q in zip(C, cl, pr) if c not in ('none', 'mark') and q > float(_os.environ.get('BASS_Q', '0.4'))]
    if len(let) < 2: return []
    yc = np.median([b[1] + b[3] / 2 for b, _ in let])
    if _os.environ.get('BASS_ALL'):   # für das Anlernen: alle buchstabengroßen Flecken der Zeile mit Merkmal
        return [(b[0] + b[2] / 2, ('?', b, f, c)) for (b, f), c in zip(C, cl) if abs(b[1] + b[3] / 2 - yc) < 0.9 * d and 0.5 * d <= b[3] <= 2.4 * d]
    row = [(b, c) for b, c in let if abs(b[1] + b[3] / 2 - yc) < 0.9 * d]   # die Buchstabenzeile
    # übereinander gesetzte Buchstaben (Grundbass und Begleiter gleichzeitig, z. B. bei Rosenzopf oder Michlbauer): ein zweiter Buchstabe
    # knapp über oder unter einem Buchstaben der Zeile, an derselben Stelle
    extra = [(b, c) for b, c in let if 0.9 * d <= abs(b[1] + b[3] / 2 - yc) <= 2.4 * d]
    marks = [b for (b, _), c in zip(C, cl) if c == 'mark']
    def name(b, c):
        x0, y0, w, h = b
        if any(x0 + w - 0.2 * d <= m[0] <= x0 + w + 0.8 * d and m[1] < y0 + 0.6 * h for m in marks): c += "'"
        return c
    out = []
    for b, c in sorted(row, key=lambda t: t[0][0]):
        cx = b[0] + b[2] / 2; n = name(b, c)
        st_ = [name(b2, c2) for b2, c2 in extra if abs(b2[0] + b2[2] / 2 - cx) < 0.6 * d]
        if st_:   # Grundbass (groß) zuerst, dann Begleiter
            n = '+'.join(sorted([n] + st_[:1], key=lambda q: q[0].islower()))
        out.append((cx, n))
    return out

def strip_long_lines(bw, d):
    """lange waagrechte Striche (Druckstriche, Hilfslinien) entfernen, damit sie nicht mit Buchstaben verschmelzen"""
    k = cv2.getStructuringElement(cv2.MORPH_RECT, (max(3, int(2.5 * d)), 1))
    lines = cv2.morphologyEx(bw, cv2.MORPH_OPEN, k)
    return cv2.subtract(bw, lines)

def band_of(st, si, H, d):
    ys = st[si]; y0 = int(ys[4] + 0.8 * d)
    nxt = st[si + 1][0] if si + 1 < len(st) else H
    y1 = int(min(ys[4] + 10 * d, nxt - 2.5 * d, H))
    return y0, y1

def feat(mask, d, ys4, box):
    x, y, w, h = box
    c = mask[y:y + h, x:x + w]
    s = max(w, h); pad = np.zeros((s, s), np.uint8); pad[(s - h) // 2:(s - h) // 2 + h, (s - w) // 2:(s - w) // 2 + w] = c
    v = cv2.resize(pad, (20, 20), interpolation=cv2.INTER_AREA).ravel().astype(np.float32) / 255.0
    return np.concatenate([v, [h / d, w / d, (y + h / 2 - ys4) / d, float((c > 0).mean())]]).astype(np.float32)

def candidates(bw, st, si, d):
    """Flecken im Band unter System si: Liste (box, Merkmal); box in Seitenpixeln (x, y, w, h)"""
    H, W = bw.shape; y0, y1 = band_of(st, si, H, d)
    if y1 - y0 < 1.5 * d: return []
    sub = strip_long_lines(bw[y0:y1], d)
    # Notenhälse, die von oben in die Buchstabenzeile ragen, abschneiden (sonst kleben sie an Buchstaben, etwa bei Rosenzopf)
    top = (sub[0] > 0) & (bw[max(0, y0 - 3):y0].min(0) > 0 if y0 >= 3 else sub[0] > 0)   # nur Striche, die schon über dem Band beginnen
    for x in (np.where(top)[0] if not _os.environ.get('NOSTEMCUT') else []):
        col = sub[:, x] > 0; L = int(np.argmin(col)) if not col.all() else len(col)
        if L >= 0.5 * d:
            sub[:L, x] = 0
    sub = cv2.morphologyEx(sub, cv2.MORPH_CLOSE, np.ones((3, 3), np.uint8))   # dünne, blasse Striche (Fotos) zusammenhalten
    n, lab, S, _ = cv2.connectedComponentsWithStats((sub > 0).astype(np.uint8))
    out = []
    for i in range(1, n):
        x, y, w, h, a = S[i]
        if not (0.15 * d <= h <= 2.6 * d and 0.1 * d <= w <= 2.6 * d): continue
        m = ((lab[y:y + h, x:x + w] == i) * 255).astype(np.uint8)
        full = np.zeros((h, w), np.uint8); full[:] = m
        box = (int(x), int(y + y0), int(w), int(h))
        out.append((box, feat_box(full, d, st[si][4], box)))
    return out

def feat_box(m, d, ys4, box):
    x, y, w, h = box
    s = max(w, h); pad = np.zeros((s, s), np.uint8); pad[(s - h) // 2:(s - h) // 2 + h, (s - w) // 2:(s - w) // 2 + w] = m
    v = cv2.resize(pad, (20, 20), interpolation=cv2.INTER_AREA).ravel().astype(np.float32) / 255.0
    return np.concatenate([v, [h / d, w / d, 0.0, float((m > 0).mean())]]).astype(np.float32)   # Höhe unter dem System zählt nicht (je nach Satz verschieden)
