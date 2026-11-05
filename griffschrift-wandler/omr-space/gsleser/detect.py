import os as _os; HERE=_os.path.dirname(_os.path.abspath(__file__))
import cv2, numpy as np, json, sys, pickle
sys.path.insert(0, HERE); import gsread as G; from train import patch, CL
import os as _os0
INK = float(_os0.environ.get('GS_INK', '0.12'))
MOD = HERE if _os.path.exists(HERE+'/heads.pkl') else _os.path.join(HERE, '..', 'gsleser')   # Forschungsordner: Modelle aus gsleser
clf = pickle.load(open(MOD+'/heads.pkl', 'rb'))
# Dreiecke (4. Reihe in der Notation von Elke Margetich, z. B. Philipp Ramoner): eigener kleiner Klassifikator auf denselben Ausschnitten
_TP = _os.environ.get('GS_TRI', MOD+'/tri.pkl')
TRI = pickle.load(open(_TP, 'rb')) if _os.path.exists(_TP) else None
TRI_THR, TRI_HEAD, TRI_SURE = 0.5, 0.3, 0.95   # Dreieck, wenn dessen Klassifikator ≥ 0,5 meldet und der Kopf-Klassifikator ≥ 0,3 (oder Dreieck ≥ 0,95)
CLX = CL + ['t_f', 't_h']
TRI_MIN = int(_os.environ.get('GS_TRIMIN', '3')); TRI_FRAC = float(_os.environ.get('GS_TRIFRAC', '0.6'))   # Seite mit Dreieck-Köpfen?
TRI_NB = float(_os.environ.get('GS_TRINB', '0.7'))   # Dreieck neben anderem Kopf eine Lage daneben: näher als so viele Linienabstände = derselbe Kopf
BEAM = _os.environ.get('GS_BEAM', '1') == '1'; BEAM_L = float(_os.environ.get('GS_BEAML', '1.8')); BEAM_F = float(_os.environ.get('GS_BEAMF', '0.5'))
import gzip as _gz
# zweites Kopfmodell für Fotos (angelernt an Wandler-Bildern, die wie Handyfotos verschlechtert wurden: blasse, dünne, unscharfe Striche);
# wird nur bei Fotos zugemischt, saubere Seiten lesen wie bisher
_H2 = _os.environ.get('GS_H2', MOD + '/heads2.pkl.gz')
H2 = pickle.load(_gz.open(_H2, 'rb')) if _os.path.exists(_H2) and not _os.environ.get('NOH2') else None
H2_W = float(_os.environ.get('GS_H2W', '0.5'))   # Mischgewicht; an einer nicht angelernten Rosenzopf-Seite am besten
# drittes Kopfmodell für saubere Scans fremder Schreibweisen (angelernt an Vektor-PDFs mit genauen Köpfen aus der Textebene,
# MC-Score und Finale, auch verschlechtert wie Scans); wird nur bei sauberen Seiten zugemischt (GS_H3W)
_H3 = _os.environ.get('GS_H3', MOD + '/heads3.pkl.gz')
H3 = pickle.load(_gz.open(_H3, 'rb')) if _os.path.exists(_H3) and not _os.environ.get('NOH3') else None
H3_W = float(_os.environ.get('GS_H3W', '0.5'))
TRI_BOOST = float(_os.environ.get('GS_TRIB', '0.25'))   # Vorrang der Dreieck-Meldungen beim Ausdünnen

# Kandidaten-Vorauswahl mit dem Zerlegungsnetz von homr (Versuch, GS_HOMR=1): homr findet Notenlinien, Hälse und runde Köpfe
# sehr verlässlich, Kreuze meldet es meist als Vorzeichen. Kopfkandidaten werden deshalb nur dort geprüft, wo das Netz einen Kopf
# oder ein Vorzeichen sieht (etwas aufgeweitet); welche Kopfform es ist, entscheidet weiter der eigene Klassifikator.
HOMR = _os.environ.get('GS_HOMR', '0') == '1'
HOMR_DIL = float(_os.environ.get('GS_HOMRDIL', '0.35'))   # Aufweitung der Masken in Linienabständen
HOMR_MODE = _os.environ.get('GS_HOMRMODE', 'round')
HOMR_ADD = float(_os.environ.get('GS_HOMRADD', '0'))   # > 0: runder Kopf auch ohne eigenen Treffer, wenn das Netz die Kopffläche mindestens so weit abdeckt (fremde Schriften)
_HCACHE = {}
def homr_masks(g):
    import zlib; key = (g.shape, zlib.crc32(g[::7, ::7].tobytes()))
    if key in _HCACHE: return _HCACHE[key]
    _dc = _os.environ.get('GS_HOMRCACHE')   # Ordner für zwischengespeicherte Masken (nur für Prüfläufe)
    _fn = _os.path.join(_dc, '%d_%d_%x.npz' % (g.shape[0], g.shape[1], key[1] & 0xffffffff)) if _dc else None
    if _fn and _os.path.exists(_fn):
        res = dict(np.load(_fn)); _HCACHE.clear(); _HCACHE[key] = res; return res
    from homr.resize import resize_image
    from homr import color_adjust
    from homr.segmentation.inference_segnet import extract
    from homr.main import download_weights
    download_weights(False, False, False)
    bgr = cv2.cvtColor(g, cv2.COLOR_GRAY2BGR); img = resize_image(bgr); pre = color_adjust.apply_clahe(img)
    r = extract(pre, '', step_size=320, use_cache=False, use_gpu_inference=False)
    H, W = g.shape
    up = lambda m: cv2.resize(m.astype(np.uint8), (W, H), interpolation=cv2.INTER_NEAREST)
    res = dict(notehead=up(r.notehead), clefs_keys=up(r.clefs_keys), stems=up(r.stems_rests), staff=up(r.staff))
    if _fn: _os.makedirs(_dc, exist_ok=True); np.savez_compressed(_fn, **res)
    _HCACHE.clear(); _HCACHE[key] = res
    return res

LAST = {}
def detect(path, scale=2, thr=0.5):
    g = cv2.imread(path, 0); _, bw = cv2.threshold(g, 200, 255, cv2.THRESH_BINARY_INV)
    out = []
    HM = homr_masks(g) if HOMR else None
    # Linien suchen mit hellerer Schwelle: gedruckte Notenlinien sind oft grau
    _, bwl = cv2.threshold(g, 240, 255, cv2.THRESH_BINARY_INV)
    staves_res = []
    for si, ys in enumerate(G.staves(bwl)):
        d = float(np.mean(np.diff(ys))); mid = ys[2]
        nol = G.clean_gray(g, d, int(_os.environ.get('GS_CLT', '60'))) if (G.FOTO_LINES and _os.environ.get('GS_CLEAN')) else G.stems_off(G.remove_lines(bw, ys, d), d)
        # Kandidaten: alle Lagen und Spalten, an denen nach Entfernen von Linien und Hälsen genug Tinte in Kopfgröße liegt
        ink = cv2.boxFilter((nol > 0).astype(np.float32), -1, (max(3, int(d * 0.9)), max(3, int(d * 0.6))))
        gate = None; cov = None
        if HM is not None:
            k_ = max(3, int(d * HOMR_DIL) * 2 + 1); se_ = cv2.getStructuringElement(cv2.MORPH_ELLIPSE, (k_, k_))
            if HOMR_MODE == 'all':   # alle Kandidaten nur an Netz-Köpfen oder -Vorzeichen
                ink = ink * cv2.dilate(((HM['notehead'] > 0) | (HM['clefs_keys'] > 0)).astype(np.uint8), se_)
            else:                    # nur runde Köpfe brauchen die Bestätigung des Netzes; Kreuze und Dreiecke kennt es nicht
                gate = cv2.dilate((HM['notehead'] > 0).astype(np.uint8), se_)
        cov = cv2.boxFilter((HM['notehead'] > 0).astype(np.float32), -1, (max(3, int(d * 0.9)), max(3, int(d * 0.6)))) if (HM is not None and HOMR_ADD > 0) else None
        cand = []; F = []
        # einmal je System rändern statt bei jedem Ausschnitt (gleiches Ergebnis wie train.patch, aber viel schneller)
        r = int(d * 0.9); bwp = cv2.copyMakeBorder(bw, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0); nolp = cv2.copyMakeBorder(nol, r * 2, r * 2, r * 2, r * 2, cv2.BORDER_CONSTANT, value=0)
        def fpatch(pad, x, y):
            x0, y0 = int(round(x)) - r, int(round(y)) - r
            return cv2.resize(pad[y0 + 2 * r:y0 + 4 * r, x0 + 2 * r:x0 + 4 * r], (20, 20), interpolation=cv2.INTER_AREA).ravel()
        for p in range(-12, 13):   # Lage -13 (sechseinhalb Linien unter dem System) kommt nicht vor, dort stehen die Bassbuchstaben
            y = int(round(mid + 0.5 - p * d / 2))
            if not (0 <= y < bw.shape[0]): continue
            row = ink[y]
            for x in range(int(d * 2), bw.shape[1] - int(d), 2):
                if row[x] >= INK:
                    F.append(np.concatenate([fpatch(bwp, x, y), fpatch(nolp, x, y)])); cand.append((x, y, p, float(row[x])))
        if not cand: continue
        F = np.array(F); pr = clf.predict_proba(F)
        if H2 is not None and H2_W > 0 and (G.FOTO_LINES or _os.environ.get('GS_H2ALL')):   # GS_H2ALL: Fotomodell auch auf Scans zumischen (Versuch)
            p2 = np.zeros_like(pr); p2[:, H2.classes_] = H2.predict_proba(F); pr = (1 - H2_W) * pr + H2_W * p2
        if H3 is not None and H3_W > 0 and not G.FOTO_LINES:
            p3 = np.zeros_like(pr); p3[:, H3.classes_] = H3.predict_proba(F); pr = (1 - H3_W) * pr + H3_W * p3
        pt = TRI.predict_proba(F) if TRI is not None else None
        def make_hits(tri_on):
            hits = []
            for i in range(len(cand)):
                tp = pt[i][1] + pt[i][2] if pt is not None else 0
                if tri_on and cand[i][2] % 2 == 0 and tp >= TRI_THR and (1 - pr[i][0] >= TRI_HEAD or tp >= TRI_SURE):   # 4. Reihe liegt immer auf geraden Lagen
                    hits.append((round(float(tp), 3), cand[i][:3], 5 if pt[i][1] >= pt[i][2] else 6))
                elif pr[i][0] < 1 - thr:
                    c_ = int(np.argmax(pr[i][1:])) + 1
                    if gate is not None and CLX[c_].startswith('o') and not gate[int(cand[i][1]), int(cand[i][0])]: continue
                    hits.append((round(1 - pr[i][0], 3), cand[i][:3], c_))
                elif cov is not None and cov[int(cand[i][1]), int(cand[i][0])] >= HOMR_ADD:
                    hits.append((round(float(cov[int(cand[i][1]), int(cand[i][0])]) * 0.6, 3), cand[i][:3], 1 if pr[i][1] >= pr[i][2] else 2))
            return hits
        # Formgüte: Übereinstimmung mit einer Kopfvorlage genau an der Stelle (trennt echte Köpfe von Treffern eine Halbstufe daneben)
        T = G.templates(d)
        def good(x, y, c):
            t = T['o_f' if CLX[c].startswith('o') else ('t_f' if CLX[c].startswith('t') else 'x_f')]; hh, ww = t.shape
            y0, x0 = int(round(y)) - hh // 2, int(round(x)) - ww // 2
            win = nol[max(0, y0):y0 + hh, max(0, x0):x0 + ww].astype(np.float32)
            if win.shape != t.shape: return 0.0
            tt = t.astype(np.float32); return float((win * tt).sum() / (np.sqrt((win * win).sum() * (tt * tt).sum()) + 1e-6))
        bm = cv2.morphologyEx(nol, cv2.MORPH_OPEN, cv2.getStructuringElement(cv2.MORPH_RECT, (max(3, int(d * BEAM_L)), 1))) if BEAM else None
        def select(hits):
            ref = sorted(((sc, good(x, y, c), (x, y, p), c) for sc, (x, y, p), c in hits), key=lambda r: (-r[0], -r[1]))
            # Dreieck-Meldungen haben Vorrang vor runden Köpfen an derselben Stelle (der Kopf-Klassifikator hält Dreiecke oft für runde Köpfe)
            ref = sorted(ref, key=lambda r: -(r[0] * 0.5 + r[1] + (TRI_BOOST if CLX[r[3]].startswith('t') else 0)))
            keep = []
            for sc, gq, (x, y, p), c in ref:
                # Dreieck neben Kreuz oder rundem Kopf eine Lage tiefer/höher (Sekunde links und rechts vom Hals) darf näher stehen
                def apart(k):
                    if abs(k[2][2] - p) > 1: return True
                    lim = d * (1.3 if k[2][2] == p else (TRI_NB if (k[3] >= 5) != (c >= 5) else 0.8))
                    return abs(x - k[2][0]) > lim
                if all(apart(k) for k in keep): keep.append((sc, gq, (x, y, p), c))
            keep = [(sc, xyp, c) for sc, gq, xyp, c in keep]
            # Balken (dicke waagrechte Striche unter/über Notengruppen) sind keine Köpfe: Meldungen, deren Kopffläche
            # überwiegend auf einem langen waagrechten Strich liegt, verwerfen
            if bm is not None:
                def on_beam(x, y):
                    hx, hy = int(d * 0.6), int(d * 0.35); win = bm[max(0, int(y) - hy):int(y) + hy + 1, max(0, int(x) - hx):int(x) + hx + 1]
                    return win.size and (win > 0).mean() > BEAM_F
                keep = [k for k in keep if not on_beam(k[1][0], k[1][1])]
            return keep
        k_tri = select(make_hits(True)); k_no = select(make_hits(False))
        staves_res.append((si, d, mid, k_tri, k_no))
    # Dreiecke gibt es nur in der Schreibweise mit Dreieck-Köpfen (z. B. Elke Margetich). Auf einer Seite ohne solche Köpfe meldet
    # der Dreieck-Klassifikator nur vereinzelt und unsicher etwas; dann gilt die Lesung ohne Dreiecke für die ganze Seite.
    tri = [sc for _, _, _, kt, _ in staves_res for sc, _, c in kt if c >= 5]
    tri_page = len(tri) >= TRI_MIN and sum(t >= 0.9 for t in tri) >= TRI_FRAC * len(tri)
    LAST['tri_page'] = bool(tri_page)   # für gs2xml: Seite in Dreieck-Schreibweise (Margetich)
    for si, d, mid, k_tri, k_no in staves_res:
        keep = k_tri if tri_page else k_no
        for sc, (x, y, p), c in keep:
            out.append({'x': x / scale, 'y': y / scale, 'p': p, 't': CLX[c][0], 'fill': CLX[c].endswith('f'), 'sc': round(float(sc), 2), 'st': si, 'd': d / scale, 'mid': mid / scale})
    return out

if __name__ == '__main__':
    files = json.load(open(HERE+'/test.json'))
    tot = dict(gt=0, det=0, hit=0, pos=0, typ=0)
    for f in files:
        gt = json.load(open(f[:-4] + '.json')); det = detect(f); mid = gt['staff'][2]; used = set()
        for h in gt['heads']:
            best = None
            for i, x in enumerate(det):
                if i in used: continue
                dd = abs(x['x'] - h['x']) + abs(x['y'] - h['y'])
                if dd < 6 and (best is None or dd < best[0]): best = (dd, i)
            if best:
                used.add(best[1]); x = det[best[1]]; tot['hit'] += 1
                tot['pos'] += x['p'] == round((mid - h['y']) / 4); tot['typ'] += x['t'] == h['t']
        tot['gt'] += len(gt['heads']); tot['det'] += len(det)
    print(tot, 'Trefferquote', round(tot['hit'] / tot['gt'], 3), 'Genauigkeit', round(tot['hit'] / max(1, tot['det']), 3))
