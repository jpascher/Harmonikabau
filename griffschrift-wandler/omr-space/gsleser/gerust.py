# Gerüst aus der Notenerkennung (homr, Versuch, GS_GERUEST=1): homr liest jedes System vollständig (Takte, Einsätze, Notenwerte,
# Pausen, Akkorde, Lage der Köpfe). Der Griffschrift-Leser liefert nur noch Kopfform, Reihe, Balgrichtung und Bassbuchstaben.
# Beides wird je System über die Lagen der Köpfe zusammengeführt; Notenwerte, Pausen und Taktstriche kommen von homr.
# Die Tonhöhe, die homr meldet, wird nur als Lage auf dem System verwendet (Schlüssel beachtet); Vorzeichen spielen keine Rolle.
import os, re
from fractions import Fraction
import cv2, numpy as np

DIAT = 'CDEFGAB'
CLEF_REF = {'G': ('G', 4), 'F': ('F', 3), 'C': ('C', 4)}
MINMATCH = float(os.environ.get('GS_GMIN', '0.5'))   # Anteil zugeordneter Akkorde, ab dem homr für ein System gilt
MAXLOSS = float(os.environ.get('GS_GLOSS', '0.15'))   # bleibt mehr als dieser Anteil der Akkorde des Lesers ohne homr-Gegenstück, gilt für das ganze System der Leser
TAKTWEISE = os.environ.get('GS_GTAKT', '1') == '1'
QUAL = float(os.environ.get('GS_GQUAL', '0.6'))   # so viele von homrs eigenen Takten im System müssen gleich lang sein, sonst gilt der Leser   # Takte mit Akkorden, die nur der Leser sieht: diese Takte vom Leser


def _diat(step, octv):
    return DIAT.index(step) + 7 * int(octv)


def _mid(clef):
    m = re.match(r'clef[-_]([GFC])(\d)', clef or 'clef_G2')
    if not m: return _diat('B', 4)
    s, o = CLEF_REF[m.group(1)]
    return _diat(s, o) + 2 * (3 - int(m.group(2)))


def _units(rh):
    """homr-Notenwert (Kern: note_8, note_4., note_24 …) -> Zwölftel-Viertel und Triole"""
    m = re.match(r'(?:note|rest)[-_](\d+)(G?)(\.*)$', rh)
    if not m or m.group(2): return None, False
    n = int(m.group(1))
    if n == 0: return None, False
    u = Fraction(48, n)
    add = u / 2
    for _ in m.group(3): u += add; add /= 2
    if u.denominator != 1: return None, False
    trip = n % 3 == 0
    return int(u), trip


_RUN = {}
def homr_systems(g, st, x_ext):
    """Jedes System des Lesers (fünf Linienhöhen, Ausdehnung) mit homr lesen; die Systeme selbst kommen vom Leser, weil homrs eigene
    Systemsuche Anker (Schlüssel, Taktstriche) braucht, die in Griffschrift oft fehlen. Ergebnis je System: Liste von Ereignissen
    ('bar',) | ('rest', u, trip) | ('note', u, trip, [Lagen]) oder None."""
    import zlib
    key = (g.shape, zlib.crc32(g[::5, ::5].tobytes()), len(st))
    if key in _RUN: return _RUN[key]
    _dc = os.environ.get('GS_GCACHE')   # Ordner für zwischengespeicherte homr-Lesungen (nur Prüfläufe)
    _fn = os.path.join(_dc, '%d_%d_%x_%d.json' % (g.shape[0], g.shape[1], key[1], len(st))) if _dc else None
    if _fn and os.path.exists(_fn):
        import json
        out = [None if e is None else [tuple(x) for x in e] for e in json.load(open(_fn))]
        _RUN.clear(); _RUN[key] = out; return out
    from homr.main import download_weights
    from homr.debug import Debug
    from homr.model import MultiStaff, Staff, StaffPoint
    from homr.staff_parsing import parse_staff_image
    from homr.staff_regions import StaffRegions
    from homr.resize import resize_image
    from homr.transformer.configs import Config
    download_weights(False, False, False)
    from homr import color_adjust
    img = color_adjust.apply_clahe(resize_image(cv2.cvtColor(g, cv2.COLOR_GRAY2BGR))); sc = img.shape[1] / g.shape[1]
    if img.ndim == 3: img = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    dbg = Debug(img, 'gerust.png', False)
    tc = Config(); tc.use_gpu_inference = False; tc.use_coreml_encoder = False
    staffs = []
    for si, ys in enumerate(st):
        x0, x1 = x_ext.get(si, (0, g.shape[1] - 1))
        pts = [StaffPoint(float(x) * sc, [float(y) * sc for y in ys], 0) for x in np.arange(x0, x1, 10)]
        staffs.append(MultiStaff([Staff(pts)], []))
    regions = StaffRegions(staffs)
    out = []
    for si, ms in enumerate(staffs):
        try: toks = parse_staff_image(dbg, si, ms.staffs[0], img, regions, tc)
        except Exception: out.append(None); continue
        ev = []; clef = 'clef_G2'; chord = False
        for s in toks:
            r = s.rhythm
            if r.startswith('clef'): clef = r; continue
            if r == 'chord': chord = True; continue
            if 'barline' in r or 'repeat' in r: ev.append(('bar',)); chord = False; continue
            if r.startswith('rest'):
                u, tr = _units(r)
                if u: ev.append(('rest', u, tr))
                chord = False; continue
            if r.startswith('note'):
                u, tr = _units(r)
                m = re.match(r'([A-G])(\d)', s.pitch or '')
                if not u or not m: chord = False; continue
                p = _diat(m.group(1), m.group(2)) - _mid(clef)
                if chord and ev and ev[-1][0] == 'note': ev[-1][3].append(p)
                else: ev.append(('note', u, tr, [p]))
                chord = False
        out.append(ev)
    if _fn:
        import json; os.makedirs(_dc, exist_ok=True); json.dump(out, open(_fn, 'w'))
    _RUN.clear(); _RUN[key] = out
    return out


def _score(R, H):
    R = set(R); H = set(H); inter = len(R & H)
    near = sum(1 for h in H - R if (h - 1) in R or (h + 1) in R)
    if not inter and not near: return -2.0
    return 2 * inter + 0.7 * near - 0.5 * (len(R - H) + len(H - R) - near)


def _align(A, B):
    n, m = len(A), len(B); S = np.zeros((n + 1, m + 1)); S[:, 0] = -np.arange(n + 1); S[0, :] = -np.arange(m + 1); bt = {}
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            c = [(S[i - 1, j - 1] + _score(A[i - 1], B[j - 1]), 'm'), (S[i - 1, j] - 1, 'a'), (S[i, j - 1] - 1, 'b')]
            S[i, j], bt[i, j] = max(c)
    i, j, pr = n, m, []
    while i > 0 and j > 0:
        k = bt[i, j]
        if k == 'm':
            if _score(A[i - 1], B[j - 1]) > 0: pr.append((i - 1, j - 1))
            i -= 1; j -= 1
        elif k == 'a': i -= 1
        else: j -= 1
    return pr[::-1]


def rebuild(g, st, x_ext, meas, medge, button, pitch, table, off, stats=None):
    """meas/medge des Lesers je System durch das homr-Gerüst ersetzen (wo homr das System lesen konnte und genug zusammenpasst)."""
    try: H = homr_systems(g, st, x_ext)
    except Exception as e:
        if stats is not None: stats['fehler'] = str(e)[:100]
        return meas, medge
    newm, newe = [], []
    for si, ys in enumerate(st):
        d = (ys[4] - ys[0]) / 4
        mine = [(m, e) for (s_, m), e in zip(meas, medge) if s_ == si]
        cand = [H[si]] if si < len(H) and H[si] else []
        groups = [gr for m, _ in mine for gr in m if gr[0]['kind'] == 'n']
        if not cand or not groups:
            newm += [(si, m) for m, _ in mine]; newe += [e for _, e in mine]; continue
        ev = cand[0]
        # Selbstprüfung: homrs eigene Takte im System (ohne ersten und letzten) sollen überwiegend gleich lang sein;
        # sonst hat homr das System nicht verstanden (z. B. dichte Akkorde mit Dreiecken bei Philipp) und der Leser gilt
        sm_ = []; c_ = 0; seen = False
        for e in ev:
            if e[0] == 'bar':
                if seen: sm_.append(c_)
                c_ = 0; seen = False
            else: c_ += e[1]; seen = True
        if seen: sm_.append(c_)
        inner = sm_[1:-1] if len(sm_) > 3 else sm_
        if inner:
            import collections as _co
            top = _co.Counter(inner).most_common(1)[0][1]
            if top < QUAL * len(inner):
                newm += [(si, m) for m, _ in mine]; newe += [e for _, e in mine]
                if stats is not None: stats['unstimmig'] = stats.get('unstimmig', 0) + 1
                continue
        A = [[p for it in gr for p in it.get('ps', [])] for gr in groups]
        hn = [k for k, e in enumerate(ev) if e[0] == 'note']
        # homr liest die Lagen in manchen Systemen durchgehend um eine oder zwei Stufen verschoben (Schlüssel, Linienlage):
        # die Verschiebung mit der besten Zuordnung gilt für das ganze System
        best = None
        for sh in (0, -1, 1, -2, 2):
            Bs = [[p + sh for p in ev[k][3]] for k in hn]; pr_ = _align(A, Bs)
            v = sum(_score(A[i], Bs[j]) for i, j in pr_) - 0.5 * abs(sh)
            if best is None or v > best[0]: best = (v, sh, pr_, Bs)
        _, hsh, pr, B = best
        if stats is not None: stats.setdefault('verschiebung', []).append(hsh)
        if stats is not None:
            stats.setdefault('systeme', 0); stats['systeme'] += 1
            stats.setdefault('paare', 0); stats['paare'] += len(pr)
            stats.setdefault('leser', 0); stats['leser'] += len(A)
            stats.setdefault('homr', 0); stats['homr'] += len(B)
        if len(pr) < MINMATCH * max(len(A), len(B)):
            newm += [(si, m) for m, _ in mine]; newe += [e for _, e in mine]
            if stats is not None: stats['abgelehnt'] = stats.get('abgelehnt', 0) + 1
            continue
        # Lücken zwischen zugeordneten Paaren: gleich viele übrige Akkorde des Lesers und Noten von homr der Reihe nach paaren
        # (homr liest Kreuze oft eine Lage daneben; Kopfform und Lage gelten dann vom Leser, Notenwert von homr)
        pr2 = list(pr); bnd = [(-1, -1)] + pr + [(len(A), len(B))]
        for (i0, j0), (i1, j1) in zip(bnd, bnd[1:]):
            if i1 - i0 == j1 - j0 and i1 - i0 > 1:
                pr2 += [(i0 + q, j0 + q) for q in range(1, i1 - i0)]
        pr = sorted(pr2)
        if len(A) - len(pr) > MAXLOSS * len(A):   # homr hat zu viele Akkorde nicht gesehen: System wie bisher vom Leser
            newm += [(si, m) for m, _ in mine]; newe += [e for _, e in mine]
            if stats is not None: stats['zuvielverlust'] = stats.get('zuvielverlust', 0) + 1
            continue
        match = {hn[j]: i for i, j in pr}
        # x je homr-Ereignis: zugeordnete Akkorde liefern x, die übrigen werden dazwischen verteilt
        xs = [None] * len(ev)
        for k, i in match.items(): xs[k] = float(np.mean([it['x'] for it in groups[i]]))
        known = [k for k in range(len(ev)) if xs[k] is not None]
        for k in range(len(ev)):
            if xs[k] is not None: continue
            L = [q for q in known if q < k]; R = [q for q in known if q > k]
            if L and R: a, b = L[-1], R[0]; xs[k] = xs[a] + (k - a) / (b - a) * (xs[b] - xs[a])
            elif L: xs[k] = xs[L[-1]] + (k - L[-1]) * 1.5 * d
            elif R: xs[k] = xs[R[0]] - (R[0] - k) * 1.5 * d
            else: xs[k] = 0.0
        def ndir(k):
            q = min(match, key=lambda q: abs(q - k)) if match else None
            return groups[match[q]][0]['dir'] if q is not None else 'Z'
        cur = []; ms_ = []
        for k, e in enumerate(ev):
            if e[0] == 'bar':
                if cur: ms_.append(cur); cur = []
                continue
            if e[0] == 'rest':
                cur.append([dict(kind='r', st=si, x=xs[k], u=e[1], hu=e[1], htrip=e[2])])
            elif k in match:
                gr = [dict(it) for it in groups[match[k]]]
                gr[0]['hu'] = e[1]; gr[0]['htrip'] = e[2]; cur.append(gr)
            else:   # Note, die nur homr kennt: runde Köpfe an seinen Lagen, Balgrichtung wie die nächste zugeordnete Note
                dr = ndir(k); hs = []
                for p in sorted(set(q + hsh for q in e[3])):
                    r_, n_ = button(p, 'o'); pi = pitch(r_, n_, dr, table)
                    if pi: hs.append((r_, n_, pi + off))
                if hs: cur.append([dict(kind='n', st=si, x=xs[k], dir=dr, heads=hs, fi=-1, hy=[], hu=e[1], htrip=e[2], nurhomr=True)])
                elif True: cur.append([dict(kind='r', st=si, x=xs[k], u=e[1], hu=e[1], htrip=e[2])])
        if cur: ms_.append(cur)
        # Taktstrich, den homr übersehen hat: Ist ein Takt länger als die übliche Taktlänge und liegt ein Taktstrich des Lesers
        # zwischen zwei seiner Ereignisse, wird dort geteilt
        sums = [sum(gr[0]['hu'] for gr in m) for m in ms_]
        if sums:
            import collections
            Lh = collections.Counter(sums).most_common(1)[0][0]
            rbars = [e[1] for _, e in mine][:-1]
            spl = []
            for m, sm in zip(ms_, sums):
                parts = [m]
                if sm > Lh:
                    for bx in rbars:
                        last = parts[-1]
                        for q in range(1, len(last)):
                            if last[q - 1][0]['x'] < bx < last[q][0]['x']:
                                parts[-1:] = [last[:q], last[q:]]; break
                    if len(parts) > 1 and stats is not None: stats['geteilt'] = stats.get('geteilt', 0) + len(parts) - 1
                spl += parts
            ms_ = spl
        if stats is not None:
            stats['nurhomr'] = stats.get('nurhomr', 0) + sum(1 for k, e in enumerate(ev) if e[0] == 'note' and k not in match)
            stats['nurleser'] = stats.get('nurleser', 0) + len(A) - len(pr)
        matched_r = {i for i, _ in pr}; spans = []
        for q, m in enumerate(ms_):
            x0 = min(it['x'] for gr in m for it in gr); x1 = max(it['x'] for gr in m for it in gr)
            a = (x0 + max(it['x'] for gr in ms_[q - 1] for it in gr)) / 2 if q else x0 - 2 * d
            b = (x1 + min(it['x'] for gr in ms_[q + 1] for it in gr)) / 2 if q + 1 < len(ms_) else x1 + 3 * d
            spans.append((m, (a, b, d)))
        # Takte, in denen der Leser Akkorde sieht, die homr fehlen: dort gelten die Takte des Lesers (Notenwerte vom Leser),
        # damit keine Töne verloren gehen; der Rest des Systems behält homrs Rhythmus
        if TAKTWEISE:
            ux = [float(np.mean([it['x'] for it in groups[i]])) for i in range(len(A)) if i not in matched_r]
            rc = [float(np.mean([it['x'] for gr in m for it in gr])) for m, _ in mine]
            used = set(); out_ = []
            for m, (a, b, d_) in spans:
                if any(a <= x < b for x in ux):
                    rm = [q for q, c in enumerate(rc) if a <= c < b and q not in used]
                    if rm:
                        for q in rm: used.add(q); out_.append(mine[q])
                        if stats is not None: stats['takte_leser'] = stats.get('takte_leser', 0) + 1
                        continue
                if any(a <= c < b for q, c in enumerate(rc) if q in used): continue   # schon durch Leser-Takte ersetzt
                out_.append((m, (a, b, d_)))
            spans = out_
        for m, e in spans: newm.append((si, m)); newe.append(e)
    return newm, newe
