# Vektor-PDFs (Notensatzprogramme wie MC-Score oder Finale): Die Notenköpfe stehen als Zeichen einer Notenschrift in der Textebene,
# mit genauer Lage und eindeutiger Form. vektor.py liest sie dort aus, statt sie im Bild zu suchen. Erkannt werden nur Schriften,
# deren Zeichen hier zugeordnet sind; andere PDFs (Scans) werden wie bisher über das Bild gelesen.
import json, os

# Zeichen -> (Kopfform, voll): 'o' runder Kopf, 'x' Kreuz; voll = False bei hohlen Köpfen (Halbe, Ganze) und beim Kreuz im Kreis
FONTS = {
    'MC6': {'E': ('o', True), 'D': ('o', False), 'O': ('x', True), 'P': ('x', False)},      # MC-Score 6 (z. B. Günther Pacher)
    'Maestro': {'œ': ('o', True), '˙': ('o', False), 'w': ('o', False),          # Finale, Schrift Maestro
                '‹': ('x', True), '⊗': ('x', False)},
}


# senkrechte Verschiebung von der Grundlinie zur Kopfmitte, in Schriftgrößen (gemessen an den Bildern)
DY = {('MC6', 'E'): -0.18, ('MC6', 'D'): -0.18, ('MC6', 'O'): -0.04, ('MC6', 'P'): -0.04}


def _font_key(name):
    for k in FONTS:
        if name == k or name.endswith('+' + k): return k
    return None


def heads_pdf_page(page, dpi):
    """Köpfe einer PDF-Seite (PyMuPDF) in Bildpunkten einer Wiedergabe mit dpi: Liste von {x, y, t, fill}; leer, wenn keine bekannte Notenschrift."""
    z = dpi / 72.0; out = []
    raw = page.get_text('rawdict')
    for b in raw.get('blocks', []):
        for l in b.get('lines', []):
            for s in l.get('spans', []):
                fk = _font_key(s.get('font', ''))
                if not fk: continue
                for ch in s.get('chars', []):
                    m = FONTS[fk].get(ch['c'])
                    if not m: continue
                    x0, y0, x1, y1 = ch['bbox']; ox, oy = ch['origin']
                    w = min(x1 - x0, s['size'] * 0.45) if x1 > x0 else s['size'] * 0.3
                    out.append(dict(x=(ox + w / 2) * z, y=(oy + DY.get((fk, ch['c']), 0.0) * s['size']) * z, t=m[0], fill=m[1], c=fk + ':' + ch['c'], size=s['size']))
    return out


BASS_OK = set('ABCDEFGHabcdefgh')


def bass_pdf_page(page, dpi):
    """Bassbuchstaben (einzelne Buchstaben A–H, a–h in einer Textschrift, nicht in der Notenschrift) in Bildpunkten; X-Bass und
    Helikon-Marken (H1, H2, H(k)) werden nicht übernommen, weil der Wandler sie nicht als Ton kennt."""
    z = dpi / 72.0; out = []
    raw = page.get_text('rawdict')
    for b in raw.get('blocks', []):
        for l in b.get('lines', []):
            for s in l.get('spans', []):
                if _font_key(s.get('font', '')): continue
                txt = ''.join(ch['c'] for ch in s.get('chars', []))
                chars = s.get('chars', [])
                for i, ch in enumerate(chars):
                    c = ch['c']
                    if c not in BASS_OK: continue
                    prev = chars[i - 1]['c'] if i else ' '; nxt = chars[i + 1]['c'] if i + 1 < len(chars) else ' '
                    if prev.strip() or nxt.strip(): continue   # nur allein stehende Buchstaben (keine Wörter, kein H1, H(k))
                    x0, y0, x1, y1 = ch['bbox']
                    out.append(dict(x=(x0 + x1) / 2 * z, y=(y0 + y1) / 2 * z, c=c))
    return out


def save_sidecars(pdf_path, pages_png, dpi):
    """Für jede wiedergegebene Seite eine Datei <png>.vk.json mit den Köpfen anlegen (nur wenn welche gefunden wurden)."""
    import fitz
    doc = fitz.open(pdf_path); n = 0
    for i, fn in enumerate(pages_png):
        if i >= len(doc): break
        hs = heads_pdf_page(doc[i], dpi)
        if len(hs) >= 10:
            json.dump(dict(dpi=dpi, heads=hs, bass=bass_pdf_page(doc[i], dpi)), open(fn + '.vk.json', 'w')); n += 1
    return n


def load(path):
    fn = path + '.vk.json'
    if os.environ.get('GS_VEKTOR', '1') != '1' or not os.path.exists(fn): return None
    try: return json.load(open(fn))['heads']
    except Exception: return None


def det_from_heads(heads, fx, st):
    """Köpfe (Bildpunkte der Originalseite) -> Einträge wie detect.detect (Lage p zum System, Form, voll), Bild um fx skaliert."""
    import numpy as np
    out = []
    mids = [(ys[2], (ys[4] - ys[0]) / 4) for ys in st]
    for h in heads:
        x, y = h['x'] * fx, h['y'] * fx
        if not mids: break
        si = int(np.argmin([abs(y - m) for m, _ in mids])); mid, d = mids[si]
        if abs(y - mid) > 9 * d: continue
        p = int(round((mid + 0.5 - y) * 2 / d))
        out.append(dict(x=x, y=y, p=p, t=h['t'], fill=bool(h['fill']), sc=1.0, st=si, d=d, mid=mid))
    return out


def load_bass(path):
    fn = path + '.vk.json'
    if os.environ.get('GS_VEKTOR', '1') != '1' or os.environ.get('GS_VEKTORBASS', '1') != '1' or not os.path.exists(fn): return None
    try: return json.load(open(fn)).get('bass')
    except Exception: return None


def letters_from_bass(bass, fx, st, si):
    """Bassbuchstaben unter System si im Format von bassread.read_letters: Liste (x, Name), übereinander stehende als 'C+c'."""
    ys = st[si]; d = (ys[4] - ys[0]) / 4; lo = ys[4] + 0.8 * d
    hi = min(ys[4] + 11 * d, st[si + 1][0] - 1.5 * d) if si + 1 < len(st) else ys[4] + 11 * d
    L = sorted([(b['x'] * fx, b['y'] * fx, b['c']) for b in (bass or []) if lo <= b['y'] * fx <= hi])
    if not L: return []
    import numpy as np
    yc = float(np.median([y for _, y, _ in L]))
    row = [(x, y, c) for x, y, c in L if abs(y - yc) < 0.9 * d]
    extra = [(x, y, c) for x, y, c in L if 0.9 * d <= abs(y - yc) <= 2.6 * d]
    out = []
    for x, y, c in row:
        st_ = [c2 for x2, y2, c2 in extra if abs(x2 - x) < 0.7 * d]
        n = '+'.join(sorted([c] + st_[:1], key=lambda q: q[0].islower())) if st_ else c
        out.append((x, n))
    return out
