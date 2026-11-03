# Bassbuchstaben samt Lage aus Philipps Heft (PDF): Die Buchstaben stehen dort als Text in Comic Sans, aber ohne Zeichenzuordnung und mit
# einer eigenen Teilschrift je Seite. Gleiche Buchstaben haben in allen Teilschriften denselben Umriss; die Umrisse werden daher einmal
# von Hand benannt (LAB), danach ist jeder Buchstabe jeder Seite bekannt. Ergebnis: ph_letters.json mit (Seite, Buchstabe, Rechteck in PDF-Punkten).
import sys, io, json, hashlib, pymupdf
from fontTools.ttLib import TTFont
from fontTools.pens.recordingPen import RecordingPen
PDF = sys.argv[1] if len(sys.argv) > 1 else '/home/claude/philipp/pdf/Griffschrift.pdf'
d = pymupdf.open(PDF); MAP = {}; U = []
for pn in range(len(d)):
    for f in d[pn].get_fonts():
        if 'Comic' not in f[3]: continue
        t = TTFont(io.BytesIO(d.extract_font(f[0])[3])); gs = t.getGlyphSet()
        for gi, gn in enumerate(t.getGlyphOrder()):
            rp = RecordingPen(); gs[gn].draw(rp); h = hashlib.md5(repr(rp.value).encode()).hexdigest()[:10]
            MAP[(pn + 1, gi)] = h
            if h not in U: U.append(h)
# Reihenfolge der Umrisse wie beim ersten Durchlauf (von Hand angesehen): 0 leer, B b A C c ' W a D d 0 Leerzeichen ( )
LAB = [None, 'B', 'b', 'A', 'C', 'c', "'", 'W', 'a', 'D', 'd', '0', ' ', '(', ')']
H2L = {h: LAB[i] for i, h in enumerate(U) if i < len(LAB)}
out = []
for pn in range(len(d)):
    for b in d[pn].get_text('rawdict')['blocks']:
        for l in b.get('lines', []):
            for s in l['spans']:
                if 'Comic' not in s['font']: continue
                for ch in s['chars']:
                    lab = H2L.get(MAP.get((pn + 1, ord(ch['c']))))
                    if lab and lab.strip(): out.append((pn + 1, lab, ch['bbox']))
json.dump(out, open('ph_letters.json', 'w')); print(len(out))
