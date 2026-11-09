# Vergleich der Kopf-Erkenner auf Blättern ohne bekannte Lösung (Fotos, Handschrift, eingesandte PDFs): gezählt werden die Takte,
# die nicht aufgehen (Summe der Notenwerte ≠ Taktlänge), und die gelesenen Griffe, über den ganzen Leser (gs2xml.read).
# Ein Lauf je Erkenner, weil das gelernte Modell je Vorgang nur einmal geladen wird:
#   python3 yolo_vergleich_fotos.py alt BILD …
#   GS_YOLOMODELL=…/modell.onnx python3 yolo_vergleich_fotos.py yolo BILD …
# PDFs werden seitenweise als Bild gelesen (pymupdf, 200 dpi).
import sys, os, tempfile
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, '..', 'gsleser'))
os.environ['GS_YOLOFOTO'] = '1'
import gs2xml as X


def seiten(p):
    if not p.lower().endswith('.pdf'): return [p]
    import fitz
    out = []
    for i, pg in enumerate(fitz.open(p)):
        f = os.path.join(tempfile.gettempdir(), f'_yvf_{os.getpid()}_{i}.png'); pg.get_pixmap(dpi=200).save(f); out.append(f)
    return out


if __name__ == '__main__':
    X.ERKENNER = sys.argv[1]; S = dict(takte=0, offen=0, griffe=0)
    for p in sys.argv[2:]:
        t = o = gr = 0
        for f in seiten(p):
            try: meas, L = X.read(f)
            except Exception as e: print(os.path.basename(p), 'Fehler', str(e)[:60]); continue
            for m, Lm in zip(meas, L):
                t += 1; o += sum(e['u'] for e in m) != Lm; gr += sum(1 for e in m if not e['rest'])
        print(f'{os.path.basename(p)}: Takte {t} offen {o} Griffe {gr}', flush=True)
        S['takte'] += t; S['offen'] += o; S['griffe'] += gr
    print('GESAMT', X.ERKENNER, os.environ.get('GS_YOLOMODELL', ''), S)
