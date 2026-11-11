# Philipps Heft (nur Messung, Seiten und Griffe nicht im Repository): Griffe UND Rhythmus gegen seine bekannten Griffe (ph_pins.json),
# je für drei Lesarten: Leser allein | Gerüst aus homr, Dreieck-Seiten ohne homr (Stand v158) | Gerüst mit runden Köpfen für homr (v159/v160).
# Rhythmus: für je zwei aufeinanderfolgende zugeordnete Akkorde wird der gelesene Abstand der Einsätze (aus den Notenwerten,
# Pausen eingerechnet) mit dem Abstand in Philipps Griffen verglichen; Maßstab ist das Verhältnis, das am häufigsten vorkommt.
# Aufruf (Arbeitsordner /tmp/philipp mit pages/ und ph_pins.json aus teil4_philipp_fuchs.zip, philipp_pruefen.py dort):
#   GS_YOLO=1 GS_YOLOMODELL=…/griffschrift_koepfe2.onnx python3 philipp_rhythmus.py ../gsleser
import sys, os, collections
GS = sys.argv[1]; sys.path.insert(0, GS); sys.path.insert(0, '/tmp/philipp')
import gs2xml as X
from align import align
from philipp_pruefen import PAGES, PINS, pin_chords
from philipp_pruefen_gs2xml import D450, key

def lesen(k, art):
    X.GERUEST = art != 'leser'; os.environ['GS_GERUESTTRI'] = '0' if art == 'bisher' else '1'
    ev = []; t = 0
    for p in PAGES[k]:
        meas, L = X.read(f'/tmp/philipp/pages/p-{p:02d}.png', -2, table=D450)
        for m in meas:
            for e in m:
                if not e['rest']: ev.append((t, e))
                t += e['u']
    return ev

ARTEN = ('leser', 'bisher', 'rund')
T = {a: collections.Counter() for a in ARTEN}
for k in PAGES:
    B = pin_chords(PINS[k]); Bs = [{(r, n) for r, n, _ in b} for b in B]; TS = sorted({round(t_, 3) for t_, *_ in PINS[k]})
    zeile = [k]
    for a in ARTEN:
        ev = lesen(k, a); A = [{(r, n) for r, n, _ in e['heads']} for _, e in ev]
        pr = align([{key(r, n) for r, n in x} for x in A], [{key(r, n) for r, n in b} for b in Bs])
        st = T[a]; st['pins'] += len(B); st['paar'] += len(pr); st['gleich'] += sum(A[i] == Bs[j] for i, j in pr)
        # Rhythmus: Paare (i, j) und (i+1, j+1) direkt hintereinander
        ver = []
        for (i, j), (i2, j2) in zip(pr, pr[1:]):
            if i2 == i + 1 and j2 == j + 1 and j2 < len(TS):
                dl = ev[i2][0] - ev[i][0]; dp = TS[j2] - TS[j]
                if dl > 0 and dp > 0: ver.append(round(dl / dp, 2))
        if ver:
            sk = collections.Counter(ver).most_common(1)[0][0]
            ok = sum(abs(v - sk) <= 0.02 * sk for v in ver)
            st['rpaare'] += len(ver); st['rok'] += ok
            zeile.append(f"{a}: Griffe {sum(A[i] == Bs[j] for i, j in pr)}/{len(B)} Rhythmus {ok}/{len(ver)}")
    print(' | '.join(zeile), flush=True)
print()
for a in ARTEN:
    s = T[a]
    print(f"{a}: Akkorde ganz gleich {s['gleich']}/{s['pins']} ({s['gleich']/max(1,s['pins']):.3f}), Rhythmus richtig {s['rok']}/{s['rpaare']} ({s['rok']/max(1,s['rpaare']):.3f})")
