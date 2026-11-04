# Lesung eines Fotos oder Scans gegen ein von Hand abgelesenes Soll: Einsätze (System, x) und Dauern (-v: falsche Dauern einzeln).
# Das Soll (soll_<name>.py mit TAKTE und S = {System: [(x, Dauer), …]}) liegt nur lokal, nicht im Repo (Blätter fremder Hefte).
import sys, importlib; sys.path.insert(0, '/home/claude/gsw/m'); sys.path.insert(0, '/home/claude/lab')
import gs2xml as X
key, img = sys.argv[1], sys.argv[2]; T = importlib.import_module('soll_' + key)
meas, L = X.read(img, 0, bass=[])
ev = [(e['dbg'][0], e['dbg'][2], e['u'], e['rest']) for m in meas for e in m]
n = hit = dur = 0; used = set()
for si, row in T.S.items():
    for t in row:
        x, u = t[0], t[1]; n += 1
        c = [i for i, e in enumerate(ev) if i not in used and e[0] == si and abs(e[1] - x) < 28]
        if not c: continue
        i = min(c, key=lambda i: abs(ev[i][1] - x)); used.add(i); hit += 1; dur += ev[i][2] == u
        if ev[i][2] != u and '-v' in sys.argv: print('  falsch', si, x, 'soll', u, 'gelesen', ev[i][2])
print(key, 'Soll-Einsätze', n, 'gefunden', hit, 'Dauer richtig', dur, '| überzählig', len(ev) - len(used), '| Takte', len(meas), 'Soll', T.TAKTE + 1)
