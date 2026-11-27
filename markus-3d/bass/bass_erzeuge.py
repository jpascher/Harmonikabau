# Bassteil-Variante erzeugen (23.11.2026)
# Aufruf: python3 bass_erzeuge.py <teilung 18|19|20> <laenge_oder_diskant> <ausgabe.step> [reihenabstand 14]
# Laenge: Endmass des Gehaeuses; mindestens 380. Seit 27.11.2026 genau wie eingegeben (auf 0,1 mm), damit Bass und
# Diskant dieselben Aussenmasse bekommen (vorher auf volle 10 mm aufgerundet). Die Rohteile der Zeichnungen haben wie
# bisher 2 mm Uebermass (Front 382 bei 380).
import sys, os, math, subprocess
Q = os.path.dirname(os.path.abspath(__file__))                 # Bass-Quellordner (dieser Ordner)
T = float(sys.argv[1]); LD = float(sys.argv[2]); OUT = os.path.abspath(sys.argv[3])
RA = sys.argv[4] if len(sys.argv) > 4 else "14"
L_END = max(380.0, round(LD, 1))
Z = os.path.join(Q, "varianten", f"T{T:g}_L{L_END:g}")   # Arbeitskopie je Variante
subprocess.run([sys.executable, os.path.dirname(os.path.abspath(__file__)) + "/bass_variante.py", Q, Z, str(T), RA], check=True)
for t in ("front", "rueck", "seite"):
    subprocess.run([sys.executable, "teile.py", "rahmen.dxf", f"Markus_{t}_v4.step", t], cwd=Z, check=True, stdout=subprocess.DEVNULL)
env = dict(os.environ, BASS_DL=str(L_END - 380.0))
print(f"Bassteil: Teilung {T:g} mm, Reihenabstand {RA} mm, Gehaeuse {L_END:g} mm (Rohmass {L_END + 2:g})")
r = subprocess.run([sys.executable, "kammern.py", "boden2.dxf", OUT], cwd=Z, env=env)
sys.exit(r.returncode)
