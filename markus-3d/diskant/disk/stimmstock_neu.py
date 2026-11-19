# Stimmstoecke Diskant Markus aus stimmstock_parameter.txt neu rechnen (ohne den ganzen Diskantteil, ca. 3 min).
# Aufruf:  python3 stimmstock_neu.py <ausgabeordner> [parameterdatei]
# Ergebnis im Ausgabeordner:
#   Stimmstock_Jakob_<n>.reihe.step          Stock + Sohle + Stimmplatten je Reihe, eingebaut (fuer Fusion)
#   Schieber_neu.step                         Schieber, Schrauben, Rastschrauben, Schraubleiste passend zur Stockhoehe
#   fraesen/Stimmstock_<reihe>_fraesen.step   aufgeklappte Fraesversion (Vorder-/Rueckwand flach, Keil)
#   fraesen/Stimmstock_<reihe>_fraesen.dxf    dazu Konturen und Taschen je Fraestiefe (Layer) fuer CamBam
import os, sys, re
aus = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else "stimmstock_aus")
if len(sys.argv) > 2: os.environ["STPARAM"] = os.path.abspath(sys.argv[2])
os.makedirs(aus, exist_ok=True)
os.environ["FRAES"] = os.path.join(aus, "fraesen")
hier = os.path.dirname(os.path.abspath(__file__)); os.chdir(hier)
src = open("diskant.py", encoding="utf-8").read()
src = src[:src.index("# ---------------- Verdeck (aussen, Z < 0)")]
sys.argv = ["diskant.py", "BESASDES", os.path.join(aus, "unbenutzt.step")]
g = {"__file__": os.path.join(hier, "diskant.py"), "__name__": "__stimmstock__"}
exec(compile(src, "diskant.py", "exec"), g)
cq = g["cq"]; teile = g["teile"]
holz, alu, stahl = cq.Color(0.85, 0.65, 0.45), cq.Color(0.80, 0.80, 0.82), cq.Color(0.6, 0.6, 0.65)
for r in ("1.reihe", "2.reihe", "3.reihe", "4.reihe"):
    a = cq.Assembly(name="Stimmstock_Jakob_" + r)
    for n, s, f in teile:
        if n.startswith("Stimmstock_" + r) or n.startswith("Stimmplatte_" + r):
            a.add(cq.Workplane().add(s), name=n, color=alu if f == "Alu" else holz)
    a.export(os.path.join(aus, f"Stimmstock_Jakob_{r}.step"))
a = cq.Assembly(name="Schieber_neu")
for n, s, f in teile:
    if re.match(r"^(Schieber_|Schraube_Schieber_|Rastschraube_|Schraubleiste_Stimmstoecke)", n):
        a.add(cq.Workplane().add(s), name=n, color=holz if f == "Holz" else stahl)
a.export(os.path.join(aus, "Schieber_neu.step"))
print("fertig:", aus)
