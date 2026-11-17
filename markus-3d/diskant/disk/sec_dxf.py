import sys, os, math, pickle
os.environ["NUR"] = "Diskantboden|Griffstock|Kunststoff|Dreieck|Hebel|Verstaerkung|Verdeck|Rahmen|Schraube"
src = open('diskant.py').read().split('# ---------------- Ausgabe')[0]
sys.argv = ['diskant.py', 'BESASDES', 'x.step']
exec(src)
import cadquery as cq
D_ = {n_: s_ for n_, s_, f_ in teile}
pickle.dump({n: s.exportBrep(f"sec_brep/{n.replace('/','__')}.brep") or f"sec_brep/{n.replace('/','__')}.brep" for n, s in D_.items()} if os.makedirs("sec_brep", exist_ok=True) is None else 0, open("sec_brep/index.pkl", "wb"))
print("DONE", len(D_))
