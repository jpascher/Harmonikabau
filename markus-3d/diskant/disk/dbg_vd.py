import sys
src=open('diskant.py').read().split('# ---------------- Ausgabe')[0]
sys.argv=['diskant.py','BESASDES','x.step']
exec(src)
import pickle
d={n_:s_ for n_,s_,f_ in teile if f_=="Verdeck" or n_ in ("Griffstock_Koerper","Griffstock_Abdeckung","Rahmen_vorne","Rahmen_hinten")}
for n_,s_ in d.items(): s_.exportBrep(f"dbg_{n_}.brep")
print("DONE")
