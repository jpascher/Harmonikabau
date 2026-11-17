import sys, os
sys.argv=['diskant.py','BESASDES',sys.argv[1]]
exec(open('diskant.py').read())
import shutil; shutil.rmtree("jp_brep", ignore_errors=True); os.makedirs("jp_brep")
for n_,s_,f_ in teile:
    if n_.startswith(("Griffstock","Hebel_","Kunststoff","Dreieck","Diskantboden","Verdeck","Schraube_Stuetze","Rahmen","Stimmstock","Stimmplatte")): s_.exportBrep("jp_brep/"+n_.replace("/","__")+".brep")
print("BREP DONE")
hoch=[(round(s_.BoundingBox().zmax,1), n_) for n_,s_,f_ in teile if s_.BoundingBox().zmax > H_RAHMEN + 10.0 + 0.01]
print("ueber Randleiste:", sorted(hoch)[-10:])
