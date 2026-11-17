import sys, os
os.environ["NUR"]="x"
src=open('diskant.py').read().split('# ---------------- Ausgabe')[0]
sys.argv=['diskant.py','BESASDES','x.step']
exec(src)
os.makedirs("jp_brep", exist_ok=True)
for n_,s_,f_ in teile:
    if n_.startswith(("Griffstock","Hebel_","Kunststoff","Dreieck","Diskantboden","Verdeck","Schraube_Stuetze","Rahmen","Stimmstock","Stimmplatte","Rahmeneinlage")): s_.exportBrep("jp_brep/"+n_.replace("/","__")+".brep")
print("DONE")
