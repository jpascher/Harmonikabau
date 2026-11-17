import sys
sys.argv=['diskant.py','BESASDES','x.step']
src=open('diskant.py').read().split('# Schraubleiste an der Wand der hohen Seite')[0]
exec(src)
zs=[(s_.BoundingBox().zmax, n_) for n_, s_, f_ in teile if n_.startswith(("Stimmstock","Stimmplatte"))]
print("max Z", max(zs)); print("Rahmen", [ (n_, round(s_.BoundingBox().zmax,1)) for n_,s_,f_ in teile if n_.startswith("Rahmen")])
