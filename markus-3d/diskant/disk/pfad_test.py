import math, ezdxf, cadquery as cq, json, sys
from ezdxf import bbox
BXO=380.0
bm=ezdxf.readfile("boden.dxf").modelspace()
src=open("hebel.py").read()
exec(src.split("T_ALU, T_HOLZ")[0])
FUEHR=eval(src.split("FUEHR = ")[1].split("\nPFADE")[0])
PF=[]
for r,x,k,kr in zuordnung:
    up,Yp=lager_von(r); v=k[0]
    if kr!=r: PF.append((r,x,[(Yp,x),(32,x),(44,v),(k[1],v)])); continue
    f=FUEHR[r]; lane=v+f["lane"]; ya,yb=f["y_lane"]; y0=min(p[1] for p in ENDE[r])
    if f["y_fin"]:
        yf0,yf1=f["y_fin"]; PF.append((r,x,[(Yp,x),(ya,x),(yb,lane),(yf0,lane),(yf1,v),(k[1],v)]))
    else:
        ye=y0+0.5; PF.append((r,x,[(Yp,x),(ya,x),(yb,lane),(ye,v),(k[1],v)]))
def xb(pf,y):
    for (a,xa),(b,xb_) in zip(pf,pf[1:]):
        if a<=y<=b: return xa+(xb_-xa)*(y-a)/max(b-a,1e-9)
for yy in [30+0.5*i for i in range(300)]:
    xs=sorted((xb(pf,yy),r,x) for r,x,pf in PF if r!=4 and xb(pf,yy) is not None)
    for a,b in zip(xs,xs[1:]):
        if b[0]-a[0]<3.5: print("Y %.1f  %s  %s  Abstand %.1f"%(yy,(a[1],a[2],round(a[0],1)),(b[1],b[2],round(b[0],1)),b[0]-a[0]))
