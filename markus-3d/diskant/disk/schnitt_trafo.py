# Transformation Griffstock-Schnitt (griffstock-alternativ1-schnitt.dxf) -> Modell (u nach aussen, Y)
import ezdxf, ezdxf.path, math
from ezdxf import bbox
SM = ezdxf.readfile("gschnitt.dxf").modelspace()
P_A, P_B = (-280.56, -10.2), (-141.66, 16.8)          # Bodenaussenseite im Schnitt (links weit weg, rechts Griffstockseite)
_L = math.hypot(P_B[0]-P_A[0], P_B[1]-P_A[1])
D = ((P_B[0]-P_A[0])/_L, (P_B[1]-P_A[1])/_L)
N = (-D[1], D[0])                                       # nach aussen (Klappenseite)
def kl_mitten():
    out = []
    for e in SM.query('LWPOLYLINE[layer=="Klappen"]'):
        pts = [(x, y) for x, y, *_ in e.get_points()]
        out.append((sum(p[0] for p in pts)/len(pts), sum(p[1] for p in pts)/len(pts)))
    return sorted(out)
def su(p):
    v = (p[0]-P_B[0], p[1]-P_B[1]); return v[0]*D[0]+v[1]*D[1], v[0]*N[0]+v[1]*N[1]
KL_S = [su(k)[0] for k in kl_mitten()]                 # von weit weg (Reihe 1) bis nahe (Reihe 4)
MODELL_Y = {1: 169.5, 2: 134.8, 3: 99.7, 4: 64.0}       # Klappenreihen-Mitten im Modell
# Y = Y0 - s ; Y0 so, dass die Klappenreihen im Mittel passen
Y0 = sum(MODELL_Y[i+1] + KL_S[i] for i in range(4)) / 4
def modell(p):
    s, u = su(p); return u, Y0 - s
if __name__ == "__main__":
    print("Richtung Boden %.2f Grad, Y0 %.2f" % (math.degrees(math.atan2(D[1], D[0])), Y0))
    for i, s in enumerate(KL_S): print(" Klappenreihe", i+1, "Schnitt-Y %.1f  Modell %.1f" % (Y0 - s, MODELL_Y[i+1]))
    for e in SM.query('CIRCLE[layer=="Hebel-4r"]'):
        u, Y = modell((e.dxf.center.x, e.dxf.center.y)); print(" Lager u %.1f Y %.1f" % (u, Y))
    for e in SM.query('INSERT'):
        if e.dxf.insert.x < -80:
            u, Y = modell((e.dxf.insert.x, e.dxf.insert.y)); print(" Knopf u %.1f Y %.1f" % (u, Y))
