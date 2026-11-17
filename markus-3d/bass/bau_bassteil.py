# Bassteil Markus - vereinfachtes 3D-Modell aus Draufsicht + Schnitt der Hauptzeichnung
# Achsen: X = Laenge (380), Y = Breite (205), Z = Hoehe ab Deckelkante (0) Richtung Balg
import sys, ezdxf, cadquery as cq

DXF = sys.argv[1]; OUT = sys.argv[2]

# ---- Parameter (aus Zeichnung gelesen bzw. angenommen -> hier korrigieren) ----
L, B        = 380.0, 205.0   # Aussenmass (GHAussenabmessungRechteck)
WAND        = 4.0            # Seitenwand (Innenmass 372 x 197)
H_RAHMEN    = 96.0           # Laengsseite im Schnitt 3420.8 -> 3516.8
BALGR_H, BALGR_B = 6.0, 13.5 # Balgrahmen
Z_BODEN, T_BODEN = 24.4, 6.0 # Boden (6 mm laut Fraesdaten, Lage angenommen)
Z_KAM_OBEN  = 81.7           # Oberkante Basskammern / Unterkante Stimmplatte
T_STPL      = 3.0            # Bassstimmplatte
KAM_Y0, KAM_Y1 = 11.0, 105.0 # Basskammern quer (Schnitt)
KAM_X0, KAM_X1 = 14.0, 366.0
KAM_WAND    = 5.0
TRENN_X     = [64.5, 117.5, 193.5, 292.5]   # Trennwaende (angenommen: zwischen den Klappenoeffnungen)
BEGL_Y0, BEGL_Y1 = 107.5, 142.0
BEGL_X0, BEGL_X1 = 70.0, 340.0
Z_BEGL_STST = 61.5           # Begleiter-Stimmstock Grundplatte
BEGL_FIN_H  = 42.5
OX, OY = 2924.3, 1521.7      # Nullpunkt der Draufsicht in der DXF

top = Z_BODEN + T_BODEN
d = ezdxf.readfile(DXF); m = d.modelspace()

def box(x0, x1, y0, y1, z0, z1):
    return cq.Workplane("XY").box(x1-x0, y1-y0, z1-z0, centered=False).translate((x0, y0, z0))

# ---- Rahmenteile aus Markus-Bass-Rahmen.dxf (Front 382 lang, Seitenteil 207, beide 76 hoch) ----
RD = ezdxf.readfile(sys.argv[3]).modelspace()
H_R, Y_TOP = 76.0, 71.8          # Rahmenhoehe; Oberkante in der Rahmenzeichnung = Deckelseite (Z=0)
H_BLENDE = 20.0                  # 20-mm-Blende zwischen Rahmen und Balgrahmen
def ausschnitte(x0, x1, mit_bohrungen):
    """2D-Ausschnitte (u, v) eines Teils: u entlang der Wand, v = Z"""
    out = []
    for p in RD.query('LWPOLYLINE[layer=="Layer2"]'):            # Ornament-Durchbrueche
        pts = list(p.get_points("xy"))
        u0, u1 = min(a for a, _ in pts) - x0, x1 - max(a for a, _ in pts)
        if u0 >= 0 and u1 >= 0 and min(u0, u1) < 40:          # nur Ornamente an den Enden des Teils
            out.append(("poly", [(a-x0, Y_TOP-b) for a, b in pts]))
    if mit_bohrungen:
        for lay in ("Hörner-Bohrung", "bohrungen-bassdrücker"):
            for c in RD.query(f'CIRCLE[layer=="{lay}"]'):
                if -4.2 <= c.dxf.center.y <= Y_TOP:
                    out.append(("kreis", (c.dxf.center.x-x0, Y_TOP-c.dxf.center.y, c.dxf.radius)))
    return out
def wandplatte(laenge, cuts, achse):
    """Platte laenge x H_R x WAND; achse 'x' = Laengswand, 'y' = Seitenwand"""
    w = cq.Workplane("XZ").rect(laenge, H_R, centered=False).extrude(-WAND)   # u->X, v->Z, Dicke +Y
    for typ, g in cuts:
        if typ == "poly":
            t = cq.Workplane("XZ").polyline(g).close().extrude(-WAND-2).translate((0, -1, 0))
        else:
            t = cq.Workplane("XZ").center(g[0], g[1]).circle(g[2]).extrude(-WAND-2).translate((0, -1, 0))
        w = w.cut(t)
    if achse == "y":
        w = w.rotate((0, 0, 0), (0, 0, 1), 90).translate((WAND, 0, 0))
    return w
front = wandplatte(L, ausschnitte(-382.0, 0.0, True), "x")                      # Y = 0
rueck = wandplatte(L, ausschnitte(-382.0, 0.0, False), "x").translate((0, B-WAND, 0))
seite = ausschnitte(-207.0, 0.0, False)
sl = wandplatte(B, seite, "y")                                                   # X = 0
sr = wandplatte(B, seite, "y").translate((L-WAND, 0, 0))
rahmen = front.union(rueck).union(sl).union(sr)
blende = box(0, L, 0, B, H_R, H_R+H_BLENDE).cut(box(WAND, L-WAND, WAND, B-WAND, H_R-1, H_R+H_BLENDE+1))
balgr  = box(0, L, 0, B, H_RAHMEN, H_RAHMEN+BALGR_H).cut(
         box(BALGR_B, L-BALGR_B, BALGR_B, B-BALGR_B, H_RAHMEN-1, H_RAHMEN+BALGR_H+1))

# Boden mit echten Oeffnungen aus der Draufsicht
boden = box(WAND, L-WAND, WAND, B-WAND, Z_BODEN, top)
n_loch = 0
for lay in ["Füllung-DH-Öffnungen2", "Begleiterstimmstock-fuellung-öffnungen"]:
    for c in m.query(f'CIRCLE[layer=="{lay}"]'):
        x, y = c.dxf.center.x-OX, c.dxf.center.y-OY
        boden = boden.cut(cq.Workplane("XY").circle(c.dxf.radius).extrude(T_BODEN+2)
                          .translate((x, y, Z_BODEN-1))); n_loch += 1
for lay in ["2Füllung-Öffnug", "Füllung-luftklappe"]:
    for p in m.query(f'LWPOLYLINE[layer=="{lay}"]'):
        pts = [(q[0]-OX, q[1]-OY) for q in p.get_points("xy")]
        if not (0 < min(a for a, _ in pts) < L): continue          # nur Teile in der Draufsicht
        boden = boden.cut(cq.Workplane("XY").polyline(pts).close().extrude(T_BODEN+2)
                          .translate((0, 0, Z_BODEN-1))); n_loch += 1

# Basskammern mit Trennwaenden
kam = box(KAM_X0, KAM_X1, KAM_Y0, KAM_Y1, top, Z_KAM_OBEN).cut(
      box(KAM_X0+KAM_WAND, KAM_X1-KAM_WAND, KAM_Y0+KAM_WAND, KAM_Y1-KAM_WAND, top-1, Z_KAM_OBEN+1))
for tx in TRENN_X:
    kam = kam.union(box(tx-KAM_WAND/2, tx+KAM_WAND/2, KAM_Y0, KAM_Y1, top, Z_KAM_OBEN))
stpl = box(KAM_X0, KAM_X1, KAM_Y0, KAM_Y1, Z_KAM_OBEN, Z_KAM_OBEN+T_STPL)

# Begleiter: Zwischenteil + Stimmstock (Grundplatte + 3 Stege laut Schnitt)
bzw = box(BEGL_X0, BEGL_X1, BEGL_Y0, BEGL_Y1, top, Z_BEGL_STST).cut(
      box(BEGL_X0+KAM_WAND, BEGL_X1-KAM_WAND, BEGL_Y0+KAM_WAND, BEGL_Y1-KAM_WAND, top-1, Z_BEGL_STST+1))
bst = box(BEGL_X0, BEGL_X1, BEGL_Y0, BEGL_Y1, Z_BEGL_STST, Z_BEGL_STST+3.6)
for y0 in (BEGL_Y0+5, BEGL_Y0+16.4, BEGL_Y0+26.4):
    bst = bst.union(box(BEGL_X0, BEGL_X1, y0, y0+3.6, Z_BEGL_STST+3.6, Z_BEGL_STST+3.6+BEGL_FIN_H))

# Luftklappen-Trennplatte (oberer Bereich im Schnitt)
lkp = box(WAND, L-WAND, 147.0, B-WAND, Z_KAM_OBEN, Z_KAM_OBEN+3)

teile = [("Front", front, (0.80,0.62,0.40)), ("Rueckteil", rueck, (0.78,0.60,0.38)),
         ("Seitenteil_links", sl, (0.72,0.55,0.35)), ("Seitenteil_rechts", sr, (0.72,0.55,0.35)),
         ("Blende_20mm", blende, (0.62,0.45,0.28)),("Balgrahmen", balgr, (0.55,0.35,0.20)),
         ("Boden", boden, (0.85,0.75,0.55)), ("Basskammern", kam, (0.60,0.70,0.85)),
         ("Bassstimmplatte", stpl, (0.75,0.75,0.78)), ("Begleiter-Zwischenteil", bzw, (0.65,0.80,0.60)),
         ("Begleiter-Stimmstock", bst, (0.45,0.65,0.45)), ("Luftklappen-Trennplatte", lkp, (0.85,0.55,0.55))]
asm = cq.Assembly(name="Bassteil_Markus")
for n, s, col in teile:
    asm.add(s, name=n, color=cq.Color(*col))
asm.save(f"{OUT}/Bassteil_Markus.step")
for n, s, _ in teile:
    cq.exporters.export(s, f"{OUT}/stl_{n}.stl")
print("Bohrungen/Oeffnungen im Boden:", n_loch)
