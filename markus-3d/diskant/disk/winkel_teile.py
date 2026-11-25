# 22.11.2026: Winkel-Griffstock - Kunststoffteil, Dreieckleiste und deren Schrauben entfallen; Alu-Winkel 46753, Achse Ø3,
# Schrauben durch die 5-mm-Umbiegung in den Boden
_weg = ("Kunststoffteil_Lagerblock", "Dreieckleiste", "Schraube_Kunststoff_", "Schraube_Boden_Dreieck_")
teile[:] = [t_ for t_ in teile if not t_[0].startswith(_weg)]
def _wpts(pg): return [(float(a_), float(b_)) for a_, b_ in list(pg.exterior.coords)[:-1]]
_wx0 = W_XOFF + W_SCHLITZ_X[0] - 5.0; _wx1 = W_XOFF + W_SCHLITZ_X[-1] + 5.0          # aussen je ein Hebel neben der Lasche
_wx0, _wx1 = max(_wx0, XM0 + 0.5), min(_wx1, XM1 - 0.5)          # Platte zwischen die massiven Griffstock-Enden gekuerzt
winkel = prisma_uy(_wpts(W_PROFIL), (_wx0 + _wx1) / 2, _wx1 - _wx0)
for i_ in range(W_SCHLITZE):
    winkel = winkel.cut(prisma_uy(_wpts(W_SCHLITZ), W_XOFF + W_SCHLITZ_X[i_], W_SCHLITZ_B))
# Laschen: aus der Platte ausgeklinkt und 90 Grad aufgestellt -> Dreiecksloch daneben in der Platte
_wtr = cq.Matrix([[1.0, 0.0, 0.0, 0.0], [0.0, W_ES[1], W_EN[1], W_O[1]], [0.0, -W_ES[0], -W_EN[0], -W_O[0]]])
_wla = []
for xl_ in W_LASCHEN:
    xa_ = W_XOFF + xl_
    if xa_ < _wx0 or xa_ + 2.0 > _wx1: print(f"   Winkel: Lasche bei X {xa_:.1f} ausserhalb der gekuerzten Platte"); continue
    _l = prisma_uy(_wpts(W_LASCHE), xa_ + 1.0, 2.0)
    _l = _l.cut(cq.Solid.makeCylinder(W_ACHSE_D / 2, 4.0, cq.Vector(xa_ - 1.0, W_ACHSE[1], -W_ACHSE[0]), cq.Vector(1, 0, 0)))
    winkel = winkel.fuse(_l)
    _dr = cq.Workplane("XY").polyline([(xa_ + 2.0, 4.0), (xa_ + 26.0, 4.0), (xa_ + 29.0, 6.5), (xa_ + 2.0, 29.0)]).close().extrude(W_T + 2.0).translate((0, 0, -1.0)).val()
    winkel = winkel.cut(_dr.transformGeometry(_wtr))
    _wla.append(xa_)
# Schrauben durch die Umbiegung in den Boden (Blechschraube 2,6): zwischen den Laschen und an den Enden
_wsx = [_wx0 + 6.0] + [(a_ + b_ + 2.0) / 2 for a_, b_ in zip(_wla, _wla[1:])] + [_wx1 - 6.0]
for i_, xs_ in enumerate(_wsx):
    winkel = winkel.cut(_bs(cq.Vector(xs_, W_LOCH_Y, -W_T), cq.Vector(0, 0, 1), W_T + 10.0, f"Schraube_Winkel_Boden_{i_+1}"))
teile.append(("Alu_Winkel_46753", winkel.clean(), "Alu"))
_wach = cq.Solid.makeCylinder(W_ACHSE_D / 2, _wx1 - _wx0 - 1.0, cq.Vector(_wx0 + 0.5, W_ACHSE[1], -W_ACHSE[0]), cq.Vector(1, 0, 0))
teile.append(("Achse_Winkel_3mm", _wach, "Stahl"))
print(f"Alu-Winkel 46753: Platte X {_wx0:.1f}-{_wx1:.1f} ({_wx1 - _wx0:.0f} mm), Laschen {len(_wla)}, Schlitze {W_SCHLITZE}, "
      f"Schrauben {len(_wsx)}, Volumen {winkel.Volume():.0f} mm3")
