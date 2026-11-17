# Animation der Bassmechanik Markus: jede Taste nacheinander druecken, zuletzt Lufttaste
import sys, math, pickle, numpy as np
sys.argv = ["kammern.py", "boden2.dxf", "/tmp/claude-0/-home-claude/cf042130-e5d9-5fc2-925a-aa4079250231/scratchpad/anim/tmp.step"]
G = {"__name__": "__main__"}
exec(open("kammern.py", encoding="utf-8").read(), G)
import cadquery as cq
mech, teile, klappen = G["mech"], G["teile"], G["klappen"]
achsen, stange_achse, zuordnung = G["achsen"], G["stange_achse"], G["zuordnung"]
L, Z_ACHSE, Z_GG, HUB = G["L"], G["Z_ACHSE"], G["Z_GG"], G["HUB"]
def mesh(s, tol=0.15):
    v, t = s.tessellate(tol, 0.3)
    return np.array([(p.x, p.y, p.z) for p in v], dtype=float), np.array(t, dtype=np.int64)
teile_d = dict(teile); mech_d = dict(mech); kl_d = dict(klappen)
def naechste_achse(x, y):
    best = None
    for ai, (a0, a1) in enumerate(achsen):
        f, d = G["fuss"]((x, y), (a0, a1))
        if best is None or d < best[1]: best = (ai, d)
    return best[0]
gruppen = []   # je Achse: bewegte Teile
for ai in sorted(stange_achse):
    drehen, schieben = [f"Achse_{ai+1}"], []
    for n, s in mech:
        if n.startswith(("Arm_", "Gelenkbolzen_", "Gabelkopf_", "Klappenhebel_", "Feder_Achse_")):
            if n.startswith("Feder_Achse_"):
                if n == f"Feder_Achse_{ai+1}": drehen.append(n)
                continue
            bb = s.BoundingBox()
            if n.startswith("Klappenhebel_"):
                v, _ = mesh(s, 0.5); top = v[v[:, 2] > v[:, 2].max() - 0.5]
                x_, y_ = top[:, 0].mean(), top[:, 1].mean()
            else:
                x_, y_ = (bb.xmin + bb.xmax) / 2, (bb.ymin + bb.ymax) / 2
            if naechste_achse(x_, y_) == ai: (schieben if n.startswith("Gabelkopf_") else drehen).append(n)
    kl = [e[0] for e in zuordnung.get(ai, [])]
    for k in kl: drehen += [f"{k}_Leder", f"{k}_Filz", f"{k}_Epoxy"]
    for r, X in stange_achse[ai]:
        schieben += [f"Knopf_{r}_{X:.0f}", f"Knopfzylinder_{r}_{X:.0f}", f"Schubstange_{r}_{X:.0f}"]
    gruppen.append((ai, drehen, schieben, stange_achse[ai]))
# Gabelkopf_i gehoert zur Gruppe von Arm_i (Mittelpunkt des Gabelkopfs liegt naeher an der Nachbarachse)
gruppen = [(ai, dr, [n for n in sc if not n.startswith("Gabelkopf_")], kn) for ai, dr, sc, kn in gruppen]
for ai, dr, sc, kn in gruppen:
    for n in dr:
        if n.startswith("Arm_"): sc.append("Gabelkopf_" + n.split("_")[1])
# Drehsinn: Gabelgelenk-Punkt (unter der Achse) soll nach -Y wandern; Winkel aus 6 mm Hub
h = Z_ACHSE - Z_GG; alpha = math.asin(HUB / h)
daten = {"gruppen": gruppen, "achsen": achsen, "Z_ACHSE": Z_ACHSE, "alpha": alpha, "HUB": HUB, "L": L,
         "lk": {k: G[k] for k in ("Y_AX_LK", "Z_AX_LK", "winkel", "HUB_LT")}}
alle = {}
zeigen = ["Bassboden", "Front", "Front_Aufdoppelung_2mm", "Front_Aufdoppelung_5mm", "Frontblende_dummy"]
alle["Bassboden"] = mesh(G["boden"], 0.3)
alle["Front"] = mesh(G["front"], 0.3)
for n in ("Front_Aufdoppelung_2mm", "Front_Aufdoppelung_5mm", "Zwischen_Deckel", "Lufttaste", "Lufttaste_Filz"):
    alle[n] = mesh(teile_d[n], 0.3)
for n, s in mech: alle[n] = mesh(s, 0.12)
for n, s in klappen: alle[n] = mesh(s, 0.3)
pickle.dump((alle, daten), open("anim/daten.pkl", "wb"))
print("Gruppen", len(gruppen), [ (g[0]+1, len(g[1]), len(g[2])) for g in gruppen])
