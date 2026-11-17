# Drehgelenke fuer alle Hebel des Diskantteils (Diskantteil_Markus_v21.step).
# Je Hebel-Komponente "Hebel_Rr_nn": wie-gebaut-Drehgelenk an der Zylinderflaeche der Achse,
# Gegenstueck = Kunststoffteil (wird fixiert). Anschlag: Ruhe (0) bis gedrueckt (Knopf 3,5 mm).
import adsk.core, adsk.fusion, traceback, math, os

LOG = r"C:\ACAD-DWG\0-Markus-3D\Hebel_Gelenke_log.txt"

def alle_koerper(occ):
    """Koerper der Hebel-Komponente samt Unterkomponenten, jeweils mit Namen (Teil, Koerper)."""
    out = [(occ.name, b) for b in occ.bRepBodies]
    for sub in occ.childOccurrences:
        out += [(sub.name + "/" + n, b) for n, b in alle_koerper(sub)]
    return out

def run(context):
    app = adsk.core.Application.get(); ui = app.userInterface
    log = []
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        root = design.rootComponent
        boden = None
        for occ in root.allOccurrences:
            if occ.name.startswith("Kunststoffteil"):
                boden = occ; break
        if boden is None:
            ui.messageBox("Kunststoffteil nicht gefunden"); return
        boden.isGrounded = True
        hebel = [o for o in root.occurrences if o.name.startswith("Hebel_R")]
        ok = 0
        for occ in hebel:
            try:
                kb = alle_koerper(occ)
                achse = [b for n, b in kb if "Achse" in n or "Achse" in b.name]
                knopf = [b for n, b in kb if "Knopf" in n or "Knopf" in b.name]
                if not achse or not knopf:
                    log.append(f"{occ.name}: Achse/Knopf nicht gefunden"); continue
                zyl = [f for f in achse[0].faces if f.geometry.surfaceType == adsk.core.SurfaceTypes.CylinderSurfaceType]
                fl = max(zyl, key=lambda f: f.area)
                cyl = adsk.core.Cylinder.cast(fl.geometry)
                ax = cyl.axis; org = cyl.origin
                kc = knopf[0].boundingBox
                kx = (kc.minPoint.y + kc.maxPoint.y) / 2 - org.y
                kz = (kc.minPoint.z + kc.maxPoint.z) / 2 - org.z
                th = 0.35 / math.hypot(kx, kz)                     # 3,5 mm am Knopf (Fusion intern cm)
                kreis = [e for e in fl.edges if e.geometry.curveType == adsk.core.Curve3DTypes.Circle3DCurveType]
                geo = adsk.fusion.JointGeometry.createByCurve(kreis[0], adsk.fusion.JointKeyPointTypes.CenterKeyPoint)
                ax = adsk.core.Circle3D.cast(kreis[0].geometry).normal
                inp = root.asBuiltJoints.createInput(occ, boden, geo)
                inp.setAsRevoluteJointMotion(adsk.fusion.JointDirections.ZAxisJointDirection)
                j = root.asBuiltJoints.add(inp)
                j.name = "Gelenk_" + occ.name.split(":")[0]
                lim = adsk.fusion.RevoluteJointMotion.cast(j.jointMotion).rotationLimits
                # gedrueckt = Drehung um -th um +X (Modell); Achsrichtung der Flaeche beachten
                lo, hi = (-th, 0.0) if ax.x > 0 else (0.0, th)
                lim.isMinimumValueEnabled = True; lim.minimumValue = lo
                lim.isMaximumValueEnabled = True; lim.maximumValue = hi
                ok += 1
                log.append(f"{occ.name}: ok, Weg {math.degrees(th):.2f} Grad")
            except:
                log.append(f"{occ.name}: FEHLER {traceback.format_exc().splitlines()[-1]}")
        try:
            with open(LOG, "w", encoding="utf-8") as f: f.write("\n".join(log))
        except: pass
        ui.messageBox(f"Drehgelenke angelegt: {ok} von {len(hebel)} Hebeln.\nKunststoffteil ist fixiert.\nProtokoll: {LOG}")
    except:
        ui.messageBox("Fehler:\n" + traceback.format_exc())
