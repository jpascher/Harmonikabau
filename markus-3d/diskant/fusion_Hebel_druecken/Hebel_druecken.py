# Schaltet einen Hebel zwischen Ruhe und gedrueckt um (Gelenk aus Hebel_Gelenke).
# Welcher Hebel: HEBEL unten aendern (z. B. "Hebel_R2_21").
import adsk.core, adsk.fusion, traceback
HEBEL = "Hebel_R1_08"

def run(context):
    app = adsk.core.Application.get(); ui = app.userInterface
    try:
        root = adsk.fusion.Design.cast(app.activeProduct).rootComponent
        j = None
        for g in list(root.asBuiltJoints) + list(root.joints):
            if g.name == "Gelenk_" + HEBEL: j = g; break
        if j is None:
            ui.messageBox("Gelenk_" + HEBEL + " nicht gefunden"); return
        m = adsk.fusion.RevoluteJointMotion.cast(j.jointMotion)
        lim = m.rotationLimits
        ziel = lim.minimumValue if abs(lim.minimumValue) > abs(lim.maximumValue) else lim.maximumValue
        m.rotationValue = 0.0 if abs(m.rotationValue) > 1e-6 else ziel
        app.activeViewport.refresh()
        import math
        ui.messageBox(f"{HEBEL}: Gelenk jetzt {math.degrees(m.rotationValue):.2f} Grad "
                      f"({'gedrueckt' if abs(m.rotationValue) > 1e-6 else 'Ruhe'}), Anschlag {math.degrees(lim.minimumValue):.2f} .. {math.degrees(lim.maximumValue):.2f} Grad")
    except:
        ui.messageBox("Fehler:\n" + traceback.format_exc())
