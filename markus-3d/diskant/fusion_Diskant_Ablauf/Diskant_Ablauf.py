# Ablauf im Fusion-Fenster (fuer Bildschirmaufnahme):
# 1. Verdeck nach unten wegfahren   2. Griffstock ausblenden
# 3. Hebel Reihe 1..4 druecken und loslassen   4. alle Hebel zusammen   5. alles zurueck
# Benoetigt die Gelenke aus "Hebel_Gelenke". Ansicht vorher selbst einstellen.
import adsk.core, adsk.fusion, traceback, time

SCHRITTE = 4          # Zwischenschritte je Bewegung (mehr = fluessiger, aber langsamer)
PAUSE = 0.05          # Sekunden zwischen den Bildern
VERDECK_WEG = 6.0     # cm nach unten (-Z)

def zeigen(app):
    app.activeViewport.refresh(); adsk.doEvents(); time.sleep(PAUSE)

def run(context):
    app = adsk.core.Application.get(); ui = app.userInterface
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        root = design.rootComponent
        occs = list(root.occurrences)
        verdeck = [o for o in occs if o.name.startswith("Verdeck:") or o.name.split(":")[0] == "Verdeck"]
        gs = [o for o in occs if o.name.startswith("Griffstock")]
        gelenke = {}
        for j in list(root.asBuiltJoints):
            if j.name.startswith("Gelenk_Hebel_R"):
                r = int(j.name[len("Gelenk_Hebel_R")])
                gelenke.setdefault(r, []).append(j)
        if not gelenke:
            ui.messageBox("Keine Hebel-Gelenke gefunden - zuerst das Skript Hebel_Gelenke ausfuehren."); return
        def ziel(j):
            lim = adsk.fusion.RevoluteJointMotion.cast(j.jointMotion).rotationLimits
            return lim.minimumValue if abs(lim.minimumValue) > abs(lim.maximumValue) else lim.maximumValue
        def setze(js, f):
            for j in js: adsk.fusion.RevoluteJointMotion.cast(j.jointMotion).rotationValue = f * ziel(j)
            zeigen(app)
        def schiebe(os_, start, dz, f):
            for o, m0 in zip(os_, start):
                m = m0.copy(); t = m.translation; t.z = t.z + f * dz; m.translation = t
                o.transform = m
            zeigen(app)
        v0 = [o.transform.copy() for o in verdeck]
        # 1. Verdeck weg
        for i in range(1, SCHRITTE + 1): schiebe(verdeck, v0, -VERDECK_WEG, i / SCHRITTE)
        # 2. Griffstock aus
        for o in gs: o.isLightBulbOn = False
        zeigen(app); time.sleep(0.5)
        # 3. reihenweise
        for r in sorted(gelenke):
            js = gelenke[r]
            for i in range(1, SCHRITTE + 1): setze(js, i / SCHRITTE)
            time.sleep(0.3)
            for i in range(SCHRITTE - 1, -1, -1): setze(js, i / SCHRITTE)
        # 4. alle zusammen
        alle = [j for r in gelenke for j in gelenke[r]]
        for i in range(1, SCHRITTE + 1): setze(alle, i / SCHRITTE)
        time.sleep(0.5)
        for i in range(SCHRITTE - 1, -1, -1): setze(alle, i / SCHRITTE)
        # 5. zurueck
        for o in gs: o.isLightBulbOn = True
        zeigen(app)
        for i in range(SCHRITTE - 1, -1, -1): schiebe(verdeck, v0, -VERDECK_WEG, i / SCHRITTE)
        for o, m0 in zip(verdeck, v0): o.transform = m0
        zeigen(app)
    except:
        ui.messageBox("Fehler:\n" + traceback.format_exc())
