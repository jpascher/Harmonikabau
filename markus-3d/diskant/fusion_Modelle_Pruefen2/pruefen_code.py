# 23.11.2026: STEP-Modelle der Reihe nach laden, pruefen (Komponenten/Koerper, offene Koerper, Bild), wieder schliessen ohne Speichern
import adsk.core, adsk.fusion, traceback, os, time
ZIEL = r"C:\ACAD-DWG\0-Markus-3D\Pruefung"
# Liste wird bei jedem Lauf frisch gelesen (Fusion behaelt das Skript im Speicher): eine STEP-Datei je Zeile
LOG = ZIEL + r"\pruefung_log.txt"

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get(); ui = app.userInterface
        os.makedirs(ZIEL, exist_ok=True)
        LISTE = [z.strip() for z in open(ZIEL + r"\liste.txt", encoding="utf-8") if z.strip() and not z.startswith("#")]
        zu = []                                   # frueher importierte, nie gespeicherte Modelle schliessen (nur Ansicht)
        for i in range(app.documents.count - 1, -1, -1):
            d = app.documents.item(i)
            try: gesp = d.isSaved and d.dataFile is not None
            except: gesp = False
            if not gesp and (d.name.startswith("Diskantteil") or d.name.startswith("Einzelhebel")):
                zu.append(d.name); d.close(False)
        with open(LOG, "w", encoding="utf-8") as f: f.write("Pruefung %s, geschlossen: %s\n" % (time.strftime("%d.%m.%Y %H:%M"), ", ".join(zu)))
        for pfad in LISTE:
            name = os.path.splitext(os.path.basename(pfad))[0]
            t0 = time.time(); zeile = name + ": "
            try:
                doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
                des = adsk.fusion.Design.cast(app.activeProduct)
                im = app.importManager
                im.importToTarget(im.createSTEPImportOptions(pfad), des.rootComponent)
                adsk.doEvents()
                root = des.rootComponent
                n_occ = root.allOccurrences.count
                n_k, n_offen, n_leer = 0, 0, 0
                for comp in des.allComponents:
                    for b in comp.bRepBodies:
                        n_k += 1
                        if not b.isSolid: n_offen += 1
                        if b.faces.count == 0: n_leer += 1
                vp = app.activeViewport; cam = vp.camera
                cam.viewOrientation = adsk.core.ViewOrientations.IsoTopRightViewOrientation
                cam.isFitView = True; vp.camera = cam; adsk.doEvents(); vp.fit(); adsk.doEvents()
                vp.saveAsImageFile(ZIEL + "\\" + name + "_iso.png", 1600, 1000)
                cam = vp.camera; cam.viewOrientation = adsk.core.ViewOrientations.FrontViewOrientation
                cam.isFitView = True; vp.camera = cam; adsk.doEvents(); vp.fit(); adsk.doEvents()
                vp.saveAsImageFile(ZIEL + "\\" + name + "_vorne.png", 1600, 1000)
                zeile += "Komponenten %d, Koerper %d, nicht geschlossen %d, leer %d, %.0f s" % (n_occ, n_k, n_offen, n_leer, time.time() - t0)
                doc.close(False)
            except:
                zeile += "FEHLER " + traceback.format_exc().replace("\n", " | ")
            with open(LOG, "a", encoding="utf-8") as f: f.write(zeile + "\n")
        with open(LOG, "a", encoding="utf-8") as f: f.write("ENDE\n")
    except:
        if ui: ui.messageBox(traceback.format_exc())
