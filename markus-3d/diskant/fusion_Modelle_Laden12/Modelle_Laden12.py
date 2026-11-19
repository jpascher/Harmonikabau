# Laedt die STEP-Dateien der Liste je in ein neues Fusion-Dokument und speichert sie im Projekt "Stimmstock"
import adsk.core, adsk.fusion, traceback
LISTE = [(r"C:\ACAD-DWG\0-Markus-3D\0-marie\Diskantteil_Marie_v16.step", "Diskantteil_Marie_v16"),
         (r"C:\ACAD-DWG\0-Markus-3D\Diskantteil_Markus_v33.step", "Diskantteil_Markus_v33")]

def run(context):
    ui = None
    try:
        app = adsk.core.Application.get(); ui = app.userInterface
        ziel = None
        for p in app.data.dataProjects:
            if p.name == "Stimmstock": ziel = p
        if ziel is None: ziel = app.data.dataProjects.item(0)
        fertig = []
        for pfad, name in LISTE:
            doc = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
            des = adsk.fusion.Design.cast(app.activeProduct)
            im = app.importManager
            im.importToTarget(im.createSTEPImportOptions(pfad), des.rootComponent)
            app.activeViewport.fit()
            doc.saveAs(name, ziel.rootFolder, "Generator 19.11.2026: Messingbuchsen, Alu-Arme verjuengt, Rastmarken Kunststoffteil", "")
            fertig.append(name)
        ui.messageBox("Geladen und gespeichert (Projekt " + ziel.name + "):\n" + "\n".join(fertig))
    except:
        if ui: ui.messageBox(traceback.format_exc())
