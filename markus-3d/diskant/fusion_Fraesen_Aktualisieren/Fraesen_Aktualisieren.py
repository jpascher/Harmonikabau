# Fraesen_Aktualisieren (Johann 19.11.2026) - ein Fusion-Skript fuer alle Instrumente, ohne feste Teilenamen.
# Liest jede C:\ACAD-DWG\0-Markus-3D\Fraesen\<Name>\fraes_liste.json (von erzeuge.py geschrieben) und macht im Dokument
# "Stimmstoecke_<Name>_Fraesen" (wird bei Bedarf geoeffnet oder neu angelegt):
#   1. alle Teile ersetzen: Stimmstoecke laut Liste (untereinander), Griffstock Teil 1/2 daneben
#   2. CAM neu: Stimmstoecke wie Jakob (je Stock Vorder-/Rueckseite, 2x adaptiv, 2,5 mm, Grenze Bodenkontur, ausserhalb,
#      Rohteil seitlich 0); Griffstock je Teil "oben"/"Deckelseite" (6 mm Schruppen/Schlichten der Stufen,
#      Kanaele einzeln 3,175/2,5 mm, Innenraum-Rechteck ausschneiden, Deckel-Ausnehmung vollflaechig)
#   3. speichern - Werkzeugwege werden NICHT berechnet.
import adsk.core, adsk.fusion, adsk.cam, traceback, json, os
BASIS = "C:\\ACAD-DWG\\0-Markus-3D\\Fraesen\\"
PROJEKT = "Stimmstock"
VORBILD = "Stimmstöcke-1-reihe-jakob"
LOG = "C:\\ACAD-DWG\\0-Markus-3D\\fraesen_aktualisieren_log.txt"
SETUP_PAR = ["job_stockMode", "job_stockOffsetMode", "job_stockOffsetSides", "job_stockOffsetTop", "job_stockOffsetBottom"]
OP_SKIP = ("stock", "surface", "connectionMove", "milling_connection", "machiningBoundary", "restMaterial", "view_",
           "tool_", "holder_", "operation_", "boundary", "useSilhouette", "silhouette", "minimumSilhouette")
WERKZEUGE = {"6": dict(d=6.0, lcf=25.0, body=28.0, oal=60.0, nr=6, txt="6 mm Schaftfraeser"),
             "3.175": dict(d=3.175, lcf=22.0, body=24.0, oal=50.0, nr=3, txt="3,175 mm Schaftfraeser"),
             "2.5": None}
log = []

def _suche(folder, name, tiefe=0):
    for i in range(folder.dataFiles.count):
        if folder.dataFiles.item(i).name == name: return folder.dataFiles.item(i)
    if tiefe < 4:
        for i in range(folder.dataFolders.count):
            r = _suche(folder.dataFolders.item(i), name, tiefe + 1)
            if r: return r
    return None

def _p(ps, n, e):
    try:
        q = ps.itemByName(n)
        if q is None: return False
        q.expression = e; return True
    except: return False

def _alle_koerper(o):
    r = [(b, o.component.name, o.name) for b in o.bRepBodies]
    for c in o.childOccurrences: r += _alle_koerper(c)
    return r

def _bodenflaeche(body):
    best = None
    for f in body.faces:
        try:
            pl = adsk.core.Plane.cast(f.geometry)
            if pl is None or abs(pl.normal.z) < 0.99: continue
            z = f.pointOnFace.z
            if best is None or z < best[0] - 1e-4 or (abs(z - best[0]) < 1e-4 and f.area > best[1].area): best = (z, f)
        except: pass
    return best[1] if best else None

def _werkzeug(basis_json, w):
    t = adsk.cam.Tool.createFromJson(basis_json)
    if w is None: return t
    for n, e in (("tool_diameter", "%g mm" % w["d"]), ("tool_shaftDiameter", "%g mm" % w["d"]), ("tool_fluteLength", "%g mm" % w["lcf"]),
                 ("tool_bodyLength", "%g mm" % w["body"]), ("tool_overallLength", "%g mm" % w["oal"]), ("tool_number", str(w["nr"])),
                 ("tool_description", "'%s'" % w["txt"])):
        if not _p(t.parameters, n, e): log.append("Werkzeug: %s nicht gesetzt" % n)
    return t

def _skizze(comp, name, z_mm=0.0):
    pl = comp.xYConstructionPlane
    if abs(z_mm) > 1e-6:
        pi = comp.constructionPlanes.createInput(); pi.setByOffset(comp.xYConstructionPlane, adsk.core.ValueInput.createByReal(z_mm / 10.0))
        pl = comp.constructionPlanes.add(pi); pl.name = name + "_Ebene"
    sk = comp.sketches.add(pl); sk.name = name
    return sk

def _rechteck(sk, r):
    x0, y0, x1, y1 = [v / 10.0 for v in r]
    return [sk.sketchCurves.sketchLines.addTwoPointRectangle(adsk.core.Point3D.create(x0, y0, 0), adsk.core.Point3D.create(x1, y1, 0)).item(0)]

def _kette(param, geo):
    """geo = Liste von Kurven (eine Kette) oder Liste von Listen (mehrere Ketten, z. B. mehrere Rechtecke)"""
    cc = param.value; sels = cc.getCurveSelections(); sels.clear()
    for g_ in (geo if geo and isinstance(geo[0], list) else [geo]):
        ch = sels.createNewChainSelection(); ch.inputGeometry = g_
    cc.applyCurveSelections(sels)

def _op(setup, strat, name, tool, par, grenze=None, kontur=None, enthalten="'inside'"):
    oi = setup.operations.createInput(strat); oi.tool = tool; oi.displayName = name
    op = setup.operations.add(oi)
    for n, e in par.items():
        if not _p(op.parameters, n, e): log.append("  %s: Parameter %s nicht gesetzt" % (name, n))
    if grenze is not None:
        _p(op.parameters, "boundaryMode", "'selection'"); _p(op.parameters, "boundaryContainment", enthalten)
        _kette(op.parameters.itemByName("machiningBoundarySel"), grenze)
    if kontur is not None:
        _kette(op.parameters.itemByName("contours"), kontur)
    return op

def _flaechen(s, gruppen, wz):
    """je Hoehe ein Schlichtlauf (Flaechen-Schlichten), Grenze = Rechtecke dieser Hoehe, seitlich verlaengert"""
    for z, rechtecke in gruppen:
        try: _op(s, "horizontal", "Flaeche Z %.2f schlichten 0,5" % z, wz, {"useStockToLeave": "false"}, grenze=rechtecke)
        except:
            log.append(s.name + ": Flaechen-Schlichten nicht moeglich, adaptiv: " + traceback.format_exc().splitlines()[-1])
            _op(s, "adaptive", "Flaeche Z %.2f schlichten 0,5" % z, wz, {"useStockToLeave": "false", "machineCavities": "true"}, grenze=rechtecke)

def _jakob(app):
    df = None
    for j in range(app.data.dataProjects.count):
        df = _suche(app.data.dataProjects.item(j).rootFolder, VORBILD)
        if df: break
    vdoc = app.documents.open(df, True); adsk.doEvents()
    vcam = adsk.cam.CAM.cast(vdoc.products.itemByProductType("CAMProductType")); vs = vcam.setups.item(0)
    spar = {n: vs.parameters.itemByName(n).expression for n in SETUP_PAR}
    vops = [vs.allOperations.item(i) for i in range(vs.allOperations.count)]
    opar = []
    for o in vops:
        d = {}
        for i in range(o.parameters.count):
            q = o.parameters.item(i)
            if q.name.startswith(OP_SKIP): continue
            try:
                if q.isEditable: d[q.name] = q.expression
            except: pass
        opar.append((o.strategy, d))
    wj = vops[0].tool.toJson()
    feeds = {n: vops[0].parameters.itemByName(n).expression for n in ("tool_spindleSpeed", "tool_feedCutting", "tool_feedPlunge", "tool_feedEntry", "tool_feedExit", "tool_feedRamp")}
    vdoc.close(False); adsk.doEvents()
    return spar, opar, wj, feeds

def _dokument(app, name):
    for i in range(app.documents.count):
        if app.documents.item(i).name.startswith(name): return app.documents.item(i)
    for j in range(app.data.dataProjects.count):
        df = _suche(app.data.dataProjects.item(j).rootFolder, name)
        if df: d = app.documents.open(df, True); adsk.doEvents(); return d
    d = app.documents.add(adsk.core.DocumentTypes.FusionDesignDocumentType)
    ziel = [p for p in app.data.dataProjects if p.name == PROJEKT] or [app.data.dataProjects.item(0)]
    d.saveAs(name, ziel[0].rootFolder, "neu angelegt von Fraesen_Aktualisieren", "")
    return d

def _ein_instrument(app, ui, ordner, L, jakob):
    spar, opar, wj, feeds = jakob
    doc = _dokument(app, L["dokument"]); doc.activate(); adsk.doEvents()
    ui.workspaces.itemById("FusionSolidEnvironment").activate(); adsk.doEvents()
    des = adsk.fusion.Design.cast(doc.products.itemByProductType("DesignProductType")); root = des.rootComponent
    for i in range(root.occurrences.count - 1, -1, -1): root.occurrences.item(i).deleteMe()
    im = app.importManager
    st_occ = []
    for i, t in enumerate(L["stimmstoecke"]["teile"]):
        occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create()); occ.component.name = t
        im.importToTarget(im.createSTEPImportOptions(os.path.join(ordner, L["stimmstoecke"]["ordner"], t + ".step")), occ.component)
        m = adsk.core.Matrix3D.create(); m.translation = adsk.core.Vector3D.create(0, -12.0 * i, 0); occ.transform2 = m
        st_occ.append(occ)
    G = L["griffstock"]; J = json.load(open(os.path.join(ordner, G["ordner"], G["json"]), encoding="utf-8"))
    gs_occ = {}
    for k, (tn, datei) in enumerate(sorted(G["teile"].items())):
        occ = root.occurrences.addNewComponent(adsk.core.Matrix3D.create()); occ.component.name = datei
        im.importToTarget(im.createSTEPImportOptions(os.path.join(ordner, G["ordner"], datei + ".step")), occ.component)
        m = adsk.core.Matrix3D.create(); m.translation = adsk.core.Vector3D.create(0, 15.0 + 15.0 * k, 0); occ.transform2 = m
        gs_occ[tn] = occ
    if des.designType == adsk.fusion.DesignTypes.ParametricDesignType and des.snapshots.hasPendingSnapshot: des.snapshots.add()
    geo = {}
    for tn, occ in gs_occ.items():
        T = J["teile"][tn]; comp = occ.component
        cx = lambda Lc, occ=occ: [c.createForAssemblyContext(occ) for c in Lc]
        sk = _skizze(comp, "Grenzen_" + tn)
        g = {"flaechen": {sd: [(fl["z"], [cx(_rechteck(sk, r_)) for r_ in fl["r"]]) for fl in T["flaechen"][sd]] for sd in ("oben", "unten")},
             "kanaele": [(i + 1, k["tiefe_von_oben"], cx(_rechteck(sk, k["r"]))) for i, k in enumerate(T["kanaele"])]}
        sk2 = _skizze(comp, "Innenraum_" + tn, -T["dicke"] - 0.3)
        g["innen"] = cx(_rechteck(sk2, T["innenraum"]))
        geo[tn] = g
    ui.workspaces.itemById("CAMEnvironment").activate(); adsk.doEvents()
    cam = adsk.cam.CAM.cast(doc.products.itemByProductType("CAMProductType"))
    for i in range(cam.setups.count - 1, -1, -1):
        if cam.setups.item(i).name.startswith(("Stimmstock ", "Griffstock ")): cam.setups.item(i).deleteMe()
    W = {k: _werkzeug(wj, v) for k, v in WERKZEUGE.items()}
    # --- Stimmstoecke (wie Jakob) ---
    for occ in st_occ:
        koerper = {}
        for b, cn, on in _alle_koerper(occ):
            for nm in (b.name, cn, on):
                if nm.startswith("Vorderseite"): koerper.setdefault("Vorderseite", b)
                if nm.startswith("Rueckseite"): koerper.setdefault("Rueckseite", b)
        reihe = occ.component.name.replace("Stimmstock_", "").replace("_fraesen", "")
        for seite in ("Vorderseite", "Rueckseite"):
            b = koerper.get(seite)
            if b is None: log.append("%s: %s fehlt" % (occ.component.name, seite)); continue
            si = cam.setups.createInput(adsk.cam.OperationTypes.MillingOperation); si.models = [b]
            si.name = "Stimmstock %s %s" % (reihe, seite); s = cam.setups.add(si)
            for n, e in spar.items(): _p(s.parameters, n, e)
            for n, e in (("job_stockMode", "'default'"), ("job_stockOffsetMode", "'simple'"), ("job_stockOffsetSides", "0 mm"),
                         ("job_programName", "'%s-%s'" % (reihe, seite))): _p(s.parameters, n, e)
            try:
                _p(s.parameters, "wcs_origin_mode", "'point'")
                s.parameters.itemByName("wcs_origin_point").value.value = [occ.component.originConstructionPoint.createForAssemblyContext(occ)]
            except: _p(s.parameters, "wcs_origin_mode", "'modelOrigin'")
            boden = _bodenflaeche(b)
            rand = [boden.loops.item(i).edges.item(0) for i in range(boden.loops.count) if boden.loops.item(i).isOuter]
            for k, (strat, par) in enumerate(opar):
                op = _op(s, strat, "Ausraeumen %d" % (k + 1), W["2.5"], {}, grenze=rand, enthalten="'outside'")
                for n, e in par.items(): _p(op.parameters, n, e)
                for n, e in feeds.items(): _p(op.parameters, n, e)
                _p(op.parameters, "boundaryMode", "'selection'"); _p(op.parameters, "boundaryContainment", "'outside'"); _p(op.parameters, "boundaryOffset", "0 mm")
    # --- Griffstock ---
    for tn, occ in sorted(gs_occ.items()):
        T = J["teile"][tn]; g = geo[tn]
        body = max([b for b, _, _ in _alle_koerper(occ)], key=lambda b: b.volume)
        wk = W["3.175"] if abs(T.get("fraeser_kanal", 2.5) - 3.175) < 0.01 else W["2.5"]
        for seite in ("oben", "Deckelseite"):
            si = cam.setups.createInput(adsk.cam.OperationTypes.MillingOperation); si.models = [body]
            si.name = "Griffstock %s %s" % (tn, seite); s = cam.setups.add(si)
            for n, e in (("job_stockMode", "'default'"), ("job_stockOffsetMode", "'simple'"), ("job_stockOffsetSides", "0 mm"),
                         ("job_stockOffsetTop", "1 mm"), ("job_stockOffsetBottom", "1 mm"),
                         ("job_programName", "'GS-%s-%s'" % (tn, seite))): _p(s.parameters, n, e)
            if seite == "Deckelseite" and not _p(s.parameters, "wcs_orientation_flipZ", "true"): log.append(s.name + ": Z nicht gewendet")
            sil = {"boundaryMode": "'silhouette'", "boundaryContainment": "'inside'", "machineCavities": "true"}
            if seite == "oben":
                par = dict(sil); par.update({"useStockToLeave": "true", "simpleStockToLeave": "false", "stockToLeave": "0 mm", "verticalStockToLeave": "0.5 mm"})
                _op(s, "adaptive", "Schruppen 6 mm (0,5 Aufmass Boeden)", W["6"], par)
                _flaechen(s, g["flaechen"]["oben"], W["6"])
                for n_, t_, gr in g["kanaele"]:
                    _op(s, "adaptive", "Kanal %d (%.1f tief)" % (n_, t_), wk, {"useStockToLeave": "false", "machineCavities": "true"}, grenze=gr)
                try: _op(s, "contour2d", "Innenraum ausschneiden (Rechteck)", W["6"], {"doMultipleDepths": "true", "maximumStepdown": "4 mm"}, kontur=g["innen"])
                except: log.append(s.name + ": 2D-Kontur Innenraum nicht moeglich: " + traceback.format_exc().splitlines()[-1])
            else:      # Rueckseite: schruppen mit 0,5 Aufmass auf den Boeden, dann jede Flaeche schlichten
                par = dict(sil); par.update({"useStockToLeave": "true", "simpleStockToLeave": "false", "stockToLeave": "0 mm", "verticalStockToLeave": "0.5 mm"})
                _op(s, "adaptive", "Schruppen 6 mm" + (" (Deckel-Ausnehmung vollflaechig)" if J["teile"][tn].get("deckel_ausnehmung") else ""), W["6"], par)
                _flaechen(s, g["flaechen"]["unten"], W["6"])
    doc.save("Fraesteile %s %s (erzeugt %s), CAM neu, Werkzeugwege noch nicht berechnet" % (L["name"], L.get("version", ""), L.get("erzeugt", "")))
    log.append("%s: %d Stimmstoecke, Griffstock %d Kanaele, %d Setups, gespeichert" % (L["dokument"], len(st_occ), sum(len(J["teile"][t]["kanaele"]) for t in J["teile"]), cam.setups.count))

def run(context):
    app = adsk.core.Application.get(); ui = app.userInterface
    try:
        listen = []
        for d in sorted(os.listdir(BASIS)):
            f = os.path.join(BASIS, d, "fraes_liste.json")
            if os.path.isfile(f): listen.append((os.path.join(BASIS, d), json.load(open(f, encoding="utf-8"))))
        if not listen: ui.messageBox("Keine fraes_liste.json unter " + BASIS); return
        jakob = _jakob(app)
        for ordner, L in listen: _ein_instrument(app, ui, ordner, L, jakob)
        open(LOG, "w", encoding="utf-8").write("\n".join(log))
        ui.messageBox("Fraesteile aktualisiert (ohne Berechnung):\n\n" + "\n".join(l for l in log if "gespeichert" in l or "nicht" in l or "fehlt" in l))
    except:
        log.append(traceback.format_exc()); open(LOG, "w", encoding="utf-8").write("\n".join(log))
        ui.messageBox("Fehler:\n" + traceback.format_exc())
