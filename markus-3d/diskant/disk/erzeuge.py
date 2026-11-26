# Ein Aufruf je Instrument (Johann 19.11.2026): alles aus der Tastenliste neu erzeugen.
#   python3 erzeuge.py <tastenliste.txt> <ausgabeordner> [version]
# Ergebnis im Ausgabeordner (auf den PC nach C:\ACAD-DWG\0-Markus-3D\Fraesen\<Name>\ kopieren):
#   Diskantteil_<Name>_<version>.step (+ .kollision.json)   ganzes Modell
#   Stimmstoecke\Stimmstock_<reihe>_<ton>_fraesen.step/.dxf  Stimmstock-Fraesteile (je Reihe Vorder-/Rueckseite)
#   Griffstock\Griffstock_Teil1/Teil2_<Name>_fraesen.step/.dxf, Griffstock_<Name>_fraesen.json
#   fraes_liste.json                                          alles, was die Fusion-Skripte brauchen (keine festen Namen dort)
# Aendert sich Tastenzahl, Knopfabstand oder Plattenlaenge, genuegt dieser Aufruf; danach in Fusion "Fraesen_Aktualisieren".
import sys, os, subprocess, json, glob, datetime

def main(liste, ordner, version="neu"):
    hier = os.path.dirname(os.path.abspath(__file__))
    sys.path.insert(0, hier)
    import instrument
    _inst = instrument.lies(os.path.abspath(liste))
    name = str(_inst["name"])
    winkel = str(_inst.get("griffstock", "")).startswith("winkel")
    os.makedirs(ordner, exist_ok=True)
    st_dir, gs_dir = os.path.join(ordner, "Stimmstoecke"), os.path.join(ordner, "Griffstock")
    for d in (st_dir, gs_dir):
        os.makedirs(d, exist_ok=True)
        for f in glob.glob(os.path.join(d, "*")): os.remove(f)          # alte Teile weg (z. B. andere Reihenzahl)
    gs_step = os.path.join(ordner, "_griffstock_roh.step")
    tmp_out = os.path.join(ordner, "_modell.step")
    teile_dir = os.path.join(ordner, "Teile")
    os.makedirs(teile_dir, exist_ok=True)
    for f in glob.glob(os.path.join(teile_dir, "*")): os.remove(f)
    dump, kst_step = os.path.join(ordner, "_teile.pkl"), os.path.join(ordner, "_kunststoff.step")
    env = dict(os.environ, INSTRUMENT=os.path.abspath(liste), KOLLISION="1", FRAES=st_dir, GS_STEP=gs_step,
               TEILE_DUMP=dump, KST_STEP=kst_step)
    print("1) Modell, Stimmstock-Fraesteile, Griffstock ...", flush=True)
    r = subprocess.run([sys.executable, "diskant.py", "BESASDES", tmp_out], cwd=hier, env=env, capture_output=True, text=True)
    open(os.path.join(ordner, "generator.log"), "w").write(r.stdout + r.stderr)
    if r.returncode != 0: raise SystemExit("Generator-Fehler, siehe generator.log")
    modell = os.path.join(ordner, f"Diskantteil_{name}_{version}.step")
    os.replace(tmp_out, modell)
    if os.path.exists(tmp_out + ".kollision.json"): os.replace(tmp_out + ".kollision.json", modell + ".kollision.json")
    print("2) Griffstock-Fraesdateien ...", flush=True)
    r2 = subprocess.run([sys.executable, "gs_fraes.py", gs_step, gs_dir, name], cwd=hier, capture_output=True, text=True)
    open(os.path.join(ordner, "gs_fraes.log"), "w").write(r2.stdout + r2.stderr)
    gs_ok = r2.returncode == 0
    if not gs_ok and not winkel: raise SystemExit("gs_fraes-Fehler, siehe gs_fraes.log")
    if not gs_ok: print("   Hinweis: Winkel-Griffstock - Griffstock-Fraesdateien werden noch nicht automatisch erzeugt (26.11.2026)", flush=True)
    if os.path.exists(gs_step): os.remove(gs_step)
    print("3) Alu-Abwicklungen, Holzteile, Kunststoffteil ...", flush=True)
    r3 = subprocess.run([sys.executable, "teile_dxf.py", dump, kst_step, teile_dir], cwd=hier, capture_output=True, text=True)
    open(os.path.join(ordner, "teile_dxf.log"), "w").write(r3.stdout + r3.stderr)
    if r3.returncode != 0: raise SystemExit("teile_dxf-Fehler, siehe teile_dxf.log")
    for f_ in (dump, kst_step, kst_step.replace('.step', '_fertig.step')):
        if os.path.exists(f_): os.remove(f_)
    stoecke = sorted(os.path.basename(f)[:-5] for f in glob.glob(os.path.join(st_dir, "Stimmstock_*_fraesen.step")))
    liste_ = {"name": name, "version": version, "erzeugt": datetime.datetime.now().isoformat(timespec="minutes"),
              "dokument": f"Stimmstoecke_{name}_Fraesen", "modell": os.path.basename(modell),
              "stimmstoecke": {"ordner": "Stimmstoecke", "teile": stoecke},
              "griffstock": {"ordner": "Griffstock", "json": f"Griffstock_{name}_fraesen.json",
                             "teile": {t: f"Griffstock_{t}_{name}_fraesen" for t in ("Teil1", "Teil2")}}}
    if not gs_ok: liste_["griffstock"] = None
    json.dump(liste_, open(os.path.join(ordner, "fraes_liste.json"), "w", encoding="utf-8"), indent=1, ensure_ascii=False)
    print(f"fertig: {name}, {len(stoecke)} Stimmstoecke" + (", Griffstock Teil 1/2" if gs_ok else ", ohne Griffstock-Fraesdateien (Winkel)") + f" -> {ordner}")

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3] if len(sys.argv) > 3 else "neu")
