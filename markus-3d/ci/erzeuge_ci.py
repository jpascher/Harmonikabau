# Harmonika-Modell automatisch erzeugen (GitHub Actions, Johann 26.11.2026)
#   python3 markus-3d/ci/erzeuge_ci.py <teil> <variante> <teilung> <laenge> <ausgabeordner> [vorschau ja|nein]
#   teil      : bass | diskant
#   variante  : Bass  -> markus, marie38, eduard, gerhard, andreas, v340
#               Diskant -> instrument:griffstock:knoepfe:knopfabstand (z. B. markus:kugellager:50:19, siehe tastenliste_bauen.py)
#                          oder Name einer vorhandenen Tastenliste ohne "tastenliste_" und ".txt", z. B. marie_winkel
#   teilung   : Knopfteilung in mm (18, 19, 20) oder "standard" (Wert des jeweiligen Skripts); beim Diskant steht sie in der Tastenliste
#   laenge    : Gehaeuselaenge (Aussenmass fertig) in mm, leer = wie gezeichnet bzw. beim Diskant wie die Knopfzahl verlangt;
#               fuer Bass und Diskant gleich einzugeben (27.11.2026), damit beide Teile dieselben Aussenmasse bekommen
#   Umgebungsvariable BREITE: Gehaeusebreite (Tiefe, Aussenmass fertig) 200-215 mm, leer = 205 (masse.py)
#   Umgebungsvariable BELEGUNG: freie Bassbelegung in Kurzform "1.1=V077, 1.2=BL152, ..." (bisher nur v340)
#   Umgebungsvariable PLATTEN: eigene Stimmplattenlaengen fuer den Diskant, z. B. "R1: 32 33 ...; R2: ..." (leer = gezeichnet)
# Ergebnis im Ausgabeordner: <Modell>_step.zip (beim Diskant samt Fraesdaten), vorschau.glb, vorschau.png, notizen.md, protokoll.txt
import sys, os, glob, shutil, subprocess, time, zipfile

HIER = os.path.dirname(os.path.abspath(__file__))
M3D = os.path.dirname(HIER)
BASS = {"markus": ("bass", None), "marie38": ("bass/marie38", "kammern38.py"), "eduard": ("bass/eduard", "eduard.py"),
        "gerhard": ("bass/gerhard", "gerhard.py"), "andreas": ("bass/andreas", "andreas.py"), "v340": ("bass/v340", "v340.py")}

NAMEN = {"marie38": "Marie38", "eduard": "Eduard", "gerhard": "Gerhard", "andreas": "Andreas", "v340": "340"}

def lauf(cmd, cwd, log, env=None):
    t0 = time.time()
    with open(log, "a", encoding="utf-8") as f:
        f.write("\n$ " + " ".join(cmd) + f"   (in {os.path.relpath(cwd, M3D)})\n"); f.flush()
        r = subprocess.run(cmd, cwd=cwd, stdout=f, stderr=subprocess.STDOUT, env=env)
        f.write(f"-> Rueckgabe {r.returncode}, {(time.time() - t0) / 60:.1f} min\n")
    print(open(log, encoding="utf-8").read()[-3000:], flush=True)
    if r.returncode != 0: raise SystemExit(f"Fehler bei {' '.join(cmd)} (siehe protokoll.txt)")

def packe(zipname, dateien, basis):
    with zipfile.ZipFile(zipname, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as z:
        for p in dateien: z.write(p, os.path.relpath(p, basis))

def masse_anpassen(step, teil, breite, laenge, log, ordner):
    """Breite (und beim Bass ausser Markus die Laenge) mit masse.py nachtraeglich aendern; gibt den neuen STEP zurueck"""
    b_ = (breite or "").strip(); l_ = (laenge or "").strip()
    if b_ in ("", "205"): b_ = ""
    if not b_ and not l_: return step
    neu = step[:-5] + (f"_B{b_}" if b_ else "") + (f"_L{l_}" if l_ else "") + ".step"
    cmd = [sys.executable, os.path.join(HIER, "masse.py"), step, neu, "--teil", teil, "--bericht", os.path.join(ordner, "masse.json")]
    if b_: cmd += ["--breite", b_]
    if l_: cmd += ["--laenge", l_]
    lauf(cmd, HIER, log)
    if os.path.abspath(step) != os.path.abspath(neu) and os.path.exists(step) and teil == "bass": os.remove(step)
    return neu

def main():
    teil, variante, teilung, laenge, aus = sys.argv[1:6]
    vorschau = (sys.argv[6] if len(sys.argv) > 6 else "ja").lower() in ("ja", "true", "1")
    aus = os.path.abspath(aus); os.makedirs(aus, exist_ok=True)
    log = os.path.join(aus, "protokoll.txt"); open(log, "w").close()
    py = sys.executable
    if teil == "bass":
        if variante not in BASS: raise SystemExit(f"unbekannte Bassvariante {variante}")
        ordner, skript = BASS[variante]
        cwd = os.path.join(M3D, ordner)
        tag = "" if teilung in ("", "standard") else f"_T{teilung}"
        if variante == "markus":                              # Laenge im Generator selbst (streck.py)
            t_ = "20" if teilung in ("", "standard") else teilung
            l_ = laenge or "380"
            step = os.path.join(aus, f"Bassteil_Markus_T{t_}_L{l_}.step")
            lauf([py, "bass_erzeuge.py", t_, l_, step], cwd, log)
        else:
            bel = os.environ.get("BELEGUNG", "").strip()
            env_ = None
            if bel:                                               # freie Belegung (Webseite "Bassbelegung zuordnen")
                mj = glob.glob(os.path.join(cwd, "belegung_modell_*.json"))
                if not mj: raise SystemExit(f"Freie Belegung gibt es fuer {variante} noch nicht (bisher nur v340)")
                bel_datei = os.path.join(aus, "belegung.txt")
                open(bel_datei, "w", encoding="utf-8").write(f"modell = {variante}\n" + "\n".join(t.strip() for t in bel.replace(";", ",").split(",") if t.strip()) + "\n")
                env_ = dict(os.environ, BASS_BELEGUNG=bel_datei)
                tag += "_Belegung"
            step = os.path.join(aus, f"Bassteil_{NAMEN[variante]}{tag}.step")
            lauf([py, skript, step] + ([] if teilung in ("", "standard") else [teilung]), cwd, log, env_)
        modell = step
        l_aendern = laenge if variante != "markus" else ""
        modell = masse_anpassen(modell, "bass", os.environ.get("BREITE", ""), l_aendern, log, aus)
        name = os.path.basename(modell)[:-5]
        packe(os.path.join(aus, name + "_step.zip"), [modell] + glob.glob(os.path.join(aus, "belegung.txt")) + glob.glob(os.path.join(aus, "masse.json")), aus)
    elif teil == "diskant":
        if ":" in variante:                                   # Auswahl der Webseite: instrument:griffstock:knoepfe:knopfabstand
            sys.path.insert(0, HIER); import tastenliste_bauen
            liste, lname, hinweise = tastenliste_bauen.baue(*variante.split(":"), os.path.join(aus, "tastenliste.txt"),
                                                            os.environ.get("PLATTEN", ""), laenge)
            with open(log, "a", encoding="utf-8") as f:
                f.write(f"Tastenliste {lname} zusammengesetzt\n" + "".join(f"Hinweis: {h}\n" for h in hinweise))
        else:
            liste = os.path.join(M3D, "diskant", "disk", f"tastenliste_{variante}.txt")
            if not os.path.exists(liste): raise SystemExit(f"Tastenliste fehlt: {liste}")
            hinweise = []
        roh = os.path.join(aus, "diskant_roh")
        lauf([py, "erzeuge.py", liste, roh, "ci"], os.path.join(M3D, "diskant", "disk"), log)
        modell = glob.glob(os.path.join(roh, "Diskantteil_*_ci.step"))[0]
        neu_ = masse_anpassen(modell, "diskant", os.environ.get("BREITE", ""), "", log, roh)
        if neu_ != modell: os.remove(modell); modell = neu_
        name = os.path.basename(modell)[:-5]
        alle = [p for p in glob.glob(os.path.join(roh, "**", "*"), recursive=True) if os.path.isfile(p)]
        shutil.copy(liste, os.path.join(roh, "tastenliste.txt"))
        alle = [p for p in glob.glob(os.path.join(roh, "**", "*"), recursive=True) if os.path.isfile(p)]
        packe(os.path.join(aus, name + "_mit_Fraesdaten.zip"), alle, roh)
    else:
        raise SystemExit("teil muss bass oder diskant sein")
    hinweise = locals().get("hinweise", [])
    if vorschau:
        lauf([py, os.path.join(HIER, "vorschau.py"), modell, os.path.join(aus, "vorschau.glb"), os.path.join(aus, "vorschau.png")], HIER, log)
    groesse = os.path.getsize(modell) / 1e6
    with open(os.path.join(aus, "notizen.md"), "w", encoding="utf-8") as f:
        f.write(f"**{name}** – Teil: {teil}, Variante: {variante}, Teilung: {teilung or 'standard'}"
                + (f", Länge: {laenge} mm" if laenge else "") + (f", Breite: {os.environ.get('BREITE')} mm" if os.environ.get("BREITE", "").strip() not in ("", "205") else "") + "\n\n"
                f"STEP {groesse:.0f} MB, gepackt im ZIP. `vorschau.glb` (grob, mit Farben) und `vorschau.png` zur schnellen Kontrolle; "
                f"`protokoll.txt` enthält die Ausgabe des Generators einschließlich Kollisionsprüfung.\n"
                + "".join(f"\n> Hinweis: {h}" for h in hinweise))
    for p in (modell, os.path.join(aus, "diskant_roh"), os.path.join(aus, "tastenliste.txt"), os.path.join(aus, "belegung.txt"), os.path.join(aus, "masse.json")):          # nur die gepackten Dateien bleiben
        if os.path.isdir(p): shutil.rmtree(p)
        elif os.path.exists(p): os.remove(p)
    print("fertig:", sorted(os.listdir(aus)))

if __name__ == "__main__":
    main()
