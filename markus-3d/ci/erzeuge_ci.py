# Harmonika-Modell automatisch erzeugen (GitHub Actions, Johann 26.11.2026)
#   python3 markus-3d/ci/erzeuge_ci.py <teil> <variante> <teilung> <laenge> <ausgabeordner> [vorschau ja|nein]
#   teil      : bass | diskant
#   variante  : Bass  -> markus, marie38, eduard, gerhard, andreas, v340
#               Diskant -> Name der Tastenliste ohne "tastenliste_" und ".txt", z. B. markus, marie_winkel, markus48_winkel
#   teilung   : Knopfteilung in mm (18, 19, 20) oder "standard" (Wert des jeweiligen Skripts); beim Diskant steht sie in der Tastenliste
#   laenge    : Gehaeuselaenge in mm, nur fuer den Markus-Bass (Diskantlaenge genuegt, wird auf volle 10 mm aufgerundet)
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
        if variante == "markus":
            t_ = "20" if teilung in ("", "standard") else teilung
            l_ = laenge or "380"
            step = os.path.join(aus, f"Bassteil_Markus_T{t_}_L{l_}.step")
            lauf([py, "bass_erzeuge.py", t_, l_, step], cwd, log)
        else:
            step = os.path.join(aus, f"Bassteil_{NAMEN[variante]}{tag}.step")
            lauf([py, skript, step] + ([] if teilung in ("", "standard") else [teilung]), cwd, log)
        modell = step
        name = os.path.basename(step)[:-5]
        packe(os.path.join(aus, name + "_step.zip"), [step] + glob.glob(step + ".*"), aus)
    elif teil == "diskant":
        liste = os.path.join(M3D, "diskant", "disk", f"tastenliste_{variante}.txt")
        if not os.path.exists(liste): raise SystemExit(f"Tastenliste fehlt: {liste}")
        roh = os.path.join(aus, "diskant_roh")
        lauf([py, "erzeuge.py", liste, roh, "ci"], os.path.join(M3D, "diskant", "disk"), log)
        modell = glob.glob(os.path.join(roh, "Diskantteil_*_ci.step"))[0]
        name = os.path.basename(modell)[:-5]
        alle = [p for p in glob.glob(os.path.join(roh, "**", "*"), recursive=True) if os.path.isfile(p)]
        packe(os.path.join(aus, name + "_mit_Fraesdaten.zip"), alle, roh)
    else:
        raise SystemExit("teil muss bass oder diskant sein")
    if vorschau:
        lauf([py, os.path.join(HIER, "vorschau.py"), modell, os.path.join(aus, "vorschau.glb"), os.path.join(aus, "vorschau.png")], HIER, log)
    groesse = os.path.getsize(modell) / 1e6
    with open(os.path.join(aus, "notizen.md"), "w", encoding="utf-8") as f:
        f.write(f"**{name}** – Teil: {teil}, Variante: {variante}, Teilung: {teilung or 'standard'}"
                + (f", Länge: {laenge}" if teil == "bass" and variante == "markus" else "") + "\n\n"
                f"STEP {groesse:.0f} MB, gepackt im ZIP. `vorschau.glb` (grob, mit Farben) und `vorschau.png` zur schnellen Kontrolle; "
                f"`protokoll.txt` enthält die Ausgabe des Generators einschließlich Kollisionsprüfung.\n")
    for p in (modell, os.path.join(aus, "diskant_roh")):          # nur die gepackten Dateien bleiben
        if os.path.isdir(p): shutil.rmtree(p)
        elif os.path.exists(p): os.remove(p)
    print("fertig:", sorted(os.listdir(aus)))

if __name__ == "__main__":
    main()
