# Modelle automatisch auf GitHub rechnen

Die Modelle lassen sich direkt auf GitHub rechnen, ohne dass am eigenen Rechner etwas läuft. Dafür gibt es den Ablauf „Harmonika erzeugen“ (`.github/workflows/harmonika_erzeugen.yml`). Weil das Repo öffentlich ist, kostet das nichts; jeder Lauf bekommt eine eigene Maschine mit 4 Kernen und 16 GB Speicher und darf bis zu zweieinhalb Stunden rechnen.

Gestartet wird auf github.com/jpascher/Harmonikabau unter **Actions → Harmonika erzeugen → Run workflow**. Dort wählt man in Auswahlfeldern, ob ein Bass- oder Diskantteil entstehen soll, welche Variante, die Knopfteilung und beim Markus-Bass die Gehäuselänge. Beim Diskant steht alles Übrige (Knopfzahl, Teilung, Griffstock, Chörigkeit) in der gewählten Tastenliste aus `diskant/disk/`. Starten darf nur, wer Schreibrecht am Repo hat.

Nach etwa 15 bis 30 Minuten liegt das Ergebnis unter **Releases** als eigener Eintrag „modell-Nummer“ mit festen Download-Links: das gepackte STEP-Modell (beim Diskant zusammen mit allen Fräsdaten), eine grobe farbige Vorschau `vorschau.glb`, ein Bild `vorschau.png` und das Protokoll des Generators samt Kollisionsprüfung. Dieselben Dateien hängen 30 Tage lang auch als Artefakt am Lauf. Schlägt ein Lauf fehl, wird wenigstens das Protokoll gesichert.

Die Übersichtsseite `docs/index.html` zeigt alle gerechneten Modelle mit Vorschaubild und Download-Links und hat einen Knopf zum Starten. Sie wird einmalig unter **Settings → Pages** mit „Deploy from a branch“, Branch `main`, Ordner `/docs` eingeschaltet und ist dann unter jpascher.github.io/Harmonikabau erreichbar. Die 3D-Vorschau `vorschau.glb` kann man dort nach dem Herunterladen auswählen und drehbar ansehen.

Technisch ruft der Ablauf `markus-3d/ci/erzeuge_ci.py` auf, das je nach Auswahl den passenden Generator startet (`bass_erzeuge.py`, `kammern38.py`, `eduard.py`, `gerhard.py`, `andreas.py`, `v340.py` oder beim Diskant `disk/erzeuge.py`), das Ergebnis packt und mit `vorschau.py` die Vorschau erzeugt. Die Python-Pakete stehen mit festen Versionen in `ci/requirements.txt`. Am eigenen Rechner funktioniert derselbe Aufruf genauso, z. B. `python3 markus-3d/ci/erzeuge_ci.py bass v340 standard 380 ergebnis`.

Hinweis zum Diskant: Die Dateien vom 23.11. (Knopfanzahl-Varianten, Griffstock-Endform ital) lagen bisher nur in `diskant/` statt in `diskant/disk/`, wo `erzeuge.py` sie sucht. Sie sind jetzt auch nach `diskant/disk/` kopiert; die Kopien in `diskant/` selbst können gelöscht werden.
