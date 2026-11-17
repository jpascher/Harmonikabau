# Steirische „Markus“ – 3D-Modell (Diskant- und Bassteil)

Dieser Ordner enthält die Quellen, aus denen das 3D-Modell der Steirischen Harmonika „Markus“ entsteht: die Python-Skripte, die die Teile rechnen, und die DXF-/CamBam-Zeichnungen, aus denen sie gelesen werden. Die fertigen Modelle selbst (STEP-Dateien mit 30 bis 200 MB) und die Videos liegen nicht im Repo, weil Git mit so großen Binärdateien schlecht umgeht; sie werden als Anhang eines GitHub-Releases bereitgestellt. Aktueller Stand ist der Diskantteil v21 (17.11.2026) und der Bassteil v46.

## Wie das Modell entsteht

Alle Teile werden nicht von Hand konstruiert, sondern aus den Werkstattzeichnungen gerechnet. Die Skripte lesen die Konturen aus den DXF-Dateien (dieselben Zeichnungen, nach denen gefräst wird), extrudieren sie mit CadQuery zu Körpern und setzen sie zum Diskant- bzw. Bassteil zusammen. Ändert man eine Zeichnung, rechnet man das Modell einfach neu. Gebraucht werden Python 3 mit cadquery 2.8, shapely 2.2, ezdxf 1.4 und pyvista 0.49, für die Videos zusätzlich ffmpeg und xvfb. Die Rechnung dauert für den ganzen Diskantteil etwa acht Minuten und läuft am besten auf einem schnellen Rechner oder in der Cloud.

## Diskantteil (`diskant/`)

Das Hauptskript ist `disk/diskant.py`. Es baut Boden, Rahmen, Stimmstöcke, Griffstock, Verdeck und über `hebel.py` die ganze Diskantmechanik mit 50 Hebeln. Als erstes Argument bekommt es die Stimmung, also den Namen eines Ordners unter `disk/var/` (z. B. `BESASDES`, `GCFB`, `ADGC`), dort liegen die Stimmstock-Zeichnungen je Reihe. Aufruf: `cd disk` und dann `python3 diskant.py BESASDES Diskantteil_Markus_vNN.step`. Die Prüfskripte `jp_check2.py` (Hebel in Ruhe und gedrückt auf Kollisionen) und `stz_check.py` (Verdeckstützen und Schrauben) arbeiten mit den BREP-Teilen, die `voll_mit_brep.py` zusätzlich ausgibt. Die Rahmenteile rechnet `rahmen_d.py`, die Videos `anim_d/anim_diskant.py` (Modi `mit`, `ohne`, `explo`). Der wichtigste Eingang für die Mechanik ist `disk/jp_schnitt.dxf`, der überarbeitete Hebelschnitt mit Holzteilen, Alu, Knöpfen und Kammern.

Die drei Ordner `fusion_*` sind Skripte für Autodesk Fusion. Man bindet sie dort über Dienstprogramme → Zusatzmodule → „+“ → „Skript vom Gerät“ ein. `Hebel_Gelenke` legt die 50 Drehgelenke an, `Hebel_druecken` schaltet einen einzelnen Hebel um (Name oben im Skript ändern), und `Diskant_Ablauf` spielt den ganzen Ablauf vor: Verdeck weg, Griffstock aus, Reihen drücken, alles zurück.

Genaue Maße und Festlegungen des aktuellen Stands stehen in `diskant/STAND_Diskant_Markus.md`. Offen ist noch, dass bei der Explosionsdarstellung Griffstock-Teil 2, Verdeckleiste und Stützen nicht mitfahren.

## Bassteil (`bass/`)

Das Hauptskript ist `kammern.py`: Es setzt den Rahmen aus Front, Rückteil und Seitenteilen zusammen, baut Boden, Basskammern, Begleiter-Stimmstöcke, Deckel und Hörner und bindet über `mechanik.py` die Bassmechanik mit Klappen, Achsen, Lagerböcken, Schubstangen und Knöpfen ein. Aufruf: `python3 kammern.py boden2.dxf Bassteil_Markus_vNN.step`. Die Rahmenteile (`Markus_front_v4.step`, `Markus_rueck_v4.step`, `Markus_seite_v4.step`), die `kammern.py` einliest, werden vorher mit `teile.py` aus der Rahmenzeichnung erzeugt. `bau_bassteil.py` ist das erste, vereinfachte Modell aus Draufsicht und Schnitt und dient nur noch als Bezug. Die Videos entstehen mit `anim.py` und `render.py` (Tasten nacheinander drücken) sowie `explo_daten.py` und `explo.py` (Explosion). Die Unterlagen im Ordner `mech_d/` sind die Diskant-Hebelzeichnungen, die beim Bassteil zum Vergleich mitgenommen wurden.

Einige Skripte enthalten noch fest eingetragene Pfade aus der Arbeitsumgebung, in der sie entstanden sind (z. B. in `anim.py`); diese vor dem Ausführen an den eigenen Ordner anpassen.

## Lizenz und Verwendung

Zeichnungen und Skripte stammen von Johann Pascher. Sie dürfen zum Nachbauen, Lernen und Weiterentwickeln frei verwendet werden; Ableitungen sollen ebenfalls frei bleiben.
