# Diskant – Dateien für GitHub (Stand 20.11.2026)

Alle ZIPs sind in Repo-Struktur aufgebaut. Jedes ZIP wird im Wurzelordner des Repos `Harmonikabau` entpackt und landet dann automatisch unter `markus-3d/diskant/`. Bestehende Dateien werden dabei überschrieben.

## 1. Diskant_Repo_20-11.zip (4,4 MB) – Quellen und Skripte

Dieses ZIP enthält die aktuellen Quellen: `disk/` mit `diskant.py`, `hebel.py` (Messingbuchse, Verjüngung des langen Arms, Rastmarken), `gs_fraes.py`, `teile_dxf.py`, `erzeuge.py`, `instrument.py`, den Tastenlisten für Marie (mit Schräglage) und Markus und `stimmstock_parameter.txt`. Dazu kommen die Fusion-Skripte `fusion_Fraesen_Aktualisieren` und `fusion_Modelle_Laden12` (neu) sowie `fusion_Diskant_Ablauf`, `fusion_Hebel_Gelenke` und `fusion_Hebel_druecken`. Die Beschreibungen `LIESMICH_Fraesen.md` und `STAND_Diskant_Markus.md` sowie `anim_d/`, `rahmen.dxf` und `rahmen_d.py` sind ebenfalls dabei. Sicherungskopien (`*_vXX_vor_*`), Debug-Ausgaben und Logs sind bewusst weggelassen.

## 2. Diskant_Fraesdaten_20-11.zip (22 MB) – erzeugte Fräsdaten

Dieses ZIP entpackt nach `markus-3d/diskant/Fraesen/Marie` und `markus-3d/diskant/Fraesen/Markus`. Es enthält Stimmstöcke, Griffstock und Teile (Alu-Hebel, Holzteile, Kunststoffteil) sowie `fraes_liste.json` und `kollision.json`. Es ist der neueste Stand mit den korrigierten Fräsflächen; seit 20.11. treffen alle 50 Rastmarken im Kunststoffteil Material (Reihe 4 um 1,25 mm neben den Hebelschlitz versetzt). Die Daten lassen sich jederzeit mit `python3 erzeuge.py <tastenliste> <ordner> <version>` neu erzeugen.

## 3. Diskantteil_Marie_v16_step.zip (29 MB) und Diskantteil_Markus_v33_step.zip (34 MB) – 3D-Modelle

Diese beiden ZIPs entpacken nach `markus-3d/diskant/Modelle/`.

**Achtung:** Entpackt sind die STEP-Dateien 178 MB bzw. 208 MB groß. GitHub lehnt jede einzelne Datei über 100 MB ab, der Push schlägt dann fehl. Dafür gibt es zwei Lösungen:

- **Einfach (empfohlen):** Die Modelle nicht entpacken, sondern die ZIP-Dateien selbst ins Repo legen, z. B. nach `markus-3d/diskant/Modelle/`. Beide liegen deutlich unter 100 MB.
- **Mit Git LFS:** Die `.step`-Dateien vor dem Hinzufügen über Git LFS verfolgen lassen:

  ```
  git lfs install
  git lfs track "*.step"
  git add .gitattributes
  ```

  Danach die STEP-Dateien normal hinzufügen. Beim kostenlosen GitHub-Konto ist der LFS-Speicher und -Datenverkehr begrenzt (aktuelles Kontingent in den GitHub-Einstellungen unter „Billing“ prüfen).

## Ablage

Alle ZIPs liegen auf dem Stick unter `D:\0-Markus-3D` und unter `C:\ACAD-DWG\0-Markus-3D`. Die Original-STEP-Modelle liegen unter `C:\ACAD-DWG\0-Markus-3D\0-marie`.
