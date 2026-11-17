# Diskantteil „Markus“ – Stand 17.11.2026

Aktueller Stand ist **Diskantteil_Markus_v21.step**. In Fusion ist er gespeichert als „Diskantteil_Markus_v21“, mit 50 Hebel-Drehgelenken.

## Erzeugen

Die Rechnung läuft in der Cloud bzw. auf einem schnellen Rechner. Benötigt werden Python 3 mit cadquery 2.8, shapely 2.2, ezdxf 1.4 und pyvista 0.49; für die Videos zusätzlich ffmpeg und xvfb.

```
cd disk
python3 diskant.py BESASDES Diskantteil_Markus_vNN.step          # ganzer Diskantteil (ca. 8 min)
python3 voll_mit_brep.py Diskantteil_Markus_vNN.step             # dasselbe + BREPs in jp_brep/ fuer die Pruefskripte
python3 jp_check2.py $(ls jp_brep | grep Knopf | sed 's/__Knopf.brep//')   # Hebel-Kollisionen Ruhe/gedrueckt
python3 stz_check.py                                             # Verdeckstuetzen + Schrauben
```

Mit `NUR="Diskantboden|Griffstock|...|Verdeck"` wird nur eine Teilmenge ausgegeben. `ANIMDIR=../anim_d` schreibt die Netze für die Animationen.

Animationen:

```
cd anim_d
python3 anim_diskant.py mit|ohne|explo Datei.mp4
```

Mit `OHNE_VD=_nv` davor entsteht die Fassung ohne Verdeck.

Rahmenteile mit Höhe 51:

```
python3 rahmen_d.py rahmen.dxf disk/Dr_front_51.step front 51
python3 rahmen_d.py rahmen.dxf disk/Dr_rueck_51.step rueck 51
python3 rahmen_d.py rahmen.dxf disk/Dr_seite_51.step seite 51
```

## Eingaben

| Datei | Inhalt |
|---|---|
| jp_schnitt.dxf | Johanns überarbeiteter Hebelschnitt: Holzteile, Alu, Knöpfe, Kammern, Teilung (Layer 0) |
| Schnitt_Hebel_4_Layer.dxf | ursprünglicher Schnitt, als Bezug für die Alu-Änderungen |
| blume1.dxf | Rosette für die Verdeck-Front |
| venus.dxf | Venus-Blume (Front-einlage) für die Rahmen-Einlage vorne |
| eduard_einlage.dxf | nicht mehr verwendet |
| verdeck.dxf, gs.dxf, gs2.dxf, gschnitt.dxf, boden.dxf, rigel_j.dxf, kst_cb.json, drehpunkte.json, var/ | Grundlagen: Verdeck, Griffstock-Schnitte, Boden, Kunststoffteil, Stimmstöcke |

## Stand v21 (wichtigste Festlegungen)

- **Rahmen:** 51 hoch, Randleiste oben bei 61 (Stimmstöcke bis 60,5). Eck-Ornamente zentriert, auf 75 % verkleinert. Vorne mittig die Venus-Blume, Ring und Muster schwarz eingelegt, 1,5 tief.
- **Verdeck:**
  - Wandstärke 3,5. Eine Gruppe „Verdeck“ mit Deckplatte, Seiten, Front, Leiste und Rändelschrauben.
  - Leiste 6 × 5 an der Kante zum Griffstock, X 55–333, Vorderseite R3, Enden R2.
  - Front mit Rosette in der Mitte.
- **Verdeckstützen:**
  - 3 mm, 25 hoch, innen an den Seitenteilen (Mitte, Y ≈ 110,8), 1,1 mm Abstand zur Wand.
  - Je Fuß eine M2,4 × 8 Sechskantschraube von außen durch den Fuß in den Boden.
  - Rändelschraube M2,4 waagrecht durch die Seitenwand (Z −14,25).
- **Griffstock:**
  - Oberer Körper als Kamm: Schlitze für Reihe 1 lang, für Reihe 2 kurz, beide unten offen.
  - Teil 2 als Kamm mit Boden: Schlitze für Reihe 3 und 4, 2 mm Boden hinter x −120.
  - Schlitze 5,4 breit. Abdeckung nur bis Schlitzende + 3 mm Auflage.
- **Hebel:** Holzteile, Alu und Knöpfe nach jp_schnitt.dxf. Alle 50 kollisionsfrei, in Ruhe und gedrückt.

## Fusion-Skripte

Einbinden über Dienstprogramme → Zusatzmodule → „+“ → Skript vom Gerät.

| Skript | Zweck |
|---|---|
| Hebel_Gelenke | legt die 50 Drehgelenke an (Achse, Anschlag Ruhe … gedrückt), fixiert das Kunststoffteil |
| Hebel_druecken | schaltet einen Hebel um (`HEBEL = "Hebel_R1_08"`) |
| Diskant_Ablauf | spielt den Ablauf vor: Verdeck weg, Griffstock aus, Reihen drücken, alles zurück; dauert ca. 5 min |

## Offen

- Explosion: Teil 2, Verdeckleiste und Stützen fahren noch nicht mit.
- Fusion-Drehbuch (Storyboard) per Hand war auf dem alten Laptop nicht machbar.
