# Diskant-Generator, Stand 22.11.2026

Neu sind die Griffstock-Varianten „kolben“ und „kolben12“ (verdeckt, Kunststoffzylinder unter den Knöpfen, `kolben_geo.py`) sowie „winkel“ (Alu-Winkel 46753 mit einer Achse Ø3 für alle Hebel, Lagerlaschen nach 4685/468401, Führung der Hebel nach ital.df.dxf, Spiralfedern und Gummi 031404; `winkel_geo.py`, `winkel_teile.py`, `lasche_4685_knick.dxf`). Die Winkel-Variante gibt es nur mit 20 mm Knopfabstand (`tastenliste_marie_winkel.txt`, `tastenliste_markus_winkel.txt`).

Die verdeckten Griffstöcke haben eine innen glatte 4-mm-Abdeckung nach dem Schnitt griffstock-alternativ DEZ 2018; die Hebel halten mindestens 2 mm Abstand zur Innenseite, Freischnitte gibt es keine mehr. `instrument.py` richtet die Knopfpyramide immer mittig am Korpus aus, abhängig von der konfigurierten Knopfanzahl.

Aufruf wie bisher: `INSTRUMENT=tastenliste_….txt python3 diskant.py BESASDES ausgabe.step`. Mit der Umgebungsvariablen `NUR=…` werden nur ausgewählte Teile ausgegeben, z. B. die Mechanik ohne Stimmstöcke.
