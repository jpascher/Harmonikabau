# Bassteil Marie 38 cm (Helikon in zwei Ebenen mit Verbindungskanal)

Erzeugen: `python3 kammern38.py Bassteil_Marie38_vN.step [Knopfteilung]` (Standard 19 mm) im Ordner markus-3d/bass/marie38.
kammern38.py baut Gehaeuse, Kammern, Etage, Formteile, Begleiterstimmstoecke, Luftklappe (Lufttaste gedrueckt dargestellt)
und bindet am Ende mech38.py ein (Achsen, Lagerboecke, Gelenke, Schubstangen, Klappenhebel, Rueckdruckfedern) samt
Kollisionspruefung in Ruhe und mit jeder gedrueckten Taste. Eingaben: bassteil-alles.dxf (Hauptzeichnung), hb_trennwaende.dxf
(Formteile), kb_lufthebel.dxf, Markus-Rahmen (Markus_*_v4.step, rahmen.dxf, verdeck.dxf, horn.dxf, fus2.dxf, blende.dxf),
lagerbeide.dxf und kb_gelenkkopf.stp (Mechanik wie Markus/V3). Stand 25.11.2026: v16.

Untergeordnet, erst bei der Fertigung bereinigen: keine Mittelstuetzen an langen Achsen, Haube auf der Abdeckplatte ueber
dem Lufthebel, Klappenabsenkung pruefen, vier Lagerboecke mit weniger als 0,6 mm Luft.
