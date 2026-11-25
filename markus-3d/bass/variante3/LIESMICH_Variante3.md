# Bassteil Variante 3 „traditionell, flacher Einbau“

Grundlage ist das Markus-Bassteil (Gehäuse, flacher Boden, Deckel mit Rosetten und Leder, Füße, Hörner, F-Löcher, Begleiter-Stimmstock, Zwischendeckel mit Luftklappe, Lufttaste). Die Markus-Kammern sind durch die Helikon-Kammern in 4/5-Anordnung ersetzt (8 Kammern 98 lang mit Längsleisten 6 mm, quer liegende Zusatzkammer vorne rechts), die Knöpfe auf 9 + 8 reduziert (Reihe 1 in der Frontfläche, Reihe 2 mit Blende). Die Bassklappen-Öffnungen haben denselben Abstand zur Rückseite wie bei Markus.

Erzeugen: die übrigen Markus-Quellen (bass/) in denselben Ordner kopieren, dann `python3 teile.py rahmen.dxf Markus_front_v4.step front` (ebenso rueck, seite) und `python3 kammern.py boden2.dxf Ausgabe.step`.

Mechanik nach Markus-Muster unter dem flachen Boden: je Klappe eine Achse quer hinter der Klappe in Lagerböcken an der Bodenunterseite (Bass drei Reihen hinter den Bassklappen, Begleiter drei Reihen zwischen Bass- und Begleiterklappen, Zusatzklappe eigene Achse), Arm am Knopf-X nach unten zum Gabelgelenk, Schubstange Ø2 zum Knopf, Klappenhebel zur Klappenmitte, Rückdruckfeder je Achse. Armlänge je Achse so, dass der Klappenhub möglichst 9 mm beträgt. Die Luftklappe ist nach der Standardlösung (Alu-Lufthebel unter dem Zwischendeckel) gebaut.


## Stand 24.11.2026 abends (v4)

Der Boden ist 3,5 mm dick (Oberkante wie Markus), dadurch liegt die Mechanik 2,5 mm höher. Begleiter-Zwischenteil, Sole, Stimmstock und Begleiterklappen sind mittig ins Gehäuse gesetzt, beidseits schließt eine 5-mm-Leiste den Außenring. Die Luftklappenöffnung liegt wie beim Standard-Bassteil relativ zum Hebeldrehpunkt.

Achsen, Schenkel, Druckstangen und Klappenhebel sind Ø2,4 rund; die Gelenke bestehen aus dem Fusion-Gelenkkopf (kb_gelenkkopf.step) an der Druckstange und einer auf 8 mm gekürzten Gabel am Schenkel. Alle Arme sind 17,4 mm lang (Gelenk 1,5 mm über Gaze/Deckel), damit die Klappen beim Öffnen nicht auf die darunter laufenden Schubstangen schwenken; die Stangen laufen unter den Klappen ganz tief und steigen erst hinter der letzten Klappe im Weg.

Lagerbuchsen nach lagerbeide.dxf (einfach mit Hohlkehle, doppelt), Reihenabstand 3,2 mm, Achsen stecken von beiden Seiten in dieselbe Buchse; hinter einem Gelenkarm steht keine Buchse, die Rückdruckfeder sitzt immer auf der Buchsenseite. Die drei Bassachsen auf der Seite der Zusatzkammer (Klappen 6, 7, 8) sind leicht schräg: Buchse am Klappenende knapp hinter der Klappe, zum Knopf hin nach hinten aufgefächert. Die Zusatztaste hat ihren Arm direkt am Knopf, der Boden ist links der Zusatzkammer für die Buchse verlängert, der Gelenkkopf-Schaft dort auf 8 mm gekürzt.

Die Zusatzachse liegt knapp hinter der Zusatzklappe; ihre Buchse im offenen Bereich links der Zusatzkammer sitzt an einem Lagerblock, der an Bodenkante und Kammerleiste geleimt ist. Beim Zusatzgelenk zeigt der Gelenkkopf nach hinten (vom Knopf weg); die Stange macht ein senkrechtes U nach oben und läuft auf Z 16 über die gedrückte Lufttaste zum Knopf. Der Schenkel ist dafür 4,5 mm seitlich versetzt, der Gelenkkopf-Schaft auf 5 mm gekürzt, der Zusatzklappenhebel liegt neben der Taste. Geprüft wird auch mit gedrückter Lufttaste und mit jeder Basstaste einzeln gedrückt (Klappe dreht mit dem Klappenhebel um die Achse): keine Überschneidung. Farbgebung wie beim Standard-Bassteil 4/5 (gleiche Teileart, gleiche Farbe).
