# Griffschrift-Wandler

Der Griffschrift-Wandler setzt Normalnoten in Griffschrift für die Steirische Harmonika um. Er läuft vollständig im Browser: `index.html` öffnen, eine MusicXML-Datei (.musicxml, .xml, .mxl) oder eine MIDI-Datei (.mid) hineinziehen, und die Seite zeigt die Griffschrift nach dem System Max Rosenzopf (ohne Schlüssel und Vorzeichen, Unterstrich = Druck, Bass mit Knopfbuchstaben) oder wahlweise eine Griff-Tabelle mit Reihe und Knopfnummer.

Die Griffwahl geht vom Bass aus. Das Programm übernimmt Akkordsymbole aus der Datei oder schätzt die Stufen aus der Melodie, bei MIDI-Dateien aus der Bassspur (mit Erkennung des Wechselbasses). Dann sucht es über das ganze Stück die Folge von Grundbass und Begleiter mit den kürzesten Wegen; ein Balgwechsel am selben Knopf zählt als besonders flüssig. Die Melodie samt Doppelgriffen folgt der Balgrichtung des Basses, Töne gegen die Bassrichtung werden rot markiert.

Wählbar sind die Stimmung (G-C-F-B, A-D-G-C, B-Es-As-Des, F-B-Es-As), der Diskant (450 Öllerer und Schaborak, 446 Normal) und der Bass (Schaborak 16 bis 22 Knöpfe, Normalbass). Die Belegungen stammen aus der Akkordsuche für Handharmonika auf hyperbox.org. Das Stück wird bei Bedarf automatisch in eine gut spielbare Lage transponiert.

Wiederholungen werden bei MIDI-Dateien erkannt und mit Wiederholungszeichen und 1./2. Schluss dargestellt; wahlweise lässt sich das Stück ausgeschrieben anzeigen. Im Bearbeiten-Modus kann man die Balgrichtung einzelner Griffe oder ganzer Bassabschnitte umschalten und Zeilenumbrüche setzen; die Zahl der Takte je Zeile ist ebenfalls wählbar. Der Drucken-Knopf gibt nur die Griffschrift auf A4 aus.

Die Beispiele im Ordner `midi/` und in der Seite sind Volkstänze von volksmusik.cc, Satz: Franz Fuchs; die Herkunft wird beim Laden angezeigt. Dazu kommen eigene Beispiele (Ländler mit Terzen, Stille Nacht, Freude schöner Götterfunken, Alle meine Entchen).

Offen sind noch die Notation der Bass-Zusatzknöpfe (derzeit Tonname mit Stern), die Notenkopfrichtung bei Doppelgriffen nach Fuchs und die Erkennung von Auftakten in MIDI-Dateien.

## English summary

The Griffschrift converter turns standard notation into tablature (Griffschrift) for the Styrian button accordion. It runs entirely in the browser: open `index.html` and drop in a MusicXML or MIDI file. The fingering is planned from the bass side. Chord degrees are taken from chord symbols, from the bass track of a MIDI file or estimated from the melody, and the bass buttons are chosen for the shortest hand movements across the whole piece. The melody, including double stops, follows the bellows direction of the bass. Tunings, treble layouts (450 Öllerer/Schaborak, 446 normal) and bass layouts (Schaborak 16–22 buttons, normal bass) are selectable. The bundled folk dance examples come from volksmusik.cc, arrangements by Franz Fuchs.
