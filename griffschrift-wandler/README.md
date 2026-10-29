# Griffschrift-Wandler

[Deutsch](#deutsch) · [English](#english)

Online: https://jpascher.github.io/Harmonikabau/griffschrift-wandler/ (English: `?lang=en`)

## Deutsch

Der Griffschrift-Wandler setzt Normalnoten in Griffschrift für die Steirische Harmonika um. Er läuft vollständig im Browser: `index.html` öffnen, eine MusicXML-Datei (.musicxml, .xml, .mxl) oder eine MIDI-Datei (.mid) hineinziehen, und die Seite zeigt die Griffschrift nach dem System Max Rosenzopf (ohne Schlüssel und Vorzeichen, Unterstrich = Druck, Bass mit Knopfbuchstaben) oder wahlweise eine Griff-Tabelle mit Reihe und Knopfnummer.

Die Seite ist zweisprachig. Der Knopf „English“ oben rechts schaltet auf Englisch um, „Deutsch“ wieder zurück; die Wahl wird im Browser gemerkt. Beim ersten Besuch richtet sich die Sprache nach der Spracheinstellung des Browsers, mit `?lang=de` oder `?lang=en` an der Adresse lässt sie sich fest vorgeben.

Die Griffwahl geht vom Bass aus. Das Programm übernimmt Akkordsymbole aus der Datei oder schätzt die Stufen aus der Melodie, bei MIDI-Dateien aus der Bassspur (mit Erkennung des Wechselbasses). Dann sucht es über das ganze Stück die Folge von Grundbass und Begleiter mit den kürzesten Wegen; ein Balgwechsel am selben Knopf zählt als besonders flüssig. Die Melodie samt Doppelgriffen folgt der Balgrichtung des Basses, Töne gegen die Bassrichtung werden rot markiert.

Wählbar sind die Stimmung (alle acht 4-Reiher-Stimmungen von F-B-Es-As bis C-F-B-Es), der Diskant (450 Öllerer und Schaborak, 446 Normal) und der Bass. Normalbass und Schaborak stammen aus der Akkordsuche für Handharmonika auf hyperbox.org; dazu kommen die Bassbelegungen aus dem Tonbelegungs-Konfigurator von Strasser (harmonika.com), darunter Übergangsbässe, Mollbässe und 1 Mollbass in mehreren Ausführungen sowie die Systeme Klassisch, Strohmeier, Glaser, Pichler, Hubmann, Senekovic, Schaborak, Kalcher, Egger, Svoboda, Zweverink, Guggi und Lenz. Das Stück wird bei Bedarf automatisch in eine gut spielbare Lage transponiert.

Für die Begleiter gibt es Mollbelegungen. Neben der Belegung, wie sie gebaut ist, legt eine Vorlage die Begleiter der Hilfsbassreihe für die Paralleltonarten auf Moll (bei G-C-F-B d-, a- und e-Moll). Weil Steirische oft individuell gebaut werden, lässt sich außerdem jeder Begleiter im Bassbild per Klick zwischen Dur und Moll umschalten, getrennt für Zug und Druck; die eigene Belegung wird je Bass-System im Browser gemerkt und die Griffwahl sofort neu gerechnet. Fehlt für eine Stufe der passende Bass oder Begleiter, zeigt die Statuszeile das an, und die automatische Transposition sucht eine Lage, in der möglichst alles vorhanden ist. Der „Moll-Walzer“ unter den Beispielen zeigt den Unterschied zwischen den Belegungen.

Wiederholungen werden bei MIDI-Dateien erkannt und mit Wiederholungszeichen und 1./2. Schluss dargestellt; wahlweise lässt sich das Stück ausgeschrieben anzeigen. Beim Abspielen werden auch im wiederholten Durchgang die gespielten Noten markiert. Im Bearbeiten-Modus kann man die Balgrichtung einzelner Griffe oder ganzer Bassabschnitte umschalten und Zeilenumbrüche setzen; die Zahl der Takte je Zeile ist ebenfalls wählbar. Der Drucken-Knopf gibt nur die Griffschrift auf A4 aus. Zum Üben lassen sich Takte über die Taktnummer (Umschalt-Klick für einen Bereich) oder über „Takt von – bis“ auswählen; „Abspielen“ spielt dann nur diese Takte, mit „Schleife“ immer wieder.

### Gescannte Noten und PDF

Gescannte Noten oder PDF-Dateien kann der Wandler nicht direkt lesen, denn dafür braucht es eine Notenerkennung (OMR), die aus dem Bild wieder Noten macht. Der Weg führt über MusicXML: Mit dem kostenlosen, quelloffenen Programm Audiveris (Windows, macOS, Linux; Installationsprogramm auf der Audiveris-Seite bei GitHub) öffnet man die PDF- oder Bilddatei, startet die Erkennung („Book“ → „Transcribe“) und exportiert das Ergebnis als MusicXML (.mxl). Diese Datei zieht man in den Wandler.

Keine Notenerkennung ist fehlerfrei. Falsche Notenwerte, vergessene Vorzeichen, Bögen, die als Haltebögen gelesen werden, und vermischte Stimmen sind typisch. Am bequemsten öffnet man die .mxl-Datei deshalb vorher im kostenlosen MuseScore, hört sie an, bessert aus und exportiert wieder als MusicXML. Gute Ergebnisse bringen Scans mit 300 dpi in Graustufen, gerade aufgelegt und ohne Schatten; saubere Druckausgaben und einstimmige Melodien mit Akkordsymbolen werden am besten erkannt. MuseScore bietet über musescore.com außerdem einen eigenen PDF-Import an, der ein Konto voraussetzt.

Die Seite enthält eine ausführliche Hilfe (Knopf „Hilfe“) mit allen Funktionen, der Erklärung der Griffschrift und diesem Ablauf.

### Beispiele und offene Punkte

Die Beispiele im Ordner `midi/` und in der Seite sind Volkstänze von volksmusik.cc, Satz: Franz Fuchs; die Herkunft wird beim Laden angezeigt. Dazu kommen eigene Beispiele (Ländler mit Terzen, Stille Nacht, Freude schöner Götterfunken, Alle meine Entchen).

Offen sind noch die Notation der Bass-Zusatzknöpfe (derzeit Tonname mit Stern), die Notenkopfrichtung bei Doppelgriffen nach Fuchs und die Erkennung von Auftakten in MIDI-Dateien.

## English

The Griffschrift-Wandler turns standard notation into Griffschrift, the button tablature for the Styrian accordion (Steirische Harmonika). It runs entirely in the browser: open `index.html`, drop in a MusicXML file (.musicxml, .xml, .mxl) or a MIDI file (.mid), and the page shows the Griffschrift in the system of Max Rosenzopf (no clef or accidentals, underline = push, bass with button letters), or alternatively a grip table with row and button number.

The page is bilingual. The “English” button at the top right switches to English, “Deutsch” switches back; the choice is remembered in the browser. On the first visit the language follows the browser's language setting, and `?lang=de` or `?lang=en` in the address fixes it.

Fingering is planned from the bass. The program takes chord symbols from the file or estimates the harmonic degrees from the melody, and for MIDI files from the bass track (alternating bass is recognised). Then it searches the whole tune for the sequence of bass notes and chord buttons with the shortest moves; a bellows change on the same button counts as especially smooth. The melody, double stops included, follows the bellows direction of the bass, and notes played against the bass direction are marked red.

You can choose the tuning (all eight 4-row tunings from F-B-Es-As to C-F-B-Es), the treble layout (450 Öllerer and Schaborak, 446 standard) and the bass layout. Standard bass and Schaborak come from the chord finder for button accordion on hyperbox.org; in addition there are the bass layouts from the Strasser key-layout configurator (harmonika.com), including transition basses, minor basses and 1 minor bass in several versions, as well as the Classic, Strohmeier, Glaser, Pichler, Hubmann, Senekovic, Schaborak, Kalcher, Egger, Svoboda, Zweverink, Guggi and Lenz systems. If needed, the tune is transposed automatically into a key that lies well on the instrument.

Minor chord layouts are available for the chord buttons. Besides the layout as built, a template switches the chord buttons of the auxiliary bass row for the relative minor keys to minor (on G-C-F-B: D, A and E minor). Because Styrian accordions are often built individually, every chord button can also be switched between major and minor by clicking it in the bass diagram, separately for pull and push; your own layout is remembered in the browser for each bass system, and the fingering is recalculated at once. If the matching bass or chord button for a degree is missing, the status line shows it, and automatic transposition looks for a key in which as much as possible is available. The “Minor waltz” among the examples shows the difference between the layouts.

In MIDI files, repeats are recognised and shown with repeat signs and first/second endings; the tune can also be shown written out in full. During playback, the notes being played are highlighted in the repeated pass as well. In edit mode you can switch the bellows direction of single grips or whole bass sections and set line breaks; the number of measures per line can also be chosen. The Print button outputs only the Griffschrift on A4. For practice, select measures by their measure number (Shift-click for a range) or with “Practice: measures from – to”; “Play” then plays only those measures, over and over with “Loop”.

### Scanned music and PDF

The converter cannot read scanned music or PDF files directly, because that needs optical music recognition (OMR) to turn the image back into notes. The way there is through MusicXML: open the PDF or image file in Audiveris, a free, open-source program (Windows, macOS, Linux; the installer is on the Audiveris page on GitHub), run recognition (“Book” → “Transcribe”) and export the result as MusicXML (.mxl). Drag that file into the converter.

No music recognition is perfect. Wrong note values, missing accidentals, slurs read as ties and mixed-up voices are typical. The easiest fix is to open the .mxl file in the free MuseScore first, listen to it, correct it and export it again as MusicXML. Scans at 300 dpi in greyscale, lying straight and without shadows give good results; clean printed music and single-line melodies with chord symbols are recognised best. MuseScore also offers its own PDF import through musescore.com, which requires an account.

The page has a detailed built-in help (the “Help” button) covering all functions, how to read Griffschrift and this workflow.

### Examples and open points

The examples in the `midi/` folder and in the page are Austrian folk dances from volksmusik.cc, arrangements by Franz Fuchs; the source is shown when a tune is loaded. There are also some own examples (Ländler in thirds, Silent Night, Ode to Joy, Alle meine Entchen).

Still open: the notation of extra bass buttons (currently the note name with a star), the note-head direction in double stops as described by Fuchs, and detecting pickups (anacruses) in MIDI files.
