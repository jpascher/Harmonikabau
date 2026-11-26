# Regeln für die freie Bassbelegung (abgeleitet aus den gebauten Modellen, Stand 26.11.2026)

Diese Regeln sind aus den fünf Bassteilen abgeleitet, deren Mechanik nach demselben Baukasten gerechnet wird: Eduard (`eduard/mech_eduard.py`), Gerhard (`gerhard/mech_gerhard.py`), Andreas (`andreas/mech_andreas.py`), 340 schmal (`v340/mech_v340.py`) und Marie 38 (`marie38/mech38.py`). Der Markus-Bass (`mechanik.py`) liest seine Mechanik noch direkt aus der Zeichnung und folgt denselben Bauteilen, aber nicht derselben Tabelle; er wird später auf die Tabelle umgestellt. Neue Bauteile kommen nicht vor; jede Belegung muss sich aus Achse, Lagerbock, Schenkel mit Gelenk, Schubstange, Klappenhebel und Feder zusammensetzen lassen.

## 1. Grundregel: Knöpfe = Klappen + Kopplungen

Jede Klappe hat genau eine Achse, und jede Achse bedient genau eine Klappe. Eine Kopplung entsteht, indem ein zweiter Knopf mit einem eigenen Schenkel auf dieselbe Achse greift. Damit ist die Knopfzahl immer die Zahl der Klappen plus die Zahl der Kopplungen. Alle fünf Modelle erfüllen das genau:

| Modell | Klappen (= Achsen) | davon Grundbass / Begleiter | Kopplungen | Knöpfe |
|---|---|---|---|---|
| Eduard | 18 | 10 / 8 | 3 | 21 |
| Gerhard | 24 | 12 / 12 | 5 | 29 |
| Andreas | 21 | 13 / 8 | 2 | 23 |
| 340 schmal | 17 | 9 / 8 | 1 | 18 |
| Marie 38 | 18 | 10 / 8 | 1 | 19 |

Eine Begleiterklappe in Tripelteilung öffnet drei Kammern, zählt hier aber als eine Klappe. Die Tonhöhe spielt für die Mechanik keine Rolle; es zählt nur, welcher Knopf welche Klappe öffnet.

## 2. Was je Modell fest ist und was frei wählbar wird

Fest sind beim gewählten Modell das Gehäuse, die Kammern und damit Zahl, Art und Lage der Klappen (Grundbass-Klappen in ihren Reihen, Begleiterklappen links und rechts), die Knopfreihen mit ihrer Knopfzahl und die Lage der Achsgassen. Frei wird die Belegung: welcher Knopf welche Klappe öffnet, welche Knöpfe gekoppelt sind, und die Knopfteilung innerhalb der Grenzen des Modells. Die Gesamtknopfzahl ergibt sich aus Klappen plus gewünschten Kopplungen; passt sie nicht zu den Knopfreihen des Modells, wird das gemeldet.

## 3. Achsgassen

Die Achsen liegen alle parallel zu den Knopfreihen (in X) unter dem Boden, in Gassen zwischen den Klappenreihen. Jedes Modell hat je nach Kammerbild 7 bis 12 solcher Gassen: hinter den Begleiterklappen bzw. am Rückteil, zwischen Begleiter- und hinterer Grundbassreihe, zwischen den beiden Grundbassreihen und vorne am Bodenrand. Benachbarte Gassen liegen meist 3,3 bis 5,5 mm auseinander (Gerhard 3,3–4,0 mm, Eduard 5 mm, Andreas 5–5,5 mm), mindestens so viel, dass Lagerböcke (8,6 mm lang, 4–5,2 mm breit) und Federn nebeneinander Platz haben; die Böcke benachbarter Gassen werden gegeneinander versetzt.

In einer Gasse dürfen mehrere Achsen hintereinander liegen, solange sich ihre X-Bereiche nicht überschneiden. Der X-Bereich einer Achse reicht vom äußersten Schenkel bis zum Klappenhebel, plus je einem Lagerbock an beiden Enden (Achsende 1 mm vor der Bockschraube). Zwei Achsen in einer Linie dürfen sich einen gemeinsamen Bock teilen; sonst brauchen ihre Enden einen Abstand von einer Bocklänge. Wo ein Bock an Eckblock, Eckwinkel, Deckelauflage oder Luftklappe stößt, wird er auf eine feste Lage gesetzt (`FEST_BOCK`).

## 4. Klappenhebel

Der Klappenhebel geht von der Achse unter die Klappenmitte. Bei Grundbass-Klappen läuft er gerade in Y, die Achse muss also in einer Gasse liegen, die direkt an die Klappenreihe grenzt (davor oder dahinter). Bei den Begleiterklappen an den Seiten läuft er schräg vom Achsende zur runden Klappe; diese Achsen liegen deshalb in den Gassen, die nahe an den Begleiterklappen vorbeiführen, und reichen bis an die Seite heran. Läuft ein Hebel über eine fremde Klappe oder ein Herz, bekommt er einen Umweg über einen Zwischenpunkt (`HEBEL_WEG`, bisher dreimal bei Andreas). Je Achse gibt es genau einen Klappenhebel.

## 5. Schenkel und Schubstangen

Jeder Knopf hat einen Knopfzylinder, eine Schubstange Ø 2,4 bis zum Gelenkkopf und einen Schenkel Ø 2,4 vom Gelenk zur Achse. Die Stange läuft vom Knopf senkrecht hoch auf Gelenkhöhe knapp über dem Deckel und dann in Y bis zur Gasse ihrer Achse; sie unterquert dabei alle Klappen. Der Schenkel sitzt normal auf der X-Lage des Knopfes. Stehen zwei Stangen in derselben X-Lage im Weg (Knopf einer anderen Reihe dahinter), wird der Schenkel um 4 bis 5 mm seitlich versetzt (in den Tabellen z. B. `K3(1) - 4.0`). Im Bereich der Lufttaste weichen Stangen seitlich oder in der Höhe aus (Luftband).

Eine Kopplung bedeutet zwei Schenkel auf einer Achse. Die gekoppelten Knöpfe können in verschiedenen Reihen und weit auseinander liegen (bei Eduard bis 150 mm); die Achse reicht dann über beide Schenkel und den Hebel.

## 6. Federn und Prüfung

Jede Achse hat eine Rückdruckfeder, normalerweise am freien Ende zwischen Gelenkhebel und Lagerbock, links- oder rechtsgewunden je nach Seite; wo dort kein Platz ist, als Torsionsfeder direkt neben dem Arm (`TORSION`) oder an einem bestimmten Arm (`FEDER_ARM`). Jede fertige Belegung wird wie bisher in Ruhe, mit jeder einzeln gedrückten Taste und mit allen gedrückten Tasten auf Kollisionen geprüft.

## 7. Wie die automatische Zuordnung arbeiten soll

Die Belegung wird in einer Datei je Instrument festgelegt (wie die Tastenliste beim Diskant): Modell, Knopfteilung, und für jeden Knopf die Klappe, die er öffnet. Zwei Knöpfe mit derselben Klappe sind gekoppelt. Daraus rechnet das Programm:

1. Prüfen, ob die Zahlen aufgehen (jede Klappe mindestens ein Knopf, Knöpfe = Klappen + Kopplungen, Knopfzahl passt zu den Reihen des Modells).
2. Für jede Achse den X-Bereich und die Gassen, aus denen der Hebel die Klappe erreicht (Regel 4).
3. Die Achsen so auf die Gassen verteilen, dass sich in keiner Gasse X-Bereiche überschneiden; bevorzugt die Gasse, die der Klappe am nächsten liegt, lange und gekoppelte Achsen zuerst.
4. Schenkel versetzen, wo Stangen in derselben X-Lage zusammentreffen (Regel 5), Böcke setzen, Federn an das freie Ende.
5. Mit der vorhandenen Mechanik-Berechnung bauen und prüfen.

Geht eine Belegung nicht auf, meldet das Programm, welche Achse in keiner Gasse Platz findet oder welche Teile kollidieren, statt ein falsches Modell zu liefern. Die bestehenden Modelle dienen als Prüfstein: ihre heutige Belegung, als Datei geschrieben, muss wieder dieselbe Mechanik ergeben.

## 8. Stand der Umsetzung

Umgesetzt in `belegung.py` (Prüfung und Verteilung mit Rücksprung, wenn die erste Verteilung nicht aufgeht) für den Bassteil 340 schmal: Klappen, Gassen und Knopflagen stehen in `v340/belegung_modell_340.json`, die Belegung wie gebaut in `v340/belegung_340.txt`. Diese ergibt genau die gebauten Gassen. Eine Probe mit anderer Kopplung (Knopf 2.1 auf H077 statt H218) und getauschten Begleitern rechts wurde vollständig gebaut: Mechanikprüfung und Kollisionsprüfung ohne neue Berührungen. Eine Belegung mit V077 auf Knopf 1.3 wird richtig als nicht baubar gemeldet. Die Webseite `docs/bass_belegung.html` macht dieselbe Prüfung im Browser. Die anderen Bassteile werden nach demselben Muster umgestellt.

