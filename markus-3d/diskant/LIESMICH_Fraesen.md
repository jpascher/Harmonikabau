# Fräs- und Fertigungsunterlagen Diskantteil (Marie / Markus) — Stand 19.11.2026

## Ablauf in zwei Schritten

Alles entsteht aus der Tastenliste (`tastenliste_marie.txt`, `tastenliste_markus.txt`). Ändern sich Tastenzahl, Knopfabstand, Plattenlängen, Versatz oder Neigung, genügt ein Aufruf des Generators: `python3 erzeuge.py <tastenliste> <Ausgabeordner> <Version>`. Er rechnet das ganze Modell samt Kollisionsprüfung und legt im Ausgabeordner alle Fertigungsdateien ab; der Ordner wird nach `C:\ACAD-DWG\0-Markus-3D\Fraesen\<Name>\` kopiert. Danach in Fusion das Skript **Fraesen_Aktualisieren** starten: es liest jede `fraes_liste.json` unter `Fraesen\`, setzt im Dokument `Stimmstoecke_<Name>_Fraesen` alle Frästeile neu ein und legt die CAM-Einrichtung neu an. Die Werkzeugwege werden dabei bewusst nicht berechnet (in Fusion: Rechtsklick auf „Setups“ → „Generieren“).

## Ordner je Instrument

`Fraesen\<Name>\` enthält die Kollisionsliste des Modells (`…kollision.json`), die Liste für Fusion (`fraes_liste.json`, Teilenamen und Ordner) und drei Unterordner. Die Unterordner `Stimmstoecke` und `Griffstock` liegen direkt vor, die übrigen Teile kommen als `Teile_<Name>.zip` und werden im Ordner des Instruments entpackt (ergibt `Teile\`).

## Stimmstöcke (`Stimmstoecke\`)

Je Reihe ein STEP und ein DXF des Stimmstocks in Fräslage: Vorder- und Rückseite flach aufgeklappt, Außenfläche oben (Z = 0), Kammern nach unten, Keil als eigenes Teil. Kammerbreite, -tiefe und Abstufung folgen Jakobs Regeln aus der Plattenlänge; die DXF hat je Frästiefe einen Layer. In Fusion entstehen je Stock zwei Setups (Vorder- und Rückseite) mit je zwei adaptiven Bearbeitungen nach Jakobs Vorlage, 2,5-mm-Schaftfräser; die Außenkontur wird nicht gefräst (Grenze = Außenkontur der Bodenfläche, Werkzeug außerhalb, Rohteil seitlich 0 mm).

## Griffstock (`Griffstock\`)

Teil 1 (Körper) und Teil 2 (Kamm mit Boden) je als STEP und DXF in Fräslage: Deckelseite unten, Oberseite bei Z = 0. Gefräst wird von oben und von der Deckelseite, die Außenkontur nie. Der 6-mm-Fräser schruppt mit 0,5 mm Aufmaß auf den Böden und schlichtet danach jede waagrechte Fläche einzeln mit 0,5 mm (Layer `Flaeche_oben_…` bzw. `Flaeche_unten_…`); seitlich läuft er um 6 mm über die Fläche hinaus, aber nur dort, wo daneben nichts höher steht, damit die Stufen erhalten bleiben. Der Innenraum zwischen den Schenkeln wird als Rechteck ausgeschnitten (Layer `Innenraum_Rechteck_durch`), die Deckel-Ausnehmung von Teil 1 vollflächig von der Deckelseite ausgeräumt. Die Kanäle für die Hebel werden einzeln als Rechtecke gefräst, zur offenen Seite um den Fräserdurchmesser verlängert, damit an der Mündung keine Rundung bleibt: Teil 1 mit dem 3,175-mm-Fräser (18,9 mm tief), Teil 2 mit dem 2,5-mm-Fräser (14,6 mm tief). Alle Maße stehen zusätzlich in `Griffstock_<Name>_fraesen.json`.

## Alu-Hebel (`Teile\Alu_Hebel_abgewickelt\`)

Je Hebel eine DXF der Abwicklung (gestreckt), Blech 2 mm. Layer `Kontur` ist die Schnittkontur, `Bohrung_D6_Messingbuchse` die Bohrung für die eingepresste Messingbuchse (außen 6, innen 2, 2 mm lang), `Biegelinie` die Biegestellen mit Winkel und Richtung im Layer `Text`. Längsbiegungen (zur Klappe hin) sind mit „zur Federseite (+X)“ bzw. „von der Federseite weg (-X)“ bezeichnet, bei den Reihen 1 und 2 kommen die beiden gegenläufigen Biegungen der Kröpfung dazu. Die Abwicklung ist für Biegen über Halbrund 5 mm gerechnet (Innenradius 2,5 mm, neutrale Faser k = 0,4). Der lange Arm läuft bis zur Klappe auf halbe Höhe zu, das flache Endstück auf der Klappe ist 2 mm hoch.

## Holzteile (`Teile\Holzteile\`)

Je Reihe eine DXF, da alle Holzteile einer Reihe gleich sind; die Stückzahl steht im Dateinamen und im Text. Neben der Kontur (`Kontur_Holz`) gibt es Ein- und Ausfahrbogen R3 an der längsten geraden Kante (`Einfahren_R3`, `Ausfahren_R3`) und, nur als Referenz, das darunterliegende Alu-Teil mit der Ø6-Bohrung (`Alu_Referenz`, `Bohrung_Referenz_D6`). Holzdicke 4,8 mm.

## Kunststoffteil (`Teile\`)

`Kunststoffteil_<Name>_fraesen.step` ist das Frästeil in Fräslage (Lagerseite oben, Z = 0): die Abschrägungen bleiben als Rohblock stehen und werden nicht gefräst, gefräst werden nur Nuten, Taschen, Durchbrüche, Schraubenlöcher und Rastmarken. `…_fertig_in_Fraeslage.step` ist das Fertigteil in derselben Lage zum Vergleich. Die DXF zeigt die Taschen je Tiefe auf eigenen Layern (`Tasche_n_bis_T…`) und die Längsbohrungen Ø2 unter den Lagernuten als Bahnen für den Kugelkopffräser Ø2 mit dünnem Hals (`Kugelfraeser_D2_Laengsbahn_T<Achstiefe>`, über die ganze Länge, 3 mm Ein-/Ausfahren an den Stirnseiten). Die Rastmarken für die Federschenkel liegen in Flucht mit jedem Hebel an der hinteren Kante der Auflagefläche, 2 mm lang, 0,5 mm tief, 0,8 mm breit; die Federschenkel werden in der Praxis bis dorthin verlängert (im Modell nicht).

## Bekannte offene Punkte

Von den 50 Rastmarken schneiden 39 ins Kunststoffteil, 11 treffen an ihrer Stelle kein Material (liegen in einer vorhandenen Öffnung); bitte im Modell prüfen. Die Auflagefläche der Federschenkel liegt etwa 5 mm unter der Lagerseite.
