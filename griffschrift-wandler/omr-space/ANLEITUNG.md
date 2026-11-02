# Eigener Erkennungsserver für „Noten scannen“ (vorbereitet, nicht aktiv)

Dieser Ordner enthält einen Erkennungsdienst, der Fotos oder PDFs gedruckter Noten mit homr in MusicXML umwandelt. Der Wandler kann ihn direkt aufrufen: Im Dialog „📷 Noten scannen …“ wählt man nach der Bestätigung „Foto oder PDF senden …“, und das Ergebnis wird gleich geladen. Solange im Wandler keine Serveradresse eingetragen ist, bleibt dieser Weg ausgeblendet, und es gibt nur die Links zu den fremden Diensten.

homr (Version 0.7) rechnet mit onnxruntime auf der CPU und braucht keine Grafikkarte. Im Test hier lag eine einzeilige Melodie auf zwei Kernen bei 6 bis 20 Sekunden, eine zweiseitige PDF bei etwa 20 Sekunden; ein sauber gesetztes Notenbild und ein absichtlich verschlechtertes Foto davon wurden fehlerfrei erkannt. homr liest nur Bilder; PDF-Seiten wandelt omr_core.py mit PyMuPDF in Bilder um und hängt die Takte aller Seiten hintereinander.

Es gibt zwei Fassungen. Für einen eigenen virtuellen Server liegt im Unterordner vps eine Fassung mit eigener Anleitung (vps/INSTALL.md). Diese Datei beschreibt die Fassung für Hugging Face.

## Hugging Face Space (kostenlos mit ZeroGPU)

Seit Juli 2026 laufen Gradio-Spaces auf der kostenlosen Hardware „CPU basic“ nur mit PRO-Abo (9 Dollar im Monat). Kostenlos bleibt „ZeroGPU“: laut Forum für Konten mit bestätigter E-Mail, die älter als 30 Tage sind, höchstens zwei solche Spaces. app.py ist dafür vorbereitet (die Funktion ist mit `spaces.GPU` gekennzeichnet, homr selbst rechnet trotzdem auf der CPU). Wie viel Rechenzeit ZeroGPU anonymen Besuchern, also Aufrufen aus dem Wandler, pro Tag gewährt, ist nicht genau bekannt; bei wenigen Benutzern sollte es reichen.

1. Auf huggingface.co anmelden (Konto mit bestätigter E-Mail, älter als 30 Tage).
2. „New Space“: Name zum Beispiel „notenerkennung“, SDK „Gradio“, Hardware „ZeroGPU“, Sichtbarkeit „Public“.
3. Unter „Files“ die vier Dateien app.py, omr_core.py, requirements.txt und README.md hochladen (README.md ersetzt die vorgegebene; sie legt Python 3.11 fest, das homr braucht).
4. Warten, bis oben „Running“ steht. Der erste Start dauert einige Minuten, weil Pakete installiert und beim ersten Aufruf die Modelle geladen werden. Im Space selbst mit einem Notenfoto ausprobieren.
5. Die Adresse des Space steht unter „⋯ → Embed this Space → Direct URL“, zum Beispiel https://NAME-notenerkennung.hf.space. Diese Adresse im Wandler in index.html bei `const OMR_URL="";` eintragen und pushen.

Ein Space, der länger nicht benutzt wurde, schläft und braucht beim nächsten Aufruf einige Minuten zum Aufwachen; der Wandler zeigt das an. Die Datei wird auch dann nur nach dem Häkchen „Verstanden: die Datei wird an einen fremden Dienst gesendet“ gesendet.

Die Anbindung im Wandler ist mit der Gradio-Schnittstelle geschrieben (Datei hochladen, Aufruf „scan“, Ergebnis abholen) und mit einem nachgebauten Server geprüft, aber noch nicht mit einem echten Space. Falls der Browser den Aufruf wegen fremder Herkunft (CORS) ablehnt oder Gradio eine andere Schnittstelle erwartet, bitte die Fehlermeldung aus dem Scan-Dialog weitergeben.
