# Erkennungsserver auf dem eigenen virtuellen Server (hyperbox.org)

Dieser Ordner enthält einen kleinen Webdienst, der Fotos oder PDFs gedruckter Noten mit homr in MusicXML umwandelt. Der Griffschrift-Wandler schickt die Datei nach dem Häkchen „Verstanden“ dorthin, fragt alle drei Sekunden nach dem Stand und lädt das Ergebnis direkt. Die Erkennung läuft als Auftrag, weil Cloudflare Antworten nach 100 Sekunden abbricht und eine Seite ein bis zwei Minuten brauchen kann. Hochgeladene Dateien liegen nur während der Erkennung in einem Temp-Ordner, Ergebnisse werden nach dem Abholen gelöscht, spätestens nach 15 Minuten. Es läuft immer nur eine Erkennung gleichzeitig, höchstens zehn Aufträge warten.

Voraussetzungen: Linux (die Befehle gelten für Debian oder Ubuntu), Python 3.11 oder 3.12, mindestens 4 GB Arbeitsspeicher (besser 8 GB; nachsehen mit `free -h`), einige GB freier Platz für Python-Pakete und Modelle. Eine Grafikkarte ist nicht nötig.

## 1. Programm einrichten (als root)

    apt install python3 python3-venv nginx
    useradd --system --home /opt/griffschrift-omr --create-home omr
    cp server.py requirements.txt /opt/griffschrift-omr/
    cd /opt/griffschrift-omr
    sudo -u omr python3 -m venv venv
    sudo -u omr venv/bin/pip install -r requirements.txt

Einmal von Hand testen; beim ersten Aufruf lädt homr seine Modelle herunter, das dauert:

    sudo -u omr venv/bin/homr /pfad/zu/einem/notenfoto.png

Danach liegt neben dem Foto eine Datei mit der Endung .musicxml.

## 2. Als Dienst starten

    cp griffschrift-omr.service /etc/systemd/system/
    systemctl daemon-reload
    systemctl enable --now griffschrift-omr
    curl http://127.0.0.1:8701/health

Die Antwort muss `{"ok":true,"queue":0}` sein. Der Dienst hört nur auf dem Rechner selbst (127.0.0.1), nach außen geht es über nginx.

## 3. nginx und Adresse

In Cloudflare einen DNS-Eintrag `scan` anlegen (Typ A, IP des virtuellen Servers, Proxy eingeschaltet). Dann:

    cp nginx-scan.conf /etc/nginx/sites-available/scan
    ln -s /etc/nginx/sites-available/scan /etc/nginx/sites-enabled/
    nginx -t && systemctl reload nginx

nginx-scan.conf hört auf Port 80. Steht Cloudflare unter „SSL/TLS“ auf „Flexibel“, genügt das. Steht es auf „Voll“ oder „Voll (strikt)“, braucht nginx zusätzlich Port 443 mit einem Zertifikat, am einfachsten ein Ursprungszertifikat aus Cloudflare (SSL/TLS → Ursprungsserver) in einem zweiten server-Block mit `listen 443 ssl;`. Prüfen im Browser: https://scan.hyperbox.org/health

## 4. Im Wandler einschalten

In griffschrift-wandler/index.html die Zeile `const OMR_URL="";` ändern in

    const OMR_URL="https://scan.hyperbox.org";

und pushen. Dann erscheint im Dialog „📷 Noten scannen …“ nach dem Häkchen die Zeile „Erkennungsserver des Wandlers (scan.hyperbox.org)“. Zugriffe erlaubt der Dienst nur von https://jpascher.github.io (Einstellung `OMR_ORIGINS` in der .service-Datei; mehrere Adressen mit Komma trennen).

## Pflege

Protokoll ansehen: `journalctl -u griffschrift-omr -f`. homr aktualisieren: `sudo -u omr /opt/griffschrift-omr/venv/bin/pip install -U homr`, danach `systemctl restart griffschrift-omr`. Grenzen lassen sich in der .service-Datei mit `Environment=` setzen: `OMR_MAX_MB` (Dateigröße, Vorgabe 20), `OMR_TIMEOUT` (Sekunden je Erkennung, Vorgabe 600), `OMR_MAX_QUEUE` (wartende Aufträge, Vorgabe 10). `MemoryMax=6G` begrenzt den Arbeitsspeicher; bei einem Server mit 4 GB auf 3G senken.

homr steht unter der AGPL-3.0: Wer den Dienst öffentlich anbietet, muss den Quelltext zugänglich machen. Das ist erfüllt, solange server.py hier im öffentlichen Repository liegt und homr unverändert verwendet wird.
