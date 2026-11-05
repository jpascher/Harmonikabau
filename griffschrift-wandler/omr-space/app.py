# Notenerkennung (OMR) für den Griffschrift-Wandler als Hugging Face Space (Gradio): Foto oder PDF -> MusicXML mit homr
# (https://github.com/liebharc/homr, AGPL-3.0). Der Wandler ruft die Funktion „scan“ über die Gradio-Schnittstelle auf
# (/gradio_api/upload, /gradio_api/call/scan). homr rechnet mit onnxruntime auf der CPU, eine Grafikkarte braucht es nicht.
import os, sys, json, tempfile
import gradio as gr
from omr_core import scan_file
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "gsleser"))

# homr rechnet auf der CPU. Auf ZeroGPU-Hardware läuft normaler Code im CPU-Teil des Space ohne Rechenzeit-Kontingent;
# nur mit @spaces.GPU gekennzeichnete Funktionen verbrauchen das (für anonyme Aufrufe sehr knappe) GPU-Kontingent.
# ZeroGPU verlangt beim Start mindestens eine solche Funktion – daher eine leere, die nie aufgerufen wird.
try:
    import spaces
    @spaces.GPU(duration=1)
    def _gpu_platzhalter():
        return None
except Exception:
    pass

def scan(file):
    if file is None:
        raise gr.Error("Keine Datei")
    src = file if isinstance(file, str) else file.name
    try:
        xml = scan_file(src, timeout=300, max_pages=4)
    except Exception as e:
        raise gr.Error("Erkennung fehlgeschlagen: " + str(e)[-400:])
    out = os.path.join(tempfile.mkdtemp(), "noten.musicxml")
    with open(out, "w", encoding="utf-8") as f:
        f.write(xml)
    return out

def scan_gs(file, shift, layout=""):
    """Griffschrift (Steirische, 3- und 4-Reiher; Stil Franz Fuchs, 4. Reihe auch als Dreieck wie bei Elke Margetich)
    -> JSON mit MusicXML und Griffen (Pins) für den Wandler. layout: Diskantbelegung des Wandlers als JSON (leer = 446 Normal)."""
    if file is None:
        raise gr.Error("Keine Datei")
    src = file if isinstance(file, str) else file.name
    try:
        import gs2xml
        pages = [src]
        if src.lower().endswith(".pdf"):
            import fitz
            doc = fitz.open(src); tmp = tempfile.mkdtemp(); pages = []
            for i, pg in enumerate(doc):
                if i >= 4: break
                fn = os.path.join(tmp, f"s{i}.png"); pg.get_pixmap(dpi=200).save(fn); pages.append(fn)
            try:   # Vektor-PDF (MC-Score, Finale): Köpfe aus der Textebene neben die Seitenbilder legen (gsleser/vektor.py)
                import vektor
                vektor.save_sidecars(src, pages, 200)
            except Exception:
                pass
        bilder = []   # Kontrollbilder je Seite (JPEG): erkannte Köpfe, Notenwerte, nicht aufgehende Takte
        # Schreibweise: steckt im Belegungs-JSON als „style“ (so bleiben ältere Wandler und die Schnittstelle mit drei Werten gültig); fehlt sie: automatisch
        style = "auto"; erkenner = ""
        try:
            lj = json.loads(layout) if layout else {}
            style = lj.get("style") or "auto"
            erkenner = lj.get("erkenner") or ""   # Kopf-Erkenner je Scan: "alt" = bisheriger, "yolo" = gelernter; leer = Space-Variable GS_YOLO
        except Exception:
            pass
        gs2xml.ERKENNER = erkenner if erkenner in ("alt", "yolo") else ""
        text = {}   # Texterkennung (Tesseract): Titel, Bearbeiter, Teilbezeichnungen, gewählte Schreibweise
        seiten = []   # Daten zum Nachbearbeiten im Wandler (entzerrtes Bild, Systeme, Ereignisse je Takt)
        xml, pins = gs2xml.convert(pages, shift=int(shift or 0), table=gs2xml.parse_table(layout), bilder=bilder, style=style, text=text, seiten_out=seiten)
        try:
            text["erkenner"] = gs2xml.D.LAST.get("erkenner", "")   # welcher Kopf-Erkenner wirklich lief (ohne Modell: bisheriger)
        except Exception:
            pass
    except Exception as e:
        raise gr.Error("Erkennung fehlgeschlagen: " + str(e)[-400:])
    import base64
    tmp = tempfile.mkdtemp(); out = os.path.join(tmp, "griffschrift.json"); imgs = []
    for i, b in enumerate(bilder):
        fn = os.path.join(tmp, f"kontrollbild_{i + 1}.jpg")
        with open(fn, "wb") as f:
            f.write(b)
        imgs.append(fn)
    with open(out, "w", encoding="utf-8") as f:
        json.dump({"xml": xml, "pins": pins, "text": text, "bild": ["data:image/jpeg;base64," + base64.b64encode(b).decode() for b in bilder], "seiten": seiten}, f)
    return out, imgs

nn = gr.Interface(
    fn=scan,
    inputs=gr.File(label="Foto oder PDF gedruckter Noten", file_types=["image", ".pdf"]),
    outputs=gr.File(label="MusicXML"),
    title="Notenerkennung für den Griffschrift-Wandler",
    description="Foto oder PDF gedruckter Noten → MusicXML (homr). Die Datei wird nur für die Erkennung verwendet und nach spätestens einer Stunde gelöscht.",
    api_name="scan",
    flagging_mode="never",
    delete_cache=(1800, 3600),   # hochgeladene Dateien und Ergebnisse nach spätestens einer Stunde löschen
)

gs = gr.Interface(
    fn=scan_gs,
    inputs=[gr.File(label="Foto oder PDF gedruckter Griffschrift", file_types=["image", ".pdf"]),
            gr.Number(label="Stimmung als Verschiebung (A-D-G-C = 0, G-C-F-B = -2, B-Es-As-Des = 1 …)", value=0, precision=0),
            gr.Textbox(label="Diskantbelegung als JSON (leer = 446 Normal; der Wandler schickt seine eingestellte Belegung mit; dazu optional \"style\": auto, fuchs, michlbauer, rosenzopf, margetich, bayerisch oder steinkellner)", value="")],
    outputs=[gr.File(label="JSON mit MusicXML und Griffen"), gr.Gallery(label="Kontrollbild", columns=1, height="auto")],
    title="Griffschrift lesen (Versuch)",
    description="Gedruckte Steirische Griffschrift (3- und 4-Reiher; Stil Franz Fuchs, 4. Reihe auch als Dreieck wie bei Elke Margetich) → Normalnoten (MusicXML) und Griffe mit Balgrichtung für den Griffschrift-Wandler. Versuchsstand. Das Kontrollbild zeigt, was gelesen wurde: Kreise um die erkannten Köpfe (Farbe = Reihe: 1 blau, 2 grün, 3 rot, 4 lila), über dem System der Notenwert (8 Achtel, 4 Viertel, 4. punktierte Viertel, P Pause), unter dem System die Taktnummer; rot hinterlegte Takte mit „!“ gehen nicht auf.",
    api_name="scan_gs",
    flagging_mode="never",
    delete_cache=(1800, 3600),
)

demo = gr.TabbedInterface([nn, gs], ["Normalnoten", "Griffschrift (Versuch)"])

if __name__ == "__main__":
    demo.queue(max_size=8).launch()
