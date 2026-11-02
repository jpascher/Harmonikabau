# Notenerkennung (OMR) für den Griffschrift-Wandler als Hugging Face Space (Gradio): Foto oder PDF -> MusicXML mit homr
# (https://github.com/liebharc/homr, AGPL-3.0). Der Wandler ruft die Funktion „scan“ über die Gradio-Schnittstelle auf
# (/gradio_api/upload, /gradio_api/call/scan). homr rechnet mit onnxruntime auf der CPU, eine Grafikkarte braucht es nicht.
import os, tempfile
import gradio as gr
from omr_core import scan_file

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

demo = gr.Interface(
    fn=scan,
    inputs=gr.File(label="Foto oder PDF gedruckter Noten", file_types=["image", ".pdf"]),
    outputs=gr.File(label="MusicXML"),
    title="Notenerkennung für den Griffschrift-Wandler",
    description="Foto oder PDF gedruckter Noten → MusicXML (homr). Die Datei wird nur für die Erkennung verwendet und nach spätestens einer Stunde gelöscht.",
    api_name="scan",
    flagging_mode="never",
    delete_cache=(1800, 3600),   # hochgeladene Dateien und Ergebnisse nach spätestens einer Stunde löschen
)

if __name__ == "__main__":
    demo.queue(max_size=8).launch()
