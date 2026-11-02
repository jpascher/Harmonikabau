# Notenerkennung (OMR) für den Griffschrift-Wandler als Hugging Face Space (Gradio): Foto oder PDF -> MusicXML mit homr
# (https://github.com/liebharc/homr, AGPL-3.0). Der Wandler ruft die Funktion „scan“ über die Gradio-Schnittstelle auf
# (/gradio_api/upload, /gradio_api/call/scan). homr rechnet mit onnxruntime auf der CPU, eine Grafikkarte braucht es nicht.
import os, tempfile
import gradio as gr
from omr_core import scan_file

try:                      # ZeroGPU-Space (kostenlos für Konten älter als 30 Tage): Funktion muss so gekennzeichnet sein
    import spaces
    gpu = spaces.GPU(duration=120)
except Exception:         # CPU-Space oder lokal
    gpu = lambda f: f

@gpu
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
    description="Foto oder PDF gedruckter Noten → MusicXML (homr). Die Datei wird nur für die Erkennung verwendet und nicht aufbewahrt.",
    api_name="scan",
    flagging_mode="never",
)

if __name__ == "__main__":
    demo.queue(max_size=8).launch()
