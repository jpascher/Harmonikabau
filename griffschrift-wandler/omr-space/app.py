# Notenerkennung (OMR) für den Griffschrift-Wandler: Foto oder PDF gedruckter Noten -> MusicXML, mit homr
# (https://github.com/liebharc/homr, AGPL-3.0). Läuft als Gradio-Space bei Hugging Face; der Wandler ruft die
# Funktion „scan“ über die Gradio-Schnittstelle auf (/gradio_api/upload, /gradio_api/call/scan).
import os, shutil, subprocess, tempfile, glob
import gradio as gr

try:                      # ZeroGPU (kostenlos für ältere, bestätigte Konten): Funktion bekommt für die Dauer des Aufrufs eine GPU
    import spaces
    gpu = spaces.GPU(duration=180)
except Exception:         # CPU-Space (PRO) oder lokal: ohne Dekorator
    gpu = lambda f: f

@gpu
def scan(file):
    if file is None:
        raise gr.Error("Keine Datei")
    src = file if isinstance(file, str) else file.name
    work = tempfile.mkdtemp()
    name = os.path.basename(src)
    dst = os.path.join(work, name)
    shutil.copy(src, dst)
    # homr schreibt die MusicXML-Datei neben das Bild (gleicher Name, Endung .musicxml)
    r = subprocess.run(["homr", dst], capture_output=True, text=True, timeout=900)
    out = sorted(glob.glob(os.path.join(work, "*.musicxml")))
    if r.returncode != 0 or not out:
        raise gr.Error("Erkennung fehlgeschlagen: " + (r.stderr or r.stdout)[-800:])
    return out[0]

demo = gr.Interface(
    fn=scan,
    inputs=gr.File(label="Foto oder PDF gedruckter Noten", file_types=["image", ".pdf"]),
    outputs=gr.File(label="MusicXML"),
    title="Notenerkennung für den Griffschrift-Wandler",
    description="Lädt ein Foto oder eine PDF gedruckter Noten und liefert MusicXML (homr). "
                "Die Datei wird auf diesem Server verarbeitet und danach nicht aufbewahrt.",
    api_name="scan",
    flagging_mode="never",
)

if __name__ == "__main__":
    demo.queue(max_size=8).launch()
