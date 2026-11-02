# Notenerkennung (OMR) für den Griffschrift-Wandler auf einem eigenen Server: Foto oder PDF -> MusicXML mit homr
# (https://github.com/liebharc/homr, AGPL-3.0). Kleiner Webdienst mit FastAPI. Weil eine Seite ein bis zwei Minuten
# brauchen kann und Cloudflare Antworten nach 100 Sekunden abbricht, läuft die Erkennung als Auftrag:
#   GET  /health          -> {"ok": true, "queue": n}
#   POST /scan            -> Formularfeld "file" (Bild oder PDF); Antwort sofort: {"id": "..."}
#                            mit kind=gs und shift=n: Griffschrift lesen (Ergebnis JSON mit MusicXML und Griffen)
#   GET  /job/{id}        -> {"state": "queued"|"running"|"done"|"error", "pos": n, "error": "..."}
#   GET  /job/{id}/xml    -> MusicXML (danach wird der Auftrag gelöscht)
# Hochgeladene Dateien liegen nur bis zum Ende der Erkennung in einem Temp-Ordner; Ergebnisse höchstens 15 Minuten.
import asyncio, glob, os, shutil, subprocess, sys, tempfile, time, uuid
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from omr_core import scan_file   # liegt im selben Ordner (beim Einrichten mitkopieren)
from fastapi import FastAPI, File, Form, HTTPException, UploadFile
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import Response

ALLOWED = os.environ.get("OMR_ORIGINS", "https://jpascher.github.io").split(",")
MAX_MB = int(os.environ.get("OMR_MAX_MB", "20"))
TIMEOUT = int(os.environ.get("OMR_TIMEOUT", "600"))
MAX_QUEUE = int(os.environ.get("OMR_MAX_QUEUE", "10"))
KEEP = 15 * 60

app = FastAPI(title="Notenerkennung Griffschrift-Wandler")
app.add_middleware(CORSMiddleware, allow_origins=ALLOWED, allow_methods=["GET", "POST"], allow_headers=["*"])
JOBS = {}                  # id -> {"state", "path", "work", "xml", "error", "t"}
ORDER = []                 # Warteschlange (ids)
lock = asyncio.Lock()      # immer nur eine Erkennung gleichzeitig (Arbeitsspeicher)

def run_homr(path):
    return scan_file(path, timeout=TIMEOUT)

def run_gs(path, shift):
    # Griffschrift-Leser (Ordner gsleser neben server.py bzw. eine Ebene höher im Repository)
    import json
    here = os.path.dirname(os.path.abspath(__file__))
    for d in (os.path.join(here, "gsleser"), os.path.join(here, "..", "gsleser")):
        if os.path.isdir(d) and d not in sys.path:
            sys.path.insert(0, d)
    import gs2xml
    pages = [path]
    if path.lower().endswith(".pdf"):
        import fitz
        doc = fitz.open(path); pages = []
        for i, pg in enumerate(doc):
            if i >= 4: break
            fn = path + f".s{i}.png"; pg.get_pixmap(dpi=200).save(fn); pages.append(fn)
    xml, pins = gs2xml.convert(pages, shift=shift)
    return json.dumps({"xml": xml, "pins": pins})

def cleanup():
    now = time.time()
    for k in [k for k, j in JOBS.items() if j["state"] in ("done", "error") and now - j["t"] > KEEP]:
        JOBS.pop(k, None)

async def work(jid):
    j = JOBS[jid]
    async with lock:
        j["state"] = "running"
        try:
            if j.get("kind") == "gs":
                j["xml"] = await asyncio.to_thread(run_gs, j["path"], j.get("shift", 0))
            else:
                j["xml"] = await asyncio.to_thread(run_homr, j["path"])
            j["state"] = "done"
        except subprocess.TimeoutExpired:
            j["state"], j["error"] = "error", "Erkennung dauerte zu lange"
        except Exception as e:
            j["state"], j["error"] = "error", "Erkennung fehlgeschlagen: " + str(e)
        finally:
            shutil.rmtree(j["work"], ignore_errors=True)
            j["t"] = time.time()
            if jid in ORDER:
                ORDER.remove(jid)

@app.get("/health")
def health():
    return {"ok": True, "queue": len(ORDER)}

@app.post("/scan")
async def scan(file: UploadFile = File(...), kind: str = Form("nn"), shift: int = Form(0)):
    cleanup()
    if len(ORDER) >= MAX_QUEUE:
        raise HTTPException(503, "Gerade zu viele Aufträge, bitte später noch einmal")
    ext = os.path.splitext(file.filename or "")[1].lower()
    if ext not in (".png", ".jpg", ".jpeg", ".webp", ".tif", ".tiff", ".bmp", ".pdf"):
        raise HTTPException(400, "Nur Bilder (PNG, JPG …) oder PDF")
    data = await file.read()
    if len(data) > MAX_MB * 1024 * 1024:
        raise HTTPException(413, f"Datei größer als {MAX_MB} MB")
    work_dir = tempfile.mkdtemp(prefix="omr-")
    path = os.path.join(work_dir, "noten" + ext)
    with open(path, "wb") as f:
        f.write(data)
    jid = uuid.uuid4().hex
    JOBS[jid] = {"state": "queued", "path": path, "work": work_dir, "xml": None, "error": "", "t": time.time(), "kind": kind, "shift": shift}
    ORDER.append(jid)
    asyncio.create_task(work(jid))
    return {"id": jid}

@app.get("/job/{jid}")
def job(jid: str):
    j = JOBS.get(jid)
    if not j:
        raise HTTPException(404, "Auftrag unbekannt oder abgelaufen")
    pos = ORDER.index(jid) if jid in ORDER else 0
    return {"state": j["state"], "pos": pos, "error": j["error"]}

@app.get("/job/{jid}/xml")
def job_xml(jid: str):
    j = JOBS.get(jid)
    if not j or j["state"] != "done":
        raise HTTPException(404, "Kein Ergebnis")
    JOBS.pop(jid, None)
    mt = "application/json" if j.get("kind") == "gs" else "application/vnd.recordare.musicxml+xml"
    return Response(content=j["xml"], media_type=mt)
