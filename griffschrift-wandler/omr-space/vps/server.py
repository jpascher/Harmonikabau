# Notenerkennung (OMR) für den Griffschrift-Wandler auf einem eigenen Server: Foto oder PDF -> MusicXML mit homr
# (https://github.com/liebharc/homr, AGPL-3.0). Kleiner Webdienst mit FastAPI. Weil eine Seite ein bis zwei Minuten
# brauchen kann und Cloudflare Antworten nach 100 Sekunden abbricht, läuft die Erkennung als Auftrag:
#   GET  /health          -> {"ok": true, "queue": n}
#   POST /scan            -> Formularfeld "file" (Bild oder PDF); Antwort sofort: {"id": "..."}
#   GET  /job/{id}        -> {"state": "queued"|"running"|"done"|"error", "pos": n, "error": "..."}
#   GET  /job/{id}/xml    -> MusicXML (danach wird der Auftrag gelöscht)
# Hochgeladene Dateien liegen nur bis zum Ende der Erkennung in einem Temp-Ordner; Ergebnisse höchstens 15 Minuten.
import asyncio, glob, os, shutil, subprocess, tempfile, time, uuid
from fastapi import FastAPI, File, HTTPException, UploadFile
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
    r = subprocess.run(["homr", path], capture_output=True, text=True, timeout=TIMEOUT)
    out = sorted(glob.glob(os.path.join(os.path.dirname(path), "*.musicxml")))
    if r.returncode != 0 or not out:
        raise RuntimeError((r.stderr or r.stdout or "homr ohne Ergebnis")[-600:])
    with open(out[0], encoding="utf-8") as f:
        return f.read()

def cleanup():
    now = time.time()
    for k in [k for k, j in JOBS.items() if j["state"] in ("done", "error") and now - j["t"] > KEEP]:
        JOBS.pop(k, None)

async def work(jid):
    j = JOBS[jid]
    async with lock:
        j["state"] = "running"
        try:
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
async def scan(file: UploadFile = File(...)):
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
    JOBS[jid] = {"state": "queued", "path": path, "work": work_dir, "xml": None, "error": "", "t": time.time()}
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
    return Response(content=j["xml"], media_type="application/vnd.recordare.musicxml+xml")
