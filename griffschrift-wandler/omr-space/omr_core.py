# Gemeinsamer Kern für beide Fassungen (Hugging Face: app.py, eigener Server: vps/server.py):
# Bild oder PDF -> MusicXML mit homr (https://github.com/liebharc/homr, AGPL-3.0). homr liest nur Bilder;
# PDF-Seiten werden mit PyMuPDF in Bilder (300 dpi) umgewandelt, einzeln erkannt und die Takte zu einem Stück zusammengefügt.
import glob, os, shutil, subprocess, sys, tempfile
import xml.etree.ElementTree as ET

IMG = (".png", ".jpg", ".jpeg", ".webp", ".tif", ".tiff", ".bmp")

def _homr(img, timeout):
    exe = shutil.which("homr") or os.path.join(os.path.dirname(sys.executable), "homr")
    r = subprocess.run([exe, img], capture_output=True, text=True, timeout=timeout)
    out = os.path.splitext(img)[0] + ".musicxml"
    if r.returncode != 0 or not os.path.exists(out):
        lines = [l for l in (r.stderr + "\n" + r.stdout).splitlines() if l.strip() and not l.startswith(("Downloaded", "Processing"))]
        raise RuntimeError(" / ".join(lines[-2:])[-400:] or "homr ohne Ergebnis")
    return out

def _merge(files):
    """Takte der 1. Stimme aller Seiten hintereinander in die erste Datei; Taktnummern fortlaufend."""
    tree = ET.parse(files[0]); root = tree.getroot()
    part = root.find("part")
    n = len(part.findall("measure"))
    for f in files[1:]:
        p2 = ET.parse(f).getroot().find("part")
        if p2 is None:
            continue
        for m in p2.findall("measure"):
            n += 1; m.set("number", str(n)); part.append(m)
    return ET.tostring(root, encoding="unicode", xml_declaration=True)

def scan_file(path, timeout=600, max_pages=8):
    ext = os.path.splitext(path)[1].lower()
    work = tempfile.mkdtemp(prefix="omrc-")
    try:
        imgs = []
        if ext == ".pdf":
            import fitz  # PyMuPDF
            doc = fitz.open(path)
            if doc.page_count > max_pages:
                raise RuntimeError(f"PDF hat {doc.page_count} Seiten, erlaubt sind höchstens {max_pages}")
            for i, page in enumerate(doc):
                fn = os.path.join(work, f"seite{i+1:02d}.png")
                page.get_pixmap(dpi=300).save(fn); imgs.append(fn)
        elif ext in IMG:
            fn = os.path.join(work, "seite01" + ext); shutil.copy(path, fn); imgs.append(fn)
        else:
            raise RuntimeError("Nur Bilder (PNG, JPG …) oder PDF")
        outs = [_homr(i, timeout) for i in imgs]
        if len(outs) == 1:
            with open(outs[0], encoding="utf-8") as f:
                return f.read()
        return _merge(outs)
    finally:
        shutil.rmtree(work, ignore_errors=True)
