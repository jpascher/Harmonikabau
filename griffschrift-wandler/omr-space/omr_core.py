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

def _ties(xml):
    """homr schreibt Haltebögen als Bindebögen (slur). Ein Bogen von einem Akkord zum unmittelbar folgenden mit genau denselben
    Tönen ist ein Haltebogen: dort werden die Bögen durch tie/tied ersetzt, sonst würde der Ton zweimal angeschlagen."""
    try: root = ET.fromstring(xml.encode("utf-8") if isinstance(xml, str) else xml)
    except ET.ParseError: return xml
    def pit(n):
        p = n.find("pitch")
        if p is None: return None
        al = p.find("alter")
        return (p.findtext("step"), p.findtext("octave"), al.text if al is not None else "0")
    for part in root.iter("part"):
        groups = []
        for m in part.iter("measure"):
            for n in m.findall("note"):
                if n.find("chord") is not None and groups: groups[-1].append(n)
                else: groups.append([n])
        def slurs(g, typ): return [(n, sl) for n in g for nt in n.findall("notations") for sl in nt.findall("slur") if sl.get("type") == typ]
        for g1, g2 in zip(groups, groups[1:]):
            st, sp = slurs(g1, "start"), slurs(g2, "stop")
            if not st or not sp: continue
            p1 = sorted(filter(None, map(pit, g1))); p2 = sorted(filter(None, map(pit, g2)))
            if not p1 or p1 != p2: continue
            for g, typ, lst in ((g1, "start", st), (g2, "stop", sp)):
                for n, sl in lst:
                    for nt in n.findall("notations"):
                        if sl in list(nt): nt.remove(sl)
                for n in g:
                    if n.find("pitch") is None: continue
                    kids = list(n); i = next((k + 1 for k, c in enumerate(kids) if c.tag == "duration"), len(kids))
                    n.insert(i, ET.Element("tie", type=typ))
                    nt = n.find("notations")
                    if nt is None: nt = ET.SubElement(n, "notations")
                    nt.append(ET.Element("tied", type=typ))
    return ET.tostring(root, encoding="unicode", xml_declaration=True)

def _title(img):
    """Titel aus der Seite lesen: die Textzeile mit der größten Schrift im oberen Viertel (homr nimmt manchmal den Untertitel)."""
    try:
        import pytesseract
        from PIL import Image
        im = Image.open(img).convert("L"); w, h = im.size; top = im.crop((0, 0, w, int(h * 0.25)))
        D = pytesseract.image_to_data(top, lang="deu+ita+eng" if "ita" in pytesseract.get_languages() else "deu", config="--psm 3", output_type=pytesseract.Output.DICT)
        lines = {}
        for i, t in enumerate(D["text"]):
            t = (t or "").strip()
            if not t or float(D["conf"][i]) < 60: continue
            k = (D["block_num"][i], D["par_num"][i], D["line_num"][i]); lines.setdefault(k, []).append((D["left"][i], t, D["height"][i]))
        best = None
        for ws in lines.values():
            txt = " ".join(t for _, t, _h in sorted(ws)); hgt = sorted(x[2] for x in ws)[len(ws) // 2]
            if sum(ch.isalpha() for ch in txt) < 3: continue
            if best is None or hgt > best[0]: best = (hgt, txt)
        return best[1] if best else None
    except Exception:
        return None

def _set_title(xml, title):
    if not title: return xml
    import html
    t = html.escape(title.title() if title.isupper() else title)
    if "<work-title>" in xml:
        import re
        return re.sub(r"<work-title>.*?</work-title>", "<work-title>" + t + "</work-title>", xml, count=1, flags=re.S)
    return xml

def _split_double(fn):
    """Querformat mit zwei Notenseiten nebeneinander (Doppelseite, z. B. Stimmenblätter von Blasmusikverlagen): an der hellsten
    senkrechten Spalte im mittleren Drittel in eine linke und eine rechte Seite teilen, weil homr ein System über die ganze Breite
    erwartet und sonst Zeilen beider Seiten vermischt (gefunden 13.11.2026 an einer Flügelhorn-Stimme). Hochformat oder ohne
    leere Mittelspalte: unverändert. OMR_SPLIT=0 schaltet das Teilen ab."""
    if os.environ.get("OMR_SPLIT", "1") != "1": return [fn]
    try:
        import numpy as np
        from PIL import Image
        im = Image.open(fn).convert("L"); w, h = im.size
        if w < 1.25 * h: return [fn]
        a = np.asarray(im, dtype=np.float32); ink = (a < 128).mean(axis=0)   # Anteil dunkler Punkte je Spalte
        k = max(3, w // 200); sm = np.convolve(ink, np.ones(k) / k, mode="same")
        lo, hi = int(w * 0.35), int(w * 0.65); x = lo + int(np.argmin(sm[lo:hi]))
        body = float(np.median(sm[int(w * 0.1):int(w * 0.9)]))
        if sm[x] > 0.15 * body or sm[x] > 0.01: return [fn]   # keine leere Mittelspalte: keine Doppelseite
        if ink[:x].mean() < 0.3 * body or ink[x:].mean() < 0.3 * body: return [fn]   # eine Hälfte fast leer
        base = os.path.splitext(fn)[0]; out = []
        for j, (x0, x1) in enumerate(((0, x), (x, w))):
            f2 = f"{base}_{'ab'[j]}.png"; im.crop((x0, 0, x1, h)).save(f2); out.append(f2)
        return out
    except Exception:
        return [fn]

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
        imgs = [g for i in imgs for g in _split_double(i)]   # Doppelseiten im Querformat in zwei Seiten teilen
        outs = [_homr(i, timeout) for i in imgs]
        if len(outs) == 1:
            with open(outs[0], encoding="utf-8") as f:
                xml = f.read()
        else:
            xml = _merge(outs)
        xml = _ties(xml) if os.environ.get("OMR_TIES", "1") == "1" else xml
        if os.environ.get("OMR_TITLE", "1") == "1": xml = _set_title(xml, _title(imgs[0]))
        return xml
    finally:
        shutil.rmtree(work, ignore_errors=True)
