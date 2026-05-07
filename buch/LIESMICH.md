# Akkordeon-Buch · Vollständige Projektmappe

Komplette Materialsammlung — alle Dateien jetzt mit ≤ 41 MB.

---

## 📕 Druckfertige PDFs (KDP-Format 6×9 Zoll)

### Deutsch
- `Akkordeon_Einband.pdf` — Einbändige Komplettausgabe (~635 Seiten)
- `Akkordeon_Band_1.pdf` — Aufbau, Bauteile, Hersteller (260 Seiten)
- `Akkordeon_Band_2.pdf` — Akkordeon-Arten und verwandte Instrumente (375 Seiten)

### Englisch
- `Accordion_Volume_1.pdf` — Construction, Components, Builders (249 Seiten)
- `Accordion_Volume_2.pdf` — Accordion Types and Related Instruments (339 Seiten)

---

## 📦 LaTeX-Quellen (jeweils zweiteilig zum Download)

| Buch | Teil 1 (Quellen + Bilder A–K) | Teil 2 (Bilder L–Z) |
|---|---|---|
| **DE Band 1** | `LaTeX_DE_Band1_Teil1.zip` (41 MB) | `LaTeX_DE_Band1_Teil2.zip` (30 MB) |
| **DE Band 2** | `LaTeX_DE_Band2_Teil1.zip` (41 MB) | `LaTeX_DE_Band2_Teil2.zip` (30 MB) |
| **EN Volume 1** | `LaTeX_EN_Volume1_Teil1.zip` (22 MB) | `LaTeX_EN_Volume1_Teil2.zip` (13 MB) |
| **EN Volume 2** | `LaTeX_EN_Volume2_Teil1.zip` (20 MB) | `LaTeX_EN_Volume2_Teil2.zip` (18 MB) |
| **Build-Skripte** | `Build_Scripts.zip` (46 KB) | — |

**Inhalt von Teil 1 (jedes Buchs):**
- `buch.tex` — Master-LaTeX-Datei mit Präambel
- `kapitel/` — alle 45–46 Kapitel-Dateien
- `bilder/` — erste Hälfte der Bilder (alphabetisch sortiert)

**Inhalt von Teil 2:**
- `bilder/` — zweite Hälfte der Bilder

**Wieder zusammenführen:** Beide ZIPs in dasselbe Verzeichnis entpacken — dann verschmelzen die `bilder/`-Ordner automatisch.

```bash
unzip LaTeX_EN_Volume1_Teil1.zip
unzip LaTeX_EN_Volume1_Teil2.zip   # in dasselbe Verzeichnis!
lualatex buch.tex                  # zweimal aufrufen für korrektes Inhaltsverzeichnis
lualatex buch.tex
```

**Inhalt von Build_Scripts.zip:**
- `verarbeite_quellen.py`, `verarbeite_quellen_en.py` — Wikitext → LaTeX
- `baue_komplett_en.py`, `baue_baende_en.py` — Build-Wrapper
- `einleitung_en.tex`, `einleitung_en_band1.tex`, `einleitung_en_band2.tex`

---

## 📰 Marketing-Texte

- `DE_KDP_Texte.md` — Deutsche KDP-Beschreibungen, Spine, Cover
- `EN_KDP_Texte_Paperback.md` — KDP-Texte für beide englischen Paperback-Bände
- `EN_Back_Cover_Texte.md` — Druckfertige Back-Cover-Texte
- `EN_Kindle_Beide_Baende.md` — Kindle-eBook-Beschreibungen
- `EN_Kindle_Volume_2_komplett.md` — Komplettes Kindle-Listing für Volume II
- `Cover_Prompts_Gemini.md` — Gemini-Bildprompts

---

## 📚 Wiki-Quellen

- `Wiki_DE_Originalartikel.zip` — 91 deutsche Wikipedia-Artikel (474 KB)
- `Wiki_EN_Uebersetzungen.zip` — 91 englische Übersetzungen (333 KB)
- `bilder_index.json` — slug → Bilder-Mapping

---

## ⚖ Lizenz

Alle 91 Artikel sind **CC-BY-SA 4.0** (Wikipedia-Lizenz).
**Bei KDP:** **KDP Select NICHT aktivieren** — sonst Verstoß.

## 🛠 Compiler

- **LuaLaTeX** (TeX Live 2023+)
- Schrift: **Inter** + **DejaVu Sans Mono**
- Build-Skripte: **Python 3.10+**
