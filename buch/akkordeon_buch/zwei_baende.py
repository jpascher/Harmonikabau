#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
zwei_baende.py
==============
Baut das Buch "Das Akkordeon und seine Geschichte" als zwei
KDP-Paperback-Bände (jeweils unter 400 Seiten).

  Band 1: Teile I-VI  (Akkordeon-Kern: Konstruktion, Personen, Bauformen)
  Band 2: Teile VII-XII (Erweiterungen: Automaten, Förderer, Glossar)

Aufruf:
  python3 zwei_baende.py [--quellen PFAD] [--ausgabe PFAD]

Erzeugt:
  buch_band1/buch.pdf
  buch_band2/buch.pdf
"""

import argparse
import json
import shutil
import subprocess
import sys
from pathlib import Path

# Eigene Module — wir importieren das Hauptskript direkt
sys.path.insert(0, str(Path(__file__).resolve().parent))
import verarbeite_quellen as vq


BAND_AUFTEILUNG = [
    # (Bandnummer, Bandtitel, Untertitel-Spezifikum, Teile-Indices)
    (
        1,
        "Das Akkordeon und seine Geschichte",
        "Band I --- Aufbau, Bauteile, Hersteller",
        [0, 1, 2, 3, 4, 5, 6],   # Teile I-VII
    ),
    (
        2,
        "Das Akkordeon und seine Geschichte",
        "Band II --- Akkordeon-Arten und verwandte Instrumente",
        [7, 8, 9, 10, 11, 12], # Teile VIII-XIII
    ),
]


def baue_band(band_nr: int, titel: str, untertitel: str,
              teile_indices: list, quellen: Path, basis_ausgabe: Path,
              kompilieren: bool) -> bool:
    print()
    print("=" * 64)
    print(f" BAND {band_nr}: {untertitel}")
    print("=" * 64)

    ausgabe = basis_ausgabe.parent / f"buch_band{band_nr}"
    kapitel_dir = ausgabe / "kapitel"
    bilder_dir = ausgabe / "bilder"
    kapitel_dir.mkdir(parents=True, exist_ok=True)
    bilder_dir.mkdir(parents=True, exist_ok=True)

    # Bilder kopieren (alle, ist einfacher und nimmt nicht viel Platz)
    quell_bilder = quellen / "bilder"
    if quell_bilder.exists():
        for src in quell_bilder.iterdir():
            if src.is_file():
                ziel = bilder_dir / src.name
                if not ziel.exists() or ziel.stat().st_size != src.stat().st_size:
                    shutil.copy2(src, ziel)

    # Bilder-Index laden
    idx_pfad = quellen / "bilder_index.json"
    bilder_index = json.loads(idx_pfad.read_text(encoding="utf-8"))

    # Nur die ausgewählten Teile bauen
    kapitel_pfade = []
    for ti in teile_indices:
        teil_titel, artikel = vq.KAPITEL[ti]
        print(f"\n=== {teil_titel} ===")
        pfade = []
        for art in artikel:
            p = vq.baue_kapitel(art, quellen, kapitel_dir, bilder_index)
            if p:
                pfade.append(p)
        kapitel_pfade.append((teil_titel, pfade))

    # Master-Datei mit angepasstem Titel + Untertitel
    schreibe_band_tex(ausgabe, titel, untertitel, kapitel_pfade, band_nr)

    if kompilieren:
        print(f"\n>>> Band {band_nr} kompilieren...")
        for lauf in (1, 2):
            print(f"  LuaLaTeX-Lauf {lauf} von 2...")
            subprocess.run(
                ["lualatex", "-interaction=nonstopmode", "buch.tex"],
                cwd=ausgabe,
                stdout=subprocess.DEVNULL,
            )
        pdf = ausgabe / "buch.pdf"
        if pdf.exists():
            seiten = ermittle_seitenzahl(pdf)
            mb = pdf.stat().st_size / (1024 * 1024)
            print(f"\n  >>> Band {band_nr} fertig: {pdf} ({seiten} Seiten, {mb:.1f} MB)")
            return True
        return False

    return True


def schreibe_band_tex(ausgabe: Path, titel: str, untertitel: str,
                       kapitel_pfade: list, band_nr: int) -> None:
    """Schreibt die Master-LaTeX-Datei mit Band-spezifischer Titelseite."""
    titelseite = TITELSEITE_BAND.replace("__BAND__", untertitel)

    lines = [
        vq.PREAMBLE,
        "",
        f"\\title{{{titel}}}",
        f"\\author{{{vq.BUCH_AUTOR}}}",
        "",
        "\\begin{document}",
        "",
        titelseite,
        "",
        vq.LIZENZSEITE,
        "",
        "\\tableofcontents",
        "\\cleardoublepage",
        "",
    ]
    for teil_titel, pfade in kapitel_pfade:
        if not pfade:
            continue
        lines.append(f"\\part{{{vq.escape_latex(teil_titel)}}}")
        for p in pfade:
            rel = p.relative_to(ausgabe).as_posix()
            lines.append(f"\\input{{{rel}}}")
        lines.append("")
    lines.append("\\end{document}")
    (ausgabe / "buch.tex").write_text("\n".join(lines), encoding="utf-8")


def ermittle_seitenzahl(pdf: Path) -> int:
    """Liest die Seitenzahl aus dem PDF (via pdfinfo)."""
    try:
        out = subprocess.check_output(["pdfinfo", str(pdf)], text=True)
        for line in out.splitlines():
            if line.startswith("Pages:"):
                return int(line.split(":")[1].strip())
    except Exception:
        pass
    return -1


# Modifizierte Titelseite mit Band-Bezeichnung
TITELSEITE_BAND = r"""\thispagestyle{empty}
\begin{titlepage}
\centering
\vspace*{2.0cm}
{\sffamily\itshape\color{t0gray}\Large
Eine Sammlung von Wikipedia-Artikeln\par}
\vspace{1.4cm}
{\sffamily\bfseries\color{t0blue}\fontsize{28}{34}\selectfont
Das Akkordeon\\[0.3em]
und seine Geschichte\par}
\vspace{1.0cm}
{\sffamily\large\color{t0blue!70!black}
__BAND__\par}
\vspace{1.6cm}
{\centering\rule{6cm}{0.4pt}\par}
\vspace{1.0cm}
{\sffamily\large
zusammengestellt von\par}
\vspace{0.5em}
{\sffamily\bfseries\Large
Johann Pascher\par}
\vfill
{\sffamily\small\color{t0gray}
Inhaltliche Quelle: deutschsprachige Wikipedia\\[0.3em]
Lizenz: Creative Commons Attribution--ShareAlike 4.0 International\\[0.3em]
\textsc{cc-by-sa 4.0}\par}
\vspace{0.6cm}
{\sffamily\small\color{t0gray}
Zusammenstellung: \today\par}
\end{titlepage}

\clearpage
\thispagestyle{empty}
~\vfill
{\footnotesize\color{t0gray}
Diese Buchausgabe ist eine private Zusammenstellung von Artikeln aus der
deutschsprachigen Wikipedia. Die Artikel stehen unter der Creative-Commons-Lizenz
\textsc{cc-by-sa 4.0}; ihre Autorenschaft ist die Gemeinschaft der jeweiligen
Wikipedia-Bearbeiter. Die vollständige Versionsgeschichte ist über die
Quellenlinks am Ende jedes Kapitels zugänglich.\par}
\clearpage
"""


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--quellen", type=Path,
                        default=Path("akkordeon_quellen"),
                        help="Pfad zum Quellenverzeichnis")
    parser.add_argument("--ausgabe", type=Path,
                        default=Path("buch_ausgabe"),
                        help="Basis-Ausgabeverzeichnis (Bände kommen daneben)")
    parser.add_argument("--nicht-kompilieren", action="store_true",
                        help="Nur LaTeX erzeugen, nicht kompilieren")
    args = parser.parse_args()

    quellen = args.quellen.resolve()
    if not quellen.exists():
        print(f"FEHLER: {quellen} nicht gefunden", file=sys.stderr)
        return 1

    ausgabe = args.ausgabe.resolve()

    # Beide Bände bauen
    for band_nr, titel, untertitel, teile in BAND_AUFTEILUNG:
        baue_band(band_nr, titel, untertitel, teile,
                  quellen, ausgabe, not args.nicht_kompilieren)

    print()
    print("=" * 64)
    print(" Beide Bände fertig.")
    print("=" * 64)
    return 0


if __name__ == "__main__":
    sys.exit(main())
