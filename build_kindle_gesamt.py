#!/usr/bin/env python3
"""
Harmonikabau — Gesamt-PDF im Kindle-Buchformat (6 × 9 Zoll)
Zweiseitiges Layout mit KDP-konformen Rändern.

KDP-Anforderungen bei ~200 Seiten:
  - Bundsteg (innen):  min. 12.700 mm (0.5")  → wir verwenden 16 mm
  - Außenrand:         min.  6.350 mm (0.25")  → wir verwenden 10 mm
  - Oben/Unten:        min.  6.350 mm (0.25")  → wir verwenden 10 mm

Ungerade Seiten (rechts): Bundsteg = links
Gerade Seiten (links):    Bundsteg = rechts
"""
import os, shutil
import fitz  # PyMuPDF — für Inhalts-Begrenzungsrahmen
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                 Table, TableStyle, PageBreak)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping
# pypdf nicht mehr benötigt (fitz übernimmt PDF-Assemblierung)

# ══════════════════════════════════════════════════════════
# Konstanten
# ══════════════════════════════════════════════════════════
W_PT = 6 * 72    # 432 pt
H_PT = 9 * 72    # 648 pt

GUTTER   = 16 * mm   # Bundsteg / innerer Rand (KDP min 12.7mm)
OUTSIDE  = 10 * mm   # Äußerer Rand (KDP min 6.35mm)
TOP      = 10 * mm   # Oben (KDP min 6.35mm)
BOTTOM   = 10 * mm   # Unten (KDP min 6.35mm)

# Für ReportLab: symmetrisch mit Gutter-Rand auf beiden Seiten (sicher)
RL_MARGIN = GUTTER  # 16mm auf beiden Seiten — sicher für jede Position
PW = W_PT - 2 * RL_MARGIN

# ══════════════════════════════════════════════════════════
# Fonts
# ══════════════════════════════════════════════════════════
pdfmetrics.registerFont(TTFont('DejaVu',  '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuB', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuI', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBI','/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
addMapping('DejaVu', 0, 0, 'DejaVu');  addMapping('DejaVu', 1, 0, 'DejaVuB')
addMapping('DejaVu', 0, 1, 'DejaVuI'); addMapping('DejaVu', 1, 1, 'DejaVuBI')

# ══════════════════════════════════════════════════════════
# Farben & Styles
# ══════════════════════════════════════════════════════════
DB  = HexColor('#16213e')
AR  = HexColor('#e94560')
KG  = HexColor('#2e7d32')
LG  = HexColor('#f5f5f5')
KBG = HexColor('#e8f5e9')

styles = getSampleStyleSheet()
for sn in styles.byName:
    s = styles.byName[sn]
    if hasattr(s, 'fontName'):
        if 'Bold' in s.fontName:   s.fontName = 'DejaVuB'
        elif 'Italic' in s.fontName: s.fontName = 'DejaVuI'
        else: s.fontName = 'DejaVu'

sT   = ParagraphStyle('T',  parent=styles['Title'],  fontSize=20, textColor=DB, spaceAfter=4, alignment=TA_CENTER, fontName='DejaVuB')
sST  = ParagraphStyle('ST', parent=styles['Normal'], fontSize=11, textColor=DB, alignment=TA_CENTER, spaceAfter=2, fontName='DejaVu')
sAb  = ParagraphStyle('Ab', parent=styles['Italic'], fontSize=8.5, textColor=HexColor('#555555'), spaceAfter=6, fontName='DejaVuI')
sCh  = ParagraphStyle('Ch', parent=styles['Heading1'],fontSize=12, textColor=DB, spaceBefore=12, spaceAfter=5, fontName='DejaVuB')
sSCh = ParagraphStyle('SCh',parent=styles['Heading2'],fontSize=10.5, textColor=DB, spaceBefore=8, spaceAfter=3, fontName='DejaVuB')
sB   = ParagraphStyle('Bo', parent=styles['Normal'], fontSize=9.5, leading=13, spaceAfter=5, alignment=TA_JUSTIFY, fontName='DejaVu')
sBI  = ParagraphStyle('BI', parent=sB, fontName='DejaVuI')
sK   = ParagraphStyle('KB', parent=sB, fontSize=9.5, backColor=KBG, borderPadding=5, borderColor=KG, borderWidth=1, spaceAfter=6, fontName='DejaVu')
sTH  = ParagraphStyle('TH', parent=sB, fontSize=8, fontName='DejaVuB', alignment=TA_CENTER, leading=10)
sTL  = ParagraphStyle('TDL',parent=sB, fontSize=8, alignment=TA_LEFT, leading=10, fontName='DejaVu')
sTrT = ParagraphStyle('TrT',parent=sT, fontSize=16, spaceBefore=0, spaceAfter=4)
sTrS = ParagraphStyle('TrS',parent=sAb, fontSize=9, alignment=TA_CENTER, spaceAfter=0)

def hr():
    return Table([['']], colWidths=[PW],
                 style=TableStyle([('LINEBELOW',(0,0),(-1,-1),1.5,AR),
                                   ('FONTNAME',(0,0),(-1,-1),'DejaVu')]))

def mk_tbl(hdr, rows, cw=None):
    data = [[Paragraph(h, sTH) for h in hdr]]
    for row in rows:
        data.append([Paragraph(str(c), sTL) for c in row])
    cw = cw or [PW / len(hdr)] * len(hdr)
    t = Table(data, colWidths=cw, repeatRows=1)
    t.setStyle(TableStyle([
        ('BACKGROUND',(0,0),(-1,0),DB), ('TEXTCOLOR',(0,0),(-1,0),white),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,LG]),
        ('GRID',(0,0),(-1,-1),0.4,HexColor('#cccccc')),
        ('VALIGN',(0,0),(-1,-1),'MIDDLE'),
        ('TOPPADDING',(0,0),(-1,-1),2), ('BOTTOMPADDING',(0,0),(-1,-1),2),
        ('FONTNAME',(0,0),(-1,-1),'DejaVu')])); return t


# ══════════════════════════════════════════════════════════
# ReportLab-Seiten
# ══════════════════════════════════════════════════════════
def build_einleitung(outpath):
    doc = SimpleDocTemplate(outpath, pagesize=(W_PT, H_PT),
                            leftMargin=RL_MARGIN, rightMargin=RL_MARGIN,
                            topMargin=TOP, bottomMargin=BOTTOM)
    story = []
    story.append(Spacer(1, 45*mm))
    story.append(Paragraph('Harmonikabau', ParagraphStyle('BT', parent=sT, fontSize=28, spaceAfter=4)))
    story.append(Paragraph('Akustik, Konstruktion und Praxis', ParagraphStyle('BST', parent=sST, fontSize=14, spaceAfter=6)))
    story.append(Spacer(1,2*mm)); story.append(hr()); story.append(Spacer(1,4*mm))
    story.append(Paragraph('Von der Stimmzunge bis zum Gehäuse', ParagraphStyle('Sub', parent=sAb, fontSize=11, alignment=TA_CENTER)))
    story.append(Spacer(1, 35*mm))
    story.append(Paragraph('Johann Pascher', ParagraphStyle('Au', parent=sST, fontSize=12, spaceAfter=3)))
    story.append(Paragraph('Linz, Österreich — 2025', ParagraphStyle('Yr', parent=sAb, fontSize=9, alignment=TA_CENTER)))
    story.append(PageBreak())

    story.append(Spacer(1, 4*mm))
    story.append(Paragraph('Einleitung', sCh))
    story.append(Paragraph(
        'Dieses Buch fasst eine Reihe technischer Dokumente zusammen, die über viele Jahre '
        'hinweg aus der praktischen Arbeit am Harmonikabau entstanden sind. Sie behandeln die '
        'Akustik der Stimmzunge und Kammer, die Strömungsmechanik im Spalt, die Kopplung von '
        'Zungen und Kammern, die Materialwahl und Konstruktion von Stimmplatten, Bälgen und '
        'Gehäusen sowie die musikalischen Themen Stimmung und Tremolo.', sB))
    story.append(Paragraph(
        'Die Dokumente richten sich an Harmonikabauer, Reparateure und alle, die verstehen '
        'wollen, warum ein Instrument klingt, wie es klingt – und was man daran ändern kann.', sB))

    story.append(Paragraph('Zum Autor', sCh))
    story.append(Paragraph(
        'Meine ersten Erfahrungen mit Harmonikas und deren Reparatur machte ich Anfang der '
        '1980er-Jahre. Was ich über den Bau und die Instandsetzung dieser Instrumente weiß, '
        'wurde mir zunächst von einem erfahrenen Harmonikabauer vermittelt – nicht in einer '
        'formalen Lehre oder Prüfung, sondern durch Zusehen, Mitmachen und viele Gespräche '
        'in der Werkstatt.', sB))
    story.append(Paragraph(
        'Von Beruf bin ich gelernter Fernsehmechaniker. Später war ich als praktischer Lehrer '
        'für Elektronik und Nachrichtentechnik tätig. Diese technische Ausbildung hat meinen '
        'Zugang zum Harmonikabau geprägt: Ich denke in Frequenzen, Impedanzen und '
        'Schwingkreisen – Begriffe, die in der Akustik der Stimmzunge unmittelbar anwendbar sind.', sB))
    story.append(Paragraph(
        'Ab dem Jahr 2000 begann ich, Harmonikas zu bauen. Heute, nach über 25 Jahren '
        'Erfahrung im Bau und in der Reparatur von Harmonikainstrumenten, möchte ich die '
        'dabei gewonnenen Erkenntnisse in dieser Dokumentensammlung weitergeben.', sB))

    story.append(Paragraph('Zur Entstehung der Dokumente', sCh))
    story.append(Paragraph(
        'Die einzelnen Dokumente sind unabhängig voneinander entstanden und behandeln jeweils '
        'ein abgegrenztes Thema. Sie wurden nicht von Anfang an als Buch geplant, sondern sind '
        'das Ergebnis konkreter Fragen, die sich in der Werkstatt stellten: Warum spricht eine '
        'Zunge schlecht an? Wie beeinflusst die Kammergeometrie den Klang? Was passiert, wenn '
        'mehrere Zungen in dieselbe Kammer klingen?', sB))
    story.append(Paragraph(
        'Die Antworten habe ich zunächst für mich selbst aufgeschrieben und mit physikalischen '
        'Modellen untermauert. Dabei kamen mir meine Kenntnisse aus der Elektronik und '
        'Nachrichtentechnik zugute – die Analogie zwischen elektrischen Schwingkreisen und '
        'akustischen Resonatoren ist direkt und fruchtbar.', sB))
    story.append(Paragraph(
        'Die Dokumente erheben nicht den Anspruch einer akademischen Abhandlung. Die angegebenen '
        'Zahlenwerte sind häufig Größenordnungsabschätzungen auf Basis vereinfachter Modelle '
        'und empirischer Kalibrierung. Sie sollen ein qualitatives Verständnis der Zusammenhänge '
        'liefern – kein Ersatz für eigene Messungen am konkreten Instrument sein.', sB))

    story.append(Paragraph('Übersicht der Dokumente', sCh))
    story.append(Paragraph('<i>Die Sammlung gliedert sich in folgende Dokumente:</i>', sB))
    tcw = [14*mm, PW - 14*mm]

    story.append(Paragraph('Akustik und Strömung', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0002','Strömungsanalyse Bass-Stimmzunge 50\u202FHz – Version 8'],
        ['0003','Impedanzvergleich: Durchschlagzunge vs. Labialpfeife'],
        ['0004','Frequenzvariation der Stimmzunge durch Kammerkopplung'],
        ['0005','Frequenzverschiebung als Indikator der Ansprache'],
        ['0006','Zeichenerklärung'],
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph('Kammer und Klang', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0007','Diskant-Stimmstock – Kammerfrequenzen D3–C6'],
        ['0008','Klangveränderung durch Kammergeometrie'],
        ['0009','Frequenzkopplung mehrerer Zungen'],
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph('Stimmzunge und Stimmplatte', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0010','Güte der Stimmzunge'],
        ['0011','Kanalgeometrie der Stimmplatte'],
        ['0012','Zungensteifigkeit'],
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph('Kopplung und Obertonmoden', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0015','Akustische Kopplung und Impedanzanpassung'],
        ['0016','Obertonmoden der Basszunge: Profilierung und Inharmonizität'],
        ['0017','Diskant-Stimmzungen: Obertonmoden F3 bis C6'],
        ['0018','Hörbarkeit der Biegemoden: Kammer-Saugkreis, Transiente, Torsion'],
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph('Stimmung und Tremolo', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0019','Stimmung und Differenztöne'],
        ['0020','Tremolo: Schwebungsphysik, Typen und Stimmungspraxis'],
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph('Konstruktion und Material', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0021','Stimmplatten: Qualität, Hersteller und Güteklassen'],
        ['0022','Balg: Querschnitt, Faltenzahl und Instrumentengröße'],
        ['0023','Gehäuse und Mechanik'],
        ['0024','Praxishinweise: Kritische Handgriffe im Harmonikabau'],
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))
    story.append(Paragraph('Spezialthemen', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0500','Leitfaden zur ästhetischen Forensik bei Akkordeon-Gehäusen'],
    ], cw=tcw))
    story.append(Spacer(1, 4*mm))

    story.append(Paragraph('Hinweise zur Lektüre', sCh))
    story.append(Paragraph(
        'Die Dokumente können grundsätzlich in beliebiger Reihenfolge gelesen werden. '
        'Querverweise am Fuß jedes Dokuments zeigen, welche anderen Dokumente verwandte '
        'Themen behandeln. Die Zeichenerklärung in Dokument\u202F0006 fasst die verwendeten '
        'Formelzeichen zusammen und dient als Nachschlagewerk.', sB))
    story.append(Paragraph(
        'Zu vielen Dokumenten existieren Berechnungsskripte in Python, die auf GitHub '
        'verfügbar sind. Diese Skripte erlauben es, die angegebenen Zahlenwerte selbst '
        'nachzurechnen und mit eigenen Messungen zu vergleichen.', sB))
    story.append(Spacer(1, 8*mm))
    story.append(Paragraph('<i>Johann Pascher — Linz, 2025</i>', sBI))
    doc.build(story)
    print(f'  ✓ Einleitung')


def build_schlusswort(outpath):
    doc = SimpleDocTemplate(outpath, pagesize=(W_PT, H_PT),
                            leftMargin=RL_MARGIN, rightMargin=RL_MARGIN,
                            topMargin=TOP, bottomMargin=BOTTOM)
    story = []
    story.append(Spacer(1, 8*mm))
    story.append(Paragraph('Schlusswort', sT))
    story.append(Spacer(1, 2*mm)); story.append(hr()); story.append(Spacer(1, 4*mm))
    story.append(Paragraph(
        'Die Stimmzunge ist eines der ältesten und zugleich am wenigsten verstandenen '
        'Klangerzeugungsprinzipien der Musikinstrumente. Sie ist mechanisch einfach – '
        'ein dünnes Metallblatt, das durch einen Schlitz schwingt – und akustisch komplex: '
        'Die Wechselwirkung zwischen Zunge, Spalt, Kammer und Balg erzeugt ein Verhalten, '
        'das sich einfachen Modellen entzieht.', sB))
    story.append(Paragraph(
        'Diese Dokumentensammlung hat versucht, die physikalischen Grundlagen dieses '
        'Zusammenspiels aufzuarbeiten – nicht mit dem Anspruch auf Vollständigkeit, sondern '
        'mit dem Ziel, die wesentlichen Mechanismen sichtbar zu machen. Von der '
        'Strömungsanalyse im Zungenspalt über die Impedanzkopplung der Kammer bis zur '
        'Hörbarkeit einzelner Oberton-Moden zieht sich ein roter Faden: Die Akustik der '
        'Harmonika lässt sich mit denselben Werkzeugen beschreiben, die in der Elektronik '
        'und Nachrichtentechnik seit Jahrzehnten bewährt sind – Resonanz, Impedanz, '
        'Güte, Kopplung.', sB))
    story.append(Paragraph(
        'Manches in diesen Dokumenten ist gesichert, manches ist Modell und Abschätzung. '
        'Die Grenze zwischen beiden habe ich nach bestem Wissen kenntlich gemacht. Wo '
        'vereinfachte Modelle an ihre Grenzen stoßen, ist das vermerkt. Die begleitenden '
        'Python-Skripte erlauben es, jede Rechnung selbst nachzuvollziehen und an eigenen '
        'Messwerten zu prüfen.', sB))
    story.append(Paragraph(
        'Der Harmonikabau lebt von der Werkstatt, nicht vom Schreibtisch. Kein Dokument '
        'ersetzt das Hören, Prüfen und Korrigieren am konkreten Instrument. Aber ein '
        'Verständnis der Physik kann helfen, die richtigen Fragen zu stellen – und die '
        'Erfahrung schneller in gute Ergebnisse zu übersetzen.', sK))
    story.append(Paragraph(
        'Ich hoffe, dass diese Sammlung für andere Harmonikabauer, Reparateure und '
        'Interessierte von Nutzen ist. Wer Fehler findet, bessere Modelle kennt oder '
        'eigene Messungen beitragen möchte, ist herzlich eingeladen, mich zu kontaktieren. '
        'Die Dokumente und Berechnungsskripte sind auf GitHub verfügbar und werden bei '
        'Bedarf aktualisiert.', sB))
    story.append(Spacer(1, 3*mm)); story.append(hr()); story.append(Spacer(1, 4*mm))
    story.append(Paragraph('<i>Johann Pascher — Linz, 2025</i>', sBI))
    doc.build(story)
    print(f'  ✓ Schlusswort')


def build_separator(outpath, dok_nr, titel, untertitel=''):
    doc = SimpleDocTemplate(outpath, pagesize=(W_PT, H_PT),
                            leftMargin=RL_MARGIN, rightMargin=RL_MARGIN,
                            topMargin=TOP, bottomMargin=BOTTOM)
    story = []
    story.append(Spacer(1, 55*mm))
    story.append(Paragraph(f'Dok.\u202F{dok_nr}', sST))
    story.append(Spacer(1, 3*mm))
    story.append(Paragraph(titel, sTrT))
    if untertitel:
        story.append(Spacer(1, 2*mm))
        story.append(Paragraph(f'<i>{untertitel}</i>', sTrS))
    story.append(Spacer(1, 4*mm)); story.append(hr())
    doc.build(story)


# ══════════════════════════════════════════════════════════
# Festes Body-Crop-Fenster — Kopf- und Fußzeilen ausblenden
# ══════════════════════════════════════════════════════════
#
# Alle LaTeX-Quelldokumente (A4, 595×842pt) haben:
#   Kopfzeile: fitz y ≈ 36–45pt (sehr nahe Seitenoberrand)
#   Fußzeile:  fitz y ≈ 807–815pt (sehr nahe Seitenunterrand)
#   Body:      fitz x ≈ 55–540pt, y ≈ 65–800pt
#
# Das Crop-Fenster schließt Kopf- und Fußzeile aus und sorgt
# für einen EINHEITLICHEN Skalierungsfaktor auf allen Seiten.

# Body-Grenzen in fitz-Koordinaten (y von oben, A4 = 842pt hoch)
BODY_X0          = 55    # linke Grenze Textbereich
BODY_X1          = 540   # rechte Grenze Textbereich
BODY_WIDTH       = BODY_X1 - BODY_X0   # 485pt → fixer Horizontalscale
HEADER_END_FITZ  = 65    # LaTeX-Kopfzeile endet hier (y von oben)
FOOTER_START_FITZ= 790   # LaTeX-Fußzeile beginnt hier (Dok 0002: y=793, andere: ≥802)

SAFETY = 8   # pt KDP-Sicherheitsabstand zur Randlinie (~3mm)


def _page_body_rect(page) -> fitz.Rect:
    """Tatsächlicher Inhaltsbereich dieser Seite (fitz-Koord.), ohne Kopf/Fußzeile.
    Breite = BODY_WIDTH (fix), Höhe = nur soweit Inhalt vorhanden."""
    rects = []
    zone = fitz.Rect(BODY_X0, HEADER_END_FITZ, BODY_X1, FOOTER_START_FITZ)

    for b in page.get_text('blocks'):
        r = fitz.Rect(b[:4]) & zone
        if not r.is_empty and r.width > 1 and r.height > 1:
            rects.append(r)
    for d in page.get_drawings():
        raw = d.get('rect')
        if raw:
            r = fitz.Rect(raw) & zone
            if not r.is_empty and r.width > 1 and r.height > 1:
                rects.append(r)
    for img in page.get_images(full=True):
        for ir in page.get_image_rects(img[0]):
            r = fitz.Rect(ir) & zone
            if not r.is_empty and r.width > 1 and r.height > 1:
                rects.append(r)

    if not rects:
        # Leere Seite: minimales Rechteck oben
        return fitz.Rect(BODY_X0, HEADER_END_FITZ, BODY_X1, HEADER_END_FITZ + 20)

    PAD = 3   # pt Randpuffer um den Inhalt
    y0 = max(HEADER_END_FITZ, min(r.y0 for r in rects) - PAD)
    y1 = min(FOOTER_START_FITZ, max(r.y1 for r in rects) + PAD)
    return fitz.Rect(BODY_X0, y0, BODY_X1, y1)


# ══════════════════════════════════════════════════════════
# Zusammenbau
# ══════════════════════════════════════════════════════════
def main():
    import subprocess
    basedir = os.path.dirname(os.path.abspath(__file__))
    tmpdir  = os.path.join(basedir, '_tmp_kindle')
    os.makedirs(tmpdir, exist_ok=True)

    einl_path = os.path.join(tmpdir, 'einleitung.pdf')
    schl_path = os.path.join(tmpdir, 'schlusswort.pdf')
    print('Erzeuge Rahmenteile ...')
    build_einleitung(einl_path)
    build_schlusswort(schl_path)

    dokumente = [
        ('bass_50hz_v8.pdf',                          '0002','Strömungsanalyse Bass-Stimmzunge 50\u202FHz','Version 8'),
        ('0003_impedanz_vergleich_De.pdf',             '0003','Impedanzvergleich','Durchschlagzunge vs. Labialpfeife'),
        ('0004_frequenzvariation_zwei_filter_De.pdf',  '0004','Frequenzvariation der Stimmzunge','durch Kammerkopplung'),
        ('0005_ansprache_frequenz_kopplung_De.pdf',    '0005','Frequenzverschiebung','als Indikator der Ansprache'),
        ('0006_zeichenerklaerung_De.pdf',              '0006','Zeichenerklärung',''),
        ('0007_diskant_kammerfrequenzen_De.pdf',       '0007','Diskant-Stimmstock','Kammerfrequenzen D3–C6'),
        ('0008_klangveraenderung_De.pdf',              '0008','Klangveränderung','durch Kammergeometrie'),
        ('0009_frequenzkopplung_mehrere_zungen_De.pdf','0009','Frequenzkopplung','mehrerer Zungen'),
        ('0010_guete_stimmplatte_De.pdf',              '0010','Güte der Stimmzunge',''),
        ('0011_kanalgeometrie_De.pdf',                 '0011','Kanalgeometrie der Stimmplatte',''),
        ('0012_steifigkeit_De.pdf',                    '0012','Zungensteifigkeit',''),
        ('0015_kopplung_De.pdf',                       '0015','Akustische Kopplung','und Impedanzanpassung'),
        ('0016_obertonmoden_De.pdf',                   '0016','Obertonmoden der Basszunge','Profilierung und Inharmonizität'),
        ('0017_diskant_De.pdf',                        '0017','Diskant-Stimmzungen','Obertonmoden F3 bis C6'),
        ('0018_hoerbarkeit_De.pdf',                    '0018','Hörbarkeit der Biegemoden','Kammer-Saugkreis, Transiente, Torsion'),
        ('0019_stimmung_De.pdf',                       '0019','Stimmung und Differenztöne',''),
        ('0020_tremolo_De.pdf',                        '0020','Tremolo','Schwebungsphysik, Typen und Stimmungspraxis'),
        ('0021_stimmplatten_De.pdf',                   '0021','Stimmplatten','Qualität, Hersteller und Güteklassen'),
        ('0022_balg_De.pdf',                           '0022','Balg','Querschnitt, Faltenzahl und Instrumentengröße'),
        ('0023_gehaeuse_De.pdf',                       '0023','Gehäuse und Mechanik',''),
        ('0024_praxishinweise_De.pdf',                 '0024','Praxishinweise','Kritische Handgriffe im Harmonikabau'),
        ('0500_forensik_De.pdf',                       '0500','Ästhetische Forensik','bei Akkordeon-Gehäusen'),
    ]

    # fitz als Haupt-Assembler
    dst = fitz.open()
    page_counter = 0

    def add_full_pages(pdf_path, label=''):
        """ReportLab-Seiten (bereits 6×9) vollflächig einbetten."""
        nonlocal page_counter
        src = fitz.open(pdf_path)
        n   = len(src)
        for i in range(n):
            dp = dst.new_page(width=W_PT, height=H_PT)
            # 1:1 – Seite ist bereits 6×9
            dp.show_pdf_page(fitz.Rect(0, 0, W_PT, H_PT), src, i)
            page_counter += 1
        src.close()
        if label:
            print(f'  + {label} ({n} S., bis S. {page_counter})')

    def add_source_pages(filepath, dok_nr):
        """Quelldokument-Seiten einbetten:
        - FESTER Horizontalscale (einheitliche Schriftgröße)
        - PRO SEITE: Clip bis zum echten Inhaltsende (kein Leerraum durch A4-Whitespace)
        - Header/Footer ausgeblendet
        """
        nonlocal page_counter
        src = fitz.open(filepath)
        n   = len(src)
        for i in range(n):
            page_counter += 1
            is_odd = (page_counter % 2 == 1)
            ml = GUTTER if is_odd else OUTSIDE
            mr = OUTSIDE if is_odd else GUTTER

            # Effektiver Zielbereich
            eff_w = (W_PT - ml - mr) - 2 * SAFETY
            eff_h = (H_PT - TOP - BOTTOM) - 2 * SAFETY

            # Tatsächlicher Inhaltsbereich dieser Seite
            body = _page_body_rect(src[i])

            # FESTER Scale (Breite): gleich für alle Seiten → einheitliche Schrift
            scale = eff_w / BODY_WIDTH        # ≈ 0.706, immer gleich
            new_w = BODY_WIDTH * scale         # füllt Breite komplett
            new_h = min(body.height * scale, eff_h)  # nur echter Inhalt, max eff_h

            # Horizontal zentrieren (geringfügig, da new_w ≈ eff_w),
            # vertikal OBEN ausrichten
            tx = ml + SAFETY + (eff_w - new_w) / 2
            ty = TOP + SAFETY                  # fitz: y von oben

            target = fitz.Rect(tx, ty, tx + new_w, ty + new_h)

            dp = dst.new_page(width=W_PT, height=H_PT)
            dp.show_pdf_page(target, src, i, clip=body)

        src.close()
        print(f'  + Dok. {dok_nr}: {n} S., bis S. {page_counter}')

    print('Zusammenbau ...')
    add_full_pages(einl_path, 'Einleitung')

    for filename, dok_nr, titel, untertitel in dokumente:
        filepath = os.path.join(basedir, filename)
        if not os.path.exists(filepath):
            print(f'  ⚠ {filename} nicht gefunden')
            continue
        sep_path = os.path.join(tmpdir, f'sep_{dok_nr}.pdf')
        build_separator(sep_path, dok_nr, titel, untertitel)
        add_full_pages(sep_path)
        add_source_pages(filepath, dok_nr)

    add_full_pages(schl_path, 'Schlusswort')

    # Gerade Seitenzahl für KDP
    if page_counter % 2 == 1:
        dst.new_page(width=W_PT, height=H_PT)
        page_counter += 1
        print(f'  + Leerseite → {page_counter} S.')

    # Speichern
    outfile  = os.path.join(basedir, 'Harmonikabau_Kindle_6x9.pdf')
    tmp_out  = outfile + '.tmp.pdf'
    dst.save(tmp_out, garbage=4, deflate=True)
    dst.close()

    # Ghostscript: alle Schriften vollständig einbetten
    gs_cmd = [
        'gs', '-dNOPAUSE', '-dBATCH', '-dQUIET',
        '-sDEVICE=pdfwrite',
        '-dCompatibilityLevel=1.4',
        '-dEmbedAllFonts=true',
        '-dSubsetFonts=true',
        '-dPDFSETTINGS=/prepress',
        f'-sOutputFile={outfile}',
        tmp_out
    ]
    result = subprocess.run(gs_cmd, capture_output=True, text=True)
    if result.returncode == 0:
        os.remove(tmp_out)
        print('  Schriften eingebettet (Ghostscript)')
    else:
        os.rename(tmp_out, outfile)
        print('  ⚠ Ghostscript nicht verfügbar')

    print(f'\n✓ {outfile}')
    print(f'  {page_counter} Seiten, 6×9 Zoll, zweiseitig')
    print(f'  Scale: {(W_PT-GUTTER-OUTSIDE-2*SAFETY)/BODY_WIDTH:.4f} (einheitlich, breitenbasiert)')
    print(f'  Body-Crop: {BODY_WIDTH:.0f}pt breit, höhenadaptiv pro Seite')

    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == '__main__':
    main()
