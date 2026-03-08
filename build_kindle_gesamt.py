#!/usr/bin/env python3
"""
Harmonikabau — Gesamt-PDF im Kindle-Buchformat (6 × 9 Zoll)

Erzeugt ein einzelnes PDF mit:
  - Einleitung (neu gesetzt in 6×9)
  - Alle Einzeldokumente (skaliert auf 6×9)
  - Schlusswort (neu gesetzt in 6×9)

Trennseiten zwischen den Dokumenten für klare Kapitelgrenzen.
"""
import os, sys
from io import BytesIO
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm, inch
from reportlab.lib.colors import HexColor, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer,
                                 Table, TableStyle, PageBreak, KeepTogether)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping
from pypdf import PdfReader, PdfWriter, Transformation, PageObject

# ── Konstanten ──
W_INCH, H_INCH = 6, 9
W_PT, H_PT = W_INCH * 72, H_INCH * 72  # 432 × 648 pt
# Zweiseitiges Layout: Bundsteg (innen) etwas größer als Außenrand
MARGIN_INSIDE  = 11 * mm   # Bundsteg / Gutter (KDP min. 12.7mm bei 200 S. – Inhalt hat eigenen Rand)
MARGIN_OUTSIDE =  6 * mm   # Außenrand (KDP min. 6.35mm)
MARGIN_TB = 10 * mm        # Oben/Unten (KDP min. 6.35mm)
# Für ReportLab-Seiten: symmetrisch setzen, nachher verschieben
MARGIN_LR = (MARGIN_INSIDE + MARGIN_OUTSIDE) / 2
PW = W_PT - MARGIN_INSIDE - MARGIN_OUTSIDE
# Verschiebung pro Seite: 5mm Unterschied → 2.5mm Shift pro Richtung
GUTTER_SHIFT = 2.5 * mm

# ── Fonts ──
pdfmetrics.registerFont(TTFont('DejaVu',  '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuB', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuI', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBI','/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
addMapping('DejaVu', 0, 0, 'DejaVu');  addMapping('DejaVu', 1, 0, 'DejaVuB')
addMapping('DejaVu', 0, 1, 'DejaVuI'); addMapping('DejaVu', 1, 1, 'DejaVuBI')

# ── Farben ──
DB  = HexColor('#16213e')
AR  = HexColor('#e94560')
KG  = HexColor('#2e7d32')
LG  = HexColor('#f5f5f5')
KBG = HexColor('#e8f5e9')

# ── Styles ──
styles = getSampleStyleSheet()
for sn in styles.byName:
    s = styles.byName[sn]
    if hasattr(s, 'fontName'):
        if 'Bold' in s.fontName:   s.fontName = 'DejaVuB'
        elif 'Italic' in s.fontName: s.fontName = 'DejaVuI'
        else: s.fontName = 'DejaVu'

sT   = ParagraphStyle('T',   parent=styles['Title'],   fontSize=20, textColor=DB,
                       spaceAfter=4, alignment=TA_CENTER, fontName='DejaVuB')
sST  = ParagraphStyle('ST',  parent=styles['Normal'],  fontSize=11, textColor=DB,
                       alignment=TA_CENTER, spaceAfter=2, fontName='DejaVu')
sAb  = ParagraphStyle('Ab',  parent=styles['Italic'],  fontSize=8.5, textColor=HexColor('#555555'),
                       spaceAfter=6, fontName='DejaVuI')
sCh  = ParagraphStyle('Ch',  parent=styles['Heading1'],fontSize=12, textColor=DB,
                       spaceBefore=12, spaceAfter=5, fontName='DejaVuB')
sSCh = ParagraphStyle('SCh', parent=styles['Heading2'],fontSize=10.5, textColor=DB,
                       spaceBefore=8, spaceAfter=3, fontName='DejaVuB')
sB   = ParagraphStyle('Bo',  parent=styles['Normal'],  fontSize=9.5, leading=13,
                       spaceAfter=5, alignment=TA_JUSTIFY, fontName='DejaVu')
sBI  = ParagraphStyle('BI',  parent=sB, fontName='DejaVuI')
sK   = ParagraphStyle('KB',  parent=sB, fontSize=9.5, backColor=KBG, borderPadding=5,
                       borderColor=KG, borderWidth=1, spaceAfter=6, fontName='DejaVu')
sTH  = ParagraphStyle('TH',  parent=sB, fontSize=8, fontName='DejaVuB', alignment=TA_CENTER, leading=10)
sTL  = ParagraphStyle('TDL', parent=sB, fontSize=8, alignment=TA_LEFT, leading=10, fontName='DejaVu')

# Trennseiten-Styles
sTrT = ParagraphStyle('TrT', parent=sT, fontSize=16, spaceBefore=0, spaceAfter=4)
sTrS = ParagraphStyle('TrS', parent=sAb, fontSize=9, alignment=TA_CENTER, spaceAfter=0)

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
        ('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
    return t


# ══════════════════════════════════════════════════════════
# TEIL 1: Einleitung als 6×9 PDF
# ══════════════════════════════════════════════════════════
def build_einleitung(outpath):
    pagesize = (W_PT, H_PT)
    doc = SimpleDocTemplate(outpath, pagesize=pagesize,
                            leftMargin=MARGIN_LR, rightMargin=MARGIN_LR,
                            topMargin=MARGIN_TB, bottomMargin=MARGIN_TB)
    story = []

    # Titelseite
    story.append(Spacer(1, 45*mm))
    story.append(Paragraph('Harmonikabau',
                 ParagraphStyle('BT', parent=sT, fontSize=28, spaceAfter=4)))
    story.append(Paragraph('Technische Dokumentation',
                 ParagraphStyle('BST', parent=sST, fontSize=14, spaceAfter=6)))
    story.append(Spacer(1,2*mm)); story.append(hr()); story.append(Spacer(1,4*mm))
    story.append(Paragraph('Akustik, Konstruktion und Praxis',
                 ParagraphStyle('Sub', parent=sAb, fontSize=11, alignment=TA_CENTER)))
    story.append(Spacer(1, 35*mm))
    story.append(Paragraph('Johann Pascher',
                 ParagraphStyle('Au', parent=sST, fontSize=12, spaceAfter=3)))
    story.append(Paragraph('Linz, Österreich — 2025',
                 ParagraphStyle('Yr', parent=sAb, fontSize=9, alignment=TA_CENTER)))
    story.append(PageBreak())

    # Einleitung
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

    # Zum Autor
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

    # Zur Entstehung
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

    # Dokumentenübersicht
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
    ], cw=tcw))
    story.append(Spacer(1, 2*mm))

    story.append(Paragraph('Spezialthemen', sSCh))
    story.append(mk_tbl(['Dok.','Titel'], [
        ['0500','Leitfaden zur ästhetischen Forensik bei Akkordeon-Gehäusen'],
    ], cw=tcw))
    story.append(Spacer(1, 4*mm))

    # Hinweise
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
    print(f'  ✓ Einleitung 6×9: {outpath}')


# ══════════════════════════════════════════════════════════
# TEIL 2: Schlusswort als 6×9 PDF
# ══════════════════════════════════════════════════════════
def build_schlusswort(outpath):
    pagesize = (W_PT, H_PT)
    doc = SimpleDocTemplate(outpath, pagesize=pagesize,
                            leftMargin=MARGIN_LR, rightMargin=MARGIN_LR,
                            topMargin=MARGIN_TB, bottomMargin=MARGIN_TB)
    story = []
    story.append(Spacer(1, 8*mm))
    story.append(Paragraph('Dok. 9999', sST))
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
    print(f'  ✓ Schlusswort 6×9: {outpath}')


# ══════════════════════════════════════════════════════════
# TEIL 3: Trennseiten-PDF erzeugen
# ══════════════════════════════════════════════════════════
def build_separator(outpath, dok_nr, titel, untertitel=''):
    """Erzeugt eine einzelne Trennseite im 6×9-Format."""
    pagesize = (W_PT, H_PT)
    doc = SimpleDocTemplate(outpath, pagesize=pagesize,
                            leftMargin=MARGIN_LR, rightMargin=MARGIN_LR,
                            topMargin=MARGIN_TB, bottomMargin=MARGIN_TB)
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
# TEIL 4: Seite auf 6×9 skalieren
# ══════════════════════════════════════════════════════════
def scale_page_to_6x9(page):
    """Skaliert eine Seite proportional auf 6×9 Zoll.
    Inhalt wird zentriert; die Bundsteg-Verschiebung erfolgt nachher global."""
    SAFE_H = 10 * 72 / 25.4   # 10mm seitlich (Verschiebung kommt danach)
    SAFE_V =  8 * 72 / 25.4   #  8mm oben/unten
    target_w = W_PT - 2 * SAFE_H
    target_h = H_PT - 2 * SAFE_V

    mb = page.mediabox
    src_w = float(mb.width)
    src_h = float(mb.height)
    scale_x = target_w / src_w
    scale_y = target_h / src_h
    scale = min(scale_x, scale_y)
    new_w = src_w * scale
    new_h = src_h * scale
    tx = (W_PT - new_w) / 2
    ty = (H_PT - new_h) / 2

    new_page = PageObject.create_blank_page(width=W_PT, height=H_PT)
    op = Transformation().scale(scale, scale).translate(tx / scale, ty / scale)
    page.add_transformation(op)
    new_page.merge_page(page)
    return new_page


def apply_gutter_shift(writer):
    """Verschiebt alle Seiten für zweiseitiges Layout.
    Ungerade Seiten (rechts): Inhalt nach rechts (Bundsteg links).
    Gerade Seiten (links): Inhalt nach links (Bundsteg rechts)."""
    shift = float(GUTTER_SHIFT)
    for i in range(len(writer.pages)):
        page = writer.pages[i]
        if i % 2 == 0:  # Seite 1,3,5... (0-indexed gerade = ungerade Seitenzahl)
            dx = shift    # nach rechts
        else:            # Seite 2,4,6...
            dx = -shift   # nach links
        op = Transformation().translate(dx, 0)
        page.add_transformation(op)


# ══════════════════════════════════════════════════════════
# TEIL 5: Zusammenbauen
# ══════════════════════════════════════════════════════════
def main():
    basedir = os.path.dirname(os.path.abspath(__file__))
    tmpdir = os.path.join(basedir, '_tmp_kindle')
    os.makedirs(tmpdir, exist_ok=True)

    # 1) Einleitung und Schlusswort im 6×9-Format erzeugen
    einl_path = os.path.join(tmpdir, 'einleitung_6x9.pdf')
    schl_path = os.path.join(tmpdir, 'schlusswort_6x9.pdf')
    build_einleitung(einl_path)
    build_schlusswort(schl_path)

    # 2) Dokumentliste mit Titeln (für Trennseiten)
    dokumente = [
        ('bass_50hz_v8.pdf',                     '0002', 'Strömungsanalyse Bass-Stimmzunge 50\u202FHz', 'Version 8'),
        ('0003_impedanz_vergleich_De.pdf',        '0003', 'Impedanzvergleich', 'Durchschlagzunge vs. Labialpfeife'),
        ('0004_frequenzvariation_zwei_filter_De.pdf','0004','Frequenzvariation der Stimmzunge', 'durch Kammerkopplung'),
        ('0005_ansprache_frequenz_kopplung_De.pdf','0005', 'Frequenzverschiebung', 'als Indikator der Ansprache'),
        ('0006_zeichenerklaerung_De.pdf',         '0006', 'Zeichenerklärung', ''),
        ('0007_diskant_kammerfrequenzen_De.pdf',  '0007', 'Diskant-Stimmstock', 'Kammerfrequenzen D3–C6'),
        ('0008_klangveraenderung_De.pdf',         '0008', 'Klangveränderung', 'durch Kammergeometrie'),
        ('0009_frequenzkopplung_mehrere_zungen_De.pdf','0009','Frequenzkopplung', 'mehrerer Zungen'),
        ('0010_guete_stimmplatte_De.pdf',         '0010', 'Güte der Stimmzunge', ''),
        ('0011_kanalgeometrie_De.pdf',            '0011', 'Kanalgeometrie der Stimmplatte', ''),
        ('0012_steifigkeit_De.pdf',               '0012', 'Zungensteifigkeit', ''),
        ('0015_kopplung_De.pdf',                  '0015', 'Akustische Kopplung', 'und Impedanzanpassung'),
        ('0016_obertonmoden_De.pdf',              '0016', 'Obertonmoden der Basszunge', 'Profilierung und Inharmonizität'),
        ('0017_diskant_De.pdf',                   '0017', 'Diskant-Stimmzungen', 'Obertonmoden F3 bis C6'),
        ('0018_hoerbarkeit_De.pdf',               '0018', 'Hörbarkeit der Biegemoden', 'Kammer-Saugkreis, Transiente, Torsion'),
        ('0019_stimmung_De.pdf',                  '0019', 'Stimmung und Differenztöne', ''),
        ('0020_tremolo_De.pdf',                   '0020', 'Tremolo', 'Schwebungsphysik, Typen und Stimmungspraxis'),
        ('0021_stimmplatten_De.pdf',              '0021', 'Stimmplatten', 'Qualität, Hersteller und Güteklassen'),
        ('0022_balg_De.pdf',                      '0022', 'Balg', 'Querschnitt, Faltenzahl und Instrumentengröße'),
        ('0023_gehaeuse_De.pdf',                  '0023', 'Gehäuse und Mechanik', ''),
        ('0500_forensik_De.pdf',                  '0500', 'Ästhetische Forensik', 'bei Akkordeon-Gehäusen'),
    ]

    # 3) Gesamt-PDF zusammenbauen
    writer = PdfWriter()

    # Einleitung (bereits 6×9)
    print('Zusammenbau ...')
    reader = PdfReader(einl_path)
    for page in reader.pages:
        writer.add_page(page)
    print(f'  + Einleitung ({len(reader.pages)} Seiten)')

    # Einzeldokumente
    for filename, dok_nr, titel, untertitel in dokumente:
        filepath = os.path.join(basedir, filename)
        if not os.path.exists(filepath):
            print(f'  ⚠ {filename} nicht gefunden — übersprungen')
            continue

        # Trennseite
        sep_path = os.path.join(tmpdir, f'sep_{dok_nr}.pdf')
        build_separator(sep_path, dok_nr, titel, untertitel)
        sep_reader = PdfReader(sep_path)
        for p in sep_reader.pages:
            writer.add_page(p)

        # Dokumentseiten skaliert einfügen
        doc_reader = PdfReader(filepath)
        n = len(doc_reader.pages)
        for page in doc_reader.pages:
            scaled = scale_page_to_6x9(page)
            writer.add_page(scaled)
        print(f'  + Dok. {dok_nr}: {n} Seiten')

    # Schlusswort (bereits 6×9)
    reader = PdfReader(schl_path)
    for page in reader.pages:
        writer.add_page(page)
    print(f'  + Schlusswort ({len(reader.pages)} Seiten)')

    # 4) Bundsteg-Verschiebung für zweiseitiges Layout
    print('Bundsteg-Verschiebung (zweiseitig) ...')
    apply_gutter_shift(writer)

    # 5) Schreiben
    outfile = os.path.join(basedir, 'Harmonikabau_Kindle_6x9.pdf')
    with open(outfile, 'wb') as f:
        writer.write(f)

    total = len(writer.pages)
    print(f'\n✓ {outfile}')
    print(f'  {total} Seiten, Format 6×9 Zoll ({W_PT}×{H_PT} pt)')

    # Aufräumen
    import shutil
    shutil.rmtree(tmpdir, ignore_errors=True)


if __name__ == '__main__':
    main()
