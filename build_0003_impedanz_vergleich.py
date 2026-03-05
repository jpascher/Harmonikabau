#!/usr/bin/env python3
"""Generate PDF: Impedanzvergleich Durchschlagzunge vs. Labialpfeife"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm, cm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor, white, black
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, KeepTogether, HRFlowable
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import os

# Colors matching v8 LaTeX style
DARKBLUE = HexColor('#16213e')
ACCENTRED = HexColor('#e94560')
KEYGREEN = HexColor('#2e7d32')
WARNRED = HexColor('#c62828')
LIGHTGREEN = HexColor('#e8f5e9')
LIGHTWARN = HexColor('#ffebee')
LIGHTGRAY = HexColor('#f5f5f5')
MIDGRAY = HexColor('#9e9e9e')

# Try to register a font with German characters
try:
    pdfmetrics.registerFont(TTFont('DejaVu', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuBold', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuItalic', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuBoldItalic', '/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
    FONT = 'DejaVu'
    FONTB = 'DejaVuBold'
    FONTI = 'DejaVuItalic'
    FONTBI = 'DejaVuBoldItalic'
except:
    FONT = 'Helvetica'
    FONTB = 'Helvetica-Bold'
    FONTI = 'Helvetica-Oblique'
    FONTBI = 'Helvetica-BoldOblique'

# PDF setup
outpath = '/mnt/user-data/outputs/0003_impedanz_vergleich_De.pdf'
doc = SimpleDocTemplate(
    outpath, pagesize=A4,
    leftMargin=25*mm, rightMargin=25*mm,
    topMargin=25*mm, bottomMargin=25*mm
)

# Styles
styles = getSampleStyleSheet()

def make_style(name, parent='Normal', **kwargs):
    base = styles[parent]
    return ParagraphStyle(name, parent=base, fontName=kwargs.get('fontName', FONT),
                          fontSize=kwargs.get('fontSize', 10),
                          leading=kwargs.get('leading', 14),
                          spaceAfter=kwargs.get('spaceAfter', 6),
                          spaceBefore=kwargs.get('spaceBefore', 0),
                          alignment=kwargs.get('alignment', TA_JUSTIFY),
                          textColor=kwargs.get('textColor', black),
                          firstLineIndent=kwargs.get('firstLineIndent', 0),
                          leftIndent=kwargs.get('leftIndent', 0),
                          rightIndent=kwargs.get('rightIndent', 0))

sTitle = make_style('sTitle', fontSize=18, fontName=FONTB, alignment=TA_CENTER,
                     textColor=DARKBLUE, spaceAfter=4, leading=22)
sSubtitle = make_style('sSubtitle', fontSize=11, alignment=TA_CENTER,
                        textColor=DARKBLUE, spaceAfter=2, leading=14)
sNote = make_style('sNote', fontSize=8.5, fontName=FONTI, alignment=TA_LEFT,
                    textColor=MIDGRAY, spaceAfter=8, leading=11)
sSection = make_style('sSection', fontSize=14, fontName=FONTB, textColor=DARKBLUE,
                       spaceBefore=14, spaceAfter=6, leading=18)
sSubsec = make_style('sSubsec', fontSize=11.5, fontName=FONTB, textColor=HexColor('#2a4a7f'),
                      spaceBefore=10, spaceAfter=4, leading=15)
sBody = make_style('sBody', fontSize=10, leading=14, spaceAfter=6)
sBodyBold = make_style('sBodyBold', fontSize=10, fontName=FONTB, leading=14, spaceAfter=6)
sKey = make_style('sKey', fontSize=9.5, leading=13, spaceAfter=4, leftIndent=4*mm, rightIndent=4*mm)
sWarn = make_style('sWarn', fontSize=9.5, leading=13, spaceAfter=4, leftIndent=4*mm, rightIndent=4*mm)
sFormula = make_style('sFormula', fontSize=10, alignment=TA_CENTER, fontName=FONTI,
                       spaceBefore=4, spaceAfter=6, leading=14)
sCaption = make_style('sCaption', fontSize=8.5, fontName=FONTI, alignment=TA_CENTER,
                       spaceAfter=8, leading=11)

story = []

# ============ TITLE ============
story.append(Paragraph('Wie die Stimmzunge die Impedanz verändert', sTitle))
story.append(Paragraph('Systematischer Vergleich: Durchschlagzunge (Akkordeon) vs. Labialpfeife (Orgel)', sSubtitle))
story.append(HRFlowable(width='80%', thickness=2, color=ACCENTRED, spaceAfter=4, spaceBefore=4))
story.append(Paragraph(
    'Ergänzungskapitel zu bass_50hz_v8. Stellt den fundamentalen Impedanz-Unterschied zwischen '
    'einer geflöteten Orgelpfeife und einer Stimmzunge in der Basskammer systematisch dar. '
    'Zahlenwerte beziehen sich auf die 50-Hz-Basszunge aus v8.', sNote))
story.append(Spacer(1, 4*mm))


def section(num, title):
    story.append(Paragraph(f'Kapitel {num}: {title}', sSection))

def subsection(title):
    story.append(Paragraph(title, sSubsec))

def body(text):
    story.append(Paragraph(text, sBody))

def formula(text):
    story.append(Paragraph(text, sFormula))

def keybox(text):
    tbl = Table([[Paragraph(text, sKey)]], colWidths=[doc.width - 8*mm])
    tbl.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), LIGHTGREEN),
        ('BOX', (0,0), (-1,-1), 1, KEYGREEN),
        ('TOPPADDING', (0,0), (-1,-1), 3*mm),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3*mm),
        ('LEFTPADDING', (0,0), (-1,-1), 4*mm),
        ('RIGHTPADDING', (0,0), (-1,-1), 4*mm),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 3*mm))

def warnbox(text):
    tbl = Table([[Paragraph(text, sWarn)]], colWidths=[doc.width - 8*mm])
    tbl.setStyle(TableStyle([
        ('BACKGROUND', (0,0), (-1,-1), LIGHTWARN),
        ('BOX', (0,0), (-1,-1), 1, WARNRED),
        ('TOPPADDING', (0,0), (-1,-1), 3*mm),
        ('BOTTOMPADDING', (0,0), (-1,-1), 3*mm),
        ('LEFTPADDING', (0,0), (-1,-1), 4*mm),
        ('RIGHTPADDING', (0,0), (-1,-1), 4*mm),
    ]))
    story.append(tbl)
    story.append(Spacer(1, 3*mm))

def make_table(data, col_widths=None, caption=None):
    if col_widths is None:
        col_widths = [doc.width / len(data[0])] * len(data[0])
    # Wrap strings in Paragraphs
    sCell = make_style('sCell', fontSize=8.5, leading=11, spaceAfter=0, alignment=TA_LEFT)
    sCellB = make_style('sCellB', fontSize=8.5, fontName=FONTB, leading=11, spaceAfter=0, alignment=TA_LEFT)
    wrapped = []
    for i, row in enumerate(data):
        wr = []
        for cell in row:
            st = sCellB if i == 0 else sCell
            wr.append(Paragraph(str(cell), st))
        wrapped.append(wr)
    
    tbl = Table(wrapped, colWidths=col_widths, repeatRows=1)
    style_cmds = [
        ('BACKGROUND', (0,0), (-1,0), HexColor('#e8eaf6')),
        ('GRID', (0,0), (-1,-1), 0.5, MIDGRAY),
        ('TOPPADDING', (0,0), (-1,-1), 2*mm),
        ('BOTTOMPADDING', (0,0), (-1,-1), 2*mm),
        ('LEFTPADDING', (0,0), (-1,-1), 2*mm),
        ('RIGHTPADDING', (0,0), (-1,-1), 2*mm),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ]
    # Alternating row colors
    for i in range(1, len(data)):
        if i % 2 == 0:
            style_cmds.append(('BACKGROUND', (0,i), (-1,i), LIGHTGRAY))
    tbl.setStyle(TableStyle(style_cmds))
    story.append(tbl)
    if caption:
        story.append(Paragraph(caption, sCaption))
    story.append(Spacer(1, 2*mm))


# ============ CHAPTER 1 ============
section(1, 'Warum dieser Vergleich?')

body('Die Basskammer des Akkordeons und eine gefaltete Orgelpfeife haben äußerlich viel gemeinsam: '
     'ein gefalteter Luftweg, ein Resonanzvolumen, eine 180°-Faltung (Gehrungsbiegung), eine Öffnung '
     'zur Außenluft. Der Coltman-Effekt (Impedanzstörung an der Faltung, v8 Abschnitt 5.2) gilt für '
     'beide Systeme. Die akustische Verkürzung Δl ≈ 0,79 · D ist geometrisch, nicht vom '
     'Anregungsmechanismus abhängig.')

body('Aber die <b>Art der Tonerzeugung</b> ist fundamental verschieden — und damit die Rolle, die die '
     'Kammerimpedanz für das Schwingverhalten spielt. Dieser Unterschied erklärt, warum Orgelbauer ihre '
     'Rohre als Resonatoren optimieren, während Akkordeonbauer ihre Kammern möglichst akustisch '
     '„unsichtbar" halten wollen (v8 Kapitel 10).')


# ============ CHAPTER 2 ============
section(2, 'Die Labialpfeife: Das Rohr bestimmt die Frequenz')

subsection('Anregungsmechanismus')
body('Bei einer geflöteten (labialen) Orgelpfeife erzeugt ein flacher Luftstrahl (Jet) aus dem Kernspalt '
     'eine Instabilität am Labium (Schneide). Der Jet schwingt seitlich hin und her — getrieben durch die '
     '<b>akustische Rückkopplung aus dem Rohr</b>. Die Frequenz des Tons wird <b>nicht</b> vom Jet bestimmt '
     '(der hat keine Eigenfrequenz im musikalischen Sinn), sondern von der Rohrresonanz.')

subsection('Eingangsimpedanz des Rohres')
body('Das Rohr wirkt als akustischer Resonator mit frequenzabhängiger Eingangsimpedanz:')
formula('Z<sub>Rohr</sub>(f) = j · ρc · (S<sub>Mund</sub>/S<sub>Rohr</sub>) · tan(2πf · L<sub>eff</sub> / c)')
body('(gedackte Pfeife; bei offener Pfeife cot statt tan). Die Impedanz hat <b>scharfe Maxima</b> '
     'bei den Rohrresonanzen (λ/4 bei gedackt, λ/2 bei offen).')

subsection('Kopplung: Jet als Geschwindigkeitsgenerator')
body('Der Jet am Labium ist ein <b>Geschwindigkeitsgenerator</b> (hohe Schnelle, niedrige Impedanz). '
     'Er „sieht" das Rohr als Last. Die Schwingung entsteht dort, wo die <b>Rohrimpedanz maximal</b> '
     'ist — weil dort die größte Druckrückkopplung auf den Jet wirkt. Der Jet passt sich der Rohrlänge an.')

keybox('<b>Schlüsselaussage Labialpfeife:</b> Das Rohr ist der frequenzbestimmende Resonator. Der Jet ist '
       'eine breitbandige Energiequelle, die bei der Rohrresonanz maximal ankoppelt. Die Impedanz am '
       'Labium ist <b>zeitlich konstant</b> (lineare Kopplung) — das Rohrende bewegt sich nicht.')


# ============ CHAPTER 3 ============
section(3, 'Die Durchschlagzunge: Die Zunge bestimmt die Frequenz')

subsection('Anregungsmechanismus')
body('Die Stimmzunge ist ein <b>mechanischer Oszillator</b> mit eigener Eigenfrequenz (Euler-Bernoulli). '
     'Die 50-Hz-Basszunge in v8 hat:')
formula('f<sub>Zunge</sub> = (1,875<super>2</super> / 2π L<super>2</super>) · √(EI / ρ<sub>s</sub>A) ≈ 50 Hz')
body('Die Frequenz wird durch Länge, Steifigkeit und Masse der Zunge bestimmt — <b>nicht</b> durch die '
     'Kammer. Die Kammer hat f<sub>H</sub> = 274 Hz — weit entfernt von 50 Hz.')

subsection('Der Spalt als nichtlinearer Impedanzwandler')
body('Hier liegt der fundamentale Unterschied. Der Spalt zwischen Zunge und Stimmplatte verwandelt '
     'die mechanische Impedanz der Zunge in eine <b>akustische Impedanz</b>, die die Kammer „sieht":')
formula('Z<sub>akust,Spalt</sub>(t) = ρ · v<sub>Spalt</sub>(t) / (2 · S<sub>Spalt</sub>(t))')
body('Da S<sub>Spalt</sub>(t) sich im Schwingzyklus um <b>Faktor 100</b> ändert (v8 Kapitel 8: '
     'von ~0 bis 58 mm²), ist die akustische Impedanz am Spalt <b>zeitvariabel und stark nichtlinear</b>.')

subsection('Mechanische Impedanz der Zunge')
body('Die Zunge selbst hat die mechanische Impedanz eines gedämpften Feder-Masse-Systems:')
formula('Z<sub>mech</sub>(f) = R<sub>Dämpfung</sub> + j(ω · m<sub>eff</sub> − k<sub>eff</sub>/ω)')
body('Bei f = f<sub>Zunge</sub>: Z<sub>mech</sub> = R<sub>Dämpfung</sub> (rein reell, Minimum). '
     'Die Güte Q = ω<sub>0</sub> · m<sub>eff</sub> / R bestimmt die Bandbreite.')

keybox('<b>Schlüsselaussage Durchschlagzunge:</b> Die Zunge ist der frequenzbestimmende Oszillator. '
       'Die Kammer ist <b>kein</b> Resonator im Sinne der Frequenzselektion, sondern ein '
       'Impedanz-Modulator, der über die Phase der Rückkopplung die Obertöne und die Einschwingzeit '
       'beeinflusst. Die Impedanz am Spalt ist <b>zeitvariabel</b> (nichtlineare Kopplung).')


# ============ CHAPTER 4 ============
section(4, 'Gegenüberstellung der Impedanzverhältnisse')

data = [
    ['Eigenschaft', 'Labialpfeife (Orgel)', 'Durchschlagzunge (Akkordeon)'],
    ['Frequenzbestimmung', 'Rohrresonanz (Impedanz-Maximum)', 'Mechanische Eigenfrequenz der Zunge'],
    ['Anregung', 'Jet am Labium (Geschwindigkeits-generator, breitbandig)',
     'Bernoulli-Sog am Spalt (selbsterregt bei f_Zunge)'],
    ['Impedanz am Anregungsort', 'Zeitlich konstant (Rohrmund bewegt sich nicht)',
     'Periodisch geschaltet: Faktor ~200 pro Zyklus'],
    ['Kopplungsart', 'Quasi-linear (kleine Jet-Auslenkung)',
     'Stark nichtlinear (S_Spalt ~ 0 ... 58 mm²)'],
    ['Kammer/Rohr-Rolle', 'Frequenzselektiv (Rohr ist das Instrument)',
     'Phasenkopplung (Kammer beeinflusst Obertöne, nicht Grundton)'],
    ['Kopplung (E-Technik)', 'Parallelresonanz (LC parallel)', 'Serienkopplung (Luft durch Kammer)'],
    ['Resonante Kammer', 'Erwünscht (Resonanz = Tonerzeugung)', 'Schädlich (entzieht Energie)'],
    ['Obertonstruktur', 'Vom Rohr erzwungen (harmonisch)', 'Von Zunge erzeugt, von Kammer gefiltert'],
    ['Coltman-Effekt', 'Korrektur der Rohrlänge → Intonation',
     'Korrektur von f_H → Phasenverschiebung der Obertöne'],
]
make_table(data, col_widths=[35*mm, 55*mm, 60*mm],
           caption='Tabelle 1: Impedanz-Vergleich: Labialpfeife vs. Durchschlagzunge')


# ============ CHAPTER 5 ============
section(5, 'Der Spalt als periodisch geschaltete Impedanz')

subsection('Elektrotechnische Analogie')
body('In der E-Technik-Analogie (v8 Kapitel 7) ist die Labialpfeife ein <b>linearer LC-Schwingkreis</b>, '
     'der von einer Konstantstromquelle (Jet) angeregt wird. Die Resonanzfrequenz '
     'f<sub>0</sub> = 1/(2π√LC) ist zeitunabhängig.')

body('Die Durchschlagzunge dagegen ist ein LC-Kreis, dessen <b>Abschlusswiderstand R(t) mit der '
     'Zungenfrequenz zwischen Null und Unendlich geschaltet wird</b>:')

formula('R<sub>Spalt</sub>(t) = ∞  wenn Spalt geschlossen;  '
        'R<sub>Spalt</sub>(t) = ρv / (2 · S<sub>Spalt</sub>(t))  wenn Spalt offen')

body('Das ist ein <b>parametrischer Oszillator</b>: ein System, bei dem ein Schaltkreis-Parameter (R) '
     'periodisch moduliert wird. Parametrische Oszillatoren können Energie bei <b>anderen Frequenzen</b> '
     'als der Anregungsfrequenz aufnehmen oder abgeben — was erklärt, warum die Kammerresonanz '
     'bei 274 Hz (weit von 50 Hz) trotzdem die Obertöne der 50-Hz-Zunge beeinflussen kann.')

subsection('Zahlenwerte für die 50-Hz-Basszunge')

data2 = [
    ['Phase im Zyklus', 'S_Spalt [mm²]', 'v_Spalt [m/s]', 'Z_akust [Pa·s/m³]'],
    ['Spalt geschlossen', '≈ 0', '—', '→ ∞'],
    ['Spalt minimal (Seitenspalt)', '7,5', '22', '1,8 × 10⁶'],
    ['Spalt Ruhelage', '4,0', '29', '4,3 × 10⁶'],
    ['Spalt maximal (Durchschwung)', '58', '2', '2,1 × 10⁴'],
]
make_table(data2, col_widths=[45*mm, 30*mm, 28*mm, 37*mm],
           caption='Tabelle 2: Akustische Spaltimpedanz in verschiedenen Phasen des Schwingzyklus (bei 1 kPa). '
                   'Faktor zwischen Minimum und Maximum: ~200.')

warnbox('Die Impedanz am Spalt ändert sich um mehr als <b>zwei Größenordnungen</b> innerhalb einer '
        'einzigen Schwingungsperiode von 20 ms. Bei einer Labialpfeife ist die entsprechende Schwankung '
        '&lt; 1 %. Dies ist der fundamentale physikalische Grund, warum lineare akustische Modelle '
        '(Helmholtz, Transmission-Line) für die Stimmzunge nur im Einschwingbereich (kleine Amplitude) '
        'zuverlässig sind.')


# ============ CHAPTER 6 ============
section(6, 'Konsequenz: Warum die Kammer anders optimiert werden muss')

subsection('Orgelpfeife: Resonanz erwünscht')
body('Der Orgelbauer <b>maximiert</b> die Resonanzgüte des Rohres. Ein höheres Q bedeutet: '
     'schärfere Frequenzselektion → reinerer Ton; stärkere Rückkopplung auf den Jet → stabilere '
     'Schwingung; weniger Energieverlust → lauterer Ton bei gleichem Winddruck. '
     'Deshalb werden Orgelpfeifen aus glattem Material (Zinn, Blei, Holz mit glatter Innenfläche) '
     'gefertigt, und die Coltman-Kompensation dient der <b>exakten Intonation</b>.')

subsection('Stimmzungenkammer: Resonanz schädlich')
body('Der Akkordeonbauer will das Gegenteil: Die Kammer soll möglichst <b>nicht resonieren</b>. '
     'Eine resonante Kammer bei einem Oberton der Zunge (v8 Kapitel 10): entzieht Energie über '
     'Serienkopplung → langsame Ansprache; erzwingt Phasenverschiebungen → verändert Klangfarbe; '
     'koppelt an die Nachbar-Zunge → Schwebung.')

body('Deshalb ist die Empfehlung aus v8: f<sub>H</sub> so wählen, dass kein Oberton der Zunge in '
     'das Resonanzband fällt. Die Coltman-Korrektur hat hier eine <b>andere Funktion</b> als im '
     'Orgelbau: Sie verschiebt nicht die Intonation (die bestimmt die Zunge), sondern die '
     '<b>Phasenlage der Impedanz bei den Obertönen</b>.')

subsection('Quantitativer Vergleich')

data3 = [
    ['Parameter', 'Orgelpfeife (C2, gedackt)', 'Basszunge (50 Hz)'],
    ['Effektive Rohrlänge', '~1300 mm (gedackt λ/4)', '169 mm (Kammerweg)'],
    ['f_Resonanz / f_Ton', '= 1,00 (Resonanz ist der Ton)', '= 5,48 (f_H / f_Zunge, weit entfernt)'],
    ['Resonanz-Güte Q', '20–50 (gewünscht: hoch)', '5–10 (gewünscht: niedrig)'],
    ['Akust. Länge der Faltung', 'Kritisch (Intonation)', 'Sekundär (f_H-Verschiebung)'],
    ['Impedanz am Mund', '~konstant (linear)', 'Faktor ~200 pro Zyklus'],
    ['Coltman-Korrektur Zweck', 'Rohrstimmung', 'Oberton-Phasenlage'],
]
make_table(data3, col_widths=[40*mm, 50*mm, 50*mm],
           caption='Tabelle 3: Quantitativer Vergleich Orgelpfeife (C2 gedackt) vs. Basszunge (50 Hz)')


# ============ CHAPTER 7 ============
section(7, 'Die drei Zeitskalen der Impedanz')

body('Die Analyse zeigt, dass bei der Durchschlagzunge <b>drei Zeitskalen</b> der Impedanz relevant '
     'sind, die bei der Labialpfeife keine Rolle spielen:')

body('<b>1. Schnelle Zeitskala — Spaltmodulation (T = 1/f<sub>Zunge</sub> = 20 ms):</b> '
     'Die Impedanz am Spalt wird 50-mal pro Sekunde zwischen fast Null und fast Unendlich geschaltet. '
     'Dieser Effekt erzeugt die nichtlineare Selbsterregung und bestimmt das Obertonspektrum. '
     'Bei der Labialpfeife gibt es kein Äquivalent — der Rohrmund hat konstanten Querschnitt.')

body('<b>2. Mittlere Zeitskala — Einschwingvorgang (τ = Q/(πf) ~ 100–600 ms):</b> '
     'Während des Amplitudenaufbaus wächst die Schwingungsamplitude, und die <b>zeitgemittelte</b> '
     'Spaltimpedanz sinkt (v8 Kapitel 8: mittlere Spaltfläche steigt von 4 mm² auf 32 mm²). '
     'Die Kammer „sieht" eine sich über hunderte Millisekunden langsam ändernde Last. '
     'In dieser Phase ist die Helmholtz-Näherung am besten gültig.')

body('<b>3. Langsame Zeitskala — Statische Spaltverengung (Balgdruck-Änderung):</b> '
     'Der Balgdruck biegt die Zunge statisch (v8 Kapitel 9). Bei steigendem Druck sinkt h<sub>eff</sub>, '
     'steigt Z<sub>Spalt,Ruhe</sub>, und der Arbeitspunkt der Schwingung verschiebt sich. '
     'Diese Zeitskala ist ~1 s (Balgbewegung des Spielers).')

body('Bei der Labialpfeife gibt es nur eine relevante Zeitskala: die (feste) akustische Periodendauer '
     'des Rohres. Der Winddruck beeinflusst die Amplitude, aber nicht die Impedanzverhältnisse am Labium.')


# ============ CHAPTER 8 ============
section(8, 'Warum lineare Modelle trotzdem nützlich sind')

body('Trotz der starken Nichtlinearität der Spaltkopplung ist das lineare Helmholtz-Modell aus v8 '
     '<b>nicht falsch</b> — es beschreibt die Kammer unter der Annahme, dass die <b>zeitgemittelte '
     'Spaltimpedanz</b> als konstanter Abschluss wirkt. Das ist gültig für:')

body('<b>Einschwingbereich (kleine Amplitude):</b> S<sub>Spalt</sub> ≈ 4 mm², wenig Modulation. '
     'Helmholtz und Phasenkopplung (v8 Kapitel 11) beschreiben korrekt, ob ein Oberton Energie gewinnt '
     'oder verliert.')

body('<b>Oberton-Analyse:</b> Die Frage, ob der 5. Oberton (250 Hz) bei f<sub>H</sub> = 274 Hz Energie '
     'gewinnt oder verliert, hängt von der <b>mittleren</b> Impedanz ab. Die schnelle Modulation mittelt '
     'sich über die viel längere Oberton-Periode heraus (50-Hz-Modulation bei 250-Hz-Abfrage → 5 Zyklen Mittelung).')

body('<b>Vergleich der Trennwandvarianten:</b> Alle drei Varianten erfahren dieselbe nichtlineare '
     'Spaltmodulation. Der <b>Unterschied</b> zwischen A, B und C liegt in der linearen Kammer-Impedanz '
     'Z(f) — und genau die beschreibt das Helmholtz-Modell.')

warnbox('<b>Wo das lineare Modell versagt:</b> Bei voller Amplitude (A ~ 20 mm) wird die zeitgemittelte '
        'Spaltfläche ~8-mal so groß (v8 Kapitel 8). Der Spalt ist nicht mehr der alleinige Engpass, '
        'die Klappe wird relevant, und die Kammer-Impedanz wirkt nicht mehr als Abschluss eines festen '
        'Halses, sondern als Abschluss eines periodisch umgeschalteten Ventils. In dieser Phase bräuchte '
        'man eine zeitaufgelöste FSI-Simulation.')


# ============ ZUSAMMENFASSUNG ============
section(9, 'Zusammenfassung')

keybox('<b>Kernaussage:</b> Die Stimmzunge verwandelt die Kammer-Impedanz von einem '
       '<b>frequenzselektierenden Resonator</b> (Orgelpfeife: Rohr bestimmt Frequenz, Impedanz zeitlich '
       'konstant, Resonanz erwünscht) in einen <b>phasenmodulierenden Energiepuffer</b> (Akkordeon: '
       'Zunge bestimmt Frequenz, Impedanz zeitvariabel um Faktor 200, Resonanz schädlich).<br/><br/>'
       'Der Spalt als nichtlinearer, zeitvariabler Impedanzwandler ist das entscheidende Element, das '
       'im Orgelpfeifen-Modell <b>kein Äquivalent</b> hat. In der E-Technik-Analogie: Die Orgelpfeife '
       'ist ein linearer Schwingkreis mit Konstantstrom-Anregung; die Stimmzunge ist ein parametrischer '
       'Oszillator, bei dem der Abschlusswiderstand mit 50 Hz zwischen Kurzschluss und Leerlauf '
       'geschaltet wird.<br/><br/>'
       'Für den Instrumentenbauer folgt: <b>Die Kammer ist kein Resonanzkörper</b> (wie bei der '
       'Orgelpfeife), sondern ein <b>Impedanz-Umfeld</b>, das die Zunge „sieht". Die Optimierung zielt '
       'nicht auf Resonanzverstärkung, sondern auf Phasenverträglichkeit: f<sub>H</sub> muss so liegen, '
       'dass kein Oberton der Zunge in das Resonanzband fällt (v8 Kapitel 12), und die Wandform '
       'moduliert die Bandbreite dieses Bandes (v8 Kapitel 11).')

story.append(Spacer(1, 6*mm))
story.append(HRFlowable(width='60%', thickness=1, color=ACCENTRED, spaceAfter=4, spaceBefore=2))
story.append(Paragraph('<i>Die Berechnung sagt, wo man suchen soll. Der praktische Test sagt, was man findet.</i>',
                        make_style('sClosing', fontSize=9, fontName=FONTI, alignment=TA_CENTER, leading=12)))


# ============ BUILD ============
def header_footer(canvas, doc):
    canvas.saveState()
    canvas.setFont(FONT, 8)
    canvas.setFillColor(MIDGRAY)
    canvas.drawString(25*mm, A4[1] - 15*mm, 'Dok. 0003 — Impedanzvergleich: Durchschlagzunge vs. Labialpfeife')
    canvas.drawRightString(A4[0] - 25*mm, A4[1] - 15*mm, f'Seite {doc.page}')
    canvas.drawCentredString(A4[0]/2, 12*mm, 'Dok. 0002 · 0004 · 0005 · 0006 · 0007 · 0008 · berechnungen_verifikation.py')
    canvas.restoreState()

doc.build(story, onFirstPage=header_footer, onLaterPages=header_footer)
print(f'PDF created: {outpath}')
