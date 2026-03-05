#!/usr/bin/env python3
"""Dok 0006: Zeichenerklärung — Formelzeichen aller Dokumente"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor, black
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, HRFlowable, KeepTogether
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

DARKBLUE = HexColor('#16213e')
ACCENTRED = HexColor('#e94560')
KEYGREEN = HexColor('#2e7d32')
WARNRED = HexColor('#c62828')
LIGHTGREEN = HexColor('#e8f5e9')
LIGHTWARN = HexColor('#ffebee')
LIGHTGRAY = HexColor('#f5f5f5')
MIDGRAY = HexColor('#9e9e9e')
LIGHTBLUE = HexColor('#e8eaf6')
LIGHTORANGE = HexColor('#fff3e0')

try:
    pdfmetrics.registerFont(TTFont('DejaVu', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuBold', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuItalic', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
    FONT='DejaVu'; FONTB='DejaVuBold'; FONTI='DejaVuItalic'
except:
    FONT='Helvetica'; FONTB='Helvetica-Bold'; FONTI='Helvetica-Oblique'

outpath = '/mnt/user-data/outputs/0006_zeichenerklaerung_De.pdf'
doc = SimpleDocTemplate(outpath, pagesize=A4,
    leftMargin=25*mm, rightMargin=25*mm, topMargin=25*mm, bottomMargin=25*mm)

styles = getSampleStyleSheet()
def ms(name, **kw):
    return ParagraphStyle(name, parent=styles['Normal'], fontName=kw.get('fn',FONT),
        fontSize=kw.get('fs',10), leading=kw.get('ld',14), spaceAfter=kw.get('sa',6),
        spaceBefore=kw.get('sb',0), alignment=kw.get('al',TA_JUSTIFY),
        textColor=kw.get('tc',black), leftIndent=kw.get('li',0), rightIndent=kw.get('ri',0))

sTitle = ms('sT',fn=FONTB,fs=17,al=TA_CENTER,tc=DARKBLUE,sa=4,ld=21)
sSubtitle = ms('sSt',fs=11,al=TA_CENTER,tc=DARKBLUE,sa=2,ld=14)
sNote = ms('sN',fs=8.5,fn=FONTI,tc=MIDGRAY,sa=8,ld=11)
sSection = ms('sSec',fs=13,fn=FONTB,tc=DARKBLUE,sb=12,sa=6,ld=17)
sSubsec = ms('sSub',fs=11,fn=FONTB,tc=HexColor('#2a4a7f'),sb=8,sa=4,ld=14)
sBody = ms('sB',fs=10,ld=14,sa=6)
sKey = ms('sK',fs=9.5,ld=13,sa=4,li=4*mm,ri=4*mm)
sWarn = ms('sW',fs=9.5,ld=13,sa=4,li=4*mm,ri=4*mm)
sCell = ms('sCell',fs=8.5,ld=11,sa=0,al=TA_LEFT)
sCellB = ms('sCellB',fs=8.5,fn=FONTB,ld=11,sa=0,al=TA_LEFT)
sCellSym = ms('sCellSym',fs=9.5,fn=FONTI,ld=12,sa=0,al=TA_LEFT)
sCaption = ms('sC',fs=8.5,fn=FONTI,al=TA_CENTER,sa=8,ld=11)

story = []

def section(t): story.append(Paragraph(t, sSection))
def subsection(t): story.append(Paragraph(t, sSubsec))
def body(t): story.append(Paragraph(t, sBody))

def keybox(t):
    tbl = Table([[Paragraph(t, sKey)]], colWidths=[doc.width-8*mm])
    tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),LIGHTGREEN),('BOX',(0,0),(-1,-1),1,KEYGREEN),
        ('TOPPADDING',(0,0),(-1,-1),3*mm),('BOTTOMPADDING',(0,0),(-1,-1),3*mm),
        ('LEFTPADDING',(0,0),(-1,-1),4*mm),('RIGHTPADDING',(0,0),(-1,-1),4*mm)]))
    story.append(tbl); story.append(Spacer(1,3*mm))

def warnbox(t):
    tbl = Table([[Paragraph(t, sWarn)]], colWidths=[doc.width-8*mm])
    tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),LIGHTWARN),('BOX',(0,0),(-1,-1),1,WARNRED),
        ('TOPPADDING',(0,0),(-1,-1),3*mm),('BOTTOMPADDING',(0,0),(-1,-1),3*mm),
        ('LEFTPADDING',(0,0),(-1,-1),4*mm),('RIGHTPADDING',(0,0),(-1,-1),4*mm)]))
    story.append(tbl); story.append(Spacer(1,3*mm))

def infobox(t, bg=LIGHTORANGE, fr=HexColor('#e65100')):
    tbl = Table([[Paragraph(t, sKey)]], colWidths=[doc.width-8*mm])
    tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),bg),('BOX',(0,0),(-1,-1),1,fr),
        ('TOPPADDING',(0,0),(-1,-1),3*mm),('BOTTOMPADDING',(0,0),(-1,-1),3*mm),
        ('LEFTPADDING',(0,0),(-1,-1),4*mm),('RIGHTPADDING',(0,0),(-1,-1),4*mm)]))
    story.append(tbl); story.append(Spacer(1,3*mm))


# Column widths for symbol tables
CW_SYM = 28*mm   # Symbol
CW_UNIT = 18*mm   # Einheit
CW_DESC = 60*mm   # Beschreibung
CW_DOK = 28*mm    # Dokumente
CW_TOTAL = CW_SYM + CW_UNIT + CW_DESC + CW_DOK

def sym_table(title, rows, caption=None):
    """Create a symbol table with header."""
    header = [
        Paragraph('<b>Zeichen</b>', sCellB),
        Paragraph('<b>Einheit</b>', sCellB),
        Paragraph('<b>Beschreibung</b>', sCellB),
        Paragraph('<b>Dok.</b>', sCellB),
    ]
    data = [header]
    for row in rows:
        data.append([
            Paragraph(row[0], sCellSym),
            Paragraph(row[1], sCell),
            Paragraph(row[2], sCell),
            Paragraph(row[3], sCell),
        ])
    tbl = Table(data, colWidths=[CW_SYM, CW_UNIT, CW_DESC, CW_DOK], repeatRows=1)
    sc = [
        ('BACKGROUND', (0,0), (-1,0), LIGHTBLUE),
        ('GRID', (0,0), (-1,-1), 0.4, MIDGRAY),
        ('TOPPADDING', (0,0), (-1,-1), 1.5*mm),
        ('BOTTOMPADDING', (0,0), (-1,-1), 1.5*mm),
        ('LEFTPADDING', (0,0), (-1,-1), 2*mm),
        ('RIGHTPADDING', (0,0), (-1,-1), 2*mm),
        ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ]
    for i in range(1, len(data)):
        if i % 2 == 0:
            sc.append(('BACKGROUND', (0,i), (-1,i), LIGHTGRAY))
    tbl.setStyle(TableStyle(sc))
    story.append(tbl)
    if caption:
        story.append(Paragraph(caption, sCaption))
    story.append(Spacer(1, 2*mm))


# ==========================================================================
# TITLE
# ==========================================================================
story.append(Paragraph('Dok. 0006: Zeichenerklärung', sTitle))
story.append(Paragraph('Formelzeichen, Indizes und Konventionen aller Dokumente', sSubtitle))
story.append(HRFlowable(width='80%', thickness=2, color=ACCENTRED, spaceAfter=4, spaceBefore=4))
story.append(Paragraph(
    'Vollständige Zeichenerklärung für alle Formelzeichen in Dok. 0001–0008. '
    'Geordnet nach physikalischem Gebiet. Jedes Zeichen verweist auf die Dokumente, '
    'in denen es erstmals definiert oder wesentlich verwendet wird.', sNote))

story.append(Spacer(1, 2*mm))
ref_data = [
    ['Dok. 0001', 'Analyse zur Instrumentenakustik und Wahrnehmungspsychologie'],
    ['Dok. 0002', 'Strömungsanalyse Bass-Stimmzunge 50 Hz — v8 (Hauptdokument)'],
    ['Dok. 0003', 'Impedanzvergleich Durchschlagzunge vs. Labialpfeife'],
    ['Dok. 0004', 'Frequenzvariation — Zwei-Filter-Modell'],
    ['Dok. 0005', 'Frequenzverschiebung als Indikator der Ansprache'],
    ['<b>Dok. 0006</b>', '<b>Zeichenerklärung (dieses Dokument)</b>'],
    ['Dok. 0007', 'Diskant-Stimmstock — Kammerfrequenzen'],
    ['Dok. 0008', 'Klangveränderung durch Kammergeometrie'],
    ['Dok. 0500', 'Leitfaden zur ästhetischen Forensik bei Akkordeon-Gehäusen'],
    ['', '<a href="berechnungen_verifikation.py">berechnungen_verifikation.py</a> — Verifikationsskript'],
]
ref_tbl = Table([[Paragraph(r[0], sCell), Paragraph(r[1], sCell)] for r in ref_data],
                colWidths=[25*mm, doc.width-25*mm-4*mm])
ref_tbl.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm),('BACKGROUND',(0,5),(-1,5),LIGHTGREEN)]))
story.append(ref_tbl)
story.append(Spacer(1, 4*mm))


# ==========================================================================
section('1. Geometrische Größen')
# ==========================================================================

sym_table('Geometrie', [
    # Symbol, Unit, Description, Docs
    ['L', 'mm', 'Freie Zungenlänge (70 mm)', '0002'],
    ['W', 'mm', 'Zungenbreite (8 mm)', '0002'],
    ['t', 'mm', 'Zungendicke (0,354 mm)', '0002'],
    ['h, h_max', 'mm', 'Aufbiegung am freien Ende (Ruhelage: 1,5 mm)', '0002, 0004'],
    ['h_eff', 'mm', 'Effektive Aufbiegung unter Balgdruck', '0002'],
    ['S_Spalt', 'mm²', 'Effektive Spaltfläche = W·h_max/3 (Cantilever-Profil)', '0002, 0004'],
    ['S_eff', 'mm²', 'Zeitabhängige effektive Durchströmfläche', '0002'],
    ['A_eff', 'm²', 'Effektive Kopplungsfläche Spalt = S_Spalt (in SI)', '0004, 0005'],
    ['S_Klappe', 'mm²', 'Effektive Klappenöffnungsfläche (244 mm²)', '0002'],
    ['W_Schlitz', 'mm', 'Schlitzbreite in der Stimmplatte (9 mm)', '0002'],
    ['W_k', 'mm', 'Kanalbreite (23,8 mm bei Furnier-Trennwand)', '0002'],
    ['H_Kam', 'mm', 'Kammerhöhe (40 mm)', '0002'],
    ['V', 'cm³', 'Kammervolumen (180 cm³)', '0002, 0004, 0005'],
    ['V_VK', 'cm³', 'Vorkammervolumen (67,2 cm³)', '0002'],
    ['gap_fold', 'mm', 'Faltungsspalt am Trennwandende (19 mm)', '0002'],
    ['L_eff', 'mm', 'Effektive akustische Länge', '0002, 0003'],
    ['Δl', 'mm', 'Akustische Verkürzung durch Faltung (≈ 15 mm)', '0002, 0003'],
    ['d_jet', 'mm', 'Äquivalenter Strahldurchmesser am Klappeneingang', '0002'],
    ['L_core', 'mm', 'Strahlkernlänge (≈ 5·d_jet)', '0002'],
    ['R', 'mm', 'Viertelkreis-Radius (Eckenverrundung)', '0002'],
    ['x', 'mm', 'Auslenkung der Zungenspitze (zeitabhängig)', '0002, 0004'],
    ['A (Amplitude)', 'mm', 'Schwingungsamplitude der Zungenspitze', '0002'],
], caption='Tabelle 1: Geometrische Größen. Zahlenwerte in Klammern beziehen sich auf die 50-Hz-Basszunge aus Dok. 0002.')


# ==========================================================================
section('2. Materialeigenschaften und Stoffgrößen')
# ==========================================================================

sym_table('Material', [
    ['E', 'Pa', 'Elastizitätsmodul der Zunge (Federstahl: ≈ 200 GPa)', '0002'],
    ['I', 'm⁴', 'Flächenträgheitsmoment = W·t³/12', '0002'],
    ['ρ_s', 'kg/m³', 'Dichte Zungenmaterial (Stahl: 7800)', '0002, 0004'],
    ['ρ', 'kg/m³', 'Dichte Luft (1,2 kg/m³ bei 20 °C)', '0002–0005'],
    ['c', 'm/s', 'Schallgeschwindigkeit in Luft (343 m/s bei 20 °C)', '0002–0005'],
    ['μ', 'Pa·s', 'Dynamische Viskosität der Luft (1,8×10⁻⁵)', '0002'],
    ['η', '—', 'Materialverlustfaktor (innere Reibung)', '0002'],
    ['m_total', 'g', 'Gesamtmasse der Zunge (= ρ_s·L·W·t)', '0004'],
    ['m_eff', 'g', 'Effektive Masse im 1. Cantilever-Mode (0,2427·m_total)', '0004'],
    ['k_mech', 'N/m', 'Mechanische Federsteifigkeit = ω₁²·m_eff', '0004, 0005'],
], caption='Tabelle 2: Materialeigenschaften und abgeleitete mechanische Größen')


# ==========================================================================
section('3. Strömungsgrößen')
# ==========================================================================

sym_table('Strömung', [
    ['Δp', 'Pa', 'Balgdruck (Druckdifferenz Kammer–Außenluft)', '0002'],
    ['Δp_min', 'Pa', 'Schwellendruck für Selbsterregung (≈ 224 Pa)', '0002'],
    ['Δp_crit', 'Pa', 'Grenzdruck für statische Spaltschließung (≈ 221 Pa)', '0002'],
    ['Δp_dyn', 'Pa', 'Praktischer Einsetzdruck (≈ 90 Pa, 30–50 % von Δp_min)', '0002'],
    ['Δp_B', 'Pa', 'Bernoulli-Druck am Spalt = ½ρ·v_tan²', '0002'],
    ['v_Spalt', 'm/s', 'Strömungsgeschwindigkeit im Spalt = √(2Δp/ρ)', '0002'],
    ['v_tan', 'm/s', 'Tangentialkomponente der Geschwindigkeit am Spalt', '0002'],
    ['v_jet', 'm/s', 'Strahlgeschwindigkeit am Klappeneingang', '0002'],
    ['q, Q (Volumenstrom)', 'ml/s', 'Volumenstrom durch den Spalt (≈ 115 ml/s bei 1 kPa)', '0002, 0004'],
    ['Re', '—', 'Reynolds-Zahl (Re_Spalt ≈ 950, Re_Kanal ≈ 239)', '0002'],
    ['f (Reibungsfaktor)', '—', 'Darcy-Reibungsfaktor (= 64/Re bei laminar)', '0002'],
    ['w_tip', 'mm', 'Statische Durchbiegung der Zungenspitze unter Balgdruck', '0002'],
], caption='Tabelle 3: Strömungsmechanische Größen')


# ==========================================================================
section('4. Verlustbeiwerte')
# ==========================================================================

sym_table('Verluste', [
    ['ζ', '—', 'Verlustbeiwert (allgemein)', '0002'],
    ['ζ_ges', '—', 'Gesamtverlustbeiwert des Luftwegs', '0002'],
    ['ζ_Klappe', '—', 'Verlustbeiwert der Klappenöffnung (Ausblas: 0,50; Umlenk: 1,70)', '0002'],
    ['ζ_entry', '—', 'Schlitz-Eintrittsverlust (0,5; neu in v8 als K3b)', '0002'],
    ['ζ_Kopplung', '—', 'Kopplungsverlust Vorkammer–Hauptkammer (4×10⁻⁴)', '0002'],
], caption='Tabelle 4: Strömungsmechanische Verlustbeiwerte')


# ==========================================================================
section('5. Richtungserhaltung und Wandform')
# ==========================================================================

sym_table('Richtung', [
    ['α', '°', 'Eintrittswinkel des Luftstrahls (gegen Horizontale)', '0002'],
    ['β (Wandwinkel)', '°', 'Halbwinkel der schrägen Trennwand (Variante B: 3,7°)', '0002'],
    ['η_ges', '—', 'Gesamter Richtungserhaltungsfaktor am Spalt', '0002'],
    ['η_Faltung', '—', 'Richtungserhaltung der 180°-Faltung (≈ 0,6)', '0002'],
    ['η_Wand', '—', 'Richtungserhaltung der Trennwand (A: 0,3; B: 0,5; C: 0,7)', '0002'],
    ['f_tan', '—', 'Tangential-Fraktion am Spalt = (1 + η_ges)/2', '0002'],
], caption='Tabelle 5: Richtungserhaltung. Variante A = gerade, B = schräg, C = parabolisch.')


# ==========================================================================
section('6. Frequenzen und Schwingungsgrößen')
# ==========================================================================

sym_table('Frequenzen', [
    ['f_1, f_Zunge, f_mech', 'Hz', 'Mechanische Eigenfrequenz der Zunge (50 Hz)', '0002–0005'],
    ['f_H', 'Hz', 'Helmholtz-Resonanzfrequenz der Kammer (274 Hz)', '0002–0005'],
    ['f_n = n·f_1', 'Hz', 'n-ter Oberton der Zunge', '0002, 0004, 0005'],
    ['f_frei', 'Hz', 'Zungenfrequenz ohne Kammer (auf Prüfblock)', '0005'],
    ['f_Kammer', 'Hz', 'Zungenfrequenz in der Kammer montiert', '0005'],
    ['f_1\', f_2\'', 'Hz', 'Gekoppelte Eigenfrequenzen (Zungen-/Kammermode)', '0004'],
    ['ω, ω_1, ω_H', 'rad/s', 'Kreisfrequenz = 2π·f', '0002–0005'],
    ['n', '—', 'Oberton-Nummer (ganzzahlig)', '0002, 0004, 0005'],
    ['λ', 'm', 'Wellenlänge = c/f', '0002, 0003'],
    ['Δf', 'Hz, Cent', 'Frequenzverschiebung durch Kammerkopplung', '0004, 0005'],
    ['Δf_n', 'Cent', 'Frequenzverschiebung des n-ten Obertons', '0005'],
    ['Δf_netto', 'Cent', 'Netto-Frequenzverschiebung (Summe über alle OT)', '0005'],
], caption='Tabelle 6: Frequenzen und Schwingungsgrößen. 1 Cent = 1/100 Halbton = 1/1200 Oktave.')


# ==========================================================================
section('7. Dämpfung und Güte')
# ==========================================================================

sym_table('Dämpfung', [
    ['Q, Q_1', '—', 'Mechanische Güte der Zunge (50–150)', '0002, 0004, 0005'],
    ['Q_H', '—', 'Akustische Güte des Helmholtz-Resonators (5–10)', '0002, 0004, 0005'],
    ['Q_ges', '—', 'Gesamtgüte (harmonische Summe aller Verluste)', '0002'],
    ['Q_Mat', '—', 'Materialdämpfung (Stahl: 1000–5000)', '0002'],
    ['Q_Luft', '—', 'Luftdämpfung (Squeeze-Film im Spalt)', '0002'],
    ['Q_Einsp', '—', 'Einspannungsdämpfung (Schrauben/Nieten)', '0002'],
    ['Q_Strahl', '—', 'Strahlungsdämpfung (Schallabstrahlung)', '0002'],
    ['ζ (Dämpfungsmaß)', '—', 'Dämpfungsverhältnis = 1/(2Q)', '0002, 0005'],
    ['Δζ_n', '—', 'Zusätzliche Dämpfung des n-ten OT durch Kammer', '0005'],
    ['τ', 'ms', 'Einschwingzeit = Q/(π·f)', '0002, 0005'],
    ['Δτ', 'ms', 'Verlängerung der Einschwingzeit durch Kammerkopplung', '0005'],
], caption='Tabelle 7: Dämpfung und Güte. Die Gesamtgüte bestimmt die Ansprache.')


# ==========================================================================
section('8. Akustische Impedanz')
# ==========================================================================

sym_table('Impedanz', [
    ['Z_H(f)', 'Pa·s/m³', 'Komplexe akustische Impedanz des Helmholtz-Resonators', '0002–0005'],
    ['Z_Rohr(f)', 'Pa·s/m³', 'Eingangsimpedanz eines Rohres (Orgelpfeife)', '0003'],
    ['Z_akust,Spalt(t)', 'Pa·s/m³', 'Zeitvariable akustische Spaltimpedanz', '0003'],
    ['Z_mech(f)', 'N·s/m', 'Mechanische Impedanz der Zunge', '0003, 0004'],
    ['Z_0', 'Pa·s/m³', 'Referenz-Impedanz = 1/(ω_H · C_a)', '0005'],
    ['R_H(f)', 'Pa·s/m³', 'Realteil von Z_H = Resistanz (→ Dämpfung/Ansprache)', '0005'],
    ['X_H(f)', 'Pa·s/m³', 'Imaginärteil von Z_H = Reaktanz (→ Frequenzverschiebung)', '0005'],
    ['C_a', 'm³/Pa', 'Akustische Compliance der Kammer = V/(ρc²)', '0004, 0005'],
    ['R_Dämpfung', 'N·s/m', 'Mechanischer Dämpfungswiderstand der Zunge', '0003, 0004'],
], caption='Tabelle 8: Akustische und mechanische Impedanzgrößen')


# ==========================================================================
section('9. Kopplungsgrößen (Zwei-Filter-Modell)')
# ==========================================================================

sym_table('Kopplung', [
    ['H_1(f)', '—', 'Übertragungsfunktion Filter 1 (Zunge, Bandpass bei f_1)', '0004, 0005'],
    ['H_2(f)', '—', 'Übertragungsfunktion Filter 2 (Kammer, Bandpass bei f_H)', '0004, 0005'],
    ['H_ges(f)', '—', 'Geschlossene Übertragungsfunktion = H_1/(1 − κ_eff·H_1·H_2)', '0004, 0005'],
    ['κ', '—', 'Akustische Kopplungskonstante (≈ 0,008)', '0004'],
    ['κ_eff', '—', 'Effektive Kopplung (mit Bernoulli-Verstärkung)', '0004, 0005'],
    ['β (Bernoulli)', '—', 'Bernoulli-Verstärkungsfaktor = δp_Kammer / p_Bernoulli (≈ 0,6 %)', '0004, 0005'],
    ['Δk', 'N/m', 'Zusätzliche Federsteifigkeit durch Kammerluft = A²_eff·ρc²/V', '0004'],
    ['G(f)', '—', 'Impedanz-Verstärkungsfaktor nahe f_H (Spitze = Q_H)', '0004, 0005'],
    ['D(f)', '—', 'Nenner der Resonanzfunktion = (1−(f/f_H)²)² + (f/(f_H·Q_H))²', '0005'],
    ['φ, φ_Zunge, φ_Kammer', '°', 'Phasenwinkel (Barkhausen-Bedingung: Σφ = 0)', '0004'],
], caption='Tabelle 9: Kopplungsgrößen des Zwei-Filter-Modells (Dok. 0004–0005)')


# ==========================================================================
section('10. Indizes und Konventionen')
# ==========================================================================

body('Die Dokumente verwenden folgende Indexkonventionen:')

idx_data = [
    ['Index', 'Bedeutung', 'Beispiel'],
    ['_Spalt', 'Am Spalt zwischen Zunge und Platte', 'v_Spalt, S_Spalt'],
    ['_Klappe, _Klap', 'An der Klappenöffnung', 'S_Klappe, ζ_Klappe'],
    ['_eff', 'Effektiver (gemittelter/gewichteter) Wert', 'S_eff, h_eff, κ_eff'],
    ['_ges', 'Gesamtwert (Summe aller Beiträge)', 'ζ_ges, Q_ges, η_ges'],
    ['_max, _min', 'Maximum / Minimum im Schwingzyklus', 'h_max, S_min'],
    ['_tip', 'Am freien Zungenende', 'w_tip, y_tip'],
    ['_tan', 'Tangentialkomponente (parallel zur Platte)', 'v_tan, f_tan'],
    ['_H', 'Helmholtz (Kammerresonanz)', 'f_H, Q_H, Z_H'],
    ['_VK', 'Vorkammer', 'V_VK'],
    ['_mech', 'Mechanisch (Zungeneigenschaft)', 'f_mech, Z_mech, k_mech'],
    ['_akust', 'Akustisch', 'Δl_akust, Z_akust'],
    ['_frei', 'Ohne Kammer (Prüfblock)', 'f_frei'],
    ['_Kammer', 'In der Kammer montiert', 'f_Kammer, p_Kammer'],
    ['_n', 'n-ter Oberton', 'f_n, Δf_n, Δτ_n, Δζ_n'],
    ['_1, _2', 'Mode 1 (Zunge) / Mode 2 (Kammer)', 'f_1, f_2, f_1\', f_2\''],
    ['(t)', 'Zeitabhängig (im Schwingzyklus)', 'S_Spalt(t), Z_akust(t), R(t)'],
    ['(f)', 'Frequenzabhängig', 'Z_H(f), H_1(f), G(f), R_H(f)'],
]
data_wrapped = [[Paragraph(c, sCellB if i==0 else sCell) for c in row] for i, row in enumerate(idx_data)]
tbl = Table(data_wrapped, colWidths=[28*mm, 52*mm, 52*mm], repeatRows=1)
tbl.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),LIGHTBLUE), ('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm), ('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm), ('VALIGN',(0,0),(-1,-1),'TOP'),
] + [('BACKGROUND',(0,i),(-1,i),LIGHTGRAY) for i in range(2,len(idx_data),2)]))
story.append(tbl)
story.append(Paragraph('Tabelle 10: Indexkonventionen', sCaption))

story.append(Spacer(1, 3*mm))

body('<b>Einheitenkonvention:</b> Die Dokumente verwenden SI-Einheiten mit praktischen Anpassungen: '
     'Längen in mm (Geometrie) oder m (Formeln), Drücke in Pa, Frequenzen in Hz, '
     'Volumenströme in ml/s, Volumina in cm³. In den Formeln von Dok. 0004–0005 stehen alle '
     'Größen in SI-Basiseinheiten (m, kg, s, Pa). Frequenzverschiebungen werden in Cent angegeben '
     '(1 Cent = 1200·log₂(f₂/f₁) für den Vergleich mit der musikalischen Wahrnehmungsschwelle '
     'von ca. ±1–2 Cent).')



# ============ BUILD ============
def hf(canvas, doc):
    canvas.saveState(); canvas.setFont(FONT, 8); canvas.setFillColor(MIDGRAY)
    canvas.drawString(25*mm, A4[1]-15*mm, 'Dok. 0006 — Zeichenerklärung')
    canvas.drawRightString(A4[0]-25*mm, A4[1]-15*mm, f'Seite {doc.page}')
    canvas.drawCentredString(A4[0]/2, 12*mm, 'Formelzeichen aller Dokumente (Dok. 0001–0008)')
    canvas.restoreState()

doc.build(story, onFirstPage=hf, onLaterPages=hf)
print(f'PDF created: {outpath}')
