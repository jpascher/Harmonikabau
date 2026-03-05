#!/usr/bin/env python3
"""Dok 0005: Frequenzverschiebung als Indikator der Ansprache"""

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor, white, black
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, HRFlowable
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
import math

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

outpath = '/mnt/user-data/outputs/0005_ansprache_frequenz_kopplung_De.pdf'
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
sSection = ms('sSec',fs=14,fn=FONTB,tc=DARKBLUE,sb=14,sa=6,ld=18)
sSubsec = ms('sSub',fs=11.5,fn=FONTB,tc=HexColor('#2a4a7f'),sb=10,sa=4,ld=15)
sBody = ms('sB',fs=10,ld=14,sa=6)
sKey = ms('sK',fs=9.5,ld=13,sa=4,li=4*mm,ri=4*mm)
sWarn = ms('sW',fs=9.5,ld=13,sa=4,li=4*mm,ri=4*mm)
sFormula = ms('sF',fs=10,al=TA_CENTER,fn=FONTI,sb=4,sa=6,ld=14)
sCaption = ms('sC',fs=8.5,fn=FONTI,al=TA_CENTER,sa=8,ld=11)
sCell = ms('sCell',fs=8.5,ld=11,sa=0,al=TA_LEFT)
sCellB = ms('sCellB',fs=8.5,fn=FONTB,ld=11,sa=0,al=TA_LEFT)

story = []

def section(n, t): story.append(Paragraph(f'Kapitel {n}: {t}', sSection))
def subsection(t): story.append(Paragraph(t, sSubsec))
def body(t): story.append(Paragraph(t, sBody))
def formula(t): story.append(Paragraph(t, sFormula))

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

def make_table(data, cw=None, cap=None):
    if cw is None: cw = [doc.width/len(data[0])]*len(data[0])
    wrapped = []
    for i, row in enumerate(data):
        wrapped.append([Paragraph(str(c), sCellB if i==0 else sCell) for c in row])
    tbl = Table(wrapped, colWidths=cw, repeatRows=1)
    sc = [('BACKGROUND',(0,0),(-1,0),LIGHTBLUE),('GRID',(0,0),(-1,-1),0.5,MIDGRAY),
          ('TOPPADDING',(0,0),(-1,-1),2*mm),('BOTTOMPADDING',(0,0),(-1,-1),2*mm),
          ('LEFTPADDING',(0,0),(-1,-1),2*mm),('RIGHTPADDING',(0,0),(-1,-1),2*mm),
          ('VALIGN',(0,0),(-1,-1),'TOP')]
    for i in range(1,len(data)):
        if i%2==0: sc.append(('BACKGROUND',(0,i),(-1,i),LIGHTGRAY))
    tbl.setStyle(TableStyle(sc))
    story.append(tbl)
    if cap: story.append(Paragraph(cap, sCaption))
    story.append(Spacer(1,2*mm))


# === Physical constants ===
f1 = 50.0; f_H = 274.0; Q_H = 7.0
rho = 1.2; c = 343.0; V = 1.8e-4
A_eff = 4e-6
rho_s = 7800.0; L_r = 0.07; W_r = 0.008; t_r = 0.000354
m_tot = rho_s*L_r*W_r*t_r; m_eff = 0.2427*m_tot
omega1 = 2*math.pi*f1
k_mech = omega1**2 * m_eff
C_a = V/(rho*c**2)


# ============ TITLE ============
story.append(Paragraph('Dok. 0005: Frequenzverschiebung als Indikator der Ansprache', sTitle))
story.append(Paragraph('Beweis des Zusammenhangs: Δf und τ als Real- und Imaginärteil derselben Kopplung', sSubtitle))
story.append(HRFlowable(width='80%', thickness=2, color=ACCENTRED, spaceAfter=4, spaceBefore=4))
story.append(Paragraph(
    'Beweist, dass Frequenzverschiebung und Ansprache-Verschlechterung physikalisch untrennbar verknüpft sind. '
    'Beide sind Projektionen derselben komplexen Kopplungsimpedanz Z_H(f). Die Frequenzverschiebung ist '
    'daher ein messbarer Indikator für die Ansprachequalität — und umgekehrt. '
    'Zahlenwerte: 50-Hz-Basszunge aus Dok. 0002.', sNote))

story.append(Spacer(1, 2*mm))
# Doc reference table
ref_data = [
    ['Dok. 0001', 'Analyse zur Instrumentenakustik und Wahrnehmungspsychologie'],
    ['Dok. 0002', 'Strömungsanalyse Bass-Stimmzunge 50 Hz — v8 (Hauptdokument)'],
    ['Dok. 0003', 'Impedanzvergleich Durchschlagzunge vs. Labialpfeife'],
    ['Dok. 0004', 'Frequenzvariation — Zwei-Filter-Modell'],
    ['<b>Dok. 0005</b>', '<b>Frequenzverschiebung als Indikator der Ansprache (dieses Dokument)</b>'],
    ['Dok. 0006', 'Zeichenerklärung'],
    ['Dok. 0007', 'Diskant-Stimmstock — Kammerfrequenzen'],
    ['Dok. 0008', 'Klangveränderung durch Kammergeometrie'],
    ['Dok. 0500', 'Leitfaden zur ästhetischen Forensik bei Akkordeon-Gehäusen'],
    ['', '<a href="berechnungen_verifikation.py">berechnungen_verifikation.py</a> — Verifikationsskript'],
]
ref_tbl = Table([[Paragraph(r[0], sCell), Paragraph(r[1], sCell)] for r in ref_data],
                colWidths=[25*mm, doc.width-25*mm-4*mm])
ref_tbl.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm),('BACKGROUND',(0,4),(-1,4),LIGHTGREEN)]))
story.append(ref_tbl)
story.append(Spacer(1, 4*mm))


# ============ CHAPTER 1 ============
section(1, 'Die zentrale These')

body('Dok. 0002 (Kap. 10) zeigt: Eine resonante Kammer verschlechtert die Ansprache, weil sie über '
     'Serienkopplung der Zunge Energie entzieht. Dok. 0004 zeigt: Die Kammer verschiebt die '
     'Zungenfrequenz um 1–5 Cent. Beide Effekte werden durch dieselbe physikalische Größe verursacht — '
     'die <b>komplexe Kopplungsimpedanz</b> Z<sub>H</sub>(f).')

body('Die These dieses Dokuments:')

keybox('<b>Frequenzverschiebung und Ansprache-Verschlechterung sind der Imaginär- und Realteil '
       'derselben komplexen Kopplungsfunktion.</b> Sie sind durch die Kramers-Kronig-Relation '
       'kausal verknüpft. An der Frequenzverschiebung lässt sich ablesen, ob die Ansprache besser '
       'oder schlechter wird — und wie weit man von der kritischen Resonanz entfernt ist.')


# ============ CHAPTER 2 ============
section(2, 'Die komplexe Kopplungsimpedanz')

subsection('Zerlegung in Real- und Imaginärteil')

body('Die Helmholtz-Impedanz der Kammer bei Frequenz f (Dok. 0002 Kap. 11, Dok. 0004 Kap. 3):')
formula('Z<sub>H</sub>(f) = R<sub>H</sub>(f) + j·X<sub>H</sub>(f)')

body('mit dem <b>Realteil</b> (Resistanz = Energiedissipation):')
formula('R<sub>H</sub>(f) = Re(Z<sub>H</sub>) = Z<sub>0</sub> · (f/f<sub>H</sub>) / (Q<sub>H</sub> · D(f))')

body('und dem <b>Imaginärteil</b> (Reaktanz = Energiespeicherung):')
formula('X<sub>H</sub>(f) = Im(Z<sub>H</sub>) = Z<sub>0</sub> · (1 − (f/f<sub>H</sub>)²) / D(f)')

body('wobei D(f) = (1 − (f/f<sub>H</sub>)²)² + (f/(f<sub>H</sub> · Q<sub>H</sub>))² und '
     'Z<sub>0</sub> = 1/(ω<sub>H</sub> · C<sub>a</sub>).')

subsection('Was Real- und Imaginärteil physikalisch bedeuten')

body('<b>Realteil R<sub>H</sub>(f) → Energieverlust → Ansprache-Verschlechterung:</b> '
     'Jede Druckschwankung, die gegen den Realteil der Kammerimpedanz arbeitet, wird in Wärme '
     'umgewandelt (Reibung an Kammerwänden, Verwirbelung). Diese Energie fehlt der Zunge. '
     'Die zusätzliche Dämpfung durch die Kammer-Kopplung ist:')
formula('Δζ<sub>n</sub> = κ<sub>eff</sub>² · R<sub>H</sub>(nf<sub>1</sub>) / Z<sub>0</sub>')

body('Die Einschwingzeit verlängert sich um:')
formula('Δτ<sub>n</sub> = Δζ<sub>n</sub> / (π · f<sub>1</sub>)')

body('<b>Imaginärteil X<sub>H</sub>(f) → Frequenzverschiebung:</b> '
     'Der Imaginärteil wirkt als reaktive Last — kapazitiv (Feder) unterhalb von f<sub>H</sub>, '
     'induktiv (Masse) oberhalb. Die Frequenzverschiebung des n-ten Obertons ist (Dok. 0004):')
formula('Δf<sub>n</sub> = −κ<sub>eff</sub>² · f<sub>n</sub> · X<sub>H</sub>(f<sub>n</sub>) / (2 · Z<sub>0</sub>)')

keybox('<b>Schlüsselaussage:</b> Δτ (Ansprache) ∝ R<sub>H</sub>(f) und Δf (Frequenzverschiebung) ∝ X<sub>H</sub>(f). '
       'Beide stammen aus <b>derselben</b> komplexen Impedanz Z<sub>H</sub> = R<sub>H</sub> + jX<sub>H</sub>. '
       'Sie verhalten sich zueinander wie der Sinus und Kosinus eines Winkels — kennt man den einen, '
       'kennt man den anderen.')


# ============ CHAPTER 3 ============
section(3, 'Quantitativer Beweis: Die Zahlentabelle')

body('Für jeden Oberton n der 50-Hz-Zunge berechnen wir R<sub>H</sub>, X<sub>H</sub>, die daraus '
     'resultierende Ansprache-Verschlechterung Δτ und die Frequenzverschiebung Δf:')

# Calculate for each overtone
data = [['OT n', 'f [Hz]', 'R_H/Z₀ (Dämpf.)', 'X_H/Z₀ (Reaktanz)', 'Δτ [ms]', 'Δf_n [Cent]', 'Diagnose']]

kappa_eff = 0.008  # from Dok 0004
beta_amp = 10.0  # Bernoulli amplification factor (effective)

for n in [1, 2, 3, 4, 5, 6, 7, 8, 10, 15, 20]:
    f = n * f1
    r = f / f_H
    D = (1 - r**2)**2 + (r/Q_H)**2
    R_norm = (r/Q_H) / D  # R_H / Z_0
    X_norm = (1 - r**2) / D  # X_H / Z_0 (note: this is actually -X/Z_0 for the reactive part)
    
    # Damping increment
    d_zeta = kappa_eff**2 * R_norm * beta_amp
    d_tau = d_zeta / (math.pi * f1) * 1000  # ms
    
    # Frequency shift (imaginary part causes shift)
    # Below f_H: X > 0 (capacitive, stiffness), f shifts up
    # Above f_H: X < 0 (inductive, mass), f shifts down  
    df_cents = -kappa_eff**2 * beta_amp * X_norm * 1200 / (2 * math.log(2))  # rough cents scaling
    # More precise: use phase-based model
    phase = math.degrees(math.atan(Q_H * (r - 1/r)))
    
    # Effective frequency shift via phase pulling
    # df/f = -tan(phase_chamber) / (2*Q_reed) * coupling
    Q_reed = 75  # mid-range estimate
    df_rel = -math.tan(math.radians(phase)) * kappa_eff**2 * beta_amp / 2
    df_cents_v2 = 1200 * math.log2(1 + abs(df_rel)) * (1 if df_rel > 0 else -1) if abs(df_rel) < 0.5 else 0
    
    # Diagnose
    if R_norm > 2.0:
        diag = '<b>KRITISCH: max. Dämpfung</b>'
    elif R_norm > 0.5:
        diag = 'Merklich — nahe Resonanz'
    elif abs(X_norm) > 1.0 and R_norm < 0.3:
        diag = 'Δf groß, Dämpfung klein'
    else:
        diag = 'Entkoppelt — gut'
    
    data.append([
        str(n), f'{f:.0f}',
        f'{R_norm:.3f}', f'{X_norm:+.3f}',
        f'{d_tau:.1f}', f'{df_cents_v2:+.1f}' if abs(df_cents_v2) > 0.01 else '≈ 0',
        diag
    ])

make_table(data, cw=[14*mm, 16*mm, 27*mm, 27*mm, 16*mm, 18*mm, 32*mm],
    cap='Tabelle 1: Dämpfung (Realteil) und Frequenzverschiebung (Imaginärteil) der Kammerimpedanz '
        'bei jedem Oberton der 50-Hz-Zunge. f_H = 274 Hz, Q_H = 7, κ_eff = 0,008 × β.')

# ============ KEY PATTERN ============
subsection('Das Muster: Drei Zonen')

body('Die Tabelle zeigt ein klares Muster, das sich in drei Zonen gliedert:')

body('<b>Zone 1 — Weit unter f<sub>H</sub> (OT 1–3, f &lt; 0,55·f<sub>H</sub>):</b> '
     'R<sub>H</sub> ≈ 0, X<sub>H</sub> ≈ +1 (reine Compliance). '
     'Minimale Dämpfung, kleine positive Frequenzverschiebung. '
     '<b>Gute Ansprache</b>, Tonhöhe leicht erhöht. Die Kammer ist „unsichtbar".')

body('<b>Zone 2 — Nahe f<sub>H</sub> (OT 5–6, f ≈ 0,9–1,1·f<sub>H</sub>):</b> '
     'R<sub>H</sub> maximal (Spitze bei f = f<sub>H</sub>), X<sub>H</sub> wechselt das Vorzeichen. '
     'Maximale Dämpfung, Frequenzverschiebung ändert Richtung. '
     '<b>Schlechteste Ansprache</b>, aber die Frequenzverschiebung des Einzelobertons '
     'geht durch Null — genau am Resonanzpunkt.')

body('<b>Zone 3 — Weit über f<sub>H</sub> (OT 8+, f &gt; 1,5·f<sub>H</sub>):</b> '
     'R<sub>H</sub> ≈ 0, X<sub>H</sub> ≈ −1 (reine Inertanz). '
     'Minimale Dämpfung, kleine negative Frequenzverschiebung. '
     '<b>Gute Ansprache</b>, Tonhöhe leicht erniedrigt.')

warnbox('<b>Das Paradoxon:</b> Genau bei optimaler Resonanz (OT exakt auf f<sub>H</sub>) '
        'ist die Dämpfung maximal — aber die Frequenzverschiebung <b>dieses</b> Obertons ist Null! '
        'Das scheint zu widersprechen. Die Auflösung: Am Resonanzpunkt wird alle Kopplungsenergie '
        'in Dissipation umgewandelt (Realteil), nichts bleibt für Reaktanz (Imaginärteil) übrig. '
        'Der Oberton wird nicht frequenzverschoben, aber <b>gedämpft</b> — das ist schlimmer.')


# ============ CHAPTER 4 ============
section(4, 'Der Vorzeichenwechsel als Warnsignal')

body('Die entscheidende praktische Erkenntnis liegt nicht im Absolutwert der Frequenzverschiebung, '
     'sondern im <b>Vorzeichenwechsel</b>:')

infobox(
    '<b>Das Vorzeichen-Kriterium:</b><br/><br/>'
    '• Δf &gt; 0 (Ton wird <b>höher</b>): Der wirksame Oberton liegt <b>unter</b> f<sub>H</sub>. '
    'Die Kammer wirkt als Feder (Compliance). Ansprache ist gut, aber verschlechtert sich, '
    'wenn Δf weiter wächst.<br/><br/>'
    '• Δf ≈ 0 mit schlechter Ansprache: Ein Oberton sitzt <b>exakt auf</b> f<sub>H</sub>. '
    '<b>Worst Case.</b> Maximale Energieextraktion, aber kein Frequenzsignal!<br/><br/>'
    '• Δf &lt; 0 (Ton wird <b>tiefer</b>): Der wirksame Oberton liegt <b>über</b> f<sub>H</sub>. '
    'Die Kammer wirkt als Masse (Inertanz). Ansprache war schlecht und wird wieder besser.<br/><br/>'
    '• Δf ≈ 0 mit guter Ansprache: Kein Oberton in der Nähe von f<sub>H</sub>. <b>Best Case.</b>',
    bg=LIGHTBLUE, fr=DARKBLUE
)

body('Der <b>Vorzeichenwechsel</b> von Δf ist das Warnsignal: Wenn sich die Frequenzverschiebung '
     'beim Variieren der Kammergeometrie von positiv nach negativ ändert (oder umgekehrt), '
     'hat man die Resonanz <b>durchfahren</b>. Am Nulldurchgang ist die Ansprache am schlechtesten.')

subsection('Nyquist-Diagramm der Kopplung')

body('In der Signaltheorie stellt man komplexe Übertragungsfunktionen als Nyquist-Diagramm dar '
     '(Realteil vs. Imaginärteil, parametrisch über die Frequenz). Für die Kammer-Zungen-Kopplung:')

body('<b>Horizontale Achse:</b> R<sub>H</sub>(f) = Ansprache-Verschlechterung (Dämpfung)')
body('<b>Vertikale Achse:</b> X<sub>H</sub>(f) = Frequenzverschiebung (Reaktanz)')

body('Die Kurve ist ein <b>Kreis</b> (für einen einfachen Helmholtz-Resonator). '
     'Jeder Oberton der Zunge liegt auf einem bestimmten Punkt dieses Kreises:')

# Calculate Nyquist circle points
data_nyq = [['OT n', 'f [Hz]', 'R_H (→ Δτ)', 'X_H (→ Δf)', 'Position auf dem Kreis']]
for n in [1, 3, 5, 6, 8, 10]:
    f = n * f1
    r = f / f_H
    D = (1 - r**2)**2 + (r/Q_H)**2
    R_n = (r/Q_H) / D
    X_n = (1 - r**2) / D
    
    if R_n > 1.5:
        pos = '⬤ Spitze des Kreises (Resonanz)'
    elif X_n > 0.3:
        pos = '↑ Obere Hälfte (kapazitiv, Δf > 0)'
    elif X_n < -0.3:
        pos = '↓ Untere Hälfte (induktiv, Δf < 0)'
    else:
        pos = '○ Nahe Ursprung (entkoppelt)'
    
    data_nyq.append([str(n), f'{f:.0f}', f'{R_n:.3f}', f'{X_n:+.3f}', pos])

make_table(data_nyq, cw=[14*mm, 18*mm, 25*mm, 25*mm, 55*mm],
    cap='Tabelle 2: Nyquist-Positionen der Obertöne. Am Kreisgipfel (OT 5–6): '
        'maximale Dämpfung, Frequenzverschiebung wechselt Vorzeichen.')

body('Der 5. Oberton (250 Hz) liegt in der oberen Kreishälfte (X &gt; 0, Δf &gt; 0), '
     'der 6. Oberton (300 Hz) in der unteren (X &lt; 0, Δf &lt; 0). Dazwischen, bei f<sub>H</sub> = 274 Hz, '
     'liegt der Kreisgipfel — maximales R, X = 0. Die Tatsache, dass kein Oberton exakt auf den Gipfel '
     'trifft, ist der Grund für die relativ gute Ansprache der 50-Hz-Zunge in dieser Kammer.')


# ============ CHAPTER 5 ============
section(5, 'Der mathematische Beweis: Kramers-Kronig')

body('Die Verknüpfung von R<sub>H</sub> und X<sub>H</sub> ist nicht zufällig, sondern '
     '<b>physikalisch zwingend</b>. Jede kausale, lineare Antwortfunktion erfüllt die '
     'Kramers-Kronig-Relationen:')

formula('X<sub>H</sub>(ω) = −(1/π) · P ∫ R<sub>H</sub>(ω\') / (ω\' − ω) dω\'')
formula('R<sub>H</sub>(ω) = +(1/π) · P ∫ X<sub>H</sub>(ω\') / (ω\' − ω) dω\'')

body('(P = Cauchyscher Hauptwert). Das bedeutet: <b>Wenn man den Frequenzverlauf der '
     'Frequenzverschiebung kennt (alle X<sub>H</sub>(f)), kann man daraus die '
     'Ansprache-Verschlechterung berechnen (alle R<sub>H</sub>(f)) — und umgekehrt.</b>')

body('Für den Helmholtz-Resonator ist die Kramers-Kronig-Relation explizit verifizierbar:')

formula('Z<sub>H</sub>(f) = Z<sub>0</sub> / (1 − (f/f<sub>H</sub>)² + j·f/(f<sub>H</sub>·Q<sub>H</sub>))')

body('Der Realteil (Lorentz-Linie, Glockenform mit Spitze bei f<sub>H</sub>) und der Imaginärteil '
     '(Dispersionskurve, S-Form mit Nulldurchgang bei f<sub>H</sub>) sind durch die '
     'Hilbert-Transformation verknüpft — genau wie Absorption und Brechungsindex in der Optik.')

keybox('<b>Physikalische Bedeutung:</b> Die Kramers-Kronig-Relation ist eine Folge der <b>Kausalität</b>. '
       'Die Kammer kann nicht „wählen", nur Frequenz zu verschieben ohne Energie zu dissipieren '
       '(oder umgekehrt). Beide Effekte sind zwei Seiten derselben Medaille. Ein System, das die '
       'Frequenz stark zieht, <b>muss</b> auch Energie dissipieren — und ein System, das Energie '
       'dissipiert, <b>muss</b> auch die Frequenz verschieben (außer exakt bei Resonanz, wo X = 0).')


# ============ CHAPTER 6 ============
section(6, 'Die praktische Messvorschrift')

body('Aus dem Beweis folgt eine praktische Methode, die Ansprache <b>ohne Spieltest</b> zu bewerten:')

subsection('Methode: Frequenzmessung mit und ohne Kammer')

body('<b>Schritt 1:</b> Zungenfrequenz <b>frei</b> messen (Stimmplatte auf Prüfblock, keine Kammer). '
     'Ergebnis: f<sub>frei</sub>.')

body('<b>Schritt 2:</b> Zungenfrequenz <b>in der Kammer</b> messen (gleicher Balgdruck). '
     'Ergebnis: f<sub>Kammer</sub>.')

body('<b>Schritt 3:</b> Δf = f<sub>Kammer</sub> − f<sub>frei</sub> berechnen.')

body('<b>Schritt 4:</b> Interpretation:')

data_interp = [
    ['Messergebnis', 'Physikalische Bedeutung', 'Ansprache-Prognose'],
    ['Δf ≈ 0, Ansprache gut', 'Keine Obertöne nahe f_H. Kammer entkoppelt.',
     '<b>Optimal</b> — nichts zu ändern'],
    ['Δf > 0 (Ton höher)', 'Nächster OT <b>unter</b> f_H. Kammer wirkt als Feder.',
     'Gut, aber wird schlechter wenn Δf wächst'],
    ['Δf maximal positiv', 'OT nähert sich f_H von unten. Starke Kopplung.',
     '<b>Warnung</b> — Resonanz nahe'],
    ['Δf → 0 (von + nach −)', '<b>OT durchfährt f_H.</b> Vorzeichenwechsel!',
     '<b>Worst Case</b> — maximale Dämpfung'],
    ['Δf < 0 (Ton tiefer)', 'Nächster OT <b>über</b> f_H. Kammer wirkt als Masse.',
     'Besserung — Resonanz überschritten'],
    ['Δf ≈ 0, Ansprache gut', 'OT weit von f_H entfernt.',
     '<b>Optimal</b> — entkoppelt'],
]
make_table(data_interp, cw=[30*mm, 50*mm, 50*mm],
    cap='Tabelle 3: Interpretationsschlüssel — Frequenzverschiebung → Ansprache-Prognose')

subsection('Methode: Kammervolumen variieren')

body('Elegant ist die <b>kontinuierliche Variation</b>: Eine verstellbare Kammer (z.B. mit verschiebbarem '
     'Boden) ermöglicht es, das Volumen V und damit f<sub>H</sub> = (c/2π)·√(S<sub>Hals</sub>/(V·l<sub>Hals</sub>)) '
     'stetig zu verändern. Man misst gleichzeitig Frequenz und Ansprache:')

body('Wenn man f<sub>H</sub> von oben nach unten durchstimmt (Volumen vergrößern), '
     'durchfährt f<sub>H</sub> nacheinander die Obertöne 8, 7, 6, 5, 4, 3, 2. Bei jedem '
     'Durchgang beobachtet man: (a) Δf wechselt das Vorzeichen — (b) die Ansprache wird kurzzeitig '
     'schlechter — (c) beides gleichzeitig. Die <b>Synchronizität</b> von (a) und (b) ist der '
     'experimentelle Beweis der Kramers-Kronig-Verknüpfung.')


# ============ CHAPTER 7 ============
section(7, 'Quantitatives Vorhersagemodell')

body('Aus dem Zwei-Filter-Modell (Dok. 0004) und der Kramers-Kronig-Zerlegung ergibt sich ein '
     'geschlossenes Vorhersagemodell. Für die 50-Hz-Zunge bei verschiedenen f<sub>H</sub>-Werten:')

# Big calculation table
data_pred = [['f_H [Hz]', 'Nächster OT', 'R_H/Z₀', 'X_H/Z₀', 'Δτ [ms]', 'Δf [Cent]', 'Ansprache']]

for fH_test in [150, 175, 200, 225, 250, 274, 300, 325, 350, 400, 500]:
    # Find most coupled overtone (closest to fH)
    n_best = round(fH_test / f1)
    if n_best < 1: n_best = 1
    f_ot = n_best * f1
    r = f_ot / fH_test
    D = (1 - r**2)**2 + (r/Q_H)**2
    R_n = (r/Q_H) / D
    X_n = (1 - r**2) / D
    
    d_zeta = kappa_eff**2 * R_n * beta_amp
    d_tau = d_zeta / (math.pi * f1) * 1000
    
    # Frequency shift sign depends on X
    phase = math.degrees(math.atan(Q_H * (r - 1/r)))
    Q_reed = 75
    df_rel = -math.tan(math.radians(phase)) * kappa_eff**2 * beta_amp / 2
    if abs(df_rel) < 0.5:
        df_c = 1200 * math.log2(1 + abs(df_rel)) * (1 if df_rel > 0 else -1)
    else:
        df_c = 600 * df_rel  # linear approx for large values
    
    if R_n > 3:
        ansp = '⬤ Schlecht'
    elif R_n > 1:
        ansp = '▲ Merklich verschlechtert'
    elif R_n > 0.3:
        ansp = '△ Leicht verschlechtert'
    else:
        ansp = '○ Gut'
    
    data_pred.append([
        f'{fH_test}', f'{n_best}. OT ({f_ot:.0f} Hz)',
        f'{R_n:.2f}', f'{X_n:+.2f}',
        f'{d_tau:.1f}', f'{df_c:+.1f}' if abs(df_c) > 0.05 else '≈ 0',
        ansp
    ])

make_table(data_pred, cw=[18*mm, 30*mm, 19*mm, 19*mm, 16*mm, 18*mm, 30*mm],
    cap='Tabelle 4: Vorhersage von Ansprache und Frequenzverschiebung bei verschiedenen f_H-Werten. '
        'Rot: f_H trifft Oberton exakt (R_H maximal, X_H ≈ 0).')

body('Das Muster ist eindeutig:')

body('<b>f<sub>H</sub> = 200, 250, 300, 350, 400 Hz</b> (exakt auf Obertönen): '
     'R<sub>H</sub> maximal (= Q<sub>H</sub>), X<sub>H</sub> ≈ 0, '
     'Δτ maximal, Δf ≈ 0. <b>Schlechteste Ansprache bei geringstem Frequenzsignal.</b>')

body('<b>f<sub>H</sub> = 225, 274, 325 Hz</b> (zwischen Obertönen): '
     'R<sub>H</sub> moderat, X<sub>H</sub> ≠ 0, '
     'Δτ mäßig, Δf merklich. <b>Bessere Ansprache, aber Frequenzverschiebung als „Warnsignal".</b>')

warnbox('<b>Das zentrale Ergebnis:</b> Die Frequenzverschiebung ist ein „Kanarienvogel im Bergwerk". '
        'Solange Δf groß ist, ist man <b>nicht</b> auf der Resonanz — die Ansprache ist mäßig '
        'belastet, aber nicht kritisch. Wenn Δf gegen Null geht <b>und gleichzeitig</b> das Vorzeichen '
        'wechselt, hat man die Resonanz durchfahren — die Ansprache ist am schlechtesten. '
        'Wenn Δf danach wieder anwächst (mit umgekehrtem Vorzeichen), bessert sich die Ansprache.')


# ============ CHAPTER 8 ============
section(8, 'Die Gesamtbilanz: Summe über alle Obertöne')

body('Die Tabellen 1 und 4 zeigen den Effekt für einzelne Obertöne. Der <b>Gesamteffekt</b> auf den '
     'Grundton ist die Summe über alle Obertöne, gewichtet mit ihrer Amplitude a<sub>n</sub>:')

formula('Δf<sub>gesamt</sub> = Σ<sub>n</sub> (a<sub>n</sub>/n) · Δf<sub>n</sub>')
formula('Δτ<sub>gesamt</sub> = Σ<sub>n</sub> a<sub>n</sub>² · Δτ<sub>n</sub>')

body('Da die Obertonamplituden a<sub>n</sub> mit n abfallen (typisch a<sub>n</sub> ~ 1/n für eine '
     'Durchschlagzunge), dominieren die niedrigen Obertöne. Aber die Obertöne nahe f<sub>H</sub> '
     'werden durch G(f) verstärkt — es gibt also einen <b>Wettbewerb</b> zwischen der abfallenden '
     'Amplitude und der steigenden Kammer-Kopplung.')

body('Für die 50-Hz-Zunge mit f<sub>H</sub> = 274 Hz: Die Obertöne 5 (250 Hz) und 6 (300 Hz) '
     'dominieren die Kopplung. Ihre Beiträge haben <b>entgegengesetztes Vorzeichen</b> '
     '(OT 5: Δf &gt; 0, OT 6: Δf &lt; 0), heben sich teilweise auf. Das Nettoergebnis ist positiv '
     '(Δf &gt; 0), weil OT 5 näher an f<sub>H</sub> liegt als OT 6 (250 Hz vs. 300 Hz zu 274 Hz: '
     'Abstand 24 vs. 26 Hz).')

keybox('<b>Praktische Faustregel:</b> Das Vorzeichen der Netto-Frequenzverschiebung zeigt an, '
       'von welcher Seite der nächste Oberton an f<sub>H</sub> herankommt:<br/><br/>'
       '• Δf<sub>netto</sub> &gt; 0 → dominanter OT liegt <b>unter</b> f<sub>H</sub> → f<sub>H</sub> <b>senken</b> verschlechtert Ansprache<br/>'
       '• Δf<sub>netto</sub> &lt; 0 → dominanter OT liegt <b>über</b> f<sub>H</sub> → f<sub>H</sub> <b>heben</b> verschlechtert Ansprache<br/>'
       '• Δf<sub>netto</sub> ≈ 0 → entweder entkoppelt (gute Ansprache) oder auf Resonanz (schlechte Ansprache) → Ansprache-Test entscheidet')


# ============ CHAPTER 9 — SETUP ============
section(9, 'Das Problem der fehlenden Referenzfrequenz')

body('Die Messvorschrift (f<sub>frei</sub> vs. f<sub>Kammer</sub>) hat eine fundamentale '
     '<b>praktische Einschränkung</b>:')

warnbox('<b>Es gibt keine absolute Referenzfrequenz.</b> Die Frequenz f<sub>frei</sub> '
        '(Zunge auf dem Prüfblock, ohne Kammer) ist <b>selbst</b> eine setup-abhängige Größe. '
        'Der Prüfblock hat ein eigenes Volumen, eine eigene akustische Impedanz. '
        'f<sub>frei</sub> ist nicht die „wahre" Eigenfrequenz der Zunge — es ist die '
        'Frequenz der Zunge <b>in einem anderen akustischen Umfeld</b>.')

body('Die physikalische Eigenfrequenz im Vakuum (f<sub>vak</sub>) ist nicht praktisch messbar — '
     'die Zunge braucht Luft zum Schwingen. Selbst eine Zunge im offenen Raum hat eine akustische '
     'Umgebung: die freie Schallabstrahlung als schwache Belastung.')

subsection('Die drei „Referenzfrequenzen" des Akkordeonbauers')

body('<b>f<sub>Prüfblock</sub>:</b> Zunge auf standardisiertem Holzblock, definierter Druck. '
     'Am nächsten an f<sub>vak</sub>, aber nicht identisch (der Block hat eine Compliance). '
     'Geschätzter Unterschied: 0,1–0,5 Cent.')

body('<b>f<sub>Kammer,offen</sub>:</b> Zunge in der Kammer montiert, Klappe offen. '
     'Die offene Klappe verändert die Kammerimpedanz gegenüber dem Spielbetrieb.')

body('<b>f<sub>Kammer,Spiel</sub>:</b> Zunge in der Kammer, Klappe geschlossen, Balg-angeblasen. '
     'Der „echte" Betriebszustand — aber hier ist die Frequenz bereits durch die volle '
     'Kammerkopplung verändert.')

subsection('Das Vorstimm-Setup als lokale Referenz')

body('Da keine der drei Referenzfrequenzen „absolut" ist, arbeitet der erfahrene Instrumentenbauer '
     'mit einer <b>lokalen Referenz</b>: Er definiert sein persönliches Vorstimm-Setup '
     '(bestimmter Prüfblock, bestimmter Blasdruck, bestimmte Raumtemperatur) als <b>Nullpunkt</b>. '
     'Das Setup ist die „Eichung" des Stimmvorgangs.')

body('Die Δf-Methode funktioniert trotzdem, weil sie eine <b>Differenz</b> misst: '
     'Δf = f<sub>Kammer</sub> − f<sub>Setup</sub>. Die absolute Frequenz spielt keine Rolle — '
     'nur die Differenz und vor allem der <b>Vorzeichenwechsel</b>. '
     'Der Vorzeichenwechsel ist setup-unabhängig, weil er eine Eigenschaft der '
     'Kammer-Zungen-Kopplung ist, nicht der absoluten Frequenz.')

subsection('Was die Erfahrung liefert, was die Berechnung nicht kann')

body('<b>1. Die Offset-Korrektur:</b> Wie viel Cent liegt f<sub>Setup</sub> über f<sub>vak</sub>? '
     'Das weiß nur, wer jahrelang mit demselben Prüfblock gestimmt hat und die resultierenden '
     'Instrumente spielen gehört hat.')

body('<b>2. Die Temperaturkompensation:</b> f<sub>Zunge</sub> ändert sich mit der Temperatur '
     '(E-Modul: ~0,5 Cent/°C). c ändert sich ebenfalls (+0,17 %/°C → verändert f<sub>H</sub> '
     'und Δf). Der erfahrene Stimmer kompensiert beides intuitiv.')

body('<b>3. Die Druckabhängigkeit:</b> Der Balgdruck verändert h<sub>eff</sub> '
     '(Dok. 0002 Kap. 9) und damit S<sub>Spalt</sub>, κ, und Δf. '
     'Das Vorstimm-Setup definiert einen Arbeitspunkt auf der Druckkurve.')

body('<b>4. Die Materialstreuung:</b> 0,1 mm Unterschied in der Aufbiegung → 7 % in S<sub>eff</sub> '
     '(Dok. 0002 Kap. 13). Die Berechnung gibt den Mittelwert — die Erfahrung sagt, '
     'wie weit der Einzelfall abweichen kann.')

keybox('<b>Die Rolle der Berechnung im Stimmvorgang:</b> Die Formeln aus Dok. 0002–0005 '
       'ersetzen nicht die Erfahrung. Sie erklären, <b>warum</b> die Erfahrungsregeln '
       'funktionieren, sagen <b>in welche Richtung</b> eine Änderung wirkt, identifizieren '
       'die <b>empfindlichen Parameter</b>, und warnen, <b>wann</b> man aufpassen muss '
       '(Vorzeichenwechsel → Resonanz). Was die Berechnung <b>nicht</b> kann: den exakten '
       'Cent-Wert vorhersagen, die Materialstreuung kompensieren, das Setup kennen.')


# ============ CHAPTER 10 — KLARSTELLUNG ============
section(10, 'Klarstellung: Ansprache vs. Lautstärke — kein Trade-off bei der Kammerabstimmung')

body('Ein naheliegendes Missverständnis muss ausgeräumt werden: Wenn optimale Ansprache bedeutet, '
     'dass die Kammer der Zunge möglichst <b>wenig</b> Energie entzieht — reduziert dann gute '
     'Ansprache die maximal erzielbare Lautstärke?')

body('<b>Nein.</b> Die Antwort ergibt sich direkt aus der Physik der Serienkopplung '
     '(Dok. 0002 Kap. 10). Zwei Effekte müssen sauber getrennt werden:')

subsection('Kammerabstimmung: Win-Win')

body('Wenn ein Oberton exakt auf f<sub>H</sub> liegt, ist R<sub>H</sub> maximal — die Kammer '
     'entzieht der Zunge Energie über Serienkopplung. Dieser Energieverlust verschlechtert '
     '<b>beides gleichzeitig</b>:')

body('<b>Ansprache:</b> Die entzogene Energie fehlt beim Amplitudenaufbau. '
     'Die Einschwingzeit τ verlängert sich um Δτ ∝ R<sub>H</sub>.')

body('<b>Lautstärke:</b> Im eingeschwungenen Zustand fließt pro Zyklus Energie in die '
     'Kammer-Dissipation ab, statt in Schall umgewandelt zu werden. '
     'Die Gleichgewichtsamplitude sinkt: A<sub>max</sub> ∝ 1/√R<sub>H</sub>.')

body('Beide Effekte stammen aus derselben Quelle — R<sub>H</sub>(f). Wenn R<sub>H</sub> sinkt '
     '(kein Oberton auf f<sub>H</sub>), verbessern sich <b>beide</b>: schnellere Ansprache '
     '<b>und</b> größere Amplitude. Bei der Kammerabstimmung gibt es keinen Trade-off.')

subsection('Zungen-Güte Q₁: Der echte Trade-off')

body('Der Trade-off zwischen Ansprache und Lautstärke existiert — aber er liegt <b>nicht</b> '
     'in der Kammerabstimmung, sondern in der Zunge selbst:')

body('<b>Hohes Q<sub>1</sub></b> (geringe Dämpfung): Viel gespeicherte Energie → <b>lauter</b>. '
     'Aber viele Zyklen zum Einschwingen → <b>langsame Ansprache</b> (τ = Q<sub>1</sub>/(πf)).')

body('<b>Niedriges Q<sub>1</sub></b> (starke Dämpfung): Wenige Zyklen → <b>schnelle Ansprache</b>. '
     'Aber geringere Gleichgewichtsamplitude → <b>leiser</b>.')

body('Dieser Trade-off wird durch Material, Einspannung und Aufbiegung bestimmt '
     '(Dok. 0002 Kap. 9) — <b>nicht</b> durch die Kammergeometrie.')

keybox('<b>Die Unterscheidung:</b><br/><br/>'
       '<b>Kammerabstimmung</b> (f<sub>H</sub> relativ zu Obertönen): Beeinflusst die '
       '<b>externe Dämpfung</b> durch die Kammer. Optimierung = R<sub>H</sub> minimieren '
       '= Win-Win (bessere Ansprache <i>und</i> mehr Lautstärke).<br/><br/>'
       '<b>Zungengüte Q<sub>1</sub></b> (Material, Einspannung): Beeinflusst die '
       '<b>interne Dämpfung</b> der Zunge. Hier gibt es den echten Trade-off: '
       'mehr Q<sub>1</sub> = lauter aber langsamer, weniger Q<sub>1</sub> = schneller aber leiser.<br/><br/>'
       'Der Instrumentenbauer optimiert beides <b>unabhängig</b>: Die Kammer so abstimmen, '
       'dass sie möglichst wenig stört (f<sub>H</sub> zwischen Obertönen). Die Zunge so wählen, '
       'dass Q<sub>1</sub> zum musikalischen Kontext passt (Bass: eher niedrig für schnelle '
       'Ansprache; Diskant: eher hoch für Lautstärke und Sustain).')


# ============ CHAPTER 11 — ZUSAMMENFASSUNG ============
section(11, 'Zusammenfassung')

keybox(
    '<b>Bewiesener Zusammenhang:</b> Frequenzverschiebung Δf und Ansprache-Verschlechterung Δτ '
    'sind der Imaginär- und Realteil derselben komplexen Kopplungsimpedanz Z<sub>H</sub>(f). '
    'Sie sind durch die Kramers-Kronig-Relation physikalisch zwingend verknüpft.<br/><br/>'
    
    '<b>Optimale Resonanz mit Obertönen ergibt die schlechteste Ansprache:</b> '
    'Wenn f<sub>H</sub> = n · f<sub>Zunge</sub> (ganzzahliges Vielfaches), ist R<sub>H</sub> maximal '
    '(= Q<sub>H</sub>), die Energieextraktion maximal, die Ansprache am schlechtesten. '
    'Gleichzeitig ist X<sub>H</sub> = 0, die Frequenzverschiebung dieses Obertons verschwindet.<br/><br/>'
    
    '<b>Die Frequenzverschiebung ist ein messbarer Indikator:</b><br/>'
    '• Großes |Δf| → starke Kopplung, aber nicht kritisch (Oberton neben f<sub>H</sub>)<br/>'
    '• Δf → 0 mit Vorzeichenwechsel → Oberton auf f<sub>H</sub> → <b>schlechteste Ansprache</b><br/>'
    '• Δf → 0 ohne Vorzeichenwechsel → Entkopplung → <b>beste Ansprache</b><br/><br/>'
    
    '<b>Messvorschrift:</b> Zungenfrequenz frei vs. in Kammer messen. Das Vorzeichen von Δf zeigt, '
    'von welcher Seite ein Oberton an f<sub>H</sub> herankommt. Der Vorzeichenwechsel beim Variieren '
    'der Kammergeometrie markiert den Resonanzdurchgang — dort ist die Ansprache am schlechtesten.<br/><br/>'
    
    '<b>Fehlende absolute Referenz:</b> f<sub>frei</sub> ist setup-abhängig. '
    'Das Vorstimm-Setup (Prüfblock, Blasdruck, Temperatur) definiert den lokalen Nullpunkt. '
    'Der Vorzeichenwechsel von Δf ist setup-unabhängig — er ist die physikalische Observable.<br/><br/>'
    
    '<b>Kein Trade-off bei der Kammerabstimmung:</b> Optimale Abstimmung (f<sub>H</sub> zwischen '
    'Obertönen) verbessert Ansprache <i>und</i> Lautstärke gleichzeitig, weil beide unter '
    'derselben Quelle (R<sub>H</sub>) leiden. Der Trade-off zwischen Ansprache und Lautstärke '
    'liegt in der Zungengüte Q<sub>1</sub> (Material, Einspannung) — nicht in der Kammer.<br/><br/>'
    
    '<b>Die Berechnung sagt, wo man suchen soll. Das Setup definiert den Nullpunkt. '
    'Die Erfahrung sagt, was man findet.</b>'
)

story.append(Spacer(1, 6*mm))
story.append(HRFlowable(width='60%', thickness=1, color=ACCENTRED, spaceAfter=4, spaceBefore=2))
story.append(Paragraph(
    '<i>Was man nicht messen kann, kann man nicht optimieren. Die Frequenzverschiebung ist die Messgröße — das Setup ist die Eichung.</i>',
    ms('sCl', fs=9, fn=FONTI, al=TA_CENTER, ld=12)))


# ============ BUILD ============
def hf(canvas, doc):
    canvas.saveState(); canvas.setFont(FONT, 8); canvas.setFillColor(MIDGRAY)
    canvas.drawString(25*mm, A4[1]-15*mm, 'Dok. 0005 — Frequenzverschiebung als Indikator der Ansprache')
    canvas.drawRightString(A4[0]-25*mm, A4[1]-15*mm, f'Seite {doc.page}')
    canvas.drawCentredString(A4[0]/2, 12*mm, 'Querverweise: Dok. 0002 (v8) · Dok. 0003 (Impedanzvergleich) · Dok. 0004 (Zwei-Filter)')
    canvas.restoreState()

doc.build(story, onFirstPage=hf, onLaterPages=hf)
print(f'PDF created: {outpath}')
