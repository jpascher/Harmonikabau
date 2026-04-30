#!/usr/bin/env python3
"""Dok. 0024 — Praxishinweise: Kritische Handgriffe im Harmonikabau"""
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuB','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuI','/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBI','/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
addMapping('DejaVu',0,0,'DejaVu'); addMapping('DejaVu',1,0,'DejaVuB')
addMapping('DejaVu',0,1,'DejaVuI'); addMapping('DejaVu',1,1,'DejaVuBI')
DB=HexColor('#16213e'); AR=HexColor('#e94560'); KG=HexColor('#2e7d32'); WR=HexColor('#c62828')
LG=HexColor('#f5f5f5'); KBG=HexColor('#e8f5e9'); WBG=HexColor('#ffebee')
WP=A4[0]; PW=WP-50*mm
styles=getSampleStyleSheet()
for sn in styles.byName:
    s=styles.byName[sn]
    if hasattr(s,'fontName'):
        if 'Bold' in s.fontName: s.fontName='DejaVuB'
        elif 'Italic' in s.fontName: s.fontName='DejaVuI'
        else: s.fontName='DejaVu'
sT=ParagraphStyle('T',parent=styles['Title'],fontSize=18,textColor=DB,spaceAfter=4,alignment=TA_CENTER,fontName='DejaVuB')
sST=ParagraphStyle('ST',parent=styles['Normal'],fontSize=12,textColor=DB,alignment=TA_CENTER,spaceAfter=2,fontName='DejaVu')
sAb=ParagraphStyle('Ab',parent=styles['Italic'],fontSize=9,textColor=HexColor('#555555'),spaceAfter=8,fontName='DejaVuI')
sCh=ParagraphStyle('Ch',parent=styles['Heading1'],fontSize=14,textColor=DB,spaceBefore=14,spaceAfter=6,fontName='DejaVuB')
sB=ParagraphStyle('Bo',parent=styles['Normal'],fontSize=10,leading=14,spaceAfter=6,alignment=TA_JUSTIFY,fontName='DejaVu')
sK=ParagraphStyle('KB',parent=sB,fontSize=10,backColor=KBG,borderPadding=6,borderColor=KG,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sW=ParagraphStyle('WB',parent=sB,fontSize=10,backColor=WBG,borderPadding=6,borderColor=WR,borderWidth=1,spaceAfter=8,fontName='DejaVu')
def hr(): return Table([['']], colWidths=[PW], style=TableStyle([('LINEBELOW',(0,0),(-1,-1),2,AR),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))

# ══ PDF ══
outfile='0024_praxishinweise_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]

story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0024',sST))
story.append(Paragraph('Praxishinweise:<br/>Kritische Handgriffe im Harmonikabau',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Kritische Tätigkeiten bei Reparatur und Neubau von Harmonikas. '
    'Handgriffe, die Spezialwissen und Übung erfordern. '
    'Referenz: Dok.\u20090023 (Gehäuse und Mechanik).',sAb))

# Kap. 1
story.append(Paragraph('1. Handwerk und Spezialwissen',sCh))
story.append(Paragraph(
    'Viele Tätigkeiten im Harmonikabau sind grundsätzlich Arbeiten, die ein guter '
    'Handwerker problemlos bewältigt. Dennoch werden dabei Handfertigkeit und Ausdauer '
    'verlangt, die man nicht unterschätzen sollte. Bei vielen Arbeitsschritten ist ein '
    'gewisses Spezialwissen erforderlich, und es ist ratsam, sich die Tätigkeiten von '
    'jemandem vorzeigen zu lassen.',sB))
story.append(Paragraph(
    'Inzwischen gibt es diverse Videos, die einzelne Arbeitsschritte zeigen \u2014 '
    'auch von mir \u2014, doch diese sind nicht systematisch zusammengestellt. '
    'Besser ist es, die Tätigkeiten von einer Person direkt vorgeführt zu bekommen. '
    'Das Grundprinzip <i>Vormachen und Nachmachen</i> kann nie übertroffen werden, '
    'da unmittelbar Rückmeldungen gegeben und umgesetzt werden können. So wie ein '
    'Lernvideo zum Erlernen eines Musikstücks nie eine Person ersetzt, die daneben '
    'sitzt und korrigiert.',sB))
story.append(Paragraph(
    'Im Folgenden führe ich kritische Aufgaben auf \u2014 sowohl bei Reparaturen '
    'als auch beim Neuanfertigen \u2014, die naturgemäß alles andere als vollständig sind.',sB))

# Kap. 2
story.append(Paragraph('2. Diskantmechanik: Achsen entfernen und einsetzen',sCh))
story.append(Paragraph(
    'Die Diskantmechanik zu zerlegen und wieder zusammenzubauen ist eine der ersten '
    'Tätigkeiten, die man beherrschen muss. Es gibt unterschiedliche Mechaniken; die '
    'häufigsten sind die <b>offene traditionelle Mechanik</b> mit einer oder mehreren '
    'Achsen und die <b>verdeckte Metallmechanik</b> mit einer Achse.',sB))
story.append(Paragraph(
    'Der wichtigste Vorgang ist das Entfernen und Einsetzen der Achse. Dabei ist '
    'äußerste Vorsicht geboten, speziell bei offenen Mechaniken, da die Achsen fest '
    'sitzen können. Man versucht, eine Achse am offenen Ende zu klemmen und mit '
    'Rotationsbewegung und gleichzeitigem Ziehen herauszuarbeiten. Dabei kann es '
    'vorkommen, dass man den Holzkamm bricht, wenn zu hohe Zugkräfte aufgewendet werden.',sB))
story.append(Paragraph(
    'Man verwendet ein <b>Spannfutter</b> zum Klemmen der Achse oder direkt einen '
    '<b>Akkuschrauber</b>, mit dem man die Achse langsam rotiert und dann zieht. '
    'Es gibt auch Vorrichtungen, die speziell für diese Aufgabe gebaut sind. '
    'Sowohl die Achse als auch der Kamm können bei diesem Vorgang beschädigt werden.',sW))
story.append(Paragraph(
    'Bei der <b>verdeckten Metallmechanik</b> ist derselbe Vorgang möglich, jedoch '
    'ist die Lage einfacher: Die Achsen sind 3\u2009mm stark und man bricht nicht so '
    'leicht etwas.',sB))

# Kap. 3
story.append(Paragraph('3. Lagerung der Diskantmechanik',sCh))
story.append(Paragraph(
    'Beim Bau einer Diskantmechanik sind die <b>Lagerungen</b> das Schwierigste.',sB))
story.append(Paragraph(
    'Bei der verdeckten Mechanik bekommt man die Hebel mit den <b>Messinglaschen</b> '
    'und Bohrungen für die Achsen halbfertig. Die Lager müssen aufgebohrt und kalibriert '
    'werden \u2014 auf eine 3-mm-Stahlachse. An die Achse ist eine hohe Genauigkeit und '
    'Gleichmäßigkeit gefordert.',sB))
story.append(Paragraph(
    'Die Achsen haben etwas <b>Überlänge</b>. Man kürzt ein Stück ab; den Abfall '
    'verwendet man zum Kalibrieren der Lager. Die komplette Achse kommt erst in einem '
    'der letzten Arbeitsgänge zum Einsatz, wenn die Lager fertig kalibriert sind.',sB))
story.append(Paragraph(
    'Man bohrt die Bohrungen der Lager auf, aber <b>etwas kleiner als die Achse</b>. '
    'Dann feilt man auf das Abfallstück eine Fase, spannt das Stück in den Akkuschrauber '
    'und dreht die Achse in die Lager.',sB))
story.append(Paragraph(
    '<b>Wichtig:</b> Man verwendet <b>Polierpaste und Wasser</b> zum Kühlen. '
    'Das ist entscheidend \u2014 wenn man die Lager verreibt, ist die Lagerbuchse '
    'bereits zu groß. Am besten führt man das Aufreiben direkt im Behälter mit '
    'Kühlwasser durch.',sW))
story.append(Paragraph(
    'Den Vorgang wiederholt man mehrmals, zum Schluss mit der langen Achse. '
    'Das <b>Lagerspiel</b> darf nicht zu groß werden. Wird das Spiel dennoch zu viel, '
    'kann man das Lager durch <b>Quetschen</b> wieder verengen \u2014 es gibt dafür '
    'spezielle Zangen.',sB))

# Kap. 4
story.append(Paragraph('4. Abschluss',sCh))
story.append(Paragraph(
    'Weitere Hinweise gebe ich hier derzeit nicht. Es gibt viele weitere Handgriffe, '
    'aber die hier beschriebenen sind die wichtigsten, die jeder kennen sollte. '
    'Viele andere Tätigkeiten sind einfacher selbst zu erlernen oder von Videos abzusehen.',sB))

story.append(Spacer(1,6*mm))
story.append(Paragraph(
    '<i>Die beste Anleitung ist ein erfahrener Harmonikabauer, der einem über die '
    'Schulter schaut.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0024 \u2014 Praxishinweise \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
