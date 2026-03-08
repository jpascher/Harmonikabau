#!/usr/bin/env python3
"""Dok. 0023 — Gehäuse und Mechanik"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image
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
sF=ParagraphStyle('Fo',parent=sB,fontSize=10,alignment=TA_CENTER,spaceAfter=8,spaceBefore=4,fontName='DejaVu')
sK=ParagraphStyle('KB',parent=sB,fontSize=10,backColor=KBG,borderPadding=6,borderColor=KG,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sW=ParagraphStyle('WB',parent=sB,fontSize=10,backColor=WBG,borderPadding=6,borderColor=WR,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sTH=ParagraphStyle('TH',parent=sB,fontSize=8,fontName='DejaVuB',alignment=TA_CENTER,leading=10)
sTD=ParagraphStyle('TD',parent=sB,fontSize=8,alignment=TA_CENTER,leading=10,fontName='DejaVu')
sTL=ParagraphStyle('TDL',parent=sB,fontSize=8,alignment=TA_LEFT,leading=10,fontName='DejaVu')
def hr(): return Table([['']], colWidths=[PW], style=TableStyle([('LINEBELOW',(0,0),(-1,-1),2,AR),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
def mk_tbl(hdr,rows,cw=None):
    data=[[Paragraph(h,sTH) for h in hdr]]
    for row in rows:
        data.append([Paragraph(str(c),sTL) if i==0 else Paragraph(str(c),sTD) for i,c in enumerate(row)])
    cw=cw or [PW/len(hdr)]*len(hdr)
    t=Table(data,colWidths=cw,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),DB),('TEXTCOLOR',(0,0),(-1,0),white),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,LG]),('GRID',(0,0),(-1,-1),0.5,HexColor('#cccccc')),
        ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),2),('BOTTOMPADDING',(0,0),(-1,-1),2),
        ('FONTNAME',(0,0),(-1,-1),'DejaVu')])); return t

# ══ PDF ══
outfile='0023_gehaeuse_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0023',sST))
story.append(Paragraph('Geh\u00e4use und Mechanik',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('Konstruktion von Geh\u00e4use, Mechanik, Lagern und Hebeln '
    'bei vierreihigen diatonischen Instrumenten. Historische Entwicklung und '
    'moderne Tendenz. Akustische Relevanz. Referenz: Dok. 0005\u20130012, 0022.',sAb))

# Kap. 1: Gehäuse
story.append(Paragraph('1. Geh\u00e4use: Aufbau und Materialien',sCh))
story.append(Paragraph(
    'Das Geh\u00e4use eines vierreihigen diatonischen Instruments besteht aus '
    'dem <b>Diskantteil</b> (rechte Hand), dem <b>Bassteil</b> (linke Hand) '
    'und dem <b>Balg</b> dazwischen. Diskant- und Bassteil sind die tragenden '
    'Strukturen, die Stimmst\u00f6cke, Mechanik und Kanzellen aufnehmen.',sB))

rows_mat=[]
rows_mat.append(['Massivholz','Traditionell. Fichte, Ahorn, Nuss. Schwerer, schwingt etwas mit.'])
rows_mat.append(['Sperrholz / Multiplex','H\u00e4ufig. Leichter, formstabil, schwingt wenig.'])
rows_mat.append(['Aluminium','Leicht, stabil. Schwingt wenig, k\u00fchle Haptik.'])
rows_mat.append(['Kunststoff / GFK','Selten bei diatonischen. Leicht, unempfindlich.'])
story.append(mk_tbl(['Material','Eigenschaft'],rows_mat,cw=[30*mm,90*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Das Geh\u00e4usematerial hat \u2014 wie bereits in den Dokumenten zu '
    'Stimmst\u00f6cken und Basskammern (Dok. 0005\u20130008) gezeigt \u2014 '
    '<b>wenig Einfluss</b> auf den Klang, auch wenn man das nicht komplett '
    'ignorieren sollte. Die richtige Analogie ist eine <b>Lautsprecherbox</b>, '
    'nicht ein Saiteninstrument mit Resonanzk\u00f6rper. Was f\u00fcr Basskammern gilt, '
    'gilt auch f\u00fcr das Geh\u00e4use.',sB))
story.append(Paragraph(
    'Das <b>innere Volumen</b> wirkt sich ebenfalls wenig aus \u2014 sonst '
    'w\u00fcrde das Instrument wesentlich anders klingen bei ge\u00f6ffnetem und '
    'geschlossenem Balg. Die <b>Masse</b> des Instruments tr\u00e4gt mehr zum '
    'Klangverhalten bei als das Material.',sB))
story.append(Paragraph(
    'Das Geh\u00e4use sollte <b>wenig schwingen</b> \u2014 was aber von manchen '
    'Spielern negativ gesehen wird, die gerne die Schwingungen des Geh\u00e4uses '
    'wahrnehmen. <b>Zuh\u00f6rer nehmen das Instrument immer anders wahr als '
    'der Spieler</b> \u2014 der Spieler sp\u00fcrt Vibrationen, die der Zuh\u00f6rer '
    'nicht h\u00f6rt.',sK))

# Kap. 2: Diskantteil
story.append(Paragraph('2. Diskantteil',sCh))
story.append(Paragraph(
    'Der Diskantteil nimmt die <b>Diskant-Stimmst\u00f6cke</b>, die <b>Mechanik</b> '
    '(Tasten, Hebel, Lager, Klappen) und die <b>Diskant-Kanzellen</b> auf. '
    'Die Gr\u00f6\u00dfe des Diskantteils wird bestimmt durch:',sB))
story.append(Paragraph(
    '\u2022 Die Anzahl der Tasten (4 Reihen, typisch 46\u201360 Tasten)<br/>'
    '\u2022 Die Stimmstock-Anordnung (gerade oder schr\u00e4g gestellt \u2014 '
    'schr\u00e4g bei dickeren Stimmplatten, da tiefe Zungen weiter ausschwingen)<br/>'
    '\u2022 Die Mechanik-H\u00f6he (Lager, Hebel, Federung)<br/>'
    '\u2022 Register-Mechanik (Umschaltung der Ch\u00f6re)',sB))
story.append(Paragraph(
    'Moderne Diskantteile sind <b>h\u00f6her</b> als historische, da die '
    'Stimmst\u00f6cke schr\u00e4g gestellt werden m\u00fcssen (tiefere Zungen schwingen '
    'weiter aus) und die Lager mehr Platz brauchen (Dok. 0022: '
    'die gr\u00f6\u00dfere Bauh\u00f6he reduziert den Balg-Bewegungsbereich).',sB))

story.append(Paragraph('<b>Diskantboden: Material und Anforderungen</b>',sCh))
story.append(Paragraph(
    'Der Diskantboden (Auflagefl\u00e4che f\u00fcr die Stimmst\u00f6cke) muss drei '
    'Anforderungen erf\u00fcllen:',sB))
story.append(Paragraph(
    '\u2022 <b>Dichtheit:</b> Die Verbindung zu den Stimmst\u00f6cken muss '
    'luftdicht sein \u2014 jede Undichtigkeit kostet Ansprache und Lautst\u00e4rke.<br/>'
    '\u2022 <b>Formstabilit\u00e4t:</b> Der Boden darf sich nicht verziehen \u2014 '
    'Verzug \u00f6ffnet Spalte zwischen Boden und Stimmst\u00f6cken.<br/>'
    '\u2022 <b>M\u00f6glichst d\u00fcnn:</b> Die Materialst\u00e4rke des Bodens '
    'erh\u00f6ht die Tonlochtiefe. Bei den hohen T\u00f6nen (kurze Kan\u00e4le) '
    'wird das problematisch \u2014 der Str\u00f6mungswiderstand steigt, '
    'die Ansprache leidet.',sB))

rows_boden=[]
rows_boden.append(['D\u00fcnnes Sperrholz','Formstabil, d\u00fcnn m\u00f6glich','Bew\u00e4hrt'])
rows_boden.append(['Massives Ahorn (ausgesucht)','Sehr gut, sorgf\u00e4ltig ausw\u00e4hlen','Bew\u00e4hrt'])
rows_boden.append(['Aluminium','Formstabil, d\u00fcnn, dicht. Etwas brillanteres Klangbild.','Bew\u00e4hrt (Akkordeon)'])
rows_boden.append(['Karbon','Sehr d\u00fcnn und steif. Schwer zu fr\u00e4sen, teuer, unangenehmes Klangbild.','Nicht empfohlen'])
story.append(mk_tbl(['Material','Eigenschaft','Status'],rows_boden,cw=[32*mm,52*mm,24*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Holzb\u00f6den sind zu bevorzugen \u2014 ob das subjektiv ist oder beweisbar, '
    'bleibt offen. D\u00fcnnes Sperrholz oder ausgesuchtes massives Ahorn '
    'haben sich bew\u00e4hrt. Entscheidend ist nicht das Material an sich, '
    'sondern dass der Boden <b>d\u00fcnn, dicht und formstabil</b> ist. '
    'Ein d\u00fcnner Boden, der sich verzieht, ist schlechter als ein '
    'etwas dickerer, der plan bleibt.',sB))

# Kap. 3: Bassteil
story.append(Paragraph('3. Bassteil',sCh))
story.append(Paragraph(
    'Der Bassteil nimmt die <b>Bass-Stimmst\u00f6cke</b>, die <b>Bass-Mechanik</b> '
    'und die <b>Begleiter-Mechanik</b> auf. Bei diatonischen Instrumenten '
    'unterscheidet man <b>Einzelton-Bass</b> (jede Taste ein Ton) und '
    '<b>Akkord-Bass</b> (Taste l\u00f6st Grundton + Terz + Quinte aus).',sB))
story.append(Paragraph(
    'Die Bass-Mechanik ist <b>nicht komplexer</b> als die Diskant-Mechanik. '
    'Speziell die Diskant-Mechanik f\u00fcr dreich\u00f6rige Instrumente oder '
    'mit Registern ist viel anspruchsvoller zu fertigen. '
    'Jedes kleinste Ger\u00e4usch in der Diskant-Mechanik ist <b>sehr st\u00f6rend</b>, '
    'w\u00e4hrend leichte Ger\u00e4usche beim Bass eher in Kauf genommen werden.',sB))
story.append(Paragraph(
    'Die <b>Helikonb\u00e4sse</b> (tiefste Oktave) brauchen die gr\u00f6\u00dften '
    'Stimmplatten und damit den meisten Platz. '
    'Die Anordnung der Bass-Stimmst\u00f6cke bestimmt '
    'wesentlich die Breite und Tiefe des Bassteils.',sB))

# Kap. 4: Mechanik
story.append(Paragraph('4. Mechanik: Hebel, Lager, Klappen',sCh))
story.append(Paragraph(
    'Die Mechanik \u00fcbertr\u00e4gt den Tastendruck auf die Klappen \u00fcber den '
    'Stimmst\u00f6cken. Die Qualit\u00e4t der Mechanik bestimmt das <b>Spielgef\u00fchl</b> '
    'und die <b>Ger\u00e4uschentwicklung</b>.',sB))

rows_mech=[]
rows_mech.append(['Hebel','Aluminium oder Stahl','Verbindung Taste\u2013Klappe. Verst\u00e4rkt gegen\u00fcber fr\u00fcher (weniger Eigenschwingung).'])
rows_mech.append(['Lager','Achsen in Lagerb\u00f6cken','Pr\u00e4ziser als fr\u00fcher, aber mehr Platzbedarf.'])
rows_mech.append(['Klappen','Holz, Leder, Filz','Verschlie\u00dfen die Kanzellen\u00f6ffnungen. Dichtung: Leder oder Filz.'])
rows_mech.append(['Federn','Stahldraht','R\u00fcckstellkraft. Bestimmt Tastengewicht und Anschlaggef\u00fchl.'])
rows_mech.append(['Tastatur','Holz, Kunststoff, Perlmutt','Spieloberfl\u00e4che. Kn\u00f6pfe bei diatonischen.'])
story.append(mk_tbl(['Bauteil','Material','Funktion'],rows_mech,cw=[18*mm,30*mm,68*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Lagerqualit\u00e4t wird \u00fcberbewertet:</b> Jedes Lager, das ger\u00e4uscharm '
    'ist und wenig Reibung verursacht, ist geeignet. Die Bewegungen der Mechanik '
    'sind <b>minimal</b> (wenige Millimeter Hub) und es liegen <b>keine hohen '
    'Geschwindigkeiten</b> vor. Der Einsatz von Kugellagern ist zwar m\u00f6glich, '
    'bringt aber objektiv <b>keine Verbesserung</b> gegen\u00fcber einfachen '
    'Gleitlagern \u2014 f\u00fcr die Haltbarkeit und Abnutzung jedoch schon '
    'g\u00fcnstiger. Einfache Achsen in sauber gebohrten Lagerb\u00f6cken '
    'funktionieren seit \u00fcber 100 Jahren einwandfrei.',sB))

story.append(Paragraph('<b>Klappen: Anforderungen und Aufbau</b>',sCh))
story.append(Paragraph(
    'Die Anforderungen an die Klappen sind vielf\u00e4ltig und teilweise '
    'widerspr\u00fcchlich: <b>Wenig Ger\u00e4usche, Dichtheit, wenig Gewicht '
    'und wenig Bauh\u00f6he.</b>',sB))
story.append(Paragraph(
    '<b>Aufbau (von unten nach oben):</b>',sB))
story.append(Paragraph(
    '\u2022 <b>Deckel:</b> Wei\u00dfblech (d\u00fcnn, aber etwas schwerer), '
    'Holz, Hartkarton oder Resopal \u2014 alles geeignet.<br/>'
    '\u2022 <b>Filz:</b> Naturfilz, nicht zu dick. Dicker Filz sitzt ein '
    'und durch Alterung vergr\u00f6\u00dfert sich der Tastenhub. '
    'Bew\u00e4hrt: <b>4 mm</b> f\u00fcr Reihe 3 und 4, <b>2 mm</b> f\u00fcr '
    'die erste Reihe.<br/>'
    '\u2022 <b>Leder:</b> Sehr weich und d\u00fcnn \u2014 dient als Dichtung und '
    'd\u00e4mpft das Aufschlagger\u00e4usch. D\u00fcnnes Leder ist hier entscheidend '
    'f\u00fcr leises Schlie\u00dfen.<br/>'
    '\u2022 <b>Kleber:</b> Zwischen Filz und Leder sowie zwischen Filz und Deckel. '
    'Eine besondere Herausforderung. '
    'Moderne Klebstoffe auf <b>Silikonbasis</b> oder \u00e4hnliche flexible, '
    'aber zugleich dichtende Materialien sind optimal.',sB))

story.append(Paragraph(
    '<b>Leder am Boden:</b> Manche Bauer bringen am Diskantboden ebenfalls '
    'eine Lederschicht an. Das ist eher <b>kontraproduktiv</b> \u2014 '
    'es erh\u00f6ht die Tonlochtiefe und damit den Str\u00f6mungswiderstand '
    'bei den hohen T\u00f6nen.',sB))

story.append(Paragraph(
    '<b>\u00dcbergang Klappe\u2013Hebel:</b> Untergeordnet, aber reduziert '
    'die Ger\u00e4usche, wenn ein <b>Gummielement</b> dazwischen sitzt. '
    'Der <b>Druckpunkt</b> ist von Bedeutung \u2014 wenn die Klappen flexibel '
    'am Hebel befestigt werden, gibt es einen sp\u00fcrbaren Druckpunkt beim Spielen.',sB))

story.append(Paragraph(
    '<b>Federdruck vs. Dichtheit:</b> Je gr\u00f6\u00dfer die Kanzellen\u00f6ffnungen, '
    'desto mehr Federkraft ist bei den Hebeln erforderlich, um die Dichtheit '
    'bei maximalem Balgdruck zu gew\u00e4hrleisten. Federdruck und Dichtheit '
    'sind <b>kontraproduktiv</b> zum Spielkomfort \u2014 mehr Feder = '
    'h\u00e4rterer Anschlag. Der Spielraum ist gering.',sW))

# Kap. 5: Historische Entwicklung
story.append(Paragraph('5. Historische Entwicklung',sCh))
story.append(Paragraph(
    'Bemerkenswert wenig ist seit den Anf\u00e4ngen wirklich neu erfunden worden \u2014 '
    '<b>Verkaufsprospekte bereits vor dem Ersten Weltkrieg</b> zeigen die gleichen '
    'Grundkonstruktionen. Was sich verbessert hat:',sB))

rows_hist=[]
rows_hist.append(['Mechanik','Einfach, direkt','Pr\u00e4ziser, ger\u00e4usch\u00e4rmer, aber mehr Platzbedarf'])
rows_hist.append(['Lager','Einfache Achsen','Pr\u00e4zisionslager, weniger Spiel, aber gr\u00f6\u00dfer'])
rows_hist.append(['Hebel','D\u00fcnn, leicht','Verst\u00e4rkt (weniger Eigenschwingung)'])
rows_hist.append(['Stimmst\u00f6cke','Gerade eingebaut','Schr\u00e4g gestellt (tiefe Zungen schwingen weiter)'])
rows_hist.append(['Stimmplatten','D\u00fcnner im Diskant','Dicker (bessere Ansprache, aber mehr Platzbedarf)'])
rows_hist.append(['Geh\u00e4use','Kompakt, niedrig','H\u00f6her (mehr Einbauten)'])
story.append(mk_tbl(['Bauteil','Historisch','Modern'],rows_hist,cw=[22*mm,32*mm,60*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Konsequenz:</b> Alle Verbesserungen zusammen f\u00fchren zu <b>gr\u00f6\u00dferen '
    'Instrumenten</b>. Historische Instrumente mit kleineren Geh\u00e4useabmessungen '
    'existieren durchaus \u2014 sie funktionieren, sind aber weniger pr\u00e4zise '
    'in der Mechanik. Das ist ein Trade-off: '
    'Pr\u00e4zision und Komfort kosten Platz und Gewicht.',sK))

# Kap. 6: Akustische Relevanz
story.append(Paragraph('6. Akustische Relevanz des Geh\u00e4uses',sCh))
story.append(Paragraph(
    'Das Geh\u00e4use ist <b>kein Resonanzk\u00f6rper</b> wie bei Saiteninstrumenten. '
    'Die korrekte Analogie ist eine <b>Lautsprecherbox</b>: Die W\u00e4nde sollen '
    'den Schall abstrahlen, aber m\u00f6glichst wenig selber schwingen. '
    'Was akustisch relevant ist:',sB))
story.append(Paragraph(
    '\u2022 <b>Masse:</b> Tr\u00e4gt mehr zum Klangverhalten bei als das Material. '
    'Schwere Instrumente schwingen weniger \u2192 neutralerer Klang.<br/>'
    '\u2022 <b>Steifigkeit:</b> Steife W\u00e4nde schwingen weniger. '
    'D\u00fcnne Holzw\u00e4nde schwingen mit \u2014 das kann erw\u00fcnscht oder '
    'unerw\u00fcnscht sein, ist aber <b>kein Qualit\u00e4tsmerkmal</b>.<br/>'
    '\u2022 <b>Mechanikraum als Kammfilter:</b> Der Hohlraum der Diskant-Mechanik '
    'bildet einen passiven Resonator, der bestimmte Frequenzen d\u00e4mpft '
    '(Saugkreis-Prinzip, Dok. 0005\u20130008).<br/>'
    '\u2022 <b>Verdeck als Filter:</b> Bestimmt, welche Frequenzen '
    'ged\u00e4mpft werden. Geschlossen = st\u00e4rkere D\u00e4mpfung der H\u00f6hen.<br/>'
    '\u2022 <b>Inneres Volumen:</b> Wirkt sich wenig aus \u2014 sonst w\u00fcrde '
    'das Instrument bei ge\u00f6ffnetem und geschlossenem Balg wesentlich '
    'anders klingen.<br/>'
    '\u2022 <b>Dichtheit:</b> Undichtigkeiten reduzieren den '
    'effektiven Balgdruck und damit Lautst\u00e4rke und Ansprache (Dok. 0010).',sB))
story.append(Paragraph(
    '<b>Spieler vs. Zuh\u00f6rer:</b> Der Spieler nimmt Vibrationen des Geh\u00e4uses '
    'wahr (taktil und \u00fcber die K\u00f6rperschallleitung). Der Zuh\u00f6rer h\u00f6rt '
    'davon nichts \u2014 er h\u00f6rt nur den abgestrahlten Schall. '
    'Manche Spieler bevorzugen ein schwingendes Geh\u00e4use, weil es sich '
    '\u201elebendiger\u201c anf\u00fchlt. Akustisch f\u00fcr den Zuh\u00f6rer macht es '
    'kaum einen Unterschied.',sB))

story.append(Paragraph('<b>Einfluss des Holzes f\u00fcr B\u00f6den und Rahmen</b>',sCh))
story.append(Paragraph(
    'Das verwendete Holz f\u00fcr B\u00f6den und Rahmen der K\u00e4sten hat '
    '<b>wenig, aber wahrnehmbare</b> Effekte. Alles ist brauchbar \u2014 Fichte, '
    'Ahorn, Nuss, Sperrholz \u2014 aber die Feinheiten, wie sich das Material '
    'auswirkt, sind <b>schwer einsch\u00e4tzbar</b> und mit heutigen Mitteln '
    '<b>nicht berechenbar</b>. Es gibt keine Formel, die den Klangunterschied '
    'zwischen einem Fichten- und einem Ahornboden vorhersagt. '
    'Man ist auf Erfahrungswerte und subjektive Beurteilung angewiesen.',sB))

story.append(Paragraph('<b>Verdeck und Diskant-Mechanikraum</b>',sCh))
story.append(Paragraph(
    'Der Raum zwischen Stimmst\u00f6cken und Verdeck \u2014 der Bereich der '
    'Diskant-Mechanik \u2014 bildet einen <b>minimalen Resonanzk\u00f6rper</b> aus, '
    'der als <b>Kammfilter</b> wirkt (Dok. 0005\u20130008: Saugkreis-Prinzip). '
    'Wie bei den Basskammern: Er <b>d\u00e4mpft</b> bei bestimmten Frequenzen, '
    'er <b>verst\u00e4rkt nicht</b>.',sB))
story.append(Paragraph(
    'Das <b>Verdeck</b> (Diskantabdeckung) wirkt als <b>Filter</b>. '
    'Ob geschlossen, mit Jalousien oder offen \u2014 es wird damit nur '
    'etwas ged\u00e4mpft. Am transparentesten klingt der Diskant entweder '
    '<b>ohne Verdeck</b> oder mit einem Verdeck mit m\u00f6glichst vielen '
    'Durchbr\u00fcchen.',sB))

story.append(Paragraph('<b>Mathematischer Beleg: Nur die Summe z\u00e4hlt</b>',sCh))
story.append(Paragraph(
    'Jede \u00d6ffnung im Verdeck hat eine akustische Impedanz, die umgekehrt '
    'proportional zur Fl\u00e4che ist. Mehrere \u00d6ffnungen wirken als '
    '<b>Parallelschaltung</b>:',sB))
story.append(Paragraph(
    'Z_gesamt = (\u03c1\u00b7c + j\u03c9\u03c1\u00b7l_eff) / A_total',sF))
story.append(Paragraph(
    'Da alle \u00d6ffnungen dieselbe Materialdicke haben, h\u00e4ngt die '
    'Gesamtimpedanz <b>nur von A_total = \u03a3 A_i</b> ab \u2014 '
    'nicht von der Form der einzelnen Durchbr\u00fcche. '
    '1 gro\u00dfe \u00d6ffnung = 100 kleine L\u00f6cher, solange die Gesamtfl\u00e4che '
    'gleich ist. Die Form (Kreis, Rechteck, Schlitz, Stern) spielt keine Rolle, '
    'weil die akustische Impedanz nur von der Fl\u00e4che abh\u00e4ngt \u2014 '
    'solange die \u00d6ffnungsdimensionen klein sind gegen\u00fcber der '
    'Wellenl\u00e4nge (\u03bb = 34 cm bei 1 kHz).',sB))
story.append(Paragraph(
    '<b>Position der Durchbr\u00fcche:</b> Die Position wirkt <b>indirekt</b> '
    '\u2014 in den Bereichen, wo sich Klappen \u00f6ffnen, gelangt der Schall '
    'am direktesten nach au\u00dfen. \u00d6ffnungen an anderen Stellen werden '
    'durch den Mechanikraum zus\u00e4tzlich gefiltert (Kammfilter). '
    'Durchbr\u00fcche \u00fcber den Klappen sind daher akustisch wirksamer '
    'als Durchbr\u00fcche am Rand.',sB))
story.append(Paragraph(
    'Beim Bass gilt dasselbe Prinzip: Die <b>Anzahl und Gr\u00f6\u00dfe der '
    '\u00d6ffnungen am Bassboden</b> bestimmt das Klangbild der B\u00e4sse \u2014 '
    'als Filter, nicht als Verst\u00e4rker. '
    'Mehr \u00d6ffnungen = weniger D\u00e4mpfung der Oberfrequenzen (brillanter). '
    'Weniger \u00d6ffnungen = st\u00e4rkere D\u00e4mpfung (w\u00e4rmerer Bassklang).',sB))
story.append(Paragraph(
    '<b>Zusammengefasst:</b> Mechanikraum, Verdeck und Bass\u00f6ffnungen '
    'wirken als <b>passive Filter</b> (Kammfilter, D\u00e4mpfung). '
    'Material der W\u00e4nde \u2192 wenig Einfluss. '
    '\u00d6ffnungsgeometrie \u2192 <b>gro\u00dfer Einfluss</b>, '
    'weil sie bestimmt, was ged\u00e4mpft wird und was nicht.',sK))

# Kap. 7: Geräuschminimierung
story.append(Paragraph('7. Ger\u00e4uschminimierung',sCh))
story.append(Paragraph(
    'Mechanische Ger\u00e4usche (Klappern, Klicken, Aufschlag) sind ein '
    'wesentliches Qualit\u00e4tsmerkmal. Quellen:',sB))
story.append(Paragraph(
    '\u2022 <b>Klappenaufschlag:</b> Filz- oder Lederd\u00e4mpfung auf den Kanzellen. '
    'Zu viel D\u00e4mpfung reduziert den Ventilhub.<br/>'
    '\u2022 <b>Lagerspiel:</b> Pr\u00e4zise Lager reduzieren seitliches Klappern. '
    'Zu eng eingepasste Lager erh\u00f6hen die Reibung.<br/>'
    '\u2022 <b>Hebelresonanz:</b> Verst\u00e4rkte Hebel haben h\u00f6here Eigenfrequenzen, '
    'die au\u00dferhalb des h\u00f6rbaren Bereichs liegen. D\u00fcnne Hebel k\u00f6nnen '
    'mitschwingen und \u201eklingeln\u201c.<br/>'
    '\u2022 <b>Tastenanschlag:</b> Der obere und untere Anschlag der Tasten '
    'braucht D\u00e4mpfung (Filz, Gummi). Ohne D\u00e4mpfung: h\u00f6rbares Klicken.',sB))
story.append(Paragraph(
    'Alle diese Ma\u00dfnahmen erh\u00f6hen den <b>Platzbedarf</b> und das '
    '<b>Gewicht</b> \u2014 das ist der Preis f\u00fcr eine leise Mechanik.',sB))

# Kap. 8: Praktische Empfehlungen
story.append(Paragraph('8. Praktische Empfehlungen f\u00fcr den Harmonikabauer',sCh))
story.append(Paragraph(
    'Grunds\u00e4tzlich ist jeder, der Harmonikas baut, gut beraten, '
    '<b>bew\u00e4hrte Instrumente eher zu kopieren als Neues auszuprobieren</b>. '
    'Die Neuerungen \u2014 wenn sie noch so innovativ wirken m\u00f6gen \u2014 sind '
    'alle mit kleinen Variationen bereits einmal gebaut worden. '
    'Fast alles funktioniert.',sB))
story.append(Paragraph(
    '\u00dcber Details mag man unterschiedlicher Meinung sein. '
    '<b>Gewohnheiten spielen viel mehr Rolle als objektive Bewertungen.</b> '
    'Der Kunde ist K\u00f6nig.',sK))

story.append(Paragraph('<b>Subjektive, aber wichtige Details:</b>',sB))
story.append(Paragraph(
    'Kein Mensch ist anatomisch komplett gleich. Daher sind die folgenden '
    'Punkte hochgradig individuell \u2014 es gibt keine \u201erichtige\u201c L\u00f6sung, '
    'nur die passende f\u00fcr den jeweiligen Spieler:',sB))

rows_subj=[]
rows_subj.append(['Knopfgr\u00f6\u00dfe und -form','Rund, oval, flach, gew\u00f6lbt \u2014 Geschmackssache'])
rows_subj.append(['Knopfabst\u00e4nde','Abh\u00e4ngig von Fingergr\u00f6\u00dfe und Spieltechnik'])
rows_subj.append(['Position der Bassriemen','Individuell anpassen, kein Standard'])
rows_subj.append(['Bassboden: Leder ja/nein','Tradition vs. Praktikabilit\u00e4t'])
rows_subj.append(['Bassriemen-Form','Breit, schmal, gepolstert, nicht gepolstert \u2014 alles machbar'])
rows_subj.append(['Riemen-Verstellungen','Position und Art der Befestigung individuell'])
rows_subj.append(['Neigung Griffbrett','Mag eine Rolle spielen, anatomisch bedingt'])
rows_subj.append(['Balgschoner','\u00dcberfl\u00fcssig in den meisten F\u00e4llen, aber heute \u00fcblich \u2014 '
    'in Varianten, die eigentlich kein Balgschoner mehr sind'])
rows_subj.append(['Luftknopf / Lufttaste','Fr\u00fcher Luftknopf (reichte, kein Problem). '
    'Heute Lufttaste im Daumenbereich \u2014 nur Gewohnheit, nicht besser.'])
story.append(mk_tbl(['Detail','Anmerkung'],rows_subj,cw=[32*mm,82*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Alle diese Punkte sind <b>machbar</b> \u2014 keiner davon ist ein '
    'technisches Problem. Es sind Kundenwu\u0308nsche, und der Instrumentenbauer '
    'sollte sie ernst nehmen. Die beste Harmonika ist die, '
    'die der Spieler nicht mehr wahrnimmt, weil alles passt.',sB))

story.append(Paragraph('<b>Tragriemen</b>',sCh))
story.append(Paragraph(
    'Das selbe subjektive Bild ergibt sich bei den Tragriemen: '
    'Ein Riemen oder zwei, Art und Breite \u2014 alles individuell. '
    'Ein gut <b>gepolsterter Riemen im Bereich der Schulter</b> ist heute '
    'bevorzugt. Zwei Riemen mit R\u00fcckenverbinder erm\u00f6glichen '
    'einigema\u00dfen stabiles Spielen im Stehen \u2014 aber auch da gilt '
    'eindeutig: Was gew\u00fcnscht wird, nicht was besser ist.',sB))
story.append(Paragraph(
    '<b>Historisch:</b> Steirische wurden mit <b>einem Riemen</b> f\u00fcrs '
    'Spielen im Sitzen gebaut. Die Riemenbefestigung war <b>mittig oben '
    'am Balg</b> \u2014 das gew\u00e4hrleistete eine stabilere Anbindung des '
    'Instruments. Heute ist die Befestigung <b>nach hinten gewandert</b> '
    '(wie beim Akkordeon) \u2014 praktischer f\u00fcr zwei Riemen und Stehspiel, '
    'aber eine andere Gewichtsverteilung.',sB))

# Kap. 9
story.append(Paragraph('9. Zusammenfassung',sCh))
pts=[
    '<b>Geh\u00e4use = Lautsprecherbox, nicht Resonanzk\u00f6rper.</b> '
    'W\u00e4nde und Material haben wenig Einfluss (Dok. 0005\u20130008). '
    'Masse tr\u00e4gt mehr bei als das Material.',
    '<b>Holz f\u00fcr B\u00f6den/Rahmen:</b> Wenig, aber wahrnehmbare Effekte. '
    'Nicht berechenbar, nur Erfahrungswerte.',
    '<b>Verdeck und \u00d6ffnungen = passive Filter.</b> '
    'Mechanikraum wirkt als Kammfilter (d\u00e4mpft, verst\u00e4rkt nicht). '
    'Verdeck und Bass\u00f6ffnungen bestimmen, was ged\u00e4mpft wird.',
    '<b>Mechanik:</b> Hebel, Lager, Klappen, Federn. Heute pr\u00e4ziser und '
    'ger\u00e4usch\u00e4rmer, aber mit mehr Platzbedarf.',
    '<b>Wenig Neues erfunden:</b> Grundkonstruktion seit \u00fcber 100 Jahren '
    'unver\u00e4ndert. Verbesserungen: Pr\u00e4zision, D\u00e4mpfung, Verst\u00e4rkung.',
    '<b>Moderne Instrumente gr\u00f6\u00dfer:</b> Schr\u00e4ge Stimmst\u00f6cke, dickere '
    'Stimmplatten, pr\u00e4zisere Lager \u2192 mehr Platz n\u00f6tig \u2192 '
    'gr\u00f6\u00dfere Bauh\u00f6he \u2192 weniger Balg-Bewegungsbereich (Dok. 0022).',
    '<b>Inneres Volumen:</b> Wenig Einfluss \u2014 sonst w\u00fcrde das Instrument '
    'bei ge\u00f6ffnetem und geschlossenem Balg anders klingen.',
    '<b>Spieler \u2260 Zuh\u00f6rer:</b> Der Spieler sp\u00fcrt Vibrationen, '
    'die der Zuh\u00f6rer nicht h\u00f6rt. Ein schwingendes Geh\u00e4use f\u00fchlt '
    'sich lebendiger an, klingt aber f\u00fcr den Zuh\u00f6rer kaum anders.',
    '<b>Trade-off:</b> Pr\u00e4zision und Komfort kosten Platz und Gewicht. '
    'Historische Instrumente waren kompakter, aber die Mechanik weniger pr\u00e4zise.',
    '<b>Bew\u00e4hrtes kopieren, nicht neu erfinden.</b> '
    'Gewohnheiten spielen mehr Rolle als objektive Bewertungen. '
    'Der Kunde ist K\u00f6nig \u2014 alle Details sind machbar.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Das Geh\u00e4use ist eine Lautsprecherbox, kein Resonanzk\u00f6rper \u2014 '
    'der Spieler sp\u00fcrt es, der Zuh\u00f6rer h\u00f6rt es kaum.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0023 \u2014 Geh\u00e4use und Mechanik \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
