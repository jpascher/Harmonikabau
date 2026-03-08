#!/usr/bin/env python3
"""Dok. 0021 — Stimmplatten: Qualität, Hersteller, Güteklassen"""
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
outfile='0021_stimmplatten_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0021',sST))
story.append(Paragraph('Stimmplatten:<br/>Qualit\u00e4t, Hersteller und G\u00fcteklassen',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('\u00dcbersicht der aktuellen Stimmplattenhersteller f\u00fcr Harmonikainstrumente. '
    'Mensuren, G\u00fcteklassen, Fertigungsmethoden und eine n\u00fcchterne Einsch\u00e4tzung '
    'der Qualit\u00e4tsunterschiede. '
    'Grundlagen: siehe Wikipedia-Artikel \u201eDurchschlagende Zunge\u201c. '
    'Referenz: Dok. 0010\u20130012.',sAb))

# Kap. 1
story.append(Paragraph('1. Hersteller: Wer fertigt Stimmplatten?',sCh))
story.append(Paragraph(
    'Derzeit gibt es weltweit nur eine Handvoll Hersteller von Stimmplatten '
    'f\u00fcr Harmonikainstrumente. Die meisten sind in Italien ans\u00e4ssig, '
    '<b>Harmonikas</b> (Eigent\u00fcmer Titlbach) in Tschechien. '
    'Alle Arbeitsg\u00e4nge erfolgen '
    'nach wie vor <b>halbmanuell</b> \u2014 es gibt keine vollautomatische '
    'Fertigung.',sB))

rows_herst=[]
rows_herst.append(['Harmonikas (Titlbach)','Tschechien','Vollsortiment, alle G\u00fcteklassen. Zus\u00e4tzlich: EDM, Langschliff (Nastrino), DIX-Nachbau','Standard + Kunden + DIX'])
rows_herst.append(['Antonelli / Salpa (AV)','Castelfidardo (IT)','Vollsortiment','Standard-Mensuren'])
rows_herst.append(['Artigianale','Castelfidardo (IT)','Handgefertigt, h\u00f6here Klassen','Standard-Mensuren'])
rows_herst.append(['Cagnoni','Castelfidardo (IT)','Vollsortiment','Standard-Mensuren'])
rows_herst.append(['Binci','Castelfidardo (IT)','Maschinell konisch gefeilt, \u00e4ltere Maschinen','Standard, eher breitere Zungen'])
story.append(mk_tbl(['Hersteller','Standort','Besonderheit','Mensuren'],rows_herst,
    cw=[28*mm,28*mm,42*mm,24*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Binci</b> arbeitet mit \u00e4lteren Maschinen und feilt die Tonzungenkan\u00e4le '
    'maschinell konisch. Bei allen anderen Herstellern ist die Stanzung oder '
    'das Erodieren mit modernen Werkzeugen so genau, dass <b>kein Feilen '
    'der Kan\u00e4le erforderlich</b> ist. Manuelles Feilen der Kan\u00e4le erfolgt '
    'bei keinem Hersteller \u2014 au\u00dfer man gibt das speziell in Auftrag. '
    '<b>Harmonikas (Titlbach)</b> ist der einzige Hersteller, der erodierte '
    'Stimmplatten (EDM = Electrical Discharge Machining) anbietet. '
    'EDM-Stimmplatten haben eine glatte Oberfl\u00e4che und sehr enge Spalte '
    '(Dok. 0010: Ra = 1,6 \u00b5m, hydraulisch glatt).',sB))

# Kap. 2
story.append(Paragraph('2. Mensuren',sCh))
story.append(Paragraph(
    'Die Mensuren (Abmessungen der Stimmplatten) sind bei den Akkordeon-Herstellern '
    '<b>sehr \u00e4hnlich</b>. Das liegt daran, dass die Stimmst\u00f6cke und Kanzellen '
    'der Akkordeons genormt sind \u2014 die Stimmplatten m\u00fcssen passen. '
    'Leichte Variationen gibt es bei der <b>Zungenbreite</b> '
    '(Dok. 0017), aber L\u00e4nge und Breite der Stimmplatten sind '
    'weitgehend einheitlich.',sB))
story.append(Paragraph(
    'Wo sich die Hersteller unterscheiden, ist die <b>Zungenbreite</b>: '
    'Binci, Cagnoni und Harmonikas (Titlbach) bieten teilweise '
    '<b>breitere Zungenmensuren</b> an. Breitere Zungen bewegen mehr Luft '
    'pro Schwingung (\u2192 mehr Lautst\u00e4rke, Dok. 0017: Volumenstrom \u221d W) '
    'und haben eine h\u00f6here Torsionssteifigkeit (\u221d W\u00b3). '
    'Die Tonh\u00f6he wird durch die Breite nicht beeinflusst (Dok. 0017). '
    'Harmonikas-Stimmplatten werden auf Bestellung auch nach '
    'Kundenmensur gefertigt.',sB))

# Kap. 3: Materialien
story.append(Paragraph('3. Materialien: Rahmen und Zungen',sCh))
story.append(Paragraph('<b>Rahmenmaterial:</b>',sB))
story.append(Paragraph(
    'Das Rahmenmaterial der Stimmplatte kann gew\u00e4hlt werden \u2014 '
    'allerdings bietet nicht jeder Hersteller alle Materialien an:',sB))

rows_mat=[]
rows_mat.append(['Messing','Traditionell, schwer, gute D\u00e4mpfung','Weit verbreitet'])
rows_mat.append(['Zink','Leichter als Messing, verschiedene H\u00e4rten','H\u00e4ufig'])
rows_mat.append(['Duralumin','Leicht, verschiedene H\u00e4rten','Gewichtsersparnis'])
rows_mat.append(['Magnesium','Sehr leicht','Selten, Spezialanfertigung'])
rows_mat.append(['Neusilber','Hart, korrosionsbest\u00e4ndig','Selten, h\u00f6here Preisklasse'])
story.append(mk_tbl(['Material','Eigenschaften','Verf\u00fcgbarkeit'],rows_mat,
    cw=[22*mm,52*mm,38*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph('<b>Zungenmaterial:</b>',sB))
story.append(Paragraph(
    'Beim Zungenmaterial gibt es bei Sonderanfertigungen zwar M\u00f6glichkeiten '
    'f\u00fcr andere Materialien, aber beim Akkordeon ist praktisch ausschlie\u00dflich '
    '<b>Federstahl</b> in Verwendung. Der Stahl stammt heute bei <b>allen</b> '
    'Herstellern vom selben Lieferanten \u2014 der <b>voestalpine</b> \u2014 in '
    'einheitlicher H\u00e4rte. Kein Hersteller verwendet anderen Stahl.',sB))
story.append(Paragraph(
    '<b>Blue Star (Antonelli):</b> Antonelli hat in letzter Zeit eine '
    '\u201eBlue Star\u201c-Serie eingef\u00fchrt. Die Blauf\u00e4rbung der Zungen '
    'hat jedoch <b>keinen nachweisbaren akustischen Effekt</b>. '
    'Vergleiche zeigen, dass es sich eher um ein Verkaufsargument handelt. '
    'Erste Muster waren qualitativ eher vergleichbar mit Tipo a mano '
    'als mit einer h\u00f6heren Klasse, was inzwischen aber verbessert '
    'werden konnte.',sW))

# Kap. 4 (was vorher 3 war)
story.append(Paragraph('4. G\u00fcteklassen',sCh))
story.append(Paragraph(
    'Die Hersteller bieten ihre Stimmplatten in verschiedenen G\u00fcteklassen an. '
    'Die Bezeichnungen variieren, aber das Prinzip ist \u00e4hnlich:',sB))

rows_guete=[]
rows_guete.append(['Tipo a mano','Handgefertigt (a mano)','H\u00f6chste Qualit\u00e4t. Zungen einzeln profiliert und eingepasst.'])
rows_guete.append(['Export / Super','Maschinell + Nacharbeit','Gute Qualit\u00e4t. Maschinelle Vorfertigung, manuelle Endbearbeitung.'])
rows_guete.append(['Standard / Durall','Maschinell','Standardqualit\u00e4t. F\u00fcr die meisten Instrumente ausreichend.'])
rows_guete.append(['Economy / Tourist','Maschinell, minimal','Einstiegsqualit\u00e4t. Funktioniert, aber gr\u00f6\u00dfere Toleranzen.'])
story.append(mk_tbl(['Klasse','Fertigung','Beschreibung'],rows_guete,
    cw=[25*mm,30*mm,68*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die G\u00fcteklassen unterscheiden sich haupts\u00e4chlich in der '
    '<b>Passgenauigkeit</b> (Spaltweite, Dok. 0010), der <b>Profilierung</b> '
    '(Dok. 0016) und der <b>Oberfl\u00e4chenqualit\u00e4t</b>. '
    'H\u00f6here Klassen haben engere Spalte (\u2192 bessere Ansprache), '
    'feinere Profilierung (\u2192 gleichm\u00e4\u00dfigere Ansprache) '
    'und glattere Oberfl\u00e4chen (\u2192 weniger Luftverluste).',sB))

# Kap. 4
story.append(Paragraph('5. N\u00fcchterne Einsch\u00e4tzung: Wie gro\u00df sind die Unterschiede?',sCh))
story.append(Paragraph(
    'Die Variationsm\u00f6glichkeiten, dass unterschiedliche Stimmplatten '
    'unterschiedlich klingen oder unterschiedlich gute Ansprache besitzen, '
    'ist <b>minimal</b>. Ein exaktes Messverfahren, um dies zu kategorisieren, '
    'gibt es nicht \u2014 man ist auf <b>Erfahrungswerte</b> angewiesen, '
    'wobei Kundenw\u00fcnsche und eine gewisse Voreingenommenheit immer mitschwingen.',sB))
story.append(Paragraph(
    '<b>Alle angebotenen Stimmplatten funktionieren \u2014 auch die g\u00fcnstigsten.</b> '
    'Der Unterschied zwischen einer Economy- und einer Tipo-a-mano-Stimmplatte '
    'liegt nicht darin, ob sie funktioniert, sondern wie <b>gleichm\u00e4\u00dfig</b> '
    'sie \u00fcber den gesamten Tonbereich und Druckbereich funktioniert. '
    'Die teuerste Stimmplatte ist die, bei der der Stimmer am wenigsten '
    'nacharbeiten muss.',sK))

# Kap. 5
story.append(Paragraph('6. Was bestimmt die Qualit\u00e4t?',sCh))
story.append(Paragraph(
    'Die Qualit\u00e4t einer Stimmplatte h\u00e4ngt von Faktoren ab, die mit heutigen '
    'Mitteln <b>nicht in Tabellen erfasst werden k\u00f6nnen</b>:',sB))

rows_qual=[]
rows_qual.append(['Spaltweite','Dok. 0010','Enger = bessere Ansprache, aber empfindlicher (K\u00e4lte, Dok. 0018)'])
rows_qual.append(['Profilierung','Dok. 0016','Gleichm\u00e4\u00dfig \u00fcber den Tonbereich = gleichm\u00e4\u00dfige Ansprache'])
rows_qual.append(['Schlitzausformung','','Konisch vs. gerade. Binci: maschinell konisch'])
rows_qual.append(['Zungensteifigkeit','Dok. 0012','Bestimmt Ansprache und Druckabh\u00e4ngigkeit'])
rows_qual.append(['Materialqualit\u00e4t','','Stahlsorte, Homogenit\u00e4t, Eigenspannungen'])
rows_qual.append(['Oberfl\u00e4chenrauheit','Dok. 0010','Bei EDM (Harmonikas): Ra = 1,6 \u00b5m, hydraulisch glatt'])
rows_qual.append(['Nietqualit\u00e4t','','Fester Sitz, kein Spiel, keine D\u00e4mpfung'])
story.append(mk_tbl(['Faktor','Referenz','Auswirkung'],rows_qual,
    cw=[25*mm,18*mm,80*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Leider stellen die europ\u00e4ischen Stimmplattenhersteller '
    '<b>kaum Daten</b> bereit \u2014 weder Abmessungen noch Toleranzen noch '
    'Materialspezifikationen. Die Instrumentenbauer sind auf eigene Tests '
    'und auf ihre Erfahrungen angewiesen.',sW))

# Kap. 6
story.append(Paragraph('7. Harmonikas (Titlbach): EDM, Langschliff und DIX',sCh))
story.append(Paragraph(
    'Die Firma <b>Harmonikas</b> (Eigent\u00fcmer Titlbach, Tschechien) ist in mehrfacher '
    'Hinsicht ein Sonderfall \u2014 neben dem Standardsortiment bietet sie drei '
    'Verfahren an, die kein anderer Hersteller im Programm hat:',sB))

story.append(Paragraph('<b>6.1 Erodierte Stimmplatten (EDM)</b>',sCh))
story.append(Paragraph(
    'Funkenerosion (EDM) schneidet den Schlitz mit Funken statt mit Fr\u00e4sern. '
    'Extrem enge Spalte (wenige Mikrometer Toleranz), glatte Oberfl\u00e4che '
    '(Ra \u2248 1,6 \u00b5m, Dok. 0010), kein Grat, konische Schlitze ohne Aufwand. '
    'Nachteil: Engere Spalte = empfindlicher gegen K\u00e4lte-Kratzen '
    '(Dok. 0018: Paradoxon der Passgenauigkeit).',sB))

story.append(Paragraph('<b>6.2 DIX-Mensur (Nachbau)</b>',sCh))
story.append(Paragraph(
    'Harmonikas (Titlbach) bietet als einziger Nachbauten der historischen <b>DIX-Mensur</b> an. '
    'DIX-Stimmplatten galten als Referenz f\u00fcr Qualit\u00e4t und werden '
    'in historischen Instrumenten h\u00e4ufig angetroffen. Die Originale werden '
    'nicht mehr hergestellt.',sB))

story.append(Paragraph('<b>6.3 Langschliff-Zungen</b>',sCh))
story.append(Paragraph(
    'Bei konventionellen Stimmzungen wird die Oberfl\u00e4che quer zur L\u00e4ngsachse '
    'geschliffen. Dabei entstehen feine <b>Querrillen</b> auf der Zungenoberfl\u00e4che. '
    'Beim <b>Langschliff</b> werden die einzelnen Zungen kalt in <b>L\u00e4ngsrichtung</b> '
    'geschliffen. Das ist aufw\u00e4ndiger (jede Zunge einzeln), ergibt aber eine '
    'Oberfl\u00e4che ohne Querrillen.',sB))
story.append(Paragraph(
    '<b>Warum ist das relevant?</b> Querrillen sind potenzielle Sollbruchstellen \u2014 '
    'die Zunge schwingt in L\u00e4ngsrichtung, und Risse beginnen bevorzugt an Querrillen. '
    'Langschliff-Zungen haben L\u00e4ngsrillen, die parallel zur Schwingungsrichtung '
    'verlaufen und daher <b>keine Sollbruchstellen</b> bilden. '
    'Es ist anzunehmen, dass Langschliff-Zungen <b>haltbarer</b> sind \u2014 '
    'ein exakter Nachweis steht allerdings aus.',sB))

story.append(Paragraph('<b>Nastrino vs. Normal</b>',sCh))
story.append(Paragraph(
    'Harmonikas bietet Langschliff-Zungen in zwei Varianten an:',sB))
rows_lz=[]
rows_lz.append(['Langschliff Normal','Zungen aus Standardblech geschnitten, dann einzeln l\u00e4ngs geschliffen'])
rows_lz.append(['Langschliff Nastrino','Zungen aus <b>schmalem Band</b> gestanzt (ital. nastrino = B\u00e4ndchen). '
    'Die Breite des Bandes entspricht der Breite des <b>Zungenfu\u00dfes</b>. '
    'Die eigentliche Zungenbreite wird trotzdem gestanzt. '
    'Aufw\u00e4ndigerer Stanzvorgang, daf\u00fcr gleichm\u00e4\u00dfigere Materialstruktur.'])
story.append(mk_tbl(['Variante','Beschreibung'],rows_lz,cw=[32*mm,92*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Der Nastrino-Vorgang ist aufw\u00e4ndiger: Das schmale Band muss exakt '
    'auf die Breite des Zungenfu\u00dfes gewalzt sein, die Zungenform wird '
    'daraus gestanzt. Zusammen mit dem kalten L\u00e4ngsschliff entsteht eine Zunge '
    'mit minimal gestresster Materialstruktur.',sB))
story.append(Paragraph(
    '<b>Preis:</b> Ein Satz Nastrino-Langschliff kostet mindestens das '
    '<b>Doppelte</b> der besten professionellen Stimmplatte (Tipo a mano). '
    'F\u00fcr h\u00f6chste Anspr\u00fcche ist dieser Preis gerechtfertigt.',sB))

# Kap. 7: Nacharbeit
story.append(Paragraph('8. Nacharbeit: Was Harmonikabauer an Stimmplatten verbessern',sCh))
story.append(Paragraph(
    'Harmonikabauer arbeiten <b>nicht generell jeden Stimmsatz</b> nach \u2014 '
    'die folgenden Arbeitsg\u00e4nge (8.1\u20138.3) werden nur bei Bedarf und '
    'bei h\u00f6heren Anspr\u00fcchen durchgef\u00fchrt. '
    'Beim <b>Vorstimmen</b> wird aber auf jeden Fall die <b>Aufbiegung</b> '
    'eingestellt (8.4) \u2014 das ist Pflicht bei jedem Instrument.',sB))

story.append(Paragraph('<b>8.1 Profilierung nachschleifen</b>',sCh))
story.append(Paragraph(
    'Wird ein Stimmsatz nachgeschliffen, kauft man einen um <b>einen Halbton '
    'h\u00f6heren</b> Stimmsatz und schleift diesen auf die Solltonh\u00f6he herunter. '
    'Dadurch steht gen\u00fcgend Material f\u00fcr die Profilierung zur Verf\u00fcgung, '
    'ohne die Zunge zu d\u00fcnn werden zu lassen.',sB))
story.append(Paragraph(
    'Die Zungendicke wird durch Querschliff reduziert \u2014 meist im vorderen '
    'Bereich der Zunge. Das verringert die Steifigkeit und hat die in '
    'Dok. 0010\u20130012 und 0016 besprochenen Effekte: '
    'Die <b>Ansprache wird besser</b> (weniger Druck zum Anschwingen), '
    'aber die <b>maximal erzielbare Lautst\u00e4rke nimmt ab</b> (weichere Zunge '
    'gibt bei hohem Druck nach) und die <b>Tonhaltigkeit sinkt</b> '
    '(st\u00e4rkere Druckabh\u00e4ngigkeit der Frequenz, Dok. 0012). '
    'Die Profilierung ist immer ein Kompromiss zwischen Ansprache und '
    'Lautst\u00e4rkereserve.',sB))

story.append(Paragraph('<b>8.2 Spaltma\u00df minimieren</b>',sCh))
story.append(Paragraph(
    'Die Einpassung der Zunge wird nachjustiert und das Spaltma\u00df '
    'auf ein Minimum verringert (Dok. 0010). Das geschieht durch:',sB))
story.append(Paragraph(
    '\u2022 <b>Nach vorne Schlagen der Stimmzunge</b> und Abfeilen am beweglichen '
    'Ende \u2014 m\u00f6glich bei Zungen, die sich vom Niet zum freien Ende '
    '<b>verj\u00fcngen</b>. Durch das Schlagen wird die Zunge breiter im Schlitz '
    'positioniert, das Abfeilen stellt die Passform wieder her.<br/>'
    '\u2022 Bei <b>parallelen Stimmzungen</b> (keine Verj\u00fcngung) bleibt nur die '
    'M\u00f6glichkeit des nachtr\u00e4glichen <b>Verstemmen</b>s \u2014 '
    'leichtes Stauchen des Rahmens neben der Zunge.',sB))

story.append(Paragraph('<b>7.3 Symmetrisches Schwingen einstellen</b>',sCh))
story.append(Paragraph(
    'Die Zunge soll m\u00f6glichst <b>symmetrisch</b> schwingen, ohne seitliches '
    'Auslenken. Asymmetrisches Schwingen f\u00fchrt zu seitlichen '
    '<b>Torsionskratzern</b> (Dok. 0018). Das Einstellen erfolgt durch '
    'gezieltes Pressen an bestimmten Punkten der Zunge und im '
    'Zungenfu\u00dfbereich \u2014 reine Erfahrungssache, mit feinen Dornen, '
    'eingespannt in Handpressen. Ziel: Die Zunge schwingt exakt mittig '
    'durch den Schlitz, ohne seitliche Drift.',sB))

story.append(Paragraph('<b>8.4 Aufbiegung einstellen (immer)</b>',sCh))
story.append(Paragraph(
    'Die Aufbiegung wird <b>bei jedem Instrument</b> beim Vorstimmen eingestellt \u2014 '
    'das ist kein optionaler Arbeitsgang. '
    'Sie bestimmt Ansprache, Dynamikbereich und '
    'Klangcharakter (Dok. 0010, 0019). Grundregeln:',sB))
story.append(Paragraph(
    '\u2022 Die Aufbiegung sollte <b>eher im letzten freien Drittel</b> der Zunge '
    'erfolgen \u2014 nicht gleichm\u00e4\u00dfig \u00fcber die gesamte L\u00e4nge.<br/>'
    '\u2022 Es schadet nicht, wenn der Bereich <b>vom Niet bis zum letzten Drittel</b> '
    'der Zunge etwas in den Kanal eintaucht \u2014 oft verbessert das die Ansprache. '
    'Nur das letzte freie Drittel steht nach oben auf.<br/>'
    '\u2022 Die Form der Aufbiegung hat auch einen <b>klanglichen Aspekt</b>: '
    'Eine Zunge, die nur am Ende aufgebogen ist, klingt anders als eine, '
    'die gleichm\u00e4\u00dfig gebogen ist.',sB))

story.append(Paragraph('<b>7.5 Lein\u00f6l-Schutzfilm</b>',sCh))
story.append(Paragraph(
    'Als letzter Arbeitsgang sollte ein <b>sehr d\u00fcnner Film aus Lein\u00f6l</b> '
    'auf die Stimmplatten aufgetragen werden. Lein\u00f6l polymerisiert an der Luft '
    'zu einer harten, d\u00fcnnen Schutzschicht. Es sch\u00fctzt vor Korrosion '
    '(Feuchtigkeit aus der Atemluft bei Mundharmonikas, Kondenswasser im Balg) '
    'und verringert die Reibung der Zunge im Schlitz. '
    'Wichtig: Nur ein <b>hauchdunner</b> Film \u2014 zu viel \u00d6l verklebt die Zunge '
    'oder ver\u00e4ndert die Masse und damit die Stimmung.',sB))

story.append(Paragraph('<b>7.6 Ventilierung</b>',sCh))
story.append(Paragraph(
    'Die Ventilierung ist eine der <b>wichtigsten</b> Arbeiten. Das Ventil verschlie\u00dft '
    'den Schlitz der gegen\u00fcberliegenden Zunge und verhindert Luftverlust. '
    'Je nach Zungenl\u00e4nge und Tonh\u00f6he muss das Ventil angepasst werden \u2014 '
    'es soll immer <b>sauber schlie\u00dfen</b>, aber <b>nicht zu steif</b> sein. '
    'Ein zu steifes Ventil d\u00e4mpft die Ansprache, ein zu weiches schlie\u00dft '
    'nicht sauber. Immer ein Kompromiss.',sB))

story.append(Paragraph(
    '<b>Eigenresonanz:</b> Ventile besitzen eine eigene Resonanzfrequenz. '
    'Diese darf <b>auf keinen Fall mit dem Ton zusammenfallen</b> \u2014 '
    'sonst entsteht eine Schnarr-Modulation des Tones. '
    'Die Schnarrt\u00f6ne sind oft auch Untert\u00f6ne. '
    'Die Ventill\u00e4nge und -steifigkeit m\u00fcssen so gew\u00e4hlt werden, '
    'dass die Eigenresonanz au\u00dferhalb des Tonbereichs der zugeordneten Zunge liegt.',sW))

story.append(Paragraph(
    '<b>Aufschlagger\u00e4usche:</b> Das Ventil soll beim Schlie\u00dfen keine '
    'h\u00f6rbaren Aufschlagger\u00e4usche erzeugen. Daher verwendet man heute im '
    'unteren Bereich (lange Ventile, tiefe T\u00f6ne) oft <b>kombinierte Ventile</b>:',sB))

rows_vent=[]
rows_vent.append(['Leder + Vinyl','Mehrere Lagen','Tiefe T\u00f6ne. Weich, ger\u00e4uscharm, gute D\u00e4mpfung.'])
rows_vent.append(['Leder + Metallfeder','Verst\u00e4rkt','Tiefe\u2013mittlere T\u00f6ne. Kontrollierter R\u00fcckschlag.'])
rows_vent.append(['Vinyl (d\u00fcnn)','Einlagig','Mittlere\u2013hohe T\u00f6ne. Leicht, schnelle Reaktion.'])
rows_vent.append(['Kein Ventil','Entf\u00e4llt','H\u00f6chste T\u00f6ne (kurze Zungen). Luftverlust vertretbar gering.'])
story.append(mk_tbl(['Material','Bauform','Einsatz'],rows_vent,
    cw=[28*mm,22*mm,68*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Bei den <b>kleinsten Stimmplatten</b> (k\u00fcrzeste Zungen, h\u00f6chste T\u00f6ne) '
    'werden Ventile <b>weggelassen</b>, da der Luftdurchsatz bei der Zunge in '
    'Gegenrichtung bereits so gering ist, dass der Verlust vertretbar wird.',sB))

story.append(Paragraph(
    '<b>Verstimmung durch Ventile:</b> Das Ventil beeinflusst die Frequenz '
    'der Stimmzunge \u2014 es muss daher <b>mit Ventilen gestimmt</b> werden. '
    'Ohne Ventil ist die Frequenz h\u00f6her (einige Cent, abh\u00e4ngig von der Tonlage). '
    'Wer erst stimmt und dann ventiliert, stimmt falsch.',sW))

story.append(Paragraph('<b>Kleber:</b>',sB))
story.append(Paragraph(
    '\u2022 <b>Pattex (Kontaktkleber):</b> Heute am h\u00e4ufigsten verwendet. '
    '<b>Achtung:</b> Bei sehr d\u00fcnnen Vinyl-Ventilen darf Pattex <b>nicht</b> '
    'verwendet werden \u2014 das enthaltene L\u00f6sungsmittel verformt den d\u00fcnnen '
    'Vinyl-Film bereits w\u00e4hrend des Austrocknens. Die Verformung wird erst '
    'nach Tagen sichtbar.<br/>'
    '\u2022 <b>Schellack (eingedickt):</b> Fr\u00fcher \u00fcblich. Guter Halt, '
    'aber aufw\u00e4ndiger in der Verarbeitung.<br/>'
    '\u2022 <b>L\u00f6sungsmittelfreie Kleber:</b> Moderne Alternative, verformt '
    'auch d\u00fcnnes Vinyl nicht. Etwas schwieriger zu verarbeiten '
    '(l\u00e4ngere Abl\u00fcftzeit, andere Konsistenz).',sB))

rows_na=[]
rows_na.append(['Profilierung nachschleifen','Ansprache \u2191, Lautst\u00e4rke \u2193, Tonhaltigkeit \u2193','Kompromiss Ansprache/Reserve'])
rows_na.append(['Spaltma\u00df minimieren','Ansprache \u2191, K\u00e4lteempfindlichkeit \u2191','Dok. 0010, 0018'])
rows_na.append(['Symmetrie einstellen','Kein Torsionskratzen','Erfahrungssache, Handpressen'])
rows_na.append(['Aufbiegung optimieren','Ansprache, Dynamik, Klang','Eher am freien Ende'])
rows_na.append(['Lein\u00f6l-Schutzfilm','Korrosionsschutz, weniger Reibung','Hauchdünn, letzter Arbeitsgang'])
rows_na.append(['Ventilierung','Sauberes Schlie\u00dfen, keine Schnarrt\u00f6ne','Material + L\u00e4nge an Ton anpassen'])
story.append(mk_tbl(['Arbeitsgang','Wirkung','Anmerkung'],rows_na,
    cw=[32*mm,40*mm,42*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Alle diese Arbeitsg\u00e4nge erfordern viel Erfahrung und sind '
    '<b>nicht standardisierbar</b>. Jede Zunge verh\u00e4lt sich etwas anders. '
    'Die Nacharbeit ist der Grund, warum ein gut eingestelltes Instrument '
    'besser klingt als die Summe seiner Teile \u2014 und warum der Harmonikabauer '
    'mindestens so wichtig ist wie die Stimmplatte.',sK))

# Kap. 8
story.append(Paragraph('9. Zusammenfassung',sCh))
pts=[
    '<b>F\u00fcnf Hersteller:</b> Harmonikas, AV (Antonelli/Salpa), Artigianale, Cagnoni, Binci '
    '\u2014 alle in Italien. Harmonikas (Titlbach, Tschechien) als Sonderfall: EDM, DIX-Nachbau, Langschliff.',
    '<b>Mensuren sind sehr \u00e4hnlich</b> bei Akkordeon-Stimmplatten. '
    'Zungenbreiten variieren: Binci, Cagnoni, Harmonikas eher breiter. Harmonikas: Kundenmensur.',
    '<b>G\u00fcteklassen:</b> Economy bis Tipo a mano. Alle funktionieren. '
    'Der Unterschied liegt in der Gleichm\u00e4\u00dfigkeit \u00fcber den Tonbereich.',
    '<b>Die Qualit\u00e4tsunterschiede sind minimal</b> und nicht exakt messbar. '
    'Kein objektives Kategorisierungsverfahren. Erfahrungswerte, Kundenw\u00fcnsche '
    'und Voreingenommenheit schwingen immer mit.',
    '<b>Was tats\u00e4chlich z\u00e4hlt:</b> Spaltweite, Profilierung, Schlitzform, '
    'Materialqualit\u00e4t \u2014 aber die Hersteller liefern kaum Daten.',
    '<b>Alle angebotenen Stimmplatten funktionieren</b> \u2014 auch die g\u00fcnstigsten. '
    'Die teuerste ist die, bei der der Stimmer am wenigsten nacharbeiten muss.',
    '<b>Langschliff (Harmonikas):</b> L\u00e4ngsschliff statt Querschliff \u2192 keine Querrillen '
    '\u2192 vermutlich haltbarer. Nastrino: aus schmalem Band gestanzt. '
    'Preis: mindestens doppelt so teuer wie Tipo a mano.',
    '<b>Nacharbeit des Harmonikabauers</b> ist entscheidend: Profilierung nachschleifen, '
    'Spaltma\u00df minimieren, symmetrisches Schwingen einstellen, Aufbiegung optimieren. '
    'Erfahrungssache \u2014 der Baumeister ist mindestens so wichtig wie die Platte.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Die Stimmplatte ist das Herz des Instruments \u2014 '
    'aber alle Herzen schlagen. Die Kunst liegt im Stimmen, nicht in der Platte.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0021 \u2014 Stimmplatten \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
