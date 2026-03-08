#!/usr/bin/env python3
"""Dok. 0022 — Balg: Querschnitt, Falten, Instrumentengröße"""
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
def fig2img(fig,w=160,dpi=200):
    buf=BytesIO(); fig.savefig(buf,format='png',dpi=dpi,bbox_inches='tight'); plt.close(fig)
    buf.seek(0); img=Image(buf); r=img.imageWidth/img.imageHeight
    img.drawWidth=w*mm; img.drawHeight=w*mm/r; return img
plt.rcParams['font.family']='DejaVu Sans'; plt.rcParams['font.size']=10

import numpy as np

# Instrumente: (Name, Breite mm, Höhe mm, Faltentiefe mm, Anzahl Falten)
instr=[
    ("Kleine Steirische (46)", 250, 180, 15, 18),
    ("Standard Steirische (48\u201350)", 290, 200, 17, 20),
    ("Gro\u00dfe Steirische (60)", 320, 220, 19, 22),
    ("Knopfakkordeon 72 Bass", 280, 190, 16, 19),
    ("Pianoakkordeon 96 Bass", 290, 190, 17, 20),
    ("Pianoakkordeon 120 Bass", 300, 190, 18, 22),
]

def calc(b,h,d,n):
    A=b*h/100; hub=d*n*0.75/10; V=A*hub/1000
    return A,hub,V

# ══ DIAGRAMME ══
def fig_volumen():
    fig,ax=plt.subplots(figsize=(10,5))
    dims=[(20,34,13),(20,34,18),(20,36,15),(20,38,15),(20,38,18),(20,40,18)]
    names=[f'{b}\u00d7{h}, {n}F' for b,h,n in dims]
    vols=[b*h*3.5*n*0.75/1000 for b,h,n in dims]
    hubs=[3.5*n*0.75 for _,_,n in dims]
    x=np.arange(len(names)); w=0.35
    ax.bar(x-w/2,vols,w,color='#1565c0',label='Nutzbares Volumen [l]')
    ax2=ax.twinx()
    ax2.bar(x+w/2,hubs,w,color='#e94560',alpha=0.7,label='Nutzbarer Hub [cm]')
    ax.set_xticks(x); ax.set_xticklabels(names,fontsize=8,rotation=15,ha='right')
    ax.set_ylabel('Volumen [Liter]',color='#1565c0'); ax2.set_ylabel('Hub [cm]',color='#e94560')
    ax.legend(loc='upper left',fontsize=8); ax2.legend(loc='upper right',fontsize=8)
    ax.set_title('Abb. 1: Nutzbares Balgvolumen (Faltentiefe 35 mm)')
    ax.grid(True,alpha=0.3,axis='y'); fig.tight_layout(); return fig

def fig_druck():
    fig,ax=plt.subplots(figsize=(10,5))
    F_vals=np.linspace(5,40,100)
    for nm,b,h in [('20\u00d734 cm',200,340),('20\u00d736 cm',200,360),
                    ('20\u00d738 cm',200,380),('20\u00d740 cm',200,400)]:
        A_m2=b*h*1e-6
        p=F_vals/A_m2
        ax.plot(F_vals,p,lw=2,label=nm)
    ax.set_xlabel('Armkraft [N]'); ax.set_ylabel('Balgdruck [Pa]')
    ax.set_title('Abb. 2: Balgdruck vs. Armkraft (vierreihig diatonisch)')
    ax.legend(fontsize=8); ax.grid(True,alpha=0.3)
    ax.axhspan(200,600,alpha=0.05,color='green')
    ax.text(8,400,'Typischer Spieldruckbereich',fontsize=8,color='green')
    fig.tight_layout(); return fig

# ══ PDF ══
outfile='0022_balg_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0022',sST))
story.append(Paragraph('Balg:<br/>Querschnitt, Faltenzahl und Instrumentengr\u00f6\u00dfe',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('Zusammenhang zwischen Balgquerschnitt, Faltenzahl, nutzbarem Volumen, '
    'Hub und Instrumentengr\u00f6\u00dfe. Luftverbrauch und Spielzeit. '
    'Druck und Armkraft. Referenz: Dok. 0010\u20130012.',sAb))

# Kap. 1
story.append(Paragraph('1. Instrumentengr\u00f6\u00dfe: Was bestimmt die Abmessungen?',sCh))
story.append(Paragraph(
    'Die Instrumentengr\u00f6\u00dfe wird von der <b>Menge an Einbauten</b> bestimmt \u2014 '
    'Stimmst\u00f6cke, Mechanik, Lager, Register. Heutige vierreihige Instrumente '
    'sind eher <b>gr\u00f6\u00dfer</b> als traditionelle, obwohl bemerkenswerterweise '
    'wenig wirklich Neues erfunden wurde \u2014 Verkaufsprospekte vor dem Ersten '
    'Weltkrieg zeigen bereits die gleichen Grundkonstruktionen.',sB))
story.append(Paragraph(
    '<b>Was sich verbessert hat:</b> Die Mechanik arbeitet heute pr\u00e4ziser und '
    'ger\u00e4usch\u00e4rmer. Jedoch wurde durch die verbesserten Lager auch der '
    '<b>Platzbedarf h\u00f6her</b>. Stimmplatten im Diskantbereich bei den tieferen '
    'T\u00f6nen wurden dicker. Stimmst\u00f6cke m\u00fcssen heute schrag gestellt werden, '
    'da die tiefen Zungen weiter ausschwingen als bei Instrumenten um 1900. '
    'Hebel wurden verst\u00e4rkt, damit die Eigenschwingungen geringer werden.',sB))
story.append(Paragraph(
    'Man findet daher durchaus <b>traditionelle historische Instrumente</b> '
    'mit kleineren Balgabmessungen als heute \u00fcblich. '
    'Die gr\u00f6\u00dfere Bauh\u00f6he im geschlossenen Zustand reduziert den '
    'verf\u00fcgbaren Bewegungsbereich des Balges.',sB))
story.append(Paragraph(
    '<b>Anmerkung:</b> Dieses Dokument behandelt ausschlie\u00dflich '
    '<b>vierreihige Instrumente</b>, da diese die \u00fcbliche Variante sind. '
    'Drei- oder f\u00fcnfreihige Instrumente ver\u00e4ndern die Situation wesentlich '
    '\u2014 insbesondere die Balgf\u00fchrung.',sW))

# Kap. 2
story.append(Paragraph('2. Balgabmessungen vierreihiger diatonischer Instrumente',sCh))
story.append(Paragraph(
    '\u00dcbliche Balgabmessungen (Breite \u00d7 H\u00f6he, Au\u00dfenma\u00df):',sB))

rows_inst=[]
rows_inst.append(['20 \u00d7 34','Kompakt, traditionell','Weniger Volumen, k\u00fcrzere Phrasen'])
rows_inst.append(['20 \u00d7 36','H\u00e4ufig','Guter Kompromiss'])
rows_inst.append(['20 \u00d7 38','Weit verbreitet','Standard bei modernen Instrumenten'])
rows_inst.append(['20 \u00d7 40','Gro\u00df','Viel Volumen, l\u00e4ngere Phrasen'])
story.append(mk_tbl(['B \u00d7 H [cm]','Einordnung','Eigenschaft'],rows_inst,
    cw=[22*mm,30*mm,60*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die Breite kann etwas mehr als 20 cm betragen (20,5 oder 21 cm). '
    'Die Faltentiefe ist bei den meisten Instrumenten <b>35 mm</b>. '
    'Wenn die Einbauten weniger Tiefe erzwingen, muss im Gegenzug '
    'die <b>Faltenzahl erh\u00f6ht</b> werden.',sB))
story.append(Paragraph(
    '<b>Geschlossene H\u00f6he:</b> 18 geschlossene Falten ergeben <b>65 mm</b> '
    'H\u00f6he (\u2248 3,6 mm pro Falte). '
    'Verwendet man 14 Falten, reduziert sich diese H\u00f6he auf \u2248 50 mm \u2014 '
    'das sind <b>15 mm weniger</b>, die direkt als zus\u00e4tzlicher Bewegungsbereich '
    'zur Verf\u00fcgung stehen. Weniger Falten = weniger geschlossene H\u00f6he = '
    'mehr nutzbarer Hub bei gleicher Arml\u00e4nge.',sB))

# Kap. 3
story.append(Paragraph('3. Faltenzahl, Hub und Bewegungsbereich',sCh))
story.append(Paragraph(
    '<b>Nicht die Faltenzahl ist ma\u00dfgeblich</b> \u2014 sondern der Bereich, '
    'in dem der Balg vom geschlossenen Zustand bis zur maximalen L\u00e4nge '
    'variiert werden kann. Bei tiefen Falten (35 mm) und einem beweglichen Balg '
    'reichen oft <b>weniger Falten</b> als bei flachen Falten.',sB))
story.append(Paragraph(
    'Dabei spielt die <b>Arml\u00e4nge</b> eine Rolle: Ein gro\u00dfer Mensch kann '
    'nat\u00fcrlich mehr Falten ausnutzen als ein Kind. '
    'Die <b>H\u00f6he des Instruments im geschlossenen Zustand</b> (Balg + Diskant + Bass) '
    'sollte minimal sein. '
    'Die geschlossene Balgh\u00f6he allein betr\u00e4gt 65 mm bei 18 Falten '
    'bzw. 50 mm bei 14 Falten \u2014 dazu kommen Diskant- und Basskorpus. '
    'Jeder Millimeter weniger geschlossene H\u00f6he ist ein Millimeter mehr Hub.',sB))

# Volumen-Tabelle mit realen Maßen
story.append(Paragraph(
    '<b>Nutzbares Volumen:</b> V = A \u00d7 nutzbarer Hub. '
    'Beispielrechnung mit 35 mm Faltentiefe, 75 % Nutzung. '
    'Faltenzahl: 13\u201318 (Maximum 18).',sB))
rows_vol=[]
for bh, b_cm, h_cm in [('20\u00d734',20,34),('20\u00d736',20,36),('20\u00d738',20,38),('20\u00d740',20,40)]:
    A=b_cm*h_cm
    for nf in [13,15,18]:
        hub_cm=3.5*nf*0.75
        V=A*hub_cm/1000
        rows_vol.append([bh,str(A),str(nf),f'{hub_cm:.0f}',f'{V:.1f}'])
story.append(mk_tbl(['B\u00d7H [cm]','A [cm\u00b2]','Falten','Hub [cm]','V [Liter]'],rows_vol,
    cw=[18*mm,16*mm,14*mm,16*mm,16*mm]))
story.append(Spacer(1,3*mm))

fig1=fig_volumen(); story.append(fig2img(fig1,150)); story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Faustregel:</b> Pro Falte bei 35 mm Tiefe \u2248 26 mm nutzbarer Hub '
    '(35 \u00d7 0,75). 18 Falten \u00d7 26 mm \u2248 47 cm nutzbarer Hub.',sB))

# Kap. 3
story.append(Paragraph('4. Luftverbrauch und Spielzeit',sCh))
story.append(Paragraph(
    'Der Luftverbrauch h\u00e4ngt ab von: Anzahl gleichzeitig klingender Zungen, '
    'Tonh\u00f6he (tiefe T\u00f6ne = mehr Luft), Lautst\u00e4rke '
    '(mehr Druck = mehr Durchfluss) und Registrierung (mehr Ch\u00f6re = mehr Verbrauch).',sB))

rows_verb=[]
rows_verb.append(['Tiefer Bass (A1, 55 Hz)','3\u20135','Gro\u00dfe Amplitude, viel Luft'])
rows_verb.append(['Mittlerer Ton (A3, 220 Hz)','1\u20132','Moderate Amplitude'])
rows_verb.append(['Hoher Diskant (A5, 880 Hz)','0,3\u20130,8','Kleine Amplitude, wenig Luft'])
story.append(mk_tbl(['Tonlage','l/min (eine Zunge)','Anmerkung'],rows_verb,cw=[32*mm,28*mm,52*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Bei <b>Akkordspiel</b> (Bass + Akkord = 6\u201312 Zungen gleichzeitig): '
    'Verbrauch \u2248 5\u201315 l/min je nach Dynamik und Registrierung. '
    'Bei 10 Liter nutzbarem Volumen und 10 l/min Verbrauch: '
    '<b>Spielzeit pro Balgzug \u2248 60 Sekunden.</b> '
    'Bei lautem Bassspiel mit vielen Ch\u00f6ren deutlich k\u00fcrzer.',sB))
story.append(Paragraph(
    '<b>Konsequenz f\u00fcr die Instrumentengr\u00f6\u00dfe:</b> '
    'Wer \u00fcberwiegend Bass und Akkordbegleitung spielt (Begleiter), '
    'braucht mehr Volumen als ein Melodiespieler. Wer laut spielt, '
    'braucht mehr Volumen als ein leiser Spieler. '
    'Die Instrumentengr\u00f6\u00dfe sollte zum Spielstil passen.',sK))

# Kap. 4
story.append(Paragraph('5. Querschnitt und Balgdruck',sCh))
story.append(Paragraph(
    'Der Balgdruck ergibt sich aus der Armkraft geteilt durch die Querschnittsfl\u00e4che:',sB))
story.append(Paragraph('p = F / A',sF))
story.append(Paragraph(
    'Gr\u00f6\u00dferer Querschnitt bei gleicher Armkraft = <b>weniger Druck</b>. '
    'Das hat zwei Seiten: F\u00fcr leises Spiel braucht man weniger Kraft, '
    'aber f\u00fcr maximale Lautst\u00e4rke reicht die Armkraft m\u00f6glicherweise nicht. '
    'Kleinere B\u00e4lge erzeugen bei gleicher Kraft mehr Druck.',sB))

fig2=fig_druck(); story.append(fig2img(fig2,150)); story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Typischer Spieldruckbereich: 200\u2013600 Pa (Dok. 0012). '
    'Die Querschnitte reichen von 20\u00d734 = 680 cm\u00b2 bis 20\u00d740 = 800 cm\u00b2 '
    '\u2014 das sind 18 % Unterschied in der Fl\u00e4che. '
    'Bei 20 N Armkraft: 680 cm\u00b2 \u2192 294 Pa, 800 cm\u00b2 \u2192 250 Pa. '
    'Der gr\u00f6\u00dfere Querschnitt braucht <b>mehr Kraft</b> f\u00fcr denselben Druck.',sB))

# Kap. 5
story.append(Paragraph('6. Faltentiefe und Material',sCh))
story.append(Paragraph(
    'Die Faltentiefe (Abstand Berg\u2013Tal) bestimmt das Volumen pro Falte. '
    'Bei den meisten vierreihigen diatonischen Instrumenten: <b>35 mm</b>. '
    'Einige haben weniger Tiefe, wenn die Einbauten das erzwingen \u2014 '
    'dann muss die Faltenzahl erh\u00f6ht werden, um dasselbe Volumen zu erhalten.',sB))

story.append(Paragraph('<b>Vergleich: 35 mm vs. 28 mm bei gleichem \u00d6ffnungswinkel</b>',sCh))
story.append(Paragraph(
    'Die Gesamtl\u00e4nge des ge\u00f6ffneten Balges muss in die verf\u00fcgbare '
    'Arml\u00e4nge passen. Pro Falte verbraucht man die geschlossene H\u00f6he '
    '(\u22483,6 mm) plus den Hub (2d\u00b7sin(\u03b8/2)):',sB))
story.append(Paragraph(
    'n_max = R_eff / [h_geschl + 2d\u00b7sin(\u03b8/2)]',sF))
story.append(Paragraph(
    'Bei gleichem \u00d6ffnungswinkel passen <b>mehr flache Falten</b> in denselben Raum '
    '(jede Falte verbraucht weniger Platz). Aber jede Falte liefert auch '
    '<b>weniger Hub</b>. Der Netto-Effekt:',sB))

rows_fv=[]
import math
A=20*38; R_eff=40; h_c_pf=6.5/18
for theta_deg in [40,50,60,70,80]:
    th=math.radians(theta_deg)
    ext35=2*3.5*math.sin(th/2); ext28=2*2.8*math.sin(th/2)
    n35=min(int(R_eff/(h_c_pf+ext35)),18); n28=min(int(R_eff/(h_c_pf+ext28)),18)
    hub35=n35*ext35; hub28=n28*ext28
    V35=A*hub35/1000; V28=A*hub28/1000
    dV=(V28/V35-1)*100 if V35>0 else 0
    rows_fv.append([f'{theta_deg}\u00b0',str(n35),f'{hub35:.1f}',f'{V35:.1f}',
        str(n28),f'{hub28:.1f}',f'{V28:.1f}',f'{dV:+.0f} %'])
story.append(mk_tbl(['\u03b8','n (35)','Hub 35','V 35 [l]','n (28)','Hub 28','V 28 [l]','\u0394V'],rows_fv,
    cw=[12*mm,12*mm,14*mm,14*mm,12*mm,14*mm,14*mm,14*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Ergebnis:</b> 35 mm Falten gewinnen bei typischen Betriebswinkeln '
    '(40\u201370\u00b0) ca. <b>3\u20137 %</b> Volumen gegen\u00fcber 28 mm \u2014 '
    'trotz weniger Falten. Bei 28 mm kompensieren die zus\u00e4tzlichen 2\u20133 Falten '
    'den geringeren Hub pro Falte teilweise, aber nicht vollst\u00e4ndig. '
    'Der Vorteil tiefer Falten ist <b>nicht dramatisch</b>, aber real: '
    'Weniger Falten, mehr Hub, weniger Verschlei\u00df, weicherer Balggang.',sK))

story.append(Paragraph('<b>Braucht man mehr als 14 Falten bei 35 mm?</b>',sCh))
story.append(Paragraph(
    'Jede zus\u00e4tzliche Falte verbraucht 3,6 mm geschlossene H\u00f6he \u2014 '
    'das wird vom verf\u00fcgbaren Hub abgezogen. <b>Mehr Falten = weniger Hub '
    '= weniger Volumen:</b>',sB))

rows_14=[]
import math
R_eff=40; hc=6.5/18; A=20*38
for n in [13,14,15,16,18]:
    hub=R_eff-n*hc; V=A*hub/1000
    sin_h=(R_eff/n-hc)/(2*3.5); th=2*math.degrees(math.asin(min(sin_h,1)))
    gain=(V/(A*(R_eff-14*hc)/1000)-1)*100
    rows_14.append([str(n),f'{n*hc:.1f}',f'{hub:.1f}',f'{V:.1f}',f'{th:.0f}\u00b0',f'{gain:+.1f} %'])
story.append(mk_tbl(['Falten','Geschl. [cm]','Hub [cm]','V [l]','\u03b8 (35 mm)','vs. 14 Falten'],rows_14,
    cw=[14*mm,20*mm,16*mm,14*mm,16*mm,20*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '14 Falten \u00d7 35 mm: Hub = 34,9 cm, V = 26,6 l, \u03b8 = 42\u00b0 (moderate \u00d6ffnung). '
    '18 Falten \u00d7 35 mm: Hub = 33,5 cm, V = 25,5 l, \u03b8 = 31\u00b0 \u2014 '
    '<b>1,1 Liter weniger</b> bei 4 Falten mehr.',sB))
story.append(Paragraph(
    '<b>M\u00f6gliche Gr\u00fcnde f\u00fcr mehr als 14 Falten:</b> '
    'Weicherer Balggang (jede Falte \u00f6ffnet weniger \u2192 weniger Widerstand), '
    'gleichm\u00e4\u00dfigerer Druck bei Richtungswechsel, kurze Arme (Kind). '
    'Aber physikalisch gewinnt man mit 14 Falten bei 35 mm <b>mehr Hub und '
    'mehr Volumen</b> als mit 18 Falten. Historisch waren weniger Falten '
    'ausreichend und \u00fcblich \u2014 der Trend zu mehr Falten ist eine '
    '<b>moderne Tendenz</b>, die physikalisch nicht begr\u00fcndet ist.',sK))
story.append(Paragraph(
    'Das Balgmaterial (Karton, Leinen, Leder, Eckverst\u00e4rkungen) bestimmt '
    'die <b>Luftdichtheit</b> und die <b>Lebensdauer</b>. Ein undichter Balg '
    'verliert Luft und reduziert die effektive Spielzeit. '
    'Die Ecken sind die kritischen Stellen \u2014 dort knickt das Material '
    'bei jeder Falte.',sB))

story.append(Paragraph(
    'Ein Balg kann <b>nie absolut dicht</b> sein. Es werden Dichtungen bei den '
    'Klappen und am Balg selbst verwendet. Die Dichtheit des Balges variiert '
    'je nach <b>Leder</b>, welches in den beweglichen Eins\u00e4tzen bei den Ecken '
    'verwendet wird. D\u00fcnneres, weicheres Leder ist beweglicher '
    '(weniger Balgwiderstand), aber weniger dicht. Dickeres Leder ist dichter, '
    'aber steifer. Auch hier ein Kompromiss.',sB))
story.append(Paragraph(
    'Undichtigkeiten reduzieren den <b>effektiven Druck</b> und das '
    '<b>nutzbare Volumen</b>. Bei einem typischen Balg geht ein Teil des '
    'Luftvolumens durch Leckagen verloren, bevor er die Stimmzungen erreicht. '
    'Je langsamer der Balgwechsel, desto gr\u00f6\u00dfer der relative Verlust.',sB))

# Diagramm: Druck vs. Querschnitt
fig_pq,ax_pq=plt.subplots(figsize=(10,5))
A_vals=np.linspace(600,850,100)  # cm², 20×30 bis 21×40
for F,col,ls in [(10,'#2e7d32','--'),(15,'#1565c0','-'),(20,'#ff9800','-'),(30,'#e94560','-')]:
    p=F/(A_vals*1e-4)  # Pa
    ax_pq.plot(A_vals,p,ls,color=col,lw=2,label=f'F = {F} N')
ax_pq.set_xlabel('Balgquerschnitt [cm\u00b2]')
ax_pq.set_ylabel('Balgdruck [Pa]')
ax_pq.set_title('Abb. 3: Balgdruck vs. Querschnitt bei verschiedenen Armkr\u00e4ften')
ax_pq.legend(fontsize=9); ax_pq.grid(True,alpha=0.3)
ax_pq.axhspan(200,600,alpha=0.05,color='green')
ax_pq.text(620,400,'Spieldruckbereich',fontsize=8,color='green')
# Markierungen für reale Querschnitte
for bh,A_real in [('20\u00d734',680),('20\u00d736',720),('20\u00d738',760),('20\u00d740',800)]:
    ax_pq.axvline(x=A_real,color='gray',ls=':',lw=0.8,alpha=0.5)
    ax_pq.text(A_real,ax_pq.get_ylim()[1]*0.95,bh,fontsize=7,ha='center',color='gray')
fig_pq.tight_layout()
story.append(fig2img(fig_pq,150)); story.append(Spacer(1,3*mm))

# Diagramm: Volumen vs. Balgweg für verschiedene Querschnitte
fig_vw,(ax_v1,ax_v2)=plt.subplots(1,2,figsize=(12,5))
# Links: V vs. Balgweg
weg=np.linspace(0,50,200)  # cm Balgweg
for bh,b,h,col in [('20\u00d734',20,34,'#2e7d32'),('20\u00d736',20,36,'#1565c0'),
                     ('20\u00d738',20,38,'#ff9800'),('20\u00d740',20,40,'#e94560')]:
    A=b*h  # cm²
    V=A*weg/1000  # Liter
    ax_v1.plot(weg,V,'-',color=col,lw=2,label=f'{bh} cm ({A} cm\u00b2)')
ax_v1.set_xlabel('Balgweg [cm]')
ax_v1.set_ylabel('Volumen [Liter]')
ax_v1.set_title('Volumen vs. Balgweg')
ax_v1.legend(fontsize=8); ax_v1.grid(True,alpha=0.3)
ax_v1.axvspan(25,45,alpha=0.05,color='blue')
ax_v1.text(35,1,'Typischer\nHubbereich',fontsize=8,color='blue',ha='center')

# Rechts: Vergleich mit Leckage
weg2=np.linspace(0,50,200)
A_38=20*38  # cm²
V_ideal=A_38*weg2/1000
# Leckage: V_eff = V_ideal - Leckrate × Zeit
# Zeit ≈ Weg / Geschwindigkeit, Geschwindigkeit ≈ 20 cm/s (schnell) oder 5 cm/s (langsam)
for v_balg,name,ls in [(20,'Schneller Balgwechsel (20 cm/s)','-'),
                        (10,'Moderater Balgwechsel (10 cm/s)','--'),
                        (5,'Langsamer Balgwechsel (5 cm/s)',':')]:
    leck_rate=0.3  # l/s typische Leckage
    zeit=weg2/(v_balg+0.001)  # s
    V_eff=V_ideal - leck_rate*zeit
    V_eff=np.maximum(V_eff,0)
    ax_v2.plot(weg2,V_eff,ls,color='#e94560',lw=2,label=name)
ax_v2.plot(weg2,V_ideal,'-',color='#1565c0',lw=2,label='Ideal (dicht)')
ax_v2.set_xlabel('Balgweg [cm]')
ax_v2.set_ylabel('Nutzbares Volumen [Liter]')
ax_v2.set_title('Effektives Volumen mit Leckage (20\u00d738)')
ax_v2.legend(fontsize=7); ax_v2.grid(True,alpha=0.3)
fig_vw.suptitle('Abb. 4: Luftvolumen vs. Balgweg',fontweight='bold')
fig_vw.tight_layout(); 
story.append(fig2img(fig_vw,158)); story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Abb. 4 links zeigt den linearen Zusammenhang zwischen Balgweg und Volumen '
    'f\u00fcr verschiedene Querschnitte. Der Unterschied zwischen 20\u00d734 und 20\u00d740 '
    'betr\u00e4gt bei 40 cm Hub: {:.1f} vs. {:.1f} Liter \u2014 <b>{:.0f} % mehr</b>. '
    'Abb. 4 rechts zeigt den Einfluss der Leckage: Bei langsamen Balgwechseln '
    'geht ein gr\u00f6\u00dferer Anteil des Volumens verloren. '
    'Ein dichter Balg ist bei <b>langsamen Phrasen</b> wichtiger als bei schnellem Spiel.'.format(
        20*34*40/1000, 20*40*40/1000, (40-34)/34*100),sK))

# Kap. 6
story.append(Paragraph('7. Zusammenfassung',sCh))
pts=[
    '<b>Instrumentengr\u00f6\u00dfe:</b> Bestimmt durch die Einbauten. '
    'Moderne Instrumente gr\u00f6\u00dfer als historische (dickere Stimmplatten, '
    'schr\u00e4ge Stimmst\u00f6cke, verst\u00e4rkte Hebel, pr\u00e4zisere Lager).',
    '<b>Balgabmessungen (vierreihig diatonisch):</b> '
    '20\u00d734 bis 20\u00d740 cm. Faltentiefe 35 mm. '
    'Breite 20\u201321 cm.',
    '<b>Nicht die Faltenzahl ist ma\u00dfgeblich</b> \u2014 sondern der nutzbare '
    'Bewegungsbereich (geschlossen bis maximal ge\u00f6ffnet). '
    'Arml\u00e4nge und Bauh\u00f6he im geschlossenen Zustand bestimmen den Hub.',
    '<b>Luftverbrauch:</b> 5\u201315 l/min bei Akkordspiel. '
    'Spielzeit pro Balgzug \u2248 30\u201390 Sekunden je nach Spielweise.',
    '<b>Druck = Kraft / Fl\u00e4che.</b> Querschnitt 20\u00d734 = 680 cm\u00b2 '
    'bis 20\u00d740 = 800 cm\u00b2 (\u224818 % Unterschied). '
    'Bei gleicher Armkraft erzeugt der kleinere Querschnitt mehr Druck.',
    '<b>Nur vierreihige Instrumente</b> behandelt. Drei- oder f\u00fcnfreihige '
    'ver\u00e4ndern die Situation wesentlich.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Der Balg ist die Lunge des Instruments \u2014 '
    'seine Gr\u00f6\u00dfe bestimmt den Atem.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0022 \u2014 Balg \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
