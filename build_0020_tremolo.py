#!/usr/bin/env python3
"""Dok. 0020 — Tremolo: Schwebungsphysik, Typen, Stimmungspraxis"""
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
pdfmetrics.registerFont(TTFont('DejaVuMono','/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))
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
sF=ParagraphStyle('Fo',parent=sB,fontSize=10,alignment=TA_CENTER,spaceAfter=8,spaceBefore=4,fontName='DejaVuMono')
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

# ══ DIAGRAMME ══
def fig_schweb():
    """Schwebungsfrequenz vs. Grundton für verschiedene Tremolo-Stärken"""
    fig,ax=plt.subplots(figsize=(10,5.5))
    freqs=np.linspace(100,1200,200)
    for dc,col,ls,name in [(1,'#2e7d32','-','1 Cent (leicht)'),(3,'#1565c0','-','3 Cent (mittel)'),
        (5,'#ff9800','-','5 Cent (kr\u00e4ftig)'),(8,'#e94560','--','8 Cent (Wiener)'),
        (12,'#9c27b0','--','12 Cent (Musette)')]:
        schweb=freqs*dc/1731.23
        ax.plot(freqs,schweb,ls,color=col,lw=2,label=name)
    ax.set_xlabel('Grundfrequenz [Hz]'); ax.set_ylabel('Schwebung [Hz]')
    ax.set_title('Abb. 1: Tremolo-Schwebung vs. Grundfrequenz')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3)
    ax.set_xticks([131,220,262,440,523,880,1047])
    ax.set_xticklabels(['C3','A3','C4','A4','C5','A5','C6'],fontsize=8)
    fig.tight_layout(); return fig

def fig_zeitverlauf():
    """Zeitverlauf: Schwebung bei verschiedenen Tremolo-Stärken"""
    fig,axes=plt.subplots(3,1,figsize=(10,7),sharex=True)
    t=np.linspace(0,1.0,20000)
    f1=440.0
    for ax,dc,name in zip(axes,[2,5,12],['Leicht (2 Cent)','Kr\u00e4ftig (5 Cent)','Musette (12 Cent)']):
        f2=f1*2**(dc/1200)
        sig=np.sin(2*np.pi*f1*t)+np.sin(2*np.pi*f2*t)
        env=2*np.abs(np.cos(np.pi*(f2-f1)*t))
        ax.fill_between(t*1000,-env,env,alpha=0.15,color='#1565c0')
        ax.plot(t*1000,env,'r-',lw=1.5,label=f'Einh\u00fcllende ({f2-f1:.2f} Hz)')
        ax.set_title(f'{name}: {dc} Cent = {f2-f1:.2f} Hz Schwebung',fontsize=10)
        ax.legend(fontsize=8); ax.set_ylabel('Ampl.'); ax.set_ylim(-2.3,2.3)
    axes[2].set_xlabel('Zeit [ms]'); axes[0].set_xlim(0,1000)
    fig.suptitle('Abb. 2: Tremolo-Zeitverlauf bei A4 (440 Hz)',fontweight='bold')
    fig.tight_layout(); return fig

def fig_wechselwirkung():
    """Wechselwirkung Tremolo × Stimmung"""
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
    # Links: Tremolo-Spreizung auf Differenztönen
    dcs=np.linspace(0,15,50)
    f1=261.63; f2=329.63
    spreads=[]
    for dc in dcs:
        df1=f1*dc/1731.23; df2=f2*dc/1731.23
        diffs=[f2-f1,(f2+df2)-f1,f2-(f1+df1),(f2+df2)-(f1+df1)]
        spreads.append(max(diffs)-min(diffs))
    ax1.plot(dcs,spreads,'r-',lw=2.5)
    ax1.axhline(y=2.6,color='blue',ls='--',lw=1,label='Terz-Schwebung gleichstufig')
    ax1.set_xlabel('Tremolo [Cent]'); ax1.set_ylabel('Spreizung Differenzton [Hz]')
    ax1.set_title('Tremolo-Spreizung auf Terz C4\u2013E4')
    ax1.legend(fontsize=9); ax1.grid(True,alpha=0.3)
    ax1.axvspan(0,3,alpha=0.05,color='green'); ax1.axvspan(5,15,alpha=0.05,color='red')
    ax1.text(1.5,3.5,'fein\nmerkbar',fontsize=8,color='green',ha='center')
    ax1.text(10,3.5,'verwischt\nStimmung',fontsize=8,color='red',ha='center')
    # Rechts: Gleichmäßigkeit bei rein auf C-Dur vs. temperiert
    tonarten=range(12)
    noten=['C','C#','D','Eb','E','F','F#','G','G#','A','Bb','B']
    schweb_temp=[]; schweb_rein=[]
    # Rein auf C-Dur: Die drei Dur-Dreiklänge I, IV, V werden rein gestimmt
    # → E, A, B jeweils 13.7 Cent abgesenkt (reine große Terz zu C, F, G)
    korr = [0]*12
    korr[4] = -13.7  # E abgesenkt (reine Terz C-E)
    korr[9] = -13.7  # A abgesenkt (reine Terz F-A)
    korr[11] = -13.7 # B abgesenkt (reine Terz G-B)
    for s in tonarten:
        f1v=440*2**((s-9)/12)
        # Temperiert: große Terz = 400 Cent
        f2t=f1v*2**(4/12)
        schweb_temp.append(abs((f2t-f1v)-f1v/4))
        # Rein auf C: Grundton und Terz haben evtl. Korrektur
        terz_semi = (s+4) % 12
        f2r=f1v*2**((400 + korr[terz_semi] - korr[s])/1200)
        schweb_rein.append(abs((f2r-f1v)-f1v/4))
    ax2.bar(np.arange(12)-0.2,schweb_temp,0.35,color='#1565c0',label='Temperiert')
    ax2.bar(np.arange(12)+0.2,schweb_rein,0.35,color='#e94560',label='Rein (auf C)')
    # Tremolo-Spreizung addiert sich: bei 5 Cent Tremolo
    dc_trem=5
    spreiz_temp=[s + s*dc_trem/100*2 for s in schweb_temp]  # vereinfacht
    spreiz_rein=[s + s*dc_trem/100*2 for s in schweb_rein]
    # Besser: exakte Spreizung berechnen
    spreiz_rein2=[]
    for s in range(12):
        f1v=440*2**((s-9)/12)
        terz_semi=(s+4)%12
        f2r=f1v*2**((400+korr[terz_semi]-korr[s])/1200)
        df1=f1v*dc_trem/1731.23; df2=f2r*dc_trem/1731.23
        diffs=[f2r-f1v,(f2r+df2)-f1v,f2r-(f1v+df1),(f2r+df2)-(f1v+df1)]
        spread=max(diffs)-min(diffs)
        base=abs((f2r-f1v)-f1v/4)
        spreiz_rein2.append(base+spread)
    spreiz_temp2=[]
    for s in range(12):
        f1v=440*2**((s-9)/12)
        f2t=f1v*2**(4/12)
        df1=f1v*dc_trem/1731.23; df2=f2t*dc_trem/1731.23
        diffs=[f2t-f1v,(f2t+df2)-f1v,f2t-(f1v+df1),(f2t+df2)-(f1v+df1)]
        spread=max(diffs)-min(diffs)
        base=abs((f2t-f1v)-f1v/4)
        spreiz_temp2.append(base+spread)
    ax2.bar(np.arange(12)-0.2,spreiz_temp2,0.35,color='#1565c0',alpha=0.3,edgecolor='#1565c0',lw=1,ls='--')
    ax2.bar(np.arange(12)+0.2,spreiz_rein2,0.35,color='#e94560',alpha=0.3,edgecolor='#e94560',lw=1,ls='--')
    ax2.text(0.5,max(spreiz_rein2)*0.95,'Gest\u00e4nderter Bereich\n= mit 5 Ct Tremolo',fontsize=7,color='gray')
    ax2.set_xticks(range(12)); ax2.set_xticklabels(noten,fontsize=8)
    ax2.set_ylabel('Differenzton-Schwebung [Hz]')
    ax2.set_title('Alle Tonarten: Rein auf C vs. temperiert\n(ausgef\u00fcllt = ohne Tremolo, transparent = mit 5 Ct Tremolo)')
    ax2.legend(fontsize=8); ax2.grid(True,alpha=0.3,axis='y')
    fig.suptitle('Abb. 3: Wechselwirkung Tremolo \u00d7 Stimmung',fontweight='bold')
    fig.tight_layout(); return fig

# ══ PDF ══
outfile='0020_tremolo_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0020',sST))
story.append(Paragraph('Tremolo:<br/>Schwebungsphysik, Typen und Stimmungspraxis',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('Physik der Schwebung, Tremolo-Typen (Null bis Musette), Wechselwirkung '
    'mit reiner Stimmung und Differenzt\u00f6nen (Dok. 0019). Warum Instrumente mit '
    'viel Tremolo temperiert gestimmt werden m\u00fcssen. Referenz: Dok. 0012, 0019.',sAb))
story.append(Paragraph('<b>Berechnungsskript:</b> <b>pruefskript_0020_tremolo.py</b> '
    'im selben Verzeichnis auf GitHub.',sAb))

# Kap. 1
story.append(Paragraph('1. Schwebungsphysik',sCh))
story.append(Paragraph(
    'Zwei T\u00f6ne f\u2081 und f\u2082 = f\u2081 + \u0394f erzeugen eine Schwebung '
    'mit der Frequenz f_schweb = |\u0394f|. Das Geh\u00f6r nimmt sie als periodische '
    'Lautst\u00e4rkeschwankung wahr. Die Umrechnung zwischen Cent und Hz:',sB))
story.append(Paragraph('\u0394f \u2248 f\u2081 \u00d7 \u0394c / 1731',sF))
story.append(Paragraph(
    'Die Schwebungsfrequenz ist <b>proportional zur Grundfrequenz</b>. '
    'Gleiche Cent-Verstimmung erzeugt im oberen Diskant eine doppelt so schnelle '
    'Schwebung wie eine Oktave tiefer. Manche Stimmer gleichen das aus, '
    'indem sie im oberen Bereich etwas weniger Tremolo geben.',sB))
fig1=fig_schweb(); story.append(fig2img(fig1,150)); story.append(Spacer(1,3*mm))

# Tabelle: Schwebungsfrequenzen
rows_sf=[]
for dc,name in [(1,'Sehr leicht'),(2,'Leicht'),(3,'Mittel'),(5,'Kr\u00e4ftig'),(8,'Wiener'),(12,'Musette')]:
    row=[f'{dc} ({name})']
    for f in [130.8,261.6,440,523.3,880,1047]:
        row.append(f'{f*dc/1731.23:.2f}')
    rows_sf.append(row)
story.append(mk_tbl(['\u0394c [Cent]','C3','C4','A4','C5','A5','C6'],rows_sf,
    cw=[28*mm,16*mm,16*mm,16*mm,16*mm,16*mm,16*mm]))

# Kap. 2
story.append(Paragraph('2. Tremolo-Typen',sCh))
fig2=fig_zeitverlauf(); story.append(fig2img(fig2,155)); story.append(Spacer(1,3*mm))
rows_typ=[]
rows_typ.append(['Null-Tremolo','< 0,5','Einch\u00f6rig. Voraussetzung f\u00fcr <b>reine Stimmung</b>.'])
rows_typ.append(['Leicht','1\u20133','Steirisch, alpenl\u00e4ndisch. W\u00e4rme ohne Vibrato.'])
rows_typ.append(['Mittel','3\u20135','Typisch f\u00fcr Begleiter. Lebendiger Klang.'])
rows_typ.append(['Kr\u00e4ftig','5\u20138','Wiener Stimmung. Deutliches Vibrato.'])
rows_typ.append(['Musette','10\u201315','Franz\u00f6sisch. Starkes Vibrato, Caf\u00e9-Klang.'])
story.append(mk_tbl(['Typ','Cent','Beschreibung'],rows_typ,
    cw=[25*mm,18*mm,82*mm]))
story.append(Paragraph(
    '<b>Feinstimmung:</b> Jedes Instrument muss feingestimmt werden \u2014 '
    'unabh\u00e4ngig vom Tremolo. Bei Instrumenten mit viel Tremolo mag es '
    'vertretbar sein, nicht ganz so sauber zu stimmen und auf das '
    'aufw\u00e4ndige Stimmen im eingebauten Zustand zu verzichten. '
    'Null-Tremolo ist nicht Voraussetzung f\u00fcr Feinstimmung, '
    'sondern f\u00fcr <b>reine Stimmung</b> (Dok. 0019).',sB))

# Kap. 3
story.append(Paragraph('3. Mehrch\u00f6rige Stimmung',sCh))
story.append(Paragraph(
    '<b>Einch\u00f6rig:</b> Kein Tremolo. Reine Stimmung m\u00f6glich und sinnvoll. '
    '<b>Zweich\u00f6rig:</b> 1 \u00d7 Grundstimmung + 1 \u00d7 Tremolo-Chor. '
    '<b>Dreich\u00f6rig:</b> Verschiedene Konfigurationen (2+1 oder 1+1+1). '
    '<b>Vierch\u00f6rig (Musette):</b> 2 \u00d7 Grundstimmung + 2 \u00d7 Tremolo (\u00b1).',sB))
story.append(Paragraph(
    'Typische Verstimmung: Steirisch +1 bis +3 Cent (nur oben), '
    'Wiener +5 bis +8 Cent (oben), Musette \u00b18 bis \u00b112 Cent (symmetrisch). '
    'Das Tremolo steigt proportional mit der Frequenz \u2014 bei gleichen Cent '
    'klingt der obere Diskant schneller als der Bass.',sB))

# Kap. 4
story.append(Paragraph('4. Wechselwirkung: Tremolo \u00d7 Stimmung',sCh))
story.append(Paragraph(
    'Tremolo erzeugt eine <b>Spreizung</b> auf den Differenzt\u00f6nen. '
    'Wenn Grundton und Terz jeweils ihren Tremolo-Chor haben, schwankt der '
    'Differenzton zwischen vier Werten (Grund+Grund, Grund+Tremolo, '
    'Tremolo+Grund, Tremolo+Tremolo). Die Spreizung w\u00e4chst mit dem Tremolo.',sB))
fig3=fig_wechselwirkung(); story.append(fig2img(fig3,158)); story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Abb. 3 rechts zeigt den doppelten Effekt: Die <b>ausgef\u00fcllten</b> Balken sind '
    'die Schwebung ohne Tremolo \u2014 bei temperierter Stimmung (blau) gleichm\u00e4\u00dfig, '
    'bei reiner Stimmung auf C (rot) mit drei Nullstellen und drei Verschlechterungen. '
    'Die <b>transparenten</b> Balken zeigen, was bei 5 Cent Tremolo passiert: '
    'Die Tremolo-Spreizung <b>addiert sich</b> zur bestehenden Schwebung. '
    'Die ohnehin schlechteren Terzen (E\u2013G#, A\u2013C#, B\u2013Eb) werden durch das '
    'Tremolo \u00fcberproportional verschlechtert \u2014 die bereits doppelt so hohen '
    'Balken wachsen noch weiter. Die guten Terzen (C\u2013E, F\u2013A, G\u2013B) waren '
    'bei Null und steigen nur auf die Tremolo-Spreizung an.',sB))

# Spreizungs-Tabelle
rows_sp=[]
f1=261.63; f2=329.63
for dc in [0,1,3,5,8,12]:
    df1=f1*dc/1731.23; df2=f2*dc/1731.23
    diffs=[f2-f1,(f2+df2)-f1,f2-(f1+df1),(f2+df2)-(f1+df1)]
    sp=max(diffs)-min(diffs)
    ct=1200*np.log2(max(diffs)/min(diffs)) if min(diffs)>0 and sp>0.001 else 0
    name={0:'kein',1:'leicht',3:'mittel',5:'kr\u00e4ftig',8:'Wiener',12:'Musette'}[dc]
    rows_sp.append([f'{dc} ({name})',f'{min(diffs):.1f}',f'{max(diffs):.1f}',f'{sp:.2f}',f'{ct:.0f}'])
story.append(mk_tbl(['Tremolo','Diff min [Hz]','Diff max [Hz]','Spreizung [Hz]','Spreiz. [Cent]'],rows_sp,
    cw=[24*mm,22*mm,22*mm,22*mm,22*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Ab 3 Cent Tremolo</b> ist die Spreizung (14,5 Cent) so gro\u00df wie die '
    'gesamte Terzen-Korrektur (14 Cent). Die reine Stimmung wird durch das '
    'Tremolo <b>wieder verwischt</b>. Bei Musette (12 Cent) betr\u00e4gt die '
    'Spreizung 57 Cent \u2014 das ist fast ein Halbton.',sW))

# Kap. 5
story.append(Paragraph('5. Warum temperiert bei Tremolo?',sCh))
story.append(Paragraph(
    '<b>Temperiert + Tremolo:</b> Alle Terzen in allen Tonarten schweben '
    'gleich schnell. Das Tremolo klingt <b>gleichm\u00e4\u00dfig</b> \u2014 '
    'kein Akkord f\u00e4llt aus dem Rahmen.',sB))
story.append(Paragraph(
    '<b>Rein + Tremolo:</b> Terzen in der Bezugstonart schweben weniger '
    '(n\u00e4her an reiner Stimmung), aber Terzen in anderen Tonarten '
    'schweben <b>mehr</b>. Das Tremolo klingt <b>ungleichm\u00e4\u00dfig</b> \u2014 '
    'manche Akkorde vibrieren schneller als andere. Das wirkt verstimmt, '
    'nicht rein.',sB))
story.append(Paragraph(
    '<b>Instrumente mit Tremolo > 3 Cent:</b> Temperiert (gleichstufig) stimmen. '
    'Die Vorteile reiner Stimmung werden durch das Tremolo aufgehoben. '
    'Gleichm\u00e4\u00dfiges Tremolo \u00fcber alle Tonarten ist wichtiger als '
    'reine Terzen in einer Tonart.',sK))

# Kap. 6: Bezeichnungen
story.append(Paragraph('6. Tremolo-Bezeichnungen (Industrie)',sCh))
story.append(Paragraph('Hersteller verwenden unterschiedliche Bezeichnungen. '
    'Die Angabe erfolgt entweder in <b>Cent</b> (tonh\u00f6henunabh\u00e4ngig) oder '
    'in <b>Hz bei A4 = 440 Hz</b> (nur f\u00fcr diesen Ton g\u00fcltig). '
    'Umrechnung: 1 Hz bei A4 \u2248 3,9 Cent.',sB))
rows_bez=[]
for hz,ct,name,genre in [(0,0,'Unison / Dry / Secco','Klassik, Balkan'),
    (0.5,2,'Concert / Swing Secco','Jazz, Tango'),
    (1,4,'Swing / Swing Mosso','Gypsy Jazz, Klezmer'),
    (2,7,'Demi-Swing','Irish'),
    (2.5,10,'American','Cajun, Qu\u00e9b\u00e9cois'),
    (3,12,'Moderate / Slovenisch','Tex-Mex, Alpin'),
    (4,15,'Standard / German','Polka, Continental'),
    (5,18,'Fast / Modern French','Modernes Musette'),
    (6,22,'Very Fast / Old French','Altes Musette'),
    (7,25,'Extremely Fast','Scottish')]:
    rows_bez.append([f'{ct}',f'{hz:.1f}',name,genre])
story.append(mk_tbl(['Cent','Hz (A4)','Bezeichnung','Stil'],rows_bez,
    cw=[14*mm,16*mm,38*mm,38*mm]))
story.append(Paragraph(
    '<b>Achtung:</b> Die Hz-Werte gelten nur bei A4 = 440 Hz. '
    'Bei C4 (262 Hz) ist 1 Hz \u2248 6,6 Cent, bei C6 (1047 Hz) nur \u2248 1,7 Cent. '
    'F\u00fcr die Stimmung ist <b>Cent</b> die sicherere Angabe '
    '\u2014 sie ist tonh\u00f6henunabh\u00e4ngig.',sB))

# Kap. 7: Sym/Asym
story.append(Paragraph('7. Symmetrisch und asymmetrisch',sCh))
story.append(Paragraph(
    '<b>Dreich\u00f6rig symmetrisch (MMM):</b> M(\u2212x) + M(0) + M(+x). '
    'Obere und untere Zunge besitzen dieselbe Verstimmung in Cent. '
    'Standard f\u00fcr franz\u00f6sisches Musette.',sB))
story.append(Paragraph(
    '<b>Asymmetrisch (z.B. slovenisch):</b> Die Unterschwebung ist z.B. '
    '<b>doppelt so viel</b> wie die Oberschwebung: M(\u22122x) + M(0) + M(+x). '
    'Beispiel: M(\u221210 Cent) + M(0) + M(+5 Cent). Das ergibt einen '
    'anderen Klangcharakter \u2014 die Unterschwebung klingt w\u00e4rmer, '
    'die Oberschwebung brillanter.',sB))
story.append(Paragraph(
    '<b>Wichtig: Symmetrisch in Hz \u2260 symmetrisch in Cent.</b> '
    'Cent ist logarithmisch. \u00b15 Cent ergibt +1,273 Hz oben, aber '
    '\u22121,269 Hz unten \u2014 eine kleine Asymmetrie in Hz. '
    'Umgekehrt: \u00b12 Hz ergibt +7,85 Cent oben, aber \u22127,89 Cent unten. '
    'Bei \u00fcblichen Tremolo-Werten (< 20 Cent) ist der Unterschied < 1 % '
    'und praktisch irrelevant. Aber das Prinzip gilt.',sB))
story.append(Paragraph(
    '<b>Zweich\u00f6rig (MM):</b> Nur eine Tremolo-Zunge: M(0) + M(+x). '
    'Man stimmt die obere <b>und</b> die untere Zunge \u2014 nicht die mittlere '
    'und die obere (die mittlere Zunge existiert nicht).',sK))

story.append(Paragraph('<b>Warum h\u00f6rt man nur 1 Hz bei 2 Hz unter + 1 Hz oben?</b>',sCh))
story.append(Paragraph(
    'Drei Zungen M(\u22122 Hz) + M(0) + M(+1 Hz) erzeugen drei Schwebungs-Paare: '
    '1 Hz, 2 Hz und 3 Hz. Aber man h\u00f6rt <b>nicht</b> drei separate Schwebungen, '
    'sondern die Einh\u00fcllende der Summe. Diese wiederholt sich mit der Frequenz '
    'GCD(1, 2) = <b>1 Hz</b>. Innerhalb jeder Sekunde variiert die Amplitude \u2014 '
    'die 2-Hz- und 3-Hz-Schwebungen erzeugen ein Muster, aber das Muster '
    'wiederholt sich exakt jede Sekunde.',sB))
story.append(Paragraph(
    '<b>Der Variationseffekt:</b> Wenn die Schwebungswerte <b>ganzzahlige</b> Vielfache '
    'sind (z.B. 2:1), wiederholt sich das Muster exakt \u2192 klingt gleichm\u00e4\u00dfig, '
    'etwas mechanisch. Wenn die Werte <b>kein</b> einfaches Verh\u00e4ltnis haben '
    '(z.B. 2,3 Hz / 1,1 Hz), wiederholt sich das Muster nie exakt \u2192 '
    'die Einh\u00fcllende <b>variiert st\u00e4ndig</b> \u2192 klingt lebendig, organisch.',sB))

rows_var=[]
rows_var.append(['\u22122,0 / +1,0','2:1 (ganzzahlig)','1,0 s','gleichm\u00e4\u00dfig, mechanisch'])
rows_var.append(['\u22122,3 / +1,1','irrational','10 s','lebendig, variierend'])
rows_var.append(['\u22122,0 / +1,4','irrational','5 s','lebendig, variierend'])
rows_var.append(['\u22121,8 / +1,2','3:2','1,7 s','langsame Variation'])
story.append(mk_tbl(['Unter/Ober [Hz]','Verh\u00e4ltnis','Wiederholung','Klangcharakter'],rows_var,
    cw=[25*mm,28*mm,22*mm,35*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph(
    '<b>Praxis:</b> F\u00fcr lebendigen Klang sollten Unter- und Oberschwebung '
    '<b>kein</b> einfaches ganzzahliges Verh\u00e4ltnis haben. '
    'Die Unregelm\u00e4\u00dfigkeit ist die Qualit\u00e4t \u2014 '
    'erfahrene Stimmer lassen die Werte bewusst \u201ekrumm\u201c.',sK))

# Kap. 8: Verläufe
story.append(Paragraph('8. Tremolo-Verl\u00e4ufe',sCh))
story.append(Paragraph(
    '<b>Konstant Hz:</b> Oben trockener. <b>Konstant Cent:</b> Oben sch\u00e4rfer. '
    '<b>S-Form</b> (die \u00fcblichste Praxis): Unten reduziert (w\u00e4rmer), '
    'in der Mitte am st\u00e4rksten, oben wieder reduziert (nicht zu scharf). '
    'Die Abweichung vom konstanten Cent-Verlauf bildet ein S: '
    'unten abgeflacht, Mitte steil ansteigend, oben wieder abgeflacht.',sB))
# S-Form Diagramm
fig_s,ax_s=plt.subplots(figsize=(10,4.5))
tx=np.arange(37); fx=[174.6*2**(i/12) for i in tx]
# Konstant Cent → Hz steigt linear mit f
c_const=5.0
hz_const=[f*c_const/1731.23 for f in fx]
# Konstant Hz
hz_flat=[2.0]*37
# S-Form: unten weniger, Mitte Maximum, oben wieder weniger
# Sigmoid-artig: c(x) = c_base + c_peak × sin(π×x/36)² aber mit Abflachung oben
x_norm=[i/36.0 for i in tx]
# S: unten 3 Cent, Mitte 6 Cent, oben 4 Cent
c_s=[3.0 + 3.0*np.sin(np.pi*xn)**0.8 * (1 - 0.4*xn) for xn in x_norm]
hz_s=[f*c/1731.23 for f,c in zip(fx,c_s)]
ax_s.plot(tx,hz_const,'b--',lw=1.5,label=f'Konstant {c_const:.0f} Cent')
ax_s.plot(tx,hz_flat,'g--',lw=1.5,label='Konstant 2,0 Hz')
ax_s.plot(tx,hz_s,'r-',lw=2.5,label='S-Form')
# Cent-Werte als zweite Info
ax_s2=ax_s.twinx()
ax_s2.plot(tx,c_s,'r:',lw=1,alpha=0.5)
ax_s2.set_ylabel('Cent (S-Form)',color='red',alpha=0.5)
ax_s2.tick_params(axis='y',labelcolor='red')
ax_s.set_xticks([0,7,12,19,24,31,36]); ax_s.set_xticklabels(['F3','C4','F4','C5','F5','C6','F6'],fontsize=8)
ax_s.set_ylabel('Schwebung [Hz]'); ax_s.set_title('Abb. 4: Tremolo-Verl\u00e4ufe \u2014 S-Form vs. konstant')
ax_s.legend(fontsize=9,loc='upper left'); ax_s.grid(True,alpha=0.3); fig_s.tight_layout()
story.append(fig2img(fig_s,148)); story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>S-Form:</b> Unten weniger Tremolo als der konstante Cent-Verlauf '
    '(\u2248 3 Cent statt 5), in der Mitte das Maximum (\u2248 6 Cent), '
    'oben wieder reduziert (\u2248 4 Cent). Das verhindert, dass der Bass '
    'matscht (zu viel Schwebung bei tiefen T\u00f6nen) und der Diskant '
    'scharf klingt (zu schnelle Schwebung bei hohen T\u00f6nen). '
    'Umgekehrt (viel unten, wenig oben) ergibt keinen musikalischen Sinn.',sK))

# Kap. 9: Akkord-Verstimmung
story.append(Paragraph('9. Akkord-Verstimmung durch Tremolo',sCh))
story.append(Paragraph(
    'Cross-Kombinationen (Grund einer Note + Tremolo einer anderen) erzeugen '
    '<b>verschiedene Terzen</b>. Die Terz \u201eschwimmt\u201c um \u00b1dc Cent:',sB))
rows_cr=[]; f1c=261.63; f3c=f1c*5/4
for dc in [0,3,5,8,12]:
    t=1200*np.log2(f3c/f1c); f3t=f3c*2**(dc/1200); f1t=f1c*2**(dc/1200)
    c1=1200*np.log2(f3t/f1c); c2=1200*np.log2(f3c/f1t)
    rows_cr.append([str(dc),f'{c1-t:+.1f}',f'{c2-t:+.1f}',f'{abs(c1-c2):.0f}'])
story.append(mk_tbl(['Trem. [Ct]','Cross\u2191','Cross\u2193','Spreizung'],rows_cr,cw=[20*mm,22*mm,22*mm,22*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Bei 5 Cent: \u00b15 Cent = 10 Cent Spreizung (\u224870 % der Terzen-Korrektur). '
    'Bei 8 Cent: 16 Cent Spreizung \u2014 <b>mehr</b> als die Korrektur. '
    'Die reine Stimmung wird durch das Tremolo komplett \u00fcberschrieben.',sW))

# Kap. 10: Feinstimmung
story.append(Paragraph('10. Feinstimmung: Werkzeug und Praxis',sCh))
story.append(Paragraph(
    'Feinstimmen der ausgebauten Stimmst\u00f6cke ist nicht zuverl\u00e4ssig \u2014 '
    'es gibt <b>Einbau-Verstimmungen</b>. Alle Verschraubungen m\u00fcssen eingebaut sein. '
    'Werkzeug: Kratzer und Feilen bevorzugen (nicht rotierende Schleifer). '
    'Hilfswerkzeug zum Herausangeln der innenliegenden Zungen: '
    'gebogene Dr\u00e4hte (\u00fcber Nachbarkanal, f\u00fcr tiefe Zungen) '
    'oder kleine Magnete + Lospl\u00e4ttchen. '
    'Bei sehr kurzen Zungen: Stimmstock ausbauen und Zunge \u00fcber das Loch herausdr\u00fccken.',sB))

# Kap. 11
story.append(Paragraph('11. Zusammenfassung',sCh))
pts=[
    '<b>Schwebung = f\u2081 \u00d7 \u0394c / 1731:</b> Proportional zu Grundfrequenz und '
    'Verstimmung. Im Diskant doppelt so schnell wie eine Oktave tiefer.',
    '<b>Tremolo-Typen:</b> Null (&lt; 0,5 Cent) bis Musette (10\u201315 Cent). '
    'Null-Tremolo ist Voraussetzung f\u00fcr <b>reine Stimmung</b> (Dok. 0019). '
    'Feingestimmt werden muss jedes Instrument.',
    '<b>Wechselwirkung:</b> Tremolo erzeugt eigene Spreizung auf Differenzt\u00f6nen. '
    'Ab 3 Cent: Reine Stimmung wird verwischt. Ab 5 Cent: Ungleichm\u00e4\u00dfigkeit h\u00f6rbar.',
    '<b>Tremolo > 3 Cent = temperiert stimmen.</b> Die Vorteile reiner Stimmung '
    'werden durch das Tremolo aufgehoben. Gleichm\u00e4\u00dfiges Tremolo ist wichtiger '
    'als reine Terzen.',
    '<b>Tremolo steigt mit der Frequenz:</b> Gleiche Cent = schnellere Schwebung oben. '
    'Manche Stimmer geben im oberen Diskant etwas weniger Cent.',
    '<b>S-Form ist Standard:</b> Unten weniger, Mitte Maximum, oben wieder reduziert. '
    'Verhindert Bass-Matsch und Diskant-Sch\u00e4rfe.',
    '<b>Akkord-Verstimmung:</b> Cross-Kombinationen erzeugen \u00b1dc Cent Spreizung. '
    'Bei 5 Cent Tremolo: 10 Cent Spreizung \u2014 fast die gesamte Terzen-Korrektur.',
    '<b>Bezeichnungen:</b> Dry (0) \u2192 Swing (4) \u2192 Demi-Swing (7) \u2192 '
    'American (10) \u2192 Slovenisch (12) \u2192 German (15) \u2192 French (18\u201322) Cent.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Das Tremolo macht den Klang lebendig \u2014 '
    'aber es verwischt die Feinheiten. Wer reinen Klang will, '
    'muss auf das Tremolo verzichten.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0020 \u2014 Tremolo \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
