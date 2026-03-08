#!/usr/bin/env python3
"""Dok. 0019 — Stimmung und Differenztöne"""
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

f_A4=440.0
noten=['C','C#','D','Eb','E','F','F#','G','G#','A','Bb','B']
def fgl(s,o=4): return f_A4*2**((s-9+(o-4)*12)/12)
rein={'kl. Terz':(6,5,3),'gr. Terz':(5,4,4),'Quarte':(4,3,5),
      'Quinte':(3,2,7),'kl. Sexte':(8,5,8),'gr. Sexte':(5,3,9)}

# ══ DIAGRAMME ══
def fig_intervalle():
    fig,ax=plt.subplots(figsize=(10,5))
    names=[]; deltas=[]; cols=[]
    for nm,(n,d,s) in rein.items():
        c_gl=s*100; c_r=1200*np.log2(n/d); delta=c_gl-c_r
        names.append(nm); deltas.append(delta)
        cols.append('#e94560' if abs(delta)>10 else '#2e7d32' if abs(delta)<5 else '#ff9800')
    ax.barh(range(len(names)),deltas,color=cols,edgecolor='black',lw=0.5)
    ax.set_yticks(range(len(names))); ax.set_yticklabels(names)
    ax.set_xlabel('Abweichung gleichstufig \u2212 rein [Cent]')
    ax.set_title('Abb. 1: Gleichstufig vs. rein \u2014 Terzen sind das Problem')
    ax.axvline(x=0,color='black',lw=1); ax.grid(True,alpha=0.3,axis='x')
    ax.axvspan(-5,5,alpha=0.1,color='green'); ax.text(0,-0.5,'\u00b15 Cent',fontsize=8,color='green',ha='center')
    fig.tight_layout(); return fig

def fig_schwebung():
    fig,ax=plt.subplots(figsize=(10,5))
    noten=['C','C#','D','Eb','E','F','F#','G','G#','A','Bb','B']
    schw_gr=[]; schw_kl=[]
    for s in range(12):
        f1=fgl(s); 
        dg_gr=f1*(2**(4/12)-1); dr_gr=f1*(5/4-1); schw_gr.append(abs(dg_gr-dr_gr))
        dg_kl=f1*(2**(3/12)-1); dr_kl=f1*(6/5-1); schw_kl.append(abs(dg_kl-dr_kl))
    x=np.arange(12); w=0.35
    ax.bar(x-w/2,schw_gr,w,color='#e94560',label='Gro\u00dfe Terz (gl. vs. rein)')
    ax.bar(x+w/2,schw_kl,w,color='#1565c0',label='Kleine Terz (gl. vs. rein)')
    ax.set_xticks(x); ax.set_xticklabels(noten); ax.set_ylabel('Schwebung [Hz]')
    ax.set_title('Abb. 2: Schwebung der Differenzt\u00f6ne \u2014 alle Tonarten')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3,axis='y')
    fig.tight_layout(); return fig

def fig_diff_ton():
    """Zeitlicher Verlauf: Differenzton bei reiner vs. gleichstufiger Terz"""
    fig,(ax1,ax2)=plt.subplots(2,1,figsize=(10,6),sharex=True)
    t=np.linspace(0,0.1,10000)
    f1=261.63
    from scipy.signal import hilbert
    f2r=f1*5/4; mix_r=np.sin(2*np.pi*f1*t)+np.sin(2*np.pi*f2r*t)
    ax1.plot(t*1000,mix_r,'b-',lw=0.5,alpha=0.5)
    ax1.plot(t*1000,np.abs(hilbert(mix_r)),'r-',lw=1.5,label='Einh\u00fcllende')
    ax1.set_title('Reine gro\u00dfe Terz (5/4): KEINE Schwebung'); ax1.legend(fontsize=9)
    ax1.set_ylabel('Amplitude'); ax1.set_ylim(-2.5,2.5)
    f2g=f1*2**(4/12); mix_g=np.sin(2*np.pi*f1*t)+np.sin(2*np.pi*f2g*t)
    ax2.plot(t*1000,mix_g,'b-',lw=0.5,alpha=0.5)
    ax2.plot(t*1000,np.abs(hilbert(mix_g)),'r-',lw=1.5,label=f'Schwebung \u2248{abs(f2g-f2r):.1f} Hz')
    ax2.set_title('Gleichstufige gro\u00dfe Terz: Schwebung sichtbar'); ax2.legend(fontsize=9)
    ax2.set_xlabel('Zeit [ms]'); ax2.set_ylabel('Amplitude'); ax2.set_ylim(-2.5,2.5)
    fig.suptitle('Abb. 3: Rein vs. gleichstufig \u2014 C4+E4',fontweight='bold')
    fig.tight_layout(); return fig

def fig_diff_uebersicht():
    """Differenztöne 1. Ordnung: Gleichstufig vs. Rein, über Tonbereich"""
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(13,5.5))
    # Links: Große Terz - Differenztöne und Schwebung
    semis=range(0,12); octs=[3,4,5,6]
    for o,col,ms in zip(octs,['#2e7d32','#1565c0','#e94560','#9c27b0'],[4,5,6,7]):
        ff=[]; dg=[]; dr=[]; sw=[]
        for s in semis:
            f1=fgl(s,o); ff.append(f1)
            dg.append(f1*(2**(4/12)-1)); dr.append(f1*(5/4-1))
            sw.append(abs(dg[-1]-dr[-1]))
        ax1.plot(ff,sw,'o-',color=col,ms=ms,lw=1.5,label=f'Oktave {o}')
    ax1.set_xlabel('Grundfrequenz f\u2081 [Hz]'); ax1.set_ylabel('Schwebung [Hz]')
    ax1.set_title('Gro\u00dfe Terz: Schwebung der Differenzt\u00f6ne')
    ax1.legend(fontsize=8); ax1.grid(True,alpha=0.3)
    ax1.axhspan(0,2,alpha=0.05,color='green'); ax1.axhspan(4,20,alpha=0.05,color='red')
    ax1.text(200,1,'wahrnehmbar',fontsize=8,color='green')
    ax1.text(200,8,'unangenehm',fontsize=8,color='red')
    # Rechts: Alle Intervalle bei C4
    intervalle=['kl.Terz','gr.Terz','Quarte','Quinte','kl.Sexte','gr.Sexte']
    verh={'kl.Terz':(6,5,3),'gr.Terz':(5,4,4),'Quarte':(4,3,5),'Quinte':(3,2,7),'kl.Sexte':(8,5,8),'gr.Sexte':(5,3,9)}
    f1=261.63; sw_int=[]; names_int=[]
    for nm in intervalle:
        n,d,s=verh[nm]; dg=f1*(2**(s/12)-1); dr=f1*(n/d-1); sw_int.append(abs(dg-dr))
        names_int.append(nm)
    cols2=['#e94560' if s>3 else '#ff9800' if s>1 else '#2e7d32' for s in sw_int]
    ax2.barh(range(len(names_int)),sw_int,color=cols2,edgecolor='black',lw=0.5)
    ax2.set_yticks(range(len(names_int))); ax2.set_yticklabels(names_int)
    ax2.set_xlabel('Schwebung [Hz] bei C4'); ax2.set_title('Alle Intervalle: Schwebung bei C4')
    ax2.grid(True,alpha=0.3,axis='x')
    ax2.axvline(x=2,color='orange',ls='--',lw=1); ax2.axvline(x=0.5,color='green',ls='--',lw=1)
    fig.suptitle('Abb. 4: Differenzt\u00f6ne \u2014 \u00dcbersicht',fontweight='bold')
    fig.tight_layout(); return fig

def fig_akkord():
    """Dur-Dreiklang: Differenztöne gleichstufig vs. rein"""
    fig,ax=plt.subplots(figsize=(10,5.5))
    octs=[3,4,5,6]
    # C-E Schwebung, E-G Schwebung, C-G Schwebung
    for pair,ratio_gl,ratio_rein,col,lab in [
        ('C-E',2**(4/12),5/4,'#e94560','gro\u00dfe Terz C-E'),
        ('E-G',2**(3/12),6/5,'#1565c0','kleine Terz E-G'),
        ('C-G',2**(7/12),3/2,'#2e7d32','Quinte C-G')]:
        ff=[]; sw=[]
        for o in [3,3.5,4,4.5,5,5.5,6]:
            f1=261.63*2**(o-4)
            dg=f1*(ratio_gl-1); dr=f1*(ratio_rein-1); sw.append(abs(dg-dr)); ff.append(f1)
        ax.plot(ff,sw,'o-',color=col,lw=2,ms=5,label=lab)
    ax.set_xlabel('Grundfrequenz [Hz]'); ax.set_ylabel('Schwebung [Hz]')
    ax.set_title('Abb. 5: Dur-Dreiklang \u2014 Schwebung pro Intervall')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3)
    ax.set_xscale('log'); ax.set_xticks([130,262,523,1047]); ax.set_xticklabels(['C3','C4','C5','C6'])
    ax.axhspan(0,2,alpha=0.05,color='green'); ax.axhspan(5,20,alpha=0.05,color='red')
    fig.tight_layout(); return fig

# ══ PDF ══
outfile='0019_stimmung_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0019',sST))
story.append(Paragraph('Stimmung und Differenzt\u00f6ne',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('Gleichstufig temperierte, reine und historische Stimmungen. '
    'Differenzt\u00f6ne bei Terzen und Quinten. Euler-Tonnetz und M\u00f6glichkeiten '
    'diatonischer Instrumente. Druckabh\u00e4ngigkeit und Tremolo als Grenzen. '
    'Referenz: Dok. 0012.',sAb))
story.append(Paragraph('<b>Berechnungsskript:</b> <b>pruefskript_0019_stimmung.py</b> '
    'im selben Verzeichnis auf GitHub.',sAb))

# Kap. 1
story.append(Paragraph('1. Stimmungssysteme im Vergleich',sCh))
story.append(Paragraph(
    'Die <b>gleichstufig temperierte Stimmung</b> (12-TET) teilt die Oktave in 12 gleiche '
    'Halbtonschritte (Faktor 2\u00b9\u2070\u00b2 = 1,05946). Kein Intervall au\u00dfer der Oktave '
    'ist rein \u2014 alle sind Kompromisse. Die <b>reine Stimmung</b> verwendet ganzzahlige '
    'Frequenzverh\u00e4ltnisse (3/2 f\u00fcr die Quinte, 5/4 f\u00fcr die gro\u00dfe Terz), '
    'ist aber nicht auf alle Tonarten \u00fcbertragbar. '
    '<b>Valotti</b> und andere historische Temperierungen verteilen die Abweichung '
    'ungleichm\u00e4\u00dfig \u2014 manche Tonarten klingen reiner, andere weniger. '
    'Valotti temperiert sechs Quinten (F\u2013C\u2013G\u2013D\u2013A\u2013E\u2013B) je um 1/6 pythagor\u00e4isches '
    'Komma und l\u00e4sst die \u00fcbrigen sechs rein.',sB))

# Valotti-Tabelle
rows_val=[]
val_c={'C':+6,'C#':0,'D':+2,'Eb':+4,'E':-2,'F':+8,'F#':-2,'G':+4,'G#':+2,'A':0,'Bb':+6,'B':-4}
for n in noten: rows_val.append([n,f'{val_c.get(n,0):+d}'])
story.append(mk_tbl(['Note','Valotti vs. gleichstufig [Cent]'],rows_val,
    cw=[25*mm,50*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Die Abweichungen sind klein (max. \u00b18 Cent). Valotti ist ein Kompromiss '
    'mit <b>entsch\u00e4rften Terzen</b> in h\u00e4ufigen Tonarten (C, G, D, F). '
    'In entfernten Tonarten (F#, C#) sind die Terzen pythagor\u00e4isch '
    '(schlechter als gleichstufig).',sB))
story.append(Paragraph(
    '<b>Wichtig:</b> Diese Tabelle gilt nur f\u00fcr <b>C-Dur als Bezugstonart</b>. '
    'F\u00fcr jede andere Grundtonart muss der Bezugston verschoben werden. '
    'Bei einem G-C-F-B Instrument (Bezug C) passt die Tabelle direkt. '
    'Bei einem A-D-G-C Instrument (Bezug D) m\u00fcssen alle Werte um 2 Halbt\u00f6ne '
    'verschoben werden: D bekommt den Wert von C (+6), E den von D (+2) usw. '
    'Die entsch\u00e4rften Terzen liegen dann um D, A, E, G \u2014 '
    'passend f\u00fcr die h\u00e4ufigsten Tonarten dieses Instruments.',sW))

# Begleiter-Praxis
story.append(Paragraph(
    '<b>Begleiter (diatonisch):</b> In der Praxis werden Begleiter oft so gestimmt, '
    'dass die <b>Dur-Terzen nur noch wenig Tremolo</b> besitzen \u2014 aber nicht komplett '
    'Null. Moll-Terzen sollten <b>temperiert</b> (gleichstufig) bleiben, damit der '
    'Moll-Charakter erhalten bleibt.',sB))

fig1=fig_intervalle(); story.append(fig2img(fig1,150)); story.append(Spacer(1,3*mm))

# Tabelle
rows_iv=[]
for nm,(n,d,s) in rein.items():
    cg=s*100; cr=1200*np.log2(n/d); delta=cg-cr
    rows_iv.append([nm,f'{cg:.0f}',f'{cr:.1f}',f'{delta:+.1f}',f'{n}/{d}'])
story.append(mk_tbl(['Intervall','Gleich [Cent]','Rein [Cent]','\u0394 [Cent]','Verh\u00e4ltnis'],rows_iv,
    cw=[28*mm,22*mm,22*mm,20*mm,22*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Das Problem sind die Terzen</b> (\u00b114 Cent), nicht die Quinten (\u00b12 Cent). '
    'Die Quarte und Quinte sind in der gleichstufigen Stimmung fast rein. '
    'Die Terzen und Sexten weichen deutlich ab.',sK))

# Kap. 2
story.append(Paragraph('2. Differenzt\u00f6ne: Was man h\u00f6rt, wenn Terzen erklingen',sCh))
story.append(Paragraph(
    'Wenn zwei T\u00f6ne f\u2081 und f\u2082 gleichzeitig erklingen, erzeugt die Nichtlinearit\u00e4t '
    'des Geh\u00f6rs einen <b>Differenzton</b> bei f_diff = |f\u2082 \u2212 f\u2081|. '
    'Bei einer <b>reinen gro\u00dfen Terz</b> (5/4) ist der Differenzton '
    'f_diff = 5/4 \u00d7 f\u2081 \u2212 f\u2081 = 1/4 \u00d7 f\u2081 \u2014 '
    'exakt zwei Oktaven unter f\u2081. Keine Schwebung, reiner Klang. '
    'Bei der <b>gleichstufigen Terz</b> weicht der Differenzton um '
    '\u22482,6 Hz ab (bei C4) \u2014 das erzeugt eine h\u00f6rbare Schwebung.',sB))
story.append(Paragraph('f_diff = |f\u2082 \u2212 f\u2081|',sF))

fig2=fig_schwebung(); story.append(fig2img(fig2,150)); story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die Schwebung w\u00e4chst proportional zur Grundfrequenz \u2014 '
    'und wird im oberen Bereich zunehmend <b>unangenehm</b>:',sB))

# Tabelle: Schwebung über den Tonbereich
rows_schw=[]
for note,s,o in [('C3',0,3),('C4',0,4),('E4',4,4),('G4',7,4),('C5',0,5),('E5',4,5),('C6',0,6)]:
    f1=fgl(s,o); f2g=f1*2**(4/12); f2r=f1*5/4
    dg=f2g-f1; dr=f2r-f1; schw=abs(dg-dr)
    terz=f'{note}\u2013{["C","C#","D","Eb","E","F","F#","G","G#","A","Bb","B"][(s+4)%12]}{o+(1 if s+4>=12 else 0)}'
    if schw<2: bew='wahrnehmbar'
    elif schw<4: bew='deutlich h\u00f6rbar'
    elif schw<7: bew='ST\u00d6REND'
    elif schw<11: bew='SEHR UNANGENEHM'
    else: bew='EXTREM'
    rows_schw.append([terz,f'{f1:.0f}',f'{dg:.1f}',f'{dr:.1f}',f'{schw:.1f}',bew])
story.append(mk_tbl(['Terz','f\u2081 [Hz]','Diff gl.','Diff rein','Schwebung','Bewertung'],rows_schw,
    cw=[22*mm,18*mm,18*mm,18*mm,20*mm,28*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Ab C5 aufw\u00e4rts werden gleichstufige Terzen als <b>Flimmern</b> wahrgenommen. '
    'Bei C6\u2013E6 betr\u00e4gt die Schwebung \u00fcber 10 Hz \u2014 das klingt wie '
    'ungewolltes Tremolo. <b>Reine Terzen machen im oberen Diskant '
    'den gr\u00f6\u00dften Unterschied.</b>',sK))

fig3=fig_diff_ton(); story.append(fig2img(fig3,150)); story.append(Spacer(1,3*mm))

fig4=fig_diff_uebersicht(); story.append(fig2img(fig4,158)); story.append(Spacer(1,3*mm))

# Große Terz: Vollständige Tabelle über den Tonbereich
story.append(Paragraph('<b>Gro\u00dfe Terz:</b> Differenzt\u00f6ne rein vs. gleichstufig \u2014 '
    'der reine Differenzton ist immer exakt 2 Oktaven unter dem Grundton:',sB))
rows_gt=[]
for s,o in [(0,3),(5,3),(0,4),(4,4),(7,4),(0,5),(4,5),(7,5),(0,6),(4,6)]:
    f1=fgl(s,o); f2g=f1*2**(4/12); f2r=f1*5/4
    dg=f2g-f1; dr=f2r-f1; schw=abs(dg-dr)
    n1=noten[s]+str(o); n2=noten[(s+4)%12]+str(o+(1 if s+4>=12 else 0))
    nd=noten[s]+str(o-2)
    bew=''
    if schw<2: bew='wahrnehmbar'
    elif schw<4: bew='deutlich'
    elif schw<7: bew='st\u00f6rend'
    elif schw<11: bew='UNANGENEHM'
    else: bew='EXTREM'
    rows_gt.append([f'{n1}\u2013{n2}',f'{f1:.0f}',f'{dg:.1f}',f'{dr:.1f}',nd,f'{schw:.1f}',bew])
story.append(mk_tbl(['Terz','f\u2081','Diff gl.','Diff rein','= Note','Schw. [Hz]','Bewertung'],rows_gt,
    cw=[18*mm,16*mm,16*mm,16*mm,14*mm,18*mm,22*mm]))
story.append(Spacer(1,3*mm))

# Dur-Dreiklang
story.append(Paragraph('<b>Dur-Dreiklang C-E-G:</b> Alle Differenzt\u00f6ne im Akkord:',sB))
fig5=fig_akkord(); story.append(fig2img(fig5,150)); story.append(Spacer(1,3*mm))

rows_akk=[]
for o in [3,4,5]:
    fC=fgl(0,o); fE=fgl(4,o); fG=fgl(7,o); fEr=fC*5/4; fGr=fC*3/2
    rows_akk.append([f'C{o}\u2013E{o}',f'{fC:.0f}',f'{fE-fC:.1f}',f'{fEr-fC:.1f}',
        f'{abs((fE-fC)-(fEr-fC)):.1f}',f'= C{o-2}'])
    rows_akk.append([f'E{o}\u2013G{o}',f'{fgl(4,o):.0f}',f'{fG-fE:.1f}',f'{fGr-fEr:.1f}',
        f'{abs((fG-fE)-(fGr-fEr)):.1f}',f'\u2248 G{o-2}'])
    rows_akk.append([f'C{o}\u2013G{o}',f'{fC:.0f}',f'{fG-fC:.1f}',f'{fGr-fC:.1f}',
        f'{abs((fG-fC)-(fGr-fC)):.1f}',f'= C{o-1}'])
story.append(mk_tbl(['Paar','f\u2081','Diff gl.','Diff rein','Schw.','Rein = '],rows_akk,
    cw=[20*mm,16*mm,18*mm,18*mm,16*mm,20*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Im reinen Dreiklang sind <b>alle</b> Differenzt\u00f6ne exakt auf Stufen des Akkords: '
    'C-E \u2192 C (2 Okt. tiefer), C-G \u2192 C (1 Okt. tiefer), E-G \u2192 C (rein). '
    'Der gesamte Akkord klingt wie eine Einheit. '
    'Im gleichstufigen Dreiklang erzeugen die Terzen Schwebungen \u2014 '
    'je h\u00f6her der Akkord, desto schlimmer.',sK))

# Kap. 3
story.append(Paragraph('3. Quinten: Fast rein, kaum Schwebung',sCh))
story.append(Paragraph(
    'Die gleichstufige Quinte weicht nur 2 Cent von der reinen ab. '
    'Der Differenzton bei C-G: f_diff = 130,4 Hz (gleichstufig) vs. '
    '130,8 Hz (rein) \u2014 Schwebung nur 0,44 Hz. Das ist im Tremolo '
    'und bei normaler Spieldynamik nicht wahrnehmbar.',sB))
story.append(Paragraph(
    '<b>Konsequenz:</b> Quinten m\u00fcssen nicht nachgestimmt werden \u2014 '
    'sie sind in der gleichstufigen Stimmung bereits fast rein. '
    'Die Verbesserung durch reine Quintenstimmung betr\u00e4gt < 0,5 Hz '
    'und geht in der Druckabh\u00e4ngigkeit unter (Dok. 0012).',sK))

# Kap. 4
story.append(Paragraph('4. Diatonische Instrumente: M\u00f6glichkeiten und Grenzen',sCh))
story.append(Paragraph(
    'Diatonische Instrumente (Steirische, Wiener) bieten eine besondere M\u00f6glichkeit: '
    'Auf Druck und Zug kommen dieselben Tonstufen vor, und gleiche T\u00f6ne erscheinen '
    'mehrfach in verschiedenen Reihen. Diese Duplikate k\u00f6nnen <b>um Nuancen '
    'verschieden gestimmt</b> werden \u2014 zum Beispiel angelehnt an das Euler-Tonnetz, '
    'wo Terzen rein gestimmt werden.',sB))
story.append(Paragraph(
    'Das Euler-Tonnetz ordnet T\u00f6ne auf zwei Achsen an: Quinten (horizontal) und '
    'Terzen (vertikal). Alle Terzen und Quinten sind rein \u2014 aber gleiche Tonnamen '
    'haben verschiedene Frequenzen (syntontisches Komma = 21,5 Cent). '
    'Ein diatonisches Instrument kann f\u00fcr bestimmte Tonarten eine Ann\u00e4herung '
    'realisieren, indem Duplikate verschieden gestimmt werden.',sB))
story.append(Paragraph(
    '<b>Praxis-Einsch\u00e4tzung:</b> Ein derartiges Instrument mit m\u00f6glichst reinen '
    'Terzen macht kaum Sinn auf einem vierreihigen Instrument. '
    'Dur-Akkorde klingen hervorragend, aber <b>Moll-Akkorde</b> und erst recht '
    '<b>Septimen</b> oder <b>\u00dcbergangst\u00f6ne</b> werden problematisch \u2014 '
    'die reine kleine Terz (6/5) weicht 16 Cent von der gleichstufigen ab, '
    'in die entgegengesetzte Richtung wie die gro\u00dfe Terz. '
    'Wer beide rein stimmt, hat 30 Cent Differenz zwischen gro\u00dfer und '
    'kleiner Terz auf demselben Ton.',sB))
story.append(Paragraph(
    'F\u00fcr <b>Chor\u00e4le</b> und einfache <b>alpenl\u00e4ndische Jodler</b> \u2014 '
    'Musik, die \u00fcberwiegend in Dur-Dreikl\u00e4ngen steht und wenige Moll-Passagen hat \u2014 '
    'ist reine Terzstimmung ein <b>besonderes Klangerlebnis</b>. '
    'Die Differenzt\u00f6ne verschwinden, der Akkord verschmilzt, '
    'das Instrument klingt wie eine Einheit. '
    'F\u00fcr anspruchsvolleres Repertoire mit Moll, Septakkorden und chromatischen '
    '\u00dcberg\u00e4ngen ist die gleichstufige Stimmung der sichere Kompromiss.',sK))

# Kap. 5
story.append(Paragraph('5. Voraussetzungen und Einschr\u00e4nkungen',sCh))
story.append(Paragraph(
    '<b>Einch\u00f6rig = Null-Tremolo:</b> Besondere Stimmpraktiken machen nur Sinn, '
    'wenn das Instrument einch\u00f6rig gestimmt ist \u2014 also <b>kein Tremolo</b>. '
    'Bereits 1 Hz Tremolo erzeugt \u00b11,7 Cent Schwankung bei 440 Hz \u2014 '
    'Feinheiten unter 3 Cent gehen darin unter. Die Terzen-Differenz (14 Cent) '
    'bleibt h\u00f6rbar, die Quinten-Differenz (2 Cent) nicht. '
    'Null-Tremolo ist die Grundvoraussetzung f\u00fcr jede Feinstimmung.',sW))
story.append(Paragraph(
    '<b>Instrumente mit viel Tremolo: Temperiert stimmen!</b> '
    'Bei Instrumenten mit kr\u00e4ftigem Tremolo (Wiener Stimmung, Musette) '
    'sollte die Grundstimmung <b>gleichstufig temperiert</b> bleiben. '
    'Der Grund: Wenn man zus\u00e4tzlich zur Tremolo-Schwebung noch die '
    'Stimmungsabweichung der reinen Terzen hinzuf\u00fcgt, addieren sich die '
    'Schwebungen. Das Tremolo klingt dann <b>ungleichm\u00e4\u00dfig</b> \u2014 '
    'manche Akkorde schweben schneller als andere. Das wirkt verstimmt, '
    'nicht rein. Temperierte Stimmung gew\u00e4hrleistet ein gleichm\u00e4\u00dfiges '
    'Tremolo \u00fcber alle Tonarten.',sW))
story.append(Paragraph(
    'Details zur Tremolo-Physik, Schwebungsberechnung und Wechselwirkung '
    'mit der Stimmung: siehe <b>Dok. 0020</b>.',sB))
story.append(Paragraph(
    '<b>Druckabh\u00e4ngigkeit (Dok. 0012):</b> Der Balgdruck verstimmt den Ton. '
    'Im Bass: \u22485 Cent bei normalem Dynamikwechsel. Das \u00fcberspielt '
    'Quinten-Feinheiten, aber nicht Terzen-Korrekturen. '
    'Im Diskant ist die Druckabh\u00e4ngigkeit gering (\u22480,4 Cent) \u2014 '
    'dort lohnt sich Feinstimmung am meisten.',sB))
story.append(Paragraph(
    '<b>Referenzton:</b> Bei diatonischen Instrumenten ist die <b>2. Reihe auf Druck</b> '
    'die Bezugsreihe \u2014 sie bestimmt die Grundtonart und damit den Referenzton. '
    'Bei G-C-F-B: 2. Reihe = C \u2192 Referenz C. '
    'Bei A-D-G-C: 2. Reihe = D \u2192 Referenz D (nicht A!). '
    'Nur bei einem E-A-D-G Instrument w\u00e4re A = 440 Hz die nat\u00fcrliche Referenz. '
    'Die Wahl einer anderen Reihe als Bezug beim Spielen ist problematisch \u2014 '
    'die 2. Reihe auf Druck ist der Standard.',sB))
story.append(Paragraph(
    '<b>Differenzt\u00f6ne werden im oberen Bereich unangenehm:</b> '
    'Die Schwebung der Differenzt\u00f6ne bei gleichstufigen Terzen verdoppelt sich '
    'pro Oktave. Bei C4\u2013E4: 2,6 Hz (deutlich h\u00f6rbar). '
    'Bei C5\u2013E5: 5,2 Hz (st\u00f6rend). '
    'Bei C6\u2013E6: 10,4 Hz (extrem unangenehm \u2014 fast wie ungewolltes Tremolo). '
    'Ab dem oberen Diskant werden gleichstufige Terzen als <b>Flimmern</b> '
    'wahrgenommen. Reine Terzen machen dort den gr\u00f6\u00dften Unterschied.',sK))
story.append(Paragraph(
    '<b>Jede Stimmung ist ein Kompromiss:</b> Sie sollte auf die Musikst\u00fccke '
    'abgestimmt sein, die tats\u00e4chlich gespielt werden. Reine Terzen in C-Dur '
    'bedeuten unreinere Terzen in anderen Tonarten. '
    'Die Wahl der Stimmung ist eine musikalische Entscheidung, keine physikalische.',sB))

# Kap. 6: Moll-Terz
story.append(Paragraph('6. Die Moll-Terz: 6/5 oder 19/16?',sCh))
story.append(Paragraph(
    'In der Literatur findet man unterschiedliche Angaben zur \u201ereinen\u201c Moll-Terz:',sB))
rows_mt=[]
rows_mt.append(['6/5 (5-Limit)','315,6','+15,6','Klassisch rein, 5. und 6. Oberton'])
rows_mt.append(['19/16 (19-Limit)','297,5','\u22122,5','19. Oberton, nahe an gleichstufig'])
rows_mt.append(['7/6 (septimal)','266,9','\u221233,1','7. Oberton, sehr eng'])
rows_mt.append(['Gleichstufig','300,0','0','Referenz'])
story.append(mk_tbl(['Verh\u00e4ltnis','Cent','vs. gleich','Herkunft'],rows_mt,
    cw=[25*mm,18*mm,18*mm,52*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Die <b>6/5-Terz</b> (315,6 Cent) ist das Gegenst\u00fcck zur gro\u00dfen Terz 5/4 \u2014 '
    'zusammen ergeben sie eine reine Quinte (5/4 \u00d7 6/5 = 3/2). '
    'Sie ist 16 Cent h\u00f6her als gleichstufig. '
    'Die <b>19/16-Terz</b> (297,5 Cent) entspricht dem 19. Oberton \u2014 '
    'nur 2,5 Cent unter der gleichstufigen und klingt f\u00fcr moderne Ohren '
    '\u201erichtig\u201c als Moll.',sB))
story.append(Paragraph(
    'In der Renaissance wurde die Moll-Terz C\u2013Eb oft als gro\u00dfe Terz '
    '<b>unter</b> der Quinte gestimmt (Eb = reine Terz unter G \u2192 316 Cent). '
    'Das erkl\u00e4rt, warum Moll-Dreikl\u00e4nge als Dissonanz galten.',sB))
story.append(Paragraph(
    '<b>Praxis f\u00fcr Begleiter:</b> Moll-Terzen <b>temperiert</b> (gleichstufig) belassen. '
    'Eine reine 6/5-Terz (\u224816 Cent h\u00f6her) ver\u00e4ndert den Moll-Charakter. '
    'Die gleichstufige 300 Cent liegt nahe an 19/16 (297,5 Cent) und wahrt den '
    'vertrauten Moll-Klang. Dur-Terzen beim Begleiter: nicht ganz Null-Tremolo, '
    'sondern mit <b>minimalem Rest-Tremolo</b> \u2014 das mildert die Sch\u00e4rfe, '
    'beh\u00e4lt aber den lebendigen Klang.',sK))

# Kap. 7: Praktische Empfehlungen
story.append(Paragraph('7. Praktische Empfehlungen f\u00fcr den Stimmer',sCh))
story.append(Paragraph(
    '<b>Absolut korrekt stimmen ist nicht m\u00f6glich</b> \u2014 es bleibt immer ein '
    'kleiner Kompromiss. Die folgenden Empfehlungen helfen, den bestm\u00f6glichen '
    'Kompromiss zu finden.',sB))

story.append(Paragraph('<b>Dur-Terzen: Nicht zum Extrem machen</b>',sCh))
story.append(Paragraph(
    'Wenn man die gro\u00dfe Terz absenkt (Richtung rein, 5/4), sollte man sie '
    '<b>nicht ganz auf Null-Schwebung</b> bringen, sondern ein leichtes Tremolo belassen \u2014 '
    'und zwar <b>etwas enger als rein</b>. Warum? '
    'Ein Instrument verstimmt sich im Laufe der Zeit. '
    'Wenn die Terz knapp unter der reinen Stimmung liegt (etwas enger), '
    'geht jede Verstimmung <b>in Richtung gleichstufig</b> \u2014 das Intervall '
    'wird breiter, n\u00e4hert sich 400 Cent und klingt dabei immer noch akzeptabel.',sB))
story.append(Paragraph(
    'Macht man die Terz dagegen <b>zu knapp</b> (zu weit abgesenkt, unter 386 Cent), '
    'besteht die Gefahr, dass Verstimmungen <b>in die falsche Richtung</b> gehen \u2014 '
    'das Intervall wird noch enger, klingt zunehmend \u201efalsch\u201c und '
    'unangenehm. Das sollte <b>unbedingt vermieden</b> werden.',sW))
story.append(Paragraph(
    'Richtwert: Die gro\u00dfe Terz ca. 10\u201312 Cent unter gleichstufig stimmen '
    '(statt der vollen 14 Cent f\u00fcr rein). '
    'Das ergibt \u22482\u20133 Cent Sicherheitsabstand unter der reinen Terz '
    'und l\u00e4sst Raum f\u00fcr Verstimmungen in die richtige Richtung.',sB))

story.append(Paragraph('<b>Quinten: Rein ist m\u00f6glich \u2014 mit Vorbehalt</b>',sCh))
story.append(Paragraph(
    'Komplett reine Quinten (702 Cent) sind nur m\u00f6glich, wenn man gleichzeitig '
    'die <b>Oktaven spreizt</b>. 12 reine Quinten ergeben 7 Oktaven + 23,5 Cent '
    '(pythagor\u00e4isches Komma). Wenn die Oktaven exakt rein bleiben (1200 Cent), '
    'm\u00fcssen die Quinten leicht enger sein (700 Cent gleichstufig). '
    'Der Kompromiss:',sB))

rows_komp=[]
rows_komp.append(['Quinten','eher enger bis rein','698\u2013702 Cent','Sicherheitsabstand nach unten'])
rows_komp.append(['Oktaven','eher rein bis etwas weiter','1200\u20131202 Cent','Kompensiert engere Quinten'])
rows_komp.append(['Gro\u00dfe Terzen','etwas enger als rein','388\u2013390 Cent','Spielraum f\u00fcr Verstimmung'])
rows_komp.append(['Kleine Terzen','temperiert (gleichstufig)','300 Cent','Moll-Charakter bewahren'])
story.append(mk_tbl(['Intervall','Richtung','Zielbereich','Begr\u00fcndung'],rows_komp,
    cw=[25*mm,32*mm,28*mm,40*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Dieses Prinzip \u2014 Quinten eher enger, Oktaven eher weiter \u2014 '
    'ist konsistent und physikalisch begr\u00fcndet. '
    'Es bietet <b>Spielraum f\u00fcr nat\u00fcrliche Verstimmung</b>, ohne dass das '
    'Instrument jemals \u201efalsch\u201c klingt. '
    'Im Gegensatz dazu f\u00fchrt eine \u00fcbertrieben reine Stimmung '
    'bei der kleinsten Verstimmung zu h\u00f6rbaren Problemen.',sK))

story.append(Paragraph('<b>Stimmprozedur nach Geh\u00f6r (ohne Stimmger\u00e4t)</b>',sCh))
story.append(Paragraph(
    'Stimmt man nach Geh\u00f6r, funktioniert das sobald man den '
    '<b>Bezugs-Gleichton in der 2. Reihe</b> festgelegt hat. '
    'Von dort ausgehend erweitert man systematisch:',sB))
story.append(Paragraph(
    '1. Bezugston festlegen \u2014 2. Reihe, Druck (z.B. D bei A-D-G-C Instrument)',sB))
story.append(Paragraph(
    '2. Quinten nach unten und oben erweitern \u2014 eher etwas enger als rein',sB))
story.append(Paragraph(
    '3. Oktaven pr\u00fcfen \u2014 eher rein bis etwas weiter',sB))
story.append(Paragraph(
    '4. Terzen pr\u00fcfen \u2014 leichtes Tremolo muss h\u00f6rbar bleiben '
    '(nicht Null-Schwebung, etwas enger als rein)',sB))
story.append(Paragraph(
    '5. Korrigieren, bis Quinten, Oktaven und Terzen im Gleichgewicht sind',sB))

story.append(Paragraph('<b>Druck rein \u2014 Zug temperiert</b>',sCh))
story.append(Paragraph(
    'Auf einem diatonischen Instrument bietet sich eine elegante Strategie an: '
    'Die <b>Druckt\u00f6ne extrem rein</b> stimmen (Terzen fast rein, Quinten rein), '
    'aber die <b>Zugt\u00f6ne temperiert</b> (gleichstufig) belassen.',sB))
story.append(Paragraph(
    '<b>Warum?</b> Auf Druck erklingen \u00fcberwiegend Dur-Dreikl\u00e4nge und Grundakkorde \u2014 '
    'dort profitiert der Klang maximal von reinen Terzen und Quinten. '
    'Auf Zug erklingen h\u00e4ufig Septimen, Sekunden und Moll-Passagen \u2014 '
    'dort sollen diese Intervalle ihren <b>charakteristischen Spannungsklang</b> '
    'behalten. Temperierte Zugt\u00f6ne gew\u00e4hrleisten das.',sB))
story.append(Paragraph(
    'Dieses Prinzip nutzt die Eigenheit des diatonischen Instruments: '
    'Druck und Zug haben verschiedene Tonvorr\u00e4te. '
    'Auf Druck dominieren Tonika und Dominante (Dur), auf Zug die Subdominante '
    'und Durchgangst\u00f6ne (oft Moll, Septimen). '
    'Zwei verschiedene Stimmungen auf demselben Instrument \u2014 '
    'je nach Balgrichtung der passende Kompromiss.',sK))

story.append(Paragraph('<b>Bass-Grundt\u00f6ne: Grund- und Oktavzunge</b>',sCh))
story.append(Paragraph(
    'Bass-Grundt\u00f6ne bestehen aus einer Grundton-Zunge und einer Oktav-Zunge. '
    'Diese beiden Zungen haben unterschiedliche Profilierungen und damit '
    '<b>unterschiedliche Druckabh\u00e4ngigkeiten</b> (Dok. 0012). '
    'Je nach Schliff und Verhalten der Zungen ist die Basszunge meist '
    'um ca. <b>5 Cent h\u00f6her</b> zu stimmen als die Oktavzunge.',sB))
story.append(Paragraph(
    'Das muss aber <b>nicht</b> sein \u2014 wenn die Profilierungen optimal sind '
    'und die T\u00f6ne bei unterschiedlichem Spieldruck <b>nicht auseinanderdriften</b>, '
    'ist keine Korrektur n\u00f6tig. Der Test: Bei pp und ff spielen und h\u00f6ren, '
    'ob die Schwebung zwischen Grund- und Oktavzunge sich \u00e4ndert. '
    'Wenn ja \u2192 die Zunge mit st\u00e4rkerer Drift etwas in Gegenrichtung vorstimmen. '
    'Wenn nein \u2192 beide auf denselben Ton stimmen.',sB))
story.append(Paragraph(
    'Der Hintergrund: Die l\u00e4ngere, tiefere Grundton-Zunge ist meist <b>weniger steif</b> '
    'als die k\u00fcrzere Oktav-Zunge und hat daher die <b>gr\u00f6\u00dfere '
    'Druckabh\u00e4ngigkeit</b> (Dok. 0012). Bei h\u00f6herem Spieldruck driftet die '
    'Grundton-Zunge st\u00e4rker nach unten als die Oktav-Zunge. '
    'Die Vorstimmung der Grundton-Zunge einige Cent nach oben kompensiert das '
    'f\u00fcr den \u00fcblichen Spieldruckbereich.',sB))
story.append(Paragraph(
    '<b>Spieleranpassung:</b> Je nach Spieler kann das angepasst werden. '
    'Wichtig ist auch die <b>Aufbiegung</b> der Zungen \u2014 sie beeinflusst '
    'die Druckabh\u00e4ngigkeit. Bei normalem Spieldruck sollte der Bass-Grundton '
    'eher <b>einige Cent zu hoch</b> sein, da er bei kr\u00e4ftigem Druck nachgibt. '
    'Ein Spieler mit kr\u00e4ftigem Balg braucht mehr Vorstimmung als ein '
    'leiser Spieler. Der Stimmer sollte den Spieler kennen \u2014 oder '
    'f\u00fcr den mittleren Spieldruck stimmen und dem Spieler die M\u00f6glichkeit '
    'zur Feinjustierung lassen.',sB))

story.append(Paragraph(
    '<b>Anmerkung:</b> Die Empfehlung \u201eMoll-Terzen temperiert\u201c '
    '(Kap. 6) gilt f\u00fcr die <b>Begleiter</b> \u2014 '
    'nicht f\u00fcr den Diskant. Im Diskant richtet sich die Stimmung '
    'nach der Druck/Zug-Strategie (oben).',sW))

# Kap. 8: Vorstimmen
story.append(Paragraph('8. Vorstimmen',sCh))
story.append(Paragraph(
    'Vorgestimmt wird entweder vom <b>Stimmplattenhersteller</b> (der auch nach '
    'Bedarf ventiliert) oder vom Harmonikabauer selbst. '
    'Wird selber vorgestimmt, erfolgt dies <b>bevor</b> die Stimmplatten '
    'aufgewachst werden.',sB))

story.append(Paragraph('<b>Ablauf:</b>',sB))
story.append(Paragraph(
    '1. <b>Sortieren:</b> Kontrollieren, welche Stimmplatten einen '
    'positiven (+) und welche einen negativen (\u2212) Offset besitzen. '
    'Entsprechend sortieren.',sB))
story.append(Paragraph(
    '2. <b>Stimmtisch:</b> Die Stimmplatten auf einem Stimmtisch mit '
    'mehreren aufgereihten Kanzellen, einer Einzelkanzelle oder drei Kanzellen '
    'mit Stimmger\u00e4t vorstimmen. Genauigkeit: ca. <b>\u00b12 Cent</b> '
    '\u2014 den Grundton und die Tremolos.',sB))
story.append(Paragraph(
    '3. <b>Erfahrungswerte als Offset:</b> Dabei sollten bereits die '
    'bekannten Erfahrungswerte als Offset mit einflie\u00dfen \u2014 '
    'Werte, die man kennt, wenn man die Stimmplatten sp\u00e4ter auf dem Stimmstock '
    'und im Instrument nachmisst. Das ist <b>stimmungsabh\u00e4ngig</b>: '
    'Ein B-Es-As-Des-Satz ist anders zu behandeln als ein F-B-Es-As-Satz. '
    'Die Einbau-Verstimmung h\u00e4ngt von der Position auf dem Stimmstock, '
    'der Kammergeometrie und den Verschraubungen ab (Dok. 0021).',sB))
story.append(Paragraph(
    '4. <b>Nicht nur das Stimmger\u00e4t:</b> Man h\u00f6rt genau hin, '
    'wie sich die Stimmplatte verh\u00e4lt \u2014 Anschwingverhalten und Klang. '
    'Immer im Hinterkopf: Eine Stimmplatte, die auf dem Stimmtisch '
    'sauber anspricht, kann im Instrument anders reagieren. '
    'Das Stimmger\u00e4t zeigt die Frequenz, aber nicht die Qualit\u00e4t '
    'der Schwingung.',sK))
story.append(Paragraph(
    '<b>Hohe T\u00f6ne bevorzugt etwas zu hoch vorstimmen:</b> '
    'Bei den kurzen Zungen der hohen T\u00f6ne kommt man beim Feinstimmen '
    'im Instrument schwer zum H\u00f6herstimmen heran (Dok. 0021: Stimmstock '
    'ausbauen oder \u00fcber das Loch dr\u00fccken). Tieferkratzen geht dagegen '
    'problemlos. Daher diese Zungen bevorzugt <b>einige Cent zu hoch</b> '
    'vorstimmen \u2014 die Feinstimmung erfolgt dann nur durch Kratzen '
    'nach unten.',sB))

# Kap. 9
story.append(Paragraph('9. Zusammenfassung',sCh))
pts=[
    '<b>Terzen sind das Problem:</b> \u00b114 Cent Abweichung in der gleichstufigen Stimmung. '
    'Quinten sind fast rein (\u00b12 Cent).',
    '<b>Differenzt\u00f6ne bei Terzen werden im oberen Bereich unangenehm:</b> '
    'Schwebung 2,6 Hz bei C4, 5,2 Hz bei C5, 10,4 Hz bei C6 \u2014 '
    'verdoppelt sich pro Oktave. Ab C5 st\u00f6rend, ab C6 extrem unangenehm.',
    '<b>Reine Terzen:</b> Der Differenzton f_diff = 1/4 \u00d7 f\u2081 '
    'ist schwebungsfrei. Reine Stimmung lohnt sich besonders im oberen Diskant.',
    '<b>Null-Tremolo ist Voraussetzung:</b> Bereits 1 Hz Tremolo verwischt '
    'Feinheiten < 3 Cent. Nur einch\u00f6rig gestimmte Instrumente profitieren.',
    '<b>Referenzton = 2. Reihe auf Druck:</b> G-C-F-B \u2192 Referenz C. '
    'A-D-G-C \u2192 Referenz D (nicht A!). Nur E-A-D-G \u2192 Referenz A.',
    '<b>Druckabh\u00e4ngigkeit begrenzt die Genauigkeit</b> im Bass (Dok. 0012). '
    'Feinstimmung lohnt sich im Diskant mehr.',
    '<b>Jede Stimmung ist ein Kompromiss</b> und sollte auf das Repertoire '
    'abgestimmt sein.',
    '<b>Praxis:</b> Terzen etwas enger als rein (Spielraum f\u00fcr Verstimmung in '
    'die richtige Richtung). Quinten enger bis rein, Oktaven rein bis weiter. '
    'Nie zum Extrem \u2014 Verstimmungen d\u00fcrfen nicht in die falsche Richtung gehen.',
    '<b>Vorstimmen:</b> Vor dem Aufwachsen, auf \u00b12 Cent genau. '
    'Erfahrungswerte f\u00fcr Einbau-Offset mit einbeziehen (stimmungsabh\u00e4ngig). '
    'Nicht nur Stimmger\u00e4t \u2014 auch H\u00f6ren: Anschwingverhalten und Klang.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Reine Terzen klingen rein \u2014 aber nur in der Tonart, '
    'f\u00fcr die sie gestimmt sind. Der kluge Stimmer l\u00e4sst Spielraum: '
    'etwas enger als rein, damit Verstimmungen in die richtige Richtung gehen.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0019 \u2014 Stimmung \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
