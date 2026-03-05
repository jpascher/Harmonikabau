#!/usr/bin/env python3
"""Dok 0008 v2: Klangveränderung durch Kammergeometrie — von Grund auf korrekt"""

import math, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.enums import TA_CENTER, TA_LEFT, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    PageBreak, HRFlowable, Image
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

DARKBLUE=HexColor('#16213e'); ACCENTRED=HexColor('#e94560')
KEYGREEN=HexColor('#2e7d32'); WARNRED=HexColor('#c62828')
LIGHTGREEN=HexColor('#e8f5e9'); LIGHTWARN=HexColor('#ffebee')
LIGHTGRAY=HexColor('#f5f5f5'); MIDGRAY=HexColor('#9e9e9e')
LIGHTBLUE=HexColor('#e8eaf6'); LIGHTORANGE=HexColor('#fff3e0')
CRITRED=HexColor('#ffcdd2'); WARNORG=HexColor('#fff3e0')

try:
    pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuBold','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuItalic','/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
    FONT='DejaVu'; FONTB='DejaVuBold'; FONTI='DejaVuItalic'
except:
    FONT='Helvetica'; FONTB='Helvetica-Bold'; FONTI='Helvetica-Oblique'

outpath='/mnt/user-data/outputs/0008_klangveraenderung_De.pdf'
doc=SimpleDocTemplate(outpath,pagesize=A4,leftMargin=18*mm,rightMargin=18*mm,topMargin=25*mm,bottomMargin=18*mm)
styles=getSampleStyleSheet()
def ms(name,**kw):
    return ParagraphStyle(name,parent=styles['Normal'],fontName=kw.get('fn',FONT),
        fontSize=kw.get('fs',10),leading=kw.get('ld',14),spaceAfter=kw.get('sa',6),
        spaceBefore=kw.get('sb',0),alignment=kw.get('al',TA_JUSTIFY),
        textColor=kw.get('tc',black),leftIndent=kw.get('li',0),rightIndent=kw.get('ri',0))
sTitle=ms('sT',fn=FONTB,fs=16,al=TA_CENTER,tc=DARKBLUE,sa=4,ld=20)
sSubtitle=ms('sSt',fs=10.5,al=TA_CENTER,tc=DARKBLUE,sa=2,ld=14)
sNote=ms('sN',fs=8.5,fn=FONTI,tc=MIDGRAY,sa=8,ld=11)
sSec=ms('sSec',fs=13,fn=FONTB,tc=DARKBLUE,sb=12,sa=6,ld=17)
sSub=ms('sSub',fs=11,fn=FONTB,tc=HexColor('#2a4a7f'),sb=8,sa=4,ld=14)
sBody=ms('sB',fs=10,ld=14,sa=6)
sKey=ms('sK',fs=9.5,ld=13,sa=4,li=4*mm,ri=4*mm)
sWarn=ms('sW',fs=9.5,ld=13,sa=4,li=4*mm,ri=4*mm)
sF=ms('sF',fs=10,al=TA_CENTER,fn=FONTI,sb=4,sa=6,ld=14)
sCap=ms('sC',fs=8.5,fn=FONTI,al=TA_CENTER,sa=8,ld=11)
sCell=ms('sCell',fs=7.5,ld=10,sa=0,al=TA_LEFT)
sCellB=ms('sCellB',fs=7.5,fn=FONTB,ld=10,sa=0,al=TA_LEFT)
sCellC=ms('sCellC',fs=7.5,ld=10,sa=0,al=TA_CENTER)
sCellCB=ms('sCellCB',fs=7.5,fn=FONTB,ld=10,sa=0,al=TA_CENTER)

story=[]
def section(n,t): story.append(Paragraph(f'Kapitel {n}: {t}',sSec))
def subsection(t): story.append(Paragraph(t,sSub))
def body(t): story.append(Paragraph(t,sBody))
def formula(t): story.append(Paragraph(t,sF))
def keybox(t):
    tbl=Table([[Paragraph(t,sKey)]],colWidths=[doc.width-8*mm])
    tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),LIGHTGREEN),('BOX',(0,0),(-1,-1),1,KEYGREEN),
        ('TOPPADDING',(0,0),(-1,-1),3*mm),('BOTTOMPADDING',(0,0),(-1,-1),3*mm),
        ('LEFTPADDING',(0,0),(-1,-1),4*mm),('RIGHTPADDING',(0,0),(-1,-1),4*mm)]))
    story.append(tbl); story.append(Spacer(1,3*mm))
def warnbox(t):
    tbl=Table([[Paragraph(t,sWarn)]],colWidths=[doc.width-8*mm])
    tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,-1),LIGHTWARN),('BOX',(0,0),(-1,-1),1,WARNRED),
        ('TOPPADDING',(0,0),(-1,-1),3*mm),('BOTTOMPADDING',(0,0),(-1,-1),3*mm),
        ('LEFTPADDING',(0,0),(-1,-1),4*mm),('RIGHTPADDING',(0,0),(-1,-1),4*mm)]))
    story.append(tbl); story.append(Spacer(1,3*mm))


# ======================= PHYSICS =======================
c=343.0; Q_H=7.0

def notch_amplitude(f, f_H, Q_H=7.0):
    r=f/f_H; D=(1-r**2)**2+(r/Q_H)**2
    R_norm=(r/Q_H)/D; kappa_eff=0.05
    return 1.0/(1.0+kappa_eff*R_norm)

def notch_phase(f, f_H, Q_H=7.0):
    r=f/f_H
    return math.degrees(math.atan(Q_H*(r-1/r)))


# === CHARTS ===
f_tone=440.0
fH_values=[1200,1700,2800]
fH_colors=['#c62828','#e65100','#2e7d32']
fH_labels=['f_H=1200 (Öffn.↓)','f_H=1700 (Referenz)','f_H=2800 (Tiefe↓)']

fig,axes=plt.subplots(2,2,figsize=(8,7.5),dpi=150)
fig.suptitle('Dok. 0008: Die Kammer als Kerbfilter im Obertonspektrum',
             fontsize=12,fontweight='bold',color='#16213e')

freqs_cont=np.linspace(100,5000,1000)

ax=axes[0,0]
for fH,col,lbl in zip(fH_values,fH_colors,fH_labels):
    amps=[notch_amplitude(f,fH) for f in freqs_cont]
    ax.plot(freqs_cont,amps,'-',color=col,lw=1.5,label=lbl,alpha=0.8)
for n in range(1,12):
    f_ot=n*f_tone
    if f_ot<5000:
        ax.axvline(x=f_ot,color='gray',ls=':',lw=0.5,alpha=0.4)
        ax.text(f_ot,1.02,f'{n}',fontsize=6,ha='center',color='gray')
ax.set_xlabel('Frequenz [Hz]',fontsize=9); ax.set_ylabel('Amplitude',fontsize=9)
ax.set_title('Kerbfilter-Wirkung auf A4-Obertöne',fontsize=10,color='#16213e')
ax.legend(fontsize=6.5,loc='lower right'); ax.grid(True,alpha=0.3)
ax.set_ylim(0.7,1.05); ax.set_xlim(100,5000)

ax=axes[0,1]
for fH,col,lbl in zip(fH_values,fH_colors,fH_labels):
    phases=[notch_phase(f,fH) for f in freqs_cont]
    ax.plot(freqs_cont,phases,'-',color=col,lw=1.5,label=lbl,alpha=0.8)
for n in range(1,12):
    f_ot=n*f_tone
    if f_ot<5000: ax.axvline(x=f_ot,color='gray',ls=':',lw=0.5,alpha=0.4)
ax.axhline(y=0,color='black',lw=0.5,alpha=0.3)
ax.set_xlabel('Frequenz [Hz]',fontsize=9); ax.set_ylabel('Phase [°]',fontsize=9)
ax.set_title('Phasenverschiebung der Obertöne',fontsize=10,color='#16213e')
ax.legend(fontsize=6.5); ax.grid(True,alpha=0.3)
ax.set_ylim(-95,95); ax.set_xlim(100,5000)

ax=axes[1,0]
n_max=10; x=np.arange(1,n_max+1); width=0.25
a_free=[1.0/n for n in x]
for j,(fH,col,lbl) in enumerate(zip(fH_values,fH_colors,fH_labels)):
    a_ch=[a_free[i]*notch_amplitude((i+1)*f_tone,fH) for i in range(n_max)]
    ax.bar(x+(j-1)*width,a_ch,width,color=col,alpha=0.7,label=lbl)
ax.plot(x,a_free,'k--',lw=1,alpha=0.4,label='Frei (ohne Kammer)')
ax.set_xlabel('Oberton n',fontsize=9); ax.set_ylabel('Amplitude a_n',fontsize=9)
ax.set_title('Obertonspektrum A4',fontsize=10,color='#16213e')
ax.legend(fontsize=6,loc='upper right'); ax.grid(True,alpha=0.3,axis='y'); ax.set_xticks(x)

ax=axes[1,1]
for j,(fH,col,lbl) in enumerate(zip(fH_values,fH_colors,fH_labels)):
    phases_ot=[notch_phase(n*f_tone,fH) for n in x]
    ax.bar(x+(j-1)*width,phases_ot,width,color=col,alpha=0.7,label=lbl)
ax.axhline(y=0,color='black',lw=0.5)
ax.set_xlabel('Oberton n',fontsize=9); ax.set_ylabel('Phase [°]',fontsize=9)
ax.set_title('Phasenverschiebung pro Oberton',fontsize=10,color='#16213e')
ax.legend(fontsize=6); ax.grid(True,alpha=0.3,axis='y'); ax.set_xticks(x)

plt.tight_layout(rect=[0,0,1,0.95])
chart1='/home/claude/chart_0008_v2_notch.png'
plt.savefig(chart1,dpi=150,bbox_inches='tight'); plt.close()

# Chart 2: Klangcharakter over range
fig,axes=plt.subplots(2,1,figsize=(8,5.5),dpi=150)
fig.suptitle('Klangcharakter über den Tonbereich — Variante A vs. C',fontsize=12,fontweight='bold',color='#16213e')

notes_freq=[440*2**((m-69)/12) for m in range(50,85)]
nn=['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
notes_name=[nn[m%12]+str(m//12-1) for m in range(50,85)]
W=15.5e-3;L_max_m=50e-3;L_min_m=30e-3;dep8=8e-3;dep3=3e-3;d10=10e-3
S_neck=math.pi*(d10/2)**2;l_eff=8e-3+0.85*d10

def klang_analysis(freq,fH):
    worst_loss=0;worst_n=0;worst_phase=0
    for n in range(1,11):
        fot=n*freq;amp=notch_amplitude(fot,fH)
        loss_dB=-20*math.log10(max(amp,0.01));phase=abs(notch_phase(fot,fH))
        if loss_dB>worst_loss: worst_loss=loss_dB;worst_n=n;worst_phase=phase
    return worst_n,worst_loss,worst_phase

losses_A=[];losses_C=[];affected_A=[];affected_C=[];phases_A=[];phases_C=[]
for i,freq in enumerate(notes_freq):
    frac=i/34;L=L_max_m-frac*(L_max_m-L_min_m)
    V_A=W*L*dep8;fH_A=(c/(2*math.pi))*math.sqrt(S_neck/(V_A*l_eff))
    nA,lA,pA=klang_analysis(freq,fH_A);losses_A.append(lA);affected_A.append(nA);phases_A.append(pA)
    dC=dep8*(dep3/dep8)**frac;V_C=W*L*dC;fH_C=(c/(2*math.pi))*math.sqrt(S_neck/(V_C*l_eff))
    nC,lC,pC=klang_analysis(freq,fH_C);losses_C.append(lC);affected_C.append(nC);phases_C.append(pC)

ax=axes[0]
ax.bar([f-4 for f in notes_freq],losses_A,width=8,color='#c62828',alpha=0.6,label='A: Ref. (8mm)')
ax.bar([f+4 for f in notes_freq],losses_C,width=8,color='#2e7d32',alpha=0.6,label='C: konvex (8→3)')
ax.set_ylabel('Max. OT-Dämpfung [dB]',fontsize=9)
ax.set_title('Stärkste Obertondämpfung pro Ton',fontsize=10,color='#16213e')
ax.legend(fontsize=7);ax.grid(True,alpha=0.3,axis='y');ax.set_xlim(notes_freq[0]-15,notes_freq[-1]+15)
ax.set_xticklabels([])
ax2=ax.twinx()
ax2.plot(notes_freq,affected_A,'r--',lw=0.8,alpha=0.5);ax2.plot(notes_freq,affected_C,'g--',lw=0.8,alpha=0.5)
ax2.set_ylabel('Betroffener OT n',fontsize=8,color='gray')

ax=axes[1]
ax.bar([f-4 for f in notes_freq],phases_A,width=8,color='#c62828',alpha=0.6,label='A')
ax.bar([f+4 for f in notes_freq],phases_C,width=8,color='#2e7d32',alpha=0.6,label='C')
ax.set_ylabel('Max. Phasenverschiebung [°]',fontsize=9);ax.set_xlabel('Tonfrequenz [Hz]',fontsize=9)
ax.set_title('Stärkste Phasenverschiebung pro Ton',fontsize=10,color='#16213e')
ax.legend(fontsize=7);ax.grid(True,alpha=0.3,axis='y');ax.set_xlim(notes_freq[0]-15,notes_freq[-1]+15)
tpos=notes_freq[::3];tlbl=notes_name[::3]
ax.set_xticks(tpos);ax.set_xticklabels(tlbl,fontsize=6,rotation=45)

plt.tight_layout(rect=[0,0,1,0.95])
chart2='/home/claude/chart_0008_v2_klang.png'
plt.savefig(chart2,dpi=150,bbox_inches='tight');plt.close()


# ======================= PDF =======================
story.append(Paragraph('Dok. 0008: Klangveränderung durch Kammergeometrie',sTitle))
story.append(Paragraph('Die Kammer als Kerbfilter im Obertonspektrum',sSubtitle))
story.append(HRFlowable(width='80%',thickness=2,color=ACCENTRED,spaceAfter=3,spaceBefore=3))
story.append(Paragraph(
    'Analysiert, wie die Kammergeometrie (Volumen, Öffnung) die Klangfarbe des Tons '
    'verändert. Die Kammer wirkt als frequenzabhängiger Kerbfilter (Notch): Sie dämpft '
    'Obertöne nahe f_H und verschiebt deren Phase. Beide Effekte verändern die Wellenform '
    'und damit den Klangeindruck.',sNote))

story.append(Spacer(1,2*mm))
ref_data=[
    ['Dok. 0002','Strömungsanalyse Bass-Stimmzunge 50 Hz — v8'],
    ['Dok. 0005','Frequenzverschiebung als Indikator der Ansprache'],
    ['Dok. 0007','Diskant-Stimmstock — Kammerfrequenzen (7 Varianten)'],
    ['<b>Dok. 0008</b>','<b>Klangveränderung durch Kammergeometrie (dieses Dokument)</b>'],
]
ref_tbl=Table([[Paragraph(r[0],sCell),Paragraph(r[1],sCell)] for r in ref_data],
              colWidths=[25*mm,doc.width-25*mm-4*mm])
ref_tbl.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm),('BACKGROUND',(0,3),(-1,3),LIGHTGREEN)]))
story.append(ref_tbl)
story.append(Spacer(1,4*mm))


# ==========================================================================
section(1,'Die Kammer als Kerbfilter — nicht als Verstärker')
# ==========================================================================

body('Ein verbreitetes Missverständnis: Die Kammer „verstärkt" den Ton, wie der Resonanzkörper '
     'einer Gitarre. <b>Das Gegenteil ist der Fall.</b> Dok. 0002 (Kap. 10) und Dok. 0003 zeigen: '
     'Wegen der <b>Serienkopplung</b> entzieht die Kammer der Zunge bei Resonanz Energie — '
     'sie dämpft, statt zu verstärken.')

body('Die Kammer wirkt auf das Obertonspektrum der Zunge wie ein <b>Kerbfilter</b> (Notch-Filter): '
     'Sie erzeugt ein „Loch" im Spektrum bei der Kammerresonanz f<sub>H</sub>. '
     'Die Übertragungsfunktion für den n-ten Oberton bei Frequenz f<sub>n</sub> = n · f:')

formula('A<sub>n</sub>(f<sub>H</sub>) = 1 / (1 + κ<sub>eff</sub> · R<sub>H</sub>(f<sub>n</sub>))')

body('Bei f<sub>n</sub> = f<sub>H</sub> ist R<sub>H</sub> maximal (= Q<sub>H</sub>) und die '
     'Dämpfung am stärksten. Weit entfernt von f<sub>H</sub> ist R<sub>H</sub> ≈ 0 '
     'und A<sub>n</sub> ≈ 1 (ungedämpft).')

keybox('<b>Die Kammer ist kein Resonanzkörper, sondern ein Kerbfilter.</b> Sie erzeugt ein Loch '
       'im Obertonspektrum bei f<sub>H</sub>. Die Zunge selbst erzeugt alle Obertöne — die Kammer '
       'filtert einige davon heraus. Welche betroffen sind, hängt von f<sub>H</sub> ab '
       '(Volumen und Öffnung, Dok. 0007). Wie stark sie betroffen sind, hängt von der '
       'Kopplungsstärke κ<sub>eff</sub> ab.')


# ==========================================================================
section(2,'Amplitudendämpfung und Phasenverschiebung')
# ==========================================================================

body('Dok. 0005 zeigt: Der Realteil R<sub>H</sub>(f) und der Imaginärteil X<sub>H</sub>(f) '
     'der Kammerimpedanz sind durch Kramers-Kronig verknüpft. Für den Klang bedeutet das: '
     '<b>Jede Amplitudendämpfung geht mit einer Phasenverschiebung einher</b> — und umgekehrt.')

subsection('Amplitudendämpfung (Realteil R_H)')

body('R<sub>H</sub>(f) hat die Form einer <b>Lorentz-Glocke</b>, zentriert bei f<sub>H</sub> '
     'mit Halbwertsbreite f<sub>H</sub>/Q<sub>H</sub>. Obertöne innerhalb dieser Breite werden '
     'in ihrer Amplitude reduziert. Musikalisch:')

data_klang=[
    ['Betroffener OT','Frequenzbereich','Klangeffekt'],
    ['n = 2 (Oktave)','2·f','Klang wird <b>dünn, hohl</b> — 25 % der Energie fehlt'],
    ['n = 3 (Duodezime)','3·f','Klang wird <b>nasal, kehlig</b> — 11 % fehlt'],
    ['n = 4 (Doppeloktave)','4·f','Klang verliert <b>Klarheit</b> — 6 % fehlt'],
    ['n = 5–6','5–6·f','Klang verliert <b>Brillanz</b> — wenig auffällig'],
    ['n ≥ 7','≥ 7·f','<b>Kaum hörbar</b> — Amplitude ohnehin < 2 %'],
]
dw=[[Paragraph(c,sCellCB if i==0 else sCellC) for c in row] for i,row in enumerate(data_klang)]
tbl=Table(dw,colWidths=[28*mm,28*mm,75*mm],repeatRows=1)
tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHTBLUE),('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('BACKGROUND',(0,1),(-1,1),CRITRED),('BACKGROUND',(0,2),(-1,2),WARNORG)]))
story.append(tbl)
story.append(Paragraph('Tabelle 1: Klangeffekt nach Oberton-Ordnung. Rot: kritisch. Orange: hörbar.',sCap))

subsection('Phasenverschiebung (Imaginärteil X_H)')

body('X<sub>H</sub>(f) hat die Form einer <b>S-Kurve</b> (Dispersionskurve): '
     'Unterhalb von f<sub>H</sub> wird die Phase <b>nacheilend</b> verschoben, '
     'oberhalb <b>voreilend</b>. Am stärksten knapp neben f<sub>H</sub>.')

formula('Φ<sub>n</sub> = arctan(Q<sub>H</sub> · (f<sub>n</sub>/f<sub>H</sub> − f<sub>H</sub>/f<sub>n</sub>))')

body('Die Phasenverschiebung verändert die <b>Wellenform</b> des Tons, ohne die '
     'Einzelamplituden wesentlich zu ändern. Das Ohr nimmt das als Veränderung '
     'der <b>Schärfe</b> oder <b>Weichheit</b> wahr. Der Effekt ist am stärksten, '
     'wenn benachbarte Obertöne <b>unterschiedliche</b> Phasenverschiebungen erfahren — '
     'das verformt die Wellenform am deutlichsten.')

warnbox('<b>Die Phasenverschiebung erklärt, warum erfahrene Instrumentenbauer '
        'Klangunterschiede hören, die ein Spektralanalysator kaum zeigt.</b> '
        'Die Amplituden ändern sich nur um wenige Prozent, aber die Phasenbeziehungen '
        'ändern sich deutlich — und das formt die Wellenform um.')


# ==========================================================================
section(3,'Diagramme')
# ==========================================================================

if os.path.exists(chart1):
    story.append(Image(chart1,width=170*mm,height=155*mm))
    story.append(Paragraph(
        'Abbildung 1: Kerbfilter-Wirkung am Beispiel A4 (440 Hz). '
        'Oben links: Amplitudendämpfung — drei verschiedene f_H. Senkrechte Linien = Obertöne. '
        'Oben rechts: Phasenverschiebung (S-Kurve mit Nulldurchgang bei f_H). '
        'Unten links: Obertonspektrum für drei Varianten. Gestrichelt: freie Zunge ohne Kammer. '
        'Unten rechts: Phasenverschiebung pro Oberton.',sCap))


# ==========================================================================
story.append(PageBreak())
section(4,'Was die Variationsfaktoren am Klang verändern')
# ==========================================================================

body('Die sieben Varianten aus Dok. 0007 verschieben f<sub>H</sub> — und damit die Position '
     'des Kerbfilters im Obertonspektrum:')

subsection('Tiefe senken (Varianten B, C, Ci): f_H steigt → Loch wandert nach oben')

body('Wenn die Tiefe von 8 auf 3 mm sinkt, steigt f<sub>H</sub> von 1500–1950 Hz auf '
     '1500–3190 Hz. Das Kerbfilter-Loch verschiebt sich von den 3.–5. Obertönen '
     'nach oben zu den 5.–8. Obertönen.')

body('<b>Klangeffekt:</b> Die energiereichen niedrigen Obertöne (2., 3.) werden <b>freigegeben</b>. '
     'Der Klang wird <b>voller, grundtöniger, wärmer</b>. Die Brillanz (hohe Obertöne ab dem 5.) '
     'wird leicht gedämpft — aber diese Obertöne haben ohnehin wenig Energie (a<sub>5</sub>² = 4 %). '
     'Der Verlust ist kaum hörbar.')

body('Innerhalb der Tiefengruppe: C (konvex) schiebt f<sub>H</sub> im Mittelbereich am weitesten '
     'nach oben → am vollsten. Ci (konkav) am wenigsten → am dünnsten der drei. '
     'Der Unterschied ist <b>subtil</b> — alle drei klingen deutlich voller als die Referenz A.')

subsection('Öffnung verkleinern (Varianten D, E, Ei): f_H sinkt → Loch wandert nach unten')

body('Wenn ∅ von 10 auf 6 mm sinkt, fällt f<sub>H</sub> auf 67 %. Das Kerbfilter-Loch wandert '
     '<b>nach unten</b> — in den Bereich der 2.–3. Obertöne.')

body('<b>Klangeffekt:</b> Die energiereichen niedrigen Obertöne werden <b>gedämpft</b>. '
     'Der Klang wird <b>dünn, hohl, leblos</b>. Der 2. Oberton (Oktave) trägt 25 % der '
     'Gesamtenergie — wenn er gedämpft wird, verliert der Ton seinen „Körper". '
     'Das ist nicht nur eine Ansprache-Verschlechterung (Dok. 0005), sondern auch '
     'eine <b>Klangverarmung</b>.')

subsection('Q_H: Die Breite des Kerbfilters')

body('Die Güte Q<sub>H</sub> bestimmt die Breite des Kerbfilter-Lochs. '
     'Die quantitative Analyse zeigt, dass Q<sub>H</sub> hauptsächlich durch die '
     '<b>Strahlungsverluste an der Öffnung</b> bestimmt wird — nicht durch die Wandform:')

data_qh=[
    ['','Q_rad (Strahlung)','Q_visc (Wandreibung)','Q_total','Strahlung dominiert?'],
    ['Diskant (5 cm³, 1691 Hz)','43','88','29','Ja (67 %)'],
    ['Bass (180 cm³, 274 Hz)','263','113','79','Nein (30 %) — Wand 70 %'],
]
dw=[[Paragraph(c,sCellCB if i==0 else sCellC) for c in row] for i,row in enumerate(data_qh)]
tbl=Table(dw,colWidths=[35*mm,28*mm,28*mm,18*mm,28*mm],repeatRows=1)
tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHTBLUE),('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
story.append(tbl)
story.append(Paragraph('Tabelle 2: Zerlegung von Q_H. Die Wandform ändert Q_visc um ±5–10 % (über die Fläche), '
    'was Q_total um maximal ±3 % (Diskant) bis ±7 % (Bass) verschiebt — nicht hörbar.',sCap))

body('Die Wandform (gerade, schräg, parabolisch) ändert die Kammer-Innenfläche um maximal '
     '±5–10 %. Da Q<sub>visc</sub> nur einen Teil der Gesamtverluste ausmacht, '
     'verschiebt das Q<sub>H</sub> um höchstens <b>±3 % (Diskant)</b> bis <b>±7 % (Bass)</b>. '
     'Die Kerbfilter-Breite ändert sich dabei um wenige Hz — <b>nicht hörbar</b>.')

body('Musikalisch wirkt Q<sub>H</sub> trotzdem, wenn auch nicht über die Wandform einstellbar:')
body('<b>Hohes Q<sub>H</sub></b> (große Öffnung, wenig Strahlung): <b>Schmales, tiefes Loch.</b> '
     'Ein einzelner Oberton wird stark gedämpft, die Nachbarn bleiben frei. '
     'Musikalisch: Eine <b>scharfe Klangfärbung</b>.')
body('<b>Niedriges Q<sub>H</sub></b> (kleine Öffnung, viel Strahlung): <b>Breites, flaches Loch.</b> '
     'Mehrere Obertöne werden leicht gedämpft. '
     'Musikalisch: Eine <b>diffuse Mattheit</b>.')

keybox('<b>Q<sub>H</sub> wird fast ausschließlich durch die Öffnungsgeometrie bestimmt</b> '
       '(Durchmesser, Halslänge), nicht durch die Wandform. Der hörbare Unterschied zwischen '
       'verschiedenen Trennwandformen (Dok. 0002 Praxisbefund) läuft nicht über Q<sub>H</sub>, '
       'sondern über die <b>instationäre Impulsantwort</b> — Reflexionsmuster, Position des '
       'akustischen Endes, Druckstoß-Transport. Diese Effekte sind real, aber nicht im '
       'Helmholtz-Modell berechenbar.')


# ==========================================================================
section(5,'Klangcharakter über den Tonbereich')
# ==========================================================================

if os.path.exists(chart2):
    story.append(Image(chart2,width=170*mm,height=120*mm))
    story.append(Paragraph(
        'Abbildung 2: Oben: Stärkste Obertondämpfung pro Ton (dB) für Referenz A (rot) und '
        'Variante C (grün). Gestrichelt: betroffener Oberton n. '
        'Unten: Maximale Phasenverschiebung. Variante C hat durchweg weniger Dämpfung.',sCap))

body('<b>Unteres Drittel (D3–C#4):</b> f<sub>H</sub> weit über den niedrigen Obertönen. '
     'Kammer „unsichtbar". Klang voll und natürlich. Kein Unterschied zwischen A und C.')

body('<b>Mittleres Drittel (D4–C#5):</b> Bei A kommt der 3.–5. Oberton in f<sub>H</sub>-Nähe. '
     'Leichte Dämpfung — Klang wird etwas „nasaler". Bei C ist die Dämpfung geringer, '
     'weil f<sub>H</sub> höher liegt.')

body('<b>Oberes Drittel (D5–C6):</b> Bei A trifft der 2. Oberton f<sub>H</sub> — das Loch '
     'sitzt bei der Oktave. Klang wird <b>deutlich dünner und hohler</b>. Bei C liegt f<sub>H</sub> '
     'beim 3.–4. Oberton — Klang bleibt voller, verliert nur etwas Brillanz.')

keybox('<b>Der Klang-Vorteil von Variante C konzentriert sich auf das obere Drittel.</b> '
       'Dort, wo A die 2. Obertöne dämpft (25 % Energie), lässt C sie frei und dämpft nur '
       'die schwächeren 4.–5. Obertöne (4–6 % Energie). In den unteren zwei Dritteln '
       'klingen beide gleich.')


# ==========================================================================
section(6,'Die drei Klangdimensionen')
# ==========================================================================

body('Drei voneinander unterscheidbare Dimensionen wirken auf den Klang:')

data_dim=[
    ['Dimension','Physik','Bestimmt durch','Klangeindruck'],
    ['Position des Lochs','f_H — Zentrum des Kerbfilters',
     'Kammervolumen V, Öffnung S (Dok. 0007)',
     'Welcher OT fehlt → dünn/voll, nasal/klar'],
    ['Breite des Lochs','Q_H — Schärfe des Kerbfilters',
     'Öffnungsgeometrie (67–70 %). Wandform < 7 %.',
     'Scharf (1 OT weg) vs. diffus (mehrere matt)'],
    ['Tiefe des Lochs','κ_eff — Stärke der Kopplung',
     'Spaltfläche A_eff, Volumen V',
     'Starke Färbung vs. kaum hörbar'],
]
dw=[[Paragraph(c,sCellCB if i==0 else sCell) for c in row] for i,row in enumerate(data_dim)]
tbl=Table(dw,colWidths=[27*mm,33*mm,38*mm,38*mm],repeatRows=1)
tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),DARKBLUE),('TEXTCOLOR',(0,0),(-1,0),white),
    ('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),2*mm),('BOTTOMPADDING',(0,0),(-1,-1),2*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm),('VALIGN',(0,0),(-1,-1),'TOP')]+
    [('BACKGROUND',(0,i),(-1,i),LIGHTGRAY) for i in range(2,4,2)]))
story.append(tbl)
story.append(Paragraph('Tabelle 3: Die drei Klangdimensionen. Position ist der dominante Hebel.',sCap))

body('Die <b>Position</b> (f<sub>H</sub>) ist der mit Abstand wirksamste Hebel und über das '
     'Kammervolumen direkt einstellbar (Dok. 0007). Die <b>Breite</b> (Q<sub>H</sub>) hängt von '
     'der Öffnungsgeometrie ab, nicht von der Wandform. Die <b>Tiefe</b> (κ<sub>eff</sub>) wird '
     'primär durch die Zunge bestimmt.')


# ==========================================================================
section(7,'Klangbilanz der Varianten aus Dok. 0007')
# ==========================================================================

data_var=[
    ['Variante','f_H-Bereich','Loch bei','Klangcharakter'],
    ['A (Ref. 8 mm)','1513–1953 Hz','OT 3–5, OT 2 bei hohen Tönen',
     'Unten/Mitte: natürlich. Oben: dünn, hohl'],
    ['C (Tiefe konvex)','1513–3189 Hz','OT 5–8, OT 3–4 bei hohen Tönen',
     'Über gesamten Bereich: voll, warm'],
    ['D (∅ lin. 10→6)','1513–1315 Hz','OT 2–3, bei hohen unter Grundton',
     'Dünn, leblos, hohl — „Körper" fehlt'],
]
dw=[[Paragraph(c,sCellCB if i==0 else sCell) for c in row] for i,row in enumerate(data_var)]
tbl=Table(dw,colWidths=[28*mm,28*mm,35*mm,45*mm],repeatRows=1)
tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHTBLUE),('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),2*mm),('BOTTOMPADDING',(0,0),(-1,-1),2*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm),('VALIGN',(0,0),(-1,-1),'TOP'),
    ('BACKGROUND',(0,2),(-1,2),LIGHTGREEN),('BACKGROUND',(0,3),(-1,3),LIGHTWARN)]))
story.append(tbl)
story.append(Paragraph('Tabelle 4: Klangcharakter der Hauptvarianten. '
    'Grün: bester Klang (C). Orange: schlechtester (D).',sCap))


# ==========================================================================
section(8,'Ansprache, Lautstärke und Klang: Dasselbe Optimum')
# ==========================================================================

body('Dok. 0005 (Kap. 10) zeigt: Optimale Kammerabstimmung verbessert Ansprache <b>und</b> '
     'Lautstärke gleichzeitig (Win-Win über R<sub>H</sub>). Dieses Dokument ergänzt: '
     '<b>Auch der Klang profitiert vom selben Optimum.</b>')

body('Der Grund ist physikalisch zwingend: Alle drei — Ansprache, Lautstärke, Klang — '
     'leiden unter derselben Quelle: R<sub>H</sub>(f). Er entzieht Energie '
     '(→ schlechtere Ansprache und Lautstärke) <b>und</b> dämpft Obertöne '
     '(→ ärmerer Klang). Wenn R<sub>H</sub> minimiert wird, verbessern sich alle drei.')

keybox('<b>Es gibt keinen Trade-off zwischen Ansprache, Lautstärke und Klang '
       'bei der Kammeroptimierung.</b> Alle drei werden durch dasselbe Optimum bedient: '
       'f<sub>H</sub>/f ≥ 3, minimales Volumen, große Öffnung.<br/><br/>'
       'Der einzige echte Trade-off liegt in der <b>Zungengüte Q<sub>1</sub></b> '
       '(Dok. 0005 Kap. 10): Höheres Q<sub>1</sub> = lauter und obertonreicher, '
       'aber langsamer. Das ist eine Zungen-Entscheidung, keine Kammer-Entscheidung.')


# ==========================================================================
section(9,'Zusammenfassung')
# ==========================================================================

keybox(
    '<b>Die Kammer ist kein Resonanzkörper, sondern ein Kerbfilter.</b> '
    'Sie erzeugt ein Loch im Obertonspektrum bei f<sub>H</sub>.<br/><br/>'
    
    '<b>Welcher Oberton betroffen ist, bestimmt den Klangcharakter:</b><br/>'
    '• 2. OT gedämpft (f<sub>H</sub>/f ≈ 2): <b>dünn, hohl</b> — 25 % Energie fehlt<br/>'
    '• 3. OT gedämpft (f<sub>H</sub>/f ≈ 3): <b>nasal</b> — 11 % fehlt<br/>'
    '• 5.+ OT gedämpft (f<sub>H</sub>/f ≥ 5): <b>kaum verändert</b> — &lt; 4 % Energie<br/><br/>'
    
    '<b>Drei Klangdimensionen:</b> Position (f<sub>H</sub> über Volumen/Öffnung — dominanter Hebel), '
    'Breite (Q<sub>H</sub> über Öffnungsgeometrie — Wandform &lt; 7 %), '
    'Tiefe (κ<sub>eff</sub> über Spaltfläche).<br/><br/>'
    
    '<b>Die Phasenverschiebung</b> (Kramers-Kronig) verändert die Wellenform der Obertöne — '
    'auch wenn die Amplituden im Spektralanalysator kaum verändert erscheinen. '
    'Das erklärt, warum geschulte Ohren Unterschiede hören, die Messgeräte nicht zeigen.<br/><br/>'
    
    '<b>Ansprache, Lautstärke und Klang profitieren vom selben Optimum:</b> '
    'f<sub>H</sub>/f ≥ 3, minimales Volumen, große Öffnung. Kein Trade-off bei der '
    'Kammeroptimierung — nur bei der Zungengüte Q<sub>1</sub>.'
)

story.append(Spacer(1,6*mm))
story.append(HRFlowable(width='60%',thickness=1,color=ACCENTRED,spaceAfter=4,spaceBefore=2))
story.append(Paragraph(
    '<i>Die Kammer erzeugt den Klang nicht — sie filtert ihn. Weniger Filter, mehr Musik.</i>',
    ms('sCl',fs=9,fn=FONTI,al=TA_CENTER,ld=12)))


def hf(canvas,doc):
    canvas.saveState();canvas.setFont(FONT,8);canvas.setFillColor(MIDGRAY)
    canvas.drawString(18*mm,A4[1]-15*mm,'Dok. 0008 — Klangveränderung durch Kammergeometrie')
    canvas.drawRightString(A4[0]-18*mm,A4[1]-15*mm,f'Seite {doc.page}')
    canvas.drawCentredString(A4[0]/2,10*mm,'Querverweise: Dok. 0002 · Dok. 0005 · Dok. 0007')
    canvas.restoreState()

doc.build(story,onFirstPage=hf,onLaterPages=hf)
print(f'PDF: {outpath}')
