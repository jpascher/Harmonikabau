#!/usr/bin/env python3
"""Dok 0007 v6: Vollständige Analyse mit praktischer Tiefenbeschränkung"""

import math, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

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
LIGHTBLUE=HexColor('#e8eaf6'); CRITRED=HexColor('#ffcdd2')
WARNORG=HexColor('#fff3e0'); LIGHTYEL=HexColor('#fff9c4')

try:
    pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuBold','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
    pdfmetrics.registerFont(TTFont('DejaVuItalic','/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
    FONT='DejaVu'; FONTB='DejaVuBold'; FONTI='DejaVuItalic'
except:
    FONT='Helvetica'; FONTB='Helvetica-Bold'; FONTI='Helvetica-Oblique'

outpath='/mnt/user-data/outputs/0007_diskant_kammerfrequenzen_De.pdf'
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
sCell=ms('sCell',fs=7.2,ld=9.5,sa=0,al=TA_LEFT)
sCellB=ms('sCellB',fs=7.2,fn=FONTB,ld=9.5,sa=0,al=TA_LEFT)
sCellC=ms('sCellC',fs=7.2,ld=9.5,sa=0,al=TA_CENTER)
sCellCB=ms('sCellCB',fs=7.2,fn=FONTB,ld=9.5,sa=0,al=TA_CENTER)

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
c=343.0
nn_names=['C','C#','D','D#','E','F','F#','G','G#','A','A#','B']
notes=[]
for midi in range(50,85):
    notes.append((midi,nn_names[midi%12]+str(midi//12-1),440.0*2**((midi-69)/12.0)))
N=len(notes)
W_ch=15.5e-3; L_max=50.0e-3; L_min=30.0e-3

def curve_linear(a,b,x): return a+x*(b-a)
def curve_convex(a,b,x): return a*(b/a)**x          # under line for a>b
def curve_concave(a,b,x): return a+b-a*(b/a)**(1-x)  # above line for a>b

def fH_calc(V,d_open):
    S=math.pi*(d_open/2)**2; l_eff=8.0e-3+0.85*d_open
    return (c/(2*math.pi))*math.sqrt(S/(V*l_eff))

def evaluate(freq,fH,Q_H=7.0):
    best_Rw=0;best_n=1;best_Gw=0;total=0
    for n in range(1,21):
        fot=n*freq;a=1.0/n;ew=a**2
        r=fot/fH;D=(1-r**2)**2+(r/Q_H)**2
        Rn=(r/Q_H)/D;G=1.0/math.sqrt(D)
        Rw=ew*Rn;Gw=a*G;total+=Rw
        if Rw>best_Rw: best_Rw=Rw;best_n=n;best_Gw=Gw
    return best_n,best_Gw,total

def calc_variant(depth_fn,diam_fn):
    res=[]
    for i,(midi,name,freq) in enumerate(notes):
        frac=i/(N-1); L=L_max-frac*(L_max-L_min)
        d=depth_fn(frac);diam=diam_fn(frac)
        V=W_ch*L*d;f_H=fH_calc(V,diam)
        n_dom,Gw,Rtot=evaluate(freq,f_H)
        res.append(dict(name=name,freq=freq,L=L*1e3,depth=d*1e3,
            diam=diam*1e3,V=V*1e6,f_H=f_H,ratio=f_H/freq,
            n_dom=n_dom,G_weighted=Gw,R_total=Rtot))
    return res

def cnt(res):
    cr=sum(1 for r in res if r['n_dom']<=2 and r['G_weighted']>0.5)
    mk=sum(1 for r in res if r['n_dom']<=3 and r['G_weighted']>0.3 and not(r['n_dom']<=2 and r['G_weighted']>0.5))
    return cr,mk,N-cr-mk

dep8=8.0e-3;dep3=3.0e-3;d10=10.0e-3;d6=6.0e-3

# Practical minimum depth: tongue must not hit floor
# Diskant tongue max tip deflection estimates:
# Long tongues (D3 ~35mm tongue): amplitude ~3-4mm below plate → need ~5mm clearance
# Short tongues (C6 ~8mm tongue): amplitude ~1.5-2mm → need ~3mm clearance
# Add 0.5mm safety → d_min follows tongue length
# Approximation: d_min ~ 2.5 + 0.04*L_tongue, or simply: d_min(D3)~5mm, d_min(C6)~3mm

# Build all variants
V={}
V['A']=('Referenz','Tiefe 8 mm, ∅ 10 mm',lambda f:dep8,lambda f:d10,'#78909c','x')
V['B']=('Tiefe lin.','Tiefe lin. 8→3, ∅ 10',lambda f:curve_linear(dep8,dep3,f),lambda f:d10,'#c62828','s')
V['C']=('Tiefe konvex','Tiefe konvex 8→3, ∅ 10',lambda f:curve_convex(dep8,dep3,f),lambda f:d10,'#2e7d32','^')
V['Ci']=('Tiefe konkav','Tiefe konkav 8→3, ∅ 10',lambda f:curve_concave(dep8,dep3,f),lambda f:d10,'#81c784','v')
V['D']=('∅ lin.','Tiefe 8, ∅ lin. 10→6',lambda f:dep8,lambda f:curve_linear(d10,d6,f),'#1565c0','o')
V['E']=('∅ konvex','Tiefe 8, ∅ konvex 10→6',lambda f:dep8,lambda f:curve_convex(d10,d6,f),'#e65100','D')
V['Ei']=('∅ konkav','Tiefe 8, ∅ konkav 10→6',lambda f:dep8,lambda f:curve_concave(d10,d6,f),'#ff9800','d')

keys_all=['A','B','C','Ci','D','E','Ei']
keys_depth=['A','B','C','Ci']; keys_open=['A','D','E','Ei']

results={};counts_v={};scores={}
for k in keys_all:
    results[k]=calc_variant(V[k][2],V[k][3])
    counts_v[k]=cnt(results[k])
    scores[k]=sum(r['R_total'] for r in results[k])

# Practical depth limits
tongue_lengths = {}
for i,(midi,name,freq) in enumerate(notes):
    # Rough estimate: diskant tongue length ~ 35mm (D3) to 8mm (C6)
    frac = i/(N-1)
    L_tongue = 35.0 - frac * 27.0  # mm
    max_deflection = 0.10 * L_tongue  # ~10% of length as max amplitude
    min_depth = max_deflection + 0.5  # +0.5mm safety
    tongue_lengths[name] = (L_tongue, max_deflection, min_depth)


# ======================= CHARTS =======================
fig,axes=plt.subplots(2,2,figsize=(8,7),dpi=150)
fig.suptitle('Dok. 0007: Sieben Varianten — Tiefe vs. Öffnung',
             fontsize=12,fontweight='bold',color='#16213e')
freqs=[r['freq'] for r in results['A']]
nms=[r['name'] for r in results['A']]

# 1: Depth + practical limit
ax=axes[0,0]
for k in keys_depth:
    ax.plot(freqs,[r['depth'] for r in results[k]],
            f'{V[k][5]}-',color=V[k][4],ms=2.5,lw=1.3,label=f'{k}: {V[k][0]}',alpha=0.85)
# practical minimum
min_depths_plot = [tongue_lengths[results['A'][i]['name']][2] for i in range(N)]
ax.fill_between(freqs, 0, min_depths_plot, alpha=0.15, color='red', label='Anschlagzone')
ax.plot(freqs, min_depths_plot, 'r--', lw=1, alpha=0.6, label='Min. Tiefe (Zungenamplitude)')
ax.set_ylabel('Kammertiefe [mm]',fontsize=9)
ax.set_title('Tiefenverläufe + praktisches Minimum',fontsize=10,color='#16213e')
ax.legend(fontsize=6,loc='upper right'); ax.grid(True,alpha=0.3)
ax.set_xlim(freqs[0]-10,freqs[-1]+10); ax.set_ylim(0,9)
ax.set_xticks(freqs[::6]); ax.set_xticklabels(nms[::6],fontsize=6,rotation=45)

# 2: Opening curves
ax=axes[0,1]
for k in keys_open:
    ax.plot(freqs,[r['diam'] for r in results[k]],
            f'{V[k][5]}-',color=V[k][4],ms=2.5,lw=1.3,label=f'{k}: {V[k][0]}',alpha=0.85)
ax.set_ylabel('∅ Öffnung [mm]',fontsize=9)
ax.set_title('Öffnungsverläufe',fontsize=10,color='#16213e')
ax.legend(fontsize=6.5); ax.grid(True,alpha=0.3)
ax.set_xlim(freqs[0]-10,freqs[-1]+10); ax.set_ylim(5,11)
ax.set_xticks(freqs[::6]); ax.set_xticklabels(nms[::6],fontsize=6,rotation=45)

# 3: f_H/f
ax=axes[1,0]
for k in keys_all:
    ax.plot(freqs,[r['ratio'] for r in results[k]],
            f'{V[k][5]}-',color=V[k][4],ms=2,lw=1,label=f'{k}',alpha=0.8)
ax.axhspan(0,2,alpha=0.12,color='red'); ax.axhspan(2,3,alpha=0.07,color='orange')
ax.axhline(y=3,color='red',ls='--',alpha=0.4,lw=0.8)
ax.set_ylabel('f_H / f',fontsize=9); ax.set_ylim(0,14)
ax.set_title('f_H/f — alle Varianten',fontsize=10,color='#16213e')
ax.legend(fontsize=5.5,ncol=2,loc='upper right'); ax.grid(True,alpha=0.3)
ax.set_xlim(freqs[0]-10,freqs[-1]+10)
ax.set_xticks(freqs[::6]); ax.set_xticklabels(nms[::6],fontsize=6,rotation=45)

# 4: Summary bars
ax=axes[1,1]
ranked=sorted(keys_all,key=lambda k:scores[k])
x_pos=range(len(ranked))
vals=[scores[k] for k in ranked]
colors=[V[k][4] for k in ranked]
bars=ax.bar(x_pos,vals,color=colors,alpha=0.7,edgecolor='black',linewidth=0.5)
ax.set_xticks(x_pos)
ax.set_xticklabels([f'{k}\n{V[k][0]}' for k in ranked],fontsize=6)
ax.set_ylabel('Σ R_eff',fontsize=9)
ax.set_title('Rangliste (niedrig = besser)',fontsize=10,color='#16213e')
ax.grid(True,alpha=0.3,axis='y')
for i,v in enumerate(vals):
    c=counts_v[ranked[i]]
    ax.text(i,v+0.2,f'{c[0]}/{c[1]}/{c[2]}',ha='center',fontsize=6)
ax.text(0.5,0.97,'krit/merkl/gut',transform=ax.transAxes,ha='center',fontsize=6,color='gray')

plt.tight_layout(rect=[0,0,1,0.95])
chart_path='/home/claude/chart_0007_v6.png'
plt.savefig(chart_path,dpi=150,bbox_inches='tight'); plt.close()


# ======================= PDF =======================
story.append(Paragraph('Dok. 0007: Diskant-Stimmstock — Kammerfrequenzen',sTitle))
story.append(Paragraph('Sieben Varianten mit praktischer Tiefenbeschränkung',sSubtitle))
story.append(HRFlowable(width='80%',thickness=2,color=ACCENTRED,spaceAfter=3,spaceBefore=3))
story.append(Paragraph(
    'Vollständige Analyse für 35 Kammern (D3–C6). Sieben Varianten: Referenz + je drei '
    'Kurvenformen (linear, konvex, konkav) für Tiefe und Öffnung. '
    'Berücksichtigt die praktische Beschränkung: Die Kammertiefe darf nicht unter die '
    'maximale Zungenamplitude fallen. Amplitudengewichtete Bewertung (a_n ~ 1/n).', sNote))

story.append(Spacer(1,2*mm))
ref_data=[
    ['Dok. 0002','Strömungsanalyse Bass-Stimmzunge 50 Hz — v8'],
    ['Dok. 0004','Frequenzvariation — Zwei-Filter-Modell'],
    ['Dok. 0005','Frequenzverschiebung als Indikator der Ansprache'],
    ['Dok. 0006','Zeichenerklärung'],
    ['<b>Dok. 0007</b>','<b>Diskant-Stimmstock — Kammerfrequenzen (dieses Dokument)</b>'],
]
ref_tbl=Table([[Paragraph(r[0],sCell),Paragraph(r[1],sCell)] for r in ref_data],
              colWidths=[25*mm,doc.width-25*mm-4*mm])
ref_tbl.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.5,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('LEFTPADDING',(0,0),(-1,-1),2*mm),('BACKGROUND',(0,4),(-1,4),LIGHTGREEN)]))
story.append(ref_tbl)
story.append(Spacer(1,4*mm))


# ==========================================================================
section(1,'Geometrie und Ausgangslage')
# ==========================================================================

body('Der Diskant-Stimmstock enthält 35 Kammern für die Töne D3 (146,8 Hz) bis C6 (1046,5 Hz). '
     'Jede Kammer ist ein Quader mit konstanter Breite W = 15,5 mm und einer Länge, die linear '
     'von 50 mm (D3) auf 30 mm (C6) abnimmt. Die Auslassöffnung ist kreisrund mit ∅ 10 mm '
     '(S<sub>Hals</sub> = 78,5 mm²) und einer physischen Halshöhe von 8 mm '
     '(effektive Halslänge l<sub>eff</sub> = 16,5 mm mit Mündungskorrektur 0,85·d).')

body('Die Helmholtz-Frequenz jeder Kammer hängt von Volumen und Öffnung ab:')
formula('f<sub>H</sub> = (c/2π) · √(S<sub>Hals</sub> / (V · l<sub>eff</sub>))')

body('Die zentrale Frage: Wie sollen Kammertiefe und Öffnungsdurchmesser über den Tonbereich '
     'verlaufen, um die Kammer-Oberton-Resonanzen (Dok. 0005) zu minimieren?')

body('Aus Dok. 0005 wissen wir: Ansprache-Probleme treten auf, wenn der 2. oder 3. Oberton '
     'einer Zunge in die Nähe der Kammerresonanz f<sub>H</sub> kommt. Das entspricht '
     'f<sub>H</sub>/f &lt; 3. Oberhalb von f<sub>H</sub>/f = 4 ist die Kammer akustisch '
     '„unsichtbar" — die Kopplung durch höhere Obertöne (n ≥ 5) ist energetisch irrelevant '
     '(a<sub>n</sub>² &lt; 4 %).')


# ==========================================================================
section(2,'Die praktische Beschränkung: Zungenanschlag')
# ==========================================================================

body('Die Berechnung zeigt, dass kleinere Kammervolumina f<sub>H</sub> nach oben schieben und '
     'damit die Ansprache verbessern. Aber die Kammertiefe kann nicht beliebig reduziert werden: '
     '<b>Die Zunge muss Platz zum Ausschwingen haben.</b>')

body('Eine Durchschlagzunge schwingt bei voller Lautstärke durch den Schlitz der Stimmplatte '
     'nach unten in die Kammer hinein. Wenn die Zungenspitze den Kammerboden berührt, '
     'wird die Schwingung abrupt gestoppt — die Zunge „schlägt an". Das erzeugt ein '
     'hässliches Klappgeräusch, beschädigt auf Dauer die Zunge, und begrenzt die Lautstärke.')

subsection('Maximale Zungenamplitude im Diskant')

body('Die maximale Auslenkung der Zungenspitze hängt von der Zungenlänge und der '
     'Spieldynamik ab. Typische Werte für Diskantzungen:')

amp_data=[['Tonbereich','Zungenlänge','Max. Amplitude','Min. Kammertiefe']]
for name,Ltong,amp,dmin in [
    ('D3 (tiefster)',35,3.5,4.0),('A3',28,2.8,3.3),('A4 (Mitte)',18,1.8,2.3),
    ('A5',10,1.0,1.5),('C6 (höchster)',8,0.8,1.3)]:
    amp_data.append([name,f'{Ltong:.0f} mm',f'~{amp:.1f} mm',f'{dmin:.1f} mm'])

dw=[[Paragraph(c,sCellCB if i==0 else sCellC) for c in row] for i,row in enumerate(amp_data)]
tbl=Table(dw,colWidths=[28*mm,25*mm,28*mm,30*mm],repeatRows=1)
tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHTBLUE),('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
story.append(tbl)
story.append(Paragraph('Tabelle 1: Geschätzte Zungenamplituden und Mindesttiefen. '
    'Amplitude ≈ 10 % der Zungenlänge + 0,5 mm Sicherheit.',sCap))

body('Die minimale Kammertiefe muss also von etwa 4 mm (D3) auf 1,3 mm (C6) sinken. '
     'Die hier berechneten Varianten mit einer Mindesttiefe von 3 mm (C6) liegen '
     '<b>deutlich über</b> dem Zungenanschlag — es gibt also Spielraum nach unten.')

keybox('<b>Praktische Regel:</b> Die Kammertiefe sollte so gering wie möglich sein '
       '(minimales Volumen → höchstes f<sub>H</sub> → beste Entkopplung), aber <b>niemals</b> '
       'unter die maximale Zungenamplitude plus Sicherheitsabstand fallen. '
       'Bei einer Mindesttiefe von 3 mm (wie in unseren Varianten B, C, Ci) ist der Anschlag '
       'für alle Töne oberhalb von D3 sicher ausgeschlossen. Für die tiefsten Diskanttöne '
       '(D3–F3) könnte die Tiefe mit 3 mm knapp werden — hier wären 4 mm sicherer.')


# ==========================================================================
section(3,'Die drei Kurvenformen')
# ==========================================================================

body('Jeder geometrische Parameter (Tiefe oder Öffnungsdurchmesser) kann von seinem '
     'Startwert a auf seinen Endwert b auf drei verschiedene Arten verlaufen:')

formula('Linear: p(x) = a + x · (b − a)')
formula('Konvex: p(x) = a · (b/a)<super>x</super>')
formula('Konkav: p(x) = a + b − a · (b/a)<super>(1−x)</super>')

body('Alle drei Kurven verbinden <b>dieselben Endpunkte</b>: p(0) = a und p(1) = b. '
     'Sie unterscheiden sich in der <b>Verteilung</b> der Reduktion über den Tonbereich:')

body('<b>Linear</b> verteilt die Reduktion gleichmäßig: In jedem Halbtonschritt ändert sich '
     'der Parameter um denselben absoluten Betrag (z.B. −0,147 mm/Halbton bei 8→3 mm).')

body('<b>Konvex</b> (liegt <b>unter</b> der Geraden) reduziert am Anfang schneller als linear '
     'und am Ende langsamer. Im Mittelbereich hat der Parameter einen <b>niedrigeren</b> Wert '
     'als bei linearem Verlauf. Beispiel Tiefe: Bei A4 hat konvex 4,9 mm, linear 5,5 mm — '
     'also 0,6 mm weniger.')

body('<b>Konkav</b> (liegt <b>über</b> der Geraden) reduziert am Anfang langsamer und am Ende '
     'schneller. Im Mittelbereich hat der Parameter einen <b>höheren</b> Wert als bei linearem '
     'Verlauf. Beispiel Tiefe: Bei A4 hat konkav ca. 6,1 mm, linear 5,5 mm — also 0,6 mm mehr.')

warnbox('<b>Welche Kurve besser ist, hängt davon ab, in welche Richtung der Parameter '
        'f<sub>H</sub> verschiebt:</b><br/><br/>'
        '• Wenn <b>weniger</b> des Parameters f<sub>H</sub> <b>hebt</b> (wie bei der Tiefe: '
        'weniger V → höheres f<sub>H</sub>), ist die Kurve besser, die im Mittelbereich '
        'weniger hat → <b>konvex gewinnt</b>.<br/>'
        '• Wenn <b>mehr</b> des Parameters f<sub>H</sub> <b>hebt</b> (wie bei der Öffnung: '
        'mehr S → höheres f<sub>H</sub>), ist die Kurve besser, die im Mittelbereich '
        'mehr hat → <b>konkav gewinnt</b>.')


# ==========================================================================
section(4,'Die sieben Varianten')
# ==========================================================================

body('Wir testen alle Kombinationen: Referenz (A) + drei Tiefenvarianten bei ∅ 10 mm (B, C, Ci) '
     '+ drei Öffnungsvarianten bei Tiefe 8 mm (D, E, Ei):')

vt=[['Var.','Kurzname','Tiefe [mm]','∅ Öffnung [mm]','Wirkung auf f_H']]
descs={
    'A':('8 konstant','10 konstant','Referenz'),
    'B':('8→3 linear','10 konstant','f_H steigt (V sinkt)'),
    'C':('8→3 konvex','10 konstant','f_H steigt stärker (V sinkt schneller im Mittel)'),
    'Ci':('8→3 konkav','10 konstant','f_H steigt schwächer (V bleibt höher im Mittel)'),
    'D':('8 konstant','10→6 linear','f_H sinkt (S sinkt)'),
    'E':('8 konstant','10→6 konvex','f_H sinkt stärker (S sinkt schneller im Mittel)'),
    'Ei':('8 konstant','10→6 konkav','f_H sinkt schwächer (S bleibt höher im Mittel)'),
}
for k in keys_all:
    vt.append([k,V[k][0],descs[k][0],descs[k][1],descs[k][2]])
dw=[[Paragraph(c,sCellCB if i==0 else sCellC) for c in row] for i,row in enumerate(vt)]
tbl=Table(dw,colWidths=[11*mm,22*mm,25*mm,25*mm,48*mm],repeatRows=1)
tbl.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),LIGHTBLUE),('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('BACKGROUND',(0,1),(-1,1),LIGHTGRAY),
    ('BACKGROUND',(0,5),(-1,7),HexColor('#e3f2fd'))]))
story.append(tbl)
story.append(Paragraph('Tabelle 2: Die sieben Varianten. Obere Gruppe (B–Ci): Tiefe variiert. '
    'Untere Gruppe (D–Ei): Öffnung variiert.',sCap))

subsection('Warum die Richtung entscheidet')

body('<b>Tiefe senken</b> verkleinert das Kammervolumen V. Da f<sub>H</sub> ∝ 1/√V, '
     '<b>steigt</b> f<sub>H</sub>. Das verschiebt f<sub>H</sub>/f nach oben — weg von den '
     'niedrigen Obertönen. Das ist die <b>richtige Richtung</b>.')

body('<b>Öffnung verkleinern</b> reduziert die Halsfläche S<sub>Hals</sub> und gleichzeitig '
     'die Mündungskorrektur (l<sub>eff</sub> = 8 + 0,85·d sinkt ebenfalls). '
     'Der Nettoeffekt: f<sub>H</sub> ∝ √(S/l<sub>eff</sub>) <b>sinkt</b>. '
     'Quantitativ: ∅ von 10 auf 6 mm senkt f<sub>H</sub> auf 67 %. '
     'Das verschiebt f<sub>H</sub>/f nach <b>unten</b> — das ist die <b>falsche Richtung</b>.')


# ==========================================================================
section(5,'Diagramme')
# ==========================================================================

if os.path.exists(chart_path):
    story.append(Image(chart_path,width=170*mm,height=155*mm))
    story.append(Paragraph(
        'Abbildung 1: Oben links: Tiefenverläufe mit Anschlagzone (rot schraffiert). '
        'Oben rechts: Öffnungsverläufe. '
        'Unten links: f_H/f — rote Zone: 2. OT kann koppeln. '
        'Unten rechts: Rangliste nach Gesamt-Kopplungslast (krit/merkl/gut über Balken).',sCap))


# ==========================================================================
story.append(PageBreak())
section(6,'Hauptergebnis')
# ==========================================================================

comp=[['','A','B','C','Ci','D','E','Ei']]
for label,fn in [
    ('f_H(D3) [Hz]',lambda r:f'{r[0]["f_H"]:.0f}'),
    ('f_H(C6) [Hz]',lambda r:f'{r[-1]["f_H"]:.0f}'),
    ('f_H/f min.',lambda r:f'{min(x["ratio"] for x in r):.1f}'),
    ('Kritisch (OT≤2)',lambda r:f'{cnt(r)[0]}'),
    ('Merklich (OT=3)',lambda r:f'{cnt(r)[1]}'),
    ('Gut (OT≥4)',lambda r:f'{cnt(r)[2]}'),
    ('Σ R_eff',lambda r:f'{sum(x["R_total"] for x in r):.1f}'),
]:
    row=[label]
    for k in keys_all: row.append(fn(results[k]))
    comp.append(row)

dw=[[Paragraph(c,sCellB if i==0 else sCell) for c in row] for i,row in enumerate(comp)]
cw=[26*mm]+[17*mm]*7
tbl=Table(dw,colWidths=cw,repeatRows=1)
sc=[('BACKGROUND',(0,0),(-1,0),DARKBLUE),('TEXTCOLOR',(0,0),(-1,0),white),
    ('GRID',(0,0),(-1,-1),0.4,MIDGRAY),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm),
    ('LEFTPADDING',(0,0),(-1,-1),1.5*mm),('VALIGN',(0,0),(-1,-1),'TOP')]
for i in range(2,len(comp),2): sc.append(('BACKGROUND',(0,i),(-1,i),LIGHTGRAY))
best=min(scores,key=scores.get)
best_col=keys_all.index(best)+1
sc.append(('BACKGROUND',(best_col,1),(best_col,-1),LIGHTGREEN))
tbl.setStyle(TableStyle(sc))
story.append(tbl)
story.append(Paragraph(f'Tabelle 3: Hauptvergleich. Grün: Variante {best} (niedrigstes Σ R_eff).',sCap))

subsection('Rangliste')

for i,k in enumerate(ranked):
    c=counts_v[k]
    marker='⬤' if c[0]>5 else '▲' if c[0]>0 else '○'
    story.append(Paragraph(
        f'<b>{i+1}. {k}</b> ({V[k][0]}): {c[0]} krit., {c[1]} merkl., {c[2]} gut — '
        f'Σ R = {scores[k]:.1f} {marker}', sBody))


# ==========================================================================
section(7,'Analyse der Tiefengruppe')
# ==========================================================================

body('Alle drei Tiefenvarianten (B, C, Ci) eliminieren die kritischen Resonanzen vollständig '
     '(0 kritische Töne). Sie unterscheiden sich nur in der Zahl der merklichen Resonanzen:')

body(f'<b>C (konvex, unter der Geraden):</b> {counts_v["C"][1]} merkliche. Die Kurve gibt den '
     f'mittleren Tönen weniger Volumen → höheres f<sub>H</sub> → bessere Entkopplung. '
     f'Beispiel A4: Tiefe 4,9 mm, V = 2,78 cm³, f<sub>H</sub> = 2258 Hz, f<sub>H</sub>/f = 5,1. '
     f'<b>Beste Tiefenvariante.</b>')

body(f'<b>B (linear):</b> {counts_v["B"][1]} merkliche. Gleichmäßige Reduktion. '
     f'A4: Tiefe 5,5 mm, V = 3,13 cm³, f<sub>H</sub> = 2128 Hz, f<sub>H</sub>/f = 4,8.')

body(f'<b>Ci (konkav, über der Geraden):</b> {counts_v["Ci"][1]} merkliche. Die Kurve erhält '
     f'im Mittelbereich mehr Volumen → niedrigeres f<sub>H</sub>. '
     f'A4: Tiefe ~6,1 mm, f<sub>H</sub>/f ~4,5. '
     f'Schlechteste Tiefenvariante — aber immer noch 0 kritische.')

body('Der Unterschied zwischen den drei Tiefenvarianten ist <b>moderat</b> '
     '(12,3 vs. 12,8 vs. 13,5 in Σ R<sub>eff</sub>). Der <b>große Sprung</b> ist von der '
     'Referenz A (17,4) zu jeder der drei Tiefenvarianten — die Tiefenreduktion an sich '
     'ist viel wichtiger als die Kurvenform.')

keybox('<b>Praktische Bedeutung:</b> Der Unterschied zwischen konvex und linear entspricht '
       '2 merklichen Resonanzen weniger — hörbar, aber nicht dramatisch. Der Unterschied '
       'zwischen Tiefenreduktion und keiner (A vs. B/C/Ci) entspricht 7 kritischen Resonanzen — '
       'das ist der entscheidende Schritt. <b>Wer die Tiefe reduziert, hat 95 % des Effekts. '
       'Die Kurvenform ist Feinabstimmung.</b>')


# ==========================================================================
section(8,'Analyse der Öffnungsgruppe')
# ==========================================================================

body('Die drei Öffnungsvarianten (D, E, Ei) verschlechtern die Situation gegenüber der '
     'Referenz — sie haben <b>mehr</b> kritische Töne (12–13 statt 7):')

S10=math.pi*5e-3**2; l10=8e-3+0.85*10e-3
S6=math.pi*3e-3**2; l6=8e-3+0.85*6e-3
ratio_fH=math.sqrt(S6/l6)/math.sqrt(S10/l10)

body(f'Der Grund: ∅ von 10 auf 6 mm senkt S<sub>Hals</sub> von {S10*1e6:.1f} auf {S6*1e6:.1f} mm² '
     f'(−64 %) und l<sub>eff</sub> von {l10*1e3:.1f} auf {l6*1e3:.1f} mm (−32 %). '
     f'Der Nettoeffekt: f<sub>H</sub> sinkt auf <b>{ratio_fH*100:.0f} %</b>. '
     f'f<sub>H</sub>/f fällt auf 1,3 bei C6 — dort liegt f<sub>H</sub> <b>unter</b> der '
     f'Zungenfrequenz selbst.')

body(f'<b>Ei (konkav):</b> {counts_v["Ei"][0]} krit. — am wenigsten schädlich, weil die Kurve '
     f'die größeren Öffnungen im Mittelbereich erhält → f<sub>H</sub> bleibt dort höher.')

body(f'<b>E (konvex):</b> {counts_v["E"][0]} krit. — am schädlichsten, weil die Kurve die '
     f'Öffnung im Mittelbereich am stärksten reduziert → f<sub>H</sub> sinkt dort am meisten.')

warnbox('<b>Fazit Öffnungsgruppe:</b> Keine der drei Öffnungsvarianten ist besser als die '
        'Referenz A. Die Öffnungsverkleinerung 10→6 mm allein ist <b>kontraproduktiv</b>. '
        'Auch die beste Öffnungsvariante (Ei) hat 12 kritische Töne — fast doppelt so viele '
        'wie die Referenz (7). Die Öffnung allein zu variieren ist der falsche Ansatz.')


# ==========================================================================
section(9,'Volumen minimal halten — die praktische Optimierungsregel')
# ==========================================================================

body('Aus der Analyse folgt eine klare Optimierungsregel, die auch der praktischen '
     'Erfahrung des Instrumentenbauers entspricht:')

keybox('<b>Optimierungsregel: Volumen so gering wie die Zunge erlaubt.</b><br/><br/>'
       'Die Kammertiefe wird für jede Kammer so gewählt, dass:<br/>'
       '(1) Die Zunge bei maximaler Lautstärke <b>nicht anschlägt</b> '
       '(Tiefe ≥ max. Amplitude + 0,5 mm Sicherheit).<br/>'
       '(2) Das Volumen <b>nicht größer</b> ist als für Bedingung (1) nötig.<br/>'
       '(3) Die Öffnung <b>so groß wie möglich</b> bleibt (10 mm), weil jede Verkleinerung '
       'f<sub>H</sub> in die falsche Richtung schiebt.<br/><br/>'
       'In der Praxis bedeutet das: Die Kammertiefe folgt der Zungenamplitude — kurze Zungen '
       '(hohe Töne) brauchen weniger Platz, bekommen flachere Kammern, und profitieren '
       'automatisch vom höheren f<sub>H</sub>/f. Der konvexe Tiefenverlauf (C) ist dabei '
       'etwas besser als der lineare (B), weil er im Mittelbereich aggressiver reduziert.')

body('Das Ergebnis stimmt mit der Instrumentenbau-Praxis überein: Erfahrene Baumeister '
     'fertigen die Kammern so, dass die Zunge gerade eben nicht anschlägt, und machen die '
     'Kammer nicht tiefer als nötig. Die Berechnung zeigt jetzt <b>warum</b> das funktioniert: '
     'Weniger Volumen → höheres f<sub>H</sub> → weniger Oberton-Resonanzen → bessere Ansprache.')


# ==========================================================================
section(10,'Zusammenfassung')
# ==========================================================================

keybox(
    f'<b>Sieben Varianten für den Diskant-Stimmstock D3–C6:</b><br/><br/>'
    f'<b>Rangliste:</b><br/>'
    + '<br/>'.join(f'{i+1}. <b>{k}</b> ({V[k][0]}): {counts_v[k][0]} krit., '
                    f'{counts_v[k][1]} merkl., {counts_v[k][2]} gut'
                    for i,k in enumerate(ranked))
    + '<br/><br/>'
    '<b>Kernaussagen:</b><br/>'
    '• <b>Tiefe reduzieren</b> (B, C, Ci) hebt f<sub>H</sub> — richtige Richtung, 0 kritische Töne<br/>'
    '• <b>Öffnung verkleinern</b> (D, E, Ei) senkt f<sub>H</sub> — falsche Richtung, 12–13 kritische<br/>'
    '• Innerhalb Tiefe: <b>konvex &gt; linear &gt; konkav</b> (12,3 vs. 12,8 vs. 13,5)<br/>'
    '• Innerhalb Öffnung: <b>konkav &gt; linear &gt; konvex</b> (17,0 vs. 17,2 vs. 17,3)<br/>'
    '• Der große Sprung liegt in der <b>Tiefenreduktion an sich</b>, nicht in der Kurvenform<br/><br/>'
    '<b>Praktische Regel:</b> Kammertiefe = maximale Zungenamplitude + 0,5 mm Sicherheit. '
    'Nicht tiefer. Öffnung so groß wie möglich (10 mm). '
    'Konvexer Tiefenverlauf ist optimal, linearer fast ebenso gut.<br/><br/>'
    '<b>f<sub>H</sub>/f ≥ 3 halten → kein 2. OT-Treffer → keine kritischen Probleme.</b>'
)

story.append(Spacer(1,6*mm))
story.append(HRFlowable(width='60%',thickness=1,color=ACCENTRED,spaceAfter=4,spaceBefore=2))
story.append(Paragraph(
    '<i>Minimales Volumen, maximale Öffnung — die Kammer so klein wie die Zunge erlaubt.</i>',
    ms('sCl',fs=9,fn=FONTI,al=TA_CENTER,ld=12)))


def hf(canvas,doc):
    canvas.saveState();canvas.setFont(FONT,8);canvas.setFillColor(MIDGRAY)
    canvas.drawString(18*mm,A4[1]-15*mm,'Dok. 0007 — Diskant-Stimmstock: Kammerfrequenzen D3–C6')
    canvas.drawRightString(A4[0]-18*mm,A4[1]-15*mm,f'Seite {doc.page}')
    canvas.drawCentredString(A4[0]/2,10*mm,'7 Varianten: A (Ref.) / B,C,Ci (Tiefe) / D,E,Ei (Öffnung)')
    canvas.restoreState()

doc.build(story,onFirstPage=hf,onLaterPages=hf)
print(f'PDF: {outpath}')
print(f'\n{"Rang":>4} {"Var":>3} {"Label":<28} {"Krit":>5} {"Merk":>5} {"Gut":>5} {"ΣR":>7} {"fH/f min":>9}')
for i,k in enumerate(ranked):
    c=counts_v[k]
    print(f'{i+1:>4}   {k:>3}  {V[k][1]:<28} {c[0]:>5} {c[1]:>5} {c[2]:>5} {scores[k]:>7.1f} {min(r["ratio"] for r in results[k]):>9.1f}')
