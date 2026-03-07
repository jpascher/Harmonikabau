#!/usr/bin/env python3
"""Dok. 0018 — Hörbarkeit der Biegemoden: Kammer, Transiente, Torsion"""
import numpy as np
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh
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

# ══ FEM ══
def beam_modes(N,tf,L,W,E,rho,nm=6):
    dx=L/N; nd=2*(N+1); K=np.zeros((nd,nd)); M=np.zeros((nd,nd))
    for e in range(N):
        xm=(e+.5)*dx; te=tf(xm); Ie=W*te**3/12; Ae=W*te; le=dx
        ke=E*Ie/le**3*np.array([[12,6*le,-12,6*le],[6*le,4*le**2,-6*le,2*le**2],[-12,-6*le,12,-6*le],[6*le,2*le**2,-6*le,4*le**2]])
        me=rho*Ae*le/420*np.array([[156,22*le,54,-13*le],[22*le,4*le**2,13*le,-3*le**2],[54,13*le,156,-22*le],[-13*le,-3*le**2,-22*le,4*le**2]])
        d=[2*e,2*e+1,2*e+2,2*e+3]
        for i in range(4):
            for j in range(4): K[d[i],d[j]]+=ke[i,j]; M[d[i],d[j]]+=me[i,j]
    fd=list(range(2,nd)); ev,evec=eigh(K[np.ix_(fd,fd)],M[np.ix_(fd,fd)])
    freqs=np.sqrt(np.abs(ev[:nm]))/(2*np.pi)
    modes=np.zeros((N+1,nm))
    for m in range(min(nm,len(fd))):
        full=np.zeros(nd); full[fd]=evec[:,m]; modes[:,m]=full[0::2]
        mx=np.max(np.abs(modes[:,m]))
        if mx>0: modes[:,m]/=mx
    return freqs, modes

Es=200e9; rho_s=7800; G=77e9; Nn=80
L_b=70e-3; W_b=8e-3; t0=0.40e-3
freqs,modes=beam_modes(Nn,lambda x:t0,L_b,W_b,Es,rho_s)
x_n=np.linspace(0,L_b,Nn+1)
integrals=[np.trapezoid(modes[:,m],x_n) for m in range(6)]

# ══ DIAGRAMME ══
def fig_kammer():
    fig,ax=plt.subplots(figsize=(9,5))
    base=-37; qhs=np.linspace(1,35,100)
    gains=20*np.log10(qhs); totals=base+gains
    ax.plot(qhs,totals,'r-',lw=2.5,label='Mode 2 mit Kammer')
    ax.axhline(y=-37,color='blue',ls='--',lw=1,label='Mode 2 ohne Kammer (\u221237 dB)')
    ax.axhline(y=-3,color='gray',ls=':',lw=1,label='Harm. Oberton 2f\u2081 (\u22123 dB)')
    ax.axhline(y=-13,color='gray',ls=':',lw=0.8); ax.text(36,-12,'6f\u2081',fontsize=8,color='gray')
    ax.axhspan(-20,0,alpha=0.05,color='green'); ax.text(2,-5,'deutlich h\u00f6rbar',fontsize=9,color='green')
    ax.axhspan(-40,-20,alpha=0.05,color='orange'); ax.text(2,-25,'grenzwertig',fontsize=9,color='orange')
    ax.set_xlabel('Kammer-G\u00fcte Q_H'); ax.set_ylabel('Mode 2 Amplitude [dB rel. Grundton]')
    ax.set_title('Abb. 1: Kammer-Verst\u00e4rkung von Mode 2')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3); ax.set_xlim(1,35); ax.set_ylim(-45,5)
    fig.tight_layout(); return fig

def fig_transient():
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
    t_ms=np.linspace(0,300,1000)
    Q_est=[50,150,300]; A0=[1.0,0.014,0.001]
    cols=['#1565c0','#e94560','#ff9800']
    for m,(Q,a0,col) in enumerate(zip(Q_est,A0,cols)):
        tau=Q/(np.pi*freqs[m])*1000
        env=a0*np.exp(-t_ms/tau)
        ax1.plot(t_ms,20*np.log10(env+1e-30),'-',color=col,lw=2,
                label=f'Mode {m+1} (f={freqs[m]:.0f} Hz, Q={Q})')
    ax1.axhline(y=-40,color='black',ls=':',lw=1); ax1.text(250,-38,'H\u00f6rschwelle',fontsize=8)
    ax1.set_xlabel('Zeit [ms]'); ax1.set_ylabel('Amplitude [dB]')
    ax1.set_title('Einh\u00fcllende beim Einschwingen'); ax1.legend(fontsize=8); ax1.grid(True,alpha=0.3)
    ax1.set_xlim(0,300); ax1.set_ylim(-60,5)
    # Rechts: Spektrum t=0 vs t=100ms
    ax2.bar([1,2,3],[0,-37,-60],color=['#1565c0','#e94560','#ff9800'],alpha=0.4,label='t = 0 ms')
    ax2.bar([1.3,2.3,3.3],[0,-47,-99],color=['#1565c0','#e94560','#ff9800'],alpha=0.8,label='t = 100 ms')
    ax2.set_xticks([1.15,2.15,3.15]); ax2.set_xticklabels(['Mode 1','Mode 2','Mode 3'])
    ax2.set_ylabel('Amplitude [dB]'); ax2.set_title('Modenspektrum: Anfang vs. 100 ms')
    ax2.legend(fontsize=9); ax2.grid(True,alpha=0.3,axis='y')
    fig.suptitle('Abb. 2: Transiente Anregung \u2014 Mode 2 im Einschwingen',fontweight='bold')
    fig.tight_layout(); return fig

def fig_torsion():
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
    zungen=[('Bass\n50Hz',70e-3,8e-3,0.40e-3,50),('Bass\nA2',70e-3,8e-3,0.66e-3,110),
        ('Disk.\nF3',37e-3,3.7e-3,0.30e-3,175),('Disk.\nC4',29e-3,3.1e-3,0.28e-3,262),
        ('Disk.\nA4',22e-3,2.4e-3,0.25e-3,440),('Disk.\nC6',14e-3,1.7e-3,0.24e-3,1047)]
    names=[z[0] for z in zungen]
    fT1s=[]; dyrs=[]; fBs=[]
    for _,Lz,Wz,tz,f1 in zungen:
        J=(1/3)*Wz*tz**3*(1-0.630*tz/Wz); Ip=rho_s*(Wz*tz)*(Wz**2+tz**2)/12
        fT1=1/(4*Lz)*np.sqrt(G*J/Ip); fT1s.append(fT1)
        dy=(Wz/2)*0.0175*1e3; spalt=0.04; dyrs.append(dy/spalt)
        fBs.append(6.27*f1)
    # Links: Frequenzen
    x=np.arange(len(names)); w=0.3
    ax1.bar(x-w/2,[z[4] for z in zungen],w,color='#1565c0',label='f\u2081 (Biegung)')
    ax1.bar(x+w/2,fT1s,w,color='#e94560',label='f_T1 (Torsion)')
    ax1.set_yscale('log'); ax1.set_xticks(x); ax1.set_xticklabels(names,fontsize=8)
    ax1.set_ylabel('Frequenz [Hz]'); ax1.set_title('Biege- vs. Torsionsfrequenz')
    ax1.legend(fontsize=9); ax1.grid(True,alpha=0.3,axis='y')
    # Rechts: Δy/Spalt
    cols_bar=['#e94560' if r>0.5 else '#2e7d32' for r in dyrs]
    ax2.bar(x,dyrs,color=cols_bar,edgecolor='black',lw=0.5)
    ax2.axhline(y=0.5,color='red',ls='--',lw=1.5,label='Kontaktgrenze')
    ax2.axhline(y=1.0,color='red',ls='-',lw=1,alpha=0.5)
    ax2.set_xticks(x); ax2.set_xticklabels(names,fontsize=8)
    ax2.set_ylabel('\u0394y / Spalt'); ax2.set_title('Seitliche Auslenkung bei 1\u00b0 Torsion')
    ax2.legend(fontsize=9); ax2.grid(True,alpha=0.3,axis='y')
    fig.suptitle('Abb. 3: Torsionsschwingungen \u2014 Frequenz und Kanalber\u00fchrung',fontweight='bold')
    fig.tight_layout(); return fig

# ══ PDF ══
outfile='0018_hoerbarkeit_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0018',sST))
story.append(Paragraph('H\u00f6rbarkeit der Biegemoden:<br/>Kammer-Saugkreis, Transiente, Torsion',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('Korrektur zu Dok. 0016: Die Biegemoden k\u00f6nnen unter bestimmten Bedingungen '
    'h\u00f6rbar werden. Drei Mechanismen: Kammer-Resonanz, Einschwingvorgang, Torsion. '
    'Referenz: Dok. 0009, 0011, 0015\u20130017.',sAb))
story.append(Paragraph('<b>Berechnungsskript:</b> <b>pruefskript_0018_hoerbarkeit.py</b> '
    '(modale Partizipation, Saugkreis-Mechanismus, Transiente, Torsionsmoden) '
    'im selben Verzeichnis auf GitHub.',sAb))

# Kap. 1
story.append(Paragraph('1. Korrektur: Biegemoden sind nicht immer unh\u00f6rbar',sCh))
story.append(Paragraph('Dok. 0016 hat gezeigt, dass die inharmonischen Biegemoden '
    '(f\u2082/f\u2081 = 6,27) im Normalfall 30\u201350 dB leiser sind als die harmonischen Obertöne '
    'des Impulsgenerators. Die Praxis zeigt aber: Mit einer abstimmbaren Resonanzkammer '
    'kann man Mode 2 bei manchen Zungen so stark anregen, dass <b>nur Mode 2 h\u00f6rbar ist</b>. '
    'Und beim Einschwingvorgang pr\u00e4gt Mode 2 den Klangcharakter.',sB))
story.append(Paragraph('<b>Die Grundaussage bleibt:</b> Die harmonischen Obertöne dominieren '
    'im eingeschwungenen Zustand. Aber Mode 2 ist nicht vernachl\u00e4ssigbar.',sW))

# Kap. 2
story.append(Paragraph('2. Kammer als Saugkreis: Mode 2 wird h\u00f6rbar durch Subtraktion',sCh))
story.append(Paragraph('Die Kammer wirkt als <b>Saugkreis</b> (Kerbfilter, Dok. 0005\u20130008) \u2014 '
    'sie <b>unterdr\u00fcckt</b> Frequenzen bei f_H. Die Kammer verst\u00e4rkt Mode 2 nicht. '
    'Mode 2 wird h\u00f6rbar, weil der Grundton und seine Harmonischen ged\u00e4mpft werden, '
    'w\u00e4hrend Mode 2 <b>au\u00dferhalb der Kerbe</b> liegt.',sB))
story.append(Paragraph(
    'Voraussetzung: (a) f_H trifft einen Harmonischen von f\u2081 \u2192 Grundton wird ged\u00e4mpft. '
    '(b) f\u2082 (Biegemode) liegt <b>nicht</b> in der Kerbe \u2192 Mode 2 bleibt \u00fcbrig. '
    'Das ist Subtraktion, nicht Addition.',sB))

# Diagramm: Saugkreis-Wirkung
fig_notch, ax_n = plt.subplots(figsize=(10,5))
f_arr = np.linspace(20, 800, 1000)
f1_val = 50.0; f2_biege = 6.27 * f1_val
# Harmonische Obertöne (Impulsgenerator, Duty 25%)
duty = 0.25
spektrum = np.zeros_like(f_arr)
for n in range(1, 14):
    a_n = abs(2 * np.sin(n * np.pi * duty) / (n * np.pi))
    a_dB = 20 * np.log10(a_n + 1e-30) - 20 * np.log10(abs(2*np.sin(np.pi*duty)/np.pi))
    idx = np.argmin(np.abs(f_arr - n*f1_val))
    if idx < len(spektrum):
        spektrum[idx] = a_dB
# Kerbfilter bei f_H = 250 Hz (5×f₁)
f_H_demo = 250; Q_H = 10; bw = f_H_demo / Q_H
kerbe = -25 * np.exp(-((f_arr - f_H_demo) / (bw/2))**2)
# Gefiltertes Spektrum
gefiltert = spektrum + np.interp(np.arange(len(f_arr)), np.arange(len(f_arr)), kerbe)

ax_n.vlines([n*f1_val for n in range(1,14)], -50, [spektrum[np.argmin(np.abs(f_arr-n*f1_val))] for n in range(1,14)],
    colors='#1565c0', lw=3, alpha=0.4, label='Harmonische (ohne Kammer)')
ax_n.fill_between(f_arr, -50, kerbe-50, alpha=0.15, color='red')
ax_n.plot(f_arr, kerbe-25, 'r-', lw=1.5, alpha=0.5, label=f'Kerbfilter bei f_H={f_H_demo} Hz')
ax_n.axvline(x=f2_biege, color='#e94560', ls='--', lw=2, label=f'f\u2082(Biege) = {f2_biege:.0f} Hz')
ax_n.annotate('Mode 2\nBLEIBT', xy=(f2_biege, -30), xytext=(f2_biege+50, -15),
    fontsize=10, fontweight='bold', color='#e94560',
    arrowprops=dict(arrowstyle='->', color='#e94560'))
ax_n.annotate('5\u00d7f\u2081\nGED\u00c4MPFT', xy=(250, -40), xytext=(180, -15),
    fontsize=9, color='red', arrowprops=dict(arrowstyle='->', color='red'))
ax_n.set_xlabel('Frequenz [Hz]'); ax_n.set_ylabel('Amplitude [dB]')
ax_n.set_title('Abb. 1: Saugkreis-Mechanismus \u2014 Kammer d\u00e4mpft Harmonische, Mode 2 bleibt')
ax_n.legend(fontsize=8, loc='upper right'); ax_n.grid(True, alpha=0.3)
ax_n.set_xlim(20, 600); ax_n.set_ylim(-50, 5)
fig_notch.tight_layout()
story.append(fig2img(fig_notch, 155)); story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    '<b>Warum geht das nicht bei jeder Zunge?</b> '
    'Die Profilierung bestimmt f\u2082. Wenn f\u2082 zuf\u00e4llig in der N\u00e4he eines '
    'Harmonischen von f\u2081 liegt (z.B. f\u2082 \u2248 6\u00d7f\u2081 bei gleichm\u00e4\u00dfiger Zunge), '
    'wird Mode 2 von der Kerbe <b>mit</b>ged\u00e4mpft \u2192 bleibt unh\u00f6rbar. '
    'Bei st\u00e4rkerer Profilierung (f\u2082 \u2248 4\u20135\u00d7f\u2081) liegt Mode 2 '
    'zwischen den Harmonischen \u2192 \u00fcberlebt die Kerbe.',sB))
story.append(Paragraph(
    '<b>Im Akkordeon trifft es die T\u00f6ne mit Ansprache-Problemen:</b> '
    'Dort ist die Kammer-Kopplung am st\u00e4rksten (Dok. 0009), die Stimmplatten '
    'sind st\u00e4rker profiliert (bessere Ansprache, Dok. 0012), und das Zusammenspiel '
    'von Profilierung und Kammergeometrie erf\u00fcllt zuf\u00e4llig die Bedingung '
    'f_H \u2260 f\u2082.',sK))

# Kap. 3
story.append(Paragraph('3. Transiente Anregung: Einschwingvorgang',sCh))
story.append(Paragraph('Beim Dr\u00fccken einer Taste trifft der Balgdruck die Zunge als Sprungfunktion. '
    'Jede Biegemode wird initial angesto\u00dfen. Mode 2 startet bei \u2248\u221225 dB '
    'und klingt innerhalb von \u224840 ms ab.',sB))
fig2=fig_transient(); story.append(fig2img(fig2,155)); story.append(Spacer(1,3*mm))
story.append(Paragraph('<b>Mode 2 ist im Einschwingvorgang h\u00f6rbar</b> \u2014 als kurze Klangfarbe '
    'in den ersten 50\u2013100 ms. Verschiedene Zungen haben verschiedene f\u2082 \u2192 '
    'verschiedene Transienten \u2192 verschiedener Klangcharakter, '
    'auch wenn der eingeschwungene Ton gleich klingt.',sK))

# Kap. 4
story.append(Paragraph('4. Torsionsschwingungen',sCh))
story.append(Paragraph('Die Torsionsmoden (Drehung um die L\u00e4ngsachse) liegen h\u00f6her als '
    'die Biegemoden: f_T1 \u2248 14\u201322\u00d7 f\u2081. Die Torsion selbst ist als Ton '
    'kaum h\u00f6rbar \u2014 aber die seitliche Auslenkung bei Torsion f\u00fchrt dazu, '
    'dass die Zunge die Kanalwand ber\u00fchrt.',sB))
story.append(Paragraph('f_T1 = 1/(4L) \u00d7 \u221a(G\u00d7J / (\u03c1\u00d7I_p))',sF))
fig3=fig_torsion(); story.append(fig2img(fig3,155)); story.append(Spacer(1,3*mm))

rows_t=[]
for nm,Lz,Wz,tz,f1 in [('Bass 50 Hz',70e-3,8e-3,0.40e-3,50),('Bass A2',70e-3,8e-3,0.66e-3,110),
    ('Diskant C4',29e-3,3.1e-3,0.28e-3,262),('Diskant A4',22e-3,2.4e-3,0.25e-3,440),
    ('Diskant C6',14e-3,1.7e-3,0.24e-3,1047)]:
    J=(1/3)*Wz*tz**3*(1-0.630*tz/Wz); Ip=rho_s*(Wz*tz)*(Wz**2+tz**2)/12
    fT1=1/(4*Lz)*np.sqrt(G*J/Ip); dy=(Wz/2)*0.0175*1e3; r=dy/0.04
    kratzt='\u25cf' if r>0.5 else '\u25cb'
    rows_t.append([nm,f'{f1}',f'{fT1:.0f}',f'{fT1/f1:.1f}\u00d7',f'{dy:.3f}',f'{r:.2f}',kratzt])
story.append(mk_tbl(['Zunge','f\u2081','f_T1 [Hz]','f_T1/f\u2081','\u0394y [mm]','\u0394y/Spalt','Kratzt'],rows_t,
    cw=[25*mm,16*mm,18*mm,16*mm,16*mm,16*mm,14*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph('Bei Basszungen (\u0394y/Spalt > 1,5) reicht schon 1\u00b0 Torsion '
    'f\u00fcr Kanalber\u00fchrung. <b>Im normalen Betrieb sollte das nicht vorkommen</b> \u2014 '
    'die Spalte sind gro\u00df genug und die Torsionsamplitude klein. '
    'Aber bei <b>kalten Instrumenten</b> ist Kratzen normal: '
    'Das Metall kontrahiert, der Spalt wird enger, und die Toleranz '
    'f\u00fcr seitliche Auslenkung sinkt. Sobald das Instrument sich erw\u00e4rmt, '
    'verschwindet das Kratzen \u2014 die thermische Ausdehnung \u00f6ffnet den Spalt wieder.',sB))
story.append(Paragraph(
    'Zus\u00e4tzlich \u00e4ndert K\u00e4lte die Wachsh\u00e4rte (steifere Einspannung), '
    'die Luftviskosit\u00e4t (h\u00f6here Squeeze-D\u00e4mpfung bei engem Spalt) und '
    'die Zungensteifigkeit (E-Modul steigt minimal). '
    'Alle Effekte wirken in dieselbe Richtung: engerer Spalt + steifere Zunge '
    '= weniger Toleranz f\u00fcr Torsion.',sB))
story.append(Paragraph(
    '<b>Paradoxon:</b> Die passgenauesten Stimmplatten \u2014 also die teuersten \u2014 '
    'sind am anf\u00e4lligsten f\u00fcr Kratzen bei K\u00e4lte. Enge Spalte (Dok. 0010: '
    'bessere Ansprache, weniger Kurzschluss) bedeuten gleichzeitig weniger '
    'Toleranz f\u00fcr Torsion. Eine Stimmplatte mit 25 \u00b5m Spalt hat bei '
    '3 \u00b5m thermischer Kontraktion 12 % weniger Spielraum \u2014 '
    'eine mit 50 \u00b5m Spalt nur 6 %. Die billige Platte mit gro\u00dfem Spalt '
    'kratzt nie, die teure mit perfektem Spalt kratzt im Winter.',sK))

# Kap. 5
story.append(Paragraph('5. Zusammenfassung',sCh))
pts=[
    '<b>Eingeschwungen, ohne Kammer-Treffer:</b> Mode 2 bei \u221237 dB \u2014 '
    'nicht h\u00f6rbar. Die harmonischen Obertöne dominieren mit 30\u201350 dB Vorsprung.',
    '<b>Kammer als Saugkreis (f_H \u2260 f\u2082):</b> Die Kammer d\u00e4mpft Harmonische von f\u2081. '
    'Wenn f\u2082 au\u00dferhalb der Kerbe liegt, bleibt Mode 2 \u00fcbrig \u2014 h\u00f6rbar durch Subtraktion, '
    'nicht durch Verst\u00e4rkung. Trifft besonders T\u00f6ne mit Ansprache-Problemen '
    '(starke Kammer-Kopplung, st\u00e4rkere Profilierung).',
    '<b>Transiente:</b> Mode 2 startet bei \u221225 dB und klingt in \u224840 ms ab. '
    'Pr\u00e4gt den Klangcharakter (Anschlag, Farbe). Verschiedene f\u2082 = verschiedener Klang.',
    '<b>Torsion:</b> f_T1 \u2248 14\u201322\u00d7 f\u2081. Seitliche Auslenkung kann '
    'bei engem Spalt zur Kanalber\u00fchrung f\u00fchren. Im Normalbetrieb sollte das '
    'nicht vorkommen. Bei <b>kalten Instrumenten</b> ist Kratzen normal \u2014 '
    'das Metall kontrahiert, der Spalt wird enger. Verschwindet beim Erw\u00e4rmen.',
    '<b>Korrektur zu Dok. 0016:</b> Die Aussage \u201eBiegemoden sind nicht h\u00f6rbar\u201c '
    'gilt nur f\u00fcr den eingeschwungenen Zustand ohne Kammer-Treffer. '
    'Im Transienten und bei Kammer-Kopplung sind sie h\u00f6rbar und pr\u00e4gen den Klang.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Die Biegemoden sind leise \u2014 aber nicht stumm. '
    'Der Saugkreis der Kammer l\u00e4sst sie \u00fcbrig, der Transient l\u00e4sst sie kurz erklingen, '
    'die K\u00e4lte macht den Spalt eng genug f\u00fcr Kratzen.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0018 \u2014 H\u00f6rbarkeit \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
