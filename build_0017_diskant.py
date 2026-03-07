#!/usr/bin/env python3
"""Dok. 0017 — Diskant-Stimmzungen: F3 bis C6"""
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
def beam_modes(N,tf,L,W,E,rho,nm=4):
    dx=L/N; nd=2*(N+1); K=np.zeros((nd,nd)); M=np.zeros((nd,nd))
    for e in range(N):
        xm=(e+.5)*dx; te=tf(xm); Ie=W*te**3/12; Ae=W*te; le=dx
        ke=E*Ie/le**3*np.array([[12,6*le,-12,6*le],[6*le,4*le**2,-6*le,2*le**2],[-12,-6*le,12,-6*le],[6*le,2*le**2,-6*le,4*le**2]])
        me=rho*Ae*le/420*np.array([[156,22*le,54,-13*le],[22*le,4*le**2,13*le,-3*le**2],[54,13*le,156,-22*le],[-13*le,-3*le**2,-22*le,4*le**2]])
        d=[2*e,2*e+1,2*e+2,2*e+3]
        for i in range(4):
            for j in range(4): K[d[i],d[j]]+=ke[i,j]; M[d[i],d[j]]+=me[i,j]
    fd=list(range(2,nd)); ev,_=eigh(K[np.ix_(fd,fd)],M[np.ix_(fd,fd)])
    return np.sqrt(np.abs(ev[:nm]))/(2*np.pi)

Es=200e9; rho_s=7800; Nn=50
nN=['C','C#','D','Eb','E','F','F#','G','G#','A','Bb','B']
f_F3=174.6; n_semi=31  # F3(0) bis C6(31)

def get_LW(s): f=s/31.0; return 37e-3*(14./37.)**f, 3.7e-3*(1.7/3.7)**f
def get_note(s):
    a=5+s; return nN[a%12]+str(3+a//12)
def make_prof(t0,vd,xm_f,L):
    xm=xm_f*L; tm=t0*vd; tt=tm+(t0-tm)*0.10
    def tf(x):
        if x<=xm:
            if xm>0: return t0-(t0-tm)*(1-np.cos(np.pi*x/xm))/2
            return tm
        else:
            r=L-xm
            if r>0: return tm+(tt-tm)*(1-np.cos(np.pi*(x-xm)/r))/2
            return tm
    return tf
def cal_t0(f_t,L,W,vd,xm_f):
    tl,th=0.02e-3,0.8e-3
    for _ in range(45):
        tm=(tl+th)/2; tf=make_prof(tm,vd,xm_f,L)
        fr=beam_modes(Nn,tf,L,W,Es,rho_s,nm=3)
        if fr[0]>f_t: th=tm
        else: tl=tm
    return (tl+th)/2

# Pre-compute all data
all_d=[]
for s in range(n_semi+1):
    ft=f_F3*2**(s/12); L,W=get_LW(s); note=get_note(s)
    frac=s/31.0; vd=0.72+0.16*frac; xm=0.55+0.40*frac
    t0=cal_t0(ft,L,W,vd,xm)
    tf=make_prof(t0,vd,xm,L); fr=beam_modes(Nn,tf,L,W,Es,rho_s,nm=4); rat=fr/fr[0]
    tt=tf(L); I=W*t0**3/12; k=3*Es*I/L**3; m=rho_s*W*t0*L*0.85*1e6
    all_d.append(dict(s=s,note=note,ft=ft,L=L,W=W,t0=t0,tt=tt,vd=vd,xm=xm,
                      f1=fr[0],r2=rat[1],f2=fr[1],f3=fr[2],k=k,m=m))

# ══ DIAGRAMME ══
def fig_mensur():
    fig,(a1,a2,a3)=plt.subplots(1,3,figsize=(13,4))
    ss=[d['s'] for d in all_d]; Ls=[d['L']*1e3 for d in all_d]
    Ws=[d['W']*1e3 for d in all_d]; ts=[d['t0']*1e3 for d in all_d]
    a1.plot(ss,Ls,'b-',lw=2); a1.set_ylabel('L [mm]'); a1.set_title('L\u00e4nge')
    a2.plot(ss,Ws,'g-',lw=2); a2.set_ylabel('W [mm]'); a2.set_title('Breite')
    a3.plot(ss,ts,'r-',lw=2,label='t\u2080 (Fu\u00df)')
    a3.plot(ss,[d['tt']*1e3 for d in all_d],'r--',lw=1.5,label='t_tip')
    a3.set_ylabel('t [mm]'); a3.set_title('Dicke (kalibriert)'); a3.legend(fontsize=8)
    ticks=[0,7,12,19,24,31]; labs=['F3','C4','F4','C5','F5','C6']
    for a in [a1,a2,a3]: a.set_xticks(ticks); a.set_xticklabels(labs,fontsize=8); a.grid(True,alpha=0.3)
    fig.suptitle('Abb. 1: Diskant-Mensur F3\u2013C6',fontweight='bold'); fig.tight_layout(); return fig

def fig_f2():
    fig,ax=plt.subplots(figsize=(10,5.5))
    ss=[d['s'] for d in all_d]; f2s=[d['f2'] for d in all_d]
    ax.plot(ss,f2s,'bo-',lw=2,ms=4,label='f\u2082 (profiliert)')
    ax.axhspan(1500,3000,alpha=0.08,color='purple'); ax.text(30,2200,'f_H Diskant',fontsize=9,color='purple',ha='right')
    ax.axhspan(3000,5000,alpha=0.08,color='blue'); ax.text(30,3800,'f_H klein',fontsize=9,color='blue',ha='right')
    ax.axvspan(8,24,alpha=0.05,color='red')
    ax.text(16,900,'KRITISCHE ZONE',ha='center',fontsize=10,color='red',fontweight='bold')
    ticks=[0,7,12,19,24,31]; labs=['F3','C4','F4','C5','F5','C6']
    ax.set_xticks(ticks); ax.set_xticklabels(labs,fontsize=8)
    ax.set_ylabel('f\u2082 [Hz]'); ax.set_title('Abb. 2: f\u2082 absolut vs. Diskantkammer f_H')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3); ax.set_ylim(500,7000)
    fig.tight_layout(); return fig

def fig_vergl():
    fig,(a1,a2)=plt.subplots(1,2,figsize=(11,5))
    ss=[d['s'] for d in all_d]; ks=[d['k'] for d in all_d]; ms=[d['m'] for d in all_d]
    a1.semilogy(ss,ks,'b-',lw=2); a1.set_ylabel('k [N/m]'); a1.set_title('Steifigkeit')
    a2.plot(ss,ms,'r-',lw=2); a2.set_ylabel('Masse [mg]'); a2.set_title('Zungenmasse')
    ticks=[0,7,12,19,24,31]; labs=['F3','C4','F4','C5','F5','C6']
    for a in [a1,a2]: a.set_xticks(ticks); a.set_xticklabels(labs,fontsize=8); a.grid(True,alpha=0.3)
    fig.suptitle('Abb. 3: Steifigkeit und Masse',fontweight='bold'); fig.tight_layout(); return fig

# ══ PDF ══
outfile='0017_diskant_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]
story.append(Spacer(1,10*mm)); story.append(Paragraph('Dok. 0017',sST))
story.append(Paragraph('Diskant-Stimmzungen:<br/>Obertonmoden F3 bis C6',sT))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph('FEM-Berechnung der Biegemoden von Diskant-Stimmzungen mit variabler '
    'Mensur (L = 37\u219214 mm, W = 3,7\u21921,7 mm, F3\u2013C6). Kalibrierte Dicke am Fu\u00df, '
    'Profilierung stimmt den Ton. Vergleich mit den Basszungen aus Dok. 0016.',sAb))
story.append(Paragraph('<b>Hinweis:</b> Dieses Dokument bietet ein Erkl\u00e4rungsmodell. '
    'Die interpolierte Mensur ist eine Sch\u00e4tzung \u2014 reale Stimmplatten variieren je nach Hersteller.',sW))
story.append(Paragraph(
    '<b>Berechnungsskript:</b> Die vollst\u00e4ndige FEM-Berechnung mit allen Tabellen '
    'ist im selben Verzeichnis auf GitHub abgelegt: '
    '<b>pruefskript_0017_diskant.py</b> (W-Unabh\u00e4ngigkeit, gleichm\u00e4\u00dfige Mindestdicke, '
    'kalibrierte Profilierung, f\u2082 vs. f_H, Vergleich Bass/Diskant).',sAb))

# Kap. 1
story.append(Paragraph('1. Die Breite W hat keinen Einfluss auf die Tonh\u00f6he',sCh))
story.append(Paragraph('f = \u03b2\u00b2/(2\u03c0L\u00b2) \u00d7 t \u00d7 \u221a(E/(12\u03c1))',sF))
story.append(Paragraph('Die Breite W k\u00fcrzt sich raus (I/A = t\u00b2/12). '
    'Die Tonh\u00f6he wird ausschlie\u00dflich durch <b>t</b> und <b>L</b> bestimmt. '
    'W beeinflusst Masse (\u221d W), Volumenstrom (\u221d W), Lautst\u00e4rke und Torsionssteifigkeit (\u221d W\u00b3). '
    'Deshalb nehmen W und L stetig ab von tief zu hoch \u2014 die tiefen T\u00f6ne brauchen '
    'mehr Volumenstrom.',sB))

# Kap. 2
story.append(Paragraph('2. Diskant-Mensur: F3 (175 Hz) bis C6 (1047 Hz)',sCh))
story.append(Paragraph('Die Diskantzungen sind immer leicht profiliert (12\u201328 % Verd\u00fcnnung) '
    'und am freien Ende am d\u00fcnnsten. Bei den tieferen T\u00f6nen liegt die d\u00fcnnste Stelle '
    'eher in der Mitte. Die Dicke am Fu\u00df (t\u2080) ist dicker als die Mindestdicke \u2014 '
    'die Profilierung hebt f\u2081 an und passt den Ton an.',sB))
fig1=fig_mensur(); story.append(fig2img(fig1,155)); story.append(Spacer(1,3*mm))

# Tabelle
sel=[0,7,12,19,24,31]  # F3, C4, F4, C5, F5, C6
rows_o=[]
for s in sel:
    if s >= len(all_d): continue
    d=all_d[s]
    rows_o.append([d['note'],f'{d["ft"]:.0f}',f'{d["L"]*1e3:.1f}',f'{d["W"]*1e3:.2f}',
        f'{d["t0"]*1e3:.3f}',f'{d["tt"]*1e3:.3f}',f'{(1-d["vd"])*100:.0f}',
        f'{d["r2"]:.2f}',f'{d["f2"]:.0f}'])
story.append(mk_tbl(['Ton','f [Hz]','L','W','t\u2080','t_tip','Vd%','f\u2082/f\u2081','f\u2082'],rows_o,
    cw=[13*mm,14*mm,12*mm,13*mm,14*mm,14*mm,11*mm,14*mm,16*mm]))
story.append(Spacer(1,3*mm))
story.append(Paragraph(f'<b>f\u2082/f\u2081 = {all_d[0]["r2"]:.2f}\u2013{all_d[-1]["r2"]:.2f}</b> (profiliert). '
    'Ohne Profilierung w\u00e4re f\u2082/f\u2081 = 6,27 (geometrische Cantilever-Konstante). '
    'Die Dicke am Fu\u00df variiert von {:.3f} mm (F3) bis {:.3f} mm (C6) \u2014 '
    'die h\u00f6heren T\u00f6ne sind <b>d\u00fcnner</b>, aber die Dicke nimmt weniger stark ab '
    'als die L\u00e4nge, weil f \u221d t/L\u00b2.'.format(all_d[0]['t0']*1e3, all_d[-1]['t0']*1e3),sK))

# Kap. 3
story.append(Paragraph('3. f\u2082 absolut: Kritische Zone',sCh))
fig2=fig_f2(); story.append(fig2img(fig2,155)); story.append(Spacer(1,3*mm))
story.append(Paragraph(f'f\u2082 steigt von {all_d[0]["f2"]:.0f} Hz (F3) auf '
    f'{all_d[-1]["f2"]:.0f} Hz (C6). Die kritische Zone \u2014 wo f\u2082 auf typische '
    'Diskantkammern (f_H = 1500\u20135000 Hz) trifft \u2014 liegt bei <b>C4 bis F5</b> '
    f'(f\u2082 \u2248 1400\u20134000 Hz). Die tiefsten Diskant-T\u00f6ne (F3\u2013B3) sind unkritisch, '
    'ebenso die h\u00f6chsten (G5\u2013C6, f\u2082 > 5000 Hz).',sB))

# Kap. 4
story.append(Paragraph('4. Steifigkeit und Masse',sCh))
fig3=fig_vergl(); story.append(fig2img(fig3,150)); story.append(Spacer(1,3*mm))
story.append(Paragraph(f'Die Steifigkeit steigt von {all_d[0]["k"]:.0f} N/m (F3) auf '
    f'{all_d[-1]["k"]:.0f} N/m (C6) \u2014 Faktor {all_d[-1]["k"]/all_d[0]["k"]:.0f}. '
    f'Die Masse sinkt von {all_d[0]["m"]:.0f} mg auf {all_d[-1]["m"]:.0f} mg.',sB))

# Kap. 5
story.append(Paragraph('5. Vergleich Bass \u2014 Diskant',sCh))
story.append(Paragraph(
    '<b>Gemeinsam:</b> f\u2082/f\u2081 = 6,27 ohne Profilierung. Cantilever-Physik identisch.',sB))
story.append(Paragraph(
    '<b>Bass:</b> Eine L\u00e4nge (70 mm), Dicke + Gewicht stimmen den Ton. '
    'Starke Profilierung (30\u201350 %), f\u2082/f\u2081 \u2248 4,0\u20135,3. '
    'Kritische Zone: A1\u2013A2 (f\u2082 = 345\u2013689 Hz \u2192 Basskammern).',sB))
story.append(Paragraph(
    '<b>Diskant:</b> L und W variieren stetig (37\u219214 mm). '
    'Leichte Profilierung (12\u201328 %), f\u2082/f\u2081 \u2248 5,3\u20135,8. '
    'Kritische Zone: C4\u2013F5 (f\u2082 \u2248 1400\u20134000 Hz \u2192 Diskantkammern). '
    'Tiefe und h\u00f6chste Diskant-T\u00f6ne sind unkritisch.',sB))

# Kap. 6
story.append(Paragraph('6. Zusammenfassung',sCh))
pts=[
    '<b>W hat keinen Einfluss auf den Ton.</b> f \u221d t/L\u00b2. W bestimmt Lautst\u00e4rke und Ansprache.',
    f'<b>Mensur F3\u2013C6:</b> L = 37\u219214 mm, W = 3,7\u21921,7 mm, '
    f't\u2080 = {all_d[0]["t0"]*1e3:.2f}\u2192{all_d[-1]["t0"]*1e3:.2f} mm (profiliert).',
    f'<b>f\u2082/f\u2081 = {all_d[0]["r2"]:.2f}\u2013{all_d[-1]["r2"]:.2f}</b> (profiliert). '
    'Die Profilierung senkt das Verh\u00e4ltnis um 0,5\u20131,0 gegen\u00fcber 6,27.',
    '<b>Kritische Zone: C4\u2013F5.</b> f\u2082 \u2248 1400\u20134000 Hz trifft die Diskantkammern. '
    'Tiefe (F3\u2013B3) und h\u00f6chste (G5\u2013C6) T\u00f6ne sind unkritisch.',
    f'<b>Steifigkeit:</b> Faktor {all_d[-1]["k"]/all_d[0]["k"]:.0f} von F3 bis C6.',
]
for i,p in enumerate(pts): story.append(Paragraph(f'{i+1}. {p}',sB))
story.append(Spacer(1,6*mm))
story.append(Paragraph('<i>Bass und Diskant gehorchen derselben Physik. '
    'Die Mensur \u00e4ndert die absoluten Frequenzen \u2014 die Verh\u00e4ltnisse bleiben.</i>',sAb))

def pn(c,d):
    c.saveState(); c.setFont('DejaVu',8); c.setFillColor(HexColor('#999999'))
    c.drawCentredString(WP/2,12*mm,f'Dok. 0017 \u2014 Diskant F3\u2013C6 \u2014 Seite {c.getPageNumber()}')
    c.restoreState()
doc.build(story,onFirstPage=pn,onLaterPages=pn)
print(f'\u2713 {outfile} erzeugt')
