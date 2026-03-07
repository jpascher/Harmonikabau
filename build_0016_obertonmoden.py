#!/usr/bin/env python3
"""
Dok. 0016 — Obertonmoden der Basszunge: Profilierung und Inharmonizität
Build-Skript: PDF mit FEM-Berechnung, Diagrammen und Tabellen.
Benötigt: reportlab, matplotlib, numpy, scipy
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.linalg import eigh
from io import BytesIO

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Image, PageBreak
)
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

pdfmetrics.registerFont(TTFont('DejaVu', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuB', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuI', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBI', '/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuMono', '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))
addMapping('DejaVu', 0, 0, 'DejaVu')
addMapping('DejaVu', 1, 0, 'DejaVuB')
addMapping('DejaVu', 0, 1, 'DejaVuI')
addMapping('DejaVu', 1, 1, 'DejaVuBI')

DARKBLUE=HexColor('#16213e'); ACCENTRED=HexColor('#e94560')
KEYGREEN=HexColor('#2e7d32'); WARNRED=HexColor('#c62828')
LIGHTGRAY=HexColor('#f5f5f5'); KEYBG=HexColor('#e8f5e9'); WARNBG=HexColor('#ffebee')
W_PAGE=A4[0]; PAGE_W=W_PAGE-50*mm

styles=getSampleStyleSheet()
for sn in styles.byName:
    s=styles.byName[sn]
    if hasattr(s,'fontName'):
        if 'Bold' in s.fontName: s.fontName='DejaVuB'
        elif 'Italic' in s.fontName: s.fontName='DejaVuI'
        else: s.fontName='DejaVu'

sTitle=ParagraphStyle('T',parent=styles['Title'],fontSize=18,textColor=DARKBLUE,spaceAfter=4,alignment=TA_CENTER,fontName='DejaVuB')
sSubtitle=ParagraphStyle('ST',parent=styles['Normal'],fontSize=12,textColor=DARKBLUE,alignment=TA_CENTER,spaceAfter=2,fontName='DejaVu')
sAbstract=ParagraphStyle('Ab',parent=styles['Italic'],fontSize=9,textColor=HexColor('#555555'),spaceAfter=8,fontName='DejaVuI')
sChapter=ParagraphStyle('Ch',parent=styles['Heading1'],fontSize=14,textColor=DARKBLUE,spaceBefore=14,spaceAfter=6,fontName='DejaVuB')
sBody=ParagraphStyle('Bo',parent=styles['Normal'],fontSize=10,leading=14,spaceAfter=6,alignment=TA_JUSTIFY,fontName='DejaVu')
sFormula=ParagraphStyle('Fo',parent=sBody,fontSize=10,alignment=TA_CENTER,spaceAfter=8,spaceBefore=4,fontName='DejaVuMono')
sKeyBox=ParagraphStyle('KB',parent=sBody,fontSize=10,backColor=KEYBG,borderPadding=6,borderColor=KEYGREEN,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sWarnBox=ParagraphStyle('WB',parent=sBody,fontSize=10,backColor=WARNBG,borderPadding=6,borderColor=WARNRED,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sTH=ParagraphStyle('TH',parent=sBody,fontSize=9,fontName='DejaVuB',alignment=TA_CENTER,leading=11)
sTD=ParagraphStyle('TD',parent=sBody,fontSize=9,alignment=TA_CENTER,leading=11,fontName='DejaVu')
sTDL=ParagraphStyle('TDL',parent=sBody,fontSize=9,alignment=TA_LEFT,leading=11,fontName='DejaVu')

def hr():
    return Table([['']], colWidths=[PAGE_W], style=TableStyle([('LINEBELOW',(0,0),(-1,-1),2,ACCENTRED),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))

def make_table(header, rows, col_widths=None):
    data=[[Paragraph(h,sTH) for h in header]]
    for row in rows:
        data.append([Paragraph(str(c),sTDL) if i==0 else Paragraph(str(c),sTD) for i,c in enumerate(row)])
    cw=col_widths or [PAGE_W/len(header)]*len(header)
    t=Table(data,colWidths=cw,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),DARKBLUE),('TEXTCOLOR',(0,0),(-1,0),white),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,LIGHTGRAY]),('GRID',(0,0),(-1,-1),0.5,HexColor('#cccccc')),
        ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),
        ('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
    return t

def fig_to_image(fig, width_mm=160, dpi=200):
    buf=BytesIO(); fig.savefig(buf,format='png',dpi=dpi,bbox_inches='tight'); plt.close(fig)
    buf.seek(0); img=Image(buf); r=img.imageWidth/img.imageHeight
    img.drawWidth=width_mm*mm; img.drawHeight=width_mm*mm/r; return img

plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10

# ══════════════════════════════════════════════════════════════
# FEM-LÖSUNG
# ══════════════════════════════════════════════════════════════
def beam_modes(N, t_func, L, W, E, rho, n_modes=6):
    dx=L/N; ndof=2*(N+1)
    K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for e in range(N):
        xm=(e+0.5)*dx; te=t_func(xm); Ie=W*te**3/12; Ae=W*te; le=dx
        ke=E*Ie/le**3*np.array([[12,6*le,-12,6*le],[6*le,4*le**2,-6*le,2*le**2],[-12,-6*le,12,-6*le],[6*le,2*le**2,-6*le,4*le**2]])
        me=rho*Ae*le/420*np.array([[156,22*le,54,-13*le],[22*le,4*le**2,13*le,-3*le**2],[54,13*le,156,-22*le],[-13*le,-3*le**2,-22*le,4*le**2]])
        d=[2*e,2*e+1,2*e+2,2*e+3]
        for i in range(4):
            for j in range(4):
                K[d[i],d[j]]+=ke[i,j]; M[d[i],d[j]]+=me[i,j]
    fd=list(range(2,ndof))
    ev,evec=eigh(K[np.ix_(fd,fd)],M[np.ix_(fd,fd)])
    freqs=np.sqrt(np.abs(ev))/(2*np.pi)
    modes=np.zeros((N+1,n_modes))
    for m in range(min(n_modes,len(fd))):
        full=np.zeros(ndof); full[fd]=evec[:,m]; modes[:,m]=full[0::2]
        mx=np.max(np.abs(modes[:,m]))
        if mx>0: modes[:,m]/=mx
    return freqs[:n_modes], modes

def beam_modes_weighted(N, t_func, L, W, E, rho, m_tip=0, n_modes=6):
    """FEM Cantilever mit Spitzengewicht m_tip [kg]"""
    dx=L/N; ndof=2*(N+1)
    K=np.zeros((ndof,ndof)); M=np.zeros((ndof,ndof))
    for e in range(N):
        xm=(e+0.5)*dx; te=t_func(xm); Ie=W*te**3/12; Ae=W*te; le=dx
        ke=E*Ie/le**3*np.array([[12,6*le,-12,6*le],[6*le,4*le**2,-6*le,2*le**2],
            [-12,-6*le,12,-6*le],[6*le,2*le**2,-6*le,4*le**2]])
        me=rho*Ae*le/420*np.array([[156,22*le,54,-13*le],[22*le,4*le**2,13*le,-3*le**2],
            [54,13*le,156,-22*le],[-13*le,-3*le**2,-22*le,4*le**2]])
        d=[2*e,2*e+1,2*e+2,2*e+3]
        for i in range(4):
            for j in range(4):
                K[d[i],d[j]]+=ke[i,j]; M[d[i],d[j]]+=me[i,j]
    if m_tip > 0:
        M[2*N, 2*N] += m_tip
    fd=list(range(2,ndof))
    ev,evec=eigh(K[np.ix_(fd,fd)],M[np.ix_(fd,fd)])
    freqs=np.sqrt(np.abs(ev))/(2*np.pi)
    modes=np.zeros((N+1,n_modes))
    for m in range(min(n_modes,len(fd))):
        full=np.zeros(ndof); full[fd]=evec[:,m]; modes[:,m]=full[0::2]
        mx=np.max(np.abs(modes[:,m]))
        if mx>0: modes[:,m]/=mx
    return freqs[:n_modes], modes

E_s=200e9; rho_s=7800; L_b=70e-3; W_b=8e-3; t0=0.40e-3
N_el=80; n_m=6; x_n=np.linspace(0,L_b,N_el+1)
m_zunge = rho_s * W_b * t0 * L_b  # Zungenmasse

profiles = [
    ('Gleichmäßig (t = const)',    lambda x: t0,                                '#1565c0'),
    ('Linear 80 % (t_tip = 0,32)', lambda x: t0*(1-0.20*x/L_b),               '#2e7d32'),
    ('Linear 50 % (t_tip = 0,20)', lambda x: t0*(1-0.50*x/L_b),               '#ff9800'),
    ('Linear 30 % (t_tip = 0,12)', lambda x: t0*(1-0.70*x/L_b),               '#e94560'),
    ('Parabolisch 50 %',           lambda x: t0*(1-0.50*(x/L_b)**2),           '#9c27b0'),
    ('Umgekehrt parabolisch',      lambda x: t0*(1-0.50*np.sqrt(x/L_b+1e-12)), '#795548'),
]

results={}
for name, tf, col in profiles:
    fr, mo = beam_modes(N_el, tf, L_b, W_b, E_s, rho_s, n_m)
    results[name]={'freqs':fr,'modes':mo,'t_func':tf,'color':col}

# ══════════════════════════════════════════════════════════════
# DIAGRAMME
# ══════════════════════════════════════════════════════════════

def make_fig1():
    fig,ax=plt.subplots(figsize=(10,4))
    for name,tf,col in profiles:
        tv=[tf(x)*1000 for x in x_n]
        ax.plot(x_n*1000,tv,'-',color=col,lw=2,label=name)
    ax.set_xlabel('Position x [mm] (0 = Einspannung)')
    ax.set_ylabel('Dicke t(x) [mm]')
    ax.set_title('Abb. 1: Dickenprofile (L = 70 mm, W = 8 mm, t\u2080 = 0,40 mm)')
    ax.legend(fontsize=8,loc='lower left'); ax.grid(True,alpha=0.3); ax.set_xlim(0,70)
    fig.tight_layout(); return fig

def make_fig2():
    fig,ax=plt.subplots(figsize=(10,5))
    x_pos=np.arange(n_m); w=0.12; off=0
    for name,tf,col in profiles:
        d=results[name]; rat=d['freqs']/d['freqs'][0]
        ax.bar(x_pos+off,rat,w,label=name,color=col,alpha=0.85); off+=w
    for n in range(n_m):
        ax.axhline(y=n+1,color='gray',ls=':',lw=0.8,alpha=0.5)
    ax.set_xticks(x_pos+w*2.5)
    ax.set_xticklabels([f'Mode {i+1}' for i in range(n_m)])
    ax.set_ylabel('f_n / f\u2081'); ax.set_title('Abb. 2: Obertonverhältnisse')
    ax.legend(fontsize=7,loc='upper left'); ax.grid(True,alpha=0.3,axis='y')
    fig.tight_layout(); return fig

def make_fig3():
    fig,axes=plt.subplots(3,1,figsize=(10,9),sharex=True)
    sel=[('Gleichmäßig (t = const)','#1565c0'),('Linear 50 % (t_tip = 0,20)','#ff9800'),('Linear 30 % (t_tip = 0,12)','#e94560')]
    for ax,( pn,pc) in zip(axes,sel):
        d=results[pn]
        for m in range(4):
            ax.plot(x_n*1000,d['modes'][:,m],'-',lw=1.5,label=f'Mode {m+1} (f/f\u2081 = {d["freqs"][m]/d["freqs"][0]:.2f})')
        ax.axhline(y=0,color='black',lw=0.5); ax.set_ylabel('Auslenkung')
        ax.set_title(pn,fontsize=11,fontweight='bold',color=pc)
        ax.legend(fontsize=8,loc='upper left'); ax.grid(True,alpha=0.3)
        tv=np.array([d['t_func'](x) for x in x_n])/t0
        ax.fill_between(x_n*1000,-1.2,-1.2+0.3*tv,alpha=0.15,color=pc)
    axes[-1].set_xlabel('Position x [mm]')
    fig.suptitle('Abb. 3: Modenformen \u2014 Profilierung verschiebt die Knoten',fontsize=12,fontweight='bold')
    fig.tight_layout(); return fig

def make_fig4():
    fig,ax=plt.subplots(figsize=(10,5))
    for name,tf,col in profiles:
        d=results[name]; rat=d['freqs']/d['freqs'][0]
        dev=1200*np.log2(rat/np.arange(1,n_m+1))
        ax.plot(range(1,n_m+1),dev,'o-',color=col,lw=1.5,ms=6,label=name)
    ax.axhline(y=0,color='black',lw=1)
    ax.set_xlabel('Oberton-Nummer n'); ax.set_ylabel('Abweichung von n \u00d7 f\u2081 [Cent]')
    ax.set_title('Abb. 4: Inharmonizität \u2014 alle Profile weit von harmonisch')
    ax.legend(fontsize=7); ax.grid(True,alpha=0.3); ax.set_xticks(range(1,n_m+1))
    fig.tight_layout(); return fig

def make_fig5():
    """Steifigkeit und Grundfrequenz vs. Profilierung"""
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(11,4.5))
    taper=np.linspace(0,0.85,50)  # t_tip/t0 Verhältnis
    f1_vals=[]; k_vals=[]
    for tp in taper:
        tf_v=lambda x,tp=tp: t0*(1-tp*x/L_b)
        fr,_=beam_modes(N_el,tf_v,L_b,W_b,E_s,rho_s,1)
        f1_vals.append(fr[0])
        # Steifigkeit: k = 3EI_eff/L³ (Näherung)
        I_tip=W_b*(t0*(1-tp))**3/12
        I_root=W_b*t0**3/12
        I_eff=(I_root+I_tip)/2  # grobe Mittelung
        k_vals.append(3*E_s*I_eff/L_b**3)
    
    ax1.plot(taper*100,f1_vals,'b-',lw=2)
    ax1.set_xlabel('Verdünnung an der Spitze [%]')
    ax1.set_ylabel('Grundfrequenz f\u2081 [Hz]')
    ax1.set_title('f\u2081 vs. Profilierung')
    ax1.grid(True,alpha=0.3)
    ax1.axhline(y=f1_vals[0],color='gray',ls=':')
    
    k0=k_vals[0]
    ax2.plot(taper*100,[k/k0*100 for k in k_vals],'r-',lw=2)
    ax2.set_xlabel('Verdünnung an der Spitze [%]')
    ax2.set_ylabel('Steifigkeit k/k\u2080 [%]')
    ax2.set_title('Steifigkeit vs. Profilierung')
    ax2.grid(True,alpha=0.3)
    ax2.axhline(y=100,color='gray',ls=':')
    
    fig.suptitle('Abb. 5: Profilierung ändert f\u2081 und k gleichzeitig',fontweight='bold')
    fig.tight_layout(); return fig

def make_fig6():
    """Mode 2 Frequenz relativ zu typischen f_H-Werten"""
    fig,ax=plt.subplots(figsize=(10,5))
    taper=np.linspace(0,0.80,40)
    f2_vals=[]
    for tp in taper:
        tf_v=lambda x,tp=tp: t0*(1-tp*x/L_b)
        fr,_=beam_modes(N_el,tf_v,L_b,W_b,E_s,rho_s,4)
        f2_vals.append(fr[:4])
    f2_vals=np.array(f2_vals)
    
    mode_labels=['Mode 1','Mode 2','Mode 3','Mode 4']
    mode_colors=['#1565c0','#2e7d32','#ff9800','#e94560']
    for m in range(4):
        ax.plot(taper*100,f2_vals[:,m],'-',color=mode_colors[m],lw=2,label=mode_labels[m])
    
    # Typische f_H-Bereiche
    ax.axhspan(200,600,alpha=0.08,color='purple')
    ax.text(82,400,'f_H Bereich\n(typisch)',fontsize=9,color='purple',ha='right')
    
    ax.set_xlabel('Verdünnung an der Spitze [%]')
    ax.set_ylabel('Frequenz [Hz]')
    ax.set_title('Abb. 6: Modenfrequenzen vs. Profilierung \u2014 Kammer-Kopplung')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3)
    ax.set_yscale('log'); ax.set_ylim(30,30000)
    fig.tight_layout(); return fig

def make_fig7():
    """Gewicht an der Spitze: f1 und f2/f1 vs. Masse"""
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
    masses=np.linspace(0,1.5e-3,30)
    t_c=lambda x: t0
    f1s=[]; r2s=[]; r3s=[]
    for mt in masses:
        fr,_=beam_modes_weighted(N_el,t_c,L_b,W_b,E_s,rho_s,m_tip=mt,n_modes=4)
        f1s.append(fr[0]); r2s.append(fr[1]/fr[0]); r3s.append(fr[2]/fr[0])
    
    ax1.plot(np.array(masses)*1000,f1s,'b-',lw=2)
    ax1.set_xlabel('Spitzengewicht [g]'); ax1.set_ylabel('f\u2081 [Hz]')
    ax1.set_title('Grundfrequenz sinkt')
    ax1.grid(True,alpha=0.3)
    ax1.axhline(y=50,color='green',ls='--',lw=1,label='Zielfrequenz 50 Hz')
    ax1.legend(fontsize=9)
    # Masse bei f1=50
    f1a=np.array(f1s); ma=np.array(masses)*1000
    idx50=np.argmin(np.abs(f1a-50))
    ax1.plot(ma[idx50],50,'go',ms=10)
    ax1.text(ma[idx50]+0.05,52,f'{ma[idx50]:.2f} g',fontsize=9,color='green')
    
    ax2.plot(np.array(masses)*1000,r2s,'r-',lw=2,label='f\u2082/f\u2081')
    ax2.plot(np.array(masses)*1000,r3s,'b--',lw=1.5,label='f\u2083/f\u2081')
    ax2.set_xlabel('Spitzengewicht [g]'); ax2.set_ylabel('Frequenzverhältnis')
    ax2.set_title('Obertöne spreizen sich!')
    ax2.grid(True,alpha=0.3); ax2.legend(fontsize=9)
    ax2.axhline(y=6.27,color='gray',ls=':',lw=1)
    ax2.text(1.2,6.0,'gleichm. ohne Gew.',fontsize=8,color='gray')
    
    fig.suptitle('Abb. 7: Spitzengewicht \u2014 gegenteiliger Effekt zur Verdünnung',fontweight='bold')
    fig.tight_layout(); return fig

def make_fig8():
    """Kombination Profilierung + Gewicht: Vergleichstabelle als Balkendiagramm"""
    fig,ax=plt.subplots(figsize=(12,5))
    
    combos=[
        ('Gleichm.\nohne Gew.',lambda x:t0,0,'#1565c0'),
        ('Gleichm.\n+0,20 g',lambda x:t0,0.20e-3,'#42a5f5'),
        ('Lin. 50%\nohne Gew.',lambda x:t0*(1-0.50*x/L_b),0,'#ff9800'),
        ('Lin. 50%\n+0,10 g',lambda x:t0*(1-0.50*x/L_b),0.10e-3,'#ffb74d'),
        ('Lin. 50%\n+0,20 g',lambda x:t0*(1-0.50*x/L_b),0.20e-3,'#ffe0b2'),
        ('Lin. 30%\nohne Gew.',lambda x:t0*(1-0.70*x/L_b),0,'#e94560'),
        ('Lin. 30%\n+0,10 g',lambda x:t0*(1-0.70*x/L_b),0.10e-3,'#ef9a9a'),
        ('Lin. 30%\n+0,20 g',lambda x:t0*(1-0.70*x/L_b),0.20e-3,'#ffcdd2'),
    ]
    
    names=[]; f1v=[]; r2v=[]
    for nm,tf,mt,co in combos:
        fr,_=beam_modes_weighted(N_el,tf,L_b,W_b,E_s,rho_s,m_tip=mt,n_modes=3)
        names.append(nm); f1v.append(fr[0]); r2v.append(fr[1]/fr[0])
    
    x_pos=np.arange(len(names))
    colors_c=[c[3] for c in combos]
    
    bars=ax.bar(x_pos,r2v,color=colors_c,edgecolor='black',lw=0.5)
    ax.set_xticks(x_pos); ax.set_xticklabels(names,fontsize=8)
    ax.set_ylabel('f\u2082 / f\u2081')
    
    # f1 als Text über den Balken
    for i,(bar,f1) in enumerate(zip(bars,f1v)):
        ax.text(bar.get_x()+bar.get_width()/2,bar.get_height()+0.05,
            f'f\u2081={f1:.0f} Hz',ha='center',fontsize=8,rotation=45)
    
    ax.axhline(y=6.27,color='gray',ls=':',lw=1)
    ax.text(7.5,6.35,'gleichm. Referenz',fontsize=8,color='gray',ha='right')
    ax.set_title('Abb. 8: Kombination Profilierung + Gewicht')
    ax.grid(True,alpha=0.3,axis='y')
    fig.tight_layout(); return fig

def make_fig9():
    """Gegenüberstellung: Verdünnung vs. Gewicht bei gleichem f1-Ziel"""
    fig,ax=plt.subplots(figsize=(10,5))
    
    # Ziel: f1 ≈ 50 Hz erreichen auf verschiedene Weisen
    # 1. Gleichmäßig + Gewicht → wie viel Gewicht für 50 Hz?
    # 2. Linear 50% + Gewicht
    # 3. Linear 30% + Gewicht
    
    configs_scan=[
        ('Gleichmäßig + Gewicht',lambda x:t0,'#1565c0'),
        ('Linear 50% + Gewicht',lambda x:t0*(1-0.50*x/L_b),'#ff9800'),
        ('Linear 30% + Gewicht',lambda x:t0*(1-0.70*x/L_b),'#e94560'),
    ]
    
    for nm,tf,co in configs_scan:
        masses_s=np.linspace(0,1.0e-3,40)
        f1s_s=[]; r2s_s=[]
        for mt in masses_s:
            fr,_=beam_modes_weighted(N_el,tf,L_b,W_b,E_s,rho_s,m_tip=mt,n_modes=3)
            f1s_s.append(fr[0]); r2s_s.append(fr[1]/fr[0])
        ax.plot(f1s_s,r2s_s,'o-',color=co,lw=2,ms=3,label=nm)
    
    # 50 Hz markieren
    ax.axvline(x=50,color='green',ls='--',lw=1.5,alpha=0.7)
    ax.text(51,4.5,'Ziel: f\u2081 = 50 Hz',fontsize=9,color='green')
    
    # Harmonische Linie
    ax.axhline(y=2,color='gray',ls=':',lw=1,alpha=0.5)
    ax.text(75,2.1,'harmonisch (2\u00d7)',fontsize=8,color='gray')
    
    ax.set_xlabel('Grundfrequenz f\u2081 [Hz]')
    ax.set_ylabel('f\u2082 / f\u2081')
    ax.set_title('Abb. 9: Gleiches f\u2081-Ziel, verschiedene Wege \u2014 Verdünnung gewinnt')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3)
    ax.set_xlim(30,80); ax.set_ylim(3,11)
    fig.tight_layout(); return fig

# ── Praxis-Profile: Glatte Mulde mit wanderndem Minimum ──
def make_smooth_profile(x_min_frac, t_min_ratio):
    xm=x_min_frac*L_b; tm=t0*t_min_ratio; tt=tm+(t0-tm)*0.15
    def tf(x):
        if x<=xm:
            if xm>0: return t0-(t0-tm)*(1-np.cos(np.pi*x/xm))/2
            return tm
        else:
            r=L_b-xm
            if r>0: return tm+(tt-tm)*(1-np.cos(np.pi*(x-xm)/r))/2
            return tm
    return tf

def make_fig10():
    """Praxis-Profile und f₂/f₁ vs. Position"""
    fig,(ax1,ax2)=plt.subplots(1,2,figsize=(12,5))
    x_plot=np.linspace(0,L_b,200)
    positions=[1.0,0.9,0.8,0.7,0.6,0.5,0.4,0.3]
    cols=plt.cm.viridis(np.linspace(0,0.9,len(positions)))
    for xf,col in zip(positions,cols):
        tf=make_smooth_profile(xf,0.50)
        t_vals=[tf(x)*1000 for x in x_plot]
        ax1.plot(x_plot*1000,t_vals,'-',color=col,lw=1.5,label=f'x_min={xf:.1f}L')
        ax1.plot(xf*L_b*1000,tf(xf*L_b)*1000,'o',color=col,ms=4)
    ax1.set_xlabel('Position x [mm]'); ax1.set_ylabel('Dicke t(x) [mm]')
    ax1.set_title('Profilformen'); ax1.legend(fontsize=7,ncol=2); ax1.grid(True,alpha=0.3); ax1.set_xlim(0,70)
    for depth,col,ls in [(0.30,'#2e7d32','--'),(0.50,'#ff9800','-'),(0.70,'#e94560','-.')]:
        r2v=[]
        for xf in np.linspace(0.15,1.0,25):
            tf=make_smooth_profile(xf,depth)
            fr,_=beam_modes_weighted(N_el,tf,L_b,W_b,E_s,rho_s,n_modes=3)
            r2v.append(fr[1]/fr[0])
        ax2.plot(np.linspace(0.15,1.0,25),r2v,ls,color=col,lw=2,label=f'Verd. {(1-depth)*100:.0f}%')
    ax2.axhline(y=6.27,color='gray',ls=':',lw=1,label='gleichm. (6,27)')
    ax2.set_xlabel('Position d\u00fcnnste Stelle x_min/L'); ax2.set_ylabel('f\u2082/f\u2081')
    ax2.set_title('f\u2082/f\u2081 vs. Position'); ax2.legend(fontsize=8); ax2.grid(True,alpha=0.3)
    fig.suptitle('Abb. 10: Praxis-Profilierung \u2014 Position vs. Tiefe',fontweight='bold')
    fig.tight_layout(); return fig

def make_fig11():
    """Verdünnungstiefe: f₂/f₁ hat ein Optimum"""
    fig,ax=plt.subplots(figsize=(9,5))
    depths=np.linspace(0.10,0.95,30)
    for xf,col,lab in [(0.7,'#ff9800','x_min=0,7L'),(0.8,'#2e7d32','x_min=0,8L'),(1.0,'#1565c0','Spitze')]:
        r2v=[]
        for dp in depths:
            if t0*dp<0.01e-3: r2v.append(np.nan); continue
            tf=make_smooth_profile(xf,dp)
            fr,_=beam_modes_weighted(N_el,tf,L_b,W_b,E_s,rho_s,n_modes=3)
            r2v.append(fr[1]/fr[0])
        ax.plot((1-depths)*100,r2v,'-',color=col,lw=2,label=lab)
    ax.axhline(y=6.27,color='gray',ls=':',lw=1)
    ax.set_xlabel('Verd\u00fcnnung [%]'); ax.set_ylabel('f\u2082/f\u2081')
    ax.set_title('Abb. 11: Optimum bei \u224870% Verd\u00fcnnung')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3)
    fig.tight_layout(); return fig

def make_fig12():
    """Chromatisch: f₂ absolut vs. f_H über 2 Oktaven"""
    fig,ax=plt.subplots(figsize=(10,5))
    f_A1=55.0; semi=np.arange(25)
    fc=[f_A1*2**(s/12) for s in semi]
    ax.plot(semi,[6.27*f for f in fc],'bo-',lw=2,ms=5,label='Gleichm. (6,27\u00d7f\u2081)')
    ax.plot(semi,[4.5*f for f in fc],'rs-',lw=2,ms=5,label='Profiliert (\u22484,5\u00d7f\u2081)')
    ax.axhspan(200,350,alpha=0.1,color='purple'); ax.axhspan(300,500,alpha=0.1,color='blue')
    ax.axhspan(500,700,alpha=0.1,color='green')
    ax.text(24.5,270,'f_H gro\u00dfe Kammer',fontsize=8,color='purple',ha='right')
    ax.text(24.5,400,'f_H mittlere',fontsize=8,color='blue',ha='right')
    ax.text(24.5,600,'f_H kleine',fontsize=8,color='green',ha='right')
    nl=['A1','','B1','C2','','D2','','E2','F2','','G2','','A2','','B2','C3','','D3','','E3','F3','','G3','','A3']
    ax.set_xticks(semi); ax.set_xticklabels(nl,fontsize=7,rotation=45)
    ax.set_ylabel('f\u2082 [Hz]'); ax.set_title('Abb. 12: f\u2082 \u00fcber 2 Oktaven \u2014 Kritische Zone')
    ax.legend(fontsize=9); ax.grid(True,alpha=0.3); ax.set_ylim(100,1500)
    ax.axvspan(0,12,alpha=0.05,color='red')
    ax.text(6,1400,'KRITISCHE ZONE',ha='center',fontsize=10,color='red',fontweight='bold')
    fig.tight_layout(); return fig

outfile='0016_obertonmoden_De.pdf'
doc=SimpleDocTemplate(outfile,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=20*mm,bottomMargin=20*mm)
story=[]

story.append(Spacer(1,10*mm))
story.append(Paragraph('Dok. 0016',sSubtitle))
story.append(Paragraph('Obertonmoden der Basszunge:<br/>Profilierung und Inharmonizität',sTitle))
story.append(Spacer(1,3*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'FEM-Berechnung der Biegemoden einer Basszunge (L = 70 mm) mit '
    'sechs verschiedenen Dickenprofilen. Modenformen, Frequenzverhältnisse, '
    'Inharmonizität und Konsequenzen für Klang und Kammerkopplung. '
    'Referenz: Dok. 0005–0009, 0012, 0015.',sAbstract))
story.append(Paragraph(
    '<b>Hinweis:</b> Dieses Dokument bietet ein Erklärungsmodell. '
    'Die FEM-Berechnung erfasst die Biegemoden eines idealisierten Balkens — '
    'Strömungskopplung, nichtlineare Amplitude und Spalteffekte sind nicht enthalten.',sWarnBox))
story.append(Paragraph(
    '<b>Berechnungsskripte:</b> Die vollständigen FEM-Berechnungen mit allen Tabellen '
    'sind in zwei separaten Python-Skripten dokumentiert, die im selben Verzeichnis '
    'auf GitHub abgelegt sind: '
    '<b>pruefskript_0016_obertonmoden.py</b> (Profile, Gewichte, Steifigkeit, Knotenpositionen) '
    'und <b>pruefskript_0016_chromatisch.py</b> (chromatische Analyse über 2 Oktaven, '
    'f₂/f₁ = const bei gleichmäßigem Profil, f₂ absolut vs. f_H-Bereiche) '
    'und <b>pruefskript_0016_praxis.py</b> (Praxis-Profilierung mit wandernder '
    'dünnster Stelle, Verdünnungstiefe, Optimum-Suche).',sAbstract))

# ══ Kap. 1 ══
story.append(Paragraph('1. Obertöne eines Cantilevers sind nicht harmonisch',sChapter))
story.append(Paragraph(
    'Die Stimmzunge ist ein einseitig eingespannter Balken (Cantilever). '
    'Die Eigenfrequenzen der Biegemoden folgen nicht der harmonischen Reihe '
    '(1 : 2 : 3 : 4 ...), sondern der Cantilever-Reihe:',sBody))
story.append(Paragraph('f_n / f₁ = (β_n / β₁)²',sFormula))
story.append(Paragraph(
    'mit β₁L = 1,8751, β₂L = 4,6941, β₃L = 7,8548, β₄L = 10,9955. '
    'Das ergibt die Verhältnisse 1 : 6,27 : 17,55 : 34,39 — '
    'der 2. Oberton liegt nicht bei der Oktave (2×), sondern bei mehr als '
    'zweieinhalb Oktaven (6,27×). Der 3. Oberton bei fast viereinhalb Oktaven (17,55×).',sBody))
story.append(Paragraph(
    '<b>Die Harmonischen, die man im Klang der Zunge hört, stammen nicht '
    'von diesen Biegemoden.</b> Sie entstehen durch die periodische '
    'Unterbrechung des Luftstroms — die Zunge wirkt als Impulsgenerator '
    '(Dok. 0013, 0015). Die Impulsform erzeugt harmonische Obertöne '
    '(1f, 2f, 3f ...), auch wenn die mechanischen Moden der Zunge '
    'ganz woanders liegen.',sKeyBox))
story.append(Paragraph(
    'Warum sind die mechanischen Moden trotzdem wichtig? Weil sie bestimmen, '
    'bei welchen Frequenzen die Zunge <b>resonant angeregt</b> werden kann — '
    'durch die Kammer (Dok. 0009), durch Nachbarzungen, oder durch '
    'transiente Druckschwankungen beim Spielen.',sBody))

# ══ Kap. 2 ══
story.append(Paragraph('2. Sechs Dickenprofile',sChapter))
story.append(Paragraph(
    'Die Profilierung (Verdünnung zur Spitze hin) verändert die Masseverteilung '
    'und die lokale Biegesteifigkeit entlang der Zunge. '
    'Sechs Profile werden verglichen — alle mit demselben Fußquerschnitt '
    '(t₀ = 0,40 mm, W = 8 mm):',sBody))

fig1=make_fig1()
story.append(fig_to_image(fig1,width_mm=155))
story.append(Spacer(1,3*mm))

# Tabelle: Profile
rows_p=[]
for name,tf,col in profiles:
    d=results[name]
    t_tip=tf(L_b)*1000
    rows_p.append([name,f'{t_tip:.2f}',f'{(1-t_tip/(t0*1000))*100:.0f}',f'{d["freqs"][0]:.1f}',
        f'{d["freqs"][1]/d["freqs"][0]:.2f}',f'{d["freqs"][2]/d["freqs"][0]:.2f}'])
story.append(make_table(['Profil','t_Spitze [mm]','Verdünnung [%]','f₁ [Hz]','f₂/f₁','f₃/f₁'],rows_p,
    col_widths=[55*mm,22*mm,22*mm,18*mm,18*mm,18*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die Grundfrequenz f₁ steigt mit zunehmender Verdünnung (weniger Masse an der Spitze). '
    'Gleichzeitig sinken die Obertonverhältnisse — die Moden rücken zusammen. '
    'Bei 70 % Verdünnung (Linear 30 %) sinkt f₂/f₁ von 6,27 auf 4,07.',sBody))

# ══ Kap. 3 ══
story.append(Paragraph('3. Obertonverhältnisse im Vergleich',sChapter))
fig2=make_fig2()
story.append(fig_to_image(fig2,width_mm=155))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die grauen horizontalen Linien markieren die harmonische Reihe (1, 2, 3, 4 ...). '
    'Kein Profil kommt auch nur in die Nähe. Selbst bei starker Verdünnung (Linear 30 %) '
    'liegt Mode 2 bei 4,07× statt 2× — mehr als eine Oktave darüber.',sBody))
story.append(Paragraph(
    'Die Inharmonizität ist eine fundamentale Eigenschaft des Cantilevers und '
    'nicht durch Profilierung zu beheben. Der Grund: Die Biegesteifigkeit '
    'erzwingt f ∝ (βL)², nicht f ∝ n. Nur ein Saitenoszillator (Zugspannung '
    'ohne Biegesteifigkeit) hätte harmonische Obertöne.',sBody))

# ══ Kap. 4: Gewichte ══
story.append(Paragraph('4. Stimmgewichte: Der gegenteilige Effekt',sChapter))
story.append(Paragraph(
    'Die Profilierung (Verdünnung zur Spitze) senkt f₂/f₁ — die Moden rücken zusammen. '
    'Stimmgewichte an der Spitze wirken <b>gegenteilig</b>: Sie senken f₁ stark '
    '(mehr Masse am Schwingungsmaximum), aber die höheren Moden weniger '
    '(deren Knoten liegen teilweise nahe der Spitze). Das Ergebnis: '
    'f₂/f₁ <b>steigt</b> mit zunehmendem Gewicht.',sBody))

fig7=make_fig7()
story.append(fig_to_image(fig7,width_mm=155))
story.append(Spacer(1,3*mm))

rows_w=[]
for wn,mt in [('Ohne',0),('0,05 g',0.05e-3),('0,10 g',0.10e-3),('0,20 g',0.20e-3),
    ('0,40 g',0.40e-3),('0,80 g',0.80e-3),('1,50 g',1.50e-3)]:
    fr,_=beam_modes_weighted(N_el,lambda x:t0,L_b,W_b,E_s,rho_s,m_tip=mt,n_modes=4)
    rat=fr/fr[0]
    rows_w.append([wn,f'{mt*1000:.2f}',f'{mt/m_zunge*100:.0f}',f'{fr[0]:.1f}',f'{rat[1]:.2f}',f'{rat[2]:.2f}'])
story.append(make_table(['Gewicht','Masse [g]','% m_Zunge','f₁ [Hz]','f₂/f₁','f₃/f₁'],rows_w,
    col_widths=[25*mm,22*mm,22*mm,20*mm,22*mm,22*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    f'Die Zungenmasse beträgt {m_zunge*1000:.2f} g. Ein Gewicht von 0,20 g (11 %) '
    'senkt f₁ von 66,8 auf 55,2 Hz, erhöht aber f₂/f₁ von 6,27 auf 6,58. '
    'Bei 1,50 g (86 %) steigt f₂/f₁ auf 9,89 — die Obertöne werden noch inharmonischer.',sBody))
story.append(Paragraph(
    '<b>Der Grund:</b> Das Gewicht sitzt am Schwingungsmaximum von Mode 1. '
    'Mode 1 wird maximal gebremst. Mode 2 hat an der Spitze noch deutliche Amplitude, '
    'wird aber weniger gebremst → das Verhältnis steigt.',sKeyBox))

# ══ Kap. 5: Kombination ══
story.append(Paragraph('5. Kombination: Profilierung + Gewicht',sChapter))
story.append(Paragraph(
    'In der Praxis wird beides kombiniert: Die Zunge wird profiliert und '
    'mit einem Stimmgewicht (Lot, Zinn) versehen.',sBody))

fig8=make_fig8()
story.append(fig_to_image(fig8,width_mm=155))
story.append(Spacer(1,3*mm))

rows_c=[]
for cn,tf,mt in [('Gleichmäßig, ohne',lambda x:t0,0),('Gleichmäßig + 0,20 g',lambda x:t0,0.20e-3),
    ('Linear 50 %, ohne',lambda x:t0*(1-0.50*x/L_b),0),('Linear 50 % + 0,20 g',lambda x:t0*(1-0.50*x/L_b),0.20e-3),
    ('Linear 30 %, ohne',lambda x:t0*(1-0.70*x/L_b),0),('Linear 30 % + 0,20 g',lambda x:t0*(1-0.70*x/L_b),0.20e-3)]:
    fr,_=beam_modes_weighted(N_el,tf,L_b,W_b,E_s,rho_s,m_tip=mt,n_modes=3)
    rat=fr/fr[0]
    rows_c.append([cn,f'{fr[0]:.1f}',f'{rat[1]:.2f}',f'{rat[2]:.2f}'])
story.append(make_table(['Konfiguration','f₁ [Hz]','f₂/f₁','f₃/f₁'],rows_c,
    col_widths=[55*mm,22*mm,25*mm,25*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Verdünnung senkt f₂/f₁, Gewicht hebt es wieder an. '
    'Linear 30 % allein: f₂/f₁ = 4,07 bei f₁ = 77,5 Hz. '
    'Linear 30 % + 0,20 g: f₂/f₁ = 4,64 bei f₁ = 50,2 Hz.',sBody))

# ══ Kap. 6: Gleiches Ziel ══
story.append(Paragraph('6. Gleiches f₁-Ziel, verschiedene Wege',sChapter))
fig9=make_fig9()
story.append(fig_to_image(fig9,width_mm=150))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Bei f₁ = 50 Hz (grüne Linie): '
    'Gleichmäßig + Gewicht: f₂/f₁ ≈ 6,7 (am schlechtesten). '
    'Linear 50 % + Gewicht: f₂/f₁ ≈ 5,3 (mittel). '
    'Linear 30 % + Gewicht: f₂/f₁ ≈ 4,6 (am besten).',sBody))
story.append(Paragraph(
    '<b>Und die Ansprache?</b> Die profilierte Zunge hat bei gleichem f₁ '
    'niedrigere Steifigkeit (Dok. 0012) → bessere Ansprache. '
    'Das Gewicht senkt f₁ ohne die Steifigkeit zu ändern '
    '(Masse ≠ Steifigkeit). Deshalb ist Profilierung + wenig Gewicht '
    'besser als keine Profilierung + viel Gewicht.',sKeyBox))

# ══ Kap. 7: Modenformen ══
story.append(Paragraph('7. Modenformen und Knotenpositionen',sChapter))
fig3=make_fig3()
story.append(fig_to_image(fig3,width_mm=155))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die Profilierung verschiebt die Schwingungsknoten zur Spitze hin. '
    'Bei der gleichmäßigen Zunge liegt der erste Knoten von Mode 2 bei x ≈ 55 % der Länge. '
    'Bei Linear 30 % wandert er auf x ≈ 65 %. '
    'Das hat Konsequenzen für die Kopplung: Der Volumenstrom pro Längeneinheit '
    '(Dok. 0015, Kap. 9) wird durch die verschobenen Knoten anders gewichtet.',sBody))
story.append(Paragraph(
    'Auch die Amplitude an der Spitze ändert sich relativ: Bei stark profilierten '
    'Zungen konzentriert sich die Schwingungsenergie stärker am leichten Ende — '
    'die Spitze schwingt relativ weiter aus. Das verbessert den Volumenstrom '
    'durch den Schlitz, aber erhöht auch die Gefahr von Kanalanschlägen '
    '(Dok. 0011, Torsion).',sBody))

# ══ Kap. 5 ══
story.append(Paragraph('8. Inharmonizität in Cent',sChapter))
fig4=make_fig4()
story.append(fig_to_image(fig4,width_mm=155))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Alle Obertöne liegen tausende Cent über der harmonischen Reihe. '
    'Mode 2 liegt bei +1200 bis +2000 Cent (1–2 Oktaven) über 2f₁. '
    'Die Profilierung reduziert die Abweichung um etwa 25–35 % — '
    'von ≈ 2000 Cent (gleichmäßig) auf ≈ 1250 Cent (Linear 30 %). '
    'Selbst bei extremer Verdünnung bleibt die Inharmonizität riesig.',sBody))
story.append(Paragraph(
    '<b>Konsequenz:</b> Die Biegemoden der Zunge fallen nie auf harmonische '
    'Obertöne des Grundtons. Wenn Mode 2 bei 6,27 × f₁ liegt und f₁ = 50 Hz, '
    'dann ist f₂ ≈ 314 Hz — das fällt nicht auf einen harmonischen Oberton '
    '(250, 300, 350 Hz), sondern irgendwo dazwischen. '
    'Die Kammer-Kopplung (Dok. 0009) betrifft daher nicht die harmonischen Obertöne, '
    'sondern die mechanischen Moden — und ob diese in der Nähe von f_H liegen.',sKeyBox))

# ══ Kap. 6 ══
story.append(Paragraph('9. Steifigkeit und Grundfrequenz',sChapter))
fig5=make_fig5()
story.append(fig_to_image(fig5,width_mm=150))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Profilierung hat zwei gleichzeitige Effekte: Die Grundfrequenz steigt '
    '(weniger Masse an der Spitze), und die Steifigkeit sinkt (dünnerer Querschnitt). '
    'Beides zusammen bedeutet: Eine profilierte Zunge für denselben Ton '
    '(z.B. 50 Hz) muss länger oder am Fuß dicker sein — und sie hat bei '
    'gleicher Frequenz eine niedrigere Steifigkeit, also bessere Ansprache '
    '(Dok. 0012).',sBody))
story.append(Paragraph(
    'Das ist der eigentliche Praxisnutzen der Profilierung: Nicht die Änderung '
    'der Obertonverhältnisse (die bleiben weit von harmonisch), sondern die '
    'Kombination aus niedrigerer Steifigkeit (bessere Ansprache) bei '
    'gleicher Grundfrequenz.',sBody))

# ══ Kap. 7 ══
story.append(Paragraph('10. Kammer-Kopplung: Wo fallen die Moden hin?',sChapter))
fig6=make_fig6()
story.append(fig_to_image(fig6,width_mm=155))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die violette Zone markiert den typischen f_H-Bereich einer Basskammer '
    '(200–600 Hz). Mode 1 liegt immer darunter (50–80 Hz). '
    'Mode 2 liegt bei der gleichmäßigen Zunge bei ≈ 420 Hz — mitten im '
    'f_H-Bereich. Profilierung verschiebt Mode 2 nach unten '
    '(Linear 30 %: ≈ 315 Hz) oder lässt ihn stabil '
    '(Parabolisch: ≈ 390 Hz).',sBody))
story.append(Paragraph(
    'Wenn Mode 2 nahe f_H liegt, kann die Kammer die Zunge auf dieser '
    'Mode resonant anregen — das erzeugt „Geistertöne" im Spektrum '
    '(Dok. 0009, Zone 2). Profilierung kann Mode 2 gezielt aus dem '
    'f_H-Bereich herausschieben — das ist ein zweiter Praxisnutzen '
    'neben der Ansprache-Verbesserung.',sBody))

# ══ Kap. 11: Praxis-Profilierung ══
story.append(Paragraph('11. Praxis-Profilierung: Position der dünnsten Stelle',sChapter))
story.append(Paragraph(
    'In der Praxis werden Zungen nicht linear verdünnt. Die Zunge ist am dicksten '
    'bei der Niete (Einspannung) und verjüngt sich glatt über die Länge. '
    'Die dünnste Stelle kann an der Spitze liegen, aber auch irgendwo '
    'zwischen Spitze und Mitte — die Spitze ist dann wieder etwas dicker. '
    'Der Verlauf ist nie abrupt.',sBody))

fig10=make_fig10()
story.append(fig_to_image(fig10,width_mm=155))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    '<b>Überraschendes Ergebnis:</b> Die Position der dünnsten Stelle hat wenig Einfluss '
    'auf f₂/f₁. Zwischen x_min = 0,7L und 1,0L (Spitze) ist das Plateau praktisch flach. '
    'Ob die dünnste Stelle bei 70 % oder an der Spitze liegt — f₂/f₁ ändert sich kaum. '
    'Erst unter x_min = 0,5L (Mitte) wird f₂/f₁ merklich schlechter (steigt gegen 6,27).',sBody))
story.append(Paragraph(
    'Die <b>Tiefe</b> der Verdünnung bestimmt f₂/f₁, nicht die Position. '
    'Die Position bestimmt f₁ — und damit, wie viel Gewicht für den Zielton nötig ist.',sKeyBox))

# ══ Kap. 12: Optimum ══
story.append(Paragraph('12. Optimale Verdünnungstiefe',sChapter))

fig11=make_fig11()
story.append(fig_to_image(fig11,width_mm=145))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'f₂/f₁ hat ein <b>Minimum bei ≈ 70 % Verdünnung</b> (t_min ≈ 0,12 mm bei t₀ = 0,40 mm). '
    'Bei stärkerer Verdünnung steigt f₂/f₁ wieder an — die Zunge wird so dünn, '
    'dass der vordere Teil praktisch keine Steifigkeit mehr hat und die '
    'Modenstruktur sich fundamental ändert.',sBody))

# Tabelle: Verdünnung bei x_min=0.7
rows_depth=[]
for dp_val in [1.0, 0.80, 0.60, 0.50, 0.40, 0.30, 0.20]:
    tf_d=make_smooth_profile(0.70, dp_val)
    fr_d,_=beam_modes_weighted(N_el,tf_d,L_b,W_b,E_s,rho_s,n_modes=4)
    rd=fr_d/fr_d[0]
    rows_depth.append([f'{(1-dp_val)*100:.0f} %',f'{t0*dp_val*1000:.2f}',f'{fr_d[0]:.1f}',f'{rd[1]:.2f}',f'{rd[2]:.2f}'])
story.append(make_table(['Verdünnung','t_min [mm]','f₁ [Hz]','f₂/f₁','f₃/f₁'],rows_depth,
    col_widths=[25*mm,25*mm,22*mm,25*mm,25*mm]))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die Praxis bestätigt dieses Optimum: Erfahrene Stimmer verdünnen '
    'die Zunge bis die Ansprache stimmt — nie darüber hinaus. '
    'Das berechnete Optimum bei 70 % entspricht dem Punkt, an dem '
    'die Verbesserung der Ansprache (Dok. 0012) und die Senkung von f₂/f₁ '
    'gleichzeitig maximal sind.',sBody))

# ══ Kap. 13: Chromatisch ══
story.append(Paragraph('13. Chromatische Analyse: 2 Oktaven, gleiche Mensur',sChapter))
story.append(Paragraph(
    'Bei einer gleichmäßigen Zunge ist f₂/f₁ = 6,27 = <b>konstant</b>, '
    'unabhängig von der Dicke. Da f ∝ t bei konstantem Profil, skalieren '
    'alle Moden gleich mit t — das Verhältnis fällt raus. '
    'Die Mensur (L, Profilform) legt f₂/f₁ für alle Töne fest.',sBody))
story.append(Paragraph(
    'Was sich ändert, ist f₂ <b>absolut</b>:',sBody))

fig12=make_fig12()
story.append(fig_to_image(fig12,width_mm=155))
story.append(Spacer(1,3*mm))

story.append(Paragraph(
    'Die untere Oktave (A1–A2) ist die <b>kritische Zone</b>: f₂ liegt bei '
    '345–689 Hz (gleichmäßig) bzw. 248–495 Hz (profiliert) — '
    'mitten im typischen f_H-Bereich der Basskammern. '
    'Die obere Oktave (A2–A3): f₂ liegt über 700 Hz — weit über f_H.',sBody))
story.append(Paragraph(
    'Das erklärt, warum Kopplung und Geistertöne (Dok. 0009) '
    'hauptsächlich bei den <b>tiefen Tönen</b> auftreten: Dort trifft f₂ '
    'auf f_H. Die Kammer muss für jeden Ton so abgestimmt werden, '
    'dass f_H nicht auf f₂ fällt.',sKeyBox))

# ══ Zusammenfassung ══
story.append(Paragraph('14. Zusammenfassung',sChapter))
pts=[
    '<b>Cantilever-Obertöne sind fundamental inharmonisch:</b> '
    '1 : 6,27 : 17,55 : 34,39 (gleichmäßig). Die hörbaren Harmonischen '
    'stammen vom Impulsgenerator (Luftstrom), nicht von den Biegemoden.',
    '<b>Profilierung komprimiert die Verhältnisse:</b> Linear 30 % senkt '
    'f₂/f₁ von 6,27 auf 4,07. Glatte Praxis-Profilierung: Optimum bei '
    '≈ 70 % Verdünnung (f₂/f₁ ≈ 3,94). Darüber hinaus steigt f₂/f₁ wieder.',
    '<b>Position der dünnsten Stelle:</b> Zwischen 70 % und 100 % der Länge '
    'kaum Unterschied. Die Tiefe der Verdünnung bestimmt f₂/f₁, '
    'die Position bestimmt f₁.',
    '<b>Stimmgewichte spreizen die Verhältnisse:</b> Gegenteiliger Effekt '
    'zur Verdünnung. Profilierung + wenig Gewicht ist besser als '
    'viel Gewicht allein.',
    '<b>f₂/f₁ = const bei gleichmäßigem Profil:</b> Die Mensur legt '
    'das Verhältnis für alle 25 chromatischen Töne fest. '
    'Nur Profilierung oder Gewicht kann es ändern.',
    '<b>Kritische Zone = untere Oktave:</b> f₂ liegt für A1–A2 '
    'bei 345–689 Hz — mitten im f_H-Bereich. Dort treten Geistertöne auf '
    '(Dok. 0009). Die Kammer muss tonnweise abgestimmt werden.',
    '<b>Hauptnutzen der Profilierung:</b> Niedrigere Steifigkeit bei '
    'gleicher Grundfrequenz → bessere Ansprache (Dok. 0012). '
    'Die Änderung der Obertonverhältnisse und die Verschiebung von Mode 2 '
    'aus dem f_H-Bereich sind weitere Vorteile.',
]
for i,p in enumerate(pts):
    story.append(Paragraph(f'{i+1}. {p}',sBody))

story.append(Spacer(1,6*mm))
story.append(Paragraph(
    '<i>Die Zunge schwingt inharmonisch. Der Klang ist harmonisch. '
    'Die Profilierung ändert die Ansprache — die Obertöne macht die Luft.</i>',sAbstract))

def add_page_number(canvas,doc):
    canvas.saveState(); canvas.setFont('DejaVu',8); canvas.setFillColor(HexColor('#999999'))
    canvas.drawCentredString(W_PAGE/2,12*mm,f'Dok. 0016 — Obertonmoden — Seite {canvas.getPageNumber()}')
    canvas.restoreState()

doc.build(story,onFirstPage=add_page_number,onLaterPages=add_page_number)
print(f'✓ {outfile} erzeugt')
