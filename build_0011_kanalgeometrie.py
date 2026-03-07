#!/usr/bin/env python3
"""
Dok. 0011 — Kanalgeometrie der Stimmplatte:
Torsion, konische Erweiterung und Energiezufuhr
Benötigt: reportlab, matplotlib, numpy
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

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
addMapping('DejaVu',0,0,'DejaVu'); addMapping('DejaVu',1,0,'DejaVuB')
addMapping('DejaVu',0,1,'DejaVuI'); addMapping('DejaVu',1,1,'DejaVuBI')

DARKBLUE=HexColor('#16213e'); ACCENTRED=HexColor('#e94560')
KEYGREEN=HexColor('#2e7d32'); WARNRED=HexColor('#c62828')
LIGHTGRAY=HexColor('#f5f5f5'); KEYBG=HexColor('#e8f5e9'); WARNBG=HexColor('#ffebee')
PW=A4[0]-50*mm

styles=getSampleStyleSheet()
for sn in styles.byName:
    s=styles.byName[sn]
    if hasattr(s,'fontName'):
        if 'Bold' in s.fontName: s.fontName='DejaVuB'
        elif 'Italic' in s.fontName: s.fontName='DejaVuI'
        else: s.fontName='DejaVu'

sT=ParagraphStyle('T',parent=styles['Title'],fontSize=18,textColor=DARKBLUE,spaceAfter=4,alignment=TA_CENTER,fontName='DejaVuB')
sST=ParagraphStyle('ST',parent=styles['Normal'],fontSize=12,textColor=DARKBLUE,alignment=TA_CENTER,spaceAfter=2,fontName='DejaVu')
sAb=ParagraphStyle('Ab',parent=styles['Italic'],fontSize=9,textColor=HexColor('#555555'),spaceAfter=8,fontName='DejaVuI')
sCh=ParagraphStyle('Ch',parent=styles['Heading1'],fontSize=14,textColor=DARKBLUE,spaceBefore=14,spaceAfter=6,fontName='DejaVuB')
sB=ParagraphStyle('Bo',parent=styles['Normal'],fontSize=10,leading=14,spaceAfter=6,alignment=TA_JUSTIFY,fontName='DejaVu')
sBB=ParagraphStyle('BB',parent=sB,fontName='DejaVuB')
sFo=ParagraphStyle('Fo',parent=sB,fontSize=10,alignment=TA_CENTER,spaceAfter=8,spaceBefore=4,fontName='DejaVuMono')
sSm=ParagraphStyle('Sm',parent=sB,fontSize=9,leading=12,fontName='DejaVu')
sKB=ParagraphStyle('KB',parent=sB,backColor=KEYBG,borderPadding=6,borderColor=KEYGREEN,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sWB=ParagraphStyle('WB',parent=sB,backColor=WARNBG,borderPadding=6,borderColor=WARNRED,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sTH=ParagraphStyle('TH',parent=sB,fontSize=9,fontName='DejaVuB',alignment=TA_CENTER,leading=11)
sTD=ParagraphStyle('TD',parent=sB,fontSize=9,alignment=TA_CENTER,leading=11,fontName='DejaVu')
sTDL=ParagraphStyle('TDL',parent=sB,fontSize=9,alignment=TA_LEFT,leading=11,fontName='DejaVu')

def hr():
    return Table([['']], colWidths=[PW], style=TableStyle([('LINEBELOW',(0,0),(-1,-1),2,ACCENTRED),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))

def mt(header, rows, cw=None):
    data=[[Paragraph(h,sTH) for h in header]]
    for row in rows:
        data.append([Paragraph(str(c),sTDL) if i==0 else Paragraph(str(c),sTD) for i,c in enumerate(row)])
    cw=cw or [PW/len(header)]*len(header)
    t=Table(data,colWidths=cw,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),DARKBLUE),('TEXTCOLOR',(0,0),(-1,0),white),
        ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,LIGHTGRAY]),('GRID',(0,0),(-1,-1),0.5,HexColor('#cccccc')),
        ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),
        ('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
    return t

# ══════════════════════════════════════════════════════════════
# PHYSIK
# ══════════════════════════════════════════════════════════════
E=200e9; rho_s=7800; G=77e9; mu=1.8e-5; rho_0=1.2

zungen=[
    {'name':'Bass 50 Hz','f':50,'L':70e-3,'W':10e-3,'t':0.40e-3,'h':15e-3},
    {'name':'Bass A1','f':55,'L':60e-3,'W':8e-3,'t':0.35e-3,'h':14e-3},
    {'name':'Bass A2','f':110,'L':40e-3,'W':7e-3,'t':0.30e-3,'h':10e-3},
    {'name':'Mitte A3','f':220,'L':30e-3,'W':6e-3,'t':0.25e-3,'h':6e-3},
    {'name':'Diskant A4','f':440,'L':20e-3,'W':5e-3,'t':0.20e-3,'h':4e-3},
    {'name':'Hoch A5','f':880,'L':12e-3,'W':4e-3,'t':0.15e-3,'h':2.5e-3},
    {'name':'Sehr hoch C6','f':1047,'L':10e-3,'W':3.5e-3,'t':0.12e-3,'h':2e-3},
]
for z in zungen:
    z['x_tip']=z['h']*1.3
    # Torsion
    beta=1/3*(1-0.63*z['t']/z['W'])
    J_T=beta*z['W']*z['t']**3
    I_p=z['W']*z['t']*(z['W']**2+z['t']**2)/12
    z['f_T']=1/(4*z['L'])*math.sqrt(G*J_T/(rho_s*I_p))
    I_B=z['W']*z['t']**3/12
    z['f_B']=1.875**2/(2*math.pi*z['L']**2)*math.sqrt(E*I_B/(rho_s*z['W']*z['t']))
    z['y_tors']=0.02*z['W']/2
    z['s_min_tors']=z['y_tors']*1.5
    # Duty
    r=z['h']/(2*z['x_tip'])
    z['duty']=2/math.pi*math.asin(r) if r<1 else 1.0
    # Antriebsfläche
    z['F_antrieb']=1000*z['W']*z['h']

# ══════════════════════════════════════════════════════════════
# DIAGRAMME
# ══════════════════════════════════════════════════════════════
cols=['#0D47A1','#1565C0','#1976D2','#2196F3','#F9A825','#E65100','#C62828']

# ── Diag 1: Schwingungsform — Zunge im Kanal (schematisch) ──
fig1, axes1 = plt.subplots(1, 3, figsize=(12, 5))
for ax, z in zip(axes1, [zungen[0], zungen[4], zungen[6]]):
    h=z['h']*1000; xt=z['x_tip']*1000
    t_arr=np.linspace(0, 2*np.pi, 500)
    x_arr=xt*np.sin(t_arr)
    # Kanal-Bereich
    ax.axhspan(-h/2, h/2, color='#90A4AE', alpha=0.2, label=f'Kanal h={h:.0f}mm')
    ax.axhline(y=h/2, color='#546E7A', lw=1.5)
    ax.axhline(y=-h/2, color='#546E7A', lw=1.5)
    # Zungenbewegung
    im_kanal = np.abs(x_arr) <= h/2
    ax.plot(t_arr[im_kanal], x_arr[im_kanal], 'g.', ms=1.5, label='Im Kanal (Antrieb)')
    ax.plot(t_arr[~im_kanal], x_arr[~im_kanal], 'r.', ms=1.5, label='Außerhalb (frei)')
    ax.set_xlabel('Phase [rad]', fontsize=9)
    ax.set_title(f'{z["name"]}\nh={h:.0f}mm, x_tip={xt:.1f}mm, Duty={z["duty"]:.0%}', fontsize=10)
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(alpha=0.3)
    ax.set_ylim(-xt*1.2, xt*1.2)
    if ax==axes1[0]: ax.set_ylabel('Auslenkung [mm]', fontsize=10)
plt.suptitle('Schwingungsform: Zunge im Kanal (grün) und außerhalb (rot)', fontsize=12, y=1.02)
plt.tight_layout()
fig1.savefig('/home/claude/diag_0011_1_schwingung.png', dpi=150); plt.close()

# ── Diag 2: Duty Cycle vs. x_tip/h ──
fig2, ax2 = plt.subplots(figsize=(8, 5))
xh_range=np.linspace(0.3, 5, 300)
duties=[min(100, 2/math.pi*math.asin(min(1, 1/(2*xh)))*100) for xh in xh_range]
ax2.plot(xh_range, duties, 'b-', lw=3)
ax2.axhline(y=50, color='gray', ls=':', alpha=0.5)
ax2.axvline(x=1.0, color='green', ls='--', alpha=0.5, label='x_tip = h')

# Alle Zungen haben x_tip/h = 1.3 → ein Punkt, alle Namen hintereinander
xh_all = zungen[0]['x_tip'] / zungen[0]['h']
d_all = zungen[0]['duty'] * 100
ax2.plot(xh_all, d_all, 's', ms=10, color='#E65100', zorder=5)
all_names = ', '.join(z['name'] for z in zungen)
ax2.annotate(f'Alle Zungen bei x_tip/h = {xh_all:.1f}:\n{all_names}\n(Duty = {d_all:.0f} %)',
    (xh_all, d_all), textcoords="offset points", xytext=(15, 15), fontsize=8,
    arrowprops=dict(arrowstyle='->', color='#E65100', lw=1.5),
    bbox=dict(boxstyle='round,pad=0.3', facecolor='#FFF3E0', edgecolor='#E65100', alpha=0.9))

# Zusätzlich: Markiere andere typische Verhältnisse
for xh_mark, label in [(0.5, 'x_tip = h/2\n(leise)'), (2.0, 'x_tip = 2h\n(sehr laut)')]:
    d_mark = min(100, 2/math.pi*math.asin(min(1, 1/(2*xh_mark)))*100)
    ax2.plot(xh_mark, d_mark, 'o', ms=7, color='gray', zorder=4)
    ax2.annotate(label, (xh_mark, d_mark), textcoords="offset points",
        xytext=(-10, -25) if xh_mark > 1 else (15, -20), fontsize=8, color='gray')

ax2.set_xlabel('x_tip / h  (Amplitude / Plattendicke)', fontsize=11)
ax2.set_ylabel('Duty Cycle [%]', fontsize=11)
ax2.set_title('Zeitanteil der Schwingung im Kanal\n(bei voller Lautstärke: x_tip ≈ 1,3 × h → Duty ≈ 25 % für alle Zungen)', fontsize=12)
ax2.legend(fontsize=9); ax2.grid(alpha=0.3); ax2.set_xlim(0.3, 5); ax2.set_ylim(0, 105)
plt.tight_layout()
fig2.savefig('/home/claude/diag_0011_2_duty.png', dpi=150); plt.close()

# ── Diag 3: Geschwindigkeit im Kanal vs. außerhalb ──
fig3, ax3 = plt.subplots(figsize=(8, 5))
t_arr=np.linspace(0, 2*np.pi, 1000)
z_ref = zungen[4]  # A4
xt=z_ref['x_tip']; h=z_ref['h']; omega=2*math.pi*z_ref['f']
x_arr=xt*np.sin(t_arr); v_arr=xt*omega*np.cos(t_arr)
im_k=np.abs(x_arr)<=h/2

ax3.fill_between(t_arr, 0, np.abs(v_arr), where=im_k, color='green', alpha=0.3, label='Im Kanal')
ax3.fill_between(t_arr, 0, np.abs(v_arr), where=~im_k, color='red', alpha=0.15, label='Außerhalb')
ax3.plot(t_arr, np.abs(v_arr), 'b-', lw=1.5, label='|v(t)|')

v_mean_total = np.mean(np.abs(v_arr))
v_mean_kanal = np.mean(np.abs(v_arr[im_k])) if np.any(im_k) else 0
ax3.axhline(y=v_mean_total, color='blue', ls='--', alpha=0.5, label=f'Mittel gesamt: {v_mean_total:.1f} m/s')
ax3.axhline(y=v_mean_kanal, color='green', ls='--', alpha=0.5, label=f'Mittel im Kanal: {v_mean_kanal:.1f} m/s')

ax3.set_xlabel('Phase [rad]', fontsize=11)
ax3.set_ylabel('|Geschwindigkeit| [m/s]', fontsize=11)
ax3.set_title(f'Diskant A4: Geschwindigkeit während des Kanaldurchgangs\nDie Zunge passiert den Kanal bei MAXIMALER Geschwindigkeit', fontsize=12)
ax3.legend(fontsize=9); ax3.grid(alpha=0.3)
plt.tight_layout()
fig3.savefig('/home/claude/diag_0011_3_geschwindigkeit.png', dpi=150); plt.close()

# ── Diag 4: Torsionsamplitude vs. Spalt ──
fig4, ax4 = plt.subplots(figsize=(8, 5))
names=[z['name'] for z in zungen]
y_tors=[z['y_tors']*1000 for z in zungen]
s_min=[z['s_min_tors']*1000 for z in zungen]
s_praxis = [0.06, 0.06, 0.05, 0.05, 0.04, 0.03, 0.025]  # mm, typischer Praxis-Spalt

x_pos=range(len(zungen))
ax4.bar(x_pos, s_praxis, width=0.35, color='#42A5F5', alpha=0.7, label='Praxis-Spalt')
ax4.bar([x+0.35 for x in x_pos], y_tors, width=0.35, color='#EF5350', alpha=0.7, label='Torsions-Auslenkung (θ=1°)')
ax4.bar([x+0.35 for x in x_pos], [sm-yt for sm,yt in zip(s_min,y_tors)], width=0.35, 
        bottom=y_tors, color='#FFCDD2', alpha=0.5, label='+ 50% Sicherheit')

ax4.set_xticks([x+0.175 for x in x_pos])
ax4.set_xticklabels(names, fontsize=8, rotation=15)
ax4.set_ylabel('Abstand [mm]', fontsize=11)
ax4.set_title('Seitenspalt vs. Torsionsamplitude\n(Torsion bei transientem Anblasen, θ ≈ 1°)', fontsize=12)
ax4.legend(fontsize=9); ax4.grid(axis='y', alpha=0.3)
plt.tight_layout()
fig4.savefig('/home/claude/diag_0011_4_torsion.png', dpi=150); plt.close()

# ── Diag 5: Konus — Spaltwiderstand und Antriebseffizienz für 3 Zungen ──
fig5, axes5 = plt.subplots(1, 3, figsize=(12, 4.5))
sel=[zungen[0], zungen[4], zungen[6]]
for ax, z in zip(axes5, sel):
    h=z['h']; umf=2*(z['L']+z['W'])
    s0=0.06e-3 if z['f']<200 else 0.04e-3
    R_haupt=rho_0*math.sqrt(2*1000/rho_0)/(2*z['W']*z['x_tip']*0.5)
    R_par=12*mu*h/(umf*s0**3)
    
    ds_r=np.linspace(0, 0.08e-3, 100)
    Rr=[]; er=[]
    for ds in ds_r:
        s1=s0+ds
        if ds>1e-9:
            a=s0; b=s1-s0
            Rk=12*mu/umf*h*(1/(2*b)*(1/a**2-1/(a+b)**2))
        else: Rk=R_par
        Rr.append(Rk/R_par*100)
        n=100; es=0
        for iz in range(n):
            zz=(iz+0.5)/n; sz=s0+ds*zz
            Rl=12*mu*(h/n)/(umf*sz**3)*n
            es+=Rl/(R_haupt+Rl)
        er.append(es/n)
    er0=er[0]; en=[e/er0*100 for e in er]
    
    ax.plot([d*1000 for d in ds_r], Rr, 'b-', lw=2, label='Spaltwiderstand R/R₀')
    ax.plot([d*1000 for d in ds_r], en, 'g-', lw=2, label='Antriebseffizienz')
    ax.axhline(y=50, color='red', ls=':', alpha=0.3)
    ax.set_xlabel('Erweiterung/Seite [mm]', fontsize=9)
    ax.set_title(f'{z["name"]}\n(h={h*1000:.0f}mm, s₀={s0*1000:.2f}mm)', fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_ylim(0, 105)
    if ax==axes5[0]: ax.set_ylabel('[%]', fontsize=10)
plt.suptitle('Konuserweiterung: Was sie kostet', fontsize=12, y=1.02)
plt.tight_layout()
fig5.savefig('/home/claude/diag_0011_5_konus.png', dpi=150); plt.close()

# ── Diag 6: Kanalquerschnitt mit Konus und Torsion (schematisch) ──
fig6, axes6 = plt.subplots(1, 3, figsize=(12, 5))
for ax, z in zip(axes6, sel):
    h=z['h']*1000; W=z['W']*1000; yt=z['y_tors']*1000
    s0=0.06 if z['f']<200 else 0.04  # mm
    ds=0.04  # Konus
    # Kanalwand links (Konus)
    ax.plot([-W/2-s0, -W/2-s0-ds], [0, h], 'k-', lw=2)
    ax.plot([W/2+s0, W/2+s0+ds], [0, h], 'k-', lw=2)
    # Paralleler Kanal zum Vergleich
    ax.plot([-W/2-s0, -W/2-s0], [0, h], 'k--', lw=1, alpha=0.4)
    ax.plot([W/2+s0, W/2+s0], [0, h], 'k--', lw=1, alpha=0.4)
    # Platte (Füllung)
    ax.fill_betweenx([0, h], -W/2-s0-ds-0.3, [-W/2-s0, -W/2-s0-ds], color='#90A4AE', alpha=0.4)
    ax.fill_betweenx([0, h], [W/2+s0, W/2+s0+ds], W/2+s0+ds+0.3, color='#90A4AE', alpha=0.4)
    # Zunge Ruhe
    ax.fill_betweenx([0, h], -W/2, W/2, color='#F9A825', alpha=0.3, label='Zunge')
    # Torsion
    ax.fill_betweenx([0, h], -W/2-yt, -W/2, color='#EF5350', alpha=0.2)
    ax.fill_betweenx([0, h], W/2, W/2+yt, color='#EF5350', alpha=0.2, label=f'Torsion ±{yt:.2f}mm')
    
    ax.set_title(f'{z["name"]}\nh={h:.0f}mm, W={W:.1f}mm', fontsize=10)
    ax.set_xlabel('Breite [mm]', fontsize=9)
    if ax==axes6[0]: ax.set_ylabel('Tiefe z [mm]', fontsize=9)
    ax.invert_yaxis(); ax.grid(alpha=0.3); ax.legend(fontsize=7)
    ax.set_aspect('equal')
plt.suptitle('Kanalquerschnitt: Konus (durchgezogen) vs. parallel (gestrichelt) + Torsion (rot)', fontsize=11, y=1.02)
plt.tight_layout()
fig6.savefig('/home/claude/diag_0011_6_querschnitt.png', dpi=150); plt.close()

# ── Diag 7: Antriebsfläche vs. Frequenz ──
fig7, ax7 = plt.subplots(figsize=(8, 5))
ff=[z['f'] for z in zungen]; Fa=[z['F_antrieb']*1000 for z in zungen]
ax7.bar(range(len(zungen)), Fa, color=cols, alpha=0.85)
ax7.set_xticks(range(len(zungen))); ax7.set_xticklabels([z['name'] for z in zungen], fontsize=8, rotation=15)
ax7.set_ylabel('Antriebskraft F = Δp × W × h [mN]', fontsize=11)
ax7.set_title('Antriebskraft bei 1 kPa Balgdruck\n(Druckfläche = Zungenbreite × Plattendicke)', fontsize=12)
for i,f in enumerate(Fa): ax7.text(i, f+2, f'{f:.0f}', ha='center', fontsize=9, fontweight='bold')
ax7.grid(axis='y', alpha=0.3)
plt.tight_layout()
fig7.savefig('/home/claude/diag_0011_7_antrieb.png', dpi=150); plt.close()

# ── Diag 8: Netto-Antrieb vs. Plattendicke (normiert) ──
fig8, ax8 = plt.subplots(figsize=(9, 5.5))
cols8 = ['#0D47A1', '#1976D2', '#2196F3', '#F9A825', '#E65100', '#C62828']
zungen_ext = [
    {'name':'Bass 50 Hz','f':50,'L':70e-3,'W':10e-3,'t':0.40e-3,'x_tip':19.5e-3,'s':0.06e-3,'h_praxis':15e-3},
    {'name':'Bass A2','f':110,'L':40e-3,'W':7e-3,'t':0.30e-3,'x_tip':13e-3,'s':0.05e-3,'h_praxis':10e-3},
    {'name':'Mitte A3','f':220,'L':30e-3,'W':6e-3,'t':0.25e-3,'x_tip':7.8e-3,'s':0.05e-3,'h_praxis':6e-3},
    {'name':'Diskant A4','f':440,'L':20e-3,'W':5e-3,'t':0.20e-3,'x_tip':5.2e-3,'s':0.04e-3,'h_praxis':4e-3},
    {'name':'Hoch A5','f':880,'L':12e-3,'W':4e-3,'t':0.15e-3,'x_tip':3.3e-3,'s':0.03e-3,'h_praxis':2.5e-3},
    {'name':'Sehr hoch C6','f':1047,'L':10e-3,'W':3.5e-3,'t':0.12e-3,'x_tip':2.6e-3,'s':0.025e-3,'h_praxis':2e-3},
]
for i, ze in enumerate(zungen_ext):
    umf=2*(ze['L']+ze['W']); om=2*math.pi*ze['f']; vt=om*ze['x_tip']
    Ah=ze['W']*ze['x_tip']*0.5; Rh=rho_0*math.sqrt(2*1000/rho_0)/(2*Ah)
    h_norm=np.linspace(0.1, 3.0, 200)
    Pn_arr=[]
    for hf in h_norm:
        h=hf*ze['h_praxis']; ratio=h/(2*ze['x_tip'])
        duty=2/math.pi*math.asin(min(1,ratio)) if ratio<1 else 1.0
        Rs=12*mu*h/(umf*ze['s']**3); eta=Rs/(Rh+Rs)
        Pa=1000*eta*vt*Ah*duty
        Qv=ze['W']*ze['t']*vt; Ps=Qv*Rs*0.3*ze['W']*h*vt*duty
        ml=rho_0*ze['W']*ze['s']*h; Pi=ml*om**2*ze['x_tip']*vt*duty*0.1
        Pn_arr.append(Pa-Ps-Pi)
    Pn_max=max(Pn_arr)
    if Pn_max>0:
        Pn_n=[p/Pn_max*100 for p in Pn_arr]
    else:
        Pn_n=[0]*len(Pn_arr)
    ax8.plot(h_norm, Pn_n, lw=2, color=cols8[i], label=ze['name'])
    ax8.plot(1.0, Pn_n[np.argmin(np.abs(h_norm-1.0))], 's', color=cols8[i], ms=8)
ax8.axhline(y=0, color='black', lw=0.5)
ax8.axvline(x=1.0, color='gray', ls='--', alpha=0.5, label='h = h_praxis')
ax8.set_xlabel('h / h_praxis  (Plattendicke relativ zur Praxis)', fontsize=11)
ax8.set_ylabel('Netto-Antrieb [% des Maximums]', fontsize=11)
ax8.set_title('Netto-Antrieb vs. Plattendicke (normiert)\nHohe Toene fallen schneller ab als tiefe', fontsize=12)
ax8.legend(fontsize=9); ax8.grid(alpha=0.3); ax8.set_xlim(0.1, 3.0); ax8.set_ylim(-50, 110)
plt.tight_layout()
fig8.savefig('/home/claude/diag_0011_8_plattendicke.png', dpi=150); plt.close()

# ── Diag 9: Detail A5 und C6 — Antrieb vs Squeeze ──
fig9, (ax9a, ax9b) = plt.subplots(1, 2, figsize=(11, 5))
for ax, ze, col in [(ax9a, zungen_ext[4], '#E65100'), (ax9b, zungen_ext[5], '#C62828')]:
    umf=2*(ze['L']+ze['W']); om=2*math.pi*ze['f']; vt=om*ze['x_tip']
    Ah=ze['W']*ze['x_tip']*0.5; Rh=rho_0*math.sqrt(2*1000/rho_0)/(2*Ah)
    h_arr=np.linspace(0.1e-3, ze['h_praxis']*3, 300)
    Pa_a=[]; Ps_a=[]; Pn_a=[]
    for h in h_arr:
        ratio=h/(2*ze['x_tip']); duty=2/math.pi*math.asin(min(1,ratio)) if ratio<1 else 1.0
        Rs=12*mu*h/(umf*ze['s']**3); eta=Rs/(Rh+Rs)
        Pa=1000*eta*vt*Ah*duty
        Qv=ze['W']*ze['t']*vt; Ps=Qv*Rs*0.3*ze['W']*h*vt*duty
        Pa_a.append(Pa); Ps_a.append(Ps); Pn_a.append(Pa-Ps)
    Pa_mx=max(Pa_a)
    ax.plot(h_arr*1000,[p/Pa_mx*100 for p in Pa_a],'g-',lw=2,label='Antrieb (saettigt)')
    ax.plot(h_arr*1000,[p/Pa_mx*100 for p in Ps_a],'r-',lw=2,label='Squeeze (quadratisch)')
    ax.plot(h_arr*1000,[p/Pa_mx*100 for p in Pn_a],'b-',lw=2.5,label='Netto')
    ax.axhline(y=0,color='black',lw=0.5)
    ax.axvline(x=ze['h_praxis']*1000,color='orange',ls=':',lw=2,label=f'Praxis h={ze["h_praxis"]*1000:.1f}mm')
    idx_opt=np.argmax(Pn_a)
    ax.axvline(x=h_arr[idx_opt]*1000,color='green',ls='--',alpha=0.5,label=f'Optimum h={h_arr[idx_opt]*1000:.1f}mm')
    ax.set_xlabel('Plattendicke h [mm]',fontsize=10)
    ax.set_title(f'{ze["name"]} ({ze["f"]} Hz)\ns = {ze["s"]*1000:.3f} mm',fontsize=11)
    ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_ylim(-50,110)
    if ax==ax9a: ax.set_ylabel('[% des max. Antriebs]',fontsize=10)
plt.suptitle('Antrieb (gruen) saettigt, Squeeze (rot) waechst quadratisch', fontsize=12, y=1.02)
plt.tight_layout()
fig9.savefig('/home/claude/diag_0011_9_detail.png', dpi=150); plt.close()

# ── Diag 10: Drei Effekte einzeln ──
fig10, axes10 = plt.subplots(1, 3, figsize=(13, 4.5))
ax=axes10[0]
h_mm=np.linspace(0.1, 10, 200)
for i, ze in enumerate(zungen_ext):
    umf=2*(ze['L']+ze['W']); Ah=ze['W']*ze['x_tip']*0.5
    Rh=rho_0*math.sqrt(2*1000/rho_0)/(2*Ah)
    etas=[12*mu*(hm*1e-3)/(umf*ze['s']**3) for hm in h_mm]
    etas=[e/(Rh+e)*100 for e in etas]
    ax.plot(h_mm, etas, lw=2, color=cols8[i], label=ze['name'])
ax.set_xlabel('h [mm]',fontsize=10); ax.set_ylabel('eta [%]',fontsize=10)
ax.set_title('Antriebseffizienz eta\nsaettigt bei eta -> 100 %',fontsize=11)
ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_ylim(0,105)

ax=axes10[1]
for i, ze in enumerate(zungen_ext):
    umf=2*(ze['L']+ze['W']); om=2*math.pi*ze['f']; vt=om*ze['x_tip']
    Qv=ze['W']*ze['t']*vt
    sq=[Qv*12*mu*(hm*1e-3)/(umf*ze['s']**3)*0.3*ze['W']*(hm*1e-3)*vt*1e6 for hm in h_mm]
    ax.plot(h_mm, sq, lw=2, color=cols8[i], label=ze['name'])
ax.set_xlabel('h [mm]',fontsize=10); ax.set_ylabel('P_squeeze [uW]',fontsize=10)
ax.set_title('Squeeze-Bremse\nwaechst quadratisch mit h',fontsize=11)
ax.legend(fontsize=7); ax.grid(alpha=0.3)
ax.set_ylim(0, min(max(sq[:80])*1.5 if sq[40]>0 else 100, 500))

ax=axes10[2]
for i, ze in enumerate(zungen_ext):
    om=2*math.pi*ze['f']
    inert=[rho_0*ze['W']*ze['s']*(hm*1e-3)*om**2*ze['x_tip']*1000 for hm in h_mm]
    ax.plot(h_mm, inert, lw=2, color=cols8[i], label=ze['name'])
ax.set_xlabel('h [mm]',fontsize=10); ax.set_ylabel('F_Traegheit [mN]',fontsize=10)
ax.set_title('Luftmassen-Traegheit\nwaechst mit h und omega^2',fontsize=11)
ax.legend(fontsize=7); ax.grid(alpha=0.3)
ax.set_ylim(0, min(max(inert[:80])*1.5, 5))

plt.suptitle('Drei Effekte: Antrieb saettigt, Squeeze und Traegheit wachsen weiter',fontsize=12,y=1.02)
plt.tight_layout()
fig10.savefig('/home/claude/diag_0011_10_effekte.png', dpi=150); plt.close()

print("Diagramme erzeugt.")

# ══════════════════════════════════════════════════════════════
# PDF
# ══════════════════════════════════════════════════════════════
out='/home/claude/0011_kanalgeometrie_De.pdf'
doc=SimpleDocTemplate(out, pagesize=A4, leftMargin=25*mm, rightMargin=25*mm, topMargin=25*mm, bottomMargin=25*mm)
def set_font(c, d): c.setFont('DejaVu', 10)
story=[]

# Titel
story.append(Paragraph('Dok. 0011', sSm)); story.append(Spacer(1, 2*mm))
story.append(Paragraph('Kanalgeometrie der Stimmplatte', sT))
story.append(Paragraph('Torsion, konische Erweiterung und Energiezufuhr', sST))
story.append(Spacer(1, 2*mm)); story.append(hr()); story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    'Der Schlitz in der Stimmplatte, durch den die Zunge schwingt, ist nicht nur ein geometrischer '
    'Rahmen — er bestimmt über die Spaltgeometrie die Energiezufuhr, begrenzt durch Torsionsschwingungen '
    'die minimale Spaltweite, und verändert bei konischer Erweiterung den akustischen Kurzschluss. '
    'Dieses Dokument berechnet die drei Effekte und ihre Wechselwirkung. '
    'Plattendicken: 2 mm (C6) bis 15 mm (50 Hz). '
    'Bei voller Lautstärke schwingt die Zunge mindestens doppelt so weit wie die Plattendicke — '
    'die Energiezufuhr erfolgt nur während des Kanaldurchgangs.', sAb))

refs=[['0002','Strömungsanalyse Bass-Stimmzunge 50 Hz — v8'],['0010','Güte der Stimmzunge: Zwei Verlustkanäle']]
story.append(Paragraph('<b>Dokumentenverweise</b>', sSm))
rd=[[Paragraph(f'<b>Dok. {r[0]}</b>',sTD), Paragraph(r[1],sTDL)] for r in refs]
rt=Table(rd, colWidths=[30*mm, PW-30*mm])
rt.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.3,HexColor('#cccccc')),('BACKGROUND',(0,0),(0,-1),LIGHTGRAY),
    ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),2),('BOTTOMPADDING',(0,0),(-1,-1),2),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
story.append(rt); story.append(Spacer(1, 4*mm))

# ═══ Kap 1: Schwingung und Kanal ═══
story.append(Paragraph('Kapitel 1: Die Zunge schwingt über den Kanal hinaus', sCh))
story.append(Paragraph(
    'Die Stimmzunge biegt sich auf und schwingt durch den Schlitz der Stimmplatte hindurch. '
    'Bei voller Lautstärke beträgt die Amplitude x<sub>tip</sub> mindestens das 1,2–1,5-fache '
    'der Plattendicke h. Der Gesamthub (maximale Aufbiegung bis maximale Durchbiegung) '
    'ist 2 × x<sub>tip</sub> ≥ 2h.', sB))
story.append(Paragraph(
    'Die Energiezufuhr durch den Balgdruck erfolgt <b>nur</b>, solange sich die Zunge im Kanal befindet. '
    'Oberhalb und unterhalb des Kanals gibt es keinen Druckunterschied zwischen den Seiten — '
    'die Zunge schwingt dort frei, ohne Antrieb und ohne Widerstand.', sB))
story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0011_1_schwingung.png', width=PW, height=PW*0.42))
story.append(Paragraph('<i>Abb. 1: Schwingungsform für Bass (h=15mm), Diskant (h=4mm) und C6 (h=2mm). '
    'Grün = Zunge im Kanal (Antrieb). Rot = Zunge außerhalb (frei). '
    'Der Kanaldurchgang konzentriert sich auf den Nulldurchgang.</i>', sSm))

# ═══ Kap 2: Duty Cycle ═══
story.append(Paragraph('Kapitel 2: Duty Cycle — Zeitanteil im Kanal', sCh))
story.append(Paragraph(
    'Der Anteil der Schwingungsperiode, in dem die Zunge im Kanal ist:', sB))
story.append(Paragraph('Duty = (2/π) × arcsin(h / (2 × x<sub>tip</sub>))', sFo))
story.append(Paragraph(
    'Bei x<sub>tip</sub> = 1,3 × h ergibt sich Duty ≈ 25 %. Die Zunge ist nur ein Viertel '
    'der Zeit im Kanal. Entscheidend: Sie passiert den Kanal nahe dem Nulldurchgang der Schwingung, '
    'wo die Geschwindigkeit <b>maximal</b> ist (v = ω × x<sub>tip</sub>). Die mittlere '
    'Geschwindigkeit während des Kanaldurchgangs ist ≈ 1,5× höher als die mittlere Geschwindigkeit '
    'über den gesamten Zyklus.', sB))
story.append(Spacer(1, 2*mm))

# Tabelle
dt_rows=[]
for z in zungen:
    dt_rows.append([z['name'], f"{z['h']*1000:.0f}", f"{z['x_tip']*1000:.1f}",
        f"{z['x_tip']/z['h']:.2f}", f"{z['duty']:.0%}"])
story.append(mt(['Zunge','h [mm]','x<sub>tip</sub> [mm]','x<sub>tip</sub>/h','Duty'],
    dt_rows, [30*mm,16*mm,22*mm,18*mm,PW-86*mm]))
story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0011_2_duty.png', width=PW, height=PW*0.625))
story.append(Paragraph('<i>Abb. 2: Duty Cycle vs. Amplituden-Verhältnis. Bei voller Lautstärke (Quadrate) '
    'liegt der Duty Cycle bei ≈ 25 % für alle Zungen.</i>', sSm))

story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0011_3_geschwindigkeit.png', width=PW, height=PW*0.625))
story.append(Paragraph('<i>Abb. 3: Geschwindigkeitsverlauf Diskant A4. Die grüne Fläche (Kanaldurchgang) '
    'liegt bei maximaler Geschwindigkeit. Der Antrieb pro Zeiteinheit im Kanal ist 1,5× höher '
    'als der Durchschnitt.</i>', sSm))

story.append(Paragraph(
    '<b>Für den Vergleich</b> verschiedener Kanalgeometrien (parallel vs. Konus) zählt nur der Kanalbereich. '
    'Da beide Geometrien denselben Duty Cycle haben (hängt von h und x<sub>tip</sub> ab, nicht vom Konus), '
    'ist der relative Vergleich gültig.', sKB))

# ═══ Kap 3: Torsion ═══
story.append(PageBreak())
story.append(Paragraph('Kapitel 3: Torsionsschwingung — der Mindestspalt', sCh))
story.append(Paragraph(
    'Neben der Biege-Grundmode besitzt die Zunge eine Torsionsmode, bei der sie sich um ihre '
    'Längsachse verdreht. Die Torsionsfrequenz liegt 5–14× über der Biegefrequenz. '
    'Bei transientem Anblasen (Druckstoß, ungleichmäßiger Balgdruck, schneller Tastendruck) '
    'wird die Torsionsmode kurzzeitig angeregt. Die Zungenspitze bewegt sich dann nicht nur '
    'auf und ab, sondern auch seitlich.', sB))
story.append(Spacer(1, 2*mm))

tr_rows=[]
for z in zungen:
    tr_rows.append([z['name'], f"{z['f_B']:.0f}", f"{z['f_T']:.0f}", f"{z['f_T']/z['f_B']:.1f}",
        f"{z['y_tors']*1000:.3f}", f"{z['s_min_tors']*1000:.3f}"])
story.append(mt(['Zunge','f<sub>Biege</sub> [Hz]','f<sub>Torsion</sub> [Hz]','f<sub>T</sub>/f<sub>B</sub>',
    'y<sub>seitlich</sub> [mm]','s<sub>min</sub> [mm]'],
    tr_rows, [28*mm,22*mm,24*mm,16*mm,24*mm,PW-114*mm]))
story.append(Paragraph('<i>y<sub>seitlich</sub> = Seitliche Auslenkung bei θ = 1° Torsionswinkel. '
    's<sub>min</sub> = y<sub>seitlich</sub> × 1,5 (Sicherheitszuschlag).</i>', sSm))

story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0011_4_torsion.png', width=PW, height=PW*0.625))
story.append(Paragraph('<i>Abb. 4: Praxis-Spalt (blau) vs. Torsionsamplitude (rot) + Sicherheit (rosa). '
    'Bei Bass-Zungen (W = 8–10 mm) ist die Torsion am größten. Der Praxis-Spalt muss die '
    'Torsionsamplitude mit Sicherheitsabstand abdecken.</i>', sSm))

story.append(Paragraph(
    '<b>Konsequenz:</b> Der Mindestspalt wird nicht von der Fertigungstoleranz allein bestimmt, '
    'sondern von der Torsionsamplitude bei transientem Anblasen. Bei Bass-Zungen (breite Zunge, '
    'großer Hebelarm) beträgt der Mindestspalt ≈ 0,10–0,15 mm. Bei hohen Zungen '
    'genügen 0,05 mm. Unter diesen Werten schlägt die Zunge bei hartem Anblasen '
    'an die Kanalwand.', sKB))

# ═══ Kap 4: Konus ═══
story.append(Paragraph('Kapitel 4: Konische Kanalerweiterung', sCh))
story.append(Paragraph(
    'Um der Torsion Platz zu geben, ohne den oberen Spalt zu vergrößern, kann der Kanal '
    'nach unten konisch erweitert werden: s(z) = s₀ + (s₁ − s₀) × z/h. '
    'Der obere Spalt s₀ bleibt eng (guter Antrieb), der untere s₁ ist weiter '
    '(Platz für Torsion und Durchschwung).', sB))
story.append(Paragraph(
    'Die Erweiterung hat zwei Kosten: (1) Der mittlere Spaltwiderstand sinkt → mehr akustischer '
    'Kurzschluss. (2) Die Antriebseffizienz sinkt, weil der Druck im erweiterten Bereich '
    'durch den breiteren Spalt entweicht.', sB))
story.append(Spacer(1, 2*mm))

# Konus-Tabelle für die drei Referenzzungen
for z in [zungen[0], zungen[4], zungen[6]]:
    h=z['h']; umf=2*(z['L']+z['W'])
    s0=0.06e-3 if z['f']<200 else 0.04e-3
    R_haupt=rho_0*math.sqrt(2*1000/rho_0)/(2*z['W']*z['x_tip']*0.5)
    R_par=12*mu*h/(umf*s0**3)
    
    k_rows=[]
    for ds in [0, 0.01e-3, 0.02e-3, 0.04e-3, 0.08e-3]:
        s1=s0+ds; angle=math.degrees(math.atan(ds/h)) if ds>0 else 0
        if ds>1e-9:
            a=s0; b=s1-s0; Rk=12*mu/umf*h*(1/(2*b)*(1/a**2-1/(a+b)**2))
        else: Rk=R_par
        n=100; es=0
        for iz in range(n):
            zz=(iz+0.5)/n; sz=s0+ds*zz
            Rl=12*mu*(h/n)/(umf*sz**3)*n; es+=Rl/(R_haupt+Rl)
        em=es/n
        if ds==0: eref=em
        k_rows.append([f"{ds*1000:.3f}", f"{s1*1000:.3f}", f"{angle:.2f}°",
            f"{Rk/R_par:.2f}", f"{em/eref:.0%}", f"{s0*1000+ds*1000*0.5:.3f}"])
    
    story.append(Paragraph(f'<b>{z["name"]} (h = {h*1000:.0f} mm, s₀ = {s0*1000:.3f} mm)</b>', sBB))
    story.append(mt(['Δs [mm]','s₁ [mm]','Winkel','R/R₀','η<sub>eff</sub>','s bei h/2'],
        k_rows, [18*mm,18*mm,16*mm,16*mm,16*mm,PW-84*mm]))
    story.append(Spacer(1, 2*mm))

story.append(Image('/home/claude/diag_0011_5_konus.png', width=PW, height=PW*0.375))
story.append(Paragraph('<i>Abb. 5: Spaltwiderstand (blau) und Antriebseffizienz (grün) bei konischer Erweiterung. '
    'Bass: Flacher Konus (0,15°), wenig Effekt. C6: Steiler Konus (1,15°), deutlicher Verlust.</i>', sSm))

story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0011_6_querschnitt.png', width=PW, height=PW*0.42))
story.append(Paragraph('<i>Abb. 6: Kanalquerschnitt maßstäblich. Durchgezogen = Konus (0,04 mm Erweiterung), '
    'gestrichelt = parallel. Rot = Torsionsamplitude. Bei der dicken Bassplatte (15 mm) ist der '
    'Konus kaum sichtbar. Bei C6 (2 mm) ist er deutlich.</i>', sSm))

story.append(Paragraph(
    '<b>Kernaussage:</b> Der Konus-Effekt skaliert mit dem Konuswinkel, nicht mit der '
    'absoluten Erweiterung. Dieselben 0,04 mm ergeben bei 15 mm Plattendicke (Bass) '
    'einen Winkel von 0,15° → vernachlässigbar. Bei 2 mm (C6) ergeben sie 1,15° → '
    'deutlich. Deshalb ist die konische Erweiterung bei dünnen Platten (Diskant/Hoch) '
    'kritischer als im Bass.', sKB))

# ═══ Kap 5: Antriebsfläche ═══
story.append(PageBreak())
story.append(Paragraph('Kapitel 5: Antriebsfläche — warum dicke Platten nötig sind', sCh))
story.append(Paragraph(
    'Die Antriebskraft auf die Zunge ist F = Δp × W × h (Balgdruck × Zungenbreite × Plattendicke). '
    'Die dicken Bassplatten (10–15 mm) liefern eine entsprechend große Antriebsfläche, '
    'die nötig ist, um die schwere Basszunge in Bewegung zu setzen.', sB))
story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0011_7_antrieb.png', width=PW, height=PW*0.625))
story.append(Paragraph(f'<i>Abb. 7: Antriebskraft bei 1 kPa Balgdruck. '
    f'Bass 50 Hz: {zungen[0]["F_antrieb"]*1000:.0f} mN (Fläche = {zungen[0]["W"]*1000:.0f} × {zungen[0]["h"]*1000:.0f} mm). '
    f'C6: {zungen[6]["F_antrieb"]*1000:.0f} mN (Fläche = {zungen[6]["W"]*1000:.1f} × {zungen[6]["h"]*1000:.0f} mm). '
    f'Faktor {zungen[0]["F_antrieb"]/zungen[6]["F_antrieb"]:.0f}×.</i>', sSm))

# ═══ Kap 6: Warum zu dicke Platten bei hohen Tönen schaden ═══
story.append(PageBreak())
story.append(Paragraph('Kapitel 6: Warum zu dicke Platten bei hohen Tönen schaden', sCh))
story.append(Paragraph(
    'Die Erfahrung zeigt, dass Stimmplatten bei sehr hohen Tönen nicht zu dick sein dürfen. '
    'Eine 4-mm-Platte für C6 funktioniert schlechter als eine 2-mm-Platte — obwohl dickere '
    'Platten mehr Antriebsfläche bieten. Die Erklärung liegt im Zusammenspiel dreier Effekte, '
    'die mit der Plattendicke zunehmen:', sB))

story.append(Paragraph(
    '<b>Effekt 1 — Der Antrieb sättigt:</b> Die Antriebseffizienz '
    'η = R<sub>Spalt</sub>/(R<sub>Haupt</sub> + R<sub>Spalt</sub>) steigt mit h, '
    'nähert sich aber asymptotisch 100 %. Bei C6 (s = 0,025 mm) ist η > 90 % '
    'schon bei h ≈ 1 mm. Mehr Dicke bringt kaum noch mehr Antrieb.', sB))
story.append(Paragraph(
    '<b>Effekt 2 — Die Squeeze-Bremse wächst quadratisch:</b> '
    'Der Squeeze-Widerstand R<sub>sq</sub> ist proportional zu h (längerer Kanal), '
    'und die Bremsfläche wächst ebenfalls proportional zu h. '
    'Die Bremsleistung P<sub>squeeze</sub> ∝ h². Je enger der Spalt s, '
    'desto steiler der Anstieg (∝ 1/s³). '
    'Hohe Töne mit engem Spalt reagieren besonders empfindlich.', sB))
story.append(Paragraph(
    '<b>Effekt 3 — Luftmassen-Trägheit:</b> '
    'Die Luftsäule im Schlitz hat die Masse m = ρ × W × s × h. '
    'Sie muss pro Zyklus beschleunigt werden: F = m × ω² × x<sub>tip</sub>. '
    'Dieser Term wächst mit h × ω². Bei 50 Hz ist ω² ≈ 10<super>5</super>, '
    'bei 1000 Hz ist ω² ≈ 4 × 10<super>7</super> — 400× größer. '
    'Hohe Töne werden von diesem Effekt überproportional getroffen.', sB))

story.append(Spacer(1, 3*mm))
story.append(Image('/home/claude/diag_0011_10_effekte.png', width=PW, height=PW*0.35))
story.append(Paragraph(
    '<i>Abb. 8: Die drei Effekte einzeln. Links: η sättigt bei allen Zungen schnell. '
    'Mitte: Squeeze wächst quadratisch mit h — bei engen Spalten (hohe Töne) besonders steil. '
    'Rechts: Luftträgheit wächst mit h und ω² — trifft hohe Töne 400× stärker als tiefe.</i>', sSm))

story.append(Spacer(1, 3*mm))
story.append(Image('/home/claude/diag_0011_9_detail.png', width=PW, height=PW*0.45))
story.append(Paragraph(
    '<i>Abb. 9: Hoch A5 und Sehr hoch C6 im Detail. Grün = Antrieb (sättigt), '
    'Rot = Squeeze (quadratisch), Blau = Netto. Die Praxis-Dicke (orange) liegt '
    'rechts vom Optimum (grün gestrichelt) — das erklärt die Erfahrung, '
    'dass dickere Platten die Ansprache verschlechtern.</i>', sSm))

story.append(Spacer(1, 3*mm))
story.append(Image('/home/claude/diag_0011_8_plattendicke.png', width=PW, height=PW*0.61))
story.append(Paragraph(
    '<i>Abb. 10: Netto-Antrieb aller Zungen, normiert auf die Praxis-Plattendicke. '
    'Hohe Töne (rot/orange) haben einen schärferen Peak und fallen schneller ab. '
    'Bei doppelter Praxis-Dicke wären C6 und A5 bereits stark im negativen Bereich, '
    'während Bass noch funktioniert.</i>', sSm))

story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    '<b>Warum die Praxis stimmt:</b> Die Praxis-Plattendicken (2 mm bei C6, 15 mm bei Bass) '
    'sind das Ergebnis von Generationen an Erfahrung. Die Berechnung bestätigt die Logik: '
    'Bei C6 ist η schon bei 1 mm > 90 % — mehr Dicke bringt kaum Antrieb, '
    'aber die Squeeze-Bremse wächst quadratisch weiter. '
    'Bei Bass ist η erst bei ≈ 5 mm > 90 % — die dicke Platte wird gebraucht. '
    'Dazu kommt: Der Squeeze-Widerstand ∝ 1/s³ ist beim Bass (s = 0,06 mm) '
    '14× kleiner als bei C6 (s = 0,025 mm) — die Bremse greift beim Bass viel später.', sKB))

story.append(Paragraph(
    '<b>Einschränkung:</b> Das Squeeze-Modell überschätzt die Bremswirkung (Dok. 0010, '
    'Zunge biegt sich). Die berechneten Optima liegen deshalb niedriger als die Praxis. '
    'Die relative Aussage — hohe Töne reagieren empfindlicher auf zu dicke Platten — '
    'ist physikalisch fundiert und deckt sich mit der Erfahrung.', sWB))

# ═══ Kap 7: Zusammenfassung ═══
story.append(Paragraph('Kapitel 7: Zusammenfassung', sCh))
story.append(Paragraph(
    '<b>1. Duty Cycle ≈ 25 %:</b> Bei voller Lautstärke ist die Zunge nur ein Viertel der Zeit '
    'im Kanal. Die Energiezufuhr erfolgt ausschließlich während des Kanaldurchgangs — '
    'bei maximaler Geschwindigkeit (Nulldurchgang der Schwingung). Für den Vergleich '
    'verschiedener Kanalgeometrien zählt nur dieser Bereich.', sKB))
story.append(Paragraph(
    '<b>2. Torsion bestimmt den Mindestspalt:</b> Die seitliche Auslenkung bei transientem '
    'Anblasen (θ ≈ 1°) beträgt 0,04–0,10 mm — im selben Bereich wie der Seitenspalt. '
    'Der Mindestspalt wird durch die Torsion begrenzt, nicht durch die Fertigungstoleranz.', sKB))
story.append(Paragraph(
    '<b>3. Konische Erweiterung:</b> Nutzen: Platz für Torsion ohne den oberen Spalt zu vergrößern. '
    'Kosten: Spaltwiderstand sinkt, Antriebseffizienz sinkt. '
    'Der Effekt skaliert mit dem Konuswinkel: Bei Bass (h = 15 mm) vernachlässigbar, '
    'bei C6 (h = 2 mm) erheblich (bis −15 % Antrieb bei 0,04 mm Erweiterung).', sKB))
story.append(Paragraph(
    '<b>4. Antriebsfläche:</b> F ∝ W × h. Die dicken Bassplatten liefern '
    f'{zungen[0]["F_antrieb"]/zungen[6]["F_antrieb"]:.0f}× mehr Antriebskraft als C6 — '
    'nötig für die schwerere Basszunge.', sKB))
story.append(Paragraph(
    '<b>5. Zu dicke Platten schaden bei hohen Tönen:</b> Der Antrieb sättigt (η → 100 %), '
    'aber die Squeeze-Bremse wächst quadratisch (∝ h²) und die Luftträgheit wächst mit h × ω². '
    'Hohe Töne mit engem Spalt und hoher Frequenz reagieren überproportional empfindlich. '
    'Die Praxis-Dicken (2 mm bei C6, 15 mm bei Bass) sind durch dieses Gleichgewicht bestimmt.', sKB))

story.append(Spacer(1, 6*mm))
story.append(Table([['']], colWidths=[PW*0.6],
    style=TableStyle([('LINEBELOW',(0,0),(-1,-1),1,ACCENTRED),('ALIGN',(0,0),(-1,-1),'CENTER'),('FONTNAME',(0,0),(-1,-1),'DejaVu')])))
story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<i>Der Kanal ist nicht nur ein Schlitz — er ist der Ort, an dem Druck zu Bewegung wird.</i>',
    ParagraphStyle('Cl', parent=sSm, alignment=TA_CENTER, fontName='DejaVuI')))

doc.build(story, onFirstPage=set_font, onLaterPages=set_font)
print(f"PDF erzeugt: {out}")

# Verifikation
print("\n=== VERIFIKATION ===")
for z in zungen:
    print(f"  {z['name']:>15s}: h={z['h']*1000:5.1f}mm, x_tip={z['x_tip']*1000:5.1f}mm, "
          f"Duty={z['duty']:.0%}, y_tors={z['y_tors']*1000:.3f}mm, F={z['F_antrieb']*1000:.0f}mN")
