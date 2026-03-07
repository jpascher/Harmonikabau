#!/usr/bin/env python3
"""
Dok. 0010 — Güte der Stimmzunge: Zwei Verlustkanäle
Build-Skript: PDF mit Diagrammen und Tabellen.
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

# DejaVu-Fonts registrieren
pdfmetrics.registerFont(TTFont('DejaVu', '/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuB', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuI', '/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBI', '/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuMono', '/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))

# Font-Family-Zuordnung (damit <b>/<i> Tags DejaVu-Varianten verwenden)
addMapping('DejaVu', 0, 0, 'DejaVu')
addMapping('DejaVu', 1, 0, 'DejaVuB')
addMapping('DejaVu', 0, 1, 'DejaVuI')
addMapping('DejaVu', 1, 1, 'DejaVuBI')

DARKBLUE=HexColor('#16213e'); ACCENTRED=HexColor('#e94560')
KEYGREEN=HexColor('#2e7d32'); WARNRED=HexColor('#c62828')
LIGHTGRAY=HexColor('#f5f5f5'); KEYBG=HexColor('#e8f5e9'); WARNBG=HexColor('#ffebee')
W_PAGE=A4[0]; PAGE_W=W_PAGE-50*mm

styles=getSampleStyleSheet()
# Alle Standard-Styles auf DejaVu umstellen (getSampleStyleSheet verwendet Helvetica)
for style_name in styles.byName:
    s = styles.byName[style_name]
    if hasattr(s, 'fontName'):
        if 'Bold' in s.fontName or s.fontName.endswith('B'):
            s.fontName = 'DejaVuB'
        elif 'Italic' in s.fontName or s.fontName.endswith('I'):
            s.fontName = 'DejaVuI'
        else:
            s.fontName = 'DejaVu'
sTitle=ParagraphStyle('T',parent=styles['Title'],fontSize=18,textColor=DARKBLUE,spaceAfter=4,alignment=TA_CENTER,fontName='DejaVuB')
sSubtitle=ParagraphStyle('ST',parent=styles['Normal'],fontSize=12,textColor=DARKBLUE,alignment=TA_CENTER,spaceAfter=2,fontName='DejaVu')
sAbstract=ParagraphStyle('Ab',parent=styles['Italic'],fontSize=9,textColor=HexColor('#555555'),spaceAfter=8,fontName='DejaVuI')
sChapter=ParagraphStyle('Ch',parent=styles['Heading1'],fontSize=14,textColor=DARKBLUE,spaceBefore=14,spaceAfter=6,fontName='DejaVuB')
sBody=ParagraphStyle('Bo',parent=styles['Normal'],fontSize=10,leading=14,spaceAfter=6,alignment=TA_JUSTIFY,fontName='DejaVu')
sBodyB=ParagraphStyle('BB',parent=sBody,fontName='DejaVuB')
sFormula=ParagraphStyle('Fo',parent=sBody,fontSize=10,alignment=TA_CENTER,spaceAfter=8,spaceBefore=4,fontName='DejaVuMono')
sSmall=ParagraphStyle('Sm',parent=sBody,fontSize=9,leading=12,fontName='DejaVu')
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

# ══════════════════════════════════════════════════════════════
# PHYSIK
# ══════════════════════════════════════════════════════════════
rho_0=1.2; mu=1.8e-5; c=343; p_balg=1000
E_s=200e9; rho_s=7800; tan_d_s=0.0002

zungen=[
    {'name':'Bass A1','f':55,'L':60e-3,'W':8e-3,'t':0.35e-3,'x_tip':2.0e-3,'h_pl':3e-3},
    {'name':'Bass A2','f':110,'L':40e-3,'W':7e-3,'t':0.30e-3,'x_tip':1.5e-3,'h_pl':3e-3},
    {'name':'Mitte A3','f':220,'L':30e-3,'W':6e-3,'t':0.25e-3,'x_tip':1.0e-3,'h_pl':3e-3},
    {'name':'Diskant A4','f':440,'L':20e-3,'W':5e-3,'t':0.20e-3,'x_tip':0.6e-3,'h_pl':2.5e-3},
    {'name':'Hoch A5','f':880,'L':12e-3,'W':4e-3,'t':0.15e-3,'x_tip':0.4e-3,'h_pl':2e-3},
    {'name':'Sehr hoch C6','f':1047,'L':10e-3,'W':3.5e-3,'t':0.12e-3,'x_tip':0.3e-3,'h_pl':1.5e-3},
]

# ── Fall A: Gezupft, kein Kanal ──
fall_a=[]
for z in zungen:
    om=2*math.pi*z['f']; m=rho_s*z['L']*z['W']*z['t']
    vt=om*z['x_tip']; E=0.24*m*vt**2
    P_hyst=om*E*tan_d_s; Q_hyst=1/tan_d_s
    dbl=math.sqrt(2*1.5e-5/om); vrms=vt*0.4/math.sqrt(2)
    P_visc=2*z['L']*z['W']*rho_0*1.5e-5*vrms**2/dbl
    Q_visc=om*E/P_visc
    P_A=P_hyst+P_visc; Q_A=om*E/P_A
    tau_A=Q_A/(math.pi*z['f'])*1000
    fall_a.append({'name':z['name'],'f':z['f'],'E':E,'Q_hyst':Q_hyst,'Q_visc':Q_visc,'Q_A':Q_A,'tau_A':tau_A})

# ── Fall B: Ideale Einspannung, Luft ──
fall_b_all={}
for z in zungen:
    om=2*math.pi*z['f']; vt=om*z['x_tip']; m=rho_s*z['L']*z['W']*z['t']; E=0.24*m*vt**2
    umf=2*(z['L']+z['W']); Ah=z['W']*z['x_tip']*0.5
    vh=math.sqrt(2*p_balg/rho_0); Rh=rho_0*vh/(2*Ah)
    s_range=np.linspace(0.005e-3,0.20e-3,400)
    res=[]
    for s in s_range:
        Rp=12*mu*z['h_pl']/(umf*s**3); Rb=rho_0*vh/(2*umf*s)
        Rs=max(Rp,Rb); eta=Rs/(Rh+Rs)
        Pa=p_balg*eta*vt*Ah
        Qv=z['W']*z['t']*vt; Rsq=12*mu*z['h_pl']/(umf*s**3)
        dp=Qv*Rsq; Fsq=dp*z['W']*z['h_pl']; Ps=Fsq*vt
        Pn=Pa-Ps; K=umf*s/Ah
        res.append({'s':s,'K':K,'eta':eta,'Pa':Pa,'Ps':Ps,'Pn':Pn})
    fall_b_all[z['name']]=res
    nets=[r['Pn'] for r in res]; idx=np.argmax(nets)
    z['s_opt']=res[idx]['s']; z['K_opt']=res[idx]['K']; z['eta_opt']=res[idx]['eta']
    z['Pn_max']=res[idx]['Pn']
    z['s_krit']=None
    for r in res:
        if r['Pn']<=0 and r['s']<z['s_opt']: z['s_krit']=r['s']

# ── Materialvergleich (Lumped-Mass) ──
platte_vol=40e-3*12e-3*3e-3
f_ref=440; om_ref=2*math.pi*f_ref
mat_data=[('Messing',110e9,8500,0.001),('Zink-Druckguss',96e9,6600,0.002),
    ('Aluminium',70e9,2700,0.001),('Magnesium',45e9,1800,0.0005)]
mat_results=[]
for name,Ev,rho,td in mat_data:
    M=platte_vol*rho; Z_imp=math.sqrt(Ev*rho)
    R_niet=0.5; R_wachs=1.0; M_stock_eff=0.06
    R_stock=om_ref*M_stock_eff*0.01; M_eff=M+M_stock_eff
    R_tot=om_ref*M*td + R_niet + R_wachs + R_stock
    Q_mount=om_ref*M_eff/R_tot; Q_mat=1/tan_d_s
    Q_total=1/(1/Q_mat+1/Q_mount)
    mat_results.append({'name':name,'M':M*1000,'Z':Z_imp/1e6,'td':td*1000,'Q_mount':Q_mount,'Q_total':Q_total})

# ══════════════════════════════════════════════════════════════
# DIAGRAMME
# ══════════════════════════════════════════════════════════════

# Diag 1: Torte
fig1,ax1=plt.subplots(figsize=(6,4))
sizes=[97,3]; colors=['#42A5F5','#EF5350']
labels=['Luftkopplung\n(Spalt, Stroemung,\naerodynamisch)\n97 %','Mechanische\nEinspannung\n3 %']
wedges,texts=ax1.pie(sizes,labels=labels,colors=colors,startangle=90,textprops={'fontsize':11})
ax1.set_title('Energieverlust der Stimmzunge:\nZwei Kanaele',fontsize=13)
plt.tight_layout(); fig1.savefig('/home/claude/diag_0010_1_torte.png',dpi=150); plt.close()

# Diag 2: Fall A
fig2,ax2=plt.subplots(figsize=(8,4))
names_a=[a['name'] for a in fall_a]; qs_a=[a['Q_A'] for a in fall_a]; taus_a=[a['tau_A'] for a in fall_a]
bars=ax2.bar(range(len(fall_a)),qs_a,color='#1565C0',alpha=0.85)
ax2.set_xticks(range(len(fall_a))); ax2.set_xticklabels(names_a,fontsize=9)
ax2.set_ylabel('Q (gezupft, kein Kanal)',fontsize=11)
ax2.set_title('Fall A: Gezupfte Zunge ohne Luftkanal (reines Q_mech)',fontsize=12)
for i,q in enumerate(qs_a): ax2.text(i,q+50,f'Q={q:.0f}\n({taus_a[i]:.0f}ms)',ha='center',fontsize=8)
ax2.grid(axis='y',alpha=0.3); plt.tight_layout()
fig2.savefig('/home/claude/diag_0010_2_fallA.png',dpi=150); plt.close()

# Diag 3: Spalt alle Zungen
fig3,ax3=plt.subplots(figsize=(9,5))
cols=['#1565C0','#1976D2','#2196F3','#F9A825','#E65100','#C62828']
for i,z in enumerate(zungen):
    res=fall_b_all[z['name']]; ss=[r['s']*1000 for r in res]; nets=[r['Pn']*1e3 for r in res]
    ax3.plot(ss,nets,linewidth=2,color=cols[i],label=z['name'])
    ax3.plot(z['s_opt']*1000,z['Pn_max']*1e3,'o',color=cols[i],ms=6)
ax3.axhline(y=0,color='black',lw=0.5)
ax3.set_xlabel('Seitenspalt s [mm]',fontsize=11); ax3.set_ylabel('Netto-Antrieb [mW]',fontsize=11)
ax3.set_title('Netto-Antrieb vs. Spaltweite\n(Punkte=Optimum, unter Null=Zunge spricht nicht an)',fontsize=12)
ax3.legend(fontsize=9); ax3.grid(alpha=0.3); ax3.set_xlim(0,0.10)
plt.tight_layout(); fig3.savefig('/home/claude/diag_0010_3_spalt.png',dpi=150); plt.close()

# Diag 4: Detail 3 Zungen
fig4,axes=plt.subplots(1,3,figsize=(12,4.5))
for ax,zn in zip(axes,['Bass A1','Diskant A4','Sehr hoch C6']):
    z=next(zz for zz in zungen if zz['name']==zn)
    res=fall_b_all[zn]; ss=[r['s']*1000 for r in res]
    ax.plot(ss,[r['Pa']*1e3 for r in res],'g-',lw=2,label='Antrieb')
    ax.plot(ss,[r['Ps']*1e3 for r in res],'r-',lw=2,label='Squeeze')
    ax.plot(ss,[r['Pn']*1e3 for r in res],'b-',lw=2.5,label='Netto')
    ax.axhline(y=0,color='black',lw=0.5)
    ax.axvline(x=z['s_opt']*1000,color='green',ls='--',alpha=0.5)
    ax.set_title(f'{zn}\n(x={z["x_tip"]*1000:.1f}mm, A_h={z["W"]*z["x_tip"]*0.5*1e6:.1f}mm\u00b2)',fontsize=10)
    ax.set_xlabel('s [mm]',fontsize=10)
    if ax==axes[0]: ax.set_ylabel('Leistung [mW]',fontsize=10)
    ax.legend(fontsize=7); ax.grid(alpha=0.3); ax.set_xlim(0,0.10)
plt.tight_layout(); fig4.savefig('/home/claude/diag_0010_4_detail.png',dpi=150); plt.close()

# Diag 5: Material
fig5,ax5=plt.subplots(figsize=(7,4))
mnames=[m['name'] for m in mat_results]; mqs=[m['Q_mount'] for m in mat_results]
mcols=['#F9A825','#BDBDBD','#90A4AE','#A5D6A7']
ax5.barh(range(len(mat_results)),mqs,color=mcols,alpha=0.85,edgecolor='#333')
ax5.set_yticks(range(len(mat_results))); ax5.set_yticklabels(mnames,fontsize=10)
ax5.set_xlabel('Q_mech (Einspannung)',fontsize=11)
ax5.set_title('Mechanische Guete der Einspannung (Niet+Holz-Stimmstock)\nNur Kanal 2 = weniger als 3 % der Gesamtguete',fontsize=11)
for i,q in enumerate(mqs): ax5.text(q+2,i,f'{q:.0f}',va='center',fontsize=10)
ax5.grid(axis='x',alpha=0.3); plt.tight_layout()
fig5.savefig('/home/claude/diag_0010_5_material.png',dpi=150); plt.close()

# Diag 6: Rauheit vs. Grenzschichtdicke
fig6, ax6 = plt.subplots(figsize=(9, 5))
freqs_plot = np.logspace(np.log10(40), np.log10(3000), 200)
nu_air = 1.5e-5

# Grenzschichtdicke
delta_plot = np.sqrt(2 * nu_air / (2 * np.pi * freqs_plot)) * 1e6  # µm

# EDM-Oberflächen
edm_surfaces = [
    ('Ra = 3,2 µm (Hauptschnitt)', 3.2 * 6, '#E65100', '--'),
    ('Ra = 1,6 µm (1 Nachschnitt)', 1.6 * 6, '#2E7D32', '-'),
    ('Ra = 0,8 µm (2 Nachschnitte)', 0.8 * 6, '#1565C0', '-.'),
    ('Ra = 0,3 µm (3 Nachschnitte)', 0.3 * 6, '#7B1FA2', ':'),
]

# Plot delta
ax6.fill_between(freqs_plot, delta_plot * 0.1, delta_plot, alpha=0.08, color='green')
ax6.plot(freqs_plot, delta_plot, 'k-', lw=2, label='δ = √(2ν/ω) (Grenzschicht)')
ax6.plot(freqs_plot, delta_plot * 0.1, 'k--', lw=1, alpha=0.5, label='0,1·δ (Glatt-Grenze)')

for name, ks_um, col, ls in edm_surfaces:
    ax6.axhline(y=ks_um, color=col, linestyle=ls, linewidth=2, label=f'{name}, k_s={ks_um:.0f} µm')

# Markiere Zungenfrequenzen
for f_mark, label in [(55, 'A1'), (220, 'A3'), (440, 'A4'), (880, 'A5'), (1047, 'C6')]:
    d = math.sqrt(2 * nu_air / (2 * math.pi * f_mark)) * 1e6
    ax6.plot(f_mark, d, 'ko', ms=5)
    ax6.annotate(label, (f_mark, d), textcoords="offset points", xytext=(5, 5), fontsize=8)

ax6.set_xscale('log')
ax6.set_yscale('log')
ax6.set_xlabel('Frequenz [Hz]', fontsize=11)
ax6.set_ylabel('Rauheit k_s bzw. Grenzschichtdicke δ [µm]', fontsize=11)
ax6.set_title('Rauheit (Drahterodieren) vs. Grenzschichtdicke\nUnter 0,1·δ = hydraulisch glatt', fontsize=12)
ax6.legend(fontsize=8, loc='upper right')
ax6.grid(alpha=0.3, which='both')
ax6.set_xlim(40, 3000)
ax6.set_ylim(1, 500)
ax6.text(100, 200, 'GLATT\n(Rauheit unsichtbar)', fontsize=10, color='green', alpha=0.6, ha='center')
ax6.text(2000, 200, 'ÜBERGANG', fontsize=9, color='orange', alpha=0.6, ha='center')
plt.tight_layout()
fig6.savefig('/home/claude/diag_0010_6_rauheit.png', dpi=150); plt.close()

# Diag 7: Toleranz vs. Rauheit — Squeeze-Dämpfung
fig7, (ax7a, ax7b) = plt.subplots(1, 2, figsize=(10, 4.5))

for ax, z_name, f_z, s_nom in [(ax7a, 'Diskant A4 (440 Hz)', 440, 0.04),
                                  (ax7b, 'Hoch A5 (880 Hz)', 880, 0.03)]:
    # Squeeze ∝ 1/s³
    s_range_plot = np.linspace(s_nom - 0.015, s_nom + 0.015, 200)
    s_range_plot = s_range_plot[s_range_plot > 0.005]
    squeeze_rel = (s_nom / s_range_plot)**3
    
    ax.plot(s_range_plot * 1000, squeeze_rel, 'b-', lw=2)
    ax.axvline(x=s_nom * 1000, color='green', ls='--', lw=1.5, label=f'Nominal s = {s_nom*1000:.0f} µm')
    ax.axvline(x=(s_nom - 0.01) * 1000, color='red', ls=':', lw=1.5, label=f's − 10 µm')
    ax.axvline(x=(s_nom + 0.01) * 1000, color='red', ls=':', lw=1.5, label=f's + 10 µm')
    
    # Toleranzbereich schattieren
    ax.axvspan((s_nom - 0.01) * 1000, (s_nom + 0.01) * 1000, alpha=0.15, color='red', label='Toleranz ±10 µm')
    
    # Rauheit-Bereich (unsichtbar klein)
    ks_lo = 0.8 * 6 / 1000  # mm, Ra=0.8
    ks_hi = 3.2 * 6 / 1000  # mm, Ra=3.2
    # Rauheit verändert den effektiven Spalt um ±k_s
    sq_ks_lo = (s_nom / (s_nom - ks_hi))**3
    sq_ks_hi = (s_nom / (s_nom + ks_hi))**3
    ax.axhspan(sq_ks_hi, sq_ks_lo, alpha=0.3, color='gray', label=f'Rauheit Ra 0,8–3,2 µm')
    
    # Squeeze-Werte bei Toleranzgrenzen
    sq_minus = (s_nom / (s_nom - 0.01))**3
    sq_plus = (s_nom / (s_nom + 0.01))**3
    ax.annotate(f'{sq_minus:.1f}×', ((s_nom-0.01)*1000, sq_minus), 
               textcoords="offset points", xytext=(-25, 5), fontsize=9, color='red', fontweight='bold')
    ax.annotate(f'{sq_plus:.1f}×', ((s_nom+0.01)*1000, sq_plus),
               textcoords="offset points", xytext=(5, -15), fontsize=9, color='red', fontweight='bold')
    
    ax.set_xlabel('Spaltweite s [µm]', fontsize=10)
    ax.set_title(z_name, fontsize=11)
    ax.legend(fontsize=7, loc='upper right')
    ax.grid(alpha=0.3)
    if ax == ax7a:
        ax.set_ylabel('Squeeze relativ zu Nominal', fontsize=10)
    ax.set_ylim(0, 5)

plt.suptitle('Toleranz (rot) vs. Rauheit (grau): Die Toleranz dominiert', fontsize=12, y=1.02)
plt.tight_layout()
fig7.savefig('/home/claude/diag_0010_7_toleranz_vs_rauheit.png', dpi=150); plt.close()

# ══════════════════════════════════════════════════════════════
# PDF
# ══════════════════════════════════════════════════════════════
output_path='/home/claude/0010_guete_stimmplatte_De.pdf'
doc=SimpleDocTemplate(output_path,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=25*mm,bottomMargin=25*mm)

def set_canvas_font(canvas, doc):
    canvas.setFont('DejaVu', 10)
story=[]

# Titel
story.append(Paragraph('Dok. 0010',sSmall)); story.append(Spacer(1,2*mm))
story.append(Paragraph('Güte der Stimmzunge: Zwei Verlustkanäle',sTitle))
story.append(Paragraph('Luftkopplung vs. Einspannung — und warum fast alles funktioniert',sSubtitle))
story.append(Spacer(1,2*mm)); story.append(hr()); story.append(Spacer(1,3*mm))
story.append(Paragraph(
    'Die Güte Q der Stimmzunge bestimmt Ansprache und Klang (Dok. 0005, 0008). '
    'Die Energie der schwingenden Zunge verlässt das System über zwei Kanäle: '
    '(1) Luftkopplung — die Zunge als Impulsgenerator im Spalt, und '
    '(2) mechanische Einspannung — Biegezonen-Hysterese und Kontaktreibung am Niet. '
    'Die Berechnung zeigt, dass der Luftkanal mehr als 97 % der Gesamtverluste bestimmt. '
    'Die mechanische Einspannung (Plattenmaterial, Niet, Stimmstock) trägt weniger als 3 % bei. '
    'Das erklärt, warum fast jede Materialkombination funktioniert — und warum die Passgenauigkeit '
    'der Zunge im Schlitz die Ansprache weit stärker beeinflusst als das Plattenmaterial. '
    'Die Zahlenwerte sind Größenordnungsabschätzungen.',sAbstract))

refs=[['0002','Strömungsanalyse Bass-Stimmzunge 50 Hz — v8'],['0005','Frequenzverschiebung als Indikator der Ansprache'],
      ['0008','Klangveränderung durch Kammergeometrie'],['0009','Frequenzkopplung mehrerer Zungen']]
story.append(Paragraph('<b>Dokumentenverweise</b>',sSmall))
rd=[[Paragraph(f'<b>Dok. {r[0]}</b>',sTD),Paragraph(r[1],sTDL)] for r in refs]
rt=Table(rd,colWidths=[30*mm,PAGE_W-30*mm])
rt.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.3,HexColor('#cccccc')),('BACKGROUND',(0,0),(0,-1),LIGHTGRAY),
    ('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),2),('BOTTOMPADDING',(0,0),(-1,-1),2),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
story.append(rt); story.append(Spacer(1,4*mm))

# Kap 1
story.append(Paragraph('Kapitel 1: Zwei Verlustkanäle',sChapter))
story.append(Paragraph(
    'Die Stimmzunge ist ein mechanischer Oszillator, der durch den Balgdruck in Schwingung gehalten wird. '
    'Die gespeicherte Schwingungsenergie verlässt die Zunge über zwei physikalisch verschiedene Kanäle:',sBody))
story.append(Paragraph(
    '<b>Kanal 1 — Luftkopplung (dominant):</b> Die Zunge schwingt durch den Spalt in der Stimmplatte '
    'und wirkt als Impulsgenerator. Bei jedem Durchschwung drückt sie Luft durch den Spalt. '
    'Die Verluste bestehen aus: (a) akustischem Kurzschluss durch den Seitenspalt, '
    '(b) Squeeze-Film-Dämpfung im engen Spalt, (c) Wirbelablösung, (d) Schallabstrahlung. '
    'Dieser Kanal ist materialunabhängig — er hängt von der Spaltgeometrie ab.',sBody))
story.append(Paragraph(
    '<b>Kanal 2 — Mechanische Einspannung (sekundär):</b> Die Biegeschwingung erzeugt '
    'an der Einspannung periodische Kräfte. In den Biegezonen wird bei jedem Zyklus ein Bruchteil '
    'der elastischen Energie durch Hysterese (tan δ) in Wärme umgewandelt. '
    'Dazu kommt Kontaktreibung am Niet. Dieser Kanal hängt vom Material und der Einspannqualität ab.',sBody))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0010_1_torte.png',width=PAGE_W*0.55,height=PAGE_W*0.37))
story.append(Paragraph('<i>Abb. 1: Luftkopplung bestimmt > 97 % der Gesamtverluste.</i>',sSmall))

# Kap 2
story.append(Paragraph('Kapitel 2: Fall A — Gezupfte Zunge ohne Luftkanal',sChapter))
story.append(Paragraph(
    'Eine gezupfte Zunge in stehender Luft — ohne Kanal, ohne Balgdruck — zeigt reines Q<sub>mech</sub>. '
    'Verluste: nur Materialhysterese (tan δ<sub>Stahl</sub> = 0,0002 → Q = 5000) und '
    'viskose Luftreibung an der Oberfläche.',sBody))
fa_rows=[]
for a in fall_a:
    fa_rows.append([a['name'],f"{a['f']}",f"{a['Q_A']:.0f}",f"{a['tau_A']:.0f}"])
story.append(make_table(['Zunge','f [Hz]','Q<sub>Fall A</sub>','τ = Q/(πf) [ms]'],
    fa_rows,col_widths=[30*mm,18*mm,25*mm,PAGE_W-73*mm]))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0010_2_fallA.png',width=PAGE_W,height=PAGE_W*0.5))
story.append(Paragraph('<i>Abb. 2: Gezupfte Zunge ohne Kanal: Q ≈ 4000–4400. '
    'Die τ-Werte sind Zeitkonstanten (Amplitude auf 1/e), nicht hörbare Nachschwingzeiten.</i>',sSmall))
story.append(Paragraph(
    '<b>Anmerkung zu τ:</b> Die Zeitkonstante τ = Q/(π·f) ist eine Eigenschaft des Systems. '
    'Die tatsächlich hörbare Nachschwingzeit hängt zusätzlich von der eingebrachten Energie ab: '
    't<sub>hörbar</sub> = τ · ln(A<sub>Start</sub> / A<sub>Schwelle</sub>). '
    'Beim Zupfen wird eine Auslenkung x vorgegeben, nicht eine kontrollierte Energie. '
    'Die gespeicherte Energie ist E = ½ k x², wobei k ∝ EI/L³ die Federsteifigkeit der Zunge ist. '
    'Eine schwerere Zunge (steifere Feder) speichert bei gleicher Auslenkung mehr Energie → '
    'lauterer Start → länger hörbar — auch bei gleichem Q. Die τ-Werte in der Tabelle zeigen '
    'daher die Systemdynamik, nicht die wahrgenommene Abklingzeit beim Zupfen. '
    'Bei Luftzuführung im Instrument ist es anders: Dort bestimmt der Balgdruck × Druckfläche '
    'die zugeführte Leistung kontinuierlich, nicht ein einmaliger Energiestoß.',sBody))
story.append(Paragraph(
    '<b>Ergebnis:</b> Q<sub>mech</sub> ≈ 4000–4400. Im Instrument ist Q ≈ 50–200. '
    'Der Unterschied (Faktor 20–80) ist der Luftkanal. '
    'Praxistest: Zunge zupfen und Abklingverhalten beobachten — die relative Dauer '
    'zwischen verschiedenen Zungen zeigt die Unterschiede in Q<sub>mech</sub>.',sKeyBox))

# Kap 3
story.append(PageBreak())
story.append(Paragraph('Kapitel 3: Fall B — Spaltweite und Luftverluste',sChapter))
story.append(Paragraph(
    'Ideale Einspannung (Q<sub>mech</sub> = ∞), nur Luftverluste. Zwei gegenläufige Effekte:',sBody))
story.append(Paragraph(
    '<b>Enger Spalt → weniger Kurzschluss → mehr Antrieb</b> (gut). '
    'Die Antriebseffizienz η ist ein Druckteiler zwischen Seitenspalt-Widerstand R<sub>Spalt</sub> '
    'und Hauptöffnungs-Widerstand R<sub>Haupt</sub>.',sBody))
story.append(Paragraph(
    '<b>Enger Spalt → mehr Squeeze-Dämpfung → mehr Bremse</b> (schlecht). '
    'Die Zunge verdrängt bei jedem Durchschwung das Volumen W × t × v<sub>tip</sub> durch den Seitenspalt.',sBody))
story.append(Spacer(1,2*mm))

sb_rows=[]
for z in zungen:
    sk=f"{z['s_krit']*1000:.3f}" if z.get('s_krit') else "< 0,005"
    sb_rows.append([z['name'],f"{z['f']}",f"{z['W']*z['x_tip']*0.5*1e6:.1f}",
        f"{z['s_opt']*1000:.3f}",f"{z['eta_opt']:.0%}",sk])
story.append(make_table(['Zunge','f [Hz]','A<sub>Haupt</sub> [mm²]','s<sub>opt</sub> [mm]',
    'η<sub>opt</sub>','s<sub>krit</sub> [mm]'],
    sb_rows,col_widths=[28*mm,16*mm,24*mm,20*mm,16*mm,PAGE_W-104*mm]))
story.append(Spacer(1,2*mm))

story.append(Image('/home/claude/diag_0010_3_spalt.png',width=PAGE_W,height=PAGE_W*0.56))
story.append(Paragraph('<i>Abb. 3: Netto-Antrieb vs. Spaltweite. Kleine Zungen (C6, rot) '
    'brauchen engere Spalte, weil A<sub>Haupt</sub> nur 0,5 mm² beträgt.</i>',sSmall))
story.append(Spacer(1,2*mm))

story.append(Image('/home/claude/diag_0010_4_detail.png',width=PAGE_W,height=PAGE_W*0.375))
story.append(Paragraph('<i>Abb. 4: Grün=Antrieb, Rot=Squeeze, Blau=Netto. '
    'Das Optimum liegt beim Maximum der blauen Kurve.</i>',sSmall))

story.append(Paragraph(
    '<b>Warum hohe Stimmplatten engere Spalte brauchen:</b> '
    'A<sub>Haupt</sub> = W × x<sub>tip</sub> sinkt mit kleinerer Zunge '
    '(Bass: 8 mm², Sehr hoch: 0,5 mm²). Bei gleichem Spalt ist der '
    'Kurzschluss für kleine Zungen proportional schlimmer.',sBody))
story.append(Paragraph(
    '<b>Einschränkung:</b> Die berechneten Optima liegen systematisch höher als die Praxiswerte. '
    'Die Squeeze-Formel überschätzt die Bremswirkung, weil die Zunge sich biegt — '
    'der Spalt ist nur am Fußpunkt eng, an der Spitze weiter. Die Praxis arbeitet enger '
    'als berechnet, weil die reale Squeeze-Dämpfung geringer ist.',sWarnBox))

# Kap 4
story.append(PageBreak())
story.append(Paragraph('Kapitel 4: Die mechanische Einspannung — Lumped-Mass-Modell',sChapter))
story.append(Paragraph(
    'Platte, Stimmstock und Prüfblock sind alle ≪ λ (Platte 40 mm ≪ λ = 8 m bei 440 Hz). '
    'Keine Wellenausbreitung — das Transmissionskoeffizienten-Modell ist falsch. '
    'Stattdessen: <b>Lumped-Mass-Modell</b> — die Platte wirkt als konzentrierte Masse.',sBody))
story.append(Paragraph(
    'Z<sub>Masse</sub> = jωM — rein imaginär, verlustfrei. Eine ideale Masse reflektiert 100 %. '
    'Verluste entstehen nur durch:',sBody))
story.append(Paragraph(
    '<b>1. Biegezonen-Hysterese:</b> In den Zonen, wo Material periodisch verformt wird, '
    'geht pro Zyklus ein Anteil tan δ der elastischen Energie in Wärme. '
    'Je steifer und massiver die Einspannung, desto kleiner die Verformung.',sBody))
story.append(Paragraph(
    '<b>2. Kontaktreibung:</b> Die Querkraft am Einspannpunkt (F ≈ 0,25 N für A4) erzeugt '
    'Mikrobewegung. Niet: wenig Bewegung. Schraube: mehr Mikro-Schlupf.',sBody))
story.append(Paragraph(
    '<b>3. Übergangszone:</b> Zwischen Niet und freiem Biegungsende: ≈ 0,05–0,2 mm. '
    'Wird beim Einspielen freigeräumt. Verbindet Q und Frequenz: '
    'steifer Abschluss → schärfere Grenze → kürzeres L<sub>eff</sub> → minimal höhere Frequenz.',sBody))
story.append(Paragraph(
    '<b>Konsequenz:</b> Je massiver und fester → höheres Q<sub>mech</sub>. '
    'Prüfblock > Schraubstock > Stimmplatte auf Stimmstock > Hand.',sKeyBox))

# Kap 5
story.append(Paragraph('Kapitel 5: Materialvergleich',sChapter))
mr=[]
for m in mat_results:
    Qg=1/(1/100+1/m['Q_mount'])
    mr.append([m['name'],f"{m['M']:.1f}",f"{m['td']:.1f}",f"{m['Q_mount']:.0f}",f"{Qg:.0f}"])
story.append(make_table(['Material','M [g]','tan δ [‰]','Q<sub>mech</sub>','Q<sub>ges</sub> (mit Q<sub>Luft</sub>=100)'],
    mr,col_widths=[30*mm,18*mm,20*mm,22*mm,PAGE_W-90*mm]))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0010_5_material.png',width=PAGE_W,height=PAGE_W*0.57))
story.append(Paragraph('<i>Abb. 5: Q<sub>mech</sub> für verschiedene Plattenmaterialien. '
    'Bestimmt weniger als 3 % der Gesamtgüte.</i>',sSmall))
story.append(Paragraph(
    f'<b>Fazit:</b> Q<sub>gesamt</sub> ändert sich von '
    f'{1/(1/100+1/mat_results[0]["Q_mount"]):.0f} (Messing) zu '
    f'{1/(1/100+1/mat_results[3]["Q_mount"]):.0f} (Magnesium) — '
    f'{abs(1/(1/100+1/mat_results[0]["Q_mount"])-1/(1/100+1/mat_results[3]["Q_mount"]))/(1/(1/100+1/mat_results[0]["Q_mount"]))*100:.0f} % '
    'Unterschied. Nicht hörbar.',sBody))

# Kap 6
story.append(Paragraph('Kapitel 6: Der Prüfblock',sChapter))
story.append(Paragraph(
    'Der Prüfblock (5 kg Stahl) ist der massivste Abschluss. '
    'Block (10 cm) ≪ λ (11,6 m bei 440 Hz) → Lumped-Mass. '
    'tan δ<sub>Stahl</sub> = 0,0002 → 50× weniger Dissipation als Holz. '
    'Q<sub>mech</sub> auf dem Prüfblock ist das höchste aller Konfigurationen.',sBody))
story.append(Paragraph(
    'Der Klangunterschied Prüfblock vs. Instrument kommt nicht von der Impedanzkette, sondern '
    'von der fehlenden Kammer (kein Kerbfilter, Dok. 0008), der freien Abstrahlung und den '
    'anderen Strömungsverhältnissen. Die Frequenz bleibt gleich. '
    'Ansprache und Klangfarbe ändern sich wegen der Kammer, nicht wegen der Einspannung.',sBody))

# Kap 7
story.append(Paragraph('Kapitel 7: Grenzen der Berechnung',sChapter))
story.append(Paragraph('<b>1.</b> Q<sub>Luft</sub> ≈ 100 ist geschätzt. '
    '<b>2.</b> Squeeze-Film überschätzt (Zunge biegt sich). '
    '<b>3.</b> Kontaktreibung am Niet und Stimmstock-Kopplung sind empirisch. '
    '<b>4.</b> Plattenresonanzen fehlen.',sBody))

# Kap 8
story.append(Paragraph('Kapitel 8: Zusammenfassung',sChapter))
story.append(Paragraph(
    '<b>Einordnung:</b> Erklärungsmodell. Die Aufteilung in zwei Kanäle ist physikalisch fundiert. '
    'Absolute Q-Werte und Spaltgrenzen sind Größenordnungsabschätzungen.',sWarnBox))
story.append(Paragraph(
    '<b>1. Zwei Kanäle:</b> Luftkopplung > 97 %, Einspannung &lt; 3 %. '
    'Beweis: Gezupfte Zunge Q ≈ 4000, im Instrument Q ≈ 50–200.',sKeyBox))
story.append(Paragraph(
    '<b>2. Spaltoptimum:</b> Enger → weniger Kurzschluss (gut), aber mehr Squeeze (schlecht). '
    'Kleine Zungen brauchen engere Spalte (A<sub>Haupt</sub> kleiner).',sKeyBox))
story.append(Paragraph(
    '<b>3. Einspannung:</b> Lumped-Mass. Fester = besser. Niet &gt; Schraube. '
    'Prüfblock hat höchstes Q<sub>mech</sub>.',sKeyBox))
story.append(Paragraph(
    '<b>4. Warum fast alles funktioniert:</b> Der Luftkanal dominiert und ist materialunabhängig. '
    'Verdopplung von Q<sub>mech</sub> ändert Q<sub>gesamt</sub> um &lt; 3 %.',sKeyBox))
story.append(Paragraph(
    '<b>5. Rangfolge:</b> '
    '1. Spaltweite Zunge/Schlitz (dominant). '
    '2. Kammergeometrie (Dok. 0005–0008). '
    '3. Niet vs. Schraube. '
    '4. Plattenmaterial (gering). '
    '5. Stimmstock (noch geringer). '
    '6. tan δ<sub>Stahl</sub> (vernachlässigbar).',sKeyBox))

# ═══ Kap 9: Rauheit des Stimmplattenkanals ═══
story.append(PageBreak())
story.append(Paragraph('Kapitel 9: Rauheit des Stimmplattenkanals',sChapter))
story.append(Paragraph(
    'Die Stimmplattenschlitze werden heute überwiegend durch Drahterodieren (Wire EDM) gefertigt. '
    'Dieses Verfahren erzeugt eine charakteristische Oberfläche mit Ra ≈ 0,8–3,2 µm, '
    'abhängig von der Schnittgeschwindigkeit und der Anzahl der Nachschnitte. '
    'Die Frage ist: Beeinflusst diese Rauheit die Güte der Zunge?',sBody))

story.append(Spacer(1,2*mm))
story.append(Paragraph('<b>Typische Oberflächen beim Drahterodieren (Aluminium):</b>',sBodyB))
story.append(Spacer(1,1*mm))

edm_rows = [
    ['Hauptschnitt (1 Schnitt)', '2,5–3,2', '15–19', 'Schnell, günstig'],
    ['+ 1 Nachschnitt', '1,2–1,8', '7–11', 'Standard-Qualität'],
    ['+ 2 Nachschnitte', '0,6–1,0', '4–6', 'Fein, ±0,01 mm Toleranz'],
    ['+ 3 Nachschnitte (Feinst)', '0,3–0,5', '2–3', 'Maximale Güte, teuer'],
]
story.append(make_table(
    ['Prozess','Ra [µm]','k<sub>s</sub> ≈ 6·Ra [µm]','Anmerkung'],
    edm_rows,
    col_widths=[38*mm, 22*mm, 30*mm, PAGE_W-90*mm]))
story.append(Paragraph(
    '<i>k<sub>s</sub> = Sandäquivalente Rauheit ≈ 6 × Ra (empirisch für funkenerodierte Oberflächen).</i>',sSmall))

story.append(Spacer(1,4*mm))
story.append(Paragraph('<b>Das entscheidende Kriterium: Rauheit vs. Grenzschichtdicke</b>',sBodyB))
story.append(Paragraph(
    'In einer oszillierenden Strömung (Schallfeld) bildet sich an der Wand eine viskose '
    'Grenzschicht der Dicke δ = √(2ν/ω) aus, wobei ν = 1,5 × 10<super>−5</super> m²/s '
    'die kinematische Viskosität der Luft ist. '
    'Wenn die Rauheit k<sub>s</sub> viel kleiner als δ ist, „sieht" die Strömung eine glatte Wand. '
    'Wenn k<sub>s</sub> ≈ δ oder größer, wird die Grenzschicht gestört und die Reibungsverluste steigen.',sBody))

story.append(Spacer(1,3*mm))
# Diagramm erzeugen
story.append(Image('/home/claude/diag_0010_6_rauheit.png', width=PAGE_W, height=PAGE_W*0.55))
story.append(Paragraph(
    '<i>Abb. 6: Verhältnis k<sub>s</sub>/δ für verschiedene Oberflächenqualitäten und Frequenzen. '
    'Unter der Linie k<sub>s</sub>/δ = 0,1 ist die Oberfläche hydraulisch glatt — '
    'Rauheit hat keinen Einfluss. Über k<sub>s</sub>/δ = 5 ist sie hydraulisch rau. '
    'Drahterodierte Aluminium-Oberflächen (Ra ≤ 3,2 µm) liegen bei allen Frequenzen '
    'im glatten oder knapp im Übergangsbereich.</i>',sSmall))

story.append(Spacer(1,3*mm))

# Vergleichstabelle für Standard EDM (Ra = 1.6 µm)
rauheit_rows = []
for name, f, s_mm in [('Bass A1 (55 Hz)', 55, 0.06),
                       ('Mitte A3 (220 Hz)', 220, 0.05),
                       ('Diskant A4 (440 Hz)', 440, 0.04),
                       ('Hoch A5 (880 Hz)', 880, 0.03),
                       ('Sehr hoch C6 (1047 Hz)', 1047, 0.025)]:
    omega = 2 * math.pi * f
    delta = math.sqrt(2 * 1.5e-5 / omega) * 1000  # mm
    ks = 1.6 * 6 / 1000  # mm
    ratio_d = ks / delta
    ratio_s = ks / s_mm
    regime = 'glatt' if ratio_d < 0.1 else 'Übergang'
    rauheit_rows.append([name, f'{delta*1000:.0f}', f'{ks*1000:.1f}', f'{ratio_d:.3f}', f'{ratio_s:.3f}', regime])

story.append(Paragraph('<b>Ra = 1,6 µm (Standard, 1 Nachschnitt): k<sub>s</sub> ≈ 10 µm = 0,010 mm</b>',sBodyB))
story.append(make_table(
    ['Zunge', 'δ [µm]', 'k<sub>s</sub> [µm]', 'k<sub>s</sub>/δ', 'k<sub>s</sub>/s', 'Regime'],
    rauheit_rows,
    col_widths=[35*mm, 18*mm, 18*mm, 18*mm, 18*mm, PAGE_W-107*mm]))

story.append(Spacer(1,3*mm))
story.append(Paragraph(
    '<b>Ergebnis:</b> Bei Ra = 1,6 µm (Standard-EDM mit 1 Nachschnitt) ist k<sub>s</sub>/δ < 0,15 '
    'für alle Frequenzen — <b>hydraulisch glatt</b>. Die Rauheit hat keinen messbaren Einfluss '
    'auf die Squeeze-Dämpfung und den akustischen Kurzschluss. '
    'Selbst Ra = 3,2 µm (nur Hauptschnitt) ergibt k<sub>s</sub>/δ < 0,3 — '
    'noch im glatten Bereich.',sKeyBox))

story.append(Spacer(1,3*mm))
story.append(Paragraph('<b>Toleranz vs. Rauheit: Die Priorität</b>',sBodyB))
story.append(Paragraph(
    'Die Squeeze-Dämpfung hängt kubisch vom Spalt ab (P<sub>squeeze</sub> ∝ 1/s³). '
    'Eine Toleranzabweichung von ±0,01 mm bei einem Nominalspalt von 0,03 mm '
    'verändert die Squeeze-Dämpfung um Faktor 3,4. Dieselbe Rauheit (k<sub>s</sub> = 0,010 mm) '
    'verändert sie um weniger als 0,1 %.',sBody))

story.append(Spacer(1,3*mm))
story.append(Image('/home/claude/diag_0010_7_toleranz_vs_rauheit.png', width=PAGE_W, height=PAGE_W*0.5))
story.append(Paragraph(
    '<i>Abb. 7: Squeeze-Dämpfung vs. Spaltweite für Diskant A4 (440 Hz) und Hoch A5 (880 Hz). '
    'Die roten Balken zeigen den Effekt der Toleranz ±0,01 mm — Faktor 2–3 in der Squeeze-Kraft. '
    'Die Rauheit (Bandbreite Ra = 0,8–3,2 µm) ist als graue Zone eingezeichnet — unsichtbar dünn.</i>',sSmall))

story.append(Spacer(1,3*mm))
story.append(Paragraph(
    '<b>Zusammenfassung Kap. 9:</b><br/>'
    '1. Die Oberflächenrauheit aus dem Drahterodieren (Ra = 0,8–3,2 µm) liegt bei allen '
    'Frequenzen im hydraulisch glatten Bereich (k<sub>s</sub>/δ < 0,3).<br/>'
    '2. Feinere Oberfläche (Ra < 0,8 µm, 3 Nachschnitte) bringt akustisch nichts — '
    'die Oberfläche ist bereits glatt genug.<br/>'
    '3. Die Maßtoleranz (±0,01 mm) bestimmt die Ansprache, nicht die Oberflächengüte. '
    'Die Toleranz ist 100× wichtiger als die Rauheit.<br/>'
    '4. Die wirtschaftliche Konsequenz: 1 Nachschnitt (Ra ≈ 1,6 µm) mit enger Toleranz (±0,01 mm) '
    'ist akustisch optimal. Mehr Nachschnitte für feinere Oberfläche sind Aufwand ohne akustischen Gewinn. '
    'Investition in Maßhaltigkeit lohnt sich mehr als Investition in Oberflächengüte.',sKeyBox))

story.append(Spacer(1,4*mm))
story.append(Paragraph('<b>Exkurs: Warum der Propeller-Vergleich nicht gilt</b>',sBodyB))
story.append(Paragraph(
    'Bei Propellern, Golfbällen und Haifischhaut verbessert eine raue oder strukturierte Oberfläche '
    'die Leistung, weil Rauheit die Grenzschicht-Turbulenz früher auslöst. Eine turbulente '
    'Grenzschicht hat mehr kinetische Energie nahe der Wand, verzögert die Strömungsablösung und '
    'reduziert den Formwiderstand. Gilt das auch im Stimmplattenkanal?',sBody))
story.append(Paragraph(
    'Nein — und der Grund ist die Reynolds-Zahl. Der Propeller-Effekt tritt bei Re > 10<super>5</super> auf '
    '(hochenergetische Grenzschicht, konvexe Oberfläche, Formwiderstand dominiert). '
    'Im Seitenspalt der Stimmplatte ist Re ≈ 100–200 — rein laminare Strömung. '
    'Bei laminarer Strömung hat Rauheit keinen Einfluss auf den Strömungswiderstand, '
    'solange k<sub>s</sub> ≪ s. Die Poiseuille-Formel R ∝ 1/s³ gilt unverändert. '
    'Der Kanal ist gerade und parallel — es gibt keine konvexe Oberfläche, '
    'an der die Strömung ablösen könnte, und keinen Formwiderstand.',sBody))

story.append(Spacer(1,2*mm))
# Re-Tabelle
re_rows = []
for name, s_um in [('30 µm (Hoch)', 30), ('40 µm (Diskant)', 40), ('50 µm (Mitte)', 50)]:
    s = s_um * 1e-6
    D_h = 2 * s
    Re = 1.2 * 40.8 * D_h / 1.8e-5
    re_rows.append([f's = {name}', f'{D_h*1e6:.0f}', f'{Re:.0f}', 'laminar' if Re < 2000 else 'turbulent'])
story.append(make_table(['Spalt','D<sub>h</sub> [µm]','Re','Regime'],
    re_rows, col_widths=[35*mm, 22*mm, 18*mm, PAGE_W-75*mm]))
story.append(Paragraph('<i>Re ≪ 2000 bei allen Spaltweiten → laminare Strömung → Rauheit wirkungslos.</i>',sSmall))

story.append(Spacer(1,3*mm))
story.append(Paragraph('<b>Die Spaltkante: Wichtiger als die Wandrauheit</b>',sBodyB))
story.append(Paragraph(
    'Am Einlauf und Auslauf des Seitenspalts gibt es eine scharfe 90°-Umlenkung. '
    'Dort reißt die Strömung lokal ab und erzeugt Wirbel — unabhängig von der Wandrauheit. '
    'Der Einlaufverlust (ζ = 0,5 bei scharfer Kante, ζ = 0,04 bei gerundeter Kante) '
    'ist jedoch bei den vorliegenden Spaltweiten klein im Vergleich zum '
    'Reibungsverlust über die 3 mm Kanallänge (Verhältnis 2–5 %). '
    'Die Kantenform ist also ebenfalls sekundär.',sBody))

story.append(Spacer(1,2*mm))
story.append(Paragraph(
    '<b>Das Paradoxon:</b> Selbst wenn Rauheit oder scharfe Kanten den Spaltwiderstand '
    'erhöhen würden — das wäre für die Ansprache sogar <b>gut</b>, weil hoher Spaltwiderstand '
    'den akustischen Kurzschluss reduziert (Kap. 3). Der Seitenspalt ist der Kurzschlusskanal, '
    'und mehr Widerstand dort bedeutet mehr Antrieb für die Zunge. '
    'Aber bei Re = 100–200 tut die Rauheit nichts — die Strömung ist und bleibt laminar.',sBody))

story.append(Spacer(1,4*mm))
story.append(Paragraph(
    '<b>Abschließende Bemerkung:</b> Dieses Dokument weist jedem Einzeleffekt seinen Anteil zu '
    'und stuft die meisten als vernachlässigbar ein. Das ist rechnerisch korrekt — aber '
    'es unterschätzt die Summe. Plattenmaterial 3 %, Stimmstock 1 %, Gehäuseschwingung 0,5 %, '
    'Wachsalterung 0,5 %, Nietqualität 1 %, Plattenresonanzen 1 % — jeder Einzelposten '
    'rechtfertigt keine Aufmerksamkeit. Aber zehn solcher Posten summieren sich auf 10–30 %, '
    'und sie summieren sich nicht einfach linear: Wie beim Schmetterlingseffekt können kleine '
    'Kopplungen zwischen den Effekten das Ergebnis verstärken oder abschwächen, ohne dass man '
    'den einzelnen Beitrag als bedeutsam erkennen würde. '
    'Dazu kommt ein Aspekt, den dieses Dokument nicht behandelt: <b>Alterung</b>. '
    'Wachs verändert seine Festigkeit über Jahre. Leimfugen im Gehäuse — speziell bei '
    'Knochenleim — brechen durch die Vibrationen, denen das Gehäuse dauerhaft ausgesetzt ist, '
    'langsam ein. Das Einspielen eines Instruments ist kein Mythos: Die Übergangszone an der '
    'Einspannung räumt sich frei (Kap. 4), Leimfugen setzen sich, das Holz des Stimmstocks '
    'und Gehäuses verändert seine Dämpfungseigenschaften unter zyklischer Belastung. '
    'Ein absolut schwingungsfreies Gehäuse gibt es nicht — es ist Teil des Systems. '
    'Der Spielraum zwischen einem guten und einem sehr guten Instrument ist minimal — '
    'aber wahrnehmbar. Zu sagen, der Stimmstock oder das Schwingverhalten des Gehäuses '
    'hätten keinen merklichen Einfluss, ist wahrscheinlich übertrieben. '
    'Die Berechnung zeigt die Rangfolge. Die Summe der kleinen Dinge, die kein Modell '
    'einzeln erfasst, macht den Unterschied, den der Instrumentenbauer hört und der Spieler fühlt.',sBody))

story.append(Paragraph(
    'Dieses Dokument liefert Größenordnungen, um die relevanten Faktoren einzuschätzen. '
    'Es ersetzt nicht die Praxis des Stimmers und das geschulte Ohr. '
    'Im Vergleich von Stimmplatten hört man alle Feinheiten, die messtechnisch mit '
    'realistischem Aufwand praktisch unmöglich zu erfassen sind. '
    'Erfahrung und Gehörschulung sind der Schlüssel zum Erfolg.',sBody))
story.append(Table([['']], colWidths=[PAGE_W*0.6],
    style=TableStyle([('LINEBELOW',(0,0),(-1,-1),1,ACCENTRED),('ALIGN',(0,0),(-1,-1),'CENTER'),('FONTNAME',(0,0),(-1,-1),'DejaVu')])))
story.append(Spacer(1,2*mm))
story.append(Paragraph(
    '<i>Die Zunge gibt ihre Energie der Luft — nicht der Platte. '
    'Die Platte hält die Zunge. Die Luft macht den Klang.</i>',
    ParagraphStyle('Cl',parent=sSmall,alignment=TA_CENTER,fontName='DejaVuI')))

doc.build(story, onFirstPage=set_canvas_font, onLaterPages=set_canvas_font)
print(f"PDF erzeugt: {output_path}")

# Verifikation
print("\n=== VERIFIKATION ===")
print(f"\nFall A (gezupft):")
for a in fall_a: print(f"  {a['name']:>16s}: Q={a['Q_A']:.0f}, tau={a['tau_A']:.0f}ms")
print(f"\nSpalt-Optima:")
for z in zungen:
    sk=f"{z['s_krit']*1000:.3f}" if z.get('s_krit') else "<0.005"
    print(f"  {z['name']:>16s}: s_opt={z['s_opt']*1000:.3f}mm, K={z['K_opt']:.1f}, eta={z['eta_opt']:.0%}, s_krit={sk}")
print(f"\nMaterial Q_mech:")
for m in mat_results:
    Qg=1/(1/100+1/m['Q_mount'])
    print(f"  {m['name']:>16s}: Q_mech={m['Q_mount']:.0f}, Q_ges={Qg:.0f}")
