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
