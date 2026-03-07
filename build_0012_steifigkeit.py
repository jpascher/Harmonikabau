#!/usr/bin/env python3
"""
Dok. 0012 — Zungensteifigkeit, Verstimmung und Ansprache
Build-Skript: PDF mit Diagrammen und Tabellen.
"""
import math, numpy as np, matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from reportlab.lib.pagesizes import A4; from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak
from reportlab.pdfbase import pdfmetrics; from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.fonts import addMapping

pdfmetrics.registerFont(TTFont('DejaVu','/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuB','/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuI','/usr/share/fonts/truetype/dejavu/DejaVuSans-Oblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuBI','/usr/share/fonts/truetype/dejavu/DejaVuSans-BoldOblique.ttf'))
pdfmetrics.registerFont(TTFont('DejaVuMono','/usr/share/fonts/truetype/dejavu/DejaVuSansMono.ttf'))
addMapping('DejaVu',0,0,'DejaVu');addMapping('DejaVu',1,0,'DejaVuB')
addMapping('DejaVu',0,1,'DejaVuI');addMapping('DejaVu',1,1,'DejaVuBI')

DB=HexColor('#16213e');AR=HexColor('#e94560');KG=HexColor('#2e7d32');WR=HexColor('#c62828')
LG=HexColor('#f5f5f5');KB=HexColor('#e8f5e9');WB=HexColor('#ffebee');PW=A4[0]-50*mm
styles=getSampleStyleSheet()
for sn in styles.byName:
    s=styles.byName[sn]
    if hasattr(s,'fontName'):
        if 'Bold' in s.fontName:s.fontName='DejaVuB'
        elif 'Italic' in s.fontName:s.fontName='DejaVuI'
        else:s.fontName='DejaVu'
sT=ParagraphStyle('T',parent=styles['Title'],fontSize=18,textColor=DB,spaceAfter=4,alignment=TA_CENTER,fontName='DejaVuB')
sST=ParagraphStyle('ST',parent=styles['Normal'],fontSize=12,textColor=DB,alignment=TA_CENTER,spaceAfter=2,fontName='DejaVu')
sAb=ParagraphStyle('Ab',parent=styles['Italic'],fontSize=9,textColor=HexColor('#555'),spaceAfter=8,fontName='DejaVuI')
sCh=ParagraphStyle('Ch',parent=styles['Heading1'],fontSize=14,textColor=DB,spaceBefore=14,spaceAfter=6,fontName='DejaVuB')
sB=ParagraphStyle('Bo',parent=styles['Normal'],fontSize=10,leading=14,spaceAfter=6,alignment=TA_JUSTIFY,fontName='DejaVu')
sBB=ParagraphStyle('BB',parent=sB,fontName='DejaVuB')
sFo=ParagraphStyle('Fo',parent=sB,fontSize=10,alignment=TA_CENTER,spaceAfter=8,spaceBefore=4,fontName='DejaVuMono')
sSm=ParagraphStyle('Sm',parent=sB,fontSize=9,leading=12,fontName='DejaVu')
sKB=ParagraphStyle('KB',parent=sB,backColor=KB,borderPadding=6,borderColor=KG,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sWB=ParagraphStyle('WB',parent=sB,backColor=WB,borderPadding=6,borderColor=WR,borderWidth=1,spaceAfter=8,fontName='DejaVu')
sTH=ParagraphStyle('TH',parent=sB,fontSize=9,fontName='DejaVuB',alignment=TA_CENTER,leading=11)
sTD=ParagraphStyle('TD',parent=sB,fontSize=9,alignment=TA_CENTER,leading=11,fontName='DejaVu')
sTDL=ParagraphStyle('TDL',parent=sB,fontSize=9,alignment=TA_LEFT,leading=11,fontName='DejaVu')
def hr():return Table([['']], colWidths=[PW], style=TableStyle([('LINEBELOW',(0,0),(-1,-1),2,AR),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
def mt(hdr,rows,cw=None):
    data=[[Paragraph(h,sTH)for h in hdr]]
    for row in rows:data.append([Paragraph(str(c),sTDL)if i==0 else Paragraph(str(c),sTD)for i,c in enumerate(row)])
    cw=cw or[PW/len(hdr)]*len(hdr)
    t=Table(data,colWidths=cw,repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),DB),('TEXTCOLOR',(0,0),(-1,0),white),('ROWBACKGROUNDS',(0,1),(-1,-1),[white,LG]),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),3),('BOTTOMPADDING',(0,0),(-1,-1),3),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
    return t

# ══════════════ PHYSIK ══════════════
E=200e9;rho_s=7800;rho_0=1.2;alpha_cal=0.001
zungen=[
    {'name':'Bass 50 Hz','f':50,'L':70e-3,'W':10e-3,'t':0.40e-3,'h':15e-3,'s':0.06e-3},
    {'name':'Bass A1','f':55,'L':60e-3,'W':8e-3,'t':0.35e-3,'h':14e-3,'s':0.06e-3},
    {'name':'Bass A2','f':110,'L':40e-3,'W':7e-3,'t':0.30e-3,'h':10e-3,'s':0.05e-3},
    {'name':'Mitte A3','f':220,'L':30e-3,'W':6e-3,'t':0.25e-3,'h':6e-3,'s':0.05e-3},
    {'name':'Diskant A4','f':440,'L':20e-3,'W':5e-3,'t':0.20e-3,'h':4e-3,'s':0.04e-3},
    {'name':'Hoch A5','f':880,'L':12e-3,'W':4e-3,'t':0.15e-3,'h':2.5e-3,'s':0.03e-3},
    {'name':'Sehr hoch C6','f':1047,'L':10e-3,'W':3.5e-3,'t':0.12e-3,'h':2e-3,'s':0.025e-3},
]
for z in zungen:
    I=z['W']*z['t']**3/12;z['k']=3*E*I/z['L']**3;z['m']=rho_s*z['L']*z['W']*z['t']
    z['f_calc']=1.875**2/(2*math.pi*z['L']**2)*math.sqrt(E*I/(rho_s*z['W']*z['t']))
    xg=z['s']*3;alpha=alpha_cal*z['W']*z['h']/xg;z['alpha']=alpha
    z['df_cents']=[]; z['p_th']=z['k']*z['s']/(z['W']*z['h'])*0.5
    for p in [500,1000,1500]:z['df_cents'].append(-alpha*p/z['k']*1200)

# ══════════════ DIAGRAMME ══════════════
cols=['#0D47A1','#1565C0','#1976D2','#2196F3','#F9A825','#E65100','#C62828']
# Diag 1: Verstimmung vs Druck
fig1,ax1=plt.subplots(figsize=(9,5.5));p_a=np.linspace(100,2000,200)
for i,z in enumerate(zungen):
    xg=z['s']*3;al=alpha_cal*z['W']*z['h']/xg
    ax1.plot(p_a,[-al*p/z['k']*1200 for p in p_a],lw=2,color=cols[i],label=f'{z["name"]} (k={z["k"]:.0f})')
ax1.axhline(y=0,color='black',lw=0.5);ax1.axhspan(-5,5,color='green',alpha=0.1,label='Hörgrenze ±5 Cent')
ax1.set_xlabel('Balgdruck [Pa]',fontsize=11);ax1.set_ylabel('Verstimmung [Cent]',fontsize=11)
ax1.set_title('Frequenzverschiebung durch Balgdruck\n(weiche Zungen verstimmen stärker)',fontsize=12)
ax1.legend(fontsize=8);ax1.grid(alpha=0.3);ax1.set_xlim(100,2000);plt.tight_layout()
fig1.savefig('/home/claude/diag_0012_1_verstimmung.png',dpi=150);plt.close()

# Diag 2: Trade-off
fig2,ax2=plt.subplots(figsize=(8,6))
for i,z in enumerate(zungen):
    ax2.plot(z['p_th'],abs(z['df_cents'][1]),'o',ms=12,color=cols[i],zorder=5)
    ax2.annotate(z['name'],(z['p_th'],abs(z['df_cents'][1])),textcoords="offset points",xytext=(8,5),fontsize=8)
p_t=np.linspace(10,5000,200);zr=zungen[4];C_=zr['p_th']*abs(zr['df_cents'][1])
ax2.plot(p_t,C_/p_t,'k--',lw=1,alpha=0.4,label='Δf × p_th ≈ const')
ax2.axhspan(0,5,color='green',alpha=0.08,label='Δf < 5 Cent (stabil)')
ax2.axvspan(0,300,color='blue',alpha=0.05,label='p_th < 300 Pa (gut)')
ax2.set_xlabel('Schwellendruck p_th [Pa]  (→ schlechtere Ansprache)',fontsize=11)
ax2.set_ylabel('|Verstimmung| bei 1 kPa [Cent]  (→ instabilerer Ton)',fontsize=11)
ax2.set_title('Der Kompromiss: Ansprache vs. Stimmstabilität\n(Links unten = ideal, aber physikalisch unmöglich)',fontsize=12)
ax2.legend(fontsize=8);ax2.grid(alpha=0.3);ax2.set_xlim(0,800);ax2.set_ylim(0,max(abs(z['df_cents'][1])for z in zungen)*1.2)
plt.tight_layout();fig2.savefig('/home/claude/diag_0012_2_tradeoff.png',dpi=150);plt.close()

# Diag 3: k ∝ t³
fig3,ax3=plt.subplots(figsize=(8,5));t_r=np.linspace(0.08,0.50,200)
Lr=20e-3;Wr=5e-3
k_t=[E*Wr*(t*1e-3)**3/(4*Lr**3)for t in t_r];f_t=[1.875**2/(2*math.pi*Lr**2)*math.sqrt(E*(t*1e-3)**2/(12*rho_s))for t in t_r]
ax3b=ax3.twinx()
ax3.plot(t_r,k_t,'b-',lw=2.5,label='k (∝ t³)');ax3b.plot(t_r,f_t,'r--',lw=2,label='f (∝ t)')
ax3.set_xlabel('Zungendicke t [mm]',fontsize=11);ax3.set_ylabel('k [N/m]',fontsize=11,color='blue')
ax3b.set_ylabel('f [Hz]',fontsize=11,color='red')
ax3.set_title('k ∝ t³, aber f ∝ t\n(10% dünner → 27% weicher, aber nur 10% tiefer)',fontsize=12)
ax3.axvline(x=0.20,color='green',ls='--',alpha=0.5,label='A4: t=0,20mm')
l1,la1=ax3.get_legend_handles_labels();l2,la2=ax3b.get_legend_handles_labels()
ax3.legend(l1+l2,la1+la2,fontsize=9);ax3.grid(alpha=0.3);plt.tight_layout()
fig3.savefig('/home/claude/diag_0012_3_kubisch.png',dpi=150);plt.close()

# Diag 4: Stimmgewicht vs Verdünnung
fig4,(ax4a,ax4b)=plt.subplots(1,2,figsize=(10,5))
zr=zungen[4];f_target=zr['f_calc']*0.9
t_w2=zr['t']*0.9;k_w2=E*zr['W']*t_w2**3/(4*zr['L']**3)
xg=zr['s']*3;al_r=alpha_cal*zr['W']*zr['h']/xg;p2=np.linspace(100,2000,100)
ax4a.plot(p2,[-al_r*p/zr['k']*1200 for p in p2],'b-',lw=2,label=f'Original (k={zr["k"]:.0f})')
ax4a.plot(p2,[-al_r*p/zr['k']*1200 for p in p2],'g--',lw=2,label=f'+Masse (k={zr["k"]:.0f})')
ax4a.plot(p2,[-al_r*p/k_w2*1200 for p in p2],'r-',lw=2,label=f'Dünner (k={k_w2:.0f})')
ax4a.axhspan(-5,5,color='green',alpha=0.1);ax4a.set_xlabel('Balgdruck [Pa]',fontsize=10)
ax4a.set_ylabel('Verstimmung [Cent]',fontsize=10)
ax4a.set_title('f um 10% senken:\nMasse vs. Verdünnung',fontsize=11)
ax4a.legend(fontsize=8);ax4a.grid(alpha=0.3)
cats=['Original','+Masse\n(gleich k)','Dünner\n(k −27%)']
p_ths=[zr['p_th'],zr['p_th'],zr['p_th']*k_w2/zr['k']]
drifts=[abs(al_r*1000/zr['k']*1200),abs(al_r*1000/zr['k']*1200),abs(al_r*1000/k_w2*1200)]
w=0.35;x_p=range(3)
ax4b.bar([x-w/2 for x in x_p],p_ths,w,color='#42A5F5',alpha=0.7,label='p_th [Pa]')
ax4b2=ax4b.twinx()
ax4b2.bar([x+w/2 for x in x_p],drifts,w,color='#EF5350',alpha=0.7,label='|Δf| [Cent]')
ax4b.set_xticks(x_p);ax4b.set_xticklabels(cats,fontsize=9)
ax4b.set_ylabel('p_th [Pa]',fontsize=10,color='blue');ax4b2.set_ylabel('Verstimmung [Cent]',fontsize=10,color='red')
ax4b.set_title('Vergleich',fontsize=11);ax4b.legend(loc='upper left',fontsize=8);ax4b2.legend(loc='upper right',fontsize=8)
plt.tight_layout();fig4.savefig('/home/claude/diag_0012_4_stimmgewicht.png',dpi=150);plt.close()

# Diag 5: Steifigkeits-Landkarte
fig5,ax5=plt.subplots(figsize=(9,5.5))
for i,z in enumerate(zungen):
    ax5.scatter(z['f'],z['k'],s=z['W']*1e3*50,color=cols[i],alpha=0.7,edgecolors='black',zorder=5)
    ax5.annotate(z['name'],(z['f'],z['k']),textcoords="offset points",xytext=(8,5),fontsize=8)
ax5.set_xlabel('Frequenz [Hz]',fontsize=11);ax5.set_ylabel('k [N/m]',fontsize=11)
ax5.set_title('Steifigkeit vs. Frequenz\n(Kreisfläche ∝ Zungenbreite)',fontsize=12)
ax5.grid(alpha=0.3);ax5.set_xlim(0,1200);ax5.set_ylim(0,450);plt.tight_layout()
fig5.savefig('/home/claude/diag_0012_5_landkarte.png',dpi=150);plt.close()
print("Diagramme erzeugt.")

# ══════════════ PDF ══════════════
out='/home/claude/0012_steifigkeit_De.pdf'
doc=SimpleDocTemplate(out,pagesize=A4,leftMargin=25*mm,rightMargin=25*mm,topMargin=25*mm,bottomMargin=25*mm)
def sf(c,d):c.setFont('DejaVu',10)
story=[]
story.append(Paragraph('Dok. 0012',sSm));story.append(Spacer(1,2*mm))
story.append(Paragraph('Zungensteifigkeit, Verstimmung und Ansprache',sT))
story.append(Paragraph('Warum der Stimmer einen Kompromiss sucht',sST))
story.append(Spacer(1,2*mm));story.append(hr());story.append(Spacer(1,3*mm))
story.append(Paragraph('Die Steifigkeit der Stimmzunge bestimmt zwei Eigenschaften, die einander widersprechen: '
    'Ansprache und Stimmstabilität. Eine weiche Zunge spricht leichter an, verstimmt aber stärker '
    'mit dem Balgdruck. Eine steife Zunge hält den Ton stabiler, braucht aber mehr Druck zum Anschwingen. '
    'Dieses Dokument berechnet den physikalischen Zusammenhang und zeigt, warum der Stimmer einen Kompromiss '
    'zwischen beiden sucht — und welche Werkzeuge ihm dafür zur Verfügung stehen.',sAb))
refs=[['0010','Güte der Stimmzunge: Zwei Verlustkanäle'],['0011','Kanalgeometrie der Stimmplatte']]
story.append(Paragraph('<b>Dokumentenverweise</b>',sSm))
rd=[[Paragraph(f'<b>Dok. {r[0]}</b>',sTD),Paragraph(r[1],sTDL)]for r in refs]
rt=Table(rd,colWidths=[30*mm,PW-30*mm])
rt.setStyle(TableStyle([('GRID',(0,0),(-1,-1),0.3,HexColor('#ccc')),('BACKGROUND',(0,0),(0,-1),LG),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('TOPPADDING',(0,0),(-1,-1),2),('BOTTOMPADDING',(0,0),(-1,-1),2),('FONTNAME',(0,0),(-1,-1),'DejaVu')]))
story.append(rt);story.append(Spacer(1,4*mm))

# Kap 1
story.append(Paragraph('Kapitel 1: Steifigkeit und Frequenz',sCh))
story.append(Paragraph('Die Stimmzunge ist ein Cantilever-Balken. Frequenz und Steifigkeit hängen beide von der Geometrie ab, aber <b>unterschiedlich stark</b>:',sB))
story.append(Paragraph('f ∝ t / L²         k ∝ W × t³ / L³',sFo))
story.append(Paragraph('Die Steifigkeit k hängt <b>kubisch</b> von der Dicke t ab, die Frequenz f nur linear. '
    '10 % dünner bedeutet 27 % weicher, aber nur 10 % tiefer. Das ist der Schlüssel zum Verständnis des Kompromisses.',sB))
story.append(Spacer(1,2*mm))
t_rows=[]
for z in zungen:
    t_rows.append([z['name'],f"{z['f']}",f"{z['L']*1000:.0f}",f"{z['W']*1000:.0f}",f"{z['t']*1000:.2f}",f"{z['k']:.0f}"])
story.append(mt(['Zunge','f [Hz]','L [mm]','W [mm]','t [mm]','k [N/m]'],t_rows,[28*mm,16*mm,16*mm,16*mm,16*mm,PW-92*mm]))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0012_3_kubisch.png',width=PW,height=PW*0.625))
story.append(Paragraph('<i>Abb. 1: Steifigkeit (blau) wächst kubisch mit der Dicke, Frequenz (rot) nur linear. '
    'Das bedeutet: Kleine Änderungen in t haben große Auswirkungen auf k.</i>',sSm))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0012_5_landkarte.png',width=PW,height=PW*0.61))
story.append(Paragraph('<i>Abb. 2: Steifigkeit vs. Frequenz für alle Zungengrößen. '
    'Kreisfläche proportional zur Zungenbreite.</i>',sSm))

# Kap 2
story.append(PageBreak())
story.append(Paragraph('Kapitel 2: Verstimmung durch Balgdruck',sCh))
story.append(Paragraph('Der Balgdruck erzeugt eine Strömung durch den Spalt. '
    'Der Bernoulli-Effekt erzeugt dabei einen Unterdruck, der die Zunge zum Schlitz hin zieht — '
    'eine „negative Feder", die die effektive Steifigkeit senkt:',sB))
story.append(Paragraph('k<sub>eff</sub> = k<sub>Zunge</sub> − k<sub>Bernoulli</sub>(p)',sFo))
story.append(Paragraph('Da f ∝ √k, sinkt die Frequenz mit steigendem Druck. Je weicher die Zunge '
    '(kleines k<sub>Zunge</sub>), desto stärker wirkt sich die druckabhängige Bernoulli-Feder relativ aus:',sB))
story.append(Paragraph('Δf/f ≈ −k<sub>Bernoulli</sub> / (2 × k<sub>Zunge</sub>)',sFo))
story.append(Spacer(1,2*mm))
v_rows=[]
for z in zungen:
    v_rows.append([z['name'],f"{z['k']:.0f}",f"{z['df_cents'][0]:+.1f}",f"{z['df_cents'][1]:+.1f}",f"{z['df_cents'][2]:+.1f}"])
story.append(mt(['Zunge','k [N/m]','Δf bei 500 Pa','Δf bei 1000 Pa','Δf bei 1500 Pa'],v_rows,
    [28*mm,18*mm,24*mm,26*mm,PW-96*mm]))
story.append(Paragraph('<i>Verstimmung in Cent (negative Werte = Frequenz sinkt). Hörgrenze ≈ 2–5 Cent.</i>',sSm))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0012_1_verstimmung.png',width=PW,height=PW*0.61))
story.append(Paragraph('<i>Abb. 3: Verstimmung vs. Balgdruck. Der grüne Bereich (±5 Cent) markiert die Hörgrenze. '
    'Bass-Zungen verlassen den stabilen Bereich schnell, Diskant-Zungen bleiben drin.</i>',sSm))
story.append(Paragraph(
    '<b>Ergebnis:</b> Bass-Zungen (k ≈ 80–150 N/m) verstimmen um 4–11 Cent pro kPa — '
    'hörbar bei Dynamikwechseln. Über den gesamten Spieldruckbereich (200–2000 Pa) beträgt die '
    'Drift bis zu 20 Cent — ein Viertelton (50 Cent) ist das absolute Limit. '
    'Diskant-Zungen (k ≈ 250–400 N/m) verstimmen um 0,3–0,8 Cent/kPa — '
    'kaum hörbar. Der Stimmer berücksichtigt diesen Effekt: Bass-Zungen werden für mittleren '
    'Spieldruck gestimmt, nicht für Extremwerte.',sKB))

# Kap 3
story.append(Paragraph('Kapitel 3: Der Kompromiss — Ansprache vs. Stabilität',sCh))
story.append(Paragraph('Der Schwellendruck (minimaler Balgdruck zum Anschwingen) ist proportional zur Steifigkeit:',sB))
story.append(Paragraph('p<sub>threshold</sub> ∝ k × s / (W × h)',sFo))
story.append(Paragraph('Weiche Zunge → niedriger Schwellendruck → gute Ansprache. '
    'Aber: Weiche Zunge → mehr Verstimmung. Es gibt kein Entkommen:',sB))
story.append(Paragraph('Δf × p<sub>threshold</sub> ≈ const',sFo))
story.append(Paragraph('Wer die Ansprache verbessert (p<sub>th</sub> senken), verschlechtert die Stimmstabilität (Δf steigt) '
    '— und umgekehrt. Das ist ein physikalisches Gesetz, kein Konstruktionsmangel.',sB))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0012_2_tradeoff.png',width=PW,height=PW*0.75))
story.append(Paragraph('<i>Abb. 4: Ansprache vs. Stimmstabilität. Die gestrichelte Kurve zeigt das physikalische Limit: '
    'Δf × p<sub>th</sub> ≈ const. Links unten (gute Ansprache UND stabiler Ton) ist unerreichbar.</i>',sSm))

# Kap 4
story.append(PageBreak())
story.append(Paragraph('Kapitel 4: Stimmgewicht vs. Verdünnung',sCh))
story.append(Paragraph('Der Stimmer hat zwei Werkzeuge, um die Frequenz zu senken:',sB))
story.append(Paragraph('<b>Masse hinzufügen</b> (Stimmgewicht, Lötzinn, Wachs an der Spitze): '
    'Senkt f, ohne k zu ändern. Die Verstimmung bleibt gleich (gleiches k<sub>B</sub>/k). '
    'Nachteil: Mehr Trägheit → langsameres Anschwingen.',sB))
story.append(Paragraph('<b>Zunge verdünnen</b> (Abschleifen, Feilen): '
    'Senkt f UND k. Weil k ∝ t³ kubisch sinkt, aber f ∝ t nur linear, '
    'verschlechtert sich das Verhältnis k<sub>B</sub>/k dramatisch. '
    '10 % verdünnen → 27 % weicher → 37 % mehr Verstimmung.',sB))
story.append(Spacer(1,2*mm))
story.append(Image('/home/claude/diag_0012_4_stimmgewicht.png',width=PW,height=PW*0.5))
story.append(Paragraph('<i>Abb. 5: Links: Verstimmungskurven. Original und +Masse liegen übereinander '
    '(gleiche Steifigkeit → gleiche Verstimmung). Die dünnere Zunge (rot) verstimmt stärker. '
    'Rechts: Masse ändert den Schwellendruck nicht, Verdünnung senkt ihn (bessere Ansprache, aber mehr Drift).</i>',sSm))
story.append(Paragraph(
    '<b>Konsequenz für die Praxis:</b> Frequenz durch Masse senken ist stimmstabiler. '
    'Frequenz durch Verdünnen senken verbessert die Ansprache, verschlechtert aber die Stabilität. '
    'Der erfahrene Stimmer verwendet beide Werkzeuge: Grobe Stimmung durch Geometrie (Länge, Dicke), '
    'Feinstimmung durch Stimmgewicht (Masse an der Spitze).',sKB))

# Kap 5
story.append(Paragraph('Kapitel 5: Zusammenfassung',sCh))
story.append(Paragraph('<b>1. Kubische Kopplung:</b> k ∝ t³ — die Steifigkeit reagiert 3× empfindlicher auf '
    'Dickenänderungen als die Frequenz (f ∝ t). 10 % dünner = 27 % weicher.',sKB))
story.append(Paragraph('<b>2. Verstimmung ∝ 1/k:</b> Die Bernoulli-Kraft wirkt als negative Feder. '
    'Je weicher die Zunge, desto stärker die Verstimmung. '
    'Bass: 4–11 Cent/kPa (hörbar). Diskant: 0,3–0,8 Cent/kPa (kaum hörbar). '
    'Ein Viertelton (50 Cent) ist das absolute Limit.',sKB))
story.append(Paragraph('<b>3. Physikalisches Limit:</b> Δf × p<sub>threshold</sub> ≈ const. '
    'Gute Ansprache und stabile Stimmung widersprechen sich — es gibt kein Entkommen, nur einen Kompromiss.',sKB))
story.append(Paragraph('<b>4. Werkzeuge:</b> Stimmgewicht senkt f bei konstantem k (stimmstabil). '
    'Verdünnung senkt f UND k (bessere Ansprache, aber mehr Drift). '
    'Profilierung (Verjüngung zur Spitze) erlaubt die Entkopplung von Masse und Steifigkeit.',sKB))
story.append(Paragraph('<b>5. Einordnung:</b> Die berechneten Cent-Werte sind Größenordnungsabschätzungen '
    '(kalibriert auf Praxiswerte). Die qualitative Aussage — weiche Zungen verstimmen stärker — '
    'ist physikalisch exakt. Die quantitativen Werte hängen von der Zungenaufbiegung, '
    'dem Spaltprofil und der Spieltechnik ab.',sWB))

story.append(Spacer(1,6*mm))
story.append(Table([['']], colWidths=[PW*0.6],style=TableStyle([('LINEBELOW',(0,0),(-1,-1),1,AR),('ALIGN',(0,0),(-1,-1),'CENTER'),('FONTNAME',(0,0),(-1,-1),'DejaVu')])))
story.append(Spacer(1,2*mm))
story.append(Paragraph('<i>Die Zunge sucht ihr Gleichgewicht zwischen Freiheit und Kontrolle — genau wie der Spieler.</i>',
    ParagraphStyle('Cl',parent=sSm,alignment=TA_CENTER,fontName='DejaVuI')))

doc.build(story,onFirstPage=sf,onLaterPages=sf)
print(f"PDF erzeugt: {out}")
for z in zungen:print(f"  {z['name']:>16s}: k={z['k']:.0f} N/m, Δf@1kPa={z['df_cents'][1]:+.1f}ct, p_th={z['p_th']:.0f}Pa")
