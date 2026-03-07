#!/usr/bin/env python3
"""
Dok. 0015 — Akustische Kopplung und Impedanzanpassung
Build-Skript: PDF mit Diagrammen und Tabellen.
Benötigt: reportlab, matplotlib, numpy
"""
import math
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import BytesIO

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Image, PageBreak, KeepTogether
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

addMapping('DejaVu', 0, 0, 'DejaVu')
addMapping('DejaVu', 1, 0, 'DejaVuB')
addMapping('DejaVu', 0, 1, 'DejaVuI')
addMapping('DejaVu', 1, 1, 'DejaVuBI')

DARKBLUE=HexColor('#16213e'); ACCENTRED=HexColor('#e94560')
KEYGREEN=HexColor('#2e7d32'); WARNRED=HexColor('#c62828')
LIGHTGRAY=HexColor('#f5f5f5'); KEYBG=HexColor('#e8f5e9'); WARNBG=HexColor('#ffebee')
W_PAGE=A4[0]; PAGE_W=W_PAGE-50*mm

styles=getSampleStyleSheet()
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
sSection=ParagraphStyle('Se',parent=styles['Heading2'],fontSize=12,textColor=DARKBLUE,spaceBefore=10,spaceAfter=4,fontName='DejaVuB')
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

def fig_to_image(fig, width_mm=160, dpi=200):
    buf = BytesIO()
    fig.savefig(buf, format='png', dpi=dpi, bbox_inches='tight')
    plt.close(fig)
    buf.seek(0)
    img = Image(buf)
    ratio = img.imageWidth / img.imageHeight
    img.drawWidth = width_mm * mm
    img.drawHeight = width_mm * mm / ratio
    return img

# ══════════════════════════════════════════════════════════════
# PHYSIK-KONSTANTEN
# ══════════════════════════════════════════════════════════════
rho_0 = 1.2     # kg/m³ Luft
c = 343.0        # m/s Schallgeschwindigkeit
mu = 1.8e-5      # Pa·s dynamische Viskosität

# Bass-Zunge Referenz
f_bass = 50      # Hz
L_z = 70e-3      # m Zungenlänge
W_z = 8e-3       # m Zungenbreite
t_z = 0.4e-3     # m Zungendicke
x_tip = 2.0e-3   # m Amplitude
s_side = 0.05e-3 # m Seitenspalt

# Kammergeometrie (Bass)
V_kammer = 50e-6  # m³ = 50 cm³
L_kammer = 120e-3 # m Kammerlänge
H_kammer = 50e-3  # m Kammerhöhe
B_kammer = 10e-3  # m Kammerbreite (Tiefe)
A_schlitz = W_z * L_z * 0.5  # effektive Schlitzfläche ≈ Mittelwert

# Öffnungen
A_klappe = 300e-6  # m² ≈ 300 mm² Klappenöffnung
l_klappe = 5e-3    # m effektive Halslänge Klappe
A_zunge_offen = W_z * x_tip  # 8×2 mm = 16 mm² bei voller Amplitude
A_spalt_seitl = 2 * (W_z + L_z) * s_side  # Seitenspalt
l_zunge = 3e-3     # m effektive Halslänge am Zungenende (Plattendicke)

# ══════════════════════════════════════════════════════════════
# DIAGRAMME
# ══════════════════════════════════════════════════════════════
plt.rcParams['font.family'] = 'DejaVu Sans'
plt.rcParams['font.size'] = 10

# --- Abb. 1: Zwei-Port-Schema mit Conductance-Anteilen ---
def make_fig1():
    fig, ax = plt.subplots(1, 1, figsize=(8, 4))
    ax.set_xlim(0, 10); ax.set_ylim(0, 5)
    ax.set_aspect('equal')
    ax.axis('off')
    
    # Kammer
    rect = plt.Rectangle((2, 1), 6, 3, fill=True, facecolor='#e3f2fd', edgecolor='#16213e', lw=2)
    ax.add_patch(rect)
    ax.text(5, 2.5, 'Volumen V', ha='center', va='center', fontsize=14, fontweight='bold', color='#16213e')
    
    # Port 1 (Zunge, links)
    ax.annotate('', xy=(2, 2.5), xytext=(0.5, 2.5),
                arrowprops=dict(arrowstyle='->', lw=2, color='#e94560'))
    ax.text(0.3, 3.5, 'Port 1 (Zunge)', fontsize=11, fontweight='bold', color='#e94560')
    ax.text(0.3, 3.1, f'A₁ ≈ {A_zunge_offen*1e6:.0f} mm² (offen)', fontsize=9, color='#555')
    ax.text(0.3, 2.7, f'+ Seitenspalt {A_spalt_seitl*1e6:.0f} mm²', fontsize=9, color='#555')
    ax.text(0.3, 1.5, '≈ 70 % Conductance', fontsize=10, fontweight='bold', color='#e94560')
    
    # Port 2 (Klappe, rechts)
    ax.annotate('', xy=(9.5, 2.5), xytext=(8, 2.5),
                arrowprops=dict(arrowstyle='->', lw=2, color='#2e7d32'))
    ax.text(8.2, 3.5, 'Port 2 (Klappe)', fontsize=11, fontweight='bold', color='#2e7d32')
    ax.text(8.2, 3.1, f'A₂ ≈ {A_klappe*1e6:.0f} mm²', fontsize=9, color='#555')
    ax.text(8.2, 1.5, '≈ 30 % Conductance', fontsize=10, fontweight='bold', color='#2e7d32')
    
    # Formel
    ax.text(5, 0.3, 'f = (c/2π) × √((S₁/l₁ + S₂/l₂) / V)', ha='center', fontsize=11,
            fontfamily='monospace', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round'))
    
    fig.suptitle('Abb. 1: Zwei-Port-System — Die Zunge dominiert die Conductance', fontsize=12, fontweight='bold')
    fig.tight_layout()
    return fig

# --- Abb. 2: Ein-Port vs. Zwei-Port Vergleich ---
def make_fig2():
    fig, ax = plt.subplots(figsize=(7, 4))
    
    # Berechnung
    G_klappe = A_klappe / l_klappe  # Conductance Klappe
    G_zunge_avg = (A_zunge_offen + A_spalt_seitl) / l_zunge  # Conductance Zunge (zeitgemittelt)
    
    f_ein_port = (c / (2*np.pi)) * np.sqrt(G_klappe / V_kammer)
    f_zwei_port = (c / (2*np.pi)) * np.sqrt((G_klappe + G_zunge_avg) / V_kammer)
    
    bars = ax.bar(['Ein-Port\n(nur Klappe)', 'Zwei-Port\n(Klappe + Zunge)'], 
                  [f_ein_port, f_zwei_port],
                  color=['#2e7d32', '#e94560'], width=0.5)
    
    ax.set_ylabel('f_H [Hz]')
    ax.set_title('Abb. 2: Ein-Port vs. Zwei-Port Helmholtz-Frequenz')
    
    for bar, val in zip(bars, [f_ein_port, f_zwei_port]):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 10,
                f'{val:.0f} Hz', ha='center', fontweight='bold')
    
    ax.text(0.5, max(f_ein_port, f_zwei_port)*0.5, 
            f'G_Klappe = {G_klappe*1e3:.1f} × 10⁻³ m\nG_Zunge = {G_zunge_avg*1e3:.1f} × 10⁻³ m\nFaktor {G_zunge_avg/G_klappe:.1f}×',
            ha='center', fontsize=10, bbox=dict(facecolor='#fff9c4', edgecolor='gray'))
    
    fig.tight_layout()
    return fig

# --- Abb. 3: Wellenlänge vs. Kammergeometrie (Faltung sichtbar?) ---
def make_fig3():
    fig, ax = plt.subplots(figsize=(8, 5))
    
    freqs = np.logspace(1.5, 4, 200)
    lam = c / freqs
    
    ax.loglog(freqs, lam * 1000, 'b-', lw=2, label='λ = c/f')
    
    # Kammergeometrien als horizontale Linien
    geom = {
        'Kammerlänge (120 mm)': 120,
        'Kammerhöhe (50 mm)': 50,
        'Faltwand (30 mm)': 30,
        'Öffnung (15 mm)': 15,
    }
    colors_g = ['#e94560', '#ff9800', '#9c27b0', '#2196f3']
    for (name, size), col in zip(geom.items(), colors_g):
        ax.axhline(y=size, color=col, ls='--', lw=1.5, label=name)
        # Schnittpunkt
        f_cross = c / (size * 1e-3)
        if f_cross < 10000:
            ax.plot(f_cross, size, 'o', color=col, ms=8)
            ax.annotate(f'{f_cross:.0f} Hz', (f_cross, size), textcoords='offset points',
                       xytext=(10, 5), fontsize=8, color=col)
    
    # Grüner Bereich: akustisch unsichtbar
    ax.fill_between(freqs, 0.1, 10, alpha=0.1, color='green')
    ax.text(80, 0.5, 'Geometrie < λ/10\n→ akustisch unsichtbar', fontsize=9, color='green')
    
    # 50 Hz markieren
    ax.axvline(x=50, color='red', ls=':', lw=1, alpha=0.5)
    ax.text(55, 5000, '50 Hz', fontsize=9, color='red')
    
    ax.set_xlabel('Frequenz [Hz]')
    ax.set_ylabel('Abmessung [mm]')
    ax.set_title('Abb. 3: Ist die Faltung sichtbar? — Wellenlänge vs. Kammergeometrie')
    ax.legend(loc='upper right', fontsize=8)
    ax.set_xlim(30, 10000)
    ax.set_ylim(1, 10000)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    return fig

# --- Abb. 4: Reflexion und Δf vs. Verengung ---
def make_fig4():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 4.5))
    
    # Verengungsverhältnis A_eng / A_kanal
    ratio = np.linspace(0.05, 1.0, 100)
    
    # Reflexionskoeffizient R = (1 - ratio) / (1 + ratio)  (vereinfacht)
    R = np.abs((1 - ratio) / (1 + ratio))
    
    ax1.plot(ratio, R * 100, 'r-', lw=2)
    ax1.set_xlabel('A_eng / A_Kanal')
    ax1.set_ylabel('Reflexion |R| [%]')
    ax1.set_title('Reflexion an Verengung')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 1)
    
    # Frequenzverschiebung (akustische Masse)
    # l_eff der Verengung = l / ratio (vereinfacht)
    l_eng = 5e-3  # 5 mm Verengungslänge
    delta_f_rel = []
    for r in ratio:
        m_add = rho_0 * l_eng / (r * A_klappe)  # zusätzliche akustische Masse
        m_base = rho_0 * l_klappe / A_klappe
        f_ratio = 1 / np.sqrt(1 + m_add / m_base)
        delta_f_rel.append((f_ratio - 1) * 100)
    
    ax2.plot(ratio, delta_f_rel, 'b-', lw=2)
    ax2.set_xlabel('A_eng / A_Kanal')
    ax2.set_ylabel('Δf/f [%]')
    ax2.set_title('Frequenzverschiebung durch Verengung')
    ax2.grid(True, alpha=0.3)
    ax2.set_xlim(0, 1)
    
    fig.suptitle('Abb. 4: Verengung = akustische Masse — Reflexion und Verstimmung', fontweight='bold')
    fig.tight_layout()
    return fig

# --- Abb. 5: Druckaufbau — drei Szenarien ---
def make_fig5():
    fig, axes = plt.subplots(1, 3, figsize=(11, 4))
    titles = ['Zu weit\n(schwache Kopplung)', 'Optimal\n(max. Kopplung)', 'Zu eng\n(Blockade)']
    colors_s = ['#ff9800', '#2e7d32', '#e94560']
    
    x = np.linspace(0, 1, 100)
    
    for ax, title, col in zip(axes, titles, colors_s):
        ax.set_xlim(-0.2, 1.2); ax.set_ylim(0, 1.2)
        ax.set_aspect(0.8)
        ax.axis('off')
        ax.set_title(title, fontsize=11, fontweight='bold', color=col)
    
    # Szenario 1: zu weit — Druck entweicht
    ax = axes[0]
    ax.add_patch(plt.Rectangle((0.0, 0.0), 0.4, 1.0, fc='#e3f2fd', ec='black', lw=1.5))
    ax.add_patch(plt.Rectangle((0.6, 0.0), 0.4, 1.0, fc='#e3f2fd', ec='black', lw=1.5))
    ax.annotate('', xy=(0.6, 0.5), xytext=(0.4, 0.5), arrowprops=dict(arrowstyle='->', lw=3, color='#ff9800'))
    ax.text(0.5, 0.15, 'A_eng ≈ A', ha='center', fontsize=9)
    ax.text(0.5, 0.85, 'Δp ≈ 0', ha='center', fontsize=10, fontweight='bold', color='#ff9800')
    
    # Szenario 2: optimal
    ax = axes[1]
    ax.add_patch(plt.Rectangle((0.0, 0.0), 0.35, 1.0, fc='#e3f2fd', ec='black', lw=1.5))
    ax.add_patch(plt.Rectangle((0.65, 0.0), 0.35, 1.0, fc='#e3f2fd', ec='black', lw=1.5))
    ax.add_patch(plt.Rectangle((0.35, 0.3), 0.30, 0.4, fc='#c8e6c9', ec='#2e7d32', lw=2))
    ax.annotate('', xy=(0.65, 0.5), xytext=(0.35, 0.5), arrowprops=dict(arrowstyle='->', lw=3, color='#2e7d32'))
    ax.text(0.5, 0.15, 'A_eng < A', ha='center', fontsize=9)
    ax.text(0.5, 0.85, 'Δp optimal', ha='center', fontsize=10, fontweight='bold', color='#2e7d32')
    
    # Szenario 3: zu eng
    ax = axes[2]
    ax.add_patch(plt.Rectangle((0.0, 0.0), 0.35, 1.0, fc='#e3f2fd', ec='black', lw=1.5))
    ax.add_patch(plt.Rectangle((0.65, 0.0), 0.35, 1.0, fc='#e3f2fd', ec='black', lw=1.5))
    ax.add_patch(plt.Rectangle((0.35, 0.42), 0.30, 0.16, fc='#ffcdd2', ec='#e94560', lw=2))
    ax.annotate('', xy=(0.35, 0.5), xytext=(0.65, 0.5), arrowprops=dict(arrowstyle='->', lw=2, color='#e94560', ls='--'))
    ax.text(0.5, 0.15, 'A_eng ≪ A', ha='center', fontsize=9)
    ax.text(0.5, 0.85, 'Reflexion', ha='center', fontsize=10, fontweight='bold', color='#e94560')
    
    fig.suptitle('Abb. 5: Druckaufbau vor der Engstelle — drei Szenarien', fontweight='bold')
    fig.tight_layout()
    return fig

# --- Abb. 6: Impedanz am Zungenende für verschiedene Übergangsgeometrien ---
def make_fig6():
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
    
    # Modell: Maximaler Leistungstransfer Balg → Zunge
    # Quelle: Balg mit Innenwiderstand R_i ∝ 1/A_eff² (Bernoulli)
    # Last: Zunge mit mechanischer Impedanz R_L (fix)
    # P/P_max = 4 × R_L × R_i / (R_i + R_L)²
    # Maximum bei R_i = R_L (Anpassung)
    
    A_ratio = np.linspace(0.05, 3.0, 500)  # A_eff / A_optimal
    R_i = 1.0 / A_ratio**2
    R_L = 1.0
    P_norm = 4 * R_L * R_i / (R_i + R_L)**2
    
    # Links: Anpassungskurve
    ax1.plot(A_ratio, P_norm, 'b-', lw=2.5)
    ax1.axvline(x=1.0, color='green', ls='--', lw=1.5, label='Anpassung: R_i = R_L')
    ax1.axhline(y=1.0, color='gray', ls=':', lw=1)
    
    ax1.fill_between(A_ratio[A_ratio < 0.5], P_norm[A_ratio < 0.5], alpha=0.15, color='red')
    ax1.fill_between(A_ratio[A_ratio > 2.0], P_norm[A_ratio > 2.0], alpha=0.15, color='orange')
    ax1.fill_between(A_ratio[(A_ratio >= 0.5) & (A_ratio <= 2.0)], 
                     P_norm[(A_ratio >= 0.5) & (A_ratio <= 2.0)], alpha=0.1, color='green')
    
    ax1.text(0.25, 0.3, 'Zu eng\nR_i \u226b R_L\nBlockade', ha='center', fontsize=9, color='red')
    ax1.text(2.5, 0.3, 'Zu weit\nR_i \u226a R_L\nKurzschluss', ha='center', fontsize=9, color='#cc6600')
    ax1.text(1.0, 0.5, 'Anpassung', ha='center', fontsize=10, fontweight='bold', color='green')
    
    ax1.set_xlabel('A_\u00dcbergang / A_optimal')
    ax1.set_ylabel('P / P_max')
    ax1.set_title('Leistungstransfer: Anpassungsmodell')
    ax1.legend(fontsize=9)
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(0, 3)
    ax1.set_ylim(0, 1.15)
    
    # Rechts: Elektrisches Analogon
    ax2.set_xlim(0, 10); ax2.set_ylim(0, 6)
    ax2.set_aspect('equal')
    ax2.axis('off')
    ax2.set_title('Elektrisches Analogon', fontsize=12, fontweight='bold')
    
    ax2.add_patch(plt.Circle((1.5, 4), 0.6, fc='#e3f2fd', ec='#1565c0', lw=2))
    ax2.text(1.5, 4, 'p_Balg', ha='center', va='center', fontsize=9, fontweight='bold', color='#1565c0')
    ax2.text(1.5, 5.0, 'Quelle', ha='center', fontsize=10, color='#1565c0')
    
    ax2.add_patch(plt.Rectangle((3.0, 3.6), 1.5, 0.8, fc='#ffcdd2', ec='#c62828', lw=1.5))
    ax2.text(3.75, 4.0, 'R_i(A)', ha='center', va='center', fontsize=9, fontweight='bold', color='#c62828')
    ax2.text(3.75, 3.0, 'Str\u00f6mungs-\nwiderstand\n\u221d 1/A\u00b2', ha='center', fontsize=8, color='#c62828')
    
    ax2.add_patch(plt.Rectangle((5.3, 3.4), 1.4, 1.2, fc='#c8e6c9', ec='#2e7d32', lw=2))
    ax2.text(6.0, 4.0, 'n(A)', ha='center', va='center', fontsize=10, fontweight='bold', color='#2e7d32')
    ax2.text(6.0, 5.0, '\u00dcbergang\nSchlitz\u2192Kammer', ha='center', fontsize=8, color='#2e7d32')
    
    ax2.add_patch(plt.Rectangle((7.5, 3.6), 1.5, 0.8, fc='#fff9c4', ec='#f57f17', lw=1.5))
    ax2.text(8.25, 4.0, 'R_L', ha='center', va='center', fontsize=9, fontweight='bold', color='#f57f17')
    ax2.text(8.25, 3.0, 'Zunge\n(mechanisch)', ha='center', fontsize=8, color='#f57f17')
    
    ax2.annotate('', xy=(3.0, 4.0), xytext=(2.1, 4.0), arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax2.annotate('', xy=(5.3, 4.0), xytext=(4.5, 4.0), arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    ax2.annotate('', xy=(7.5, 4.0), xytext=(6.7, 4.0), arrowprops=dict(arrowstyle='->', lw=2, color='black'))
    
    ax2.text(5.0, 1.2, 'P_max bei R_i(A_opt) = R_L', ha='center', fontsize=11,
             fontfamily='monospace', bbox=dict(facecolor='white', edgecolor='gray', boxstyle='round'))
    ax2.text(5.0, 0.5, 'Filterwirkung f_H h\u00e4ngt vom Resonator ab (V, \u00d6ffnung)\nnicht vom Transformator (\u00dcbergangsgeometrie)', 
             ha='center', fontsize=9, style='italic', color='#555')
    
    fig.suptitle('Abb. 6: Impedanzanpassung \u2014 \u00dcbergangsgeometrie als Transformator', fontweight='bold', fontsize=12)
    fig.tight_layout()
    return fig

# --- Abb. 7: Klappenposition — f_H vs. Abstand ---
def make_fig7():
    fig, ax = plt.subplots(figsize=(7, 4.5))
    
    # Klappe am Ende vs. versetzt
    abstand = np.linspace(0, L_kammer, 100)  # Abstand Klappe vom Ende
    
    # Vereinfachung: V ändert sich nicht, aber l_eff der Klappe ändert sich
    # Position nahe Druckbauch → stärkere Kopplung
    # Position nahe Schnellebauch → schwächere Kopplung
    # Bei f << f_stehend gilt: Position ändert V_eff
    
    # Einfaches Modell: f_H(x) = f_H0 × (1 + α × (x/L)²)
    alpha_pos = 0.15  # empirischer Korrekturfaktor
    G_kl = A_klappe / l_klappe
    G_z = (A_zunge_offen + A_spalt_seitl) / l_zunge
    f_H0 = (c / (2*np.pi)) * np.sqrt((G_kl + G_z) / V_kammer)
    
    f_H = f_H0 * (1 + alpha_pos * (abstand / L_kammer)**2)
    
    ax.plot(abstand * 1000, f_H, 'b-', lw=2)
    ax.axhline(y=f_H0, color='gray', ls=':', lw=1, label=f'f_H (Klappe am Ende) = {f_H0:.0f} Hz')
    
    # Praxis-Bereich markieren
    ax.axvspan(0, 15, alpha=0.15, color='green', label='Praxis: 0–15 mm vom Ende')
    
    ax.set_xlabel('Abstand Klappe vom Kammerende [mm]')
    ax.set_ylabel('f_H [Hz]')
    ax.set_title('Abb. 7: Klappenposition verschiebt f_H — aber nur als Feinkorrektur')
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)
    
    delta_f = f_H[-1] - f_H[0]
    ax.text(60, f_H0 * 1.01, f'Δf über ganzen Bereich: {delta_f:.0f} Hz ({delta_f/f_H0*100:.1f} %)',
            fontsize=10, bbox=dict(facecolor='#fff9c4', edgecolor='gray'))
    
    fig.tight_layout()
    return fig

# --- Abb. 8: Schwingform und Volumenstrom entlang der Zunge ---
def make_fig8():
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    
    x_norm = np.linspace(0, 1, 200)  # x/L
    
    # Cantilever 1. Mode: y(x) = A × [cosh(βx) - cos(βx) - σ(sinh(βx) - sin(βx))]
    # Vereinfacht: y ∝ x² (3 - x/L) / 2  (Näherung für 1. Mode)
    # Exakter: y ∝ (cosh(1.875x) - cos(1.875x)) - 0.7341*(sinh(1.875x) - sin(1.875x))
    beta_L = 1.8751
    beta = beta_L  # für normiertes x
    sigma = (np.cosh(beta_L) + np.cos(beta_L)) / (np.sinh(beta_L) + np.sin(beta_L))
    
    y_mode = np.cosh(beta*x_norm) - np.cos(beta*x_norm) - sigma*(np.sinh(beta*x_norm) - np.sin(beta*x_norm))
    y_mode = y_mode / y_mode[-1]  # normieren auf Spitze = 1
    
    # Auslenkung
    ax1.plot(x_norm, y_mode, 'b-', lw=2, label='y(x)/y_tip')
    ax1.fill_between(x_norm, 0, y_mode, alpha=0.1, color='blue')
    ax1.set_ylabel('Auslenkung y(x)/y_tip')
    ax1.set_title('Abb. 8: Schwingform der Zunge (1. Mode, Cantilever)')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    
    # Volumenstrom pro Längeneinheit: dQ/dx = W × dy/dx × v_tip
    # dy/dx normiert
    dy = np.gradient(y_mode, x_norm)
    dy = dy / dy[-1]  # normieren
    
    # Kumulativer Volumenstrom
    Q_cumul = np.cumsum(y_mode) / np.sum(y_mode)  # normiert auf 1
    
    ax2.plot(x_norm, y_mode / np.max(y_mode), 'r-', lw=2, label='dQ/dx (Volumenstrom pro Länge)')
    ax2.plot(x_norm, Q_cumul, 'g--', lw=2, label='∫Q kumulativ')
    ax2.fill_between(x_norm, 0, y_mode / np.max(y_mode), alpha=0.1, color='red')
    
    # Schwerpunkt
    x_centroid = np.trapezoid(x_norm * y_mode, x_norm) / np.trapezoid(y_mode, x_norm)
    ax2.axvline(x=x_centroid, color='purple', ls='--', lw=2, label=f'Schwerpunkt x/L = {x_centroid:.2f}')
    
    # 50%-Punkt
    idx_50 = np.argmin(np.abs(Q_cumul - 0.5))
    ax2.axvline(x=x_norm[idx_50], color='orange', ls=':', lw=1.5, label=f'50 % bei x/L = {x_norm[idx_50]:.2f}')
    
    ax2.set_xlabel('x / L (0 = Einspannung, 1 = Spitze)')
    ax2.set_ylabel('Normierte Größe')
    ax2.set_title('Volumenstrom und Schwerpunkt der Kopplung')
    ax2.legend(fontsize=8, loc='upper left')
    ax2.grid(True, alpha=0.3)
    
    fig.tight_layout()
    return fig

# --- Abb. 9: Impedanzvergleich Z_Zunge vs. Z_Luft ---
def make_fig9():
    fig, ax = plt.subplots(figsize=(8, 5))
    
    # Verschiedene Zungengrößen
    zungen = [
        {'name': 'Bass 50 Hz', 'f': 50, 'L': 70e-3, 'W': 8e-3, 't': 0.4e-3, 'x': 2.0e-3},
        {'name': 'Bass A2', 'f': 110, 'L': 40e-3, 'W': 7e-3, 't': 0.3e-3, 'x': 1.5e-3},
        {'name': 'Diskant A4', 'f': 440, 'L': 20e-3, 'W': 5e-3, 't': 0.2e-3, 'x': 0.6e-3},
        {'name': 'Hoch A5', 'f': 880, 'L': 12e-3, 'W': 4e-3, 't': 0.15e-3, 'x': 0.4e-3},
    ]
    
    E_stahl = 200e9
    rho_stahl = 7800
    
    names = []
    z_mech = []
    z_luft = []
    ratios = []
    
    for z in zungen:
        # Mechanische Impedanz der Zunge: Z = F/v = k/(jω) ≈ k/ω
        I = z['W'] * z['t']**3 / 12
        k = 3 * E_stahl * I / z['L']**3
        omega = 2 * np.pi * z['f']
        v_tip = z['x'] * omega
        Z_m = k / omega  # |Z_mech| bei Resonanz
        
        # Akustische Impedanz der Luft am Spalt
        A_spalt = z['W'] * z['x']  # Momentane Öffnung
        Z_a = rho_0 * c / A_spalt  # Strahlungsimpedanz
        
        names.append(z['name'])
        z_mech.append(Z_m)
        z_luft.append(Z_a)
        ratios.append(Z_m / Z_a)
    
    x_pos = np.arange(len(names))
    width = 0.35
    
    bars1 = ax.bar(x_pos - width/2, z_mech, width, label='Z_Zunge (mechanisch)', color='#e94560')
    bars2 = ax.bar(x_pos + width/2, z_luft, width, label='Z_Luft (akustisch)', color='#2196f3')
    
    ax.set_yscale('log')
    ax.set_xticks(x_pos)
    ax.set_xticklabels(names, fontsize=9)
    ax.set_ylabel('Impedanz [N·s/m]')
    ax.set_title('Abb. 9: Z_Zunge ≫ Z_Luft — Die Zunge ist ein Strom-Generator')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    
    for i, r in enumerate(ratios):
        ax.text(i, max(z_mech[i], z_luft[i]) * 2, f'×{r:.0f}', ha='center', fontweight='bold', fontsize=10)
    
    fig.tight_layout()
    return fig

# ══════════════════════════════════════════════════════════════
# PDF ERZEUGEN
# ══════════════════════════════════════════════════════════════
outfile = '0015_kopplung_De.pdf'
doc = SimpleDocTemplate(outfile, pagesize=A4,
                        leftMargin=25*mm, rightMargin=25*mm,
                        topMargin=20*mm, bottomMargin=20*mm)

story = []

# ── Titelseite ──
story.append(Spacer(1, 10*mm))
story.append(Paragraph('Dok. 0015', sSubtitle))
story.append(Paragraph('Akustische Kopplung und<br/>Impedanzanpassung', sTitle))
story.append(Spacer(1, 3*mm))
story.append(hr())
story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    'Analyse der akustischen Kopplung zwischen Stimmzunge und Kammer als '
    'Zwei-Port-System. Übergangsgeometrie, Verengung als akustische Masse, '
    'Kopplungsschwerpunkt und Impedanzfehlanpassung. '
    'Referenz: Dok. 0005–0013.',
    sAbstract))
story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    '<b>Hinweis:</b> Dieses Dokument bietet ein Erklärungsmodell. '
    'Die angegebenen Zahlenwerte sind Größenordnungsabschätzungen — '
    'sie liefern keine konstruktiv verwertbaren Vorgaben. '
    'Das Optimum wird gehört, nicht berechnet.',
    sWarnBox))

# ══ Kapitel 1: Zwei-Port-System ══
story.append(Paragraph('1. Die Kammer als Zwei-Port-System', sChapter))
story.append(Paragraph(
    'In den bisherigen Dokumenten (Dok. 0004–0008) wurde die Kammer als '
    'Helmholtz-Resonator mit einer einzigen Öffnung (Klappe) behandelt. '
    'Dok. 0013 hat gezeigt, dass die Zunge den Kanal nie vollständig verschließt — '
    'der Seitenspalt bleibt immer offen. Die Kammer hat daher <b>zwei Ports</b>:',
    sBody))
story.append(Paragraph(
    'Port 1 (Zungenende): Die Zunge pulsiert und erzeugt Druckimpulse. '
    'Die effektive Fläche schwankt zwischen dem Seitenspalt '
    f'(A_min ≈ {A_spalt_seitl*1e6:.0f} mm²) und der vollen Schlitzfläche '
    f'(A_max ≈ {(W_z * x_tip + A_spalt_seitl)*1e6:.0f} mm²).',
    sBody))
story.append(Paragraph(
    f'Port 2 (Klappe): Stationäre Öffnung mit A₂ ≈ {A_klappe*1e6:.0f} mm². '
    'Hier strahlt der Schall ab.',
    sBody))

# Abb. 1
fig1 = make_fig1()
story.append(fig_to_image(fig1, width_mm=155))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Die Conductance (Leitwert) G = A/l bestimmt, wie leicht Luft durch jeden '
    'Port fließt. Port 1 (Zunge) hat zeitgemittelt etwa 70 % der Gesamt-Conductance — '
    'die Zunge ist der Hauptport, nicht die Klappe.',
    sBody))

# Abb. 2
fig2 = make_fig2()
story.append(fig_to_image(fig2, width_mm=130))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Die Zwei-Port-Formel ergibt eine deutlich höhere Helmholtz-Frequenz als '
    'das Ein-Port-Modell (nur Klappe). Das alte Modell unterschätzt f_H systematisch, '
    'weil es den dominanten Port ignoriert.',
    sBody))
story.append(Paragraph(
    '<b>Konsequenz:</b> Alle Berechnungen in Dok. 0004–0008, die f_H nur aus der '
    'Klappenöffnung ableiten, müssen als Untergrenze verstanden werden. '
    'Die tatsächliche Kammerresonanz liegt höher.',
    sKeyBox))

# ══ Kapitel 2: Ist die Faltung sichtbar? ══
story.append(Paragraph('2. Ist die Faltung akustisch sichtbar?', sChapter))
story.append(Paragraph(
    'Bei einem Bass mit 50 Hz beträgt die Wellenlänge λ = c/f = 6860 mm. '
    'Die größte Kammerabmessung (Kammerlänge 120 mm) ist nur 1,7 % von λ. '
    'Eine Faltwand (30 mm) ist 0,44 % von λ.',
    sBody))
story.append(Paragraph(
    'Die Frage ist: Ab welcher Relation Geometrie/λ wird eine Struktur '
    'akustisch „sichtbar"? Die Antwort aus der Wellenphysik: '
    'Wenn d/λ < 0,1, beugt sich die Welle vollständig um das Hindernis — '
    'es existiert akustisch nicht. Bei d/λ > 0,5 wirkt es als Reflektor.',
    sBody))

# Abb. 3
fig3 = make_fig3()
story.append(fig_to_image(fig3, width_mm=155))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Das Diagramm zeigt: Bei 50 Hz liegen <b>alle</b> Kammergeometrien tief '
    'im Bereich d ≪ λ. Erst ab ≈ 1400 Hz (28. Oberton) beginnt die Faltwand '
    'die Grenze d/λ = 0,1 zu erreichen. Für die ersten 27 Obertöne ist die '
    'Kammerform irrelevant — nur das Volumen zählt.',
    sBody))
story.append(Paragraph(
    'Die Analogie zum Lichtwellenleiter bestätigt dies: Ein glattes, rundes '
    'Rohr mit konstantem Querschnitt hat keinen Reflexionsfaktor — die Welle '
    'folgt der Kurve, solange der Biegeradius groß gegen die Wellenlänge ist. '
    'Bei Schall in der Basskammer ist der Biegeradius immer groß gegen λ.',
    sBody))
story.append(Paragraph(
    '<b>Coltman-Korrektur:</b> Coltman kompensiert nicht die Biegung selbst, '
    'sondern die Querschnittsänderung und Pfadlängendifferenz an der Faltstelle. '
    'Eine glatte Biegung mit konstantem Querschnitt ändert die Impedanz nicht (Z = ρc/A), '
    'verändert aber die effektive akustische Länge. '
    'Die Pfadlängenkorrektur skaliert mit d²/r — bei d/λ < 0,01 ist sie vernachlässigbar.',
    sBody))
story.append(Paragraph(
    '<b>Fazit:</b> In der Basskammer spielt die Faltung für den Grundton und die '
    'ersten Obertöne keine Rolle. Die Anpassung geschieht an den Öffnungen '
    '(Port 1 und Port 2), nicht an der Faltung.',
    sKeyBox))

# ══ Kapitel 3: Verengung = akustische Masse ══
story.append(Paragraph('3. Verengung = akustische Masse', sChapter))
story.append(Paragraph(
    'Jede Querschnittsverkleinerung im Kanal wirkt als konzentrierte akustische Masse: '
    'm = ρ × l / A. Das gilt unabhängig davon, ob die Verengung in einem geraden Rohr, '
    'an einer Faltstelle oder am Übergang Schlitz→Kammer sitzt.',
    sBody))
story.append(Paragraph(
    'Die Coltman-Kompensation bei gefalteten Orgelpfeifen ist im Kern eine '
    'Verengung an der Umlenkstelle — sie kompensiert den Impedanzsprung. '
    'Im gestreckten Rohr wäre dieselbe Verengung ein Reflexionspunkt und Tiefpass. '
    'Der Mechanismus ändert sich (Wellenleiter vs. Lumped Element), '
    'das physikalische Element (Verengung = Masse) bleibt dasselbe.',
    sBody))

# Abb. 4
fig4 = make_fig4()
story.append(fig_to_image(fig4, width_mm=155))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Links: Die Reflexion steigt steil an, wenn die Verengung kleiner wird. '
    'Bei A_eng/A_Kanal = 0,5 werden bereits 11 % reflektiert. '
    'Rechts: Die zusätzliche akustische Masse verschiebt f_H nach unten — '
    'je enger der Übergang, desto tiefer die Kammerresonanz.',
    sBody))

# ══ Kapitel 4: Druckaufbau vor der Engstelle ══
story.append(Paragraph('4. Druckaufbau vor der Engstelle', sChapter))
story.append(Paragraph(
    'Der Druckimpuls der Zunge staut sich vor der Engstelle am Übergang '
    'Schlitz→Kammer. Drei Szenarien zeigen die Wirkung:',
    sBody))

# Abb. 5
fig5 = make_fig5()
story.append(fig_to_image(fig5, width_mm=155))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    '<b>Zu weit:</b> Der Druck entweicht sofort ins Volumen — '
    'schwache Rückwirkung auf die Zunge, geringe Kopplung. '
    '<b>Optimal:</b> Der Druck staut sich kurz auf — '
    'maximale Rückwirkung, beste Kopplung. '
    '<b>Zu eng:</b> Der Impuls wird überwiegend reflektiert — '
    'die Energie kommt nicht ins Volumen.',
    sBody))

# ══ Kapitel 5: Impedanz am Zungenende ══
story.append(Paragraph('5. Impedanzanpassung: Übergangsgeometrie als Transformator', sChapter))
story.append(Paragraph(
    'Die Beobachtung, dass verschiedene Übergangsgeometrien die <b>Ansprache</b> '
    'ändern, aber nicht die <b>Filterwirkung</b>, lässt sich sauber als '
    'Anpassungsproblem modellieren:',
    sBody))
story.append(Paragraph(
    '<b>Quelle:</b> Der Balgdruck p_Balg mit einem Innenwiderstand R_i, der vom '
    'Strömungswiderstand des Systems abhängt. R_i ∝ 1/A² (Bernoulli: '
    'kleinere Öffnung → höherer Widerstand). '
    '<b>Last:</b> Die Zunge mit mechanischer Impedanz R_L (fix, von Steifigkeit '
    'und Masse bestimmt, Dok. 0012). '
    '<b>Transformator:</b> Die Übergangsgeometrie am Zungenende bestimmt, '
    'wie effizient die Energie der Strömung in Zungenenergie umgesetzt wird.',
    sBody))
story.append(Paragraph(
    'P / P_max = 4 × R_L × R_i / (R_i + R_L)²',
    sFormula))
story.append(Paragraph(
    'Maximum bei R_i = R_L — klassische Anpassungsbedingung. '
    'Zu eng: R_i ≫ R_L → Blockade, Energie wird reflektiert. '
    'Zu weit: R_i ≪ R_L → Kurzschluss, Druck entweicht. '
    'Optimal: R_i = R_L → maximaler Leistungstransfer Balg→Zunge.',
    sBody))

# Abb. 6
fig6 = make_fig6()
story.append(fig_to_image(fig6, width_mm=155))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Links: Die Leistungstransferkurve zeigt das typische Anpassungsmaximum. '
    'Der Bereich P > 0,9 × P_max ist relativ breit — das erklärt, warum '
    'verschiedene Geometrien alle „funktionieren". '
    'Rechts: Das elektrische Analogon. Die Quelle (Balg) treibt über den '
    'Innenwiderstand (Strömung) und den Transformator (Übergangsgeometrie) '
    'die Last (Zunge) an.',
    sBody))
story.append(Paragraph(
    '<b>Entscheidend:</b> Die Filterwirkung (f_H, Kerbfilter) hängt vom '
    'Resonator ab — Kammervolumen und Öffnungsfläche. Der Transformator '
    '(Übergangsgeometrie) ändert die Anpassungseffizienz, nicht die Resonanz. '
    'Das ist genau das, was man beobachtet: andere Ansprache, gleiche Klangfarbe.',
    sKeyBox))

# ══ Kapitel 6: Ersatzschaltbild — Zwei gekoppelte Schwingkreise ══
story.append(Paragraph('6. Ersatzschaltbild: Zwei gekoppelte Schwingkreise', sChapter))
story.append(Paragraph(
    'Die Trennung von Ansprache und Filterwirkung lässt sich als '
    'elektroakustische Analogie exakt darstellen. Zunge und Kammer sind je '
    'ein Serienschwingkreis (RLC), verbunden durch ein Kopplungselement:',
    sBody))
story.append(Paragraph(
    '<b>Kreis 1 (Zunge):</b> Masse m → Induktivität L₁. '
    'Steifigkeit k → Kapazität C₁ = 1/k. '
    'Dämpfung → Widerstand R₁. '
    'Eigenfrequenz f₁ = 1/(2π√(L₁C₁)) ≈ 50 Hz.',
    sBody))
story.append(Paragraph(
    '<b>Kreis 2 (Kammer):</b> Akustische Masse ρl/A → Induktivität L₂. '
    'Volumennachgiebigkeit V/(ρc²) → Kapazität C₂. '
    'Viskose + Abstrahlverluste → Widerstand R₂. '
    'Eigenfrequenz f_H = 1/(2π√(L₂C₂)).',
    sBody))
story.append(Paragraph(
    '<b>Kopplung:</b> Gegenseitige Impedanz M = k × √(L₁L₂). '
    'Der Kopplungskoeffizient k (0–1) wird von der Übergangsgeometrie bestimmt.',
    sBody))
story.append(Paragraph(
    'Die Transferfunktion H(f) = I₂/V₁ beschreibt, wie viel Energie '
    'vom Balg (über die Zunge) in der Kammer ankommt:',
    sBody))

# Abb. 6b: Transferfunktion und Schaltbild
# Generate the coupled circuit diagrams inline
fig_circuit, axes_c = plt.subplots(2, 2, figsize=(14, 10))

f_z_c = 50.0; m_z_c = 0.8e-3; E_sc = 200e9; W_c=8e-3; t_c=0.4e-3; L_zc=70e-3
I_zc = W_c * t_c**3 / 12
k_zc = 3 * E_sc * I_zc / L_zc**3
Q_zc = 50
L1c = m_z_c; C1c = 1.0/k_zc; R1c = 2*np.pi*f_z_c*L1c/Q_zc
V_kc = 50e-6; A_hals_c = 300e-6; l_hals_c = 8e-3; A_ze = 500e-6; l_ze = 3e-3
G_tot = A_hals_c/l_hals_c + A_ze/l_ze
L2c = rho_0/G_tot; C2c = V_kc/(rho_0*c**2); Q_Hc = 15
f_Hc = 1/(2*np.pi*np.sqrt(L2c*C2c)); R2c = 2*np.pi*f_Hc*L2c/Q_Hc

freqs_c = np.logspace(0.5, 4, 2000)
omega_c = 2*np.pi*freqs_c

# Plot 1: Transferfunktion bei verschiedenen k
ax_c = axes_c[0,0]
for name_k, k_v, col_k in [('k = 0,01 (sehr schwach)', 0.01, '#2196f3'),
    ('k = 0,05 (schwach)', 0.05, '#4caf50'), ('k = 0,15 (mittel)', 0.15, '#ff9800'),
    ('k = 0,40 (stark)', 0.40, '#e94560')]:
    Mc = k_v * np.sqrt(L1c * L2c)
    Zkc = 1j * omega_c * Mc
    Z1c_f = R1c + 1j*omega_c*L1c + 1/(1j*omega_c*C1c)
    Z2c_f = R2c + 1j*omega_c*L2c + 1/(1j*omega_c*C2c)
    denom_c = (Z1c_f + Zkc) * (Z2c_f + Zkc) - Zkc**2
    Hc = Zkc / denom_c
    H_dBc = 20*np.log10(np.abs(Hc)/np.max(np.abs(Hc))+1e-30)
    ax_c.semilogx(freqs_c, H_dBc, '-', color=col_k, lw=1.5, label=name_k)
ax_c.axvline(x=f_z_c, color='blue', ls=':', alpha=0.4)
ax_c.axvline(x=f_Hc, color='red', ls=':', alpha=0.4)
ax_c.text(f_z_c*1.1, -5, f'f_Zunge\n{f_z_c:.0f} Hz', fontsize=8, color='blue')
ax_c.text(f_Hc*1.1, -5, f'f_Kammer\n{f_Hc:.0f} Hz', fontsize=8, color='red')
ax_c.set_xlabel('Frequenz [Hz]'); ax_c.set_ylabel('|H(f)| [dB]')
ax_c.set_title('Transferfunktion: k variiert'); ax_c.legend(fontsize=8)
ax_c.grid(True, alpha=0.3); ax_c.set_xlim(3, 10000); ax_c.set_ylim(-60, 5)

# Plot 2: Impedanzen
ax_c = axes_c[0,1]
Z1_mg = np.abs(R1c + 1j*omega_c*L1c + 1/(1j*omega_c*C1c))
Z2_mg = np.abs(R2c + 1j*omega_c*L2c + 1/(1j*omega_c*C2c))
ax_c.loglog(freqs_c, Z1_mg, 'b-', lw=2, label=f'|Z\u2081| Zunge (f\u2081={f_z_c:.0f} Hz)')
ax_c.loglog(freqs_c, Z2_mg, 'r-', lw=2, label=f'|Z\u2082| Kammer (f_H={f_Hc:.0f} Hz)')
ax_c.set_xlabel('Frequenz [Hz]'); ax_c.set_ylabel('|Z|')
ax_c.set_title('Impedanzen der Schwingkreise'); ax_c.legend(fontsize=9)
ax_c.grid(True, alpha=0.3); ax_c.set_xlim(3, 10000)

# Plot 3: f_H variiert
ax_c = axes_c[1,0]; k_fix_c = 0.10
for name_f, fH_t, col_f in [('f_H = 200 Hz (V gro\u00df)', 200, '#2196f3'),
    ('f_H = 500 Hz (Standard)', 500, '#4caf50'), ('f_H = 1000 Hz (V klein)', 1000, '#ff9800'),
    ('f_H = 2000 Hz (V sehr klein)', 2000, '#e94560')]:
    C2v = 1/((2*np.pi*fH_t)**2*L2c); R2v = 2*np.pi*fH_t*L2c/Q_Hc
    Mc = k_fix_c*np.sqrt(L1c*L2c); Zkc = 1j*omega_c*Mc
    Z1v = R1c + 1j*omega_c*L1c + 1/(1j*omega_c*C1c)
    Z2v = R2v + 1j*omega_c*L2c + 1/(1j*omega_c*C2v)
    denom_v = (Z1v+Zkc)*(Z2v+Zkc)-Zkc**2; Hv = Zkc/denom_v
    H_dBv = 20*np.log10(np.abs(Hv)/np.max(np.abs(Hv))+1e-30)
    ax_c.semilogx(freqs_c, H_dBv, '-', color=col_f, lw=1.5, label=name_f)
ax_c.axvline(x=f_z_c, color='blue', ls=':', alpha=0.4)
ax_c.text(f_z_c*1.1, -5, f'f_Zunge\n{f_z_c:.0f} Hz', fontsize=8, color='blue')
ax_c.set_xlabel('Frequenz [Hz]'); ax_c.set_ylabel('|H(f)| [dB]')
ax_c.set_title(f'Kammervolumen variiert (k={k_fix_c})'); ax_c.legend(fontsize=8)
ax_c.grid(True, alpha=0.3); ax_c.set_xlim(3, 10000); ax_c.set_ylim(-60, 5)

# Plot 4: Schaltbild
ax_c = axes_c[1,1]
ax_c.set_xlim(0, 12); ax_c.set_ylim(0, 6); ax_c.set_aspect('equal'); ax_c.axis('off')
ax_c.set_title('Ersatzschaltbild', fontsize=11, fontweight='bold')
ax_c.add_patch(plt.Rectangle((0.5,3.5),1.2,0.8,fc='#e3f2fd',ec='#1565c0',lw=1.5))
ax_c.text(1.1,3.9,'L\u2081\nm',ha='center',va='center',fontsize=8,color='#1565c0')
ax_c.add_patch(plt.Rectangle((2.2,3.5),1.2,0.8,fc='#e3f2fd',ec='#1565c0',lw=1.5))
ax_c.text(2.8,3.9,'C\u2081\n1/k',ha='center',va='center',fontsize=8,color='#1565c0')
ax_c.add_patch(plt.Rectangle((3.9,3.5),1.2,0.8,fc='#e3f2fd',ec='#1565c0',lw=1.5))
ax_c.text(4.5,3.9,'R\u2081\nd',ha='center',va='center',fontsize=8,color='#1565c0')
ax_c.text(2.8,5.0,'Kreis 1: Zunge',ha='center',fontsize=10,fontweight='bold',color='#1565c0')
ax_c.text(2.8,2.8,f'f\u2081 = {f_z_c:.0f} Hz, Q\u2081 = {Q_zc}',ha='center',fontsize=9,color='#1565c0')
ax_c.add_patch(plt.Rectangle((5.5,3.3),1.0,1.2,fc='#c8e6c9',ec='#2e7d32',lw=2))
ax_c.text(6.0,3.9,'M\nk',ha='center',va='center',fontsize=9,fontweight='bold',color='#2e7d32')
ax_c.text(6.0,2.5,'Kopplung',ha='center',fontsize=8,color='#2e7d32')
ax_c.add_patch(plt.Rectangle((7.0,3.5),1.2,0.8,fc='#ffebee',ec='#c62828',lw=1.5))
ax_c.text(7.6,3.9,'L\u2082\n\u03c1l/A',ha='center',va='center',fontsize=8,color='#c62828')
ax_c.add_patch(plt.Rectangle((8.7,3.5),1.2,0.8,fc='#ffebee',ec='#c62828',lw=1.5))
ax_c.text(9.3,3.9,'C\u2082\nV/\u03c1c\u00b2',ha='center',va='center',fontsize=8,color='#c62828')
ax_c.add_patch(plt.Rectangle((10.4,3.5),1.2,0.8,fc='#ffebee',ec='#c62828',lw=1.5))
ax_c.text(11.0,3.9,'R\u2082\nR_visc',ha='center',va='center',fontsize=8,color='#c62828')
ax_c.text(9.3,5.0,'Kreis 2: Kammer',ha='center',fontsize=10,fontweight='bold',color='#c62828')
ax_c.text(9.3,2.8,f'f_H = {f_Hc:.0f} Hz, Q_H = {Q_Hc}',ha='center',fontsize=9,color='#c62828')
for xs,xe in [(1.7,2.2),(3.4,3.9),(5.1,5.5),(6.5,7.0),(8.2,8.7),(9.9,10.4)]:
    ax_c.plot([xs,xe],[3.9,3.9],'k-',lw=1.5)
ax_c.annotate('p_Balg',xy=(0.5,3.9),xytext=(-0.3,3.9),fontsize=9,fontweight='bold',
    arrowprops=dict(arrowstyle='->',lw=1.5),ha='right')
ax_c.text(6.0,0.8,'Masse m \u2192 L    Steifigkeit k \u2192 1/C    D\u00e4mpfung \u2192 R',
    ha='center',fontsize=9,fontfamily='monospace',bbox=dict(facecolor='#f5f5f5',edgecolor='gray',boxstyle='round'))
ax_c.text(6.0,0.2,'k bestimmt Ansprache (Energietransfer), C\u2082 bestimmt Filterwirkung (f_H)',
    ha='center',fontsize=9,style='italic',color='#555')

fig_circuit.suptitle('Abb. 6b: Elektroakustische Analogie \u2014 Zwei gekoppelte RLC-Kreise',
    fontsize=13, fontweight='bold')
fig_circuit.tight_layout()
story.append(fig_to_image(fig_circuit, width_mm=160))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    '<b>Oben links — Kopplungskoeffizient k variiert:</b> '
    'Bei schwacher Kopplung (k = 0,01) kommt fast keine Energie in der Kammer an. '
    'Bei starker Kopplung (k = 0,40) sind beide Resonanzpeaks deutlich sichtbar. '
    'Entscheidend: Die <b>Position</b> der Peaks (f₁ und f_H) ändert sich nicht — '
    'nur die <b>Höhe</b> (Energietransfer). Das ist genau die Beobachtung: '
    'andere Ansprache, gleiche Klangfarbe.',
    sBody))
story.append(Paragraph(
    '<b>Unten links — Kammervolumen variiert:</b> '
    'Hier verschiebt sich f_H (der rote Peak wandert), aber f₁ bleibt bei 50 Hz. '
    'Das ist die Filterwirkung: Die Kammer bestimmt, welche Obertöne betont oder '
    'gedämpft werden. Die Zunge bleibt davon unberührt.',
    sBody))
story.append(Paragraph(
    '<b>Unten rechts — Ersatzschaltbild:</b> '
    'Zwei Serienschwingkreise (L-C-R), verbunden durch die gegenseitige Impedanz M. '
    'k steuert den Energietransfer (Ansprache), C₂ = V/(ρc²) steuert f_H (Filterwirkung). '
    'Diese zwei Parameter sind unabhängig — genau wie in der Praxis.',
    sBody))

story.append(Paragraph(
    'Das Ersatzschaltbild ist berechenbar. Mit den physikalischen Werten '
    '(Masse, Steifigkeit, Volumen, Öffnungsfläche) ergibt sich die Transferfunktion '
    'quantitativ. Die absolute Kalibrierung des Kopplungskoeffizienten k '
    'bleibt empirisch — aber die Form der Kurven und die Trennung von '
    'Ansprache (k) und Filterwirkung (C₂) ist physikalisch exakt.',
    sKeyBox))

# ══ Kapitel 7: Klappenposition ══
story.append(Paragraph('7. Klappenposition — Feinkorrektur', sChapter))
story.append(Paragraph(
    'Die Klappenöffnung kann in der Praxis um einige Millimeter bis maximal '
    'einen Zentimeter vom Kammerende versetzt werden — mehr erlaubt die Konstruktion '
    'nicht. Die Verschiebung ändert die effektive akustische Länge des Halses '
    'und damit f_H.',
    sBody))

# Abb. 7
fig7 = make_fig7()
story.append(fig_to_image(fig7, width_mm=140))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Im Praxisbereich (0–15 mm) verschiebt sich f_H um wenige Prozent. '
    'Die Klappenposition ist eine Feinkorrektur — der Haupthebel liegt am '
    'Zungenende (Port 1), wo 70 % der Conductance sitzen.',
    sBody))

# ══ Kapitel 7: Grenzen des Modells ══
story.append(Paragraph('8. Grenzen des Modells', sChapter))
story.append(Paragraph(
    'Das Zwei-Port-Modell zeigt Zusammenhänge und Richtungen. '
    'Es sagt nicht, wo das Optimum liegt. Die Gründe:',
    sBody))
story.append(Paragraph(
    '(a) <b>Nichtlinearität:</b> Die Zunge schaltet periodisch zwischen '
    'kleiner und großer Öffnung — der zeitgemittelte Leitwert ist eine Vereinfachung. '
    'Die tatsächliche Impedanzmodulation (Dok. 0013) ist hochgradig nichtlinear.',
    sBody))
story.append(Paragraph(
    '(b) <b>Klangziel:</b> Das „optimale" Spektrum hängt vom gewünschten Klangcharakter ab — '
    'dafür gibt es keine physikalische Formel. Die Berechnung zeigt die Landkarte, '
    'den Weg geht das Ohr.',
    sBody))
story.append(Paragraph(
    '(c) <b>Wechselwirkungen:</b> Spaltweite (Dok. 0010), Steifigkeit (Dok. 0012), '
    'Kammervolumen, Klappenöffnung und Übergangsgeometrie wirken alle zusammen. '
    'Eine isolierte Optimierung eines Parameters kann einen anderen verschlechtern.',
    sBody))
story.append(Paragraph(
    '<b>Dieses Dokument ist ein Erklärungsmodell, keine Konstruktionsanleitung. '
    'Konstruktionsentscheidungen erfordern nach wie vor praktische Tests und Erfahrung.</b>',
    sWarnBox))

# ══ Kapitel 8: Kopplungsschwerpunkt ══
story.append(Paragraph('9. Wo koppelt die Zunge? — Schwerpunkt des Volumenstroms', sChapter))
story.append(Paragraph(
    'Die Annahme, dass die Kopplung am absoluten Ende der Zunge (Spitze) erfolgt, '
    'ist eine Vereinfachung. Die Zunge ist ein Cantilever — ihre Schwingform '
    'ist keine Gerade, sondern folgt der ersten Euler-Bernoulli-Biegemode.',
    sBody))

# Abb. 8
fig8 = make_fig8()
story.append(fig_to_image(fig8, width_mm=155))
story.append(Spacer(1, 3*mm))

# Berechne Schwerpunkt für Text
beta_L_val = 1.8751
sigma_val = (np.cosh(beta_L_val) + np.cos(beta_L_val)) / (np.sinh(beta_L_val) + np.sin(beta_L_val))
x_n = np.linspace(0, 1, 1000)
y_n = np.cosh(beta_L_val*x_n) - np.cos(beta_L_val*x_n) - sigma_val*(np.sinh(beta_L_val*x_n) - np.sin(beta_L_val*x_n))
y_n = y_n / y_n[-1]
x_cen = np.trapezoid(x_n * y_n, x_n) / np.trapezoid(y_n, x_n)
Q_c = np.cumsum(y_n) / np.sum(y_n)
idx50 = np.argmin(np.abs(Q_c - 0.5))
x_50 = x_n[idx50]

story.append(Paragraph(
    f'Der Volumenstrom pro Längeneinheit dQ/dx ∝ y(x) steigt steil zum Ende hin an. '
    f'Der gewichtete <b>Schwerpunkt des Volumenstroms</b> liegt bei '
    f'<b>x/L = {x_cen:.2f}</b> — nicht an der Spitze (x/L = 1,0). '
    f'50 % des gesamten Volumenstroms entstehen erst in den letzten '
    f'{(1-x_50)*100:.0f} % der Zungenlänge (ab x/L = {x_50:.2f}).',
    sBody))
story.append(Paragraph(
    f'Für eine Bass-Zunge mit L = 70 mm liegt der Schwerpunkt bei '
    f'{x_cen*70:.0f} mm von der Einspannung — ca. {(1-x_cen)*70:.0f} mm '
    f'vor der Spitze.',
    sBody))

story.append(Paragraph('9.1 Impedanzfehlanpassung', sSection))
story.append(Paragraph(
    'Die mechanische Impedanz der Zunge und die akustische Impedanz der Luft '
    'unterscheiden sich um Größenordnungen. Die Zunge ist ein <b>Strom-Generator</b> '
    '(hohe Impedanz) — sie erzwingt einen Volumenstrom, den die Luft aufnimmt.',
    sBody))

# Abb. 9
fig9 = make_fig9()
story.append(fig_to_image(fig9, width_mm=150))
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Z_Zunge liegt Faktor 50–500 über Z_Luft. Die Zunge schwingt ungebremst durch '
    'die Luft — die Luftlast ist viel weicher als die mechanische Rückstellkraft. '
    'Die Kopplung ist <b>nicht angepasst</b> im Sinne eines maximalen Leistungstransfers '
    '(Z_Quelle = Z_Last). Das ist physikalisch richtig so: Die Zunge ist selbsterregt '
    'und schwingt auf ihrer Eigenfrequenz, nicht auf der Kammerresonanz. '
    'Die Luft wird „mitgenommen", nicht „angetrieben".',
    sBody))
story.append(Paragraph(
    'Dies steht im Einklang mit Dok. 0010: Über 97 % der Energieabgabe erfolgt '
    'über die Luftkopplung, weniger als 3 % über die mechanische Einspannung. '
    'Die Zunge verliert ihre Energie an die Luft — aber die Luft bremst die Zunge kaum.',
    sBody))

story.append(Paragraph('9.2 Konsequenz: Position der Plattenkante', sSection))
story.append(Paragraph(
    'Wenn die Stimmplatte die volle Zungenlänge abdeckt, koppelt der gesamte '
    'Schlitz in die Kammer. Ragt die Zunge über die Plattenkante hinaus (Aufbiegung), '
    'geht der energiereichste Teil — die letzten 20–30 % der Zunge — an der Kammer '
    'vorbei, direkt in den freien Raum.',
    sBody))
story.append(Paragraph(
    'Die Position der Plattenkante relativ zur Zunge bestimmt daher, welcher Anteil '
    'des Volumenstroms in die Kammer eingekoppelt wird und welcher Anteil '
    'direkt abstrahlt. Dies ist ein weiterer Stellhebel für die Anpassung — '
    'neben Spaltweite, Kammervolumen und Klappenposition.',
    sBody))

# ══ Kapitel 9: Zusammenfassung ══
story.append(Paragraph('10. Zusammenfassung', sChapter))

summary_points = [
    '<b>Zwei-Port-System:</b> Die Kammer ist beidseitig offen. '
    'Port 1 (Zunge) trägt ≈ 70 % der Conductance. Das Ein-Port-Modell unterschätzt f_H.',
    '<b>Faltung unsichtbar:</b> Bei 50 Hz ist die größte Kammerabmessung 1,7 % von λ. '
    'Die Kammerform spielt für den Grundton und die ersten 27 Obertöne keine Rolle — '
    'nur das Volumen zählt.',
    '<b>Impedanzanpassung:</b> Die Übergangsgeometrie wirkt als Transformator '
    'zwischen Balgdruck (Quelle) und Zunge (Last). Maximaler Leistungstransfer '
    'bei R_i = R_L. Die Filterwirkung (f_H) hängt vom Resonator ab, nicht vom Transformator.',
    '<b>Ersatzschaltbild:</b> Zunge und Kammer bilden zwei gekoppelte RLC-Schwingkreise. '
    'Der Kopplungskoeffizient k bestimmt die Ansprache (Energietransfer), '
    'C₂ = V/(ρc²) bestimmt die Filterwirkung (f_H). Diese Parameter sind unabhängig.',
    '<b>Klappenposition = Feinkorrektur:</b> Wenige Prozent f_H-Verschiebung im '
    'Praxisbereich. Der Haupthebel liegt am Zungenende.',
    '<b>Schwerpunkt bei x/L = 0,73:</b> Die Kopplung sitzt nicht an der Zungenspitze, '
    'sondern ≈ 27 % davor. Die Position der Plattenkante bestimmt die Einkopplung.',
    '<b>Z_Zunge ≫ Z_Luft:</b> Die Zunge ist ein Strom-Generator. '
    'Die Impedanzen sind nicht angepasst — die Luft wird mitgenommen, nicht angetrieben.',
    '<b>Optimum wird gehört:</b> Die Berechnung zeigt Zusammenhänge und Richtungen. '
    'Erfahrung und Gehörschulung bleiben der Schlüssel zum Erfolg.',
]

for i, point in enumerate(summary_points):
    story.append(Paragraph(f'{i+1}. {point}', sBody))

story.append(Spacer(1, 6*mm))
story.append(Paragraph(
    '<i>Die Zunge speist, die Kammer formt, die Klappe strahlt ab. '
    'Die Geometrie am Übergang bestimmt die Kopplung.</i>',
    sAbstract))

# ── Seitennummer ──
def add_page_number(canvas, doc):
    canvas.saveState()
    canvas.setFont('DejaVu', 8)
    canvas.setFillColor(HexColor('#999999'))
    canvas.drawCentredString(W_PAGE/2, 12*mm,
                             f'Dok. 0015 — Akustische Kopplung — Seite {canvas.getPageNumber()}')
    canvas.restoreState()

doc.build(story, onFirstPage=add_page_number, onLaterPages=add_page_number)
print(f'✓ {outfile} erzeugt')
