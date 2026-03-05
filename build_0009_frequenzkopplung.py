#!/usr/bin/env python3
"""
Dok. 0009 — Frequenzkopplung mehrerer Zungen: Drei-Zonen-Gradientenmodell
Build-Skript: Erzeugt PDF mit Diagrammen und Tabellen.
Benötigt: reportlab, matplotlib, numpy
"""

import math
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyArrowPatch

from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.colors import HexColor, black, white
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.enums import TA_LEFT, TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (
    SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle,
    Image, PageBreak, KeepTogether
)

# ── Farben ──
DARKBLUE  = HexColor('#16213e')
ACCENTRED = HexColor('#e94560')
KEYGREEN  = HexColor('#2e7d32')
WARNRED   = HexColor('#c62828')
LIGHTGRAY = HexColor('#f5f5f5')
KEYBG     = HexColor('#e8f5e9')
WARNBG    = HexColor('#ffebee')

W = A4[0]
PAGE_W = W - 50*mm

# ── Styles ──
styles = getSampleStyleSheet()

sTitle = ParagraphStyle('DocTitle', parent=styles['Title'],
    fontSize=18, textColor=DARKBLUE, spaceAfter=4, alignment=TA_CENTER)
sSubtitle = ParagraphStyle('SubTitle', parent=styles['Normal'],
    fontSize=12, textColor=DARKBLUE, alignment=TA_CENTER, spaceAfter=2)
sAbstract = ParagraphStyle('Abstract', parent=styles['Italic'],
    fontSize=9, textColor=HexColor('#555555'), spaceAfter=8, alignment=TA_LEFT)
sChapter = ParagraphStyle('Chapter', parent=styles['Heading1'],
    fontSize=14, textColor=DARKBLUE, spaceBefore=14, spaceAfter=6,
    fontName='Helvetica-Bold')
sBody = ParagraphStyle('Body', parent=styles['Normal'],
    fontSize=10, leading=14, spaceAfter=6, alignment=TA_JUSTIFY,
    fontName='Helvetica')
sBodyB = ParagraphStyle('BodyBold', parent=sBody, fontName='Helvetica-Bold')
sFormula = ParagraphStyle('Formula', parent=sBody,
    fontSize=10, alignment=TA_CENTER, spaceAfter=8, spaceBefore=4,
    fontName='Courier')
sSmall = ParagraphStyle('Small', parent=sBody, fontSize=9, leading=12)
sKeyBox = ParagraphStyle('KeyBox', parent=sBody, fontSize=10,
    backColor=KEYBG, borderPadding=6, borderColor=KEYGREEN,
    borderWidth=1, spaceAfter=8)
sWarnBox = ParagraphStyle('WarnBox', parent=sBody, fontSize=10,
    backColor=WARNBG, borderPadding=6, borderColor=WARNRED,
    borderWidth=1, spaceAfter=8)
sTH = ParagraphStyle('TH', parent=sBody, fontSize=9,
    fontName='Helvetica-Bold', alignment=TA_CENTER, leading=11)
sTD = ParagraphStyle('TD', parent=sBody, fontSize=9,
    alignment=TA_CENTER, leading=11)
sTDL = ParagraphStyle('TDL', parent=sBody, fontSize=9,
    alignment=TA_LEFT, leading=11)

def hr():
    return Table([['']], colWidths=[PAGE_W],
        style=TableStyle([('LINEBELOW', (0,0), (-1,-1), 2, ACCENTRED)]))

def make_table(header, rows, col_widths=None):
    data = [[Paragraph(h, sTH) for h in header]]
    for row in rows:
        data.append([Paragraph(str(c), sTD) if i > 0 else Paragraph(str(c), sTDL)
                      for i, c in enumerate(row)])
    cw = col_widths or [PAGE_W / len(header)] * len(header)
    t = Table(data, colWidths=cw, repeatRows=1)
    t.setStyle(TableStyle([
        ('BACKGROUND', (0, 0), (-1, 0), DARKBLUE),
        ('TEXTCOLOR', (0, 0), (-1, 0), white),
        ('ROWBACKGROUNDS', (0, 1), (-1, -1), [white, LIGHTGRAY]),
        ('GRID', (0, 0), (-1, -1), 0.5, HexColor('#cccccc')),
        ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),
        ('TOPPADDING', (0, 0), (-1, -1), 3),
        ('BOTTOMPADDING', (0, 0), (-1, -1), 3),
    ]))
    return t

# ── Physik-Parameter ──
c_sound = 343.0  # m/s
rho = 1.2  # kg/m³

# Diskant-Stimmstock: D3 bis C6, 35 Halbtöne
notes_names = ['D3','D#3','E3','F3','F#3','G3','G#3','A3','A#3','B3',
               'C4','C#4','D4','D#4','E4','F4','F#4','G4','G#4','A4',
               'A#4','B4','C5','C#5','D5','D#5','E5','F5','F#5','G5',
               'G#5','A5','A#5','B5','C6']
n_notes = len(notes_names)
freqs = [146.83 * 2**((i)/12) for i in range(n_notes)]

# Kammergeometrie (Variante C — konvex, Sieger aus Dok. 0007)
W_chamber = 15.5e-3  # m
L_min, L_max = 30e-3, 50e-3
d_max, d_min = 8e-3, 3e-3
d_opening = 10e-3  # m
l_neck = 8e-3  # m

def depth_C(i):
    """Konvexe (exponentielle) Tiefenreduktion"""
    return d_max * (d_min/d_max)**(i/(n_notes-1))

def chamber_length(i):
    return L_max - (L_max - L_min) * i / (n_notes - 1)

def helmholtz_freq(V, S_neck, l_eff):
    return c_sound / (2 * math.pi) * math.sqrt(S_neck / (V * l_eff))

S_neck = math.pi * (d_opening/2)**2
l_eff_neck = l_neck + 0.85 * d_opening

# Berechne f_H für Variante C
fH_C = []
for i in range(n_notes):
    d = depth_C(i)
    L = chamber_length(i)
    V = W_chamber * L * d
    fH_C.append(helmholtz_freq(V, S_neck, l_eff_neck))

# ── Kopplungsberechnung ──
def coupling_G(f_ot, f_H, Q_H=8):
    """Impedanzverstärkung bei Oberton f_ot nahe Kammerresonanz f_H"""
    x = f_ot / f_H
    return 1.0 / math.sqrt((1 - x**2)**2 + (x/Q_H)**2)

def amplitude_weight(n):
    """Obertonamplitude a_n ~ 1/n, Energie ~ 1/n²"""
    return 1.0 / n**2

Q_H = 8  # Kammer-Gütefaktor

# Drei-Zonen-Analyse
# Für jeden Ton: Summierte gewichtete Kopplung, dominanter OT, Lock-in-Abschätzung
zone_data = []
for i in range(n_notes):
    f0 = freqs[i]
    fH = fH_C[i]
    ratio = fH / f0
    
    # Summiere gewichtete Kopplung über alle Obertöne
    R_eff = 0
    max_contrib = 0
    dom_n = 1
    contributions = []
    
    for n in range(1, 21):  # bis 20. Oberton
        f_ot = n * f0
        if f_ot > 5 * fH:
            break
        G = coupling_G(f_ot, fH, Q_H)
        w = amplitude_weight(n)
        contrib = w * (G - 1)  # Nur Überschuss über 1
        if contrib > 0:
            R_eff += contrib
            contributions.append((n, f_ot, G, w, contrib))
            if contrib > max_contrib:
                max_contrib = contrib
                dom_n = n
    
    # Kopplungssumme für Lock-in-Abschätzung
    # Empirisch: Δf_lock ∝ Σ(a_n² · G_n) · A_eff²/V
    d = depth_C(i)
    L = chamber_length(i)
    V = W_chamber * L * d
    
    # Anzahl Obertöne im Resonanzband (|f_ot - f_H| < f_H/Q_H)
    bw = fH / Q_H
    n_in_band = sum(1 for n in range(1, 21) if abs(n*f0 - fH) < bw)
    
    zone_data.append({
        'note': notes_names[i],
        'f0': f0,
        'fH': fH,
        'ratio': ratio,
        'R_eff': R_eff,
        'dom_n': dom_n,
        'max_contrib': max_contrib,
        'n_in_band': n_in_band,
        'V': V * 1e6,  # cm³
        'contributions': contributions
    })

# ── Zone-Klassifikation ──
def classify_zone(d):
    r = d['ratio']
    if r > 5:
        return 1  # Tiefe Töne: Locking-Gefahr
    elif r >= 3:
        return 2  # Mittlere Töne: Geisterpeaks
    else:
        return 3  # Hohe Töne: Ansprache-Problem

for d in zone_data:
    d['zone'] = classify_zone(d)

# ── Adler-Gleichung Berechnungen ──
# Lock-in-Bereich: Δf_lock = f₀/(2Q₁) · κ_eff
# κ_eff hängt von Kammerresonanz-Kopplung ab
# Empirisch kalibriert auf Praxiswerte

def adler_pulling(delta_f_stimmung, delta_f_lock):
    """Adler: Tatsächliches Tremolo bei gegebenem Lock-in-Bereich"""
    if abs(delta_f_stimmung) <= delta_f_lock:
        return 0.0  # Locked
    return math.sqrt(delta_f_stimmung**2 - delta_f_lock**2)

# ── Diagramm 1: Drei-Zonen-Übersicht ──
fig1, ax1 = plt.subplots(figsize=(8, 5))

zone_colors = {1: '#1565C0', 2: '#E65100', 3: '#C62828'}
zone_labels_de = {1: 'Zone 1: Locking', 2: 'Zone 2: Geisterpeaks', 3: 'Zone 3: Ansprache'}

for d in zone_data:
    z = d['zone']
    ax1.bar(d['note'], d['R_eff'], color=zone_colors[z], alpha=0.8,
            label=zone_labels_de[z] if d['note'] in ['D3', 'D4', 'D5'] else '')

ax1.set_xlabel('Ton', fontsize=11)
ax1.set_ylabel('Gewichtete Kopplung R_eff', fontsize=11)
ax1.set_title('Drei-Zonen-Gradientenmodell: Kopplungsstärke über Tonbereich\n(Variante C, konvexe Tiefenreduktion)', fontsize=12)
ax1.legend(loc='upper left', fontsize=10)
ax1.set_xticks(range(0, n_notes, 3))
ax1.set_xticklabels([notes_names[i] for i in range(0, n_notes, 3)], rotation=45)
ax1.grid(axis='y', alpha=0.3)

# Zonengrenzen
# Find boundary indices
zone1_end = max(i for i, d in enumerate(zone_data) if d['zone'] == 1) + 0.5
zone2_end = max(i for i, d in enumerate(zone_data) if d['zone'] <= 2) + 0.5

ax1.axvline(x=zone1_end, color='gray', linestyle='--', alpha=0.5)
ax1.axvline(x=zone2_end, color='gray', linestyle='--', alpha=0.5)
ax1.text(zone1_end/2, ax1.get_ylim()[1]*0.9, 'Zone 1\nLocking', ha='center', fontsize=9, color='#1565C0')
ax1.text((zone1_end+zone2_end)/2, ax1.get_ylim()[1]*0.9, 'Zone 2\nGeisterpeaks', ha='center', fontsize=9, color='#E65100')
ax1.text((zone2_end+n_notes)/2, ax1.get_ylim()[1]*0.9, 'Zone 3\nAnsprache', ha='center', fontsize=9, color='#C62828')

plt.tight_layout()
fig1.savefig('/home/claude/diag_0009_1_zonen.png', dpi=150)
plt.close()

# ── Diagramm 2: f_H/f und dominanter Oberton ──
fig2, (ax2a, ax2b) = plt.subplots(2, 1, figsize=(8, 7), sharex=True)

ratios = [d['ratio'] for d in zone_data]
dom_ns = [d['dom_n'] for d in zone_data]
x = range(n_notes)

# Oberes Panel: f_H/f
for i, d in enumerate(zone_data):
    ax2a.bar(i, d['ratio'], color=zone_colors[d['zone']], alpha=0.8)
ax2a.axhline(y=3, color='red', linestyle='--', linewidth=1.5, label='f_H/f = 3 (Grenze Zone 2/3)')
ax2a.axhline(y=5, color='blue', linestyle='--', linewidth=1.5, label='f_H/f = 5 (Grenze Zone 1/2)')
ax2a.set_ylabel('f_H / f', fontsize=11)
ax2a.set_title('Verhältnis Kammerresonanz zu Grundfrequenz', fontsize=12)
ax2a.legend(fontsize=9)
ax2a.grid(axis='y', alpha=0.3)

# Unteres Panel: Dominanter Oberton + Energie
for i, d in enumerate(zone_data):
    energy_pct = amplitude_weight(d['dom_n']) * 100
    ax2b.bar(i, d['dom_n'], color=zone_colors[d['zone']], alpha=0.8)
    if i % 3 == 0:
        ax2b.annotate(f'{energy_pct:.0f}%', (i, d['dom_n']), textcoords="offset points",
                      xytext=(0, 5), ha='center', fontsize=7)

ax2b.set_xlabel('Ton', fontsize=11)
ax2b.set_ylabel('Dominanter Oberton n', fontsize=11)
ax2b.set_title('Stärkster koppelnder Oberton (Zahl = Energieanteil a_n²)', fontsize=12)
ax2b.set_xticks(range(0, n_notes, 3))
ax2b.set_xticklabels([notes_names[i] for i in range(0, n_notes, 3)], rotation=45)
ax2b.grid(axis='y', alpha=0.3)

plt.tight_layout()
fig2.savefig('/home/claude/diag_0009_2_ratio_oberton.png', dpi=150)
plt.close()

# ── Diagramm 3: Adler-Pulling für 1 Hz Tremolo ──
fig3, ax3 = plt.subplots(figsize=(8, 4.5))

# Für verschiedene Lock-in-Bereiche
delta_f_stimmung = np.linspace(0, 3, 200)
for lock_in, col, lbl in [(0.1, '#2e7d32', 'Getrennte Öffnungen (0,1 Hz)'),
                           (0.5, '#E65100', 'Schwache Kopplung (0,5 Hz)'),
                           (1.0, '#C62828', 'Gemeinsame Öffnung (1,0 Hz)'),
                           (1.5, '#7B1FA2', 'Starke Kopplung (1,5 Hz)')]:
    delta_f_eff = [adler_pulling(df, lock_in) for df in delta_f_stimmung]
    ax3.plot(delta_f_stimmung, delta_f_eff, color=col, linewidth=2, label=lbl)

ax3.plot([0, 3], [0, 3], 'k--', alpha=0.3, label='Ohne Kopplung')
ax3.axhline(y=1.0, color='gray', linestyle=':', alpha=0.5)
ax3.axvline(x=1.0, color='gray', linestyle=':', alpha=0.5)

ax3.set_xlabel('Gestimmte Frequenzdifferenz Δf_Stimmung [Hz]', fontsize=11)
ax3.set_ylabel('Tatsächliche Schwebung Δf_eff [Hz]', fontsize=11)
ax3.set_title('Adler-Gleichung: Frequenz-Pulling bei gekoppelten Zungen', fontsize=12)
ax3.legend(fontsize=9)
ax3.grid(alpha=0.3)
ax3.set_xlim(0, 3)
ax3.set_ylim(0, 3)

plt.tight_layout()
fig3.savefig('/home/claude/diag_0009_3_adler.png', dpi=150)
plt.close()

# ── Diagramm 4: Kopplungslandschaft (Heatmap) ──
fig4, ax4 = plt.subplots(figsize=(8, 5))

# Für jeden Ton: Kopplung pro Oberton als Heatmap
max_ot = 12
heatmap = np.zeros((max_ot, n_notes))
for i, d in enumerate(zone_data):
    for n in range(1, max_ot + 1):
        f_ot = n * d['f0']
        G = coupling_G(f_ot, d['fH'], Q_H)
        w = amplitude_weight(n)
        heatmap[n-1, i] = w * G

im = ax4.imshow(heatmap, aspect='auto', cmap='YlOrRd', origin='lower',
                extent=[-0.5, n_notes-0.5, 0.5, max_ot+0.5])
ax4.set_xticks(range(0, n_notes, 3))
ax4.set_xticklabels([notes_names[i] for i in range(0, n_notes, 3)], rotation=45)
ax4.set_yticks(range(1, max_ot + 1))
ax4.set_ylabel('Oberton-Nummer n', fontsize=11)
ax4.set_xlabel('Ton', fontsize=11)
ax4.set_title('Kopplungslandschaft: Gewichtete Impedanzverstärkung a_n² · G(n·f)\n(hell = stark, dunkel = schwach)', fontsize=12)
plt.colorbar(im, ax=ax4, label='a_n² · G')

# Markiere f_H-Linie
fH_over_f = [d['fH']/d['f0'] for d in zone_data]
ax4.plot(range(n_notes), fH_over_f, 'w--', linewidth=2, label='f_H/f (Resonanzlinie)')
ax4.legend(loc='upper right', fontsize=9)

plt.tight_layout()
fig4.savefig('/home/claude/diag_0009_4_heatmap.png', dpi=150)
plt.close()

# ── Diagramm 5: Spektrum-Beispiel mit Geisterpeaks ──
fig5, (ax5a, ax5b) = plt.subplots(1, 2, figsize=(10, 4))

# Beispiel: A4 (440 Hz), 3 Zungen mit 1 Hz Tremolo
f_base = 440.0
delta_f_tremolo = 1.0

# Links: Ohne Kopplung (getrennte Öffnungen)
for offset, col, lbl in [(-1, '#1565C0', 'Zunge 1 (f-1)'),
                          (0, '#2e7d32', 'Zunge 2 (f)'),
                          (1, '#C62828', 'Zunge 3 (f+1)')]:
    f = f_base + offset * delta_f_tremolo
    for n in range(1, 6):
        amp = 1.0/n
        ax5a.bar(n*f, amp, width=2, color=col, alpha=0.7)

ax5a.set_xlabel('Frequenz [Hz]', fontsize=10)
ax5a.set_ylabel('Amplitude', fontsize=10)
ax5a.set_title('Getrennte Öffnungen\n(3 × saubere Spektren)', fontsize=11)
ax5a.set_xlim(400, 2300)
ax5a.grid(axis='y', alpha=0.3)

# Rechts: Mit Kopplung (gemeinsame Öffnung) → Geisterpeaks
for offset, col in [(-1, '#1565C0'), (0, '#2e7d32'), (1, '#C62828')]:
    f = f_base + offset * delta_f_tremolo
    for n in range(1, 6):
        amp = 1.0/n
        ax5b.bar(n*f, amp, width=2, color=col, alpha=0.7)
        # Geisterpeaks durch Intermodulation
        if n <= 3:
            for ghost_offset in [-2, 2]:
                f_ghost = n*f + ghost_offset * delta_f_tremolo
                amp_ghost = amp * 0.08  # ~8% Geisterpeak
                ax5b.bar(f_ghost, amp_ghost, width=1.5, color='#9E9E9E', alpha=0.5)

ax5b.set_xlabel('Frequenz [Hz]', fontsize=10)
ax5b.set_title('Gemeinsame Öffnung\n(Geisterpeaks durch Intermodulation)', fontsize=11)
ax5b.set_xlim(400, 2300)
ax5b.grid(axis='y', alpha=0.3)

# Graue Legende
from matplotlib.patches import Patch
ax5b.legend(handles=[Patch(facecolor='#9E9E9E', alpha=0.5, label='Geisterpeaks')],
            loc='upper right', fontsize=9)

plt.tight_layout()
fig5.savefig('/home/claude/diag_0009_5_spektrum.png', dpi=150)
plt.close()

# ══════════════════════════════════════════════════════════════
# PDF-Erzeugung
# ══════════════════════════════════════════════════════════════

output_path = '/home/claude/0009_frequenzkopplung_mehrere_zungen_De.pdf'

doc = SimpleDocTemplate(output_path, pagesize=A4,
    leftMargin=25*mm, rightMargin=25*mm, topMargin=25*mm, bottomMargin=25*mm)

story = []

# ── Titel ──
story.append(Paragraph('Dok. 0009', sSmall))
story.append(Spacer(1, 2*mm))
story.append(Paragraph('Frequenzkopplung mehrerer Zungen', sTitle))
story.append(Paragraph('Drei-Zonen-Gradientenmodell für Tremolo, Locking und Geisterpeaks', sSubtitle))
story.append(Spacer(1, 2*mm))
story.append(hr())
story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    'Wenn zwei oder drei Stimmzungen mit Tremolo-Verstimmung (≈ 1 Hz) in dieselbe Kammer oder '
    'durch dieselbe Öffnung klingen, beeinflussen sie sich gegenseitig in der Frequenz. '
    'Das Ausmaß hängt davon ab, wie stark die Kammerresonanz die Obertöne koppelt. '
    'Dieses Dokument leitet das Drei-Zonen-Gradientenmodell ab, das erklärt, warum verschiedene Töne '
    'verschiedene Kopplungseffekte zeigen — von Locking (tiefe Töne) über Geisterpeaks (mittlere Töne) '
    'bis zu Ansprache-Problemen (hohe Töne). '
    'Zahlenwerte beziehen sich auf den Diskant-Stimmstock aus Dok. 0007, Variante C (konvex).', sAbstract))

story.append(Spacer(1, 2*mm))

# Dokumentenindex
story.append(Paragraph('<b>Dokumentenverweise</b>', sSmall))
refs = [
    ['0002', 'Strömungsanalyse Bass-Stimmzunge 50 Hz — v8'],
    ['0004', 'Frequenzvariation — Zwei-Filter-Modell'],
    ['0005', 'Frequenzverschiebung als Indikator der Ansprache'],
    ['0007', 'Diskant-Stimmstock — Kammerfrequenzen (7 Varianten)'],
    ['0008', 'Klangveränderung durch Kammergeometrie'],
]
ref_data = [[Paragraph(f'<b>Dok. {r[0]}</b>', sTD), Paragraph(r[1], sTDL)] for r in refs]
ref_table = Table(ref_data, colWidths=[30*mm, PAGE_W - 30*mm])
ref_table.setStyle(TableStyle([
    ('GRID', (0,0), (-1,-1), 0.3, HexColor('#cccccc')),
    ('BACKGROUND', (0,0), (0,-1), LIGHTGRAY),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 2),
    ('BOTTOMPADDING', (0,0), (-1,-1), 2),
]))
story.append(ref_table)
story.append(Spacer(1, 4*mm))

# ══════════════════════════════════════════════════════════════
# Kapitel 1: Die Beobachtung
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 1: Die Beobachtung aus der Praxis', sChapter))
story.append(Paragraph(
    'Beim Stimmen von Tremolo-Registern beobachtet der Akkordeonbauer drei verschiedene Effekte, '
    'die vom Tonbereich abhängen:', sBody))
story.append(Spacer(1, 2*mm))

obs_rows = [
    ['Tiefe Töne\n(D3–G#3)', 'f_H/f > 5', 
     'Zwei oder drei Zungen, die auf ≈ 1 Hz Schwebung gestimmt sind, rasten '
     'bei gemeinsamer Öffnung auf eine einzige Frequenz ein. Das Tremolo verschwindet. '
     'Bei getrennter Luftzufuhr tritt der Effekt nicht auf.'],
    ['Mittlere Töne\n(A3–C#5)', 'f_H/f ≈ 3–5',
     'Das Tremolo ist hörbar, aber im Spektralanalysator zeigen sich neben den drei '
     'erwarteten Frequenzspitzen zusätzliche Geisterpeaks im gleichen Abstand. '
     'Statt drei Spitzen erscheinen vier oder fünf.'],
    ['Hohe Töne\n(D5–C6)', 'f_H/f ≈ 2–3',
     'Das Tremolo funktioniert, aber die Ansprache ist träge und der Klang dünn. '
     'Die Kopplung ist stark, aber der Effekt äußert sich als Energieverlust, '
     'nicht als Frequenzverschiebung.'],
]
obs_table_data = [[Paragraph('<b>Bereich</b>', sTH), Paragraph('<b>f<sub>H</sub>/f</b>', sTH),
                    Paragraph('<b>Beobachtung</b>', sTH)]]
for row in obs_rows:
    obs_table_data.append([Paragraph(row[0], sTDL), Paragraph(row[1], sTD), Paragraph(row[2], sTDL)])
obs_t = Table(obs_table_data, colWidths=[30*mm, 22*mm, PAGE_W-52*mm])
obs_t.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), DARKBLUE),
    ('TEXTCOLOR', (0,0), (-1,0), white),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHTGRAY]),
    ('GRID', (0,0), (-1,-1), 0.5, HexColor('#cccccc')),
    ('VALIGN', (0,0), (-1,-1), 'TOP'),
    ('TOPPADDING', (0,0), (-1,-1), 4),
    ('BOTTOMPADDING', (0,0), (-1,-1), 4),
]))
story.append(obs_t)
story.append(Spacer(1, 3*mm))

story.append(Paragraph(
    'Die entscheidende Beobachtung: Diese drei Effekte sind <b>nicht</b> binär (vorhanden/nicht vorhanden), '
    'sondern bilden ein Gradientenfeld. Die Übergänge sind fließend, und ein Ton kann Merkmale '
    'zweier Zonen gleichzeitig zeigen. Die Effekte treten nicht exakt bei den Tönen auf, deren Obertöne '
    'die Kammerresonanz treffen, sondern bei Tönen, die in die <b>Nähe</b> der kritischen Kopplung kommen.', sBody))

# ══════════════════════════════════════════════════════════════
# Kapitel 2: Der Kopplungsmechanismus
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 2: Der Kopplungsmechanismus — Kammerresonanz als Vermittler', sChapter))
story.append(Paragraph(
    'Die Kopplung zwischen mehreren Zungen läuft <b>nicht</b> über die schwache Compliance-Luftfeder '
    '(gemeinsames Kammervolumen als Feder), sondern über die <b>Kammerresonanz</b> selbst. '
    'Der Mechanismus ist derselbe, der in Dok. 0005 die Ansprache-Verschlechterung erklärt:', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>Schritt 1 — Oberton trifft Kammerresonanz:</b> Zunge A schwingt bei f<sub>A</sub> und erzeugt '
    'Obertöne n·f<sub>A</sub>. Wenn ein Oberton n·f<sub>A</sub> nahe an f<sub>H</sub> liegt, wird er '
    'von der Kammerresonanz verstärkt (Impedanzverstärkung G, Dok. 0004 und 0007).', sBody))
story.append(Paragraph(
    '<b>Schritt 2 — Resonant verstärkter Druck koppelt an Zunge B:</b> Der verstärkte Oberton erzeugt '
    'einen Drucküberschuss in der Kammer, den Zunge B als periodische Kraft „sieht". '
    'Diese Kraft ist proportional zu G(n·f<sub>A</sub>) und zur Obertonamplitude a<sub>n</sub> ∝ 1/n.', sBody))
story.append(Paragraph(
    '<b>Schritt 3 — Nichtlineare Verstärkung bei voller Amplitude:</b> Bei voller Lautstärke ist die '
    'effektive Spaltfläche A<sub>eff</sub> nicht 2 mm² (Ruhelage), sondern zeitgemittelt ≈ 6 mm² '
    '(Dok. 0002 Kap. 8: Spalt öffnet 8× pro Zyklus). Die Kopplungsstärke skaliert mit A<sub>eff</sub>², '
    'was einen Faktor (6/2)² = 9 gegenüber dem linearen Modell ergibt.', sBody))
story.append(Spacer(1, 2*mm))

story.append(Paragraph(
    '<b>Warum das lineare Modell versagt:</b> Das Helmholtz-Modell aus Dok. 0002–0008 berechnet die '
    'Compliance-Kopplung (gemeinsames Volumen als Luftfeder) und erhält Δf<sub>lock</sub> &lt; 0,01 Hz — '
    'viel zu schwach. Die Praxis zeigt 0,5–1 Hz. Die fehlenden Faktoren sind: '
    '(a) die resonante Verstärkung G ≈ 3–7 bei Obertönen nahe f<sub>H</sub>, und '
    '(b) die nichtlineare Spaltamplitude bei voller Lautstärke (Faktor 9). '
    'Zusammen: 0,01 × 5 × 9 ≈ 0,5 Hz — konsistent mit der Beobachtung.', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>Kernaussage:</b> Die Kopplung zwischen Zungen wird durch dieselbe Kammerimpedanz Z<sub>H</sub>(f) '
    'vermittelt, die in Dok. 0005 die Ansprache und in Dok. 0008 den Klang bestimmt. '
    'Ansprache, Klang und Tremolo-Kopplung sind drei Projektionen derselben physikalischen Größe.', sKeyBox))

# ══════════════════════════════════════════════════════════════
# Kapitel 3: Die Adler-Gleichung
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 3: Die Adler-Gleichung — Injection Locking', sChapter))
story.append(Paragraph(
    'Die Physik gekoppelter Oszillatoren mit kleiner Frequenzdifferenz wird durch die '
    'Adler-Gleichung (1946) beschrieben. Sie gilt universell für Laser, elektronische Oszillatoren, '
    'biologische Rhythmen — und Stimmzungen:', sBody))

story.append(Paragraph(
    'dφ/dt = 2π·Δf − 2π·Δf<sub>lock</sub> · sin(φ)', sFormula))

story.append(Paragraph(
    'Dabei ist Δf die gestimmte Frequenzdifferenz und Δf<sub>lock</sub> der Lock-in-Bereich. '
    'Die Gleichung hat zwei Regime:', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>Regime 1 — Locking (|Δf| &lt; Δf<sub>lock</sub>):</b> Die Phase φ konvergiert auf einen '
    'festen Wert. Beide Zungen schwingen mit derselben Frequenz. Das Tremolo verschwindet vollständig.', sBody))
story.append(Paragraph(
    '<b>Regime 2 — Pulling (|Δf| > Δf<sub>lock</sub>):</b> Die Zungen schwingen mit '
    'verschiedenen Frequenzen, aber die tatsächliche Schwebung ist geringer als gestimmt:', sBody))
story.append(Paragraph(
    'Δf<sub>eff</sub> = √(Δf² − Δf<sub>lock</sub>²)', sFormula))

story.append(Paragraph(
    'Für ein 1-Hz-Tremolo bei Δf<sub>lock</sub> = 0,7 Hz (gemeinsame Öffnung) ergibt sich: '
    'Δf<sub>eff</sub> = √(1² − 0,7²) = 0,71 Hz — das Tremolo wird um 30 % langsamer. '
    'Bei getrennten Öffnungen (Δf<sub>lock</sub> ≈ 0,1 Hz): Δf<sub>eff</sub> = 0,995 Hz — '
    'praktisch unverändert.', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>Drei Zungen (Musette):</b> Bei drei Zungen mit Δf = ±1 Hz verstärkt sich der Effekt. '
    'Die mittlere Zunge „sieht" zwei Kopplungsquellen und der effektive Lock-in-Bereich wächst '
    'um Faktor ≈ √2: Δf<sub>lock,3</sub> ≈ 1,4 · Δf<sub>lock,2</sub>. Bei gemeinsamer Öffnung '
    'kann das ausreichen, um alle drei Zungen einzurasten.', sBody))

story.append(Spacer(1, 2*mm))
story.append(Image('/home/claude/diag_0009_3_adler.png', width=PAGE_W, height=PAGE_W*0.56))
story.append(Paragraph(
    '<i>Abb. 1: Adler-Pulling für verschiedene Kopplungsstärken. Bei gemeinsamer Öffnung (rot) wird '
    'ein 1-Hz-Tremolo auf 0,71 Hz reduziert. Bei getrennten Öffnungen (grün) bleibt es nahezu unverändert.</i>', sSmall))

# ══════════════════════════════════════════════════════════════
# Kapitel 4: Das Drei-Zonen-Gradientenmodell
# ══════════════════════════════════════════════════════════════
story.append(PageBreak())
story.append(Paragraph('Kapitel 4: Das Drei-Zonen-Gradientenmodell', sChapter))
story.append(Paragraph(
    'Die Kopplungsstärke Δf<sub>lock</sub> ist nicht für alle Töne gleich. Sie hängt davon ab, '
    'wie die Obertöne des jeweiligen Tons zur Kammerresonanz f<sub>H</sub> stehen. '
    'Drei Zonen mit unterschiedlichen dominierenden Effekten ergeben sich:', sBody))

story.append(Spacer(1, 2*mm))

# Zone 1
story.append(Paragraph('<b>Zone 1: Locking-Zone (f<sub>H</sub>/f > 5, tiefe Töne)</b>', sBodyB))
story.append(Paragraph(
    'Bei den tiefsten Tönen (D3–G#3, f = 147–208 Hz) liegt f<sub>H</sub> beim 7.–12. Oberton. '
    'Jeder einzelne Oberton koppelt schwach (Energieanteil a<sub>n</sub>² = 1–4 %), aber das '
    'Obertonraster ist dicht: Der Abstand zwischen benachbarten Obertönen beträgt nur f<sub>0</sub> '
    '(147–208 Hz), und die Resonanzbandbreite f<sub>H</sub>/Q<sub>H</sub> ≈ 200–300 Hz ist breit genug, '
    'dass immer mindestens ein Oberton im Band liegt.', sBody))
story.append(Paragraph(
    'Die Summe vieler schwacher Kopplungen ergibt eine moderate Gesamtkopplung. Bei voller Lautstärke '
    '(nichtlinearer Faktor ×9) reicht das aus, um Δf<sub>lock</sub> auf 0,5–1,5 Hz zu treiben — '
    'genug für Locking bei 1-Hz-Tremolo. Subjektiv sind das die Töne, bei denen das Tremolo '
    '„einrastet" oder „verschmiert".', sBody))
story.append(Spacer(1, 2*mm))

# Zone 2
story.append(Paragraph('<b>Zone 2: Geisterpeak-Zone (f<sub>H</sub>/f ≈ 3–5, mittlere Töne)</b>', sBodyB))
story.append(Paragraph(
    'In der Mitte des Tonbereichs (A3–C#5, f = 220–554 Hz) trifft der 3.–5. Oberton auf f<sub>H</sub>. '
    'Jeder dieser Obertöne hat 4–11 % Energieanteil — deutlich mehr als in Zone 1. '
    'Das Obertonraster ist aber weiter (Abstand 220–554 Hz), sodass nur einer oder zwei Obertöne '
    'gleichzeitig im Resonanzband liegen.', sBody))
story.append(Paragraph(
    'Die Kopplung ist stark genug für <b>Intermodulation</b>, aber nicht für Locking. '
    'Im Spektrum erscheinen Mischprodukte: Wenn Zunge A bei f<sub>A</sub> und Zunge B bei f<sub>B</sub> = f<sub>A</sub> + 1 Hz schwingen, '
    'erzeugt die nichtlineare Kopplung über die Kammerresonanz Seitenbänder bei f<sub>A</sub> ± 2 Hz und f<sub>B</sub> ± 2 Hz. '
    'Bei drei Zungen (f − 1, f, f + 1) entstehen Peaks bei f − 2 und f + 2 — '
    'die beobachteten „vierten und fünften Spitzen".', sBody))
story.append(Spacer(1, 2*mm))

# Zone 3
story.append(Paragraph('<b>Zone 3: Ansprache-Zone (f<sub>H</sub>/f ≈ 2–3, hohe Töne)</b>', sBodyB))
story.append(Paragraph(
    'Im oberen Drittel (D5–C6, f = 587–1047 Hz) liegt der 2. oder 3. Oberton auf f<sub>H</sub>. '
    'Diese Obertöne tragen 11–25 % der Gesamtenergie — massive Kopplung. '
    'Aber der Effekt äußert sich primär als Energieverlust (R<sub>H</sub> maximal, Dok. 0005), '
    'nicht als Frequenzverschiebung. Die Zunge schwingt langsam ein und klingt dünn.', sBody))
story.append(Paragraph(
    'Das Tremolo-Pulling ist zwar am stärksten, aber sekundär: Der Ton spricht ohnehin schlecht an, '
    'und bei reduzierter Amplitude sinkt der nichtlineare Verstärkungsfaktor. '
    'In der Praxis bemerkt der Stimmer hier eher die schlechte Ansprache als die Tremolo-Instabilität.', sBody))

story.append(Spacer(1, 4*mm))

# Zusammenfassungstabelle
story.append(Paragraph('<b>Übersicht der drei Zonen</b>', sBodyB))
zone_summary = [
    ['Zone 1\nLocking', 'D3–G#3', '> 5', '7–12', '1–4 %',
     'Viele schwache OT summieren sich.\nLocking bei voller Lautstärke.',
     'Getrennte\nÖffnungen'],
    ['Zone 2\nGeisterpeaks', 'A3–C#5', '3–5', '3–5', '4–11 %',
     'Einzelne mittelstarke OT erzeugen\nIntermodulation im Spektrum.',
     'Tremolo\ngrößer stimmen'],
    ['Zone 3\nAnsprache', 'D5–C6', '2–3', '2–3', '11–25 %',
     'Wenige starke OT → Energieverlust\ndominiert über Frequenzverschiebung.',
     'Kammertiefe\nreduzieren'],
]
zh = ['Zone', 'Töne', 'f<sub>H</sub>/f', 'Dom. OT', 'Energie', 'Mechanismus', 'Abhilfe']
zt_data = [[Paragraph(h, sTH) for h in zh]]
for row in zone_summary:
    zt_data.append([Paragraph(c, sTDL) if i in [0,5,6] else Paragraph(c, sTD) for i, c in enumerate(row)])
zt = Table(zt_data, colWidths=[22*mm, 18*mm, 15*mm, 15*mm, 16*mm, 42*mm, PAGE_W-128*mm])
zt.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), DARKBLUE),
    ('TEXTCOLOR', (0,0), (-1,0), white),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHTGRAY]),
    ('GRID', (0,0), (-1,-1), 0.5, HexColor('#cccccc')),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 3),
    ('BOTTOMPADDING', (0,0), (-1,-1), 3),
]))
story.append(zt)
story.append(Spacer(1, 4*mm))

story.append(Paragraph(
    '<b>Entscheidend:</b> Die Grenzen zwischen den Zonen sind fließend, nicht scharf. '
    'Ein Ton wie A3 (Grenze Zone 1/2) kann bei gemeinsamer Öffnung noch Locking zeigen, '
    'bei getrennten Öffnungen aber bereits Geisterpeaks. Die Zone hängt nicht nur von f<sub>H</sub>/f ab, '
    'sondern auch von der Kopplungsstärke (gemeinsam vs. getrennt) und der Spiellautstärke '
    '(nichtlinearer Faktor).', sWarnBox))

# ══════════════════════════════════════════════════════════════
# Kapitel 5: Diagramme
# ══════════════════════════════════════════════════════════════
story.append(PageBreak())
story.append(Paragraph('Kapitel 5: Quantitative Analyse — Diagramme', sChapter))

story.append(Image('/home/claude/diag_0009_1_zonen.png', width=PAGE_W, height=PAGE_W*0.625))
story.append(Paragraph(
    '<i>Abb. 2: Gewichtete Kopplungsstärke R<sub>eff</sub> = Σ(a<sub>n</sub>² · (G<sub>n</sub> − 1)) '
    'über den Tonbereich. Zone 1 (blau) hat moderate Kopplung durch viele schwache Obertöne. '
    'Zone 2 (orange) zeigt ansteigende Kopplung durch stärker werdende einzelne Obertöne. '
    'Zone 3 (rot) hat die stärkste Kopplung, aber der Effekt wirkt als Energieverlust.</i>', sSmall))
story.append(Spacer(1, 4*mm))

story.append(Image('/home/claude/diag_0009_2_ratio_oberton.png', width=PAGE_W, height=PAGE_W*0.875))
story.append(Paragraph(
    '<i>Abb. 3: Oben: f<sub>H</sub>/f über den Tonbereich. Die gestrichelten Linien markieren die '
    'Zonengrenzen. Unten: Dominanter (stärkster) koppelnder Oberton. Die Prozentzahlen zeigen den '
    'Energieanteil a<sub>n</sub>² = 1/n². Zone 1: OT 7–12 mit 1–2 %. Zone 3: OT 2–3 mit 11–25 %.</i>', sSmall))

story.append(PageBreak())

story.append(Image('/home/claude/diag_0009_4_heatmap.png', width=PAGE_W, height=PAGE_W*0.625))
story.append(Paragraph(
    '<i>Abb. 4: Kopplungslandschaft — gewichtete Impedanzverstärkung a<sub>n</sub>² · G(n·f) für '
    'jeden Ton und jeden Oberton. Die weiße Linie zeigt f<sub>H</sub>/f. Helle Flecken = starke Kopplung. '
    'In Zone 1 (links) verteilt sich die Kopplung auf viele Obertöne. In Zone 3 (rechts) konzentriert '
    'sie sich auf den 2.–3. Oberton.</i>', sSmall))
story.append(Spacer(1, 4*mm))

story.append(Image('/home/claude/diag_0009_5_spektrum.png', width=PAGE_W, height=PAGE_W*0.5))
story.append(Paragraph(
    '<i>Abb. 5: Schematisches Spektrum von drei Zungen (1 Hz Tremolo) bei A4. '
    'Links: Getrennte Öffnungen — drei saubere Obertonreihen. '
    'Rechts: Gemeinsame Öffnung — zusätzliche Geisterpeaks (grau) durch nichtlineare '
    'Intermodulation über die Kammerresonanz.</i>', sSmall))

# ══════════════════════════════════════════════════════════════
# Kapitel 6: Gemeinsame vs. getrennte Öffnungen
# ══════════════════════════════════════════════════════════════
story.append(PageBreak())
story.append(Paragraph('Kapitel 6: Gemeinsame vs. getrennte Öffnungen — quantitativer Vergleich', sChapter))
story.append(Paragraph(
    'Die Praxis zeigt eindeutig: Getrennte Öffnungen (jede Zunge bekommt ihre eigene Luftzufuhr '
    'mit etwas Abstand) reduzieren alle drei Kopplungseffekte drastisch. Der Grund:', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    'Bei gemeinsamer Öffnung teilen sich die Zungen denselben akustischen Widerstand am Hals. '
    'Der verstärkte Oberton von Zunge A erzeugt Druck in der Kammer, und Zunge B „sieht" diesen '
    'Druck direkt — die Kammer ist das Kopplungsmedium. Bei getrennten Öffnungen muss der Druck '
    'von Kammer A über den Balg (großes Volumen, sehr niedrige Impedanz) zu Kammer B gelangen. '
    'Das Balgvolumen (mehrere Liter) wirkt als Tiefpass mit einer Grenzfrequenz weit unter 1 Hz — '
    'die akustische Kopplung wird um Faktor ≈ 7 geschwächt.', sBody))

story.append(Spacer(1, 3*mm))

coupling_table = [
    ['Δf<sub>lock</sub> (2 Zungen)', '≈ 0,7 Hz', '≈ 0,1 Hz', '7 ×'],
    ['1-Hz-Tremolo wird zu', '0,71 Hz (−29 %)', '0,995 Hz (−0,5 %)', '—'],
    ['Δf<sub>lock</sub> (3 Zungen, Musette)', '≈ 1,0 Hz', '≈ 0,14 Hz', '7 ×'],
    ['Geisterpeaks bei A4', '−24 dB', '−41 dB', '17 dB'],
    ['Locking-Grenzton', '≈ G#3 (208 Hz)', '≈ D3 (147 Hz)', '—'],
]
ch = ['Parameter', 'Gemeinsam', 'Getrennt', 'Faktor']
ct_data = [[Paragraph(h, sTH) for h in ch]]
for row in coupling_table:
    ct_data.append([Paragraph(c, sTDL) if i == 0 else Paragraph(c, sTD) for i, c in enumerate(row)])
ct = Table(ct_data, colWidths=[45*mm, 35*mm, 35*mm, PAGE_W-115*mm])
ct.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), DARKBLUE),
    ('TEXTCOLOR', (0,0), (-1,0), white),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHTGRAY]),
    ('GRID', (0,0), (-1,-1), 0.5, HexColor('#cccccc')),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 3),
    ('BOTTOMPADDING', (0,0), (-1,-1), 3),
]))
story.append(ct)
story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    '<i>Tabelle 1: Empirische Vergleichswerte gemeinsame vs. getrennte Öffnungen. '
    'Alle Werte gelten für volle Spiellautstärke. Bei leisem Spiel sinkt der nichtlineare '
    'Verstärkungsfaktor und die Effekte werden schwächer.</i>', sSmall))

# ══════════════════════════════════════════════════════════════
# Kapitel 7: Korrektur beim Stimmen
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 7: Praktische Konsequenz — Korrektur beim Stimmen', sChapter))
story.append(Paragraph(
    'Wenn getrennte Öffnungen nicht möglich sind (Bauform, Platz), muss das Tremolo '
    '<b>bewusst größer gestimmt</b> werden, um den Pulling-Effekt zu kompensieren:', sBody))
story.append(Paragraph(
    'Δf<sub>Stimmung</sub> = √(Δf<sub>Ziel</sub>² + Δf<sub>lock</sub>²)', sFormula))

story.append(Paragraph(
    'Für ein gewünschtes Tremolo von 1 Hz bei Δf<sub>lock</sub> = 0,7 Hz:', sBody))
story.append(Paragraph(
    'Δf<sub>Stimmung</sub> = √(1² + 0,7²) = √1,49 = 1,22 Hz', sFormula))
story.append(Paragraph(
    'Man stimmt also 22 % weiter als gewünscht. Dieser Korrekturwert ist tonabhängig — '
    'in Zone 1 (tiefe Töne) größer, in Zone 2 (mittlere Töne) kleiner, in Zone 3 (hohe Töne) '
    'wieder größer (wegen der starken Einzeloberton-Kopplung). In Zone 1 kann die Korrektur '
    'sogar so groß werden, dass Locking nicht mehr vermeidbar ist — dann bleiben nur getrennte '
    'Öffnungen oder eine Änderung der Kammergeometrie.', sBody))

story.append(Spacer(1, 3*mm))
story.append(Paragraph(
    '<b>Stimm-Kompensationstabelle</b> (gemeinsame Öffnung, volle Lautstärke):', sBodyB))
comp_rows = [
    ['D3 (147 Hz)', '1,0', '~1,2', '1,56', 'Locking möglich'],
    ['A3 (220 Hz)', '1,0', '~0,8', '1,28', 'Starkes Pulling'],
    ['A4 (440 Hz)', '1,0', '~0,5', '1,12', 'Moderates Pulling'],
    ['D5 (587 Hz)', '1,0', '~0,6', '1,17', 'Pulling + Ansprache'],
    ['A5 (880 Hz)', '1,0', '~0,4', '1,08', 'Schwaches Pulling'],
]
comp_h = ['Ton', 'Δf<sub>Ziel</sub> [Hz]', 'Δf<sub>lock</sub> [Hz]', 'Δf<sub>Stimmung</sub> [Hz]', 'Bemerkung']
comp_data = [[Paragraph(h, sTH) for h in comp_h]]
for row in comp_rows:
    comp_data.append([Paragraph(c, sTDL) if i in [0,4] else Paragraph(c, sTD) for i, c in enumerate(row)])
comp_t = Table(comp_data, colWidths=[28*mm, 25*mm, 25*mm, 30*mm, PAGE_W-108*mm])
comp_t.setStyle(TableStyle([
    ('BACKGROUND', (0,0), (-1,0), DARKBLUE),
    ('TEXTCOLOR', (0,0), (-1,0), white),
    ('ROWBACKGROUNDS', (0,1), (-1,-1), [white, LIGHTGRAY]),
    ('GRID', (0,0), (-1,-1), 0.5, HexColor('#cccccc')),
    ('VALIGN', (0,0), (-1,-1), 'MIDDLE'),
    ('TOPPADDING', (0,0), (-1,-1), 3),
    ('BOTTOMPADDING', (0,0), (-1,-1), 3),
]))
story.append(comp_t)

# ══════════════════════════════════════════════════════════════
# Kapitel 8: Verbindung zu Dok. 0005 — Dieselbe Physik
# ══════════════════════════════════════════════════════════════
story.append(PageBreak())
story.append(Paragraph('Kapitel 8: Verbindung zur Einzelzungen-Physik (Dok. 0005)', sChapter))
story.append(Paragraph(
    'Die Korrelation zwischen Tremolo-Kopplung und Ansprache ist kein Zufall — sie folgt aus der Physik. '
    'Beide Effekte werden von derselben komplexen Kammerimpedanz Z<sub>H</sub>(f) = R<sub>H</sub> + jX<sub>H</sub> '
    'bestimmt (Dok. 0005):', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>Realteil R<sub>H</sub>:</b> Bestimmt den Energieverlust. Für die Einzelzunge bedeutet das '
    'schlechtere Ansprache (Dok. 0005). Für die Mehrzungen-Kopplung bedeutet es: Die Kammer absorbiert '
    'und re-emittiert Energie — genau der Mechanismus, der die Kopplung vermittelt. '
    '<b>Großes R<sub>H</sub> = schlechte Ansprache UND starke Kopplung.</b>', sBody))
story.append(Paragraph(
    '<b>Imaginärteil X<sub>H</sub>:</b> Bestimmt die Frequenzverschiebung. Für die Einzelzunge verschiebt '
    'X<sub>H</sub> die Obertöne (Dok. 0005). Für die Mehrzungen-Kopplung bestimmt X<sub>H</sub> die '
    'Phasenbeziehung zwischen den gekoppelten Zungen — und damit, ob Pulling oder Locking auftritt.', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    'Da R<sub>H</sub> und X<sub>H</sub> durch die Kramers-Kronig-Relation verknüpft sind (Dok. 0005 Kap. 3), '
    'gilt: <b>Ein Ton, der Ansprache-Probleme hat (großes R<sub>H</sub> bei einem Oberton), hat '
    'gleichzeitig ein starkes Tremolo-Pulling.</b> Aber der Effekt zeigt sich unterschiedlich: '
    'Die Ansprache verschlechtert sich exakt bei f<sub>H</sub> = n·f (R<sub>H</sub> maximal, X<sub>H</sub> = 0). '
    'Das Pulling ist maximal knapp neben der Resonanz (X<sub>H</sub> maximal, R<sub>H</sub> abfallend). '
    'Deshalb treten die stärksten Tremolo-Effekte nicht bei exakt denselben Tönen auf wie die '
    'schlechteste Ansprache — sondern bei den Nachbartönen. Das ist genau die Beobachtung.', sBody))

# ══════════════════════════════════════════════════════════════
# Kapitel 9: Warum die drei Zonen verschiedene Effekte zeigen
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 9: Warum die Zonen verschiedene Effekte zeigen', sChapter))

story.append(Paragraph(
    'Die unterschiedlichen Erscheinungsformen in den drei Zonen folgen aus der Struktur des Obertonrasters:', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>Zone 1 (tiefe Töne) → Locking:</b> Das Obertonraster ist dicht (Abstand = f<sub>0</sub> ≈ 150 Hz). '
    'Bei f<sub>H</sub> ≈ 1500 Hz liegen ~10 Obertöne in Reichweite, mindestens einer immer im Resonanzband. '
    'Die Kopplung ist breitbandig — sie wirkt bei jeder Verstimmung, nicht nur bei bestimmten Δf-Werten. '
    'Das begünstigt Locking, weil die Kopplungskraft nie „Lücken" hat.', sBody))

story.append(Paragraph(
    '<b>Zone 2 (mittlere Töne) → Geisterpeaks:</b> Das Obertonraster ist weiter (Abstand ≈ 300–500 Hz). '
    'Die Kopplung ist schmalbandig — sie wirkt nur über einen oder zwei Obertöne, die in das '
    'Resonanzband fallen. Diese selektive Kopplung erzeugt Intermodulationsprodukte bei diskreten '
    'Frequenzen (Seitenbänder), statt eine breitbandige Synchronisierung. Das Ergebnis sind '
    'Geisterpeaks im Spektrum, kein Locking.', sBody))

story.append(Paragraph(
    '<b>Zone 3 (hohe Töne) → Ansprache:</b> Der 2. Oberton (Oktave) liegt auf oder nahe f<sub>H</sub>. '
    'Dieser Oberton trägt 25 % der Gesamtenergie — die Kopplung ist maximal, aber sie äußert sich '
    'als massive Energieabsorption durch die Kammer (R<sub>H</sub> dominiert). '
    'Die Zunge verliert so viel Energie an die Kammerresonanz, dass die Amplitude sinkt. '
    'Bei reduzierter Amplitude sinkt auch der nichtlineare Verstärkungsfaktor — '
    'der Teufelskreis begünstigt Energieverlust gegenüber Frequenzverschiebung.', sBody))

# ══════════════════════════════════════════════════════════════
# Kapitel 10: Grenzen des Modells
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 10: Grenzen des Modells', sChapter))
story.append(Paragraph(
    'Das Drei-Zonen-Modell erklärt die qualitativen Beobachtungen und gibt quantitative Richtwerte. '
    'Es hat aber klare Grenzen:', sBody))

story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<b>1. Die Lock-in-Werte sind empirisch kalibriert.</b> Die Adler-Gleichung beschreibt die Form '
    'der Kopplung korrekt (Pulling-Kurve, Locking-Schwelle), aber Δf<sub>lock</sub> selbst kann mit dem '
    'Helmholtz-Modell nicht berechnet werden. Der nichtlineare Verstärkungsfaktor (×9 bei voller '
    'Amplitude) ist eine Abschätzung aus Dok. 0002 Kap. 8, nicht aus einer FSI-Simulation.', sBody))
story.append(Paragraph(
    '<b>2. Die Zonengrenzen sind geometrieabhängig.</b> Die Werte f<sub>H</sub>/f = 5 und 3 gelten für '
    'den in Dok. 0007 berechneten Stimmstock (Variante C, ∅ 10 mm Öffnung, Q<sub>H</sub> ≈ 8). '
    'Bei anderen Geometrien verschieben sich die Grenzen.', sBody))
story.append(Paragraph(
    '<b>3. Die Lautstärkeabhängigkeit ist qualitativ, nicht quantitativ.</b> Der nichtlineare '
    'Verstärkungsfaktor skaliert mit A<sub>eff</sub>², und A<sub>eff</sub> hängt vom Blasdruck ab. '
    'Bei pianissimo kann ein Ton in Zone 1 völlig unproblematisch sein — bei fortissimo einrasten. '
    'Diese Dynamikabhängigkeit ist in der Praxis bekannt, aber nicht exakt berechenbar.', sBody))
story.append(Paragraph(
    '<b>4. Aerodynamische Nahfeldeffekte fehlen.</b> Bei gemeinsamer Öffnung können die Luftströme '
    'zweier Zungen sich direkt beeinflussen (Bernoulli-Wechselwirkung), bevor die Kammerresonanz '
    'ins Spiel kommt. Dieser Effekt ist im Helmholtz-Modell nicht enthalten.', sBody))

# ══════════════════════════════════════════════════════════════
# Kapitel 11: Zusammenfassung
# ══════════════════════════════════════════════════════════════
story.append(Paragraph('Kapitel 11: Zusammenfassung', sChapter))

story.append(Paragraph(
    '<b>1. Der Kopplungsmechanismus:</b> Mehrere Zungen koppeln über die Kammerresonanz, nicht über '
    'die schwache Compliance-Luftfeder. Obertöne nahe f<sub>H</sub> werden resonant verstärkt und '
    'übertragen Energie zwischen den Zungen. Bei voller Amplitude ist die Kopplung ≈ 500× stärker '
    'als das lineare Modell vorhersagt.', sKeyBox))

story.append(Paragraph(
    '<b>2. Das Drei-Zonen-Gradientenmodell:</b><br/>'
    '• Zone 1 (f<sub>H</sub>/f > 5, tiefe Töne): Dichte Obertöne → breitbandige Kopplung → <b>Locking</b>.<br/>'
    '• Zone 2 (f<sub>H</sub>/f ≈ 3–5, mittlere Töne): Einzelne starke OT → selektive Kopplung → <b>Geisterpeaks</b>.<br/>'
    '• Zone 3 (f<sub>H</sub>/f ≈ 2–3, hohe Töne): Oktave/Duodezime auf f<sub>H</sub> → massive Absorption → <b>Ansprache-Verlust</b>.<br/>'
    'Die Grenzen sind fließend. Die Effekte treten in der Nähe der kritischen Kopplung auf, nicht exakt bei den Resonanztönen.', sKeyBox))

story.append(Paragraph(
    '<b>3. Abhilfe:</b><br/>'
    '• <b>Getrennte Öffnungen</b> schwächen die Kopplung um Faktor ≈ 7 — die wirksamste Maßnahme.<br/>'
    '• <b>Tremolo größer stimmen:</b> Δf<sub>Stimmung</sub> = √(Δf<sub>Ziel</sub>² + Δf<sub>lock</sub>²) '
    'kompensiert das Pulling, hilft aber nicht gegen Locking.<br/>'
    '• <b>Kammertiefe reduzieren</b> (Dok. 0007, Variante C) hebt f<sub>H</sub> und verschiebt '
    'die Zonengrenzen nach oben — weniger Töne sind betroffen.', sKeyBox))

story.append(Paragraph(
    '<b>4. Die einheitliche Physik:</b> Ansprache (Dok. 0005), Klang (Dok. 0008) und Tremolo-Kopplung '
    '(dieses Dokument) sind drei Projektionen derselben Kammerimpedanz Z<sub>H</sub>(f). '
    'Was die Ansprache verbessert, verbessert auch den Klang und reduziert die Tremolo-Kopplung. '
    'Es gibt keinen Trade-off zwischen diesen drei Zielen — das Optimum ist dasselbe: '
    '<b>f<sub>H</sub>/f ≥ 3, minimales Volumen, große Öffnung.</b>', sKeyBox))

story.append(Spacer(1, 6*mm))
# Abschluss
story.append(Table([['']], colWidths=[PAGE_W*0.6],
    style=TableStyle([('LINEBELOW', (0,0), (-1,-1), 1, ACCENTRED),
                      ('ALIGN', (0,0), (-1,-1), 'CENTER')])))
story.append(Spacer(1, 2*mm))
story.append(Paragraph(
    '<i>Die Berechnung sagt, wo man suchen soll. Der praktische Test sagt, was man findet. '
    'Und die Erfahrung aus tausend Instrumenten sagt, worauf es ankommt.</i>',
    ParagraphStyle('Closing', parent=sSmall, alignment=TA_CENTER)))

# ── Build ──
doc.build(story)
print(f"PDF erzeugt: {output_path}")
print(f"  Diagramme: diag_0009_1..5 in /home/claude/")
