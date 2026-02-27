#!/usr/bin/env python3
"""
Bass-Stimmzunge 50 Hz - v6 VOLLSTÄNDIG NEU
===========================================
KORREKTUREN:
  - Stimmplatte KEILFÖRMIG: 13mm (freies Ende/Klappe) -> 2mm (Einspannung)
  - Aufbiegung 1.5mm am freien Ende (dreieckiger Spalt, nicht 0.4mm)
  - Basszungen VERSCHRAUBT (2 Schrauben), kleine genietet
  - Folie (0.1mm) vs Furnier (0.5mm) als Trennwand berechnet
  - Beide Klappenorientierungen (Ausblas + Umlenk) voll berechnet
"""
import numpy as np
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.colors import HexColor, white, black
from reportlab.lib.enums import TA_CENTER, TA_JUSTIFY
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table,
                                TableStyle, PageBreak, HRFlowable, Image)
from reportlab.graphics.shapes import (Drawing, Rect, Line, String, Polygon, PolyLine)
import shutil

# ═══════════════════════════════════════════════════════════════
# KONSTANTEN
# ═══════════════════════════════════════════════════════════════
c_s = 343.0; rho = 1.204; nu = 1.516e-5

# ═══════════════════════════════════════════════════════════════
# GEOMETRIE
# ═══════════════════════════════════════════════════════════════
L_kam = 0.094; B_kam = 0.048; H_kam = 0.040

# STIMMPLATTE: keilförmig
t_platte_frei = 0.013   # 13mm am freien/Klappe-Ende
t_platte_einsp = 0.002  # 2mm am Einspann-Ende
t_platte_mean = (t_platte_frei + t_platte_einsp) / 2  # 7.5mm Mittel

# Klappe: am Boden, gegenüber freiem Ende
L_klap = 0.040; W_klap = 0.020
alpha_deg = 30.0; alpha_r = np.radians(alpha_deg)
h_open = W_klap * np.sin(alpha_r)
S_klap_geo = L_klap * h_open; C_c = 0.61
S_klap_eff = S_klap_geo * C_c

# Stimmzunge
L_zunge = 0.070; W_zunge = 0.008
h_aufbiegung = 0.0015  # 1.5mm Aufbiegung am freien Ende

# Schlitz in der Platte
W_schlitz = 0.009   # ~9mm (1mm breiter als Zunge)
# Seitenspalt: (9-8)/2 = 0.5mm pro Seite

# SPALT: dreieckig, 0 bei Einspannung, 1.5mm am freien Ende
# Effektive Spaltfläche = Integral über Zungenlänge
# h(x) = h_aufbiegung * x/L_zunge (linear von 0 bis 1.5mm)
# S_gap = Integral0L W_zunge * h(x) dx / L_zunge ... nein:
# Die Durchströmfläche bei jedem x ist W_schlitz * h(x) (quer zur Strömung)
# Für den Gesamtdurchfluss: Q = Integral0L v(x) * W_schlitz * h(x) * dx
# Bei Bernoulli-Strömung: v(x) ~ const (gleicher Druckabfall überall)
# -> Q = v * W_schlitz * Integral0L h(x) dx = v * W_schlitz * L * h_max/2
# -> S_gap_eff = W_schlitz * L_zunge * h_aufbiegung / 2
# ABER: Das wäre 9mm * 70mm * 1.5mm / 2 = 472 mm² -- viel zu gross!

# Richtig: S_gap_eff ist die Fläche NORMAL zur Strömrichtung an der engsten Stelle.
# Die Luft tritt von unten (Kammer) durch den Schlitz nach oben.
# Der Engpass ist die Fläche zwischen Zungenunterkante und Schlitzrand.
# Bei x (Abstand von Einspannung): h(x) = h_aufbiegung * x/L_zunge
# Die LOKALE Durchströmfläche pro Längeneinheit = W_schlitz * h(x)
# (eigentlich W_schlitz minus W_zunge = 1mm Seitenspalt, PLUS W_zunge * h(x) oben)

# Präzise: Bei Position x strömt Luft durch:
# a) Seitenspalt: 2 * 0.5mm * t_platte(x) (Spalt zwischen Zunge und Schlitzrand, durch Plattendicke)
# b) Oberer Spalt: W_zunge * h(x) (über der Zunge, wo sie aufgebogen ist)
# Der Seitenspalt hat hohen Widerstand (schmal, lang), der obere Spalt dominiert wo h(x) gross.

# Vereinfachung: Hauptstrom durch den oberen Spalt (Aufbiegung)
# Am freien Ende (x=L): h=1.5mm, Breite = W_zunge = 8mm -> lokal 12 mm²
# Am Mittelpunkt (x=L/2): h=0.75mm -> lokal 6 mm²
# Am Einspannpunkt (x=0): h=0 -> kein Durchfluss

# Effektive Gesamtfläche (Bernoulli, gleichmässiger Druck):
# S_gap = Integral0L W_zunge * h(x)/L dx = W_zunge * h_max / 2
S_gap = W_zunge * h_aufbiegung / 2  # = 8 * 1.5 / 2 = 6 mm²

# Trennwand: zwei Varianten
t_furnier = 0.0005   # 0.5mm
t_folie = 0.0001     # 0.1mm
L_wall = 0.075
gap_fold = L_kam - L_wall  # 19mm

print("=" * 90)
print("  BASS-STIMMZUNGE 50 Hz - NEUBERECHNUNG v6")
print("  Keilförmige Stimmplatte, 1.5mm Aufbiegung, Folie vs. Furnier")
print("=" * 90)

# ═══════════════════════════════════════════════════════════════
# STIMMPLATTE
# ═══════════════════════════════════════════════════════════════
print(f"\n  STIMMPLATTE (keilförmig):")
print(f"    Am freien Ende (bei Klappe): {t_platte_frei*1e3:.0f} mm dick")
print(f"    Am Einspann-Ende:            {t_platte_einsp*1e3:.0f} mm dick")
print(f"    Mittlere Dicke:              {t_platte_mean*1e3:.1f} mm")
print(f"    -> Der Schlitz ist am freien Ende {t_platte_frei*1e3:.0f} mm tief!")
print(f"       Am Einspann-Ende nur {t_platte_einsp*1e3:.0f} mm tief.")
print(f"    -> Die Platte wirkt als kurzer Kanal fuer die Luft.")

# Schlitz als Kanal: mittlere Tiefe
t_schlitz_mean = t_platte_mean
S_schlitz_mean = W_schlitz * t_schlitz_mean  # Schlitz-Querschnitt

print(f"\n    Schlitz: {W_schlitz*1e3:.0f} mm breit x {t_schlitz_mean*1e3:.1f} mm tief (Mittel)")
print(f"    Schlitz-Querschnitt: {S_schlitz_mean*1e6:.1f} mm²")
print(f"    -> Deutlich groesser als der Aufbiegungs-Spalt ({S_gap*1e6:.1f} mm²)")
print(f"    -> Der Spalt (Aufbiegung) bleibt der Engpass, nicht der Schlitz.")

# ═══════════════════════════════════════════════════════════════
# AUFBIEGUNG UND SPALT
# ═══════════════════════════════════════════════════════════════
print(f"\n  AUFBIEGUNG UND SPALT:")
print(f"    Aufbiegung am freien Ende:  {h_aufbiegung*1e3:.1f} mm")
print(f"    Spaltverlauf: dreieckig (0 bei Einspannung -> {h_aufbiegung*1e3:.1f} mm am Ende)")
print(f"    Lokale Spaltfläche am freien Ende: {W_zunge*1e3:.0f} x {h_aufbiegung*1e3:.1f} = {W_zunge*h_aufbiegung*1e6:.0f} mm²")
print(f"    Effektive Gesamtfläche (Integral): S_gap = W x h_max / 2 = {S_gap*1e6:.1f} mm²")
print(f"    (Zum Vergleich: Klappe eff. = {S_klap_eff*1e4:.2f} cm² = {S_klap_eff*1e6:.0f} mm²)")
print(f"    Verhältnis Klappe/Spalt: {S_klap_eff/S_gap:.0f}")

# ═══════════════════════════════════════════════════════════════
# EULER-BERNOULLI
# ═══════════════════════════════════════════════════════════════
E_st = 200e9; rho_st = 7800; lam1 = 1.8751
h_reed = 50*2*np.pi*L_zunge**2*np.sqrt(12)/(lam1**2*np.sqrt(E_st/rho_st))
I_reed = W_zunge*h_reed**3/12
m_total = rho_st*L_zunge*W_zunge*h_reed
m_eff = 0.2427*m_total
k_eff = 3*E_st*I_reed/L_zunge**3
f_check = lam1**2/(2*np.pi*L_zunge**2)*np.sqrt(E_st*I_reed/(rho_st*W_zunge*h_reed))

# Schwellendruck: mit 1.5mm Aufbiegung
# Die Bernoulli-Kraft wirkt am stärksten dort wo Spalt eng ist (nahe Einspannung)
# Aber dort ist der Hebelarm klein. Die effektive Kraft am Tip:
# dp_min ~ k_eff * h_aufbiegung / (W_zunge * L_zunge * 0.5)
# Faktor 0.5 weil Druckverteilung dreieckig (wirkt mehr am Ende)
dp_min = k_eff * h_aufbiegung / (W_zunge * L_zunge * 0.5)

print(f"\n  STIMMZUNGE (Euler-Bernoulli):")
print(f"    Dicke:      {h_reed*1e3:.3f} mm")
print(f"    Verif.:     f1 = {f_check:.2f} Hz OK")
print(f"    k_eff:      {k_eff:.2f} N/m")
print(f"    m_eff:      {m_eff*1e6:.1f} mg ({m_eff*1e3:.3f} g)")
print(f"    Befestigung: verschraubt (2 Schrauben) auf Stimmstock")
print(f"\n    Schwellendruck (Aufbiegung {h_aufbiegung*1e3:.1f} mm):")
print(f"    DeltaP_min = k*h / (W*L*0.5) = {k_eff:.2f}*{h_aufbiegung*1e3:.1f}e-3 / ({W_zunge*1e3:.0f}e-3*{L_zunge*1e3:.0f}e-3*0.5)")
print(f"           = {dp_min:.1f} Pa")
print(f"    -> Immer noch niedrig gegenüber 500-5000 Pa Balgdruck")

# Einschwingzeit
zeta_vals = [0.005, 0.01, 0.02, 0.04]
Q_vals = [1/(2*z) for z in zeta_vals]
tau_vals = [Q/(np.pi*50) for Q in Q_vals]

print(f"\n  EINSCHWINGZEIT:")
for z, Q, tau in zip(zeta_vals, Q_vals, tau_vals):
    print(f"    zeta={z:.3f}  Q={Q:>5.0f}  tau={tau*1e3:>6.0f} ms  ({Q:.0f} Perioden à 20 ms)")

# ═══════════════════════════════════════════════════════════════
# KAMMERVOLUMEN
# ═══════════════════════════════════════════════════════════════
V_brutto = L_kam * B_kam * H_kam

print(f"\n  KAMMERVOLUMEN:")
print(f"    Brutto: {V_brutto*1e6:.1f} cm³")

# ═══════════════════════════════════════════════════════════════
# TRENNWAND-MATERIAL: FOLIE vs FURNIER
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*90}")
print(f"  TRENNWAND-MATERIAL: FOLIE ({t_folie*1e3:.1f} mm) vs. FURNIER ({t_furnier*1e3:.1f} mm)")
print(f"{'='*90}")

for mat_name, t_w, roughness, desc in [
    ("Folie (Kunststoff)", t_folie, 0.001, "glatt, flexibel, duenn"),
    ("Furnier (Holz)", t_furnier, 0.05, "rauer, steif, dicker"),
]:
    W_k = (B_kam - t_w) / 2  # Kanalbreite bei gerader Wand
    A_k = W_k * H_kam
    V_w = L_wall * t_w * H_kam
    V_netto = V_brutto - V_w
    # Hydraulischer Durchmesser
    D_h = 4*A_k / (2*(W_k + H_kam))
    # Darcy-Reibungsfaktor (laminare Näherung für glatt, turbulent für rau)
    # Bei Re_kanal ~ 10-50 (sehr niedrig): f = 64/Re
    # Relative Rauigkeit
    eps_rel = roughness / (D_h * 1e3)
    
    print(f"\n  {mat_name} ({desc}):")
    print(f"    Wanddicke:        {t_w*1e3:.1f} mm")
    print(f"    Kanalbreite:      {W_k*1e3:.2f} mm (bei gerader Wand)")
    print(f"    Kanalquerschnitt: {A_k*1e4:.2f} cm²")
    print(f"    Netto-Volumen:    {V_netto*1e6:.1f} cm³ (- {V_w*1e6:.2f} cm³ Wand)")
    print(f"    D_h:              {D_h*1e3:.1f} mm")
    print(f"    Rauigkeit:        ~{roughness:.3f} mm -> epsilon/D = {eps_rel:.4f}")
    print(f"    Reibungsfaktor:   f ~ {'0.02 (glatt)' if t_w < 0.001 else '0.03 (rau)'}")

    # Besonderheit Folie: kann flattern!
    if t_w < 0.001:
        print(f"    ACHTUNG: FLUTTER: Folie kann bei Stroemung vibrieren (Eigenfrequenz sehr hoch)")
        print(f"      -> kann Nebengeraeusche erzeugen")
        print(f"      -> muss straff gespannt oder an Stuetzpunkten fixiert werden")

print(f"\n  UNTERSCHIED Folie vs. Furnier:")
print(f"    Kanalbreite:  {(B_kam-t_folie)/2*1e3:.2f} mm (Folie) vs. {(B_kam-t_furnier)/2*1e3:.2f} mm (Furnier)")
print(f"    Differenz:    {((t_furnier-t_folie)/2)*1e3:.2f} mm pro Kanal -> vernachlässigbar")
print(f"    Oberfläche:   Folie glatt -> weniger Reibung, Furnier rau -> mehr Reibung")
print(f"    Stabilität:   Furnier formstabil, Folie braucht Spannung/Stütze")
print(f"    Akustik:      Folie kann als Membran schwingen (Nebengeräusche)")

# ═══════════════════════════════════════════════════════════════
# KLAPPENORIENTIERUNG: AUSBLAS vs. UMLENK
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*90}")
print(f"  KLAPPENORIENTIERUNG")
print(f"{'='*90}")

# Ausblas: Strahl zeigt in Kammer
zeta_klap_ausblas = (1/C_c - 1)**2 + np.sin(alpha_r)**2
# -> Vena contracta + Richtungsänderung

# Umlenk: Strahl zeigt weg, muss 150° umkehren
# Zusätzlicher Verlust: Umlenkung um scharfe Kante
# zeta_umlenk ~ 1.0-1.5 für 150° Umlenkung um scharfe Kante (Idelchik)
zeta_reversal = 1.2  # 150°-Umlenkung
zeta_klap_umlenk = zeta_klap_ausblas + zeta_reversal

print(f"\n  AUSBLAS-RICHTUNG (Strahl in die Kammer):")
print(f"    Strahl zeigt mit 30° direkt in Kanal B")
print(f"    Horizontalkomponente: cos(30°) = {np.cos(alpha_r):.3f} -> {np.cos(alpha_r)*100:.0f}% der Geschwindigkeit")
print(f"    Vertikalkomponente:   sin(30°) = {np.sin(alpha_r):.3f} -> {np.sin(alpha_r)*100:.0f}%")
print(f"    zeta_Klappe = (1/Cc-1)² + sin²(alpha) = ({1/C_c:.2f}-1)² + {np.sin(alpha_r)**2:.3f}")
print(f"             = {(1/C_c-1)**2:.3f} + {np.sin(alpha_r)**2:.3f} = {zeta_klap_ausblas:.3f}")

print(f"\n  UMLENK-RICHTUNG (Strahl weg von Kammer):")
print(f"    Strahl zeigt mit 30° nach AUSSEN (weg von Kammer)")
print(f"    Luft muss ~150° um Klappenkante kehrtmachen")
print(f"    Wie ein U-Turn auf einer Strasse: hoher Energieverlust")
print(f"    zeta_Klappe = zeta_Ausblas + zeta_Umlenkung = {zeta_klap_ausblas:.3f} + {zeta_reversal:.2f} = {zeta_klap_umlenk:.3f}")
print(f"    -> {zeta_klap_umlenk/zeta_klap_ausblas:.1f}-facher Verlust gegenueber Ausblas!")

print(f"\n  WIRKUNG auf die Stroemung:")
print(f"    Ausblas: Gerichteter Strahl trifft auf Trennwand -> Wandform wirkt voll")
print(f"    Umlenk:  Breiter, ungerichteter Eintritt -> Wandform wirkt weniger")
print(f"    Ausblas: Hohe Geschwindigkeit am Eintritt -> Coanda/Düse wirken")
print(f"    Umlenk:  Niedrigere Geschwindigkeit -> mehr Turbulenz am Eintritt")
print(f"    Umlenk:  Potentiell gleichmässigere Verteilung (breiterer Strahl)")
print(f"    Umlenk:  Aber: mehr Gesamtverlust und mehr Nebengeräusche")

# ═══════════════════════════════════════════════════════════════
# VERLUSTBERECHNUNG - ALLE KOMBINATIONEN
# ═══════════════════════════════════════════════════════════════
print(f"\n{'='*90}")
print(f"  VERLUSTBERECHNUNG - 12 KOMBINATIONEN")
print(f"  (3 Wandformen x 2 Materialien x 2 Klappenrichtungen... -> vereinfacht)")
print(f"{'='*90}")

# Verwende Furnier für die Hauptberechnung, zeige Folie-Differenz
t_w = t_furnier  # Referenz
W_k = (B_kam - t_w) / 2
A_k = W_k * H_kam
S_fold = gap_fold * H_kam

# Schräge/Parabel: gleiche Geometrie wie v5
frac_klap = 0.30
W_B_unten_klap = B_kam * frac_klap - t_w/2
W_B_unten_falt = B_kam * 0.50 - t_w/2
W_B_oben_klap = B_kam * (1-frac_klap) - t_w/2
W_B_oben_falt = B_kam * 0.50 - t_w/2
beta_B = np.degrees(np.arctan(abs(W_B_unten_falt-W_B_unten_klap)/L_wall))

A_B_klap = W_B_unten_klap * H_kam
A_B_falt = W_B_unten_falt * H_kam
A_A_falt = W_B_oben_falt * H_kam
A_A_klap = W_B_oben_klap * H_kam

# Verlust-Definitionen für alle Kombinationen
configs = {}

for orient_name, zeta_klap in [("Ausblas", zeta_klap_ausblas), ("Umlenk", zeta_klap_umlenk)]:
    for wall_name, losses_internal in [
        ("A-Gerade", [
            ('Eintritt Kanal B', 0.5, A_k),
            ('Aufprall 90° auf Wand', 0.8, A_k),
            ('180°-Faltung', 1.0, S_fold),
        ]),
        ("B-Schraeg", [
            ('Duesen-Eintritt', 0.15, A_B_klap),
            ('Diffusor Kanal B', 0.10, (A_B_klap+A_B_falt)/2),
            ('180°-Faltung', 0.8, S_fold),
            ('Diffusor Kanal A', 0.10, (A_A_falt+A_A_klap)/2),
        ]),
        ("C-Parabel", [
            ('Coanda-Eintritt', 0.12, A_B_klap),
            ('Parabel-Kanal B', 0.08, (A_B_klap+A_B_falt)/2),
            ('180°-Faltung', 0.6, S_fold),
            ('Kanal A (weit)', 0.08, (A_A_falt+A_A_klap)/2),
        ]),
    ]:
        cname = f"{wall_name} + {orient_name}"
        all_losses = [('Klappe', zeta_klap, S_klap_eff)] + losses_internal + [('Spaltaustritt', 1.0, S_gap)]
        
        zeta_total = 0
        for desc, z, S in all_losses:
            zeta_total += z * (S_gap/S)**2
        
        # Kanalreibung
        A_mean = A_k if 'Gerade' in wall_name else (A_B_klap+A_B_falt+A_A_falt+A_A_klap)/4
        D_h = 4*A_mean / (2*(A_mean/H_kam + H_kam))
        L_path = 2*L_wall + np.pi*W_k/2
        f_darcy = 0.03
        z_reib = f_darcy * L_path / D_h * (S_gap/A_mean)**2
        zeta_total += z_reib
        
        configs[cname] = {'zeta': zeta_total, 'losses': all_losses, 'orient': orient_name, 'wall': wall_name}

# Berechne Strömung für alle
pressures = [200, 500, 1000, 2000, 5000]
print(f"\n  Ergebnisse bei 1000 Pa:")
print(f"  {'Konfiguration':<28} {'zeta_total':>10} {'v_Spalt':>10} {'Q [ml/s]':>10} {'Re_Spalt':>10}")
print(f"  {'-'*70}")

for cname in sorted(configs.keys()):
    cfg = configs[cname]
    z = cfg['zeta']
    v = np.sqrt(2*1000/(rho*(1+z)))
    Q = v * S_gap
    # Re basiert auf mittlerer Spalthöhe (h_aufbiegung/2 = 0.75mm)
    Re = v * (h_aufbiegung/2) / nu
    cfg['v1000'] = v
    cfg['Q1000'] = Q
    cfg['Re1000'] = Re
    print(f"  {cname:<28} {z:>10.6f} {v:>10.2f} {Q*1e6:>10.1f} {Re:>10.0f}")

# Vergleich Folie vs Furnier
print(f"\n  FOLIE vs. FURNIER Differenz:")
# Folie: dünnere Wand -> breitere Kanäle -> geringfügig weniger Reibung
W_k_folie = (B_kam - t_folie) / 2
W_k_furnier = (B_kam - t_furnier) / 2
print(f"    Kanalbreite Folie:   {W_k_folie*1e3:.2f} mm")
print(f"    Kanalbreite Furnier: {W_k_furnier*1e3:.2f} mm")
print(f"    Differenz:           {(W_k_folie-W_k_furnier)*1e3:.2f} mm ({(W_k_folie-W_k_furnier)/W_k_furnier*100:.1f}%)")
print(f"    -> Stroemungstechnisch vernachlässigbar")
print(f"    -> Hauptunterschied: Oberfläche (Reibung) und Stabilität (Flutter)")
print(f"    Reibungsfaktor Folie:   f ~ 0.02 (glatt)")
print(f"    Reibungsfaktor Furnier: f ~ 0.03 (rau)")
print(f"    Kanalreibung ist aber nur {z_reib:.2e} des Gesamtverlusts -> irrelevant")

# Helmholtz
V_eff = V_brutto - L_wall*t_furnier*H_kam
r_eff = np.sqrt(S_klap_eff/np.pi); L_eff_h = 0.005+1.6*r_eff
f_H = (c_s/(2*np.pi))*np.sqrt(S_klap_eff/(V_eff*L_eff_h))
# Schlitz als zusätzlicher Hals: mittlere Tiefe des Schlitzes
# Effektive Halslänge für den Spalt = Plattendicke (mittlere)
L_eff_spalt = t_platte_mean + 0.8*np.sqrt(S_gap/np.pi)

print(f"\n  HELMHOLTZ-RESONANZ:")
print(f"    f_H = {f_H:.0f} Hz (Klappe als Hals)")
print(f"    Schlitz-Hals (Plattendicke): mittl. {t_platte_mean*1e3:.1f} mm Tiefe")
print(f"    -> Plattendicke wirkt als zusätzliche Halslänge für den Spalt")
print(f"    -> Erhöht den Strömungswiderstand des Spalts leicht")

print(f"\n  WIRKUNG DER PLATTENDICKE:")
print(f"    Am freien Ende: 13 mm -> Luft durchquert 13 mm langen Kanal im Schlitz")
print(f"    Am Einspann-Ende: 2 mm -> Kanal nur 2 mm lang")
print(f"    Zusätzlicher Kanalwiderstand: zeta_Schlitz ~ f*L/D_h")
D_h_schlitz = 2*W_schlitz*h_aufbiegung/(W_schlitz+h_aufbiegung)  # am Tip
z_schlitz = 0.03 * t_platte_frei / D_h_schlitz
print(f"    Am freien Ende: zeta ~ 0.03 * {t_platte_frei*1e3:.0f}mm / {D_h_schlitz*1e3:.2f}mm = {z_schlitz:.3f}")
print(f"    -> Gering, aber nicht vernachlässigbar (addiert sich zum Spaltverlust)")

# ═══════════════════════════════════════════════════════════════
# PDF
# ═══════════════════════════════════════════════════════════════
print(f"\n  Generiere PDF...")

C_WALL_D=HexColor('#8B4513'); C_ALU_D=HexColor('#B0B0B0'); C_STEEL_D=HexColor('#4169E1')
C_AIR_D=HexColor('#2196F3'); C_ARROW_D=HexColor('#E94560'); C_DIM_D=HexColor('#666')
C_CHAMBER_D=HexColor('#FFF8E1'); C_TRENN_D=HexColor('#2E7D32')

def make_fig_schnitt_keil():
    """Längsschnitt mit keilförmiger Stimmplatte und 1.5mm Aufbiegung"""
    w, h = 170*mm, 110*mm
    d = Drawing(w, h)
    sc = 1.3; ox, oy = 18*mm, 22*mm
    kl = 94*sc; kh = 40*sc

    # Kammer
    d.add(Rect(ox, oy, kl, kh, fillColor=C_CHAMBER_D, strokeColor=C_WALL_D, strokeWidth=1.5))

    # Stimmplatte KEILFÖRMIG oben
    # Links (Einspannung): 2mm, Rechts (frei/Klappe): 13mm
    t_l = 2*sc; t_r = 13*sc
    pts_plate = [ox, oy+kh, ox+kl, oy+kh, ox+kl, oy+kh+t_r, ox, oy+kh+t_l]
    d.add(Polygon(pts_plate, fillColor=C_ALU_D, strokeColor=black, strokeWidth=1))
    d.add(String(ox+20*sc, oy+kh+t_l+2*sc, 'Alu-Stimmplatte (keilfoermig: 2->13 mm)',
                 fontSize=7, fillColor=black))

    # Schlitz (gestrichelt) - wird tiefer nach rechts
    d.add(Line(ox+12*sc, oy+kh, ox+12*sc, oy+kh+t_l+0.5*sc,
               strokeColor=C_STEEL_D, strokeWidth=0.8, strokeDashArray=[2,2]))
    d.add(Line(ox+84*sc, oy+kh, ox+84*sc, oy+kh+t_r-1*sc,
               strokeColor=C_STEEL_D, strokeWidth=0.8, strokeDashArray=[2,2]))

    # Zunge auf Platte mit Aufbiegung
    # Einspannung links (flach auf Platte)
    zunge_y_l = oy + kh + t_l  # bündig mit Plattenoberfläche links
    # Freies Ende rechts (1.5mm aufgebogen)
    aufb_sc = 1.5*sc * 3  # übertrieben für Sichtbarkeit
    zunge_y_r = oy + kh + t_r + aufb_sc

    # Zunge als gekrümmte Linie (Aufbiegung)
    pts_z = []
    for i in range(41):
        t = i/40
        x = ox + 12*sc + t*72*sc
        # Plattenoberfläche bei x
        plate_top = oy + kh + t_l + t*(t_r - t_l)
        # Aufbiegung: quadratisch von 0 bis aufb_sc
        bend = aufb_sc * t**2
        pts_z.extend([x, plate_top + bend])
    d.add(PolyLine(pts_z, strokeColor=C_STEEL_D, strokeWidth=2))

    # Einspannung markieren
    d.add(Rect(ox+12*sc, zunge_y_l-1*sc, 8*sc, 2*sc, fillColor=None,
               strokeColor=HexColor('#1a237e'), strokeWidth=1))
    d.add(String(ox+12*sc, zunge_y_l+2*sc, 'Einsp. (Schrauben)', fontSize=5.5, fillColor=black))

    # Aufbiegung bemaßen
    tip_x = ox+84*sc
    tip_plate = oy+kh+t_r
    tip_zunge = tip_plate + aufb_sc
    d.add(Line(tip_x+3*sc, tip_plate, tip_x+3*sc, tip_zunge,
               strokeColor=C_ARROW_D, strokeWidth=0.8))
    d.add(Line(tip_x+1*sc, tip_plate, tip_x+5*sc, tip_plate,
               strokeColor=C_ARROW_D, strokeWidth=0.5))
    d.add(Line(tip_x+1*sc, tip_zunge, tip_x+5*sc, tip_zunge,
               strokeColor=C_ARROW_D, strokeWidth=0.5))
    d.add(String(tip_x+6*sc, tip_plate+aufb_sc/2, '1,5 mm', fontSize=6.5,
                 fillColor=C_ARROW_D, fontName='Helvetica-Bold'))

    # Plattendicke bemaßen (rechts)
    d.add(Line(ox+kl+3*sc, oy+kh, ox+kl+3*sc, oy+kh+t_r,
               strokeColor=C_DIM_D, strokeWidth=0.5))
    d.add(String(ox+kl+5*sc, oy+kh+t_r/2-1*sc, '13 mm', fontSize=6, fillColor=C_DIM_D))
    # Plattendicke links
    d.add(Line(ox-3*sc, oy+kh, ox-3*sc, oy+kh+t_l,
               strokeColor=C_DIM_D, strokeWidth=0.5))
    d.add(String(ox-12*sc, oy+kh+t_l/2, '2mm', fontSize=6, fillColor=C_DIM_D))

    # Klappe am Boden rechts
    klap_x = ox+kl-25*sc
    d.add(Rect(klap_x, oy-4*sc, 20*sc, 3*sc, fillColor=HexColor('#FFE0B2'),
               strokeColor=C_ARROW_D, strokeWidth=1))
    kl_len=16*sc; kl_ex=klap_x+kl_len*np.cos(alpha_r)
    kl_ey=oy-4*sc-kl_len*np.sin(alpha_r)
    d.add(Line(klap_x, oy-4*sc, kl_ex, kl_ey, strokeColor=C_ARROW_D, strokeWidth=2))
    d.add(String(klap_x-3*sc, oy-10*sc, 'Klappe 30°', fontSize=6, fillColor=C_ARROW_D, fontName='Helvetica-Bold'))

    # Trennwand
    tw_y = oy+kh/2-0.5*sc; tw_x = ox+gap_fold/0.001*sc
    d.add(Rect(tw_x, tw_y, L_wall/0.001*sc, 1*sc, fillColor=C_TRENN_D, strokeColor=black, strokeWidth=0.8))
    d.add(String(tw_x+12*sc, tw_y+2.5*sc, 'Trennwand', fontSize=6, fillColor=C_TRENN_D))
    d.add(String(ox+40*sc, tw_y+10*sc, 'Kanal A (oben)', fontSize=7, fillColor=HexColor('#1565C0'), fontName='Helvetica-Bold'))
    d.add(String(ox+40*sc, tw_y-7*sc, 'Kanal B (unten)', fontSize=7, fillColor=HexColor('#C62828'), fontName='Helvetica-Bold'))

    # Luftpfeile
    by = tw_y-8*sc; ay = tw_y+10*sc
    d.add(Line(klap_x+10*sc, oy+3*sc, klap_x+10*sc, tw_y-3*sc, strokeColor=C_AIR_D, strokeWidth=2))
    d.add(Polygon([klap_x+10*sc, oy+3*sc, klap_x+8*sc, oy+6*sc, klap_x+12*sc, oy+6*sc], fillColor=C_AIR_D))
    d.add(Line(klap_x+5*sc, by, ox+8*sc, by, strokeColor=C_AIR_D, strokeWidth=2))
    d.add(Polygon([ox+8*sc, by, ox+11*sc, by+2*sc, ox+11*sc, by-2*sc], fillColor=C_AIR_D))
    fx = ox+5*sc
    d.add(PolyLine([fx+3*sc, by, fx, by, fx, ay, fx+3*sc, ay], strokeColor=C_AIR_D, strokeWidth=2))
    d.add(Line(fx+3*sc, ay, ox+80*sc, ay, strokeColor=C_AIR_D, strokeWidth=2))
    d.add(Polygon([ox+80*sc, ay, ox+77*sc, ay+2*sc, ox+77*sc, ay-2*sc], fillColor=C_AIR_D))
    d.add(Line(ox+75*sc, oy+kh-2*sc, ox+75*sc, oy+kh+t_r-2*sc, strokeColor=C_AIR_D, strokeWidth=2))
    d.add(Polygon([ox+75*sc, oy+kh+t_r-2*sc, ox+73*sc, oy+kh+t_r-5*sc, ox+77*sc, oy+kh+t_r-5*sc], fillColor=C_AIR_D))

    # Kammermaße
    d.add(Line(ox, oy-12*sc, ox+kl, oy-12*sc, strokeColor=C_DIM_D, strokeWidth=0.5))
    d.add(String(ox+38*sc, oy-15*sc, '94 mm', fontSize=7, fillColor=C_DIM_D))

    d.add(String(ox, oy+kh+t_r+aufb_sc+4*sc,
        'Abb. 1: Laengsschnitt. Keilfoermige Platte (2->13 mm), Aufbiegung 1,5 mm am freien Ende.',
        fontSize=7, fillColor=black, fontName='Helvetica-Bold'))
    return d

def make_fig_varianten():
    w, h = 170*mm, 65*mm
    d = Drawing(w, h)
    labels = ['A: Gerade', 'B: Schraeg', 'C: Parabel']
    offsets = [5*mm, 60*mm, 115*mm]
    bw=45*mm; bh=42*mm
    for idx, (label, ox) in enumerate(zip(labels, offsets)):
        oy = 12*mm
        d.add(Rect(ox, oy, bw, bh, fillColor=C_CHAMBER_D, strokeColor=C_WALL_D, strokeWidth=1))
        d.add(Rect(ox+bw-12*mm, oy-2.5*mm, 10*mm, 2.5*mm, fillColor=HexColor('#FFE0B2'), strokeColor=C_ARROW_D, strokeWidth=0.6))
        if idx == 0:
            d.add(Line(ox+8*mm, oy+bh/2, ox+bw-3*mm, oy+bh/2, strokeColor=C_TRENN_D, strokeWidth=2))
        elif idx == 1:
            d.add(Line(ox+8*mm, oy+bh/2, ox+bw-3*mm, oy+bh*0.35, strokeColor=C_TRENN_D, strokeWidth=2))
        else:
            y_l=oy+bh/2; y_r=oy+bh*0.35; pts=[]
            for i in range(31):
                t=i/30; x=ox+8*mm+t*(bw-11*mm)
                y=y_l+t*(y_r-y_l)+4*mm*4*t*(1-t)
                pts.extend([x,y])
            d.add(PolyLine(pts, strokeColor=C_TRENN_D, strokeWidth=2))
        d.add(Line(ox+bw-7*mm, oy+2*mm, ox+bw-7*mm, oy+bh*0.3, strokeColor=C_AIR_D, strokeWidth=1.2))
        d.add(Line(ox+bw-8*mm, oy+bh*0.3, ox+12*mm, oy+bh*0.3, strokeColor=C_AIR_D, strokeWidth=1.2))
        d.add(Polygon([ox+12*mm, oy+bh*0.3, ox+14*mm, oy+bh*0.3+1.5*mm, ox+14*mm, oy+bh*0.3-1.5*mm], fillColor=C_AIR_D))
        d.add(PolyLine([ox+6*mm, oy+bh*0.3, ox+4*mm, oy+bh*0.3, ox+4*mm, oy+bh*0.7, ox+6*mm, oy+bh*0.7], strokeColor=C_AIR_D, strokeWidth=1.2))
        d.add(Line(ox+6*mm, oy+bh*0.7, ox+bw-8*mm, oy+bh*0.7, strokeColor=C_AIR_D, strokeWidth=1.2))
        d.add(Polygon([ox+bw-8*mm, oy+bh*0.7, ox+bw-10*mm, oy+bh*0.7+1.5*mm, ox+bw-10*mm, oy+bh*0.7-1.5*mm], fillColor=C_AIR_D))
        d.add(String(ox+2*mm, oy+bh+2*mm, label, fontSize=8, fillColor=black, fontName='Helvetica-Bold'))
    d.add(String(5*mm, 4*mm, 'Abb. 2: Trennwandvarianten. Bei B und C: eng bei Klappe (rechts), weit bei Faltung (links).', fontSize=6.5, fillColor=black, fontName='Helvetica-Bold'))
    return d

def make_fig_bernoulli():
    w, h = 170*mm, 50*mm
    d = Drawing(w, h)
    phases = [('1: Spalt enger', -0.15, 'v steigt, p sinkt\n-> Sog auf Zunge'), ('2: Fast zu', -0.35, 'Stroemung stoppt\n-> Feder zurueck'), ('3: Ueber Ruhelage', 0.2, 'Spalt weit offen'), ('4: Ruecklauf', 0.05, 'Zyklus neu')]
    for i, (title, dy, desc) in enumerate(phases):
        ox = 3*mm + i*43*mm; oy = 16*mm
        # Platte im Querschnitt (quer zur Zunge = Schmalseite -> gleichmässig dick)
        d.add(Rect(ox, oy+10*mm, 36*mm, 2.5*mm, fillColor=C_ALU_D, strokeColor=black, strokeWidth=0.7))
        d.add(Rect(ox+7*mm, oy+10*mm, 20*mm, 2.5*mm, fillColor=white, strokeColor=black, strokeWidth=0.4))
        zy = oy+11.2*mm + dy*8*mm
        d.add(Rect(ox+8*mm, zy, 18*mm, 0.7*mm, fillColor=C_STEEL_D, strokeColor=HexColor('#1a237e'), strokeWidth=0.5))
        if abs(dy) > 0.1:
            d.add(Line(ox+17*mm, oy+7*mm, ox+17*mm, oy+9.5*mm, strokeColor=C_AIR_D, strokeWidth=1.2))
            d.add(Polygon([ox+17*mm, oy+9.5*mm, ox+16*mm, oy+8*mm, ox+18*mm, oy+8*mm], fillColor=C_AIR_D))
        d.add(String(ox+1*mm, oy+16*mm, title, fontSize=6.5, fillColor=black, fontName='Helvetica-Bold'))
        for j, line in enumerate(desc.split('\n')):
            d.add(String(ox+1*mm, oy+1*mm-j*3*mm, line, fontSize=5.5, fillColor=C_DIM_D))
    d.add(String(3*mm, 3*mm, 'Abb. 3: Bernoulli-Kreislauf (Querschnitt quer zur Zunge - Keil nicht sichtbar)', fontSize=6.5, fillColor=black, fontName='Helvetica-Bold'))
    return d

# ── PDF bauen ──
path = "/home/claude/bass_v6.pdf"
doc = SimpleDocTemplate(path, pagesize=A4, leftMargin=16*mm, rightMargin=16*mm, topMargin=18*mm, bottomMargin=14*mm)
sty = getSampleStyleSheet()
sty.add(ParagraphStyle('MT', parent=sty['Title'], fontSize=14, spaceAfter=3*mm, textColor=HexColor('#1a1a2e'), alignment=TA_CENTER, leading=18))
sty.add(ParagraphStyle('Sub', parent=sty['Normal'], fontSize=9, alignment=TA_CENTER, textColor=HexColor('#666'), spaceAfter=4*mm))
sty.add(ParagraphStyle('Ch', parent=sty['Heading1'], fontSize=12, spaceBefore=6*mm, spaceAfter=3*mm, textColor=HexColor('#16213e'), keepWithNext=True))
sty.add(ParagraphStyle('Sec', parent=sty['Heading2'], fontSize=10, spaceBefore=4*mm, spaceAfter=2*mm, textColor=HexColor('#533483'), keepWithNext=True))
sty.add(ParagraphStyle('B', parent=sty['Normal'], fontSize=8.5, leading=12.5, alignment=TA_JUSTIFY, spaceAfter=2*mm))
sty.add(ParagraphStyle('Key', parent=sty['Normal'], fontSize=8.5, leading=12, spaceAfter=2.5*mm, backColor=HexColor('#fff3e0'), borderWidth=0.5, borderColor=HexColor('#e65100'), borderPadding=4))
sty.add(ParagraphStyle('Good', parent=sty['Normal'], fontSize=8.5, leading=12, spaceAfter=2.5*mm, backColor=HexColor('#e8f5e9'), borderWidth=0.5, borderColor=HexColor('#2e7d32'), borderPadding=4))
sty.add(ParagraphStyle('Warn', parent=sty['Normal'], fontSize=8.5, leading=12, spaceAfter=2.5*mm, backColor=HexColor('#fce4ec'), borderWidth=0.5, borderColor=HexColor('#c62828'), borderPadding=4))

st = []
st.append(Paragraph("Stroemungsanalyse einer 50-Hz-Bass-Stimmzunge", sty['MT']))
st.append(Paragraph("Keilfoermige Stimmplatte, 1,5 mm Aufbiegung, Folie vs. Furnier, Ausblas vs. Umlenk", sty['Sub']))
st.append(HRFlowable(width="100%", thickness=2, color=HexColor('#e94560')))
st.append(Spacer(1, 3*mm))

# Kap 1
st.append(Paragraph("Kapitel 1: Aufbau", sty['Ch']))
st.append(Paragraph("Im Bassteil eines Akkordeons sitzen Dutzende Stimmzungen in eigenen Kammern, alle am selben Balg. Der Balg erzeugt den Druck, die Tasten steuern ueber mechanische Hebel die Klappen. Die Klappe wird <b>nicht</b> vom Balgdruck geoeffnet - der Druck steht an allen geschlossenen Klappen an. Erst der Tastendruck oeffnet die Klappe mechanisch.", sty['Key']))

# Kap 2
st.append(Paragraph("Kapitel 2: Kammer, Stimmplatte und Aufbiegung", sty['Ch']))
st.append(Paragraph(f"Die Kammer ist ein Holzquader von {L_kam*1e3:.0f} x {B_kam*1e3:.0f} x {H_kam*1e3:.0f} mm. Die gesamte Oberseite wird von einer Aluminium-Stimmplatte verschlossen. Diese Platte ist <b>keilfoermig</b>: am freien Zungenende (ueber der Klappe) {t_platte_frei*1e3:.0f} mm dick, am Einspannende nur {t_platte_einsp*1e3:.0f} mm. Der Schlitz in der Platte ist daher am freien Ende {t_platte_frei*1e3:.0f} mm tief - die Luft durchquert einen kurzen Kanal, bevor sie die Zunge erreicht (Abbildung 1).", sty['B']))
st.append(Paragraph(f"Die Stahlzunge ({L_zunge*1e3:.0f} x {W_zunge*1e3:.0f} x {h_reed*1e3:.3f} mm) ist mit <b>zwei Schrauben</b> auf der Platte befestigt (kleinere Stimmplatten werden genietet). Das freie Ende ist um <b>{h_aufbiegung*1e3:.1f} mm aufgebogen</b> - es steht also im Ruhezustand {h_aufbiegung*1e3:.1f} mm ueber der Plattenoberflaeche. Dadurch entsteht ein dreieckiger Spalt: 0 mm an der Einspannung, {h_aufbiegung*1e3:.1f} mm am freien Ende.", sty['B']))
st.append(Paragraph(f"Die effektive Spaltflaeche (Integral ueber die Zungenlaenge) betraegt S<sub>Spalt</sub> = W x h<sub>max</sub>/2 = {W_zunge*1e3:.0f} x {h_aufbiegung*1e3:.1f}/2 = <b>{S_gap*1e6:.0f} mm<super>2</super></b>. Zum Vergleich: Die Klappenoeffnung hat {S_klap_eff*1e6:.0f} mm<super>2</super> - rund {S_klap_eff/S_gap:.0f}-mal groesser. Der Spalt bleibt der Engpass, aber weniger extrem als bei kleinerer Aufbiegung.", sty['B']))
st.append(make_fig_schnitt_keil())
st.append(Spacer(1, 2*mm))

# Kap 3
st.append(PageBreak())
st.append(Paragraph("Kapitel 3: Der gefaltete Luftweg", sty['Ch']))
st.append(Paragraph(f"Eine Trennwand teilt die Kammer in zwei Kanaele. Der Luftweg: Klappe (Boden rechts) -> Kanal B nach links -> 180°-Faltung -> Kanal A nach rechts -> Schlitz -> Spalt. Gesamtweg ueber 200 mm.", sty['B']))

# Kap 4: Folie vs Furnier
st.append(Paragraph("Kapitel 4: Trennwand - Folie oder Furnier?", sty['Ch']))
st.append(Paragraph(f"Fuer die Trennwand kommen zwei Materialien in Frage: <b>Kunststofffolie</b> (~{t_folie*1e3:.1f} mm) oder <b>Holzfurnier</b> (~{t_furnier*1e3:.1f} mm). Die Wanddicke beeinflusst die Kanalbreite: Bei gerader Wand ergibt Folie {W_k_folie*1e3:.2f} mm pro Kanal, Furnier {W_k_furnier*1e3:.2f} mm - eine Differenz von {(W_k_folie-W_k_furnier)*1e3:.2f} mm ({(W_k_folie-W_k_furnier)/W_k_furnier*100:.1f}%). Stroemungstechnisch ist das vernachlaessigbar.", sty['B']))
st.append(Paragraph("Der eigentliche Unterschied liegt in drei Eigenschaften:", sty['B']))

d_mat = [['', 'Folie (0,1 mm)', 'Furnier (0,5 mm)'],
    ['Oberflaeche', 'Glatt -> weniger Reibung\n(f ~ 0,02)', 'Rau -> mehr Reibung\n(f ~ 0,03)'],
    ['Stabilitaet', 'Flexibel -> kann flattern,\nbraucht Stuetzpunkte', 'Formstabil -> bleibt\nin Position'],
    ['Akustik', 'Kann als Membran schwingen\n-> Nebengeraeusche', 'Akustisch inert'],
    ['Verarbeitung', 'Kleben, spannen', 'Kleben, klemmen, einfacher']]
t_mat = Table(d_mat, colWidths=[25*mm, 50*mm, 50*mm])
t_mat.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),('BACKGROUND',(0,1),(0,-1),HexColor('#e8e8e8')),('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,0),(0,-1),'Helvetica-Bold'),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_mat)
st.append(Paragraph(f"Da die Kanalreibung nur {z_reib:.2e} des Gesamtverlusts ausmacht, ist der Reibungsunterschied zwischen Folie und Furnier fuer den Durchfluss bedeutungslos. Der Unterschied liegt ausschliesslich in Stabilitaet und Akustik. Furnier ist die sicherere Wahl; Folie nur sinnvoll, wenn der Platz extrem knapp ist.", sty['B']))

# Kap 5: Klappenorientierung
st.append(Paragraph("Kapitel 5: Klappenorientierung - Ausblas oder Umlenk?", sty['Ch']))
st.append(Paragraph(f"<b>Ausblas-Richtung:</b> Die Klappe oeffnet so, dass der 30°-Strahl direkt in Kanal B zeigt. Der Strahl hat eine horizontale Komponente von {np.cos(alpha_r)*100:.0f}% und eine vertikale von {np.sin(alpha_r)*100:.0f}%. Verlustbeiwert: zeta<sub>Klappe</sub> = {zeta_klap_ausblas:.3f} (Vena contracta + Richtungsaenderung).", sty['B']))
st.append(Paragraph(f"<b>Umlenk-Richtung:</b> Die Klappe oeffnet nach aussen. Die Luft muss ~150° um die Klappenkante kehrtmachen, bevor sie in die Kammer eintritt. Verlustbeiwert: zeta<sub>Klappe</sub> = {zeta_klap_umlenk:.3f} - das <b>{zeta_klap_umlenk/zeta_klap_ausblas:.1f}-fache</b> der Ausblas-Richtung.", sty['B']))

d_orient = [['', 'Ausblas', 'Umlenk'],
    ['Strahl', 'Gerichtet in Kammer', 'Erst weg, dann 150° Kurve'],
    ['zeta Klappe', f'{zeta_klap_ausblas:.3f}', f'{zeta_klap_umlenk:.3f}'],
    ['Wandform-Nutzung', 'Voll (Duese, Coanda)', 'Reduziert (breiter Eintritt)'],
    ['Nebengeraeusche', 'Weniger', 'Mehr (Abloesung an Kante)'],
    ['Anwendung', 'Standardwahl', 'Nur wenn Bauform es erfordert']]
t_orient = Table(d_orient, colWidths=[32*mm, 45*mm, 45*mm])
t_orient.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),('BACKGROUND',(0,1),(0,-1),HexColor('#e8e8e8')),('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,0),(0,-1),'Helvetica-Bold'),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_orient)
st.append(Paragraph("Der hoehere Verlust der Umlenk-Richtung wirkt sich auf die Spaltgeschwindigkeit kaum aus (Spalt dominiert), aber auf die Stroemungsqualitaet: Der Strahl tritt breiter und langsamer in die Kammer ein, die Turbulenz am Eintritt ist hoeher, und die Wandform (Duese bei B, Coanda bei C) kann nicht voll wirken.", sty['B']))

# --- Analyse: Klappenöffnungswinkel ---
st.append(Paragraph("Klappenoeffnung und Lautstaerke", sty['Sec']))

angles_deg = [30, 25, 20, 15, 10, 7, 5, 3, 2, 1]
valve_rows = [['Winkel', 'S<sub>Klap,eff</sub> [mm2]', 'S<sub>Klap</sub>/S<sub>Spalt</sub>',
               'zeta<sub>ges</sub>', 'Q [ml/s]', 'Q/Q<sub>30°</sub>']]
dP_valve = 1000  # Pa

# Reference at 30°
zeta_ref_cfg = configs['A-Gerade + Ausblas']['zeta']
v_ref = np.sqrt(2*dP_valve/(rho*zeta_ref_cfg))
Q_ref = v_ref * S_gap

for ang in angles_deg:
    ar = np.radians(ang)
    h_o = W_klap * np.sin(ar)
    S_kl = L_klap * h_o * C_c
    zeta_k_ang = (1/C_c - 1)**2 + np.sin(ar)**2  # Ausblas zeta at this angle
    
    # Recalculate total zeta with this valve area
    z_tot = zeta_k_ang * (S_gap/S_kl)**2  # valve contribution
    # Add internal losses (use A-Gerade as reference)
    for desc, z, S in [('Eintritt Kanal B', 0.5, A_k), ('Aufprall 90° auf Wand', 0.8, A_k),
                        ('180°-Faltung', 1.0, S_fold)]:
        z_tot += z * (S_gap/S)**2
    z_tot += 1.0  # Spaltaustritt
    # Kanalreibung
    z_tot += f_darcy * L_path / D_h * (S_gap/A_k)**2
    
    v_ang = np.sqrt(2*dP_valve/(rho*z_tot))
    Q_ang = v_ang * S_gap
    ratio = Q_ang / Q_ref
    
    valve_rows.append([
        f'{ang}°',
        f'{S_kl*1e6:.0f}',
        f'{S_kl/S_gap:.1f}x',
        f'{z_tot:.4f}',
        f'{Q_ang*1e6:.0f}',
        f'{ratio*100:.1f}%'
    ])

t_valve = Table(valve_rows, colWidths=[16*mm, 28*mm, 28*mm, 22*mm, 20*mm, 20*mm])
t_valve.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))

st.append(Paragraph(
    "Praktische Tests zeigen: Bei Verringerung der Klappenoeffnung treten "
    "Lautstaerke-Einbussen auf. Die stationaere Berechnung zeigt, warum:", sty['B']))

st.append(t_valve)
st.append(Spacer(1, 2*mm))

# Find threshold where Q drops significantly
ang_thresh_10 = 1  # default
for ang in angles_deg:
    ar = np.radians(ang)
    S_kl = L_klap * W_klap * np.sin(ar) * C_c
    zeta_k_ang = (1/C_c - 1)**2 + np.sin(ar)**2
    z_tot = zeta_k_ang * (S_gap/S_kl)**2 + 1.0
    for desc, z, S in [('E', 0.5, A_k), ('A', 0.8, A_k), ('F', 1.0, S_fold)]:
        z_tot += z * (S_gap/S)**2
    z_tot += f_darcy * L_path / D_h * (S_gap/A_k)**2
    v_a = np.sqrt(2*dP_valve/(rho*z_tot))
    Q_a = v_a * S_gap
    if Q_a/Q_ref < 0.95:
        ang_thresh_5pct = ang
        break
else:
    ang_thresh_5pct = 1

st.append(Paragraph(
    f"<b>Ergebnis:</b> Die stationaere Stroemung ist erstaunlich unempfindlich gegen "
    f"die Klappenoeffnung. Selbst bei 5° sind noch ueber 99% des maximalen Durchflusses "
    f"verfuegbar. Erst unter ~{ang_thresh_5pct}° sinkt Q messbar (ueber 5%).", sty['B']))

st.append(Paragraph(
    "Das die praktisch beobachteten Lautstaerke-Einbussen <b>nicht durch stationaeren "
    "Durchflussverlust</b> erklaerbar sind (der Spalt dominiert zu stark), deutet auf "
    "einen <b>instationaeren Mechanismus</b> hin:", sty['Warn']))

st.append(Paragraph(
    "<b>1. Langsamerer Druckaufbau:</b> Eine kleinere Klappenoeffnung verzoegert den "
    "Druckanstieg in der Kammer. Die Zunge braucht laenger, um die volle Amplitude "
    "zu erreichen. Bei schnellen Passagen wird die Maximalamplitude nie erreicht - "
    "das klingt leiser.", sty['B']))

st.append(Paragraph(
    "<b>2. Veraenderte Jet-Kopplung:</b> Bei kleinerem Winkel ist der Strahl flacher "
    "und schneller (gleicher Volumenstrom durch kleinere Oeffnung). Die Jet-Richtung "
    "aendert sich: Bei 30° hat der Strahl 50% vertikale Komponente, bei 10° nur noch "
    "17%. Weniger Vertikalkomponente bedeutet, dass der Strahl weniger direkt in den "
    "Kanal B gerichtet wird - aehnlich wie beim Uebergang von Ausblas zu Umlenk.", sty['B']))

st.append(Paragraph(
    "<b>3. Akustische Impedanzaenderung:</b> Die Klappenoeffnung ist der akustische "
    "Eingang der Kammer. Eine kleinere Oeffnung erhoet die akustische Impedanz am "
    "Kammereingang, was die Resonanzeigenschaften der Kammer veraendert. Bei bestimmten "
    "Oeffnungswinkeln koennen destruktive Interferenzen entstehen, die die "
    "Zungenamplitude reduzieren.", sty['B']))

st.append(Paragraph(
    "Alle drei Effekte wirken zusammen und erklaeren, warum die Lautstaerke-Einbussen "
    "in der Praxis deutlicher ausfallen als es die stationaere Berechnung erwarten "
    "laesst. Auch hier gilt: <b>Die Berechnung zeigt, dass der Effekt nicht im "
    "Durchfluss liegt - also muss er im instationaeren Verhalten liegen.</b>", sty['Key']))

# --- Kap 5b: Viertelkreis an der Faltung ---
st.append(Paragraph("Viertelkreis an der 180°-Faltung", sty['Sec']))

# Calculate bend losses with quarter-circle radii
D_h_fold = 2 * gap_fold * H_kam / (gap_fold + H_kam)  # hydraulic diameter of fold passage

# Idelchik: 90° bend with inner radius R
# Sharp (R=0): zeta_90 ~ 1.1
# R/D_h = 0.25: zeta_90 ~ 0.8
# R/D_h = 0.5: zeta_90 ~ 0.5  
# R/D_h = 1.0: zeta_90 ~ 0.2-0.3
# R/D_h = 2.0: zeta_90 ~ 0.1

# 180° turn = 2 x 90° plus additional mixing loss
# Sharp 180°: zeta_180 ~ 2.0-2.5 (used: 1.0 for A because S_fold >> S_gap makes it nearly zero)
# With radius R: zeta_180 ~ 2 * zeta_90(R/Dh) + 0.1 (residual)

radii_mm = [0, 3, 5, 8, 10, 15]
fold_rows = [['Radius R', 'R/D<sub>h</sub>', 'zeta<sub>Faltung</sub>',
              'zeta auf S<sub>Spalt</sub> bezogen', 'Anteil an zeta<sub>ges</sub>']]

zeta_total_ref = configs['A-Gerade + Ausblas']['zeta']

for R_mm in radii_mm:
    R = R_mm / 1000.0
    r_dh = R / D_h_fold
    
    # Idelchik approximation for 90° bend
    if r_dh < 0.01:
        z_90 = 1.1  # sharp
    else:
        z_90 = 0.21 * (1 + 0.7/r_dh**0.5)  # Idelchik interpolation
    z_180 = 2 * z_90 + 0.15  # two bends + mixing
    z_180 = min(z_180, 2.5)  # cap at sharp value
    
    # Scale to gap area
    z_scaled = z_180 * (S_gap / S_fold)**2
    anteil = z_scaled / zeta_total_ref * 100
    
    fold_rows.append([
        f'{R_mm} mm' if R_mm > 0 else 'scharf',
        f'{r_dh:.2f}' if R_mm > 0 else '0',
        f'{z_180:.2f}',
        f'{z_scaled:.6f}',
        f'{anteil:.4f}%'
    ])

t_fold_r = Table(fold_rows, colWidths=[18*mm, 18*mm, 28*mm, 36*mm, 30*mm])
t_fold_r.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))

st.append(Paragraph(
    f"Die 180°-Faltung am Ende der Trennwand hat scharfe Ecken. Kann man durch "
    f"Einsetzen eines Viertelkreis-Profils (Radius R) an den beiden 90°-Ecken den "
    f"Verlust reduzieren? Der hydraulische Durchmesser der Faltung ist "
    f"D<sub>h</sub> = {D_h_fold*1e3:.1f} mm (Spalt {gap_fold*1e3:.0f} mm x Hoehe "
    f"{H_kam*1e3:.0f} mm).", sty['B']))

st.append(t_fold_r)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    f"<b>Stationaeres Ergebnis:</b> Die Faltung traegt bei <b>allen Radien weniger als "
    f"0,01% zum Gesamtwiderstand</b> bei. Der Grund: Die Faltungsquerschnittsflaeche "
    f"S<sub>Faltung</sub> = {S_fold*1e6:.0f} mm2 ist {S_fold/S_gap:.0f}-mal groesser als "
    f"der Spalt ({S_gap*1e6:.0f} mm2). Der Verlustbeiwert wird mit (S<sub>Spalt</sub>/"
    f"S<sub>Faltung</sub>)2 = ({S_gap/S_fold:.4f})2 = {(S_gap/S_fold)**2:.2e} skaliert - "
    f"praktisch Null.", sty['B']))

st.append(Paragraph(
    "<b>Instationaeres Bild:</b> Fuer den Druckstoss-Durchgang ist die Frage "
    "differenzierter als beim stationaeren Durchfluss. Ein Viertelkreis-Profil "
    "mit R = 8-10 mm haette zwei Effekte:", sty['B']))

st.append(Paragraph(
    "<b>a) Reflexionen reduzieren:</b> An einer scharfen 90°-Ecke wird ein "
    "Druckpuls teilweise reflektiert (akustische Impedanzsprung). Je steiler "
    "die Geometrieaenderung, desto staerker die Reflexion. Ein Viertelkreis "
    "glaettet den Impedanzuebergang und laesst mehr Pulsenergie passieren.", sty['B']))

st.append(Paragraph(
    "<b>b) Wirbelbildung reduzieren:</b> Bei scharfen Ecken loest die Stroemung "
    "bei jedem Druckstoss sofort ab und bildet Wirbel, die Energie dissipieren. "
    "Mit Radius R > 5 mm bleibt die Stroemung laenger angelegt.", sty['B']))

st.append(Paragraph(
    "<b>Aber Vorsicht:</b> Eine <b>zu kohaerente</b> Stroemung am Spalt ist "
    "nicht wuenschenswert! Die Bernoulli-Selbsterregung braucht eine "
    "Anfangsstoerung, um die Schwingung auszuloesen. Eine perfekt gleichmaessige, "
    "laminare Stroemung wuerde die Zunge nur statisch auslenken, ohne sie in "
    "Schwingung zu versetzen. Ein gewisses Mass an Turbulenz und raeumlicher "
    "Ungleichmaessigkeit ist notwendig:", sty['Warn']))

st.append(Paragraph(
    "<b>Lokale Druckschwankungen</b> an der Zungenspitze stoessen die erste "
    "Auslenkung an. <b>Wirbelabloesungen</b> an den Spaltkanten liefern "
    "breitbandige Stoerungen, die die Zunge zum Schwingen anregen. "
    "Waere die Stroemung perfekt laminar und gleichfoermig, muesste die Zunge "
    "auf eine infinitesimale Instabilitaet warten - der Einschwingvorgang "
    "wuerde sich deutlich verzoegern.", sty['B']))

st.append(Paragraph(
    "Das Optimum liegt also nicht bei maximaler Kohaerenz, sondern bei "
    "<b>schnellem Druckaufbau mit ausreichend Stoerungen</b>. Konkret: "
    "Der Druckstoss soll schnell am Spalt ankommen (-> wenig Reflexion in "
    "der Faltung, gute Wandfuehrung), aber dort soll genuegend Turbulenz "
    "entstehen, um die Selbsterregung sofort auszuloesen. Die Turbulenz "
    "entsteht ohnehin am Spalt selbst (Re ~ 1400, Abloesungen an den "
    "Spaltkanten) - sie muss nicht zusaetzlich aus der Kammer angeliefert "
    "werden.", sty['B']))

st.append(Paragraph(
    "<b>Empfehlung:</b> Ein Viertelkreis-Profil R = 8-10 mm an den Ecken der "
    "Faltung reduziert Reflexionen und Energieverluste im Druckstoss-Transport, "
    "ohne die Turbulenz am Spalt zu beeintraechtigen (der Spalt erzeugt seine "
    "eigene Turbulenz). Es ist baulich einfach (z.B. eingeleimter Rundstab, "
    "Wachshohlkehle oder Epoxy-Verrundung). Ob der Effekt hoerbar ist, "
    "haengt davon ab, wie stark die Faltungsreflexionen den Anlauf tatsaechlich "
    "bremsen - das laesst sich nur durch Versuch klaeren.", sty['Key']))

# Kap 5b: Düsenwinkel-Experiment
st.append(PageBreak())
st.append(Paragraph("Kapitel 5b: Duesenwinkel-Experiment - warum der Anstroemwinkel die Amplitude bestimmt", sty['Ch']))

st.append(Paragraph(
    "Ein aufschlussreicher Praxistest: Stimmplatte montiert, Klappe entfernt, "
    "Anblasen mit einer Druckluftduese aus verschiedenen Winkeln. Ergebnis: "
    "Die Zungenamplitude aendert sich <b>von Null bis Maximum</b> allein durch "
    "Aenderung des Anstroemwinkels. Wie laesst sich das erklaeren?", sty['B']))

st.append(Paragraph("Die zwei Komponenten des Jets", sty['Sec']))

st.append(Paragraph(
    "Ein Luftstrahl, der unter dem Winkel alpha auf den Spalt trifft, hat zwei "
    "Komponenten:", sty['B']))

nozzle_rows = [['Winkel alpha', 'Normalkomp. sin(alpha)', 'Tangentialkomp. cos(alpha)', 'Charakter']]
for ang in [0, 10, 20, 30, 45, 60, 75, 90]:
    ar = np.radians(ang)
    char = ('Rein tangential - kein Eintritt' if ang == 0
            else 'Rein normal - statischer Druck' if ang == 90
            else 'Flach - viel Ueberstroemung' if ang <= 15
            else 'Optimal-Bereich?' if 20 <= ang <= 45
            else 'Steil - wenig Ueberstroemung')
    nozzle_rows.append([f'{ang}°', f'{np.sin(ar):.2f}', f'{np.cos(ar):.2f}', char])

t_nozzle = Table(nozzle_rows, colWidths=[22*mm, 32*mm, 34*mm, 45*mm])
t_nozzle.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_nozzle)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "<b>Normalkomponente</b> (senkrecht zur Platte, v*sin(alpha)): Treibt Luft durch "
    "den Spalt in den Schlitz hinein. Erzeugt statischen Ueberdruck unter der Zunge. "
    "Dieser Druck biegt die Zunge nach oben aus - aber er schwingt nicht von selbst. "
    "Ohne weitere Kopplung haelt er die Zunge nur statisch in einer Gleichgewichtslage.", sty['B']))

st.append(Paragraph(
    "<b>Tangentialkomponente</b> (parallel zur Platte, v*cos(alpha)): Stroemt ueber "
    "die Spalt<b>oeffnung</b> hinweg. Das ist der entscheidende Mechanismus:", sty['B']))

st.append(Paragraph("Der Floeten-Mechanismus am Zungenspalt", sty['Sec']))

st.append(Paragraph(
    "Wenn Luft ueber eine Oeffnung streicht, entsteht an der Oeffnung ein Unterdruck "
    "(Bernoulli). Genau das passiert bei der Tangentialkomponente des Jets: "
    "Sie stroemt ueber den Spalt und erzeugt einen Sog, der Luft aus dem Schlitz "
    "herauszieht - also die Zunge nach <b>unten</b> zieht. "
    "Gleichzeitig drueckt die Normalkomponente die Zunge nach <b>oben</b>.", sty['B']))

st.append(Paragraph(
    "Entscheidend ist: Beide Kraeft aendern sich, wenn die Zunge schwingt. "
    "Bewegt sich die Zunge nach unten (Spalt wird groesser), kann mehr Luft "
    "tangential ueber den Spalt stroemen -> staerkerer Sog -> Zunge wird weiter "
    "nach unten gezogen. Bewegt sie sich nach oben (Spalt wird enger), verengt "
    "sich der Kanal fuer die Tangentialstroemung -> weniger Sog -> Zunge federt "
    "zurueck. Das ist eine <b>positive Rueckkopplung</b>, die die Schwingung antreibt.", sty['B']))

st.append(Paragraph(
    "Dieser Mechanismus ist physikalisch identisch mit der Tonerzeugung bei "
    "Floeteninstrumenten: Dort stroemt ein flacher Luftstrahl ueber eine scharfe "
    "Kante (Labium), und die Wechselwirkung zwischen Jet und Resonator erzeugt "
    "den Ton. Die Zungenkante uebernimmt hier die Rolle des Labiums.", sty['B']))

st.append(Paragraph("Warum der Winkel die Amplitude bestimmt", sty['Sec']))

st.append(Paragraph(
    "<b>alpha = 0° (rein tangential):</b> Maximale Ueberstroemung, aber kein Druck "
    "im Schlitz. Die Zunge hat keine Vorauslenkung - der Bernoulli-Sog muss gegen "
    "die volle Federsteifigkeit arbeiten. Ausserdem fehlt der Durchfluss durch den "
    "Spalt, der fuer die klassische Bernoulli-Selbsterregung noetig ist. "
    "Ergebnis: Keine oder minimale Schwingung.", sty['B']))

st.append(Paragraph(
    "<b>alpha = 90° (rein normal):</b> Maximaler statischer Druck im Schlitz, "
    "aber keine Tangentialstroemung. Die Zunge wird nach oben gedrueckt und bleibt "
    "dort. Es fehlt der Bernoulli-Unterdruck an der Spaltoberseite, der die "
    "Zunge zurueckziehen wuerde. Es fehlt auch die Wechselwirkung zwischen "
    "Spaltweite und Saugkraft. Ergebnis: Statische Auslenkung, keine Schwingung.", sty['B']))

st.append(Paragraph(
    "<b>Optimaler Winkel (~ 20-40°):</b> Beide Komponenten sind in der richtigen "
    "Balance. Die Normalkomponente liefert den Grunddruck, der die Zunge "
    "in den Arbeitsbereich vorspannt. Die Tangentialkomponente liefert den "
    "geschwindigkeitsabhaengigen Bernoulli-Sog, der die Schwingung antreibt. "
    "Das Verhaeltnis entscheidet ueber die Phasenlage zwischen Druck und "
    "Zungenposition - und damit ueber die Energiezufuhr pro Zyklus.", sty['B']))

st.append(Paragraph(
    "Die Amplitude aendert sich stufenlos mit dem Winkel, weil die "
    "Energiezufuhr pro Schwingungszyklus vom Verhaeltnis der beiden Komponenten "
    "abhaengt. Am optimalen Winkel ist die Netto-Energiezufuhr maximal - die "
    "Phasenbeziehung zwischen Druckschwankung und Zungengeschwindigkeit ist "
    "so, dass in jedem Zyklus maximale Energie in die Schwingung gepumpt wird.", sty['B']))

st.append(Paragraph("Konsequenz fuer die Kammerauslegung", sty['Sec']))

st.append(Paragraph(
    "Dieses Experiment erklaert unmittelbar, warum die Klappenorientierung und die "
    "Trennwandform so wichtig sind: Sie bestimmen den <b>effektiven Anstroemwinkel</b> "
    "am Spalt.", sty['B']))

st.append(Paragraph(
    "Im normalen Betrieb (mit Klappe und Kammer) wird der Anstroemwinkel nicht "
    "direkt vom Balgdruck bestimmt, sondern von der Geometrie des Luftwegs. "
    "Die Luft kommt aus Kanal A (parallel zur Platte) und muss um 90° in den "
    "Schlitz einbiegen. Der effektive Anstroemwinkel haengt davon ab, wie die "
    "Stroemung in Kanal A verteilt ist:", sty['B']))

st.append(Paragraph(
    "<b>Schnelle, laminare Kanalstroemung</b> (gute Wandfuehrung, wenig Wirbel): "
    "Die Luft biegt gleichmaessig in den Spalt ein. Es entsteht eine definierte "
    "Tangentialkomponente, die den Bernoulli-Sog zuverlaessig erzeugt. "
    "<b>Langsame, verwirbelte Kanalstroemung</b> (schlechte Wandfuehrung, viele "
    "Reflexionen): Die Richtungsinformation ist verloren. Die Luft naehert sich "
    "dem Spalt aus allen Richtungen - die Tangentialkomponente mittelt sich "
    "teilweise heraus. Der Bernoulli-Antrieb wird schwaecher.", sty['B']))

st.append(Paragraph(
    "Das Duesenexperiment zeigt also <b>direkt</b>, warum die Trennwand das "
    "Einschwingverhalten beeinflusst: Nicht ueber den Durchfluss (der ist immer "
    "gleich), sondern ueber den Anstroemwinkel am Spalt. Eine saubere "
    "Stroemungsfuehrung in der Kammer erhaelt die Richtungsinformation des "
    "Druckstosses - und damit den optimalen Anstroemwinkel fuer die "
    "Bernoulli-Rueckkopplung. Das erklaert auch, warum Variante C (Parabel) "
    "in der Praxis die beste Ansprache liefert: Sie erhaelt die "
    "Stroemungsrichtung am besten.", sty['Key']))

# Kap 6: Trennwandvarianten
st.append(Paragraph("Kapitel 6: Drei Trennwandvarianten", sty['Ch']))
st.append(make_fig_varianten())
st.append(Spacer(1, 2*mm))
st.append(Paragraph(f"<b>A - Gerade:</b> Beide Kanaele {W_k*1e3:.1f} mm. 30°-Strahl prallt ~90° auf Wand. Stagnationspunkt, Abloesung, hoechste Turbulenz.", sty['B']))
st.append(Paragraph(f"<b>B - Schraeg (beta={beta_B:.1f}°):</b> Kanal B bei Klappe {W_B_unten_klap*1e3:.1f} mm (Duese), bei Faltung {W_B_unten_falt*1e3:.1f} mm. Kanal A bei Spalt {W_B_oben_klap*1e3:.1f} mm (weit). Diffusor-Halbwinkel nur {np.degrees(np.arctan((W_B_unten_falt-W_B_unten_klap)/(2*L_wall))):.1f}° -> sicher abloesungsfrei.", sty['B']))
st.append(Paragraph(f"<b>C - Parabel:</b> Gleiche Endpunkte wie B, aber stetige Kruemmung. Coanda-Effekt fuehrt Strahl entlang konvexer Wand. Sauberste Stroemung.", sty['B']))

# Ergebnistabelle
v1k = np.sqrt(2*1000/(rho*2))
Q1k = v1k * S_gap
d_res = [['(bei 1000 Pa, Furnier, Ausblas)', 'A: Gerade', 'B: Schraeg', 'C: Parabel']]
for cname_list in [['A-Gerade + Ausblas', 'B-Schraeg + Ausblas', 'C-Parabel + Ausblas']]:
    d_res.append(['zeta total'] + [f"{configs[c]['zeta']:.6f}" for c in cname_list])
    d_res.append(['v Spalt [m/s]'] + [f"{configs[c]['v1000']:.2f}" for c in cname_list])
    d_res.append(['Q [ml/s]'] + [f"{configs[c]['Q1000']*1e6:.1f}" for c in cname_list])
d_res.append(['Kanal B Klappe [mm]', f'{W_k*1e3:.1f}', f'{W_B_unten_klap*1e3:.1f}', f'{W_B_unten_klap*1e3:.1f}'])
d_res.append(['Kanal A Spalt [mm]', f'{W_k*1e3:.1f}', f'{W_B_oben_klap*1e3:.1f}', f'{W_B_oben_klap*1e3:.1f}'])
t_res = Table(d_res, colWidths=[40*mm, 28*mm, 28*mm, 28*mm])
t_res.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),('BACKGROUND',(0,1),(0,-1),HexColor('#e8e8e8')),('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,1),(0,-1),'Helvetica-Bold'),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_res)

st.append(Paragraph(
    "<b>Achtung: Diese Tabelle zeigt nur den stationaeren Zustand.</b> "
    "Die nahezu identischen Volumenströme taueschen darueber hinweg, dass "
    "die Trennwand das Einschwingverhalten in der Praxis entscheidend "
    "beeinflusst. Praktische Messungen zeigen deutliche Unterschiede in der "
    "Ansprache zwischen verschiedenen Wandformen - obwohl der stationaere "
    "Durchfluss bei allen fast gleich ist.", sty['Warn']))

st.append(Paragraph("Warum die Trennwand trotzdem entscheidend ist:", sty['Sec']))

st.append(Paragraph(
    "Der Widerspruch loest sich auf, wenn man den <b>instationaeren Anlauf</b> "
    "betrachtet statt den Gleichgewichtszustand. In den ersten Millisekunden nach "
    "Klappenoeffnung passiert folgendes: Ein Druckstoss laeuft vom Ventil durch "
    "Kanal B, um die Faltung, durch Kanal A zum Spalt. Die Trennwand ist der "
    "Wellenleiter fuer diesen Druckstoss. Ihre Form bestimmt, wie schnell "
    "der Druck am Spalt ankommt und wie viel Energie auf dem Weg durch "
    "Reflexionen und Wirbelbildung verlorengeht.", sty['B']))

st.append(Paragraph(
    "<b>Gerade Wand (A):</b> Der 30-Grad-Strahl prallt auf die Wand, erzeugt einen "
    "Stagnationspunkt und Wirbel. Der Druckstoss wird breit gestreut und "
    "teilweise reflektiert. Am Spalt kommt er verspaetet und mit reduzierter "
    "Amplitude an. Der Druckaufbau ist langsam - die Zunge erreicht den "
    "Schwellendruck fuer die Bernoulli-Rueckkopplung spaeter.", sty['B']))

st.append(Paragraph(
    "<b>Schraege Wand (B):</b> Die Duese am Eintritt beschleunigt den Druckstoss. "
    "Weniger Reflexionen auf dem Weg. Am Spalt kommt der Druck schneller auf "
    "den Schwellenwert. Der Anlauf ist kuerzer.", sty['B']))

st.append(Paragraph(
    "<b>Parabolische Wand (C):</b> Die stetige Kruemmung fuehrt den Druckstoss "
    "ohne Abloesung entlang der Wand (Coanda-Effekt). Minimale Reflexionen, "
    "schnellster Druckaufbau am Spalt. Die Turbulenz, die die Selbsterregung "
    "ausloest, entsteht ohnehin am Spalt selbst (Abloesungen an den Kanten bei "
    "Re ~ 1400) und muss nicht aus der Kammer angeliefert werden.", sty['B']))

st.append(Paragraph(
    "In der Elektrotechnik-Analogie (Kapitel 10): Die stationaere Berechnung "
    "entspricht dem Gleichstromwiderstand eines Kabels - fuer alle drei Varianten "
    "fast gleich. Aber die <b>Impulsantwort</b> (wie schnell ein Spannungssprung "
    "am Ende ankommt) haengt von der Wellenimpedanz, Reflexionen und Dispersion "
    "ab - und die sind voellig verschieden. Fuer die Ansprache zaehlt nicht der "
    "Gleichstromwiderstand, sondern die Impulsantwort.", sty['Key']))

# Kap 7: Bernoulli + Einschwingzeit
st.append(PageBreak())
st.append(Paragraph("Kapitel 7: Bernoulli-Selbsterregung und Einschwingzeit", sty['Ch']))
st.append(Paragraph("Die Selbsterregung erfolgt ueber den Bernoulli-Mechanismus (Abbildung 3): Verengt sich der Spalt, steigt v, sinkt p, Unterdruck verstaerkt Auslenkung. Am Umkehrpunkt bricht die Stroemung zusammen, Feder treibt Zunge zurueck.", sty['B']))
st.append(make_fig_bernoulli())
st.append(Spacer(1, 2*mm))
st.append(Paragraph(f"Schwellendruck: DeltaP<sub>min</sub> = {dp_min:.1f} Pa (mit {h_aufbiegung*1e3:.1f} mm Aufbiegung). Niedrig gegenueber 500-5000 Pa Balgdruck.", sty['B']))
st.append(Paragraph("Die Einschwingzeit tau = Q/(pi*f) ist der eigentliche Engpass:", sty['B']))
rows = [['zeta', 'Q', 'tau [ms]', 'Perioden']]
for z, Q, tau in zip(zeta_vals, Q_vals, tau_vals):
    rows.append([f'{z}', f'{Q:.0f}', f'{tau*1e3:.0f}', f'{Q:.0f} x 20 ms'])
t_tau = Table(rows, colWidths=[22*mm, 22*mm, 22*mm, 30*mm])
t_tau.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),('FONTSIZE',(0,0),(-1,-1),8),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(1,1),(-1,-1),'Courier'),('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_tau)
st.append(Paragraph("Q haengt von Material, Geometrie und Einspannung ab und muss gemessen werden.", sty['B']))

# ═══════════════════════════════════════════════════════════════
# Kap 8: Dämpfung - was sie beeinflusst und wie
# ═══════════════════════════════════════════════════════════════
st.append(PageBreak())
st.append(Paragraph("Kapitel 8: Daempfung - was sie beeinflusst und wie", sty['Ch']))

st.append(Paragraph(
    "Die Guete Q bestimmt die Einschwingzeit und damit die Ansprache. Q ist keine "
    "Materialkonstante, sondern die Summe aller Energieverluste pro Schwingungszyklus. "
    "Diese lassen sich in vier physikalisch unterscheidbare Mechanismen zerlegen, "
    "deren Beitraege sich addieren:", sty['B']))

st.append(Paragraph(
    "<b>1. Materialdaempfung (innere Reibung):</b> "
    "Jedes Metall wandelt bei Verformung einen Bruchteil der elastischen Energie in "
    "Waerme um. Der Verlustfaktor eta ist eine Materialeigenschaft. Fuer Federstahl "
    "liegt eta bei 0,0002-0,001, fuer Messing bei 0,001-0,003 - also bis zu zehnmal "
    "hoeher. In der Guete-Zerlegung: Q<sub>Material</sub> = 1/eta. Federstahl kommt "
    "allein auf Q<sub>Material</sub> ~ 1000-5000, Messing auf 300-1000. Dieser "
    "Beitrag ist frequenzunabhaengig.", sty['B']))

st.append(Paragraph(
    "<b>2. Luftdaempfung (viskose Verluste im Spalt):</b> "
    "Die schwingende Zunge verdraengt Luft durch den engen Spalt. Dabei entsteht ein "
    "Squeeze-Film-Effekt: Die Luft wird bei jedem Zyklus durch die Engstelle gepresst "
    "und wieder zurueckgesaugt. Die dissipierte Leistung skaliert mit "
    "P<sub>visc</sub> ~ mu*W*L<super>3</super>*omega<super>2</super>*a<super>2</super>/h<super>3</super>, "
    "wobei h die Spalthoehe ist. Mit " + f"{h_aufbiegung*1e3:.1f}" + " mm Aufbiegung ist h deutlich "
    "groesser als bei einem 0,4-mm-Spalt - der Squeeze-Film-Beitrag sinkt mit "
    "h<super>-3</super>, also um Faktor (1,5/0,4)<super>3</super> ~ 53. "
    "Bei der grossen Aufbiegung einer Basszunge ist die Luftdaempfung daher "
    "ein untergeordneter Beitrag.", sty['B']))

st.append(Paragraph(
    "<b>3. Einspannungsdaempfung:</b> "
    "An der Verschraubung (zwei Schrauben bei Basszungen) wird bei jeder Schwingung "
    "mikroskopisch Energie in Reibung umgewandelt. Die Guete der Einspannung haengt "
    "vom Anpressdruck, der Oberflaechenguete und dem Plattenmaterial ab. Aluminium "
    "ist weicher als Stahl - bei jeder Schwingung verformt sich die Kontaktzone "
    "minimal. Dieser Beitrag ist schwer zu berechnen, aber erfahrungsgemaess einer "
    "der groessten Einzelposten. Eine sauber gefraeste Auflageflaeche und "
    "gleichmaessig angezogene Schrauben minimieren ihn.", sty['B']))

st.append(Paragraph(
    "<b>4. Schallabstrahlung:</b> "
    "Die schwingende Zunge strahlt Schall ab - das ist ja der Zweck. Die abgestrahlte "
    "Leistung geht dem Schwinger als Daempfung verloren. Bei 50 Hz und einer Flaeche "
    "von 70 x 8 mm ist der Strahlungswiderstand gering (die Zunge ist viel kleiner als "
    "die Wellenlaenge von 6,9 m). Dieser Beitrag ist bei tiefen Basszungen klein, "
    "waechst aber mit f<super>2</super> und wird bei hohen Diskant-Frequenzen relevant.", sty['B']))

st.append(Paragraph(
    "Die Gesamt-Guete ergibt sich aus den Einzelbeitraegen als harmonische Summe: "
    "1/Q<sub>ges</sub> = 1/Q<sub>Mat</sub> + 1/Q<sub>Luft</sub> + 1/Q<sub>Einsp</sub> + "
    "1/Q<sub>Strahl</sub>. Der kleinste Einzelwert dominiert. Bei einer typischen "
    "Basszunge aus Federstahl mit Schraubbefestigung und 1,5 mm Aufbiegung duerfte "
    "Q<sub>ges</sub> zwischen 50 und 150 liegen, wobei die Einspannung und "
    "Materialdaempfung die groessten Beitraege liefern.", sty['Key']))

st.append(Paragraph(
    "Was folgt daraus fuer die Praxis? Q laesst sich <b>senken</b> (schnellere Ansprache) "
    "durch: Materialwahl (Messing statt Stahl), Massebeladung am Tip mit weichem Metall "
    "(erhoehte Materialdaempfung des Gewichts), geringere Aufbiegung (mehr Squeeze-Film), "
    "oder weichere Einspannung. Q laesst sich <b>erhoehen</b> (laengeres Nachklingen, "
    "lauterer Ton) durch: haerteren Stahl, praezisere Einspannung, groessere Aufbiegung. "
    "Der Instrumentenbauer navigiert zwischen diesen Polen - Messung (z.B. Ausschwing-"
    "versuch, Logarithmisches Dekrement) ist der einzige zuverlaessige Weg, Q einer "
    "konkreten Zunge zu bestimmen.", sty['B']))

# ═══════════════════════════════════════════════════════════════
# Kap 9: Propeller, Windräder und Hinterkanten-Optimierung
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 9: Lehren aus der Propeller- und Windkraft-Aerodynamik", sty['Ch']))

st.append(Paragraph(
    "In der Windkraft- und Propellertechnik sind in den letzten Jahrzehnten erhebliche "
    "Fortschritte bei der Optimierung von Stroemungskanten erzielt worden. Viele dieser "
    "Erkenntnisse betreffen dieselben physikalischen Phaenomene, die auch in der "
    "Stimmzungenkammer auftreten.", sty['B']))

st.append(Paragraph("<b>Hinterkanten-Optimierung (Trailing Edge):</b> "
    "Windkraftrotoren verwenden gezackte oder gekaemmte Hinterkanten (Serrations), um "
    "tonale Geraeusche zu reduzieren. Der Mechanismus: An einer scharfen Kante loesen "
    "sich Wirbel kohaerenter in einer einzigen Frequenz ab (Karman-Wirbel). Gezackte "
    "Kanten brechen diese Kohaerenz auf - die Wirbel loesen sich an verschiedenen "
    "Positionen zu verschiedenen Zeiten, das Geraeusch wird breitbandig statt tonal. "
    "Dasselbe Prinzip gilt fuer die Trennwand-Hinterkante am Faltspalt. Eine "
    "ausgefranste oder gezackte Wandkante wuerde den 180°-Umlenkungsverlust nicht "
    "senken, aber die Geraeuschqualitaet aendern.", sty['B']))

st.append(Paragraph("<b>Vorderkanten-Optimierung (Leading Edge):</b> "
    "Buckelwal-inspirierte Vorderkanten (Tubercles) verhindern den gleichzeitigen "
    "Stroemungsabriss ueber die gesamte Spannweite. Auf die Kammer uebertragen: "
    "Der 30°-Strahl aus der Klappe trifft auf die Trennwand wie auf eine "
    "Fluegelvorderkante. Eine wellenfoermige Wandkante (anstatt der geraden Kante "
    "bei Variante A) koennte die grossflaechige Abloesung in viele kleine, "
    "oertlich begrenzte Abloesungen aufteilen. Das Ergebnis waere weniger "
    "niederfrequente Druckschwankung am Zungenspalt.", sty['B']))

st.append(Paragraph("<b>Diffusor-Geometrie (Windkanal-Forschung):</b> "
    "Bei Windkraft-Diffusoren (Shrouded Turbines) ist bekannt, dass der "
    "Diffusor-Halbwinkel unter 7° bleiben muss, um Abloesung zu vermeiden. "
    "Die schraege Trennwand (Variante B) mit 3,7° liegt sicher darunter. "
    "Aber die Windkraftforschung zeigt auch: Stetige Kruemmung (wie bei Variante C) "
    "ist einem geraden Diffusor bei gleichem Oeffnungswinkel ueberlegen, weil "
    "der Druckgradient gleichmaessiger verteilt wird. Die Analogie stuetzt die "
    "Erwartung, dass C besser als B abschneidet.", sty['B']))

st.append(Paragraph(
    "Die Uebertragbarkeit hat Grenzen: Propeller und Windraeder arbeiten bei "
    "Reynolds-Zahlen von 10<super>5</super>-10<super>7</super>, die Stimmzungenkammer bei "
    "Re ~ 50-1500. Bei niedrigem Re sind Viskositaetseffekte staerker, "
    "Grenzschichten dicker, und Uebergaenge zwischen laminar und turbulent "
    "verlaufen anders. Dennoch sind die geometrischen Prinzipien "
    "(stetige Kruemmung, Vermeidung scharfer Kanten, Diffusor-Winkel) "
    "Reynolds-unabhaengig gueltig.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 10: Kammer-Zunge Kopplung: Serien- statt Parallelkreis
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 10: Warum resonante Kammern die Ansprache verschlechtern", sty['Ch']))

st.append(Paragraph(
    "Ein verbreiteter Gedanke ist, dass eine resonante Kammer die Zunge unterstuetzen "
    "koennte - aehnlich wie ein Resonanzkoerper einer Gitarre. Praktische Messungen "
    "zeigen das Gegenteil: Kammern, deren Helmholtz-Frequenz ein ganzzahliges "
    "Vielfaches der Zungenfrequenz ist, bringen ein schlechteres Ansprechverhalten "
    "als nicht-resonante Kammern.", sty['B']))

st.append(Paragraph(
    "Die Erklaerung liegt in der Art der Kopplung. In der Elektrotechnik gibt es "
    "zwei grundlegend verschiedene Resonanzkreis-Topologien:", sty['B']))

st.append(Paragraph(
    "<b>Parallelkopplung</b> (Schwingkreis mit L und C parallel): Beide Elemente "
    "werden gleichzeitig erregt, die Energie pendelt zwischen ihnen. Der Strom "
    "durch die Quelle ist minimal bei Resonanz - die Impedanz ist maximal. "
    "Ein Gitarrenkoerper arbeitet so: Die Saitenenergie koppelt ueber den Steg "
    "gleichzeitig in den Koerper, Decke und Luft schwingen gemeinsam, die "
    "Abstrahlung wird effizienter.", sty['B']))

st.append(Paragraph(
    "<b>Serienkopplung</b> (L und C in Reihe): Die Energie muss <i>durch</i> "
    "beide Elemente hindurch. Bei Resonanz sinkt die Gesamtimpedanz - der Strom "
    "steigt. Aber jedes Element hat seine eigene Zeitkonstante, und die Gesamtantwort "
    "ist langsamer als die des schnelleren Elements allein.", sty['B']))

st.append(Paragraph(
    "Die Stimmzungenkammer ist eine <b>Serienkopplung</b>. Die Luft muss den Weg "
    "Balg -> Klappe -> Kanal B -> Faltung -> Kanal A -> Schlitz -> Spalt -> Zunge "
    "in genau dieser Reihenfolge durchlaufen. Die Kammer sitzt als akustische "
    "Impedanz <i>im Signalweg</i> zwischen Druckquelle (Balg) und Last (Zunge). "
    "Jede Kammer-Resonanz fuegt dem System eine eigene Zeitkonstante hinzu.", sty['Key']))

st.append(Paragraph(
    "Was passiert bei einer resonanten Kammer? Wenn die Helmholtz-Frequenz ein "
    "Vielfaches der Zungenfrequenz ist, kann die Kammer in einem der Zugngen-Obertöne "
    "mitschwingen. Dabei wird Energie aus der Hauptschwingung der Zunge in die "
    "Kammer-Eigenschwingung umgeleitet. Die Kammer speichert diese Energie und gibt "
    "sie zeitversetzt zurueck. Fuer den Amplitudenaufbau der Zunge wirkt das wie "
    "eine zusaetzliche Daempfung: Ein Teil der Bernoulli-Energie, die die Zunge "
    "aufbauen sollte, verschwindet in der Kammerresonanz.", sty['B']))

st.append(Paragraph(
    "In der Elektrotechnik-Analogie: Der Serienkreis hat bei Resonanz minimale "
    "Impedanz. Das bedeutet maximaler Stromfluss - aber auch maximale Verluste "
    "in den parasitaeren Widerstaenden des Kreises. Uebertragen: Bei Kammerresonanz "
    "ist der Volumenstrom durch den Spalt zwar maximal, aber die Phasenbeziehung "
    "zwischen Druck und Zungenposition verschiebt sich. Die Bernoulli-Rueckkopplung "
    "braucht eine bestimmte Phasenlage, um Energie in die Zunge zu pumpen. Wenn die "
    "Kammer die Phase verschiebt, wird die Rueckkopplung weniger effizient - "
    "oder im schlimmsten Fall sogar bremsend.", sty['B']))

st.append(Paragraph(
    "Die praktische Konsequenz: <b>Die Kammer soll akustisch moeglichst 'unsichtbar' "
    "sein.</b> Das heisst: Ihre Eigenfrequenzen sollen weit von der Zungenfrequenz "
    "und deren Obertoenen entfernt liegen. Mit f<sub>H</sub> = " + f"{f_H:.0f}" + " Hz bei einer "
    "50-Hz-Zunge ist das Verhaeltnis 9,2 - kein ganzzahliges Vielfaches. Das ist "
    "guenstig. Eine Kammer mit f<sub>H</sub> = 500 Hz (10x) waere problematisch, "
    "weil der 10. Oberton der Zunge (10 x 50 = 500 Hz) in die Kammerresonanz "
    "fallen wuerde.", sty['Warn']))

st.append(Paragraph(
    "Praktische Erfahrung zeigt, dass die Trennwand das Einschwingverhalten "
    "entscheidend beeinflusst - obwohl alle Varianten im stationaeren Zustand "
    "nahezu identische Volumenströme liefern. Die Erklaerung liegt in der "
    "Serienkopplung: Die Kammer ist kein Verstaerker, sondern ein Wellenleiter. "
    "Ihre Form bestimmt nicht den Durchfluss, sondern die Impulsantwort - "
    "wie schnell und mit wie viel Amplitude ein Druckstoss am Spalt ankommt. "
    "Reflexionen, Wirbelabloesungen und Phasenverschiebungen in der Kammer "
    "bremsen den Anlauf der Bernoulli-Rueckkopplung. Je sauberer die "
    "Wandgeometrie den Druckstoss fuehrt, desto schneller spricht die Zunge an.", sty['B']))

st.append(Paragraph(
    "Die Optimierung hat damit zwei Ebenen: Die Zunge selbst (Q, Aufbiegung, "
    "Material) bestimmt die Grundansprache. Die Kammergeometrie (Trennwandform, "
    "Klappenrichtung) bestimmt, wie viel davon im Instrument tatsaechlich ankommt. "
    "Beides muss stimmen.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 11: Resonanzanalyse - welche Töne sind betroffen?
# ═══════════════════════════════════════════════════════════════
# Erzeuge Diagramm mit matplotlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

note_names_de = ['C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'B', 'H']
def note_freq(ni, oc):
    return 440.0 * 2**((oc-4)*12 + (ni-9))/12.0

# Korrektur: die Formel braucht Klammern
def note_freq(ni, oc):
    semitones = (oc - 4) * 12 + (ni - 9)
    return 440.0 * 2**(semitones / 12.0)

Q_H_lo, Q_H_hi = 5, 10
bw_wide = f_H / Q_H_lo; bw_narrow = f_H / Q_H_hi

bass_notes = []
for oc in range(1, 5):
    for ni in range(12):
        f = note_freq(ni, oc)
        if 30 <= f <= 300:
            bass_notes.append((f"{note_names_de[ni]}{oc}", f, ni, oc))

res = []
for name, f, ni, oc in bass_notes:
    n_c = max(1, round(f_H / f))
    best_n, best_d = n_c, abs(n_c*f - f_H)
    for nt in [n_c-1, n_c+1]:
        if nt >= 1 and abs(nt*f - f_H) < best_d:
            best_n, best_d = nt, abs(nt*f - f_H)
    abst = best_n*f - f_H
    sev = 2 if abs(abst) < bw_narrow/2 else (1 if abs(abst) < bw_wide/2 else 0)
    res.append({'name': name, 'f': f, 'n': best_n, 'f_ot': best_n*f,
                'abstand': abst, 'severity': sev})

# Diagramm
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(11, 8.5), gridspec_kw={'height_ratios': [2, 1.2]})
names_r = [r['name'] for r in res]
abstande_r = [r['abstand'] for r in res]
sevs = [r['severity'] for r in res]
obs = [r['n'] for r in res]
cols = ['#4CAF50' if s==0 else '#FF9800' if s==1 else '#F44336' for s in sevs]

bars = ax1.bar(range(len(res)), abstande_r, color=cols, edgecolor='#333', linewidth=0.5, width=0.8)
ax1.axhspan(-bw_narrow/2, bw_narrow/2, color='#F44336', alpha=0.12, label=f'Kritisch (Q_H=10, +/-{bw_narrow/2:.0f} Hz)')
ax1.axhspan(-bw_wide/2, -bw_narrow/2, color='#FF9800', alpha=0.10, label=f'Ungünstig (Q_H=5, +/-{bw_wide/2:.0f} Hz)')
ax1.axhspan(bw_narrow/2, bw_wide/2, color='#FF9800', alpha=0.10)
ax1.axhline(y=0, color='red', linewidth=1.5, label=f'f_H = {f_H:.0f} Hz')
for i, (bar, n) in enumerate(zip(bars, obs)):
    y = bar.get_height()
    ax1.text(i, y+(2 if y>=0 else -5), f'n={n}', ha='center', va='bottom' if y>=0 else 'top',
             fontsize=5.5, color='#333', rotation=90)
ax1.set_xticks(range(len(res)))
ax1.set_xticklabels(names_r, rotation=90, fontsize=7, fontweight='bold')
ax1.set_ylabel('Abstand naechster Oberton zu f_H [Hz]', fontsize=10)
ax1.set_title(f'Welche Basstoene werden durch die Kammerresonanz (f_H = {f_H:.0f} Hz) benachteiligt?\n'
              f'Temperierte Stimmung, Kammer {L_kam*1e3:.0f} x {B_kam*1e3:.0f} x {H_kam*1e3:.0f} mm', fontsize=11, fontweight='bold')
ax1.legend(loc='upper right', fontsize=8); ax1.set_ylim(-60, 60); ax1.grid(axis='y', alpha=0.3)
ax1t = ax1.twiny(); ax1t.set_xlim(ax1.get_xlim())
tp = list(range(0, len(res), 4))
ax1t.set_xticks(tp); ax1t.set_xticklabels([f'{res[i]["f"]:.1f} Hz' for i in tp], fontsize=7)
ax1t.set_xlabel('Grundfrequenz', fontsize=8)

# Oberton-Landkarte
ax2.axhspan(f_H-bw_wide/2, f_H+bw_wide/2, color='#FF9800', alpha=0.15, label='Resonanzband Q=5')
ax2.axhspan(f_H-bw_narrow/2, f_H+bw_narrow/2, color='#F44336', alpha=0.15, label='Resonanzband Q=10')
ax2.axhline(y=f_H, color='red', linewidth=1.5, linestyle='--', alpha=0.8)
for i, r in enumerate(res):
    for n in range(1, int(600/r['f'])+2):
        fot = n*r['f']
        if fot > 600: break
        d2 = abs(fot - f_H)
        if d2 < bw_narrow/2: c,s,m = '#F44336',25,'D'
        elif d2 < bw_wide/2: c,s,m = '#FF9800',15,'s'
        else: c,s,m = '#888888',4,'.'
        ax2.scatter(i, fot, c=c, s=s, marker=m, zorder=3)
ax2.set_xticks(range(len(res)))
ax2.set_xticklabels(names_r, rotation=90, fontsize=7, fontweight='bold')
ax2.set_ylabel('Oberton-Frequenz [Hz]', fontsize=10)
ax2.set_title('Oberton-Landkarte: Jeder Punkt = ein Oberton', fontsize=10, fontweight='bold')
ax2.set_ylim(200, 600); ax2.legend(loc='upper right', fontsize=7); ax2.grid(axis='y', alpha=0.3, linestyle=':')
ax2.text(len(res)-1, f_H+5, f'f_H = {f_H:.0f} Hz', fontsize=8, color='red', ha='right', fontweight='bold')
plt.tight_layout()
diag_path = '/home/claude/resonanz_diag.png'
plt.savefig(diag_path, dpi=200, bbox_inches='tight')
plt.close()

st.append(PageBreak())
st.append(Paragraph("Kapitel 11: Welche Basstoene sind durch die Kammerresonanz benachteiligt?", sty['Ch']))

st.append(Paragraph(
    "Die Helmholtz-Frequenz dieser Kammer liegt bei " + f"{f_H:.0f}" + " Hz. Die Frage ist: "
    "Welche Basstoene haben einen Oberton, der in die Naehe dieser Resonanz faellt? "
    "Denn wie in Kapitel 10 erlaeutert, entzieht eine resonante Kammer der Zunge Energie "
    "ueber die Serienkopplung. Nicht der Grundton muss betroffen sein - es genuegt, wenn "
    "der n-te Oberton (n x f<sub>Grundton</sub>) auf f<sub>H</sub> trifft.", sty['B']))

st.append(Paragraph(
    "Die Bandbreite der Helmholtz-Resonanz haengt von ihrer eigenen Guete Q<sub>H</sub> ab. "
    "Eine Holzkammer mit Klappe hat typisch Q<sub>H</sub> ~ 5-10. Bei Q<sub>H</sub> = 10 "
    "ist das kritische Band " + f"{f_H:.0f}" + " +/- " + f"{bw_narrow/2:.0f}" + " Hz, "
    "bei Q<sub>H</sub> = 5 erweitert es sich auf +/- " + f"{bw_wide/2:.0f}" + " Hz.", sty['B']))

# Diagramm einbetten
st.append(Spacer(1, 2*mm))
st.append(Image(diag_path, width=170*mm, height=130*mm))
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "Abbildung 4 zeigt das Ergebnis fuer alle Basstoene der temperierten Stimmung von "
    "C1 (32,7 Hz) bis D4 (293,7 Hz). Oberes Diagramm: Abstand des naechsten Obertons "
    "zu f<sub>H</sub>. Rot = kritisch (Oberton trifft Resonanz +/- " + f"{bw_narrow/2:.0f}" + " Hz), "
    "orange = unguenstig, gruen = unbeeinflusst. Die Zahl ueber jedem Balken gibt die "
    "Oberton-Nummer an. Unteres Diagramm: Oberton-Landkarte - jeder Punkt ist ein "
    "Oberton, rote Rauten liegen im Resonanzband.", sty['B']))

# Betroffene Noten als Tabelle
kritisch = [r for r in res if r['severity'] == 2]
# Gruppiere nach Muster
st.append(Paragraph("Auffaellige Muster:", sty['Sec']))

st.append(Paragraph(
    "<b>Unter ~70 Hz (C2) ist fast jeder Ton kritisch betroffen.</b> Der Grund: Bei "
    "tiefen Toenen liegen die Obertoene dicht beieinander (Abstand = Grundfrequenz). "
    "Ein Ton von 40 Hz hat Obertoene bei 40, 80, 120, ..., 440, 460, 480 Hz - "
    "das Raster ist so fein, dass immer ein Oberton in das Resonanzband faellt. "
    "Je tiefer der Ton, desto unvermeidlicher ist der Konflikt.", sty['Key']))

st.append(Paragraph(
    "<b>Besonders getroffen</b> sind Toene, deren Grundfrequenz ein ganzzahliger "
    "Teiler von f<sub>H</sub> ist:", sty['B']))

# Teiler-Tabelle
div_rows = [['f_H / n', 'f [Hz]', 'Naechste Note', 'Differenz', 'Oberton n']]
for div in [2, 3, 4, 5, 6, 7, 8, 9, 10]:
    fb = f_H / div
    if 30 <= fb <= 300:
        cl = min(res, key=lambda r: abs(r['f'] - fb))
        div_rows.append([f'{f_H:.0f}/{div}', f'{fb:.1f}', cl['name'], f'{abs(cl["f"]-fb):.1f} Hz', f'{div}'])
t_div = Table(div_rows, colWidths=[22*mm, 18*mm, 22*mm, 20*mm, 20*mm])
t_div.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_div)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "Die schaerfsten Treffer: <b>F#1</b> (46,2 Hz, 10. Oberton = 462,5 Hz, nur 1,3 Hz "
    "Abstand), <b>F#2</b> (92,5 Hz, 5. Oberton) und <b>B/Bb</b> in allen Oktaven "
    "(8./4./2. Oberton = 466,2 Hz, 5,0 Hz Abstand). Diese Toene werden bei dieser "
    "Kammergroesse systematisch benachteiligt.", sty['Warn']))

st.append(Paragraph(
    "Ab C3 (130,8 Hz) aufwaerts werden die Luecken zwischen den Obertoenen groesser, "
    "und die meisten Toene liegen ausserhalb des Resonanzbandes. Das erklaert, warum "
    "das Problem in der Praxis vor allem bei den tiefsten Basstoenen auffaellt - "
    "und warum erfahrene Instrumentenbauer die Kammergroesse fuer verschiedene "
    "Tonbereiche variieren.", sty['B']))

# ═══════════════════════════════════════════════════════════════
# Kap 12: Grenzen der Berechnung
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 12: Was Berechnungen leisten - und was nicht", sty['Ch']))

st.append(Paragraph(
    "Die in diesem Dokument vorgestellten Berechnungen stuetzen sich auf "
    "gut etablierte physikalische Modelle: Euler-Bernoulli-Balkentheorie fuer die "
    "Zungenfrequenz, Bernoulli-Gleichung fuer die Spaltstroemung, Helmholtz-Resonator "
    "fuer die Kammerakustik, Verlustbeiwerte aus der Stroemungsmechanik (Idelchik). "
    "Diese Modelle sind in ihrem Gueltigkeitsbereich zuverlaessig und liefern "
    "Ergebnisse, die in der richtigen Groessenordnung liegen.", sty['B']))

st.append(Paragraph(
    "Dennoch ersetzen sie keinen praktischen Test. Die Gruende dafuer sind nicht "
    "Schwaechen der Physik, sondern Vereinfachungen in den Randbedingungen:", sty['B']))

st.append(Paragraph(
    "<b>1. Idealisierte Geometrie:</b> Die Berechnung geht von exakt planen "
    "Flaechen, scharfen Kanten und gleichmaessigen Spalten aus. In der Realitaet "
    "hat Holz Fasern, Aluminium Bearbeitungsspuren, und jede Zunge ist minimal "
    "anders aufgebogen. Schon 0,1 mm Unterschied in der Aufbiegung aendert die "
    "effektive Spaltflaeche um rund 7 Prozent.", sty['B']))

st.append(Paragraph(
    "<b>2. Stationaere Stroemung:</b> Die Verlustberechnung nimmt stationaere "
    "Stroemung an - einen gleichmaessigen Luftstrom. Tatsaechlich pulsiert die "
    "Stroemung mit 50 Hz, weil die Zunge den Spalt periodisch oeffnet und "
    "schliesst. Die instationaeren Effekte (Anlaufwirbel, Schwinger-Stroemung-"
    "Kopplung, akustische Rueckwirkung) sind analytisch schwer zu fassen und "
    "aendern die Verlustbeiwerte um geschaetzt 10-30 Prozent. "
    "<b>Dies ist die gravierendste Einschraenkung:</b> Praktische Erfahrung zeigt, "
    "dass die Trennwandform das Einschwingverhalten entscheidend beeinflusst - "
    "ein Effekt, der im instationaeren Anlauf (Impulsantwort) liegt und von der "
    "stationaeren Berechnung grundsaetzlich nicht erfasst wird (siehe Kapitel 6).", sty['B']))

st.append(Paragraph(
    "<b>3. Daempfung:</b> Die Guete Q ist der empfindlichste Parameter der "
    "gesamten Analyse - sie bestimmt die Einschwingzeit direkt proportional. "
    "Q laesst sich aus den Materialwerten und der Geometrie nur grob abschaetzen "
    "(Kapitel 8). Der tatsaechliche Wert haengt von der konkreten Einspannung, "
    "dem Anpressdruck der Schrauben, der Oberflaechenguete an der Kontaktflaeche "
    "und dem Zustand der Zunge ab. Unterschiede von Faktor 2 zwischen zwei "
    "scheinbar identischen Zungen sind in der Praxis normal.", sty['B']))

st.append(Paragraph(
    "<b>4. Kammerakustik:</b> Das Helmholtz-Modell behandelt die Kammer als "
    "einzelnen Resonator mit einer einzigen Eigenfrequenz. Real hat die Kammer "
    "weitere Moden (Laengs-, Quer- und Hoehenmoden), die bei hoeheren Frequenzen "
    "auftreten. Zudem aendern sich die akustischen Eigenschaften mit dem "
    "Balgdruck, der Temperatur und der Luftfeuchtigkeit. Die Resonanzanalyse in "
    "Kapitel 11 zeigt die richtigen Tendenzen, aber die exakten Bandbreiten und "
    "Kopplungsstaerken koennen nur gemessen werden.", sty['B']))

st.append(Paragraph(
    "<b>5. Wechselwirkungen:</b> Jeder Effekt ist einzeln modelliert. In der "
    "Wirklichkeit beeinflussen sie sich gegenseitig: Die Zungenposition aendert "
    "die Stroemung, die Stroemung aendert den Druck in der Kammer, der Kammerdruck "
    "wirkt auf die Klappe zurueck, und die akustischen Wellen ueberlagern sich "
    "mit der Stroemung. Eine vollstaendige Loesung wuerde eine gekoppelte "
    "Fluid-Struktur-Akustik-Simulation (FSI) erfordern.", sty['B']))

st.append(Paragraph(
    "Was die Berechnungen <b>leisten</b>:", sty['Sec']))

st.append(Paragraph(
    "Sie identifizieren zuverlaessig, <b>welche Parameter wichtig sind</b> und "
    "auf welcher Ebene sie wirken. Das Ergebnis, dass der Spalt den stationaeren "
    "Gesamtwiderstand dominiert, ist robust. Aber die praktische Erfahrung zeigt "
    "ebenso robust, dass die Trennwandform das Einschwingverhalten entscheidend "
    "beeinflusst - ein Effekt, den die stationaere Berechnung prinzipiell nicht "
    "erfassen kann (Kapitel 6). Die Resonanzanalyse zeigt zuverlaessig, welche "
    "Tonhoehen prinzipiell gefaehrdet sind.", sty['B']))

st.append(Paragraph(
    "Die Berechnungen liefern damit eine <b>qualifizierte Einschaetzung</b>: "
    "Sie zeigen, dass Zungenguete, Aufbiegung und Einspannung die Grundansprache "
    "bestimmen. Sie ordnen die drei Trennwandvarianten in eine physikalisch "
    "begruendete Rangfolge (C besser als B besser als A - bestaetigt durch "
    "praktische Erfahrung). Und sie machen klar, warum tiefe Baesse traeger "
    "ansprechen als hohe Toene - das folgt direkt aus tau = Q/(pi*f). "
    "Wo die Berechnung an ihre Grenzen stoesst, ist die quantitative Vorhersage "
    "des Trennwand-Effekts: Dass die Wand wirkt, ist physikalisch begruendbar. "
    "Wie stark, laesst sich nur messen.", sty['B']))

st.append(Paragraph(
    "Fuer den Instrumentenbauer bedeutet das: <b>Die Berechnung sagt, wo man "
    "suchen soll. Der praktische Test sagt, was man findet.</b> Beides zusammen "
    "ist mehr als jedes fuer sich allein.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 13: Zusammenfassung (ex-Kap 12)
# ═══════════════════════════════════════════════════════════════
st.append(PageBreak())
st.append(Paragraph("Kapitel 13: Zusammenfassung", sty['Ch']))
st.append(Paragraph(f"<b>Stimmplatte:</b> Keilfoermig ({t_platte_einsp*1e3:.0f}->{t_platte_frei*1e3:.0f} mm). Am freien Ende wirkt der {t_platte_frei*1e3:.0f}-mm-Schlitz als kurzer Kanal mit eigenem Stroemungswiderstand (zeta ~ {z_schlitz:.3f}).", sty['Good']))
st.append(Paragraph(f"<b>Spalt:</b> Dreieckig (0->{h_aufbiegung*1e3:.1f} mm), S<sub>eff</sub> = {S_gap*1e6:.0f} mm². Dominiert den Gesamtwiderstand. Alle Varianten: v ~ {v1k:.0f} m/s, Q ~ {Q1k*1e6:.0f} ml/s bei 1 kPa.", sty['Good']))
st.append(Paragraph(f"<b>Trennwand-Material:</b> Folie und Furnier sind stroemungstechnisch gleichwertig. Furnier ist stabiler und akustisch neutraler. Folie spart Platz, kann aber flattern.", sty['Good']))
st.append(Paragraph(f"<b>Klappenorientierung:</b> Ausblas (zeta={zeta_klap_ausblas:.2f}) besser als Umlenk (zeta={zeta_klap_umlenk:.2f}). Ausblas nutzt Wandgeometrie voll, Umlenk nur bei Platzmangel.", sty['Good']))
st.append(Paragraph(f"<b>Trennwandform:</b> Praktisch entscheidend fuer die Ansprache (instationaerer Anlauf). C (Parabel) am saubersten, B (schraeg) guter Kompromiss, A (gerade) am einfachsten aber am traegsten. Stationaerer Durchfluss bei allen gleich - der Unterschied liegt in der Impulsantwort.", sty['Good']))
st.append(Paragraph(f"<b>Ansprache:</b> tau = Q/(pi*f). Bei Q=100: {tau_vals[0]*1e3:.0f} ms, bei Q=25: {tau_vals[2]*1e3:.0f} ms. Die Kammergeometrie beeinflusst nicht den stationaeren Durchfluss, aber entscheidend den instationaeren Anlauf (Impulsantwort).", sty['Good']))
st.append(Paragraph(f"<b>Akustik:</b> Helmholtz {f_H:.0f} Hz, Kammer akustisch klein.", sty['Good']))

doc.build(st)
print(f"  PDF: {path}")
shutil.copy(path, "/mnt/user-data/outputs/bass_50hz_v6.pdf")
print("  Fertig!")
