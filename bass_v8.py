#!/usr/bin/env python3
"""
Bass-Stimmzunge 50 Hz - v8 KORRIGIERT
======================================
KORREKTUREN gegenüber v6 (bereits in v7):
  1. Spaltfläche: Faktor 1/3 statt 1/2 (Cantilever-Biegeprofil quadratisch)
  2. Laminare Strömung: f = 64/Re für alle Kanalstrecken (Re << 2300)
  3. Schwellendruck: korrekte Gleichlast auf Cantilever -> 8EI*h/(W*L^4)
  4. Helmholtz: Plattendicke als Halslänge (Schlitz = akustischer Hals)
  5. Reynolds-Zahlen: ueberall ausgewiesen

NEU in v8 gegenüber v7:
  K3b. Schlitz-Eintrittsverlust: zeta_entry = 0.5 (scharfkantig) addiert
  K6.  Dynamische Spaltänderung: statische Durchbiegung unter Blasdruck
       -> iterative Lösung h_eff(dp), Grenzdruck dp_crit
       -> erklärt nichtlineares Verhalten bei hohen Drücken
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
# Cantilever-Biegeprofil: y(x) ~ y_tip * (x/L)²  (Näherung)
# S_gap = Integral[W * y(x) dx, 0, L] / L = W * h_max * (1/3)
# (v6 hatte Faktor 1/2 = lineares Profil, aber Cantilever biegt quadratisch)
S_gap = W_zunge * h_aufbiegung / 3  # v7: Faktor 1/3 statt 1/2
# Exakt: Cantilever y=(3Lx²-x³)/(2L³)*y_tip -> Integral/L = 3/8*y_tip
# Aber (x/L)² ist die bessere Näherung für den Spalt-Engpass -> 1/3

# Trennwand: zwei Varianten
t_furnier = 0.0005   # 0.5mm
t_folie = 0.0003     # 0.3mm
L_wall = 0.075
gap_fold = L_kam - L_wall  # 19mm

print("=" * 90)
print("  BASS-STIMMZUNGE 50 Hz - NEUBERECHNUNG v8 (korrigiert)")
print("  v7-Basis + Schlitz-Eintrittsverlust + dynamische Spaltänderung")
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
print(f"    -> Die Platte wirkt als kurzer Kanal für die Luft.")

# Schlitz als Kanal: mittlere Tiefe
t_schlitz_mean = t_platte_mean
S_schlitz_mean = W_schlitz * t_schlitz_mean  # Schlitz-Querschnitt

print(f"\n    Schlitz: {W_schlitz*1e3:.0f} mm breit x {t_schlitz_mean*1e3:.1f} mm tief (Mittel)")
print(f"    Schlitz-Querschnitt: {S_schlitz_mean*1e6:.1f} mm²")
print(f"    -> Deutlich größer als der Aufbiegungs-Spalt ({S_gap*1e6:.1f} mm²)")
print(f"    -> Der Spalt (Aufbiegung) bleibt der Engpass, nicht der Schlitz.")

# ═══════════════════════════════════════════════════════════════
# AUFBIEGUNG UND SPALT
# ═══════════════════════════════════════════════════════════════
print(f"\n  AUFBIEGUNG UND SPALT (v7: Faktor 1/3):")
print(f"    Aufbiegung am freien Ende:  {h_aufbiegung*1e3:.1f} mm")
print(f"    Spaltverlauf: Cantilever-Profil y(x) ~ y_tip*(x/L)²")
print(f"    Lokale Spaltfläche am freien Ende: {W_zunge*1e3:.0f} x {h_aufbiegung*1e3:.1f} = {W_zunge*h_aufbiegung*1e6:.0f} mm²")
print(f"    Effektive Gesamtfläche: S_gap = W x h_max / 3 = {S_gap*1e6:.1f} mm²")
print(f"    (v6 hatte: W x h_max / 2 = {W_zunge*h_aufbiegung/2*1e6:.1f} mm² -- zu gross)")
print(f"    (Zum Vergleich: Klappe eff. = {S_klap_eff*1e6:.0f} mm²)")
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

# Schwellendruck: v7 KORRIGIERT
# Gleichlast auf Cantilever: delta_tip = q*L^4 / (8*E*I), q = p*W
# -> p_min = 8*E*I*h_aufbiegung / (W*L^4)
# (v6 hatte: k_eff*h / (W*L*0.5) = 6EI*h/(W*L^4) -- falscher Hebelarm)
dp_min = 8 * E_st * I_reed * h_aufbiegung / (W_zunge * L_zunge**4)

print(f"\n  STIMMZUNGE (Euler-Bernoulli):")
print(f"    Dicke:      {h_reed*1e3:.3f} mm")
print(f"    Verif.:     f1 = {f_check:.2f} Hz OK")
print(f"    k_eff:      {k_eff:.2f} N/m")
print(f"    m_eff:      {m_eff*1e6:.1f} mg ({m_eff*1e3:.3f} g)")
print(f"    Befestigung: verschraubt (2 Schrauben) auf Stimmstock")
print(f"\n    Schwellendruck (v7: korrekte Gleichlast auf Cantilever):")
print(f"    DeltaP_min = 8*E*I*h / (W*L^4)")
print(f"           = 8*{E_st:.0e}*{I_reed:.3e}*{h_aufbiegung:.4f} / ({W_zunge}*{L_zunge}^4)")
print(f"           = {dp_min:.1f} Pa")
dp_min_v6 = k_eff * h_aufbiegung / (W_zunge * L_zunge * 0.5)
print(f"    (v6 hatte: {dp_min_v6:.1f} Pa -- Faktor {dp_min/dp_min_v6:.2f})")
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
        print(f"    ACHTUNG: FLUTTER: Folie kann bei Strömung vibrieren (Eigenfrequenz sehr hoch)")
        print(f"      -> kann Nebengeräusche erzeugen")
        print(f"      -> muss straff gespannt oder an Stützpunkten fixiert werden")

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
print(f"    Wie ein U-Turn auf einer Straße: hoher Energieverlust")
print(f"    zeta_Klappe = zeta_Ausblas + zeta_Umlenkung = {zeta_klap_ausblas:.3f} + {zeta_reversal:.2f} = {zeta_klap_umlenk:.3f}")
print(f"    -> {zeta_klap_umlenk/zeta_klap_ausblas:.1f}-facher Verlust gegenüber Ausblas!")

print(f"\n  WIRKUNG auf die Strömung:")
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
        ("B-Schräg", [
            ('Düsen-Eintritt', 0.15, A_B_klap),
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
        
        # Kanalreibung - v7: LAMINARE Berechnung f = 64/Re
        A_mean = A_k if 'Gerade' in wall_name else (A_B_klap+A_B_falt+A_A_falt+A_A_klap)/4
        D_h = 4*A_mean / (2*(A_mean/H_kam + H_kam))
        L_path = 2*L_wall + np.pi*W_k/2
        # Geschwindigkeit im Kanal bei 1000 Pa (Vorabschätzung via Spalt-Dominanz)
        v_gap_est = np.sqrt(2*1000/(rho*2.0))  # zeta_ges ~ 2
        v_kanal_est = v_gap_est * S_gap / A_mean
        Re_kanal = v_kanal_est * D_h / nu
        # Laminar: f = 64/Re (Rechteckkanal: Korrekturfaktor ~0.9, vernachlässigt)
        f_darcy = 64.0 / max(Re_kanal, 1)  # Sicherheit: min Re=1
        z_reib = f_darcy * L_path / D_h * (S_gap/A_mean)**2
        zeta_total += z_reib
        
        configs[cname] = {'zeta': zeta_total, 'losses': all_losses, 'orient': orient_name, 
                          'wall': wall_name, 'Re_kanal': Re_kanal, 'f_darcy': f_darcy,
                          'z_reib': z_reib}

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
    # Re im Spalt: charakteristische Länge = mittlere Spalthöhe
    h_mean_gap = h_aufbiegung / 3  # v7: bei quadrat. Profil ist h_mean = h_max/3
    Re_gap = v * h_mean_gap / nu
    # Re im Schlitz: char. Länge = D_h des Schlitzes
    v_schlitz = v * S_gap / S_schlitz_mean
    D_h_schlitz_loc = 2*W_schlitz*t_platte_mean/(W_schlitz+t_platte_mean)
    Re_schlitz = v_schlitz * D_h_schlitz_loc / nu
    cfg['v1000'] = v
    cfg['Q1000'] = Q
    cfg['Re_gap'] = Re_gap
    cfg['Re_schlitz'] = Re_schlitz
    print(f"  {cname:<28} {z:>10.6f} {v:>8.2f} {Q*1e6:>8.1f}  Re_Sp={Re_gap:>5.0f}  Re_Schl={Re_schlitz:>5.0f}  Re_Kan={cfg['Re_kanal']:>4.0f}  f={cfg['f_darcy']:.3f}")

# Vergleich Folie vs Furnier
print(f"\n  FOLIE vs. FURNIER Differenz:")
# Folie: dünnere Wand -> breitere Kanäle -> geringfügig weniger Reibung
W_k_folie = (B_kam - t_folie) / 2
W_k_furnier = (B_kam - t_furnier) / 2
print(f"    Kanalbreite Folie:   {W_k_folie*1e3:.2f} mm")
print(f"    Kanalbreite Furnier: {W_k_furnier*1e3:.2f} mm")
print(f"    Differenz:           {(W_k_folie-W_k_furnier)*1e3:.2f} mm ({(W_k_folie-W_k_furnier)/W_k_furnier*100:.1f}%)")
print(f"    -> Strömungstechnisch vernachlässigbar")
print(f"    -> Hauptunterschied: Oberfläche (Reibung) und Stabilität (Flutter)")
print(f"    Reibungsfaktor Folie:   f ~ 0.02 (glatt)")
print(f"    Reibungsfaktor Furnier: f ~ 0.03 (rau)")
print(f"    Kanalreibung ist aber nur {z_reib:.2e} des Gesamtverlusts -> irrelevant")

# Helmholtz - v7: Plattendicke als Halslänge
# Der akustische Hals ist der Schlitz in der Stimmplatte
# Halsquerschnitt = Schlitz (W_schlitz * t_platte) oder Spalt (S_gap), je nach Betrachtung
# Halslänge = mittlere Plattendicke + Endkorrektur
V_eff = V_brutto - L_wall*t_furnier*H_kam
# Variante 1: Klappe als Hals (wie v6)
r_klap = np.sqrt(S_klap_eff/np.pi)
L_eff_klap = 0.005 + 1.6*r_klap
f_H_klap = (c_s/(2*np.pi))*np.sqrt(S_klap_eff/(V_eff*L_eff_klap))
# Variante 2: Schlitz als Hals (v7 NEU - physikalisch korrekter)
# Halsquerschnitt = Schlitz-Querschnitt (W_schlitz * t_platte_mean)
# Halslänge = Plattendicke (der Kanal, den die Luft durchqueren muss)
# Endkorrektur: 0.8 * sqrt(S/pi) beidseitig
S_hals = S_schlitz_mean  # W_schlitz * t_platte_mean = 67.5 mm²
r_hals = np.sqrt(S_hals/np.pi)
L_hals = t_platte_mean + 2 * 0.8 * r_hals  # Plattendicke + 2x Endkorrektur
f_H = (c_s/(2*np.pi))*np.sqrt(S_hals/(V_eff*L_hals))

print(f"\n  HELMHOLTZ-RESONANZ (v7 korrigiert):")
print(f"    Variante 1 - Klappe als Hals:  f_H = {f_H_klap:.0f} Hz")
print(f"      S_Hals = {S_klap_eff*1e6:.0f} mm², L_Hals = {L_eff_klap*1e3:.1f} mm")
print(f"    Variante 2 - Schlitz als Hals: f_H = {f_H:.0f} Hz  <-- v7")
print(f"      S_Hals = {S_hals*1e6:.1f} mm², L_Hals = {L_hals*1e3:.1f} mm")
print(f"      (Plattendicke {t_platte_mean*1e3:.1f} mm + 2x Endkorrektur {0.8*r_hals*1e3:.1f} mm)")
print(f"    Der Schlitz ist der engere Hals und bestimmt die Resonanz")

# ═══════════════════════════════════════════════════════════════
# VORKAMMER-ANALYSE (v8 neu)
# ═══════════════════════════════════════════════════════════════
# Zusätzliche Kammer ZWISCHEN Klappe und Hauptkammer
# Abmessungen: 48 x 20 mm Querschnitt, 70 mm lang
# Kopplung zur Hauptkammer über die Klappenöffnung

print(f"\n{'='*90}")
print(f"  VORKAMMER-ANALYSE (v8 neu)")
print(f"{'='*90}")

# Vorkammer-Geometrie
L_vk = 0.070   # 70 mm Länge
B_vk = 0.048   # 48 mm Breite (= Kammerbreite)
H_vk = 0.020   # 20 mm Höhe
V_vk = L_vk * B_vk * H_vk  # Volumen
S_vk_quer = B_vk * H_vk     # Querschnitt

# Kopplung Vorkammer -> Hauptkammer:
# Die Klappe öffnet jetzt in die Vorkammer statt direkt in die Hauptkammer.
# Zwischen Vorkammer und Hauptkammer braucht es eine Öffnung.
# Annahme: Die Öffnung hat dieselben Abmessungen wie die Klappenöffnung
# (die Klappe sitzt jetzt am EINGANG der Vorkammer, am Ausgang ist eine
# feste Öffnung zur Hauptkammer)
S_kopplung = S_klap_eff  # ~244 mm² (gleiche Abmessungen wie Klappe)

# Kopplungshals: Wanddicke zwischen Vorkammer und Hauptkammer
t_kopplung = 0.005  # 5 mm Wandstärke (Holz)
r_kopplung = np.sqrt(S_kopplung / np.pi)
L_eff_kopplung = t_kopplung + 2 * 0.8 * r_kopplung  # mit Endkorrektur

print(f"\n  Vorkammer-Geometrie:")
print(f"    Breite x Höhe x Länge: {B_vk*1e3:.0f} x {H_vk*1e3:.0f} x {L_vk*1e3:.0f} mm")
print(f"    Volumen V_vk: {V_vk*1e6:.1f} cm³")
print(f"    Querschnitt: {S_vk_quer*1e6:.0f} mm²")
print(f"    Kopplung zur Hauptkammer: S = {S_kopplung*1e6:.0f} mm² (= Klappenöffnung)")
print(f"    Kopplungshals-Länge: {L_eff_kopplung*1e3:.1f} mm (inkl. Endkorrektur)")
print(f"\n  Zum Vergleich:")
print(f"    Hauptkammer V_hk: {V_eff*1e6:.1f} cm³")
print(f"    Gesamtvolumen: {(V_vk+V_eff)*1e6:.1f} cm³ ({V_vk/V_eff*100:.0f}% Zunahme)")

# ═══════════════════════════════════════════════════════════════
# MODELL 1: Einzelkammer f_H mit Gesamtvolumen
# ═══════════════════════════════════════════════════════════════
# Wenn Vorkammer und Hauptkammer akustisch stark gekoppelt sind,
# wirken sie zusammen als eine grössere Kammer

V_gesamt = V_vk + V_eff
f_H_gesamt = (c_s/(2*np.pi))*np.sqrt(S_hals/(V_gesamt*L_hals))

print(f"\n  MODELL 1: Gesamtvolumen (stark gekoppelt)")
print(f"    V_gesamt = {V_gesamt*1e6:.1f} cm³")
print(f"    f_H (Schlitz als Hals) = {f_H_gesamt:.0f} Hz  (ohne VK: {f_H:.0f} Hz)")
print(f"    Absenkung: {(1-f_H_gesamt/f_H)*100:.1f}%")

# ═══════════════════════════════════════════════════════════════
# MODELL 2: Gekoppelte Helmholtz-Resonatoren
# ═══════════════════════════════════════════════════════════════
# System: Balg -> VK (V1) -> Kopplungshals (S_k, L_k) -> HK (V2) -> Schlitz -> Zunge
#
# Akustische Massen: M_k = rho * L_k / S_k (Kopplungshals)
#                    M_s = rho * L_s / S_s (Schlitz)
# Akustische Compliances: C1 = V1 / (rho * c²) (Vorkammer)
#                          C2 = V2 / (rho * c²) (Hauptkammer)
#
# Eigenfrequenzen des gekoppelten Systems:
# | M_k  0  | |q1_dd| + |1/C1+1/C2  -1/C2| |q1| = 0
# | 0   M_s | |q2_dd|   |-1/C2       1/C2 | |q2|
#
# Das ergibt eine 2x2-Eigenwertgleichung

M_k = rho * L_eff_kopplung / S_kopplung  # akustische Masse Kopplungshals
M_s = rho * L_hals / S_hals              # akustische Masse Schlitz
C1 = V_vk / (rho * c_s**2)              # Compliance Vorkammer
C2 = V_eff / (rho * c_s**2)             # Compliance Hauptkammer

print(f"\n  MODELL 2: Gekoppelte Helmholtz-Resonatoren")
print(f"    Akustische Massen:")
print(f"      M_Kopplung = rho*L/S = {M_k:.2f} kg/m4")
print(f"      M_Schlitz  = rho*L/S = {M_s:.2f} kg/m4")
print(f"    Akustische Compliances:")
print(f"      C_VK = V/(rho*c²) = {C1:.3e} m³/Pa")
print(f"      C_HK = V/(rho*c²) = {C2:.3e} m³/Pa")

# Eigenfrequenzen: det(K - omega² M) = 0
# K11 = 1/C1 + 1/C2, K12 = K21 = -1/C2, K22 = 1/C2
# M11 = M_k, M22 = M_s (diagonal)
#
# (K11 - w²M_k)(K22 - w²M_s) - K12² = 0
# M_k*M_s * w⁴ - (K11*M_s + K22*M_k) * w² + (K11*K22 - K12²) = 0

K11 = 1/C1 + 1/C2
K12 = -1/C2
K22 = 1/C2

a_coeff = M_k * M_s
b_coeff = -(K11 * M_s + K22 * M_k)
c_coeff = K11 * K22 - K12**2  # = 1/(C1*C2)

discriminant = b_coeff**2 - 4*a_coeff*c_coeff
w2_1 = (-b_coeff - np.sqrt(discriminant)) / (2*a_coeff)
w2_2 = (-b_coeff + np.sqrt(discriminant)) / (2*a_coeff)

f_mode1 = np.sqrt(w2_1) / (2*np.pi)
f_mode2 = np.sqrt(w2_2) / (2*np.pi)

# Einzelresonanzen (zum Vergleich)
# VK allein mit Kopplungshals: f_vk = c/(2pi) * sqrt(S_k/(V1*L_k))
f_vk_allein = (c_s/(2*np.pi)) * np.sqrt(S_kopplung / (V_vk * L_eff_kopplung))
# HK allein mit Schlitz: = f_H (bereits berechnet)

print(f"\n    Einzelresonanzen (ungekoppelt):")
print(f"      VK mit Kopplungshals: f_VK = {f_vk_allein:.0f} Hz")
print(f"      HK mit Schlitz:       f_HK = {f_H:.0f} Hz")

print(f"\n    Gekoppelte Eigenfrequenzen:")
print(f"      Mode 1 (niedrig):  f_1 = {f_mode1:.0f} Hz")
print(f"      Mode 2 (hoch):     f_2 = {f_mode2:.0f} Hz")

# Modenformen: Wer schwingt wie?
# Eigenvektor für Mode 1: (K - w1² M) * v = 0
# (K11 - w1²*M_k) * v1 + K12 * v2 = 0
# -> v2/v1 = -(K11 - w1²*M_k) / K12
ratio_1 = -(K11 - w2_1*M_k) / K12
ratio_2 = -(K11 - w2_2*M_k) / K12

print(f"\n    Modenformen (Druckverhältnis VK/HK):")
print(f"      Mode 1: p_VK/p_HK = {ratio_1:.2f}  ", end='')
if ratio_1 > 0:
    print("(gleichphasig - beide Kammern schwingen zusammen)")
else:
    print("(gegenphasig - Kammern schwingen gegeneinander)")
print(f"      Mode 2: p_VK/p_HK = {ratio_2:.2f}  ", end='')
if ratio_2 > 0:
    print("(gleichphasig)")
else:
    print("(gegenphasig - Kammern schwingen gegeneinander)")

# ═══════════════════════════════════════════════════════════════
# WIRKUNG AUF DIE ZUNGE
# ═══════════════════════════════════════════════════════════════
# Phase der Impedanz am Schlitz für beide Konfigurationen
print(f"\n  WIRKUNG AUF DIE ZUNGE (Phasenverschiebung):")
print(f"  {'Oberton':>8} {'f [Hz]':>8} {'OHNE VK':>16} {'MIT VK Mode 1':>16} {'MIT VK Mode 2':>16}")
print(f"  {'':>8} {'':>8} {'Phase [°]':>16} {'Phase [°]':>16} {'Phase [°]':>16}")
print(f"  {'-'*70}")

Q_H_est = 7  # Kammer-Güte
phase_vk_results = []

for n in [1, 2, 3, 5, 6, 8, 10, 15, 20]:
    f_n = n * 50
    
    # OHNE Vorkammer: einfacher Helmholtz
    x_ohne = f_n / f_H
    phi_ohne = np.degrees(np.arctan(Q_H_est * (x_ohne - 1/x_ohne)))
    
    # MIT Vorkammer: Impedanz am Schlitz mit Dämpfung
    omega = 2 * np.pi * f_n
    
    # Dämpfungswiderstände (akustisch): R = omega_res * M / Q
    R_s = 2*np.pi*f_H * M_s / Q_H_est        # Schlitz-Dämpfung
    R_k = 2*np.pi*f_vk_allein * M_k / Q_H_est # Kopplungshals-Dämpfung
    
    # Impedanz am Schlitz:
    # Z_Schlitz = (R_s + j*w*M_s) + 1/(j*w*C2)  ← HK-Seite
    # Z_VK_arm = (R_k + j*w*M_k) + 1/(j*w*C1)  ← VK+Kopplungshals
    # Parallel: Z_HK || Z_VK_arm, dann Schlitz-Impedanz in Serie
    
    Z_schlitz_arm = R_s + 1j*omega*M_s
    Z_HK_cap = 1/(1j*omega*C2)
    Z_vk_arm = R_k + 1j*omega*M_k + 1/(1j*omega*C1)
    
    # Am Schlitz: Zunge sieht M_s + (C2 parallel mit VK-Arm)
    Y_parallel = 1j*omega*C2 + 1/Z_vk_arm
    Z_parallel = 1 / Y_parallel
    Z_total = Z_schlitz_arm + Z_parallel
    
    phi_mit = np.degrees(np.angle(Z_total))
    
    phase_vk_results.append((n, f_n, phi_ohne, phi_mit))
    
    # Nächste Mode-Resonanzen
    dist1 = abs(f_n - f_mode1)
    dist2 = abs(f_n - f_mode2)
    
    print(f"  {n:>8} {f_n:>8} {phi_ohne:>16.1f} {phi_mit:>16.1f}")

# Anti-Resonanz: wo Y_parallel = 0, d.h. 1/(j*w*C2) + 1/Z_vk = 0
# -> j*w*C2 + 1/(j*w*M_k + 1/(j*w*C1)) = 0
# Bei Anti-Resonanz wird die Impedanz am Schlitz maximal -> Druck am Schlitz maximal
# Das passiert wenn die VK als Tilger wirkt
# f_anti ≈ f_vk_allein (VK-Eigenresonanz)
f_anti = f_vk_allein  # Näherung
print(f"\n    Anti-Resonanz (VK als Tilger): ~{f_anti:.0f} Hz")
print(f"    Bei dieser Frequenz blockiert die VK den Durchfluss -> maximaler Druck in HK")

# Zusammenfassung
print(f"\n  ZUSAMMENFASSUNG VORKAMMER-WIRKUNG:")
print(f"    OHNE VK: Ein Helmholtz bei f_H = {f_H:.0f} Hz")
print(f"    MIT  VK: Zwei Moden bei f_1 = {f_mode1:.0f} Hz und f_2 = {f_mode2:.0f} Hz")
print(f"             Anti-Resonanz bei ~{f_anti:.0f} Hz (VK als Tilger)")
print(f"    Verstimmung: f_1 liegt {(1-f_mode1/f_H)*100:.0f}% unter dem Original")
print(f"    Neues Volumen: {V_gesamt*1e6:.1f} cm³ (+{V_vk/V_eff*100:.0f}%)")
print(f"")
print(f"    Für 50-Hz-Zunge (Obertöne 250/300 Hz nahe f_H = {f_H:.0f} Hz):")

# Prüfe ob die Modenaufspaltung die 5./6. OT-Kopplung verbessert oder verschlechtert
for n_check in [5, 6]:
    f_check_val = n_check * 50
    x_o = f_check_val / f_H
    phi_o = np.degrees(np.arctan(Q_H_est * (x_o - 1/x_o)))
    # Phase mit VK (gedämpft)
    omega_c = 2 * np.pi * f_check_val
    R_s_c = 2*np.pi*f_H * M_s / Q_H_est
    R_k_c = 2*np.pi*f_vk_allein * M_k / Q_H_est
    Z_vk_c = R_k_c + 1j*omega_c*M_k + 1/(1j*omega_c*C1)
    Y_par_c = 1j*omega_c*C2 + 1/Z_vk_c
    Z_par_c = 1 / Y_par_c
    Z_tot_c = R_s_c + 1j*omega_c*M_s + Z_par_c
    phi_m = np.degrees(np.angle(Z_tot_c))
    print(f"      {n_check}. OT ({f_check_val} Hz): Phase ohne VK = {phi_o:.1f}°, mit VK = {phi_m:.1f}°")
    if abs(phi_m) < abs(phi_o):
        print(f"        -> VK VERBESSERT die Kopplung (Phase näher an 0°)")
    else:
        print(f"        -> VK VERSCHLECHTERT die Kopplung (Phase weiter von 0°)")

# Stationärer Durchfluss: zusätzlicher Widerstand
# Kopplungshals als zusätzlicher Verlust
zeta_kopplung = 1.5  # Eintritt + Austritt (scharfkantig)
z_kopplung_scaled = zeta_kopplung * (S_gap / S_kopplung)**2
print(f"\n    Stationärer Durchfluss:")
print(f"      Zusätzlicher Verlust: zeta_Kopplung = {zeta_kopplung} (Ein+Austritt)")
print(f"      Auf Spalt bezogen: {z_kopplung_scaled:.2e}")
print(f"      -> VERNACHLÄSSIGBAR (S_Kopplung/S_Spalt = {S_kopplung/S_gap:.0f})")

# Akustische Laufzeit in der Vorkammer
t_vk = L_vk / c_s
print(f"\n    Akustische Laufzeit VK: {t_vk*1e3:.2f} ms (L={L_vk*1e3:.0f} mm / c={c_s:.0f} m/s)")
print(f"    Rundlauf VK: {2*t_vk*1e3:.2f} ms")
print(f"    Bei 50 Hz: {1/(50*2*t_vk):.0f} Rundläufe/Periode -> quasi-statisch")

# Schlitz als zusätzlicher Hals: mittlere Tiefe des Schlitzes
L_eff_spalt = t_platte_mean + 0.8*np.sqrt(S_gap/np.pi)

print(f"\n  WIRKUNG DER PLATTENDICKE (v8: laminar + Eintrittsverlust):")
print(f"    Am freien Ende: 13 mm -> Luft durchquert 13 mm langen Kanal im Schlitz")
print(f"    Am Einspann-Ende: 2 mm -> Kanal nur 2 mm lang")
D_h_schlitz = 2*W_schlitz*h_aufbiegung/(W_schlitz+h_aufbiegung)  # am Tip
# v8: Geschwindigkeit im Schlitz = Gap-Geschwindigkeit * (S_gap / S_schlitz)
v_gap_est = np.sqrt(2*1000/(rho*2.0))  # ca. 29 m/s am Spalt
v_schlitz_est = v_gap_est * S_gap / S_schlitz_mean  # viel langsamer (großer Querschnitt)
Re_schlitz_tip = v_schlitz_est * D_h_schlitz / nu
f_schlitz = 64.0 / max(Re_schlitz_tip, 1)  # laminar (Re << 2300)
# [K3b] Eintrittsverlust: scharfkantiger Eintritt in kurzen Kanal
zeta_schlitz_entry = 0.5  # scharfkantiger Eintritt (Idelchik)
z_schlitz_friction = f_schlitz * t_platte_frei / D_h_schlitz
z_schlitz = z_schlitz_friction + zeta_schlitz_entry  # v8: Gesamt = Reibung + Eintritt
print(f"    D_h am Tip: {D_h_schlitz*1e3:.2f} mm")
print(f"    v_Schlitz: {v_schlitz_est:.2f} m/s  (v_Spalt * S_gap/S_schlitz = {v_gap_est:.0f} * {S_gap*1e6:.0f}/{S_schlitz_mean*1e6:.0f})")
print(f"    Re_Schlitz am Tip: {Re_schlitz_tip:.0f}  {'(laminar)' if Re_schlitz_tip < 2300 else '(TURBULENT!)'}")
print(f"    f_Darcy (laminar): {f_schlitz:.3f}  (v6 hatte: 0.030 turbulent)")
print(f"    [K3b] Schlitz-Eintrittsverlust: zeta_entry = {zeta_schlitz_entry}")
print(f"    L/D_h = {t_platte_frei*1e3:.0f}/{D_h_schlitz*1e3:.2f} = {t_platte_frei/D_h_schlitz:.1f} (kurzer Kanal!)")
print(f"    zeta_Schlitz = zeta_entry + f*L/D_h = {zeta_schlitz_entry} + {z_schlitz_friction:.3f} = {z_schlitz:.3f}")
print(f"    (v7 hatte: {z_schlitz_friction:.3f} -- ohne Eintrittsverlust)")

# ═══════════════════════════════════════════════════════════════
# [K6] DYNAMISCHE SPALTAENDERUNG UNTER BLASDRUCK
# ═══════════════════════════════════════════════════════════════
# VOR Schwingungsbeginn: Der statische Blasdruck biegt die Zunge nach unten.
# Dadurch verkleinert sich der Spalt -> weniger Durchfluss -> höhere
# Geschwindigkeit nötig -> stärkerer Bernoulli-Sog -> weitere Biegung.
# Das ist ein selbstverstärkender Effekt mit Grenzdruck.
#
# Cantilever unter Gleichlast q = dp * W_zunge:
#   w_tip = q * L^4 / (8 * E * I) = dp * W * L^4 / (8*E*I)
# Die Zunge biegt sich NACH UNTEN (in den Schlitz hinein), d.h.
# der Spalt (Aufbiegung über der Platte) wird KLEINER:
#   h_eff = h_aufbiegung - w_tip(dp)
#
# Grenzdruck: wo h_eff -> 0 (Zunge liegt flach auf Platte)
#   dp_crit = 8*E*I*h_aufbiegung / (W*L^4) = dp_min  (!)
# Das ist identisch mit dem Schwellendruck -- physikalisch logisch:
# Der Druck, der die Zunge gerade flach drückt, ist derselbe,
# der nötig ist, um die Aufbiegung statisch zu ueberwinden.

print(f"\n{'='*90}")
print(f"  [K6] DYNAMISCHE SPALTAENDERUNG UNTER BLASDRUCK (NEU in v8)")
print(f"{'='*90}")

# Statische Durchbiegung als Funktion des Drucks
def w_static_tip(dp_val):
    """Cantilever-Tip-Auslenkung unter Gleichlast dp*W"""
    return dp_val * W_zunge * L_zunge**4 / (8 * E_st * I_reed)

# Effektive Spalthöhe
h_gap_min = 0.02e-3  # 20 um Minimalspalt (Oberflächenrauigkeit)

def h_eff_static(dp_val):
    """Effektive Aufbiegung nach statischer Druckverformung"""
    w = w_static_tip(dp_val)
    return max(h_aufbiegung - w, h_gap_min)

def S_gap_eff(dp_val):
    """Effektive Spaltfläche unter Druck (mit 1/3-Faktor)"""
    return W_zunge * h_eff_static(dp_val) / 3

# Grenzdruck: h_eff -> h_gap_min
dp_crit = 8 * E_st * I_reed * (h_aufbiegung - h_gap_min) / (W_zunge * L_zunge**4)

print(f"\n  Statische Durchbiegung unter Blasdruck:")
print(f"    w_tip(dp) = dp * W * L^4 / (8*E*I)")
print(f"    w_tip(dp) = dp * {W_zunge*1e3:.0f}e-3 * ({L_zunge*1e3:.0f}e-3)^4 / (8 * {E_st:.0e} * {I_reed:.3e})")
print(f"    w_tip(1000 Pa) = {w_static_tip(1000)*1e3:.4f} mm")
print(f"    w_tip(dp_min={dp_min:.0f} Pa) = {w_static_tip(dp_min)*1e3:.4f} mm = h_aufbiegung  (Identitaet!)")
print(f"")
print(f"  Grenzdruck (Spalt schließt):")
print(f"    dp_crit = {dp_crit:.1f} Pa ({dp_crit/100:.1f} mbar)")
print(f"    dp_min  = {dp_min:.1f} Pa ({dp_min/100:.1f} mbar)  (Schwellendruck)")
print(f"    -> dp_crit ≈ dp_min: Der Spalt schließt genau beim Schwellendruck!")
print(f"    -> Das ist physikalisch konsistent: Beide berechnen denselben Effekt.")

# Iterative Lösung: Bei gegebenem dp, wie gross ist der effektive Spalt
# unter Berücksichtigung der Strömungs-Rückwirkung?
# Einfaches Modell: statischer Druck biegt Zunge, Bernoulli-Sog verstärkt
# Für die Grundanalyse genügt die statische Lösung (ohne Bernoulli-Feedback).

# Druckscan: h_eff, S_gap, Q als Funktion von dp
print(f"\n  Druckscan (statische Durchbiegung):")
print(f"  {'dp [Pa]':>10} {'dp [mbar]':>10} {'w_tip [mm]':>10} {'h_eff [mm]':>10} {'S_gap [mm2]':>10} {'S/S_ruhe':>8}")
print(f"  {'-'*62}")

dp_scan_vals = [100, 200, 500, 1000, 1500, 2000, 3000, 4000, 5000, dp_min, dp_crit]
dp_scan_vals = sorted(set([int(round(v)) for v in dp_scan_vals]))

for dp_val in dp_scan_vals:
    w = w_static_tip(dp_val)
    h_e = h_eff_static(dp_val)
    S_e = S_gap_eff(dp_val)
    ratio = S_e / S_gap if S_gap > 0 else 0
    marker = "  <-- dp_min=dp_crit" if abs(dp_val - dp_min) < 1 else ""
    print(f"  {dp_val:>10.0f} {dp_val/100:>10.1f} {w*1e3:>10.4f} {h_e*1e3:>10.4f} {S_e*1e6:>10.2f} {ratio:>8.3f}{marker}")

# Effektiver Volumenstrom mit dynamischer Spaltänderung
print(f"\n  Volumenstrom mit dynamischer Spaltfläche:")
print(f"  {'dp [Pa]':>10} {'S_stat [mm2]':>12} {'S_dyn [mm2]':>12} {'Q_stat [ml/s]':>14} {'Q_dyn [ml/s]':>14} {'Q_dyn/Q_stat':>12}")
print(f"  {'-'*78}")

zeta_ref = configs['C-Parabel + Ausblas']['zeta']  # Referenz-Verlustbeiwert

for dp_val in [200, 500, 1000, 2000, 3000, 5000]:
    # Statisch (konstanter Spalt)
    v_stat = np.sqrt(2*dp_val/(rho*(1+zeta_ref)))
    Q_stat = v_stat * S_gap
    
    # Dynamisch (Spalt verkleinert sich)
    S_dyn = S_gap_eff(dp_val)
    # Verlustbeiwert skaliert: zeta_total bleibt gleich bezogen auf v_Spalt
    # Aber S ist jetzt kleiner -> Q = v * S_dyn
    v_dyn = np.sqrt(2*dp_val/(rho*(1+zeta_ref)))  # v ändert sich kaum
    Q_dyn = v_dyn * S_dyn
    
    print(f"  {dp_val:>10.0f} {S_gap*1e6:>12.2f} {S_dyn*1e6:>12.2f} {Q_stat*1e6:>14.1f} {Q_dyn*1e6:>14.1f} {Q_dyn/Q_stat if Q_stat>0 else 0:>12.3f}")

# Physikalische Interpretation
print(f"\n  PHYSIKALISCHE INTERPRETATION:")
print(f"    1. Bei niedrigem Druck (dp << dp_crit): Spalt fast unverändert, S ≈ S_ruhe")
print(f"    2. Bei dp -> dp_crit: Spalt schließt sich rapide (nichtlinear)")
print(f"    3. Bei dp > dp_crit: Zunge liegt auf Platte, kein statischer Durchfluss")
print(f"       -> Aber: Die Zunge beginnt zu SCHWINGEN, bevor dp_crit erreicht wird!")
print(f"       -> Dynamische Bernoulli-Rückkopplung setzt bei ca. 0.3-0.5 x dp_min ein")
print(f"    4. Typischer Balgdruck 500-5000 Pa -> bei 50 Hz liegt dp_min bei {dp_min:.0f} Pa")
dp_dynamic_est = 0.4 * dp_min
print(f"    5. Geschätzter dynamischer Einsetzpunkt: ~{dp_dynamic_est:.0f} Pa ({dp_dynamic_est/100:.0f} mbar)")
print(f"       (Bernoulli-Instabilität setzt ein, bevor Spalt statisch schließt)")

# ═══════════════════════════════════════════════════════════════
# [NEU] AKUSTISCHE PHASENANALYSE DER KAMMER
# ═══════════════════════════════════════════════════════════════
# Die entscheidende Frage ist nicht "saubere Strömung", sondern:
# Kommt der Energienachschub mit der richtigen PHASENLAGE an der Zunge an?
#
# Modell: Transmission Line (TL) - jeder Kammerabschnitt hat eine
# akustische Impedanz Z = rho*c / A. An Querschnittsänderungen
# entstehen Reflexionen mit Koeffizient Gamma = (Z2-Z1)/(Z2+Z1).
# Der reflektierte Anteil kommt mit Phasenverschiebung zurück.

print(f"\n{'='*90}")
print(f"  AKUSTISCHE PHASENANALYSE DER KAMMER")
print(f"{'='*90}")

# Akustische Laufzeiten
# Weg: Spalt -> Kanal A -> Faltung -> Kanal B -> Klappe -> zurück
L_kanal_A = L_wall  # 75 mm
L_kanal_B = L_wall  # 75 mm
L_faltung = np.pi * gap_fold / 2  # Halbkreis um die Faltung ~ 30 mm
L_gesamt_einfach = L_kanal_A + L_faltung + L_kanal_B  # ~ 180 mm
L_gesamt_rund = 2 * L_gesamt_einfach  # hin und zurück

t_einfach = L_gesamt_einfach / c_s
t_rund = L_gesamt_rund / c_s
T_50 = 1.0 / 50.0  # Schwingungsperiode bei 50 Hz

print(f"\n  Akustische Laufzeiten:")
print(f"    Kanal A:          {L_kanal_A*1e3:.0f} mm  -> {L_kanal_A/c_s*1e3:.3f} ms")
print(f"    Faltung (pi*r/2): {L_faltung*1e3:.0f} mm  -> {L_faltung/c_s*1e3:.3f} ms")
print(f"    Kanal B:          {L_kanal_B*1e3:.0f} mm  -> {L_kanal_B/c_s*1e3:.3f} ms")
print(f"    Einfach (Spalt->Klappe): {L_gesamt_einfach*1e3:.0f} mm  -> {t_einfach*1e3:.3f} ms")
print(f"    Rundlauf (Spalt->Klappe->Spalt): {L_gesamt_rund*1e3:.0f} mm  -> {t_rund*1e3:.3f} ms")
print(f"    Schwingungsperiode 50 Hz: {T_50*1e3:.0f} ms")
print(f"    Rundläufe pro Periode: {T_50/t_rund:.1f}")
print(f"    -> Kammer ist akustisch SEHR klein gegenüber Grundton")

# Phasenverschiebung als Funktion der Frequenz
# phi = 2*pi*f * t_rund = omega * t_rund
# Konstruktive Verstärkung wenn phi = n * 2*pi (n = 0, 1, 2, ...)
# d.h. f = n / t_rund = n * c_s / (2*L_gesamt)
# Destruktive Dämpfung wenn phi = (n+0.5) * 2*pi

f_konstruktiv_1 = c_s / (2 * L_gesamt_einfach)  # erste konstruktive Frequenz
f_destruktiv_1 = c_s / (4 * L_gesamt_einfach)   # erste destruktive Frequenz

print(f"\n  Phasenbeziehung (Stehwellen in der Kammer):")
print(f"    Erste konstruktive Frequenz: f_k1 = c/(2L) = {f_konstruktiv_1:.0f} Hz")
print(f"    Erste destruktive Frequenz:  f_d1 = c/(4L) = {f_destruktiv_1:.0f} Hz")
print(f"    Helmholtz-Resonanz:          f_H = {f_H:.0f} Hz")

# Transmission-Line Analyse: Reflexionskoeffizienten
# Z = rho*c / A (akustische Impedanz eines Kanalabschnitts)
# An Querschnittsänderung A1 -> A2:
#   Gamma = (Z2-Z1)/(Z2+Z1) = (A1-A2)/(A1+A2)  (Impedanz ~ 1/A)
#   |Gamma|^2 = reflektierte Leistung
#   Transmission = 1 - |Gamma|^2

# Querschnitte für jede Variante
# WICHTIG: Das TL-Modell mit diskreten Impedanzsprüngen ist hier FALSCH.
# Der Spalt-zu-Kanal-Übergang (4 -> 950 mm²) reflektiert 98.8% der Energie
# für ALLE Varianten gleich. Der Puls schafft es nie als kohärenter Puls
# durch diesen Impedanzsprung.
#
# Die physikalisch korrekte Sichtweise:
# Die Bernoulli-Rückkopplung ist ein LOKALER Effekt am Spalt.
# Die Kammer wirkt als Druckreservoir / Compliance.
# Die Phase wird bestimmt durch die LOKALE Druck-Geschwindigkeits-Beziehung
# am Spalt, die von der Kammerimpedanz abhängt.
#
# Kammerimpedanz bei Frequenz f:
#   Z_Kammer(f) = dp/dQ = j*omega*rho/V_eff (Compliance bei niedrigen f)
#                       + Helmholtz-Resonanz bei f_H
#
# Die Phase der Kammerimpedanz bestimmt, ob Energie in die Zunge gepumpt wird.

# Kammer als akustischer Resonator mit Güte Q_H
Q_H_est = 7  # geschätzte Kammer-Güte (Holz, mit Klappe)

print(f"\n  KORREKTES MODELL: Kammerimpedanz und Phasenlage")
print(f"  (Transmission-Line-Modell NICHT anwendbar: Spalt reflektiert 98.8%)")
print(f"  -> Kammer wirkt als Druckreservoir/Compliance, nicht als Wellenleiter")

print(f"\n  Kammerimpedanz Z_Kammer(f) = dp/dQ:")
print(f"    Tiefe Frequenzen (f << f_H): Z ~ j*omega*rho*c²/V -> reine Compliance")
print(f"    Bei f_H = {f_H:.0f} Hz: Z wird real -> Resonanz, maximaler Energieaustausch")
print(f"    Hohe Frequenzen (f >> f_H): Z ~ j*omega*M_Hals -> Masse-dominiert")

# Phasenlage des Kammerdrucks relativ zur Spaltgeschwindigkeit
# Bei Helmholtz-Resonator:
#   dp_Kammer / v_Spalt = Z_Kammer(omega)
#   Phase von Z bestimmt, ob dp im Takt mit der Zungenbewegung ist
#
# Für Bernoulli-Verstärkung brauchen wir:
#   dp_Kammer maximal positiv wenn Zunge IN den Spalt bewegt (Spalt wird enger)
#   Das heißt: dp muss IN PHASE mit der Spaltverengung sein
#   -> dp muss ~90° VORAUS der Zungenauslenkung sein
#      (weil Geschwindigkeit = Ableitung der Auslenkung)

print(f"\n  Phase der Kammerimpedanz bei Zungenfrequenz und Obertönen:")
print(f"  Q_H (Kammer-Güte, geschätzt): {Q_H_est}")
print(f"  {'Oberton n':>10} {'f [Hz]':>8} {'f/f_H':>8} {'Phase Z [°]':>12} {'Wirkung auf Zunge':>25}")
print(f"  {'-'*68}")

phase_results = []
for n in [1, 2, 3, 5, 6, 8, 10, 15, 20, 25, 30]:
    f_n = n * 50
    # Phase des Helmholtz-Resonators
    # Z(omega) = (j*omega*M + R + 1/(j*omega*C))
    # Normiert: Z/Z_0 = 1 + j*Q*(f/f_H - f_H/f)
    # Phase = arctan(Q*(f/f_H - f_H/f))
    x = f_n / f_H
    phi_Z = np.degrees(np.arctan(Q_H_est * (x - 1/x)))
    
    # Für Bernoulli-Verstärkung:
    # Optimale Phase: ~0° (Druck IN PHASE mit Volumenstrom)
    # Bei phi > 0: Druck eilt VORAUS -> Compliance-dominiert -> schwache Kopplung
    # Bei phi < 0: Druck eilt NACH -> Masse-dominiert -> kann bremsend wirken
    # Bei phi = 0: Resonanz -> maximale Kopplung
    
    if abs(phi_Z) < 15:
        wirkung = 'STARKE Kopplung (Resonanz)'
    elif abs(phi_Z) < 45:
        wirkung = 'Gute Kopplung'
    elif abs(phi_Z) < 75:
        wirkung = 'Schwache Kopplung'
    else:
        wirkung = 'Entkoppelt (90° Phase)'
    
    phase_results.append((n, f_n, x, phi_Z, wirkung))
    print(f"  {n:>10} {f_n:>8} {x:>8.2f} {phi_Z:>12.1f} {wirkung:>25}")

print(f"\n  INTERPRETATION:")
print(f"    Die Kammerimpedanz hat bei f_H = {f_H:.0f} Hz Phase 0° (Resonanz).")
print(f"    Nur Obertöne NAHE f_H koppeln stark an die Kammer.")
print(f"    Der Grundton (50 Hz, f/f_H = {50/f_H:.2f}) ist weit unter der Resonanz")
print(f"    -> Phase ~-90° -> fast reine Compliance -> schwache Kopplung.")
print(f"    Das ist GUENSTIG: Die Kammer stört den Grundton nicht.")

print(f"\n  WAS DIE WANDFORM WIRKLICH BEEINFLUSST:")
print(f"    Nicht die Phase (die wird durch f_H und Q_H bestimmt, nicht durch die Wand),")
print(f"    sondern die EFFEKTIVE KAMMER-GUETE Q_H.")
print(f"    - Gerade Wand (A): Stagnationspunkt erzeugt Wirbel -> Energie wird")
print(f"      in Turbulenz dissipiert -> Q_H wird GESENKT -> Resonanzband BREITER")
print(f"      -> Mehr Obertöne werden beeinflusst (teils positiv, teils negativ)")
print(f"    - Parabolische Wand (C): Wenig Dissipation -> Q_H bleibt HOEHER")
print(f"      -> Resonanzband SCHMALER -> Weniger Obertöne betroffen")
print(f"      -> Aber die betroffenen werden STAERKER beeinflusst")
print(f"    Der Unterschied ist also: A beeinflusst viele OT schwach,")
print(f"    C beeinflusst wenige OT stark. Welches besser klingt: Geschmackssache.")
print(f"")
print(f"    Die Phase selbst wird primär durch die Kammergeometrie bestimmt")
print(f"    (Volumen, Halsquerschnitt, Halslänge -> f_H), NICHT durch die Wandform.")
print(f"    Die Wandform moduliert die Güte Q_H und damit die Bandbreite.")

print(f"\n  SCHLUSSFOLGERUNG (korrigiertes Modell):")
print(f"    1. Die Kammer ist kein Wellenleiter, sondern ein Helmholtz-Resonator.")
print(f"       Der Spalt reflektiert 98.8% -> TL-Modell nicht anwendbar.")
print(f"    2. Die PHASE des Energienachschubs wird durch f_H und Q_H bestimmt,")
print(f"       nicht durch die Wandform. Gleichmäßige Strömung hilft NICHTS")
print(f"       wenn die Phase nicht stimmt!")
print(f"    3. Die Wandform beeinflusst Q_H (effektive Kammer-Güte):")
print(f"       A senkt Q_H (breitere Kopplung, schwächer pro OT)")
print(f"       C erhält Q_H (schmalere Kopplung, stärker pro OT)")
print(f"    4. Der Instrumentenbauer stimmt f_H ab (Kammergröße),")
print(f"       nicht die Wandform. Die Wandform ist Feinabstimmung.")

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
    d.add(String(ox+20*sc, oy+kh+t_l+2*sc, 'Alu-Stimmplatte (keilförmig: 2->13 mm)',
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
        'Abb. 1: Längsschnitt. Keilförmige Platte (2->13 mm), Aufbiegung 1,5 mm am freien Ende.',
        fontSize=7, fillColor=black, fontName='Helvetica-Bold'))
    return d

def make_fig_varianten():
    w, h = 170*mm, 65*mm
    d = Drawing(w, h)
    labels = ['A: Gerade', 'B: Schräg', 'C: Parabel']
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
    phases = [('1: Spalt enger', -0.15, 'v steigt, p sinkt\n-> Sog auf Zunge'), ('2: Fast zu', -0.35, 'Strömung stoppt\n-> Feder zurück'), ('3: Über Ruhelage', 0.2, 'Spalt weit offen'), ('4: Rücklauf', 0.05, 'Zyklus neu')]
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
path = "/home/claude/bass_v8.pdf"
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

# Hilfsfunktion: Tabellenzellen mit HTML-Tags als Paragraph rendern
_tbl_sty = ParagraphStyle('TblCell', parent=sty['Normal'], fontSize=7, leading=9, alignment=TA_CENTER)
_tbl_sty_hdr = ParagraphStyle('TblHdr', parent=sty['Normal'], fontSize=7, leading=9, alignment=TA_CENTER, textColor=white, fontName='Helvetica-Bold')
def P(txt, hdr=False):
    """Wrap text in Paragraph if it contains HTML tags, else return as-is."""
    if isinstance(txt, str) and ('<sub>' in txt or '<super>' in txt or '<sup>' in txt or '<b>' in txt):
        return Paragraph(txt, _tbl_sty_hdr if hdr else _tbl_sty)
    return txt
def Prow(row, hdr=False):
    """Wrap all cells in a row."""
    return [P(c, hdr=hdr) for c in row]

st = []
st.append(Paragraph("Strömungsanalyse einer 50-Hz-Bass-Stimmzunge", sty['MT']))
st.append(Paragraph("v8 -- v7-Basis + Schlitz-Eintrittsverlust (K3b) + dynamische Spaltänderung (K6)", sty['Sub']))
st.append(HRFlowable(width="100%", thickness=2, color=HexColor('#e94560')))
st.append(Spacer(1, 3*mm))

# Kap 1
st.append(Paragraph("Kapitel 1: Aufbau", sty['Ch']))
st.append(Paragraph("Im Bassteil eines Akkordeons sitzen Dutzende Stimmzungen in eigenen Kammern, alle am selben Balg. Der Balg erzeugt den Druck, die Tasten steuern über mechanische Hebel die Klappen. Die Klappe wird <b>nicht</b> vom Balgdruck geöffnet - der Druck steht an allen geschlossenen Klappen an. Erst der Tastendruck öffnet die Klappe mechanisch.", sty['Key']))

# Kap 2
st.append(Paragraph("Kapitel 2: Kammer, Stimmplatte und Aufbiegung", sty['Ch']))
st.append(Paragraph(f"Die Kammer ist ein Holzquader von {L_kam*1e3:.0f} x {B_kam*1e3:.0f} x {H_kam*1e3:.0f} mm. Die gesamte Oberseite wird von einer Aluminium-Stimmplatte verschlossen. Diese Platte ist <b>keilförmig</b>: am freien Zungenende (über der Klappe) {t_platte_frei*1e3:.0f} mm dick, am Einspannende nur {t_platte_einsp*1e3:.0f} mm. Der Schlitz in der Platte ist daher am freien Ende {t_platte_frei*1e3:.0f} mm tief - die Luft durchquert einen kurzen Kanal, bevor sie die Zunge erreicht (Abbildung 1).", sty['B']))
st.append(Paragraph(f"Die Stahlzunge ({L_zunge*1e3:.0f} x {W_zunge*1e3:.0f} x {h_reed*1e3:.3f} mm) ist mit <b>zwei Schrauben</b> auf der Platte befestigt (kleinere Stimmplatten werden genietet). Das freie Ende ist um <b>{h_aufbiegung*1e3:.1f} mm aufgebogen</b> - es steht also im Ruhezustand {h_aufbiegung*1e3:.1f} mm über der Plattenoberfläche. Dadurch entsteht ein dreieckiger Spalt: 0 mm an der Einspannung, {h_aufbiegung*1e3:.1f} mm am freien Ende.", sty['B']))
st.append(Paragraph(f"Die effektive Spaltfläche berücksichtigt das Cantilever-Biegeprofil y(x) ~ y<sub>tip</sub> x (x/L)2. Der Spalt ist am Tip am größten, fällt aber quadratisch ab. Das Integral ergibt S<sub>Spalt</sub> = W x h<sub>max</sub>/3 = {W_zunge*1e3:.0f} x {h_aufbiegung*1e3:.1f}/3 = <b>{S_gap*1e6:.0f} mm<super>2</super></b> (v6 hatte Faktor 1/2 = {W_zunge*h_aufbiegung/2*1e6:.0f} mm2 -- lineares Profil). Zum Vergleich: Die Klappenöffnung hat {S_klap_eff*1e6:.0f} mm<super>2</super> - rund {S_klap_eff/S_gap:.0f}-mal größer.", sty['B']))
st.append(make_fig_schnitt_keil())
st.append(Spacer(1, 2*mm))

# Kap 3
st.append(PageBreak())
st.append(Paragraph("Kapitel 3: Der gefaltete Luftweg", sty['Ch']))
st.append(Paragraph(f"Eine Trennwand teilt die Kammer in zwei Kanäle. Der Luftweg: Klappe (Boden rechts) -> Kanal B nach links -> 180°-Faltung -> Kanal A nach rechts -> Schlitz -> Spalt. Gesamtweg über 200 mm.", sty['B']))

# Kap 4: Folie vs Furnier
st.append(Paragraph("Kapitel 4: Trennwand - Folie oder Furnier?", sty['Ch']))
st.append(Paragraph(f"Für die Trennwand kommen zwei Materialien in Frage: <b>Kunststofffolie</b> (~{t_folie*1e3:.1f} mm) oder <b>Holzfurnier</b> (~{t_furnier*1e3:.1f} mm). Die Wanddicke beeinflusst die Kanalbreite: Bei gerader Wand ergibt Folie {W_k_folie*1e3:.2f} mm pro Kanal, Furnier {W_k_furnier*1e3:.2f} mm - eine Differenz von {(W_k_folie-W_k_furnier)*1e3:.2f} mm ({(W_k_folie-W_k_furnier)/W_k_furnier*100:.1f}%). Strömungstechnisch ist das vernachlässigbar.", sty['B']))
st.append(Paragraph("Der eigentliche Unterschied liegt in drei Eigenschaften:", sty['B']))

d_mat = [['', 'Folie (0,1 mm)', 'Furnier (0,5 mm)'],
    ['Oberfläche', 'Glatt -> weniger Reibung\n(f ~ 0,02)', 'Rau -> mehr Reibung\n(f ~ 0,03)'],
    ['Stabilität', 'Flexibel -> kann flattern,\nbraucht Stützpunkte', 'Formstabil -> bleibt\nin Position'],
    ['Akustik', 'Kann als Membran schwingen\n-> Nebengeräusche', 'Akustisch inert'],
    ['Verarbeitung', 'Kleben, spannen', 'Kleben, klemmen, einfacher']]
t_mat = Table(d_mat, colWidths=[25*mm, 50*mm, 50*mm])
t_mat.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),('BACKGROUND',(0,1),(0,-1),HexColor('#e8e8e8')),('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,0),(0,-1),'Helvetica-Bold'),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_mat)
st.append(Paragraph(f"Da die Kanalreibung nur {z_reib:.2e} des Gesamtverlusts ausmacht, ist der Reibungsunterschied zwischen Folie und Furnier für den Durchfluss bedeutungslos. Der Unterschied liegt ausschließlich in Stabilität und Akustik. Furnier ist die sicherere Wahl; Folie nur sinnvoll, wenn der Platz extrem knapp ist.", sty['B']))

# Kap 5: Klappenorientierung
st.append(Paragraph("Kapitel 5: Klappenorientierung - Ausblas oder Umlenk?", sty['Ch']))
st.append(Paragraph(f"<b>Ausblas-Richtung:</b> Die Klappe öffnet so, dass der 30°-Strahl direkt in Kanal B zeigt. Der Strahl hat eine horizontale Komponente von {np.cos(alpha_r)*100:.0f}% und eine vertikale von {np.sin(alpha_r)*100:.0f}%. Verlustbeiwert: zeta<sub>Klappe</sub> = {zeta_klap_ausblas:.3f} (Vena contracta + Richtungsänderung).", sty['B']))
st.append(Paragraph(f"<b>Umlenk-Richtung:</b> Die Klappe öffnet nach außen. Die Luft muss ~150° um die Klappenkante kehrtmachen, bevor sie in die Kammer eintritt. Verlustbeiwert: zeta<sub>Klappe</sub> = {zeta_klap_umlenk:.3f} - das <b>{zeta_klap_umlenk/zeta_klap_ausblas:.1f}-fache</b> der Ausblas-Richtung.", sty['B']))

d_orient = [['', 'Ausblas', 'Umlenk'],
    ['Strahl', 'Gerichtet in Kammer', 'Erst weg, dann 150° Kurve'],
    ['zeta Klappe', f'{zeta_klap_ausblas:.3f}', f'{zeta_klap_umlenk:.3f}'],
    ['Wandform-Nutzung', 'Voll (Düse, Coanda)', 'Reduziert (breiter Eintritt)'],
    ['Nebengeräusche', 'Weniger', 'Mehr (Ablösung an Kante)'],
    ['Anwendung', 'Standardwahl', 'Nur wenn Bauform es erfordert']]
t_orient = Table(d_orient, colWidths=[32*mm, 45*mm, 45*mm])
t_orient.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),('BACKGROUND',(0,1),(0,-1),HexColor('#e8e8e8')),('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,0),(0,-1),'Helvetica-Bold'),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_orient)
st.append(Paragraph("Der höhere Verlust der Umlenk-Richtung wirkt sich auf die Spaltgeschwindigkeit kaum aus (Spalt dominiert), aber auf die Strömungsqualität: Der Strahl tritt breiter und langsamer in die Kammer ein, die Turbulenz am Eintritt ist höher, und die Wandform (Düse bei B, Coanda bei C) kann nicht voll wirken.", sty['B']))

# --- Analyse: Klappenöffnungswinkel ---
st.append(Paragraph("Klappenöffnung und Lautstärke", sty['Sec']))

angles_deg = [30, 25, 20, 15, 10, 7, 5, 3, 2, 1]
valve_rows = [Prow(['Winkel', 'S<sub>Klap,eff</sub> [mm2]', 'S<sub>Klap</sub>/S<sub>Spalt</sub>',
               'zeta<sub>ges</sub>', 'Q [ml/s]', 'Q/Q<sub>30°</sub>'], hdr=True)]
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
    "Praktische Tests zeigen: Bei Verringerung der Klappenöffnung treten "
    "Lautstärke-Einbußen auf. Die stationäre Berechnung zeigt, warum:", sty['B']))

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
    f"<b>Ergebnis:</b> Die stationäre Strömung ist erstaunlich unempfindlich gegen "
    f"die Klappenöffnung. Selbst bei 5° sind noch über 99% des maximalen Durchflusses "
    f"verfügbar. Erst unter ~{ang_thresh_5pct}° sinkt Q messbar (über 5%).", sty['B']))

st.append(Paragraph(
    "Das die praktisch beobachteten Lautstärke-Einbußen <b>nicht durch stationären "
    "Durchflussverlust</b> erklärbar sind (der Spalt dominiert zu stark), deutet auf "
    "einen <b>instationären Mechanismus</b> hin:", sty['Warn']))

st.append(Paragraph(
    "<b>1. Langsamerer Druckaufbau:</b> Eine kleinere Klappenöffnung verzögert den "
    "Druckanstieg in der Kammer. Die Zunge braucht länger, um die volle Amplitude "
    "zu erreichen. Bei schnellen Passagen wird die Maximalamplitude nie erreicht - "
    "das klingt leiser.", sty['B']))

st.append(Paragraph(
    "<b>2. Veränderte Jet-Kopplung:</b> Bei kleinerem Winkel ist der Strahl flacher "
    "und schneller (gleicher Volumenstrom durch kleinere Öffnung). Die Jet-Richtung "
    "ändert sich: Bei 30° hat der Strahl 50% vertikale Komponente, bei 10° nur noch "
    "17%. Weniger Vertikalkomponente bedeutet, dass der Strahl weniger direkt in den "
    "Kanal B gerichtet wird - ähnlich wie beim Übergang von Ausblas zu Umlenk.", sty['B']))

st.append(Paragraph(
    "<b>3. Akustische Impedanzänderung:</b> Die Klappenöffnung ist der akustische "
    "Eingang der Kammer. Eine kleinere Öffnung erhöht die akustische Impedanz am "
    "Kammereingang, was die Resonanzeigenschaften der Kammer verändert. Bei bestimmten "
    "Öffnungswinkeln können destruktive Interferenzen entstehen, die die "
    "Zungenamplitude reduzieren.", sty['B']))

st.append(Paragraph(
    "Alle drei Effekte wirken zusammen und erklären, warum die Lautstärke-Einbußen "
    "in der Praxis deutlicher ausfallen als es die stationäre Berechnung erwarten "
    "lässt. Auch hier gilt: <b>Die Berechnung zeigt, dass der Effekt nicht im "
    "Durchfluss liegt - also muss er im instationären Verhalten liegen.</b>", sty['Key']))

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
fold_rows = [Prow(['Radius R', 'R/D<sub>h</sub>', 'zeta<sub>Faltung</sub>',
              'zeta auf S<sub>Spalt</sub> bezogen', 'Anteil an zeta<sub>ges</sub>'], hdr=True)]

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
    f"D<sub>h</sub> = {D_h_fold*1e3:.1f} mm (Spalt {gap_fold*1e3:.0f} mm x Höhe "
    f"{H_kam*1e3:.0f} mm).", sty['B']))

st.append(t_fold_r)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    f"<b>Stationäres Ergebnis:</b> Die Faltung trägt bei <b>allen Radien weniger als "
    f"0,01% zum Gesamtwiderstand</b> bei. Der Grund: Die Faltungsquerschnittsfläche "
    f"S<sub>Faltung</sub> = {S_fold*1e6:.0f} mm2 ist {S_fold/S_gap:.0f}-mal größer als "
    f"der Spalt ({S_gap*1e6:.0f} mm2). Der Verlustbeiwert wird mit (S<sub>Spalt</sub>/"
    f"S<sub>Faltung</sub>)2 = ({S_gap/S_fold:.4f})2 = {(S_gap/S_fold)**2:.2e} skaliert - "
    f"praktisch Null.", sty['B']))

st.append(Paragraph(
    "<b>Instationäres Bild:</b> Für den Druckstoß-Durchgang ist die Frage "
    "differenzierter als beim stationären Durchfluss. Ein Viertelkreis-Profil "
    "mit R = 8-10 mm hätte zwei Effekte:", sty['B']))

st.append(Paragraph(
    "<b>a) Reflexionen reduzieren:</b> An einer scharfen 90°-Ecke wird ein "
    "Druckpuls teilweise reflektiert (akustische Impedanzsprung). Je steiler "
    "die Geometrieänderung, desto stärker die Reflexion. Ein Viertelkreis "
    "glättet den Impedanzübergang und lässt mehr Pulsenergie passieren.", sty['B']))

st.append(Paragraph(
    "<b>b) Wirbelbildung reduzieren:</b> Bei scharfen Ecken löst die Strömung "
    "bei jedem Druckstoß sofort ab und bildet Wirbel, die Energie dissipieren. "
    "Mit Radius R > 5 mm bleibt die Strömung länger angelegt.", sty['B']))

st.append(Paragraph(
    "<b>Aber Vorsicht:</b> Eine <b>zu kohärente</b> Strömung am Spalt ist "
    "nicht wünschenswert! Die Bernoulli-Selbsterregung braucht eine "
    "Anfangsstörung, um die Schwingung auszulösen. Eine perfekt gleichmäßige, "
    "laminare Strömung würde die Zunge nur statisch auslenken, ohne sie in "
    "Schwingung zu versetzen. Ein gewisses Maß an Turbulenz und räumlicher "
    "Ungleichmäßigkeit ist notwendig:", sty['Warn']))

st.append(Paragraph(
    "<b>Lokale Druckschwankungen</b> an der Zungenspitze stößen die erste "
    "Auslenkung an. <b>Wirbelablösungen</b> an den Spaltkanten liefern "
    "breitbandige Störungen, die die Zunge zum Schwingen anregen. "
    "Wäre die Strömung perfekt laminar und gleichförmig, müsste die Zunge "
    "auf eine infinitesimale Instabilität warten - der Einschwingvorgang "
    "würde sich deutlich verzögern.", sty['B']))

st.append(Paragraph(
    "Das Optimum liegt also nicht bei maximaler Kohärenz, sondern bei "
    "<b>schnellem Druckaufbau mit ausreichend Störungen</b>. Konkret: "
    "Der Druckstoß soll schnell am Spalt ankommen (-> wenig Reflexion in "
    "der Faltung, gute Wandführung), aber dort soll genügend Turbulenz "
    "entstehen, um die Selbsterregung sofort auszulösen. Die Turbulenz "
    "entsteht ohnehin am Spalt selbst (Re ~ 1400, Ablösungen an den "
    "Spaltkanten) - sie muss nicht zusätzlich aus der Kammer angeliefert "
    "werden.", sty['B']))

st.append(Paragraph(
    "<b>Empfehlung:</b> Ein Viertelkreis-Profil R = 8-10 mm an den Ecken der "
    "Faltung reduziert Reflexionen und Energieverluste im Druckstoß-Transport, "
    "ohne die Turbulenz am Spalt zu beeinträchtigen (der Spalt erzeugt seine "
    "eigene Turbulenz). Ob der Effekt hörbar ist, "
    "hängt davon ab, wie stark die Faltungsreflexionen den Anlauf tatsächlich "
    "bremsen - das lässt sich nur durch Versuch klären.", sty['Key']))

# Kap 5b: Düsenwinkel-Experiment
st.append(PageBreak())
st.append(Paragraph("Kapitel 5b: Düsenwinkel-Experiment - warum der Anströmwinkel die Amplitude bestimmt", sty['Ch']))

st.append(Paragraph(
    "Ein aufschlussreicher Praxistest: Stimmplatte montiert, Klappe entfernt, "
    "Anblasen mit einer Druckluftdüse aus verschiedenen Winkeln. Ergebnis: "
    "Die Zungenamplitude ändert sich <b>von Null bis Maximum</b> allein durch "
    "Änderung des Anströmwinkels. Wie lässt sich das erklären?", sty['B']))

st.append(Paragraph("Die zwei Komponenten des Jets", sty['Sec']))

st.append(Paragraph(
    "Ein Luftstrahl, der unter dem Winkel alpha auf den Spalt trifft, hat zwei "
    "Komponenten:", sty['B']))

nozzle_rows = [['Winkel alpha', 'Normalkomp. sin(alpha)', 'Tangentialkomp. cos(alpha)', 'Charakter']]
for ang in [0, 10, 20, 30, 45, 60, 75, 90]:
    ar = np.radians(ang)
    char = ('Rein tangential - kein Eintritt' if ang == 0
            else 'Rein normal - statischer Druck' if ang == 90
            else 'Flach - viel Überströmung' if ang <= 15
            else 'Optimal-Bereich?' if 20 <= ang <= 45
            else 'Steil - wenig Überströmung')
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
    "den Spalt in den Schlitz hinein. Erzeugt statischen Überdruck unter der Zunge. "
    "Dieser Druck biegt die Zunge nach oben aus - aber er schwingt nicht von selbst. "
    "Ohne weitere Kopplung hält er die Zunge nur statisch in einer Gleichgewichtslage.", sty['B']))

st.append(Paragraph(
    "<b>Tangentialkomponente</b> (parallel zur Platte, v*cos(alpha)): Strömt über "
    "die Spalt<b>öffnung</b> hinweg. Das ist der entscheidende Mechanismus:", sty['B']))

st.append(Paragraph("Der Flöten-Mechanismus am Zungenspalt", sty['Sec']))

st.append(Paragraph(
    "Wenn Luft über eine Öffnung streicht, entsteht an der Öffnung ein Unterdruck "
    "(Bernoulli). Genau das passiert bei der Tangentialkomponente des Jets: "
    "Sie strömt über den Spalt und erzeugt einen Sog, der Luft aus dem Schlitz "
    "herauszieht - also die Zunge nach <b>unten</b> zieht. "
    "Gleichzeitig drückt die Normalkomponente die Zunge nach <b>oben</b>.", sty['B']))

st.append(Paragraph(
    "Entscheidend ist: Beide Kräft ändern sich, wenn die Zunge schwingt. "
    "Bewegt sich die Zunge nach unten (Spalt wird größer), kann mehr Luft "
    "tangential über den Spalt strömen -> stärkerer Sog -> Zunge wird weiter "
    "nach unten gezogen. Bewegt sie sich nach oben (Spalt wird enger), verengt "
    "sich der Kanal für die Tangentialströmung -> weniger Sog -> Zunge federt "
    "zurück. Das ist eine <b>positive Rückkopplung</b>, die die Schwingung antreibt.", sty['B']))

st.append(Paragraph(
    "Dieser Mechanismus ist physikalisch identisch mit der Tonerzeugung bei "
    "Flöteninstrumenten: Dort strömt ein flacher Luftstrahl über eine scharfe "
    "Kante (Labium), und die Wechselwirkung zwischen Jet und Resonator erzeugt "
    "den Ton. Die Zungenkante übernimmt hier die Rolle des Labiums.", sty['B']))

st.append(Paragraph("Warum der Winkel die Amplitude bestimmt", sty['Sec']))

st.append(Paragraph(
    "<b>alpha = 0° (rein tangential):</b> Maximale Überströmung, aber kein Druck "
    "im Schlitz. Die Zunge hat keine Vorauslenkung - der Bernoulli-Sog muss gegen "
    "die volle Federsteifigkeit arbeiten. Außerdem fehlt der Durchfluss durch den "
    "Spalt, der für die klassische Bernoulli-Selbsterregung nötig ist. "
    "Ergebnis: Keine oder minimale Schwingung.", sty['B']))

st.append(Paragraph(
    "<b>alpha = 90° (rein normal):</b> Maximaler statischer Druck im Schlitz, "
    "aber keine Tangentialströmung. Die Zunge wird nach oben gedrückt und bleibt "
    "dort. Es fehlt der Bernoulli-Unterdruck an der Spaltoberseite, der die "
    "Zunge zurückziehen würde. Es fehlt auch die Wechselwirkung zwischen "
    "Spaltweite und Saugkraft. Ergebnis: Statische Auslenkung, keine Schwingung.", sty['B']))

st.append(Paragraph(
    "<b>Optimaler Winkel (~ 20-40°):</b> Beide Komponenten sind in der richtigen "
    "Balance. Die Normalkomponente liefert den Grunddruck, der die Zunge "
    "in den Arbeitsbereich vorspannt. Die Tangentialkomponente liefert den "
    "geschwindigkeitsabhängigen Bernoulli-Sog, der die Schwingung antreibt. "
    "Das Verhältnis entscheidet über die Phasenlage zwischen Druck und "
    "Zungenposition - und damit über die Energiezufuhr pro Zyklus.", sty['B']))

st.append(Paragraph(
    "Die Amplitude ändert sich stufenlos mit dem Winkel, weil die "
    "Energiezufuhr pro Schwingungszyklus vom Verhältnis der beiden Komponenten "
    "abhängt. Am optimalen Winkel ist die Netto-Energiezufuhr maximal - die "
    "Phasenbeziehung zwischen Druckschwankung und Zungengeschwindigkeit ist "
    "so, dass in jedem Zyklus maximale Energie in die Schwingung gepumpt wird.", sty['B']))

st.append(Paragraph("Konsequenz für die Kammerauslegung", sty['Sec']))

st.append(Paragraph(
    "Dieses Experiment erklärt unmittelbar, warum die Klappenorientierung und die "
    "Trennwandform so wichtig sind: Sie bestimmen den <b>effektiven Anströmwinkel</b> "
    "am Spalt.", sty['B']))

st.append(Paragraph(
    "Im normalen Betrieb (mit Klappe und Kammer) wird der Anströmwinkel nicht "
    "direkt vom Balgdruck bestimmt, sondern von der Geometrie des Luftwegs. "
    "Die Luft kommt aus Kanal A (parallel zur Platte) und muss um 90° in den "
    "Schlitz einbiegen. Der effektive Anströmwinkel hängt davon ab, wie die "
    "Strömung in Kanal A verteilt ist:", sty['B']))

st.append(Paragraph(
    "<b>Gerichtete Kanalströmung</b> (Wandführung erhält Strömungsrichtung): "
    "Die Luft biegt mit definierter Richtung in den Spalt ein. Es entsteht eine "
    "Tangentialkomponente, die den Bernoulli-Sog zuverlässig erzeugt. "
    "<b>Richtungslose, großskalig verwirbelte Kanalströmung</b> (viele "
    "Reflexionen, Stagnationspunkte): Die Richtungsinformation ist verloren. "
    "Die Luft nähert sich dem Spalt aus allen Richtungen - die "
    "Tangentialkomponente mittelt sich teilweise heraus. Der "
    "Bernoulli-Antrieb wird schwächer.", sty['B']))

st.append(Paragraph(
    "Das Düsenexperiment zeigt also <b>direkt</b>, warum die Trennwand das "
    "Einschwingverhalten beeinflusst: Nicht über den Durchfluss (der ist immer "
    "gleich), sondern über den <b>Anströmwinkel</b> am Spalt. Die Wandform "
    "bestimmt, wie viel von der ursprünglichen Richtungsinformation des "
    "Druckstoßes am Spalt noch ankommt. Die Turbulenz, die für die "
    "Selbsterregung nötig ist, entsteht am Spalt selbst (Ablösungen an den "
    "Spaltkanten bei Re ~ 1400, vgl. Kap. 5b) - sie muss nicht aus der Kammer "
    "mitgeliefert werden. Variante C (Parabel) liefert die beste Ansprache, weil "
    "sie die Richtungsinformation am besten erhält - nicht weil sie die "
    "<b>sauberste</b> Strömung liefert.", sty['Key']))

# Kap 6: Trennwandvarianten
st.append(Paragraph("Kapitel 6: Drei Trennwandvarianten", sty['Ch']))
st.append(make_fig_varianten())
st.append(Spacer(1, 2*mm))
st.append(Paragraph(f"<b>A - Gerade:</b> Beide Kanäle {W_k*1e3:.1f} mm. 30°-Strahl prallt ~90° auf Wand. Stagnationspunkt, Ablösung, stärkster Richtungsverlust und Reflexion.", sty['B']))
st.append(Paragraph(f"<b>B - Schräg (beta={beta_B:.1f}°):</b> Kanal B bei Klappe {W_B_unten_klap*1e3:.1f} mm (Düse), bei Faltung {W_B_unten_falt*1e3:.1f} mm. Kanal A bei Spalt {W_B_oben_klap*1e3:.1f} mm (weit). Diffusor-Halbwinkel nur {np.degrees(np.arctan((W_B_unten_falt-W_B_unten_klap)/(2*L_wall))):.1f}° -> sicher ablösungsfrei.", sty['B']))
st.append(Paragraph(f"<b>C - Parabel:</b> Gleiche Endpunkte wie B, aber stetige Krümmung. Coanda-Effekt führt Strahl entlang konvexer Wand. Beste Richtungserhaltung des Druckstoßes.", sty['B']))

# Ergebnistabelle
v1k = np.sqrt(2*1000/(rho*2))
Q1k = v1k * S_gap
d_res = [['(bei 1000 Pa, Furnier, Ausblas)', 'A: Gerade', 'B: Schräg', 'C: Parabel']]
for cname_list in [['A-Gerade + Ausblas', 'B-Schräg + Ausblas', 'C-Parabel + Ausblas']]:
    d_res.append(['zeta total'] + [f"{configs[c]['zeta']:.6f}" for c in cname_list])
    d_res.append(['v Spalt [m/s]'] + [f"{configs[c]['v1000']:.2f}" for c in cname_list])
    d_res.append(['Q [ml/s]'] + [f"{configs[c]['Q1000']*1e6:.1f}" for c in cname_list])
d_res.append(['Kanal B Klappe [mm]', f'{W_k*1e3:.1f}', f'{W_B_unten_klap*1e3:.1f}', f'{W_B_unten_klap*1e3:.1f}'])
d_res.append(['Kanal A Spalt [mm]', f'{W_k*1e3:.1f}', f'{W_B_oben_klap*1e3:.1f}', f'{W_B_oben_klap*1e3:.1f}'])
t_res = Table(d_res, colWidths=[40*mm, 28*mm, 28*mm, 28*mm])
t_res.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),('BACKGROUND',(0,1),(0,-1),HexColor('#e8e8e8')),('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(0,1),(0,-1),'Helvetica-Bold'),('ALIGN',(0,0),(-1,-1),'CENTER'),('VALIGN',(0,0),(-1,-1),'MIDDLE'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_res)

st.append(Paragraph(
    "<b>Achtung: Diese Tabelle zeigt nur den stationären Zustand.</b> "
    "Die nahezu identischen Volumenströme täuschen darüber hinweg, dass "
    "die Trennwand das Einschwingverhalten in der Praxis entscheidend "
    "beeinflusst. Praktische Messungen zeigen deutliche Unterschiede in der "
    "Ansprache zwischen verschiedenen Wandformen - obwohl der stationäre "
    "Durchfluss bei allen fast gleich ist.", sty['Warn']))

st.append(Paragraph("Warum die Trennwand trotzdem entscheidend ist:", sty['Sec']))

st.append(Paragraph(
    "Der Widerspruch löst sich auf, wenn man den <b>instationären Anlauf</b> "
    "betrachtet statt den Gleichgewichtszustand. In den ersten Millisekunden nach "
    "Klappenöffnung passiert folgendes: Ein Druckstoß läuft vom Ventil durch "
    "Kanal B, um die Faltung, durch Kanal A zum Spalt. In dieser initialen Phase "
    "wirkt die Trennwand als Strömungsführung für den Druckaufbau. "
    "Ihre Form bestimmt, wie schnell der Druck am Spalt ankommt. "
    "Im eingeschwungenen Zustand (nach wenigen ms) wirkt die Kammer dann als "
    "Helmholtz-Resonator, und die Phasenlage der Kammerimpedanz bestimmt die "
    "Oberton-Kopplung (siehe Kap. 10b).", sty['B']))

st.append(Paragraph(
    "<b>Gerade Wand (A):</b> Der 30-Grad-Strahl prallt auf die Wand, erzeugt einen "
    "Stagnationspunkt und Wirbel. Der Druckstoß wird breit gestreut und "
    "teilweise reflektiert. Am Spalt kommt er verspätet und mit reduzierter "
    "Amplitude an. Der Druckaufbau ist langsam - die Zunge erreicht den "
    "Schwellendruck für die Bernoulli-Rückkopplung später.", sty['B']))

st.append(Paragraph(
    "<b>Schräge Wand (B):</b> Die Düse am Eintritt beschleunigt den Druckstoß. "
    "Weniger Reflexionen auf dem Weg. Am Spalt kommt der Druck schneller auf "
    "den Schwellenwert. Der Anlauf ist kürzer.", sty['B']))

st.append(Paragraph(
    "<b>Parabolische Wand (C):</b> Die stetige Krümmung führt den Druckstoß "
    "entlang der Wand (Coanda-Effekt). Minimale Reflexionen, daher schnellster "
    "Druckaufbau am Spalt. Entscheidend ist nicht, dass die Strömung 'sauber' "
    "oder laminar ankommt - entscheidend ist, dass der Druckstoß seine "
    "Richtungsinformation behält und wenig Energie durch Reflexionen verliert. "
    "Die Turbulenz, die die Selbsterregung auslöst, entsteht am Spalt selbst "
    "(Ablösungen an den Kanten bei Re ~ 1400, vgl. Kap. 5b).", sty['B']))

st.append(Paragraph(
    "In der Elektrotechnik-Analogie (Kapitel 10): Die stationäre Berechnung "
    "entspricht dem Gleichstromwiderstand eines Kabels - für alle drei Varianten "
    "fast gleich. Aber die <b>Impulsantwort</b> (wie schnell ein Spannungssprung "
    "am Ende ankommt) hängt von der Wellenimpedanz, Reflexionen und Dispersion "
    "ab - und die sind völlig verschieden. Für die Ansprache zählt nicht der "
    "Gleichstromwiderstand, sondern die Impulsantwort.", sty['Key']))

# ===============================================================
# Kap 6b: Dynamische Spaltfläche
# ===============================================================
st.append(PageBreak())
st.append(Paragraph("Kapitel 6b: Die dynamische Spaltfläche - warum 6 mm2 nur die halbe Wahrheit ist", sty['Ch']))

st.append(Paragraph(
    "Alle bisherigen Berechnungen verwenden die Ruhelage-Spaltfläche "
    "S<sub>Spalt</sub> = " + f"{S_gap*1e6:.0f}" + " mm2. Aber die Zunge schwingt - und bei voller "
    "Lautstärke schwingt die Zungenspitze ca. 20 mm nach unten durch den "
    "Schlitz. Dabei gibt die Zunge den <b>gesamten Schlitzquerschnitt</b> frei.", sty['B']))

st.append(Paragraph("Geometrie des Durchschwingvorgangs", sty['Sec']))

st.append(Paragraph(
    "Der Schlitz hat die Breite W<sub>Schlitz</sub> = 9 mm und verläuft durch "
    "die gesamte keilförmige Platte (2-13 mm dick). Der maximale "
    "Schlitzquerschnitt (wenn die Zunge komplett entfernt wäre) beträgt "
    "W<sub>Schlitz</sub> x t<sub>mittel</sub> = 9 x 7,5 = "
    "<b>67,5 mm2</b> - das 11-fache der Ruhelage-Spaltfläche.", sty['B']))

# Berechne dynamische Spaltflächen
dyn_N = 500
dyn_x = np.linspace(0, L_zunge, dyn_N)
dyn_dx = dyn_x[1] - dyn_x[0]
side_gap_w = (W_schlitz - W_zunge) / 2  # 0.5mm

def cantilever_shape(xi, y_tip):
    return y_tip * (3*L_zunge*xi**2 - xi**3) / (2*L_zunge**3)

def local_plate_t(xi):
    return t_platte_einsp + (t_platte_frei - t_platte_einsp) * xi / L_zunge

def calc_S_eff(y_tip_val):
    S_sum = 0
    L_frei = 0
    for xi in dyn_x:
        y = cantilever_shape(xi, y_tip_val)
        tp = local_plate_t(xi)
        if y > 0:
            local_S = W_zunge * y
        elif y > -tp:
            local_S = 2 * side_gap_w * tp
        else:
            local_S = W_schlitz * tp
            L_frei += dyn_dx
        S_sum += local_S * dyn_dx
    return S_sum / L_zunge, L_frei

# Statische Tabelle: verschiedene Tip-Positionen
dyn_tips = [1.5, 0.5, 0, -2, -5, -10, -15, -20, -25, -30]
dyn_rows = [Prow(['Tip-Position', 'Zustand', 'S<sub>eff</sub> [mm2]', 
             'S/S<sub>Ruhe</sub>', 'Freigegebene Länge'], hdr=True)]

for yt_mm in dyn_tips:
    yt = yt_mm / 1000.0
    S_e, L_f = calc_S_eff(yt)
    S_mm2 = S_e * 1e6
    ratio = S_mm2 / (S_gap * 1e6)
    
    if yt_mm > 0: zust = "Über Platte"
    elif yt_mm > -2: zust = "Schließen/Eintauchen"
    elif yt_mm > -13: zust = "Teilw. durchgeschwungen"
    else: zust = "Großteils frei"
    
    dyn_rows.append([f'{yt_mm:+.0f} mm', zust, f'{S_mm2:.1f}', 
                     f'{ratio:.1f}x', f'{L_f*1e3:.0f} mm'])

t_dyn = Table(dyn_rows, colWidths=[22*mm, 36*mm, 24*mm, 20*mm, 30*mm])
t_dyn.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))

st.append(Paragraph(
    "Die Tabelle zeigt den effektiven Durchströmquerschnitt (Mittelwert über "
    "die Zungenlänge) für verschiedene Tip-Positionen. Drei Phasen:", sty['B']))

st.append(t_dyn)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "<b>Phase 1 - Zunge über Platte (+1,5 bis 0 mm):</b> "
    "Normaler Spaltbetrieb. S sinkt von 4,5 auf 0 mm2. "
    "In dieser Phase gilt die bisherige Berechnung.", sty['B']))

st.append(Paragraph(
    "<b>Phase 2 - Zunge im Schlitz (0 bis ~-13 mm):</b> "
    "Die Zunge taucht in den Schlitz ein. Der Hauptquerschnitt ist blockiert, "
    "nur die schmalen Seitenspalte (je 0,5 mm) bleiben offen. "
    "S ist auf ca. 7,5 mm2 begrenzt (Seitenspalt-Querschnitt). "
    "Diese Phase ist der <b>Engpass im Schwingzyklus</b>.", sty['B']))

st.append(Paragraph(
    "<b>Phase 3 - Zunge unter der Platte (unter -13 mm):</b> "
    "Die Zungenspitze hat den Schlitz verlassen. Von der Spitze her wird der "
    "volle Schlitzquerschnitt freigegeben (W<sub>Schlitz</sub> x t<sub>Platte</sub>). "
    "Bei -20 mm Tip: S<sub>eff</sub> = 42 mm2 (<b>7-fach</b>), "
    "28 mm der Zungenlänge sind frei. "
    "Bei -30 mm: S<sub>eff</sub> = 54 mm2 (<b>9-fach</b>).", sty['B']))

# Zeitlicher Verlauf
st.append(Paragraph("Zeitlicher Verlauf im Schwingzyklus", sty['Sec']))

amp_rows = [Prow(['Amplitude', 'Tip<sub>min</sub>', 'S<sub>min</sub> [mm2]',
             'S<sub>max</sub> [mm2]', 'S<sub>mittel</sub> [mm2]', 'Mittel/Ruhe'], hdr=True)]

for A_mm in [1, 2, 5, 10, 15, 20, 25]:
    A = A_mm / 1000.0
    N_t = 200
    S_list = []
    for phase in np.linspace(0, 2*np.pi, N_t):
        yt = h_aufbiegung - A + A * np.cos(phase)
        S_e, _ = calc_S_eff(yt)
        S_list.append(S_e * 1e6)
    
    tip_min = (h_aufbiegung - 2*A) * 1e3
    amp_rows.append([
        f'{A_mm} mm', f'{tip_min:+.1f} mm', f'{min(S_list):.1f}',
        f'{max(S_list):.1f}', f'{np.mean(S_list):.1f}',
        f'{np.mean(S_list)/(S_gap*1e6):.1f}x'])

t_amp = Table(amp_rows, colWidths=[18*mm, 20*mm, 22*mm, 22*mm, 24*mm, 20*mm])
t_amp.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))

st.append(t_amp)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "<b>Entscheidend:</b> Bei voller Lautstärke (A ~ 20 mm) oszilliert der "
    "Durchströmquerschnitt zwischen fast Null und 58 mm2 - ein Faktor von "
    "über <b>100</b> innerhalb einer einzigen Schwingungsperiode von 20 ms. "
    "Der zeitgemittelte Querschnitt liegt bei 32 mm2 - mehr als "
    "<b>5-mal</b> so gross wie in der Ruhelage.", sty['Warn']))

st.append(Paragraph("Konsequenzen für die Strömungsanalyse", sty['Sec']))

st.append(Paragraph(
    "<b>1. Der Spalt dominiert nur beim Einschwingen:</b> "
    "Zu Beginn der Schwingung (kleine Amplitude) ist der Spalt mit 6 mm2 "
    "tatsächlich der Engpass - die Kammer mit ~950 mm2 Kanalquerschnitt "
    "ist irrelevant. Aber sobald die volle Amplitude erreicht ist, ist "
    "der Spalt in jeder Halbperiode für einen Bruchteil der Zeit weit offen. "
    "In dieser Phase werden die Kammer, die Klappe (244 mm2) und die "
    "Trennwand zu den relativ engsten Stellen.", sty['B']))

st.append(Paragraph(
    "<b>2. Der Volumenstrom pulsiert massiv:</b> "
    "Während eines Zyklus wechselt der Durchfluss zwischen praktisch Null "
    "(Zunge schließt Spalt) und einem Maximum, das weit über dem stationären "
    "Wert liegt (weil S so viel größer ist). Der instationäre Volumenstrom "
    "ist kein sanftes Sinussignal - er hat extrem steile Flanken. "
    "Diese Pulse müssen durch die Kammer.", sty['B']))

st.append(Paragraph(
    "<b>3. Die Kammer wird zum Puffer:</b> "
    "Die Kammer mit ihrem Volumen von " + f"{L_kam*B_kam*H_kam*1e6:.0f}" + " cm3 "
    "muss die extremen Volumenstrom-Schwankungen abpuffern. Wenn der Spalt "
    "kurzzeitig 10x so viel Luft durchlässt wie im stationären Fall, muss "
    "diese Luft irgendwo herkommen - aus dem Kammervolumen. Der Druck in "
    "der Kammer fällt bei jedem Durchschwing-Puls kurzzeitig ab und wird "
    "durch die Klappe nachgespeist. Die Klappenöffnung (244 mm2) begrenzt "
    "die Nachspeisung.", sty['B']))

st.append(Paragraph(
    "<b>4. Die Klappenfläche wird bei voller Lautstärke relevant:</b> "
    "Im stationären Bild ist die Klappe mit S<sub>Klap</sub>/S<sub>Spalt</sub> = "
    f"{S_klap_eff/S_gap:.0f}x völlig überdimensioniert. Aber wenn der Spalt "
    "kurzzeitig 58 mm2 hat, ist das Verhältnis nur noch "
    f"S<sub>Klap</sub>/S<sub>Spalt,dyn</sub> = {S_klap_eff*1e6/58:.1f}x. "
    "Damit wird die Klappenöffnung relevant - und erklärt direkt, warum "
    "eine kleinere Klappenöffnung bei voller Lautstärke die Amplitude "
    "begrenzt (Kapitel 5).", sty['Key']))

st.append(Paragraph(
    "<b>5. Die Trennwand sieht pulsierenden Gegendruck:</b> "
    "Bei jeder Durchschwing-Phase entsteht ein Unterdruck-Puls in der Kammer "
    "(Luft wird durch den offenen Spalt abgesaugt). Dieser Puls läuft "
    "durch die Kammer zurück zur Klappe. Wie die Kammer diesen Puls "
    "transportiert (Reflexionen, Dämpfung, Dispersion) hängt von der "
    "Trennwandform ab - genau wie der Druckstoß beim Einschwingen. "
    "Die Trennwandform beeinflusst also nicht nur das Einschwingen, "
    "sondern auch den <b>eingeschwungenen Betrieb bei voller Lautstärke</b>.", sty['B']))

st.append(Paragraph(
    "<b>Zusammenfassung:</b> Die stationäre Berechnung mit S<sub>Spalt</sub> = "
    f"{S_gap*1e6:.0f} mm2 beschreibt zuverlässig das Einschwingen (kleine "
    "Amplitude) und gibt die richtigen Relationen zwischen den Varianten. "
    "Bei voller Lautstärke ändert sich das Bild fundamental: Der Spalt "
    "dominiert nicht mehr allein, die Kammer und die Klappe werden zu "
    "Mitspielerern, und der Volumenstrom wird extrem pulsierend. "
    "In dieser Phase wird die Kammergeometrie auch für den <b>stationären</b> "
    "(zeitgemittelten) Widerstand relevant - nicht nur für die Impulsantwort.", sty['Key']))

# Kap 7: Bernoulli + Einschwingzeit
st.append(PageBreak())
st.append(Paragraph("Kapitel 7: Bernoulli-Selbsterregung und Einschwingzeit", sty['Ch']))
st.append(Paragraph("Die Selbsterregung erfolgt über den Bernoulli-Mechanismus (Abbildung 3): Verengt sich der Spalt, steigt v, sinkt p, Unterdruck verstärkt Auslenkung. Am Umkehrpunkt bricht die Strömung zusammen, Feder treibt Zunge zurück.", sty['B']))
st.append(make_fig_bernoulli())
st.append(Spacer(1, 2*mm))
st.append(Paragraph(f"Schwellendruck: DeltaP<sub>min</sub> = {dp_min:.1f} Pa (mit {h_aufbiegung*1e3:.1f} mm Aufbiegung). Niedrig gegenüber 500-5000 Pa Balgdruck.", sty['B']))
st.append(Paragraph("Die Einschwingzeit tau = Q/(pi*f) ist der eigentliche Engpass:", sty['B']))
rows = [['zeta', 'Q', 'tau [ms]', 'Perioden']]
for z, Q, tau in zip(zeta_vals, Q_vals, tau_vals):
    rows.append([f'{z}', f'{Q:.0f}', f'{tau*1e3:.0f}', f'{Q:.0f} x 20 ms'])
t_tau = Table(rows, colWidths=[22*mm, 22*mm, 22*mm, 30*mm])
t_tau.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),('FONTSIZE',(0,0),(-1,-1),8),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),('FONTNAME',(1,1),(-1,-1),'Courier'),('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_tau)
st.append(Paragraph("Q hängt von Material, Geometrie und Einspannung ab und muss gemessen werden.", sty['B']))

# ═══════════════════════════════════════════════════════════════
# Kap 7b: [K6] Statische Spaltverengung unter Blasdruck
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 7b: Statische Spaltverengung unter Blasdruck (v8 neu)", sty['Ch']))

st.append(Paragraph(
    "Die bisherige Berechnung verwendet die Ruhelage-Aufbiegung h<sub>max</sub> = "
    f"{h_aufbiegung*1e3:.1f} mm als festen Wert. Aber <b>bevor</b> die Zunge zu schwingen "
    "beginnt, liegt bereits Balgdruck an. Dieser Druck wirkt als Gleichlast auf den "
    "Cantilever und biegt die Zungenspitze nach unten - der Spalt wird enger.", sty['B']))

st.append(Paragraph(
    "Cantilever unter Gleichlast q = DeltaP x W<sub>Zunge</sub>: "
    "w<sub>tip</sub> = q x L<super>4</super> / (8 x E x I) = "
    "DeltaP x W x L<super>4</super> / (8EI). "
    "Die effektive Aufbiegung wird h<sub>eff</sub> = h<sub>max</sub> - w<sub>tip</sub>(DeltaP).", sty['B']))

# Druckscan-Tabelle im PDF
k6_rows = [Prow(['DeltaP [Pa]', 'DeltaP [mbar]', 'w<sub>tip</sub> [mm]',
            'h<sub>eff</sub> [mm]', 'S<sub>Spalt</sub> [mm2]', 'S/S<sub>Ruhe</sub>'], hdr=True)]
for dp_val in [100, 200, 500, 1000, 2000, 3000, int(round(dp_min))]:
    w = w_static_tip(dp_val)
    h_e = h_eff_static(dp_val)
    S_e = S_gap_eff(dp_val)
    ratio = S_e / S_gap
    k6_rows.append([f'{dp_val}', f'{dp_val/100:.1f}', f'{w*1e3:.4f}',
                    f'{h_e*1e3:.4f}', f'{S_e*1e6:.2f}', f'{ratio:.3f}'])

t_k6 = Table(k6_rows, colWidths=[22*mm, 22*mm, 22*mm, 22*mm, 24*mm, 20*mm])
t_k6.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_k6)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    f"<b>Grenzdruck:</b> Der Spalt schließt sich bei "
    f"DeltaP<sub>crit</sub> = 8EI x h<sub>max</sub> / (W x L<super>4</super>) = "
    f"<b>{dp_crit:.0f} Pa ({dp_crit/100:.1f} mbar)</b>. Das ist identisch mit dem "
    f"Schwellendruck DeltaP<sub>min</sub> = {dp_min:.0f} Pa - und das ist kein "
    f"Zufall: Beide Formeln berechnen denselben physikalischen Effekt (Gleichlast, "
    f"die die Aufbiegung kompensiert).", sty['Key']))

st.append(Paragraph(
    "Das bedeutet: Bei genau dem Druck, der die Schwingung auslösen soll, ist der "
    "Spalt bereits fast geschlossen. In der Praxis beginnt die Schwingung aber "
    "deutlich früher, weil die Bernoulli-Instabilität nicht erst bei voller "
    "statischer Schließung einsetzt:", sty['B']))

st.append(Paragraph(
    "<b>1. Dynamischer Bernoulli-Effekt:</b> Schon bei ca. 30-50% von DeltaP<sub>min</sub> "
    "erzeugt die Strömung durch den enger werdenden Spalt genügend Bernoulli-Sog, "
    "um die erste Schwingungsauslenkung auszulösen. Der Spalt muss nicht Null werden - "
    "es genügt, dass die Druck-Geschwindigkeits-Kopplung instabil wird.", sty['B']))

st.append(Paragraph(
    "<b>2. Nichtlineare Verstärkung:</b> Die Spaltverengung ist selbstverstärkend: "
    "Weniger Spalt -> höhere Geschwindigkeit -> stärkerer Bernoulli-Sog -> "
    "weitere Verengung. Dieser positive Rückkopplungskreis macht das System "
    "schon vor der statischen Schließung instabil.", sty['B']))

st.append(Paragraph(
    f"<b>3. Geschätzter Einsetzpunkt:</b> Erfahrungswerte und analytische Modelle "
    f"(Fletcher, Cottingham) zeigen, dass durchschlagende Zungen typisch bei "
    f"30-50% des statischen Schwellendrucks zu schwingen beginnen. Für diese Zunge: "
    f"DeltaP<sub>dyn</sub> ~ {dp_dynamic_est:.0f} Pa ({dp_dynamic_est/100:.0f} mbar). "
    f"Bei diesem Druck ist der Spalt noch auf "
    f"{h_eff_static(dp_dynamic_est)*1e3:.2f} mm geöffnet "
    f"({h_eff_static(dp_dynamic_est)/h_aufbiegung*100:.0f}% der Ruhelage).", sty['B']))

st.append(Paragraph(
    "<b>Konsequenz für die Praxis:</b> Die statische Spaltverengung erklärt, warum "
    "der effektive Balgdruck-Bereich einer Zunge begrenzt ist. Unterhalb von "
    f"~{dp_dynamic_est:.0f} Pa spricht die Zunge nicht an (zu wenig Anregung). "
    f"Oberhalb von ~{dp_crit:.0f} Pa schließt der Spalt vollständig, die Zunge "
    "kann nicht mehr frei schwingen und wird in den Schlitz gedrückt. "
    "Dazwischen liegt der <b>Arbeitsbereich</b>, in dem Amplitude und Lautstärke "
    "mit dem Druck steigen. Dieser Bereich kann durch Änderung der Aufbiegung "
    "verschoben werden: Mehr Aufbiegung -> höherer dp_crit -> breiterer Bereich, "
    "aber auch mehr Luftverbrauch im Leerlauf.", sty['Warn']))

# ═══════════════════════════════════════════════════════════════
# Kap 8: Dämpfung - was sie beeinflusst und wie
# ═══════════════════════════════════════════════════════════════
st.append(PageBreak())
st.append(Paragraph("Kapitel 8: Dämpfung - was sie beeinflusst und wie", sty['Ch']))

st.append(Paragraph(
    "Die Güte Q bestimmt die Einschwingzeit und damit die Ansprache. Q ist keine "
    "Materialkonstante, sondern die Summe aller Energieverluste pro Schwingungszyklus. "
    "Diese lassen sich in vier physikalisch unterscheidbare Mechanismen zerlegen, "
    "deren Beiträge sich addieren:", sty['B']))

st.append(Paragraph(
    "<b>1. Materialdämpfung (innere Reibung):</b> "
    "Jedes Metall wandelt bei Verformung einen Bruchteil der elastischen Energie in "
    "Wärme um. Der Verlustfaktor eta ist eine Materialeigenschaft. Für Federstahl "
    "liegt eta bei 0,0002-0,001, für Messing bei 0,001-0,003 - also bis zu zehnmal "
    "höher. In der Güte-Zerlegung: Q<sub>Material</sub> = 1/eta. Federstahl kommt "
    "allein auf Q<sub>Material</sub> ~ 1000-5000, Messing auf 300-1000. Dieser "
    "Beitrag ist frequenzunabhängig.", sty['B']))

st.append(Paragraph(
    "<b>2. Luftdämpfung (viskose Verluste im Spalt):</b> "
    "Die schwingende Zunge verdrängt Luft durch den engen Spalt. Dabei entsteht ein "
    "Squeeze-Film-Effekt: Die Luft wird bei jedem Zyklus durch die Engstelle gepresst "
    "und wieder zurückgesaugt. Die dissipierte Leistung skaliert mit "
    "P<sub>visc</sub> ~ mu*W*L<super>3</super>*omega<super>2</super>*a<super>2</super>/h<super>3</super>, "
    "wobei h die Spalthöhe ist. Mit " + f"{h_aufbiegung*1e3:.1f}" + " mm Aufbiegung ist h deutlich "
    "größer als bei einem 0,4-mm-Spalt - der Squeeze-Film-Beitrag sinkt mit "
    "h<super>-3</super>, also um Faktor (1,5/0,4)<super>3</super> ~ 53. "
    "Bei der großen Aufbiegung einer Basszunge ist die Luftdämpfung daher "
    "ein untergeordneter Beitrag.", sty['B']))

st.append(Paragraph(
    "<b>3. Einspannungsdämpfung:</b> "
    "An der Verschraubung (zwei Schrauben bei Basszungen) wird bei jeder Schwingung "
    "mikroskopisch Energie in Reibung umgewandelt. Die Güte der Einspannung hängt "
    "vom Anpressdruck, der Oberflächengüte und dem Plattenmaterial ab. Aluminium "
    "ist weicher als Stahl - bei jeder Schwingung verformt sich die Kontaktzone "
    "minimal. Dieser Beitrag ist schwer zu berechnen, aber erfahrungsgemäß einer "
    "der größten Einzelposten. Eine sauber gefräste Auflagefläche und "
    "gleichmäßig angezogene Schrauben minimieren ihn.", sty['B']))

st.append(Paragraph(
    "<b>4. Schallabstrahlung:</b> "
    "Die schwingende Zunge strahlt Schall ab - das ist ja der Zweck. Die abgestrahlte "
    "Leistung geht dem Schwinger als Dämpfung verloren. Bei 50 Hz und einer Fläche "
    "von 70 x 8 mm ist der Strahlungswiderstand gering (die Zunge ist viel kleiner als "
    "die Wellenlänge von 6,9 m). Dieser Beitrag ist bei tiefen Basszungen klein, "
    "wächst aber mit f<super>2</super> und wird bei hohen Diskant-Frequenzen relevant.", sty['B']))

st.append(Paragraph(
    "Die Gesamt-Güte ergibt sich aus den Einzelbeiträgen als harmonische Summe: "
    "1/Q<sub>ges</sub> = 1/Q<sub>Mat</sub> + 1/Q<sub>Luft</sub> + 1/Q<sub>Einsp</sub> + "
    "1/Q<sub>Strahl</sub>. Der kleinste Einzelwert dominiert. Bei einer typischen "
    "Basszunge aus Federstahl mit Schraubbefestigung und 1,5 mm Aufbiegung dürfte "
    "Q<sub>ges</sub> zwischen 50 und 150 liegen, wobei die Einspannung und "
    "Materialdämpfung die größten Beiträge liefern.", sty['Key']))

st.append(Paragraph(
    "Was folgt daraus für die Praxis? Q lässt sich <b>senken</b> (schnellere Ansprache) "
    "durch: Materialwahl (Messing statt Stahl), Massebeladung am Tip mit weichem Metall "
    "(erhöhte Materialdämpfung des Gewichts), geringere Aufbiegung (mehr Squeeze-Film), "
    "oder weichere Einspannung. Q lässt sich <b>erhöhen</b> (längeres Nachklingen, "
    "lauterer Ton) durch: härteren Stahl, präzisere Einspannung, größere Aufbiegung. "
    "Der Instrumentenbauer navigiert zwischen diesen Polen - Messung (z.B. Ausschwing-"
    "versuch, Logarithmisches Dekrement) ist der einzige zuverlässige Weg, Q einer "
    "konkreten Zunge zu bestimmen.", sty['B']))

# ═══════════════════════════════════════════════════════════════
# Kap 9: Propeller, Windräder und Hinterkanten-Optimierung
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 9: Lehren aus der Propeller- und Windkraft-Aerodynamik", sty['Ch']))

st.append(Paragraph(
    "In der Windkraft- und Propellertechnik sind in den letzten Jahrzehnten erhebliche "
    "Fortschritte bei der Optimierung von Strömungskanten erzielt worden. Viele dieser "
    "Erkenntnisse betreffen dieselben physikalischen Phänomene, die auch in der "
    "Stimmzungenkammer auftreten.", sty['B']))

st.append(Paragraph("<b>Hinterkanten-Optimierung (Trailing Edge):</b> "
    "Windkraftrotoren verwenden gezackte oder gekämmte Hinterkanten (Serrations), um "
    "tonale Geräusche zu reduzieren. Der Mechanismus: An einer scharfen Kante lösen "
    "sich Wirbel kohärenter in einer einzigen Frequenz ab (Karman-Wirbel). Gezackte "
    "Kanten brechen diese Kohärenz auf - die Wirbel lösen sich an verschiedenen "
    "Positionen zu verschiedenen Zeiten, das Geräusch wird breitbandig statt tonal. "
    "Dasselbe Prinzip gilt für die Trennwand-Hinterkante am Faltspalt. Eine "
    "ausgefranste oder gezackte Wandkante würde den 180°-Umlenkungsverlust nicht "
    "senken, aber die Geräuschqualität ändern.", sty['B']))

st.append(Paragraph("<b>Vorderkanten-Optimierung (Leading Edge):</b> "
    "Buckelwal-inspirierte Vorderkanten (Tubercles) verhindern den gleichzeitigen "
    "Strömungsabriss über die gesamte Spannweite. Auf die Kammer übertragen: "
    "Der 30°-Strahl aus der Klappe trifft auf die Trennwand wie auf eine "
    "Flügelvorderkante. Eine wellenförmige Wandkante (anstatt der geraden Kante "
    "bei Variante A) könnte die großflächige Ablösung in viele kleine, "
    "örtlich begrenzte Ablösungen aufteilen. Das Ergebnis wäre weniger "
    "niederfrequente Druckschwankung am Zungenspalt.", sty['B']))

st.append(Paragraph("<b>Diffusor-Geometrie (Windkanal-Forschung):</b> "
    "Bei Windkraft-Diffusoren (Shrouded Turbines) ist bekannt, dass der "
    "Diffusor-Halbwinkel unter 7° bleiben muss, um Ablösung zu vermeiden. "
    "Die schräge Trennwand (Variante B) mit 3,7° liegt sicher darunter. "
    "Aber die Windkraftforschung zeigt auch: Stetige Krümmung (wie bei Variante C) "
    "ist einem geraden Diffusor bei gleichem Öffnungswinkel überlegen, weil "
    "der Druckgradient gleichmäßiger verteilt wird. Die Analogie stützt die "
    "Erwartung, dass C besser als B abschneidet.", sty['B']))

st.append(Paragraph(
    "Die Übertragbarkeit hat Grenzen: Propeller und Windräder arbeiten bei "
    "Reynolds-Zahlen von 10<super>5</super>-10<super>7</super>, die Stimmzungenkammer bei "
    "Re ~ 50-1500. Bei niedrigem Re sind Viskositätseffekte stärker, "
    "Grenzschichten dicker, und Übergänge zwischen laminar und turbulent "
    "verlaufen anders. Dennoch sind die geometrischen Prinzipien "
    "(stetige Krümmung, Vermeidung scharfer Kanten, Diffusor-Winkel) "
    "Reynolds-unabhängig gültig.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 10: Kammer-Zunge Kopplung: Serien- statt Parallelkreis
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 10: Warum resonante Kammern die Ansprache verschlechtern", sty['Ch']))

st.append(Paragraph(
    "Ein verbreiteter Gedanke ist, dass eine resonante Kammer die Zunge unterstützen "
    "könnte - ähnlich wie ein Resonanzkörper einer Gitarre. Praktische Messungen "
    "zeigen das Gegenteil: Kammern, deren Helmholtz-Frequenz ein ganzzahliges "
    "Vielfaches der Zungenfrequenz ist, bringen ein schlechteres Ansprechverhalten "
    "als nicht-resonante Kammern.", sty['B']))

st.append(Paragraph(
    "Die Erklärung liegt in der Art der Kopplung. In der Elektrotechnik gibt es "
    "zwei grundlegend verschiedene Resonanzkreis-Topologien:", sty['B']))

st.append(Paragraph(
    "<b>Parallelkopplung</b> (Schwingkreis mit L und C parallel): Beide Elemente "
    "werden gleichzeitig erregt, die Energie pendelt zwischen ihnen. Der Strom "
    "durch die Quelle ist minimal bei Resonanz - die Impedanz ist maximal. "
    "Ein Gitarrenkörper arbeitet so: Die Saitenenergie koppelt über den Steg "
    "gleichzeitig in den Körper, Decke und Luft schwingen gemeinsam, die "
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
    "Jede Kammer-Resonanz fügt dem System eine eigene Zeitkonstante hinzu.", sty['Key']))

st.append(Paragraph(
    "Was passiert bei einer resonanten Kammer? Wenn die Helmholtz-Frequenz ein "
    "Vielfaches der Zungenfrequenz ist, kann die Kammer in einem der Zugngen-Obertöne "
    "mitschwingen. Dabei wird Energie aus der Hauptschwingung der Zunge in die "
    "Kammer-Eigenschwingung umgeleitet. Die Kammer speichert diese Energie und gibt "
    "sie zeitversetzt zurück. Für den Amplitudenaufbau der Zunge wirkt das wie "
    "eine zusätzliche Dämpfung: Ein Teil der Bernoulli-Energie, die die Zunge "
    "aufbauen sollte, verschwindet in der Kammerresonanz.", sty['B']))

st.append(Paragraph(
    "In der Elektrotechnik-Analogie: Der Serienkreis hat bei Resonanz minimale "
    "Impedanz. Das bedeutet maximaler Stromfluss - aber auch maximale Verluste "
    "in den parasitären Widerständen des Kreises. Übertragen: Bei Kammerresonanz "
    "ist der Volumenstrom durch den Spalt zwar maximal, aber die Phasenbeziehung "
    "zwischen Druck und Zungenposition verschiebt sich. Die Bernoulli-Rückkopplung "
    "braucht eine bestimmte Phasenlage, um Energie in die Zunge zu pumpen. Wenn die "
    "Kammer die Phase verschiebt, wird die Rückkopplung weniger effizient - "
    "oder im schlimmsten Fall sogar bremsend.", sty['B']))

st.append(Paragraph(
    "Die praktische Konsequenz: <b>Die Kammer soll akustisch möglichst unsichtbar "
    "sein.</b> Das heißt: Ihre Eigenfrequenzen sollen weit von der Zungenfrequenz "
    "und deren Obertönen entfernt liegen. Mit f<sub>H</sub> = " + f"{f_H:.0f}" + " Hz bei einer "
    "50-Hz-Zunge ist das Verhältnis 9,2 - kein ganzzahliges Vielfaches. Das ist "
    "günstig. Eine Kammer mit f<sub>H</sub> = 500 Hz (10x) wäre problematisch, "
    "weil der 10. Oberton der Zunge (10 x 50 = 500 Hz) in die Kammerresonanz "
    "fallen würde.", sty['Warn']))

st.append(Paragraph(
    "Praktische Erfahrung zeigt, dass die Trennwand das Einschwingverhalten "
    "entscheidend beeinflusst - obwohl alle Varianten im stationären Zustand "
    "nahezu identische Volumenströme liefern. Die Erklärung liegt in der "
    "Serienkopplung: Die Kammer bestimmt nicht den Durchfluss, sondern die "
    "akustische Impedanz am Spalt. Die Wandform beeinflusst die effektive "
    "Kammer-Güte Q<sub>H</sub> und damit, wie schmalbandig oder breitbandig "
    "die Kammer an die Obertöne der Zunge koppelt (siehe Kap. 10b).", sty['B']))

st.append(Paragraph(
    "Die Optimierung hat damit zwei Ebenen: Die Zunge selbst (Q, Aufbiegung, "
    "Material) bestimmt die Grundansprache. Die Kammergeometrie (Volumen, "
    "Halsquerschnitt -> f<sub>H</sub>) bestimmt die Phasenlage des "
    "Energienachschubs. Die Wandform moduliert die Kammer-Güte Q<sub>H</sub>. "
    "Alles muss zusammenpassen.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 10b: Akustische Phasenanalyse
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 10b: Warum die Phasenlage entscheidet - nicht die Strömungsqualität (v8 neu)", sty['Ch']))

st.append(Paragraph(
    "Die bisherige Argumentation (<i>saubere Strömungsführung = bessere Ansprache</i>) "
    "greift zu kurz. Entscheidend ist nicht, wie gleichmäßig die Strömung am Spalt "
    "ankommt, sondern ob der <b>Energienachschub mit der richtigen Phasenlage</b> an "
    "der Zunge eintrifft. Eine perfekt gleichmäßige Strömung hilft nichts, wenn die "
    "Phase nicht stimmt.", sty['Warn']))

st.append(Paragraph(
    "Die Kammer ist ein akustischer Helmholtz-Resonator. Jeder Druckpuls, den die schwingende "
    "Zunge erzeugt (Spalt verengt sich -> Druckstoß), ändert den Kammerdruck. "
    "Damit dieser Druckaufbau die Schwingung "
    "<b>verstärkt</b>, muss er genau dann am Spalt wirken, wenn die Zunge gerade "
    "Energie aufnehmen kann - er muss <b>phasenrichtig</b> sein.", sty['B']))

st.append(Paragraph("Akustische Laufzeiten", sty['Sec']))

lauf_rows = [['Abschnitt', 'Länge [mm]', 'Laufzeit [ms]'],
    ['Kanal A', f'{L_kanal_A*1e3:.0f}', f'{L_kanal_A/c_s*1e3:.3f}'],
    ['Faltung', f'{L_faltung*1e3:.0f}', f'{L_faltung/c_s*1e3:.3f}'],
    ['Kanal B', f'{L_kanal_B*1e3:.0f}', f'{L_kanal_B/c_s*1e3:.3f}'],
    [Paragraph('<b>Einfach (Spalt->Klappe)</b>', _tbl_sty),
     f'{L_gesamt_einfach*1e3:.0f}', f'{t_einfach*1e3:.3f}'],
    [Paragraph('<b>Rundlauf (hin+zurück)</b>', _tbl_sty),
     f'{L_gesamt_rund*1e3:.0f}', f'{t_rund*1e3:.3f}']]
t_lauf = Table(lauf_rows, colWidths=[45*mm, 30*mm, 25*mm])
t_lauf.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_lauf)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    f"Rundlaufzeit: <b>{t_rund*1e3:.2f} ms</b>. "
    f"Schwingungsperiode bei 50 Hz: <b>{T_50*1e3:.0f} ms</b>. "
    f"Verhältnis: {T_50/t_rund:.0f} Rundläufe pro Periode. "
    "Die Kammer ist akustisch <b>sehr klein</b> gegenüber dem Grundton.", sty['B']))

st.append(Paragraph("Phasen-Relevanz: Grundton vs. Obertöne", sty['Sec']))

phase_rows = [Prow(['Oberton n', 'f [Hz]', 'T [ms]', 'Rundl./T',
               'phi [Grad]', 'Wirkung'], hdr=True)]
for n in [1, 2, 5, 10, 15, 20]:
    f_n = n * 50
    T_n = 1.0 / f_n
    rundl = T_n / t_rund
    phi = (360 * t_rund * f_n) % 360
    if phi > 180: phi = 360 - phi
    if phi < 30: wirkung = 'KONSTRUKTIV'
    elif phi > 150: wirkung = 'DESTRUKTIV'
    else: wirkung = f'neutral'
    phase_rows.append([f'{n}', f'{f_n}', f'{T_n*1e3:.2f}',
                       f'{rundl:.1f}', f'{phi:.0f}', wirkung])

t_phase = Table(phase_rows, colWidths=[18*mm, 18*mm, 18*mm, 18*mm, 20*mm, 28*mm])
t_phase.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_phase)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "<b>Für den Grundton (50 Hz):</b> "
    f"Mit {T_50/t_rund:.0f} Rundläufen pro Periode mittelt sich die Phase. "
    "Der Druckpuls ist quasi-statisch - die Wandform beeinflusst den Grundton "
    "kaum über die Phase. Dies erklärt, warum alle drei Varianten im stationären "
    "Zustand identischen Durchfluss liefern.", sty['B']))

st.append(Paragraph(
    "<b>Für Obertöne (ab ~500 Hz):</b> "
    "Ab dem 10. Oberton macht der Puls nur noch 2 Rundläufe pro Periode - die "
    "Phase wird entscheidend. Ein Puls, der mit der falschen Phase ankommt, "
    "bremst die Schwingung statt sie anzutreiben. Die Wandform bestimmt "
    "nun nicht den Durchfluss, sondern wie viel Energie als kohärenter Puls "
    "mit definierter Phase am Spalt ankommt.", sty['Key']))

st.append(Paragraph("Reflexionen und Kohärenz: A vs. B vs. C", sty['Sec']))

st.append(Paragraph(
    "Ein Transmission-Line-Modell mit diskreten Impedanzsprüngen zeigt: "
    "Der Spalt-zu-Kanal-Übergang (4 mm<super>2</super> -> 950 mm<super>2</super>) "
    "reflektiert für <b>alle drei Varianten gleich</b> 98,8% der Wellenenergie. "
    "Ein Druckpuls schafft es nie als kohärente Welle durch diesen Sprung. "
    "Das TL-Modell kann die Wandformen nicht unterscheiden.", sty['Warn']))

st.append(Paragraph(
    "Die physikalisch korrekte Sichtweise: Die Bernoulli-Rückkopplung ist ein "
    "<b>lokaler Effekt am Spalt</b>. Die Kammer wirkt nicht als Wellenleiter, "
    "sondern als akustisches Reservoir mit Impedanz Z<sub>Kammer</sub>(f). "
    "Die Phase dieser Impedanz bestimmt, ob Energie in die Zunge gepumpt wird.", sty['B']))

st.append(Paragraph("Kammerimpedanz und Oberton-Kopplung", sty['Sec']))

st.append(Paragraph(
    "Die Kammerimpedanz Z(f) hat die Phase eines Helmholtz-Resonators: "
    "Phase(Z) = arctan(Q<sub>H</sub> x (f/f<sub>H</sub> - f<sub>H</sub>/f)). "
    "Bei f = f<sub>H</sub> ist die Phase Null (Resonanz, maximaler Energieaustausch). "
    "Weit unterhalb: Phase ~ -90° (Compliance-dominiert). "
    "Weit oberhalb: Phase ~ +90° (Masse-dominiert).", sty['B']))

# Phase-Tabelle im PDF
phase_pdf_rows = [Prow(['Oberton n', 'f [Hz]', 'f/f<sub>H</sub>',
                   'Phase Z [Grad]', 'Kopplung'], hdr=True)]
for n, f_n, x, phi_Z, wirkung in phase_results:
    if n in [1, 2, 5, 6, 8, 10, 15, 20]:
        phase_pdf_rows.append([f'{n}', f'{f_n}', f'{x:.2f}',
                               f'{phi_Z:.1f}', wirkung.replace('STARKE Kopplung (Resonanz)', 'STARK (Resonanz)')])

t_ph = Table(phase_pdf_rows, colWidths=[18*mm, 18*mm, 18*mm, 24*mm, 40*mm])
t_ph.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_ph)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    f"Der Grundton (50 Hz, f/f<sub>H</sub> = {50/f_H:.2f}) liegt weit unter der "
    f"Kammerresonanz. Phase ~ -90° -> fast reine Compliance -> <b>schwache Kopplung</b>. "
    f"Das ist günstig: Die Kammer stört den Grundton nicht. Erst Obertöne nahe "
    f"f<sub>H</sub> = {f_H:.0f} Hz koppeln stark an die Kammer.", sty['B']))

st.append(Paragraph("Was die Wandform wirklich beeinflusst", sty['Sec']))

st.append(Paragraph(
    "Die <b>Phase</b> der Kammerimpedanz wird durch f<sub>H</sub> und Q<sub>H</sub> "
    "bestimmt - also durch Volumen, Halsquerschnitt und Halslänge. Die Wandform "
    "ändert diese Größen kaum (alle drei Varianten haben dasselbe Volumen und "
    "denselben Schlitz als Hals).", sty['B']))

st.append(Paragraph(
    "Was die Wandform beeinflusst ist die <b>effektive Kammer-Güte</b> Q<sub>H</sub>:", sty['B']))

st.append(Paragraph(
    "<b>Gerade Wand (A):</b> Der Stagnationspunkt erzeugt großskalige Wirbel. "
    "Diese dissipieren Energie -> Q<sub>H</sub> wird gesenkt -> das Resonanzband wird "
    "<b>breiter</b>. Mehr Obertöne werden beeinflusst, aber jeder einzelne schwächer.", sty['B']))

st.append(Paragraph(
    "<b>Parabolische Wand (C):</b> Weniger Dissipation in der Kammer -> Q<sub>H</sub> "
    "bleibt höher -> das Resonanzband wird <b>schmaler</b>. Weniger Obertöne betroffen, "
    "aber die betroffenen werden stärker beeinflusst.", sty['B']))

st.append(Paragraph(
    "<b>Die Konsequenz für die Praxis:</b> Die Wandform ist kein Strömungsoptimierungs-"
    "Problem ('sauberer = besser'), sondern ein <b>Abstimmungsproblem</b>. "
    "C verstärkt die Phasenkopplung in <b>beide</b> Richtungen: Wenn die Kammer-"
    "Obertöne günstig zum Ton liegen, profitiert C mehr als A. Wenn sie ungünstig "
    "liegen (destruktive Phase), kann A mit seinem breiteren, aber schwächeren "
    "Resonanzband sogar besser abschneiden. Das erklärt, warum erfahrene "
    "Instrumentenbauer die Kammergröße für verschiedene Tonhöhen variieren: "
    "Sie stimmen f<sub>H</sub> ab, nicht die Strömung.", sty['Key']))

st.append(Paragraph("Warum die Faltung die Phase ändert", sty['Sec']))

st.append(Paragraph(
    "Die 180°-Faltung ist eine <b>akustische Diskontinuität</b>. An der scharfen "
    "Umlenkung ändert sich die effektive Querschnittsfläche abrupt und die "
    "Strömungsrichtung kehrt sich um. Ein Druckpuls, der hier ankommt, wird "
    "teilweise reflektiert. Der reflektierte Anteil läuft zurück zum Spalt und "
    "überlagert sich dort mit dem direkten Kammerdruck.", sty['B']))

st.append(Paragraph(
    "Diese Überlagerung ist der Mechanismus: Der reflektierte Puls kommt mit "
    "einer Zeitverzögerung an (Doppelter Faltungs-Abstand / Schallgeschwindigkeit). "
    "Je nach Frequenz ist die Überlagerung konstruktiv oder destruktiv - "
    "die <b>effektive Phase</b> des Kammerdrucks am Spalt wird verschoben. "
    "Zusätzlich ändert die Faltung die <b>effektive akustische Masse</b> des Halses "
    "(der Luftpfad Schlitz-Kammer wird länger und gewundener), was f<sub>H</sub> "
    "geringfügig nach unten verschiebt. Beides zusammen bewirkt: Unterschiedliche "
    "Faltungsgeometrien erzeugen unterschiedliche Phasenlagen am Spalt.", sty['B']))

st.append(Paragraph(
    "Die Abrundung (Viertelkreis R = 8-10 mm) reduziert die Stärke der Reflexion "
    "an der Faltung: Ein stetiger Querschnittsübergang erzeugt weniger Reflexion "
    "als ein abrupter Knick. Das bedeutet weniger Phasenverzerrung am Spalt - "
    "nicht unbedingt eine <b>bessere</b> Phase, aber eine <b>berechenbarere</b>: "
    "Die Phase kommt näher an den idealen Helmholtz-Wert heran, der durch "
    "V, S<sub>Hals</sub> und L<sub>Hals</sub> bestimmt ist.", sty['B']))

st.append(Paragraph("Hierarchie der Effekte für die Ansprache", sty['Sec']))

st.append(Paragraph(
    "Die Kammergeometrie beeinflusst die Ansprache über <b>mehrere Mechanismen</b>. "
    "Diese sind physikalisch verschieden und unterschiedlich wichtig:", sty['B']))

st.append(Paragraph(
    "<b>1. Akustische Phasenkopplung (dominant):</b> "
    "Die Phase der Kammerimpedanz Z(f) bestimmt, ob Obertöne der Zunge Energie "
    "aus der Kammer aufnehmen oder an sie abgeben. Wird primär durch f<sub>H</sub> "
    "bestimmt (Kammervolumen, Halsquerschnitt, Halslänge). Die Wandform hat darauf "
    "nur geringen Einfluss (über Q<sub>H</sub> und Reflexionen). "
    "Dies ist der Hauptgrund, warum die <b>Kammergröße</b> für verschiedene "
    "Tonhöhen variiert wird.", sty['B']))

st.append(Paragraph(
    "<b>2. Kammer-Güte Q<sub>H</sub> (sekundär):</b> "
    "Bestimmt die Bandbreite der Kammer-Zungen-Kopplung. Höheres Q<sub>H</sub> "
    "(weniger Dissipation, z.B. Variante C) -> schmalere, stärkere Kopplung. "
    "Niedrigeres Q<sub>H</sub> (mehr Dissipation, z.B. Variante A) -> breitere, "
    "schwächere Kopplung. Die Wandform beeinflusst dies direkt über die "
    "Strömungsverluste in der Kammer.", sty['B']))

st.append(Paragraph(
    "<b>3. Anströmwinkel am Spalt (untergeordnet, aber praktisch relevant):</b> "
    "Das Düsenexperiment (Kap. 5b) zeigt eindeutig: Der Winkel, unter dem "
    "die Luft auf die Zunge trifft, bestimmt das Verhältnis zwischen "
    "Normalkomponente (statischer Druck) und Tangentialkomponente (Bernoulli-Sog). "
    "Nur bei ca. 20-40° stimmt die Balance für maximale Energiezufuhr. "
    "Die Wandform beeinflusst diesen Winkel, weil sie die Strömungsrichtung "
    "in Kanal A lenkt, bevor die Luft in den Spalt einbiegt.", sty['B']))

st.append(Paragraph(
    "Dieser dritte Effekt erklärt die <b>praktische Erfahrung</b> der "
    "Instrumentenbauer, dass die Trennwandform hörbar wirkt, obwohl "
    "der stationäre Durchfluss identisch ist und die Phasenkopplung "
    "primär durch die Kammergröße bestimmt wird. Die gleichmäßige "
    "Luftführung selbst ist untergeordnet - aber der <b>Winkel</b>, "
    "unter dem die Luft am Spalt ankommt, ist es nicht. Das ist ein "
    "aerodynamischer Effekt, der vom akustischen Phasenmechanismus "
    "unabhängig ist und vor allem in den ersten Schwingungszyklen "
    "(dem Einschwingvorgang) eine Rolle spielt.", sty['B']))

st.append(Paragraph(
    "<b>Zusammenfassung:</b> Die Kammergeometrie beeinflusst die Ansprache "
    "über drei verschiedene Kanäle: (1) Akustische Phase (dominant, durch "
    "Kammergröße), (2) Resonanzbandbreite (Q<sub>H</sub>, durch Wandform), "
    "(3) Anströmwinkel (durch Wandform, besonders beim Einschwingvorgang). "
    "Dass alle drei Kanäle existieren, erklärt, warum weder eine rein "
    "akustische Optimierung noch eine rein strömungstechnische allein zum "
    "besten Ergebnis führt - der erfahrene Instrumentenbauer optimiert "
    "intuitiv alle drei gleichzeitig.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 10c: Vorkammer-Analyse
# ═══════════════════════════════════════════════════════════════
st.append(PageBreak())
st.append(Paragraph("Kapitel 10c: Vorkammer zwischen Klappe und Hauptkammer (v8 neu)", sty['Ch']))

st.append(Paragraph(
    "Was passiert, wenn zwischen Klappe und Hauptkammer eine zusätzliche Kammer "
    "eingebaut wird? Abmessungen: 48 x 20 x 70 mm "
    f"(Volumen {V_vk*1e6:.1f} cm<super>3</super>). "
    "Die Kopplung zur Hauptkammer erfolgt über eine Öffnung mit den Abmessungen "
    "der Klappenöffnung.", sty['B']))

st.append(Paragraph("Geometrie", sty['Sec']))

vk_geom_rows = [
    ['Parameter', 'Vorkammer', 'Hauptkammer'],
    ['Volumen', f'{V_vk*1e6:.1f} cm3', f'{V_eff*1e6:.1f} cm3'],
    ['Querschnitt', f'{S_vk_quer*1e6:.0f} mm2', f'{B_kam*H_kam*1e6:.0f} mm2'],
    ['Kopplung/Hals', f'{S_kopplung*1e6:.0f} mm2 (Klappenoffng.)', f'{S_hals*1e6:.0f} mm2 (Schlitz)'],
    [Paragraph('<b>f<sub>H</sub> (einzeln)</b>', _tbl_sty),
     f'{f_vk_allein:.0f} Hz', f'{f_H:.0f} Hz'],
]
t_vk_geom = Table(vk_geom_rows, colWidths=[35*mm, 40*mm, 40*mm])
t_vk_geom.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_vk_geom)
st.append(Spacer(1, 2*mm))

st.append(Paragraph("Gekoppelte Helmholtz-Resonatoren", sty['Sec']))

st.append(Paragraph(
    "Zwei Kammern mit verbindendem Hals bilden ein <b>gekoppeltes Zweimassensystem</b>. "
    "Die akustischen Massen (M = rho x L/S in den Hälsen) und Compliances "
    "(C = V/(rho x c<super>2</super>) der Kammern) ergeben eine 2x2-Eigenwertgleichung. "
    "Die Eigenfrequenzen spalten sich auf:", sty['B']))

vk_mode_rows = [
    ['', 'Ohne VK', 'Mit VK'],
    ['Mode 1 (gleichphasig)', f'{f_H:.0f} Hz', f'{f_mode1:.0f} Hz'],
    ['Mode 2 (gegenphasig)', '--', f'{f_mode2:.0f} Hz'],
    ['Anti-Resonanz', '--', f'~{f_anti:.0f} Hz'],
    [Paragraph('<b>Einfaches Modell</b> (V<sub>ges</sub>)', _tbl_sty),
     f'{f_H:.0f} Hz', f'{f_H_gesamt:.0f} Hz'],
]
t_vk_mode = Table(vk_mode_rows, colWidths=[40*mm, 35*mm, 35*mm])
t_vk_mode.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#533483')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_vk_mode)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    f"<b>Mode 1 ({f_mode1:.0f} Hz):</b> Beide Kammern schwingen gleichphasig "
    f"(p<sub>VK</sub>/p<sub>HK</sub> = {ratio_1:.1f}). "
    f"Gegenüber dem Original ({f_H:.0f} Hz) um {(1-f_mode1/f_H)*100:.0f}% abgesenkt - "
    f"das zusätzliche Volumen vergrößert die effektive Compliance.", sty['B']))

st.append(Paragraph(
    f"<b>Mode 2 ({f_mode2:.0f} Hz):</b> Kammern schwingen gegenphasig "
    f"(p<sub>VK</sub>/p<sub>HK</sub> = {ratio_2:.2f}). "
    f"Liegt weit über dem relevanten Frequenzbereich - für die 50-Hz-Zunge "
    f"spielt dieser Mode keine Rolle.", sty['B']))

st.append(Paragraph(
    f"<b>Anti-Resonanz (~{f_anti:.0f} Hz):</b> Bei dieser Frequenz wirkt die VK als "
    f"akustischer Tilger: Die Luft im Kopplungshals schwingt so, dass sie den "
    f"Durchfluss blockiert. In der Hauptkammer entsteht maximaler Druck.", sty['B']))

st.append(Paragraph("Phasenverschiebung am Spalt", sty['Sec']))

vk_phase_rows = [Prow(['Oberton n', 'f [Hz]',
                   'Phase ohne VK [Grad]', 'Phase mit VK [Grad]', 'Änderung'], hdr=True)]
for n, f_n, phi_o, phi_m in phase_vk_results:
    if n in [1, 3, 5, 6, 8, 10, 15]:
        delta = phi_m - phi_o
        if abs(phi_m) < abs(phi_o):
            aend = 'BESSER'
        elif abs(phi_m) > abs(phi_o) + 5:
            aend = 'schlechter'
        else:
            aend = '~gleich'
        vk_phase_rows.append([f'{n}', f'{f_n}', f'{phi_o:.1f}', f'{phi_m:.1f}', aend])

t_vk_phase = Table(vk_phase_rows, colWidths=[18*mm, 18*mm, 30*mm, 30*mm, 22*mm])
t_vk_phase.setStyle(TableStyle([
    ('BACKGROUND',(0,0),(-1,0),HexColor('#16213e')),('TEXTCOLOR',(0,0),(-1,0),white),
    ('FONTSIZE',(0,0),(-1,-1),7.5),('FONTNAME',(0,0),(-1,0),'Helvetica-Bold'),
    ('ALIGN',(0,0),(-1,-1),'CENTER'),('GRID',(0,0),(-1,-1),0.5,HexColor('#ccc')),
    ('ROWBACKGROUNDS',(0,1),(-1,-1),[white,HexColor('#f5f5f5')]),
    ('TOPPADDING',(0,0),(-1,-1),1.5*mm),('BOTTOMPADDING',(0,0),(-1,-1),1.5*mm)]))
st.append(t_vk_phase)
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "Der entscheidende Effekt: Die VK verschiebt Mode 1 von 274 auf "
    f"{f_mode1:.0f} Hz. Dadurch ändert sich die Phasenlandschaft für die "
    "Obertöne der 50-Hz-Zunge:", sty['B']))

st.append(Paragraph(
    f"<b>5. Oberton (250 Hz):</b> Lag ohne VK knapp UNTER der Resonanz "
    f"(Phase -53°). Liegt jetzt knapp UEBER Mode 1 ({f_mode1:.0f} Hz) -> "
    f"Phase springt auf +42°. Die Phase ist näher an 0° -> "
    f"<b>bessere Kopplung</b>.", sty['B']))

st.append(Paragraph(
    f"<b>6. Oberton (300 Hz):</b> Lag ohne VK knapp UEBER der Resonanz "
    f"(Phase +51°). Liegt jetzt noch weiter oberhalb -> Phase +71° -> "
    f"<b>schlechtere Kopplung</b>.", sty['B']))

st.append(Paragraph(
    "<b>Grundton und tiefe Obertöne (50-150 Hz):</b> Kaum verändert. "
    "Beide liegen weit unter jeder Kammerresonanz -> reine Compliance -> "
    "Phase ~ -90° mit und ohne VK.", sty['B']))

st.append(Paragraph("Bewertung", sty['Sec']))

st.append(Paragraph(
    f"<b>Stationärer Durchfluss:</b> Der zusätzliche Kopplungshals "
    f"(zeta = 1,5 auf {S_kopplung*1e6:.0f} mm<super>2</super>) ergibt "
    f"bezogen auf den Spalt nur zeta = {z_kopplung_scaled:.2e} -> "
    f"vernachlässigbar. Die VK ändert den Volumenstrom nicht.", sty['B']))

st.append(Paragraph(
    "<b>Akustisch:</b> Die VK senkt die Hauptresonanz um 16% (von 274 auf "
    f"{f_mode1:.0f} Hz). Das kann gezielt genutzt werden, um die Resonanz "
    "von bestimmten Obertönen wegzuschieben - oder heranzuschieben. "
    "Die Wirkung hängt kritisch davon ab, welche Obertöne in die Nähe "
    "der neuen Mode 1 fallen.", sty['B']))

st.append(Paragraph(
    f"<b>Für die 50-Hz-Zunge konkret:</b> Gemischtes Ergebnis. Der 5. OT "
    f"(250 Hz) koppelt besser (Phase von -53° auf +42°, näher an 0°), "
    f"aber der 6. OT (300 Hz) koppelt schlechter (von +51° auf +71°, "
    f"weiter von 0°). Ob die Gesamtwirkung positiv oder negativ ist, "
    f"hängt von der relativen Stärke dieser Obertöne ab und muss "
    f"am Instrument erprobt werden.", sty['B']))

st.append(Paragraph(
    "<b>Praktische Empfehlung:</b> Die VK ist ein <b>Abstimmungswerkzeug</b>: "
    "Durch Variation der VK-Länge (und damit des Volumens) lässt sich "
    "Mode 1 kontinuierlich verschieben. Optimal wäre eine VK-Länge, bei der "
    "ein wichtiger Oberton genau auf Mode 1 fällt (Phase = 0°, maximale "
    "Kopplung). Da die VK den stationären Durchfluss nicht beeinflusst, "
    "ist sie ein reines Akustik-Tuning ohne Nebenwirkungen auf die "
    "Grundansprache.", sty['Key']))

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
ax1.set_ylabel('Abstand nächster Oberton zu f_H [Hz]', fontsize=10)
ax1.set_title(f'Welche Basstöne werden durch die Kammerresonanz (f_H = {f_H:.0f} Hz) benachteiligt?\n'
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
st.append(Paragraph("Kapitel 11: Welche Basstöne sind durch die Kammerresonanz benachteiligt?", sty['Ch']))

st.append(Paragraph(
    "Die Helmholtz-Frequenz dieser Kammer liegt bei " + f"{f_H:.0f}" + " Hz. Die Frage ist: "
    "Welche Basstöne haben einen Oberton, der in die Nähe dieser Resonanz fällt? "
    "Denn wie in Kapitel 10 erläutert, entzieht eine resonante Kammer der Zunge Energie "
    "über die Serienkopplung. Nicht der Grundton muss betroffen sein - es genügt, wenn "
    "der n-te Oberton (n x f<sub>Grundton</sub>) auf f<sub>H</sub> trifft.", sty['B']))

st.append(Paragraph(
    "Die Bandbreite der Helmholtz-Resonanz hängt von ihrer eigenen Güte Q<sub>H</sub> ab. "
    "Eine Holzkammer mit Klappe hat typisch Q<sub>H</sub> ~ 5-10. Bei Q<sub>H</sub> = 10 "
    "ist das kritische Band " + f"{f_H:.0f}" + " +/- " + f"{bw_narrow/2:.0f}" + " Hz, "
    "bei Q<sub>H</sub> = 5 erweitert es sich auf +/- " + f"{bw_wide/2:.0f}" + " Hz.", sty['B']))

# Diagramm einbetten
st.append(Spacer(1, 2*mm))
st.append(Image(diag_path, width=170*mm, height=130*mm))
st.append(Spacer(1, 2*mm))

st.append(Paragraph(
    "Abbildung 4 zeigt das Ergebnis für alle Basstöne der temperierten Stimmung von "
    "C1 (32,7 Hz) bis D4 (293,7 Hz). Oberes Diagramm: Abstand des nächsten Obertons "
    "zu f<sub>H</sub>. Rot = kritisch (Oberton trifft Resonanz +/- " + f"{bw_narrow/2:.0f}" + " Hz), "
    "orange = ungünstig, grün = unbeeinflusst. Die Zahl über jedem Balken gibt die "
    "Oberton-Nummer an. Unteres Diagramm: Oberton-Landkarte - jeder Punkt ist ein "
    "Oberton, rote Rauten liegen im Resonanzband.", sty['B']))

# Betroffene Noten als Tabelle
kritisch = [r for r in res if r['severity'] == 2]
# Gruppiere nach Muster
st.append(Paragraph("Auffällige Muster:", sty['Sec']))

st.append(Paragraph(
    "<b>Unter ~70 Hz (C2) ist fast jeder Ton kritisch betroffen.</b> Der Grund: Bei "
    "tiefen Tönen liegen die Obertöne dicht beieinander (Abstand = Grundfrequenz). "
    "Ein Ton von 40 Hz hat Obertöne bei 40, 80, 120, ..., 440, 460, 480 Hz - "
    "das Raster ist so fein, dass immer ein Oberton in das Resonanzband fällt. "
    "Je tiefer der Ton, desto unvermeidlicher ist der Konflikt.", sty['Key']))

st.append(Paragraph(
    "<b>Besonders getroffen</b> sind Töne, deren Grundfrequenz ein ganzzahliger "
    "Teiler von f<sub>H</sub> ist:", sty['B']))

# Teiler-Tabelle
div_rows = [['f_H / n', 'f [Hz]', 'Nächste Note', 'Differenz', 'Oberton n']]
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
    "Die schärfsten Treffer: <b>F#1</b> (46,2 Hz, 10. Oberton = 462,5 Hz, nur 1,3 Hz "
    "Abstand), <b>F#2</b> (92,5 Hz, 5. Oberton) und <b>B/Bb</b> in allen Oktaven "
    "(8./4./2. Oberton = 466,2 Hz, 5,0 Hz Abstand). Diese Töne werden bei dieser "
    "Kammergröße systematisch benachteiligt.", sty['Warn']))

st.append(Paragraph(
    "Ab C3 (130,8 Hz) aufwärts werden die Lücken zwischen den Obertönen größer, "
    "und die meisten Töne liegen ausserhalb des Resonanzbandes. Das erklärt, warum "
    "das Problem in der Praxis vor allem bei den tiefsten Basstönen auffällt - "
    "und warum erfahrene Instrumentenbauer die Kammergröße für verschiedene "
    "Tonbereiche variieren.", sty['B']))

# ═══════════════════════════════════════════════════════════════
# Kap 12: Grenzen der Berechnung
# ═══════════════════════════════════════════════════════════════
st.append(Paragraph("Kapitel 12: Was Berechnungen leisten - und was nicht", sty['Ch']))

st.append(Paragraph(
    "Die in diesem Dokument vorgestellten Berechnungen stützen sich auf "
    "gut etablierte physikalische Modelle: Euler-Bernoulli-Balkentheorie für die "
    "Zungenfrequenz, Bernoulli-Gleichung für die Spaltströmung, Helmholtz-Resonator "
    "für die Kammerakustik, Verlustbeiwerte aus der Strömungsmechanik (Idelchik). "
    "Diese Modelle sind in ihrem Gültigkeitsbereich zuverlässig und liefern "
    "Ergebnisse, die in der richtigen Größenordnung liegen.", sty['B']))

st.append(Paragraph(
    "Dennoch ersetzen sie keinen praktischen Test. Die Gründe dafür sind nicht "
    "Schwächen der Physik, sondern Vereinfachungen in den Randbedingungen:", sty['B']))

st.append(Paragraph(
    "<b>1. Idealisierte Geometrie:</b> Die Berechnung geht von exakt planen "
    "Flächen, scharfen Kanten und gleichmäßigen Spalten aus. In der Realität "
    "hat Holz Fasern, Aluminium Bearbeitungsspuren, und jede Zunge ist minimal "
    "anders aufgebogen. Schon 0,1 mm Unterschied in der Aufbiegung ändert die "
    "effektive Spaltfläche um rund 7 Prozent.", sty['B']))

st.append(Paragraph(
    "<b>2. Stationäre Strömung:</b> Die Verlustberechnung nimmt stationäre "
    "Strömung an - einen gleichmäßigen Luftstrom. Tatsächlich pulsiert die "
    "Strömung mit 50 Hz, weil die Zunge den Spalt periodisch öffnet und "
    "schließt. Die instationären Effekte (Anlaufwirbel, Schwinger-Strömung-"
    "Kopplung, akustische Rückwirkung) sind analytisch schwer zu fassen und "
    "ändern die Verlustbeiwerte um geschätzt 10-30 Prozent. "
    "<b>Dies ist die gravierendste Einschränkung:</b> Praktische Erfahrung zeigt, "
    "dass die Trennwandform das Einschwingverhalten entscheidend beeinflusst - "
    "ein Effekt, der im instationären Anlauf (Impulsantwort) liegt und von der "
    "stationären Berechnung grundsätzlich nicht erfasst wird (siehe Kapitel 6).", sty['B']))

st.append(Paragraph(
    "<b>3. Dämpfung:</b> Die Güte Q ist der empfindlichste Parameter der "
    "gesamten Analyse - sie bestimmt die Einschwingzeit direkt proportional. "
    "Q lässt sich aus den Materialwerten und der Geometrie nur grob abschätzen "
    "(Kapitel 8). Der tatsächliche Wert hängt von der konkreten Einspannung, "
    "dem Anpressdruck der Schrauben, der Oberflächengüte an der Kontaktfläche "
    "und dem Zustand der Zunge ab. Unterschiede von Faktor 2 zwischen zwei "
    "scheinbar identischen Zungen sind in der Praxis normal.", sty['B']))

st.append(Paragraph(
    "<b>4. Kammerakustik:</b> Das Helmholtz-Modell behandelt die Kammer als "
    "einzelnen Resonator mit einer einzigen Eigenfrequenz. Real hat die Kammer "
    "weitere Moden (Längs-, Quer- und Höhenmoden), die bei höheren Frequenzen "
    "auftreten. Zudem ändern sich die akustischen Eigenschaften mit dem "
    "Balgdruck, der Temperatur und der Luftfeuchtigkeit. Die Resonanzanalyse in "
    "Kapitel 11 zeigt die richtigen Tendenzen, aber die exakten Bandbreiten und "
    "Kopplungsstärken können nur gemessen werden.", sty['B']))

st.append(Paragraph(
    "<b>5. Wechselwirkungen:</b> Jeder Effekt ist einzeln modelliert. In der "
    "Wirklichkeit beeinflussen sie sich gegenseitig: Die Zungenposition ändert "
    "die Strömung, die Strömung ändert den Druck in der Kammer, der Kammerdruck "
    "wirkt auf die Klappe zurück, und die akustischen Wellen überlagern sich "
    "mit der Strömung. Eine vollständige Lösung würde eine gekoppelte "
    "Fluid-Struktur-Akustik-Simulation (FSI) erfordern.", sty['B']))

st.append(Paragraph(
    "Was die Berechnungen <b>leisten</b>:", sty['Sec']))

st.append(Paragraph(
    "Sie identifizieren zuverlässig, <b>welche Parameter wichtig sind</b> und "
    "auf welcher Ebene sie wirken. Das Ergebnis, dass der Spalt den stationären "
    "Gesamtwiderstand dominiert, ist robust. Aber die praktische Erfahrung zeigt "
    "ebenso robust, dass die Trennwandform das Einschwingverhalten entscheidend "
    "beeinflusst - ein Effekt, den die stationäre Berechnung prinzipiell nicht "
    "erfassen kann (Kapitel 6). Die Resonanzanalyse zeigt zuverlässig, welche "
    "Tonhöhen prinzipiell gefährdet sind.", sty['B']))

st.append(Paragraph(
    "Die Berechnungen liefern damit eine <b>qualifizierte Einschätzung</b>: "
    "Sie zeigen, dass Zungengüte, Aufbiegung und Einspannung die Grundansprache "
    "bestimmen. Sie ordnen die drei Trennwandvarianten in eine physikalisch "
    "begründete Rangfolge (C besser als B besser als A - bestätigt durch "
    "praktische Erfahrung). Und sie machen klar, warum tiefe Bässe träger "
    "ansprechen als hohe Töne - das folgt direkt aus tau = Q/(pi*f). "
    "Wo die Berechnung an ihre Grenzen stößt, ist die quantitative Vorhersage "
    "des Trennwand-Effekts: Dass die Wand wirkt, ist physikalisch begründbar. "
    "Wie stark, lässt sich nur messen.", sty['B']))

st.append(Paragraph(
    "Für den Instrumentenbauer bedeutet das: <b>Die Berechnung sagt, wo man "
    "suchen soll. Der praktische Test sagt, was man findet.</b> Beides zusammen "
    "ist mehr als jedes für sich allein.", sty['Key']))

# ═══════════════════════════════════════════════════════════════
# Kap 13: Zusammenfassung (ex-Kap 12)
# ═══════════════════════════════════════════════════════════════
st.append(PageBreak())
st.append(Paragraph("Kapitel 13: Zusammenfassung", sty['Ch']))
st.append(Paragraph(f"<b>Stimmplatte:</b> Keilförmig ({t_platte_einsp*1e3:.0f}->{t_platte_frei*1e3:.0f} mm). Am freien Ende wirkt der {t_platte_frei*1e3:.0f}-mm-Schlitz als kurzer Kanal mit eigenem Strömungswiderstand (zeta ~ {z_schlitz:.3f}, davon {zeta_schlitz_entry} Eintrittsverlust [K3b]).", sty['Good']))
st.append(Paragraph(f"<b>Spalt:</b> Dreieckig (0->{h_aufbiegung*1e3:.1f} mm), S<sub>eff</sub> = {S_gap*1e6:.0f} mm². Dominiert den Gesamtwiderstand. Alle Varianten: v ~ {v1k:.0f} m/s, Q ~ {Q1k*1e6:.0f} ml/s bei 1 kPa.", sty['Good']))
st.append(Paragraph(f"<b>Dynamische Spaltänderung [K6]:</b> Der Balgdruck biegt die Zunge statisch nach unten. Grenzdruck dp<sub>crit</sub> = {dp_crit:.0f} Pa identisch mit Schwellendruck. Effektiver Schwingungsbeginn bei ~{dp_dynamic_est:.0f} Pa ({dp_dynamic_est/100:.0f} mbar). Arbeitsbereich {dp_dynamic_est:.0f}-{dp_crit:.0f} Pa.", sty['Good']))
st.append(Paragraph(f"<b>Trennwand-Material:</b> Folie und Furnier sind strömungstechnisch gleichwertig. Furnier ist stabiler und akustisch neutraler. Folie spart Platz, kann aber flattern.", sty['Good']))
st.append(Paragraph(f"<b>Klappenorientierung:</b> Ausblas (zeta={zeta_klap_ausblas:.2f}) besser als Umlenk (zeta={zeta_klap_umlenk:.2f}). Ausblas nutzt Wandgeometrie voll, Umlenk nur bei Platzmangel.", sty['Good']))
st.append(Paragraph(f"<b>Trennwandform:</b> Praktisch entscheidend für die Ansprache (instationärer Anlauf). C (Parabel) beste Richtungserhaltung, B (schräg) guter Kompromiss, A (gerade) am einfachsten aber am trägsten. Stationärer Durchfluss bei allen gleich - der Unterschied liegt in der Impulsantwort (Reflexionen, Richtungsverlust). Die Turbulenz für die Selbsterregung erzeugt der Spalt selbst.", sty['Good']))
st.append(Paragraph(f"<b>Ansprache:</b> tau = Q/(pi*f). Bei Q=100: {tau_vals[0]*1e3:.0f} ms, bei Q=25: {tau_vals[2]*1e3:.0f} ms. Die Kammergeometrie beeinflusst nicht den stationären Durchfluss, aber entscheidend den instationären Anlauf (Impulsantwort).", sty['Good']))
st.append(Paragraph(f"<b>Akustik:</b> Helmholtz {f_H:.0f} Hz, Kammer akustisch klein.", sty['Good']))
st.append(Paragraph(f"<b>Phasenkopplung [v8 neu]:</b> Die Ansprache wird durch drei Mechanismen bestimmt: (1) akustische Phase (dominant, durch f<sub>H</sub>), (2) Kammer-Güte Q<sub>H</sub> (durch Wandform), (3) Anströmwinkel am Spalt (durch Wandform, beim Einschwingen). Gleichmäßige Strömung allein reicht nicht - die Phase muss stimmen.", sty['Good']))
st.append(Paragraph(f"<b>Vorkammer [v8 neu]:</b> 48x20x70 mm ({V_vk*1e6:.0f} cm<super>3</super>) senkt Mode 1 von {f_H:.0f} auf {f_mode1:.0f} Hz (-16%). 5. OT (250 Hz) koppelt besser, 6. OT (300 Hz) schlechter. VK-Länge ist Abstimmwerkzeug für die Oberton-Phasenlage ohne Einfluss auf den stationären Durchfluss.", sty['Good']))

doc.build(st)
print(f"  PDF: {path}")
shutil.copy(path, "/mnt/user-data/outputs/bass_50hz_v8.pdf")
print("  Fertig!")
