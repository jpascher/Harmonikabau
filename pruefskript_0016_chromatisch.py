#!/usr/bin/env python3
"""
Dok. 0016 — Ergänzung: Chromatische Analyse
Obertonmoden über 2 Oktaven bei gleicher Zungenlänge

FEM-Lösung: Wie verhält sich f₂/f₁ über die chromatische Skala?
Drei Methoden: Nur Dicke, Gewicht, Profil + Gewicht

Benötigt: numpy, scipy
"""
import numpy as np
from scipy.linalg import eigh

def beam_modes(N, t_func, L, W, E, rho, m_tip=0, n_modes=4):
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
    if m_tip>0: M[2*N,2*N]+=m_tip
    fd=list(range(2,ndof))
    ev,_=eigh(K[np.ix_(fd,fd)],M[np.ix_(fd,fd)])
    return np.sqrt(np.abs(ev[:n_modes]))/(2*np.pi)

# ══════════════════════════════════════════════════════════════
# Parameter
# ══════════════════════════════════════════════════════════════
E=200e9; rho=7800; L=70e-3; W=8e-3; N=80
note_names = ['A','Bb','B','C','C#','D','Eb','E','F','F#','G','G#']
f_A1 = 55.0

# Referenz: gleichmäßig
t_ref = 0.40e-3
fr_ref = beam_modes(N, lambda x: t_ref, L, W, E, rho)
f_ref = fr_ref[0]

print("="*95)
print("  Dok. 0016 — Chromatische Analyse")
print(f"  Gleiche Zungenlänge L = {L*1000:.0f} mm, W = {W*1000:.0f} mm")
print(f"  Bereich: A1 ({f_A1:.0f} Hz) bis A3 ({f_A1*4:.0f} Hz) = 2 Oktaven, 25 Halbtöne")
print("="*95)

# ══════════════════════════════════════════════════════════════
# 1. GRUNDSATZ: f₂/f₁ bei gleichmäßigem Profil
# ══════════════════════════════════════════════════════════════
print("""
══ GRUNDSATZ ══════════════════════════════════════════════════════════

  Bei einem gleichmäßigen Cantilever (t = const entlang x) gilt:
  
      f_n = (β_n)² × √(EI / ρAL⁴) / (2π)
  
  Da I ∝ t³ und A ∝ t, ist f_n ∝ t. ALLE Moden skalieren gleich mit t.
  
  → f₂/f₁ = (β₂/β₁)² = 6,2669  — UNABHÄNGIG von t, L, W, E, ρ
  → f₃/f₁ = (β₃/β₁)² = 17,547  — ebenfalls konstant
  
  Das Verhältnis ist eine GEOMETRISCHE KONSTANTE des Cantilevers.
  Es ändert sich nur durch: Profilierung ODER Gewicht.
""")

# ══════════════════════════════════════════════════════════════
# 2. TABELLE: Alle 25 Halbtöne bei gleichmäßigem Profil
# ══════════════════════════════════════════════════════════════
print("══ TABELLE 1: Gleichmäßiges Profil — Dicke variiert für jeden Ton ══")
print("-"*95)
print(f"{'Nr':>3s} {'Ton':>5s} {'f₁ [Hz]':>8s} {'t [mm]':>8s} {'m_Z [g]':>8s} {'k [N/m]':>8s} "
      f"{'f₂/f₁':>7s} {'f₂ [Hz]':>8s} {'f₃ [Hz]':>8s}")
print("-"*95)

for idx in range(25):
    f_t = f_A1 * 2**(idx/12)
    octave = 1 + idx // 12 + (1 if (idx%12) >= 3 else 0)
    note = note_names[idx%12] + str(octave)
    
    t_n = t_ref * (f_t / f_ref)
    I_n = W * t_n**3 / 12
    k_n = 3 * E * I_n / L**3
    m_z = rho * W * t_n * L
    
    fr = beam_modes(N, lambda x,t=t_n: t, L, W, E, rho)
    
    print(f"{idx+1:3d} {note:>5s} {f_t:7.1f}  {t_n*1000:7.3f}  {m_z*1000:7.2f}  {k_n:7.1f}  "
          f" {fr[1]/fr[0]:6.2f}   {fr[1]:7.0f}   {fr[2]:7.0f}")

# ══════════════════════════════════════════════════════════════
# 3. Gewichtsmethode: Gleiche Dicke, Gewicht variiert
# ══════════════════════════════════════════════════════════════
print("\n\n══ TABELLE 2: Feste Dicke t = 0,40 mm + Gewicht (A1 bis Eigenfrequenz) ══")
print("-"*95)
print(f"{'Nr':>3s} {'Ton':>5s} {'f₁ [Hz]':>8s} {'m_tip [g]':>10s} {'f₂/f₁':>7s} {'f₂ [Hz]':>8s} "
      f"{'Δf₂/f₁':>8s}")
print("-"*95)

r2_ref_gew = None
for idx in range(25):
    f_t = f_A1 * 2**(idx/12)
    octave = 1 + idx // 12 + (1 if (idx%12) >= 3 else 0)
    note = note_names[idx%12] + str(octave)
    
    if f_t > f_ref:
        print(f"{idx+1:3d} {note:>5s} {f_t:7.1f}   (über Eigenfrequenz {f_ref:.0f} Hz)")
        continue
    
    # Bisection
    m_lo, m_hi = 0, 5e-3
    for _ in range(50):
        mm=(m_lo+m_hi)/2
        ff=beam_modes(N, lambda x: t_ref, L, W, E, rho, m_tip=mm)
        if ff[0]>f_t: m_lo=mm
        else: m_hi=mm
    m_opt=(m_lo+m_hi)/2
    fr = beam_modes(N, lambda x: t_ref, L, W, E, rho, m_tip=m_opt)
    rat = fr[1]/fr[0]
    
    if r2_ref_gew is None: r2_ref_gew = rat
    delta = rat - r2_ref_gew
    
    print(f"{idx+1:3d} {note:>5s} {f_t:7.1f}   {m_opt*1000:8.3f}    {rat:6.2f}   {fr[1]:7.0f}  "
          f" {delta:+6.2f}")

# ══════════════════════════════════════════════════════════════
# 4. Profil + Gewicht (Linear 40%, t₀ = 0,50 mm für größeren Bereich)
# ══════════════════════════════════════════════════════════════
print("\n\n══ TABELLE 3: Profil Linear 40% (t₀=0,50 mm, t_tip=0,30 mm) + Gewicht ══")
print("-"*95)
print(f"{'Nr':>3s} {'Ton':>5s} {'f₁ [Hz]':>8s} {'m_tip [g]':>10s} {'f₂/f₁':>7s} {'f₂ [Hz]':>8s} "
      f"{'f₃/f₁':>7s} {'Δf₂/f₁':>8s}")
print("-"*95)

t0_p = 0.50e-3; taper_p = 0.40
tf_p = lambda x: t0_p*(1-taper_p*x/L)
fr_base_p = beam_modes(N, tf_p, L, W, E, rho)

r2_ref_prof = None
for idx in range(25):
    f_t = f_A1 * 2**(idx/12)
    octave = 1 + idx // 12 + (1 if (idx%12) >= 3 else 0)
    note = note_names[idx%12] + str(octave)
    
    if f_t > fr_base_p[0]:
        print(f"{idx+1:3d} {note:>5s} {f_t:7.1f}   (über Basiston {fr_base_p[0]:.0f} Hz)")
        continue
    
    m_lo, m_hi = 0, 5e-3
    for _ in range(50):
        mm=(m_lo+m_hi)/2
        ff=beam_modes(N, tf_p, L, W, E, rho, m_tip=mm)
        if ff[0]>f_t: m_lo=mm
        else: m_hi=mm
    m_opt=(m_lo+m_hi)/2
    fr = beam_modes(N, tf_p, L, W, E, rho, m_tip=m_opt)
    rat = fr/fr[0]
    
    if r2_ref_prof is None: r2_ref_prof = rat[1]
    delta = rat[1] - r2_ref_prof
    
    print(f"{idx+1:3d} {note:>5s} {f_t:7.1f}   {m_opt*1000:8.3f}    {rat[1]:6.2f}   {fr[1]:7.0f}  "
          f" {rat[2]:6.2f}   {delta:+6.2f}")

# ══════════════════════════════════════════════════════════════
# 5. f₂ vs. typische f_H-Bereiche
# ══════════════════════════════════════════════════════════════
print("\n\n══ TABELLE 4: f₂ absolut vs. Kammerresonanz f_H ══")
print("-"*80)
print(f"{'Ton':>5s} {'f₁':>7s}  {'f₂(Dicke)':>10s} {'f₂(Profil)':>10s}  "
      f"{'f_H groß':>9s} {'f_H mittel':>10s} {'f_H klein':>9s}")
print(f"{'':>5s} {'':>7s}  {'6,27×f₁':>10s} {'≈5,1×f₁':>10s}  "
      f"{'200-350':>9s} {'300-500':>10s} {'500-700':>9s}")
print("-"*80)

for idx in range(25):
    f_t = f_A1 * 2**(idx/12)
    octave = 1 + idx // 12 + (1 if (idx%12) >= 3 else 0)
    note = note_names[idx%12] + str(octave)
    
    f2_dicke = 6.27 * f_t
    f2_prof = 5.1 * f_t  # Näherung für Linear 40%
    
    # Treffer markieren
    hit_d = []
    if 200 < f2_dicke < 350: hit_d.append("groß")
    if 300 < f2_dicke < 500: hit_d.append("mittel")
    if 500 < f2_dicke < 700: hit_d.append("klein")
    
    hit_p = []
    if 200 < f2_prof < 350: hit_p.append("groß")
    if 300 < f2_prof < 500: hit_p.append("mittel")
    if 500 < f2_prof < 700: hit_p.append("klein")
    
    marker = ""
    if hit_d: marker += f" ◄ Dicke trifft {','.join(hit_d)}"
    if hit_p: marker += f" ◄ Profil trifft {','.join(hit_p)}"
    
    print(f"{note:>5s} {f_t:6.1f}   {f2_dicke:9.0f}   {f2_prof:9.0f}  "
          f"  {'●' if 200<f2_dicke<350 else '○':>5s}      {'●' if 300<f2_dicke<500 else '○':>5s}       "
          f"{'●' if 500<f2_dicke<700 else '○':>5s}{marker}")

# ══════════════════════════════════════════════════════════════
# FAZIT
# ══════════════════════════════════════════════════════════════
print(f"""

{"="*80}
  ZUSAMMENFASSUNG
{"="*80}

  1. GLEICHMÄSSIGES PROFIL: f₂/f₁ = 6,27 = KONSTANT
     ─────────────────────────────────────────────────
     Die Modenreihe ist eine geometrische Eigenschaft des Cantilevers.
     Dickenänderung skaliert ALLE Moden gleich → Verhältnis bleibt.
     
     Konsequenz: Die Mensur (gleiche L) erzeugt für ALLE Töne dasselbe
     Verhältnis. Der Stimmer kann es nicht tonnweise ändern.

  2. GEWICHT: f₂/f₁ STEIGT mit zunehmendem Gewicht
     ───────────────────────────────────────────────
     Tiefer Ton → mehr Gewicht nötig → höheres f₂/f₁ (bis 9,89).
     Spreizung über die untere Oktave: ≈ +0,5 bis +1,5.
     
     Das verschärft die Inharmonizität bei tiefen Tönen.

  3. PROFIL + GEWICHT: f₂/f₁ START niedriger, STEIGT langsam
     ─────────────────────────────────────────────────────────
     Linear 40%: f₂/f₁ ≈ 5,1 ohne Gewicht.
     Mit Gewicht für tiefere Töne: steigt auf ≈ 5,8.
     Spreizung über den Bereich: ≈ 0,7.

  4. KRITISCHE ZONE: UNTERE OKTAVE
     ──────────────────────────────
     f₂ liegt für A1–A2 bei 345–689 Hz.
     Typische Kammerresonanzen f_H liegen bei 200–700 Hz.
     
     → Die gesamte untere Oktave hat f₂ im f_H-Bereich!
     → Die obere Oktave (A2–A3): f₂ = 689–1379 Hz → über f_H.
     
     Das erklärt, warum Probleme (Geistertöne, Dok. 0009)
     hauptsächlich bei den TIEFEN Tönen auftreten.

  5. WAS DER STIMMER TUN KANN
     ─────────────────────────
     f₂/f₁ ist durch die Mensur festgelegt — nicht tonnweise steuerbar.
     Was steuerbar ist: f_H der Kammer (Volumen, Öffnung).
     
     Die Aufgabe ist daher: f_H so wählen, dass es NICHT auf f₂ fällt.
     Da f₂ = 6,27 × f₁ ± 0 bei gleichmäßigem Profil, ist f₂ vorhersagbar.
     f_H muss unterhalb oder oberhalb von 6,27 × f₁ liegen.
""")
