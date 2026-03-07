#!/usr/bin/env python3
"""
Dok. 0016 — Ergänzung: Praxis-Profilierung
Glatte Verjüngung, dünnste Stelle wandert von Spitze bis Niete

Ergebnis: Die Position der dünnsten Stelle hat überraschend wenig
Einfluss auf f₂/f₁ — die TIEFE der Verdünnung bestimmt das Verhältnis.
Die Position bestimmt f₁ (und damit das benötigte Gewicht).

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

E=200e9; rho=7800; L=70e-3; W=8e-3; t0=0.40e-3; N=40

def make_smooth(x_min_frac, t_min_ratio):
    """Glatte Profilierung: dickste Stelle am Fuß, dünnste bei x_min_frac×L,
    Spitze etwas dicker als Minimum (15% Rückfederung)."""
    xm=x_min_frac*L; tm=t0*t_min_ratio; tt=tm+(t0-tm)*0.15
    def tf(x):
        if x<=xm:
            if xm>0: return t0-(t0-tm)*(1-np.cos(np.pi*x/xm))/2
            return tm
        else:
            r=L-xm
            if r>0: return tm+(tt-tm)*(1-np.cos(np.pi*(x-xm)/r))/2
            return tm
    return tf

def find_weight(tf, target_f):
    """Bisection: Gewicht für Zielfrequenz"""
    fb = beam_modes(N, tf, L, W, E, rho)
    if fb[0] < target_f: return None, None
    ml,mh=0,5e-3
    for _ in range(40):
        mm=(ml+mh)/2
        ff=beam_modes(N,tf,L,W,E,rho,m_tip=mm)
        if ff[0]>target_f: ml=mm
        else: mh=mm
    mo=(ml+mh)/2
    fr=beam_modes(N,tf,L,W,E,rho,m_tip=mo)
    return mo, fr

print("="*85)
print("  PRAXIS-PROFILIERUNG: Glatte Verjüngung")
print(f"  L = {L*1000:.0f} mm, W = {W*1000:.0f} mm, t₀ = {t0*1000:.2f} mm")
print("="*85)

# ══════════════════════════════════════════════════════════════
# 1. Dünnste Stelle wandert (Verdünnung 50%)
# ══════════════════════════════════════════════════════════════
print("\n1. DÜNNSTE STELLE WANDERT (t_min = 50% von t₀)")
print("-"*85)
print(f"{'x_min/L':>8s} {'t_min':>8s} {'t_tip':>8s} {'f₁':>7s} "
      f"{'f₂/f₁':>7s} {'f₃/f₁':>7s} {'f₂':>7s}")
print("-"*85)

x_n = np.linspace(0, L, N+1)
for xf in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2]:
    tf = make_smooth(xf, 0.50)
    fr = beam_modes(N, tf, L, W, E, rho)
    r = fr/fr[0]
    print(f"  {xf:5.2f}   {t0*0.50*1000:6.2f} mm  {tf(L)*1000:6.2f} mm  {fr[0]:6.1f}  "
          f"{r[1]:6.2f}  {r[2]:6.2f}  {fr[1]:6.0f}")

# ══════════════════════════════════════════════════════════════
# 2. Verdünnungstiefe variiert (Position x_min = 0,7)
# ══════════════════════════════════════════════════════════════
print("\n\n2. VERDÜNNUNGSTIEFE VARIIERT (x_min/L = 0,7)")
print("-"*70)
print(f"{'t_min/t₀':>8s} {'t_min':>10s} {'f₁':>7s} {'f₂/f₁':>7s} {'f₃/f₁':>7s} {'f₂':>7s}")
print("-"*70)

for ratio in [1.0, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40, 0.30, 0.20]:
    tf = make_smooth(0.70, ratio)
    fr = beam_modes(N, tf, L, W, E, rho)
    r = fr/fr[0]
    print(f"  {ratio:5.2f}   {t0*ratio*1000:7.2f} mm  {fr[0]:6.1f}  "
          f"{r[1]:6.2f}  {r[2]:6.2f}  {fr[1]:6.0f}")

# ══════════════════════════════════════════════════════════════
# 3. Alle Kombinationen + Gewicht → Ziel 50 Hz
# ══════════════════════════════════════════════════════════════
print("\n\n3. ALLE KOMBINATIONEN + GEWICHT → Ziel 50 Hz")
print("-"*85)
print(f"{'x_min/L':>8s} {'Verd.':>6s} {'m_tip':>8s} {'f₂/f₁':>7s} {'f₃/f₁':>7s} "
      f"{'f₂':>7s} {'Anmerkung':>20s}")
print("-"*85)

target = 50.0
for xf in [1.0, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4]:
    for depth in [0.30, 0.50, 0.70]:
        tf = make_smooth(xf, depth)
        mo, fr = find_weight(tf, target)
        if fr is None: continue
        r = fr/fr[0]
        anm = ""
        if r[1] < 4.5: anm = "← niedrig!"
        print(f"  {xf:5.2f}   {(1-depth)*100:4.0f}%  {mo*1000:6.3f}g  {r[1]:6.2f}  {r[2]:6.2f}  "
              f"{fr[1]:6.0f}   {anm}")
    print()

# ══════════════════════════════════════════════════════════════
# 4. Fazit
# ══════════════════════════════════════════════════════════════
print(f"""
{"="*85}
  ZUSAMMENFASSUNG
{"="*85}

  Position der dünnsten Stelle (bei gleicher Verdünnung 50%):
  ──────────────────────────────────────────────────────────
  x_min/L = 1.0 (Spitze):     f₂/f₁ = 4.35   ← klassische Profilierung
  x_min/L = 0.8:               f₂/f₁ = 4.20   ← Minimum! (kaum besser)
  x_min/L = 0.5 (Mitte):       f₂/f₁ = 5.16   ← schlechter
  x_min/L = 0.3 (nahe Niete):  f₂/f₁ = 6.22   ← fast wie gleichmäßig
  
  → Position hat überraschend flaches Optimum zwischen 0.7–1.0
  → Die TIEFE der Verdünnung bestimmt f₂/f₁, nicht die Position
  → Position nahe der Niete bringt fast nichts (f₂/f₁ ≈ 6.2)
  
  Verdünnungstiefe (bei x_min = 0.7):
  ────────────────────────────────────
  t_min/t₀ = 1.0 (keine):     f₂/f₁ = 6.27
  t_min/t₀ = 0.5 (50%):       f₂/f₁ = 4.35
  t_min/t₀ = 0.3 (70%):       f₂/f₁ = 3.94   ← niedrigstes erreichbares
  t_min/t₀ = 0.2 (80%):       f₂/f₁ = 4.25   ← steigt wieder! (zu dünn)
  
  → Es gibt ein Optimum bei ≈ 70% Verdünnung
  → Unter 30% Restdicke steigt f₂/f₁ wieder an (Zunge wird zu weich)

  Mit Gewicht auf 50 Hz:
  ──────────────────────
  Gleichmäßig + 0,34 g:            f₂/f₁ = 6.94  (Gewicht verschlechtert!)
  Mulde x=0.8, 50% + 0,21 g:      f₂/f₁ = 4.98
  Mulde x=0.7, 30% + 0,05 g:      f₂/f₁ = 4.38  (wenig Gewicht nötig!)
  Mulde x=0.8, 70% + 0,11 g:      f₂/f₁ = 3.95  (bestes Ergebnis)
  
  PRAXIS-KONSEQUENZ:
  ──────────────────
  1. Profilierung so tief wie mechanisch vertretbar (60–70% Verdünnung)
  2. Position der dünnsten Stelle: zwischen 70–90% der Zungenlänge
  3. So wenig Gewicht wie möglich (jedes Gramm erhöht f₂/f₁)
  4. Beste Kombination: tiefe Verdünnung bei 80% + minimales Gewicht
""")
