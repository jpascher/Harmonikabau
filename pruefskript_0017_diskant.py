#!/usr/bin/env python3
"""
Dok. 0017 — Berechnungsskript
Diskant-Stimmzungen: Obertonmoden F3 (175 Hz) bis C6 (1047 Hz)

FEM-Lösung: Kalibrierte Mensur, Profilierung, Modenanalyse
Beweis: W hat keinen Einfluss auf die Tonhöhe

Benötigt: numpy, scipy
"""
import numpy as np
from scipy.linalg import eigh

def beam_modes(N, t_func, L, W, E, rho, n_modes=4):
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
    fd=list(range(2,ndof))
    ev,_=eigh(K[np.ix_(fd,fd)],M[np.ix_(fd,fd)])
    return np.sqrt(np.abs(ev[:n_modes]))/(2*np.pi)

# ══════════════════════════════════════════════════════════════
# Parameter
# ══════════════════════════════════════════════════════════════
E = 200e9       # Pa E-Modul Stahl
rho = 7800      # kg/m³ Dichte Stahl
N = 50          # FEM-Elemente
beta1 = 1.8751  # 1. Cantilever-Eigenwert

nN = ['C','C#','D','Eb','E','F','F#','G','G#','A','Bb','B']

# Bereich: F3 (175 Hz) bis C6 (1047 Hz) = 31 Halbtöne
f_F3 = 174.6
n_semi = 31  # 0=F3, 30=B5, 31=C6

# Mensur (Johann): F3: L=37mm, W=3.7mm; C6: L=14mm, W=1.7mm
def get_LW(semi):
    frac = semi / 31.0
    L = 37e-3 * (14.0/37.0)**frac
    W = 3.7e-3 * (1.7/3.7)**frac
    return L, W

def get_note(semi):
    abs_semi = 5 + semi  # F = Index 5 ab C
    return nN[abs_semi % 12] + str(3 + abs_semi // 12)

def make_profile(t0, verd, xmin_frac, L):
    """Glatte Profilierung: t0 am Fuß, Minimum bei xmin_frac×L,
    Spitze etwas dicker als Minimum (10% Rückfederung)"""
    xm = xmin_frac * L
    tm = t0 * verd
    tt = tm + (t0 - tm) * 0.10
    def tf(x):
        if x <= xm:
            if xm > 0:
                return t0 - (t0 - tm) * (1 - np.cos(np.pi * x / xm)) / 2
            return tm
        else:
            r = L - xm
            if r > 0:
                return tm + (tt - tm) * (1 - np.cos(np.pi * (x - xm) / r)) / 2
            return tm
    return tf

def calibrate_t0(f_target, L, W, verd, xmin_frac):
    """Bisection: t₀ finden so dass f₁(profiliert) = f_Ziel"""
    tl, th = 0.02e-3, 0.8e-3
    for _ in range(45):
        tm = (tl + th) / 2
        tf = make_profile(tm, verd, xmin_frac, L)
        fr = beam_modes(N, tf, L, W, E, rho, n_modes=2)
        if fr[0] > f_target:
            th = tm
        else:
            tl = tm
    return (tl + th) / 2


print("=" * 100)
print(f"  Dok. 0017 — Diskant-Stimmzungen: F3 ({f_F3:.0f} Hz) bis C6 (1047 Hz)")
print(f"  Mensur: L = 37→14 mm, W = 3,7→1,7 mm")
print(f"  FEM: {N} Elemente, Euler-Bernoulli, Hermite")
print("=" * 100)

# ══════════════════════════════════════════════════════════════
# 1. BEWEIS: W hat keinen Einfluss auf f₁
# ══════════════════════════════════════════════════════════════
print("\n\n1. BEWEIS: Breite W hat keinen Einfluss auf die Tonhöhe")
print("-" * 55)
L_test = 30e-3; t_test = 0.15e-3
print(f"   L = {L_test*1e3:.0f} mm, t = {t_test*1e3:.2f} mm, W variiert:\n")
print(f"   {'W [mm]':>8s} {'f₁ [Hz]':>10s} {'f₂ [Hz]':>10s} {'f₂/f₁':>8s}")
print("   " + "-" * 40)
for W_v in [1.0e-3, 2.0e-3, 3.0e-3, 5.0e-3, 8.0e-3]:
    fr = beam_modes(N, lambda x: t_test, L_test, W_v, E, rho)
    print(f"   {W_v*1e3:6.1f}   {fr[0]:9.2f}   {fr[1]:9.1f}   {fr[1]/fr[0]:7.3f}")
print(f"\n   → f₁ und f₂/f₁ IDENTISCH bei jeder Breite.")
print(f"   → f = β²/(2πL²) × t × √(E/(12ρ)) — W fällt raus.")
print(f"\n   W bestimmt:")
print(f"     Masse:            m = ρ × W × t × L        → ∝ W")
print(f"     Volumenstrom:     Q ∝ W × x_tip            → ∝ W (= Lautstärke)")
print(f"     Torsionssteifig.: ∝ W³                      → breitere = steifer")
print(f"     Spaltverhältnis:  Haupt/Nebenöffnung ∝ W   → breitere = weniger Kurzschluss")

# ══════════════════════════════════════════════════════════════
# 2. GLEICHMÄSSIG: t = const (Mindestdicke für Zielton)
# ══════════════════════════════════════════════════════════════
print("\n\n2. GLEICHMÄSSIG (t = const) — Mindestdicke für jeden Ton")
print("-" * 95)
print(f"{'Nr':>3s} {'Ton':>5s} {'f [Hz]':>7s} {'L [mm]':>7s} {'W [mm]':>7s} "
      f"{'t_min':>8s} {'m [mg]':>8s} {'k [N/m]':>8s} {'f₂/f₁':>7s} {'f₂ [Hz]':>8s}")
print("-" * 95)

for semi in range(n_semi + 1):
    f_t = f_F3 * 2**(semi/12)
    L_r, W_r = get_LW(semi)
    note = get_note(semi)
    
    # Analytisch: Mindestdicke für gleichmäßige Zunge
    t_min = f_t * 2*np.pi * L_r**2 / beta1**2 * np.sqrt(12*rho/E)
    
    fr = beam_modes(N, lambda x, t=t_min: t, L_r, W_r, E, rho, n_modes=4)
    rat = fr / fr[0]
    m_z = rho * W_r * t_min * L_r * 1e6  # mg
    I_r = W_r * t_min**3 / 12
    k_r = 3 * E * I_r / L_r**3
    
    if semi % 3 == 0 or semi % 12 == 0:
        print(f"{semi+1:3d} {note:>5s} {f_t:6.0f}  {L_r*1e3:6.2f}  {W_r*1e3:5.2f}  "
              f"{t_min*1e3:7.3f}  {m_z:7.1f}   {k_r:7.1f}  {rat[1]:6.2f}   {fr[1]:7.0f}")

print(f"\n   f₂/f₁ = 6,27 = KONSTANT — identisch mit Bass (Dok. 0016)")

# ══════════════════════════════════════════════════════════════
# 3. PROFILIERT: Kalibrierte Dicke am Fuß
# ══════════════════════════════════════════════════════════════
print("\n\n3. PROFILIERT — t₀ kalibriert per Bisection: f₁(profiliert) = f_Ziel")
print("-" * 105)
print(f"{'Nr':>3s} {'Ton':>5s} {'f [Hz]':>7s} {'L':>6s} {'W':>5s} {'t₀':>7s} "
      f"{'t_tip':>7s} {'Vd%':>4s} {'x_min':>5s} {'f₁':>7s} {'f₂/f₁':>7s} "
      f"{'f₂':>7s} {'f₃':>7s} {'k':>7s} {'m[mg]':>6s}")
print("-" * 105)

all_data = []
for semi in range(n_semi + 1):
    f_t = f_F3 * 2**(semi/12)
    L_r, W_r = get_LW(semi)
    note = get_note(semi)
    
    # Profilparameter: tiefe → stärker, Mitte; hohe → leichter, Spitze
    frac = semi / 31.0
    verd = 0.72 + 0.16 * frac   # t_min/t₀ = 0.72 → 0.88
    xmin = 0.55 + 0.40 * frac   # 0.55L → 0.95L
    
    t0_cal = calibrate_t0(f_t, L_r, W_r, verd, xmin)
    tf = make_profile(t0_cal, verd, xmin, L_r)
    fr = beam_modes(N, tf, L_r, W_r, E, rho, n_modes=4)
    rat = fr / fr[0]
    t_tip = tf(L_r)
    I_r = W_r * t0_cal**3 / 12
    k_r = 3 * E * I_r / L_r**3
    m_z = rho * W_r * t0_cal * L_r * 0.85 * 1e6  # mg (≈85% wegen Profilierung)
    
    all_data.append({
        'semi': semi, 'note': note, 'ft': f_t, 'L': L_r, 'W': W_r,
        't0': t0_cal, 'tt': t_tip, 'vd': verd, 'xm': xmin,
        'f1': fr[0], 'r2': rat[1], 'f2': fr[1], 'f3': fr[2],
        'k': k_r, 'm': m_z
    })
    
    if semi % 2 == 0 or semi % 12 == 0:
        print(f"{semi+1:3d} {note:>5s} {f_t:5.0f} {L_r*1e3:5.1f} {W_r*1e3:4.2f} "
              f"{t0_cal*1e3:6.3f} {t_tip*1e3:6.3f} {(1-verd)*100:3.0f}  {xmin:4.2f} "
              f"{fr[0]:6.1f}  {rat[1]:6.2f} {fr[1]:6.0f} {fr[2]:6.0f} "
              f"{k_r:6.1f} {m_z:5.0f}")

# ══════════════════════════════════════════════════════════════
# 4. f₂ vs. Diskantkammer f_H
# ══════════════════════════════════════════════════════════════
print("\n\n4. f₂ ABSOLUT vs. DISKANTKAMMER f_H")
print("   Diskantkammern: V ≈ 2–10 cm³, f_H ≈ 1500–5000 Hz")
print("-" * 80)
print(f"{'Ton':>5s} {'f₁':>7s} {'f₂':>7s} {'f₂/f₁':>7s}  {'Bewertung':>30s}")
print("-" * 80)

for d in all_data:
    hit = ""
    if d['f2'] < 1500:
        hit = "unter f_H → unkritisch"
    elif 1500 <= d['f2'] < 3000:
        hit = "◄ IM f_H-BEREICH"
    elif 3000 <= d['f2'] < 5000:
        hit = "◄ kleine Kammer"
    else:
        hit = "über f_H → unkritisch"
    
    if d['semi'] % 3 == 0 or d['semi'] % 12 == 0:
        print(f"  {d['note']:>5s} {d['ft']:6.0f}  {d['f2']:6.0f}  {d['r2']:6.2f}   {hit}")

# ══════════════════════════════════════════════════════════════
# 5. VERGLEICH Bass vs. Diskant
# ══════════════════════════════════════════════════════════════
print(f"""

{"=" * 80}
  VERGLEICH: BASS (Dok. 0016) vs. DISKANT
{"=" * 80}

                         Bass                    Diskant
  ──────────────────────────────────────────────────────────
  Bereich:               A1–A3 (55–220 Hz)      F3–C6 (175–1047 Hz)
  Zungenlänge:           70 mm (eine Mensur)     37→14 mm (variiert)
  Zungenbreite:          8 mm (fix)              3,7→1,7 mm (variiert)
  Dicke am Fuß:          0,33–1,32 mm            {all_data[0]['t0']*1e3:.2f}–{all_data[-1]['t0']*1e3:.2f} mm
  Profilierung:          30–50 %                 12–28 %
  
  f₂/f₁ (gleichmäßig):  6,27                    6,27  (identisch!)
  f₂/f₁ (profiliert):   4,0–5,3                 {all_data[0]['r2']:.1f}–{all_data[-1]['r2']:.1f}
  
  Kritische Zone:        A1–A2                   C4–F5
                         f₂ = 345–689 Hz         f₂ ≈ 1400–4000 Hz
                         → Basskammern            → Diskantkammern
  
  W-Einfluss auf Ton:    KEINER                  KEINER
  W bestimmt:            Lautstärke, Ansprache    Lautstärke, Ansprache

  GEMEINSAMKEIT:
  ─────────────
  f₂/f₁ = 6,27 (gleichmäßig) ist eine universelle Cantilever-Konstante.
  Profilierung senkt sie — im Bass stärker (mehr Verdünnung möglich).
  Die Harmonischen im Klang stammen vom Impulsgenerator, nicht von den Moden.
""")
