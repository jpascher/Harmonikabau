#!/usr/bin/env python3
"""
Dok. 0016 — Berechnungsskript
Obertonmoden der Basszunge: Profilierung, Gewichte, Kombination

FEM-Lösung (Euler-Bernoulli, Hermite-Elemente) für Cantilever-Balken
mit variablem Dickenprofil und optionalem Spitzengewicht.

Benötigt: numpy, scipy
"""
import numpy as np
from scipy.linalg import eigh

# ══════════════════════════════════════════════════════════════
# FEM-LÖSUNG
# ══════════════════════════════════════════════════════════════

def beam_modes(N, t_func, L, W, E, rho, m_tip=0, n_modes=6):
    """
    FEM für Cantilever-Balken mit variablem Dickenprofil t(x)
    und optionalem Spitzengewicht.

    Parameter:
        N        : Anzahl Elemente
        t_func   : Funktion t(x) → Dicke [m] an Position x [m]
        L        : Balkenlänge [m]
        W        : Balkenbreite [m]
        E        : E-Modul [Pa]
        rho      : Dichte [kg/m³]
        m_tip    : Spitzengewicht [kg] (default: 0)
        n_modes  : Anzahl berechneter Moden

    Rückgabe:
        freqs    : Eigenfrequenzen [Hz], Array der Länge n_modes
        modes    : Modenformen [N+1 × n_modes], normiert auf max=1
    """
    dx = L / N
    ndof = 2 * (N + 1)  # je 2 DOF pro Knoten (w, θ)

    K = np.zeros((ndof, ndof))
    M = np.zeros((ndof, ndof))

    for e in range(N):
        x_mid = (e + 0.5) * dx
        t_e = t_func(x_mid)
        I_e = W * t_e**3 / 12       # Flächenträgheitsmoment
        A_e = W * t_e               # Querschnittsfläche
        le = dx

        # Hermite-Balkenelemente (Steifigkeitsmatrix)
        k_e = E * I_e / le**3 * np.array([
            [ 12,    6*le,  -12,    6*le  ],
            [ 6*le,  4*le**2, -6*le, 2*le**2],
            [-12,   -6*le,   12,   -6*le  ],
            [ 6*le,  2*le**2, -6*le, 4*le**2]
        ])

        # Konsistente Massenmatrix
        m_e = rho * A_e * le / 420 * np.array([
            [ 156,    22*le,   54,   -13*le  ],
            [ 22*le,  4*le**2, 13*le, -3*le**2],
            [ 54,     13*le,  156,   -22*le  ],
            [-13*le, -3*le**2,-22*le,  4*le**2]
        ])

        dofs = [2*e, 2*e+1, 2*e+2, 2*e+3]
        for i in range(4):
            for j in range(4):
                K[dofs[i], dofs[j]] += k_e[i, j]
                M[dofs[i], dofs[j]] += m_e[i, j]

    # Spitzengewicht: Punktmasse am letzten Knoten
    if m_tip > 0:
        tip_dof = 2 * N   # w-DOF des Spitzenknotens
        M[tip_dof, tip_dof] += m_tip

    # Randbedingungen: Einspannung bei x=0 (w=0, θ=0)
    free_dofs = list(range(2, ndof))
    K_free = K[np.ix_(free_dofs, free_dofs)]
    M_free = M[np.ix_(free_dofs, free_dofs)]

    eigenvalues, eigenvectors = eigh(K_free, M_free)

    freqs = np.sqrt(np.abs(eigenvalues)) / (2 * np.pi)

    # Modenformen rekonstruieren
    modes = np.zeros((N + 1, n_modes))
    for m in range(min(n_modes, len(free_dofs))):
        full_mode = np.zeros(ndof)
        full_mode[free_dofs] = eigenvectors[:, m]
        modes[:, m] = full_mode[0::2]   # nur Verschiebungen
        mx = np.max(np.abs(modes[:, m]))
        if mx > 0:
            modes[:, m] /= mx

    return freqs[:n_modes], modes


# ══════════════════════════════════════════════════════════════
# PARAMETER
# ══════════════════════════════════════════════════════════════

E_stahl  = 200e9    # Pa    E-Modul Stahl
rho_stahl = 7800    # kg/m³ Dichte Stahl
L_zunge  = 70e-3    # m     Zungenlänge
W_zunge  = 8e-3     # m     Zungenbreite
t0       = 0.40e-3  # m     Dicke am Fuß

N_elem   = 80       # Anzahl FEM-Elemente
n_modes  = 6        # berechnete Moden

m_zunge  = rho_stahl * W_zunge * t0 * L_zunge

print("=" * 80)
print(f"  Dok. 0016 — Obertonmoden der Basszunge")
print(f"  L = {L_zunge*1000:.0f} mm, W = {W_zunge*1000:.0f} mm, "
      f"t₀ = {t0*1000:.2f} mm, E = {E_stahl/1e9:.0f} GPa, ρ = {rho_stahl} kg/m³")
print(f"  Zungenmasse: {m_zunge*1000:.2f} g, FEM: {N_elem} Elemente")
print("=" * 80)

# ══════════════════════════════════════════════════════════════
# 1. ANALYTISCHE REFERENZ (gleichmäßiger Cantilever)
# ══════════════════════════════════════════════════════════════
print("\n1. ANALYTISCHE REFERENZ (gleichmäßiger Cantilever)")
print("-" * 60)
beta_n = [1.8751, 4.6941, 7.8548, 10.9955, 14.1372, 17.2788]
ratios_analyt = [(b / beta_n[0])**2 for b in beta_n]
I0 = W_zunge * t0**3 / 12
A0 = W_zunge * t0
f1_analyt = beta_n[0]**2 / (2 * np.pi * L_zunge**2) * np.sqrt(E_stahl * I0 / (rho_stahl * A0))
print(f"  f₁ (analytisch) = {f1_analyt:.1f} Hz")
print(f"  Verhältnisse: {' : '.join([f'{r:.2f}' for r in ratios_analyt])}")
print(f"  Harmonisch:   {' : '.join([f'{n+1:.2f}' for n in range(6)])}")

# ══════════════════════════════════════════════════════════════
# 2. DICKENPROFILE (ohne Gewicht)
# ══════════════════════════════════════════════════════════════
print("\n\n2. DICKENPROFILE (ohne Gewicht)")
print("-" * 80)
print(f"{'Profil':35s} {'t_tip':>8s} {'Verd.':>6s} {'f₁':>8s} "
      f"{'f₂/f₁':>7s} {'f₃/f₁':>7s} {'f₄/f₁':>7s} {'f₅/f₁':>7s} {'f₆/f₁':>7s}")
print("-" * 80)

profiles = [
    ('Gleichmäßig (t = const)',    lambda x: t0),
    ('Linear 80 % (t_tip=0,32)',   lambda x: t0 * (1 - 0.20 * x / L_zunge)),
    ('Linear 50 % (t_tip=0,20)',   lambda x: t0 * (1 - 0.50 * x / L_zunge)),
    ('Linear 30 % (t_tip=0,12)',   lambda x: t0 * (1 - 0.70 * x / L_zunge)),
    ('Parabolisch 50 %',           lambda x: t0 * (1 - 0.50 * (x / L_zunge)**2)),
    ('Umgekehrt parabolisch',      lambda x: t0 * (1 - 0.50 * np.sqrt(x / L_zunge + 1e-12))),
]

for name, t_func in profiles:
    freqs, modes = beam_modes(N_elem, t_func, L_zunge, W_zunge, E_stahl, rho_stahl, n_modes=n_modes)
    ratios = freqs / freqs[0]
    t_tip = t_func(L_zunge) * 1000
    verd = (1 - t_tip / (t0 * 1000)) * 100
    print(f"{name:35s} {t_tip:7.2f}mm {verd:5.0f}% {freqs[0]:7.1f} Hz "
          f"{ratios[1]:6.2f}  {ratios[2]:6.2f}  {ratios[3]:6.2f}  {ratios[4]:6.2f}  {ratios[5]:6.2f}")

# ══════════════════════════════════════════════════════════════
# 3. STIMMGEWICHTE AN DER SPITZE (gleichmäßiges Profil)
# ══════════════════════════════════════════════════════════════
print("\n\n3. STIMMGEWICHTE AN DER SPITZE (gleichmäßiges Profil)")
print("-" * 80)
print(f"{'Gewicht':20s} {'m_tip':>8s} {'%m_Z':>6s} {'f₁':>8s} "
      f"{'f₂/f₁':>7s} {'f₃/f₁':>7s} {'f₄/f₁':>7s} {'f₅/f₁':>7s} {'f₆/f₁':>7s}")
print("-" * 80)

t_const = lambda x: t0
weight_configs = [
    ('Ohne Gewicht',           0),
    ('0,02 g (sehr leicht)',   0.02e-3),
    ('0,05 g (leicht)',        0.05e-3),
    ('0,10 g',                 0.10e-3),
    ('0,15 g',                 0.15e-3),
    ('0,20 g',                 0.20e-3),
    ('0,30 g',                 0.30e-3),
    ('0,40 g (≈ m_Z/4)',       0.40e-3),
    ('0,60 g',                 0.60e-3),
    ('0,80 g (≈ m_Z/2)',       0.80e-3),
    ('1,00 g',                 1.00e-3),
    ('1,50 g (≈ m_Z)',         1.50e-3),
]

for name, m_t in weight_configs:
    freqs, _ = beam_modes(N_elem, t_const, L_zunge, W_zunge, E_stahl, rho_stahl,
                          m_tip=m_t, n_modes=n_modes)
    ratios = freqs / freqs[0]
    pct = m_t / m_zunge * 100
    print(f"{name:20s} {m_t*1000:7.2f} g {pct:5.1f}% {freqs[0]:7.1f} Hz "
          f"{ratios[1]:6.2f}  {ratios[2]:6.2f}  {ratios[3]:6.2f}  {ratios[4]:6.2f}  {ratios[5]:6.2f}")

# ══════════════════════════════════════════════════════════════
# 4. KOMBINATION: Profilierung + Gewicht
# ══════════════════════════════════════════════════════════════
print("\n\n4. KOMBINATION: Profilierung + Gewicht")
print("-" * 80)
print(f"{'Konfiguration':40s} {'f₁':>7s} {'f₂/f₁':>7s} {'f₃/f₁':>7s} {'f₄/f₁':>7s}")
print("-" * 80)

L = L_zunge  # Kurzname für Lambdas
combos = [
    ('Gleichmäßig, ohne Gewicht',          lambda x: t0,                        0),
    ('Gleichmäßig + 0,10 g',               lambda x: t0,                        0.10e-3),
    ('Gleichmäßig + 0,20 g',               lambda x: t0,                        0.20e-3),
    ('Gleichmäßig + 0,40 g',               lambda x: t0,                        0.40e-3),
    ('',                                    None,                                 0),
    ('Linear 80 %, ohne Gewicht',           lambda x: t0*(1-0.20*x/L),          0),
    ('Linear 80 % + 0,10 g',               lambda x: t0*(1-0.20*x/L),          0.10e-3),
    ('Linear 80 % + 0,20 g',               lambda x: t0*(1-0.20*x/L),          0.20e-3),
    ('',                                    None,                                 0),
    ('Linear 50 %, ohne Gewicht',           lambda x: t0*(1-0.50*x/L),          0),
    ('Linear 50 % + 0,10 g',               lambda x: t0*(1-0.50*x/L),          0.10e-3),
    ('Linear 50 % + 0,20 g',               lambda x: t0*(1-0.50*x/L),          0.20e-3),
    ('Linear 50 % + 0,40 g',               lambda x: t0*(1-0.50*x/L),          0.40e-3),
    ('',                                    None,                                 0),
    ('Linear 30 %, ohne Gewicht',           lambda x: t0*(1-0.70*x/L),          0),
    ('Linear 30 % + 0,05 g',               lambda x: t0*(1-0.70*x/L),          0.05e-3),
    ('Linear 30 % + 0,10 g',               lambda x: t0*(1-0.70*x/L),          0.10e-3),
    ('Linear 30 % + 0,20 g',               lambda x: t0*(1-0.70*x/L),          0.20e-3),
    ('Linear 30 % + 0,40 g',               lambda x: t0*(1-0.70*x/L),          0.40e-3),
    ('',                                    None,                                 0),
    ('Parabolisch 50 %, ohne Gewicht',      lambda x: t0*(1-0.50*(x/L)**2),     0),
    ('Parabolisch 50 % + 0,10 g',          lambda x: t0*(1-0.50*(x/L)**2),     0.10e-3),
    ('Parabolisch 50 % + 0,20 g',          lambda x: t0*(1-0.50*(x/L)**2),     0.20e-3),
]

for name, tf, mt in combos:
    if tf is None:
        print()
        continue
    freqs, _ = beam_modes(N_elem, tf, L_zunge, W_zunge, E_stahl, rho_stahl,
                          m_tip=mt, n_modes=4)
    ratios = freqs / freqs[0]
    print(f"{name:40s} {freqs[0]:6.1f} Hz {ratios[1]:6.2f}  {ratios[2]:6.2f}  {ratios[3]:6.2f}")

# ══════════════════════════════════════════════════════════════
# 5. ZIELFREQUENZ 50 Hz: Welches Gewicht braucht jedes Profil?
# ══════════════════════════════════════════════════════════════
print("\n\n5. ZIELFREQUENZ 50 Hz: Benötigtes Gewicht pro Profil")
print("-" * 70)

target_profiles = [
    ('Gleichmäßig',           lambda x: t0),
    ('Linear 80 %',           lambda x: t0*(1-0.20*x/L)),
    ('Linear 50 %',           lambda x: t0*(1-0.50*x/L)),
    ('Linear 30 %',           lambda x: t0*(1-0.70*x/L)),
    ('Parabolisch 50 %',      lambda x: t0*(1-0.50*(x/L)**2)),
]

f_target = 50.0

print(f"{'Profil':25s} {'m_tip für 50Hz':>15s} {'f₂/f₁':>8s} {'f₃/f₁':>8s}")
print("-" * 60)

for name, tf in target_profiles:
    # Bisection: Gewicht finden, das f1 = 50 Hz ergibt
    m_lo, m_hi = 0, 2e-3
    for _ in range(50):
        m_mid = (m_lo + m_hi) / 2
        fr, _ = beam_modes(N_elem, tf, L_zunge, W_zunge, E_stahl, rho_stahl,
                           m_tip=m_mid, n_modes=3)
        if fr[0] > f_target:
            m_lo = m_mid
        else:
            m_hi = m_mid

    m_opt = (m_lo + m_hi) / 2
    fr_opt, _ = beam_modes(N_elem, tf, L_zunge, W_zunge, E_stahl, rho_stahl,
                           m_tip=m_opt, n_modes=4)
    rat = fr_opt / fr_opt[0]
    print(f"{name:25s} {m_opt*1000:12.3f} g  {rat[1]:7.2f}  {rat[2]:7.2f}")

# ══════════════════════════════════════════════════════════════
# 6. KNOTENPOSITIONEN (Mode 2, erste Nullstelle)
# ══════════════════════════════════════════════════════════════
print("\n\n6. KNOTENPOSITIONEN von Mode 2 (erste Nullstelle)")
print("-" * 50)

x_nodes = np.linspace(0, L_zunge, N_elem + 1)

check_profiles = [
    ('Gleichmäßig',       lambda x: t0,                  0),
    ('Linear 50 %',       lambda x: t0*(1-0.50*x/L),    0),
    ('Linear 30 %',       lambda x: t0*(1-0.70*x/L),    0),
    ('Gleichm. + 0,20 g', lambda x: t0,                  0.20e-3),
    ('Lin. 30% + 0,20 g', lambda x: t0*(1-0.70*x/L),    0.20e-3),
]

print(f"{'Konfiguration':25s} {'Knoten x/L':>10s} {'x [mm]':>8s}")
print("-" * 50)

for name, tf, mt in check_profiles:
    _, modes = beam_modes(N_elem, tf, L_zunge, W_zunge, E_stahl, rho_stahl,
                          m_tip=mt, n_modes=3)
    mode2 = modes[:, 1]
    # Nullstelle finden (Vorzeichenwechsel)
    for i in range(1, len(mode2)):
        if mode2[i-1] * mode2[i] < 0:
            # Lineare Interpolation
            x_zero = x_nodes[i-1] - mode2[i-1] * (x_nodes[i] - x_nodes[i-1]) / (mode2[i] - mode2[i-1])
            print(f"{name:25s} {x_zero/L_zunge:9.3f}   {x_zero*1000:7.1f}")
            break

# ══════════════════════════════════════════════════════════════
# 7. STEIFIGKEIT (statische Federrate k = F/x_tip)
# ══════════════════════════════════════════════════════════════
print("\n\n7. STEIFIGKEIT (effektive Federrate)")
print("-" * 60)

print(f"{'Konfiguration':30s} {'k [N/m]':>10s} {'k/k₀ [%]':>10s}")
print("-" * 55)

k0 = 3 * E_stahl * (W_zunge * t0**3 / 12) / L_zunge**3  # Referenz (gleichmäßig)

for name, tf in [('Gleichmäßig',lambda x:t0),('Linear 80%',lambda x:t0*(1-0.20*x/L)),
    ('Linear 50%',lambda x:t0*(1-0.50*x/L)),('Linear 30%',lambda x:t0*(1-0.70*x/L)),
    ('Parabolisch 50%',lambda x:t0*(1-0.50*(x/L)**2))]:
    # Numerische Steifigkeit: Einheitslast an Spitze
    dx = L_zunge / N_elem
    ndof = 2 * (N_elem + 1)
    K = np.zeros((ndof, ndof))
    for e in range(N_elem):
        xm = (e + 0.5) * dx
        te = tf(xm)
        Ie = W_zunge * te**3 / 12
        le = dx
        ke = E_stahl * Ie / le**3 * np.array([
            [12,6*le,-12,6*le],[6*le,4*le**2,-6*le,2*le**2],
            [-12,-6*le,12,-6*le],[6*le,2*le**2,-6*le,4*le**2]])
        d = [2*e, 2*e+1, 2*e+2, 2*e+3]
        for i in range(4):
            for j in range(4):
                K[d[i], d[j]] += ke[i, j]
    fd = list(range(2, ndof))
    K_free = K[np.ix_(fd, fd)]
    # Einheitskraft an Spitze (w-DOF)
    F = np.zeros(len(fd))
    tip_idx = 2 * N_elem - 2  # Index in freien DOFs
    F[tip_idx] = 1.0
    u = np.linalg.solve(K_free, F)
    w_tip = u[tip_idx]
    k_eff = 1.0 / w_tip
    print(f"{name:30s} {k_eff:9.1f}   {k_eff/k0*100:8.1f}%")

print("\n" + "=" * 80)
print("  Berechnung abgeschlossen.")
print("=" * 80)
