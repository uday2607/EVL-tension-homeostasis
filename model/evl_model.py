"""
Core ODE simulation for the EVL minimal vertex model.

Paper: Hino, Kapoor, Gubbala, Hannezo & Heisenberg
       "Apical domain mechanosensation regulates tissue tension homeostasis"
       Institute of Science and Technology Austria (ISTA)

Model summary
-------------
Each EVL cell is a prism with apical radius r(t) and height h(t).
Volume is conserved: V = V0 (fixed at t = 6 hpf), so h = V0 / (π r²).

Tissue tension:
    γ = 2 T_a − Γ_l · h / r

Apical contractility via Hill function:
    T_a = Σ · M^n / (M^n + M_H^n)

Equation of motion (overdamped):
    η · dr/dt = σ_ring(t) − γ

Mechanosensitive myosin dynamics:
    τ_M · dM/dt = M_H / (1 + β · ε⁺) − M
    ε = π r² / A₀ − 1   (fractional apical strain)
    ε⁺ = max(ε, 0)       (only tensile strain activates aPKC)
"""

import numpy as np
from experimental_data import V0, R0, A0, GAMMA_0

# ── Default parameters ────────────────────────────────────────────────────────

PARAMS = dict(
    ETA    = 14.0,     # lateral drag (nN·min/µm)
    KF     = 0.003,    # ring ramp rate (nN/µm/min)
    DELTA  = 0.40,     # ring offset above equilibrium for Control (nN/µm)
    BETA   = 12.0,     # aPKC mechanosensitivity
    TAU_M  = 90.0,     # mechanosensitive feedback timescale (min)
    N_HILL = 0.50,     # Hill exponent
    M_HIGH = 0.80,     # Myosin saturation level (model units)
    GL     = 1.50,     # lateral junction tension Γ_l (nN/µm)
    G0_WT  = GAMMA_0,  # WT tissue tension at t = 6 hpf (nN/µm)
    FRAC_DN = 0.65,    # DN-aPKC ring force as fraction of Control
)

T_START = 6.0    # simulation start (hpf)
T_END   = 9.0    # simulation end (hpf)
DT      = 0.5    # time step (min)
STEPS   = int((T_END - T_START) * 60 / DT)


def hill(M, M_HIGH, N_HILL):
    return M**N_HILL / (M**N_HILL + M_HIGH**N_HILL)


def simulate(sring0, SIG, p=None,
             is_dn=False, r_floor=R0):
    """
    Run the EVL model for one condition.

    Parameters
    ----------
    sring0  : float — initial ring force σ₀ (nN/µm)
    SIG     : float — junction tension gain Σ (sets T_a scale)
    p       : dict  — parameters (defaults to PARAMS if None)
    is_dn   : bool  — if True, M is frozen at M_HIGH (DN-aPKC: no feedback)
    r_floor : float — minimum allowed radius (µm); use R0 for normal epiboly

    Returns
    -------
    dict with keys 'hpf', 'r', 'M', 'gam'
        hpf : time array (hpf)
        r   : apical radius (µm)
        M   : Myosin II level (model units)
        gam : tissue tension γ (nN/µm)
    """
    if p is None:
        p = PARAMS

    ETA    = p['ETA'];   KF    = p['KF'];   BETA  = p['BETA']
    TAU_M  = p['TAU_M']; GL    = p['GL'];   M_HIGH = p['M_HIGH']
    N_HILL = p['N_HILL']

    r = R0; M = M_HIGH
    ts = []; rs = []; Ms = []; gs = []

    for i in range(STEPS + 1):
        t_min = i * DT
        t_hpf = T_START + t_min / 60.0
        M_    = M_HIGH if is_dn else float(np.clip(M, 0, M_HIGH))
        h     = V0 / (np.pi * r**2)
        gam   = 2 * SIG * hill(M_, M_HIGH, N_HILL) - GL * h / r

        ts.append(t_hpf); rs.append(r); Ms.append(M_); gs.append(gam)
        if i == STEPS:
            break

        def deriv(t2, r2, M2):
            r2  = max(r2, r_floor)
            M2  = float(np.clip(M2, 0, M_HIGH))
            ep2 = np.pi * r2**2 / A0 - 1
            h2  = V0 / (np.pi * r2**2)
            g2  = 2 * SIG * hill(M2, M_HIGH, N_HILL) - GL * h2 / r2
            sr2 = sring0 + KF * t2
            Meq = M_HIGH / (1 + BETA * max(ep2, 0))
            drdt = (sr2 - g2) / ETA
            dMdt = 0.0 if is_dn else (Meq - M2) / TAU_M
            return drdt, dMdt

        # 4th-order Runge-Kutta
        k1r, k1M = deriv(t_min, r, M)
        k2r, k2M = deriv(t_min + DT/2, r + DT/2*k1r, M + DT/2*k1M)
        k3r, k3M = deriv(t_min + DT/2, r + DT/2*k2r, M + DT/2*k2M)
        k4r, k4M = deriv(t_min + DT,   r + DT*k3r,   M + DT*k3M)

        r += DT/6 * (k1r + 2*k2r + 2*k3r + k4r)
        r  = max(r, r_floor)
        M += DT/6 * (k1M + 2*k2M + 2*k3M + k4M)
        M  = float(np.clip(M, 0, M_HIGH))

    return dict(hpf=np.array(ts), r=np.array(rs),
                M=np.array(Ms), gam=np.array(gs))


def build_conditions(p=None):
    """
    Return simulated outputs for all four experimental conditions.

    Returns dict with keys: 'control', 'ca', 'dn_area', 'dn_tens', 'dnca_tens'
        control   — Control (full mechanosensing, σ₀ = γ₀ + Δ)
        ca        — CA-Mypt YSL (σ₀ = γ₀, mechanosensing intact)
        dn_area   — DN-aPKC apical area (β=0, reduced ring force)
        dn_tens   — DN-aPKC tissue tension (β=0, WT ring force, elevated SIG)
        dnca_tens — DN-aPKC + CA-Mypt tension (β=0, σ₀=γ₀, additive SIG)
    """
    if p is None:
        p = PARAMS

    GL     = p['GL'];   G0_WT = p['G0_WT']
    DELTA  = p['DELTA']; FRAC_DN = p['FRAC_DN']

    h6   = V0 / (np.pi * R0**2)
    SIG     = G0_WT + GL * h6 / R0          # calibrated: γ(M_H, R0) = G0_WT
    SIG_DN  = 1.337 + GL * h6 / R0          # calibrated to DN tension at t=6
    SIG_DNCA = (SIG + SIG_DN) / 2           # additive combination

    sring0_wt = G0_WT + DELTA
    sring0_ca = G0_WT                        # ring at mechanical equilibrium
    sring0_dn = FRAC_DN * sring0_wt         # reduced ring force (EVL-only perturbation)

    return {
        'control'   : simulate(sring0_wt, SIG,     p=p),
        'ca'        : simulate(sring0_ca, SIG,     p=p),
        'dn_area'   : simulate(sring0_dn, SIG,     p=p, is_dn=True),
        'dn_tens'   : simulate(sring0_wt, SIG_DN,  p=p, is_dn=True),
        'dnca_tens' : simulate(sring0_ca, SIG_DNCA,p=p, is_dn=True),
    }


def at(sol, key, t_hpf):
    """Interpolate a solution variable at a given time."""
    return sol[key][np.argmin(np.abs(sol['hpf'] - t_hpf))]
