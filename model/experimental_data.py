"""
Experimental data used for model calibration and comparison.

Paper: Hino, Kapoor, Gubbala, Hannezo & Heisenberg
       "Apical domain mechanosensation regulates tissue tension homeostasis"
       Institute of Science and Technology Austria (ISTA)

All values are taken directly from the paper.
No raw data files are needed to run the model — everything is hardcoded here.

Units:
  Tension : nN/µm
  Area    : µm²  (normalized values are dimensionless)
  Volume  : µm³
  Radius  : µm
  Time    : hpf
"""

import numpy as np

# ── Tissue tension (micropipette aspiration, nN/µm) ──────────────────────────
# Normalized to Control at t = 6 hpf

TENSION = {
    "WT_6hpf":  1.000,   # reference (= γ₀ = 0.983 nN/µm absolute)
    "WT_9hpf":  0.983,   # stays approximately flat (homeostasis)
    "DN_6hpf":  1.337,   # elevated due to ~6 h of high Myosin prior to window
    "DN_9hpf":  2.030,   # rises substantially without mechanosensing
    "CA_9hpf":  0.890,   # slightly reduced
}

GAMMA_0 = 0.983          # absolute WT tension at 6 hpf (nN/µm), from MPA

# ── Apical area (normalized to WT at 6 hpf) ──────────────────────────────────

AREA = {
    "WT_6hpf":  1.000,
    "WT_9hpf":  2.220,   # ~2.2-fold expansion
    "DN_9hpf":  1.120,   # barely expands (high tension resists spreading)
    "CA_8hpf":  0.940,   # slightly smaller than WT (reduced ring force)
}

# ── Apical Myosin II (A.U., normalized to display scale 1-6) ─────────────────
# Measured as mean apical fluorescence intensity (Tg(actb2:Myl12.1-EGFP))

MYOSIN_AU = {
    "WT_6hpf":  6.32,    # absolute A.U. at 6 hpf
    "WT_9hpf":  2.20,    # drops ~3-fold during epiboly
}

# ── EVL cell volume (µm³, from Supp. Fig. 1c) ────────────────────────────────
# Same 10 cells tracked at 5, 7, 9 hpf (paired measurements)
# Mean values; volume decreases ~17% from 5 to 9 hpf

VOLUME_HPF = np.array([5.0, 7.0, 9.0])          # time points (hpf)
VOLUME_MEAN = np.array([7205.9, 6291.7, 5968.8]) # mean cell volume (µm³)

# Volume used in model: V0 = V at t = 6 hpf (interpolated), held constant
V0 = float(np.interp(6.0, VOLUME_HPF, VOLUME_MEAN))   # = 6748.8 µm³

# ── Initial cell geometry ─────────────────────────────────────────────────────

R0 = 16.0                          # initial apical radius (µm)
A0 = np.pi * R0**2                 # initial apical area (µm²) ≈ 804 µm²
H0 = V0 / (np.pi * R0**2)         # initial height (µm), from volume conservation


if __name__ == "__main__":
    print("=== Experimental data summary ===")
    print(f"  γ₀ (WT tension @ 6 hpf) = {GAMMA_0:.3f} nN/µm")
    print(f"  V₀ (cell volume @ 6 hpf) = {V0:.1f} µm³")
    print(f"  r₀ = {R0} µm,  h₀ = {H0:.2f} µm,  A₀ = {A0:.1f} µm²")
    print(f"  WT area expansion 6→9 hpf: {AREA['WT_9hpf']:.2f}x")
    print(f"  WT Myosin drop 6→9 hpf: {MYOSIN_AU['WT_6hpf']:.2f} → {MYOSIN_AU['WT_9hpf']:.2f} A.U.")
