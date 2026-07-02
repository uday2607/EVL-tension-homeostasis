"""
Run all four EVL model conditions and compare outputs to experimental targets.

Paper: Hino, Kapoor, Gubbala, Hannezo & Heisenberg
       "Apical domain mechanosensation regulates tissue tension homeostasis"
       Institute of Science and Technology Austria (ISTA)

Usage
-----
    python run_model.py
"""

import numpy as np
from evl_model import build_conditions, at, PARAMS, A0
from experimental_data import GAMMA_0

G_REF = GAMMA_0
A_REF = A0
M_HIGH = PARAMS['M_HIGH']
M_JUNC = 1.0; M_APIC = 5.0   # display scale: M=0 → 1.0, M=M_HIGH → 6.0


def norm_area(sol, t):
    return np.pi * at(sol, 'r', t)**2 / A_REF

def norm_tens(sol, t):
    return at(sol, 'gam', t) / G_REF

def norm_myo(sol, t):
    return (M_JUNC + at(sol, 'M', t) / M_HIGH * M_APIC) / (M_JUNC + M_APIC)


if __name__ == "__main__":
    conds = build_conditions()
    wt  = conds['control']
    ca  = conds['ca']
    da  = conds['dn_area']
    dt_ = conds['dn_tens']

    print("=" * 60)
    print("EVL Minimal Vertex Model — Output Summary")
    print("=" * 60)
    print(f"\n{'Quantity':<30} {'Model':>8}  {'Experiment':>12}")
    print("-" * 55)
    print(f"{'WT area @ 9 hpf':<30} {norm_area(wt, 9):>8.3f}  {'2.220':>12}")
    print(f"{'DN-aPKC area @ 9 hpf':<30} {norm_area(da, 9):>8.3f}  {'1.120':>12}")
    print(f"{'CA-Mypt area @ 8 hpf':<30} {norm_area(ca, 8):>8.3f}  {'0.940':>12}")
    print(f"{'WT tension @ 9 hpf':<30} {norm_tens(wt, 9):>8.3f}  {'0.983':>12}")
    print(f"{'DN-aPKC tension @ 9 hpf':<30} {norm_tens(dt_, 9):>8.3f}  {'2.030':>12}")
    print(f"{'WT myosin @ 9 hpf (norm)':<30} {norm_myo(wt, 9):>8.3f}  {'0.350':>12}")
    print("-" * 55)
    print("\nAll values normalized to Control at t = 6 hpf.")
    print(f"Parameters: {PARAMS}")
