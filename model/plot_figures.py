"""
Generate all model figures for the paper.

Paper: Hino, Kapoor, Gubbala, Hannezo & Heisenberg
       "Apical domain mechanosensation regulates tissue tension homeostasis"
       Institute of Science and Technology Austria (ISTA)

Figures produced
----------------
  fig2p_apical_area.svg      — Fig 2p: apical area vs time
  fig2q_tissue_tension.svg   — Fig 2q: tissue tension vs time
  fig2r_apical_myosin.svg    — Fig 2r: apical Myosin II vs time
  suppfig3b_sensitivity.svg  — Supp Fig 3b: tornado sensitivity chart
  suppfig1d_vol_sensitivity.svg — Supp Fig 1d: volume assumption sensitivity

Usage
-----
    python plot_figures.py

Outputs saved to ../../figures/
"""

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os, sys

matplotlib.rcParams['svg.fonttype'] = 'none'
matplotlib.rcParams['font.family'] = 'sans-serif'
matplotlib.rcParams['font.sans-serif'] = ['Arial', 'Helvetica', 'DejaVu Sans']

sys.path.insert(0, os.path.dirname(__file__))
from evl_model import build_conditions, at, simulate, PARAMS, A0, STEPS, DT, T_START
from experimental_data import GAMMA_0, V0, R0

OUT = os.path.join(os.path.dirname(__file__), '..', 'figures')
os.makedirs(OUT, exist_ok=True)

G_REF  = GAMMA_0
A_REF  = A0
M_HIGH = PARAMS['M_HIGH']
M_JUNC = 1.0; M_APIC = 5.0
M_REF  = M_JUNC + M_APIC

# Color palette (Tableau Dark2 — no pure RGB)
COL = dict(
    control = '#B07AA1',   # mauve
    ca      = '#76B7B2',   # slate teal
    dn      = '#F28E2B',   # vivid orange
    dnca    = '#EDC948',   # gold
)
LS = dict(control='-', ca='--', dn='-.', dnca=':')
LW = 1.3; FS = 7; FSL = 6


def style_ax(ax, ylabel, ylim, yticks, legend_loc='best'):
    ax.set_xlabel('Time (hpf)', fontsize=FS)
    ax.set_ylabel(ylabel, fontsize=FS)
    ax.set_xlim(6, 9); ax.set_xticks([6, 7, 8, 9])
    ax.set_ylim(*ylim); ax.set_yticks(yticks)
    ax.tick_params(labelsize=FS, length=2, pad=1.5, width=0.7)
    for sp in ['top', 'right']:
        ax.spines[sp].set_visible(False)
    ax.spines['left'].set_linewidth(0.7)
    ax.spines['bottom'].set_linewidth(0.7)
    ax.legend(fontsize=FSL, frameon=False, loc=legend_loc,
              handlelength=1.5, handletextpad=0.4, labelspacing=0.28)


# ── Build all conditions ──────────────────────────────────────────────────────
conds = build_conditions()
wt   = conds['control']
ca   = conds['ca']
dn   = conds['dn_area']
dt_  = conds['dn_tens']
dnca = conds['dnca_tens']

mask = wt['hpf'] >= 6


# ── Fig 2p: Apical area ───────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(1.55, 1.70))
for sol, lbl, key in [(dn, 'DN-aPKC', 'dn'),
                       (wt, 'Control', 'control'),
                       (ca, 'CA-Mypt (YSL)', 'ca')]:
    ax.plot(sol['hpf'][mask], np.pi*sol['r'][mask]**2/A_REF,
            color=COL[key], lw=LW, ls=LS[key], label=lbl)
style_ax(ax, 'Apical area (rel. to WT, 6 hpf)',
         (0, 3.2), [0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0], 'upper left')
plt.tight_layout(pad=0.4)
plt.savefig(os.path.join(OUT, 'fig2p_apical_area.svg'), bbox_inches='tight')
plt.close()

# ── Fig 2q: Tissue tension ────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(1.55, 1.70))
for sol, lbl, key in [(dt_,  'DN-aPKC',                 'dn'),
                       (dnca, 'DN-aPKC + CA-Mypt (YSL)', 'dnca'),
                       (wt,   'Control',                  'control'),
                       (ca,   'CA-Mypt (YSL)',            'ca')]:
    ax.plot(sol['hpf'][mask], sol['gam'][mask]/G_REF,
            color=COL[key], lw=LW, ls=LS[key], label=lbl)
style_ax(ax, 'Tissue tension (rel. to WT, 6 hpf)',
         (0.5, 2.2), [0.5, 1.0, 1.5, 2.0], 'upper left')
plt.tight_layout(pad=0.4)
plt.savefig(os.path.join(OUT, 'fig2q_tissue_tension.svg'), bbox_inches='tight')
plt.close()

# ── Fig 2r: Apical Myosin II ──────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(1.55, 1.70))
for sol, lbl, key in [(dn, 'DN-aPKC', 'dn'),
                       (wt, 'Control', 'control'),
                       (ca, 'CA-Mypt (YSL)', 'ca')]:
    yy = (M_JUNC + sol['M'][mask]/M_HIGH*M_APIC) / M_REF
    ax.plot(sol['hpf'][mask], yy, color=COL[key], lw=LW, ls=LS[key], label=lbl)
style_ax(ax, 'Apical myosin (rel. to WT, 6 hpf)',
         (0, 1.5), [0, 0.5, 1.0, 1.5], 'lower left')
plt.tight_layout(pad=0.4)
plt.savefig(os.path.join(OUT, 'fig2r_apical_myosin.svg'), bbox_inches='tight')
plt.close()

print("Saved Fig 2p, 2q, 2r.")


# ── Supp Fig 3b: Sensitivity tornado ─────────────────────────────────────────
def run_wt(p):
    GL = p['GL']; G0_WT = p['G0_WT']; DELTA = p['DELTA']
    h6 = V0 / (np.pi * R0**2)
    SIG = G0_WT + GL * h6 / R0
    return simulate(G0_WT + DELTA, SIG, p=p)

def at_t(sol, key, t):
    return sol[key][np.argmin(np.abs(sol['hpf'] - t))]

nom = run_wt(PARAMS)
nom_g9 = at_t(nom, 'gam', 9)
nom_m9 = M_JUNC + at_t(nom, 'M', 9) / M_HIGH * M_APIC
nom_a9 = np.pi * at_t(nom, 'r', 9)**2 / 1e3

PARAM_KEYS  = ['BETA', 'KF', 'DELTA', 'TAU_M', 'ETA']
PARAM_LABELS = [r'$\beta$ (aPKC sens.)', r'$k_f$ (ring ramp)',
                r'$\Delta$ (ring offset)', r'$\tau_M$ (myosin $\tau$)', r'$\eta$ (drag)']
OUT_LABELS  = ['Tissue tension@9 hpf', 'Apical myosin@9 hpf', 'Apical area@9 hpf']
OUT_COLORS  = ['#d62728', '#2ca02c', '#1f77b4']
PERT = 0.30

sens_pos = np.zeros((3, len(PARAM_KEYS)))
sens_neg = np.zeros((3, len(PARAM_KEYS)))
noms9 = [nom_g9, nom_m9, nom_a9]

for ci, pkey in enumerate(PARAM_KEYS):
    for sign, arr in [(+PERT, sens_pos), (-PERT, sens_neg)]:
        p = dict(PARAMS); p[pkey] = PARAMS[pkey] * (1 + sign)
        s = run_wt(p)
        vals = [at_t(s, 'gam', 9),
                M_JUNC + at_t(s, 'M', 9)/M_HIGH*M_APIC,
                np.pi * at_t(s, 'r', 9)**2 / 1e3]
        for ri in range(3):
            arr[ri, ci] = 100 * (vals[ri] - noms9[ri]) / noms9[ri]

fig, axes = plt.subplots(1, 3, figsize=(7.0, 2.2))
fig.subplots_adjust(wspace=0.18, left=0.23, right=0.97, top=0.85, bottom=0.14)
y = np.arange(len(PARAM_KEYS)); bh = 0.30

for col, (ax, out_lbl, oc) in enumerate(zip(axes, OUT_LABELS, OUT_COLORS)):
    pos = sens_pos[col]; neg = sens_neg[col]
    ax.barh(y+bh/2, pos, height=bh, color=oc, alpha=0.90, label=f'+{int(PERT*100)}%')
    ax.barh(y-bh/2, neg, height=bh, color=oc, alpha=0.38, label=f'-{int(PERT*100)}%')
    ax.axvline(0, color='k', lw=0.7)
    ax.set_yticks(y)
    ax.set_yticklabels(PARAM_LABELS if col == 0 else [], fontsize=FSL)
    ax.set_xlabel('% change in output', fontsize=FSL)
    ax.set_title(out_lbl, fontsize=FSL, pad=3, color=oc, fontweight='bold')
    ax.tick_params(labelsize=FSL, length=2, width=0.6, pad=1.5)
    for sp in ['top', 'right']: ax.spines[sp].set_visible(False)
    xmax = max(np.max(np.abs(pos)), np.max(np.abs(neg))) * 1.35
    ax.set_xlim(-xmax, xmax)
    thresh = xmax * 0.06
    for i, (p_, n_) in enumerate(zip(pos, neg)):
        pad_ = xmax * 0.03
        if abs(p_) > thresh:
            ax.text(p_+np.sign(p_)*pad_, i+bh/2, f'{p_:+.1f}%',
                    va='center', ha='left' if p_>0 else 'right', fontsize=4.8)
        if abs(n_) > thresh:
            ax.text(n_+np.sign(n_)*pad_, i-bh/2, f'{n_:+.1f}%',
                    va='center', ha='left' if n_>0 else 'right', fontsize=4.8)

axes[0].legend(fontsize=5.5, frameon=False, loc='lower right',
               handlelength=1.2, handletextpad=0.3, labelspacing=0.2)
plt.savefig(os.path.join(OUT, 'suppfig3b_sensitivity.svg'), bbox_inches='tight')
plt.close()
print("Saved Supp Fig 3b (sensitivity tornado).")


# ── Supp Fig 1d: Volume assumption sensitivity ────────────────────────────────
def V_const(t):  return V0
def V_dec30(t):  return V0 * (1 - 0.30 * (np.clip(t, 6, 9) - 6) / 3)
def V_inc30(t):  return V0 * (1 + 0.30 * (np.clip(t, 6, 9) - 6) / 3)

def sim_vol(V_func):
    p = PARAMS
    GL = p['GL']; G0_WT = p['G0_WT']; DELTA = p['DELTA']
    h6 = V_func(6.0) / (np.pi * R0**2)
    SIG = G0_WT + GL * h6 / R0
    sring0 = G0_WT + DELTA
    r = R0; M = M_HIGH; ts = []; gs = []; rs = []
    for i in range(STEPS + 1):
        t_min = i * DT; t_hpf = T_START + t_min / 60.0
        M_ = float(np.clip(M, 0, M_HIGH))
        V  = V_func(t_hpf); h = V / (np.pi * r**2)
        gam = 2*SIG*((M_**p['N_HILL'])/(M_**p['N_HILL']+M_HIGH**p['N_HILL'])) - GL*h/r
        ts.append(t_hpf); gs.append(gam); rs.append(r)
        if i == STEPS: break
        def deriv(t2, r2, M2):
            r2 = max(r2, R0); M2 = float(np.clip(M2, 0, M_HIGH))
            ep2 = np.pi*r2**2/A0 - 1
            V2 = V_func(T_START + t2/60); h2 = V2/(np.pi*r2**2)
            g2 = 2*SIG*((M2**p['N_HILL'])/(M2**p['N_HILL']+M_HIGH**p['N_HILL'])) - GL*h2/r2
            Meq = M_HIGH / (1 + p['BETA']*max(ep2, 0))
            return (sring0 + p['KF']*t2 - g2)/p['ETA'], (Meq - M2)/p['TAU_M']
        k1r,k1M=deriv(t_min,r,M); k2r,k2M=deriv(t_min+DT/2,r+DT/2*k1r,M+DT/2*k1M)
        k3r,k3M=deriv(t_min+DT/2,r+DT/2*k2r,M+DT/2*k2M); k4r,k4M=deriv(t_min+DT,r+DT*k3r,M+DT*k3M)
        r += DT/6*(k1r+2*k2r+2*k3r+k4r); r = max(r, R0)
        M += DT/6*(k1M+2*k2M+2*k3M+k4M); M = float(np.clip(M, 0, M_HIGH))
    return dict(hpf=np.array(ts), gam=np.array(gs), r=np.array(rs))

runs   = [sim_vol(V_const), sim_vol(V_dec30), sim_vol(V_inc30)]
labels = ['V = const', 'V -30%', 'V +30%']
colors = ['#B07AA1', '#F28E2B', '#76B7B2']
lstyle = ['-', '--', '--']

H_IN = 120/72
fig, axes = plt.subplots(1, 2, figsize=(3.2, H_IN))
fig.patch.set_alpha(0)
specs = [
    ('Tissue tension (rel. to WT, 6 hpf)',
     [s['gam']/G_REF for s in runs], (0.5, 2.0), [0.5, 1.0, 1.5, 2.0]),
    ('Apical area (rel. to WT, 6 hpf)',
     [np.pi*s['r']**2/A_REF for s in runs], (0, 3.2), [0, 1.0, 2.0, 3.0]),
]
for ax, (ylabel, yarrs, ylim, yticks) in zip(axes, specs):
    for yy, lbl, col, ls in zip(yarrs, labels, colors, lstyle):
        ax.plot(runs[0]['hpf'][mask], yy[mask], color=col, lw=LW, ls=ls, label=lbl)
    ax.set_xlabel('Time (hpf)', fontsize=FS); ax.set_ylabel(ylabel, fontsize=FS)
    ax.set_xlim(6, 9); ax.set_xticks([6, 7, 8, 9])
    ax.set_ylim(*ylim); ax.set_yticks(yticks)
    ax.tick_params(labelsize=FS, length=2, pad=1.5, width=0.7)
    for sp in ['top','right']: ax.spines[sp].set_visible(False)
    ax.set_facecolor('none')
axes[0].legend(fontsize=FSL, frameon=False, loc='upper left',
               handlelength=1.5, handletextpad=0.4, labelspacing=0.28)
plt.tight_layout(pad=0.4)
plt.savefig(os.path.join(OUT, 'suppfig1d_vol_sensitivity.svg'),
            bbox_inches='tight', transparent=True)
plt.close()
print("Saved Supp Fig 1d (volume sensitivity).")
print(f"\nAll figures saved to: {os.path.abspath(OUT)}")
