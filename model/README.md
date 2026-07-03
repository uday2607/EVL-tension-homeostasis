# Minimal Vertex Model — EVL Tension Homeostasis

## Biological context

During zebrafish epiboly (6-9 hpf), the EVL is pulled vegetally by the YSL actomyosin ring. Without regulation, this stretching would passively increase tissue tension. Instead, tension is observed to remain approximately constant — a phenomenon we term **tension homeostasis**. The model shows this arises from a mechanosensitive negative-feedback circuit: stretch activates aPKC, which downregulates apical Myosin II, reducing contractility and counteracting the passive tension increase.

---

## Model description

### Cell geometry

Each EVL cell is represented as a prism with:
- Apical radius `r(t)` (evolves in time)
- Height `h(t)`
- **Volume conservation**: V = V(t=6 hpf) = constant throughout the simulation window (6-9 hpf). This gives `h = V / (π r²)`.

### Tissue tension

From the cell energy and the Young-Laplace relation:

```
γ(t) = 2 T_a(t) − Γ_l · h(t) / r(t)
```

- `T_a`: active apical junction tension (set by Myosin II)
- `Γ_l = 1.50 nN/µm`: passive lateral junction tension
- The geometric term `Γ_l · h/r` decreases as cells spread, partially compensating the tension change

### Equation of motion

In the overdamped limit:

```
η · dr/dt = σ_ring(t) − γ(t)
```

The YSL ring force ramps linearly:

```
σ_ring(t) = σ₀ + k_f · t        (t in minutes from 6 hpf)
```

### Myosin II and contractility

Apical tension scales with Myosin II via a Hill function:

```
T_a = Σ · M^n / (M^n + M_H^n)        n = 0.5
```

The sublinear (n < 1) response means a moderate myosin drop gives substantial tension relief, making homeostasis achievable at realistic parameter values.

### Mechanosensitive feedback

Myosin II relaxes toward a strain-dependent equilibrium:

```
τ_M · dM/dt = M_H / (1 + β · ε⁺) − M
```

- `ε = π r² / A₀ − 1`: fractional apical strain
- `ε⁺ = max(ε, 0)`: only tensile strain activates aPKC
- `β`: aPKC mechanosensitivity
- `τ_M`: timescale of mechanosensitive cascade (~90 min)

---

## Experimental conditions → model parameters

| Condition | Biology | Model setting |
|---|---|---|
| **Control** | Normal EVL + YSL | σ₀ = γ₀ + Δ, β = 12, M evolves |
| **CA-Mypt (YSL)** | Reduced YSL ring force | σ₀ = γ₀ (ring at equilibrium), β = 12, M evolves |
| **DN-aPKC** | aPKC blocked in EVL | σ₀ = γ₀ + Δ (YSL unchanged), β = 0, M frozen at M_H |
| **DN-aPKC + CA-Mypt (YSL)** | Both perturbations | σ₀ = γ₀, β = 0, M frozen at M_H |

For DN-aPKC: Myosin II is elevated for ~6 h prior to the model window, so Σ_DN is recalibrated to match the elevated tension measured at t = 6 hpf. The YSL ring force is unchanged (DN-aPKC is EVL-only).

---

## Parameters

| Symbol | Value | Description | How determined |
|---|---|---|---|
| η | 14.0 nN·min/µm | Lateral drag | Calibrated to Control area timecourse |
| k_f | 0.003 nN/µm/min | Ring ramp rate | Calibrated to Control area at t = 9 hpf |
| Δ | 0.40 nN/µm | Ring offset (Control) | Calibrated to Control spreading rate |
| β | 12.0 | aPKC mechanosensitivity | Calibrated to Control Myosin drop |
| τ_M | 90 min | Feedback timescale | Calibrated to Myosin decrease timescale |
| n | 0.50 | Hill exponent | Chosen for sublinear Myosin-tension coupling |
| M_H | 0.80 | Myosin saturation level | Sets Myosin scale (model units) |
| Γ_l | 1.50 nN/µm | Lateral junction tension | Calibrated from geometry and tension at t = 6 hpf |
| γ₀ | 0.983 nN/µm | EVL tension at t = 6 hpf | Directly from micropipette aspiration |
| r₀ | 16 µm | Initial cell radius | From experimental data |

---

## Files

| File | Description |
|---|---|
| `experimental_data.py` | All hardcoded experimental values used for calibration and comparison. No raw data files needed. |
| `evl_model.py` | Core ODE solver. Returns time series of r(t), M(t), γ(t) for any condition. |
| `run_model.py` | Runs all four conditions and prints key outputs vs experimental targets. |
| `plot_figures.py` | Generates all model figures: apical area (Fig 2p), tissue tension (Fig 2q), apical myosin (Fig 2r), sensitivity tornado (Supp Fig 3b), volume sensitivity (Supp Fig 1d). |
