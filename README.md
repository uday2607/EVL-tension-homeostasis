# EVL Tension Homeostasis — Code Repository

Code accompanying the paper:

> **Apical domain mechanosensation regulates tissue tension homeostasis**
> Naoya Hino, Tushna Kapoor, Uday Ram Gubbala, Edouard Hannezo, and Carl-Philipp Heisenberg
> Institute of Science and Technology Austria (ISTA)

---

## Overview

During zebrafish epiboly, the enveloping layer (EVL) undergoes ~2.2-fold apical expansion while maintaining nearly constant tissue tension. This repository contains:

1. **Minimal vertex model** (`model/`) — Python implementation of the mechanosensitive EVL cell model, including all simulation conditions, figure generation, and sensitivity analysis.
2. **Image analysis scripts** (`image_analysis/`) — MATLAB scripts used for EVL apical surface projection and Actin–Kibra spatial cross-correlation analysis.

---

## Repository Structure

```
EVL-tension-homeostasis/
├── README.md                    # this file
├── requirements.txt             # Python dependencies
├── model/
│   ├── README.md                # detailed model description
│   ├── experimental_data.py     # hardcoded experimental values (no raw data needed)
│   ├── evl_model.py             # core ODE simulation
│   ├── run_model.py             # run all conditions and print summary
│   └── plot_figures.py          # reproduce all model figures (Fig 2 o-r, Supp Fig 3)
└── image_analysis/
    ├── README.md                # usage instructions
    ├── EVL_surface_projection.m # apical surface projection from z-stack (MATLAB)
    └── Actin_kibra_xcorr.m      # Actin-Kibra spatial cross-correlation (MATLAB)
```

---

## Quick Start (Model)

```bash
# Install dependencies
pip install -r requirements.txt

# Run the model for all conditions and print key outputs
python model/run_model.py

# Reproduce all model figures
python model/plot_figures.py
```

Figures are saved to `figures/`.

---

## Requirements

- Python >= 3.9
- See `requirements.txt` for full list

MATLAB scripts require:
- MATLAB R2020b or later
- Image Processing Toolbox
- Signal Processing Toolbox (for `Actin_kibra_xcorr.m`)
