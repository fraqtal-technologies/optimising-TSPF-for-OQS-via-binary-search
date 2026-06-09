# Optimising Trotter-Suzuki Product Formulas for Open Quantum Systems via Binary Search

Code accompanying the manuscript:

> **Optimizing Gate Complexities Using Classical Search**  
> — numerical reproduction of the Trotter-step and gate-complexity bounds for open
> quantum system simulation using deterministic and randomised Trotter-Suzuki product
> formulas (TS-PF).

## Overview

This repository contains the complete numerical workflow for two benchmark models:

- **XX-spin chain** with boundary driving and local dephasing
- **Transverse-field Ising model (TFIM)** with local dephasing on a 2D lattice

For each model and each of four simulation methods (first/second-order, deterministic/randomised),
the code computes:

1. **λ = max_k ‖L_k‖◇** — the largest diamond norm among the local Liouvillian terms,
   using the Nechita et al. Choi-Jamiolkowski upper bound (Proposition 1,
   J. Math. Phys. 59, 052201, 2018). Each term is 1- or 2-local, so the Choi matrix
   is constructed only on the local subsystem, making the computation independent of
   system size.
2. **Analytic Trotter-step bounds** (N_analytic) and corresponding gate complexities.
3. **Empirically minimal Trotter-step counts** (N_min) via binary search on the
   precision function, and corresponding gate complexities.

## Code Flow

The pipeline runs in four stages:

```
models/diamond_norm_bound.py          models/tfim_lattice.py
         |                                     |
         └──────────────┬──────────────────────┘
                        ↓
               compute_lambda.py
               (λ = max_k ‖L_k‖◇)
                        |
                        └──────────────────────┐
                                               ↓
computing_bounds/analytic_bounds.py   compute_bounds.py
computing_bounds/error_functions.py   (N_analytic, N_min,
computing_bounds/binary_search.py  →   gate complexities)
                                               |
                              ┌────────────────┼────────────────┐
                              ↓                ↓                ↓
                  results/data/          results/data/    (read by plots)
                  results_xx_chain.json  results_tfim.json
                              |                |
              ┌───────────────┼────────────────┘
              ↓               ↓               ↓
       plot_tfim_couplings  plot_fig2.py   plot_fig3.py
       (Fig1.png)           (Fig2.png)     (Fig3.png)
```

**Fig1** (lattice architecture) is λ-independent — it reads only from
`models/tfim_lattice.py` and does not require running `compute_bounds.py` first.

**Fig2** and **Fig3** read from the JSON datasets produced by `compute_bounds.py`.

## Repository Structure

```
.
├── models/
│   ├── diamond_norm_bound.py    # Nechita et al. CJ-matrix diamond-norm bound
│   ├── xx_spin_chain.py         # XX-chain full Liouvillian builder
│   └── tfim_lattice.py          # TFIM Liouvillian builder and lattice connectivity
├── computing_bounds/
│   ├── analytic_bounds.py       # Analytic N upper bounds (Propositions 2–5)
│   ├── error_functions.py       # Precision functions epsilon_hat (Table 1)
│   └── binary_search_methods.py # Minimum-N search (Algorithm 1)
├── compute_lambda.py            # Compute λ = max_k ‖L_k‖◇ for each model
├── compute_bounds.py            # Generate N_analytic, N_min, gate complexities
├── plot_tfim_couplings.py       # Figure 1: TFIM lattice architectures
├── plot_fig2.py                 # Figure 2: XX-spin chain bounds
├── plot_fig3.py                 # Figure 3: TFIM bounds
└── results/
    ├── data/                    # JSON datasets (output of compute_bounds.py)
    └── figures/                 # Rendered figures (Fig1.png, Fig2.png, Fig3.png)
```

## Environment and Dependencies

- Python `>=3.12`
- `numpy==1.26.0`, `scipy==1.16.0`, `matplotlib==3.10.3`, `qutip==5.2.0`
- `cvxpy==1.6.6`, `cvxopt==1.3.2` (required by qutip's SDP solver)

Dependencies are pinned in `pyproject.toml` / `uv.lock` and mirrored in `requirements.txt`.

## Quick Start

### 1) Install dependencies

Using `uv` (recommended):

```bash
uv sync
```

Or with pip:

```bash
pip install -r requirements.txt
```

### 2) Compute λ (optional diagnostic)

```bash
uv run python compute_lambda.py
```

Prints ‖L_k‖◇ for each local term in both models and reports λ.

### 3) Compute bound datasets

```bash
uv run python compute_bounds.py
```

Outputs:

- `results/data/results_xx_chain.json`
- `results/data/results_tfim.json`

### 4) Generate figures

```bash
uv run python plot_tfim_couplings.py   # Figure 1 (no compute_bounds.py needed)
uv run python plot_fig2.py             # Figure 2 (requires step 3)
uv run python plot_fig3.py             # Figure 3 (requires step 3)
```

Outputs: `results/figures/Fig1.png`, `Fig2.png`, `Fig3.png`.

## Reproducibility Notes

- λ is computed from the model at runtime — not hard-coded. `compute_bounds.py`
  always recomputes λ before generating results.
- Gate-complexity formulas follow Table 1 of the manuscript: `M*N` for first-order
  methods and `2*M*N` for second-order methods.
- TFIM lattice connectivity is defined in `models/tfim_lattice.py`; M is derived
  from connectivity rather than hard-coded.
- Plot scripts print parameter values (λ, ε, t) read directly from the JSON output,
  so printed captions always match the data actually plotted.

## Citation

If you use this code in your research, please cite the corresponding manuscript.
A DOI will be added here upon publication.
