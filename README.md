# Optimising TS-PF for OQS via Binary Search

Code accompanying research on resource-efficient simulation of Markovian open quantum systems using deterministic and randomised Trotter-Suzuki product formulas (TS-PF), with analytic bounds and binary-search-optimised empirical bounds on the required Trotter step count.

## Scope

This repository reproduces the numerical workflows used for:

- XX-spin chain with boundary driving and local dephasing
- Transverse-field Ising model (TFIM) with local dephasing
- Comparison of analytic versus empirically optimised bounds for:
  - first-order deterministic TS-PF
  - first-order randomised TS-PF
  - second-order deterministic TS-PF
  - second-order randomised TS-PF

## Repository Structure

- `models/` - model-specific Liouvillian builders and diamond-norm bound routines
- `computing-bounds/` - error functions, analytic bounds, and binary search
- `compute_lambda.py` - lambda computation (`lambda = max_k ||L_k||_diamond`)
- `compute_bounds.py` - generates bound datasets for figures
- `plot_fig2.py`, `plot_fig3.py` - generates publication figures
- `results/`
  - `results/data/` - JSON datasets used for plotting
  - `results/figures/` - rendered figure outputs

## Environment and Dependencies

- Python `>=3.12`
- `numpy==1.26.0`
- `scipy==1.16.0`
- `matplotlib==3.10.3`
- `qutip==5.2.0`

Dependencies are pinned in `pyproject.toml` / `uv.lock` (also available in `requirements.txt`).

## Quick Start

### 1) Install dependencies

Using `uv` (recommended):

```bash
uv sync
```

### 2) Compute bound datasets

```bash
uv run python compute_bounds.py
```

Outputs:

- `results/data/results_xx_chain.json`
- `results/data/results_tfim.json`

### 3) Generate figures

```bash
uv run python plot_fig2.py
uv run python plot_fig3.py
uv run python plot_tfim_couplings.py
```

Outputs:

- `results/figures/fig2.png`
- `results/figures/fig3.png`
- `results/figures/tfim_couplings.png`

### 4) (Optional) Recompute lambda diagnostics

```bash
uv run python compute_lambda.py
```

## Reproducibility Notes

- Gate-complexity formulas are implemented exactly as in Table 1 of the manuscript (`MN` for first-order methods and `2MN` for second-order methods).
- Plot scripts print figure-caption text with the parameter values used in the runs.
- TFIM lattice connectivity is encoded in `models/tfim_lattice.py`; `compute_bounds.py` derives TFIM `M` values directly from connectivity to avoid hard-coded inconsistencies.

## Citation

If you use this codebase in academic work, please cite the corresponding manuscript and any linked repository release/DOI once published.