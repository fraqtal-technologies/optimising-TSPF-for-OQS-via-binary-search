"""
Compute Trotter-step and gate-complexity bounds for Figures 2 and 3.

For each model and simulation method, computes:
  - N_analytic : analytic upper bound on Trotter steps (Propositions 2-5)
  - N_min      : empirically minimal Trotter steps via binary search (Algorithm 1)
  - g_analytic : gate complexity from N_analytic
  - g_min      : gate complexity from N_min

Gate complexity is M*N for first-order methods and 2*M*N for second-order methods,
where M is the number of Liouvillian terms (Table 1 of the manuscript).

Outputs JSON to results/data/.
"""
import json
import os

from computing_bounds.analytic_bounds import (
    N_analytic_det1,
    N_analytic_ran1,
    N_analytic_det2,
    N_analytic_ran2,
)
from computing_bounds.binary_search_methods import binary_search
from computing_bounds.error_functions import (
    epsilon_hat_det1,
    epsilon_hat_ran1,
    epsilon_hat_det2,
    epsilon_hat_ran2,
)
from models.tfim_lattice import _LATTICE_PAIRS
from compute_lambda import (
    build_local_matrices_xx,
    build_local_matrices_tfim,
    compute_lambda,
)

ROOT = os.path.dirname(__file__)
RESULTS_DATA_DIR = os.path.join(ROOT, "results", "data")


def _gate_complexity(M, N, second_order):
    """Gate complexity: M*N for first-order, 2*M*N for second-order (Table 1)."""
    return (2 if second_order else 1) * M * N


def _run_all_methods(t, lam, M, epsilon):
    """
    Compute analytic and minimal Trotter-step counts and gate complexities for all
    four simulation methods at fixed (t, lam, M, epsilon).
    """
    return {
        "det1": {
            "N_analytic": N_analytic_det1(t, lam, M, epsilon),
            "N_min":      binary_search(t, lam, M, epsilon, epsilon_hat_det1),
            "g_analytic": _gate_complexity(M, N_analytic_det1(t, lam, M, epsilon), second_order=False),
            "g_min":      _gate_complexity(M, binary_search(t, lam, M, epsilon, epsilon_hat_det1), second_order=False),
        },
        "ran1": {
            "N_analytic": N_analytic_ran1(t, lam, M, epsilon),
            "N_min":      binary_search(t, lam, M, epsilon, epsilon_hat_ran1),
            "g_analytic": _gate_complexity(M, N_analytic_ran1(t, lam, M, epsilon), second_order=False),
            "g_min":      _gate_complexity(M, binary_search(t, lam, M, epsilon, epsilon_hat_ran1), second_order=False),
        },
        "det2": {
            "N_analytic": N_analytic_det2(t, lam, M, epsilon),
            "N_min":      binary_search(t, lam, M, epsilon, epsilon_hat_det2),
            "g_analytic": _gate_complexity(M, N_analytic_det2(t, lam, M, epsilon), second_order=True),
            "g_min":      _gate_complexity(M, binary_search(t, lam, M, epsilon, epsilon_hat_det2), second_order=True),
        },
        "ran2": {
            "N_analytic": N_analytic_ran2(t, lam, M, epsilon),
            "N_min":      binary_search(t, lam, M, epsilon, epsilon_hat_ran2),
            "g_analytic": _gate_complexity(M, N_analytic_ran2(t, lam, M, epsilon), second_order=True),
            "g_min":      _gate_complexity(M, binary_search(t, lam, M, epsilon, epsilon_hat_ran2), second_order=True),
        },
    }


def compute_xx_results():
    """
    XX-spin chain: P in {2,...,8}, M = 2P+3 (P-1 bonds + 4 boundary dissipators + P dephasing).
    Parameters: t=2, epsilon=1e-3, Omega=3.94, gamma=0.31.
    """
    t = 2.0
    epsilon = 1e-3
    Omega, gamma = 3.94, 0.31

    # lambda = max_k ||L_k||_diamond via the Nechita et al. CJ-matrix bound.
    # All terms are 1- or 2-local, so lambda is independent of chain length P.
    lam = compute_lambda(build_local_matrices_xx(P=2, Omega=Omega, gamma=gamma), verbose=False)

    p_values = [2, 3, 4, 5, 6, 7, 8]
    rows = [
        {"P": p, "M": 2 * p + 3, "methods": _run_all_methods(t, lam, 2 * p + 3, epsilon)}
        for p in p_values
    ]
    return {
        "model": "xx_spin_chain",
        "params": {"t": t, "epsilon": epsilon, "lam": lam, "Omega": Omega, "gamma": gamma},
        "rows": rows,
    }


def compute_tfim_results():
    """
    TFIM: n_spins in {2,3,4,5,6}, M = n_pairs + 2*n_spins (ZZ + X-field + Z-dephasing).
    Parameters: t=5, epsilon=1e-5, J=1.0, h=0.5, gamma=0.1.
    """
    t = 5.0
    epsilon = 1e-5
    J, h, gamma = 1.0, 0.5, 0.1

    # lambda = max_k ||L_k||_diamond via the Nechita et al. CJ-matrix bound.
    # All terms are 1- or 2-local, so lambda is independent of system size.
    lam = compute_lambda(build_local_matrices_tfim(n_spins=2, J=J, h=h, gamma=gamma), verbose=False)

    n_spins_values = [2, 3, 4, 5, 6]
    rows = [
        {
            "n_spins": n,
            "M": len(_LATTICE_PAIRS[n]) + 2 * n,
            "methods": _run_all_methods(t, lam, len(_LATTICE_PAIRS[n]) + 2 * n, epsilon),
        }
        for n in n_spins_values
    ]
    return {
        "model": "tfim_lattice",
        "params": {"t": t, "epsilon": epsilon, "lam": lam, "J": J, "h": h, "gamma": gamma},
        "rows": rows,
    }


if __name__ == "__main__":
    os.makedirs(RESULTS_DATA_DIR, exist_ok=True)

    for payload, filename in [
        (compute_xx_results(),  "results_xx_chain.json"),
        (compute_tfim_results(), "results_tfim.json"),
    ]:
        path = os.path.join(RESULTS_DATA_DIR, filename)
        with open(path, "w", encoding="utf-8") as f:
            json.dump(payload, f, indent=2)
        print(f"Wrote {path}")
