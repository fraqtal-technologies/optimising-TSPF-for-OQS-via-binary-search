"""
Compute analytic/minimal Trotter-step and gate-complexity bounds for Figures 2 and 3.
"""
import json
import os
import sys

ROOT = os.path.dirname(__file__)
COMPUTING_BOUNDS_DIR = os.path.join(ROOT, "computing-bounds")
RESULTS_DATA_DIR = os.path.join(ROOT, "results", "data")
sys.path.insert(0, COMPUTING_BOUNDS_DIR)

from analytic_bounds import (  # noqa: E402
    N_analytic_det1,
    N_analytic_ran1,
    N_analytic_det2,
    N_analytic_ran2,
)
from binary_search_methods import binary_search  # noqa: E402
from error_functions import (  # noqa: E402
    epsilon_hat_det1,
    epsilon_hat_ran1,
    epsilon_hat_det2,
    epsilon_hat_ran2,
)
from models.tfim_lattice import _LATTICE_PAIRS  # noqa: E402


def _gate_complexities(M, N, second_order):
    factor = 2 if second_order else 1
    return factor * M * N


def _record_method_results(t, lam, M, epsilon):
    n_analytic_det1 = N_analytic_det1(t, lam, M, epsilon)
    n_analytic_ran1 = N_analytic_ran1(t, lam, M, epsilon)
    n_analytic_det2 = N_analytic_det2(t, lam, M, epsilon)
    n_analytic_ran2 = N_analytic_ran2(t, lam, M, epsilon)

    n_min_det1 = binary_search(t, lam, M, epsilon, epsilon_hat_det1)
    n_min_ran1 = binary_search(t, lam, M, epsilon, epsilon_hat_ran1)
    n_min_det2 = binary_search(t, lam, M, epsilon, epsilon_hat_det2)
    n_min_ran2 = binary_search(t, lam, M, epsilon, epsilon_hat_ran2)

    return {
        "det1": {
            "N_analytic": n_analytic_det1,
            "N_min": n_min_det1,
            "g_analytic": _gate_complexities(M, n_analytic_det1, second_order=False),
            "g_min": _gate_complexities(M, n_min_det1, second_order=False),
        },
        "ran1": {
            "N_analytic": n_analytic_ran1,
            "N_min": n_min_ran1,
            "g_analytic": _gate_complexities(M, n_analytic_ran1, second_order=False),
            "g_min": _gate_complexities(M, n_min_ran1, second_order=False),
        },
        "det2": {
            "N_analytic": n_analytic_det2,
            "N_min": n_min_det2,
            "g_analytic": _gate_complexities(M, n_analytic_det2, second_order=True),
            "g_min": _gate_complexities(M, n_min_det2, second_order=True),
        },
        "ran2": {
            "N_analytic": n_analytic_ran2,
            "N_min": n_min_ran2,
            "g_analytic": _gate_complexities(M, n_analytic_ran2, second_order=True),
            "g_min": _gate_complexities(M, n_min_ran2, second_order=True),
        },
    }


def compute_xx_results():
    t = 2.0
    epsilon = 1e-3
    lam = 8.119
    p_values = [2, 3, 4, 5, 6, 7, 8]
    m_values = [2 * p + 3 for p in p_values]

    rows = []
    for p, m in zip(p_values, m_values):
        rows.append(
            {
                "P": p,
                "M": m,
                "methods": _record_method_results(t=t, lam=lam, M=m, epsilon=epsilon),
            }
        )

    return {
        "model": "xx_spin_chain",
        "params": {"t": t, "epsilon": epsilon, "lam": lam, "Omega": 3.94, "gamma": 0.31},
        "rows": rows,
    }


def compute_tfim_results():
    t = 5.0
    epsilon = 1e-5
    lam = 4.00
    n_spins_values = [2, 3, 4, 5, 6]
    m_values = [len(_LATTICE_PAIRS[n]) + 2 * n for n in n_spins_values]

    rows = []
    for n_spins, m in zip(n_spins_values, m_values):
        rows.append(
            {
                "n_spins": n_spins,
                "M": m,
                "methods": _record_method_results(t=t, lam=lam, M=m, epsilon=epsilon),
            }
        )

    return {
        "model": "tfim_lattice",
        "params": {"t": t, "epsilon": epsilon, "lam": lam, "J": 1.0, "h": 0.5, "gamma": 0.1},
        "rows": rows,
    }


def _write_json(path, payload):
    with open(path, "w", encoding="utf-8") as f:
        json.dump(payload, f, indent=2)


if __name__ == "__main__":
    xx_results = compute_xx_results()
    tfim_results = compute_tfim_results()

    os.makedirs(RESULTS_DATA_DIR, exist_ok=True)
    xx_path = os.path.join(RESULTS_DATA_DIR, "results_xx_chain.json")
    tfim_path = os.path.join(RESULTS_DATA_DIR, "results_tfim.json")
    _write_json(xx_path, xx_results)
    _write_json(tfim_path, tfim_results)

    print(f"Wrote {xx_path}")
    print(f"Wrote {tfim_path}")
