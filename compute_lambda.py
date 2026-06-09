"""
Compute lambda = max_k ||L_k||_diamond for the XX-Spin Chain and TFIM models.

Uses the CJ-matrix bound (Nechita et al.) on LOCAL superoperator matrices.
Each term is 1- or 2-local: d=2 (4x4 matrix) or d=4 (16x16 matrix).
Cost is independent of total system size P.
"""
import sys
import os
import numpy as np
import qutip as qt

sys.path.insert(0, os.path.dirname(__file__))
from models.diamond_norm_bound import diamond_norm_bound
from models.tfim_lattice import _LATTICE_PAIRS


def build_local_matrices_xx(P, Omega=3.94, gamma=0.31):
    """
    Return list of (L_local_matrix, d_local, label) for each term in the XX-chain.
    XX and YY coupling are combined into a single term per bond (as in M = 2P+3).
    """
    terms = []
    zero1 = 0 * qt.qeye(2)

    # Combined XX+YY coupling per bond (2-site, d=4)
    for j in range(P - 1):
        H = qt.tensor(qt.sigmax(), qt.sigmax()) + qt.tensor(qt.sigmay(), qt.sigmay())
        mat = qt.liouvillian(H=H, c_ops=[]).full()
        terms.append((mat, 4, f"XX+YY coupling bond {j}"))

    # Boundary dissipators: separate sigma+/- per site (1-site, d=2)
    prefactor = np.sqrt(Omega / 2)
    for site_label in ["site 0 sigma+", "site 0 sigma-", "site P-1 sigma+", "site P-1 sigma-"]:
        jump = qt.create(2) if "sigma+" in site_label else qt.destroy(2)
        c_op = prefactor * jump
        mat = qt.liouvillian(H=zero1, c_ops=[c_op]).full()
        terms.append((mat, 2, f"boundary {site_label}"))

    # Local dephasing (1-site, d=2)
    dephasing_prefactor = np.sqrt(gamma / 2)
    for j in range(P):
        c_op = dephasing_prefactor * qt.sigmaz()
        mat = qt.liouvillian(H=zero1, c_ops=[c_op]).full()
        terms.append((mat, 2, f"dephasing site {j}"))

    return terms


def build_local_matrices_tfim(n_spins, J=1.0, h=0.5, gamma=0.1):
    """
    Return list of (L_local_matrix, d_local, label) for each term in the TFIM.
    """
    terms = []
    zero1 = 0 * qt.qeye(2)
    pairs = _LATTICE_PAIRS[n_spins]

    # ZZ coupling (2-site, d=4)
    for (i, k) in pairs:
        H = 1j * J * qt.tensor(qt.sigmaz(), qt.sigmaz())
        mat = qt.liouvillian(H=H, c_ops=[]).full()
        terms.append((mat, 4, f"ZZ coupling sites ({i},{k})"))

    # X-field (1-site, d=2)
    for j in range(n_spins):
        H = 1j * h * qt.sigmax()
        mat = qt.liouvillian(H=H, c_ops=[]).full()
        terms.append((mat, 2, f"X-field site {j}"))

    # Z-dephasing (1-site, d=2)
    sqrt_gamma = np.sqrt(gamma)
    for j in range(n_spins):
        c_op = sqrt_gamma * qt.sigmaz()
        mat = qt.liouvillian(H=zero1, c_ops=[c_op]).full()
        terms.append((mat, 2, f"Z-dephasing site {j}"))

    return terms


def compute_lambda(local_terms, verbose=True):
    lam = 0.0
    for mat, d_local, label in local_terms:
        bound, is_exact = diamond_norm_bound(mat, d_local)
        if verbose:
            exact_str = "exact" if is_exact else "upper bound"
            print(f"  {label:40s} ||L_k||_diamond = {bound:.6f} ({exact_str})")
        lam = max(lam, bound)
    if verbose:
        print(f"  lambda = {lam:.6f}")
    return lam


if __name__ == "__main__":
    print("=" * 60)
    print("XX-Spin Chain (P=2, Omega=3.94, gamma=0.31)")
    print("=" * 60)
    lam_xx = compute_lambda(build_local_matrices_xx(P=2, Omega=3.94, gamma=0.31))
    print()

    print("=" * 60)
    print("TFIM (n_spins=2, J=1, h=0.5, gamma=0.1)")
    print("=" * 60)
    lam_tfim = compute_lambda(build_local_matrices_tfim(n_spins=2, J=1.0, h=0.5, gamma=0.1))
