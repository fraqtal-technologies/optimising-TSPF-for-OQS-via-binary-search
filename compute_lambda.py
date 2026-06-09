"""
Compute lambda = max_k ||L_k||_diamond for the XX-spin chain and TFIM models.

For each local Liouvillian term L_k, the diamond norm is bounded using the
Nechita et al. Choi-Jamiolkowski matrix formula (Proposition 1 of Nechita,
Puchala, Pawela, Zyczkowski, J. Math. Phys. 59, 052201, 2018). Because every
term is 1- or 2-local, the Choi matrix is constructed only on the local subsystem
(4x4 for 1-site terms, 16x16 for 2-site terms), making the cost independent of
total system size.

The bound is tight (equals the exact diamond norm) when the relevant partial-trace
matrices are scalar multiples of the identity; the output labels these cases as "exact".
"""
import numpy as np
import qutip as qt

from models.diamond_norm_bound import diamond_norm_bound
from models.tfim_lattice import _LATTICE_PAIRS


def build_local_matrices_xx(P, Omega=3.94, gamma=0.31):
    """
    Build local superoperator matrices for all terms in the XX-spin chain Liouvillian.

    The XX-chain has M = 2P+3 terms:
      - (P-1) combined XX+YY coupling terms per bond  (2-site, d=4)
      -  4    boundary dissipators: sqrt(Omega/2) * sigma+/- on sites 0 and P-1  (1-site, d=2)
      -  P    local dephasing terms: sqrt(gamma/2) * Z_j  (1-site, d=2)

    Returns:
        List of (L_local_matrix, d_local, label) tuples, one per Liouvillian term.
    """
    terms = []
    zero1 = 0 * qt.qeye(2)
    prefactor = np.sqrt(Omega / 2)
    dephasing_prefactor = np.sqrt(gamma / 2)

    # Combined XX+YY coupling per bond (same for all bonds, so only geometry-independent
    # local matrix needed; bond index is kept for labelling only)
    for j in range(P - 1):
        H = qt.tensor(qt.sigmax(), qt.sigmax()) + qt.tensor(qt.sigmay(), qt.sigmay())
        mat = qt.liouvillian(H=H, c_ops=[]).full()
        terms.append((mat, 4, f"XX+YY bond {j}"))

    # Boundary dissipators: sigma+ and sigma- on site 0 and site P-1
    for site_label in ["site 0 sigma+", "site 0 sigma-", "site P-1 sigma+", "site P-1 sigma-"]:
        jump = qt.create(2) if "sigma+" in site_label else qt.destroy(2)
        mat = qt.liouvillian(H=zero1, c_ops=[prefactor * jump]).full()
        terms.append((mat, 2, f"boundary {site_label}"))

    # Local dephasing on each site
    for j in range(P):
        mat = qt.liouvillian(H=zero1, c_ops=[dephasing_prefactor * qt.sigmaz()]).full()
        terms.append((mat, 2, f"dephasing site {j}"))

    return terms


def build_local_matrices_tfim(n_spins, J=1.0, h=0.5, gamma=0.1):
    """
    Build local superoperator matrices for all terms in the TFIM Liouvillian.

    The TFIM has M = n_pairs + 2*n_spins terms:
      - n_pairs  ZZ coupling terms: commutator with J*Z_i Z_j  (2-site, d=4)
      - n_spins  X-field terms: commutator with h*X_j  (1-site, d=2)
      - n_spins  Z-dephasing terms: jump operator sqrt(gamma)*Z_j  (1-site, d=2)

    Returns:
        List of (L_local_matrix, d_local, label) tuples, one per Liouvillian term.
    """
    terms = []
    zero1 = 0 * qt.qeye(2)
    pairs = _LATTICE_PAIRS[n_spins]

    # ZZ coupling (local 2-site Liouvillian)
    for (i, k) in pairs:
        H = 1j * J * qt.tensor(qt.sigmaz(), qt.sigmaz())
        mat = qt.liouvillian(H=H, c_ops=[]).full()
        terms.append((mat, 4, f"ZZ coupling sites ({i},{k})"))

    # X-field (local 1-site Liouvillian)
    for j in range(n_spins):
        H = 1j * h * qt.sigmax()
        mat = qt.liouvillian(H=H, c_ops=[]).full()
        terms.append((mat, 2, f"X-field site {j}"))

    # Z-dephasing (local 1-site Liouvillian)
    for j in range(n_spins):
        mat = qt.liouvillian(H=zero1, c_ops=[np.sqrt(gamma) * qt.sigmaz()]).full()
        terms.append((mat, 2, f"Z-dephasing site {j}"))

    return terms


def compute_lambda(local_terms, verbose=True):
    """
    Compute lambda = max_k ||L_k||_diamond over a list of local Liouvillian terms.

    Args:
        local_terms : list of (L_local_matrix, d_local, label) as returned by
                      build_local_matrices_xx or build_local_matrices_tfim.
        verbose     : if True, print each term's bound and the final lambda.

    Returns:
        lam (float): the maximum diamond-norm bound across all terms.
    """
    lam = 0.0
    for mat, d_local, label in local_terms:
        bound, is_exact = diamond_norm_bound(mat, d_local)
        if verbose:
            tag = "exact" if is_exact else "upper bound"
            print(f"  {label:40s} ||L_k||_diamond = {bound:.6f} ({tag})")
        lam = max(lam, bound)
    if verbose:
        print(f"  lambda = {lam:.6f}")
    return lam


if __name__ == "__main__":
    print("=" * 60)
    print("XX-Spin Chain (P=2, Omega=3.94, gamma=0.31)")
    print("=" * 60)
    compute_lambda(build_local_matrices_xx(P=2, Omega=3.94, gamma=0.31))

    print()
    print("=" * 60)
    print("TFIM (n_spins=2, J=1.0, h=0.5, gamma=0.1)")
    print("=" * 60)
    compute_lambda(build_local_matrices_tfim(n_spins=2, J=1.0, h=0.5, gamma=0.1))
