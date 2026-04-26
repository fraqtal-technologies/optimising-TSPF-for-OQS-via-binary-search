import qutip as qt
import numpy as np


# Neighbour pairs for each system size (0-indexed)
_LATTICE_PAIRS = {
    2: [(0, 1)],
    3: [(0, 1), (1, 2)],
    4: [(0, 1), (2, 3), (0, 2), (1, 3)],        # 2x2 square
    # 5-spin ring gives n_pairs=5, matching the revision target M=15.
    5: [(0, 1), (1, 2), (2, 3), (3, 4), (4, 0)],
    6: [(0, 1), (1, 2), (3, 4), (4, 5), (0, 3), (1, 4), (2, 5)],  # 2x3 rectangle
}

# Expected M values from the paper
_EXPECTED_M = {2: 5, 3: 8, 4: 12, 5: 15, 6: 19}


def build_tfim_liouvillian(n_spins, J=1.0, h=0.5, gamma=0.1):
    """
    Build the TFIM Liouvillian for n_spins in {2,3,4,5,6}.

    H = -J * sum_{<i,j>} Z_i Z_j  -  h * sum_j X_j

    Liouvillian terms:
      - n_pairs ZZ coupling terms: L_k(rho) = i*J*[Z_i Z_j, rho]
      - n_spins X-field terms:     L_k(rho) = i*h*[X_j, rho]
      - n_spins Z-dephasing terms: jump operator sqrt(gamma) * Z_j

    Returns dict with keys: 'L', 'L_as_list', 'M', 'n_pairs'
    """
    if n_spins not in _LATTICE_PAIRS:
        raise ValueError(f"n_spins must be one of {list(_LATTICE_PAIRS.keys())}")

    pairs = _LATTICE_PAIRS[n_spins]
    n_pairs = len(pairs)
    zero_H = 0 * qt.tensor([qt.qeye(2)] * n_spins)
    terms = []

    # ZZ coupling terms: L_k(rho) = iJ[Z_i Z_j, rho]
    for (i, j) in pairs:
        ops = [qt.qeye(2)] * n_spins
        ops[i] = qt.sigmaz()
        ops[j] = qt.sigmaz()
        H_zz = 1j * J * qt.tensor(ops)
        terms.append(qt.liouvillian(H=H_zz, c_ops=[]))

    # X-field terms: L_k(rho) = ih[X_j, rho]
    for j in range(n_spins):
        ops = [qt.qeye(2)] * n_spins
        ops[j] = qt.sigmax()
        H_x = 1j * h * qt.tensor(ops)
        terms.append(qt.liouvillian(H=H_x, c_ops=[]))

    # Z-dephasing terms: jump operator sqrt(gamma) * Z_j
    sqrt_gamma = np.sqrt(gamma)
    for j in range(n_spins):
        ops = [qt.qeye(2)] * n_spins
        ops[j] = sqrt_gamma * qt.sigmaz()
        c_op = qt.tensor(ops)
        terms.append(qt.liouvillian(H=zero_H, c_ops=[c_op]))

    M = n_pairs + n_spins + n_spins
    assert len(terms) == M, f"Expected M={M} terms, got {len(terms)}"

    expected = _EXPECTED_M.get(n_spins)
    if expected is not None:
        assert M == expected, f"M={M} does not match expected M={expected} for n_spins={n_spins}"

    L_total = sum(terms)
    return {'L': L_total, 'L_as_list': terms, 'M': M, 'n_pairs': n_pairs}
