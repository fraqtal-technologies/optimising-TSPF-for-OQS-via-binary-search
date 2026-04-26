import qutip as qt
import numpy as np


def build_xx_liouvillian(P, Omega=3.94, gamma=0.31):
    """
    Build the XX-Spin Chain Liouvillian for P sites.

    H = sum_{j=0}^{P-2} (X_j X_{j+1} + Y_j Y_{j+1})

    Liouvillian terms (M = 2P + 3):
      - (P-1) combined XX+YY coupling terms per bond
      - 4 boundary dissipators: sqrt(Omega/2) sigma+/- on sites 0 and P-1
      - P local dephasing terms: sqrt(gamma/2) Z_j per site

    Returns dict with keys: 'L', 'L_as_list', 'M'
    """
    zero_H = 0 * qt.tensor([qt.qeye(2)] * P)
    terms = []

    # Combined XX+YY coupling per bond (single term per bond)
    for j in range(P - 1):
        ops_xx = [qt.qeye(2)] * P
        ops_xx[j] = qt.sigmax()
        ops_xx[j + 1] = qt.sigmax()
        ops_yy = [qt.qeye(2)] * P
        ops_yy[j] = qt.sigmay()
        ops_yy[j + 1] = qt.sigmay()
        H_bond = qt.tensor(ops_xx) + qt.tensor(ops_yy)
        terms.append(qt.liouvillian(H=H_bond, c_ops=[]))

    # Boundary dissipators: sigma+/- on site 0 and site P-1 (4 separate terms)
    prefactor = np.sqrt(Omega / 2)
    for site in [0, P - 1]:
        for jump in [qt.create(2), qt.destroy(2)]:
            ops = [qt.qeye(2)] * P
            ops[site] = prefactor * jump
            terms.append(qt.liouvillian(H=zero_H, c_ops=[qt.tensor(ops)]))

    # Local dephasing: sqrt(gamma/2) Z_j on each site
    dephasing_prefactor = np.sqrt(gamma / 2)
    for j in range(P):
        ops = [qt.qeye(2)] * P
        ops[j] = dephasing_prefactor * qt.sigmaz()
        terms.append(qt.liouvillian(H=zero_H, c_ops=[qt.tensor(ops)]))

    M = 2 * P + 3
    assert len(terms) == M, f"Expected M={M} terms, got {len(terms)}"

    return {'L': sum(terms), 'L_as_list': terms, 'M': M}
