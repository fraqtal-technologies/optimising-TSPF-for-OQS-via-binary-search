"""
Diamond norm upper bound via the Nechita et al. CJ-matrix formula (Eq. 32 of the paper).

Uses locality: each L_k is k-local, so the CJ matrix is computed only on the local
subsystem (d=2 for 1-site, d=4 for 2-site). This makes the cost independent of total
system size P. For O(P)-local terms the approach becomes classically intractable beyond
~10-15 qubits, but all terms in the current models are 1- or 2-local.

Formula:
    ||Phi||_diamond <= (||tr2(sqrt(J†J))||_inf + ||tr2(sqrt(JJ†))||_inf) / 2

where J is the Choi-Jamiolkowski matrix of Phi, tr2 is the partial trace over the second
tensor factor, and ||.||_inf is the spectral norm (largest singular value).

In the standard basis, J = L_mat (the superoperator matrix itself).
"""
import numpy as np
from scipy.linalg import sqrtm


def choi_matrix(L_local_matrix, d_local):
    """
    Build Choi matrix explicitly by applying the superoperator to basis matrices.

    This provides an explicit CJ map constructor for validation and convention checks.
    """
    d2 = d_local * d_local
    J = np.zeros((d2, d2), dtype=complex)

    for a in range(d_local):
        for b in range(d_local):
            E_ab = np.zeros((d_local, d_local), dtype=complex)
            E_ab[a, b] = 1.0

            vec_E_ab = E_ab.reshape(d2, order="F")
            vec_out = L_local_matrix @ vec_E_ab
            out = vec_out.reshape((d_local, d_local), order="F")

            J[a * d_local:(a + 1) * d_local, b * d_local:(b + 1) * d_local] = out

    return J


def partial_trace_second(M, d):
    """Partial trace over second factor of a (d^2 x d^2) matrix viewed as C^d ⊗ C^d."""
    return np.einsum('ijkj->ik', M.reshape(d, d, d, d))


def diamond_norm_bound(L_local_matrix, d_local):
    """
    Compute the Nechita et al. CJ-matrix upper bound on the diamond norm.

    Args:
        L_local_matrix: (d^2 x d^2) superoperator matrix on the local subsystem only.
        d_local: local Hilbert space dimension (2 for 1-site, 4 for 2-site).

    Returns:
        (bound, is_exact): bound value and whether equality holds (both partial traces
                           are scalar multiples of identity).
    """
    # In the basis convention used throughout this codebase, the Liouville-space matrix
    # directly represents the Choi matrix for the CJ bound evaluation.
    J = L_local_matrix

    sqrt_JdJ = sqrtm(J.conj().T @ J)
    sqrt_JJd = sqrtm(J @ J.conj().T)

    # Positive-semidefinite inputs can still yield tiny imaginary numerical noise.
    if np.max(np.abs(np.imag(sqrt_JdJ))) < 1e-10:
        sqrt_JdJ = np.real(sqrt_JdJ)
    if np.max(np.abs(np.imag(sqrt_JJd))) < 1e-10:
        sqrt_JJd = np.real(sqrt_JJd)

    A = partial_trace_second(sqrt_JdJ, d_local)
    B = partial_trace_second(sqrt_JJd, d_local)

    bound = (np.linalg.norm(A, ord=2) + np.linalg.norm(B, ord=2)) / 2.0

    def _is_scalar_id(M, tol=1e-8):
        d = M.shape[0]
        diag_vals = np.diag(M)
        return (np.max(np.abs(M - np.diag(diag_vals))) < tol and
                np.max(np.abs(diag_vals - diag_vals[0])) < tol)

    is_exact = _is_scalar_id(A) and _is_scalar_id(B)
    return bound, is_exact
