"""
Diamond norm upper bound via the Nechita et al. CJ-matrix formula (Proposition 1,
inequality (7) of Nechita, Puchala, Pawela, Zyczkowski, J. Math. Phys. 59, 052201 (2018)).

Uses locality: each L_k is k-local, so the Choi matrix is computed only on the local
subsystem (d=2 for 1-site, d=4 for 2-site), making the cost independent of total
system size P.

Formula:
    ||Phi||_diamond <= (||tr2(sqrt(J(Phi)^dag J(Phi)))||_inf + ||tr2(sqrt(J(Phi) J(Phi)^dag))||_inf) / 2

where J(Phi) is the Choi-Jamiolkowski matrix of Phi, tr2 is the partial trace over the
second tensor factor, and ||.||_inf is the operator norm (largest singular value).

The bound is tight when both partial-trace matrices are scalar multiples of the identity.
"""
import numpy as np
from scipy.linalg import sqrtm


def choi_matrix(L_local_matrix, d_local):
    """
    Construct the Choi-Jamiolkowski matrix of a superoperator from its matrix representation.

    The superoperator is applied to each matrix basis element E_{ab} and the output is
    assembled into the Choi matrix using the column-stacking (Fortran-order) convention.

    Args:
        L_local_matrix: (d^2 x d^2) superoperator matrix on the local subsystem.
        d_local: local Hilbert space dimension.

    Returns:
        J: (d^2 x d^2) Choi-Jamiolkowski matrix.
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
    """Partial trace over the second factor of a (d^2 x d^2) matrix viewed as C^d x C^d."""
    return np.einsum('ijkj->ik', M.reshape(d, d, d, d))


def diamond_norm_bound(L_local_matrix, d_local):
    """
    Evaluate the Nechita et al. upper bound on the diamond norm of a local superoperator.

    Args:
        L_local_matrix: (d^2 x d^2) superoperator matrix on the local subsystem only.
        d_local: local Hilbert space dimension (2 for 1-site, 4 for 2-site).

    Returns:
        (bound, is_exact): the upper bound value, and a boolean indicating whether the
                           bound is tight (both partial-trace matrices are scalar multiples
                           of the identity, which is the equality condition in Proposition 1).
    """
    J = choi_matrix(L_local_matrix, d_local)

    sqrt_JdJ = sqrtm(J.conj().T @ J)
    sqrt_JJd = sqrtm(J @ J.conj().T)

    # sqrtm can introduce small imaginary parts for near-PSD inputs; discard if negligible.
    if np.max(np.abs(np.imag(sqrt_JdJ))) < 1e-10:
        sqrt_JdJ = np.real(sqrt_JdJ)
    if np.max(np.abs(np.imag(sqrt_JJd))) < 1e-10:
        sqrt_JJd = np.real(sqrt_JJd)

    A = partial_trace_second(sqrt_JdJ, d_local)
    B = partial_trace_second(sqrt_JJd, d_local)

    # ||.||_inf in Nechita et al. is the operator norm (largest singular value).
    # In numpy this is ord=2 for a 2-D array; ord=inf gives the max-row-sum norm.
    bound = (np.linalg.norm(A, ord=2) + np.linalg.norm(B, ord=2)) / 2.0

    def _is_scalar_id(M, tol=1e-8):
        diag_vals = np.diag(M)
        return (np.max(np.abs(M - np.diag(diag_vals))) < tol and
                np.max(np.abs(diag_vals - diag_vals[0])) < tol)

    is_exact = _is_scalar_id(A) and _is_scalar_id(B)
    return bound, is_exact
