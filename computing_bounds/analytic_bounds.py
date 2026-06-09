"""
Analytic upper bounds on the minimum number of Trotter steps N.

Each function implements one of Propositions 2-5 of the manuscript. The bounds are
derived by controlling the exponential factor in the precision function: imposing
t*lam*M / N <= 1 gives exp(...) <= e, which yields a closed-form N sufficient to
guarantee epsilon_hat(t, lam, M, N) <= epsilon.

Args common to all functions:
    t       : total simulation time
    lam     : lambda = max_k ||L_k||_diamond
    M       : number of Liouvillian terms
    epsilon : target precision (diamond-norm error tolerance)

Returns:
    N_analytic (int): smallest N satisfying the analytic sufficient condition.
"""
import numpy as np
from math import ceil, e


def N_analytic_det1(t, lam, M, epsilon):
    """Analytic bound for first-order deterministic TS-PF (Proposition 2)."""
    return ceil(max(t * lam * M,
                    e * (t * lam * M)**2 / epsilon))


def N_analytic_ran1(t, lam, M, epsilon):
    """Analytic bound for first-order randomised TS-PF (Proposition 3)."""
    return ceil(max(t * lam * M,
                    np.sqrt(e * (t * lam * M)**3 / (3 * epsilon))))


def N_analytic_det2(t, lam, M, epsilon):
    """Analytic bound for second-order deterministic TS-PF (Proposition 4)."""
    return ceil(max(t * lam * M,
                    np.sqrt(e * (t * lam * M)**3 / (3 * epsilon))))


def N_analytic_ran2(t, lam, M, epsilon):
    """Analytic bound for second-order randomised TS-PF (Proposition 5)."""
    return ceil(max(t * lam * M,
                    np.sqrt(e * (t * lam)**3 * M**2 / epsilon)))
