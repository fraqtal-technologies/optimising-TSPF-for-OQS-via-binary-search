"""
Binary search for the minimum Trotter-step count N (Algorithm 1 of the manuscript).

The search finds the smallest integer N in [1, N_upper] such that
    error_function(t, lam, M, N) <= epsilon.

Phase 1 doubles N_upper until the error condition is satisfied, establishing a
valid upper bound. Phase 2 applies standard binary search within [1, N_upper]
to find the minimum such N. Total evaluations of error_function: O(log N_upper).
"""


def binary_search(t, lam, M, epsilon, error_function):
    """
    Find the minimum integer N >= 1 such that error_function(t, lam, M, N) <= epsilon.

    Args:
        t              : simulation time
        lam            : lambda = max_k ||L_k||_diamond
        M              : number of Liouvillian terms
        epsilon        : target precision
        error_function : one of epsilon_hat_{det1,ran1,det2,ran2} from error_functions.py

    Returns:
        N_min (int): minimum Trotter-step count achieving the target precision.
    """
    N_lower = 1
    N_upper = 1

    # Phase 1: double N_upper until epsilon_hat(N_upper) <= epsilon
    while error_function(t, lam, M, N_upper) > epsilon:
        N_upper *= 2

    # Phase 2: binary search on [N_lower, N_upper] for the minimum valid N
    while N_lower < N_upper:
        N_mid = (N_upper + N_lower) // 2
        if error_function(t, lam, M, N_mid) > epsilon:
            N_lower = N_mid + 1
        else:
            N_upper = N_mid

    return N_lower
