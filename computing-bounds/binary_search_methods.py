import numpy as np

def binary_search(t, lam, M, epsilon, error_function):
    """
    Find minimum integer N such that error_function(t, lam, M, N) <= epsilon.
    Implements Algorithm 1 from Appendix A of the paper.
    """
    N_lower = 1
    N_upper = 1

    # Phase 1: double N_upper until error is below epsilon
    while error_function(t, lam, M, N_upper) > epsilon:
        N_upper *= 2

    # Phase 2: lower-bound search for the minimum N with error <= epsilon.
    while N_lower < N_upper:
        N_mid = (N_upper + N_lower) // 2
        current_error = error_function(t, lam, M, N_mid)
        if current_error > epsilon:
            N_lower = N_mid + 1
        else:
            N_upper = N_mid

    return N_lower
