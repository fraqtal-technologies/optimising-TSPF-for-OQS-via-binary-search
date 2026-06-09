"""
Precision functions epsilon_hat for the four Trotter-Suzuki product-formula methods.

Each function upper-bounds the diamond-norm simulation error:

    epsilon_hat(t, lam, M, N) >= || Lambda(t) - (Lambda_tilde(t/N))^N ||_diamond

where t is simulation time, lam = max_k ||L_k||_diamond, M is the number of
Liouvillian terms, and N is the number of Trotter steps.

These are the specific forms of epsilon_hat used in the analytic bounds
(Propositions 2-5) and the binary search (Algorithm 1). Sources:
  det1, det2 : Suzuki (1990), Table 1 of the manuscript.
  ran1, ran2 : David et al. (2024), Table 1 of the manuscript.
"""
import numpy as np


def epsilon_hat_det1(t, lam, M, N):
    """First-order deterministic TS-PF error bound (Table 1, row 1)."""
    return (t * lam * M)**2 / N * np.exp(t * lam * M / N)


def epsilon_hat_ran1(t, lam, M, N):
    """First-order randomised TS-PF error bound (Table 1, row 2)."""
    return (t * lam * M)**3 / (3 * N**2) * np.exp(t * lam * M / N)


def epsilon_hat_det2(t, lam, M, N):
    """Second-order deterministic TS-PF error bound (Table 1, row 3)."""
    return (t * lam * M)**3 / (3 * N**2) * np.exp(t * lam * M / N)


def epsilon_hat_ran2(t, lam, M, N):
    """Second-order randomised TS-PF error bound (Table 1, row 4)."""
    return (t * lam)**3 * M**2 / N**2 * np.exp(t * lam * M / N)
