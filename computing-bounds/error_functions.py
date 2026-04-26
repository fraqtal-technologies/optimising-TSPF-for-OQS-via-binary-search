import numpy as np

def epsilon_hat_det1(t, lam, M, N):
    return (t * lam * M)**2 / N * np.exp(t * lam * M / N)

def epsilon_hat_ran1(t, lam, M, N):
    return (t * lam * M)**3 / (3 * N**2) * np.exp(t * lam * M / N)

def epsilon_hat_det2(t, lam, M, N):
    return (t * lam * M)**3 / (3 * N**2) * np.exp(t * lam * M / N)

def epsilon_hat_ran2(t, lam, M, N):
    return (t * lam)**3 * M**2 / N**2 * np.exp(t * lam * M / N)
