import numpy as np
from math import ceil, e

def N_analytic_det1(t, lam, M, epsilon):
    return ceil(max(t * lam * M, e * (t * lam * M)**2 / epsilon))

def N_analytic_ran1(t, lam, M, epsilon):
    return ceil(max(t * lam * M, np.sqrt(e * (t * lam * M)**3 / (3 * epsilon))))

def N_analytic_det2(t, lam, M, epsilon):
    return ceil(max(t * lam * M, np.sqrt(e * (t * lam * M)**3 / (3 * epsilon))))

def N_analytic_ran2(t, lam, M, epsilon):
    return ceil(max(t * lam * M, np.sqrt(e * (t * lam)**3 * M**2 / epsilon)))
