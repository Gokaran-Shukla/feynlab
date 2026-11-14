"""Tools for perturbative series  f(x) = sum_k a_k x^k  produced by DiagMC / CDet.

* domb_sykes(a)         ratio-test estimate of the radius of convergence R and the
                        singularity exponent (Domb & Sykes 1957), with 1/k extrapolation
* pade(a, L, M)         [L/M] Pade approximant (numerator, denominator coefficients)
* pade_eval             evaluate [L/M] at x (works beyond the radius of convergence)
* conformal_map_sum     resummation through the map x = 4 R w/(1-w)^2 for a branch cut on
                        (-inf, -R] (Rossi et al., PRL 121, 130405 (2018) style)
* misleading_convergence_check  Van Houcke et al., SciPost Phys. 16, 133 (2024): a converged
                        bold series is trustworthy only if sum_n xi^n Sigma_n also converges
                        for some xi > 1; returns the partial sums for a list of xi.
"""
from __future__ import annotations

import numpy as np


def domb_sykes(a, step=1):
    """Ratios r_k = a_k / a_{k-step}; fit r_k = (1/R^step)(1 + c/k).  Returns (R, c)."""
    a = np.asarray(a, float)
    ks = np.arange(step, len(a))
    m = (a[ks - step] != 0) & (a[ks] != 0)
    k = ks[m]
    r = a[k] / a[k - step]
    if len(k) < 2:
        return np.nan, np.nan
    p = np.polyfit(1.0 / k, r, 1)
    inv = p[1]
    R = abs(inv) ** (-1.0 / step)
    return float(R), float(p[0] / p[1])


def pade(a, L, M):
    a = np.asarray(a, float)
    A = np.zeros((M, M))
    rhs = np.zeros(M)
    for i in range(M):
        for j in range(M):
            idx = L + i - j
            A[i, j] = a[idx] if 0 <= idx < len(a) else 0.0
        rhs[i] = -a[L + 1 + i]
    q = np.concatenate([[1.0], np.linalg.solve(A, rhs)]) if M else np.array([1.0])
    p = np.array([sum(q[j] * a[i - j] for j in range(min(i, M) + 1)) for i in range(L + 1)])
    return p, q


def pade_eval(p, q, x):
    x = np.asarray(x, dtype=complex)
    return np.polyval(p[::-1], x) / np.polyval(q[::-1], x)


def conformal_map_sum(a, R, x):
    """Re-expand f(x) in w with x = 4 R w/(1 - w)^2 (maps the cut plane to the unit disk)."""
    a = np.asarray(a, float)
    N = len(a)
    # coefficients of x^k as series in w up to order N-1
    b = np.zeros(N)
    xw = np.zeros(N); xw[1:] = [4 * R * n for n in range(1, N)]   # x(w) = 4R sum n w^n
    pw = np.zeros(N); pw[0] = 1.0
    for k in range(N):
        b += a[k] * pw
        pw = np.convolve(pw, xw)[:N]
    x = np.asarray(x, float)
    # invert x = 4R w/(1-w)^2:  w = (sqrt(1+x/R) - 1)/(sqrt(1+x/R) + 1)
    s = np.sqrt(1 + x / R + 0j)
    w = (s - 1) / (s + 1)
    return np.real(sum(b[k] * w**k for k in range(N)))


def misleading_convergence_check(terms, xis=(1.0, 1.5, 2.0)):
    """terms[n] = Sigma_n (n-th order contribution).  Returns {xi: partial sums}; a series
    converging at xi = 1 but diverging for all xi > 1 is flagged as suspicious."""
    terms = np.asarray(terms)
    out = {}
    for xi in xis:
        out[xi] = np.cumsum([t * xi**n for n, t in enumerate(terms)], axis=0)
    return out
