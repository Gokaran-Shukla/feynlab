"""Connected-determinant diagrammatic Monte Carlo (CDet).

R. Rossi, PRL 119, 045701 (2017): the sum of all *connected* Feynman diagrams with a given
set of vertices is obtained from determinants by the recursion

    C(r u S) = D(r u S) - sum_{T subset S, T != S} C(r u T) * D(S \\ T),

where D(V) = prod_sigma det[G0_sigma(tau_i - tau_j) - alpha delta_ij] is the sum of all
(connected and disconnected) diagrams and r the external ('root') vertex.  The cost is
O(3^n) per configuration instead of the factorial number of diagrams, and the
disconnected (volume-divergent) pieces cancel exactly -- the basis of DiagMC directly in
the thermodynamic limit.

Here: perturbative coefficients of the local observable
    X = <(n_up - 1/2)(n_dn - 1/2)>  =  sum_k a_k U^k,
    a_k = (-1)^k / k!  int_0^beta dtau_1..dtau_k  C(r; tau_1..tau_k)
for an impurity with Hartree-shifted bath propagator G0 (any hybridization).  The
vertices are U (n_up - 1/2)(n_dn - 1/2).  Benchmark: Hubbard atom,
X = -(1/4) tanh(beta U / 4) exactly.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid


@njit(cache=True)
def _g(t, beta, tab, n):
    s = 1.0
    if t < 0.0:
        t += beta
        s = -1.0
    x = t / beta * (n - 1)
    i = int(x)
    if i >= n - 1:
        i = n - 2
    f = x - i
    return s * ((1 - f) * tab[i] + f * tab[i + 1])


@njit(cache=True)
def _det_subset(mask, times, beta, tab_u, tab_d, ntab, diag_u, diag_d):
    idx = np.empty(times.shape[0], dtype=np.int64)
    m = 0
    for i in range(times.shape[0]):
        if (mask >> i) & 1:
            idx[m] = i
            m += 1
    if m == 0:
        return 1.0
    Mu = np.empty((m, m)); Md = np.empty((m, m))
    for a in range(m):
        for b in range(m):
            if a == b:
                Mu[a, b] = diag_u; Md[a, b] = diag_d
            else:
                dt = times[idx[a]] - times[idx[b]]
                Mu[a, b] = _g(dt, beta, tab_u, ntab)
                Md[a, b] = _g(dt, beta, tab_d, ntab)
    return np.linalg.det(Mu) * np.linalg.det(Md)


@njit(cache=True)
def _sample(order, n_samples, beta, tab_u, tab_d, ntab, diag_u, diag_d, seed):
    np.random.seed(seed)
    nv = order + 1                      # vertex 0 = root at tau = 0
    full = (1 << nv) - 1
    D = np.empty(1 << nv)
    C = np.empty(1 << nv)
    times = np.zeros(nv)
    acc = 0.0
    acc2 = 0.0
    for s in range(n_samples):
        for i in range(1, nv):
            times[i] = np.random.random() * beta
        for mask in range(1 << nv):
            D[mask] = _det_subset(mask, times, beta, tab_u, tab_d, ntab, diag_u, diag_d)
        # C only for masks containing the root (bit 0)
        for mask in range(1, 1 << nv, 2):
            S = mask & ~1
            val = D[mask]
            if S == 0:                  # root alone: no proper subsets
                C[mask] = val
                continue
            # proper subsets T of S: iterate sub = (sub - 1) & S, excluding S itself
            sub = (S - 1) & S
            while True:
                val -= C[sub | 1] * D[S & ~sub]
                if sub == 0:
                    break
                sub = (sub - 1) & S
            C[mask] = val
        c = C[full]
        acc += c
        acc2 += c * c
    mean = acc / n_samples
    err = np.sqrt(max(acc2 / n_samples - mean * mean, 0.0) / n_samples)
    return mean, err


def cdet_series(grid: MatsubaraGrid, g0_iw, max_order=5, n_samples=200000, seed=3, alpha=0.5, n_tab=20001):
    """Coefficients a_k (k = 0..max_order) of X(U) = sum a_k U^k with error bars."""
    from math import factorial
    tau = np.linspace(0, grid.beta, n_tab)
    c = MatsubaraGrid.tail_coefficients(grid.iw, g0_iw, n_fit=min(60, grid.n_iw // 4))
    tab = grid.iw_to_tau(g0_iw, 1.0, c[1], c[2], tau=tau)
    n0 = -tab[-1]
    diag = n0 - alpha
    coef, err = [], []
    for k in range(0, max_order + 1):
        m, e = _sample(k, n_samples if k > 0 else 1, grid.beta, tab, tab, n_tab, diag, diag, seed + k)
        pref = (-1) ** k * grid.beta**k / factorial(k)
        coef.append(pref * m)
        err.append(abs(pref) * e)
    return np.array(coef), np.array(err)


def atom_exact_coefficients(beta, max_order):
    """Taylor coefficients of -(1/4) tanh(beta U/4) in U."""
    # tanh x = x - x^3/3 + 2x^5/15 - 17x^7/315 + 62 x^9/2835
    t = {1: 1.0, 3: -1 / 3, 5: 2 / 15, 7: -17 / 315, 9: 62 / 2835}
    return np.array([-0.25 * t.get(k, 0.0) * (beta / 4) ** k for k in range(max_order + 1)])
