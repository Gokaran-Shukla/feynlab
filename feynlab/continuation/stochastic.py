"""Stochastic analytic continuation (SAC) by spectrum averaging.

Sandvik, PRB 57, 10287 (1998); Beach, arXiv:cond-mat/0403055; Shao & Sandvik, Phys. Rep. 1003, 1 (2023).

A(w) = sum_i a delta(w - w_i) with N_delta equal-amplitude delta functions at *continuous* positions
(a kernel table on a very fine frequency grid is used for evaluation).  Configurations are sampled
with weight exp(-chi^2/(2 Theta)) by single-delta and pair moves whose step sizes adapt to keep the
acceptance near 40%.  Theta is lowered by simulated annealing (the configuration is carried over),
and the final spectrum is the average at the Theta that satisfies the Shao-Sandvik criterion
    <chi^2>(Theta) = chi^2_min + a sqrt(2 chi^2_min),   a = 0.5,
sampled with many sweeps.  The averaged spectrum is histogrammed on a fine grid and smoothed by a
Gaussian whose width is the output grid spacing.
"""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _sweeps(x, Ktab, wmin, dwt, Gr, amp, Theta, n_sweeps, step, seed, accumulate, hist, hmin, hdw):
    np.random.seed(seed)
    N = len(x)
    nt = Ktab.shape[1]
    nd = Ktab.shape[0]
    idx = np.empty(N, np.int64)
    model = np.zeros(nd)
    for i in range(N):
        k = int((x[i] - wmin) / dwt + 0.5)
        idx[i] = k
        model += amp * Ktab[:, k]
    r = model - Gr
    chi2 = (r * r).sum()
    acc1 = 0
    tr1 = 0
    acc2 = 0
    tr2 = 0
    chi_acc = 0.0
    nmeas = 0
    for s in range(n_sweeps):
        for _ in range(N):
            if np.random.random() < 0.5:
                i = np.random.randint(N)
                xn = x[i] + (2 * np.random.random() - 1) * step[0]
                k = int((xn - wmin) / dwt + 0.5)
                if k < 0 or k >= nt:
                    continue
                tr1 += 1
                r2 = r + amp * (Ktab[:, k] - Ktab[:, idx[i]])
                c2 = (r2 * r2).sum()
                if c2 <= chi2 or np.random.random() < np.exp(-(c2 - chi2) / (2 * Theta)):
                    x[i] = xn; idx[i] = k; r = r2; chi2 = c2; acc1 += 1
            else:
                i = np.random.randint(N)
                j = np.random.randint(N)
                if i == j:
                    continue
                d = (2 * np.random.random() - 1) * step[1]
                xi, xj = x[i] + d, x[j] - d          # conserves the first moment
                ki = int((xi - wmin) / dwt + 0.5)
                kj = int((xj - wmin) / dwt + 0.5)
                if ki < 0 or ki >= nt or kj < 0 or kj >= nt:
                    continue
                tr2 += 1
                r2 = r + amp * (Ktab[:, ki] - Ktab[:, idx[i]] + Ktab[:, kj] - Ktab[:, idx[j]])
                c2 = (r2 * r2).sum()
                if c2 <= chi2 or np.random.random() < np.exp(-(c2 - chi2) / (2 * Theta)):
                    x[i] = xi; x[j] = xj; idx[i] = ki; idx[j] = kj; r = r2; chi2 = c2; acc2 += 1
        # adapt step sizes towards ~40% acceptance
        if tr1 > 0:
            step[0] *= 1.25 if acc1 / tr1 > 0.45 else (0.8 if acc1 / tr1 < 0.35 else 1.0)
        if tr2 > 0:
            step[1] *= 1.25 if acc2 / tr2 > 0.45 else (0.8 if acc2 / tr2 < 0.35 else 1.0)
        step[0] = min(max(step[0], 2 * dwt), 4.0)
        step[1] = min(max(step[1], 2 * dwt), 4.0)
        acc1 = 0; tr1 = 0; acc2 = 0; tr2 = 0
        if s >= n_sweeps // 5:
            chi_acc += chi2
            nmeas += 1
            if accumulate:
                for i in range(N):
                    b = int((x[i] - hmin) / hdw)
                    if 0 <= b < len(hist):
                        hist[b] += amp
    return chi_acc / max(nmeas, 1), nmeas


def stochastic_continuation(iw, G, sigma, omega, n_delta=1000, thetas=None, n_sweeps=400, n_final=4000,
                            a_criterion=0.5, seed=1, n_table=20001, fine=8):
    """Returns dict(A on 'omega', Theta, chi2, chi2_per_point, all=[(Theta, <chi2>)])."""
    iw = np.asarray(iw); G = np.asarray(G); w = np.asarray(omega, float)
    wmin, wmax = w[0], w[-1]
    wt = np.linspace(wmin, wmax, n_table)
    dwt = wt[1] - wt[0]
    K = 1.0 / (iw[:, None] - wt[None, :])
    sig = np.broadcast_to(np.asarray(sigma, float), G.shape)
    Ktab = np.ascontiguousarray(np.vstack([K.real / sig[:, None], K.imag / sig[:, None]]))
    Gr = np.concatenate([G.real / sig, G.imag / sig])
    rng = np.random.default_rng(seed)
    x = rng.uniform(wmin + 0.25 * (wmax - wmin), wmax - 0.25 * (wmax - wmin), n_delta)
    amp = 1.0 / n_delta
    step = np.array([0.5, 0.5])
    thetas = np.logspace(3, -2, 26) if thetas is None else np.asarray(thetas)
    nh = fine * (len(w) - 1)
    hdw = (wmax - wmin) / nh
    dummy = np.zeros(1)
    res, confs = [], []
    for k, th in enumerate(thetas):
        c, _ = _sweeps(x, Ktab, wmin, dwt, Gr, amp, th, n_sweeps, step, seed + k, False, dummy, wmin, hdw)
        res.append((th, c))
        confs.append((x.copy(), step.copy()))
    chi_min = min(c for _, c in res)
    target = chi_min + a_criterion * np.sqrt(2 * chi_min)
    # largest Theta (most entropy) whose <chi2> is below the target
    ok = [i for i, (_, c) in enumerate(res) if c <= target]
    i = ok[0] if ok else int(np.argmin([c for _, c in res]))
    th = res[i][0]
    x, step = confs[i][0].copy(), confs[i][1].copy()
    hist = np.zeros(nh)
    c, nmeas = _sweeps(x, Ktab, wmin, dwt, Gr, amp, th, n_final, step, seed + 999, True, hist, wmin, hdw)
    hc = wmin + (np.arange(nh) + 0.5) * hdw
    Ah = hist / max(nmeas, 1) / hdw
    # Gaussian smoothing with the output spacing
    s = (w[1] - w[0]) if len(w) > 1 else hdw
    A = np.array([np.sum(Ah * np.exp(-0.5 * ((hc - wi) / s) ** 2)) for wi in w]) * hdw / (np.sqrt(2 * np.pi) * s)
    return dict(A=A, Theta=float(th), chi2=float(c), chi2_per_point=float(c / len(Gr)), all=res,
                omega_fine=hc, A_fine=Ah)
