"""Stochastic analytic continuation by spectrum averaging (Sandvik, PRB 57, 10287 (1998);
Beach, arXiv:cond-mat/0403055).  A(w) is represented by K delta functions on a fine grid;
configurations are sampled with weight exp(-chi^2 / (2 Theta)); the continuation is the
average spectrum at the 'temperature' Theta where <chi^2> = chi2_min + sqrt(2 N_data)."""
from __future__ import annotations

import numpy as np
from numba import njit


@njit(cache=True)
def _sweep(pos, amp, Kr, Gr, Theta, n_sweeps, seed, nw):
    np.random.seed(seed)
    K = len(pos)
    model = Kr[:, pos] @ amp
    r = model - Gr
    chi2 = (r * r).sum()
    acc = np.zeros(nw)
    nacc = 0
    for s in range(n_sweeps):
        for _ in range(K):
            i = np.random.randint(K)
            if np.random.random() < 0.5:
                # move one delta
                new = pos[i] + np.random.randint(-25, 26)
                if new < 0 or new >= nw:
                    continue
                dr = amp[i] * (Kr[:, new] - Kr[:, pos[i]])
                r2 = r + dr
                c2 = (r2 * r2).sum()
                if np.random.random() < np.exp(-(c2 - chi2) / (2 * Theta)):
                    pos[i] = new; r = r2; chi2 = c2
            else:
                # transfer weight between two deltas (conserves normalisation)
                j = np.random.randint(K)
                if i == j:
                    continue
                d = (np.random.random() - 0.5) * 0.5 * min(amp[i], amp[j]) * 2
                if amp[i] + d <= 0 or amp[j] - d <= 0:
                    continue
                dr = d * (Kr[:, pos[i]] - Kr[:, pos[j]])
                r2 = r + dr
                c2 = (r2 * r2).sum()
                if np.random.random() < np.exp(-(c2 - chi2) / (2 * Theta)):
                    amp[i] += d; amp[j] -= d; r = r2; chi2 = c2
        if s > n_sweeps // 4:
            for k in range(K):
                acc[pos[k]] += amp[k]
            nacc += 1
    return acc / max(nacc, 1), chi2


def stochastic_continuation(iw, G, sigma, omega, n_delta=4000, thetas=None, n_sweeps=1500, seed=1):
    iw = np.asarray(iw); w = np.asarray(omega)
    dw = np.gradient(w)
    K = 1.0 / (iw[:, None] - w[None, :])
    Kr = np.ascontiguousarray(np.vstack([K.real, K.imag]) / sigma)
    Gr = np.concatenate([G.real, G.imag]) / sigma
    thetas = np.logspace(2, -1, 10) if thetas is None else thetas
    rng = np.random.default_rng(seed)
    pos = rng.integers(len(w) // 4, 3 * len(w) // 4, n_delta)
    amp = np.full(n_delta, 1.0 / n_delta)
    res = []
    for th in thetas:
        A, chi2 = _sweep(pos, amp, Kr, Gr, th, n_sweeps, seed, len(w))
        res.append((th, A / dw, chi2))
    chi_min = min(c for _, _, c in res)
    target = chi_min + np.sqrt(2 * len(Gr))
    i = int(np.argmin([abs(c - target) for _, _, c in res]))
    th, A, chi2 = res[i]
    return dict(A=A, Theta=th, chi2=chi2, chi2_per_point=chi2 / len(Gr), all=res)
