"""Statistical error analysis for Monte Carlo time series.

* binning_analysis : error of the mean vs bin size (plateau = correct error)
* gamma_method     : integrated autocorrelation time tau_int with automatic
                     windowing (Wolff, Comput. Phys. Commun. 156, 143 (2004))
* jackknife        : bias-corrected estimate and error of a nonlinear function
                     of means (e.g. a ratio of sign-weighted observables)
"""
from __future__ import annotations

import numpy as np


def binning_analysis(x, min_bins: int = 16):
    x = np.asarray(x, float)
    out = []
    b = 1
    while len(x) // b >= min_bins:
        n = len(x) // b
        xb = x[: n * b].reshape(n, b).mean(1)
        out.append((b, xb.std(ddof=1) / np.sqrt(n)))
        b *= 2
    return np.array(out)


def gamma_method(x, S: float = 1.5):
    """Returns (mean, error, tau_int, window W)."""
    x = np.asarray(x, float)
    N = len(x)
    d = x - x.mean()
    f = np.fft.rfft(d, n=2 * N)
    acf = np.fft.irfft(f * np.conj(f))[:N] / np.arange(N, 0, -1)
    gamma0 = acf[0]
    if gamma0 == 0:
        return x.mean(), 0.0, 0.5, 0
    rho = acf / gamma0
    tau = 0.5
    W = 1
    for W in range(1, N // 2):
        tau = 0.5 + rho[1:W + 1].sum()
        if tau <= 0:
            tau = 1e-9
        g = np.exp(-W / (S * tau / np.log((2 * tau + 1) / max(2 * tau - 1, 1e-12)))) - (S * tau) / np.sqrt(W * N) \
            if tau > 0.5 else -1
        if g < 0:
            break
    err = np.sqrt(2 * tau * gamma0 / N)
    return float(x.mean()), float(err), float(tau), int(W)


def jackknife(func, *series, n_blocks: int = 32):
    """func(*means) -> value. Returns (estimate, error)."""
    arrs = [np.asarray(s, float) for s in series]
    N = len(arrs[0])
    bs = N // n_blocks
    arrs = [a[: bs * n_blocks].reshape(n_blocks, bs) for a in arrs]
    tot = [a.sum(1) for a in arrs]
    full = func(*[t.sum() / (bs * n_blocks) for t in tot])
    jk = np.array([func(*[(t.sum() - t[i]) / (bs * (n_blocks - 1)) for t in tot]) for i in range(n_blocks)])
    est = n_blocks * full - (n_blocks - 1) * jk.mean()
    err = np.sqrt((n_blocks - 1) / n_blocks * np.sum((jk - jk.mean()) ** 2))
    return float(est), float(err)


def ar1_series(n, phi, seed=0):
    """AR(1) test process with exactly known tau_int = (1+phi)/(2(1-phi))."""
    rng = np.random.default_rng(seed)
    x = np.empty(n)
    x[0] = rng.standard_normal()
    e = rng.standard_normal(n) * np.sqrt(1 - phi**2)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    return x
