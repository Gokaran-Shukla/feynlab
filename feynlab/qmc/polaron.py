"""Diagrammatic Monte Carlo for the 3D Froehlich polaron
(Prokof'ev & Svistunov, PRL 81, 2514 (1998); Mishchenko et al., PRB 62, 6317 (2000)).

Units hbar = m = omega_ph = 1.  H = p^2/2 + sum_q b_q^+ b_q + sum_q V_q (b_q^+ + b_-q) e^{iqr},
|V_q|^2 = 2 sqrt(2) pi alpha / q^2.

The Monte Carlo samples *all* Feynman diagrams of G(p=0, tau) (bare series, any
order): a configuration is the diagram length tau and a set of phonon arcs
(t1, t2, q).  Because |V_q|^2 d^3q ~ dq dOmega, a diagram's weight per arc is
(alpha / (2 sqrt 2 pi^2)) e^{-(t2-t1)} times the electron propagators
exp[-(k^2/2 - mu) dt] along the line.  The polaron energy follows from the
exponential decay G(0, tau) ~ Z exp[-(E - mu) tau].

Also provided: Feynman's variational energy (Phys. Rev. 97, 660 (1955)) and the
exact weak-coupling series E = -a - 0.0159196220 a^2 - 0.000806070048 a^3.
"""
from __future__ import annotations

import numpy as np
from numba import njit
from scipy import integrate, optimize


@njit(cache=True)
def _overlap(a1, a2, b1, b2):
    lo = a1 if a1 > b1 else b1
    hi = a2 if a2 < b2 else b2
    return hi - lo if hi > lo else 0.0


@njit(cache=True)
def _run(alpha, mu, tau_max, n_steps, n_bins, seed, nmax):
    np.random.seed(seed)
    t1 = np.zeros(nmax)
    t2 = np.zeros(nmax)
    q = np.zeros((nmax, 3))
    n = 0
    tau = 0.5
    hist = np.zeros(n_bins)
    order_hist = np.zeros(nmax)
    pref = np.sqrt(2.0) * alpha / np.pi   # (alpha/(2 sqrt2 pi^2)) * 4 pi
    for step in range(n_steps):
        r = np.random.random()
        if r < 0.2:
            # ---- resample diagram length (exact heat-bath on the last free segment)
            tlast = 0.0
            for i in range(n):
                if t2[i] > tlast:
                    tlast = t2[i]
            new = tlast - np.log(np.random.random()) / (-mu)
            if new < tau_max:
                tau = new
        elif r < 0.6:
            # ---- add an arc
            ta = np.random.random() * tau
            dt = -np.log(np.random.random())
            tb = ta + dt
            if tb < tau and n < nmax - 1:
                sig = 1.0 / np.sqrt(dt)
                # |q| from half-normal(sig): density P(q) = 2/(sig sqrt(2pi)) exp(-q^2/2sig^2)
                qm = abs(np.random.normal()) * sig
                Pq = 2.0 / (sig * np.sqrt(2 * np.pi)) * np.exp(-qm * qm / (2 * sig * sig))
                ct = 2 * np.random.random() - 1
                ph = 2 * np.pi * np.random.random()
                st = np.sqrt(1 - ct * ct)
                qx, qy, qz = qm * st * np.cos(ph), qm * st * np.sin(ph), qm * ct
                dS = 0.5 * qm * qm * dt
                for j in range(n):
                    ov = _overlap(ta, tb, t1[j], t2[j])
                    if ov > 0:
                        dS += (qx * q[j, 0] + qy * q[j, 1] + qz * q[j, 2]) * ov
                R = pref * tau * np.exp(-dS) / (Pq * (n + 1))
                if np.random.random() < R:
                    t1[n], t2[n] = ta, tb
                    q[n, 0], q[n, 1], q[n, 2] = qx, qy, qz
                    n += 1
        else:
            # ---- remove an arc
            if n > 0:
                i = np.random.randint(n)
                dt = t2[i] - t1[i]
                qx, qy, qz = q[i, 0], q[i, 1], q[i, 2]
                qm = np.sqrt(qx * qx + qy * qy + qz * qz)
                sig = 1.0 / np.sqrt(dt)
                Pq = 2.0 / (sig * np.sqrt(2 * np.pi)) * np.exp(-qm * qm / (2 * sig * sig))
                dS = 0.5 * qm * qm * dt   # action of the arc relative to configuration without it
                for j in range(n):
                    if j != i:
                        ov = _overlap(t1[i], t2[i], t1[j], t2[j])
                        if ov > 0:
                            dS += (qx * q[j, 0] + qy * q[j, 1] + qz * q[j, 2]) * ov
                R = (Pq * n) / (pref * tau * np.exp(-dS))
                if np.random.random() < R:
                    n -= 1
                    t1[i], t2[i] = t1[n], t2[n]
                    q[i, 0], q[i, 1], q[i, 2] = q[n, 0], q[n, 1], q[n, 2]
        b = int(tau / tau_max * n_bins)
        if b < n_bins:
            hist[b] += 1.0
        order_hist[n] += 1.0
    return hist, order_hist


def polaron_diagmc(alpha: float, mu: float | None = None, tau_max: float = 30.0,
                   n_steps: int = 50_000_000, n_bins: int = 300, n_chains: int = 4,
                   fit_window=(0.5, 0.95), seed: int = 11):
    """Ground-state energy from the exponential tail of G(p=0, tau).
    Returns dict(E, E_err, Z_rel, order_mean, tau, G)."""
    if mu is None:
        mu = feynman_energy(alpha) - 0.3
    tau_grid = (np.arange(n_bins) + 0.5) * tau_max / n_bins
    Es, hs, ords = [], [], []
    for c in range(n_chains):
        h, o = _run(alpha, mu, tau_max, n_steps, n_bins, seed + 101 * c, 4000)
        hs.append(h)
        ords.append(o)
        lo, hi = int(fit_window[0] * n_bins), int(fit_window[1] * n_bins)
        m = h[lo:hi] > 0
        p = np.polyfit(tau_grid[lo:hi][m], np.log(h[lo:hi][m]), 1, w=np.sqrt(h[lo:hi][m]))
        Es.append(mu - p[0])
    Es = np.array(Es)
    o = np.sum(ords, axis=0)
    return dict(E=float(Es.mean()), E_err=float(Es.std(ddof=1) / np.sqrt(len(Es))) if len(Es) > 1 else np.nan,
                tau=tau_grid, G=np.sum(hs, axis=0), order_mean=float((np.arange(len(o)) * o).sum() / o.sum()),
                mu=mu)


# ----------------------------------------------------------------------------
def weak_coupling_energy(alpha):
    """Exact series to O(alpha^3) (Smondyrev 1986; Höhler & Müllensiefen)."""
    return -alpha - 0.0159196220 * alpha**2 - 0.000806070048 * alpha**3


def feynman_energy(alpha: float) -> float:
    """Feynman's variational upper bound (Phys. Rev. 97, 660 (1955)):
    E(v,w) = 3(v-w)^2/(4v) - (alpha v/sqrt(pi)) int_0^inf e^{-t} / sqrt(w^2 t + (v^2-w^2)(1-e^{-vt})/v) dt."""
    def E(x):
        w = np.exp(x[0])
        v = w + np.exp(x[1])
        f = lambda t: np.exp(-t) / np.sqrt(w**2 * t + (v**2 - w**2) * (1 - np.exp(-v * t)) / v)
        # substitute t = s^2 to remove the 1/sqrt(t) singularity
        g = lambda s: 2 * s * f(s * s)
        I, _ = integrate.quad(g, 0, np.inf, limit=400)
        return 3 * (v - w) ** 2 / (4 * v) - alpha * v / np.sqrt(np.pi) * I
    best = None
    for x0 in ([np.log(3.0), np.log(0.5)], [np.log(1.0), np.log(alpha**2)], [np.log(2.0), np.log(3.0)]):
        r = optimize.minimize(E, x0, method="Nelder-Mead", options=dict(xatol=1e-10, fatol=1e-12, maxiter=4000))
        if best is None or r.fun < best.fun:
            best = r
    return float(best.fun)
