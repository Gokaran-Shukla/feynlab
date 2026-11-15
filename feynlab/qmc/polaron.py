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
def _run(alpha, mu, tau_max, n_steps, n_bins, seed, nmax, p=0.0):
    np.random.seed(seed)
    t1 = np.zeros(nmax)
    t2 = np.zeros(nmax)
    q = np.zeros((nmax, 3))
    n = 0
    tau = 0.5
    hist = np.zeros(n_bins)
    hist0 = np.zeros(n_bins)
    order_hist = np.zeros(nmax)
    hx2 = np.zeros(n_bins)          # accumulates X.X, X = sum_arcs q (t2 - t1)  (m* estimator)
    X = np.zeros(3)
    e_p = 0.5 * p * p
    pref = np.sqrt(2.0) * alpha / np.pi   # (alpha/(2 sqrt2 pi^2)) * 4 pi
    for step in range(n_steps):
        r = np.random.random()
        if r < 0.2:
            # ---- resample diagram length (exact heat-bath on the last free segment)
            tlast = 0.0
            for i in range(n):
                if t2[i] > tlast:
                    tlast = t2[i]
            new = tlast - np.log(np.random.random()) / (e_p - mu)
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
                dS = 0.5 * qm * qm * dt - p * qz * dt     # external momentum p along z
                for j in range(n):
                    ov = _overlap(ta, tb, t1[j], t2[j])
                    if ov > 0:
                        dS += (qx * q[j, 0] + qy * q[j, 1] + qz * q[j, 2]) * ov
                R = pref * tau * np.exp(-dS) / (Pq * (n + 1))
                if np.random.random() < R:
                    t1[n], t2[n] = ta, tb
                    q[n, 0], q[n, 1], q[n, 2] = qx, qy, qz
                    n += 1
                    X[0] += qx * dt
                    X[1] += qy * dt
                    X[2] += qz * dt
        else:
            # ---- remove an arc
            if n > 0:
                i = np.random.randint(n)
                dt = t2[i] - t1[i]
                qx, qy, qz = q[i, 0], q[i, 1], q[i, 2]
                qm = np.sqrt(qx * qx + qy * qy + qz * qz)
                sig = 1.0 / np.sqrt(dt)
                Pq = 2.0 / (sig * np.sqrt(2 * np.pi)) * np.exp(-qm * qm / (2 * sig * sig))
                dS = 0.5 * qm * qm * dt - p * qz * dt   # action of the arc relative to configuration without it
                for j in range(n):
                    if j != i:
                        ov = _overlap(t1[i], t2[i], t1[j], t2[j])
                        if ov > 0:
                            dS += (qx * q[j, 0] + qy * q[j, 1] + qz * q[j, 2]) * ov
                R = (Pq * n) / (pref * tau * np.exp(-dS))
                if np.random.random() < R:
                    X[0] -= qx * dt
                    X[1] -= qy * dt
                    X[2] -= qz * dt
                    n -= 1
                    t1[i], t2[i] = t1[n], t2[n]
                    q[i, 0], q[i, 1], q[i, 2] = q[n, 0], q[n, 1], q[n, 2]
        b = int(tau / tau_max * n_bins)
        if b < n_bins:
            hist[b] += 1.0
            if n == 0:
                hist0[b] += 1.0
            hx2[b] += X[0] * X[0] + X[1] * X[1] + X[2] * X[2]
        order_hist[n] += 1.0
    return hist, order_hist, hist0, hx2


def polaron_diagmc(alpha: float, mu: float | None = None, tau_max: float = 30.0,
                   n_steps: int = 50_000_000, n_bins: int = 300, n_chains: int = 4,
                   fit_window=(0.5, 0.95), seed: int = 11, p: float = 0.0):
    """Energy E(p), quasiparticle weight Z(p) from the tail G(p, tau) ~ Z exp[-(E - mu) tau].
    The absolute normalisation of the histogram is fixed by the order-0 diagram, whose weight
    exp[-(p^2/2 - mu) tau] is known exactly."""
    if mu is None:
        mu = feynman_energy(alpha) - 0.3
    tau_grid = (np.arange(n_bins) + 0.5) * tau_max / n_bins
    Es, Zs, hs, ords, Ms = [], [], [], [], []
    for c in range(n_chains):
        h, o, h0, x2 = _run(alpha, mu, tau_max, n_steps, n_bins, seed + 101 * c, 4000, p)
        hs.append(h)
        ords.append(o)
        lo, hi = int(fit_window[0] * n_bins), int(fit_window[1] * n_bins)
        m = h[lo:hi] > 0
        fit = np.polyfit(tau_grid[lo:hi][m], np.log(h[lo:hi][m]), 1, w=np.sqrt(h[lo:hi][m]))
        Es.append(mu - fit[0])
        # normalisation C: h0(tau) = C exp[-(p^2/2 - mu) tau]  (fit on the short-time region)
        m0 = h0 > 20
        C = np.exp(np.mean(np.log(h0[m0]) + (0.5 * p * p - mu) * tau_grid[m0])) if m0.any() else np.nan
        Zs.append(np.exp(fit[1]) / C)
        # m/m* = 1 - lim d<X^2>/dtau / 3   (G(p,tau)/G(0,tau) = e^{-p^2 tau/2} <e^{p.X}>)
        mm = h[lo:hi] > 50
        sl = np.polyfit(tau_grid[lo:hi][mm], x2[lo:hi][mm] / h[lo:hi][mm], 1, w=np.sqrt(h[lo:hi][mm]))[0]
        Ms.append(1.0 / (1.0 - sl / 3.0))
    Es, Zs, Ms = np.array(Es), np.array(Zs), np.array(Ms)
    o = np.sum(ords, axis=0)
    se = lambda x: float(x.std(ddof=1) / np.sqrt(len(x))) if len(x) > 1 else np.nan
    return dict(E=float(Es.mean()), E_err=se(Es), Z=float(Zs.mean()), Z_err=se(Zs),
                m_star=float(Ms.mean()), m_star_err=se(Ms),
                tau=tau_grid, G=np.sum(hs, axis=0), order_mean=float((np.arange(len(o)) * o).sum() / o.sum()),
                mu=mu, p=p)


def polaron_effective_mass(alpha, p=0.4, **kw):
    """m* from E(p) - E(0) = p^2 / (2 m*) (finite-difference route, O(p^4) bias).
    The preferred, unbiased estimator is polaron_diagmc(...)['m_star'] (momentum-fluctuation
    estimator m/m* = 1 - lim <X^2>/(3 tau), X = sum over phonon arcs of q (t2 - t1))."""
    r0 = polaron_diagmc(alpha, p=0.0, **kw)
    r1 = polaron_diagmc(alpha, p=p, **kw)
    dE = r1["E"] - r0["E"]
    err = np.hypot(r0["E_err"], r1["E_err"])
    m = p * p / (2 * dE)
    return dict(m_star=m, m_err=m * err / dE, E0=r0, Ep=r1)


def weak_coupling_mass(alpha):
    """m*/m = 1 + alpha/6 + 0.02362763 alpha^2 (Roseler 1968; Smondyrev 1986)."""
    return 1 + alpha / 6 + 0.02362763 * alpha**2


def weak_coupling_Z(alpha):
    """Z = 1 - alpha/2 + O(alpha^2) (first-order Froehlich)."""
    return 1 - alpha / 2


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
