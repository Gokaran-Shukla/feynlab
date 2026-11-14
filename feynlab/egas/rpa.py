"""Homogeneous electron gas in the random-phase approximation (Mattuck Chs. 10.4,
12.3, 13.6).  Hartree atomic units, spin-unpolarised, density n = 3/(4 pi rs^3).

* lindhard_iw(q, nu)         : Pi0(q, i nu) >= 0, both spins (T = 0)
* lindhard_real(q, w)        : Re Pi0(q, w) on the real axis
* rpa_correlation_energy(rs) : ring-diagram sum  (Gell-Mann & Brueckner 1957)
      e_c = 1/(2n) int d^3q/(2pi)^3 int dnu/2pi [ln(1 + v_q Pi0) - v_q Pi0]
  High-density limit:  e_c -> (1 - ln 2)/pi^2 * ln rs + const  (0.0311 ln rs Ha)
* plasmon_dispersion(rs, q)  : zero of eps = 1 + v_q Re Pi0;
      w^2 = w_p^2 + (3/5) v_F^2 q^2 + ...   (Mattuck Ex. 13.6)
"""
from __future__ import annotations

import numpy as np
from scipy import integrate, optimize


def kF_of_rs(rs):
    return (9 * np.pi / 4) ** (1 / 3) / rs


def lindhard_iw(q, nu, kF):
    """Lindhard function at imaginary frequency, units of states/(Ha bohr^3)."""
    q = np.asarray(q, float)
    nu = np.asarray(nu, float)
    N0 = kF / np.pi**2  # m=1, both spins
    z = q / (2 * kF)
    u = nu / (q * kF)
    with np.errstate(divide="ignore", invalid="ignore"):
        t1 = (1 - z**2 + u**2) / (8 * z) * np.log(((z + 1) ** 2 + u**2) / ((z - 1) ** 2 + u**2))
        t2 = 0.5 * u * (np.arctan((1 + z) / u) + np.arctan((1 - z) / u))
    t2 = np.where(u == 0, 0.0, t2)
    res = N0 * (0.5 + t1 - t2)
    # large-u / small-z asymptotics (avoids catastrophic cancellation):
    # Pi0/N0 = 1/(3 u^2) - (1/(15 u^4))(3 + 5 z^2 u^2 ... ) ~ leading terms
    big = (u > 50 * (1 + z)) | ~np.isfinite(res)
    if np.any(big):
        ub, zb = np.broadcast_arrays(u, z)
        ub, zb = ub[big], zb[big]
        res = np.array(res, dtype=float)
        res[big] = N0 * (1.0 / (3 * ub**2) - (1.0 / 15 + zb**2 / 3) / ub**4 * 1.0)
    return res


def lindhard_real(q, w, kF):
    q = np.asarray(q, float)
    w = np.asarray(w, float)
    N0 = kF / np.pi**2
    z = q / (2 * kF)
    u = w / (q * kF)
    def L(a):
        return (1 - a**2) * np.log(np.abs((a + 1) / (a - 1)))
    return N0 * (0.5 + (L(z - u) + L(z + u)) / (8 * z))


def rpa_correlation_energy(rs: float, n_x: int = 400, n_y: int = 200) -> float:
    """RPA correlation energy per electron (Hartree): Gauss-Legendre quadrature on
    the maps q = 2 kF x, nu = q kF y with x = t/(1-t), y = s/(1-s)."""
    kF = kF_of_rs(rs)
    n = kF**3 / (3 * np.pi**2)
    t, wt = np.polynomial.legendre.leggauss(n_x)
    t = 0.5 * (t + 1); wt = 0.5 * wt
    s_, ws = np.polynomial.legendre.leggauss(n_y)
    s_ = 0.5 * (s_ + 1); ws = 0.5 * ws
    x = t / (1 - t); dx = wt / (1 - t) ** 2
    y = s_ / (1 - s_); dy = ws / (1 - s_) ** 2
    X, Yv = np.meshgrid(x, y, indexing="ij")
    q = 2 * kF * X
    nu = q * kF * Yv
    p = 4 * np.pi / q**2 * lindhard_iw(q, nu, kF)
    f = np.log1p(p) - p
    # int d^3q/(2pi)^3 = int 4 pi q^2 dq/(2pi)^3 ; int_{-inf}^{inf} dnu/2pi = 2 int_0^inf dnu/2pi
    integrand = f * (4 * np.pi * q**2) / (2 * np.pi) ** 3 * (2 * kF) * (q * kF) / np.pi
    tot = np.sum(integrand * dx[:, None] * dy[None, :])
    return float(tot / (2 * n))


def plasma_frequency(rs):
    return np.sqrt(3.0 / rs**3)


def plasmon_dispersion(rs: float, q: np.ndarray):
    kF = kF_of_rs(rs)
    out = []
    for qq in np.atleast_1d(q):
        vq = 4 * np.pi / qq**2
        f = lambda w: 1.0 + vq * lindhard_real(qq, w, kF)
        lo = qq * kF + qq**2 / 2 + 1e-9  # upper edge of particle-hole continuum
        hi = 10 * plasma_frequency(rs) + 10 * qq * kF
        if f(lo) * f(hi) > 0:
            out.append(np.nan)  # plasmon has entered the continuum (Landau damping)
            continue
        out.append(optimize.brentq(f, lo, hi, xtol=1e-14))
    return np.array(out)
