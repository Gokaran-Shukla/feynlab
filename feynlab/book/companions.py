"""Small calculators that reproduce worked examples of Mattuck's
'A Guide to Feynman Diagrams in the Many-Body Problem' (2nd ed., Dover).

ch2  pinball_propagator          classical pinball / 'monkey + lion' partial sum (Ch. 2.3)
ch3  quantum_pinball             T-matrix sum for one level coupled to scatterers (Ch. 3.2-3.3)
ch3  impurity_self_energy        energy shift & lifetime in a random impurity metal (Ch. 3.5)
ch4  geometric_dyson             Dyson sum G = G0 / (1 - G0 Sigma) vs truncated series
ch10 friedel_oscillation         static Lindhard screening, 2 kF oscillations (Ch. 10.7)
ch18 kondo_log                   2nd-order Kondo log and T_K (Ch. 18.2)
"""
from __future__ import annotations

import numpy as np


# --------------------------------------------------------------- Chapter 2
def pinball_propagator(P0: float, p_scatter: float, n_max: int = 50):
    """Probability to go from start to finish when each of the scatterers adds a
    detour of relative weight p (geometric series, Mattuck Eq. 2.25-type sum).
    Returns (exact partial sum P0/(1-p), list of partial sums by order)."""
    partial = np.cumsum([P0 * p_scatter**n for n in range(n_max)])
    return P0 / (1 - p_scatter), partial


# --------------------------------------------------------------- Chapter 3
def quantum_pinball(omega, eps0: float, V: float, eta: float = 1e-3, n_order: int | None = None):
    """Single level eps0 with a constant scattering amplitude V (diagonal):
    G = G0 + G0 V G0 + ... = 1/(w - eps0 - V + i eta).  n_order truncates the series."""
    w = np.asarray(omega) + 1j * eta
    G0 = 1.0 / (w - eps0)
    if n_order is None:
        return 1.0 / (1.0 / G0 - V)
    return sum(G0 * (V * G0) ** n for n in range(n_order + 1))


def impurity_self_energy(k, n_imp: float, W0: float, m: float = 1.0):
    """Born-approximation self-energy of an electron in randomly placed delta-
    function impurities (3D, Mattuck Ch. 3.5):
        Re Sigma = n W0   (first order: average potential),
        1/tau    = 2 pi n W0^2 rho(e_k),  rho(e) = m^{3/2} sqrt(2 e)/(2 pi^2)  (per spin).
    Returns (quasiparticle energy, inverse lifetime)."""
    e = k**2 / (2 * m)
    rho = m**1.5 * np.sqrt(2 * e) / (2 * np.pi**2)
    return e + n_imp * W0, 2 * np.pi * n_imp * W0**2 * rho


# --------------------------------------------------------------- Chapter 4/10
def geometric_dyson(G0: complex, Sigma: complex, n_terms: int):
    """Partial sums of G0 + G0 S G0 + ... and the Dyson resummation."""
    terms = np.cumsum([G0 * (Sigma * G0) ** n for n in range(n_terms)])
    return G0 / (1 - G0 * Sigma), terms


# --------------------------------------------------------------- Chapter 10.7
def friedel_oscillation(r, rs: float = 2.0, Z: float = 1.0, nq: int = 20000, qmax_kF: float = 12.0):
    """Induced density around a point charge Z screened by the static RPA dielectric
    function; returns delta n(r).  Asymptotically delta n ~ cos(2 kF r)/r^3."""
    from ..egas.rpa import kF_of_rs, lindhard_iw
    kF = kF_of_rs(rs)
    q = np.linspace(1e-4, qmax_kF * kF, nq)
    vq = 4 * np.pi / q**2
    P = lindhard_iw(q, np.zeros_like(q) + 1e-12, kF)
    chi = -P / (1 + vq * P)  # density response
    dn_q = chi * Z * vq      # external potential of +Z charge: -Z v_q acting on electrons -> sign
    r = np.atleast_1d(r)
    # 3D Fourier transform of radial function
    out = np.array([np.trapezoid(dn_q * q * np.sin(q * rr), q) / (2 * np.pi**2 * rr) for rr in r])
    return -out


# --------------------------------------------------------------- Chapter 18
def kondo_log(T, J_rho: float, D: float = 1.0):
    """Second-order (Born + Kondo) resistivity factor  1 + 2 J rho ln(D/T)  (J<0 AF in
    Mattuck's sign convention, here J_rho>0 means antiferromagnetic) and the Kondo
    temperature T_K = D exp(-1/(2 J rho)) where perturbation theory breaks down."""
    T = np.asarray(T, float)
    factor = 1.0 + 2 * J_rho * np.log(D / T)
    T_K = D * np.exp(-1.0 / (2 * J_rho))
    return factor, T_K
