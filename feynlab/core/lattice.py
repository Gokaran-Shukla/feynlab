"""Densities of states, lattice Hilbert transforms and analytic continuation."""
from __future__ import annotations

import numpy as np


# ----------------------------------------------------------------------------
# Bethe lattice (semicircular DOS, half bandwidth D)
# ----------------------------------------------------------------------------
def bethe_dos(e, D: float = 1.0):
    e = np.asarray(e, dtype=float)
    out = np.zeros_like(e)
    m = np.abs(e) < D
    out[m] = 2.0 / (np.pi * D**2) * np.sqrt(D**2 - e[m] ** 2)
    return out


def bethe_hilbert(z, D: float = 1.0):
    """G(z) = int de rho(e)/(z-e) = 2(z - s sqrt(z^2-D^2))/D^2, branch with Im G<=0 for Im z>0."""
    z = np.asarray(z, dtype=complex)
    sq = np.sqrt(z**2 - D**2)
    # choose branch so that G ~ 1/z at large |z| and Im G has opposite sign to Im z
    sq = np.where((sq.imag * z.imag) < 0, -sq, sq)
    sq = np.where((z.imag == 0) & (np.real(z) * np.real(sq) < 0), -sq, sq)
    return 2.0 * (z - sq) / D**2


# ----------------------------------------------------------------------------
# Square / cubic tight-binding
# ----------------------------------------------------------------------------
def square_lattice_eps(nk: int, t: float = 1.0, tp: float = 0.0):
    k = 2 * np.pi * np.arange(nk) / nk
    kx, ky = np.meshgrid(k, k, indexing="ij")
    return -2 * t * (np.cos(kx) + np.cos(ky)) - 4 * tp * np.cos(kx) * np.cos(ky)


def lattice_hilbert(z, eps_k):
    """G_loc(z) = mean_k 1/(z - eps_k)."""
    z = np.asarray(z, dtype=complex)
    e = eps_k.ravel()
    return np.array([np.mean(1.0 / (zz - e)) for zz in np.atleast_1d(z)])


# ----------------------------------------------------------------------------
# Pade analytic continuation (Vidberg & Serene, J. Low Temp. Phys. 29, 179 (1977))
# ----------------------------------------------------------------------------
def pade_coefficients(z, u):
    z = np.asarray(z, dtype=np.clongdouble)
    u = np.asarray(u, dtype=np.clongdouble)
    N = len(z)
    g = np.zeros((N, N), dtype=np.clongdouble)
    g[0] = u
    for p in range(1, N):
        g[p, p:] = (g[p - 1, p - 1] - g[p - 1, p:]) / ((z[p:] - z[p - 1]) * g[p - 1, p:])
    return np.diag(g)


def pade_eval(zeval, z, a):
    zeval = np.asarray(zeval, dtype=np.clongdouble)
    N = len(a)
    A0, A1 = np.zeros_like(zeval), a[0] * np.ones_like(zeval)
    B0, B1 = np.ones_like(zeval), np.ones_like(zeval)
    for n in range(1, N):
        A2 = A1 + (zeval - z[n - 1]) * a[n] * A0
        B2 = B1 + (zeval - z[n - 1]) * a[n] * B0
        A0, A1, B0, B1 = A1, A2, B1, B2
    return np.asarray(A1 / B1, dtype=complex)


def pade_continue(iw, g_iw, omega, eta: float = 0.01, n_use: int = 64):
    """Continue G(iw_n) (first n_use frequencies) to G(omega + i eta)."""
    z = iw[:n_use]
    a = pade_coefficients(z, g_iw[:n_use])
    return pade_eval(np.asarray(omega) + 1j * eta, z, a)
