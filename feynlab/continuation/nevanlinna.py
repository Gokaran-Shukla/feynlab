"""Nevanlinna analytic continuation (Fei, Yeh & Gull, PRL 126, 056402 (2021)).

-G(z) is a Nevanlinna function (maps the upper half-plane into itself).  With the Moebius
map h(w) = (w - i)/(w + i) the data become Schur-class values theta_n = h(-G(iw_n)) in the
unit disk; the Schur algorithm interpolates them exactly, and the free function
theta_{M+1}(z) is set to 0 (no Hardy-basis optimisation here).  Pick's criterion selects
the largest set of points that admits a Nevanlinna interpolant.  Arithmetic in mpmath.
"""
from __future__ import annotations

import numpy as np

try:
    import mpmath as mp
except ImportError:  # pragma: no cover
    mp = None


def _mob(w):
    return (w - 1j) / (w + 1j)


def pick_ok(Y, theta, tol=1e-12):
    n = len(Y)
    P = np.empty((n, n), complex)
    for i in range(n):
        for j in range(n):
            P[i, j] = (1 - theta[i] * np.conj(theta[j])) / (1 - _mob(Y[i]) * np.conj(_mob(Y[j])))
    return np.linalg.eigvalsh(0.5 * (P + P.conj().T)).min() > -tol


class Nevanlinna:
    def __init__(self, iw, G, n_max=None, dps=64):
        if mp is None:
            raise ImportError("mpmath is required")
        mp.mp.dps = dps
        iw = np.asarray(iw); G = np.asarray(G)
        theta = _mob(-G)
        # Pick criterion: largest prefix admitting an interpolant
        n = len(iw) if n_max is None else min(n_max, len(iw))
        while n > 2 and not pick_ok(iw[:n], theta[:n]):
            n -= 1
        self.n_used = n
        self.Y = [mp.mpc(complex(y)) for y in iw[:n]]
        th = [mp.mpc(complex(t)) for t in theta[:n]]
        self.phis = []
        for j in range(n):
            abcd = mp.eye(2)
            for k in range(j):
                f = (self.Y[j] - self.Y[k]) / (self.Y[j] - mp.conj(self.Y[k]))
                abcd = abcd * mp.matrix([[f, self.phis[k]], [mp.conj(self.phis[k]) * f, 1]])
            num = -abcd[1, 1] * th[j] + abcd[0, 1]
            den = abcd[1, 0] * th[j] - abcd[0, 0]
            self.phis.append(num / den)

    def __call__(self, z):
        out = []
        for zz in np.atleast_1d(z):
            zz = mp.mpc(complex(zz))
            abcd = mp.eye(2)
            for k in range(self.n_used):
                f = (zz - self.Y[k]) / (zz - mp.conj(self.Y[k]))
                abcd = abcd * mp.matrix([[f, self.phis[k]], [mp.conj(self.phis[k]) * f, 1]])
            theta = abcd[0, 1] / abcd[1, 1]          # theta_{M+1} = 0
            nev = 1j * (1 + theta) / (1 - theta)
            out.append(complex(-nev))
        return np.array(out)

    def spectral(self, omega, eta=0.01):
        return -self(np.asarray(omega) + 1j * eta).imag / np.pi


def nevanlinna_auto(iw, G, omega, eta=0.02, n_max=None, tol_neg=1e-3, tol_norm=0.05):
    """Largest number of Matsubara points (<= n_max) whose interpolant gives a non-negative,
    normalised spectrum on the requested real grid (guards against loss of precision)."""
    n = len(iw) if n_max is None else min(n_max, len(iw))
    w = np.asarray(omega)
    dw = np.gradient(w)
    while n >= 4:
        nv = Nevanlinna(iw, G, n_max=n)
        A = nv.spectral(w, eta)
        if A.min() > -tol_neg and abs(np.sum(A * dw) - 1) < tol_norm:
            return dict(A=A, n_used=nv.n_used, model=nv)
        n = int(n * 0.8)
    raise RuntimeError("no Nevanlinna-compatible subset found")
