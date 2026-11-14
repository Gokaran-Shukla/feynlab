"""Fixed-point accelerators for self-consistency loops (DMFT, GF2, FLEX, SCBA).

x_{n+1} = F(x_n):  residual r = F(x) - x.
* LinearMixer   : x <- x + a r
* AndersonMixer : Anderson / Pulay (DIIS) with history m (Walker & Ni, SIAM J. Numer. Anal. 49, 1715 (2011))
* BroydenMixer  : modified Broyden (Johnson, PRB 38, 12807 (1988))
All work on complex arrays of any shape (flattened internally, real/imag stacked).
"""
from __future__ import annotations

import numpy as np


def _flat(x):
    x = np.asarray(x)
    return np.concatenate([x.real.ravel(), x.imag.ravel()]) if np.iscomplexobj(x) else x.ravel().astype(float)


def _unflat(v, like):
    like = np.asarray(like)
    if np.iscomplexobj(like):
        n = like.size
        return (v[:n] + 1j * v[n:]).reshape(like.shape)
    return v.reshape(like.shape)


class LinearMixer:
    def __init__(self, alpha=0.5):
        self.alpha = alpha

    def __call__(self, x, fx):
        return x + self.alpha * (np.asarray(fx) - x)


class AndersonMixer:
    def __init__(self, alpha=0.5, m=6):
        self.alpha, self.m = alpha, m
        self.X, self.R = [], []

    def __call__(self, x, fx):
        xv, rv = _flat(x), _flat(fx) - _flat(x)
        self.X.append(xv); self.R.append(rv)
        self.X, self.R = self.X[-self.m - 1:], self.R[-self.m - 1:]
        if len(self.X) == 1:
            return _unflat(xv + self.alpha * rv, x)
        dX = np.array([self.X[i + 1] - self.X[i] for i in range(len(self.X) - 1)]).T
        dR = np.array([self.R[i + 1] - self.R[i] for i in range(len(self.R) - 1)]).T
        gamma, *_ = np.linalg.lstsq(dR, rv, rcond=None)
        new = xv + self.alpha * rv - (dX + self.alpha * dR) @ gamma
        return _unflat(new, x)


class BroydenMixer:
    """Modified Broyden with uniform weights."""

    def __init__(self, alpha=0.5, m=8, w0=0.01):
        self.alpha, self.m, self.w0 = alpha, m, w0
        self.dF, self.U = [], []
        self.last = None

    def __call__(self, x, fx):
        xv, F = _flat(x), _flat(fx) - _flat(x)
        if self.last is not None:
            xo, Fo = self.last
            dF = F - Fo
            nrm = np.linalg.norm(dF) or 1.0
            dF /= nrm
            dx = (xv - xo) / nrm
            self.dF.append(dF)
            self.U.append(self.alpha * dF + dx)
            self.dF, self.U = self.dF[-self.m:], self.U[-self.m:]
        self.last = (xv.copy(), F.copy())
        new = xv + self.alpha * F
        if self.dF:
            k = len(self.dF)
            A = np.array([[self.dF[i] @ self.dF[j] for j in range(k)] for i in range(k)]) + self.w0**2 * np.eye(k)
            c = np.array([dF @ F for dF in self.dF])
            gam = np.linalg.solve(A, c)
            new -= sum(g * u for g, u in zip(gam, self.U))
        return _unflat(new, x)


def solve_fixed_point(F, x0, mixer=None, tol=1e-10, maxiter=500):
    mixer = mixer or AndersonMixer()
    x = np.asarray(x0)
    for it in range(maxiter):
        fx = F(x)
        res = np.max(np.abs(fx - x))
        if res < tol:
            return fx, it, res
        x = mixer(x, fx)
    return x, maxiter, res
