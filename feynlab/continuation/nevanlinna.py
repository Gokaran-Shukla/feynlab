"""Nevanlinna analytic continuation (Fei, Yeh & Gull, PRL 126, 056402 (2021)).

-G(z) is a Nevanlinna function (maps the upper half-plane into itself).  With the Moebius
map h(w) = (w - i)/(w + i) the data become Schur-class values theta_n = h(-G(iw_n)) in the
unit disk; the Schur algorithm interpolates them exactly.  The free function theta_{M+1}(z) is
either 0 or optimised in the Hardy basis f_k(z) = (z-i)^k / [sqrt(pi) (z+i)^{k+1}] to minimise
(1 - int A)^2 + lambda int (A'')^2 (Fei, Yeh, Gull, Sec. III; `hardy_optimize`).  Pick's criterion
(evaluated in multiprecision) selects the largest set of points that admits a Nevanlinna interpolant.
The working precision grows with the number of points (dps = max(64, 3 n)) because the Schur
recursion loses about one digit per point.
"""
from __future__ import annotations

import numpy as np

try:
    import mpmath as mp
except ImportError:  # pragma: no cover
    mp = None


def _mob(w):
    return (w - 1j) / (w + 1j)


def pick_ok(Y, theta, tol=1e-12, multiprecision=True):
    n = len(Y)
    if multiprecision and mp is not None:
        Ym = [mp.mpc(complex(y)) for y in Y]
        tm = [mp.mpc(complex(t)) for t in theta]
        hm = [(y - 1j) / (y + 1j) for y in Ym]
        P = mp.matrix(n, n)
        for i in range(n):
            for j in range(n):
                P[i, j] = (1 - tm[i] * mp.conj(tm[j])) / (1 - hm[i] * mp.conj(hm[j]))
        ev = mp.eighe(P, eigvals_only=True)
        return float(min(ev)) > -tol
    P = np.empty((n, n), complex)
    for i in range(n):
        for j in range(n):
            P[i, j] = (1 - theta[i] * np.conj(theta[j])) / (1 - _mob(Y[i]) * np.conj(_mob(Y[j])))
    return np.linalg.eigvalsh(0.5 * (P + P.conj().T)).min() > -tol


class Nevanlinna:
    def __init__(self, iw, G, n_max=None, dps=None):
        if mp is None:
            raise ImportError("mpmath is required")
        iw = np.asarray(iw); G = np.asarray(G)
        n = len(iw) if n_max is None else min(n_max, len(iw))
        mp.mp.dps = dps if dps is not None else max(64, 3 * n)
        theta = _mob(-G)
        # Pick criterion: largest prefix admitting an interpolant
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

    # ------------------------------------------------------------- Hardy-basis optimisation
    def _abcd(self, z):
        """Moebius coefficients theta(z) = (a t + b)/(c t + d), t = theta_{M+1}(z), normalised by d."""
        out = np.empty((len(z), 4), complex)
        for i, zz in enumerate(z):
            zz = mp.mpc(complex(zz))
            abcd = mp.eye(2)
            for k in range(self.n_used):
                f = (zz - self.Y[k]) / (zz - mp.conj(self.Y[k]))
                abcd = abcd * mp.matrix([[f, self.phis[k]], [mp.conj(self.phis[k]) * f, 1]])
            d = abcd[1, 1]
            out[i] = [complex(abcd[0, 0] / d), complex(abcd[0, 1] / d), complex(abcd[1, 0] / d), 1.0]
        return out

    @staticmethod
    def _hardy(z, H):
        k = np.arange(H)
        f = ((z[:, None] - 1j) / (z[:, None] + 1j)) ** k[None, :] / (np.sqrt(np.pi) * (z[:, None] + 1j))
        zc = np.conj(z)
        fs = np.conj(((zc[:, None] - 1j) / (zc[:, None] + 1j)) ** k[None, :] / (np.sqrt(np.pi) * (zc[:, None] + 1j)))
        return f, fs

    def spectral_theta(self, omega, eta, coeffs, H, cache=None):
        z = np.asarray(omega) + 1j * eta
        C = self._abcd(z) if cache is None else cache
        f, fs = self._hardy(z, H)
        a = coeffs[:H] + 1j * coeffs[H:2 * H]
        b = coeffs[2 * H:3 * H] + 1j * coeffs[3 * H:]
        t = f @ a + fs @ b
        th = (C[:, 0] * t + C[:, 1]) / (C[:, 2] * t + C[:, 3])
        nev = 1j * (1 + th) / (1 - th)
        return (nev).imag / np.pi        # A = -Im G / pi = Im(nev)/pi since G = -nev

    def hardy_optimize(self, omega, eta=0.02, H=20, lam=1e-3, maxiter=5000, mu_norm=100.0):
        """Optimise theta_{M+1} in the Hardy basis for a smooth, normalised, non-negative spectrum:
        minimise mu_norm (1 - int A)^2 + lam int (A'')^2 + 10 mu_norm (int max(-A,0))^2, analytic gradient."""
        from scipy.optimize import minimize
        w = np.asarray(omega, float)
        dw = np.gradient(w)
        z = w + 1j * eta
        C = self._abcd(z)
        f, fs = self._hardy(z, H)
        B = np.hstack([f, 1j * f, fs, 1j * fs])          # t = B @ c (c real, length 4H)
        h = w[1] - w[0]
        D2 = (np.eye(len(w), k=1)[1:-1] - 2 * np.eye(len(w))[1:-1] + np.eye(len(w), k=-1)[1:-1] * 0)
        D2 = np.zeros((len(w) - 2, len(w)))
        for i in range(len(w) - 2):
            D2[i, i], D2[i, i + 1], D2[i, i + 2] = 1, -2, 1
        D2 /= h**2

        def fg(c):
            t = B @ c
            den = C[:, 2] * t + C[:, 3]
            th = (C[:, 0] * t + C[:, 1]) / den
            nev = 1j * (1 + th) / (1 - th)
            A = nev.imag / np.pi
            dth = (C[:, 0] * C[:, 3] - C[:, 1] * C[:, 2]) / den**2
            dnev = 2j / (1 - th) ** 2
            J = (dnev * dth)[:, None] * B                  # d nev / d c
            dA = J.imag / np.pi                              # (nw, 4H)
            norm = np.sum(A * dw)
            d2 = D2 @ A
            neg = np.sum(np.clip(-A, 0, None) * dw)
            val = mu_norm * (1 - norm) ** 2 + lam * np.sum(d2**2) * h + 10 * mu_norm * neg**2
            g = (-2 * mu_norm * (1 - norm) * (dw @ dA) + lam * 2 * h * (d2 @ (D2 @ dA))
                 + 20 * mu_norm * neg * (-(dw * (A < 0)) @ dA))
            return val, g

        # continuation in the basis size (Fei et al.): H = 1, 2, ... with warm starts
        c0 = np.zeros(4 * H)
        cost0 = float(fg(c0)[0])
        best = None
        for Hk in range(1, H + 1):
            sel = np.concatenate([np.arange(Hk) + j * H for j in range(4)])

            def fk(ck):
                c = np.zeros(4 * H)
                c[sel] = ck
                v, g = fg(c)
                return v, g[sel]
            start = (best[sel] if best is not None else np.zeros(4 * Hk))
            r = minimize(fk, start, jac=True, method="L-BFGS-B", options=dict(maxiter=maxiter))
            c = np.zeros(4 * H)
            c[sel] = r.x
            if best is None or r.fun <= fg(best)[0]:
                best = c
        A = self.spectral_theta(w, eta, best, H, cache=C)
        return dict(A=A, coeffs=best, cost0=cost0, cost=float(fg(best)[0]), success=True)


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


def hardy_auto(nv, omega, eta=0.02, lams=(1e-4, 3e-4, 1e-3), Hs=(10, 25), norm_tol=0.02, neg_tol=1e-2):
    """Run hardy_optimize on a small (lambda, H) grid (the landscape has local minima) and return the
    smoothest spectrum (minimal int A''^2) among those that are normalised and non-negative."""
    w = np.asarray(omega, float)
    dw = np.gradient(w)
    h = w[1] - w[0]
    best, cands = None, []
    for lam in lams:
        for H in Hs:
            r = nv.hardy_optimize(w, eta, H=H, lam=lam)
            A = r["A"]
            rough = float(np.sum(((A[2:] - 2 * A[1:-1] + A[:-2]) / h**2) ** 2) * h)
            ok = abs(np.sum(A * dw) - 1) < norm_tol and A.min() > -neg_tol
            cands.append(dict(lam=lam, H=H, rough=rough, ok=ok, A=A))
            if ok and (best is None or rough < best["rough"]):
                best = cands[-1]
    if best is None:
        best = min(cands, key=lambda c: abs(np.sum(c["A"] * dw) - 1))
    return dict(A=best["A"], lam=best["lam"], H=best["H"], candidates=cands)
