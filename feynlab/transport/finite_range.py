"""Kubo-Streda Hall response with finite-range (correlated) disorder: k-dependent Born self-energy and ladder vertex.

Disorder V(r) = sum_i V_i u(r - R_i) M with <V^2> = V0^2, so the impurity line carries
    W(q) = n_i V0^2 |u(q)|^2,     Gaussian impurities  u(q) = exp(-q^2 d^2 / 4)  (range d)
(the same form factor as transport.boltzmann.SemiclassicalHall(range_d=d)).  With d = 0 everything reduces to the
delta-correlated DisorderedHall (k-independent self-energy and vertex).  Non-crossing approximation:

    Sigma^R(k)  = int_k' W(k - k') M G^R(k') M^+                     (anti-Hermitian part, self-consistent Born)
    Gam_y(k)    = v_y(k) + int_k' W(k - k') M G^R(k') Gam_y(k') G^A(k') M^+      (ladder, Bethe-Salpeter)
    sigma^I     = (1/2) int_k Tr[2 J_x G^R Gam_y G^A - J_x G^R v_y G^R - J_x G^A v_y G^A]      (e^2/h)
    sigma^II    = intrinsic - sigma^I(clean)  (disorder independent at leading order, as in DisorderedHall)

The convolutions int_k' W(k - k') f(k') are evaluated on a uniform Cartesian k grid by zero-padded FFTs (linear,
not cyclic, convolution: exact for the truncated domain), so one ladder application costs O(N^2 log N); the vertex
equation is solved with GMRES.  Any 2D band model with H(kx, ky) and dH(kx, ky) is accepted.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator, gmres

from .intrinsic import intrinsic_hall


class FiniteRangeHall:
    def __init__(self, model, mu: float, n_i: float, V0: float, d: float = 0.0, K: float = 2.0, N: int = 512,
                 center=(0.0, 0.0), form_factor=None, imp_matrix=None, scba_iter: int = 60, tol: float = 1e-11):
        """K, N: square grid [-K, K]^2 (+ center) with N^2 points.  form_factor(q) overrides the Gaussian u(q)."""
        self.model, self.mu = model, float(mu)
        self.n_i, self.V0, self.d = float(n_i), float(V0), float(d)
        h = 2 * K / N
        x = -K + h * (np.arange(N) + 0.5)
        self.kx, self.ky = np.meshgrid(x + center[0], x + center[1], indexing="ij")
        self.h, self.N = h, N
        self.wk = h * h / (2 * np.pi) ** 2
        self.H = model.H(self.kx, self.ky)
        self.vx, self.vy = model.dH(self.kx, self.ky)
        nb = self.H.shape[-1]
        self.nb = nb
        self.I = np.eye(nb)
        self.M = np.eye(nb) if imp_matrix is None else np.asarray(imp_matrix, complex)
        u = form_factor or (lambda q: np.exp(-q**2 * self.d**2 / 4.0))
        # kernel on the (2N-1)^2 grid of differences, zero padded to 2N x 2N
        dq = h * np.arange(-(N - 1), N)
        QX, QY = np.meshgrid(dq, dq, indexing="ij")
        ker = np.zeros((2 * N, 2 * N))
        ker[: 2 * N - 1, : 2 * N - 1] = self.n_i * self.V0**2 * np.abs(u(np.hypot(QX, QY))) ** 2
        self.kerF = np.fft.rfft2(ker)
        self.tol = tol
        self.sigma_R = self._scba(scba_iter)
        self.GR = np.linalg.inv(self.mu * self.I - self.H - self.sigma_R)
        self.GA = np.conj(np.swapaxes(self.GR, -1, -2))

    # ------------------------------------------------------------------ convolution
    def conv(self, F):
        """C(k) = int_k' W(k - k') M F(k') M^+  for F of shape (N, N, nb, nb)."""
        N, nb = self.N, self.nb
        F = np.einsum("ab,xybc,dc->xyad", self.M, F, np.conj(self.M))
        out = np.empty((N, N, nb, nb), complex)
        for a in range(nb):
            for b in range(nb):
                pad = np.zeros((2 * N, 2 * N), complex)
                pad[:N, :N] = F[:, :, a, b]
                re = np.fft.irfft2(np.fft.rfft2(pad.real) * self.kerF, s=(2 * N, 2 * N))
                im = np.fft.irfft2(np.fft.rfft2(pad.imag) * self.kerF, s=(2 * N, 2 * N))
                out[:, :, a, b] = (re + 1j * im)[N - 1: 2 * N - 1, N - 1: 2 * N - 1]
        return out * self.wk

    def _scba(self, n_iter):
        S = -1e-3j * np.broadcast_to(self.I, self.H.shape).copy()
        for it in range(n_iter):
            G = np.linalg.inv(self.mu * self.I - self.H - S)
            Sn = self.conv(G)
            Sn = 0.5 * (Sn - np.conj(np.swapaxes(Sn, -1, -2)))
            diff = np.max(np.abs(Sn - S))
            S = Sn if it > 3 else 0.5 * (S + Sn)
            if diff < 1e-12:
                break
        self.scba_residual = float(diff)
        return S

    # ------------------------------------------------------------------ ladder
    def vertex(self, v, ladder=True):
        """Gam(k) = v(k) + X(k),  X = conv(G^R (v + X) G^A)  solved with GMRES."""
        if not ladder:
            return v
        GR, GA = self.GR, self.GA
        shp = v.shape
        rhs = self.conv(GR @ v @ GA).ravel()

        def mv(x):
            X = x.reshape(shp)
            return x - self.conv(GR @ X @ GA).ravel()
        A = LinearOperator((rhs.size, rhs.size), matvec=mv, dtype=complex)
        X, info = gmres(A, rhs, rtol=self.tol, atol=0.0, restart=60, maxiter=400)
        self.gmres_info = info
        return v + X.reshape(shp)

    def sigma_I(self, op=None, ladder=True):
        GR, GA, vx, vy = self.GR, self.GA, self.vx, self.vy
        Jx = vx if op is None else 0.5 * (op @ vx + vx @ op)
        Gam = self.vertex(vy, ladder)
        tr = np.trace(2 * Jx @ GR @ Gam @ GA - Jx @ GR @ vy @ GR - Jx @ GA @ vy @ GA, axis1=-2, axis2=-1)
        return float(np.real(tr.sum()) * self.wk / 2.0)

    def sigma_xx(self, ladder=True):
        GR, GA, vx = self.GR, self.GA, self.vx
        Gam = self.vertex(vx, ladder)
        tr = np.trace(2 * vx @ GR @ Gam @ GA - 2 * np.real(vx @ GR @ vx @ GR), axis1=-2, axis2=-1)
        return float(np.real(tr.sum()) * self.wk / 2.0)

    def sigma_I_clean(self, op=None, eta=None):
        if eta is None:
            eta = float(np.max(np.abs(self.sigma_R)))
        GR = np.linalg.inv((self.mu + 1j * eta) * self.I - self.H)
        GA = np.conj(np.swapaxes(GR, -1, -2))
        vx, vy = self.vx, self.vy
        Jx = vx if op is None else 0.5 * (op @ vx + vx @ op)
        tr = np.trace(2 * Jx @ GR @ vy @ GA - Jx @ GR @ vy @ GR - Jx @ GA @ vy @ GA, axis1=-2, axis2=-1)
        return float(np.real(tr.sum()) * self.wk / 2.0)

    def decompose(self, op=None, intrinsic_value=None, intrinsic_grid=None):
        """intrinsic, extrinsic_gaussian (= side jump + anomalous distribution + Gaussian skew, NCA), total."""
        if intrinsic_value is None:
            intrinsic_value = intrinsic_hall(self.model, self.mu, op=op,
                                             grid=intrinsic_grid or (self.kx, self.ky, np.full(self.kx.shape, self.wk)))
        clean = self.sigma_I_clean(op)
        lad = self.sigma_I(op, ladder=True)
        return dict(intrinsic=intrinsic_value, sigma_I_clean=clean, sigma_I_ladder=lad,
                    sigma_I_bare=self.sigma_I(op, ladder=False), extrinsic_gaussian=lad - clean,
                    total=intrinsic_value - clean + lad, Gamma_broadening=float(np.max(np.abs(self.sigma_R))))
