"""Crossed impurity diagrams (X and Psi) for the Fermi-surface Hall response.

Ado, Dmitriev, Ostrovsky & Titov, EPL 111, 37004 (2015) and PRL 117, 046601 (2016)
showed that diagrams with two *crossing* impurity lines contribute at the same
order as side jump (skew scattering on pairs of close impurities).  In real space
(r = vector between the two impurities):

    sigma^X   = (n V0^2)^2 int d^2r Tr[ J^x_r G^R_{-r} J^y_r G^A_{-r} ]
    sigma^Psi = (n V0^2)^2 int d^2r Tr[ J^x_r G^R_{-r} G^R_r J^y_{-r} + h.c. ]

with G_r = int d^2p/(2pi)^2 G_p e^{ipr} and the dressed vertex
J^a_r = int d^2p/(2pi)^2 G^A_p Gam^a G^R_p e^{ipr}, Gam^a the ladder-dressed
current (A-R ordering).  Units: e^2/h for charge (e/2pi for an operator current
on the x vertex).  The (n V0^2)^2 is cancelled by the 1/(n V0^2) of each J, so
the result is O(disorder^0) and is extrapolated to vanishing broadening.

Everything is done with FFTs on an N x N momentum grid of a box of size L, so
any continuum band model (H(kx,ky), dH) can be used.
"""
from __future__ import annotations

import numpy as np


class CrossedDiagrams:
    def __init__(self, model, mu: float, nV2: float, L: float = 400.0, N: int = 1024,
                 born_iter: int = 40, lattice: bool = False):
        """lattice=True: model is a square-lattice tight-binding model (a = 1); the momentum grid
        is then the full Brillouin zone (L = N sites, periodic), r runs over lattice vectors and
        no Nyquist masking is needed (H(k) is periodic)."""
        if lattice:
            L = float(N)
        self.model, self.mu, self.nV2 = model, float(mu), float(nV2)
        self.L, self.N = float(L), int(N)
        p = 2 * np.pi * np.fft.fftfreq(N, d=L / N)
        self.px, self.py = np.meshgrid(p, p, indexing="ij")
        # drop the Nyquist row/column so that the grid is inversion symmetric
        # (otherwise odd functions such as v_x leak into the singular density channel)
        m1 = np.ones(N)
        if not lattice:
            m1[N // 2] = 0.0
        self.mask = np.outer(m1, m1)
        H = model.H(self.px, self.py)
        self.vx, self.vy = model.dH(self.px, self.py)
        nb = H.shape[-1]
        self.nb = nb
        I = np.eye(nb)
        w = 1.0 / L**2  # int d^2p/(2pi)^2 -> (1/L^2) sum_p
        z = self.mu * I - H
        S = -1e-3j * I
        for _ in range(born_iter):  # self-consistent Born, anti-Hermitian part only
            G = np.linalg.inv(z - S) * self.mask[..., None, None]
            g = w * G.reshape(-1, nb, nb).sum(0)
            Sn = nV2 * g
            Sn = 0.5 * (Sn - Sn.conj().T)
            if np.max(np.abs(Sn - S)) < 1e-12:
                S = Sn
                break
            S = 0.5 * (S + Sn)
        self.Sigma = S
        self.GR = np.linalg.inv(z - S) * self.mask[..., None, None]
        self.GA = np.conj(np.swapaxes(self.GR, -1, -2))
        self.w = w

    # --------------------------------------------------------------
    def _ar_vertex(self, v):
        """Gam = v + nV2 int G^A Gam G^R  (k-independent correction X), returns Gam(p)."""
        nb = self.nb
        GA = self.GA.reshape(-1, nb, nb)
        GR = self.GR.reshape(-1, nb, nb)
        vv = v.reshape(-1, nb, nb)
        M0 = self.w * np.einsum("kia,kab,kbj->ij", GA, vv, GR)
        Lop = self.w * np.einsum("kia,kbj->ijab", GA, GR).reshape(nb * nb, nb * nb)
        X = np.linalg.solve(np.eye(nb * nb) - self.nV2 * Lop, self.nV2 * M0.reshape(-1)).reshape(nb, nb)
        return v + X

    def _to_r(self, Ap):
        """int d^2p/(2pi)^2 A_p e^{ipr} on the conjugate real-space grid."""
        return np.fft.ifft2(Ap, axes=(0, 1)) * (self.N**2) / self.L**2

    @staticmethod
    def _minus_r(Ar):
        N = Ar.shape[0]
        idx = (-np.arange(N)) % N
        return Ar[idx][:, idx]

    def compute(self, op=None):
        vx, vy = self.vx, self.vy
        Jx_bare = vx if op is None else 0.5 * (op @ vx + vx @ op)
        Gx = self._ar_vertex(Jx_bare) if op is None else Jx_bare + (self._ar_vertex(vx) - vx) * 0.0
        if op is not None:
            # dress the operator current with the same A-R ladder
            Gx = self._ar_vertex(Jx_bare)
        Gy = self._ar_vertex(vy)
        Jx = self._to_r(self.GA @ Gx @ self.GR)
        Jy = self._to_r(self.GA @ Gy @ self.GR)
        GRr = self._to_r(self.GR)
        GAr = self._to_r(self.GA)
        GRm, GAm = self._minus_r(GRr), self._minus_r(GAr)
        Jym = self._minus_r(Jy)
        dA = (self.L / self.N) ** 2
        X = np.trace(Jx @ GRm @ Jy @ GAm, axis1=-2, axis2=-1).sum() * dA
        P = Jx @ GRm @ GRr @ Jym
        Psi = (np.trace(P, axis1=-2, axis2=-1).sum() + np.conj(np.trace(P, axis1=-2, axis2=-1).sum())) * dA
        f = self.nV2**2
        return dict(X=float(np.real(f * X)), Psi=float(np.real(f * Psi)),
                    Gamma_broadening=float(np.max(np.abs(self.Sigma))))


def crossed_hall(model, mu, op=None, runs=((0.015, 900.0, 1536), (0.0075, 1600.0, 2048)),
                 lattice: bool = False):
    """X and Psi contributions extrapolated linearly to vanishing broadening.
    runs = sequence of (n V0^2, box L, grid N); keep L >> v/Gamma and pi N/L >~ 4 k_F.
    lattice=True: tight-binding model on the full BZ (L is ignored, L = N)."""
    xs, ps, gs = [], [], []
    for nv, L, N in runs:
        r = CrossedDiagrams(model, mu, nv, L=L, N=N, lattice=lattice).compute(op)
        xs.append(r["X"]); ps.append(r["Psi"]); gs.append(r["Gamma_broadening"])
    if len(gs) > 1:
        X0 = float(np.polyval(np.polyfit(gs, xs, 1), 0.0))
        P0 = float(np.polyval(np.polyfit(gs, ps, 1), 0.0))
    else:
        X0, P0 = xs[0], ps[0]
    return dict(X=X0, Psi=P0, raw=dict(Gamma=gs, X=xs, Psi=ps))
