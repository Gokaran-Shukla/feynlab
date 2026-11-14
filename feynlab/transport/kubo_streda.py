"""Disorder-averaged Kubo-Streda Hall response with intrinsic, side-jump and
skew-scattering diagrams, for any band model and any 'operator current'.

Disorder: V(r) = sum_i V_i delta(r - R_i), identical in every orbital,
<V_i^2> = V0^2 (Gaussian part), <V_i^3> = V1^3 (third cumulant -> skew).

Diagrams (non-crossing approximation, Sinitsyn et al. PRB 75, 045315 (2007), Fig. 1):
* Born self-energy        Sigma^R = n_i V0^2 int G^R           (anti-Hermitian part)
* Fermi-surface bubble    sigma^I  = (1/2) int Tr[2 J_x G^R Gam_y G^A - J_x G^R v_y G^R - J_x G^A v_y G^A]
                          (units e^2/h for charge, e/2pi for other operators)
* ladder vertex           Gam_y = v_y + n_i V0^2 int G^R Gam_y G^A          -> side jump +
                                                                         'intrinsic skew'
* third-order vertex      + n_i V1^3 [g^R M + M g^A],  M = int G^R Gam G^A, g = int G
                                                                         -> conventional skew
The Fermi-sea term sigma^II is disorder-independent at leading order and equals
(intrinsic - bare-bubble) in the clean limit, so we report

    intrinsic       : Berry-curvature (Kubo) result, exact clean limit
    gaussian_extr.  : sigma^I(ladder) - sigma^I(bare)   (= side jump + adist + intrinsic skew)
    skew_V1         : sigma^I(ladder + V1^3) - sigma^I(ladder)
    total           : sum

Crossed (X and Psi) diagrams of the same order (Ado et al., EPL 111, 37004 (2015))
are NOT included; see docs/software_tutorials for this limitation.
"""
from __future__ import annotations

import numpy as np

from .intrinsic import intrinsic_hall, k_grid_polar, k_grid_lattice


class DisorderedHall:
    def __init__(self, model, mu: float, n_i: float, V0: float, V1: float = 0.0,
                 grid=None, scba: bool = False, scba_iter: int = 50):
        self.model = model
        self.mu = float(mu)
        self.n_i, self.V0, self.V1 = float(n_i), float(V0), float(V1)
        self.grid = grid
        kx, ky, w = grid
        self.w = w
        self.H = model.H(kx, ky)
        self.vx, self.vy = model.dH(kx, ky)
        nb = self.H.shape[-1]
        self.nb = nb
        self.I = np.eye(nb)
        self.sigma_R = self._born(scba, scba_iter)
        z = self.mu * self.I - self.H
        self.GR = np.linalg.inv(z - self.sigma_R)
        self.GA = np.conj(np.swapaxes(self.GR, -1, -2))
        self.gR_loc = self.sigma_R / (self.n_i * self.V0**2)

    # ------------------------------------------------------------------
    def _integrate(self, A):
        nb = A.shape[-1]
        return np.tensordot(self.w.ravel(), A.reshape(-1, nb, nb), axes=(0, 0))

    def _born(self, scba, n_iter):
        nb = self.nb
        z = self.mu * self.I - self.H
        G0 = np.linalg.inv(z + 1e-12j * self.I) if False else None
        # first Born with an infinitesimal broadening is ill-defined numerically on a grid;
        # iterate SCBA on the anti-Hermitian part starting from a small guess.
        gam = 1e-3
        S = -1j * gam * self.I
        for _ in range(n_iter if scba else 30):
            G = np.linalg.inv(z - S)
            g = self._integrate(G)
            S_new = self.n_i * self.V0**2 * g
            S_new = 0.5 * (S_new - np.conj(S_new.T))  # anti-Hermitian part only
            if np.max(np.abs(S_new - S)) < 1e-12:
                S = S_new
                break
            S = S_new if scba else 0.5 * (S + S_new)
        return S

    # ------------------------------------------------------------------
    def _vertex(self, vy, ladder=True, skew=False):
        """k-independent dressing of the RA vertex: Gam(k) = vy(k) + X."""
        nb = self.nb
        GR, GA = self.GR, self.GA
        M0 = self._integrate(GR @ vy @ GA)
        if not ladder and not skew:
            return np.zeros((nb, nb), complex)
        # superoperator L(X) = int G^R X G^A  as nb^2 x nb^2 matrix
        GRf = GR.reshape(-1, nb, nb) * self.w.reshape(-1, 1, 1)
        GAf = GA.reshape(-1, nb, nb)
        # (G^R X G^A)_{ij} = sum_ab GR_ia X_ab GA_bj  ->  L[(i,j),(a,b)]
        L = np.einsum("kia,kbj->ijab", GRf, GAf).reshape(nb * nb, nb * nb)
        a = self.n_i * self.V0**2 if ladder else 0.0
        b = self.n_i * self.V1**3 if skew else 0.0
        gR, gA = self.gR_loc, np.conj(self.gR_loc.T)
        # X = a M + b (gR M + M gA),  M = M0 + L X
        Id = np.eye(nb * nb)
        lgR = np.kron(gR, np.eye(nb))        # vec(gR M) = (gR (x) I) vec(M) for row-major vec
        rgA = np.kron(np.eye(nb), gA.T)      # vec(M gA)
        K = a * Id + b * (lgR + rgA)
        # X = K (M0 + L X)  ->  (I - K L) X = K M0   (row-major vec, L acts on vec(X))
        X = np.linalg.solve(Id - K @ L, K @ M0.reshape(-1))
        return X.reshape(nb, nb)

    def sigma_I(self, op=None, ladder=True, skew=False):
        """Fermi-surface Streda term (units e^2/h or e/2pi)."""
        GR, GA = self.GR, self.GA
        vx, vy = self.vx, self.vy
        Jx = vx if op is None else 0.5 * (op @ vx + vx @ op)
        X = self._vertex(vy, ladder, skew)
        Gam = vy + X
        RA = Jx @ GR @ Gam @ GA
        RR = Jx @ GR @ vy @ GR
        AA = Jx @ GA @ vy @ GA
        tr = np.trace(2 * RA - RR - AA, axis1=-2, axis2=-1)
        return float(np.real(np.sum(self.w * tr)) / 2.0)

    def sigma_xx(self, ladder=True):
        """Longitudinal conductivity (e^2/h), RA part with ladder vertex."""
        GR, GA = self.GR, self.GA
        vx = self.vx
        X = self._vertex(vx, ladder, False)
        RA = vx @ GR @ (vx + X) @ GA
        RR = vx @ GR @ vx @ GR
        tr = np.trace(2 * RA - 2 * np.real(RR), axis1=-2, axis2=-1)
        return float(np.real(np.sum(self.w * tr)) / 2.0)

    # ------------------------------------------------------------------
    def sigma_I_clean(self, op=None, eta: float | None = None):
        """Bare Fermi-surface bubble with an infinitesimal *scalar* broadening
        (clean limit); used to obtain the Fermi-sea term sigma^II = intrinsic -
        sigma^I_clean."""
        if eta is None:
            eta = float(np.max(np.abs(self.sigma_R)))
        z = (self.mu + 1j * eta) * self.I - self.H
        GR = np.linalg.inv(z)
        GA = np.conj(np.swapaxes(GR, -1, -2))
        vx, vy = self.vx, self.vy
        Jx = vx if op is None else 0.5 * (op @ vx + vx @ op)
        tr = np.trace(2 * (Jx @ GR @ vy @ GA) - Jx @ GR @ vy @ GR - Jx @ GA @ vy @ GA,
                      axis1=-2, axis2=-1)
        return float(np.real(np.sum(self.w * tr)) / 2.0)

    def decompose(self, op=None, intrinsic_grid=None, intrinsic_value=None, tail: float = 0.0):
        """Returns intrinsic, extrinsic (Gaussian: side jump + anomalous distribution +
        Gaussian skew, non-crossing approximation), skew_V1 and total.

        total = sigma^II + sigma^I(ladder [+ V1^3]),  sigma^II = intrinsic - sigma^I(clean)."""
        if intrinsic_value is None:
            intrinsic_value = intrinsic_hall(self.model, self.mu, op=op,
                                             grid=intrinsic_grid or self.grid) + tail
        clean = self.sigma_I_clean(op)
        sigma_II = intrinsic_value - clean
        lad = self.sigma_I(op, ladder=True)
        out = dict(intrinsic=intrinsic_value, sigma_II=sigma_II, sigma_I_clean=clean,
                   sigma_I_bare_born=self.sigma_I(op, ladder=False), sigma_I_ladder=lad)
        out["extrinsic_gaussian"] = lad - clean
        out["skew_V1"] = (self.sigma_I(op, ladder=True, skew=True) - lad) if self.V1 != 0.0 else 0.0
        out["total"] = sigma_II + lad + out["skew_V1"]
        out["Gamma_broadening"] = float(np.max(np.abs(self.sigma_R)))
        return out
