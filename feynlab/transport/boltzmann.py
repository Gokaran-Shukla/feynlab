"""Semiclassical (Boltzmann) theory of the anomalous / spin Hall effect with all
mechanisms resolved separately (Sinitsyn, J. Phys.: Condens. Matter 20, 023201
(2008); Sinitsyn et al. PRB 75, 045315 (2007)):

    intrinsic        -q^2 sum_n int f Omega_n                     (Berry curvature)
    side jump        q sum_l g0_l v^sj_l,   v^sj_l = sum_l' w_l'l dr_l'l
    anomalous distr. q sum_l g^adist_l v_l
    skew (V1^3)      third-order asymmetric rate  w^(3a)

The Gaussian skew term is NOT computed here but by combining this module with the
Kubo-Streda NCA total (see transport.mechanisms, validated to 1e-4 against
Sinitsyn Eqs. 65, 67, 72, B8).  Direct attempt at the n_i-independent Gaussian ('intrinsic') skew term from the fourth-
order rate w^(4a) (Sinitsyn Eq. 71).  An attempt that keeps interband self-energy
pairings and |T2|^2 interference reproduced exactly 1/2 of Eq. 71, so it is not
shipped; for the massive Dirac model the analytic value is available as
MassiveDirac.ref_intrinsic_skew().  Crossed X/Psi diagrams (Ado et al. EPL 111,
37004 (2015)) are likewise not included.

Coordinate shift (gauge invariant):
    dr_l'l = A_l' - A_l - (d_k + d_k') arg<u_l'|u_l>,   A = i<u|d_k u>.

Disorder: V(r) = sum_i V_i delta(r - R_i) (orbital-independent),
<V^2> = V0^2, <V^3> = V1^3.   Single Fermi contour in band 'band', T = 0.
Conductivities in units of e^2/h (hbar = 1).
"""
from __future__ import annotations

import numpy as np

from .intrinsic import intrinsic_hall, k_grid_polar


class FermiContour:
    def __init__(self, model, mu: float, band: int, n_theta: int = 128, kmax: float = 50.0):
        self.model, self.mu, self.band = model, float(mu), band
        th = 2 * np.pi * (np.arange(n_theta) + 0.5) / n_theta
        self.theta = th
        # bisection for k_F(theta) on the chosen band
        lo = np.zeros(n_theta)
        hi = np.full(n_theta, kmax)
        def e(k):
            return np.linalg.eigvalsh(model.H(k * np.cos(th), k * np.sin(th)))[:, band]
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            above = e(mid) > mu
            hi = np.where(above, mid, hi)
            lo = np.where(above, lo, mid)
        k = 0.5 * (lo + hi)
        self.k = k
        self.kx, self.ky = k * np.cos(th), k * np.sin(th)
        H = model.H(self.kx, self.ky)
        E, U = np.linalg.eigh(H)
        self.E, self.U = E, U
        self.u = U[:, :, band]
        hx, hy = model.dH(self.kx, self.ky)
        self.hx, self.hy = hx, hy
        u = self.u
        self.v = np.stack([np.real(np.einsum("ti,tij,tj->t", u.conj(), hx, u)),
                           np.real(np.einsum("ti,tij,tj->t", u.conj(), hy, u))], -1)
        rhat = np.stack([np.cos(th), np.sin(th)], -1)
        vr = np.sum(self.v * rhat, -1)
        # contour measure: int d^2k/(2pi)^2 delta(e - mu) F = sum_t rho_t F_t
        self.rho = k / np.abs(vr) * (2 * np.pi / n_theta) / (2 * np.pi) ** 2
        # |d u> in parallel-transport gauge: sum_{m != n} |m><m|dH|n>/(e_n - e_m)
        nb = E.shape[-1]
        du = []
        for h in (hx, hy):
            hm = np.einsum("tim,tij,tj->tm", U.conj(), h, u)  # <m|dH|n>
            de = E[:, band][:, None] - E
            fac = np.where(np.arange(nb)[None, :] == band, 0.0, 1.0 / np.where(de == 0, 1, de))
            du.append(np.einsum("tim,tm->ti", U, hm * fac))
        self.du = np.stack(du, 0)  # (2, t, nb)

    # ------------------------------------------------------------------
    def overlaps(self):
        return np.einsum("ai,bi->ab", self.u.conj(), self.u)  # <u_a|u_b>

    def coordinate_shift(self):
        """dr[l', l, :] for scattering l -> l'."""
        S = self.overlaps()  # S[a,b] = <u_a|u_b>
        Sc = np.where(np.abs(S) < 1e-12, 1e-12, S)
        dr = np.zeros(S.shape + (2,))
        for c in range(2):
            # <u_l'|d u_l>  -> M[l', l]
            M = np.einsum("ai,bi->ab", self.u.conj(), self.du[c])
            # <u_l|d u_l'>  -> N[l, l'] ; need conj(<u_l|du_l'>) arranged as [l', l]
            N = np.einsum("bi,ai->ab", self.u.conj(), self.du[c])  # [l', l] = <u_l|du_l'>
            term_k = -np.imag(M / Sc)
            term_kp = -np.imag(np.conj(N) / Sc)
            dr[..., c] = term_k + term_kp
        return dr


class SemiclassicalHall:
    def __init__(self, model, mu: float, band: int, n_i: float = 1.0, V0: float = 0.1,
                 V1: float = 0.0, n_theta: int = 128, kmax: float = 50.0):
        self.c = FermiContour(model, mu, band, n_theta, kmax)
        self.model, self.mu = model, mu
        self.n_i, self.V0, self.V1 = n_i, V0, V1
        S = self.c.overlaps()
        self.S = S
        rho = self.c.rho
        # rate l -> l' (per unit contour weight rho_l'):  W2[l', l]
        self.W2 = 2 * np.pi * n_i * V0**2 * np.abs(S) ** 2

    # ------------------------------------------------------------------
    def _collision(self, WS, WA=None):
        """Matrix C with (C phi)_l = sum_l' rho_l' [ WS[l',l](phi_l - phi_l') - WA[l,l'] phi_l' ]."""
        rho = self.c.rho
        C = np.diag((WS * rho[:, None]).sum(0)) - WS.T * rho[None, :]
        if WA is not None:
            C = C - WA * rho[None, :]
        return C

    def _solve(self, C, rhs):
        # remove zero mode (particle conservation) by adding projector onto constant
        rho = self.c.rho
        A = C + np.outer(np.ones_like(rho), rho) * C.diagonal().mean()
        return np.linalg.solve(A, rhs)

    def w3_skew(self):
        """Antisymmetric third-order rate W3[l', l] (l -> l')."""
        S = self.S
        rho = self.c.rho
        # Im[ S[l,l'] S[l',l''] S[l'',l] ] summed over l'' with rho_l''
        T = np.einsum("ab,bc,ca,c->ba", S, S, S, rho)  # index [l', l] with a=l, b=l'
        return (2 * np.pi) ** 2 * self.n_i * self.V1**3 * np.imag(T)

    # ------------------------------------------------------------------
    def compute(self, intrinsic_grid=None):
        c = self.c
        rho = c.rho
        WS = self.W2
        C2 = self._collision(WS)
        # field along y, q = 1 ;  C phi = q E . v
        phi0 = self._solve(C2, c.v[:, 1])
        phix = self._solve(C2, c.v[:, 0])
        sxx = np.sum(rho * phix * c.v[:, 0])
        dr = c.coordinate_shift()  # [l', l, :]
        vsj = np.einsum("pl,plk,p->lk", WS, dr, rho)  # v^sj_l = sum_l' w_l'l dr_l'l
        sig_sj = np.sum(rho * phi0 * vsj[:, 0])
        phi_ad = self._solve(C2, -vsj[:, 1])
        sig_ad = np.sum(rho * phi_ad * c.v[:, 0])
        out = dict(sigma_xx=2 * np.pi * sxx, side_jump=2 * np.pi * sig_sj,
                   anomalous_distribution=2 * np.pi * sig_ad)
        WA = np.zeros_like(WS)
        if self.V1 != 0.0:
            W3 = self.w3_skew()
            # WA[l, l'] = rate l' -> l  =  W3[l, l']
            # first order in W3:  C2 phi_sk = sum_l' rho_l' W3[l, l'] phi0_l'
            phi_sk = self._solve(C2, (W3 * rho[None, :]) @ phi0)
            out["skew_V1"] = 2 * np.pi * np.sum(rho * phi_sk * c.v[:, 0])
        else:
            out["skew_V1"] = 0.0
        out["intrinsic"] = intrinsic_hall(self.model, self.mu, grid=intrinsic_grid) if intrinsic_grid else np.nan
        return out
