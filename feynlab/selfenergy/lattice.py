"""Momentum-dependent self-energies of the 2D Hubbard model (weak -> moderate coupling).

Diagrams evaluated
------------------
* 'hf'   : Hartree (tadpole) only.
* '2nd'  : second-order 'pair-bubble' diagram (Mattuck Sec. 9.8):
              Sigma(r,tau) = U^2 chi0(r,tau) G(r,tau),
              chi0(r,tau)  = G(r,tau) G(-r,beta-tau)    (single bubble)
* 'rpa'  : spin+charge fluctuation exchange (Berk-Schrieffer / one-shot FLEX):
              W(q,iv) = U^2 [ 3/2 chi0/(1-U chi0) + 1/2 chi0/(1+U chi0) - chi0 ]
           = sum of bubble chains (RPA, Mattuck Ch. 10.4) and ladders in the
           particle-hole channel, minus the double-counted 2nd order piece.
           Expanding to O(U^2) recovers '2nd' exactly (tested).

Both are evaluated with the (Hartree-shifted) non-interacting G0 ('one-shot') or
self-consistently ('self_consistent=True', i.e. GF2 / FLEX).
All convolutions use FFTs in r and tau; tails of G are treated analytically.
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import MatsubaraGrid
from ..core.lattice import square_lattice_eps


class HubbardLattice2D:
    def __init__(self, nk: int = 32, t: float = 1.0, tp: float = 0.0, U: float = 2.0,
                 beta: float = 5.0, n_iw: int = 256, mu: float | None = None):
        self.nk = nk
        self.U = float(U)
        self.eps = square_lattice_eps(nk, t, tp)
        self.grid = MatsubaraGrid(beta, n_iw)
        self.mu = U / 2.0 if mu is None else float(mu)  # half filling for tp = 0

    # ------------------------------------------------------------------
    def _g_tau_r(self, g_k_iw, shift):
        """G(r,tau) from G(k,iw); shift = (eps_k - mu_eff) for the 1/iw^2 tail."""
        gk_tau = self.grid.iw_to_tau(g_k_iw, 1.0, shift[..., None] if np.ndim(shift) else shift)
        return np.fft.ifft2(gk_tau, axes=(0, 1)).real  # r-space (real for inversion-symm.)

    @staticmethod
    def _reverse(g_r_tau):
        """G(-r, beta - tau)."""
        nk = g_r_tau.shape[0]
        idx = (-np.arange(nk)) % nk
        return g_r_tau[idx][:, idx][..., ::-1]

    def chi0(self, g_r_tau):
        """Bubble chi0(q, iv) >= 0 (per spin)."""
        chi_r_tau = g_r_tau * self._reverse(g_r_tau)
        chi_q_tau = np.fft.fft2(chi_r_tau, axes=(0, 1)).real
        return self.grid.tau_to_iv(chi_q_tau), chi_r_tau

    # ------------------------------------------------------------------
    def sigma(self, kind: str = "2nd", self_consistent: bool = False, n_iter: int = 30,
              mix: float = 0.5, tol: float = 1e-6, n_target: float | None = None):
        grid = self.grid
        iw = grid.iw
        U = self.U
        nk = self.nk
        n = 0.5 if n_target is None else n_target
        hartree = U * n
        shift = self.eps - self.mu + hartree
        g = 1.0 / (iw[None, None, :] - shift[..., None])
        sig = np.zeros_like(g)
        info = {}
        for it in range(n_iter if self_consistent else 1):
            g_r_tau = self._g_tau_r(g, shift)
            chi_q_iv, chi_r_tau = self.chi0(g_r_tau)
            if kind == "hf":
                s_r_tau = np.zeros_like(g_r_tau)
            elif kind == "2nd":
                s_r_tau = U**2 * chi_r_tau * g_r_tau
            elif kind == "rpa":
                c = chi_q_iv
                if np.max(U * c[..., 0].real) >= 1.0:
                    raise RuntimeError("RPA spin susceptibility diverged (U chi0 >= 1): magnetic instability")
                w_rest = U**2 * (1.5 * c / (1 - U * c) + 0.5 * c / (1 + U * c) - 2.0 * c)
                w_rest_r_tau = np.fft.ifft2(grid.iv_to_tau(w_rest), axes=(0, 1)).real
                s_r_tau = (U**2 * chi_r_tau + w_rest_r_tau) * g_r_tau
                info["stoner_max"] = float(np.max(U * c[..., 0].real))
            else:
                raise ValueError(kind)
            s_k_tau = np.fft.fft2(s_r_tau, axes=(0, 1)).real
            sig_new = grid.tau_to_iw(s_k_tau)
            diff = float(np.max(np.abs(sig_new - sig)))
            sig = mix * sig_new + (1 - mix) * sig if self_consistent else sig_new
            g = 1.0 / (iw[None, None, :] - shift[..., None] - sig)
            if self_consistent and diff < tol:
                break
        self.sigma_k = sig
        self.g_k = g
        info.update(dict(hartree=hartree, chi0_q_iv=chi_q_iv, iterations=it + 1))
        return sig, info

    # ------------------------------------------------------------------
    def local_sigma(self):
        return self.sigma_k.mean(axis=(0, 1))

    def k_point(self, ix, iy):
        return self.sigma_k[ix, iy]
