"""Dynamical cluster approximation (DCA) for the 2D Hubbard model.

Hettler et al., PRB 58, R7475 (1998); Maier, Jarrell, Pruschke, Hettler, RMP 77, 1027 (2005).
The Brillouin zone is tiled into Nc = Lc x Lc patches around cluster momenta K;
the self-energy is constant within a patch, Sigma(k) = Sigma(K).

    G(K, iw)  = (Nc/N) sum_{k in patch K} [iw + mu - eps_k - Sigma(K, iw)]^{-1}
    G0(K, iw) = [G(K)^{-1} + Sigma(K)]^{-1}            (cluster Weiss field)

Cluster solvers:
  'second_order' : second-order perturbation theory on the cluster (cluster IPT),
                   Sigma_ij(tau) = U^2 G0~_ij(tau)^2 G0~_ji(beta - tau); Nc = 1 reduces exactly to DMFT+IPT.
  'ctint'        : numerically exact cluster CT-INT (qmc.ctint_cluster), validated against ED of a
                   cluster+bath; high frequencies use the Hartree tail Sigma -> U n/2 + c/(iw).
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import MatsubaraGrid


class DCA2D:
    def __init__(self, Lc=2, nk=64, t=0.25, tp=0.0, U=2.0, beta=20.0, n_iw=512, mu=None):
        if nk % Lc:
            raise ValueError("nk must be a multiple of Lc")
        self.Lc, self.Nc, self.nk = Lc, Lc * Lc, nk
        self.U = float(U)
        self.grid = MatsubaraGrid(beta, n_iw)
        self.mu = U / 2.0 if mu is None else mu
        k = 2 * np.pi * np.arange(nk) / nk
        kx, ky = np.meshgrid(k, k, indexing="ij")
        self.eps = -2 * t * (np.cos(kx) + np.cos(ky)) - 4 * tp * np.cos(kx) * np.cos(ky)
        # patch index: nearest cluster momentum K = 2 pi (m, n)/Lc
        step = nk // Lc
        ix = ((np.arange(nk) + step // 2) // step) % Lc
        self.patch = (ix[:, None] * Lc + ix[None, :])
        Km = 2 * np.pi * np.arange(Lc) / Lc
        self.K = np.array([(a, b) for a in Km for b in Km])

    def coarse_grain(self, sigma_K):
        iw = self.grid.iw
        G = np.zeros((self.Nc, len(iw)), complex)
        for p in range(self.Nc):
            e = self.eps[self.patch == p]
            G[p] = np.mean(1.0 / (iw[None, :] + self.mu - e[:, None] - sigma_K[p][None, :]), axis=0)
        return G

    def _K_to_r(self, AK):
        Lc = self.Lc
        return np.fft.ifft2(AK.reshape(Lc, Lc, -1), axes=(0, 1))   # A(r) = (1/Nc) sum_K e^{iKr} A(K)

    def _r_to_K(self, Ar):
        Lc = self.Lc
        return np.fft.fft2(Ar, axes=(0, 1)).reshape(self.Nc, -1)

    def cluster_sigma(self, G0K):
        g = self.grid
        Lc = self.Lc
        G0t = 1.0 / (1.0 / G0K - self.U / 2.0)             # Hartree shift (half filling)
        if abs(self.mu - self.U / 2) < 1e-12:
            G0t = 1j * G0t.imag if self.Nc == 1 else G0t
        Gr_iw = self._K_to_r(G0t)                            # (Lc, Lc, n_iw)
        # G(r, tau): r-space functions are not real in general -> transform real and imaginary parts
        Gr_tau = np.empty((Lc, Lc, g.n_tau))
        for a in range(Lc):
            for b in range(Lc):
                c1 = 1.0 if (a == 0 and b == 0) else 0.0
                Gr_tau[a, b] = self._iw_to_tau_general(Gr_iw[a, b], c1)
        idx = (-np.arange(Lc)) % Lc
        Gm = Gr_tau[idx][:, idx][..., ::-1]                  # G(-r, beta - tau)
        S_tau = self.U**2 * Gr_tau**2 * Gm
        S_iw = np.empty((Lc, Lc, g.n_iw), complex)
        for a in range(Lc):
            for b in range(Lc):
                S_iw[a, b] = g.tau_to_iw(S_tau[a, b])
        return self.U / 2.0 + self._r_to_K(S_iw)

    def _iw_to_tau_general(self, f_iw, c1):
        """G(r,tau) for cluster-translation-invariant, inversion-symmetric clusters: G(r,iw) obeys
        G(r,-iw) = conj(G(r,iw)) so G(r,tau) is real."""
        return self.grid.iw_to_tau(f_iw, c1, 0.0, 0.0)

    def cluster_sigma_ctint(self, G0K, **kw):
        from ..qmc.ctint_cluster import ClusterCTINT
        Lc = self.Lc
        sites = [(x, y) for x in range(Lc) for y in range(Lc)]
        kw = dict(kw)
        m = kw.pop("n_iw_meas", 48)
        solver = ClusterCTINT(self.grid, self.U, sites=sites, L=[Lc, Lc], K=self.K, n_iw_meas=m)
        r = solver.solve(G0K, **kw)
        self.ctint_last = r
        iw = self.grid.iw
        out = np.empty((self.Nc, len(iw)), complex)
        n = 1.0 if abs(self.mu - self.U / 2) < 1e-12 else None
        for p in range(self.Nc):
            sl = r["Sigma_K"][p]
            er = np.abs(r["err_Sigma"][p]) if r["err_Sigma"] is not None else np.zeros(m)
            sinf = self.U / 2.0 if n is not None else float(np.real(sl[-1]))
            bad = np.where(er > 0.05 + 0.05 * np.abs(sl - sinf))[0]
            mc = max(int(bad[0]) if len(bad) else m, 2)
            out[p, :mc] = sl[:mc]
            out[p, mc:] = sinf + (sl[mc - 1] - sinf) * iw[mc - 1] / iw[mc:]
        return out

    def solve(self, n_iter=100, mix=0.5, tol=1e-6, solver="second_order", ctint_kwargs=None):
        iw = self.grid.iw
        sig = np.full((self.Nc, len(iw)), self.U / 2.0, complex)
        hist = []
        for it in range(n_iter):
            G = self.coarse_grain(sig)
            G0 = 1.0 / (1.0 / G + sig)
            new = self.cluster_sigma(G0) if solver == "second_order" else self.cluster_sigma_ctint(G0, **(ctint_kwargs or {}))
            diff = float(np.max(np.abs(new - sig)))
            sig = mix * new + (1 - mix) * sig
            hist.append(diff)
            if diff < tol:
                break
        self.sigma_K = sig
        G = self.coarse_grain(sig)
        Z = [float(1 / (1 - sig[p, 0].imag / iw[0].imag)) for p in range(self.Nc)]
        return dict(iw=iw, Sigma_K=sig, G_K=G, K=self.K, Z=Z, iterations=len(hist), history=hist)
