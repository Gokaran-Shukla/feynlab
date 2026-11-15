"""Cluster CT-INT: weak-coupling continuous-time QMC for a Hubbard cluster in a bath (DCA / cluster DMFT).

Rubtsov, Savkin, Lichtenstein PRB 72, 035122 (2005); Gull et al. RMP 83, 349 (2011) Sec. III.A;
for DCA: Gull, Staar, Fuchs, Nukala, Summers, Pruschke, Schulthess, Maier, PRB 83, 075122 (2011).

The interaction U sum_i n_i,up n_i,dn on N_c cluster sites is expanded with auxiliary-field shifted vertices
(tau, site, s) exactly as in the single-site solver; the determinant matrices are
    D_sigma[p, q] = G0~_sigma(r_p - r_q, tau_p - tau_q) - alpha_sigma(s_p) delta_pq,
with the cluster-translation-invariant Weiss field G0~(r, tau) = (1/N_c) sum_K e^{iK.r} G0~(K, tau)
(G0~ includes the Hartree shift mu -> mu - U/2).  The Green function is measured in K space,
    G(K, iw) = G0~(K) - G0~(K)^2 S(K, iw) / beta,
    S(K, iw) = (1/N_c) < s sum_pq e^{iw tau_p - iK.r_p} M_pq e^{-iw tau_q + iK.r_q} >,   M = D^{-1}.
Numerically exact (no truncation); the sign problem appears away from half filling for N_c >= 4.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid
from .ctint import _g0_interp, _ins_update, _rem_update


@njit(cache=True)
def _run(beta, U, delta, gt, n0, diffidx, Nc, w, phase, n_warm, n_meas, meas_every, seed, kmax):
    """gt: (Nr, ntab) tables of G0~(r, tau) (spin symmetric); phase[K, site] = e^{-iK.r_site}."""
    np.random.seed(seed)
    ntab = gt.shape[1]
    nw = w.shape[0]
    nK = phase.shape[0]
    taus = np.zeros(kmax)
    sites = np.zeros(kmax, dtype=np.int64)
    ss = np.zeros(kmax, dtype=np.int64)
    Mu = np.zeros((kmax, kmax))
    Md = np.zeros((kmax, kmax))
    k = 0
    S = np.zeros((2, nK, nw), dtype=np.complex128)
    hist = np.zeros(kmax)
    sign = 1.0
    sign_acc = 0.0
    nmeas = 0
    Qu = np.zeros(kmax); Ru = np.zeros(kmax); Qd = np.zeros(kmax); Rd = np.zeros(kmax)
    for step in range(n_warm + n_meas):
        if np.random.random() < 0.5:
            if k < kmax - 1:
                tn = np.random.random() * beta
                xn = np.random.randint(Nc)
                sn = 1 if np.random.random() < 0.5 else -1
                au = 0.5 + sn * delta
                ad = 0.5 - sn * delta
                for i in range(k):
                    a = diffidx[sites[i], xn]
                    b = diffidx[xn, sites[i]]
                    Qu[i] = _g0_interp(taus[i] - tn, beta, gt[a], ntab)
                    Ru[i] = _g0_interp(tn - taus[i], beta, gt[b], ntab)
                    Qd[i] = Qu[i]
                    Rd[i] = Ru[i]
                lu = n0 - au
                ld = n0 - ad
                MQu = np.ascontiguousarray(Mu[:k, :k]) @ Qu[:k]
                RMu = Ru[:k] @ np.ascontiguousarray(Mu[:k, :k])
                MQd = np.ascontiguousarray(Md[:k, :k]) @ Qd[:k]
                RMd = Rd[:k] @ np.ascontiguousarray(Md[:k, :k])
                lu -= Ru[:k] @ MQu
                ld -= Rd[:k] @ MQd
                A = -U * beta * Nc / (k + 1) * lu * ld
                if np.random.random() < abs(A):
                    if A < 0:
                        sign = -sign
                    _ins_update(Mu, MQu, RMu, lu, k)
                    _ins_update(Md, MQd, RMd, ld, k)
                    taus[k] = tn
                    sites[k] = xn
                    ss[k] = sn
                    k += 1
        else:
            if k > 0:
                p = np.random.randint(k)
                A = k / (-U * beta * Nc) * Mu[p, p] * Md[p, p]
                if np.random.random() < abs(A):
                    if A < 0:
                        sign = -sign
                    last = k - 1
                    _rem_update(Mu, p, k)
                    _rem_update(Md, p, k)
                    if p != last:
                        taus[p] = taus[last]
                        sites[p] = sites[last]
                        ss[p] = ss[last]
                    k -= 1
        if step >= n_warm and (step - n_warm) % meas_every == 0:
            nmeas += 1
            sign_acc += sign
            hist[k] += 1.0
            if k > 0:
                # X[i, (K, w)] = e^{i w tau_i} e^{-iK.r_i};  S = sum_ij X_i M_ij conj(X_j)  (vectorised)
                X = np.empty((k, nK * nw), dtype=np.complex128)
                for i in range(k):
                    for iwi in range(nw):
                        e = np.exp(1j * w[iwi] * taus[i])
                        for KK in range(nK):
                            X[i, KK * nw + iwi] = e * phase[KK, sites[i]]
                Xc = np.conj(X)
                Mu_c = np.ascontiguousarray(Mu[:k, :k]).astype(np.complex128)
                Md_c = np.ascontiguousarray(Md[:k, :k]).astype(np.complex128)
                Bu = Mu_c @ Xc
                Bd = Md_c @ Xc
                for KK in range(nK):
                    for iwi in range(nw):
                        col = KK * nw + iwi
                        su = 0.0 + 0.0j
                        sd = 0.0 + 0.0j
                        for i in range(k):
                            su += X[i, col] * Bu[i, col]
                            sd += X[i, col] * Bd[i, col]
                        S[0, KK, iwi] += sign * su / Nc
                        S[1, KK, iwi] += sign * sd / Nc
    return S, hist, sign_acc, nmeas


class ClusterCTINT:
    """sites: (Nc, d) integer coordinates; L: periodicity per dimension; K: (Nc, d) cluster momenta.
    solve(G0_K) with G0_K[K, iw] the Weiss field of the physical problem (spin symmetric)."""

    def __init__(self, grid: MatsubaraGrid, U: float, sites, L, K, delta: float = 0.51, n_iw_meas: int = 64,
                 n_tau_table: int = 8001):
        self.grid, self.U, self.delta = grid, float(U), float(delta)
        self.sites = np.atleast_2d(np.asarray(sites, int))
        self.L = np.atleast_1d(np.asarray(L, int))
        self.K = np.atleast_2d(np.asarray(K, float))
        self.Nc = len(self.sites)
        self.n_iw_meas = min(n_iw_meas, grid.n_iw)
        self.n_tau_table = n_tau_table
        Nc = self.Nc
        lookup = {tuple(r % self.L): i for i, r in enumerate(self.sites)}
        self.diffidx = np.array([[lookup[tuple((self.sites[i] - self.sites[j]) % self.L)] for j in range(Nc)]
                                 for i in range(Nc)], dtype=np.int64)
        self.phase = np.exp(-1j * self.K @ self.sites.T)          # [K, site] = e^{-iK.r}

    def _tables(self, G0t_K):
        g = self.grid
        tau = np.linspace(0, g.beta, self.n_tau_table)
        # G0~(r, iw) = (1/Nc) sum_K e^{iK.r} G0~(K, iw)
        Gr = np.einsum("Kr,Kw->rw", np.conj(self.phase), G0t_K) / self.Nc
        tabs = []
        for i, r in enumerate(self.sites):
            c1 = 1.0 if not np.any(r % self.L) else 0.0
            cf = MatsubaraGrid.tail_coefficients(g.iw, Gr[i], n_fit=min(60, g.n_iw // 4))
            tabs.append(g.iw_to_tau(Gr[i], c1, float(np.real(cf[1])), float(np.real(cf[2])), tau=tau))
        return np.array(tabs)

    def solve(self, G0_K, n_warm=20000, n_meas=400000, meas_every=20, seed=99, n_chains=4, kmax=600, parallel="serial"):
        g = self.grid
        U = self.U
        G0t = 1.0 / (1.0 / np.asarray(G0_K) - U / 2.0)
        gt = self._tables(G0t)
        n0 = -gt[0, -1]
        w = g.iw[: self.n_iw_meas].imag.copy()
        Gs, hist, signs = [], 0, []
        m = self.n_iw_meas
        from ..utils.parallel import map_chains
        args = [(g.beta, U, self.delta, gt, n0, self.diffidx, self.Nc, w, self.phase, n_warm, n_meas, meas_every,
                 seed + 7919 * c, kmax) for c in range(n_chains)]
        for S, h, sa, nm in map_chains(_run, args, backend=parallel):
            Sm = 0.5 * (S[0] + S[1]) / sa
            Gs.append(G0t[:, :m] - G0t[:, :m] ** 2 * Sm / g.beta)
            hist = hist + h
            signs.append(sa / nm)
        Gs = np.array(Gs)
        G = Gs.mean(0)
        G0 = np.asarray(G0_K)[:, :m]
        sig = 1.0 / G0 - 1.0 / G
        sc = 1.0 / G0[None] - 1.0 / Gs
        err = (sc.real.std(0, ddof=1) + 1j * sc.imag.std(0, ddof=1)) / np.sqrt(n_chains) if n_chains > 1 else None
        hn = hist / hist.sum()
        return dict(G_K=G, Sigma_K=sig, err_Sigma=err, sign=float(np.mean(signs)),
                    k_mean=float((np.arange(len(hn)) * hn).sum()), iw=g.iw[:m])
