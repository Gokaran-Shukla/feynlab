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
The sign problem cannot be removed in general (it is NP-hard, Troyer & Wiese PRL 94, 170201 (2005)); this module
measures and manages it: pooled sign-reweighted estimators with jackknife errors over bins (`solve`), a pilot-run
optimisation of the auxiliary shift delta (`optimize_delta`), and an exponential cost model (`sign_cost_model`).
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid
from .ctint import _g0_interp, _ins_update, _rem_update


@njit(cache=True)
def _run(beta, U, delta, gt, n0, diffidx, Nc, w, phase, n_warm, n_meas, meas_every, seed, kmax, n_bins):
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
    n_tot = n_meas // meas_every + 1
    Sb = np.zeros((n_bins, 2, nK, nw), dtype=np.complex128)
    sb = np.zeros(n_bins)
    cb = np.zeros(n_bins)
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
            ib = min((nmeas - 1) * n_bins // n_tot, n_bins - 1)
            sb[ib] += sign
            cb[ib] += 1.0
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
                        Sb[ib, 0, KK, iwi] += sign * su / Nc
                        Sb[ib, 1, KK, iwi] += sign * sd / Nc
    return S, hist, sign_acc, nmeas, Sb, sb, cb


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

    def solve(self, G0_K, n_warm=20000, n_meas=400000, meas_every=20, seed=99, n_chains=4, kmax=600, parallel="serial",
              n_bins=10):
        """Sign-aware estimate: G = G0~ - G0~^2 <s S>/<s> / beta with <s S> and <s> pooled over all chains, and
        jackknife errors over n_bins x n_chains bins of the ratio (the sign fluctuations and their correlation with
        S are propagated; per-chain ratios would be biased when <s> is small)."""
        g = self.grid
        U = self.U
        G0t = 1.0 / (1.0 / np.asarray(G0_K) - U / 2.0)
        gt = self._tables(G0t)
        n0 = -gt[0, -1]
        w = g.iw[: self.n_iw_meas].imag.copy()
        hist = 0
        m = self.n_iw_meas
        from ..utils.parallel import map_chains
        args = [(g.beta, U, self.delta, gt, n0, self.diffidx, self.Nc, w, self.phase, n_warm, n_meas, meas_every,
                 seed + 7919 * c, kmax, n_bins) for c in range(n_chains)]
        Sbins, sbins, cbins = [], [], []
        for S, h, sa, nm, Sb, sb, cb in map_chains(_run, args, backend=parallel):
            hist = hist + h
            Sbins.append(Sb); sbins.append(sb); cbins.append(cb)
        Sb = np.concatenate(Sbins); sb = np.concatenate(sbins); cb = np.concatenate(cbins)
        keep = cb > 0
        Sb, sb, cb = Sb[keep], sb[keep], cb[keep]
        G0 = np.asarray(G0_K)[:, :m]

        def estimate(Ssum, ssum, csum):
            Sm = 0.5 * (Ssum[0] + Ssum[1]) / ssum
            Gk = G0t[:, :m] - G0t[:, :m] ** 2 * Sm / g.beta
            return Gk, 1.0 / G0 - 1.0 / Gk, ssum / csum
        G, sig, sgn = estimate(Sb.sum(0), sb.sum(), cb.sum())
        nb = len(sb)
        jG, jS, js = [], [], []
        for b in range(nb):
            a, bb, c = estimate(Sb.sum(0) - Sb[b], sb.sum() - sb[b], cb.sum() - cb[b])
            jG.append(a); jS.append(bb); js.append(c)
        jS, js = np.array(jS), np.array(js)
        jk = lambda x: np.sqrt((nb - 1) * np.mean(np.abs(x - x.mean(0)) ** 2, axis=0))
        err = jk(jS.real) + 1j * jk(jS.imag)
        hn = hist / hist.sum()
        k_mean = float((np.arange(len(hn)) * hn).sum())
        return dict(G_K=G, Sigma_K=sig, err_Sigma=err, sign=float(sgn), sign_err=float(jk(js)),
                    k_mean=k_mean, iw=g.iw[:m], n_bins_total=nb)


def optimize_delta(solver_factory, G0_K, deltas=(0.51, 0.6, 0.75, 1.0), pilot=dict(n_warm=5000, n_meas=100000),
                   **kw):
    """Pilot runs over the auxiliary-field shift delta (alpha = 1/2 +- delta). The statistical error of a
    sign-reweighted estimate scales as 1 / (<s> sqrt(N_eff)); with the CPU time t per pilot run the figure of merit
    is efficiency = <s>^2 / t (error^2 x time is what the user pays). Returns (best_delta, table)."""
    import time
    table = []
    for d in deltas:
        s = solver_factory(d)
        t0 = time.time()
        r = s.solve(G0_K, **pilot, **kw)
        dt = time.time() - t0
        table.append(dict(delta=d, sign=r["sign"], sign_err=r["sign_err"], k_mean=r["k_mean"], time=dt,
                          efficiency=r["sign"] ** 2 / dt))
    best = max(table, key=lambda x: x["efficiency"])
    return best["delta"], table


def sign_cost_model(betas, signs, Nc):
    """Fit <s>(beta) = exp(-beta Nc df) (Troyer & Wiese PRL 94, 170201 (2005): the sign is the ratio of two
    partition functions, exponentially small in beta x volume). Returns (df, predict) with predict(beta) giving
    the expected sign and the relative CPU cost ~ 1/<s>^2 for fixed error."""
    b = np.asarray(betas, float)
    y = -np.log(np.clip(np.asarray(signs, float), 1e-300, None))
    df = float(np.sum(b * y) / np.sum(b * b) / Nc)

    def predict(beta):
        s = np.exp(-beta * Nc * df)
        return dict(sign=s, relative_cost=1.0 / s**2)
    return df, predict
