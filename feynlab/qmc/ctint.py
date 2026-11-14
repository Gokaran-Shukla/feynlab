"""Continuous-time interaction-expansion QMC (CT-INT) for a single-orbital Anderson
impurity.  This is a *diagrammatic* Monte Carlo: the partition function is expanded
in powers of U,

    Z/Z0 = sum_k (-U)^k / k!  int dtau_1..dtau_k  prod_sigma det D_sigma^(k),

and the sum over all Feynman diagrams of order k with fixed vertex times is
performed exactly by the determinant (Wick's theorem).  Configurations (k,
{tau_i}, {s_i}) are sampled with a Metropolis random walk.

References: A. N. Rubtsov et al., PRB 72, 035122 (2005);
            E. Gull et al., Rev. Mod. Phys. 83, 349 (2011), Sec. III.A.

Interaction:  U n_up n_dn = (U/2) sum_{s=+-1} (n_up - a_up(s))(n_dn - a_dn(s)) + ...,
              a_up(s) = 1/2 + s*delta, a_dn(s) = 1/2 - s*delta.
The -U/2 shift of the chemical potential this generates is applied internally.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid


@njit(cache=True)
def _g0_interp(t, beta, gt_tab, ntab):
    # antiperiodic extension; table covers [0, beta]
    sgn = 1.0
    if t < 0.0:
        t += beta
        sgn = -1.0
    x = t / beta * (ntab - 1)
    i = int(x)
    if i >= ntab - 1:
        i = ntab - 2
    f = x - i
    return sgn * ((1.0 - f) * gt_tab[i] + f * gt_tab[i + 1])


@njit(cache=True)
def _ins_update(M, MQ, RM, lam, k):
    for i in range(k):
        for j in range(k):
            M[i, j] += MQ[i] * RM[j] / lam
    for i in range(k):
        M[i, k] = -MQ[i] / lam
        M[k, i] = -RM[i] / lam
    M[k, k] = 1.0 / lam


@njit(cache=True)
def _rem_update(M, p, k):
    last = k - 1
    if p != last:
        for j in range(k):
            tmp = M[p, j]
            M[p, j] = M[last, j]
            M[last, j] = tmp
        for i in range(k):
            tmp = M[i, p]
            M[i, p] = M[i, last]
            M[i, last] = tmp
    piv = M[last, last]
    for i in range(last):
        for j in range(last):
            M[i, j] -= M[i, last] * M[last, j] / piv
    for i in range(k):
        M[i, last] = 0.0
        M[last, i] = 0.0


@njit(cache=True)
def _run(beta, U, delta, gt_up, gt_dn, n0_up, n0_dn, w, n_warm, n_meas, meas_every,
         seed, kmax):
    np.random.seed(seed)
    ntab = gt_up.shape[0]
    nw = w.shape[0]
    taus = np.zeros(kmax)
    ss = np.zeros(kmax, dtype=np.int64)
    Mu = np.zeros((kmax, kmax))
    Md = np.zeros((kmax, kmax))
    k = 0
    S_up = np.zeros(nw, dtype=np.complex128)
    S_dn = np.zeros(nw, dtype=np.complex128)
    hist = np.zeros(kmax)
    sign_acc = 0.0
    sign = 1.0
    n_acc_ins = 0
    n_acc_rem = 0
    nmeas_done = 0
    total = n_warm + n_meas
    Qu = np.zeros(kmax)
    Ru = np.zeros(kmax)
    Qd = np.zeros(kmax)
    Rd = np.zeros(kmax)
    for step in range(total):
        if np.random.random() < 0.5:
            # ---------------- insertion ----------------
            if k < kmax - 1:
                tn = np.random.random() * beta
                sn = 1 if np.random.random() < 0.5 else -1
                au = 0.5 + sn * delta
                ad = 0.5 - sn * delta
                for i in range(k):
                    Qu[i] = _g0_interp(taus[i] - tn, beta, gt_up, ntab)
                    Ru[i] = _g0_interp(tn - taus[i], beta, gt_up, ntab)
                    Qd[i] = _g0_interp(taus[i] - tn, beta, gt_dn, ntab)
                    Rd[i] = _g0_interp(tn - taus[i], beta, gt_dn, ntab)
                lu = n0_up - au
                ld = n0_dn - ad
                MQu = np.ascontiguousarray(Mu[:k, :k]) @ Qu[:k]
                RMu = Ru[:k] @ np.ascontiguousarray(Mu[:k, :k])
                MQd = np.ascontiguousarray(Md[:k, :k]) @ Qd[:k]
                RMd = Rd[:k] @ np.ascontiguousarray(Md[:k, :k])
                lu -= Ru[:k] @ MQu
                ld -= Rd[:k] @ MQd
                A = -U * beta / (k + 1) * lu * ld
                if np.random.random() < abs(A):
                    n_acc_ins += 1
                    if A < 0:
                        sign = -sign
                    # update inverses
                    _ins_update(Mu, MQu, RMu, lu, k)
                    _ins_update(Md, MQd, RMd, ld, k)
                    taus[k] = tn
                    ss[k] = sn
                    k += 1
        else:
            # ---------------- removal ----------------
            if k > 0:
                p = np.random.randint(k)
                A = k / (-U * beta) * Mu[p, p] * Md[p, p]
                if np.random.random() < abs(A):
                    n_acc_rem += 1
                    if A < 0:
                        sign = -sign
                    last = k - 1
                    # move p to last position (swap rows/cols) then shrink
                    _rem_update(Mu, p, k)
                    _rem_update(Md, p, k)
                    if p != last:
                        taus[p] = taus[last]
                        ss[p] = ss[last]
                    k -= 1
        # ---------------- measurement ----------------
        if step >= n_warm and (step - n_warm) % meas_every == 0:
            nmeas_done += 1
            sign_acc += sign
            hist[k] += 1.0
            if k > 0:
                for iw_i in range(nw):
                    wn = w[iw_i]
                    su = 0.0 + 0.0j
                    sd = 0.0 + 0.0j
                    for i in range(k):
                        ei = np.exp(1j * wn * taus[i])
                        for j in range(k):
                            ej = np.exp(-1j * wn * taus[j])
                            su += ei * Mu[i, j] * ej
                            sd += ei * Md[i, j] * ej
                    S_up[iw_i] += sign * su
                    S_dn[iw_i] += sign * sd
    return S_up, S_dn, hist, sign_acc, nmeas_done, n_acc_ins, n_acc_rem


class CTINTSolver:
    """Single-orbital CT-INT impurity solver.

    Parameters
    ----------
    grid  : MatsubaraGrid
    U     : Hubbard interaction
    delta : auxiliary-field shift (0.5 avoids the sign problem at half filling for U>0;
            values 0.5-0.6 are standard)
    """

    def __init__(self, grid: MatsubaraGrid, U: float, delta: float = 0.51,
                 n_iw_meas: int = 128, n_tau_table: int = 8001):
        self.grid = grid
        self.U = float(U)
        self.delta = float(delta)
        self.n_iw_meas = min(n_iw_meas, grid.n_iw)
        self.n_tau_table = n_tau_table

    def _g0_tau_table(self, g0t_iw):
        g = self.grid
        tau = np.linspace(0, g.beta, self.n_tau_table)
        c = MatsubaraGrid.tail_coefficients(g.iw, g0t_iw, n_fit=min(60, g.n_iw // 4))
        gt = g.iw_to_tau(g0t_iw, 1.0, c[1], c[2], tau=tau)
        return gt

    def solve(self, g0_up_iw, g0_dn_iw=None, n_warm=20000, n_meas=400000,
              meas_every=20, seed=1234, n_chains: int = 1, kmax: int = 400):
        """g0_*_iw: Weiss field of the *physical* problem 1/(iw + mu - Delta(iw)).
        Returns dict with G(iw), Sigma(iw), average sign, <k>, order histogram."""
        g = self.grid
        if g0_dn_iw is None:
            g0_dn_iw = g0_up_iw
        U = self.U
        # shifted bare propagators (mu -> mu - U/2)
        g0t_up = 1.0 / (1.0 / g0_up_iw - U / 2.0)
        g0t_dn = 1.0 / (1.0 / g0_dn_iw - U / 2.0)
        gt_up = self._g0_tau_table(g0t_up)
        gt_dn = self._g0_tau_table(g0t_dn)
        n0_up = -gt_up[-1]   # G0(0^-) = -G0(beta^-) = n0
        n0_dn = -gt_dn[-1]
        w = g.iw[: self.n_iw_meas].imag.copy()
        beta = g.beta
        n_iw = self.n_iw_meas
        hist = np.zeros(kmax)
        sgn_tot = 0.0
        nm = 0
        acc = 0
        G_chains_up, G_chains_dn = [], []
        for c in range(n_chains):
            a, b, h, s, n, ai, ar = _run(beta, U, self.delta, gt_up, gt_dn, n0_up, n0_dn,
                                          w, n_warm, n_meas, meas_every, seed + 7919 * c, kmax)
            hist += h
            sgn_tot += s
            nm += n
            acc += ai + ar
            G_chains_up.append(g0t_up[:n_iw] - g0t_up[:n_iw] ** 2 * (a / s) / beta)
            G_chains_dn.append(g0t_dn[:n_iw] - g0t_dn[:n_iw] ** 2 * (b / s) / beta)
        Gu = np.array(G_chains_up)
        Gd = np.array(G_chains_dn)
        G_up = Gu.mean(0)
        G_dn = Gd.mean(0)
        sigma_up = 1.0 / g0_up_iw[:n_iw] - 1.0 / G_up
        sigma_dn = 1.0 / g0_dn_iw[:n_iw] - 1.0 / G_dn
        if n_chains > 1:
            sig_chains = 1.0 / g0_up_iw[:n_iw] - 1.0 / Gu
            err_sigma = (sig_chains.real.std(0, ddof=1) + 1j * sig_chains.imag.std(0, ddof=1)) / np.sqrt(n_chains)
            err_G = (Gu.real.std(0, ddof=1) + 1j * Gu.imag.std(0, ddof=1)) / np.sqrt(n_chains)
        else:
            err_sigma = np.full(n_iw, np.nan + 1j * np.nan)
            err_G = err_sigma.copy()
        avg_sign = sgn_tot / nm
        orders = np.arange(kmax)
        hist_n = hist / hist.sum()
        return dict(G_up=G_up, G_dn=G_dn, Sigma_up=sigma_up, Sigma_dn=sigma_dn,
                    iw=g.iw[:n_iw], sign=avg_sign, k_mean=float((orders * hist_n).sum()),
                    order_histogram=hist_n, n_meas=nm, err_Sigma=err_sigma, err_G=err_G,
                    acceptance=acc / max(1, n_chains * (n_warm + n_meas)))
