"""Continuous-time hybridization-expansion QMC (CT-HYB), segment picture.

Werner, Comanac, de' Medici, Troyer & Millis, PRL 97, 076405 (2006);
Gull et al., RMP 83, 349 (2011), Sec. III.B.

Expansion of Z in the hybridization: for each spin a configuration is a set of
'segments' [t_s, t_e) during which the impurity is occupied.  Weight

    w = prod_sigma det F_sigma  *  exp( mu (L_up + L_dn) - U O_up,dn )

with F_ij = F(t_e,i - t_s,j),  F(tau) = -Delta(-tau) = Delta(beta - tau) for
0 < tau < beta  (antiperiodic), L the total segment length, O the up/down
overlap.  Updates: insert/remove segments and anti-segments (gaps), plus the
'full line' state.  G(iw) is measured directly from the inverse hybridization
matrix M = F^{-1}:  G(iw) = (1/beta) <sum_ij M_ji e^{iw(t_e,i - t_s,j)}>  (with this F sign).
Exact at all U (no truncation); efficient at strong coupling where CT-INT is not.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid


@njit(cache=True)
def _F(t, beta, tab, n):
    s = 1.0
    if t < 0.0:
        t += beta
        s = -1.0
    x = t / beta * (n - 1)
    i = int(x)
    if i >= n - 1:
        i = n - 2
    f = x - i
    return s * ((1 - f) * tab[i] + f * tab[i + 1])


@njit(cache=True)
def _seg_overlap(a0, a1, b0, b1, beta):
    """overlap of two (possibly wrapping) intervals [a0,a1) and [b0,b1) on the circle"""
    def pieces(x0, x1):
        if x1 >= x0:
            return ((x0, x1), (0.0, 0.0))
        return ((x0, beta), (0.0, x1))
    pa = pieces(a0, a1)
    pb = pieces(b0, b1)
    tot = 0.0
    for i in range(2):
        for j in range(2):
            lo = max(pa[i][0], pb[j][0])
            hi = min(pa[i][1], pb[j][1])
            if hi > lo:
                tot += hi - lo
    return tot


@njit(cache=True)
def _overlap_with(ts, te, k, full, a0, a1, beta):
    if full:
        return a1 - a0 if a1 >= a0 else beta - a0 + a1
    tot = 0.0
    for i in range(k):
        tot += _seg_overlap(a0, a1, ts[i], te[i], beta)
    return tot


@njit(cache=True)
def _contains(ts, te, k, t):
    """index of segment containing t, else -1"""
    for i in range(k):
        if te[i] >= ts[i]:
            if ts[i] <= t < te[i]:
                return i
        else:
            if t >= ts[i] or t < te[i]:
                return i
    return -1


@njit(cache=True)
def _next_start_dist(ts, k, t, beta):
    best = beta
    for i in range(k):
        d = ts[i] - t
        if d <= 0:
            d += beta
        if d < best:
            best = d
    return best


@njit(cache=True)
def _next_end_dist(te, k, t, beta):
    best = beta
    for i in range(k):
        d = te[i] - t
        if d <= 0:
            d += beta
        if d < best:
            best = d
    return best


@njit(cache=True)
def _ins(M, k, Q, R, lam):
    MQ = M[:k, :k] @ Q[:k]
    RM = R[:k] @ M[:k, :k]
    for i in range(k):
        for j in range(k):
            M[i, j] += MQ[i] * RM[j] / lam
    for i in range(k):
        M[i, k] = -MQ[i] / lam
        M[k, i] = -RM[i] / lam
    M[k, k] = 1.0 / lam


@njit(cache=True)
def _rem(M, p, k):
    last = k - 1
    if p != last:
        for j in range(k):
            tmp = M[p, j]; M[p, j] = M[last, j]; M[last, j] = tmp
        for i in range(k):
            tmp = M[i, p]; M[i, p] = M[i, last]; M[i, last] = tmp
    piv = M[last, last]
    for i in range(last):
        for j in range(last):
            M[i, j] -= M[i, last] * M[last, j] / piv


@njit(cache=True)
def _run(beta, U, mu_up, mu_dn, Ftab_up, Ftab_dn, w, n_warm, n_meas, meas_every, seed, kmax):
    np.random.seed(seed)
    ntab = Ftab_up.shape[0]
    nw = w.shape[0]
    TS = np.zeros((2, kmax)); TE = np.zeros((2, kmax))
    Ms = np.zeros((2, kmax, kmax))
    K = np.zeros(2, dtype=np.int64)
    FULL = np.zeros(2, dtype=np.bool_)
    Gacc = np.zeros((2, nw), dtype=np.complex128)
    nacc = np.zeros(2); dacc = 0.0
    sgn_acc = 0.0
    hist = np.zeros(kmax)
    nmeas = 0
    Q = np.zeros(kmax); R = np.zeros(kmax)
    mus = np.array([mu_up, mu_dn])
    for step in range(n_warm + n_meas):
        s = 0 if np.random.random() < 0.5 else 1
        o = 1 - s
        tab = Ftab_up if s == 0 else Ftab_dn
        k = K[s]
        ts = TS[s]; te = TE[s]; M = Ms[s]
        r = np.random.random()
        if r < 0.25:
            # ---------------- insert segment
            if k < kmax - 1 and not FULL[s]:
                t0 = np.random.random() * beta
                if _contains(ts, te, k, t0) < 0:
                    lmax = _next_start_dist(ts, k, t0, beta) if k > 0 else beta
                    l = np.random.random() * lmax
                    t1 = t0 + l
                    if t1 >= beta:
                        t1 -= beta
                    ov = _overlap_with(TS[o], TE[o], K[o], FULL[o], t0, t1, beta)
                    for i in range(k):
                        Q[i] = _F(te[i] - t0, beta, tab, ntab)
                        R[i] = _F(t1 - ts[i], beta, tab, ntab)
                    lam = _F(t1 - t0, beta, tab, ntab)
                    if k > 0:
                        lam -= R[:k] @ (M[:k, :k] @ Q[:k])
                    A = beta * lmax / (k + 1) * np.exp(mus[s] * l - U * ov) * lam
                    if np.random.random() < abs(A):
                        _ins(M, k, Q, R, lam)
                        ts[k] = t0; te[k] = t1; K[s] = k + 1
        elif r < 0.5:
            # ---------------- remove segment
            if k > 0:
                p = np.random.randint(k)
                t0, t1 = ts[p], te[p]
                l = t1 - t0 if t1 >= t0 else beta - t0 + t1
                # lmax in the reverse move = distance from t0 to next start (other than p)
                best = beta
                for i in range(k):
                    if i != p:
                        d = ts[i] - t0
                        if d <= 0:
                            d += beta
                        if d < best:
                            best = d
                lmax = best
                ov = _overlap_with(TS[o], TE[o], K[o], FULL[o], t0, t1, beta)
                A = k / (beta * lmax) * np.exp(-mus[s] * l + U * ov) * M[p, p]
                if np.random.random() < abs(A):
                    _rem(M, p, k)
                    ts[p] = ts[k - 1]; te[p] = te[k - 1]; K[s] = k - 1
        elif r < 0.75:
            # ---------------- insert anti-segment (gap) inside an occupied region
            if k < kmax - 1:
                t0 = np.random.random() * beta
                a_idx = -1
                inside = FULL[s]
                if not inside:
                    a_idx = _contains(ts, te, k, t0)
                    inside = a_idx >= 0
                if inside:
                    lmax = beta if FULL[s] else _next_end_dist(te, k, t0, beta)
                    l = np.random.random() * lmax
                    t1 = t0 + l
                    if t1 >= beta:
                        t1 -= beta
                    ov = _overlap_with(TS[o], TE[o], K[o], FULL[o], t0, t1, beta)
                    # new end operator at t0 (row), new start operator at t1 (column)
                    for i in range(k):
                        Q[i] = _F(te[i] - t1, beta, tab, ntab)
                        R[i] = _F(t0 - ts[i], beta, tab, ntab)
                    lam = _F(t0 - t1, beta, tab, ntab)
                    if k > 0:
                        lam -= R[:k] @ (M[:k, :k] @ Q[:k])
                    A = beta * lmax / (k + 1) * np.exp(-mus[s] * l + U * ov) * lam
                    if np.random.random() < abs(A):
                        _ins(M, k, Q, R, lam)
                        if FULL[s]:
                            ts[k] = t1; te[k] = t0          # single wrapping segment
                            FULL[s] = False
                        else:
                            # segment a=[ts_a, te_a) -> [ts_a, t0) and [t1, te_a)
                            # F rows are ends: row a holds te_a, row k holds t0 -> swap so that
                            # index a pairs (ts_a, t0) and index k pairs (t1, te_a).
                            # swapping rows of F == swapping columns of M
                            for i in range(k + 1):
                                tmp = M[i, a_idx]; M[i, a_idx] = M[i, k]; M[i, k] = tmp
                            ts[k] = t1; te[k] = te[a_idx]; te[a_idx] = t0
                        K[s] = k + 1
        else:
            # ---------------- remove anti-segment: merge segment p with the next segment q
            if k > 0:
                p = np.random.randint(k)
                t0 = te[p]
                best = beta; q = -1
                for i in range(k):
                    d = ts[i] - t0
                    if d <= 0:
                        d += beta
                    if d < best:
                        best = d; q = i
                l = best
                t1 = ts[q]
                if k == 1:
                    lmax = beta
                else:
                    lmax = te[q] - t0
                    if lmax <= 0:
                        lmax += beta
                ov = _overlap_with(TS[o], TE[o], K[o], FULL[o], t0, t1, beta)
                # det ratio: remove row p (end t0) and column q (start t1): |M[q, p]|
                A = k / (beta * lmax) * np.exp(mus[s] * l - U * ov) * M[q, p]
                if np.random.random() < abs(A):
                    if k == 1:
                        K[s] = 0
                        FULL[s] = True
                    else:
                        newstart = ts[p]
                        newend = te[q]
                        nts = np.empty(k - 1); nte = np.empty(k - 1)
                        c = 0
                        for i in range(k):
                            if i == p or i == q:
                                continue
                            nts[c] = ts[i]; nte[c] = te[i]; c += 1
                        nts[c] = newstart; nte[c] = newend
                        for i in range(k - 1):
                            ts[i] = nts[i]; te[i] = nte[i]
                        Fm = np.empty((k - 1, k - 1))
                        for i in range(k - 1):
                            for j in range(k - 1):
                                Fm[i, j] = _F(te[i] - ts[j], beta, tab, ntab)
                        Mi = np.linalg.inv(Fm)
                        for i in range(k - 1):
                            for j in range(k - 1):
                                M[i, j] = Mi[i, j]
                        K[s] = k - 1
        # ---------------- measurement
        if step >= n_warm and (step - n_warm) % meas_every == 0:
            nmeas += 1
            for ss in range(2):
                kk = K[ss]
                hist[kk] += 1
                L = 0.0
                if FULL[ss]:
                    L = beta
                for i in range(kk):
                    l = TE[ss, i] - TS[ss, i]
                    if l < 0:
                        l += beta
                    L += l
                nacc[ss] += L / beta
                for iw in range(nw):
                    acc = 0.0 + 0.0j
                    for i in range(kk):
                        for j in range(kk):
                            dt = TE[ss, i] - TS[ss, j]
                            acc += Ms[ss, j, i] * np.exp(1j * w[iw] * dt)
                    Gacc[ss, iw] += acc / beta   # sign: F = Delta(beta - tau) convention
            dacc += _overlap_total(TS, TE, K, FULL, beta) / beta
    return Gacc / nmeas, nacc / nmeas, dacc / nmeas, hist / max(nmeas, 1)


@njit(cache=True)
def _overlap_total(TS, TE, K, FULL, beta):
    if FULL[0] and FULL[1]:
        return beta
    if FULL[0]:
        tot = 0.0
        for i in range(K[1]):
            l = TE[1, i] - TS[1, i]
            tot += l if l >= 0 else l + beta
        return tot
    if FULL[1]:
        tot = 0.0
        for i in range(K[0]):
            l = TE[0, i] - TS[0, i]
            tot += l if l >= 0 else l + beta
        return tot
    tot = 0.0
    for i in range(K[0]):
        for j in range(K[1]):
            tot += _seg_overlap(TS[0, i], TE[0, i], TS[1, j], TE[1, j], beta)
    return tot


class CTHYBSolver:
    """Single-orbital CT-HYB (segment) solver.

    solve(delta_iw, mu) with Delta(iw) the hybridization function and mu the
    chemical potential (impurity level eps_d = -mu).  Returns G, Sigma, <n>, <n_up n_dn>.
    """

    def __init__(self, grid: MatsubaraGrid, U: float, n_iw_meas: int = 64, n_tab: int = 10001):
        self.grid, self.U = grid, float(U)
        self.n_iw_meas = min(n_iw_meas, grid.n_iw)
        self.n_tab = n_tab

    def _Ftab(self, delta_iw):
        g = self.grid
        tau = np.linspace(0, g.beta, self.n_tab)
        c = MatsubaraGrid.tail_coefficients(g.iw, delta_iw, n_fit=min(60, g.n_iw // 4))
        dtau = g.iw_to_tau(delta_iw, c[0], c[1], c[2], tau=tau)
        # F(tau) = Delta(beta - tau) for 0<tau<beta
        return dtau[::-1].copy()

    def solve(self, delta_iw, mu, delta_dn=None, n_warm=20000, n_meas=2_000_000, meas_every=20,
              n_chains=4, seed=7, kmax=300):
        g = self.grid
        Fu = self._Ftab(delta_iw)
        Fd = self._Ftab(delta_iw if delta_dn is None else delta_dn)
        w = g.iw[: self.n_iw_meas].imag.copy()
        Gs, ns, ds = [], [], []
        hist = 0
        for c in range(n_chains):
            G, n, d, h = _run(g.beta, self.U, mu, mu, Fu, Fd, w, n_warm, n_meas, meas_every, seed + 977 * c, kmax)
            Gs.append(G); ns.append(n); ds.append(d); hist = hist + h
        Gs = np.array(Gs)
        G = Gs.mean(0)
        iw = g.iw[: self.n_iw_meas]
        G0inv = iw + mu - delta_iw[: self.n_iw_meas]
        sig = G0inv[None, :] - 1.0 / G
        Gsym = 0.5 * (G[0] + G[1])
        sig_chains = G0inv[None, :] - 1.0 / (0.5 * (Gs[:, 0] + Gs[:, 1]))
        err = (sig_chains.real.std(0, ddof=1) + 1j * sig_chains.imag.std(0, ddof=1)) / np.sqrt(n_chains) \
            if n_chains > 1 else np.full(len(iw), np.nan)
        return dict(iw=iw, G=Gsym, G_up=G[0], G_dn=G[1], Sigma=G0inv - 1.0 / Gsym, err_Sigma=err,
                    n=np.mean(ns, 0), docc=float(np.mean(ds)), order_histogram=hist / n_chains)
