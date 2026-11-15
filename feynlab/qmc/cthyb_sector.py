"""Sector-blocked, truncated matrix CT-HYB for up to 5 orbitals (10 flavours, 1024 atomic states).

Same configuration space, weights and estimators as `cthyb_matrix.MatrixCTHYB` (Werner & Millis PRB 74, 155107
(2006); Gull et al. RMP 83, 349 (2011)); the local trace is evaluated the way Haule (PRB 75, 155113 (2007)) does it:

* Quantum-number sectors.  H_loc is diagonalised inside blocks of fixed (N_up, N_dn); each c_f / c+_f maps one
  block to exactly one other block.  For an operator string the sequence of blocks is fixed by the starting
  block, so the trace is a sum over starting blocks of products of small dense blocks; strings that leave the
  Fock space (c on an empty flavour) are discarded before any arithmetic.  Cost per trace
  sum_{viable s0} K d_s0 d^2 instead of K 4^{n_orb} ... 4^{2 n_orb} for the dense trace.
* Energy truncation.  Atomic eigenstates with E - E_min > e_cut are dropped (Haule's truncation).  Convergence in
  e_cut is checked in the benchmarks (5 orbitals: error below the QMC noise for e_cut >= 2.5 U).
* Lazy trace bound (Semon et al. PRB 90, 075149 (2014)).  Before the exact trace, |Tr| is bounded by
  sum_{s0} d_s0 prod ||B|| e^{-dt E_min}; if the Metropolis ratio with the bound is below the random number the
  move is rejected without computing the trace.

Kanamori interactions (spin flip, pair hopping) conserve N_up and N_dn, so the blocking is exact for them.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid
from .multiorbital import _jw_ops, kanamori_hamiltonian
from .cthyb_matrix import _collect, _detF


class SectorStructure:
    """Atomic eigenbasis in (N_up, N_dn) blocks with truncation; flat arrays for the numba kernels."""

    def __init__(self, H, c_dense, nf, e_cut=None):
        dim = H.shape[0]
        states = np.arange(dim)
        nup = sum(((states >> f) & 1) for f in range(0, nf, 2))
        ndn = sum(((states >> f) & 1) for f in range(1, nf, 2))
        labels = nup * (nf + 1) + ndn
        keys = sorted(set(labels.tolist()))
        Hd = np.asarray(H)
        blk_of = {k: np.nonzero(labels == k)[0] for k in keys}
        off = np.abs(Hd[labels[:, None] != labels[None, :]]).max() if len(keys) > 1 else 0.0
        if off > 1e-12:
            raise ValueError("H_loc does not conserve (N_up, N_dn); use MatrixCTHYB")
        Es, Vs = [], []
        for k in keys:
            idx = blk_of[k]
            e, v = np.linalg.eigh(Hd[np.ix_(idx, idx)])
            Es.append(e); Vs.append(v)
        self.E0 = min(e.min() for e in Es)
        keep = []
        for e in Es:
            keep.append(np.nonzero(e - self.E0 <= (np.inf if e_cut is None else e_cut))[0])
        # drop empty blocks
        sec = [i for i in range(len(keys)) if len(keep[i]) > 0]
        self.keys = [keys[i] for i in sec]
        self.n_states_total = dim
        self.n_kept = int(sum(len(keep[i]) for i in sec))
        nsec = len(sec)
        self.nsec = nsec
        self.dims = np.array([len(keep[i]) for i in sec], np.int64)
        self.e_off = np.zeros(nsec + 1, np.int64)
        self.e_off[1:] = np.cumsum(self.dims)
        self.E = np.concatenate([Es[i][keep[i]] - self.E0 for i in sec])
        self.Emin = np.array([Es[i][keep[i]].min() - self.E0 for i in sec])
        self.U = [Vs[i][:, keep[i]] for i in sec]               # Fock -> kept eigenvectors of each block
        self.idx = [blk_of[keys[i]] for i in sec]
        pos = {keys[i]: n for n, i in enumerate(sec)}
        # operator blocks: kind 0 = c_f, 1 = c+_f
        tgt = -np.ones((nf, 2, nsec), np.int64)
        boff = -np.ones((nf, 2, nsec), np.int64)
        bnorm = np.zeros((nf, 2, nsec))
        flat = []
        cur = 0
        for f in range(nf):
            for kind in (0, 1):
                C = c_dense[f] if kind == 0 else c_dense[f].T
                dup = (1 if f % 2 == 0 else 0) * (1 if kind == 1 else -1)
                ddn = (1 if f % 2 == 1 else 0) * (1 if kind == 1 else -1)
                for s in range(nsec):
                    k = self.keys[s]
                    kt = k + dup * (nf + 1) + ddn
                    if kt not in pos:
                        continue
                    t = pos[kt]
                    B = self.U[t].T @ C[np.ix_(self.idx[t], self.idx[s])] @ self.U[s]
                    if np.abs(B).max() < 1e-14:
                        continue
                    tgt[f, kind, s] = t
                    boff[f, kind, s] = cur
                    bnorm[f, kind, s] = np.linalg.norm(B, 2)
                    flat.append(np.ascontiguousarray(B).ravel())
                    cur += B.size
        self.tgt, self.boff, self.bnorm = tgt, boff, bnorm
        self.blocks = np.concatenate(flat) if flat else np.zeros(0)

    def observable(self, Op):
        """Fock operator (block diagonal in (N_up, N_dn)) -> flat array of kept blocks, offsets o_off."""
        Op = np.asarray(Op.toarray() if hasattr(Op, "toarray") else Op)
        out = []
        for s in range(self.nsec):
            out.append((self.U[s].T @ Op[np.ix_(self.idx[s], self.idx[s])] @ self.U[s]).ravel())
        return np.concatenate(out)

    @property
    def o_off(self):
        o = np.zeros(self.nsec + 1, np.int64)
        o[1:] = np.cumsum(self.dims**2)
        return o


@njit(cache=True)
def _order_parity(times):
    K = times.shape[0]
    order = np.argsort(-times)
    inv = 0
    for a in range(K):
        for b in range(a + 1, K):
            if order[a] > order[b]:
                inv += 1
    return order, (-1.0 if inv % 2 == 1 else 1.0)


@njit(cache=True)
def _sector_bounds(order, times, kinds, flav, beta, dims, Emin, tgt, bnorm):
    """Upper bound of |tr P_s0| for every starting block s0 (0 if the string annihilates s0)."""
    K = times.shape[0]
    nsec = dims.shape[0]
    bd = np.zeros(nsec)
    for s0 in range(nsec):
        s = s0
        b = 1.0 * dims[s0]
        tprev = 0.0
        ok = True
        for a in range(K - 1, -1, -1):
            o = order[a]
            b *= np.exp(-(times[o] - tprev) * Emin[s])
            b *= bnorm[flav[o], kinds[o], s]
            s = tgt[flav[o], kinds[o], s]
            if s < 0:
                ok = False
                break
            tprev = times[o]
        if ok and s == s0:
            bd[s0] = b * np.exp(-(beta - tprev) * Emin[s])
    return bd


@njit(cache=True)
def _propagate(order, times, kinds, flav, beta, E, e_off, dims, tgt, boff, blocks, s0, ncol):
    """P = e^{-(beta - t_1) H} O_1 ... O_K e^{-t_K H} restricted to block s0, first ncol (lowest) columns."""
    K = times.shape[0]
    d0 = dims[s0]
    s = s0
    ds = d0
    P = np.zeros((d0, ncol))
    for i in range(ncol):
        P[i, i] = 1.0
    tprev = 0.0
    for a in range(K - 1, -1, -1):
        o = order[a]
        dt = times[o] - tprev
        for i in range(ds):
            fac = np.exp(-dt * E[e_off[s] + i])
            for j in range(ncol):
                P[i, j] *= fac
        f = flav[o]; kd = kinds[o]
        t = tgt[f, kd, s]
        dt_ = dims[t]
        B = blocks[boff[f, kd, s]: boff[f, kd, s] + dt_ * ds].reshape(dt_, ds)
        P = B @ P
        s = t
        ds = dt_
        tprev = times[o]
    dt = beta - tprev
    for i in range(ds):
        fac = np.exp(-dt * E[e_off[s] + i])
        for j in range(ncol):
            P[i, j] *= fac
    return P


@njit(cache=True)
def _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, obs, o_off, use_obs,
               bound_only, tol):
    """Sector trace (sign of the time-ordering permutation included). bound_only: return the upper bound of |Tr|.
    Screening, all with rigorous bounds (relative error <= n_states x tol):
      blocks with bound_s0 < tol x total are skipped;
      inside a block only the outer eigenstates n with bound_s0 e^{-l (E_n - E_min)} / d_s0 >= tol x total are
      propagated, l = beta - t_first + t_last being the wrap-around interval (matrix-vector instead of
      matrix-matrix products)."""
    K = times.shape[0]
    order, sgn = _order_parity(times)
    bd = _sector_bounds(order, times, kinds, flav, beta, dims, Emin, tgt, bnorm)
    btot = bd.sum()
    if bound_only:
        return btot
    ell = beta
    if K > 0:
        ell = beta - times[order[0]] + times[order[K - 1]]
    nsec = dims.shape[0]
    total = 0.0
    for s0 in range(nsec):
        if bd[s0] == 0.0 or bd[s0] < tol * btot:
            continue
        d0 = dims[s0]
        ncol = 0
        for n in range(d0):
            if bd[s0] / d0 * np.exp(-ell * (E[e_off[s0] + n] - Emin[s0])) >= tol * btot:
                ncol = n + 1
        if ncol == 0:
            continue
        P = _propagate(order, times, kinds, flav, beta, E, e_off, dims, tgt, boff, blocks, s0, ncol)
        for i in range(ncol):
            total += P[i, i]
    return sgn * total


@njit(cache=True)
def _trace_obs(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, obs, o_off, tol):
    """Tr[P O_k] / sign for all observables O_k inserted at tau = 0 (full columns, block screening only)."""
    K = times.shape[0]
    order, sgn = _order_parity(times)
    bd = _sector_bounds(order, times, kinds, flav, beta, dims, Emin, tgt, bnorm)
    btot = bd.sum()
    nobs = obs.shape[0]
    out = np.zeros(nobs)
    for s0 in range(dims.shape[0]):
        if bd[s0] == 0.0 or bd[s0] < tol * btot:
            continue
        d0 = dims[s0]
        P = _propagate(order, times, kinds, flav, beta, E, e_off, dims, tgt, boff, blocks, s0, d0)
        for k in range(nobs):
            Ob = obs[k, o_off[s0]: o_off[s0 + 1]].reshape(d0, d0)
            acc = 0.0
            for i in range(d0):
                for j in range(d0):
                    acc += P[i, j] * Ob[j, i]
            out[k] += acc
    return sgn * out


@njit(cache=True)
def _run_sec(beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, Ftab, w, n_warm, n_meas, meas_every, seed, kmax,
             obs, o_off, lazy, tol):
    np.random.seed(seed)
    nf = Ftab.shape[0]
    nw = w.shape[0]
    nobs = obs.shape[0]
    TS = np.zeros((nf, kmax)); TE = np.zeros((nf, kmax))
    K = np.zeros(nf, np.int64)
    dets = np.ones(nf)
    dummy = np.zeros(1)
    times, kinds, flav = _collect(TS, TE, K, nf)
    tr = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, dummy, o_off,
                    False, False, tol)
    Gacc = np.zeros((nf, nw), dtype=np.complex128)
    oacc = np.zeros(nobs)
    sacc = 0.0
    hist = np.zeros(kmax + 1)
    nmeas = 0
    n_lazy = 0
    n_trace = 0
    for step in range(n_warm + n_meas):
        f = np.random.randint(nf)
        k = K[f]
        u = np.random.random()
        trn = tr
        if np.random.random() < 0.5 and k < kmax - 1:
            TS[f, k] = np.random.random() * beta
            TE[f, k] = np.random.random() * beta
            K[f] = k + 1
            dnew, _ = _detF(TS[f], TE[f], k + 1, beta, Ftab[f])
            times, kinds, flav = _collect(TS, TE, K, nf)
            pre = (beta * beta / ((k + 1) * (k + 1))) * abs(dnew) / abs(dets[f] * tr)
            acc = False
            if dnew != 0.0:
                if lazy:
                    bd = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, dummy,
                                    o_off, False, True, tol)
                    if pre * bd < u:
                        n_lazy += 1
                    else:
                        trn = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks,
                                         dummy, o_off, False, False, tol)
                        n_trace += 1
                        acc = u < pre * abs(trn)
                else:
                    trn = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, dummy,
                                     o_off, False, False, tol)
                    n_trace += 1
                    acc = u < pre * abs(trn)
            if acc:
                dets[f] = dnew; tr = trn
            else:
                K[f] = k
        elif k > 0:
            i = np.random.randint(k)
            j = np.random.randint(k)
            sts = TS[f, i]; ste = TE[f, j]
            TS[f, i] = TS[f, k - 1]; TS[f, k - 1] = sts
            TE[f, j] = TE[f, k - 1]; TE[f, k - 1] = ste
            K[f] = k - 1
            dnew, _ = _detF(TS[f], TE[f], k - 1, beta, Ftab[f])
            times, kinds, flav = _collect(TS, TE, K, nf)
            pre = (k * k / (beta * beta)) * abs(dnew) / abs(dets[f] * tr)
            acc = False
            if dnew != 0.0:
                if lazy:
                    bd = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, dummy,
                                    o_off, False, True, tol)
                    if pre * bd < u:
                        n_lazy += 1
                    else:
                        trn = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks,
                                         dummy, o_off, False, False, tol)
                        n_trace += 1
                        acc = u < pre * abs(trn)
                else:
                    trn = _trace_sec(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, dummy,
                                     o_off, False, False, tol)
                    n_trace += 1
                    acc = u < pre * abs(trn)
            if acc:
                dets[f] = dnew; tr = trn
            else:
                K[f] = k
                TS[f, k - 1] = TS[f, i]; TS[f, i] = sts
                TE[f, k - 1] = TE[f, j]; TE[f, j] = ste
        if step >= n_warm and (step - n_warm) % meas_every == 0:
            nmeas += 1
            s = 1.0
            for ff in range(nf):
                if dets[ff] < 0:
                    s = -s
            if tr < 0:
                s = -s
            sacc += s
            for ff in range(nf):
                kk = K[ff]
                hist[kk] += 1.0 / nf
                if kk > 0:
                    _, Fm = _detF(TS[ff], TE[ff], kk, beta, Ftab[ff])
                    M = np.linalg.inv(Fm)
                    for iw in range(nw):
                        a_ = 0.0 + 0.0j
                        for a in range(kk):
                            for b in range(kk):
                                a_ += M[b, a] * np.exp(1j * w[iw] * (TE[ff, a] - TS[ff, b]))
                        Gacc[ff, iw] -= s * a_ / beta
            if nobs > 0:
                times, kinds, flav = _collect(TS, TE, K, nf)
                ov = _trace_obs(times, kinds, flav, beta, E, e_off, dims, Emin, tgt, boff, bnorm, blocks, obs, o_off,
                                tol)
                for o in range(nobs):
                    oacc[o] += s * ov[o] / tr
    return Gacc / nmeas, oacc / nmeas, sacc / nmeas, hist / max(nmeas, 1), n_lazy, n_trace


class SectorCTHYB:
    """Matrix CT-HYB with (N_up, N_dn)-blocked, energy-truncated local trace and lazy bounds; n_orb <= 5.

    Same interface and estimators as MatrixCTHYB:
        SectorCTHYB(grid, n_orb, U, J, eps, mu, Up=None, H_loc=None, e_cut=None).solve(delta_iw, ...)
    e_cut: keep atomic eigenstates with E - E_min <= e_cut (None = all)."""

    def __init__(self, grid: MatsubaraGrid, n_orb: int, U: float, J: float, eps=None, mu: float = 0.0, Up=None,
                 H_loc=None, n_iw_meas: int = 32, n_tab: int = 10001, e_cut=None):
        self.grid, self.n_orb = grid, n_orb
        self.nf = 2 * n_orb
        c = _jw_ops(self.nf)
        if H_loc is None:
            H = kanamori_hamiltonian(c, n_orb, U, J, Up)
            eps = np.zeros(n_orb) if eps is None else np.asarray(eps, float)
            for a in range(n_orb):
                for s in (0, 1):
                    k = 2 * a + s
                    H = H + (eps[a] - mu) * c[k].T @ c[k]
            H = H.toarray()
        else:
            H = np.asarray(H_loc, float)
        self.c = c
        self.c_dense = [ci.toarray() for ci in c]
        self.S = SectorStructure(H, self.c_dense, self.nf, e_cut)
        self.n_iw_meas = min(n_iw_meas, grid.n_iw)
        self.n_tab = n_tab

    def default_observables(self):
        n = [ci.T @ ci for ci in self.c]
        obs, names = [], []
        for f in range(self.nf):
            obs.append(n[f]); names.append(f"n_{f}")
        for a in range(self.n_orb):
            obs.append(n[2 * a] @ n[2 * a + 1]); names.append(f"docc_{a}")
        Sz = 0.5 * sum(n[2 * a] - n[2 * a + 1] for a in range(self.n_orb))
        obs.append(Sz @ Sz); names.append("Sz^2")
        return obs, names

    def _Ftab(self, delta_iw):
        g = self.grid
        tau = np.linspace(0, g.beta, self.n_tab)
        cf = MatsubaraGrid.tail_coefficients(g.iw, delta_iw, n_fit=min(60, g.n_iw // 4))
        dtau = g.iw_to_tau(delta_iw, cf[0], cf[1], cf[2], tau=tau)
        return -dtau[::-1].copy()

    def solve(self, delta_iw, n_warm=20000, n_meas=400_000, meas_every=25, n_chains=4, seed=11, kmax=80,
              observables=None, parallel="serial", lazy=True, tol=1e-12):
        g = self.grid
        D = np.asarray(delta_iw)
        if D.ndim == 1:
            D = np.tile(D, (self.nf, 1))
        Ftab = np.ascontiguousarray(np.array([self._Ftab(D[f]) for f in range(self.nf)]))
        ops, names = self.default_observables() if observables is None else (observables[0], list(observables[1]))
        S = self.S
        obs = np.ascontiguousarray(np.array([S.observable(o) for o in ops])) if ops else np.zeros((0, 1))
        w = g.iw[: self.n_iw_meas].imag.copy()
        from ..utils.parallel import map_chains
        args = [(g.beta, S.E, S.e_off, S.dims, S.Emin, S.tgt, S.boff, S.bnorm, S.blocks, Ftab, w, n_warm, n_meas,
                 meas_every, seed + 7919 * ch, kmax, obs, S.o_off, lazy, tol) for ch in range(n_chains)]
        Gs, Os, Ss, H, nl, nt = [], [], [], 0, 0, 0
        for Gc, oc, sc, h, a, b in map_chains(_run_sec, args, backend=parallel):
            Gs.append(Gc / sc); Os.append(oc / sc); Ss.append(sc); H = H + h; nl += a; nt += b
        Gs, Os = np.array(Gs), np.array(Os)
        se = lambda x: x.std(0, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.full(x.shape[1:], np.nan)
        out = dict(iw=g.iw[: self.n_iw_meas], G=Gs.mean(0), G_err=se(Gs.real) + 1j * se(Gs.imag),
                   sign=float(np.mean(Ss)), order_histogram=H / n_chains,
                   lazy_rejected=nl, traces_computed=nt, n_states_kept=S.n_kept, n_sectors=S.nsec)
        for k, nm in enumerate(names):
            out[nm] = float(Os[:, k].mean())
            out[nm + "_err"] = float(se(Os[:, k]))
        return out
