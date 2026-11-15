"""Matrix (general local Hamiltonian) CT-HYB: hybridization expansion with the full local trace.

Werner & Millis, PRB 74, 155107 (2006); Haule, PRB 75, 155113 (2007); Gull et al. RMP 83, 349 (2011) Sec. III.

Removes the density-density restriction of the segment solver: any local Hamiltonian H_loc on n_f flavours
(spin-orbitals) is allowed, in particular the full rotationally invariant Slater-Kanamori interaction with
spin-flip and pair-hopping terms (Hund's metals).  Hybridization is flavour-diagonal, Delta_f(iw).

Configuration: for each flavour f, k_f creators c+_f(ts) and annihilators c_f(te).
Weight
    w(C) = prod_f det F_f  *  T(C),        F_f[i, j] = F_f(te_i - ts_j),   F(tau) = -Delta(beta - tau) = Delta(-tau)
    T(C) = s_P Tr[ e^{-(beta - t_1) H} O_1 e^{-(t_1 - t_2) H} O_2 ... O_K e^{-t_K H} ] / Tr e^{-beta H}
with all operators sorted by decreasing time and s_P the parity of the permutation from the reference order
(f = 0, 1, ...; pairs c_f(te_i) c+_f(ts_i)) to time order; with this convention T reduces to the segment weight
exp(mu L - U O) for density-density interactions.  Moves: insertion/removal of a (c+, c) pair of one flavour,
and global flavour-swap moves for degenerate flavours (optional).  The sign of w is tracked.
Estimators:  G_f(iw) = -(1/beta) <s sum_ij M_ji e^{iw (te_i - ts_j)}> / <s>   (M = F^{-1}),
             <O> for local operators by insertion at tau = 0 (densities, double occupancies, Hund's moments).
Cost: dimension of the local Fock space (2^n_f = 16 for two orbitals, 64 for three) x expansion order.
"""
from __future__ import annotations

import numpy as np
from numba import njit

from ..core.matsubara import MatsubaraGrid
from .multiorbital import _jw_ops, kanamori_hamiltonian


@njit(cache=True)
def _F(t, beta, tab):
    n = tab.shape[0]
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
def _trace(times, kinds, flav, E, Cop, Cdag, beta, insert_op):
    """s_P * Tr[...] / Tr e^{-beta H} in the eigenbasis (E shifted so that min E = 0).
    kinds: 0 = annihilator, 1 = creator; reference order = as given in the arrays."""
    K = times.shape[0]
    order = np.argsort(-times)
    # permutation parity (inversions)
    inv = 0
    for a in range(K):
        for b in range(a + 1, K):
            if order[a] > order[b]:
                inv += 1
    sgn = -1.0 if inv % 2 == 1 else 1.0
    D = E.shape[0]
    # start with e^{-t_K H} (rightmost) applied to identity, build P = product
    P = np.eye(D)
    tprev = 0.0
    for a in range(K - 1, -1, -1):
        o = order[a]
        t = times[o]
        dt = t - tprev
        for i in range(D):
            fac = np.exp(-dt * E[i])
            for j in range(D):
                P[i, j] *= fac
        if kinds[o] == 0:
            P = Cop[flav[o]] @ P
        else:
            P = Cdag[flav[o]] @ P
        tprev = t
    dt = beta - tprev
    for i in range(D):
        fac = np.exp(-dt * E[i])
        for j in range(D):
            P[i, j] *= fac
    tr = 0.0
    if insert_op.shape[0] == D:
        # Tr[ P O ] with O inserted at tau = 0 (rightmost)
        Q = P @ insert_op
        for i in range(D):
            tr += Q[i, i]
    else:
        for i in range(D):
            tr += P[i, i]
    return sgn * tr


@njit(cache=True)
def _collect(TS, TE, K, nf):
    tot = 0
    for f in range(nf):
        tot += 2 * K[f]
    times = np.empty(tot)
    kinds = np.empty(tot, np.int64)
    flav = np.empty(tot, np.int64)
    c = 0
    for f in range(nf):
        for i in range(K[f]):
            times[c] = TE[f, i]; kinds[c] = 0; flav[c] = f; c += 1
            times[c] = TS[f, i]; kinds[c] = 1; flav[c] = f; c += 1
    return times, kinds, flav


@njit(cache=True)
def _detF(TS, TE, k, beta, tab):
    if k == 0:
        return 1.0, np.zeros((0, 0))
    Fm = np.empty((k, k))
    for i in range(k):
        for j in range(k):
            Fm[i, j] = _F(TE[i] - TS[j], beta, tab)
    d = np.linalg.det(Fm)
    return d, Fm


@njit(cache=True)
def _run(beta, E, Cop, Cdag, Ftab, w, n_warm, n_meas, meas_every, seed, kmax, obs, Z0):
    np.random.seed(seed)
    nf = Ftab.shape[0]
    nw = w.shape[0]
    nobs = obs.shape[0]
    TS = np.zeros((nf, kmax)); TE = np.zeros((nf, kmax))
    K = np.zeros(nf, np.int64)
    dets = np.ones(nf)
    empty = np.zeros((0, 0))
    times, kinds, flav = _collect(TS, TE, K, nf)
    tr = _trace(times, kinds, flav, E, Cop, Cdag, beta, empty)
    Gacc = np.zeros((nf, nw), dtype=np.complex128)
    oacc = np.zeros(nobs)
    sacc = 0.0
    hist = np.zeros(kmax + 1)
    nmeas = 0
    total = n_warm + n_meas
    for step in range(total):
        f = np.random.randint(nf)
        k = K[f]
        if np.random.random() < 0.5 and k < kmax - 1:
            ts = np.random.random() * beta
            te = np.random.random() * beta
            TS[f, k] = ts; TE[f, k] = te
            K[f] = k + 1
            dnew, _ = _detF(TS[f], TE[f], k + 1, beta, Ftab[f])
            times, kinds, flav = _collect(TS, TE, K, nf)
            trn = _trace(times, kinds, flav, E, Cop, Cdag, beta, empty)
            R = (beta * beta / ((k + 1) * (k + 1))) * abs(dnew * trn) / abs(dets[f] * tr)
            if np.random.random() < R:
                dets[f] = dnew; tr = trn
            else:
                K[f] = k
        elif k > 0:
            # (no 'continue' anywhere: every step must reach the measurement, otherwise configurations
            #  where a move is impossible, e.g. removal at k = 0, would be under-sampled)
            i = np.random.randint(k)   # creator to remove
            j = np.random.randint(k)   # annihilator to remove
            sts = TS[f, i]; ste = TE[f, j]
            TS[f, i] = TS[f, k - 1]; TS[f, k - 1] = sts
            TE[f, j] = TE[f, k - 1]; TE[f, k - 1] = ste
            K[f] = k - 1
            dnew, _ = _detF(TS[f], TE[f], k - 1, beta, Ftab[f])
            times, kinds, flav = _collect(TS, TE, K, nf)
            trn = _trace(times, kinds, flav, E, Cop, Cdag, beta, empty)
            R = (k * k / (beta * beta)) * abs(dnew * trn) / abs(dets[f] * tr)
            if np.random.random() < R:
                dets[f] = dnew; tr = trn
            else:
                # restore the original ordering: stored det and trace signs refer to it
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
                        acc = 0.0 + 0.0j
                        for a in range(kk):
                            for b in range(kk):
                                acc += M[b, a] * np.exp(1j * w[iw] * (TE[ff, a] - TS[ff, b]))
                        Gacc[ff, iw] -= s * acc / beta
            times, kinds, flav = _collect(TS, TE, K, nf)
            for o in range(nobs):
                oacc[o] += s * _trace(times, kinds, flav, E, Cop, Cdag, beta, obs[o]) / tr
    return Gacc / nmeas, oacc / nmeas, sacc / nmeas, hist / max(nmeas, 1)


class MatrixCTHYB:
    """General-H_loc CT-HYB for n_orb orbitals x 2 spins (flavour index 2a + s).

    H_loc = Kanamori(U, J, U') + sum_a (eps_a - mu) n_a  (or a user matrix via H_loc=...).
    solve(delta_iw) with delta_iw of shape (n_flavour, n_iw) or (n_iw,) (same for all flavours)."""

    def __init__(self, grid: MatsubaraGrid, n_orb: int, U: float, J: float, eps=None, mu: float = 0.0, Up=None,
                 H_loc=None, n_iw_meas: int = 32, n_tab: int = 10001):
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
        E, V = np.linalg.eigh(H)
        self.E0 = E.min()
        self.E = E - self.E0
        self.V = V
        self.c_dense = [ci.toarray() for ci in c]
        self.Cop = np.ascontiguousarray(np.array([V.T @ ci @ V for ci in self.c_dense]))
        self.Cdag = np.ascontiguousarray(np.array([V.T @ ci.T @ V for ci in self.c_dense]))
        self.n_iw_meas = min(n_iw_meas, grid.n_iw)
        self.n_tab = n_tab

    def local_operator(self, Op):
        """Fock-basis operator -> eigenbasis."""
        return self.V.T @ Op @ self.V

    def default_observables(self):
        n = [ci.T @ ci for ci in self.c_dense]
        obs, names = [], []
        for f in range(self.nf):
            obs.append(n[f]); names.append(f"n_{f}")
        for a in range(self.n_orb):
            obs.append(n[2 * a] @ n[2 * a + 1]); names.append(f"docc_{a}")
        Sz = 0.5 * sum(n[2 * a] - n[2 * a + 1] for a in range(self.n_orb))
        obs.append(Sz @ Sz); names.append("Sz^2")
        return [self.local_operator(o) for o in obs], names

    def _Ftab(self, delta_iw):
        g = self.grid
        tau = np.linspace(0, g.beta, self.n_tab)
        cf = MatsubaraGrid.tail_coefficients(g.iw, delta_iw, n_fit=min(60, g.n_iw // 4))
        dtau = g.iw_to_tau(delta_iw, cf[0], cf[1], cf[2], tau=tau)
        return -dtau[::-1].copy()          # F(tau) = -Delta(beta - tau) = Delta(-tau): positive weights

    def solve(self, delta_iw, n_warm=20000, n_meas=400_000, meas_every=25, n_chains=4, seed=11, kmax=60,
              observables=None, parallel="serial"):
        """parallel: 'serial' | 'process' | 'mpi' | 'auto' (feynlab.utils.parallel.map_chains)."""
        g = self.grid
        D = np.asarray(delta_iw)
        if D.ndim == 1:
            D = np.tile(D, (self.nf, 1))
        Ftab = np.ascontiguousarray(np.array([self._Ftab(D[f]) for f in range(self.nf)]))
        if observables is None:
            obs, names = self.default_observables()
        else:
            obs, names = [self.local_operator(o) for o in observables[0]], list(observables[1])
        obs = np.ascontiguousarray(np.array(obs))
        w = g.iw[: self.n_iw_meas].imag.copy()
        from ..utils.parallel import map_chains
        args = [(g.beta, self.E, self.Cop, self.Cdag, Ftab, w, n_warm, n_meas, meas_every, seed + 7919 * ch, kmax, obs, 1.0)
                for ch in range(n_chains)]
        Gs, Os, Ss, H = [], [], [], 0
        for Gc, oc, sc, h in map_chains(_run, args, backend=parallel):
            Gs.append(Gc / sc); Os.append(oc / sc); Ss.append(sc); H = H + h
        Gs, Os = np.array(Gs), np.array(Os)
        se = lambda x: x.std(0, ddof=1) / np.sqrt(len(x)) if len(x) > 1 else np.full(x.shape[1:], np.nan)
        G = Gs.mean(0)
        iw = g.iw[: self.n_iw_meas]
        out = dict(iw=iw, G=G, G_err=se(Gs.real) + 1j * se(Gs.imag), sign=float(np.mean(Ss)),
                   order_histogram=H / n_chains)
        for k, nm in enumerate(names):
            out[nm] = float(Os[:, k].mean())
            out[nm + "_err"] = float(se(Os[:, k]))
        return out
