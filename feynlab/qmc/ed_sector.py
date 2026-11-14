"""Sector-resolved exact diagonalization of the Anderson impurity model and an
ED-DMFT impurity solver (Caffarel & Krauth, PRL 72, 1545 (1994)).

The Hamiltonian conserves N_up and N_dn separately, so it is diagonalised block
by block; with 1 impurity + n_b bath sites the largest block is
C(n_s, n_s/2)^2 (n_s = 1 + n_b), e.g. 4900 for n_b = 7.

Bath parameters are obtained by fitting the hybridization function on a
fictitious Matsubara grid (beta_fit), as in Caffarel-Krauth.
"""
from __future__ import annotations

import itertools
from math import comb

import numpy as np
from scipy.optimize import minimize


def _states(ns, n):
    return [sum(1 << i for i in c) for c in itertools.combinations(range(ns), n)]


class SectorED:
    """H = eps_d n_d + U n_du n_dd + sum_l eps_l n_l + sum_l V_l (d^+ c_l + h.c.), per spin."""

    def __init__(self, eps_d, U, eps_bath, V_bath):
        self.eps_d, self.U = float(eps_d), float(U)
        self.eb = np.asarray(eps_bath, float)
        self.Vb = np.asarray(V_bath, float)
        self.ns = 1 + len(self.eb)
        self.onsite = np.concatenate([[self.eps_d], self.eb])
        self._diag_all()

    # ---------------------------------------------------------- one spin
    def _hop_matrix(self, states):
        """One-spin-species kinetic matrix in the given occupation basis."""
        idx = {s: i for i, s in enumerate(states)}
        n = len(states)
        T = np.zeros((n, n))
        for i, s in enumerate(states):
            T[i, i] = sum(self.onsite[k] for k in range(self.ns) if (s >> k) & 1)
            for l, V in enumerate(self.Vb, start=1):
                # d^+ c_l : site l occupied, site 0 empty
                if (s >> l) & 1 and not (s & 1):
                    t = s ^ (1 << l)
                    sgn = (-1) ** bin(s & ((1 << l) - 1)).count("1")  # remove c_l
                    t2 = t | 1  # add at 0: no sites below 0
                    j = idx[t2]
                    T[j, i] += V * sgn
                    T[i, j] += V * sgn
        return T

    def _diag_all(self):
        ns = self.ns
        self.sectors = {}
        st = [_states(ns, n) for n in range(ns + 1)]
        self.st = st
        hops = [self._hop_matrix(s) for s in st]
        for nu in range(ns + 1):
            for nd in range(ns + 1):
                du, dd = len(st[nu]), len(st[nd])
                H = np.kron(hops[nu], np.eye(dd)) + np.kron(np.eye(du), hops[nd])
                occ_u = np.array([s & 1 for s in st[nu]], float)
                occ_d = np.array([s & 1 for s in st[nd]], float)
                H += self.U * np.diag(np.kron(occ_u, occ_d))
                E, Vv = np.linalg.eigh(H)
                self.sectors[(nu, nd)] = (E, Vv)
        self.E0 = min(E.min() for E, _ in self.sectors.values())

    def _cdag_up_matrix(self, nu, nd):
        """<(nu+1, nd)| d_up^+ |(nu, nd)> in occupation basis."""
        a, b = self.st[nu], self.st[nu + 1]
        idx = {s: i for i, s in enumerate(b)}
        m = np.zeros((len(b), len(a)))
        for i, s in enumerate(a):
            if not (s & 1):
                m[idx[s | 1], i] = 1.0  # impurity is orbital 0 -> no sign within up species
        dd = len(self.st[nd])
        return np.kron(m, np.eye(dd))  # up-major ordering; up operators come first -> no extra sign

    # ---------------------------------------------------------- observables
    def green(self, z, beta, cutoff=1e-12):
        z = np.asarray(z, complex)
        G = np.zeros(len(z), complex)
        Zp = 0.0
        for (nu, nd), (E, V) in self.sectors.items():
            Zp += np.exp(-beta * (E - self.E0)).sum()
        for (nu, nd), (E, V) in self.sectors.items():
            if nu == self.ns:
                continue
            E2, V2 = self.sectors[(nu + 1, nd)]
            wn = np.exp(-beta * (E - self.E0))
            wm = np.exp(-beta * (E2 - self.E0))
            keep_n = wn > cutoff
            keep_m = wm > cutoff
            if not (keep_n.any() or keep_m.any()):
                continue
            M = V2.T @ self._cdag_up_matrix(nu, nd) @ V  # <m|d^+|n>
            A2 = np.abs(M) ** 2
            W = A2 * (wn[None, :] + wm[:, None])
            mask = W > cutoff * 1e-3
            dE = (E2[:, None] - E[None, :])[mask]
            Wm = W[mask]
            G += (Wm[None, :] / (z[:, None] - dE[None, :])).sum(1)
        return G / Zp

    def double_occupancy(self, beta):
        num = Zp = 0.0
        for (nu, nd), (E, V) in self.sectors.items():
            w = np.exp(-beta * (E - self.E0))
            occ_u = np.array([s & 1 for s in self.st[nu]], float)
            occ_d = np.array([s & 1 for s in self.st[nd]], float)
            dop = np.kron(occ_u, occ_d)
            num += np.sum(w * np.einsum("ik,i,ik->k", V, dop, V))
            Zp += w.sum()
        return num / Zp

    def hybridization(self, z):
        z = np.asarray(z, complex)
        return np.sum(self.Vb[None, :] ** 2 / (z[:, None] - self.eb[None, :]), axis=1)


# ----------------------------------------------------------------------------
def fit_bath(iw, delta, n_bath, x0=None, symmetric=True, n_fit=None, sum_rule=None, e_max=3.0, n_starts=8):
    """Least-squares fit Delta(iw) ~ sum_l V_l^2/(iw - e_l) (Caffarel-Krauth)."""
    n_fit = n_fit or len(iw)
    iw = iw[:n_fit]
    delta = delta[:n_fit]
    wts = 1.0 / iw.imag  # emphasise low frequencies

    if symmetric:
        npair = n_bath // 2
        has_zero = n_bath % 2 == 1

        def unpack(x):
            e = e_max * np.abs(np.tanh(x[:npair]))  # bounded bath levels
            V = x[npair:2 * npair]
            eb = np.concatenate([e, -e] + ([[0.0]] if has_zero else []))
            Vb = np.concatenate([V, V] + ([[x[-1]]] if has_zero else []))
            return eb, Vb
        if x0 is None:
            x0 = np.concatenate([np.arctanh(np.linspace(0.2, 1.5, npair) / e_max), np.full(npair, 0.3)]
                                + ([[0.3]] if has_zero else []))
    else:
        def unpack(x):
            return e_max * np.tanh(x[:n_bath]), x[n_bath:]
        if x0 is None:
            x0 = np.concatenate([np.arctanh(np.linspace(-1.5, 1.5, n_bath) / e_max), np.full(n_bath, 0.3)])

    def cost(x):
        eb, Vb = unpack(x)
        d = np.sum(Vb[None, :] ** 2 / (iw[:, None] - eb[None, :]), axis=1)
        return np.sum(wts * np.abs(d - delta) ** 2)

    if sum_rule is not None:
        # enforce sum_l V_l^2 = sum_rule (exact high-frequency weight of Delta)
        raw_unpack = unpack

        def unpack(x):  # noqa: F811
            eb, Vb = raw_unpack(x)
            return eb, Vb * np.sqrt(sum_rule / max(np.sum(Vb**2), 1e-300))

    # deterministic multi-start: the cost landscape has many local minima
    rng = np.random.default_rng(12345)
    starts = [np.asarray(x0, float)]
    for _ in range(n_starts - 1):
        starts.append(np.asarray(x0, float) + 0.5 * rng.standard_normal(len(x0)))
    best = None
    for xs in starts:
        r = minimize(cost, xs, method="BFGS", options=dict(gtol=1e-10, maxiter=5000))
        if best is None or r.fun < best.fun:
            best = r
    r = best
    eb, Vb = unpack(r.x)
    return eb, Vb, r.x, r.fun


# ============================================================================
# T = 0 Lanczos solver (sparse, ground state + continued fraction)
# ============================================================================
class LanczosED:
    """Zero-temperature Anderson impurity solver: Lanczos ground state in the
    (N_up, N_dn) sectors around half filling and the impurity Green function

        G(z) = <0| d (z - (H - E0))^{-1} d^+ |0> + <0| d^+ (z + (H - E0))^{-1} d |0>

    as Lanczos continued fractions (Dagotto, RMP 66, 763 (1994)), averaged over
    degenerate ground states.  Handles 1 + n_b sites with n_b up to ~9.
    """

    def __init__(self, eps_d, U, eps_bath, V_bath, n_lanczos: int = 300, sector_window: int = 1):
        import scipy.sparse as sp
        self.sp = sp
        self.eps_d, self.U = float(eps_d), float(U)
        self.eb = np.asarray(eps_bath, float)
        self.Vb = np.asarray(V_bath, float)
        self.ns = 1 + len(self.eb)
        self.onsite = np.concatenate([[self.eps_d], self.eb])
        self.nl = n_lanczos
        ns = self.ns
        self.st = [_states(ns, n) for n in range(ns + 1)]
        self._hop = {}
        self._H = {}
        # candidate ground-state sectors near half filling
        h = ns // 2
        cands = {(a, b) for a in range(max(0, h - sector_window), min(ns, h + sector_window + 1) + 1)
                 for b in range(max(0, h - sector_window), min(ns, h + sector_window + 1) + 1)}
        E = {}
        vecs = {}
        for s in cands:
            e, v = self._lowest(s)
            E[s], vecs[s] = e, v
        self.E0 = min(E.values())
        self.gs = [(s, vecs[s]) for s in E if E[s] - self.E0 < 1e-8]

    # ---------------------------------------------------------------- matrices
    def _hop_matrix(self, n):
        if n not in self._hop:
            self._hop[n] = SectorED._hop_matrix(self, self.st[n])
        return self._hop[n]

    def H(self, sec):
        if sec not in self._H:
            sp = self.sp
            nu, nd = sec
            du, dd = len(self.st[nu]), len(self.st[nd])
            Tu = sp.csr_matrix(self._hop_matrix(nu))
            Td = sp.csr_matrix(self._hop_matrix(nd))
            occ_u = np.array([s & 1 for s in self.st[nu]], float)
            occ_d = np.array([s & 1 for s in self.st[nd]], float)
            Hm = sp.kron(Tu, sp.identity(dd)) + sp.kron(sp.identity(du), Td) + \
                self.U * sp.diags(np.kron(occ_u, occ_d))
            self._H[sec] = Hm.tocsr()
        return self._H[sec]

    def _lowest(self, sec):
        from scipy.sparse.linalg import eigsh
        Hm = self.H(sec)
        if Hm.shape[0] <= 400:
            e, v = np.linalg.eigh(Hm.toarray())
            return e[0], v[:, 0]
        e, v = eigsh(Hm, k=1, which="SA", tol=1e-12)
        return e[0], v[:, 0]

    def _cdag(self, nu, nd):
        """sparse <(nu+1,nd)| d_up^+ |(nu,nd)>"""
        a, b = self.st[nu], self.st[nu + 1]
        idx = {s: i for i, s in enumerate(b)}
        rows, cols = [], []
        for i, s in enumerate(a):
            if not (s & 1):
                rows.append(idx[s | 1]); cols.append(i)
        m = self.sp.csr_matrix((np.ones(len(rows)), (rows, cols)), shape=(len(b), len(a)))
        return self.sp.kron(m, self.sp.identity(len(self.st[nd]))).tocsr()

    # ---------------------------------------------------------------- Lanczos
    def _cf(self, Hm, v0, n):
        """Lanczos coefficients (a_k, b_k) of Hm starting from v0 (norm returned)."""
        nrm = np.linalg.norm(v0)
        if nrm < 1e-14:
            return nrm, np.zeros(0), np.zeros(0)
        v = v0 / nrm
        vold = np.zeros_like(v)
        a, b = [], []
        beta = 0.0
        for k in range(min(n, Hm.shape[0])):
            w = Hm @ v
            alpha = float(np.dot(v, w))
            w = w - alpha * v - beta * vold
            a.append(alpha)
            beta = float(np.linalg.norm(w))
            if beta < 1e-12:
                break
            b.append(beta)
            vold, v = v, w / beta
        return nrm, np.array(a), np.array(b[: len(a) - 1])

    @staticmethod
    def _eval_cf(z, a, b):
        g = np.zeros_like(z)
        for k in range(len(a) - 1, -1, -1):
            g = 1.0 / (z - a[k] - (b[k] ** 2 * g if k < len(b) else 0.0))
        return g

    def green(self, z, beta=None):
        z = np.asarray(z, complex)
        G = np.zeros_like(z)
        for (nu, nd), psi in self.gs:
            if nu < self.ns:  # particle part: d^+ |0>
                cd = self._cdag(nu, nd)
                n2, a, b = self._cf(self.H((nu + 1, nd)), cd @ psi, self.nl)
                if n2 > 0:
                    G += n2**2 * self._eval_cf(z + self.E0, a, b)
            if nu > 0:        # hole part: d |0>
                c = self._cdag(nu - 1, nd).T
                n2, a, b = self._cf(self.H((nu - 1, nd)), c @ psi, self.nl)
                if n2 > 0:
                    G += n2**2 * (-self._eval_cf(-z + self.E0, a, b))
        return G / len(self.gs)

    def hybridization(self, z):
        z = np.asarray(z, complex)
        return np.sum(self.Vb[None, :] ** 2 / (z[:, None] - self.eb[None, :]), axis=1)
