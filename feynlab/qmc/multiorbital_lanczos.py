"""Zero-temperature Lanczos exact diagonalization for multi-orbital Kanamori impurities (ED-DMFT solver).

Impurity: n_orb orbitals with the full Kanamori interaction (spin flip and pair hopping included); each orbital
couples to its own bath levels (diagonal hybridization).  The Hamiltonian conserves N_up and N_dn, so the ground
state is found sector by sector with ARPACK (scipy eigsh).  Green functions are Lanczos continued fractions
(Dagotto, RMP 66, 763 (1994)), averaged over the (orbitally and spin) degenerate ground-state manifold.

    G_a(z) = < c_a (z - (H - E0))^{-1} c+_a > + < c+_a (z + (H - E0))^{-1} c_a >

With 3 orbitals and 2 bath sites per orbital (18 spin-orbitals) the largest sector has 15876 states, so a full
three-orbital t2g DMFT step takes seconds; T = 0 is appropriate when k_B T is far below the coherence scale.
"""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp
from scipy.sparse.linalg import eigsh

from .multiorbital import _jw_ops, kanamori_hamiltonian


def _continued_fraction(H, v0, z, sign, E0, n_lanczos=200):
    """<v0| (z - sign (H - E0))^{-1} |v0> by Lanczos."""
    nrm = np.linalg.norm(v0)
    if nrm < 1e-14:
        return np.zeros(len(z), complex)
    v = v0 / nrm
    a, b = [], [0.0]
    vprev = np.zeros_like(v)
    for _ in range(min(n_lanczos, H.shape[0])):
        w = H @ v
        al = float(np.real(np.vdot(v, w)))
        w = w - al * v - b[-1] * vprev
        a.append(al)
        bn = np.linalg.norm(w)
        if bn < 1e-12:
            break
        b.append(bn)
        vprev, v = v, w / bn
    a = np.array(a) - E0
    b = np.array(b[1:len(a)])
    g = np.zeros(len(z), complex)
    for i in range(len(a) - 1, -1, -1):
        bb = b[i] ** 2 if i < len(b) else 0.0
        g = 1.0 / (z - sign * a[i] - bb * g)
    return nrm**2 * g


class MultiOrbitalLanczos:
    def __init__(self, n_orb, U, J, eps_d, eps_bath, V_bath, Up=None, n_lanczos=200, degeneracy_tol=1e-7,
                 n_eig=6, beta=None, boltzmann_cut=1e-5):
        """beta=None: T = 0 (degenerate ground-state manifold).  beta given: low-temperature Lanczos, all states
        with Boltzmann weight > boltzmann_cut among the n_eig lowest of every sector are included (continuous in
        the bath parameters, which stabilises DMFT loops with degenerate orbitals)."""
        self.beta, self.bcut = beta, boltzmann_cut
        self.n_orb = n_orb
        eb = [np.atleast_1d(np.asarray(e, float)) for e in eps_bath]
        Vb = [np.atleast_1d(np.asarray(v, float)) for v in V_bath]
        nsites = n_orb + sum(len(e) for e in eb)
        self.n_so = 2 * nsites
        if self.n_so > 22:
            raise ValueError("MultiOrbitalLanczos: at most 11 sites")
        c = _jw_ops(self.n_so)
        self.c = c
        H = kanamori_hamiltonian(c, n_orb, U, J, Up)
        site = n_orb
        for a in range(n_orb):
            for s in (0, 1):
                H = H + eps_d[a] * (c[2 * a + s].T @ c[2 * a + s])
            for e, V in zip(eb[a], Vb[a]):
                for s in (0, 1):
                    k = 2 * site + s
                    H = H + e * (c[k].T @ c[k])
                    hop = V * (c[2 * a + s].T @ c[k])
                    H = H + hop + hop.T
                site += 1
        self.H = sp.csr_matrix(H)
        dim = 2**self.n_so
        states = np.arange(dim)
        nup = np.zeros(dim, int)
        ndn = np.zeros(dim, int)
        for i in range(self.n_so):
            bit = (states >> i) & 1
            if i % 2 == 0:
                nup += bit
            else:
                ndn += bit
        self.nup, self.ndn = nup, ndn
        self.n_lanczos = n_lanczos
        self.eb, self.Vb = eb, Vb
        self._ground_state(degeneracy_tol, n_eig)

    def _sector(self, nu, nd):
        return np.where((self.nup == nu) & (self.ndn == nd))[0]

    def _ground_state(self, tol, n_eig):
        ns = self.n_so // 2
        best = []
        for nu in range(ns + 1):
            for nd in range(ns + 1):
                if abs(nu - nd) > 2:
                    continue
                sel = self._sector(nu, nd)
                Hs = self.H[sel][:, sel]
                k = min(n_eig, len(sel) - 1) if len(sel) > 1 else 1
                if len(sel) <= 400:
                    E, V = np.linalg.eigh(Hs.toarray())
                    E, V = E[:k], V[:, :k]
                else:
                    E, V = eigsh(Hs, k=k, which="SA", tol=1e-10)
                best.append(((nu, nd), sel, E, V))
        E0 = min(E.min() for _, _, E, _ in best)
        self.E0 = E0
        self.gs = []
        for sec, sel, E, V in best:
            if self.beta is None:
                idx = np.where(E - E0 < tol)[0]
                w = np.ones(len(idx))
            else:
                w_all = np.exp(-self.beta * (E - E0))
                idx = np.where(w_all > self.bcut)[0]
                w = w_all[idx]
            for i, wi in zip(idx, w):
                self.gs.append((sec, sel, V[:, i], float(E[i]), float(wi)))
        self.Zw = sum(g[4] for g in self.gs)

    def green(self, z, a=0, s=0):
        z = np.atleast_1d(np.asarray(z, complex))
        dim = 2**self.n_so
        cop = self.c[2 * a + s].tocsr()
        G = np.zeros(len(z), complex)
        for (nu, nd), sel, v, Ei, wi in self.gs:
            psi = np.zeros(dim)
            psi[sel] = v
            up = cop.T @ psi
            sp_sec = self._sector(nu + (s == 0), nd + (s == 1))
            if len(sp_sec):
                G += wi * _continued_fraction(self.H[sp_sec][:, sp_sec], up[sp_sec], z, 1.0, Ei, self.n_lanczos)
            dn = cop @ psi
            sh_sec = self._sector(nu - (s == 0), nd - (s == 1))
            if len(sh_sec) and nu - (s == 0) >= 0 and nd - (s == 1) >= 0:
                G += wi * _continued_fraction(self.H[sh_sec][:, sh_sec], dn[sh_sec], z, -1.0, Ei, self.n_lanczos)
        return G / self.Zw

    def hybridization(self, z, a):
        z = np.asarray(z, complex)
        return np.sum(self.Vb[a][None, :] ** 2 / (z[:, None] - self.eb[a][None, :]), axis=1)

    def occupation(self, a=0, s=0):
        dim = 2**self.n_so
        n = self.c[2 * a + s].T @ self.c[2 * a + s]
        tot = 0.0
        for (nu, nd), sel, v, Ei, wi in self.gs:
            psi = np.zeros(dim)
            psi[sel] = v
            tot += wi * (psi @ (n @ psi))
        return tot / self.Zw
