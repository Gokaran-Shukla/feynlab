"""Multi-orbital impurities with the Slater-Kanamori interaction (Hund's metals,
t2g / eg shells, kagome d-bands).

H_int = U sum_a n_a,up n_a,dn + U' sum_{a!=b} n_a,up n_b,dn + (U'-J) sum_{a<b,s} n_a,s n_b,s
        - J sum_{a!=b} ( c+_a,up c_a,dn c+_b,dn c_b,up  +  c+_a,up c+_a,dn c_b,dn c_b,up )
with U' = U - 2J (rotational invariance).  Kanamori, Prog. Theor. Phys. 30, 275 (1963);
Georges, de' Medici, Mravlje, Annu. Rev. CMP 4, 137 (2013).

* KanamoriAtom      : exact multiplets of the isolated shell (+ crystal field)
* MultiOrbitalED    : impurity + bath sites, sector (N, Sz) diagonalisation, G_ab(iw)
* hubbard_one_sigma : multi-orbital Hubbard-I self-energy from the atomic G
"""
from __future__ import annotations

import itertools

import numpy as np


def _jw_ops(n_so, sparse=True):
    """Jordan-Wigner annihilation operators for n_so spin-orbitals (dim 2^n_so), sparse CSR."""
    import scipy.sparse as sp
    dim = 2**n_so
    ops = []
    states = np.arange(dim)
    for j in range(n_so):
        occ = (states >> j) & 1
        src = states[occ == 1]
        tgt = src ^ (1 << j)
        sgn = np.array([(-1) ** bin(int(s) & ((1 << j) - 1)).count("1") for s in src], float)
        c = sp.csr_matrix((sgn, (tgt, src)), shape=(dim, dim))
        ops.append(c if sparse else c.toarray())
    return ops


def kanamori_hamiltonian(c, n_orb, U, J, Up=None, idx=None):
    """H_int on spin-orbital operators c; idx(a, s) -> index (default 2a+s)."""
    Up = U - 2 * J if Up is None else Up
    idx = idx or (lambda a, s: 2 * a + s)
    n = lambda a, s: c[idx(a, s)].T @ c[idx(a, s)]
    cd = lambda a, s: c[idx(a, s)].T
    cc = lambda a, s: c[idx(a, s)]
    H = 0
    for a in range(n_orb):
        H = H + U * n(a, 0) @ n(a, 1)
    for a in range(n_orb):
        for b in range(n_orb):
            if a == b:
                continue
            H = H + Up * n(a, 0) @ n(b, 1)
            H = H - J * (cd(a, 0) @ cc(a, 1) @ cd(b, 1) @ cc(b, 0))          # spin flip
            H = H - J * (cd(a, 0) @ cd(a, 1) @ cc(b, 1) @ cc(b, 0))          # pair hopping
    for a in range(n_orb):
        for b in range(a + 1, n_orb):
            for s in (0, 1):
                H = H + (Up - J) * n(a, s) @ n(b, s)
    return H


class KanamoriAtom:
    def __init__(self, n_orb, U, J, eps=None, Up=None, mu=0.0):
        self.n_orb = n_orb
        self.c = _jw_ops(2 * n_orb)
        eps = np.zeros(n_orb) if eps is None else np.asarray(eps, float)
        H = kanamori_hamiltonian(self.c, n_orb, U, J, Up)
        for a in range(n_orb):
            for s in (0, 1):
                k = 2 * a + s
                H = H + (eps[a] - mu) * self.c[k].T @ self.c[k]
        self.c = [ci.toarray() for ci in self.c]
        self.H = H.toarray()
        self.N = sum(ci.T @ ci for ci in self.c)
        self.E, self.V = np.linalg.eigh(self.H)

    def multiplets(self, N):
        """Eigenvalues within the N-electron sector."""
        Nd = np.diag(self.N)
        sel = np.where(np.abs(Nd - N) < 1e-9)[0]
        return np.sort(np.linalg.eigvalsh(self.H[np.ix_(sel, sel)]))

    def green(self, z, beta):
        """Atomic G_aa(z) (diagonal, orbital a=0..n_orb-1, spin up)."""
        E = self.E - self.E.min()
        w = np.exp(-beta * E)
        Z = w.sum()
        out = []
        for a in range(self.n_orb):
            A = self.V.T @ self.c[2 * a] @ self.V
            num = np.abs(A) ** 2 * (w[:, None] + w[None, :])
            dE = E[None, :] - E[:, None]
            m = num > 1e-14
            out.append(np.array([np.sum(num[m] / (zz - dE[m])) for zz in np.atleast_1d(z)]) / Z)
        return np.array(out)

    def occupation(self, beta):
        E = self.E - self.E.min()
        w = np.exp(-beta * E)
        return float(np.sum(w * np.einsum("ik,ij,jk->k", self.V, self.N, self.V)) / w.sum())


def hubbard_one_sigma(z, atom: KanamoriAtom, beta, eps=None, mu=0.0):
    """Sigma_a(z) = z + mu - eps_a - G_at,a(z)^-1   (multi-orbital Hubbard-I)."""
    G = atom.green(z, beta)
    eps = np.zeros(atom.n_orb) if eps is None else np.asarray(eps)
    return np.array([np.asarray(z) + mu - eps[a] - 1.0 / G[a] for a in range(atom.n_orb)])


class MultiOrbitalED:
    """Impurity with n_orb Kanamori orbitals, each coupled to its own bath levels
    (diagonal hybridization).  Full diagonalisation per particle-number sector."""

    def __init__(self, n_orb, U, J, eps_d, eps_bath, V_bath, Up=None):
        self.n_orb = n_orb
        eb = [np.atleast_1d(e) for e in eps_bath]
        Vb = [np.atleast_1d(v) for v in V_bath]
        self.eb, self.Vb = eb, Vb
        nsites = n_orb + sum(len(e) for e in eb)
        n_so = 2 * nsites
        if n_so > 16:
            raise ValueError("MultiOrbitalED: at most 8 sites (16 spin-orbitals)")
        c = _jw_ops(n_so)
        self.c = c
        H = kanamori_hamiltonian(c, n_orb, U, J, Up)
        site = n_orb
        for a in range(n_orb):
            for s in (0, 1):
                H = H + eps_d[a] * c[2 * a + s].T @ c[2 * a + s]
            for e, V in zip(eb[a], Vb[a]):
                for s in (0, 1):
                    k = 2 * site + s
                    H = H + e * c[k].T @ c[k]
                    hop = V * c[2 * a + s].T @ c[k]
                    H = H + hop + hop.T
                site += 1
        H = H.tocsr()
        dim = H.shape[0]
        Nd = np.array([bin(i).count("1") for i in range(dim)])
        self.E = np.zeros(dim)
        self.blocks = []
        for N in np.unique(Nd):
            sel = np.where(Nd == N)[0]
            e, v = np.linalg.eigh(H[sel][:, sel].toarray())
            self.E[sel] = e
            self.blocks.append((sel, v))
        self.H = H

    def expectation(self, op, beta):
        """Thermal <op> for a Fock-space operator (sparse or dense) on the impurity+bath system."""
        E0 = self.E.min()
        tot, Z = 0.0, 0.0
        op = op.tocsr() if hasattr(op, "tocsr") else op
        for sel, v in self.blocks:
            wts = np.exp(-beta * (self.E[sel] - E0))
            if wts.max() < 1e-14:
                continue
            sub = op[sel][:, sel]
            sub = sub.toarray() if hasattr(sub, "toarray") else sub
            tot += np.sum(wts * np.einsum("ia,ij,ja->a", v, sub, v))
            Z += wts.sum()
        return float(tot / Z)

    def n_op(self, a, s):
        return self.c[2 * a + s].T @ self.c[2 * a + s]

    def hybridization(self, z, a):
        z = np.asarray(z, complex)
        return np.sum(self.Vb[a][None, :] ** 2 / (z[:, None] - self.eb[a][None, :]), axis=1)

    def green(self, z, beta, a=0, cutoff=1e-12):
        """Lehmann sum over sector pairs (N -> N-1) with Boltzmann cutoff."""
        z = np.atleast_1d(np.asarray(z, complex))
        E0 = self.E.min()
        Zp = np.exp(-beta * (self.E - E0)).sum()
        c = self.c[2 * a].tocsr()
        G = np.zeros(len(z), complex)
        byN = {len(sel) and bin(int(sel[0])).count("1"): (sel, v) for sel, v in self.blocks}
        for N, (sel, v) in byN.items():
            if N - 1 not in byN:
                continue
            sel2, v2 = byN[N - 1]
            En = self.E[sel] - E0
            Em = self.E[sel2] - E0
            wn, wm = np.exp(-beta * En), np.exp(-beta * Em)
            if wn.max() < cutoff and wm.max() < cutoff:
                continue
            A = v2.T @ (c[sel2][:, sel].toarray() @ v)        # <m, N-1| c |n, N>
            W = np.abs(A) ** 2 * (wn[None, :] + wm[:, None])
            msk = W > cutoff * 1e-3
            dE = (En[None, :] - Em[:, None])[msk]
            G += (W[msk][None, :] / (z[:, None] - dE[None, :])).sum(1)
        return G / Zp
