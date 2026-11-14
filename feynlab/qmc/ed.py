"""Exact diagonalization of a single-impurity Anderson model with a few bath sites.

H = sum_s eps_d n_ds + U n_du n_dd + sum_{l,s} eps_l n_ls + sum_{l,s} V_l (d_s^+ c_ls + h.c.)
(energies measured from the chemical potential).  Used as an exact reference for
CT-INT, IPT and second-order perturbation theory: the corresponding hybridization
function is Delta(iw) = sum_l V_l^2 / (iw - eps_l).
"""
from __future__ import annotations

import numpy as np


def _build_ops(n_orb):
    """Jordan-Wigner annihilation operators (dense) for n_orb spin-orbitals."""
    dim = 2**n_orb
    ops = []
    states = np.arange(dim)
    for j in range(n_orb):
        c = np.zeros((dim, dim))
        for s in states:
            if (s >> j) & 1:
                t = s ^ (1 << j)
                # sign = (-1)^{number of occupied orbitals with index < j}
                sgn = (-1) ** bin(s & ((1 << j) - 1)).count("1")
                c[t, s] = sgn
        ops.append(c)
    return ops


class AndersonED:
    def __init__(self, eps_d: float, U: float, eps_bath, V_bath):
        self.eps_d = float(eps_d)
        self.U = float(U)
        self.eps_bath = np.asarray(eps_bath, float)
        self.V_bath = np.asarray(V_bath, float)
        nb = len(self.eps_bath)
        self.n_orb = 2 * (1 + nb)  # orbital index: 2*site + spin, site 0 = impurity
        c = _build_ops(self.n_orb)
        self.c = c
        idx = lambda site, s: 2 * site + s
        dim = 2**self.n_orb
        H = np.zeros((dim, dim))
        n = [ci.T @ ci for ci in c]
        for s in (0, 1):
            H += self.eps_d * n[idx(0, s)]
            for l in range(nb):
                H += self.eps_bath[l] * n[idx(l + 1, s)]
                hop = self.V_bath[l] * c[idx(0, s)].T @ c[idx(l + 1, s)]
                H += hop + hop.T
        H += self.U * n[idx(0, 0)] @ n[idx(0, 1)]
        self.H = H
        self.E, self.vec = np.linalg.eigh(H)
        self.d_up = self.vec.T @ c[idx(0, 0)] @ self.vec
        self.n_imp = [self.vec.T @ n[idx(0, s)] @ self.vec for s in (0, 1)]
        self.docc_op = self.vec.T @ (n[idx(0, 0)] @ n[idx(0, 1)]) @ self.vec

    def hybridization(self, z):
        z = np.asarray(z, complex)
        return np.sum(self.V_bath[None, :] ** 2 / (z[:, None] - self.eps_bath[None, :]), axis=1)

    def green(self, z, beta: float):
        """Lehmann representation G(z) = (1/Z) sum_{mn} |<m|d^+|n>|^2 (e^{-bE_n}+e^{-bE_m})/(z + E_n - E_m)."""
        z = np.asarray(z, complex)
        E = self.E - self.E.min()
        w = np.exp(-beta * E)
        Z = w.sum()
        A = self.d_up  # <m|d|n>
        num = np.abs(A) ** 2 * (w[:, None] + w[None, :])  # indices (m, n) for <m|d|n>
        dE = E[None, :] - E[:, None]  # E_n - E_m
        mask = num > 1e-14
        nums = num[mask]
        dEs = dE[mask]
        return np.array([np.sum(nums / (zz - dEs)) for zz in z]) / Z

    def thermal(self, op_diag_basis, beta):
        E = self.E - self.E.min()
        w = np.exp(-beta * E)
        return float(np.sum(np.diag(op_diag_basis) * w) / w.sum())

    def occupation(self, beta):
        return self.thermal(self.n_imp[0], beta)

    def double_occupancy(self, beta):
        return self.thermal(self.docc_op, beta)

    def sigma(self, iw, beta):
        """Sigma = G0^{-1} - G^{-1} with G0^{-1} = iw - eps_d - Delta(iw)."""
        g = self.green(iw, beta)
        return iw - self.eps_d - self.hybridization(iw) - 1.0 / g
