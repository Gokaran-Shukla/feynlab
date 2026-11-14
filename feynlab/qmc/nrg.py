"""Numerical renormalization group (Wilson NRG) for the single-impurity Anderson model:
thermodynamics (impurity entropy, susceptibility, Kondo temperature).

Wilson, RMP 47, 773 (1975); Krishna-murthy, Wilkins & Wilson, PRB 21, 1003 & 1044 (1980);
Bulla, Costi & Pruschke, RMP 80, 395 (2008).

Flat band [-D, D] (D = 1), hybridization Gamma = pi rho V^2.  Logarithmic discretisation
with parameter Lambda -> Wilson chain with hoppings
    t_n = (1 + 1/Lambda)/2 * Lambda^(-n/2) * xi_n ,
    xi_n = (1 - Lambda^-(n+1)) / sqrt((1 - Lambda^-(2n+1))(1 - Lambda^-(2n+3)))
and impurity-chain coupling V0 = sqrt(2 Gamma A_Lambda / pi), A_Lambda the standard
discretisation correction.  Iterative diagonalisation keeps N_keep states; at iteration N
thermodynamics is evaluated at T_N = omega_N / beta_bar from the kept spectrum.
Impurity quantities = (impurity + chain) - (chain only).
"""
from __future__ import annotations

import numpy as np


def _site_ops():
    # basis |0>, |up>, |dn>, |updn>  (up before dn in Jordan-Wigner order)
    au = np.zeros((4, 4)); ad = np.zeros((4, 4))
    au[0, 1] = 1; au[2, 3] = 1          # a_up |up> = |0>, a_up |up dn> = |dn>
    ad[0, 2] = 1; ad[1, 3] = -1         # a_dn |dn> = |0>, a_dn |up dn> = -|up>
    P = np.diag([1, -1, -1, 1.0])
    sz = np.diag([0, 0.5, -0.5, 0])
    n = au.T @ au + ad.T @ ad
    return au, ad, P, sz, n


def wilson_chain(Lam, n_sites):
    n = np.arange(n_sites)
    xi = (1 - Lam ** (-(n + 1))) / np.sqrt((1 - Lam ** (-(2 * n + 1))) * (1 - Lam ** (-(2 * n + 3))))
    return 0.5 * (1 + 1 / Lam) * Lam ** (-n / 2) * xi


def _A_Lambda(Lam):
    return 0.5 * (1 + 1 / Lam) / (1 - 1 / Lam) * np.log(Lam)


class NRG:
    def __init__(self, U=0.0, eps_d=None, Gamma=0.01, Lam=2.5, n_iter=40, n_keep=400, beta_bar=2.0,
                 impurity=True, A_correction=True):
        self.U, self.Gamma, self.Lam = U, Gamma, Lam
        self.eps_d = -U / 2 if eps_d is None else eps_d
        self.n_iter, self.n_keep, self.bb = n_iter, n_keep, beta_bar
        self.impurity = impurity
        self.V0 = np.sqrt(2 * Gamma * (_A_Lambda(Lam) if A_correction else 1.0) / np.pi)
        self.t = wilson_chain(Lam, n_iter + 2)

    def run(self):
        au, ad, P, sz, nsite = _site_ops()
        I4 = np.eye(4)
        if self.impurity:
            Hd = self.eps_d * nsite + self.U * np.diag([0, 0, 0, 1.0])
            # impurity (block) + f0 (site)
            H = np.kron(Hd, I4)
            for a in (au, ad):
                f_site = np.kron(P, a)                  # f0 in combined space
                d_blk = np.kron(a, I4)
                hop = self.V0 * d_blk.T @ f_site
                H = H + hop + hop.T
            Sz = np.kron(sz, I4) + np.kron(np.eye(4), sz)
            Par = np.kron(P, P)
            f_last = [np.kron(P, au), np.kron(P, ad)]
        else:
            H = np.zeros((4, 4))
            Sz, Par = sz.copy(), P.copy()
            f_last = [au.copy(), ad.copy()]
        out = []
        E, U_ = np.linalg.eigh(H)
        for N in range(self.n_iter):
            # rotate to eigenbasis and truncate
            keep = min(self.n_keep, len(E))
            E0 = E[0]
            E = E[:keep] - E0
            Uk = U_[:, :keep]
            Sz_e = Uk.T @ Sz @ Uk
            Sz2 = Uk.T @ (Sz @ Sz) @ Uk
            Par_e = Uk.T @ Par @ Uk
            f_e = [Uk.T @ f @ Uk for f in f_last]
            # thermodynamics at T_N
            omega_N = 0.5 * (1 + 1 / self.Lam) * self.Lam ** (-(N - 1) / 2 if self.impurity else -(N - 1) / 2)
            T = omega_N / self.bb
            w = np.exp(-(E - E.min()) / T)
            Z = w.sum()
            Eavg = np.sum(w * E) / Z
            S = (Eavg - E.min()) / T + np.log(Z)
            Tchi = float(np.sum(w * np.diag(Sz2)) / Z)
            out.append((T, S, Tchi))
            # add next chain site
            tN = self.t[N] if self.impurity else self.t[N]
            dim = keep
            Hn = np.kron(np.diag(E), I4)
            for s_, a in enumerate((au, ad)):
                f_site = np.kron(Par_e, a)
                hop = tN * np.kron(f_e[s_], I4).T @ f_site
                Hn = Hn + hop + hop.T
            Sz = np.kron(Sz_e, I4) + np.kron(np.eye(dim), sz)
            Par = np.kron(Par_e, P)
            f_last = [np.kron(Par_e, au), np.kron(Par_e, ad)]
            Hn = 0.5 * (Hn + Hn.T)
            E, U_ = np.linalg.eigh(Hn)
        return np.array(out)


def impurity_thermodynamics(U=0.0, Gamma=0.01, eps_d=None, Lam=2.5, n_iter=40, n_keep=400, beta_bar=2.0):
    """Returns arrays T, S_imp, T chi_imp (impurity minus bare chain, matched iteration by iteration).
    The bare chain has one site fewer at the same iteration index, so the even/odd structure of the
    two runs is aligned by averaging consecutive iterations (standard practice)."""
    a = NRG(U, eps_d, Gamma, Lam, n_iter, n_keep, beta_bar, impurity=True).run()
    b = NRG(U, eps_d, Gamma, Lam, n_iter, n_keep, beta_bar, impurity=False).run()
    # both runs contain chain sites 0..N at iteration N (same T_N)
    T = a[:, 0]
    S = a[:, 1] - b[:, 1]
    X = a[:, 2] - b[:, 2]
    # even/odd averaging
    S2 = 0.5 * (S[1:] + S[:-1]); X2 = 0.5 * (X[1:] + X[:-1]); T2 = np.sqrt(T[1:] * T[:-1])
    return T2, S2, X2


def resonant_level_exact(T, Gamma, eps_d=0.0):
    """Exact U=0 impurity entropy and T*chi of a wide-band resonant level (2 spins), from the
    impurity spectral function rho_d(e) = (Gamma/pi)/((e-eps_d)^2 + Gamma^2):
        S_imp = 2 int de rho_d(e) s(e/T) ... (wide-band, T << D) is NOT exact: the correct
    formula uses the phase-shift derivative, F_imp = -2T int de (rho_d(e)) ln(1 + e^{-e/T}),
    S = -dF/dT (evaluated by a centred difference), chi_imp = (1/2) int rho_d (-df/de)."""
    from scipy import integrate

    def F(t):
        g = lambda x: (Gamma / np.pi) / ((x - eps_d) ** 2 + Gamma**2) * (-t) * np.logaddexp(0.0, -x / t)
        pts = sorted({eps_d, 0.0})
        v = 0.0
        edges = [-np.inf] + pts + [np.inf]
        for lo, hi in zip(edges[:-1], edges[1:]):
            v += integrate.quad(g, lo, hi, limit=400, epsabs=1e-14, epsrel=1e-11)[0]
        return 2 * v

    S, X = [], []
    for t in np.atleast_1d(T):
        h = 1e-3 * t
        S.append(-(F(t + h) - F(t - h)) / (2 * h))
        chi = integrate.quad(lambda x: (Gamma / np.pi) / ((x - eps_d) ** 2 + Gamma**2) /
                             (4 * t * np.cosh(x / (2 * t)) ** 2), -np.inf, np.inf, limit=400)[0]
        X.append(0.5 * t * chi)
    return np.array(S), np.array(X)


def kondo_temperature(T, Tchi, threshold=0.07):
    """Wilson's convention T_K: T chi_imp(T_K) = 0.07 (log interpolation)."""
    T = np.asarray(T); X = np.asarray(Tchi)
    for i in range(len(T) - 1):
        if (X[i] - threshold) * (X[i + 1] - threshold) <= 0 and X[i] != X[i + 1]:
            f = (threshold - X[i]) / (X[i + 1] - X[i])
            return float(np.exp(np.log(T[i]) + f * (np.log(T[i + 1]) - np.log(T[i]))))
    return np.nan


def haldane_TK(U, Gamma):
    """Haldane, PRL 40, 416 (1978):  T_K ~ sqrt(U Gamma / 2) exp(-pi U/(8 Gamma) + pi Gamma/(2 U))."""
    return np.sqrt(U * Gamma / 2) * np.exp(-np.pi * U / (8 * Gamma) + np.pi * Gamma / (2 * U))
