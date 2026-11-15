"""Two-particle vertex functions, DMFT susceptibilities and ladder dynamical vertex approximation (DGA).

Impurity two-particle Green function from exact diagonalization (Lehmann representation)
-----------------------------------------------------------------------------------------
    G2_{ss'}(nu, nu', w) = int_0^beta dt1 dt2 dt3 e^{-i nu t1} e^{i(nu+w) t2} e^{-i(nu'+w) t3}
                           < T c+_s(t1) c_s(t2) c+_s'(t3) c_s'(0) >
    chi_{ss'}(nu, nu', w) = G2 - beta G_s(nu) G_s'(nu') delta_{w,0}      (Rohringer et al. PRB 86, 125114 (2012))
For every time ordering t_a > t_b > t_c the integral over the simplex of exp(sum x t) is the divided difference
of e^{beta z} at the nodes (-E_i, -E_j + i phi_a, -E_k + i(phi_a + phi_b), -E_l + i(phi_a + phi_b + phi_c)),
evaluated stably (also for degenerate nodes) as the corner element of exp(beta M), M bidiagonal (Opitz 1964).
Chains whose four energies all exceed the Boltzmann cutoff are pruned.

Channels and vertices (matrices in nu, nu' at fixed w):
    chi_d = chi_uu + chi_ud,  chi_m = chi_uu - chi_ud,  chi0 = -beta G(nu) G(nu+w) delta_{nu nu'}
    full vertex       chi_r = chi0 - chi0 F_r chi0 / beta^2
    irreducible       chi_r^{-1} = chi0^{-1} + Gamma_r / beta^2
    F_ud = (F_d - F_m)/2.   Equation of motion (Schwinger-Dyson):
    Sigma(nu) = U n/2 - (U/beta^2) sum_{nu' w} F_ud(nu, nu', w) G(nu') G(nu'+w) G(nu+w)

Lattice (DMFT / ladder DGA, Toschi, Katanin, Held PRB 75, 045118 (2007); Katanin et al. PRB 80, 075104 (2009);
Rohringer et al. RMP 90, 025003 (2018)):
    chi_r(q)^{-1} = chi0(q)^{-1} + Gamma_r / beta^2,   chi0(q)_{nu} = -(beta/N) sum_k G(k,nu) G(k+q,nu+w)
    physical susceptibility chi_r(q, w) = (1/beta^2) sum_{nu nu'} chi_r(q)_{nu nu'}   (+ asymptotic correction)
    ladder DGA self-energy (both ph and transverse ph-bar ladders, double counting removed):
    Sigma(k,nu) = Sigma_DMFT(nu) - (U/beta^2) sum_{nu' w q} [ (F_d(q)/2 - 3 F_m(q)/2) - (F_d/2 - 3 F_m/2)_loc ]
                  x  G(k',nu') G(k'+q, nu'+w) G(k+q, nu+w)  summed over k'.
"""
from __future__ import annotations

import itertools

import numpy as np
from numba import njit

from .qmc.multiorbital import _jw_ops


# =============================================================================== divided differences
@njit(cache=True)
def _dd4w(y0, y1, y2, y3, beta, M, E, T, P):
    """Divided difference of exp(beta z) at four complex nodes = (exp(beta J))[0, 3], J bidiagonal
    (Opitz 1964); scaling-and-squaring + Taylor, stable for (near-)degenerate nodes. M, E, T, P: 4x4 work arrays."""
    for i in range(4):
        for j in range(4):
            M[i, j] = 0.0
    M[0, 0] = beta * y0; M[1, 1] = beta * y1; M[2, 2] = beta * y2; M[3, 3] = beta * y3
    M[0, 1] = beta; M[1, 2] = beta; M[2, 3] = beta
    nrm = 0.0
    for i in range(4):
        r = 0.0
        for j in range(i, 4):
            r += abs(M[i, j])
        if r > nrm:
            nrm = r
    s = 0
    while nrm > 0.5:
        nrm *= 0.5
        s += 1
    sc = 0.5**s
    for i in range(4):
        for j in range(i, 4):
            M[i, j] *= sc
    for i in range(4):
        for j in range(4):
            E[i, j] = 1.0 if i == j else 0.0
            T[i, j] = E[i, j]
    for n in range(1, 18):
        # T = T @ M / n  (upper triangular)
        for i in range(4):
            for j in range(3, i - 1, -1):
                acc = 0.0 + 0.0j
                for l in range(i, j + 1):
                    acc += T[i, l] * M[l, j]
                P[i, j] = acc / n
        for i in range(4):
            for j in range(i, 4):
                T[i, j] = P[i, j]
                E[i, j] += P[i, j]
    for _ in range(s):
        for i in range(4):
            for j in range(i, 4):
                acc = 0.0 + 0.0j
                for l in range(i, j + 1):
                    acc += E[i, l] * E[l, j]
                P[i, j] = acc
        for i in range(4):
            for j in range(i, 4):
                E[i, j] = P[i, j]
    return E[0, 3]


@njit(cache=True)
def _dd4(y0, y1, y2, y3, beta):
    M = np.zeros((4, 4), dtype=np.complex128); E = np.zeros((4, 4), dtype=np.complex128)
    T = np.zeros((4, 4), dtype=np.complex128); P = np.zeros((4, 4), dtype=np.complex128)
    return _dd4w(y0, y1, y2, y3, beta, M, E, T, P)


@njit(cache=True)
def _g2_kernel(chains, amps, E, phis, beta, nu, nup, om, out_shape):
    """chains: (nc, 4) state indices (i,j,k,l) for one ordering, amps: products of matrix elements,
    phis: (3, nnu, nnup, nom) phases of the ordered operators (a, b, c).  Returns sum over chains."""
    nnu, nnp, nom = out_shape
    res = np.zeros((nnu, nnp, nom), dtype=np.complex128)
    M = np.zeros((4, 4), dtype=np.complex128); Ew = np.zeros((4, 4), dtype=np.complex128)
    T = np.zeros((4, 4), dtype=np.complex128); P = np.zeros((4, 4), dtype=np.complex128)
    for c in range(chains.shape[0]):
        i, j, k, l = chains[c, 0], chains[c, 1], chains[c, 2], chains[c, 3]
        a = amps[c]
        for x in range(nnu):
            for y in range(nnp):
                for z in range(nom):
                    pa = phis[0, x, y, z]
                    pb = phis[1, x, y, z]
                    pc = phis[2, x, y, z]
                    res[x, y, z] += a * _dd4w(-E[i] + 0j, -E[j] + 1j * pa, -E[k] + 1j * (pa + pb),
                                              -E[l] + 1j * (pa + pb + pc), beta, M, Ew, T, P)
    return res


# =============================================================================== impurity ED
class ImpurityED:
    """Single-orbital Anderson impurity + discrete bath (spin-resolved bath allowed, magnetic field h).
    H = sum_s (eps_d - s h) n_s + U n_up n_dn + sum_{l s} e_ls b+ b + V_ls (c+ b + h.c.)."""

    def __init__(self, eps_d, U, eb, V, eb_dn=None, V_dn=None, h=0.0):
        eb, V = np.atleast_1d(np.asarray(eb, float)), np.atleast_1d(np.asarray(V, float))
        ebd = eb if eb_dn is None else np.atleast_1d(np.asarray(eb_dn, float))
        Vd = V if V_dn is None else np.atleast_1d(np.asarray(V_dn, float))
        nb = len(eb)
        self.nb = nb
        self.U = float(U)
        n_so = 2 * (1 + nb)
        c = [x.toarray() for x in _jw_ops(n_so)]
        self.c = c
        H = np.zeros((2**n_so, 2**n_so))
        for s, sg in ((0, 1.0), (1, -1.0)):
            H += (eps_d - sg * h) * c[s].T @ c[s]
            e_, v_ = (eb, V) if s == 0 else (ebd, Vd)
            for l in range(nb):
                b = c[2 * (l + 1) + s]
                H += e_[l] * b.T @ b + v_[l] * (c[s].T @ b + b.T @ c[s])
        H += U * (c[0].T @ c[0]) @ (c[1].T @ c[1])
        self.H = H
        E, Vv = np.linalg.eigh(H)
        self.E0 = E.min()
        self.E = E - self.E0
        self.Vv = Vv
        self.eb, self.V, self.ebd, self.Vd = eb, V, ebd, Vd
        self.eps_d, self.h = eps_d, h

    def _op(self, s, dag=False):
        o = self.c[s].T if dag else self.c[s]
        return self.Vv.T @ o @ self.Vv

    def Z(self, beta):
        return np.exp(-beta * self.E).sum()

    def expect(self, O, beta):
        Ob = self.Vv.T @ O @ self.Vv
        w = np.exp(-beta * self.E)
        return float(np.sum(w * np.diag(Ob)) / w.sum())

    def occupation(self, beta, s=0):
        return self.expect(self.c[s].T @ self.c[s], beta)

    def double_occupancy(self, beta):
        return self.expect((self.c[0].T @ self.c[0]) @ (self.c[1].T @ self.c[1]), beta)

    def hybridization(self, z, s=0):
        e_, v_ = (self.eb, self.V) if s == 0 else (self.ebd, self.Vd)
        return np.sum(v_[None, :] ** 2 / (np.asarray(z)[:, None] - e_[None, :]), axis=1)

    def green(self, z, beta, s=0):
        z = np.atleast_1d(np.asarray(z, complex))
        A = self._op(s)
        w = np.exp(-beta * self.E)
        W = np.abs(A) ** 2 * (w[None, :] + w[:, None])     # A[m, n] = <m|c|n>
        dE = self.E[None, :] - self.E[:, None]
        m = W > 1e-14
        return (W[m][None, :] / (z[:, None] - dE[m][None, :])).sum(1) / w.sum()

    def sigma(self, z, beta, s=0):
        z = np.atleast_1d(np.asarray(z, complex))
        sg = 1.0 if s == 0 else -1.0
        return z - (self.eps_d - sg * self.h) - self.hybridization(z, s) - 1.0 / self.green(z, beta, s)

    def static_susceptibility(self, beta, which="m"):
        """chi = int_0^beta <O(tau) O> - beta <O>^2 with O = (n_up - n_dn) (m) or (n_up + n_dn) (d); Lehmann."""
        nu, nd = self.c[0].T @ self.c[0], self.c[1].T @ self.c[1]
        O = (nu - nd) if which == "m" else (nu + nd)
        Ob = self.Vv.T @ O @ self.Vv
        E = self.E
        w = np.exp(-beta * E)
        Zs = w.sum()
        dE = E[None, :] - E[:, None]                                # E_n - E_m  [m, n]
        with np.errstate(divide="ignore", invalid="ignore"):
            K = np.where(np.abs(dE) > 1e-10, (w[:, None] - w[None, :]) / np.where(np.abs(dE) > 1e-10, dE, 1.0),
                         beta * w[:, None])
        val = np.sum(np.abs(Ob) ** 2 * K) / Zs
        mean = np.sum(w * np.diag(Ob)) / Zs
        return float(val - beta * mean**2)

    # ------------------------------------------------------------------ two-particle
    def g2(self, beta, n_nu=8, n_om=0, s1=0, s2=0, cutoff=30.0, tol=1e-12):
        """G2_{s1 s2}(nu, nu', w) on nu, nu' = (2n+1) pi/beta for n in [-n_nu, n_nu), w = 2 m pi/beta, |m| <= n_om.
        Returns (G2[w, nu, nu'], nu_values, om_values)."""
        A = [self._op(s1, True), self._op(s1), self._op(s2, True), self._op(s2)]
        nus = (2 * np.arange(-n_nu, n_nu) + 1) * np.pi / beta
        oms = 2 * np.arange(-n_om, n_om + 1) * np.pi / beta
        NU, NUP, OM = np.meshgrid(nus, nus, oms, indexing="ij")
        ph = np.array([-NU, NU + OM, -(NUP + OM)])               # phases of operators 1, 2, 3
        E = self.E
        low = E < cutoff / beta
        res = np.zeros(NU.shape, complex)
        for perm in itertools.permutations(range(3)):
            sign = -1.0 if sum(1 for x in range(3) for y in range(x + 1, 3) if perm[x] > perm[y]) % 2 else 1.0
            Aa, Ab, Ac, A4 = A[perm[0]], A[perm[1]], A[perm[2]], A[3]
            chains, amps = [], []
            nz_a = [np.nonzero(np.abs(Aa[i]) > tol)[0] for i in range(len(E))]
            nz_b = [np.nonzero(np.abs(Ab[j]) > tol)[0] for j in range(len(E))]
            nz_c = [np.nonzero(np.abs(Ac[k]) > tol)[0] for k in range(len(E))]
            for i in range(len(E)):
                for j in nz_a[i]:
                    for k in nz_b[j]:
                        for l in nz_c[k]:
                            amp = Aa[i, j] * Ab[j, k] * Ac[k, l] * A4[l, i]
                            if abs(amp) < tol:
                                continue
                            if not (low[i] or low[j] or low[k] or low[l]):
                                continue
                            chains.append((i, j, k, l)); amps.append(amp)
            if not chains:
                continue
            phis = np.ascontiguousarray(ph[list(perm)])
            r = _g2_kernel(np.array(chains, dtype=np.int64), np.array(amps, dtype=np.complex128), E.astype(float),
                           phis, float(beta), nus, nus, oms, NU.shape)
            res += sign * r
        res /= self.Z(beta)
        return np.transpose(res, (2, 0, 1)), nus, oms

    def generalized_susceptibilities(self, beta, n_nu=8, n_om=0, **kw):
        """chi_uu, chi_ud [w, nu, nu'] and the bubble chi0 [w, nu] (paramagnetic)."""
        Guu, nus, oms = self.g2(beta, n_nu, n_om, 0, 0, **kw)
        Gud, _, _ = self.g2(beta, n_nu, n_om, 0, 1, **kw)
        nmax = n_nu + n_om + 1
        allnu = (2 * np.arange(-nmax, nmax) + 1) * np.pi / beta
        Gall = self.green(1j * allnu, beta, 0)
        Gf = lambda x: Gall[np.searchsorted(allnu, x - 1e-9)]
        Gnu = Gf(nus)
        chi_uu, chi_ud = Guu.copy(), Gud.copy()
        i0 = int(np.argmin(np.abs(oms)))
        chi_uu[i0] -= beta * np.outer(Gnu, Gnu)
        chi_ud[i0] -= beta * np.outer(Gnu, Gnu)
        chi0 = np.array([-beta * Gnu * Gf(nus + w) for w in oms])
        return dict(chi_uu=chi_uu, chi_ud=chi_ud, chi0=chi0, nu=nus, om=oms, G=Gnu, G_of=Gf)


# =============================================================================== vertex algebra
def vertices(gs, beta):
    """F_d, F_m, Gamma_d, Gamma_m [w, nu, nu'] from generalized susceptibilities (local)."""
    chi_d = gs["chi_uu"] + gs["chi_ud"]
    chi_m = gs["chi_uu"] - gs["chi_ud"]
    out = {}
    for name, chi in (("d", chi_d), ("m", chi_m)):
        F, Gm = [], []
        for w in range(chi.shape[0]):
            c0 = gs["chi0"][w]
            c0i = np.diag(1.0 / c0)
            F.append(beta**2 * c0i @ (np.diag(c0) - chi[w]) @ c0i)
            Gm.append(beta**2 * (np.linalg.inv(chi[w]) - c0i))
        out["F_" + name] = np.array(F)
        out["Gamma_" + name] = np.array(Gm)
        out["chi_" + name] = chi
    out["F_ud"] = 0.5 * (out["F_d"] - out["F_m"])
    return out


def eom_sigma(F_ud, gs, beta, U, n_total):
    """Local Schwinger-Dyson equation Sigma(nu) = U n/2 - (U/beta^2) sum F_ud G(nu') G(nu'+w) G(nu+w)."""
    nus, oms, Gf = gs["nu"], gs["om"], gs["G_of"]
    S = np.full(len(nus), U * n_total / 2, complex)
    for iw, w in enumerate(oms):
        Gp = Gf(nus) * Gf(nus + w)                 # G(nu') G(nu'+w)
        S -= (U / beta**2) * (F_ud[iw] @ Gp) * Gf(nus + w)
    return S


def physical_susceptibility(chi_r, beta):
    """chi_r(w) = (1/beta^2) sum_{nu nu'} chi_r(nu, nu', w)."""
    return np.real(chi_r.sum(axis=(-1, -2))) / beta**2


# =============================================================================== lattice ladders
def lattice_bubble(G_k_of, nus, oms, beta, q_shift_index):
    """chi0(q)_{w, nu} = -(beta/N) sum_k G(k, nu) G(k+q, nu+w).
    G_k_of(nu_array) -> array (n_k1, n_k2, len(nu)) on a periodic k mesh; q_shift_index = integer mesh shift."""
    out = []
    Gk = G_k_of(nus)
    for w in oms:
        Gkq = np.roll(G_k_of(nus + w), shift=[-q_shift_index[0], -q_shift_index[1]], axis=(0, 1))
        out.append(-beta * np.mean(Gk * Gkq, axis=(0, 1)))
    return np.array(out)


def lattice_susceptibility(Gamma_r, chi0_q, beta):
    """chi_r(q)_{nu nu'} = [chi0(q)^{-1} + Gamma_r / beta^2]^{-1}  per bosonic frequency; returns matrices."""
    return np.array([np.linalg.inv(np.diag(1.0 / chi0_q[w]) + Gamma_r[w] / beta**2) for w in range(len(chi0_q))])


def ladder_F(Gamma_r, chi0_q, beta):
    """F_r(q) from chi_r(q) = chi0(q) - chi0(q) F_r(q) chi0(q)/beta^2."""
    chi = lattice_susceptibility(Gamma_r, chi0_q, beta)
    out = []
    for w in range(len(chi0_q)):
        c0i = np.diag(1.0 / chi0_q[w])
        out.append(beta**2 * c0i @ (np.diag(chi0_q[w]) - chi[w]) @ c0i)
    return np.array(out)


def ladder_dga_sigma(eps_k, mu, sigma_dmft_of, vert, gs, beta, U, n_total, q_stride=1, G_loc_of=None):
    """Ladder DGA self-energy (Katanin, Toschi, Held PRB 80, 075104 (2009); Rohringer et al. RMP 90, 025003):
        Sigma_EOM(k,nu) = U n/2 - (U/beta^2) sum_{q w nu'} K(q)_{nu nu'} [sum_k' G(k',nu') G(k'+q,nu'+w)] G(k+q,nu+w)
        K(q) = F_d(q)/2 - 3 F_m(q)/2 - F_d^loc/2 + F_m^loc/2      (ph + transverse ph-bar ladders, the doubly counted
                                                                    local ph-bar part removed; K -> U at lowest order)
        Sigma_DGA(k,nu) = Sigma_EOM(k,nu) - Sigma_EOM^loc(nu) + Sigma_DMFT(nu)
    where Sigma_EOM^loc is the same expression with local propagators and K_loc = -F_m^loc (this cancels the
    frequency-truncation error of the local part).  To O(U^2) Sigma_DGA equals the lattice second-order self-energy.
    eps_k: (N, N) mesh; returns Sigma[kx, ky, nu] on nu = gs['nu']."""
    nus, oms = gs["nu"], gs["om"]
    N = eps_k.shape[0]

    def Gk_of(x):
        return 1.0 / (1j * x[None, None, :] + mu - eps_k[..., None] - sigma_dmft_of(x)[None, None, :])
    Floc_d, Floc_m = vert["F_d"], vert["F_m"]
    Gk_shift = {w: Gk_of(nus + w) for w in oms}
    S = np.zeros((N, N, len(nus)), complex)
    qs = [(a, b) for a in range(0, N, q_stride) for b in range(0, N, q_stride)]
    wq = 1.0 / len(qs)
    for (a, b) in qs:
        c0q = lattice_bubble(Gk_of, nus, oms, beta, (a, b))
        Fd = ladder_F(vert["Gamma_d"], c0q, beta)
        Fm = ladder_F(vert["Gamma_m"], c0q, beta)
        Kq = 0.5 * Fd - 1.5 * Fm - 0.5 * Floc_d + 0.5 * Floc_m
        for iw, w in enumerate(oms):
            gg = -c0q[iw] / beta                                   # sum_k' G(k',nu') G(k'+q,nu'+w)
            vec = Kq[iw] @ gg
            Gkq = np.roll(Gk_shift[w], shift=[-a, -b], axis=(0, 1))
            S -= (U / beta**2) * wq * vec[None, None, :] * Gkq
    # local counterpart with impurity propagators
    Gf = gs["G_of"] if G_loc_of is None else G_loc_of
    Sloc = np.zeros(len(nus), complex)
    for iw, w in enumerate(oms):
        gg = Gf(nus) * Gf(nus + w)
        Sloc -= (U / beta**2) * ((-Floc_m[iw]) @ gg) * Gf(nus + w)
    return S - Sloc[None, None, :] + sigma_dmft_of(nus)[None, None, :]


# =============================================================================== TPSC
def tpsc(eps_k, beta, U, n=1.0, n_iw=512, tol=1e-10):
    """Two-particle self-consistent approach (Vilk & Tremblay, J. Phys. I France 7, 1309 (1997)).
    Spin and charge vertices from the local sum rules
        (T/N) sum_q chi_sp(q) = n - 2 <n_up n_dn>,   (T/N) sum_q chi_ch(q) = n + 2 <n_up n_dn> - n^2,
    with the ansatz U_sp = U <n_up n_dn> / (<n_up><n_dn>) and RPA forms chi = chi0/(1 -+ U chi0/2).
    Self-energy: Sigma(k) = U n/2 + (U/8)(T/N) sum_q [3 U_sp chi_sp(q) + U_ch chi_ch(q)] G0(k+q).
    Returns dict(U_sp, U_ch, docc, chi_sp, chi_ch, Sigma, satisfies sum rules)."""
    from .twoparticle import bubble
    from .core.matsubara import MatsubaraGrid
    g = MatsubaraGrid(beta, n_iw)
    N = eps_k.shape[0]
    # chemical potential of the non-interacting system at filling n
    from .core.matsubara import fermi
    lo, hi = eps_k.min() - 5, eps_k.max() + 5
    for _ in range(100):
        mu = 0.5 * (lo + hi)
        nn = 2 * np.mean(fermi(eps_k - mu, beta))
        lo, hi = (mu, hi) if nn < n else (lo, mu)
    chi0 = 2 * bubble(eps_k, np.zeros(n_iw), g, mu)        # convention: chi0 includes spin factor 2 (TPSC)
    chi0 = np.real(chi0)
    nb = chi0.shape[-1]

    from scipy.special import polygamma
    M = nb - 1
    Om_M = 2 * np.pi * M / beta
    psi1 = float(polygamma(1, M + 1))            # sum_{m > M} 1/m^2

    def local_sum(chi):
        # (T/N) sum_{q, Omega} with Omega >= 0 stored: chi(-Omega) = chi(Omega); the 1/Omega^2 tail beyond
        # the last frequency is added analytically:  sum_{m>M} C/Omega_m^2 = C (beta/2pi)^2 psi'(M+1)
        loc = chi.mean(axis=(0, 1))
        C = float(loc[-1]) * Om_M**2
        s = loc[0] + 2 * loc[1:].sum() + 2 * C * (beta / (2 * np.pi)) ** 2 * psi1
        return s / beta

    def chi_sp(Usp):
        return chi0 / (1 - 0.5 * Usp * chi0)

    def chi_ch(Uch):
        return chi0 / (1 + 0.5 * Uch * chi0)

    Umax = 2.0 / chi0[..., 0].max()
    lo, hi = 0.0, Umax * (1 - 1e-12)
    for _ in range(200):                         # self-consistency U_sp: sum rule with docc = U_sp n^2/(4U)
        Usp = 0.5 * (lo + hi)
        docc = Usp * (n / 2) ** 2 / U
        f = local_sum(chi_sp(Usp)) - (n - 2 * docc)
        lo, hi = (Usp, hi) if f < 0 else (lo, Usp)
        if hi - lo < tol:
            break
    Usp = 0.5 * (lo + hi)
    docc = Usp * (n / 2) ** 2 / U
    target = n + 2 * docc - n**2
    lo, hi = 0.0, 1e3
    for _ in range(200):
        Uch = 0.5 * (lo + hi)
        f = local_sum(chi_ch(Uch)) - target
        lo, hi = (lo, Uch) if f < 0 else (Uch, hi)
        if hi - lo < tol:
            break
    Uch = 0.5 * (lo + hi)
    return dict(U_sp=Usp, U_ch=Uch, docc=docc, chi_sp=chi_sp(Usp), chi_ch=chi_ch(Uch), chi0=chi0, mu=mu,
                sum_sp=local_sum(chi_sp(Usp)), sum_ch=local_sum(chi_ch(Uch)), target_sp=n - 2 * docc, target_ch=target)


# =============================================================================== ED-DMFT driver for vertex work
def ed_dmft(beta, U, n_bath=2, D=1.0, eps_k=None, h=0.0, n_iw=256, n_iter=60, mix=0.5, tol=1e-7, mu=None,
            beta_fit=None, verbose=False):
    """Finite-temperature ED-DMFT returning the converged ImpurityED (spin-resolved bath, field h).
    Bethe lattice of half bandwidth D (Delta = D^2/4 G) or a lattice given by eps_k (any shape).
    Used to supply the local vertex of DMFT consistently with the impurity model it came from."""
    from .core.matsubara import MatsubaraGrid
    from .qmc.ed_sector import fit_bath
    g = MatsubaraGrid(beta, n_iw)
    iw = g.iw
    mu = U / 2 if mu is None else mu
    eps_d = -mu
    Ds = [D**2 / 4 * (1.0 / (iw + 0.3)), D**2 / 4 * (1.0 / (iw + 0.3))]
    ed = None
    hist = []
    nf = min(n_iw, 120)
    for it in range(n_iter):
        baths = []
        for s in (0, 1):
            e_, v_, *_ = fit_bath(iw[:nf], Ds[s][:nf], n_bath, symmetric=(h == 0.0))
            baths.append((e_, v_))
        ed = ImpurityED(eps_d, U, baths[0][0], baths[0][1], baths[1][0], baths[1][1], h=h)
        newD = []
        for s in (0, 1):
            Gs = ed.green(iw, beta, s)
            if eps_k is None:
                newD.append(D**2 / 4 * Gs)
            else:
                Sg = ed.sigma(iw, beta, s)
                sg = 1.0 if s == 0 else -1.0
                Gl = np.mean(1.0 / (iw[None, :] + mu + sg * h - np.ravel(eps_k)[:, None] - Sg[None, :]), axis=0)
                newD.append(iw + mu + sg * h - Sg - 1.0 / Gl)
        diff = max(float(np.max(np.abs(newD[s] - Ds[s]))) for s in (0, 1))
        Ds = [mix * newD[s] + (1 - mix) * Ds[s] for s in (0, 1)]
        hist.append(diff)
        if verbose:
            print(it, diff)
        if diff < tol:
            break
    ed.history = hist
    ed.mu = mu
    return ed


def tpsc_sigma(res, eps_k, beta, U, n=1.0, n_nu=None):
    """TPSC self-energy (Vilk & Tremblay, J. Phys. I 7, 1309 (1997), crossing-symmetric form)
        Sigma(k, nu) = U n/2 + (U/8) (T/N) sum_q [3 U_sp chi_sp(q) + U_ch chi_ch(q)] G0(k+q, nu+Omega)
    evaluated by FFT convolution in k for every bosonic frequency.  Returns (nu, Sigma[kx, ky, nu]).
    Exact internal check:  (T/N) sum_k Sigma(k) G0(k) e^{i nu 0+} = U <n_up n_dn>   (tpsc_trace_check)."""
    chi_sp, chi_ch = res["chi_sp"], res["chi_ch"]
    nb = chi_sp.shape[-1]
    n_nu = nb if n_nu is None else n_nu
    mu = res["mu"]
    V = 3 * res["U_sp"] * chi_sp + res["U_ch"] * chi_ch            # (N, N, nb), Omega_m >= 0, even in Omega
    nus = (2 * np.arange(-n_nu, n_nu) + 1) * np.pi / beta
    Vf = np.fft.fft2(V, axes=(0, 1))
    S = np.zeros(eps_k.shape + (len(nus),), complex)
    N2 = eps_k.size
    for m in range(-nb + 1, nb):
        Vm = Vf[..., abs(m)]
        Om = 2 * np.pi * m / beta
        G0 = 1.0 / (1j * (nus[None, None, :] + Om) + mu - eps_k[..., None])
        # sum_q V(q) G0(k+q) = sum_q V(-q) G0(k+q) (V even) = (V * G0)(k)  circular convolution
        conv = np.fft.ifft2(Vm[..., None] * np.fft.fft2(G0, axes=(0, 1)), axes=(0, 1))
        S += conv / N2 / beta
    return nus, U * n / 2 + (U / 8) * S


def tpsc_trace_check(res, nus, Sigma, eps_k, beta, U, n=1.0):
    """Returns ((T/N) sum_k (Sigma - U n/2) G0 + U n^2/4,  U <n_up n_dn>_TPSC)."""
    mu = res["mu"]
    G0 = 1.0 / (1j * nus[None, None, :] + mu - eps_k[..., None])
    corr = np.real(np.sum((Sigma - U * n / 2) * G0)) / eps_k.size / beta
    return float(corr + U * n**2 / 4), float(U * res["docc"])
