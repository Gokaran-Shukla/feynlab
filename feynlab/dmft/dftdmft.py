"""DFT+DMFT interface for Wannier90 tight-binding Hamiltonians (one-shot, i.e.
without charge self-consistency).

    G_loc(iw) = (1/N_k) sum_k [ (iw + mu) 1 - H(k) - Sigma(iw) + Sigma_dc ]^{-1}
    Delta_a(iw) = iw + mu - eps_a - Sigma_a - 1/G_loc,aa          (diagonal hybridization)

* Wannier input: seedname_hr.dat (any number of orbitals), correlated subspace = a list of
  orbital indices; uncorrelated orbitals stay in the k-sum.
* Double counting: 'FLL' (fully localised limit), 'AMF' (around mean field), or 'occupation'
  (Hubbard-I only: level shift fixed so the atom carries the lattice occupation) for the
  Kanamori shell (Anisimov et al. PRB 48, 16929 (1993); Czyzyk & Sawatzky PRB 49, 14211 (1994)).
* Chemical potential: bisection to a target total electron count per unit cell.
* Impurity solvers: multi-orbital 'hubbard1' (Lichtenstein & Katsnelson PRB 57, 6884 (1998))
  or multi-orbital 'ed' (Kanamori + bath sites, diagonal hybridization).
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import MatsubaraGrid
from ..qmc.multiorbital import KanamoriAtom, MultiOrbitalED
from ..qmc.ed_sector import fit_bath


def read_hr(path):
    """Wannier90 seedname_hr.dat -> (R vectors (n_R,3) int, H_R (n_R, nw, nw))."""
    with open(path) as f:
        f.readline()
        nw = int(f.readline()); nr = int(f.readline())
        deg = []
        while len(deg) < nr:
            deg += [int(x) for x in f.readline().split()]
        data = np.loadtxt(f)
    R = data[:, :3].astype(int)
    Rs, inv = np.unique(R, axis=0, return_inverse=True)
    HR = np.zeros((len(Rs), nw, nw), complex)
    degR = np.zeros(len(Rs))
    for row, ir in zip(data, inv.ravel()):
        a, b = int(row[3]) - 1, int(row[4]) - 1
        HR[ir, a, b] = row[5] + 1j * row[6]
    # degeneracies listed in the order R first appears
    seen = []
    for r in R:
        t = tuple(r)
        if t not in seen:
            seen.append(t)
    for k, t in enumerate(seen):
        ir = np.where((Rs == np.array(t)).all(1))[0][0]
        degR[ir] = deg[k]
    return Rs, HR / degR[:, None, None]


def write_hr(path, Rs, HR, deg=None):
    nR, nw, _ = HR.shape
    deg = np.ones(nR, int) if deg is None else deg
    with open(path, "w") as f:
        f.write("written by feynlab\n%d\n%d\n" % (nw, nR))
        for i in range(0, nR, 15):
            f.write(" ".join(str(d) for d in deg[i:i + 15]) + "\n")
        for ir in range(nR):
            for b in range(nw):
                for a in range(nw):
                    h = HR[ir, a, b] * deg[ir]
                    f.write("%5d%5d%5d%5d%5d %12.6f %12.6f\n" % (*Rs[ir], a + 1, b + 1, h.real, h.imag))


def hk_on_mesh(Rs, HR, nk=(16, 16, 1)):
    k = [np.arange(n) / n for n in nk]
    K = np.stack(np.meshgrid(*k, indexing="ij"), -1).reshape(-1, 3)
    ph = np.exp(2j * np.pi * K @ Rs.T)
    return np.einsum("kr,rab->kab", ph, HR)


class DFTDMFT:
    def __init__(self, hr_file=None, Hk=None, corr=None, U=4.0, J=0.6, n_elec=2.0, beta=40.0, nk=(16, 16, 1),
                 n_iw=512, dc="FLL"):
        if Hk is None:
            Rs, HR = read_hr(hr_file)
            Hk = hk_on_mesh(Rs, HR, nk)
        self.Hk = Hk
        self.nw = Hk.shape[-1]
        self.corr = list(range(self.nw)) if corr is None else list(corr)
        self.nc = len(self.corr)
        self.U, self.J, self.dc_kind = U, J, dc
        self.n_elec = n_elec
        self.grid = MatsubaraGrid(beta, n_iw)
        self.eps_loc = np.real(np.diag(Hk.mean(0)))[self.corr]
        self.mu = 0.0

    # --------------------------------------------------------------
    def gloc(self, sigma_c, mu):
        """sigma_c: (n_corr, n_iw) correlated self-energy (with dc already subtracted)."""
        iw = self.grid.iw
        nw = self.nw
        S = np.zeros((len(iw), nw, nw), complex)
        for i, a in enumerate(self.corr):
            S[:, a, a] = sigma_c[i]
        G = np.zeros((len(iw), nw, nw), complex)
        eye = np.eye(nw)
        for k in range(self.Hk.shape[0]):
            G += np.linalg.inv((iw[:, None, None] + mu) * eye - self.Hk[k][None] - S)
        return G / self.Hk.shape[0]

    def density(self, G):
        g = self.grid
        n = 0.0
        for a in range(self.nw):
            n += 2 * g.density(G[:, a, a], c2=None)
        return n

    def find_mu(self, sigma_c, lo=-20.0, hi=20.0, tol=1e-6):
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            n = self.density(self.gloc(sigma_c, mid))
            if n > self.n_elec:
                hi = mid
            else:
                lo = mid
            if hi - lo < tol:
                break
        return 0.5 * (lo + hi)

    def double_counting(self, n_corr):
        U, J, N = self.U, self.J, n_corr
        M = self.nc
        Ubar = (U + 2 * (M - 1) * (U - 2 * J) + 2 * (M - 1) * (U - 3 * J)) / (2 * M - 1) if M > 1 else U
        Jbar = U - Ubar + (U - 2 * J - (U - 3 * J)) if M > 1 else 0.0
        if self.dc_kind == "FLL":
            return Ubar * (N - 0.5) - Jbar * (N / 2 - 0.5)
        if self.dc_kind == "AMF":
            return Ubar * N * (2 * M - 1) / (2 * M) - Jbar * N * (M - 1) / (2 * M) if M > 1 else U * N / 2
        return 0.0

    # --------------------------------------------------------------
    def solve(self, solver="hubbard1", n_iter=20, mix=0.5, n_bath=1, tol=1e-5, verbose=False):
        iw = self.grid.iw
        sig = np.zeros((self.nc, len(iw)), complex)
        self.mu = self.find_mu(sig)
        hist = []
        for it in range(n_iter):
            G = self.gloc(sig, self.mu)
            Gc = np.array([G[:, a, a] for a in self.corr])
            n_corr = sum(2 * self.grid.density(Gc[i]) for i in range(self.nc))
            vdc = self.double_counting(n_corr)
            eps = self.eps_loc - vdc
            if solver == "hubbard1":
                if self.dc_kind == "occupation":
                    # level shift s such that the atom holds the lattice occupation n_corr.
                    # At low T, N(s) has Mott plateaus: take the midpoint of the shifts giving
                    # N = n_corr -/+ 1/2, i.e. the plateau centre for integer fillings.
                    def shift_for(target):
                        lo, hi = -6 * self.U - 10, 6 * self.U + 10
                        for _ in range(70):
                            sh = 0.5 * (lo + hi)
                            nat = KanamoriAtom(self.nc, self.U, self.J, eps=self.eps_loc + sh,
                                               mu=self.mu).occupation(self.grid.beta)
                            lo, hi = (sh, hi) if nat > target else (lo, sh)
                        return 0.5 * (lo + hi)
                    sh = 0.5 * (shift_for(n_corr - 0.5) + shift_for(n_corr + 0.5))
                    vdc = -sh
                    eps = self.eps_loc - vdc
                atom = KanamoriAtom(self.nc, self.U, self.J, eps=eps, mu=self.mu)
                Gat = atom.green(iw, self.grid.beta)
                sig_new = np.array([iw + self.mu - eps[i] - 1.0 / Gat[i] for i in range(self.nc)])
            elif solver == "ed":
                eb, Vb = [], []
                for i in range(self.nc):
                    delta = iw + self.mu - self.eps_loc[i] - (sig[i] - 0) - 1.0 / Gc[i]
                    e_, v_, _, _ = fit_bath(iw[:120], delta[:120], n_bath, symmetric=False)
                    eb.append(e_ - self.mu); Vb.append(v_)
                ed = MultiOrbitalED(self.nc, self.U, self.J, eps - self.mu, eb, Vb)
                sig_new = np.array([iw - (eps[i] - self.mu) - ed.hybridization(iw, i) - 1.0 / ed.green(iw, self.grid.beta, i)
                                    for i in range(self.nc)])
            else:
                raise ValueError(solver)
            sig_new = sig_new - vdc          # Sigma_imp - Sigma_dc enters the lattice G
            diff = float(np.max(np.abs(sig_new - sig)))
            sig = mix * sig_new + (1 - mix) * sig
            self.mu = self.find_mu(sig)
            hist.append(diff)
            if verbose:
                print(it, diff, self.mu, n_corr)
            if diff < tol:
                break
        self.sigma = sig
        G = self.gloc(sig, self.mu)
        return dict(iw=iw, sigma=sig, G=G, mu=self.mu, n_total=self.density(G),
                    n_corr=sum(2 * self.grid.density(G[:, a, a]) for a in self.corr), history=hist,
                    Z=[float(1 / (1 - sig[i, 0].imag / iw[0].imag)) for i in range(self.nc)])
