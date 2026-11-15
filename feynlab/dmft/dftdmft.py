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

    def density(self, G, mu=None, sigma_c=None):
        """Total electron number (both spins) from G_loc with the exact 1/(iw)^2 tail coefficient
        c2_a = <H_aa>_k - mu + Re Sigma_a(i infinity)."""
        g = self.grid
        mu = self.mu if mu is None else mu
        hloc = np.real(np.diag(self.Hk.mean(0)))
        sinf = np.zeros(self.nw)
        if sigma_c is not None:
            for i, a in enumerate(self.corr):
                sinf[a] = float(np.real(sigma_c[i][-1]))
        n = 0.0
        for a in range(self.nw):
            n += 2 * g.density(G[:, a, a], c2=hloc[a] - mu + sinf[a])
        return n

    def find_mu(self, sigma_c, lo=-20.0, hi=20.0, tol=1e-6):
        for _ in range(80):
            mid = 0.5 * (lo + hi)
            n = self.density(self.gloc(sigma_c, mid), mid, sigma_c)
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
    def solve(self, solver="hubbard1", n_iter=20, mix=0.5, n_bath=1, tol=1e-5, verbose=False, cthyb_kwargs=None,
              sigma0=None):
        """solver: 'hubbard1' | 'ed' (bath fit, n_bath per orbital) | 'cthyb' (matrix CT-HYB, full Kanamori)."""
        iw = self.grid.iw
        sig = np.zeros((self.nc, len(iw)), complex) if sigma0 is None else np.array(sigma0, complex)
        self.mu = self.find_mu(sig)
        hist = []
        for it in range(n_iter):
            G = self.gloc(sig, self.mu)
            Gc = np.array([G[:, a, a] for a in self.corr])
            hl = np.real(np.diag(self.Hk.mean(0)))
            n_corr = sum(2 * self.grid.density(Gc[i], c2=hl[a] - self.mu + float(np.real(sig[i][-1])))
                         for i, a in enumerate(self.corr))
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
            elif solver == "cthyb":
                from ..qmc.cthyb_matrix import MatrixCTHYB
                kw = dict(n_meas=200_000, n_chains=4, meas_every=20)
                kw.update(cthyb_kwargs or {})
                nmeas_iw = kw.pop("n_iw_meas", 48)
                Delta = []
                for i in range(self.nc):
                    d_ = iw + self.mu - self.eps_loc[i] - sig[i] - 1.0 / Gc[i]
                    Delta += [d_, d_]
                solverq = MatrixCTHYB(self.grid, self.nc, self.U, self.J, eps=eps - self.mu, n_iw_meas=nmeas_iw)
                rq = solverq.solve(np.array(Delta), **kw)
                self.cthyb_last = rq
                m = len(rq["iw"])
                # high-frequency tail: Hartree-Fock constant of the Kanamori interaction from measured densities
                nfl = [rq[f"n_{f}"] for f in range(2 * self.nc)]
                Up = self.U - 2 * self.J
                sig_new = np.empty((self.nc, len(iw)), complex)
                for i in range(self.nc):
                    ss = 0.5 * (rq["G"][2 * i] + rq["G"][2 * i + 1])
                    se = np.abs(rq["G_err"][2 * i]) / np.sqrt(2)
                    s_low = iw[:m] - (eps[i] - self.mu) - Delta[2 * i][:m] - 1.0 / ss
                    s_inf = self.U * nfl[2 * i + 1] + sum(Up * nfl[2 * j + 1] + (Up - self.J) * nfl[2 * j]
                                                         for j in range(self.nc) if j != i)
                    # keep frequencies whose propagated QMC error |dG|/|G|^2 is small; beyond them use the
                    # analytic high-frequency form Sigma_inf + c/(iw) matched at the cutoff
                    err = se / np.abs(ss) ** 2
                    good = np.where(err > 0.05 + 0.05 * np.abs(s_low - s_inf))[0]
                    mc = int(good[0]) if len(good) else m
                    mc = max(mc, 2)
                    sig_new[i, :mc] = s_low[:mc]
                    sig_new[i, mc:] = s_inf + (s_low[mc - 1] - s_inf) * (iw[mc - 1] / iw[mc:])
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
        hl = np.real(np.diag(self.Hk.mean(0)))
        return dict(iw=iw, sigma=sig, G=G, mu=self.mu, n_total=self.density(G, self.mu, sig),
                    n_corr=sum(2 * self.grid.density(G[:, a, a], c2=hl[a] - self.mu + float(np.real(sig[i][-1])))
                               for i, a in enumerate(self.corr)), history=hist,
                    Z=[float(1 / (1 - sig[i, 0].imag / iw[0].imag)) for i in range(self.nc)])


    # ================================================================== charge self-consistency
    def density_matrix_k(self, sigma_c=None, mu=None):
        """Spin-resolved Wannier density matrices N(k)_ab = T sum_n G_ab(k, iw_n) e^{iw_n 0+} (per spin),
        with the analytic 1/(iw)^2 tail.  This is the quantity a DFT code needs for charge self-consistency
        (rho(r) = 2 sum_k sum_ab w_a,k(r) N(k)_ab w*_b,k(r)); see write_density_matrix."""
        sigma_c = getattr(self, "sigma", np.zeros((self.nc, self.grid.n_iw), complex)) if sigma_c is None else sigma_c
        mu = self.mu if mu is None else mu
        iw = self.grid.iw
        beta = self.grid.beta
        nw = self.nw
        S = np.zeros((len(iw), nw, nw), complex)
        for i, a in enumerate(self.corr):
            S[:, a, a] = sigma_c[i]
        eye = np.eye(nw)
        Sinf = np.real(S[-1])
        wn = iw.imag
        tail_sum = np.sum(1.0 / wn**2)
        out = np.empty((self.Hk.shape[0], nw, nw), complex)
        for k in range(self.Hk.shape[0]):
            G = np.linalg.inv((iw[:, None, None] + mu) * eye - self.Hk[k][None] - S)
            c2 = self.Hk[k] - mu * eye + Sinf
            Gs = G + np.conj(np.swapaxes(G, -1, -2))
            # (G + G^dagger) = -2 c2/w^2 + O(w^-4): subtract partial tail, add the exact sum beta^2/8 * (-2 c2)
            N = 0.5 * eye + (Gs.sum(0) + 2 * c2 * tail_sum - c2 * beta**2 / 4) / beta
            out[k] = 0.5 * (N + N.conj().T)
        return out

    def write_density_matrix(self, path, sigma_c=None, mu=None):
        """HDF5 file with N(k) (per spin), the k mesh index, mu and the correlated orbitals, for a DFT code
        (interfaces such as TRIQS/DFTTools or the VASP ICHARG=5 route read equivalent data)."""
        from ..utils.io import save
        N = self.density_matrix_k(sigma_c, mu)
        save(path, dict(N_k=N, mu=self.mu if mu is None else mu, corr=np.array(self.corr), n_total=float(2 * np.trace(N.mean(0)).real)))
        return path

    def solve_csc(self, functional, n_outer=20, mix_rho=0.5, tol_rho=1e-5, solver="hubbard1", dmft_iter=20,
                  verbose=False, **solve_kw):
        """Charge self-consistent DFT+DMFT at the level of the Wannier Hamiltonian:
            H(k) = H0(k) + V[n],   n_a = sum_k N(k)_aa (DMFT density matrix, both spins)
        where functional(n_orbital, N_local) -> local one-body potential matrix V (nw x nw), e.g. the Hartree /
        crystal-field response of the DFT functional projected on the Wannier orbitals (hartree_functional).
        The double counting is re-evaluated with the self-consistent occupations at every outer iteration.
        For a full ab-initio cycle, replace 'functional' by a call to the DFT code fed with write_density_matrix."""
        H0 = self.Hk.copy()
        eps0 = self.eps_loc.copy()
        n_prev = None
        hist = []
        sig = None
        res = None
        for it in range(n_outer):
            res = self.solve(solver=solver, n_iter=dmft_iter, sigma0=sig, **solve_kw)
            sig = res["sigma"]
            N = self.density_matrix_k(sig, res["mu"]).mean(0)
            n_orb = 2 * np.real(np.diag(N))
            if n_prev is not None:
                n_orb = mix_rho * n_orb + (1 - mix_rho) * n_prev
            V = functional(n_orb, 2 * N)
            self.Hk = H0 + V[None]
            self.eps_loc = np.real(np.diag(self.Hk.mean(0)))[self.corr]
            d = float(np.max(np.abs(n_orb - n_prev))) if n_prev is not None else np.inf
            hist.append(d)
            if verbose:
                print("CSC", it, d, n_orb)
            n_prev = n_orb
            if d < tol_rho:
                break
        res = self.solve(solver=solver, n_iter=dmft_iter, sigma0=sig, **solve_kw)
        res["csc_history"] = hist
        res["n_orbital"] = n_prev
        res["H_shift"] = self.Hk.mean(0) - H0.mean(0)
        self._H0, self._eps0 = H0, eps0
        return res


def hartree_functional(U_H):
    """Model DFT functional: orbital-resolved Hartree potential V_ab = delta_ab sum_c U_H[a, c] n_c
    (U_H symmetric matrix or scalar times (1 - delta)).  Captures the charge redistribution between
    correlated and uncorrelated orbitals that charge self-consistency corrects."""
    U_H = np.asarray(U_H, float)

    def f(n_orb, N_loc):
        nw = len(n_orb)
        M = U_H if U_H.ndim == 2 else U_H * (np.ones((nw, nw)) - np.eye(nw))
        return np.diag(M @ n_orb)
    return f


def hartree_fock_reference(Hk0, U_H, n_elec, beta, n_iter=500, mix=0.3):
    """Independent static mean-field solution H(k) = H0(k) + diag(U_H n) with Fermi occupations
    (benchmark for solve_csc with U_DMFT = 0)."""
    nw = Hk0.shape[-1]
    M = U_H if np.ndim(U_H) == 2 else U_H * (np.ones((nw, nw)) - np.eye(nw))
    n = np.full(nw, n_elec / nw)
    from ..core.matsubara import fermi
    for _ in range(n_iter):
        E, V = np.linalg.eigh(Hk0 + np.diag(M @ n)[None])
        lo, hi = E.min() - 5, E.max() + 5
        for _ in range(100):
            mu = 0.5 * (lo + hi)
            ntot = 2 * fermi(E - mu, beta).sum() / Hk0.shape[0]
            lo, hi = (mu, hi) if ntot < n_elec else (lo, mu)
        f = fermi(E - mu, beta)
        n_new = 2 * np.einsum("kan,kn->a", np.abs(V) ** 2, f) / Hk0.shape[0]
        if np.max(np.abs(n_new - n)) < 1e-10:
            n = n_new
            break
        n = mix * n_new + (1 - mix) * n
    return n, mu
