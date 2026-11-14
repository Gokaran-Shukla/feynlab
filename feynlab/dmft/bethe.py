"""Single-site DMFT for the Hubbard model on the Bethe lattice (half bandwidth D).

Self-consistency (Georges et al., RMP 68, 13 (1996), Eq. 32):
    G0^{-1}(iw) = iw + mu - t^2 G(iw),   t = D/2.

Solvers
-------
'ipt'       : iterated 2nd-order perturbation theory (half filling)  -> weak..moderate,
              qualitatively captures the Mott transition
'hubbard1'  : atomic self-energy                                       -> strong coupling
'ctint'     : numerically exact CT-INT QMC                              -> any coupling
'2opt'      : plain second order around the *non-interacting* bath (no iteration)
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import MatsubaraGrid
from ..core.lattice import bethe_hilbert
from ..selfenergy import local


class BetheDMFT:
    def __init__(self, beta: float, U: float, D: float = 1.0, mu: float | None = None,
                 n_iw: int = 1024):
        self.grid = MatsubaraGrid(beta, n_iw)
        self.U = float(U)
        self.D = float(D)
        self.t = D / 2.0
        self.mu = U / 2.0 if mu is None else float(mu)
        self.history = []

    # ------------------------------------------------------------------
    # ---- lattice-specific pieces (overridden by LatticeDMFT) ----------
    def hilbert(self, z):
        return bethe_hilbert(z, self.D)

    @property
    def second_moment(self):
        """int e^2 rho(e) de  (= t^2 = D^2/4 for the Bethe lattice): exact weight of Delta."""
        return self.t**2

    def g_lattice(self, sigma):
        return self.hilbert(self.grid.iw + self.mu - sigma)

    def weiss(self, g):
        """G0^{-1} = G^{-1} + Sigma  (for the Bethe lattice identical to iw + mu - t^2 G)."""
        sig = getattr(self, "_sigma_current", None)
        if sig is None:
            return 1.0 / (self.grid.iw + self.mu - self.t**2 * g)
        return 1.0 / (1.0 / g + sig)

    # ------------------------------------------------------------------
    def _sigma_ipt(self, g0):
        """Half-filled IPT: Sigma = U/2 + U^2 G0~(tau)^2 G0~(beta - tau) with the
        Hartree-shifted Weiss field G0~^{-1} = G0^{-1} - U/2."""
        grid = self.grid
        g0t = 1.0 / (1.0 / g0 - self.U / 2.0)
        if abs(self.mu - self.U / 2.0) < 1e-12:
            g0t = 1j * g0t.imag
        s2 = local.second_order_sigma(grid, g0t, self.U)
        return self.U / 2.0 + s2

    def _sigma_hubbard1(self, g):
        n = self.grid.density(g, c2=None)
        return local.hubbard_one_sigma(self.grid.iw, self.U, self.mu, n)

    # ------------------------------------------------------------------
    def solve(self, solver: str = "ipt", n_iter: int = 200, mix: float = 0.5, tol: float = 1e-7,
              sigma_init=None, ctint_kwargs: dict | None = None, insulating_start: bool = False,
              verbose: bool = False):
        iw = self.grid.iw
        if sigma_init is not None:
            sigma = np.asarray(sigma_init, complex)
        elif insulating_start:
            sigma = self.U / 2.0 + self.U**2 / (4.0 * iw)
        else:
            sigma = np.full_like(iw, self.U / 2.0)
        g = self.g_lattice(sigma)
        self.ctint_last = None
        for it in range(n_iter):
            self._sigma_current = sigma
            g0 = self.weiss(g)
            if solver == "ipt":
                sigma_new = self._sigma_ipt(g0)
            elif solver == "hubbard1":
                sigma_new = self._sigma_hubbard1(g)
            elif solver == "ctint":
                from ..qmc.ctint import CTINTSolver
                kw = dict(n_warm=10000, n_meas=300000, meas_every=10, n_chains=2)
                kw.update(ctint_kwargs or {})
                n_meas_iw = kw.pop("n_iw_meas", 64)
                sol = CTINTSolver(self.grid, self.U, n_iw_meas=n_meas_iw)
                res = sol.solve(g0, seed=1000 + it, **kw)
                self.ctint_last = res
                # high-frequency part: atomic-like tail Sigma = U/2 + U^2 n(1-n)/iw
                s_tail = self.U / 2.0 + self.U**2 * 0.25 / iw
                sigma_new = s_tail.copy()
                m = len(res["Sigma_up"])
                sym = 0.5 * (res["Sigma_up"] + res["Sigma_dn"])
                if abs(self.mu - self.U / 2) < 1e-12:  # enforce p-h symmetry exactly
                    sym = self.U / 2.0 + 1j * sym.imag
                sigma_new[:m] = sym
            elif solver == "ed":
                from ..qmc.ed_sector import SectorED, fit_bath
                kw = dict(n_bath=5, beta_fit=200.0, n_fit=200, beta_ed=1e4, ed_method="auto")
                kw.update(ctint_kwargs or {})
                iwf = 1j * np.pi * (2 * np.arange(kw["n_fit"]) + 1) / kw["beta_fit"]
                # hybridization on the fitting grid from the current lattice G
                sig_fit = self._sigma_on(iwf, sigma)
                g_fit = self.hilbert(iwf + self.mu - sig_fit)
                delta_fit = iwf + self.mu - sig_fit - 1.0 / g_fit
                eb, Vb, self._bath_x, _ = fit_bath(iwf, delta_fit, kw["n_bath"],
                                                    x0=getattr(self, "_bath_x", None),
                                                    symmetric=abs(self.mu - self.U / 2) < 1e-12,
                                                    sum_rule=self.second_moment)
                from ..qmc.ed_sector import LanczosED
                use_lanczos = kw["ed_method"] == "lanczos" or (kw["ed_method"] == "auto" and kw["beta_ed"] >= 1e4)
                ed = LanczosED(-self.mu, self.U, eb, Vb) if use_lanczos else SectorED(-self.mu, self.U, eb, Vb)
                self.ed_last = ed
                self._beta_ed = kw["beta_ed"]   # cluster temperature (1e4 -> ground state, T=0 DMFT)
                Gimp = ed.green(iw, self._beta_ed)
                sigma_new = iw + self.mu - ed.hybridization(iw) - 1.0 / Gimp
                self._ed_sigma_fn = (ed, )
            elif solver == "2opt":
                g00 = 1.0 / (iw + self.mu - self.U / 2 - self.t**2 * bethe_hilbert(iw, self.D))
                sigma_new = self.U / 2.0 + local.second_order_sigma(self.grid, g00, self.U)
                sigma = sigma_new
                g = self.g_lattice(sigma)
                break
            else:
                raise ValueError(f"unknown solver {solver}")
            if solver == "ed" and insulating_start and it == 0:
                pass
            if abs(self.mu - self.U / 2.0) < 1e-12:
                # particle-hole symmetry at half filling: Re Sigma = U/2 exactly
                sigma_new = self.U / 2.0 + 1j * sigma_new.imag
            diff = float(np.max(np.abs(sigma_new - sigma)))
            sigma = mix * sigma_new + (1 - mix) * sigma
            g = self.g_lattice(sigma)
            self.history.append(diff)
            if verbose:
                print(f"iter {it:3d}  |dSigma| = {diff:.3e}")
            if diff < tol:
                break
        self.sigma = sigma
        self.g = g
        return dict(iw=iw, G=g, Sigma=sigma, iterations=len(self.history),
                    converged=self.history[-1] < tol if self.history else True,
                    Z=local.z_from_matsubara(iw, sigma),
                    n=self.grid.density(g))

    def _sigma_on(self, z, sigma):
        """Sigma on another frequency set: from the last ED solution if available,
        otherwise by nearest-frequency interpolation of the current Sigma."""
        if getattr(self, "_ed_sigma_fn", None):
            ed = self._ed_sigma_fn[0]
            return z + self.mu - ed.hybridization(z) - 1.0 / ed.green(z, getattr(self, "_beta_ed", self.grid.beta))
        w = self.grid.iw.imag
        return np.interp(z.imag, w, sigma.real) + 1j * np.interp(z.imag, w, sigma.imag)

    # ------------------------------------------------------------------
    def spectral_function(self, omega, eta: float = 0.02, n_use: int = 64):
        from ..core.lattice import pade_continue
        sig = pade_continue(self.grid.iw, self.sigma, omega, eta=eta, n_use=n_use)
        g = self.hilbert(np.asarray(omega) + 1j * eta + self.mu - sig)
        return -g.imag / np.pi


class LatticeDMFT(BetheDMFT):
    """Single-site DMFT for any lattice, given either a k-mesh dispersion eps_k or a
    tabulated DOS (energies, weights).  Examples:
        LatticeDMFT.square(beta, U, t=0.25, nk=64)      # 2D square, W = 8t
        LatticeDMFT.cubic(beta, U, t=1/6, nk=32)        # 3D simple cubic, W = 12t
        LatticeDMFT.hypercubic(beta, U, t=1/sqrt(2))    # Gaussian DOS (d -> inf)
    """

    def __init__(self, beta, U, energies, weights, mu=None, n_iw=1024):
        super().__init__(beta, U, D=1.0, mu=mu, n_iw=n_iw)
        self.e = np.asarray(energies, float)
        self.wt = np.asarray(weights, float) / np.sum(weights)

    def hilbert(self, z):
        z = np.asarray(z, complex)
        return (self.wt[None, :] / (z[:, None] - self.e[None, :])).sum(1)

    @property
    def second_moment(self):
        return float(np.sum(self.wt * self.e**2) - np.sum(self.wt * self.e) ** 2)

    def weiss(self, g):
        sig = getattr(self, "_sigma_current", None)
        if sig is None:
            sig = np.full_like(g, self.U / 2.0)
        return 1.0 / (1.0 / g + sig)

    @classmethod
    def from_eps_k(cls, beta, U, eps_k, n_bins: int = 4000, **kw):
        h, edges = np.histogram(np.ravel(eps_k), bins=n_bins)
        c = 0.5 * (edges[1:] + edges[:-1])
        return cls(beta, U, c[h > 0], h[h > 0], **kw)

    @classmethod
    def square(cls, beta, U, t=0.25, nk=128, **kw):
        k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
        kx, ky = np.meshgrid(k, k, indexing="ij")
        return cls.from_eps_k(beta, U, -2 * t * (np.cos(kx) + np.cos(ky)), **kw)

    @classmethod
    def cubic(cls, beta, U, t=1 / 6, nk=48, **kw):
        k = 2 * np.pi * (np.arange(nk) + 0.5) / nk
        kx, ky, kz = np.meshgrid(k, k, k, indexing="ij")
        return cls.from_eps_k(beta, U, -2 * t * (np.cos(kx) + np.cos(ky) + np.cos(kz)), **kw)

    @classmethod
    def hypercubic(cls, beta, U, t=1 / np.sqrt(2), n: int = 4001, **kw):
        e = np.linspace(-6, 6, n)
        return cls(beta, U, e, np.exp(-e**2 / (2 * t**2)), **kw)

    @classmethod
    def semicircle(cls, beta, U, D=1.0, n: int = 4001, **kw):
        e = np.linspace(-D, D, n)
        return cls(beta, U, e, np.sqrt(np.clip(D**2 - e**2, 0, None)), **kw)
