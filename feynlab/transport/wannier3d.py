"""Three-dimensional Wannier tight-binding models and the intrinsic Hall tensor with adaptive k refinement.

H(k) = sum_R e^{i k.R} H_R from a Wannier90 seedname_hr.dat ("tight-binding convention", R in Cartesian
Angstrom); velocities dH/dk_a = sum_R i R_a e^{ik.R} H_R.  The Berry-curvature pseudovector Omega_c
(c = x, y, z  <->  (yz, zx, xy)) and the operator-resolved (spin/orbital) generalisation are computed with the
Kubo sum-over-states formula, as in intrinsic.py.

Adaptive refinement (Wang, Yates, Souza, Vanderbilt PRB 74, 195118 (2006); Wannier90 'berry_curv_adpt_kmesh'):
points of a uniform Monkhorst-Pack mesh where |Omega| exceeds a threshold, or whose cell is crossed by the
Fermi surface (|E_n - mu| < |v_n| * half cell diagonal + 4 k_B T), are replaced by an n_sub^3 sub-mesh.  The result is in units of e^2/h per Angstrom; multiply by E2H_PER_A_TO_S_PER_CM for S/cm.
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import fermi
from ..dmft.dftdmft import read_hr

E2H_PER_A_TO_S_PER_CM = 3.874045865e-5 * 1e8   # (e^2/h)/Angstrom in S/cm


class WannierTB3D:
    def __init__(self, hr_file: str | None = None, lat=None, Rs=None, HR=None, operators: dict | None = None):
        self.lat = np.eye(3) if lat is None else np.asarray(lat, float)
        if hr_file is not None:
            Rs, HR = read_hr(hr_file)
        self.Rint = np.asarray(Rs, int)
        self.R = self.Rint @ self.lat
        self.HR = np.asarray(HR, complex)
        self.nb = self.HR.shape[-1]
        self.recip = 2 * np.pi * np.linalg.inv(self.lat).T
        self.vol = abs(np.linalg.det(self.lat))
        self.operators = operators or {}

    def frac_to_cart(self, kf):
        return kf @ self.recip

    def H(self, k):
        ph = np.exp(1j * k @ self.R.T)
        return np.einsum("kr,rab->kab", ph, self.HR)

    def dH(self, k):
        ph = np.exp(1j * k @ self.R.T)
        return [np.einsum("kr,rab->kab", 1j * self.R[:, a] * ph, self.HR) for a in range(3)]

    # ------------------------------------------------------------------
    def curvature(self, k, op=None, return_velocity=False):
        """Returns E (nk, nb) and Omega (nk, nb, 3) with Omega_c = -2 Im sum_m J^a_nm v^b_mn / dE^2,
        (a, b) = (y, z), (z, x), (x, y) for c = x, y, z;  J^a = {op, v_a}/2 (op=None -> charge)."""
        E, U = np.linalg.eigh(self.H(k))
        Ud = np.conj(np.swapaxes(U, -1, -2))
        v = [Ud @ h @ U for h in self.dH(k)]
        if op is None:
            J = v
        else:
            Ob = Ud @ op @ U
            J = [0.5 * (Ob @ va + va @ Ob) for va in v]
        dE = E[:, :, None] - E[:, None, :]
        with np.errstate(divide="ignore", invalid="ignore"):
            inv = np.where(np.abs(dE) > 1e-8, 1.0 / dE**2, 0.0)
        Om = np.empty(E.shape + (3,))
        for c, (a, b) in enumerate(((1, 2), (2, 0), (0, 1))):
            Om[..., c] = -2 * np.imag(np.einsum("knm,kmn,knm->kn", J[a], v[b], inv))
        if return_velocity:
            vn = np.sqrt(sum(np.real(np.einsum("knn->kn", va)) ** 2 for va in v))
            return E, Om, vn
        return E, Om

    def _integrand(self, kf, mu, beta, op, half_diag=0.0):
        E, Om, vn = self.curvature(self.frac_to_cart(kf), op, return_velocity=True)
        f = fermi(E - mu, beta)
        # a cell is cut by the Fermi surface if some band has |E - mu| < |v| * (half cell diagonal) + 4 k_B T
        self._last_gap = np.min(np.abs(E - mu) - vn * half_diag, axis=1) - 4.0 / beta
        return np.einsum("kn,knc->kc", f, Om)

    def intrinsic_hall(self, mu: float, nk=(20, 20, 20), op=None, beta: float = 1e4,
                       adaptive: bool = True, threshold: float | None = None, n_sub: int = 5, chunk: int = 20000,
                       window: float | None = None):
        """sigma_c (c = x,y,z  <->  sigma_yz, sigma_zx, sigma_xy) in e^2/h per Angstrom
        (e/2pi per Angstrom for an operator current).  Returns dict with sigma, the number of refined
        points and the total number of k points used."""
        nk = np.asarray(nk)
        axes = [(np.arange(n) + 0.5) / n for n in nk]
        K = np.stack(np.meshgrid(*axes, indexing="ij"), -1).reshape(-1, 3)
        hd = 0.5 * np.sqrt(np.sum((np.linalg.norm(self.recip, axis=1) / nk) ** 2))   # half cell diagonal
        vals, gaps = [], []
        for i in range(0, len(K), chunk):
            vals.append(self._integrand(K[i:i + chunk], mu, beta, op, hd))
            gaps.append(self._last_gap)
        vals, gaps = np.concatenate(vals), np.concatenate(gaps)
        w0 = 1.0 / nk.prod()
        nref = 0
        if adaptive:
            mag = np.linalg.norm(vals, axis=1)
            # hot-spot refinement is opt-in (threshold): for gapped regions the uniform midpoint rule is already
            # spectrally accurate, and mixing resolutions there only adds O(h^2) quadrature error
            thr = threshold if threshold is not None else np.inf
            # refine (i) optional Berry-curvature hot spots and (ii) cells actually cut by the Fermi surface
            # (|E - mu| < |v| * half cell diagonal + 4 kT); 'window' adds an extra energy margin
            sel = np.where((mag > thr) | (gaps < (0.0 if window is None else window)))[0]
            nref = len(sel)
            if nref:
                s = (np.arange(n_sub) + 0.5) / n_sub - 0.5
                sub = np.stack(np.meshgrid(s, s, s, indexing="ij"), -1).reshape(-1, 3) / nk
                for i in range(0, nref, max(1, chunk // n_sub**3)):
                    idx = sel[i:i + max(1, chunk // n_sub**3)]
                    Ks = (K[idx][:, None, :] + sub[None]).reshape(-1, 3)
                    vs = self._integrand(Ks, mu, beta, op).reshape(len(idx), -1, 3).mean(1)
                    vals[idx] = vs
        # sigma = -2pi int d^3k/(2pi)^3 f Omega  ->  -(2pi) * (1/V_cell) * (1/(2pi)^3 * (2pi)^3) mean(...)
        sigma = -2 * np.pi * w0 * vals.sum(0) / self.vol
        return dict(sigma=sigma, sigma_S_per_cm=sigma * E2H_PER_A_TO_S_PER_CM, n_refined=nref,
                    n_k=int(nk.prod() + nref * (n_sub**3 - 1)))


def stacked_qwz_hr(m: float = 1.0, tz: float = 0.2):
    """Rs, HR of QWZ layers stacked along z with spin-independent interlayer hopping tz:
    H = sin kx sx + sin ky sy + (m + cos kx + cos ky) sz + 2 tz cos kz.
    Each layer is a Chern insulator (|C| = 1 for 0 < |m| < 2) for |2tz| below the gap, so the 3D
    quantised Hall conductivity is sigma_xy = -C e^2/(h c)  (Halperin 1987)."""
    sx = np.array([[0, 1], [1, 0]], complex); sy = np.array([[0, -1j], [1j, 0]]); sz = np.diag([1.0 + 0j, -1])
    terms = {(0, 0, 0): m * sz,
             (1, 0, 0): 0.5 * sz - 0.5j * sx, (-1, 0, 0): 0.5 * sz + 0.5j * sx,
             (0, 1, 0): 0.5 * sz - 0.5j * sy, (0, -1, 0): 0.5 * sz + 0.5j * sy,
             (0, 0, 1): tz * np.eye(2), (0, 0, -1): tz * np.eye(2)}
    Rs = np.array(list(terms.keys()))
    HR = np.array(list(terms.values()))
    return Rs, HR
