"""Further linear and nonlinear responses built on the band-model interface.

* orbital_moment      - modern-theory orbital magnetic moment m_n(k) (self-rotation part)
                        Xiao, Chang, Niu RMP 82, 1959 (2010) Eq. (3.17); valley moment of a gapped
                        Dirac cone: Xiao, Yao, Niu PRL 99, 236809 (2007).
* berry_curvature_dipole - D_ab = int f0 d_a Omega_b  (Fermi-sea) and
                        D_ab = int (-df/dE) v_a Omega_b (Fermi-surface; integration by parts).
                        Sodemann & Fu, PRL 115, 216806 (2015).
* edelstein_boltzmann - current-induced spin/orbital density, intraband (tau) term
                        Edelstein, Solid State Commun. 73, 233 (1990).
* edelstein_kubo      - constant-broadening Kubo formula, intraband (odd) + interband (even) parts
                        (Freimuth, Bluegel, Mokrousov PRB 90, 174423 (2014); Li et al. PRB 91, 134402 (2015)).
* kubo_bastin_T       - finite-temperature Hall response by Sommerfeld convolution of the T=0
                        energy-resolved Kubo-Bastin conductivity, sigma(mu,T) = int dE (-f') sigma(E,0).

Units: hbar = e = 1.  Hall conductivities in e^2/h when multiplied by 2 pi (as in intrinsic.py);
response densities per unit area per (e E tau) (Boltzmann) or per e E (Kubo).
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import fermi


def _band_matrices(model, kx, ky):
    H = model.H(kx, ky)
    hx, hy = model.dH(kx, ky)
    E, U = np.linalg.eigh(H)
    Ud = np.conj(np.swapaxes(U, -1, -2))
    return E, U, Ud, Ud @ hx @ U, Ud @ hy @ U


def _mfermi_prime(x, beta):
    """-df/dE (normalised to 1)."""
    y = np.clip(beta * x, -200, 200)
    return beta / (4 * np.cosh(0.5 * y) ** 2)


# ------------------------------------------------------------------------------------------
def orbital_moment(model, kx, ky):
    """m^z_n(k) = Im sum_{m != n} v^x_nm v^y_mn / (E_n - E_m)   (units |e|/hbar = 1, electron charge -|e|).

    Equals -(|e|/hbar) Im <d_x u_n|(H - E_n)|d_y u_n> = -i(e/2hbar)<grad u|x(H-E)|grad u>_z with e = -|e|,
    the self-rotation ('modern theory') orbital moment of the Bloch state.  Returns (E, m), shape (..., nb)."""
    E, U, Ud, vx, vy = _band_matrices(model, kx, ky)
    dE = E[..., :, None] - E[..., None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(np.abs(dE) > 1e-10, 1.0 / dE, 0.0)
    m = np.imag(np.einsum("...nm,...mn,...nm->...n", vx, vy, inv))
    return E, m


def dirac_valley_moment(v, Delta, k):
    """Analytic m(k) = Delta v^2 / [2 (Delta^2 + v^2 k^2)] for H = v k.sigma + Delta sz, identical for both
    bands (Xiao, Yao & Niu, PRL 99, 236809 (2007), with their gap = 2 Delta and v = sqrt(3) a t / 2)."""
    return Delta * v**2 / (2 * (Delta**2 + v**2 * k**2))


def orbital_magnetization(model, mu, grid, beta=1e4):
    """Modern-theory orbital magnetisation at T ~ 0 (electron charge -|e|, |e|/hbar = 1):
    M = int sum_n f_n [ m_n - (mu - E_n) Omega_n ]  (Ceresoli, Thonhauser, Vanderbilt, Resta PRB 74, 024408)."""
    from .intrinsic import operator_berry_curvature
    kx, ky, w = grid
    E, m = orbital_moment(model, kx, ky)
    _, Om = operator_berry_curvature(model, kx, ky)
    f = fermi(E - mu, beta)
    return float(np.sum(w[..., None] * f * (m - (mu - E) * Om)))


# ------------------------------------------------------------------------------------------
def berry_curvature_dipole(model, mu, grid, beta=200.0, dk=1e-4, form="surface"):
    """Berry curvature dipole components (D_x, D_y) = int d^2k/(2pi)^2 sum_n ... of Omega_z.

    form='surface':  D_a =  sum_n int (-df/dE) v^a_nn Omega_n
    form='sea'    :  D_a =  sum_n int f_n d_a Omega_n   (central finite differences, step dk)
    The two agree for a Fermi surface lying inside the grid (integration by parts); the
    nonlinear Hall current is j_a = -(e^3 tau/2) eps_{adc} D_{bd} E_b E_c (Sodemann-Fu)."""
    from .intrinsic import operator_berry_curvature
    kx, ky, w = grid
    if form == "surface":
        E, U, Ud, vx, vy = _band_matrices(model, kx, ky)
        _, Om = operator_berry_curvature(model, kx, ky)
        fp = _mfermi_prime(E - mu, beta)
        vxx = np.real(np.einsum("...nn->...n", vx))
        vyy = np.real(np.einsum("...nn->...n", vy))
        return (float(np.sum(w[..., None] * fp * vxx * Om)), float(np.sum(w[..., None] * fp * vyy * Om)))
    E, Om = operator_berry_curvature(model, kx, ky)
    f = fermi(E - mu, beta)
    out = []
    for sx_, sy_ in ((dk, 0.0), (0.0, dk)):
        _, Op = operator_berry_curvature(model, kx + sx_, ky + sy_)
        _, Om_ = operator_berry_curvature(model, kx - sx_, ky - sy_)
        out.append(float(np.sum(w[..., None] * f * (Op - Om_) / (2 * dk))))
    return tuple(out)


# ------------------------------------------------------------------------------------------
def edelstein_boltzmann(model, mu, op, grid, beta=500.0, direction="x"):
    """delta<O> per (e E tau), field along 'direction':  chi = - sum_n int (-df/dE) v^a_nn O_nn
    (electron charge -e; for E along x the drift shifts the distribution by +e E tau / hbar along -x)."""
    E, U, Ud, vx, vy = _band_matrices(model, *grid[:2])
    v = vx if direction == "x" else vy
    Ob = Ud @ op @ U
    fp = _mfermi_prime(E - mu, beta)
    vnn = np.real(np.einsum("...nn->...n", v))
    Onn = np.real(np.einsum("...nn->...n", Ob))
    return float(-np.sum(grid[2][..., None] * fp * vnn * Onn))


def rashba_edelstein_reference(m, alpha):
    """|delta S_y| = m alpha / (4 pi) per e E tau (spin s = sigma/2; both Rashba bands occupied)."""
    return m * alpha / (4 * np.pi)


def edelstein_kubo(model, mu, op, grid, Gamma=0.01, direction="x"):
    """Constant-broadening Kubo response chi = d<O>/d(eE) split into
      odd  (intraband, Fermi surface, ~1/Gamma):  - sum_n int v_nn O_nn  A_n^2 / (2 pi) ... -> Boltzmann with tau = 1/(2 Gamma)
      even (interband, intrinsic, Gamma -> 0):    sum_n f_n 2 Im sum_m O_nm v_mn / (E_n - E_m)^2
    The even part is the intrinsic (Berry-curvature-like) spin-orbit torque response."""
    kx, ky, w = grid
    E, U, Ud, vx, vy = _band_matrices(model, kx, ky)
    v = vx if direction == "x" else vy
    Ob = Ud @ op @ U
    # odd: (1/pi) int Re Tr[O G^R v G^A] restricted to n = m terms with Lorentzians
    A = (Gamma / np.pi) / ((mu - E) ** 2 + Gamma**2)
    vnn = np.real(np.einsum("...nn->...n", v))
    Onn = np.real(np.einsum("...nn->...n", Ob))
    odd = -float(np.sum(w[..., None] * Onn * vnn * A / (2 * Gamma)))
    dE = E[..., :, None] - E[..., None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(np.abs(dE) > 1e-10, dE**2 / (dE**2 + Gamma**2) ** 2, 0.0)   # regularised 1/dE^2
    f = (E < mu).astype(float)
    even = float(np.sum(w[..., None] * f * 2 * np.imag(np.einsum("...nm,...mn,...nm->...n", Ob, v, inv))))
    return dict(odd=odd, even=even, tau=1 / (2 * Gamma))


# ------------------------------------------------------------------------------------------
def kubo_bastin_T(sigma_of_E, mu, T, n_nodes=201, width=20.0):
    """sigma(mu, T) = int dE (-df/dE) sigma(E, T=0) by Gauss-Legendre quadrature on
    |E - mu| < width k_B T.  sigma_of_E: callable on an array of energies (e.g. a Kubo-Bastin or
    Kubo-Streda T=0 result evaluated at Fermi energy E)."""
    if T <= 0:
        return float(sigma_of_E(np.array([mu]))[0])
    x, wq = np.polynomial.legendre.leggauss(n_nodes)
    Es = mu + width * T * x
    weights = wq * width * T * _mfermi_prime(Es - mu, 1.0 / T)
    return float(np.sum(weights * sigma_of_E(Es)))
