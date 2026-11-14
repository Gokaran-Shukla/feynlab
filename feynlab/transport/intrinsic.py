"""Intrinsic (Berry-curvature / Karplus-Luttinger) Hall responses for any operator.

For an operator O (1 -> charge, S_z -> spin, L_z -> orbital, Q_ab -> orbital
quadrupole ...) with current J^O_x = {O, v_x}/2:

    sigma^O_xy = -e * sum_n int d^2k/(2pi)^2 f_n(k) Omega^O_n(k),
    Omega^O_n  = -2 Im sum_{m != n} <n|J^O_x|m><m|v_y|n> / (e_n - e_m)^2 .

For O = e (charge) this is the TKNN / Karplus-Luttinger formula and Omega is the
Berry curvature (Omega = -2 Im <d_x u|d_y u>), giving sigma_xy = -e^2 sum f Omega
(Sinitsyn et al. 2007 convention).  Results are returned in units of
    e^2/h  for charge   and   e/(2 pi)  for other operators (hbar = 1),
so a Chern insulator gives sigma_xy = -C (with C = (1/2pi) int Omega).
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import fermi


def k_grid_lattice(model, nk: int = 200, shift: bool = True):
    """Uniform Monkhorst-Pack grid over the BZ. Returns kx, ky, weight (sum = A_BZ/(2pi)^2)."""
    if hasattr(model, "recip"):
        b = np.asarray(model.recip)
    else:
        b = 2 * np.pi * np.eye(2)
    s = (np.arange(nk) + (0.5 if shift else 0.0)) / nk
    u, v = np.meshgrid(s, s, indexing="ij")
    kx = u * b[0, 0] + v * b[1, 0]
    ky = u * b[0, 1] + v * b[1, 1]
    area = abs(np.linalg.det(b))
    w = np.full(kx.shape, area / nk**2 / (2 * np.pi) ** 2)
    return kx, ky, w


def k_grid_polar(kmax: float, nr: int = 600, nth: int = 64, dense: list | None = None):
    """Polar grid for continuum models; 'dense' = list of (k0, width, n) windows
    (e.g. around Fermi momenta).  Returns kx, ky, weights for int d^2k/(2pi)^2."""
    edges = list(np.linspace(0, kmax, nr + 1))
    for k0, wdt, n in dense or []:
        edges += list(np.linspace(max(0, k0 - wdt), min(kmax, k0 + wdt), n))
    edges = np.unique(np.clip(edges, 0, kmax))
    kr = 0.5 * (edges[1:] + edges[:-1])
    dk = np.diff(edges)
    th = 2 * np.pi * (np.arange(nth) + 0.5) / nth
    K, T = np.meshgrid(kr, th, indexing="ij")
    W = (K * dk[:, None]) * (2 * np.pi / nth) / (2 * np.pi) ** 2
    return K * np.cos(T), K * np.sin(T), W


def operator_berry_curvature(model, kx, ky, op=None):
    """Omega^O_n(k) for all bands; op=None -> charge (ordinary Berry curvature)."""
    H = model.H(kx, ky)
    hx, hy = model.dH(kx, ky)
    E, U = np.linalg.eigh(H)
    Ud = np.conj(np.swapaxes(U, -1, -2))
    vx = Ud @ hx @ U
    vy = Ud @ hy @ U
    if op is None:
        jx = vx
    else:
        Ob = Ud @ op @ U
        jx = 0.5 * (Ob @ vx + vx @ Ob)
    dE = E[..., :, None] - E[..., None, :]
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(np.abs(dE) > 1e-10, 1.0 / dE**2, 0.0)
    # Omega_n = -2 Im sum_m jx[n,m] vy[m,n] / (E_n-E_m)^2
    omega = -2.0 * np.imag(np.einsum("...nm,...mn,...nm->...n", jx, vy, inv))
    return E, omega


def intrinsic_hall(model, mu: float, op=None, beta: float = 1e4, grid=None, nk: int = 300,
                   kmax: float | None = None):
    """sigma^O_xy (units e^2/h for charge, e/2pi otherwise). grid=(kx,ky,w) optional."""
    if grid is None:
        if model.kind == "lattice":
            grid = k_grid_lattice(model, nk)
        else:
            grid = k_grid_polar(kmax or 50.0, nr=nk * 4, nth=32)
    kx, ky, w = grid
    E, om = operator_berry_curvature(model, kx, ky, op)
    f = fermi(E - mu, beta)
    sig = -np.sum(w[..., None] * f * om)
    return float(2 * np.pi * sig)


def chern_numbers(model, nk: int = 200):
    """Band Chern numbers by the Fukui-Hatsugai-Suzuki lattice method
    (J. Phys. Soc. Jpn. 74, 1674 (2005)); exact integers on a finite grid."""
    b = np.asarray(getattr(model, "recip", 2 * np.pi * np.eye(2)))
    s = np.arange(nk) / nk
    u, v = np.meshgrid(s, s, indexing="ij")
    kx = u * b[0, 0] + v * b[1, 0]
    ky = u * b[0, 1] + v * b[1, 1]
    _, U = np.linalg.eigh(model.H(kx, ky))
    nb = U.shape[-1]
    C = []
    for n in range(nb):
        un = U[..., :, n]
        def link(a, b_):
            z = np.sum(np.conj(a) * b_, axis=-1)
            return z / np.abs(z)
        u1 = link(un, np.roll(un, -1, axis=0))
        u2 = link(un, np.roll(un, -1, axis=1))
        u2s = link(np.roll(un, -1, axis=0), np.roll(np.roll(un, -1, axis=0), -1, axis=1))
        u1s = link(np.roll(un, -1, axis=1), np.roll(np.roll(un, -1, axis=0), -1, axis=1))
        F = np.angle(u1 * u2s / (u1s * u2))
        # FHS link variables use A = <u|du>; convert to Omega = -2 Im<d_x u|d_y u>
        C.append(-F.sum() / (2 * np.pi))
    return np.array(C)
