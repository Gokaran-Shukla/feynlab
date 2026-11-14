"""Local (impurity / atomic) self-energy approximations, ordered by coupling strength.

weak      : Hartree-Fock, second-order perturbation theory (2OPT) in U
moderate  : iterated perturbation theory (IPT, Georges & Kotliar 1992) inside DMFT
strong    : Hubbard-I (atomic self-energy), exact atomic limit,
            numerically exact CT-INT QMC (see feynlab.qmc.ctint)

All functions take/return arrays on positive fermionic Matsubara frequencies of a
feynlab.core.matsubara.MatsubaraGrid and use G(tau) = -<T c c^+>.
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import MatsubaraGrid


# ----------------------------------------------------------------------------
# exact atomic limit (single Hubbard site)
# ----------------------------------------------------------------------------
def atomic_g(iw, U: float, mu: float, n_other: float):
    """Exact one-site Green's function at fixed occupation n_other of the opposite spin
    (Hubbard-I form, exact for the isolated atom when n_other is its thermal value)."""
    return (1 - n_other) / (iw + mu) + n_other / (iw + mu - U)


def atomic_occupation(U: float, mu: float, beta: float) -> float:
    """Thermal <n_sigma> of an isolated Hubbard site."""
    w = np.array([1.0, np.exp(beta * mu), np.exp(beta * mu), np.exp(beta * (2 * mu - U))])
    w = w / w.sum()
    return float(w[1] + w[3])


def atomic_sigma(iw, U: float, mu: float, beta: float):
    """Exact atomic self-energy. At half filling (mu=U/2): Sigma = U/2 + U^2/(4 iw)."""
    n = atomic_occupation(U, mu, beta)
    g = atomic_g(iw, U, mu, n)
    return iw + mu - 1.0 / g


# ----------------------------------------------------------------------------
# Hartree-Fock
# ----------------------------------------------------------------------------
def hartree_sigma(U: float, n_other: float):
    """First-order (tadpole) diagram: Sigma_sigma = U <n_{-sigma}> (Hubbard contact
    interaction has no Fock term between equal spins)."""
    return U * n_other


# ----------------------------------------------------------------------------
# second-order diagram (the 'single pair-bubble' of Mattuck Sec. 9.8)
# ----------------------------------------------------------------------------
def second_order_sigma_tau(g0_tau_up, g0_tau_dn_rev, g0_tau_dn, U: float):
    """Sigma_up(tau) = -U^2 G0_up(tau) G0_dn(tau) G0_dn(-tau)
                     =  U^2 G0_up(tau) G0_dn(tau) G0_dn(beta - tau).
    g0_tau_dn_rev must be G0_dn(beta - tau) on the same grid."""
    return U**2 * g0_tau_up * g0_tau_dn * g0_tau_dn_rev


def second_order_sigma(grid: MatsubaraGrid, g0_iw, U: float, n0: float | None = None):
    """Paramagnetic second-order self-energy of an impurity with (Hartree-shifted)
    bath propagator g0_iw. Returns Sigma^(2)(iw) (Hartree part NOT included)."""
    c2 = 0.0 if n0 is None else 0.0
    g0_tau = grid.iw_to_tau(g0_iw, 1.0, c2)
    s_tau = second_order_sigma_tau(g0_tau, g0_tau[::-1], g0_tau, U)
    return grid.tau_to_iw(s_tau)


# ----------------------------------------------------------------------------
# Hubbard-I
# ----------------------------------------------------------------------------
def hubbard_one_sigma(iw, U: float, mu: float, n_other: float):
    """Sigma_HI(iw) = U n + U^2 n(1-n)/(iw + mu - U(1-n)), the self-energy of the atom
    with occupation n of the opposite spin."""
    return U * n_other + U**2 * n_other * (1 - n_other) / (iw + mu - U * (1 - n_other))


# ----------------------------------------------------------------------------
# Yamada-Yosida (exact weak-coupling result for the symmetric Anderson model)
# ----------------------------------------------------------------------------
def yamada_yosida_Z(U: float, Delta: float) -> float:
    """Quasiparticle weight of the symmetric Anderson impurity to O(U^2):
    1/Z = 1 + (3 - pi^2/4) (U / (pi Delta))^2   (Yamada, Prog. Theor. Phys. 53, 970 (1975))."""
    return 1.0 / (1.0 + (3.0 - np.pi**2 / 4.0) * (U / (np.pi * Delta)) ** 2)


def z_from_matsubara(iw, sigma_iw) -> float:
    """Z = [1 - Im Sigma(iw_0)/w_0]^{-1} (first-Matsubara estimate)."""
    return float(1.0 / (1.0 - sigma_iw[0].imag / iw[0].imag))


def z_from_matsubara_extrap(iw, sigma_iw, n: int = 3) -> float:
    """Z from a polynomial extrapolation of Im Sigma(iw_n)/w_n to w->0."""
    w = iw[:n].imag
    y = sigma_iw[:n].imag / w
    p = np.polyfit(w, y, n - 1)
    return float(1.0 / (1.0 - np.polyval(p, 0.0)))
