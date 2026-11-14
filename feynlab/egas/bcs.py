"""BCS gap equation (Mattuck Ch. 15; Nambu ladder sum in the Cooper channel).

    1 = lambda int_0^{w_D} dxi tanh(E/2T)/E,   E = sqrt(xi^2 + Delta^2)
T_c (Delta -> 0): 1 = lambda int_0^{w_D} tanh(xi/2T)/xi  ->  T_c = 1.134 w_D e^{-1/lambda}
Universal weak-coupling ratio Delta(0)/T_c = pi e^{-gamma_E} = 1.7639.
"""
from __future__ import annotations

import numpy as np
from scipy import integrate, optimize


def gap_equation_residual(Delta, T, lam, wD):
    f = lambda xi: np.tanh(np.sqrt(xi**2 + Delta**2) / (2 * T)) / np.sqrt(xi**2 + Delta**2) if T > 0 \
        else 1.0 / np.sqrt(xi**2 + max(Delta, 1e-300) ** 2)
    val, _ = integrate.quad(f, 0, wD, limit=400, points=[min(wD, 5 * max(Delta, T))])
    return lam * val - 1.0


def gap_T0(lam, wD=1.0):
    return optimize.brentq(lambda D: gap_equation_residual(D, 0.0, lam, wD), 1e-300, wD)


def Tc(lam, wD=1.0):
    return optimize.brentq(lambda T: gap_equation_residual(1e-14, T, lam, wD), 1e-12 * wD, wD)


def gap_T(T, lam, wD=1.0):
    if T >= Tc(lam, wD):
        return 0.0
    return optimize.brentq(lambda D: gap_equation_residual(D, T, lam, wD), 1e-14, 2 * wD)
