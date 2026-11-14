"""Evaluation of individual Hubbard-family self-energy diagrams for an impurity
(Anderson model / DMFT bath / isolated atom) by Monte Carlo integration over the
internal imaginary times — i.e. a minimal *bare-series diagrammatic Monte Carlo*.

Feynman rules (G = -<T c c^+>, Sigma defined through G = G0 + G0 Sigma G):
    Sigma_D(tau_out - tau_in) = (-U)^n (-1)^F  prod_lines G0_sigma(tau_end - tau_start)
with F the number of closed fermion loops.  G0 must be the Hartree-shifted
propagator (tadpole diagrams excluded).  Checks: n=2 reproduces
Sigma = -U^2 G0(tau)^2 G0(-tau) exactly (IPT kernel).
"""
from __future__ import annotations

import numpy as np

from .enumerate import X, Y, Diagram, unique_diagrams
from ..core.matsubara import MatsubaraGrid


class G0Table:
    def __init__(self, grid: MatsubaraGrid, g0_iw, n_tab: int = 20001):
        self.beta = grid.beta
        tau = np.linspace(0, grid.beta, n_tab)
        c = MatsubaraGrid.tail_coefficients(grid.iw, g0_iw, n_fit=min(60, grid.n_iw // 4))
        self.tab = grid.iw_to_tau(g0_iw, 1.0, c[1], c[2], tau=tau)
        self.n = n_tab

    def __call__(self, t):
        t = np.asarray(t, float)
        sgn = np.where(t < 0, -1.0, 1.0)
        tt = np.where(t < 0, t + self.beta, t)
        x = tt / self.beta * (self.n - 1)
        i = np.clip(x.astype(int), 0, self.n - 2)
        f = x - i
        return sgn * ((1 - f) * self.tab[i] + f * self.tab[i + 1])


def evaluate_sigma_diagram(d: Diagram, g0_up: G0Table, g0_dn: G0Table, U: float,
                           iw: np.ndarray, n_samples: int = 400000, seed: int = 0,
                           batch: int = 100000):
    """Monte Carlo estimate of Sigma_D(iw_m) for the given frequencies.
    Returns (mean, stderr)."""
    rng = np.random.default_rng(seed)
    n = d.order
    beta = g0_up.beta
    first = d.succ[X]
    last = [v for v in d.vertices if d.succ[v] == Y][0]
    line_of = lambda v: v // 2
    l_in, l_out = line_of(first), line_of(last)
    lines = [(a, b) for a, b in d.succ.items() if a != X and b != Y]
    pref = (-U) ** n * d.sign
    w = np.asarray(iw).imag
    acc = []
    done = 0
    while done < n_samples:
        m = min(batch, n_samples - done)
        t = rng.random((m, n)) * beta
        t[:, l_in] = 0.0
        val = np.full(m, pref * beta ** (n - 1))
        for a, b in lines:
            ta = t[:, line_of(a)]
            tb = t[:, line_of(b)]
            g = g0_up if (d.spin is None or d.spin[a] == 0) else g0_dn
            val *= g(tb - ta)
        tau = t[:, l_out]
        acc.append(val[:, None] * np.exp(1j * np.outer(tau, w)))
        done += m
    acc = np.concatenate(acc)
    return acc.mean(0), acc.std(0) / np.sqrt(len(acc))


def sigma_series(grid: MatsubaraGrid, g0_iw, U: float, max_order: int = 4, n_iw: int = 8,
                 n_samples: int = 400000, seed: int = 0):
    """Order-by-order bare perturbation series of the impurity self-energy (beyond
    Hartree).  Returns {order: (Sigma_n(iw), err, n_diagrams)}."""
    tab = G0Table(grid, g0_iw)
    iw = grid.iw[:n_iw]
    out = {}
    for n in range(2, max_order + 1):
        diags = unique_diagrams(n, family="hubbard", kind="proper", exclude_tadpoles=True)
        tot = np.zeros(n_iw, complex)
        var = np.zeros(n_iw)
        for i, d in enumerate(diags):
            m, e = evaluate_sigma_diagram(d, tab, tab, U, iw, n_samples, seed + 101 * i + 7 * n)
            tot += m
            var += np.abs(e) ** 2
        out[n] = (tot, np.sqrt(var), len(diags))
    return out
