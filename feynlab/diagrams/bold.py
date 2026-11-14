"""Bold-line (skeleton) diagrammatic series for an impurity with Hubbard interaction.

Sigma[G] = U n + sum_{n=2}^{N} (all skeleton diagrams of order n evaluated with the
*dressed* G), solved self-consistently with G^{-1} = G0^{-1} - Sigma[G]
(the scheme of bold diagrammatic Monte Carlo, Prokof'ev & Svistunov PRB 77,
125101 (2008); Van Houcke et al. Nature Phys. 8, 366 (2012)).

Also provides `skeleton_series_at(G)`: the skeleton series evaluated on a given
G — the quantity studied by Kozik, Ferrero & Georges, PRL 114, 156402 (2015),
who showed that for the Hubbard atom at strong coupling the skeleton series on
G_exact converges to an *unphysical* self-energy.
"""
from __future__ import annotations

import numpy as np

from ..core.matsubara import MatsubaraGrid
from .enumerate import unique_diagrams
from .evaluate import G0Table, evaluate_sigma_diagram

_SKEL_CACHE: dict = {}


def skeletons(n):
    if n not in _SKEL_CACHE:
        _SKEL_CACHE[n] = unique_diagrams(n, family="hubbard", kind="skeleton", exclude_tadpoles=True)
    return _SKEL_CACHE[n]


def skeleton_series_at(grid: MatsubaraGrid, G_iw, U: float, max_order: int = 4, n_iw: int = 4,
                       n_samples: int = 400000, seed: int = 3):
    """Order-by-order skeleton contributions Sigma_n[G] (n >= 2) evaluated with G."""
    tab = G0Table(grid, G_iw)
    iw = grid.iw[:n_iw]
    out = {}
    for n in range(2, max_order + 1):
        tot = np.zeros(n_iw, complex)
        var = np.zeros(n_iw)
        for i, d in enumerate(skeletons(n)):
            m, e = evaluate_sigma_diagram(d, tab, tab, U, iw, n_samples, seed + 31 * i + n)
            tot += m
            var += np.abs(e) ** 2
        out[n] = (tot, np.sqrt(var), len(skeletons(n)))
    return out


def bold_selfconsistent(grid: MatsubaraGrid, g0_iw, U: float, max_order: int = 3, n_iw_mc: int = 32,
                        n_iter: int = 15, mix: float = 0.5, n_samples: int = 100000, half_filling=True):
    """Self-consistent skeleton (bold) series truncated at max_order.
    Sigma on the first n_iw_mc frequencies is from MC; above, the 2nd-order bold
    tail is used (exact to O(U^2))."""
    from ..selfenergy.local import second_order_sigma
    iw = grid.iw
    G = g0_iw.copy()
    hist = []
    for it in range(n_iter):
        n = 0.5 if half_filling else grid.density(G)
        s2 = second_order_sigma(grid, G, U)  # bold 2nd order on all frequencies (FFT)
        sig = U * n + s2
        if max_order >= 3:
            ser = skeleton_series_at(grid, G, U, max_order, n_iw_mc, n_samples, seed=17 + it)
            for k in range(3, max_order + 1):
                sig[:n_iw_mc] += ser[k][0]
        if half_filling:
            sig = U / 2 + 1j * sig.imag
        Gn = 1.0 / (1.0 / g0_iw - sig)
        diff = float(np.max(np.abs(Gn - G)))
        G = mix * Gn + (1 - mix) * G
        hist.append(diff)
    return dict(G=G, Sigma=1.0 / g0_iw - 1.0 / G, history=hist)
