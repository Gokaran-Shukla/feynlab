"""One-call decomposition of Hall responses into microscopic mechanisms.

    intrinsic                      Kubo / Berry curvature (clean limit)
    side_jump                      semiclassical, coordinate shift x Born rate
    anomalous_distribution         semiclassical
    skew_V1 (conventional skew)    semiclassical third-order rate  (~ V1^3 / (n_i V0^4))
    skew_gaussian (intrinsic skew) Kubo-Streda NCA total - (intrinsic + side jump + adist),
                                   extrapolated to vanishing disorder broadening
    X, Psi (optional)              crossed impurity-line diagrams (skew scattering on pairs of
                                   close impurities), Ado et al. EPL 111, 37004 (2015)
    total                          sum

For operators other than charge (spin, orbital, quadrupole currents) only the
Kubo-Streda route is available: intrinsic + 'extrinsic_gaussian' (side jump +
Gaussian skew, NCA) + skew_V1.
"""
from __future__ import annotations

import numpy as np

from .intrinsic import intrinsic_hall, k_grid_polar, k_grid_lattice
from .boltzmann import SemiclassicalHall
from .kubo_streda import DisorderedHall


def _default_grid(model, mu, band, kmax, dense_width=0.4):
    if model.kind == "lattice":
        return k_grid_lattice(model, 400)
    from .boltzmann import FermiContour
    c = FermiContour(model, mu, band, 32, kmax)
    kF = float(np.mean(c.k))
    spread = float(np.ptp(c.k))
    w = dense_width + spread
    return k_grid_polar(kmax, nr=300, nth=48, dense=[(kF, w, 1500), (kF, 0.05 + spread, 1500)])


def anomalous_hall_mechanisms(model, mu: float, band: int, V1_over_V0: float = 0.0,
                              n_i: float = 1.0, nV2_list=(0.01, 0.005), kmax: float = 8.0,
                              intrinsic_value: float | None = None, grid=None,
                              n_theta: int = 128, include_crossed: bool = False, crossed_runs=None):
    """Full AHE mechanism decomposition for a charge current (units e^2/h).
    V1_over_V0 = V1/V0 sets the third cumulant; skew_V1 scales as 1/(n_i V0^2)
    and is reported for the *largest-disorder* entry of nV2_list's V0 and for V0 such
    that n_i V0^2 = nV2_list[-1]."""
    grid = grid or _default_grid(model, mu, band, kmax)
    if intrinsic_value is None:
        intrinsic_value = intrinsic_hall(model, mu, grid=grid)
        if hasattr(model, "lower_band_tail") and model.kind == "continuum":
            intrinsic_value += model.lower_band_tail(kmax)
    nV2 = nV2_list[-1]
    V0 = np.sqrt(nV2 / n_i)
    V1 = V1_over_V0 * V0
    sc = SemiclassicalHall(model, mu, band, n_i=n_i, V0=V0, V1=V1, n_theta=n_theta).compute()
    totals, gammas = [], []
    for x in nV2_list:
        dh = DisorderedHall(model, mu, n_i=n_i, V0=np.sqrt(x / n_i), V1=0.0, grid=grid)
        r = dh.decompose(intrinsic_value=intrinsic_value)
        totals.append(r["total"])
        gammas.append(r["Gamma_broadening"])
    if len(totals) >= 2:
        p = np.polyfit(gammas, totals, 1)
        total_nca = float(np.polyval(p, 0.0))
    else:
        total_nca = totals[0]
    skew_g = total_nca - intrinsic_value - sc["side_jump"] - sc["anomalous_distribution"]
    out = dict(intrinsic=intrinsic_value, side_jump=sc["side_jump"],
               anomalous_distribution=sc["anomalous_distribution"],
               skew_gaussian=skew_g, skew_V1=sc["skew_V1"], sigma_xx=sc["sigma_xx"],
               total_gaussian_nca=total_nca)
    out["total"] = total_nca + sc["skew_V1"]
    if include_crossed:
        from .crossed import crossed_hall
        cr = crossed_hall(model, mu, runs=crossed_runs) if crossed_runs else crossed_hall(model, mu)
        out["X"], out["Psi"] = cr["X"], cr["Psi"]
        out["total_with_crossed"] = out["total"] + cr["X"] + cr["Psi"]
    return out


def operator_hall_mechanisms(model, mu: float, op, n_i: float = 1.0, nV2: float = 0.005,
                             V1_over_V0: float = 0.0, grid=None, kmax: float = 8.0,
                             band_for_grid: int = 0, intrinsic_value=None):
    """Spin / orbital / quadrupole Hall conductivity with disorder (units e/2pi):
    intrinsic + extrinsic_gaussian (side jump + Gaussian skew, NCA ladder) + skew_V1."""
    grid = grid or _default_grid(model, mu, band_for_grid, kmax)
    V0 = np.sqrt(nV2 / n_i)
    dh = DisorderedHall(model, mu, n_i=n_i, V0=V0, V1=V1_over_V0 * V0, grid=grid)
    return dh.decompose(op=op, intrinsic_value=intrinsic_value)
