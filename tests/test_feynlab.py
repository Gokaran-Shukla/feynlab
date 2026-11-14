"""Fast regression tests (pytest -q, ~1-2 min).  Full validation: benchmarks/run_benchmarks.py"""
import numpy as np
import pytest

from feynlab.core.matsubara import MatsubaraGrid, g0_tau_level
from feynlab.diagrams.enumerate import count_diagrams, unique_diagrams
from feynlab.selfenergy import local
from feynlab.qmc.ed import AndersonED
from feynlab.qmc.ctint import CTINTSolver
from feynlab.dmft.bethe import BetheDMFT
from feynlab.transport.models import MassiveDirac, Rashba2DEG, QWZ
from feynlab.transport.intrinsic import intrinsic_hall, chern_numbers, k_grid_polar
from feynlab.transport.boltzmann import SemiclassicalHall
from feynlab.egas import bcs, rpa


def test_fourier_roundtrip():
    g = MatsubaraGrid(10.0, 256)
    gt = g0_tau_level(0.4, g.tau, g.beta)
    assert np.abs(g.tau_to_iw(gt) - 1 / (g.iw - 0.4)).max() < 1e-4
    assert np.abs(g.iw_to_tau(1 / (g.iw - 0.4), 1, 0.4, 0.16) - gt).max() < 1e-6


@pytest.mark.parametrize("n,G,proper", [(1, 2, 2), (2, 10, 6), (3, 74, 42)])
def test_diagram_counts(n, G, proper):
    assert count_diagrams(n, "generic", "G") == G
    assert count_diagrams(n, "generic", "proper") == proper


def test_hubbard_diagram_counts_and_signs():
    d2 = unique_diagrams(2, "hubbard")
    assert len(d2) == 1 and d2[0].n_loops == 1 and d2[0].sign == -1
    assert count_diagrams(4, "hubbard", "proper", exclude_tadpoles=True) == 12


def test_atomic_sigma_half_filling():
    g = MatsubaraGrid(5.0, 64)
    s = local.atomic_sigma(g.iw, 2.0, 1.0, 5.0)
    assert np.allclose(s, 1.0 + 1.0 / g.iw)


def test_second_order_matches_atom_kernel():
    g = MatsubaraGrid(10.0, 512)
    s = local.second_order_sigma(g, 1 / g.iw, 1.0)
    assert np.allclose(s[:10], 0.25 / g.iw[:10], atol=2e-4)


def test_ed_free_limit():
    ed = AndersonED(0.2, 0.0, [-0.5, 0.7], [0.3, 0.4])
    g = MatsubaraGrid(10.0, 32)
    G = ed.green(g.iw, 10.0)
    assert np.allclose(G, 1 / (g.iw - 0.2 - ed.hybridization(g.iw)), atol=1e-10)


def test_ctint_vs_ed_quick():
    beta, U = 8.0, 2.0
    ed = AndersonED(-U / 2, U, [-1.0, 1.0], [0.6, 0.6])
    g = MatsubaraGrid(beta, 512)
    g0 = 1 / (g.iw - ed.eps_d - ed.hybridization(g.iw))
    r = CTINTSolver(g, U, n_iw_meas=4).solve(g0, n_warm=5000, n_meas=150000, meas_every=10, n_chains=3)
    ex = ed.sigma(r["iw"], beta)
    assert abs(r["Sigma_up"][0].imag - ex[0].imag) < 4 * abs(r["err_Sigma"][0].imag) + 0.02


def test_ipt_metal_and_insulator():
    m = BetheDMFT(beta=50, U=1.0, n_iw=1024).solve("ipt")
    i = BetheDMFT(beta=50, U=4.0, n_iw=1024).solve("ipt", insulating_start=True)
    assert 0.6 < m["Z"] < 0.9 and i["Z"] < 0.01


def test_chern_and_rashba_shc():
    assert np.allclose(chern_numbers(QWZ(-1.0), 40), [-1, 1], atol=1e-8) or np.allclose(chern_numbers(QWZ(-1.0), 40), [1, -1], atol=1e-8)
    r = Rashba2DEG(1.0, 0.3)
    gg = k_grid_polar(8.0, nr=3000, nth=8, dense=[(1.0, 0.5, 3000)])
    assert abs(intrinsic_hall(r, 0.5, op=r.operators["sz/2"], grid=gg) - 0.25) < 3e-3


def test_side_jump_and_skew_dirac():
    d = MassiveDirac(1.0, 0.5)
    r = SemiclassicalHall(d, 1.0, band=1, n_i=1.0, V0=0.1, V1=0.05, n_theta=64).compute()
    assert abs(r["side_jump"] - d.ref_side_jump(1.0)) < 1e-4
    assert abs(r["anomalous_distribution"] - d.ref_anomalous_distribution(1.0)) < 1e-4
    assert abs(r["skew_V1"] - d.ref_skew_V1(1.0, 1.0, 0.1, 0.05)) < 1e-4


def test_bcs_ratio():
    assert abs(bcs.gap_T0(0.2) / bcs.Tc(0.2) - 1.7639) < 3e-3


def test_plasmon_three_fifths():
    rs = 2.0
    kF = rpa.kF_of_rs(rs)
    q = np.array([0.02, 0.04]) * kF
    w = rpa.plasmon_dispersion(rs, q)
    assert abs(np.polyfit(q**2, w**2 - rpa.plasma_frequency(rs) ** 2, 1)[0] / kF**2 - 0.6) < 0.01


# ---------------------------------------------------------------- v0.2
def test_sector_ed_matches_dense():
    from feynlab.qmc.ed_sector import SectorED
    a = AndersonED(-1.2, 3.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
    b = SectorED(-1.2, 3.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
    z = 1j * np.pi * (2 * np.arange(4) + 1) / 10
    assert np.abs(a.green(z, 10.0) - b.green(z, 10.0)).max() < 1e-10


def test_lattice_dmft_equals_bethe():
    from feynlab.dmft.bethe import LatticeDMFT
    a = BetheDMFT(50, 1.5).solve("ipt")
    b = LatticeDMFT.semicircle(50, 1.5).solve("ipt")
    assert abs(a["Z"] - b["Z"]) < 1e-4


def test_feynman_polaron_energy():
    from feynlab.qmc.polaron import feynman_energy
    assert abs(feynman_energy(1.0) + 1.0130308) < 2e-6


def test_polaron_diagmc_weak_coupling():
    from feynlab.qmc.polaron import polaron_diagmc, weak_coupling_energy
    r = polaron_diagmc(0.3, mu=-0.45, tau_max=30, n_steps=4_000_000, n_chains=3, fit_window=(0.4, 0.97))
    assert abs(r["E"] - weak_coupling_energy(0.3)) < 4 * r["E_err"] + 3e-3


def test_crossed_psi_small_and_x_positive():
    from feynlab.transport.crossed import CrossedDiagrams
    r = CrossedDiagrams(MassiveDirac(1.0, 0.5), 1.0, 0.06, L=300.0, N=512).compute()
    assert 0.3 < r["X"] < 0.6 and abs(r["Psi"]) < 0.15


def test_skeleton_series_weak_coupling():
    from feynlab.diagrams.bold import skeleton_series_at
    g = MatsubaraGrid(1.0, 256)
    U = 0.5
    Gex = 0.5 * (1 / (g.iw + U / 2) + 1 / (g.iw - U / 2))
    ser = skeleton_series_at(g, Gex, U, 3, n_iw=1, n_samples=50000)
    tot = local.second_order_sigma(g, Gex, U)[0].imag + ser[3][0][0].imag
    assert abs(tot / (U**2 / (4 * g.iw[0])).imag - 1) < 0.03


def test_maxent_normalisation():
    from feynlab.core.maxent import maxent
    from feynlab.core.lattice import bethe_dos
    w = np.linspace(-3, 3, 151)
    iw = 1j * np.pi * (2 * np.arange(60) + 1) / 20
    A0 = bethe_dos(w)
    G = (A0[None] * (w[1] - w[0]) / (iw[:, None] - w[None])).sum(1)
    r = maxent(iw, G, 1e-3, w, alphas=np.logspace(2, -2, 12))
    assert abs(np.sum(r["A"]) * (w[1] - w[0]) - 1) < 0.02


# ---------------------------------------------------------------- v0.3: T=0 Lanczos ED-DMFT
def test_lanczos_matches_full_ed():
    from feynlab.qmc.ed_sector import SectorED, LanczosED
    args = (-1.5, 3.0, [-1.0, -0.2, 0.2, 1.0], [0.5, 0.3, 0.3, 0.5])
    z = 1j * np.pi * (2 * np.arange(5) + 1) / 200
    assert np.abs(SectorED(*args).green(z, 1e5) - LanczosED(*args).green(z)).max() < 1e-8


def test_lanczos_ground_state_parity():
    """odd n_bath -> even sites -> singlet in (h,h); even n_bath -> doublet (h,h+1),(h+1,h)."""
    from feynlab.qmc.ed_sector import LanczosED
    odd = LanczosED(-1.5, 3.0, [-1, -0.3, 0.0, 0.3, 1], [0.4, 0.3, 0.2, 0.3, 0.4])
    even = LanczosED(-1.5, 3.0, [-1, -0.3, 0.3, 1], [0.4, 0.3, 0.3, 0.4])
    assert len(odd.gs) == 1 and odd.gs[0][0] == (3, 3)
    assert sorted(s for s, _ in even.gs) == [(2, 3), (3, 2)]


def test_ed_dmft_metal_7bath_quick():
    d = BetheDMFT(beta=200, U=2.0, n_iw=512)
    r = d.solve("ed", n_iter=40, mix=0.5, tol=1e-5, ctint_kwargs=dict(n_bath=7))
    assert 0.25 < r["Z"] < 0.31   # 7-bath T=0 value 0.279; NRG ~0.28-0.30
