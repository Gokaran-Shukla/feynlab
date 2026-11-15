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


# ============================================================================ v1.0 (ROADMAP modules)
def test_gamma_method_and_jackknife():
    from feynlab.utils.errors import ar1_series, gamma_method, jackknife
    _, _, tau, _ = gamma_method(ar1_series(200000, 0.8))
    assert abs(tau - 4.5) / 4.5 < 0.1
    rng = np.random.default_rng(0)
    m, e = jackknife(lambda a, b: a / b, 1 + 0.1 * rng.standard_normal(5000), 2 + 0.1 * rng.standard_normal(5000))[:2]
    assert abs(m - 0.5) < 5 * e + 1e-3


def test_mixers_reach_same_fixed_point():
    from feynlab.utils.mixing import LinearMixer, AndersonMixer, BroydenMixer, solve_fixed_point
    F = lambda x: 0.9 * np.cos(x) + 0.1 * np.roll(x, 1)
    sols = [solve_fixed_point(F, np.zeros(10), mx)[0] for mx in (LinearMixer(0.5), AndersonMixer(), BroydenMixer())]
    assert np.abs(sols[0] - sols[1]).max() < 1e-7 and np.abs(sols[0] - sols[2]).max() < 1e-7


def test_hdf5_checkpoint(tmp_path):
    from feynlab.utils.io import Checkpoint
    ck = Checkpoint(str(tmp_path / "c.h5"), {"U": 2.0})
    ck.save(1, sigma=np.ones(3) * 1j)
    ck.save(2, sigma=np.ones(3) * 2j)
    it, d = ck.latest()
    assert it == 2 and np.allclose(d["sigma"], 2j)


def test_cthyb_vs_ed_quick():
    from feynlab.qmc.cthyb import CTHYBSolver
    beta = 10.0
    g = MatsubaraGrid(beta, 512)
    ed = AndersonED(-1.2, 3.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
    r = CTHYBSolver(g, 3.0, n_iw_meas=4).solve(ed.hybridization(g.iw), -ed.eps_d, n_meas=200_000, n_chains=2)
    assert abs(np.mean(r["n"]) - ed.occupation(beta)) < 0.02
    assert abs(r["docc"] - ed.double_occupancy(beta)) < 0.01


def test_kanamori_multiplets():
    from feynlab.qmc.multiorbital import KanamoriAtom
    m = KanamoriAtom(2, 4.0, 0.6).multiplets(2)
    assert abs(m[0] - 2.2) < 1e-10 and np.sum(np.abs(m - 2.2) < 1e-9) == 3


def test_dca_nc1_equals_dmft():
    from feynlab.dmft.dca import DCA2D
    from feynlab.dmft.bethe import LatticeDMFT
    r = DCA2D(Lc=1, nk=32, U=1.5, beta=10.0, n_iw=256).solve()
    k = 2 * np.pi * np.arange(32) / 32
    kx, ky = np.meshgrid(k, k, indexing="ij")
    ref = LatticeDMFT.from_eps_k(10.0, 1.5, -0.5 * (np.cos(kx) + np.cos(ky)), n_iw=256).solve("ipt", tol=1e-8)
    assert abs(r["Z"][0] - ref["Z"]) < 1e-4


def test_compressibility_sum_rule():
    from feynlab.twoparticle import bubble
    from feynlab.core.matsubara import fermi
    from feynlab.core.lattice import square_lattice_eps
    g = MatsubaraGrid(5.0, 256)
    eps = square_lattice_eps(16)
    chi = bubble(eps, np.zeros(256), g, 0.3)
    f = fermi(eps - 0.3, 5.0)
    assert abs(chi[0, 0, 0].real - np.mean(5.0 * f * (1 - f))) < 1e-6


def test_nevanlinna_interpolates():
    from feynlab.continuation.nevanlinna import Nevanlinna
    w = np.linspace(-4, 4, 401)
    dw = w[1] - w[0]
    A = np.exp(-w**2 / 0.5) / np.sqrt(0.5 * np.pi)
    iw = 1j * np.pi * (2 * np.arange(12) + 1) / 30
    G = (A[None] * dw / (iw[:, None] - w[None])).sum(1)
    assert np.abs(Nevanlinna(iw, G)(iw) - G).max() < 1e-10


def test_resummation():
    from scipy.special import binom
    from feynlab.diagrams.resum import domb_sykes, pade, pade_eval, conformal_map_sum
    a = np.array([binom(-0.5, k) for k in range(20)])
    assert abs(domb_sykes(a)[0] - 1) < 1e-6
    p, q = pade(a, 8, 8)
    assert abs(pade_eval(p, q, 3.0).real - 0.5) < 1e-6
    assert abs(conformal_map_sum(a, 1.0, 3.0) - 0.5) < 1e-6


def test_cdet_atom():
    from feynlab.qmc.cdet import cdet_series, atom_exact_coefficients
    g = MatsubaraGrid(5.0, 1024)
    c, e = cdet_series(g, 1 / g.iw, max_order=3, n_samples=50000)
    assert np.allclose(c, atom_exact_coefficients(5.0, 3), atol=1e-6 + 5 * np.max(e))


def test_nrg_resonant_level():
    from feynlab.qmc.nrg import NRG, resonant_level_exact
    a = NRG(0.0, None, 0.01, 2.5, 20, 200, 2.0, True).run()
    b = NRG(0.0, None, 0.01, 2.5, 20, 200, 2.0, False).run()
    S = 0.5 * ((a[1:, 1] - b[1:, 1]) + (a[:-1, 1] - b[:-1, 1]))
    T = np.sqrt(a[1:, 0] * a[:-1, 0])
    Se, _ = resonant_level_exact(T[[12]], 0.01)
    assert abs(S[12] / Se[0] - 1) < 0.08


def test_orbital_moment_and_bcd():
    from feynlab.transport.models import MassiveDirac, TiltedDirac
    from feynlab.transport.response import orbital_moment, dirac_valley_moment, berry_curvature_dipole
    from feynlab.transport.intrinsic import k_grid_polar
    k = np.array([0.0, 0.4])
    assert np.allclose(orbital_moment(MassiveDirac(1, 0.3), k, 0 * k)[1][:, 0], dirac_valley_moment(1, 0.3, k))
    g = k_grid_polar(2.0, nr=200, nth=128, dense=[(0.283, 0.2, 800)])
    a = berry_curvature_dipole(TiltedDirac(1, 0.1, 0.3), 0.3, g, beta=300)[0]
    b = berry_curvature_dipole(TiltedDirac(1, 0.1, 0.3), 0.3, g, beta=300, form="sea")[0]
    assert abs(a - b) < 1e-4 * abs(a) and abs(a) > 1e-3


def test_rashba_edelstein():
    from feynlab.transport.models import Rashba2DEG
    from feynlab.transport.response import edelstein_boltzmann, rashba_edelstein_reference
    from feynlab.transport.intrinsic import k_grid_polar
    g = k_grid_polar(4.0, nr=200, nth=64, dense=[(np.sqrt(2.0), 0.5, 1500)])
    s = edelstein_boltzmann(Rashba2DEG(1.0, 0.3, 0.0), 1.0, np.array([[0, -1j], [1j, 0]]) / 2, g)
    assert abs(abs(s) / rashba_edelstein_reference(1.0, 0.3) - 1) < 0.01


def test_tmatrix_identity_limit():
    from feynlab.transport.models import MassiveDirac
    from feynlab.transport.kubo_streda import DisorderedHall
    from feynlab.transport.intrinsic import k_grid_polar
    d = MassiveDirac(1.0, 0.5)
    g = k_grid_polar(6.0, nr=150, nth=32, dense=[(d.kF(1.0), 0.3, 600)])
    a = DisorderedHall(d, 1.0, n_i=1.0, V0=0.07, grid=g).sigma_I()
    b = DisorderedHall(d, 1.0, n_i=1.0, V0=0.07, grid=g, imp_matrix=np.eye(2)).sigma_I()
    assert abs(a - b) < 1e-12


def test_wannier3d_quantised():
    from feynlab.transport.wannier3d import WannierTB3D, stacked_qwz_hr
    R, H = stacked_qwz_hr(1.0, 0.2)
    w = WannierTB3D(Rs=R, HR=H, lat=np.diag([1.0, 1.0, 3.0]))
    assert abs(w.intrinsic_hall(0.0, nk=(16, 16, 4), adaptive=False)["sigma"][2] - 1 / 3) < 2e-3


def test_cli_parser():
    from feynlab.cli.main import build_parser
    p = build_parser()
    assert p.parse_args(["hall", "--model", "dirac"]).command == "hall"


# ============================================================================ v1.1 (closing the v1.0 limitations)
def test_polaron_mstar_estimator():
    from feynlab.qmc.polaron import polaron_diagmc, feynman_energy, weak_coupling_mass
    r = polaron_diagmc(0.2, mu=feynman_energy(0.2) - 0.15, tau_max=40, n_steps=8_000_000, n_chains=2,
                       fit_window=(0.4, 0.97))
    assert abs(r["m_star"] - weak_coupling_mass(0.2)) < 0.01


def test_sac_beats_threshold():
    from feynlab.continuation.stochastic import stochastic_continuation
    w = np.linspace(-4, 4, 801); dw = w[1] - w[0]
    A0 = 0.5 * (np.exp(-(w - 1.5) ** 2 / 0.18) + np.exp(-(w + 1.5) ** 2 / 0.18)) / np.sqrt(0.18 * np.pi)
    iw = 1j * np.pi * (2 * np.arange(30) + 1) / 50
    G = (A0[None] * dw / (iw[:, None] - w[None])).sum(1) + 1e-4 * np.random.default_rng(0).standard_normal(30)
    wc = np.linspace(-4, 4, 161)
    r = stochastic_continuation(iw, G, 1e-4, wc, n_sweeps=150, n_final=800)
    assert abs(np.sum(r["A"]) * (wc[1] - wc[0]) - 1) < 0.02
    assert np.sum(np.abs(r["A"] - np.interp(wc, w, A0))) * (wc[1] - wc[0]) < 0.3


def test_finite_range_sigma_xx():
    from feynlab.transport.models import MassiveDirac
    from feynlab.transport.boltzmann import SemiclassicalHall, dirac_finite_range_sigma_xx
    r = SemiclassicalHall(MassiveDirac(1.0, 0.5), 1.0, 1, n_i=1.0, V0=0.1, n_theta=64, range_d=1.0).compute()
    assert abs(r["sigma_xx"] / dirac_finite_range_sigma_xx(1, 0.5, 1.0, 1.0, 0.1, 1.0) - 1) < 1e-6


def test_semiclassical_gaussian_skew():
    from feynlab.transport.models import MassiveDirac
    from feynlab.transport.boltzmann import SemiclassicalHall
    d = MassiveDirac(1.0, 0.5)
    r = SemiclassicalHall(d, 1.0, 1, n_i=1.0, V0=0.1, n_theta=64).compute()
    assert abs(r["skew_gaussian"] - d.ref_intrinsic_skew(1.0)) < 1e-6


def test_generic_lattice_crossed_matches_square():
    from feynlab.transport.models import QWZ, WannierTB
    from feynlab.transport.crossed import CrossedDiagrams, LatticeCrossed
    from feynlab.transport.wannier3d import stacked_qwz_hr
    from feynlab.dmft.dftdmft import write_hr
    import tempfile, os
    a = CrossedDiagrams(QWZ(-1.7), 0.6, 0.03, N=128, lattice=True).compute()
    R, H = stacked_qwz_hr(-1.7, 0.0); m = R[:, 2] == 0
    f = os.path.join(tempfile.mkdtemp(), "q_hr.dat"); write_hr(f, R[m], H[m])
    b = LatticeCrossed(WannierTB(f, np.eye(3)), 0.6, 0.03, N=(128, 128)).compute()
    assert abs(a["X"] - b["X"]) < 1e-4 and abs(a["Psi"] - b["Psi"]) < 1e-4


def test_matrix_cthyb_single_orbital():
    from feynlab.qmc.cthyb_matrix import MatrixCTHYB
    g = MatsubaraGrid(5.0, 256)
    ed = AndersonED(-1.2, 3.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
    r = MatrixCTHYB(g, 1, 3.0, 0.0, eps=[-1.2], n_iw_meas=2).solve(ed.hybridization(g.iw), n_meas=150000, n_chains=2)
    assert r["sign"] > 0.99
    assert abs(r["n_0"] - ed.occupation(5.0)) < 0.01 and abs(r["docc_0"] - ed.double_occupancy(5.0)) < 0.01


def test_cluster_ctint_nc1_matches_single_site():
    from feynlab.qmc.ctint_cluster import ClusterCTINT
    from feynlab.qmc.ctint import CTINTSolver
    g = MatsubaraGrid(5.0, 256)
    ed = AndersonED(-1.0, 2.0, [-0.8, 0.8], [0.5, 0.5])
    G0 = 1.0 / (g.iw + 1.0 - ed.hybridization(g.iw))
    a = ClusterCTINT(g, 2.0, sites=[[0]], L=[1], K=[[0.0]], n_iw_meas=2).solve(G0[None], n_meas=100000, n_chains=2)
    ex = ed.sigma(g.iw[:2], 5.0)
    assert abs(a["Sigma_K"][0, 0] - ex[0]) < 0.03


def test_dftdmft_density_matrix_and_csc():
    from feynlab.dmft.dftdmft import DFTDMFT, write_hr, hartree_functional, hartree_fock_reference, read_hr, hk_on_mesh
    from feynlab.core.matsubara import fermi
    import tempfile, os
    f = os.path.join(tempfile.mkdtemp(), "c_hr.dat")
    Rs = np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)])
    HR = np.array([np.diag([0.0, 0.6]) + 0j] + [(-0.25 * np.diag([1.0, 0.6]) - 0.05 * np.array([[0, 1], [1, 0]])) + 0j] * 4)
    write_hr(f, Rs, HR)
    m = DFTDMFT(f, corr=[0], U=0.0, J=0.0, n_elec=1.4, beta=40, nk=(16, 16, 1))
    N = m.density_matrix_k(np.zeros((1, m.grid.n_iw), complex), 0.1)
    E, V = np.linalg.eigh(m.Hk)
    assert np.abs(N - np.einsum("kan,kn,kbn->kab", V, fermi(E - 0.1, 40), V.conj())).max() < 1e-5
    r = m.solve_csc(hartree_functional(1.5), n_outer=60, mix_rho=0.5, dmft_iter=5, tol_rho=1e-8)
    ref, _ = hartree_fock_reference(hk_on_mesh(*read_hr(f), (16, 16, 1)), 1.5, 1.4, 40)
    assert np.abs(r["n_orbital"] - ref).max() < 1e-4


def test_two_particle_ed_wick_and_eom():
    from feynlab.vertex import ImpurityED, vertices, eom_sigma
    beta = 5.0
    ed0 = ImpurityED(-0.2, 0.0, [0.5], [0.6])
    gs = ed0.generalized_susceptibilities(beta, n_nu=3, n_om=1)
    assert np.abs(gs["chi_uu"][1] - np.diag(gs["chi0"][1])).max() < 1e-12 and np.abs(gs["chi_ud"]).max() < 1e-12
    ed = ImpurityED(-1.0, 2.0, [0.0], [0.5])
    gs = ed.generalized_susceptibilities(beta, n_nu=8, n_om=8)
    S = eom_sigma(vertices(gs, beta)["F_ud"], gs, beta, 2.0, 2 * ed.occupation(beta))
    assert abs(S[8] - ed.sigma(1j * gs["nu"][8:9], beta)[0]) < 5e-3


def test_tpsc_sum_rules():
    from feynlab.vertex import tpsc
    from feynlab.core.lattice import square_lattice_eps
    r = tpsc(square_lattice_eps(16), 5.0, 2.0, 1.0, n_iw=128)
    assert abs(r["sum_sp"] - r["target_sp"]) < 1e-6 and abs(r["sum_ch"] - r["target_ch"]) < 1e-6
    assert r["U_sp"] < 2.0 < r["U_ch"]


def test_holstein_ved():
    from feynlab.qmc.lattice_polaron import LatticePolaronVED, holstein_weak_coupling
    r = LatticePolaronVED("holstein", 1.0, 1.0, 1.0, n_gen=12).ground_state()
    assert abs(r["E"] + 2.469684723932870) < 1e-8
    w = LatticePolaronVED("holstein", 1.0, 1.0, 0.1, n_gen=8).ground_state()
    assert abs(w["E"] - holstein_weak_coupling(1.0, 1.0, 0.1)) < 1e-5


def test_parallel_backends_identical():
    from feynlab.qmc.cthyb import CTHYBSolver
    g = MatsubaraGrid(5.0, 256)
    ed = AndersonED(-1.0, 2.0, [-0.5, 0.5], [0.5, 0.5])
    s = CTHYBSolver(g, 2.0, n_iw_meas=2)
    a = s.solve(ed.hybridization(g.iw), 1.0, n_meas=40000, n_chains=2, parallel="serial")
    b = s.solve(ed.hybridization(g.iw), 1.0, n_meas=40000, n_chains=2, parallel="process")
    assert abs(a["docc"] - b["docc"]) < 1e-12


# ============================================================================ v1.2: QE charge self-consistency
def test_qe_dmft_h5_layout(tmp_path):
    """The archive must match QE's add_dmft_occ.f90: Fortran delta_n(2, nb, nb, nk) <-> h5py (nk, nb, nb, 2) with
    D[k, i, j] = Delta N_ij (QE transposes on read), band_window (1, nk, 2)."""
    import h5py
    from feynlab.dft.qe_csc import write_dmft_h5
    rng = np.random.default_rng(0)
    dN = rng.standard_normal((5, 3, 3)) + 1j * rng.standard_normal((5, 3, 3))
    p = str(tmp_path / "x_dmft.h5")
    write_dmft_h5(p, (21, 23), dN)
    with h5py.File(p) as h:
        bw = h["dft_misc_input/band_window"][...]
        d = h["dft_update/delta_N"][...]
    assert bw.shape == (1, 5, 2) and bw[0, 0, 0] == 21 and bw[0, 0, 1] == 23
    F = np.transpose(d, (3, 2, 1, 0))                     # what Fortran sees: delta_n(re/im, i, j, k)
    n_dmft = F[0] + 1j * F[1]                              # n_dmft(i, j, k)
    assert np.allclose(np.transpose(n_dmft, (1, 0, 2)), np.transpose(dN, (1, 2, 0)))   # after QE's TRANSPOSE


def test_qe_window_density_noninteracting():
    """Sigma = 0: N^DMFT(k) = diag(f(eps_k - mu)), so Delta N = 0 when QE uses Fermi-Dirac smearing at the same T."""
    from feynlab.dft.qe_csc import window_density_matrix, find_mu_window
    from feynlab.core.matsubara import fermi
    rng = np.random.default_rng(1)
    eps = np.sort(rng.uniform(-1, 1, (40, 3)), axis=1)
    V = np.array([np.linalg.qr(rng.standard_normal((3, 3)) + 1j * rng.standard_normal((3, 3)))[0] for _ in range(40)])
    sig = np.zeros((3, 512), complex)
    mu = find_mu_window(eps, V, sig, 40.0, 512, 1.0)
    N = window_density_matrix(eps, V, sig, mu, 40.0, 512)
    f = fermi(eps - mu, 40.0)
    assert np.abs(N - np.einsum("km,mn->kmn", f, np.eye(3))).max() < 1e-5
    assert abs(2 * np.real(np.trace(N, axis1=1, axis2=2)).mean() - 1.0) < 1e-6


def test_w90_u_matrix_reader(tmp_path):
    from feynlab.dft.qe_csc import read_u_matrix
    U = np.array([[[1 + 2j, 3 + 4j], [5 + 6j, 7 + 8j]]])
    lines = ["header", "1 2 2", "", " 0.0 0.0 0.0"]
    for n in range(2):                 # Wannier90 writes ((u(i, n), i = 1, num_bands), n = 1, num_wann)
        for i in range(2):
            lines.append(f"{U[0, i, n].real} {U[0, i, n].imag}")
    p = tmp_path / "s_u.mat"
    p.write_text("\n".join(lines) + "\n")
    assert np.allclose(read_u_matrix(str(p)), U)


def test_dftdmft_bath_levels_measured_from_mu(tmp_path):
    """Regression (v1.2): fitted bath levels are already measured from mu; shifting them again by -mu emptied the
    bath away from mu = 0.  At n = 0.8 the ED and Lanczos solvers must agree with each other and with matrix
    CT-HYB (no bath fit, Z = 0.497 +- 0.005 from benchmarks/qe_csc)."""
    from feynlab.dmft.dftdmft import DFTDMFT, write_hr
    f = str(tmp_path / "one_hr.dat")
    write_hr(f, np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]), np.array([[[0j]]] + [[[-0.25 + 0j]]] * 4))
    Z = []
    for solver in ("ed", "lanczos"):
        m = DFTDMFT(f, U=2.0, J=0.0, n_elec=0.8, beta=10.0, nk=(16, 16, 1), n_iw=256)
        Z.append(m.solve(solver, n_iter=40, mix=0.5, n_bath=3, tol=1e-6)["Z"][0])
    assert abs(Z[0] - Z[1]) < 2e-3
    assert abs(Z[0] - 0.497) < 0.02


def test_multiorbital_lanczos_matches_full_ed():
    from feynlab.qmc.multiorbital_lanczos import MultiOrbitalLanczos
    from feynlab.qmc.multiorbital import MultiOrbitalED
    U, J = 2.5, 0.4
    e = [-(3 * U - 5 * J) / 2 + 0.2, -(3 * U - 5 * J) / 2 + 0.5]
    eb, Vb = [[-0.5, 0.6], [-0.4, 0.7]], [[0.5, 0.4], [0.45, 0.5]]
    z = 1j * np.pi * (2 * np.arange(3) + 1) / 20.0
    a = MultiOrbitalED(2, U, J, e, eb, Vb).green(z, 20.0, a=0)
    b = MultiOrbitalLanczos(2, U, J, e, eb, Vb, beta=20.0, n_eig=40).green(z, a=0)
    assert np.abs(a - b).max() < 1e-5


def _rpa_vertex_setup(U=0.6, beta=4.0, nn=6, N=8):
    from feynlab.core.lattice import square_lattice_eps
    eps = square_lattice_eps(N, 0.25)
    nus = (2 * np.arange(-nn, nn) + 1) * np.pi / beta
    oms = 2 * np.arange(-nn + 1, nn) * np.pi / beta
    nO, nF = len(oms), len(nus)
    v = dict(F_d=np.full((nO, nF, nF), U + 0j), F_m=np.full((nO, nF, nF), -U + 0j),
             Gamma_d=np.full((nO, nF, nF), U + 0j), Gamma_m=np.full((nO, nF, nF), -U + 0j))
    sd = lambda x: np.full(len(x), U / 2 + 0j)
    Gl = lambda x: np.mean(1 / (1j * x[None, None, :] - eps[..., None]), axis=(0, 1))
    return eps, v, dict(nu=nus, om=oms, G_of=Gl), sd


def test_lambda_dga_reduces_to_ladder_and_fulfils_sum_rule():
    from feynlab.vertex import ladder_dga_sigma, lambda_dga
    U, beta = 0.6, 4.0
    eps, v, gs, sd = _rpa_vertex_setup(U, beta)
    S0 = ladder_dga_sigma(eps, U / 2, sd, v, gs, beta, U, 1.0)
    r0 = lambda_dga(eps, U / 2, sd, v, gs, beta, U, 1.0, lambda_mode=None)
    assert np.abs(r0["Sigma"] - S0).max() < 1e-12
    # gamma = 1 for a constant (RPA) vertex
    target = 0.8 * r0["sum_lattice"]
    r = lambda_dga(eps, U / 2, sd, v, gs, beta, U, 1.0, lambda_mode="m", target_m=target)
    assert r["lam_m"] > 0 and abs(r["sum_lattice"] - target) < 1e-10
    w0 = len(gs["om"]) // 2
    assert np.allclose(r["chi_m_lambda"][:, w0], 1 / (1 / r["chi_m"][:, w0] + r["lam_m"]))
    rf = lambda_dga(eps, U / 2, sd, v, gs, beta, U, 1.0, lambda_mode="fixed", lam_m=r["lam_m"])
    assert np.abs(rf["Sigma"] - r["Sigma"]).max() < 1e-12
