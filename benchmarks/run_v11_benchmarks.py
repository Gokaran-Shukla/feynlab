"""feynlab v1.1 benchmark suite: the checks that close the v1.0 limitations.
Usage:   python benchmarks/run_v11_benchmarks.py
Writes benchmarks/RESULTS_v11.md and benchmarks/results_v11.json.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
results = []


def record(group, name, value, reference, tol, source, kind="abs"):
    value, reference = float(np.real(value)), float(reference)
    err = abs(value - reference) if kind == "abs" else abs(value - reference) / max(abs(reference), 1e-300)
    ok = bool(err <= tol) if kind != "info" else None
    results.append(dict(group=group, name=name, value=value, reference=reference, error=err, tol=tol, kind=kind,
                        passed=ok, source=source))
    print(f"[{'INFO' if ok is None else 'PASS' if ok else 'FAIL'}] {group:11s} {name:62s} {value: .8g}  ref {reference: .8g}  err {err:.2e}",
          flush=True)


t0 = time.time()
from feynlab.core.matsubara import MatsubaraGrid  # noqa: E402
from feynlab.qmc.ed import AndersonED  # noqa: E402

# ------------------------------------------------------------------ 1. full-Kanamori CT-HYB (matrix trace)
from feynlab.qmc.cthyb_matrix import MatrixCTHYB  # noqa: E402
from feynlab.qmc.multiorbital import MultiOrbitalED  # noqa: E402

beta = 5.0
g = MatsubaraGrid(beta, 512)
U, J = 2.5, 0.4
lvl = -(3 * U - 5 * J) / 2 + 0.2
ed = MultiOrbitalED(2, U, J, [lvl, lvl + 0.3], [[-0.5, 0.6], [-0.4, 0.7]], [[0.5, 0.4], [0.45, 0.5]])
D = np.array([ed.hybridization(g.iw, a) for a in (0, 0, 1, 1)])
r = MatrixCTHYB(g, 2, U, J, eps=[lvl, lvl + 0.3], n_iw_meas=3).solve(D, n_meas=1_500_000, n_chains=8, meas_every=20)
src = "Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED"
n = lambda a, s: ed.n_op(a, s)
record("CT-HYB-mat", "2-orbital Kanamori (spin flip + pair hopping): average sign", r["sign"], 1.0, 1e-12, src)
for a in (0, 1):
    record("CT-HYB-mat", f"orbital {a}: n_up vs ED", r[f"n_{2 * a}"], ed.expectation(n(a, 0), beta), max(4 * r[f"n_{2 * a}_err"], 2e-3), src)
    record("CT-HYB-mat", f"orbital {a}: <n_up n_dn> vs ED", r[f"docc_{a}"], ed.expectation(n(a, 0) @ n(a, 1), beta),
           max(4 * r[f"docc_{a}_err"], 2e-3), src)
Sz = 0.5 * sum(n(a, 0) - n(a, 1) for a in (0, 1))
record("CT-HYB-mat", "<S_z^2> (Hund's moment) vs ED", r["Sz^2"], ed.expectation(Sz @ Sz, beta), max(4 * r["Sz^2_err"], 2e-3), src)
Ge = ed.green(g.iw[:1], beta, a=0)
record("CT-HYB-mat", "Im G(i w_0), orbital 0, vs ED", r["G"][0][0].imag, Ge[0].imag, 0.01, src)

# ------------------------------------------------------------------ 2. cluster CT-INT vs ED (dimer in a bath)
from feynlab.qmc.ctint_cluster import ClusterCTINT  # noqa: E402
from feynlab.qmc.multiorbital import _jw_ops  # noqa: E402

Uc, tc, muc = 2.0, -0.3, 1.15
eb, Vb = np.array([0.4, -0.6]), np.array([0.5, 0.6])
c = [x.toarray() for x in _jw_ops(8)]
cd = lambda i, s: c[2 * i + s].T
cc = lambda i, s: c[2 * i + s]
H = np.zeros((256, 256))
for s in (0, 1):
    for i in (0, 1):
        H += -muc * cd(i, s) @ cc(i, s)
    H += tc * (cd(0, s) @ cc(1, s) + cd(1, s) @ cc(0, s))
    for K, sg in enumerate([1, -1]):
        cK = (cc(0, s) + sg * cc(1, s)) / np.sqrt(2)
        H += eb[K] * cd(2 + K, s) @ cc(2 + K, s) + Vb[K] * (cK.T @ cc(2 + K, s) + cc(2 + K, s).T @ cK)
for i in (0, 1):
    H += Uc * cd(i, 0) @ cc(i, 0) @ cd(i, 1) @ cc(i, 1)
E, V = np.linalg.eigh(H); E -= E.min(); wB = np.exp(-beta * E)


def GKed(sg, z):
    cK = V.T @ ((cc(0, 0) + sg * cc(1, 0)) / np.sqrt(2)) @ V
    A = np.abs(cK) ** 2 * (wB[None, :] + wB[:, None]); dE = E[None, :] - E[:, None]
    return np.array([np.sum(A / (zz - dE)) for zz in z]) / wB.sum()


epsK = np.array([tc, -tc])
G0 = np.array([1 / (g.iw + muc - epsK[K] - Vb[K] ** 2 / (g.iw - eb[K])) for K in range(2)])
rc = ClusterCTINT(g, Uc, sites=[[0], [1]], L=[2], K=[[0.0], [np.pi]], n_iw_meas=3).solve(G0, n_meas=1_500_000, n_chains=8)
for K, sg in enumerate([1, -1]):
    se = 1 / G0[K, 0] - 1 / GKed(sg, g.iw[:1])[0]
    record("ClusterCTINT", f"dimer+bath, K={'0' if K == 0 else 'pi'}: Im Sigma(i w_0) vs ED", rc["Sigma_K"][K, 0].imag, se.imag,
           max(4 * abs(rc["err_Sigma"][K, 0].imag), 3e-3), "Gull et al. PRB 83, 075122 (2011) vs ED")
    record("ClusterCTINT", f"dimer+bath, K={'0' if K == 0 else 'pi'}: Re Sigma(i w_0) vs ED", rc["Sigma_K"][K, 0].real, se.real,
           max(4 * abs(rc["err_Sigma"][K, 0].real), 3e-3), "vs ED")

# ------------------------------------------------------------------ 3. DFT+DMFT: density matrix, CSC, matrix CT-HYB solver
from feynlab.dmft.dftdmft import DFTDMFT, write_hr, read_hr, hk_on_mesh, hartree_functional, hartree_fock_reference  # noqa: E402
from feynlab.core.matsubara import fermi  # noqa: E402

f = "/tmp/feynlab_csc_hr.dat"
Rs = np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)])
HR = np.array([np.diag([0.0, 0.6]) + 0j] + [(-0.25 * np.diag([1.0, 0.6]) - 0.05 * np.array([[0, 1], [1, 0]])) + 0j] * 4)
write_hr(f, Rs, HR)
m = DFTDMFT(f, corr=[0], U=0.0, J=0.0, n_elec=1.4, beta=40, nk=(24, 24, 1))
Nk = m.density_matrix_k(np.zeros((1, m.grid.n_iw), complex), 0.1)
Ee, Ve = np.linalg.eigh(m.Hk)
record("DFT+DMFT", "Wannier density matrix N(k) vs exact Fermi occupations, max error",
       np.abs(Nk - np.einsum("kan,kn,kbn->kab", Ve, fermi(Ee - 0.1, 40), Ve.conj())).max(), 0.0, 1e-5, "exact (U = 0)")
rr = m.solve_csc(hartree_functional(1.5), n_outer=80, mix_rho=0.5, dmft_iter=5, tol_rho=1e-9)
ref, _ = hartree_fock_reference(hk_on_mesh(*read_hr(f), (24, 24, 1)), 1.5, 1.4, 40)
record("DFT+DMFT", "charge self-consistency (U_DMFT = 0): n_orbital 0 vs Hartree mean field", rr["n_orbital"][0], ref[0], 1e-4,
       "independent SCF reference")
record("DFT+DMFT", "charge self-consistency: total electron number", rr["n_orbital"].sum(), 1.4, 1e-4, "electron count")
t = 0.25
write_hr("/tmp/feynlab_one_hr.dat", np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]),
         np.array([[[0j]]] + [[[-t + 0j]]] * 4))
mm = DFTDMFT("/tmp/feynlab_one_hr.dat", U=2.0, J=0.0, n_elec=1.0, beta=10.0, nk=(32, 32, 1), n_iw=256)
r1 = mm.solve("cthyb", n_iter=12, mix=0.6, cthyb_kwargs=dict(n_meas=300000, n_chains=4, n_iw_meas=32))
from feynlab.dmft.bethe import LatticeDMFT  # noqa: E402
k = 2 * np.pi * np.arange(32) / 32
kx, ky = np.meshgrid(k, k, indexing="ij")
L = LatticeDMFT.from_eps_k(10.0, 2.0, -2 * t * (np.cos(kx) + np.cos(ky)), n_iw=256)
r2 = L.solve("cthyb", n_iter=12, mix=0.6, ctint_kwargs=dict(n_meas=600000, n_chains=4, n_iw_meas=32))
record("DFT+DMFT", "one-band lattice: Z (matrix CT-HYB in DFT+DMFT) vs segment CT-HYB", r1["Z"][0], r2["Z"], 0.02,
       "two independent CT-HYB implementations")

# ------------------------------------------------------------------ 4. DCA with cluster CT-INT
from feynlab.dmft.dca import DCA2D  # noqa: E402

d1 = DCA2D(Lc=1, nk=64, U=1.5, beta=10.0, n_iw=256).solve(n_iter=8, mix=0.7, solver="ctint",
                                                         ctint_kwargs=dict(n_meas=400000, n_chains=4, meas_every=40, n_iw_meas=24))
Lb = LatticeDMFT.from_eps_k(10.0, 1.5, -0.5 * (np.cos(np.meshgrid(2 * np.pi * np.arange(64) / 64, 2 * np.pi * np.arange(64) / 64, indexing="ij")[0])
                                               + np.cos(np.meshgrid(2 * np.pi * np.arange(64) / 64, 2 * np.pi * np.arange(64) / 64, indexing="ij")[1])), n_iw=256)
rcs = Lb.solve("ctint", n_iter=8, mix=0.7, ctint_kwargs=dict(n_meas=400000, n_chains=4))
record("DCA", "N_c = 1, cluster CT-INT vs single-site CT-INT DMFT (Z)", d1["Z"][0], rcs["Z"], 0.01, "Hettler et al. (1998)")
d4 = DCA2D(Lc=2, nk=64, U=1.5, beta=10.0, n_iw=256).solve(n_iter=6, mix=0.7, solver="ctint",
                                                         ctint_kwargs=dict(n_meas=300000, n_chains=2, meas_every=60, n_iw_meas=16))
Zk = dict(zip([tuple(np.round(K / np.pi, 1)) for K in d4["K"]], d4["Z"]))
record("DCA", "N_c = 4 (exact cluster solver): Z(pi,0) < Z(0,0)", float(Zk[(1.0, 0.0)] < Zk[(0.0, 0.0)]), 1.0, 0,
       "Maier et al. RMP 77, 1027 (2005)")
record("DCA", "N_c = 4 CT-INT: Z(pi,0)", Zk[(1.0, 0.0)], 0, 0, "informational", kind="info")

# ------------------------------------------------------------------ 5. two-particle vertex, DMFT susceptibility, TPSC
from feynlab.vertex import ImpurityED, vertices, eom_sigma, ed_dmft, lattice_susceptibility, tpsc, tpsc_sigma, tpsc_trace_check  # noqa: E402

ed0 = ImpurityED(-0.2, 0.0, [0.5], [0.6])
gs0 = ed0.generalized_susceptibilities(5.0, n_nu=4, n_om=1)
record("vertex", "U = 0: chi_uu = -beta G G delta (Wick), max error", np.abs(gs0["chi_uu"][1] - np.diag(gs0["chi0"][1])).max(), 0.0, 1e-12,
       "Wick's theorem")
edh = ImpurityED(-1.0, 2.0, [0.0], [0.5])
gsh = edh.generalized_susceptibilities(5.0, n_nu=16, n_om=16)
Seom = eom_sigma(vertices(gsh, 5.0)["F_ud"], gsh, 5.0, 2.0, 2 * edh.occupation(5.0))
record("vertex", "Schwinger-Dyson Sigma(i nu_0) from the ED vertex vs ED Sigma", Seom[16].imag, edh.sigma(1j * gsh["nu"][16:17], 5.0)[0].imag,
       2e-3, "Rohringer et al. RMP 90, 025003 (2018)")
bt, Ub, Db = 8.0, 1.5, 1.0
edm = ed_dmft(bt, Ub, n_bath=2, D=Db)
hh = 0.004
mag = lambda e: e.occupation(bt, 0) - e.occupation(bt, 1)
chi_th = (mag(ed_dmft(bt, Ub, n_bath=2, D=Db, h=hh, mix=0.6)) - mag(ed_dmft(bt, Ub, n_bath=2, D=Db, h=-hh, mix=0.6))) / (2 * hh)
vals = []
for nn in (10, 20):
    gs = edm.generalized_susceptibilities(bt, n_nu=nn, n_om=0)
    v = vertices(gs, bt)
    G = gs["G"]
    c0q = np.array([-bt * G**2 / (1 - Db**2 / 4 * G**2)])
    vals.append(2 * np.real(lattice_susceptibility(v["Gamma_m"], c0q, bt)[0].sum()) / bt**2)
chi_bse = vals[1] + (vals[1] - vals[0])
record("vertex", "DMFT chi_m(q=0) from Bethe-Salpeter (1/N extrapolated) vs dm/dh", chi_bse, chi_th, 0.005,
       "thermodynamic consistency of DMFT (Baym-Kadanoff; Georges et al. RMP 68, 13 (1996))", kind="rel")
from feynlab.selfenergy.lattice import HubbardLattice2D  # noqa: E402
from feynlab.vertex import ladder_dga_sigma  # noqa: E402
Nk, bq, nn = 16, 5.0, 12
kk = 2 * np.pi * np.arange(Nk) / Nk
KX, KY = np.meshgrid(kk, kk, indexing="ij")
epsq = -0.5 * (np.cos(KX) + np.cos(KY))
ratios = []
for Uq in (0.1, 0.05):
    sig0 = (lambda H: (lambda x: np.full(len(x), H + 0j)))(Uq / 2)
    Gloc = (lambda Uq_: (lambda x: np.mean(1 / (1j * x[None, None, :] - epsq[..., None]), axis=(0, 1))))(Uq)
    nus = (2 * np.arange(-nn, nn) + 1) * np.pi / bq
    oms = 2 * np.arange(-nn, nn + 1) * np.pi / bq
    nO, nF = len(oms), len(nus)
    vb = dict(F_d=np.full((nO, nF, nF), Uq + 0j), F_m=np.full((nO, nF, nF), -Uq + 0j),
              Gamma_d=np.full((nO, nF, nF), Uq + 0j), Gamma_m=np.full((nO, nF, nF), -Uq + 0j))
    Sd = ladder_dga_sigma(epsq, Uq / 2, sig0, vb, dict(nu=nus, om=oms, G_of=Gloc), bq, Uq, 1.0)
    dd = Sd[..., nn] - Sd[..., nn].mean()
    s2, _ = HubbardLattice2D(nk=Nk, t=0.25, U=Uq, beta=bq, n_iw=512).sigma("2nd")
    ratios.append(dd[8, 0].imag / (s2[8, 0, 0] - s2[..., 0].mean()).imag)
record("DGA", "ladder DGA / 2nd-order non-local Sigma(pi,0), U=0.05 (weak-coupling limit)", ratios[1], 1.0, 0.1,
       "Katanin, Toschi, Held PRB 80, 075104 (2009): DGA -> 2nd order as U -> 0")
record("DGA", "(ratio - 1) scales linearly in U: [r(0.1)-1]/[r(0.05)-1]", (ratios[0] - 1) / (ratios[1] - 1), 2.0, 0.3,
       "ladder resummation is O(U^3)")
from feynlab.core.lattice import square_lattice_eps  # noqa: E402
eps = square_lattice_eps(32)
rt = tpsc(eps, 5.0, 4.0, 1.0, n_iw=128)
record("TPSC", "U=4t: local spin sum rule", rt["sum_sp"], rt["target_sp"], 1e-8, "Vilk & Tremblay J. Phys. I 7, 1309 (1997)")
record("TPSC", "U=4t: local charge sum rule", rt["sum_ch"], rt["target_ch"], 1e-8, "Vilk & Tremblay (1997)")
r0 = tpsc(eps, 5.0, 1e-6, 1.0, n_iw=256)
record("TPSC", "U -> 0: double occupancy = n^2/4 (Pauli + bosonic tail)", r0["docc"], 0.25, 2e-4, "exact")
diffs = []
for nw in (64, 128):
    rr_ = tpsc(eps, 5.0, 4.0, 1.0, n_iw=nw)
    nus, S = tpsc_sigma(rr_, eps, 5.0, 4.0, 1.0, n_nu=nw)
    a, b = tpsc_trace_check(rr_, nus, S, eps, 5.0, 4.0, 1.0)
    diffs.append(a - b)
record("TPSC", "Tr(Sigma G0) - U<n_up n_dn>, extrapolated (1/N) to infinite cutoff", 2 * diffs[1] - diffs[0], 0.0, 1e-3,
       "exact TPSC identity (Vilk & Tremblay 1997)")

# ------------------------------------------------------------------ 6. transport
from feynlab.transport.models import MassiveDirac, DiracQuadraticMass, QWZ, WannierTB  # noqa: E402
from feynlab.transport.boltzmann import SemiclassicalHall, dirac_finite_range_sigma_xx  # noqa: E402
from feynlab.transport.mechanisms import anomalous_hall_mechanisms  # noqa: E402

d = MassiveDirac(1.0, 0.5)
for mu in (0.7, 1.0, 2.0):
    rs = SemiclassicalHall(d, mu, 1, n_i=1.0, V0=0.1).compute()
    record("Gauss skew", f"Dirac mu={mu}: semiclassical w(4a) skew vs Sinitsyn", rs["skew_gaussian"], d.ref_intrinsic_skew(mu), 1e-6,
           "Sinitsyn et al. PRB 75, 045315 (2007) Eq. (72)")
mq = DiracQuadraticMass(1.0, 0.3, -0.5)
km = anomalous_hall_mechanisms(mq, 0.5, 1, kmax=3.0)
record("Gauss skew", "quadratic-mass Dirac: semiclassical w(4a) vs Kubo-Streda residual", km["skew_gaussian_semiclassical"],
       km["skew_gaussian"], 0.002, "two independent routes", kind="rel")
for rd in (0.5, 2.0):
    rs = SemiclassicalHall(d, 1.0, 1, n_i=1.0, V0=0.1, range_d=rd).compute()
    record("range", f"finite-range impurities d={rd}: sigma_xx vs Bessel formula", rs["sigma_xx"],
           dirac_finite_range_sigma_xx(1, 0.5, 1.0, 1.0, 0.1, rd), 1e-8, "analytic (this work)", kind="rel")
from feynlab.transport.crossed import CrossedDiagrams, LatticeCrossed  # noqa: E402
from feynlab.transport.wannier3d import WannierTB3D, stacked_qwz_hr  # noqa: E402
a_ = CrossedDiagrams(QWZ(-1.7), 0.6, 0.03, N=512, lattice=True).compute()
R_, H_ = stacked_qwz_hr(-1.7, 0.0)
mk = R_[:, 2] == 0
Rn = np.stack([R_[mk, 0] - R_[mk, 1], R_[mk, 1], R_[mk, 2]], 1)
write_hr("/tmp/feynlab_obl_hr.dat", Rn, H_[mk])
b_ = LatticeCrossed(WannierTB("/tmp/feynlab_obl_hr.dat", np.array([[1, 0, 0], [1, 1, 0], [0, 0, 1.0]])), 0.6, 0.03, N=(512, 512)).compute()
record("crossed", "oblique-cell Wannier model vs square-lattice code: X", b_["X"], a_["X"], 1e-5, "same model, different unit cell", kind="rel")
c3 = LatticeCrossed(WannierTB3D(Rs=R_, HR=H_, lat=np.diag([1, 1, 3.0])), 0.6, 0.09, N=(512, 512, 2)).compute()
record("crossed", "3D stacked layers (x c): X vs 2D", 3 * c3["X"], a_["X"], 1e-5, "layer decoupling", kind="rel")

# ------------------------------------------------------------------ 7. continuation, polaron, lattice polarons
from feynlab.continuation.stochastic import stochastic_continuation  # noqa: E402
from feynlab.continuation.nevanlinna import Nevanlinna, hardy_auto  # noqa: E402
from feynlab.core.maxent import maxent  # noqa: E402

w = np.linspace(-4, 4, 801); dw = w[1] - w[0]
A0 = 0.5 * (np.exp(-(w - 1.5) ** 2 / 0.18) + np.exp(-(w + 1.5) ** 2 / 0.18)) / np.sqrt(0.18 * np.pi)
iw = 1j * np.pi * (2 * np.arange(40) + 1) / 50
G = (A0[None] * dw / (iw[:, None] - w[None])).sum(1)
rng = np.random.default_rng(0)
Gn = G + 1e-4 * (rng.standard_normal(40) + 1j * rng.standard_normal(40))
wc = np.linspace(-4, 4, 201); wcd = wc[1] - wc[0]
Aex = np.interp(wc, w, A0)
sac = stochastic_continuation(iw, Gn, 1e-4, wc)
me = maxent(iw, Gn, 1e-4, wc)
Am = me["A"] if isinstance(me, dict) else me
record("continuation", "SAC L1 error (two-Gaussian, sigma=1e-4)", np.sum(np.abs(sac["A"] - Aex)) * wcd, 0.0, 0.15,
       "Shao & Sandvik, Phys. Rep. 1003, 1 (2023)")
record("continuation", "MaxEnt L1 error (same data, for comparison)", np.sum(np.abs(Am - Aex)) * wcd, 0, 0, "informational", kind="info")
wn = w[::4]; dwn = wn[1] - wn[0]; eta = 0.02
Aeta = np.array([np.sum(A0 * dw * (eta / np.pi) / ((x - w) ** 2 + eta**2)) for x in wn])
for nn in (20, 30, 40):
    iwn = 1j * np.pi * (2 * np.arange(nn) + 1) / 50
    Gnn = (A0[None] * dw / (iwn[:, None] - w[None])).sum(1)
    nv = Nevanlinna(iwn, Gnn)
    hr = hardy_auto(nv, wn, eta)
    record("continuation", f"Nevanlinna + Hardy optimisation, {nn} points: L1 error", np.sum(np.abs(hr["A"] - Aeta)) * dwn, 0.0, 0.1,
           "Fei, Yeh, Gull PRL 126, 056402 (2021)")
from feynlab.qmc.polaron import polaron_diagmc, feynman_energy, weak_coupling_mass  # noqa: E402
for al, tol in ((0.2, 0.002), (0.5, 0.005)):
    rp = polaron_diagmc(al, mu=feynman_energy(al) - 0.15, tau_max=40, n_steps=30_000_000, n_chains=4, fit_window=(0.4, 0.97))
    record("polaron", f"Froehlich alpha={al}: m*/m (momentum-fluctuation estimator)", rp["m_star"], weak_coupling_mass(al), tol,
           "Roseler 1968; Smondyrev 1986", kind="abs")
from feynlab.qmc.lattice_polaron import LatticePolaronVED, holstein_weak_coupling, ssh_weak_coupling  # noqa: E402
rv = LatticePolaronVED("holstein", 1.0, 1.0, 1.0, n_gen=16).ground_state()
record("lattice pol.", "1D Holstein t=w=g=1: E0 (VED)", rv["E"], -2.46968472393287, 1e-10,
       "Bonca, Trugman, Batistic PRB 60, 1633 (1999); Ku, Trugman, Bonca PRB 65, 174306 (2002)")
rw = LatticePolaronVED("holstein", 1.0, 1.0, 0.1, n_gen=8).ground_state()
record("lattice pol.", "1D Holstein g=0.1: E0 vs second-order perturbation theory", rw["E"], holstein_weak_coupling(1, 1, 0.1), 5e-6,
       "Rayleigh-Schroedinger")
rs_ = LatticePolaronVED("ssh", 1.0, 1.0, 0.1, n_gen=8).ground_state()
record("lattice pol.", "1D SSH g=0.1: E0 vs second-order perturbation theory", rs_["E"], ssh_weak_coupling(1, 1, 0.1), 1e-4,
       "Marchand et al. PRL 105, 266605 (2010) model")
record("lattice pol.", "1D Holstein t=w=g=1: effective mass m*/m", LatticePolaronVED("holstein", 1.0, 1.0, 1.0, n_gen=14).effective_mass(), 0, 0,
       "informational", kind="info")

# ------------------------------------------------------------------ report
npass = sum(bool(r["passed"]) for r in results)
ninfo = sum(r["passed"] is None for r in results)
rt_ = time.time() - t0
json.dump(dict(results=results, runtime_s=rt_), open(os.path.join(HERE, "results_v11.json"), "w"), indent=1)
lines = ["# feynlab v1.1 benchmark results (closing the v1.0 limitations)", "",
         f"{npass}/{len(results) - ninfo} checks passed, plus {ninfo} informational rows (runtime {rt_:.0f} s)", "",
         "| group | check | feynlab | reference | error | tol | pass | source |", "|---|---|---|---|---|---|---|---|"]
for r in results:
    lines.append(f"| {r['group']} | {r['name']} | {r['value']:.8g} | {r['reference']:.8g} | {r['error']:.2e} "
                 f"| {r['tol']:.1e}{' (rel)' if r['kind']=='rel' else ''} | {'info' if r['passed'] is None else 'yes' if r['passed'] else '**NO**'} | {r['source']} |")
open(os.path.join(HERE, "RESULTS_v11.md"), "w").write("\n".join(lines) + "\n")
print(f"\n{npass}/{len(results) - ninfo} passed (+{ninfo} info) in {rt_:.0f}s")
