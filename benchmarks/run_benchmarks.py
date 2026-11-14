"""feynlab benchmark suite: every number is compared with an exact or published
reference.  Usage:   python benchmarks/run_benchmarks.py [--quick]
Writes benchmarks/RESULTS.md, benchmarks/results.json and figures in benchmarks/figures/.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))

from feynlab.core.matsubara import MatsubaraGrid, g0_tau_level  # noqa: E402
from feynlab.diagrams.enumerate import count_diagrams  # noqa: E402
from feynlab.diagrams.evaluate import sigma_series  # noqa: E402
from feynlab.qmc.ctint import CTINTSolver  # noqa: E402
from feynlab.qmc.ed import AndersonED  # noqa: E402
from feynlab.selfenergy import local  # noqa: E402
from feynlab.selfenergy.lattice import HubbardLattice2D  # noqa: E402
from feynlab.dmft.bethe import BetheDMFT  # noqa: E402
from feynlab.egas import rpa, bcs  # noqa: E402
from feynlab.transport.models import MassiveDirac, Rashba2DEG, QWZ, Haldane  # noqa: E402
from feynlab.transport.intrinsic import intrinsic_hall, chern_numbers, k_grid_polar  # noqa: E402
from feynlab.transport.mechanisms import anomalous_hall_mechanisms  # noqa: E402
from feynlab.transport.kubo_streda import DisorderedHall  # noqa: E402

QUICK = "--quick" in sys.argv
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]  # validated categorical slots 1-4
INK, MUTED = "#1f1f1e", "#6b6a64"

results = []


def record(group, name, value, reference, tol, source, kind="abs"):
    value = float(np.real(value))
    reference = float(reference)
    err = abs(value - reference) if kind == "abs" else abs(value - reference) / max(abs(reference), 1e-300)
    ok = bool(err <= tol) if kind != "info" else None
    results.append(dict(group=group, name=name, value=value, reference=reference, error=err,
                        tol=tol, kind=kind, passed=ok, source=source))
    print(f"[{'INFO' if ok is None else 'PASS' if ok else 'FAIL'}] {group:12s} {name:58s} {value: .6g}  ref {reference: .6g}  err {err:.2e}")


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)


t0 = time.time()
# ============================================================== 1. Matsubara
g = MatsubaraGrid(20.0, 512)
gt = g0_tau_level(0.7, g.tau, g.beta)
record("matsubara", "tau->iw transform, max |error| (free level)", np.abs(g.tau_to_iw(gt) - 1 / (g.iw - 0.7)).max(),
       0.0, 1e-4, "exact 1/(iw - e)")
record("matsubara", "iw->tau transform, max |error| (free level)", np.abs(g.iw_to_tau(1 / (g.iw - 0.7), 1, 0.7, 0.49) - gt).max(),
       0.0, 1e-6, "exact -e^{-e tau}/(1+e^{-beta e})")

# ============================================================== 2. diagram counts
pub_G = {1: 2, 2: 10, 3: 74}
pub_proper = {1: 2, 2: 6, 3: 42}
pub_improper = {1: 0, 2: 4, 3: 32}
for n in (1, 2, 3):
    record("diagrams", f"# connected G diagrams, order {n}", count_diagrams(n, "generic", "G"), pub_G[n], 0,
           "Mohr, physics/0512022; Molinari PRB 71, 113102 (2005)")
    record("diagrams", f"# proper self-energy diagrams, order {n}", count_diagrams(n, "generic", "proper"), pub_proper[n], 0,
           "Mohr, physics/0512022")
    record("diagrams", f"# improper diagrams, order {n}", count_diagrams(n, "generic", "improper"), pub_improper[n], 0,
           "Mohr, physics/0512022")
mol = {2: {0: 2, 1: 1}, 3: {0: 10, 1: 9, 2: 1}}
if not QUICK:
    mol[4] = {0: 74, 1: 91, 2: 23, 3: 1}
for n, ref in mol.items():
    got = count_diagrams(n, "generic", "proper", exclude_tadpoles=True, by_loops=True)
    for l, v in ref.items():
        record("diagrams", f"# Sigma diagrams (Hartree G), order {n}, {l} loops", got.get(l, 0), v, 0,
               "Molinari PRB 71, 113102 (2005), Eq. (11)")
if not QUICK:
    record("diagrams", "# connected G diagrams, order 4", count_diagrams(4, "generic", "G"), 706, 0,
           "Mohr, physics/0512022")
for n, v in {2: 1, 3: 2, 4: 12}.items():
    record("diagrams", f"# Hubbard Sigma diagrams (HF-shifted G0), order {n}",
           count_diagrams(n, "hubbard", "proper", exclude_tadpoles=True), v, 0,
           "standard (e.g. Yamada 1975: 1, 2, 12)")

# ============================================================== 3. diagrammatic MC
beta = 10.0
g = MatsubaraGrid(beta, 512)
ser = sigma_series(g, 1 / g.iw, 1.0, max_order=4, n_iw=1, n_samples=200000 if QUICK else 800000)
record("diagMC", "Hubbard atom: sum of 3rd-order Sigma diagrams (Im, iw0)", ser[3][0][0].imag, 0.0,
       4 * ser[3][1][0], "exact atom: Sigma = U/2 + U^2/(4iw) (orders>2 cancel)")
record("diagMC", "Hubbard atom: sum of 4th-order Sigma diagrams (Im, iw0)", ser[4][0][0].imag, 0.0,
       4 * ser[4][1][0], "exact atom")
record("diagMC", "Hubbard atom: 2nd-order Sigma(iw0) Im", ser[2][0][0].imag, -0.25 / g.iw[0].imag,
       4 * ser[2][1][0], "exact U^2/(4 i w0)")
U = 1.0
ed = AndersonED(-U / 2, U, [-1.0, 0.0, 1.0], [0.5, 0.6, 0.5])
ser = sigma_series(g, 1 / (g.iw - ed.hybridization(g.iw)), U, max_order=4, n_iw=3,
                   n_samples=200000 if QUICK else 600000)
s_tot = U / 2 + sum(v for v, _, _ in ser.values())
err = np.sqrt(sum(e**2 for _, e, _ in ser.values()))
s_ed = ed.sigma(g.iw[:3], beta)
for m in range(3):
    record("diagMC", f"Anderson U=1: Sigma series to 4th order, Im Sigma(iw{m})", s_tot[m].imag, s_ed[m].imag,
           4 * err[m] + 2e-3, "exact diagonalization (this package)")

# ============================================================== 4. CT-INT
g = MatsubaraGrid(10.0, 1024)
U = 3.0
ed = AndersonED(-U / 2 + 0.3, U, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
g0 = 1 / (g.iw - ed.eps_d - ed.hybridization(g.iw))
res = CTINTSolver(g, U, n_iw_meas=10).solve(g0, n_warm=20000, n_meas=300000 if QUICK else 1000000,
                                             meas_every=10, n_chains=4)
ex = ed.sigma(res["iw"], 10.0)
for m in range(4):
    e = abs(res["err_Sigma"][m].imag)
    record("CT-INT", f"Anderson U=3 (no p-h symm.): Im Sigma(iw{m}) vs ED", res["Sigma_up"][m].imag, ex[m].imag,
           3 * e + 3e-3, "exact diagonalization")
record("CT-INT", "average sign (repulsive U, delta=0.51)", res["sign"], 1.0, 1e-9, "Gull et al. RMP 83, 349 (2011)")
fig_ctint = (res, ex)

# ============================================================== 5. weak coupling exact
Gam = 1.0
gY = MatsubaraGrid(400.0, 4096)
g0 = 1 / (gY.iw + 1j * Gam)
for Uy in (0.5, 1.0):
    s2 = local.second_order_sigma(gY, g0, Uy)
    Zinv = 1 / local.z_from_matsubara_extrap(gY.iw, s2, 4)
    ref = (3 - np.pi**2 / 4) * (Uy / (np.pi * Gam)) ** 2
    record("weak", f"Yamada-Yosida: 1/Z - 1 from 2nd-order Sigma, U/Gamma={Uy}", Zinv - 1, ref, 0.03,
           "Yamada, Prog. Theor. Phys. 53, 970 (1975)", kind="rel")
L = HubbardLattice2D(nk=32, U=1.0, beta=5.0, n_iw=256)
L.sigma("2nd")
c = MatsubaraGrid.tail_coefficients(L.grid.iw, L.local_sigma(), 40)
record("weak", "2D Hubbard 2nd-order Sigma tail coefficient U^2 n(1-n)", c[0], 0.25, 1e-3,
       "exact moment sum rule")
Ls = HubbardLattice2D(nk=16, U=0.05, beta=5.0, n_iw=128)
s2, _ = Ls.sigma("2nd")
s3, _ = Ls.sigma("rpa")
record("weak", "spin-fluctuation (RPA/FLEX) -> 2nd order as U->0 (rel. diff, U=0.05)",
       np.abs(s3 - s2).max() / np.abs(s2).max(), 0.0, 0.02, "diagram expansion O(U^3)")

# ============================================================== 6. DMFT
Us = np.arange(1.0, 4.01, 0.2)
Z_met, Z_ins = [], []
for Ub in Us:
    d = BetheDMFT(beta=100, U=Ub, n_iw=2048)
    Z_met.append(d.solve("ipt", n_iter=400, mix=0.6)["Z"])
    d = BetheDMFT(beta=100, U=Ub, n_iw=2048)
    Z_ins.append(d.solve("ipt", n_iter=400, mix=0.6, insulating_start=True)["Z"])
Z_met, Z_ins = np.array(Z_met), np.array(Z_ins)
uc2 = Us[np.argmax(Z_met < 0.01)] if np.any(Z_met < 0.01) else np.nan
uc1 = Us[np.argmax(Z_ins < 0.01)] if np.any(Z_ins < 0.01) else np.nan  # smallest U where insulator survives
info_dmft = dict(Uc1_ipt=float(uc1), Uc2_ipt=float(uc2))
record("DMFT", "IPT Uc2 (metal disappears), beta=100 [IPT approximation; not a pass/fail check]", uc2, 2.94, 0,
       "exact 2.92-3.0D; IPT is known to overestimate U_c2 -> use ED/CT-INT", kind="info")
record("DMFT", "IPT Uc1 (insulator disappears), beta=100 [informational]", uc1, 2.39, 0,
       "exact 2.39(2)D", kind="info")
d = BetheDMFT(beta=10, U=1.0, n_iw=1024)
r_ipt = d.solve("ipt")
d2 = BetheDMFT(beta=10, U=1.0, n_iw=1024)
r_qmc = d2.solve("ctint", n_iter=6 if QUICK else 10, mix=0.7,
                 ctint_kwargs=dict(n_meas=150000 if QUICK else 400000, n_chains=2, n_iw_meas=64))
record("DMFT", "U=D, beta=10: Z(CT-INT DMFT) vs Z(IPT)", r_qmc["Z"], r_ipt["Z"], 0.03,
       "IPT exact to O(U^2); agreement expected at weak coupling", kind="abs")
d3 = BetheDMFT(beta=50, U=6.0, n_iw=1024)
r_h1 = d3.solve("hubbard1", n_iter=100)
gap_proxy = r_h1["Sigma"][0].imag
record("DMFT", "U=6D Hubbard-I: Sigma(iw0) -> atomic U^2/(4 iw0)", gap_proxy,
       (36 / 4) / -d3.grid.iw[0].imag, 0.01, "atomic limit (exact at t=0)", kind="rel")

# ============================================================== 7. electron gas / BCS
rs = 2.0
kF = rpa.kF_of_rs(rs)
q = np.array([0.02, 0.04, 0.06]) * kF
w = rpa.plasmon_dispersion(rs, q)
coef = np.polyfit(q**2, w**2 - rpa.plasma_frequency(rs) ** 2, 1)[0] / kF**2
record("e-gas", "RPA plasmon dispersion coefficient alpha in w^2=wp^2+alpha vF^2 q^2", coef, 0.6, 0.01,
       "Mattuck Ex. 13.6 (3/5)", kind="abs")
e3, e4 = rpa.rpa_correlation_energy(1e-3), rpa.rpa_correlation_energy(1e-4)
record("e-gas", "RPA correlation energy slope d e_c / d ln rs (rs->0)", (e3 - e4) / np.log(10),
       (1 - np.log(2)) / np.pi**2, 5e-4, "Gell-Mann & Brueckner, Phys. Rev. 106, 364 (1957)", kind="rel")
record("e-gas", "RPA correlation energy constant (Ha)", e3 - (1 - np.log(2)) / np.pi**2 * np.log(1e-3), -0.0711, 0.0005,
       "RPA high-density constant -0.0711 Ha", kind="abs")
for lam in (0.2, 0.25):
    record("BCS", f"Delta(0)/Tc, lambda={lam}", bcs.gap_T0(lam) / bcs.Tc(lam), np.pi * np.exp(-np.euler_gamma), 2e-3,
           "BCS weak coupling pi e^-gamma = 1.7639", kind="rel")
record("BCS", "Tc / (w_D e^{-1/lambda}), lambda=0.2", bcs.Tc(0.2) / np.exp(-1 / 0.2), 2 * np.exp(np.euler_gamma) / np.pi,
       1e-3, "BCS 2e^gamma/pi = 1.134", kind="rel")

# ============================================================== 8. intrinsic Hall
for m, Cref in ((-2.5, 0), (-1.0, 1), (1.0, 1), (2.5, 0)):
    sxy = intrinsic_hall(QWZ(m), 0.0, nk=200)
    record("Berry", f"QWZ m={m}: |sigma_xy| (e^2/h) at mu=0", abs(sxy), Cref, 1e-3,
           "Qi, Wu, Zhang PRB 74, 085308 (2006): |C|=1 for 0<|m|<2, 0 otherwise")
    record("Berry", f"QWZ m={m}: Kubo sigma_xy = -C(lower band, FHS lattice method)", sxy,
           -chern_numbers(QWZ(m), 60)[0], 1e-3, "TKNN; Fukui-Hatsugai-Suzuki JPSJ 74, 1674 (2005)")
record("Berry", "Haldane (t2=0.1,phi=pi/2,M=0): |sigma_xy| (e^2/h)", abs(intrinsic_hall(Haldane(), 0.0, nk=240)), 1.0, 1e-3,
       "Haldane PRL 61, 2015 (1988)")
record("Berry", "Haldane trivial (M=1): lowest band Chern number", chern_numbers(Haldane(M=1.0), 90)[0], 0.0, 1e-6,
       "Haldane PRL 61, 2015 (1988)")
r2 = Rashba2DEG(1.0, 0.3)
for mu in (0.5, 1.0):
    kF = np.sqrt(2 * mu)
    gg = k_grid_polar(8.0, nr=4000, nth=8, dense=[(kF, 0.5, 4000)])
    record("Berry", f"Rashba intrinsic SHC (e/2pi) mu={mu}", intrinsic_hall(r2, mu, op=r2.operators["sz/2"], grid=gg),
           0.25, 2e-3, "Sinova et al. PRL 92, 126603 (2004): e/8pi")

# ============================================================== 9. AHE mechanisms
dirac = MassiveDirac(1.0, 0.5)
mus = [0.6, 0.7, 0.85, 1.0, 1.3, 1.6, 2.0, 2.5] if not QUICK else [0.7, 1.0, 2.0]
mech = []
for mu in mus:
    kmax = 8.0
    gg = k_grid_polar(kmax, nr=6000, nth=8, dense=[(dirac.kF(mu), 0.05, 4000)])
    intr = intrinsic_hall(dirac, mu, grid=gg) + dirac.lower_band_tail(kmax)
    r = anomalous_hall_mechanisms(dirac, mu, band=1, V1_over_V0=0.3, nV2_list=(0.01, 0.005),
                                  intrinsic_value=intr)
    mech.append((mu, r))
    V0 = np.sqrt(0.005)
    if mu in (0.7, 1.0, 2.0):
        src = "Sinitsyn et al. PRB 75, 045315 (2007)"
        record("AHE", f"Dirac mu={mu}: intrinsic (e^2/h)", r["intrinsic"], dirac.ref_intrinsic(mu), 2e-3, src + " Eq. 58")
        record("AHE", f"Dirac mu={mu}: side jump", r["side_jump"], dirac.ref_side_jump(mu), 2e-3, src + " Eq. 65")
        record("AHE", f"Dirac mu={mu}: anomalous distribution", r["anomalous_distribution"],
               dirac.ref_anomalous_distribution(mu), 2e-3, src + " Eq. 67")
        record("AHE", f"Dirac mu={mu}: Gaussian (intrinsic) skew", r["skew_gaussian"], dirac.ref_intrinsic_skew(mu), 3e-3,
               src + " Eq. 72, 1st term")
        record("AHE", f"Dirac mu={mu}: conventional skew (V1/V0=0.3)", r["skew_V1"],
               dirac.ref_skew_V1(mu, 1.0, V0, 0.3 * V0), 2e-3, src + " Eq. 72, 2nd term")
        record("AHE", f"Dirac mu={mu}: total, Gaussian disorder (NCA)", r["total_gaussian_nca"],
               dirac.ref_total_gaussian(mu), 3e-3, src + " Eq. B8")
        D = 0.5
        k2 = mu**2 - D**2
        inv_ttr = (k2 + 4 * D**2) / (4 * (k2 + D**2)) * 0.005 * mu
        record("AHE", f"Dirac mu={mu}: sigma_xx Boltzmann vs Drude (e^2/h)", r["sigma_xx"],
               np.pi * mu / (2 * np.pi) * (k2 / mu**2) / inv_ttr, 1e-3, src + " Eq. 54 (tau_tr)", kind="rel")
for mu in (0.5, 1.0):
    kF = np.sqrt(2 * mu)
    gg = k_grid_polar(6.0, nr=300, nth=48, dense=[(kF, 0.6, 1500), (kF, 0.1, 1500)])
    dh = DisorderedHall(r2, mu, n_i=1.0, V0=np.sqrt(0.002), grid=gg)
    res = dh.decompose(op=r2.operators["sz/2"])
    record("SHE", f"Rashba SHC with ladder vertex, total (e/2pi), mu={mu}", res["total"], 0.0, 0.015,
           "Inoue, Bauer, Molenkamp PRB 70, 041303(R) (2004): exact cancellation")


# ============================================================== 10. v0.2: crossed diagrams (Ado et al.)
from feynlab.transport.crossed import crossed_hall  # noqa: E402
eps_, m_ = 1.0, 0.5
cr = crossed_hall(dirac, eps_)
srcA = "Ado, Dmitriev, Ostrovsky, Titov, EPL 111, 37004 (2015)"
record("AHE-X", "Dirac eps=1: X diagram (e^2/h)", cr["X"], 4 * eps_ * m_ * (eps_**2 - m_**2) / (eps_**2 + 3 * m_**2) ** 2,
       0.015, srcA + " Eq. 17", kind="rel")
record("AHE-X", "Dirac eps=1: Psi diagram (e^2/h)", cr["Psi"], 0.0, 0.01, srcA + " Eq. 17 (Psi = 0)")
nca = [r for m0, r in mech if abs(m0 - 1.0) < 1e-12][0]["total_gaussian_nca"]
record("AHE-X", "Dirac eps=1: NCA total vs Ado Eq. 3", nca, -4 * eps_ * m_ * (eps_**2 + m_**2) / (eps_**2 + 3 * m_**2) ** 2,
       3e-3, srcA + " Eq. 3")
record("AHE-X", "Dirac eps=1: full sigma_xy (NCA + X + Psi)", nca + cr["X"] + cr["Psi"],
       -8 * eps_ * m_**3 / (eps_**2 + 3 * m_**2) ** 2, 0.01, srcA + " Eq. 4")

# ============================================================== 11. v0.2: polaron DiagMC
from feynlab.qmc.polaron import polaron_diagmc, feynman_energy, weak_coupling_energy  # noqa: E402
srcH = "Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I"
for a, ref in ((1, -1.0130308), (3, -3.1333335), (5, -5.4401445)):
    record("polaron", f"Feynman variational E0, alpha={a}", feynman_energy(a), ref, 2e-6,
           "Feynman, Phys. Rev. 97, 660 (1955); values from " + srcH)
pol_ref = {0.5: (weak_coupling_energy(0.5), 1e-4, "exact weak-coupling series (Smondyrev 1986)"),
           1.0: (-1.01662, 0.00047, srcH), 2.0: (-2.06957, 0.00084, srcH), 3.0: (-3.16829, 0.00136, srcH)}
for a, (ref, eref, src) in pol_ref.items():
    r = polaron_diagmc(a, mu=feynman_energy(a) - 0.15, tau_max=40, n_steps=40_000_000 if QUICK else 150_000_000,
                       n_chains=8, fit_window=(0.4, 0.97), seed=7)
    record("polaron", f"DiagMC E0 (bare series, all orders), alpha={a}  [+-{r['E_err']:.4f}]", r["E"], ref,
           3 * np.hypot(r["E_err"], eref), src)

# ============================================================== 12. v0.2: bold-line / skeleton series (Kozik et al.)
from feynlab.diagrams.bold import skeleton_series_at  # noqa: E402
from feynlab.selfenergy.local import second_order_sigma as _s2  # noqa: E402
gk = MatsubaraGrid(1.0, 512)
srcK = "Kozik, Ferrero, Georges, PRL 114, 156402 (2015)"
for bU, expect_phys in ((0.5, True), (6.0, False)):
    Uk = bU
    Gex = 0.5 * (1 / (gk.iw + Uk / 2) + 1 / (gk.iw - Uk / 2))
    ser = skeleton_series_at(gk, Gex, Uk, 4, n_iw=1, n_samples=200000 if QUICK else 400000)
    tot = _s2(gk, Gex, Uk)[0].imag + ser[3][0][0].imag + ser[4][0][0].imag
    ratio = tot / (Uk**2 / (4 * gk.iw[0])).imag
    if expect_phys:
        record("bold", f"atom betaU={bU}: skeleton series(G_exact)/Sigma_exact", ratio, 1.0, 0.005,
               srcK + ": physical branch at weak coupling")
    else:
        record("bold", f"atom betaU={bU}: skeleton series(G_exact) misses Sigma_exact (1-ratio > 0.5)",
               float(1 - ratio > 0.5), 1.0, 0, srcK + ": convergence to unphysical branch at strong coupling")

# ============================================================== 13. v0.2: ED-DMFT Mott transition (T = 0)
from feynlab.dmft.bethe import LatticeDMFT  # noqa: E402
from feynlab.qmc.ed_sector import SectorED  # noqa: E402
a_ = AndersonED(-1.2, 3.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
b_ = SectorED(-1.2, 3.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
z_ = 1j * np.pi * (2 * np.arange(5) + 1) / 10
record("ED", "sector-resolved ED vs dense ED, max|dG|", np.abs(a_.green(z_, 10.) - b_.green(z_, 10.)).max(), 0.0, 1e-10, "exact")


def _ed_branch(Ulist, insulating, nb):
    x = S = Uo = None
    out = []
    for U in Ulist:
        d = BetheDMFT(beta=200, U=U, n_iw=1024)
        d._ed_sigma_fn = None
        if x is not None:
            d._bath_x = x
        r = d.solve("ed", n_iter=150, mix=0.5, tol=1e-6, ctint_kwargs=dict(n_bath=nb),
                    insulating_start=(insulating and S is None),
                    sigma_init=None if S is None else S - Uo / 2 + U / 2)
        x, S, Uo = d._bath_x, r["Sigma"], U
        out.append((U, float(r["Z"]), r["iterations"]))
    return out


sys.path.insert(0, os.path.join(HERE, "ed_dmft"))
from analyze import uc1 as _uc1, uc2_linear as _uc2lin  # noqa: E402
srcU = ("consensus U_c2 = 2.92-3.0D, U_c1 = 2.39(2)D: Karski, Raas, Uhrig PRB 72, 113110 (2005) and refs; "
        "Bulla PRL 83, 136 (1999)")
met = _ed_branch(list(np.round(np.arange(2.0, 2.91, 0.1), 3)), False, 7)         # odd bath -> metal
uc2_ed = float(_uc2lin(met)[0])
record("DMFT-ED", "Lanczos ED-DMFT (7 bath = odd, T=0) U_c2, linear Z->0 (window 0.04<=Z<=0.15)", uc2_ed,
       2.96, 0.04, srcU)
ins = _ed_branch(list(np.round(np.arange(3.0, 2.24, -0.05), 3)), True, 6)        # even bath -> insulator
uc1_ed = float(_uc1(ins))
record("DMFT-ED", "Lanczos ED-DMFT (6 bath = even, T=0) U_c1 (insulator lost)", uc1_ed, 2.39, 0.05, srcU)
info_dmft.update(Uc2_ed=uc2_ed, Uc1_ed=uc1_ed, ed_metal=met, ed_insulator=ins)
z2 = [z for u, z, _ in met if abs(u - 2.0) < 1e-9][0]
from feynlab.qmc.ed_sector import LanczosED as _LZ  # noqa: E402
_args = (-1.5, 3.0, [-1.0, -0.2, 0.2, 1.0], [0.5, 0.3, 0.3, 0.5])
_zz = 1j * np.pi * (2 * np.arange(6) + 1) / 200
record("DMFT-ED", "T=0 Lanczos ED vs full sector ED (4 bath), max|dG|",
       np.abs(_LZ(*_args).green(_zz) - SectorED(*_args).green(_zz, 1e5)).max(), 0.0, 1e-8, "exact")

a1 = BetheDMFT(50, 2.5).solve("ipt")
b1 = LatticeDMFT.semicircle(50, 2.5).solve("ipt")
record("DMFT", "LatticeDMFT(semicircle DOS) Z == BetheDMFT Z (U=2.5)", b1["Z"], a1["Z"], 1e-4, "identity")

# ============================================================== 14. v0.2: MaxEnt
from feynlab.core.maxent import maxent  # noqa: E402
from feynlab.core.lattice import bethe_dos  # noqa: E402
wme = np.linspace(-4, 4, 401)
iwm = 1j * np.pi * (2 * np.arange(200) + 1) / 50
A0 = 0.5 * np.exp(-(wme - 1.5) ** 2 / 0.18) / np.sqrt(0.18 * np.pi) + 0.5 * np.exp(-(wme + 1.5) ** 2 / 0.18) / np.sqrt(0.18 * np.pi)
Gme = (A0[None, :] * (wme[1] - wme[0]) / (iwm[:, None] - wme[None, :])).sum(1)
rng_ = np.random.default_rng(0)
rme = maxent(iwm, Gme + 1e-4 * (rng_.standard_normal(200) + 1j * rng_.standard_normal(200)), 1e-4, wme)
record("MaxEnt", "two-peak spectrum: chi^2 per point", rme["chi2_per_point"], 1.0, 0.2, "Jarrell & Gubernatis Phys. Rep. 269, 133 (1996)")
record("MaxEnt", "two-peak spectrum: L1 error of A(w)", np.sum(np.abs(rme["A"] - A0)) * (wme[1] - wme[0]), 0.0, 0.3,
       "synthetic test (noise 1e-4)")

# ============================================================== figures
import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

plt.rcParams.update({"font.size": 10, "axes.edgecolor": MUTED, "text.color": INK, "axes.labelcolor": INK})

# (a) AHE mechanisms vs mu, lines = analytic, markers = feynlab
fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
mm = np.linspace(0.52, 2.6, 200)
curves = [("intrinsic", dirac.ref_intrinsic, "intrinsic"),
          ("side jump + adist", lambda m: 2 * dirac.ref_side_jump(m), None),
          ("Gaussian skew", dirac.ref_intrinsic_skew, "skew_gaussian"),
          ("total (V1=0)", dirac.ref_total_gaussian, "total_gaussian_nca")]
for i, (lab, f, key) in enumerate(curves):
    ax.plot(mm, f(mm), color=C[i], lw=2, label=lab + " (analytic)")
    xs = [m for m, _ in mech]
    if key is None:
        ys = [r["side_jump"] + r["anomalous_distribution"] for _, r in mech]
    else:
        ys = [r[key] for _, r in mech]
    ax.plot(xs, ys, "o", ms=6, color=C[i], mec="white", mew=1.5)
    ax.annotate(lab, (mm[-1], f(mm[-1])), xytext=(4, 0), textcoords="offset points", color=INK,
                fontsize=8, va="center")
ax.set_xlabel("Fermi energy  mu / v (Delta = 0.5)")
ax.set_ylabel("sigma_xy  (e^2/h)")
ax.set_title("Massive Dirac AHE: feynlab (dots) vs Sinitsyn et al. 2007 (lines)", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8, loc="lower left")
ax.set_xlim(0.5, 3.1)
style(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "ahe_dirac_mechanisms.png"))

# (b) IPT Z(U) hysteresis
fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
ax.plot(Us, Z_met, "-o", color=C[0], lw=2, ms=5, mec="white", label="from metallic start")
ax.plot(Us, Z_ins, "-s", color=C[1], lw=2, ms=5, mec="white", label="from insulating start")
ax.axvline(2.94, color=MUTED, lw=1, ls="--")
ax.annotate("NRG U_c2 = 2.94D", (2.94, 0.6), xytext=(4, 0), textcoords="offset points", fontsize=8, color=MUTED)
ax.set_xlabel("U / D")
ax.set_ylabel("quasiparticle weight Z")
ax.set_title("Bethe-lattice DMFT (IPT), T = 0.01D: Mott hysteresis", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8)
style(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "ipt_mott_hysteresis.png"))

# (c) CT-INT vs ED
res, ex = fig_ctint
fig, ax = plt.subplots(figsize=(6.4, 4.0), dpi=150)
wn = res["iw"].imag
ax.plot(wn, ex.imag, color=C[0], lw=2, label="exact diagonalization")
ax.errorbar(wn, res["Sigma_up"].imag, yerr=3 * np.abs(res["err_Sigma"].imag), fmt="o", color=C[1], ms=6,
            mec="white", label="CT-INT (3 sigma)")
ax.set_xlabel("Matsubara frequency  w_n")
ax.set_ylabel("Im Sigma(i w_n)")
ax.set_title("CT-INT vs ED, Anderson impurity U=3, beta=10", fontsize=10, loc="left")
ax.legend(frameon=False, fontsize=8)
style(ax)
fig.tight_layout()
fig.savefig(os.path.join(FIG, "ctint_vs_ed.png"))

# ============================================================== report
npass = sum(bool(r["passed"]) for r in results)
ninfo = sum(r["passed"] is None for r in results)
with open(os.path.join(HERE, "results.json"), "w") as f:
    json.dump(dict(results=results, dmft=info_dmft, runtime_s=time.time() - t0), f, indent=1)
lines = ["# feynlab benchmark results", "",
         f"{npass}/{len(results) - ninfo} checks passed, plus {ninfo} informational rows  (runtime {time.time()-t0:.0f} s, {'quick' if QUICK else 'full'} mode)", "",
         f"IPT coexistence (beta=100): U_c1 ~ {info_dmft['Uc1_ipt']:.1f}D, U_c2 ~ {info_dmft['Uc2_ipt']:.1f}D; "
         f"Lanczos ED-DMFT (6/7 bath sites, T=0): U_c1 ~ {info_dmft.get('Uc1_ed', float('nan')):.2f}D, U_c2 ~ {info_dmft.get('Uc2_ed', float('nan')):.2f}D "
         "(consensus 2.39D / 2.92-3.0D).", "",
         "| group | check | feynlab | reference | error | tol | pass | source |",
         "|---|---|---|---|---|---|---|---|"]
for r in results:
    lines.append(f"| {r['group']} | {r['name']} | {r['value']:.6g} | {r['reference']:.6g} | {r['error']:.2e} "
                 f"| {r['tol']:.1e}{' (rel)' if r['kind']=='rel' else ''} | {'info' if r['passed'] is None else 'yes' if r['passed'] else '**NO**'} | {r['source']} |")
lines += ["", "Figures: figures/ahe_dirac_mechanisms.png, figures/ipt_mott_hysteresis.png, figures/ctint_vs_ed.png"]
open(os.path.join(HERE, "RESULTS.md"), "w").write("\n".join(lines))
print(f"\n{npass}/{len(results) - ninfo} passed (+{ninfo} info) in {time.time()-t0:.0f}s")
