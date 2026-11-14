"""feynlab v1.0 benchmark suite for the ROADMAP modules (A1-D22).
Usage:   python benchmarks/run_v1_benchmarks.py
Writes benchmarks/RESULTS_v1.md, benchmarks/results_v1.json and figures in benchmarks/figures/.
Every row is compared with an exact result, an independent solver, or a published formula.
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
FIG = os.path.join(HERE, "figures")
os.makedirs(FIG, exist_ok=True)
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MUTED = "#6b6a64"
results = []


def record(group, name, value, reference, tol, source, kind="abs"):
    value, reference = float(np.real(value)), float(reference)
    err = abs(value - reference) if kind == "abs" else abs(value - reference) / max(abs(reference), 1e-300)
    ok = bool(err <= tol) if kind != "info" else None
    results.append(dict(group=group, name=name, value=value, reference=reference, error=err, tol=tol, kind=kind,
                        passed=ok, source=source))
    print(f"[{'INFO' if ok is None else 'PASS' if ok else 'FAIL'}] {group:11s} {name:60s} {value: .6g}  ref {reference: .6g}  err {err:.2e}",
          flush=True)


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.grid(True, color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)


import matplotlib  # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

t0 = time.time()
# ===================================================================== D17-D19 infrastructure
from feynlab.utils.errors import ar1_series, gamma_method, jackknife  # noqa: E402
from feynlab.utils.mixing import LinearMixer, AndersonMixer, BroydenMixer, solve_fixed_point  # noqa: E402
from feynlab.dmft.bethe import BetheDMFT  # noqa: E402

x = ar1_series(400000, 0.9, seed=1) if "seed" in ar1_series.__code__.co_varnames else ar1_series(400000, 0.9)
_, err, tau, W = gamma_method(x)
record("errors", "Gamma method tau_int, AR(1) rho=0.9", tau, (1 + 0.9) / (2 * (1 - 0.9)), 0.08,
       "Wolff, CPC 156, 143 (2004); exact (1+rho)/(2(1-rho))", kind="rel")
rng = np.random.default_rng(2)
a_, b_ = 1 + 0.1 * rng.standard_normal(20000), 2 + 0.1 * rng.standard_normal(20000)
jm, je = jackknife(lambda p, q: p / q, a_, b_)[:2]
record("errors", "jackknife <a>/<b> (a~1, b~2)", jm, 0.5, 4 * je + 1e-4, "exact ratio of means")

its, Zs = {}, {}
for m in ("linear", "anderson", "broyden"):
    r = BetheDMFT(100, 2.6).solve("ipt", mixer=m, mix=0.5, tol=1e-9, n_iter=500)
    its[m], Zs[m] = r["iterations"], r["Z"]
record("mixing", "IPT U=2.6 beta=100: Z(Anderson) - Z(linear)", Zs["anderson"] - Zs["linear"], 0.0, 1e-6,
       "same fixed point (Anderson, J. ACM 12, 547 (1965))")
record("mixing", "IPT U=2.6: iterations linear / Anderson", its["linear"] / its["anderson"], 1.0, 0,
       "speed-up factor (informational)", kind="info")

# ===================================================================== A1 CT-HYB
from feynlab.core.matsubara import MatsubaraGrid  # noqa: E402
from feynlab.qmc.cthyb import CTHYBSolver  # noqa: E402
from feynlab.qmc.ed import AndersonED  # noqa: E402

beta = 10.0
g = MatsubaraGrid(beta, 1024)
for U, shift in ((3.0, 0.3), (4.0, 0.0)):
    ed = AndersonED(-U / 2 + shift, U, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
    r = CTHYBSolver(g, U, n_iw_meas=6).solve(ed.hybridization(g.iw), -ed.eps_d, n_meas=1_000_000, n_chains=4)
    record("CT-HYB", f"U={U}: occupation per spin vs ED", np.mean(r["n"]), ed.occupation(beta), 0.01,
           "Werner et al. PRL 97, 076405 (2006) vs exact diagonalization")
    record("CT-HYB", f"U={U}: double occupancy vs ED", r["docc"], ed.double_occupancy(beta), 0.006,
           "Werner et al. PRL 97, 076405 (2006) vs exact diagonalization")
    ex = ed.sigma(r["iw"], beta)
    record("CT-HYB", f"U={U}: Im Sigma(i w_0) vs ED", r["Sigma"][0].imag, ex[0].imag, max(0.03, 4 * abs(r["err_Sigma"][0])),
           "exact diagonalization")

# ===================================================================== A2/B11 multi-orbital
from feynlab.qmc.multiorbital import KanamoriAtom, MultiOrbitalED  # noqa: E402

U, J = 4.0, 0.6
mlt = KanamoriAtom(2, U, J).multiplets(2)
record("Kanamori", "2 orbitals, N=2: lowest multiplet (3-fold) = U - 3J", mlt[0], U - 3 * J, 1e-10,
       "Georges, de' Medici, Mravlje, Annu. Rev. CMP 4, 137 (2013)")
record("Kanamori", "2 orbitals, N=2: highest multiplet = U + J", mlt[-1], U + J, 1e-10, "Georges et al. (2013)")
m3 = KanamoriAtom(3, U, J).multiplets(3)
record("Kanamori", "t2g^3 ground state 4S energy = 3U - 9J", m3[0], 3 * U - 9 * J, 1e-10, "Georges et al. (2013)")
record("Kanamori", "t2g^3 ground-state degeneracy (S=3/2)", np.sum(np.abs(m3 - m3[0]) < 1e-9), 4, 0, "Hund's rule")
z = 1j * np.pi * (2 * np.arange(4) + 1) / 10
a = AndersonED(-1.0, 2.0, [-0.5, 0.7], [0.4, 0.5]).green(z, 10.0)
b = MultiOrbitalED(1, 2.0, 0.0, [-1.0], [[-0.5, 0.7]], [[0.4, 0.5]]).green(z, 10.0)
record("Kanamori", "MultiOrbitalED one-orbital limit, max|G - G_ED|", np.abs(np.asarray(a) - np.asarray(b)).max(), 0, 1e-10,
       "single-orbital AndersonED")

# ===================================================================== A3 DFT+DMFT
from feynlab.dmft.dftdmft import DFTDMFT, write_hr, read_hr  # noqa: E402

t_ = 0.25
Rs = np.array([(0, 0, 0), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)])
HR = np.array([np.diag([0.0, 0.1]) + 0j] + [-t_ * np.diag([1.0, 0.6]) + 0j] * 4)
hr = os.path.join("/tmp", "feynlab_bench_hr.dat")
write_hr(hr, Rs, HR)
R2, H2 = read_hr(hr)
record("DFT+DMFT", "hr.dat write/read round trip, sum|H_R|", np.abs(H2).sum(), np.abs(HR).sum(), 1e-5, "Wannier90 format")
m0 = DFTDMFT(hr, U=0.0, J=0.0, n_elec=2.0, beta=40, nk=(24, 24, 1)).solve("hubbard1", n_iter=30, mix=0.6, verbose=False)
record("DFT+DMFT", "U=0: total occupation (target 2)", m0["n_total"], 2.0, 2e-3, "electron count")
record("DFT+DMFT", "U=0: Z of orbital 1", np.asarray(m0["Z"])[0], 1.0, 1e-3, "non-interacting")
m4 = DFTDMFT(hr, U=4.0, J=0.5, n_elec=2.0, beta=40, nk=(24, 24, 1), dc="occupation").solve("hubbard1", n_iter=80, mix=0.3,
                                                                                         verbose=False)
record("DFT+DMFT", "U=4 J=0.5 half filling, Hubbard-I: max Z (Mott)", np.max(m4["Z"]), 0.0, 0.02,
       "Mott insulator (Lichtenstein & Katsnelson PRB 57, 6884 (1998))")

# ===================================================================== A4 DCA
from feynlab.dmft.dca import DCA2D  # noqa: E402
from feynlab.dmft.bethe import LatticeDMFT  # noqa: E402

U, beta = 1.5, 10.0
r1 = DCA2D(Lc=1, nk=64, U=U, beta=beta, n_iw=256).solve()
k = 2 * np.pi * np.arange(64) / 64
kx, ky = np.meshgrid(k, k, indexing="ij")
refd = LatticeDMFT.from_eps_k(beta, U, -0.5 * (np.cos(kx) + np.cos(ky)), n_iw=256).solve("ipt", tol=1e-8)
record("DCA", "N_c = 1 DCA = single-site DMFT (Z)", r1["Z"][0], refd["Z"], 1e-4, "Hettler et al. PRB 61, 12739 (2000)")
r4 = DCA2D(Lc=2, nk=64, U=U, beta=beta, n_iw=256).solve()
Kl = [tuple(np.round(K / np.pi, 1)) for K in r4["K"]]
Zk = dict(zip(Kl, r4["Z"]))
record("DCA", "N_c = 4: Z(pi,0) < Z(0,0) (antinodal suppression)", float(Zk[(1.0, 0.0)] < Zk[(0.0, 0.0)]), 1.0, 0,
       "Maier et al. RMP 77, 1027 (2005)")

# ===================================================================== B10 two-particle
from feynlab.twoparticle import bubble, rpa  # noqa: E402
from feynlab.core.matsubara import fermi  # noqa: E402
from feynlab.core.lattice import square_lattice_eps  # noqa: E402

g5 = MatsubaraGrid(5.0, 256)
eps = square_lattice_eps(32)
chi = bubble(eps, np.zeros(256), g5, 0.3)
dndmu = np.mean(g5.beta * fermi(eps - 0.3, g5.beta) * (1 - fermi(eps - 0.3, g5.beta)))
record("two-part.", "chi0(q=0, Omega=0) = dn/dmu (compressibility sum rule)", chi[0, 0, 0].real, dndmu, 1e-6, "exact identity")
rp = rpa(chi, 1.0)
record("two-part.", "RPA chi_s(Q) = chi0/(1 - U chi0) at Q=(pi,pi)", rp["spin"][16, 16, 0].real,
       chi[16, 16, 0].real / (1 - chi[16, 16, 0].real), 1e-12, "RPA (Berk-Schrieffer)")

# ===================================================================== B9 NRG
from feynlab.qmc.nrg import NRG, resonant_level_exact, impurity_thermodynamics, kondo_temperature  # noqa: E402

a = NRG(0.0, None, 0.01, 2.5, 26, 400, 2.0, True).run()
b = NRG(0.0, None, 0.01, 2.5, 26, 400, 2.0, False).run()
S = a[:, 1] - b[:, 1]
T = a[:, 0]
S2, T2 = 0.5 * (S[1:] + S[:-1]), np.sqrt(T[1:] * T[:-1])
Se, _ = resonant_level_exact(T2[[12, 16]], 0.01)
record("NRG", "U=0 resonant level: S_imp(T ~ Gamma) / exact", S2[12] / Se[0], 1.0, 0.06, "exact resonant-level entropy")
record("NRG", "U=0 resonant level: S_imp(T << Gamma) / exact", S2[16] / Se[1], 1.0, 0.06, "exact resonant-level entropy")
Gm = 0.01
TKs, rs, curves = [], [6.0, 8.0, 10.0], []
for rr in rs:
    Tt, St, Xt = impurity_thermodynamics(U=rr * Gm, Gamma=Gm, n_iter=60, n_keep=300)
    TKs.append(kondo_temperature(Tt, Xt))
    curves.append((Tt, St, Xt))
Tt, St, Xt = curves[-1]
record("NRG", "Kondo U=10 Gamma: S_imp(T -> 0)", St[-1], 0.0, 0.03, "Wilson, RMP 47, 773 (1975); Krishna-murthy et al. PRB 21, 1003")
im = np.argmin(np.abs(np.log(Tt) - np.log(np.sqrt(TKs[-1] * 10 * Gm))))
record("NRG", "Kondo U=10 Gamma: local-moment plateau S_imp = ln 2", St[im], np.log(2), 0.06, "Krishna-murthy, Wilkins, Wilson (1980)")
xr = np.pi * np.array(rs) / 8 - np.pi / (2 * np.array(rs))
yr = np.log(np.array(TKs) / np.sqrt(np.array(rs) * Gm * Gm / 2))
record("NRG", "Haldane slope d ln T_K / d[pi U/8G - pi G/2U]", np.polyfit(xr, yr, 1)[0], -1.0, 0.05,
       "Haldane, PRL 40, 416 (1978)")
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.3))
for i, (rr, (Tt_, St_, Xt_)) in enumerate(zip(rs, curves)):
    ax[0].semilogx(Tt_ / Gm, St_, color=C[i], lw=2, label=f"U = {rr:.0f} Γ")
    ax[1].semilogx(Tt_ / TKs[i], Xt_, color=C[i], lw=2, label=f"U = {rr:.0f} Γ")
ax[0].axhline(np.log(2), color=MUTED, ls="--", lw=1)
ax[0].axhline(np.log(4), color=MUTED, ls=":", lw=1)
ax[0].set_xlabel("T / Γ"); ax[0].set_ylabel("S_imp"); ax[0].legend(frameon=False, fontsize=8)
ax[1].set_xlabel("T / T_K"); ax[1].set_ylabel("T χ_imp"); ax[1].legend(frameon=False, fontsize=8)
ax[0].set_title("NRG impurity entropy", fontsize=10, loc="left"); ax[1].set_title("Kondo scaling collapse", fontsize=10, loc="left")
for a_x in ax:
    style(a_x)
fig.tight_layout(); fig.savefig(os.path.join(FIG, "nrg_kondo.png"), dpi=160); plt.close(fig)

# ===================================================================== B8 continuation
from feynlab.continuation.nevanlinna import Nevanlinna, nevanlinna_auto  # noqa: E402
from feynlab.continuation.stochastic import stochastic_continuation  # noqa: E402
from feynlab.core.maxent import maxent  # noqa: E402

w = np.linspace(-4, 4, 801)
dw = w[1] - w[0]
A0 = 0.5 * (np.exp(-(w - 1.5) ** 2 / 0.18) + np.exp(-(w + 1.5) ** 2 / 0.18)) / np.sqrt(0.18 * np.pi)
iw = 1j * np.pi * (2 * np.arange(20) + 1) / 50
G = (A0[None, :] * dw / (iw[:, None] - w[None, :])).sum(1)
nv = Nevanlinna(iw, G)
record("continuation", "Nevanlinna interpolant reproduces input, max|dG|", np.abs(nv(iw) - G).max(), 0.0, 1e-10,
       "Fei, Yeh, Gull PRL 126, 056402 (2021)")
record("continuation", "Nevanlinna G(2i) vs exact (analytic continuation off axis)", abs(nv([2j])[0] - (A0 * dw / (2j - w)).sum()), 0.0,
       1e-4, "Fei, Yeh, Gull (2021)")
wc = w[::8]
An = nevanlinna_auto(iw, G, wc, eta=0.05)
An = An[0] if isinstance(An, tuple) else An
An = np.asarray(An["A"] if isinstance(An, dict) else An)
record("continuation", "Nevanlinna spectral normalisation", np.sum(An) * (wc[1] - wc[0]), 1.0, 0.03, "sum rule")
iw40 = 1j * np.pi * (2 * np.arange(40) + 1) / 50
G40 = (A0[None, :] * dw / (iw40[:, None] - w[None, :])).sum(1)
Gn = G40 + 1e-4 * (np.random.default_rng(0).standard_normal(40) + 1j * np.random.default_rng(1).standard_normal(40))
wsc = np.linspace(-4, 4, 201)
sc = stochastic_continuation(iw40, Gn, 1e-4, wsc)
record("continuation", "stochastic continuation: normalisation", np.sum(sc["A"]) * (wsc[1] - wsc[0]), 1.0, 0.02,
       "Mishchenko et al. PRB 62, 6317 (2000); Sandvik PRB 57, 10287 (1998)")
fig, ax = plt.subplots(figsize=(5.4, 3.3))
ax.plot(w, A0, color="#1f1f1e", lw=2.4, label="exact")
ax.plot(wc, An, color=C[0], lw=2, label="Nevanlinna (η = 0.05)")
ax.plot(wsc, sc["A"], color=C[1], lw=1.6, label="stochastic (σ = 1e-4)")
try:
    me = maxent(iw40, Gn, 1e-4, wsc)
    ax.plot(wsc, me["A"] if isinstance(me, dict) else me, color=C[2], lw=1.6, label="MaxEnt")
except Exception:  # noqa: BLE001
    pass
ax.set_xlabel("ω"); ax.set_ylabel("A(ω)"); ax.legend(frameon=False, fontsize=8)
ax.set_title("Analytic continuation of a two-peak spectrum", fontsize=10, loc="left")
style(ax); fig.tight_layout(); fig.savefig(os.path.join(FIG, "continuation.png"), dpi=160); plt.close(fig)

# ===================================================================== C13 CDet / resummation
from feynlab.qmc.cdet import cdet_series, atom_exact_coefficients  # noqa: E402
from feynlab.diagrams.resum import domb_sykes, pade, pade_eval, conformal_map_sum  # noqa: E402

g5 = MatsubaraGrid(5.0, 1024)
ca, ea = cdet_series(g5, 1 / g5.iw, max_order=4, n_samples=200000)
ex = atom_exact_coefficients(5.0, 4)
record("CDet", "Hubbard atom a_1 (double occupancy series) vs exact", ca[1], ex[1], 1e-6 + 4 * ea[1],
       "Rossi PRL 119, 045701 (2017); exact -(1/4) tanh(beta U/4)")
record("CDet", "Hubbard atom a_3 vs exact", ca[3], ex[3], 1e-6 + 4 * ea[3], "exact")
edb = AndersonED(0.0, 0.0, [-1.0, 0.0, 1.0], [0.5, 0.4, 0.5])
c, e = cdet_series(g5, 1 / (g5.iw - edb.hybridization(g5.iw)), max_order=3, n_samples=400000)
Us = np.linspace(-0.3, 0.3, 13)
p = np.polyfit(Us, [AndersonED(-u / 2, u, [-1.0, 0.0, 1.0], [0.5, 0.4, 0.5]).double_occupancy(5.0) - 0.25 for u in Us], 6)[::-1]
record("CDet", "Anderson model a_1 vs ED Taylor fit", c[1], p[1], max(1e-3, 4 * e[1]), "exact diagonalization")
from scipy.special import binom  # noqa: E402
aa = np.array([binom(-0.5, kk) for kk in range(20)])
record("resum", "Domb-Sykes radius of (1+x)^(-1/2)", domb_sykes(aa)[0], 1.0, 1e-6, "Domb & Sykes Proc. R. Soc. A 240, 214 (1957)")
pp, qq = pade(aa, 8, 8)
record("resum", "Pade[8/8] beyond radius, x = 3", pade_eval(pp, qq, 3.0).real, 0.5, 1e-6, "Baker & Graves-Morris (1996)")
record("resum", "conformal-map resummation, x = 3", conformal_map_sum(aa, 1.0, 3.0), 0.5, 1e-6,
       "Rossi, Ohgoe, Van Houcke, Werner PRL 121, 130405 (2018)")

# ===================================================================== C15 polaron Z, m*
from feynlab.qmc.polaron import polaron_diagmc, polaron_effective_mass, feynman_energy, weak_coupling_Z, weak_coupling_mass  # noqa: E402

al = 0.2
kw = dict(mu=feynman_energy(al) - 0.15, tau_max=40, n_steps=30_000_000, n_chains=4, fit_window=(0.4, 0.97))
rz = polaron_diagmc(al, **kw)
record("polaron", "Froehlich alpha=0.2: Z vs 1 - alpha/2", rz["Z"], weak_coupling_Z(al), 0.035,
       "Mishchenko, Prokof'ev, Sakamoto, Svistunov PRB 62, 6317 (2000)")
mm = polaron_effective_mass(al, p=0.5, **kw)
record("polaron", "Froehlich alpha=0.2: m*/m vs 1 + alpha/6 + 0.0236 alpha^2", mm["m_star"], weak_coupling_mass(al), 0.06,
       "Roseler (1968); Smondyrev, Theor. Math. Phys. 68, 653 (1986)")

# ===================================================================== C12/C14 transport extensions
from feynlab.transport.models import MassiveDirac, TiltedDirac, Rashba2DEG, QWZ, DiracQuadraticMass  # noqa: E402
from feynlab.transport.kubo_streda import DisorderedHall  # noqa: E402
from feynlab.transport.intrinsic import k_grid_polar, intrinsic_hall  # noqa: E402
from feynlab.transport.crossed import crossed_hall  # noqa: E402
from feynlab.transport.response import (orbital_moment, dirac_valley_moment, berry_curvature_dipole,  # noqa: E402
                                        edelstein_boltzmann, edelstein_kubo, rashba_edelstein_reference, kubo_bastin_T)
from feynlab.transport.wannier3d import WannierTB3D, stacked_qwz_hr  # noqa: E402

d = MassiveDirac(1.0, 0.5)
mu = 1.0
gp = k_grid_polar(8.0, nr=300, nth=48, dense=[(d.kF(mu), 0.4, 1500), (d.kF(mu), 0.05, 1500)])
V0 = 0.05
ni = 0.005 / V0**2
born = DisorderedHall(d, mu, n_i=ni, V0=V0, grid=gp).sigma_I()
tm = DisorderedHall(d, mu, n_i=ni, V0=V0, grid=gp, t_matrix=True).sigma_I()
skw = DisorderedHall(d, mu, n_i=ni, V0=V0, V1=V0, grid=gp).sigma_I(skew=True)
record("T-matrix", "T-matrix sigma_I = Born + third-order skew (n_i V0^2 fixed)", tm, skw, 0.01,
       "Sinitsyn et al. PRB 75, 045315 (2007) Eq. (72)", kind="rel")
tp = DisorderedHall(d, mu, n_i=0.005 / 0.25, V0=0.5, grid=gp, t_matrix=True).sigma_I()
tn = DisorderedHall(d, mu, n_i=0.005 / 0.25, V0=-0.5, grid=gp, t_matrix=True).sigma_I()
record("T-matrix", "strong scatterers V0=+-0.5: even part ~ Born", 0.5 * (tp + tn), born, 0.03, "Born limit", kind="rel")
record("T-matrix", "identity impurity matrix reproduces scalar disorder",
       DisorderedHall(d, mu, n_i=ni, V0=V0, grid=gp, imp_matrix=np.eye(2)).sigma_I(), born, 1e-12, "identity")

qwz = QWZ(-2 + 0.3)
lat = crossed_hall(qwz, 0.6, lattice=True, runs=((0.03, 0, 1024), (0.015, 0, 2048)))
con = crossed_hall(DiracQuadraticMass(1.0, 0.3, -0.5), 0.6, runs=((0.03, 1300.0, 1024), (0.015, 2600.0, 2048)))
record("crossed", "lattice (BZ-FFT) QWZ X vs continuum k.p limit", lat["X"], con["X"], 0.02,
       "Ado et al. EPL 111, 37004 (2015), k.p limit of QWZ", kind="rel")
record("crossed", "lattice QWZ Psi vs continuum k.p limit", lat["Psi"], con["Psi"], 0.08,
       "Ado et al. (2015); Psi != 0 from the k^2 mass", kind="rel")

kk = np.array([0.0, 0.3, 0.7, 1.5])
_, mk = orbital_moment(MassiveDirac(1.0, 0.3), kk, 0 * kk)
record("OAM", "Dirac valley orbital moment m(k), max deviation", np.abs(mk[:, 0] - dirac_valley_moment(1, 0.3, kk)).max(), 0,
       1e-10, "Xiao, Yao, Niu PRL 99, 236809 (2007)")
td = TiltedDirac(1.0, 0.1, 0.3)
gb = k_grid_polar(3.0, nr=400, nth=256, dense=[(np.sqrt(0.09 - 0.01), 0.2, 1500)])
Ds = berry_curvature_dipole(td, 0.3, gb, beta=300)
Dse = berry_curvature_dipole(td, 0.3, gb, beta=300, form="sea")
record("BCD", "tilted Dirac: D_x Fermi-surface = Fermi-sea", Ds[0], Dse[0], 1e-6, "Sodemann & Fu PRL 115, 216806 (2015)",
       kind="rel")
record("BCD", "untilted Dirac: D_x = 0 (C3/mirror)", berry_curvature_dipole(TiltedDirac(1, 0.1, 0.0), 0.3, gb, beta=300)[0], 0.0,
       1e-10, "Sodemann & Fu (2015)")
r2 = Rashba2DEG(1.0, 0.3, 0.0)
SY = np.array([[0, -1j], [1j, 0]]) / 2
gr = k_grid_polar(4.0, nr=400, nth=128, dense=[(np.sqrt(2.0), 0.5, 2000)])
eb = edelstein_boltzmann(r2, 1.0, SY, gr)
record("Edelstein", "Rashba |delta S_y| per eE tau = m alpha / 4 pi", abs(eb), rashba_edelstein_reference(1.0, 0.3), 0.005,
       "Edelstein, Solid State Commun. 73, 233 (1990)", kind="rel")
ek = edelstein_kubo(r2, 1.0, SY, gr, Gamma=0.005)
record("Edelstein", "Kubo intraband (odd) / tau = Boltzmann", ek["odd"] / ek["tau"], eb, 0.01,
       "Freimuth, Bluegel, Mokrousov PRB 90, 174423 (2014)", kind="rel")
record("Edelstein", "Rashba, no magnetisation: intrinsic (even) S_y = 0", ek["even"], 0.0, 1e-10, "time-reversal symmetry")
gT = k_grid_polar(30.0, nr=3000, nth=16, dense=[(d.kF(mu), 0.5, 3000)])
for Tk in (0.05, 0.1):
    sea = intrinsic_hall(d, mu, beta=1 / Tk, grid=gT) + d.lower_band_tail(30.0)
    conv = kubo_bastin_T(lambda Es: np.array([d.ref_intrinsic(max(e_, 0.5)) for e_ in Es]), mu, Tk)
    record("finite-T", f"Dirac sigma_int(mu=1, T={Tk}): Fermi-sea vs Bastin convolution", sea, conv, 2e-5,
           "Sommerfeld/Bastin; Sinitsyn et al. Eq. (58) at T=0")

Rq, Hq = stacked_qwz_hr(1.0, 0.2)
w3 = WannierTB3D(Rs=Rq, HR=Hq, lat=np.diag([1.0, 1.0, 3.0]))
record("Wannier3D", "stacked Chern layers: sigma_xy = C e^2/(h c)", w3.intrinsic_hall(0.0, nk=(24, 24, 8))["sigma"][2], 1 / 3.0, 1e-4,
       "Halperin, Jpn. J. Appl. Phys. 26, 1913 (1987)")
ref3 = w3.intrinsic_hall(-0.9, nk=(200, 200, 160), adaptive=False, beta=1e3)["sigma"][2]
ad3 = w3.intrinsic_hall(-0.9, nk=(40, 40, 40), beta=1e3)["sigma"][2]
record("Wannier3D", "metal mu=-0.9: adaptive 40^3 vs converged 200x200x160", ad3, ref3, 5e-3,
       "Wang, Yates, Souza, Vanderbilt PRB 74, 195118 (2006)", kind="rel")

# ===================================================================== report
npass = sum(bool(r["passed"]) for r in results)
ninfo = sum(r["passed"] is None for r in results)
rt = time.time() - t0
json.dump(dict(results=results, runtime_s=rt), open(os.path.join(HERE, "results_v1.json"), "w"), indent=1)
lines = ["# feynlab v1.0 benchmark results (ROADMAP modules)", "",
         f"{npass}/{len(results) - ninfo} checks passed, plus {ninfo} informational rows (runtime {rt:.0f} s)", "",
         "| group | check | feynlab | reference | error | tol | pass | source |", "|---|---|---|---|---|---|---|---|"]
for r in results:
    lines.append(f"| {r['group']} | {r['name']} | {r['value']:.6g} | {r['reference']:.6g} | {r['error']:.2e} "
                 f"| {r['tol']:.1e}{' (rel)' if r['kind']=='rel' else ''} | {'info' if r['passed'] is None else 'yes' if r['passed'] else '**NO**'} | {r['source']} |")
lines += ["", "Figures: figures/nrg_kondo.png, figures/continuation.png"]
open(os.path.join(HERE, "RESULTS_v1.md"), "w").write("\n".join(lines) + "\n")
print(f"\n{npass}/{len(results) - ninfo} passed (+{ninfo} info) in {rt:.0f}s")
