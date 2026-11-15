"""Figures for the v1.1 slides (run from the repository root)."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd())
OUT = os.path.join("docs", "slides", "figures")
C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100"]
MUTED = "#6b6a64"


def style(ax):
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    for s in ("left", "bottom"):
        ax.spines[s].set_color(MUTED)
    ax.tick_params(colors=MUTED, labelsize=8)
    ax.grid(True, color="#e6e5e0", lw=0.6)
    ax.set_axisbelow(True)


# ---------- continuation comparison
from feynlab.continuation.stochastic import stochastic_continuation
from feynlab.continuation.nevanlinna import Nevanlinna, hardy_auto
from feynlab.core.maxent import maxent
w = np.linspace(-4, 4, 801); dw = w[1] - w[0]
A0 = 0.5 * (np.exp(-(w - 1.5) ** 2 / 0.18) + np.exp(-(w + 1.5) ** 2 / 0.18)) / np.sqrt(0.18 * np.pi)
iw = 1j * np.pi * (2 * np.arange(40) + 1) / 50
G = (A0[None] * dw / (iw[:, None] - w[None])).sum(1)
rng = np.random.default_rng(0)
Gn = G + 1e-4 * (rng.standard_normal(40) + 1j * rng.standard_normal(40))
wc = np.linspace(-4, 4, 201)
sac = stochastic_continuation(iw, Gn, 1e-4, wc)
me = maxent(iw, Gn, 1e-4, wc); Am = me["A"] if isinstance(me, dict) else me
iw30 = iw[:30]; nv = Nevanlinna(iw30, G[:30]); wn = w[::4]
raw = nv.spectral(wn, 0.02); hd = hardy_auto(nv, wn, 0.02)
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.2))
ax[0].plot(w, A0, color="#1f1f1e", lw=2.6, label="exact")
ax[0].plot(wc, Am, color=C[0], lw=1.6, label="MaxEnt")
ax[0].plot(wc, sac["A"], color=C[1], lw=1.6, label="SAC (v1.1)")
ax[0].set_title("Noisy data (σ = 1e-4)", fontsize=9, loc="left"); ax[0].legend(frameon=False, fontsize=7)
ax[1].plot(w, A0, color="#1f1f1e", lw=2.6, label="exact")
ax[1].plot(wn, raw, color=C[3], lw=1.2, label="Nevanlinna, θ = 0")
ax[1].plot(wn, hd["A"], color=C[2], lw=1.8, label="Nevanlinna + Hardy")
ax[1].set_ylim(-0.05, 1.2)
ax[1].set_title("Noise-free data, 30 Matsubara points", fontsize=9, loc="left"); ax[1].legend(frameon=False, fontsize=7)
for a in ax:
    a.set_xlabel("ω"); a.set_ylabel("A(ω)"); style(a)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "continuation_v11.png"), dpi=170); plt.close(fig)

# ---------- finite-range transport + DGA weak-coupling ratio
from feynlab.transport.models import MassiveDirac
from feynlab.transport.boltzmann import SemiclassicalHall
d = MassiveDirac(1.0, 0.5)
ds = np.linspace(0, 2.5, 11)
rows = [SemiclassicalHall(d, 1.0, 1, n_i=1.0, V0=0.1, n_theta=96, range_d=x).compute() for x in ds]
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.2))
ax[0].plot(ds, [r["side_jump"] for r in rows], "-o", color=C[0], ms=3, label="side jump")
ax[0].plot(ds, [r["anomalous_distribution"] for r in rows], "--", color=C[1], label="anomalous distribution")
ax[0].plot(ds, [r["skew_gaussian"] for r in rows], "-s", color=C[2], ms=3, label="Gaussian skew (w⁽⁴ᵃ⁾)")
ax[0].set_xlabel("impurity range d (k_F d ≈ 0.87 d)"); ax[0].set_ylabel("σ_xy (e²/h)")
ax[0].set_title("Massive Dirac, μ = 1: smooth disorder", fontsize=9, loc="left"); ax[0].legend(frameon=False, fontsize=7)
from feynlab.vertex import ladder_dga_sigma
from feynlab.selfenergy.lattice import HubbardLattice2D
beta, N, nn = 5.0, 16, 12
k = 2 * np.pi * np.arange(N) / N; KX, KY = np.meshgrid(k, k, indexing="ij"); eps = -0.5 * (np.cos(KX) + np.cos(KY))
Us = [0.05, 0.1, 0.2, 0.3, 0.4]; rat = []
for U in Us:
    nus = (2 * np.arange(-nn, nn) + 1) * np.pi / beta; oms = 2 * np.arange(-nn, nn + 1) * np.pi / beta
    nO, nF = len(oms), len(nus)
    vb = dict(F_d=np.full((nO, nF, nF), U + 0j), F_m=np.full((nO, nF, nF), -U + 0j), Gamma_d=np.full((nO, nF, nF), U + 0j), Gamma_m=np.full((nO, nF, nF), -U + 0j))
    S = ladder_dga_sigma(eps, U / 2, lambda x: np.full(len(x), U / 2 + 0j), vb,
                         dict(nu=nus, om=oms, G_of=lambda x: np.mean(1 / (1j * x[None, None, :] - eps[..., None]), axis=(0, 1))), beta, U, 1.0)
    dd = S[..., nn] - S[..., nn].mean()
    s2, _ = HubbardLattice2D(nk=N, t=0.25, U=U, beta=beta, n_iw=512).sigma("2nd")
    rat.append(dd[8, 0].imag / (s2[8, 0, 0] - s2[..., 0].mean()).imag)
ax[1].plot(Us, rat, "-o", color=C[0], ms=4)
ax[1].axhline(1, color=MUTED, ls="--", lw=1)
ax[1].set_xlabel("U (D = 1)"); ax[1].set_ylabel("Σ_DΓA / Σ⁽²⁾  (non-local, k = (π,0))")
ax[1].set_title("Ladder DΓA → 2nd order as U → 0", fontsize=9, loc="left")
for a in ax:
    style(a)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "transport_dga_v11.png"), dpi=170); plt.close(fig)

# ---------- Holstein convergence and dispersion
from feynlab.qmc.lattice_polaron import LatticePolaronVED
gens = list(range(3, 15))
E = [LatticePolaronVED("holstein", 1.0, 1.0, 1.0, n_gen=n).ground_state()["E"] for n in gens]
fig, ax = plt.subplots(1, 2, figsize=(8.6, 3.0))
ax[0].semilogy(gens, np.abs(np.array(E) + 2.46968472393287071561) + 1e-16, "-o", color=C[0], ms=4)
ax[0].set_xlabel("generations N_gen"); ax[0].set_ylabel("|E₀ − E_exact|")
ax[0].set_title("Holstein t = ω = g = 1: exponential convergence", fontsize=9, loc="left")
v = LatticePolaronVED("holstein", 1.0, 1.0, 1.0, n_gen=11)
ks = np.linspace(0, np.pi, 17)
Ek = [v.ground_state(kk)["E"] for kk in ks]
ax[1].plot(ks, Ek, "-o", color=C[1], ms=3, label="polaron band E(k)")
ax[1].plot(ks, -2 * np.cos(ks), "--", color=MUTED, label="bare −2t cos k")
ax[1].set_xlabel("k"); ax[1].set_ylabel("E"); ax[1].legend(frameon=False, fontsize=7)
ax[1].set_title("Polaron band flattening", fontsize=9, loc="left")
for a in ax:
    style(a)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "holstein_v11.png"), dpi=170); plt.close(fig)
print("v1.1 figures written")
