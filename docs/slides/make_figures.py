"""Extra figures for the feynlab v1.0 slide deck (run from the repository root)."""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
sys.path.insert(0, os.getcwd())
OUT = os.path.join("docs", "slides", "figures")
os.makedirs(OUT, exist_ok=True)
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


from feynlab.transport.models import MassiveDirac, TiltedDirac, Rashba2DEG
from feynlab.transport.response import orbital_moment, dirac_valley_moment, berry_curvature_dipole, edelstein_boltzmann, rashba_edelstein_reference, kubo_bastin_T
from feynlab.transport.intrinsic import k_grid_polar, intrinsic_hall

fig, ax = plt.subplots(2, 2, figsize=(8.4, 5.6))
k = np.linspace(0, 2, 60)
_, m = orbital_moment(MassiveDirac(1.0, 0.3), k, 0 * k)
ax[0, 0].plot(k, dirac_valley_moment(1, 0.3, k), color=MUTED, lw=3, label="Xiao–Yao–Niu (2007)")
ax[0, 0].plot(k[::3], m[::3, 0], "o", color=C[0], ms=4, label="feynlab, lower band")
ax[0, 0].plot(k[1::3], m[1::3, 1], "s", color=C[1], ms=4, label="feynlab, upper band")
ax[0, 0].set_xlabel("k"); ax[0, 0].set_ylabel("m(k)"); ax[0, 0].set_title("Valley orbital moment (Δ = 0.3)", fontsize=9, loc="left")
ax[0, 0].legend(frameon=False, fontsize=7)
mus = np.linspace(0.12, 0.6, 14)
for i, t in enumerate((0.1, 0.2, 0.3)):
    D = []
    for mu in mus:
        kF = np.sqrt(mu**2 - 0.01)
        g = k_grid_polar(3.0, nr=200, nth=192, dense=[(kF, 0.15, 800)])
        D.append(berry_curvature_dipole(TiltedDirac(1, 0.1, t), mu, g, beta=400)[0])
    ax[0, 1].plot(mus, D, "-o", color=C[i], ms=3, lw=1.6, label=f"tilt t = {t}")
ax[0, 1].axhline(0, color=MUTED, lw=0.8)
ax[0, 1].set_xlabel("μ"); ax[0, 1].set_ylabel("D_x"); ax[0, 1].set_title("Berry-curvature dipole, tilted Dirac (Δ = 0.1)", fontsize=9, loc="left")
ax[0, 1].legend(frameon=False, fontsize=7)
r = Rashba2DEG(1.0, 0.3, 0.0)
SY = np.array([[0, -1j], [1j, 0]]) / 2
mus = np.linspace(-0.03, 1.2, 25)
S = []
for mu in mus:
    kF = np.sqrt(2 * max(mu, 0.0) + 0.09) + 0.3
    g = k_grid_polar(4.0, nr=300, nth=96, dense=[(kF, 0.7, 1500)])
    S.append(abs(edelstein_boltzmann(r, mu, SY, g)))
ax[1, 0].plot(mus, S, "-", color=C[2], lw=2, label="Boltzmann (feynlab)")
ax[1, 0].axhline(rashba_edelstein_reference(1, 0.3), color=MUTED, ls="--", lw=1, label="m α / 4π")
ax[1, 0].axvline(0, color=MUTED, lw=0.6, ls=":")
ax[1, 0].set_xlabel("μ (band crossing at 0)"); ax[1, 0].set_ylabel("|δS_y| / eEτ"); ax[1, 0].set_title("Rashba–Edelstein effect", fontsize=9, loc="left")
ax[1, 0].legend(frameon=False, fontsize=7)
d = MassiveDirac(1.0, 0.5)
gT = k_grid_polar(30.0, nr=2000, nth=12, dense=[(d.kF(1.0), 0.6, 2000)])
Ts = np.linspace(0.01, 0.2, 10)
sea = [intrinsic_hall(d, 1.0, beta=1 / T, grid=gT) + d.lower_band_tail(30.0) for T in Ts]
conv = [kubo_bastin_T(lambda Es: np.array([d.ref_intrinsic(max(e, 0.5)) for e in Es]), 1.0, T) for T in Ts]
ax[1, 1].plot(Ts, conv, color=MUTED, lw=3, label="Bastin convolution")
ax[1, 1].plot(Ts, sea, "o", color=C[3], ms=4, label="finite-T Fermi sea")
ax[1, 1].set_xlabel("k_B T"); ax[1, 1].set_ylabel("σ_int (e²/h)"); ax[1, 1].set_title("Intrinsic AHE at finite T (μ = 1, Δ = 0.5)", fontsize=9, loc="left")
ax[1, 1].legend(frameon=False, fontsize=7)
for a in ax.ravel():
    style(a)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "transport_v1.png"), dpi=170); plt.close(fig)

from feynlab.dmft.bethe import BetheDMFT
fig, ax = plt.subplots(figsize=(4.6, 3.0))
for i, mname in enumerate(("linear", "anderson", "broyden")):
    b = BetheDMFT(100, 2.6)
    rr = b.solve("ipt", mixer=mname, mix=0.5, tol=1e-9, n_iter=500)
    h = np.asarray(b.history, float)
    ax.semilogy(np.arange(1, len(h) + 1), np.abs(h) + 1e-16, "-", color=C[i], lw=1.8, label=f"{mname} ({rr['iterations']} it.)")
ax.set_xlabel("iteration"); ax.set_ylabel("‖Σ_new − Σ_old‖")
ax.set_title("IPT-DMFT, U = 2.6D, β = 100", fontsize=9, loc="left")
ax.legend(frameon=False, fontsize=7); style(ax)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "mixers.png"), dpi=170); plt.close(fig)

from feynlab.qmc.multiorbital import KanamoriAtom
fig, ax = plt.subplots(figsize=(4.6, 3.0))
Js = np.linspace(0, 0.8, 30)
for J in Js:
    e = KanamoriAtom(2, 4.0, J).multiplets(2)
    ax.plot([J] * len(e), e, ".", color=C[0], ms=3)
ax.plot(Js, 4 - 3 * Js, color=C[1], lw=1, label="U − 3J (triplet)")
ax.plot(Js, 4 + Js, color=C[2], lw=1, label="U + J")
ax.plot(Js, 4 - Js, color=C[3], lw=1, label="U − J")
ax.set_xlabel("J"); ax.set_ylabel("E (N = 2)"); ax.set_title("Kanamori two-orbital multiplets, U = 4", fontsize=9, loc="left")
ax.legend(frameon=False, fontsize=7); style(ax)
fig.tight_layout(); fig.savefig(os.path.join(OUT, "kanamori.png"), dpi=170); plt.close(fig)
print("figures written")
