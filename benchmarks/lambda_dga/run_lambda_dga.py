"""lambda-corrected ladder DGA for the 2D Hubbard model at half filling (t = 0.25, 16x16 k mesh).

DMFT vertex from finite-T ED (2 bath sites, spin-resolved), box of 16 fermionic x 15 bosonic frequencies.
Usage: python benchmarks/lambda_dga/run_lambda_dga.py     (about 35 min on 2 cores; vertices cached in ./cache)
Writes RESULTS.md, results.json and lambda_dga.png in this directory."""
import json
import os
import pickle
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from feynlab.core.lattice import square_lattice_eps  # noqa: E402
from feynlab.vertex import ed_dmft, vertices, ladder_dga_sigma, lambda_dga, tpsc  # noqa: E402

N, T_HOP, NB, NNU, NOM = 16, 0.25, 2, 8, 7
eps = square_lattice_eps(N, T_HOP)
CASES = [(10.0, 1.0), (10.0, 2.0), (13.0, 2.0)]           # (beta, U):  U = 4t and 8t


def vertex(beta, U):
    os.makedirs(os.path.join(HERE, "cache"), exist_ok=True)
    f = os.path.join(HERE, "cache", f"v_b{beta}_U{U}.pkl")
    if os.path.exists(f):
        return pickle.load(open(f, "rb"))
    t0 = time.time()
    ed = ed_dmft(beta, U, n_bath=NB, eps_k=eps, n_iter=100)
    gs = ed.generalized_susceptibilities(beta, n_nu=NNU, n_om=NOM)
    v = vertices(gs, beta)
    allnu = (2 * np.arange(-40, 40) + 1) * np.pi / beta
    d = dict(v=v, nu=gs["nu"], om=gs["om"], G=gs["G"], allnu=allnu, sig=ed.sigma(1j * allnu, beta, 0),
             docc=ed.double_occupancy(beta), dmft_iter=len(ed.history), time=time.time() - t0)
    pickle.dump(d, open(f, "wb"))
    return d


def run(beta, U):
    d = vertex(beta, U)
    allnu, sig = d["allnu"], d["sig"]

    def sd(x):
        x = np.asarray(x)
        out = U / 2 + U**2 / 4 / (1j * x)
        m = np.abs(x) < allnu.max() + 1e-9
        out[m] = sig[np.searchsorted(allnu, x[m] - 1e-9)]
        return out
    Gl = lambda x: np.mean(1 / (1j * x[None, None, :] + U / 2 - eps[..., None] - sd(x)[None, None, :]), axis=(0, 1))
    gs = dict(nu=d["nu"], om=d["om"], G_of=Gl)
    nn, w0 = len(d["nu"]) // 2, len(d["om"]) // 2
    S_lad = ladder_dga_sigma(eps, U / 2, sd, d["v"], gs, beta, U, 1.0)
    r0 = lambda_dga(eps, U / 2, sd, d["v"], gs, beta, U, 1.0, lambda_mode=None, n_tail=256)
    rl = lambda_dga(eps, U / 2, sd, d["v"], gs, beta, U, 1.0, n_tail=256)
    rx = lambda_dga(eps, U / 2, sd, d["v"], gs, beta, U, 1.0, n_tail=256, docc=d["docc"])
    iq = rl["qs"].index((N // 2, N // 2))
    an, nd = (N // 2, 0), (N // 4, N // 4)
    out = dict(beta=beta, U=U, U_over_t=U / T_HOP, T_over_t=1 / (beta * T_HOP),
               ladder_reproduced=float(np.abs(r0["Sigma"] - S_lad).max()),
               lam=rl["lam_m"], lam_min=rl["lam_min"], target=rl["sum_target"], sum_plain=r0["sum_lattice"],
               sum_lambda=rl["sum_lattice"], chiQ_plain=float(r0["chi_m"][iq, w0].real),
               chiQ_lambda=float(rl["chi_m_lambda"][iq, w0].real),
               lam_exact=rx["lam_m"], target_exact=rx["sum_target"], chiQ_exact=float(rx["chi_m_lambda"][iq, w0].real),
               ImS_dmft=float(sd(d["nu"])[nn].imag),
               ImS_plain_AN=float(r0["Sigma"][an + (nn,)].imag), ImS_plain_N=float(r0["Sigma"][nd + (nn,)].imag),
               ImS_lam_AN=float(rl["Sigma"][an + (nn,)].imag), ImS_lam_N=float(rl["Sigma"][nd + (nn,)].imag),
               ImS_exact_AN=float(rx["Sigma"][an + (nn,)].imag), ImS_exact_N=float(rx["Sigma"][nd + (nn,)].imag),
               docc=d["docc"])
    if U / T_HOP <= 4.0:                                  # TPSC is controlled at weak to intermediate coupling
        rt = tpsc(eps, beta, U, 1.0, n_iw=512)
        out["chiQ_tpsc"] = float(rt["chi_sp"][N // 2, N // 2, 0].real) / 2   # chi_sp = 2 chi_m
        out["docc_tpsc"] = float(rt["docc"])
    chi_l = np.real(rl["chi_m_lambda"][:, w0]).reshape(N, N)
    chi_p = np.real(r0["chi_m"][:, w0]).reshape(N, N)
    return out, chi_l, chi_p, rl["Sigma"][..., nn].imag


if __name__ == "__main__":
    res, maps = [], {}
    for beta, U in CASES:
        o, cl, cp, sk = run(beta, U)
        res.append(o)
        maps[(beta, U)] = (cl, cp, sk)
        print(json.dumps(o, indent=1), flush=True)
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    path = [(i, 0) for i in range(N // 2)] + [(N // 2, j) for j in range(N // 2)] + \
        [(N // 2 - i, N // 2 - i) for i in range(N // 2 + 1)]
    fig, ax = plt.subplots(1, 3, figsize=(10.5, 3.0))
    cols = {(10.0, 2.0): "#2a78d6", (13.0, 2.0): "#eb6834", (10.0, 1.0): "#1baf7a"}
    for key, (cl, cp, sk) in maps.items():
        lab = f"U={key[1] / T_HOP:.0f}t, T={1 / (key[0] * T_HOP):.2f}t"
        ax[0].plot([cp[p] for p in path], "--", color=cols[key], lw=1.5)
        ax[0].plot([cl[p] for p in path], "-", color=cols[key], lw=2, label=lab)
    ax[0].set_yscale("symlog", linthresh=1)
    ax[0].set_xticks([0, N // 2, N, len(path) - 1], [r"$\Gamma$", "X", "M", r"$\Gamma$"])
    ax[0].set_title(r"$\chi_m(\mathbf{q},\omega=0)$: $\lambda$ (solid), plain ladder (dashed)", fontsize=8)
    ax[0].legend(fontsize=7, frameon=False)
    for i, key in enumerate([(10.0, 2.0), (13.0, 2.0)]):
        im = ax[i + 1].imshow(np.fft.fftshift(maps[key][2]).T, origin="lower", cmap="Blues_r",
                              extent=[-np.pi, np.pi, -np.pi, np.pi])
        ax[i + 1].set_title(rf"Im$\Sigma_{{\lambda}}(\mathbf{{k}},i\pi T)$, U=8t, T={1 / (key[0] * T_HOP):.2f}t", fontsize=8)
        plt.colorbar(im, ax=ax[i + 1], shrink=0.8)
        ax[i + 1].set_xticks([]); ax[i + 1].set_yticks([])
    for a in ax:
        a.tick_params(labelsize=7)
    fig.tight_layout()
    fig.savefig(os.path.join(HERE, "lambda_dga.png"), dpi=200)
