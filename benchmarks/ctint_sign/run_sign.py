"""Sign problem of cluster CT-INT on a 2x2 plaquette (DCA Weiss field of the non-interacting lattice, t = 1, U = 6,
mu = U/2 - 1, i.e. hole doped): sign vs beta and cost model, auxiliary-shift optimisation, and a check that the
jackknife error bars of the pooled sign-reweighted estimator match the spread of independent runs.
Usage: python benchmarks/ctint_sign/run_sign.py   (about 15 min on 2 cores); writes results.json."""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from feynlab.core.matsubara import MatsubaraGrid  # noqa: E402
from feynlab.qmc.ctint_cluster import ClusterCTINT, optimize_delta, sign_cost_model  # noqa: E402

T_HOP, U, DMU = 1.0, 6.0, -1.0
SITES = [[0, 0], [1, 0], [0, 1], [1, 1]]
KS = [[0, 0], [np.pi, 0], [0, np.pi], [np.pi, np.pi]]


def weiss(beta):
    g = MatsubaraGrid(beta, 512)
    kk = 2 * np.pi * np.arange(16) / 16 / 2 - np.pi / 2
    G0 = []
    for Kx, Ky in KS:
        kx, ky = np.meshgrid(Kx + kk, Ky + kk, indexing="ij")
        e = -2 * T_HOP * (np.cos(kx) + np.cos(ky))
        G0.append(np.mean(1 / (g.iw[:, None] + U / 2 + DMU - e.ravel()[None, :]), axis=1))
    return g, np.array(G0)


def solver(g, delta):
    return ClusterCTINT(g, U, SITES, [2, 2], KS, delta=delta, n_iw_meas=8)


if __name__ == "__main__":
    out = {}
    betas, signs, rows = [2.0, 4.0, 6.0, 8.0, 10.0], [], []
    for b in betas:
        g, G0 = weiss(b)
        t0 = time.time()
        r = solver(g, 0.51).solve(G0, n_warm=5000, n_meas=200000, meas_every=10, n_chains=2)
        signs.append(r["sign"])
        rows.append(dict(beta=b, sign=r["sign"], sign_err=r["sign_err"], k_mean=r["k_mean"], time=time.time() - t0,
                         Sigma_pi0=[r["Sigma_K"][1, 0].real, r["Sigma_K"][1, 0].imag],
                         err=[r["err_Sigma"][1, 0].real, r["err_Sigma"][1, 0].imag]))
        print(rows[-1], flush=True)
    df, pred = sign_cost_model(betas, signs, 4)
    out["beta_scan"] = rows
    out["cost_model"] = dict(df=df, predicted={b: pred(b) for b in (12.0, 16.0, 20.0)})
    print("df", df, out["cost_model"]["predicted"], flush=True)
    g, G0 = weiss(8.0)
    best, table = optimize_delta(lambda d: solver(g, d), G0, deltas=(0.3, 0.51, 0.75, 1.0),
                                 pilot=dict(n_warm=5000, n_meas=100000, meas_every=10, n_chains=2))
    out["delta_scan"] = dict(best=best, table=table)
    print("delta scan", best, table, flush=True)
    # error-bar check: 8 independent runs vs jackknife
    vals, errs = [], []
    for sd in range(8):
        r = solver(g, 0.51).solve(G0, n_warm=5000, n_meas=100000, meas_every=10, n_chains=2, seed=1000 + 37 * sd)
        vals.append(r["Sigma_K"][1, 0]); errs.append(r["err_Sigma"][1, 0])
    vals = np.array(vals)
    out["error_check"] = dict(spread_im=float(vals.imag.std(ddof=1)), mean_jackknife_im=float(np.mean(np.abs(np.imag(errs)))),
                              spread_re=float(vals.real.std(ddof=1)), mean_jackknife_re=float(np.mean(np.abs(np.real(errs)))))
    print(out["error_check"], flush=True)
    json.dump(out, open(os.path.join(HERE, "results.json"), "w"), indent=1, default=float)
