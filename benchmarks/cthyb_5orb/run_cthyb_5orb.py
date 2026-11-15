"""5-orbital full-Kanamori CT-HYB (SectorCTHYB) against exact Lanczos ED, and speed vs the dense MatrixCTHYB.

Model: 5 orbitals, U = 4, J = 0.6, U' = U - 2J (spin flip + pair hopping), beta = 10, mu = 16.3 (atomic N = 6,
S = 2 Hund's ground state), one bath level per orbital, Delta(iw) = V^2/(iw - e_b), so ED (10 sites, 2^20 states)
is exact for the same impurity model.
Usage: python benchmarks/cthyb_5orb/run_cthyb_5orb.py   (about 1.5 h on 2 cores); writes results.json."""
import json
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", ".."))
from feynlab.core.matsubara import MatsubaraGrid  # noqa: E402
from feynlab.qmc.cthyb_matrix import MatrixCTHYB  # noqa: E402
from feynlab.qmc.cthyb_sector import SectorCTHYB  # noqa: E402
from feynlab.qmc.multiorbital_lanczos import MultiOrbitalLanczos  # noqa: E402

U, J, BETA, MU = 4.0, 0.6, 10.0, 16.3


def ed(eb, V):
    m = MultiOrbitalLanczos(5, U, J, [-MU] * 5, [[eb]] * 5, [[V]] * 5, beta=BETA, n_eig=40, sz_max=10,
                            boltzmann_cut=1e-7)
    wn = np.pi * (2 * np.arange(32) + 1) / BETA
    G = m.green(1j * wn, 0, 0)
    c = m.c
    n = [c[f].T @ c[f] for f in range(10)]
    Sz = 0.5 * sum(n[2 * a] - n[2 * a + 1] for a in range(5))
    out = dict(n=0.0, docc=0.0, Sz2=0.0)
    dim = 2**m.n_so
    for sec, sel, v, E, w in m.gs:
        psi = np.zeros(dim); psi[sel] = v
        out["n"] += w * psi @ (n[0] @ psi); out["docc"] += w * psi @ (n[0] @ (n[1] @ psi))
        x = Sz @ psi; out["Sz2"] += w * x @ x
    out = {k: float(v / m.Zw) for k, v in out.items()}
    out["G"] = G
    return out


def qmc(eb, V, e_cut=None, n_meas=300_000):
    g = MatsubaraGrid(BETA, 512)
    s = SectorCTHYB(g, 5, U, J, mu=MU, e_cut=e_cut)
    t0 = time.time()
    r = s.solve(V**2 / (g.iw - eb), n_warm=20000, n_meas=n_meas, meas_every=20, n_chains=4, seed=1)
    r["time"] = time.time() - t0
    return r


if __name__ == "__main__":
    res = {}
    for eb, V in ((-2.0, 0.5), (-1.0, 1.0)):
        e, q = ed(eb, V), qmc(eb, V)
        res[f"eb={eb},V={V}"] = dict(
            ed=dict(n=e["n"], docc=e["docc"], Sz2=e["Sz2"], G=[[z.real, z.imag] for z in e["G"][:6]]),
            qmc=dict(n=float(np.mean([q[f"n_{f}"] for f in range(10)])),
                     docc=float(np.mean([q[f"docc_{a}"] for a in range(5)])), docc_err=q["docc_0_err"],
                     Sz2=q["Sz^2"], Sz2_err=q["Sz^2_err"], G=[[z.real, z.imag] for z in q["G"].mean(0)[:6]],
                     G_err=[float(abs(x)) for x in q["G_err"][0][:6]], time=q["time"],
                     mean_order=float((np.arange(len(q["order_histogram"])) * q["order_histogram"]).sum())))
        print(json.dumps(res[f"eb={eb},V={V}"], indent=1), flush=True)
    # 3-orbital speed: sector vs dense, same statistics
    g = MatsubaraGrid(BETA, 256)
    mu3 = 2.0 / 2 + 2 * (2 * 2.0 - 5 * 0.4) / 2
    sp = {}
    for name, cls in (("sector", SectorCTHYB), ("dense", MatrixCTHYB)):
        t0 = time.time()
        r = cls(g, 3, 2.0, 0.4, mu=mu3).solve(0.25 / g.iw, n_warm=5000, n_meas=300000, meas_every=20, n_chains=4, seed=5)
        sp[name] = dict(time=time.time() - t0, n0=r["n_0"], n0_err=r["n_0_err"], Sz2=r["Sz^2"], Sz2_err=r["Sz^2_err"])
    res["3orb_speed"] = sp
    print(json.dumps(sp, indent=1))
    json.dump(res, open(os.path.join(HERE, "results.json"), "w"), indent=1)
