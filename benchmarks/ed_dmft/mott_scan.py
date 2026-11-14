"""T = 0 ED-DMFT (Lanczos) Mott transition on the Bethe lattice vs bath size.
Usage: python benchmarks/ed_dmft/mott_scan.py N_BATH  ->  benchmarks/ed_dmft/scan_nbN.json"""
import json, os, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from feynlab.dmft.bethe import BetheDMFT

nb = int(sys.argv[1])
which = sys.argv[2] if len(sys.argv) > 2 else "both"   # 'metal', 'insulator' or 'both' 
def branch(Ulist, insulating):
    x = S = Uo = None; out = []
    for U in Ulist:
        if insulating and out and out[-1][1] > 0.005:
            break  # insulator lost: rest of the branch is the metal
        if not insulating and out and out[-1][1] < 0.002:
            break  # metal lost
        t = time.time()
        d = BetheDMFT(beta=200, U=U, n_iw=1024); d._ed_sigma_fn = None
        if x is not None: d._bath_x = x
        r = d.solve("ed", n_iter=150, mix=0.5, tol=1e-6, ctint_kwargs=dict(n_bath=nb),
                    insulating_start=(insulating and S is None),
                    sigma_init=None if S is None else S - Uo / 2 + U / 2)
        x, S, Uo = d._bath_x, r["Sigma"], U
        out.append((U, float(r["Z"]), r["iterations"]))
        print(nb, "ins" if insulating else "met", U, round(r["Z"], 5), r["iterations"], f"{time.time()-t:.0f}s", flush=True)
    return out
met = branch(list(np.round(np.arange(2.0, 3.61, 0.05), 3)), False) if which != "insulator" else [(np.nan, np.nan, 0)]
ins = branch(list(np.round(np.arange(3.6, 2.19, -0.05), 3)), True) if which != "metal" else [(np.nan, np.nan, 0)]
Um = np.array([m[0] for m in met]); Zm = np.array([m[1] for m in met])
Ui = np.array([m[0] for m in ins]); Zi = np.array([m[1] for m in ins])
km = int(np.argmax(Zm < 0.005)) if np.any(Zm < 0.005) else None
ki = int(np.argmax(Zi > 0.005)) if np.any(Zi > 0.005) else None
uc2 = 0.5 * (Um[km] + Um[km - 1]) if km else float("nan")    # metal last seen / first lost
# linear-regime extrapolation (window fixed in advance: 0.04 <= Z <= 0.15)
w = (Zm >= 0.04) & (Zm <= 0.15)
uc2_lin = float(-np.polyfit(Um[w], Zm[w], 1)[1] / np.polyfit(Um[w], Zm[w], 1)[0]) if w.sum() >= 3 else float("nan")
uc1 = 0.5 * (Ui[ki] + Ui[ki - 1]) if ki else float("nan")
json.dump(dict(n_bath=nb, metal=met, insulator=ins, Uc2=uc2, Uc2_linear=uc2_lin, Uc1=uc1),
          open(os.path.join(os.path.dirname(os.path.abspath(__file__)), f"scan_nb{nb}{'' if which == 'both' else '_' + which}.json"), "w"), indent=1)
print("RESULT", nb, "Uc1", uc1, "Uc2(collapse)", uc2, "Uc2(linear extrap.)", uc2_lin)
