"""Collect ED-DMFT Mott scans (scan_nb*.json) -> table, figure, verdicts.

Physics of the bath-size parity (T = 0, particle-hole symmetric bath):
* odd n_bath  -> even number of cluster sites, singlet ground state with a zero-energy bath level:
                 describes the METAL (Kondo quasiparticle) well -> use for U_c2.
                 Its insulator acquires a spurious Kondo resonance on the zero-energy level.
* even n_bath -> odd number of sites, spin-doublet ground state, no zero mode:
                 describes the INSULATOR well -> use for U_c1; its metal collapses too early.
U_c2: linear extrapolation of Z(U) in the window 0.04 <= Z <= 0.15 (fixed in advance).
U_c1: midpoint between the last U with Z_ins < 0.005 and the first with Z > 0.005.
"""
import glob, json, os
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REF_UC1, REF_UC1_ERR = 2.39, 0.02          # Karski, Raas, Uhrig PRB 72, 113110 (2005) and refs.
REF_UC2_LO, REF_UC2_HI = 2.92, 3.00        # range of published values (NRG, ED, DMRG, PT)


def load():
    data = {}
    for f in sorted(glob.glob(os.path.join(HERE, "scan_nb*.json"))):
        d = json.load(open(f))
        nb = d["n_bath"]
        e = data.setdefault(nb, dict(metal=[], insulator=[]))
        for k in ("metal", "insulator"):
            pts = [p for p in d[k] if p[0] == p[0]]
            if pts:
                e[k] = pts
    return data


def uc2_linear(met):
    U, Z = np.array([m[0] for m in met]), np.array([m[1] for m in met])
    w = (Z >= 0.04) & (Z <= 0.15)
    if w.sum() < 3:
        return np.nan, np.nan
    p = np.polyfit(U[w], Z[w], 1)
    return -p[1] / p[0], p[0]


def uc2_collapse(met):
    U, Z = np.array([m[0] for m in met]), np.array([m[1] for m in met])
    k = np.argmax(Z < 0.002) if np.any(Z < 0.002) else None
    return 0.5 * (U[k] + U[k - 1]) if k else np.nan


def uc1(ins):
    U, Z = np.array([m[0] for m in ins]), np.array([m[1] for m in ins])
    k = np.argmax(Z > 0.005) if np.any(Z > 0.005) else None
    return 0.5 * (U[k] + U[k - 1]) if k else np.nan


if __name__ == "__main__":
    data = load()
    rows = []
    for nb, e in sorted(data.items()):
        r = dict(n_bath=nb, parity="odd" if nb % 2 else "even")
        r["Uc2_linear"], r["slope"] = uc2_linear(e["metal"]) if e["metal"] else (np.nan, np.nan)
        r["Uc2_collapse"] = uc2_collapse(e["metal"]) if e["metal"] else np.nan
        r["Uc1"] = uc1(e["insulator"]) if e["insulator"] else np.nan
        rows.append(r)
    lines = ["| n_bath | parity | U_c2 (linear Z->0) | U_c2 (metal collapse) | U_c1 (insulator lost) |", "|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['n_bath']} | {r['parity']} | {r['Uc2_linear']:.3f} | {r['Uc2_collapse']:.3f} | {r['Uc1']:.3f} |")
    lines.append(f"| ref | | {REF_UC2_LO}-{REF_UC2_HI} | | {REF_UC1} ± {REF_UC1_ERR} |")
    txt = "\n".join(lines)
    print(txt)
    open(os.path.join(HERE, "TABLE.md"), "w").write(txt + "\n")
    json.dump(rows, open(os.path.join(HERE, "summary.json"), "w"), indent=1, default=float)

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    C = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
    fig, ax = plt.subplots(figsize=(7.2, 4.4), dpi=150)
    for i, (nb, e) in enumerate(sorted(data.items())):
        if e["metal"]:
            U, Z = np.array([m[0] for m in e["metal"]]), np.array([m[1] for m in e["metal"]])
            ax.plot(U, Z, "-o", ms=3.5, lw=1.8, color=C[i % 5], label=f"{nb} bath, metal start")
        if e["insulator"]:
            U, Z = np.array([m[0] for m in e["insulator"]]), np.array([m[1] for m in e["insulator"]])
            ax.plot(U, Z, "--s", ms=3.5, lw=1.2, color=C[i % 5], alpha=0.8, label=f"{nb} bath, insulator start")
    ax.axvspan(REF_UC2_LO, REF_UC2_HI, color="#9aa3ad", alpha=0.25, lw=0)
    ax.axvspan(REF_UC1 - REF_UC1_ERR, REF_UC1 + REF_UC1_ERR, color="#9aa3ad", alpha=0.25, lw=0)
    ax.annotate("published U_c1", (REF_UC1, 0.26), ha="center", fontsize=8, color="#6b6a64")
    ax.annotate("published U_c2", ((REF_UC2_LO + REF_UC2_HI) / 2, 0.26), ha="center", fontsize=8, color="#6b6a64")
    ax.set_xlabel("U / D"); ax.set_ylabel("quasiparticle weight Z")
    ax.set_title("T = 0 ED-DMFT (Lanczos), Bethe lattice: bath-size dependence", fontsize=10, loc="left")
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.grid(True, color="#e6e5e0", lw=0.6); ax.set_axisbelow(True)
    ax.legend(frameon=False, fontsize=7, ncol=2); ax.set_ylim(-0.01, 0.3)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "ed_dmft_mott_bath_size.png"))
