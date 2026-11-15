"""Turn a QEDMFT CSC log (csc.log) into a convergence figure and a markdown table.

usage: python benchmarks/qe_csc/make_report.py csc.log [figure.png]
"""
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PAT = re.compile(r"CSC (\d+): n_win=([\d.]+).*?Z=\[([\d. ]+)\].*?\|dN\|max=([\de.+-]+)\s+E_F=([\d.]+) eV"
                 r"\s+scf_acc=([\de.+-]+) Ry\s+E=([-\d.]+) Ry")


def parse(path):
    rows = []
    for line in open(path):
        m = PAT.search(line)
        if m:
            rows.append(dict(it=int(m[1]), n_win=float(m[2]), Z=float(m[3].split()[0]), dN=float(m[4]),
                             ef=float(m[5]), acc=float(m[6]), E=float(m[7])))
    return rows


def table(rows):
    out = ["| cycle | Z (t2g) | E_F (eV) | QE scf accuracy (Ry) | E (Ry) | ΔE (Ry) | window charge |",
           "|---|---|---|---|---|---|---|"]
    prev = None
    for r in rows:
        dE = "—" if prev is None else f"{r['E'] - prev:+.1e}"
        out.append(f"| {r['it']} | {r['Z']:.3f} | {r['ef']:.4f} | {r['acc']:.2e} | {r['E']:.6f} | {dE} | {r['n_win']:.5f} |")
        prev = r["E"]
    return "\n".join(out)


def figure(rows, path):
    C, MUTED = ["#2a78d6", "#eb6834", "#1baf7a"], "#6b6a64"
    it = [r["it"] for r in rows]
    fig, ax = plt.subplots(1, 3, figsize=(9.6, 2.8))
    ax[0].semilogy(it, [r["acc"] for r in rows], "o-", color=C[0], lw=2, ms=6)
    ax[0].set_title("QE scf accuracy (Ry)", fontsize=9, color="#1f1f1e")
    E0 = rows[-1]["E"]
    ax[1].semilogy(it[:-1], [abs(r["E"] - E0) * 1e3 for r in rows[:-1]], "o-", color=C[1], lw=2, ms=6)
    ax[1].set_title(r"$|E - E_{\rm last}|$ (mRy)", fontsize=9, color="#1f1f1e")
    ax[2].plot(it, [r["Z"] for r in rows], "o-", color=C[2], lw=2, ms=6)
    ax[2].set_title(r"quasiparticle weight $Z$ (t$_{2g}$)", fontsize=9, color="#1f1f1e")
    ax[2].set_ylim(0.40, 0.45)
    for a in ax:
        for s in ("top", "right"):
            a.spines[s].set_visible(False)
        for s in ("left", "bottom"):
            a.spines[s].set_color(MUTED)
        a.tick_params(colors=MUTED, labelsize=8)
        a.grid(True, color="#e6e5e0", lw=0.6)
        a.set_axisbelow(True)
        a.set_xlabel("CSC cycle", fontsize=8, color=MUTED)
    fig.tight_layout()
    fig.savefig(path, dpi=200)


if __name__ == "__main__":
    rows = parse(sys.argv[1])
    print(table(rows))
    if len(sys.argv) > 2:
        figure(rows, sys.argv[2])
