"""feynlab command-line interface driven by TOML input files.

    feynlab dmft  input.toml      # DMFT (Bethe or lattice; ipt/ed/ctint/hubbard1/ctint)
    feynlab hall  input.toml      # intrinsic / disordered Hall responses
    feynlab polaron input.toml    # Froehlich polaron DiagMC
    feynlab bench [--quick]       # benchmark suite
    feynlab nrg   input.toml      # NRG impurity thermodynamics, Kondo temperature
    feynlab bench --suite core|v1|all
    feynlab example dmft|hall|polaron|nrg   # print an example input file
    (without an input file the example input is used)

Results go to the HDF5 file named in [output] (default feynlab_out.h5).
"""
from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

try:
    import tomllib
except ImportError:  # pragma: no cover
    import tomli as tomllib

import numpy as np

EXAMPLES = {
    "dmft": """# feynlab dmft example
[model]
lattice = "bethe"        # bethe | square | cubic | hypercubic
D = 1.0
U = 2.6
beta = 200.0
[solver]
name = "ed"              # ipt | ed | ctint | hubbard1
n_bath = 7
n_iter = 150
mix = 0.5
[output]
file = "dmft_U2.6.h5"
""",
    "hall": """# feynlab hall example
[model]
name = "MassiveDirac"    # MassiveDirac | Rashba2DEG | QWZ | Haldane | POrbitalSquare
params = {v = 1.0, Delta = 0.5}
[calc]
mu = 1.0
band = 1
mechanisms = true        # intrinsic + side jump + skew (+ crossed if include_crossed)
include_crossed = false
V1_over_V0 = 0.3
[output]
file = "hall.h5"
""",
    "polaron": """[model]
alpha = 1.0
[mc]
n_steps = 20000000
n_chains = 4
tau_max = 40
[output]
file = "polaron.h5"
""",
}


def run_dmft(cfg):
    from ..dmft.bethe import BetheDMFT, LatticeDMFT
    m, s = cfg["model"], cfg.get("solver", {})
    lat = m.get("lattice", "bethe")
    if lat == "bethe":
        d = BetheDMFT(m["beta"], m["U"], D=m.get("D", 1.0), n_iw=m.get("n_iw", 1024))
    else:
        d = getattr(LatticeDMFT, lat)(m["beta"], m["U"])
    kw = {k: v for k, v in s.items() if k in ("n_bath", "beta_fit", "n_fit", "n_meas", "n_chains", "n_iw_meas")}
    r = d.solve(s.get("name", "ipt"), n_iter=s.get("n_iter", 200), mix=s.get("mix", 0.5),
                ctint_kwargs=kw or None, insulating_start=s.get("insulating_start", False))
    return dict(iw=r["iw"], G=r["G"], Sigma=r["Sigma"], Z=r["Z"], n=r["n"], iterations=r["iterations"])


def run_hall(cfg):
    from ..transport import models
    from ..transport.mechanisms import anomalous_hall_mechanisms
    from ..transport.intrinsic import intrinsic_hall
    m, c = cfg["model"], cfg["calc"]
    model = getattr(models, m["name"])(**m.get("params", {}))
    if c.get("mechanisms", False):
        out = anomalous_hall_mechanisms(model, c["mu"], c.get("band", 1), V1_over_V0=c.get("V1_over_V0", 0.0),
                                        include_crossed=c.get("include_crossed", False))
    else:
        op = model.operators.get(c["operator"]) if c.get("operator") else None
        out = dict(sigma=intrinsic_hall(model, c["mu"], op=op, nk=c.get("nk", 300)))
    return out


def run_polaron(cfg):
    from ..qmc.polaron import polaron_diagmc, feynman_energy
    a = cfg["model"]["alpha"]
    mc = cfg.get("mc", {})
    r = polaron_diagmc(a, mu=feynman_energy(a) - 0.15, tau_max=mc.get("tau_max", 40),
                       n_steps=mc.get("n_steps", 20_000_000), n_chains=mc.get("n_chains", 4), fit_window=(0.4, 0.97))
    return dict(E=r["E"], E_err=r["E_err"], Z=r["Z"], Z_err=r["Z_err"], tau=r["tau"], G=r["G"], feynman=feynman_energy(a))


def run_nrg(cfg):
    from ..qmc.nrg import impurity_thermodynamics, kondo_temperature
    m, c = cfg["model"], cfg.get("nrg", {})
    T, S, X = impurity_thermodynamics(U=m["U"], Gamma=m["Gamma"], eps_d=m.get("eps_d"), Lam=c.get("Lambda", 2.5),
                                      n_iter=c.get("n_iter", 50), n_keep=c.get("n_keep", 300))
    return dict(T=T, S_imp=S, Tchi_imp=X, T_K=kondo_temperature(T, X) if m["U"] > 0 else np.nan)


EXAMPLES["nrg"] = """[model]
U = 0.1
Gamma = 0.01
[nrg]
Lambda = 2.5
n_iter = 60
n_keep = 300
[output]
file = "nrg.h5"
"""


def build_parser():
    p = argparse.ArgumentParser(prog="feynlab", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="command", required=True)
    for name in ("dmft", "hall", "polaron", "nrg"):
        sp = sub.add_parser(name)
        sp.add_argument("input", nargs="?", default=None)
        sp.add_argument("--model", default=None, help="(hall) shortcut: dirac | rashba | qwz")
    b = sub.add_parser("bench")
    b.add_argument("--quick", action="store_true")
    b.add_argument("--suite", choices=("core", "v1", "all"), default="core")
    e = sub.add_parser("example")
    e.add_argument("what", choices=sorted(EXAMPLES))
    return p


def main(argv=None):
    a = build_parser().parse_args(argv)
    a.cmd = a.command
    if a.cmd == "example":
        print(EXAMPLES[a.what])
        return 0
    if a.cmd == "bench":
        bdir = os.path.join(os.path.dirname(__file__), "..", "..", "benchmarks")
        rc = 0
        if a.suite in ("core", "all"):
            rc |= subprocess.call([sys.executable, os.path.join(bdir, "run_benchmarks.py")] + (["--quick"] if a.quick else []))
        if a.suite in ("v1", "all"):
            rc |= subprocess.call([sys.executable, os.path.join(bdir, "run_v1_benchmarks.py")])
        return rc
    if a.input is None:
        if a.cmd == "hall" and a.model:
            name = {"dirac": "MassiveDirac", "rashba": "Rashba2DEG", "qwz": "QWZ"}[a.model]
            cfg = tomllib.loads(EXAMPLES["hall"].replace('"MassiveDirac"', f'"{name}"'))
            if name != "MassiveDirac":
                cfg["model"]["params"] = {}
                cfg["calc"]["mechanisms"] = False
        else:
            cfg = tomllib.loads(EXAMPLES[a.cmd])
    else:
        with open(a.input, "rb") as f:
            cfg = tomllib.load(f)
    res = {"dmft": run_dmft, "hall": run_hall, "polaron": run_polaron, "nrg": run_nrg}[a.cmd](cfg)
    from ..utils.io import save
    out = cfg.get("output", {}).get("file", "feynlab_out.h5")
    save(out, {k: v for k, v in res.items()} | {"input_json": json.dumps(cfg, default=str)})
    print(json.dumps({k: (float(v) if np.ndim(v) == 0 and not isinstance(v, str) else f"array{np.shape(v)}")
                      for k, v in res.items()}, indent=1))
    print("results written to", out)
    return 0


if __name__ == "__main__":
    sys.exit(main())
