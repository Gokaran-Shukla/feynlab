"""Charge self-consistent DFT+DMFT for SrVO3 with Quantum ESPRESSO + Wannier90 + feynlab (requires pw.x >= 7.0
built with HDF5, pw2wannier90.x, wannier90.x; e.g. micromamba create -p /opt/qe -c conda-forge qe wannier90).

Pseudopotentials: PBE ultrasoft from PSlibrary generated with ld1.x (github.com/dalcorso/pslibrary):
Sr.pbe-spn-rrkjus_psl.0.2.3, V.pbe-spnl-rrkjus_psl.1.0.0, O.pbe-n-rrkjus_psl.0.1 (copy the UPF files here).

    python run_srvo3.py [U] [J]       (defaults U = 5 eV, J = 0.75 eV, beta = 40 1/eV)
Steps: (1) U = 0 plumbing test (Delta N must vanish and pw.x must reproduce the DFT energy), (2) one-shot DMFT,
(3) charge self-consistent cycle.  Results -> RESULTS.md, history.json.
"""
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", ".."))
from feynlab.dft.qe_csc import QEDMFT  # noqa: E402

U = float(sys.argv[1]) if len(sys.argv) > 1 else 5.0
J = float(sys.argv[2]) if len(sys.argv) > 2 else 0.75
beta = 40.0
degauss_Ry = 1.0 / beta / 13.605693122994
a = 7.2612 * 0.529177210903
SCF = f""" &control
    calculation='scf', prefix='svo', pseudo_dir='./', outdir='./tmp/'
 /
 &system
    ibrav=1, celldm(1)=7.2612, nat=5, ntyp=3, ecutwfc=40.0, ecutrho=320.0,
    occupations='smearing', smearing='fermi-dirac', degauss={degauss_Ry:.8f}
 /
 &electrons
    mixing_beta=0.5, conv_thr=1e-10
 /
ATOMIC_SPECIES
 Sr 87.62  Sr.pbe-spn-rrkjus_psl.0.2.3.UPF
 V  50.9415 V.pbe-spnl-rrkjus_psl.1.0.0.UPF
 O  15.9994 O.pbe-n-rrkjus_psl.0.1.UPF
ATOMIC_POSITIONS (crystal)
 Sr 0.0 0.0 0.0
 V  0.5 0.5 0.5
 O  0.5 0.5 0.0
 O  0.5 0.0 0.5
 O  0.0 0.5 0.5
K_POINTS automatic
 8 8 8 0 0 0
"""
atoms = [("Sr", (0, 0, 0)), ("V", (.5, .5, .5)), ("O", (.5, .5, 0)), ("O", (.5, 0, .5)), ("O", (0, .5, .5))]
np_ = os.environ.get("NP", "2")
kw = dict(nbnd=30, n_iw=1024, mixing_beta=0.3, pw_cmd=f"mpirun -np {np_} pw.x -nk {np_}",
          pw2wan_cmd=f"mpirun -np {np_} pw2wannier90.x")
q0 = QEDMFT("u0", SCF, (8, 8, 8), (21, 23), ["V:dxz;dyz;dxy"], U=0.0, J=0.0, beta=beta, solver="hubbard1", **kw)
for f in os.listdir("."):
    if f.endswith(".UPF"):
        os.makedirs("u0", exist_ok=True); os.makedirs("csc", exist_ok=True)
        os.system(f"cp {f} u0/; cp {f} csc/")
q0.write_inputs(); q0.write_wannier(np.eye(3) * a, atoms); q0.run_scf(); q0.run_nscf(); q0.run_wannier()
i0 = q0.dmft_step(2)
s0 = q0.run_dmft_scf()
e_dft = float(open("u0/scf.out").read().split("!    total energy")[1].split("=")[1].split()[0])
q = QEDMFT("csc", SCF, (8, 8, 8), (21, 23), ["V:dxz;dyz;dxy"], U=U, J=J, beta=beta, solver="lanczos",
           solver_kwargs=dict(n_bath=2), **kw)
h = q.cycle(n_csc=12, n_dmft_first=45, n_dmft_per=10, lattice_rows_A=np.eye(3) * a, atoms_frac=atoms)
lines = ["# SrVO3: charge self-consistent DFT+DMFT with Quantum ESPRESSO 7.x + Wannier90 + feynlab", "",
         f"U = {U} eV, J = {J} eV (Kanamori), beta = {beta} 1/eV, t2g window = bands 21-23, 8x8x8 k mesh, "
         "solver: finite-T Lanczos ED, 2 bath sites per orbital", "",
         "## U = 0 plumbing test", "",
         f"* max |Delta N| = {i0['dN_max']:.2e}, mu_DMFT - E_F = {i0['mu_dmft']:.2e} eV",
         f"* total energy after the DMFT density update: {s0['etot_Ry']:.8f} Ry (DFT: {e_dft:.8f} Ry), "
         f"scf accuracy {s0['scf_accuracy_Ry']:.1e} Ry", "",
         "## CSC cycle", "", "| cycle | Z (t2g) | QE E_F (eV) | scf accuracy (Ry) | total energy (Ry) | max abs(Delta N) |",
         "|---|---|---|---|---|---|"]
for i, x in enumerate(h):
    lines.append(f"| {i} | {np.mean(x['Z']):.4f} | {x['scf_ef']:.4f} | {x['scf_scf_accuracy_Ry']:.2e} | {x['scf_etot_Ry']:.8f} | {x['dN_max']:.3f} |")
open("RESULTS.md", "w").write("\n".join(lines) + "\n")
json.dump([{k: (np.asarray(v).tolist() if not isinstance(v, (str, bool)) else v) for k, v in d.items()} for d in h],
          open("history.json", "w"), indent=1)
