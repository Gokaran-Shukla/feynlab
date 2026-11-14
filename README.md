# feynlab — a Feynman-diagram laboratory

Companion software to R. D. Mattuck, *A Guide to Feynman Diagrams in the Many-Body Problem* (2nd ed.),
extended to research-level tools for correlated electrons and Hall transport.

**Status v1.0.0** — every check in `benchmarks/RESULTS.md` (102 core rows) and `benchmarks/RESULTS_v1.md` (54 rows, 54/54 pass) compares against an exact result, an independent solver or a published formula (PRL, PRB, RMP, EPL, Phys. Rep.).

📑 **Slides:** [`docs/slides/feynlab_slides.pdf`](docs/slides/feynlab_slides.pdf) — a 46-page Beamer deck covering every module, equation (colour-coded by physical quantity), benchmark and limitation.

## What it does

| capability | method | validated against |
|---|---|---|
| Diagram enumeration & analysis | bijection representation; connected / proper / skeleton / tadpole / loops / sign / topology / drawing | Mohr (2, 10, 74, 706); Molinari PRB 71, 113102 (loop-resolved) |
| Diagrammatic Monte Carlo (bare series) | MC over internal imaginary times, Feynman rules (−U)ⁿ(−1)^F∏G₀ | Hubbard atom (orders 3, 4 cancel); ED |
| **Polaron DiagMC** (v0.2) | Prokof'ev–Svistunov all-order sampling of G(p=0,τ); Feynman variational energy | Hahn et al. PRB 97, 134305 (2018); Feynman 1955 |
| **Bold-line / skeleton series** (v0.2) | skeleton diagrams with dressed G, self-consistent bold scheme | Kozik, Ferrero, Georges PRL 114, 156402 (2015) |
| CT-INT quantum Monte Carlo | determinant (all-diagrams-per-order) sampling, Rubtsov/Gull | exact diagonalization, sign = 1 |
| Self-energy, **weak** | Hartree–Fock, 2nd order (impurity & 2D lattice), RPA/FLEX spin+charge fluctuations, GF2 | Yamada–Yosida; exact tail moment |
| Self-energy, **moderate** | DMFT + IPT (Mott hysteresis); **T = 0 Lanczos ED-DMFT** (5–9 bath sites; U_c2 = 2.955D, U_c1 = 2.375D); any lattice (`LatticeDMFT`); MaxEnt | U_c1, U_c2 consensus (Bulla PRL 1999; Karski et al. PRB 2005) |
| Self-energy, **strong** | Hubbard-I, atomic limit, DMFT + CT-INT | atomic limit; IPT at weak U |
| Electron gas | Lindhard (iω & real ω), plasmons, RPA ring-sum correlation energy, Friedel | Mattuck Ex. 13.6 (3/5); Gell-Mann–Brueckner |
| Superconductivity | BCS gap equation, T_c | 1.764, 1.134 |
| **Hall transport** (charge, spin, orbital, orbital-quadrupole) | Berry / operator curvature; Kubo–Středa NCA (Born Σ + ladder + V₁³ skew vertex); **crossed X/Ψ diagrams** (v0.2); Boltzmann (side jump, anomalous distribution, skew) | Sinitsyn et al. PRB 75, 045315 (all 6 mechanism terms, ≤1e-5); Ado et al. EPL 111, 37004 (X, Ψ, full σ_xy); Inoue–Bauer–Molenkamp; QWZ/Haldane Chern; Rashba e/8π |
| Models | massive Dirac (+tilt, +quadratic mass), Rashba(+exchange), QWZ, Haldane, p-orbital square lattice with SOC/orbital Rashba/exchange, **Wannier90 hr.dat (2D and 3D)** | — |

### New in v1.0 (all ROADMAP items)
| capability | module | validated against |
|---|---|---|
| **CT-HYB** (segment picture) | `qmc.cthyb`, `BetheDMFT.solve('cthyb')` | ED: n, ⟨n↑n↓⟩, Σ(iω₀) at U = 3, 4 |
| **Kanamori multi-orbital** atom, Hubbard-I, sparse ED | `qmc.multiorbital` | multiplets U−3J, U+J, t₂g³ 3U−9J (Georges et al. 2013) |
| **DFT+DMFT** from Wannier90 (FLL / AMF / occupation dc) | `dmft.dftdmft` | U = 0 limit; Mott insulator |
| **DCA** (2×2 clusters) | `dmft.dca` | N_c = 1 ≡ DMFT to 6e-6; antinodal suppression |
| **Two-particle**: bubble, RPA χ_s, χ_c | `twoparticle` | compressibility sum rule 1e-8 |
| **NRG** (Wilson chain), Kondo thermodynamics | `qmc.nrg` | resonant level; ln 2 plateau; Haldane slope −0.993 |
| **Nevanlinna** and **stochastic** continuation | `continuation` | exact interpolation; sum rules (Fei, Yeh, Gull 2021) |
| **Connected-determinant DiagMC** | `qmc.cdet` | Hubbard atom (exact), Anderson vs ED |
| **Resummation** (Domb–Sykes, Padé, conformal), misleading-convergence check | `diagrams.resum` | exact series |
| **Polaron Z and m\*** | `qmc.polaron` | 1 − α/2; 1 + α/6 + 0.0236α² |
| **Matrix impurities and T-matrix** | `transport.kubo_streda` | Born + skew (Sinitsyn 2007) |
| **Lattice crossed X/Ψ** (BZ FFT) | `transport.crossed` | k·p continuum limit |
| **Modern orbital moment**, orbital magnetisation, **Berry-curvature dipole**, **Edelstein/SOT**, **finite-T Kubo–Bastin** | `transport.response` | Xiao–Yao–Niu 2007; Sodemann–Fu 2015; Edelstein 1990 |
| **3D Wannier** AHC/SHC/OHC with adaptive k refinement | `transport.wannier3d` | quantised Chern stacks (Halperin 1987) |
| Γ-method / jackknife, Anderson/Broyden mixing, HDF5 checkpoints, TOML CLI, CuPy backend, units, weekly CI, CITATION.cff | `utils`, `cli` | AR(1) τ_int; same fixed point |

## Quick start
```bash
pip install -e .
pytest -q tests
python benchmarks/run_benchmarks.py          # core suite
python benchmarks/run_v1_benchmarks.py       # v1.0 modules (~10 min)
feynlab bench --suite all                    # both, via the CLI
feynlab hall --model dirac                   # AHE mechanisms from the command line
open quizzes/quiz.html                        # 36 quizzes, 187 questions
open docs/slides/feynlab_slides.pdf           # presentation
```

```python
from feynlab.transport.models import MassiveDirac
from feynlab.transport.mechanisms import anomalous_hall_mechanisms
print(anomalous_hall_mechanisms(MassiveDirac(1.0, 0.5), mu=1.0, band=1, V1_over_V0=0.3))
```

## Documentation
* `docs/software_tutorials/01…14` — using the package (09–14 cover the v1.0 modules)
* `docs/slides/feynlab_slides.pdf` — presentation slides (source `feynlab_slides.tex`, figures via `make_figures.py`)
* `docs/book_tutorials/ch00…ch19, appendices` — one tutorial per chapter of Mattuck: goals, key ideas, guide to the book's exercises, feynlab hands-on, self-check
* `quizzes/quiz.html` (interactive), `quizzes/QUIZZES.md` (printable with answer keys), `quizzes/quizzes.json`

## Roadmap
`docs/ROADMAP.md` lists the 22 roadmap items with their v1.0 status (all implemented; partial items are marked) and the next steps.

## Known limitations
See `docs/software_tutorials/07_benchmarks_and_limitations.md`: CT-HYB is density–density only; DFT+DMFT has no charge
self-consistency; the DCA cluster solver is second order; Nevanlinna lacks Hardy optimisation; stochastic continuation is a
cross-check only; two-particle = bubble/RPA (no vertex ladders); polaron m* has a few-% statistical spread; the semiclassical
Gaussian-skew rate decomposition is partial (the Kubo route is exact).

## Citation
If you use feynlab, please cite the original papers of the methods you use (listed in each module docstring).
