# feynlab — a Feynman-diagram laboratory

Companion software to R. D. Mattuck, *A Guide to Feynman Diagrams in the Many-Body Problem* (2nd ed.),
extended to research-level tools for correlated electrons and Hall transport.

**Status v0.3.0** — every check in `benchmarks/RESULTS.md` compares against an exact or published result (PRL, PRB, RMP, EPL, Phys. Rep.).

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
| Models | massive Dirac, Rashba(+exchange), QWZ, Haldane, p-orbital square lattice with SOC/orbital Rashba/exchange, **Wannier90 hr.dat** | — |

## Quick start
```bash
pip install -e .
pytest -q tests
python benchmarks/run_benchmarks.py          # ~2 min with numba cache
open quizzes/quiz.html                        # 28 quizzes, 140 questions
```

```python
from feynlab.transport.models import MassiveDirac
from feynlab.transport.mechanisms import anomalous_hall_mechanisms
print(anomalous_hall_mechanisms(MassiveDirac(1.0, 0.5), mu=1.0, band=1, V1_over_V0=0.3))
```

## Documentation
* `docs/software_tutorials/01…07` — using the package (every code block is executed in CI-style checks)
* `docs/book_tutorials/ch00…ch19, appendices` — one tutorial per chapter of Mattuck: goals, key ideas, guide to the book's exercises, feynlab hands-on, self-check
* `quizzes/quiz.html` (interactive), `quizzes/QUIZZES.md` (printable with answer keys), `quizzes/quizzes.json`

## Roadmap
See `docs/ROADMAP.md` for the planned extensions (CT-HYB, multi-orbital Kanamori, cluster DMFT, CDet DiagMC, Wannier 3D orbital transport with crossed diagrams, GPU/HDF5/CI).

## Known limitations
See `docs/software_tutorials/07_benchmarks_and_limitations.md` (what v0.2 fixed, and what remains:
the semiclassical rate-by-rate Gaussian skew decomposition, slow bath-size convergence of U_c1 in ED,
~0.2% polaron precision, single-site DMFT, bath-parity rules for ED-DMFT).

## Citation
If you use feynlab, please cite the original papers of the methods you use (listed in each module docstring).
