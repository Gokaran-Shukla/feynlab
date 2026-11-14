# Tutorial 1 — Getting started with feynlab

## Install

```bash
cd feynlab
pip install -e .            # numpy, scipy, numba, matplotlib
pytest -q tests             # ~15 fast regression tests
python benchmarks/run_benchmarks.py --quick   # validation against published results (~5 min)
python benchmarks/run_benchmarks.py           # full suite (~15-20 min)
```

## Package map

| subpackage | what it does | book chapters |
|---|---|---|
| `core.matsubara` | Matsubara grids, exact FFT-based τ↔iωₙ transforms with analytic tails | 14 |
| `core.lattice` | Bethe/square-lattice DOS, Hilbert transforms, Padé continuation | 9, App. H/L |
| `diagrams` | enumerate / classify / draw / Monte-Carlo-evaluate Feynman diagrams | 4, 5, 9, 10, 11 |
| `selfenergy` | HF, 2nd order, RPA/FLEX (lattice); Hubbard-I, atomic, Yamada–Yosida (local) | 4.7, 9.8, 10 |
| `qmc` | CT-INT diagrammatic QMC; exact diagonalization reference | — |
| `dmft` | Bethe-lattice DMFT with IPT / Hubbard-I / CT-INT solvers | — |
| `egas` | RPA electron gas (Lindhard, plasmons, correlation energy), BCS | 10, 12, 13, 15 |
| `transport` | Berry curvature, Kubo–Středa (ladder + skew), Boltzmann (side jump, skew) | — |
| `book` | calculators for the book's worked examples | 2, 3, 10, 18 |

## Conventions (memorise these — every module uses them)

* ħ = k_B = 1. Lattice constant 1. Electron-gas module: Hartree atomic units.
* Green's function **G(τ) = −⟨T c(τ) c†(0)⟩**, fermionic ωₙ = (2n+1)π/β, G(iω) = ∫₀^β e^{iωτ}G(τ)dτ.
* Dyson: **G = G₀ + G₀ Σ G**, i.e. G⁻¹ = G₀⁻¹ − Σ.
* Feynman rule for an order-n diagram (instantaneous interaction): factor (−v)ⁿ (−1)^F ∏ G₀, F = number of closed fermion loops.
* Hall conductivities: **e²/h** for charge currents, **e/2π** (i.e. (ħ/e)·e²/h) for spin/orbital/quadrupole currents with operator O dimensionless (S_z = σ_z/2, L_z in units of ħ).

## First calculation: a free level

```python
import numpy as np
from feynlab.core.matsubara import MatsubaraGrid, g0_tau_level

g = MatsubaraGrid(beta=20.0, n_iw=512)
eps = 0.7
G_iw = 1 / (g.iw - eps)
G_tau = g.iw_to_tau(G_iw, c1=1.0, c2=eps, c3=eps**2)   # tail coefficients of 1/(iw-eps)
print(np.abs(G_tau - g0_tau_level(eps, g.tau, g.beta)).max())   # ~1e-8
print(g.density(G_iw, c2=eps), 1/(np.exp(g.beta*eps)+1))         # occupation
```

Why the tail matters: G(iω) ~ 1/iω decays so slowly that a naive Matsubara sum has an O(1) error at τ→0⁺. feynlab subtracts c₁/iω + c₂/(iω)² + c₃/(iω)³, transforms the smooth remainder by FFT, and adds the analytic transforms of the tail back. The reverse transform integrates the piecewise-linear interpolant **exactly**, so the error is O(Δτ²) with no aliasing.

## Where next
* Tutorial 2 — Feynman diagrams as data structures
* Tutorial 3 — self-energy from weak to strong coupling
* Tutorial 4 — CT-INT quantum Monte Carlo
* Tutorial 5 — Hall transport: intrinsic, side jump, skew
* Tutorial 6 — electron gas, BCS and the book companions
* Tutorial 7 — benchmarks, accuracy and known limitations
