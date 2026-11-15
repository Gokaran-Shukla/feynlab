# Tutorial 20 — Two-particle vertex for continuous baths (CT-HYB)

The Lehmann ED vertex (`vertex.ImpurityED.g2`) is exact, but it lives in a finite Fock space. The cost grows with
the number of 4-state chains: 2 bath sites (64 states) take minutes, and 3 bath sites (256 states) exhaust
memory. A DMFT bath discretized with 2 sites also shifts the lattice susceptibility. Segment CT-HYB samples the
continuous bath exactly, so the vertex comes for free from the same Markov chain.

## 1. Estimator

Removing two rows and two columns from the hybridization matrix gives the Wick-like estimator

    G2_{σσ'}(ν,ν',ω) = (1/β) ⟨ M̂_σ(ν+ω,ν) M̂_σ'(ν',ν'+ω) − δ_σσ' M̂_σ(ν+ω,ν'+ω) M̂_σ(ν',ν) ⟩

with M̂(a,b) = Σ_ij M_ji e^{iν_a τ^e_i} e^{−iν_b τ^s_j}. M̂ is built as A·Mᵀ·B with A = e^{iντ^e} and
B = e^{−iντ^s}, then the two products are accumulated for all (ν, ν', ω) in the box. Up and down spins are
symmetrized (paramagnet).

At U = 0 the estimator reduces exactly to the bubble: χ_↑↓ = 0 and χ_↑↑ = χ0 δ_νν'.

## 2. Usage

```python
from feynlab.qmc.cthyb import CTHYBSolver
from feynlab.vertex import vertices, eom_sigma, lambda_dga
r = CTHYBSolver(grid, U, n_iw_meas=64).solve(delta_iw, mu, n_meas=2_000_000, meas_every=5, g2=(16, 15))
gs = r["gs"]                      # chi_uu, chi_ud, chi0, nu, om, G, G_of (+ chain errors)
v = vertices(gs, beta)            # F_d, F_m, Gamma_d, Gamma_m
eom_sigma(v["F_ud"], gs, beta, U, 1.0)    # Schwinger-Dyson check against r["Sigma"]
```

For a lattice, set `gs["G_of"]` to the DMFT lattice G_loc and call `lambda_dga` (tutorial 18).

## 3. Results (`benchmarks/vertex_cthyb/RESULTS.md`)

* **Discrete bath vs ED.** The deviations follow the statistical errors (median 0.7–1.0 σ). CT-HYB takes 16 s
  where ED takes 292 s.
* **2D Hubbard, U = 4t, T = 0.4t.** Σ from the equation of motion equals Σ from Dyson within 1–2%.
* **Box convergence of λ-DΓA.** Box 32 is out of reach for ED. Going from box 16 to box 32 moves the local sum
  towards the exact value (0.295 → 0.326, exact 0.354). The gap between the two sum-rule targets shrinks from a
  factor 2 to 1.4.

## 4. Limits

* **Spin structure.** Single orbital and paramagnetic (symmetrized). The multi-orbital `SectorCTHYB` has no G2
  measurement yet.
* **Estimator.** The removal estimator is noisy at high frequency and low expansion order. Worm sampling with
  improved estimators would reduce the noise.
