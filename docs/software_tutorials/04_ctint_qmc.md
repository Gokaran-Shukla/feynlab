# Tutorial 4 — CT-INT: quantum Monte Carlo as a sum of Feynman diagrams

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


## The idea

Expand the partition function in the interaction:

Z/Z₀ = Σₖ (−U)ᵏ/k! ∫dτ₁…dτₖ ⟨T Hᵢₙₜ(τ₁)…Hᵢₙₜ(τₖ)⟩₀.

By Wick's theorem the thermal average at fixed vertex times is a sum over all k!-contractions =
**a determinant** (one per spin). So each Monte Carlo configuration (k, {τᵢ}, {sᵢ}) represents the
*sum of all Feynman diagrams* (connected and disconnected) with those vertices. The Metropolis walk
samples the perturbation order k and the times, i.e. it samples the diagrammatic series to all orders
(Rubtsov et al. PRB 72, 035122 (2005); Gull et al. RMP 83, 349 (2011)).

The sign problem is removed for repulsive U by the auxiliary shift
U n↑n↓ → (U/2) Σ_{s=±} (n↑ − α↑(s))(n↓ − α↓(s)), α_σ(s) = ½ + σ s δ, δ ≳ ½ (default 0.51).

## Using the solver

```python
import numpy as np
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.qmc.ctint import CTINTSolver
from feynlab.qmc.ed import AndersonED

beta, U = 10.0, 3.0
ed = AndersonED(eps_d=-U/2 + 0.3, U=U, eps_bath=[-1.0, 0.0, 1.2], V_bath=[0.5, 0.4, 0.6])
g = MatsubaraGrid(beta, 1024)
G0 = 1/(g.iw - ed.eps_d - ed.hybridization(g.iw))          # Weiss field of the *physical* problem

res = CTINTSolver(g, U, n_iw_meas=10).solve(G0, n_meas=1_000_000, n_chains=4)
print(res['sign'], res['k_mean'])                          # 1.0,  ~12
print(res['Sigma_up'][:3]); print(res['err_Sigma'][:3])
print(ed.sigma(res['iw'][:3], beta))                        # exact reference
```

Returned: `G_up/G_dn`, `Sigma_up/dn`, error bars, average sign, ⟨k⟩ and the **order histogram**
`order_histogram[k]` — the probability that the sampled diagram order is k.

## Reading the order histogram
* ⟨k⟩ grows ~ βU: lower temperature or larger U → higher orders → cost ∝ ⟨k⟩³.
* A histogram touching `kmax` (default 400) means you must increase it.
* Diagnostic identity: ⟨k⟩ = −β⟨H_int,shifted⟩, so ⟨k⟩ also measures the potential energy.

## Validation built into feynlab
| test | expected |
|---|---|
| Hubbard atom (no bath) | Σ = U/2 + U²/(4iω) |
| Anderson impurity with 3 bath sites, no particle–hole symmetry | ED within error bars (benchmarks) |
| DMFT at U = D | Z agrees with IPT (exact to O(U²)) |

## Practical tips
* Use `n_chains ≥ 4` to get honest error bars (they are inter-chain standard errors).
* The measurement is done in Matsubara space directly (no τ binning) for the first `n_iw_meas` frequencies.
* The Hubbard atom is the *hardest* case (no hybridization → long autocorrelation); real baths mix faster.

## v0.2: two more diagrammatic Monte Carlo methods

### Fröhlich polaron (Prokof'ev–Svistunov DiagMC)
```python
from feynlab.qmc.polaron import polaron_diagmc, feynman_energy
r = polaron_diagmc(alpha=1.0, mu=feynman_energy(1.0) - 0.15, tau_max=40,
                   n_steps=20_000_000, n_chains=4, fit_window=(0.4, 0.97))
print(r['E'], r['E_err'], r['order_mean'])    # ~ -1.017 (Hahn et al. PRB 2018: -1.01662(47))
```
Every configuration is one Feynman diagram of G(p=0,τ) with an arbitrary number of phonon arcs; updates add/
remove arcs and resample the diagram length. E₀ comes from the decay G ~ Z e^{−(E−μ)τ}.

### Bold-line (skeleton) series
```python
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.diagrams.bold import skeleton_series_at
g = MatsubaraGrid(1.0, 512); U = 6.0
Gex = 0.5*(1/(g.iw + U/2) + 1/(g.iw - U/2))         # exact atomic G
ser = skeleton_series_at(g, Gex, U, max_order=4, n_iw=1)
```
For βU = 6 the skeleton series converges — but to a self-energy far from the exact U²/(4iω): the
"misleading convergence" of Kozik, Ferrero & Georges (PRL 114, 156402 (2015)). For βU = 0.5 it converges to the
exact result. `bold_selfconsistent(...)` solves G = [G₀⁻¹ − Σ_skeleton[G]]⁻¹ up to a chosen order.
