# Tutorial 13 — Connected determinants, resummation, bold-series diagnostics and polaron Z / m*

Roadmap items B8–B11.

## 1. Connected-determinant DiagMC (`feynlab.qmc.cdet`)

The sum of **all** connected diagrams of order n with fixed vertices is obtained from determinants of
non-interacting Green functions by inclusion–exclusion over subsets (Rossi, PRL 119, 045701 (2017)):

  c(V) = det(V) − Σ_{S ⊊ V, S ∋ v₁} c(S) det(V \ S)

* The cost is 3ⁿ instead of the n! of explicit diagrams.
* The fermionic sign cancellation between diagrams is done exactly inside each determinant.
* Vertex times are Monte Carlo sampled.

```python
from feynlab.qmc.cdet import cdet_series, atom_exact_coefficients
coef, err = cdet_series(grid, g0_iw, max_order=4, n_samples=400000)   # X(U) = <n_up n_dn> - n^2 = sum_k a_k U^k
```
Checks:

* Hubbard atom: matches the exact Taylor series of −¼ tanh(βU/4) at every order.
* Anderson impurity: a₁ = −0.06813 (CDet) vs −0.06792 (ED fit), and a₃ = 0.002225 vs 0.002246.

**Implementation note.** c(∅) must be defined as 0 for the recursion to start correctly.

## 2. Resummation (`feynlab.diagrams.resum`)

A DiagMC series Σ a_n Uⁿ is only useful if it can be summed beyond its radius of convergence R.

| tool | purpose | test result |
|---|---|---|
| `domb_sykes(a)` | ratio plot a_n/a_{n−1} vs 1/n → 1/R and the singularity exponent | R = 1.0000 for (1+x)^(−1/2) |
| `pade(a, L, M)` | rational approximant; continues past R if the singularity is a pole or branch cut | x = 3: 0.50000001 (exact 0.5) |
| `conformal_map_sum(a, R, x)` | maps the cut plane to the disc (Rossi et al. PRL 121, 130405 (2018)) | x = 3: 0.4999999996 |
| `misleading_convergence_check(terms)` | rescaling ξ → series that converges at ξ=1 but diverges for all ξ>1 is suspect | Kozik, Ferrero, Georges PRL 114, 156402 (2015) |

For alternating-sign coefficients with a singularity off the real axis (e.g. tanh), use `step = 2`. This gives
radius 1.52 vs the exact π/2.

## 3. Bold (skeleton) series — the misleading-convergence trap

`feynlab.diagrams.bold.skeleton_series_at` evaluates skeleton diagrams with the **exact** G of the Hubbard atom.
The skeleton series then converges smoothly, but beyond βU ≈ 5 it converges to the **wrong** Σ. The Luttinger–Ward
functional is multivalued (Kozik et al. 2015). Always compare with an independent route: CDet, or the bare series
plus resummation.

## 4. Fröhlich polaron: energy, Z and effective mass (`feynlab.qmc.polaron`)

The DiagMC samples G(p, τ) = Σ diagrams. At long τ:

  G(p, τ) → Z(p) exp[−(E(p) − μ) τ]

* **E(p)** comes from the slope.
* **Z** comes from the intercept. The absolute normalisation uses the order-0 diagram, whose weight
  exp[−(p²/2 − μ)τ] is known exactly.
* **m\*** comes from E(p) − E(0) = p²/2m\*. The phonon-line action carries −p·q τ for external momentum p.

```python
from feynlab.qmc.polaron import polaron_diagmc, polaron_effective_mass, feynman_energy
r = polaron_diagmc(0.2, mu=feynman_energy(0.2) - 0.15, tau_max=40)   # E = -0.2015, Z = 0.883 +- 0.016
m = polaron_effective_mass(0.2, p=0.5, mu=feynman_energy(0.2) - 0.15, tau_max=40)   # m* = 1.05 +- 0.02
```
References:

* Weak coupling: Z = 1 − α/2 and m*/m = 1 + α/6 + 0.0236α² (Rösler 1968; Smondyrev 1986).
* DiagMC: Mishchenko et al. PRB 62, 6317 (2000).

**Limitation.** The finite p used for m* includes O(p⁴) corrections. Use p ≤ 0.5 and check stability.

Quiz: **T14**.
