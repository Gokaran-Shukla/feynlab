# Tutorial 12 — Analytic continuation: MaxEnt, Nevanlinna and stochastic

Roadmap item A5. The task is to find A(ω) ≥ 0 from G(iω_n) = ∫ dω A(ω)/(iω_n − ω). The problem is
exponentially ill-conditioned: very different A give the same G within the noise. Each method adds a different
piece of prior information.

| method | prior / principle | strengths | weaknesses |
|---|---|---|---|
| MaxEnt (`core.maxent`) | maximise αS[A] − χ²/2 with an entropy relative to a default model | stable, error-aware | over-smooth, depends on the model |
| Nevanlinna (`continuation.nevanlinna`) | G is a Nevanlinna function, so Schur interpolation is **exact** on the data | resolves sharp peaks, exact causality | noise-free data only; precision loss for many points |
| stochastic (`continuation.stochastic`) | average over all A with χ² ≈ N (Sandvik, Mishchenko) | no default model | slower, spiky individual samples |

## 1. Nevanlinna (Fei, Yeh, Gull PRL 126, 056402 (2021))

* The Möbius map h(z) = (z − i)/(z + i) turns −G(iω_n) into Schur values θ_n with |θ_n| < 1.
* The **Pick criterion** (`pick_ok`) checks that a Nevanlinna interpolant exists. If it fails, the data are too
  noisy or inconsistent.
* The Schur algorithm builds θ(z) recursively. It is done in `mpmath` with `dps = 64` digits, because the
  recursion loses about one digit per point.
* `nevanlinna_auto` picks the largest number of points whose spectrum is non-negative and normalised.

```python
from feynlab.continuation.nevanlinna import Nevanlinna, nevanlinna_auto
nv = Nevanlinna(iw, G)             # nv(iw) == G to 1e-12; nv(2j) matches the exact value to 1e-6
A = nevanlinna_auto(iw, G, omega, eta=0.05)
```
**Honest limitation.** The Hardy-basis optimisation of the free function (Fei et al., Sec. III) is not
implemented. With 60 or more points the precision loss forces n ≤ 40.

## 2. Stochastic continuation

`stochastic_continuation(iw, G, sigma, omega, n_delta=4000)` represents A as a sum of δ-functions. They are
sampled with a Metropolis walk at fictitious temperature Θ, and the samples are averaged at the Θ where
χ²/N ≈ 1. This is Beach's/Sandvik's criterion, and the sweep is accelerated with numba.

Results:

* Normalisation is exact.
* χ²/N ≈ 1 is reached.
* The L1 error of 0.8 is larger than MaxEnt's 0.21 on the two-Gaussian test.

Use it as an **independent cross-check** of MaxEnt peak positions. Do not use it as the primary method.

## 3. Practical protocol

1. Check the data: high-frequency tail 1/iω, the Pick criterion, and the covariance (Tutorial 9).
2. Run MaxEnt with two default models (flat and Gaussian). Features present in both are robust.
3. Cross-check peak positions with stochastic continuation. For converged, near-noise-free data
   (e.g. from ED), also use Nevanlinna.
4. Never interpret features narrower than ~πT (2πT is the Matsubara spacing).

Figure: `benchmarks/figures/continuation.png`. Quiz: **T13**.
