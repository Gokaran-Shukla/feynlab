# Tutorial 11 — Beyond single-site DMFT: DCA, two-particle response, NRG and the Kondo effect

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


Roadmap items A3, A6 and A7.

## 1. Dynamical cluster approximation (`feynlab.dmft.dca`)

The Brillouin zone is tiled into N_c = L_c² patches around cluster momenta K. The self-energy is constant on
each patch: Σ(k, iω) → Σ(K, iω) (Hettler et al. PRB 58, R7475 (1998); Maier et al. RMP 77, 1027 (2005)).

  Ḡ(K) = (N_c/N) Σ_{k∈K} [iω + μ − ε_k − Σ(K)]⁻¹,  𝒢₀(K)⁻¹ = Ḡ(K)⁻¹ + Σ(K)

The cluster solver here is second-order perturbation theory in K space (exact up to U²). For N_c = 1 this is
exactly DMFT with IPT.

```python
from feynlab.dmft.dca import DCA2D
r1 = DCA2D(Lc=1, nk=64, U=1.5, beta=10).solve()   # Z = 0.486644  (single-site DMFT: 0.486638)
r4 = DCA2D(Lc=2, nk=64, U=1.5, beta=10).solve()   # Z(0,0) = 0.565, Z(pi,0) = 0.475
```
**Lesson.** Momentum differentiation appears as soon as N_c > 1. The antinodal patch (π,0) is more strongly
renormalised than (0,0), because it couples to (π,π) spin fluctuations. This is the seed of the pseudogap.

## 2. Two-particle response (`feynlab.twoparticle`)

  χ₀(q, iΩ) = −(T/N) Σ_k,ν G(k, iν) G(k+q, iν+iΩ),  computed by FFT in (r, τ).
  RPA: χ_s = χ₀/(1 − Uχ₀),  χ_c = χ₀/(1 + Uχ₀),  Stoner criterion U χ₀(Q) = 1.

The **compressibility sum rule** χ₀(q→0, Ω=0) = ∂n/∂μ holds to 1e-8. The nesting peak at Q = (π,π) dominates at
half filling. `bubble(eps_k, sigma, grid, mu)` accepts a self-energy, so the same routine gives the dressed
bubble of DMFT/DCA (the first step towards DΓA / dual fermions).

## 3. Numerical renormalisation group (`feynlab.qmc.nrg`)

Wilson's NRG (RMP 47, 773 (1975); Bulla, Costi, Pruschke RMP 80, 395 (2008)) works as follows.

* The band is discretised logarithmically, ±Λ^{−n}, and mapped to a **Wilson chain** with hoppings
  t_n ∝ Λ^{−n/2}. The prefactor A_Λ = ½(1+1/Λ) ln Λ/(1−1/Λ) corrects the discretisation.
* Sites are added one at a time. The N_keep lowest states are kept, and charge and S_z are used as quantum numbers.
* Thermodynamics are evaluated at T_N = Λ^{−(N−1)/2}/β̄ with β̄ = 2. Smaller β̄ (e.g. 0.75) under-counts the
  entropy because of truncation; this effect was tested.

```python
from feynlab.qmc.nrg import impurity_thermodynamics, kondo_temperature
T, S, Tchi = impurity_thermodynamics(U=0.1, Gamma=0.01, n_iter=60, n_keep=300)
TK = kondo_temperature(T, Tchi)             # Wilson's definition T_K chi(T_K) = 0.07
```

The output shows the three fixed points of the symmetric Anderson model:

| regime | S_imp | T χ_imp |
|---|---|---|
| free orbital (T ≫ U) | ln 4 | 1/8 |
| local moment (T_K ≪ T ≪ U) | ln 2 | 1/4 |
| strong-coupling Fermi liquid (T ≪ T_K) | 0 | 0 |

Benchmarks:

* U = 0 entropy matches the exact resonant-level result within 2–4 %.
* T_K follows **Haldane's formula** T_K ∝ √(UΓ/2) exp(−πU/8Γ + πΓ/2U) with fitted slope −1.00
  (Haldane PRL 40, 416 (1978)).
* Tχ curves for different U collapse onto one universal function of T/T_K (figure `benchmarks/figures/nrg_kondo.png`).

**Why it matters for ED-DMFT.** NRG has no bath-size parity problem (Tutorial 8). Its logarithmic resolution
captures the Kondo scale at any U.

Quiz: **T12**.
