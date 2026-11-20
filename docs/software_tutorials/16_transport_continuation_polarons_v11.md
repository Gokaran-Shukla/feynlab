# Tutorial 16 — Gaussian skew from rates, finite-range disorder, crossed diagrams for Wannier lattices, better continuation, polaron masses

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


## 1. The semiclassical Gaussian ("intrinsic") skew scattering is complete

`SemiclassicalHall(...).compute()['skew_gaussian']` builds the antisymmetric fourth-order rate w⁽⁴ᵃ⁾ for any band
model from two pieces:

1. T₁T₃* + c.c. with adjacent impurity pairings: the on-shell self-energy −iπ n V₀² Σ_FS |u⟩⟨u| dressed by
   the off-shell interband propagator;
2. the |T₂|² interference between one on-shell intraband state and one off-shell interband state.

Results:

* Massive Dirac model: Sinitsyn's Gaussian skew is reproduced to 1e-6 at μ = 0.7, 1, 2.
* Dirac cone with a quadratic mass (no analytic answer): it matches the Kubo–Středa NCA residual to 1e-4.

**What went wrong in v1.0, and the lesson.** The very same terms were computed in v1.0 and judged to give "½ of
Sinitsyn's Eq. 71". In fact the reference formula had been transcribed with an extra factor of 2. Feeding the
transcribed rate into the Boltzmann solver gave exactly twice the known conductivity, which exposed the error.

> Always validate the reference expression through an independent observable before concluding that the
> physics is missing.

## 2. Finite-range (smooth) impurities

`SemiclassicalHall(model, mu, band, range_d=d)` gives Gaussian-shaped impurities with form factor
u(q) = exp(−q²d²/4).

* Born rate ∝ |u(q)|²; third-order rate ∝ u(q₁)u(q₂)u(q₃).
* The k-dependent on-shell self-energy enters w⁽⁴ᵃ⁾.
* For real, even u(q) the coordinate shift is unchanged.

Check: σ_xx agrees with the closed-form Bessel-function result `dirac_finite_range_sigma_xx` to 1e-10.

Physics: side jump and anomalous distribution decrease with d. Smooth disorder suppresses large-momentum
transfer, which carries the interband coherence (compare Ado et al. PRB 96, 235148 (2017)).

| d | σ_xx | side jump | Gaussian skew |
|---|---|---|---|
| 0 | 85.71 | −0.2143 | −0.1378 |
| 0.5 | 108.48 | −0.2056 | −0.1348 |
| 1.0 | 205.18 | −0.1802 | −0.1197 |

## 3. Crossed X/Ψ diagrams for any Bravais lattice in 2D and 3D

`LatticeCrossed(model, mu, nV2, N=(N1, N2[, N3])).compute(a, b, op)` works for `WannierTB` and `WannierTB3D`
models:

* momenta k = Σ(n_a/N_a) b_a;
* r runs over lattice vectors;
* ∫d^dk/(2π)^d → (1/(N V_cell)) Σ_k.

Checks:

* The same QWZ model written in an **oblique** unit cell reproduces the square-lattice result to 1e-6.
* Decoupled 3D layers give the 2D value per layer.
* `op=Sz` or `op=Lz` gives the spin or orbital crossed diagrams.

## 4. Analytic continuation, upgraded

* **SAC** (Shao–Sandvik): continuous-position δ functions with adaptive single and pair moves, annealing, and
  the χ²_min + ½√(2χ²_min) criterion. L1 error 0.075 on the two-Gaussian test, against 0.25 for MaxEnt.
* **Nevanlinna + Hardy optimisation**:
  * The free Schur function is expanded in Hardy functions and optimised for smoothness and normalisation,
    with an analytic gradient.
  * The basis is grown H = 1 → H_max with warm starts.
  * The smoothest normalised result over a small (λ, H) grid is selected.
  * L1 error 0.04–0.06 for 20–40 points, against about 1 without optimisation.
  * Precision grows with the number of points (dps = 3n), and the Pick test runs in multiprecision.

## 5. Polaron effective mass without finite-difference bias

  G(p,τ)/G(0,τ) = e^{−p²τ/2}⟨e^{p·X}⟩,  X = Σ_arcs q (t₂ − t₁)  ⇒  m/m* = 1 − lim ⟨X²⟩/(3τ)

The fluctuation of the phonon momentum carried along the line gives m* from a single p = 0 run:

| α | m*/m (DiagMC) | weak coupling |
|---|---|---|
| 0.2 | 1.0340 ± 0.0002 | 1.0343 |
| 0.5 | 1.0889 ± 0.0008 | 1.0892 |

## 6. Lattice polarons: Holstein and SSH (`qmc.lattice_polaron`)

Variational exact diagonalization (Bonča–Trugman–Batistić) in the translation-invariant phonon basis:

* 1D Holstein at t = ω = g = 1: E₀ = −2.46968472393287, the published 21-digit value. Converged with
  12 generations.
* Weak coupling: E₀ matches second-order perturbation theory for both Holstein and SSH coupling.
* Strong coupling: E₀ approaches the Lang–Firsov limit.
* Z(k) and m*/m are also returned.

## 7. Cheaper adaptive 3D refinement

`WannierTB3D.intrinsic_hall` now refines only cells actually crossed by the Fermi surface:
|E_n − μ| < |v_n| × (half cell diagonal) + 4k_BT. Hot-spot refinement is opt-in, because in gapped regions the
uniform midpoint rule is already spectrally accurate.

Quiz: **T17**.
