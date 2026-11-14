# Tutorial 14 — Advanced transport: T-matrix and matrix disorder, lattice crossed diagrams, orbital moments, Berry-curvature dipole, Edelstein/SOT, finite T and 3D Wannier

Roadmap items C12–C16. Everything uses the same band-model interface: `H(kx, ky)`, `dH(kx, ky)` and `operators`.

## 1. Beyond scalar Gaussian disorder (`DisorderedHall`)

```python
DisorderedHall(model, mu, n_i, V0, imp_matrix=M, t_matrix=False)
```

* **`imp_matrix=M`** sets the impurity potential to V₀M. Use σ_z for magnetic impurities, an orbital projector
  for orbital-selective disorder, or L_z for orbital-dependent scattering. The self-energy is Σ = n_i V₀² M ḡ M†
  and the ladder rung is n_i V₀² M ⊗ M*. M = 1 reproduces scalar disorder exactly.
* **`t_matrix=True`** gives the self-consistent T-matrix approximation (SCTMA):
  T = V(1 − ḡV)⁻¹ with Σ = n_i(T − V) and the ladder rung n_i T^R ⊗ T^A.
  * Expanding to third order gives n_i V³ ḡ: **conventional skew scattering emerges automatically**, with no
    need to put in a V₁ by hand.
  * Check: at fixed n_iV₀², σ_I(T-matrix) = σ_I(Born) + σ_skew(V₁ = V₀) to 0.1 %. In the Sinitsyn model the skew
    term scales as V₀ (relative to n_iV₀²) and dominates in clean samples.
  * For strong scatterers (V₀ = ±0.5) the part even in V₀ stays within 1 % of Born, while the odd part is the
    resonant skew scattering.

## 2. Crossed X and Ψ diagrams on a lattice

```python
crossed_hall(model, mu, lattice=True, runs=((0.03, 0, 1024), (0.015, 0, 2048)))
```
* The momentum grid is the full Brillouin zone, so r runs over lattice vectors.
* No Nyquist masking is needed, because H(k) is periodic.

Benchmark:

* QWZ with m = −2 + 0.3 is a massive Dirac cone at Γ with a quadratic mass correction, Δ(k) = Δ − k²/2.
  Its crossed diagrams agree with the continuum k·p model `DiracQuadraticMass(v=1, Δ=0.3, b=−0.5)`:
  X = 0.795 vs 0.790, Ψ = −0.50 vs −0.52.
* Ψ vanishes only for the linear Dirac model of Ado et al. (EPL 111, 37004 (2015)). The k² mass of real
  materials makes it finite. This is easy to miss with continuum-only codes.

## 3. Modern orbital moment and orbital magnetisation (`transport.response`)

  m_n(k) = Im Σ_{m≠n} v^x_nm v^y_mn / (E_n − E_m)   (|e|/ħ = 1, electron charge −|e|)

For the gapped Dirac cone this gives m = Δv²/[2(Δ² + v²k²)], identical in both bands. That is the valley moment
of Xiao, Yao, Niu (PRL 99, 236809 (2007)), reproduced to machine precision.

`orbital_magnetization` adds the Berry-curvature term,

  M = Σ_n f_n [m_n − (μ − E_n) Ω_n]

(Ceresoli et al. PRB 74, 024408 (2006)).

The **modern** moment includes inter-atomic (itinerant) contributions that the atom-centred approximation (ACA)
misses. Using it as the operator O = m in `intrinsic_hall(model, mu, op=...)` gives the modern-theory orbital
Hall conductivity.

## 4. Berry-curvature dipole and nonlinear Hall effect

  D_a = Σ_n ∫ f_n ∂_a Ω_n  (Fermi sea)  = Σ_n ∫ (−∂f/∂E) v^a_nn Ω_n  (Fermi surface),
  j_a = −(e³τ/2ħ²) ε_{adc} D_{bd} E_b E_c   (Sodemann & Fu PRL 115, 216806 (2015)).

`berry_curvature_dipole(model, mu, grid, form='surface'|'sea')`:

* Both forms agree to 1e-8 on the tilted Dirac cone (`TiltedDirac`).
* Without tilt D = 0 by rotational symmetry.
* For bilayer WSe₂, CIPS etc., use a `WannierTB` model and all components (D_xx, D_xy, …).

## 5. Edelstein effect and spin–orbit torque

* `edelstein_boltzmann(model, mu, op, grid)` gives δ⟨O⟩ per eEτ.
  * Rashba: |δS_y| = mα/4π, independent of μ when both bands are occupied (Edelstein 1990). Reproduced to 0.01 %.
* `edelstein_kubo(..., Gamma)` gives the constant-broadening Kubo formula, split into two parts:
  * **odd** (intraband, ∝ 1/Γ = 2τ), which matches Boltzmann with τ = 1/2Γ;
  * **even** (interband, intrinsic, Berry-curvature-like), which vanishes for non-magnetic Rashba by
    time reversal. With exchange (Rashba2DEG(Δ_ex ≠ 0)) it gives the intrinsic **damping-like** torque.
* The torque on the magnetisation is T = J_ex m × δS.
* Use op = L/2 for the **orbital** Edelstein effect.
* Reference: Freimuth, Blügel, Mokrousov PRB 90, 174423 (2014).

## 6. Finite temperature (Kubo–Bastin)

`kubo_bastin_T(sigma_of_E, mu, T)` computes σ(μ,T) = ∫ dE (−f′) σ(E, 0) by Gauss–Legendre quadrature.
`intrinsic_hall(..., beta=1/T)` gives the Fermi-sea result directly. Both agree to 1e-6 for the Dirac cone at
T = 0.05 and 0.1, once the analytic lower-band tail (`lower_band_tail`) is included. This is a stringent test of
k-grid quality.

## 7. Three-dimensional Wannier models with adaptive k refinement (`transport.wannier3d`)

```python
from feynlab.transport.wannier3d import WannierTB3D
w = WannierTB3D("seedname_hr.dat", lat=lattice_vectors_in_A, operators={"Sz": Sz, "Lz": Lz})
r = w.intrinsic_hall(mu=E_F, nk=(40, 40, 40), op=None)       # sigma_x,y,z = (sigma_yz, sigma_zx, sigma_xy)
print(r["sigma_S_per_cm"])
```
* Refinement replaces points with large |Ω| **or** a band within a window of μ by an n_sub³ sub-mesh
  (Wang, Yates, Souza, Vanderbilt PRB 74, 195118 (2006)).
* Stacked Chern layers give σ_xy = C e²/(hc) to 1e-7 (Halperin 1987).
* In the metallic case, adaptive 40³ is within 0.2 % of a converged 200×200×160 mesh, while plain 40³ is off
  by 0.6 %.
* Pass `op=Sz` (spin Hall) or `op=Lz` (orbital Hall) to get the operator-resolved tensors.

Quiz: **T15**.
