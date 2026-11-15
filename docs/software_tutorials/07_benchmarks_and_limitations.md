# Tutorial 7 — Benchmarks, accuracy and limitations

Run `python benchmarks/run_benchmarks.py` (≈30 min; `--quick` ≈15 min); bath-size study: `benchmarks/ed_dmft/` →
`benchmarks/RESULTS.md`, `results.json`, figures.

## Published results reproduced

| area | check | reference |
|---|---|---|
| diagram enumeration | 2, 10, 74, 706 G diagrams; 2, 6, 42 proper; 0, 4, 32 improper | Mohr, physics/0512022 |
| | loop-resolved Σ counts 2+ℓ, 10+9ℓ+ℓ², 74+91ℓ+23ℓ²+ℓ³ | Molinari, **PRB** 71, 113102 (2005) |
| diagrammatic MC | atom: orders 3, 4 cancel; Anderson series → ED | exact |
| **polaron DiagMC** | E₀(α) for α = 0.5, 1, 2, 3 (all-order bare series) | Hahn et al., **PRB** 97, 134305 (2018), Table I; weak-coupling series |
| | Feynman variational E₀ for α = 1, 3, 5 (to 1e-6) | Feynman, Phys. Rev. 97, 660 (1955) |
| **bold / skeleton series** | atom: physical branch at βU = 0.5, unphysical at βU = 6 | Kozik, Ferrero, Georges, **PRL** 114, 156402 (2015) |
| CT-INT | Anderson impurity vs ED, sign = 1 | Gull et al., **RMP** 83, 349 (2011) (method) |
| weak coupling | Yamada–Yosida 1/Z − 1; lattice tail U²n(1−n); RPA → 2nd order | Prog. Theor. Phys. 53, 970 (1975) |
| **DMFT, Mott** | Lanczos ED-DMFT (T = 0, 6/7 bath) U_c2, U_c1; bath-size study 5–9; IPT hysteresis (informational) | Bulla, **PRL** 83, 136 (1999); Karski et al. **PRB** 72, 113110 (2005); Georges et al. **RMP** 68, 13 (1996) |
| | LatticeDMFT(semicircle) ≡ BetheDMFT | identity |
| MaxEnt | χ²/N ≈ 1, L1 error on a two-peak spectrum | Jarrell & Gubernatis, Phys. Rep. 269, 133 (1996) |
| electron gas | plasmon 3/5; RPA e_c slope and constant | Mattuck Ex. 13.6; Gell-Mann & Brueckner, Phys. Rev. 106, 364 (1957) |
| BCS | Δ(0)/T_c = 1.764, T_c prefactor 1.134 | BCS |
| Berry | QWZ, Haldane Chern numbers; Rashba SHC e/8π | **PRB** 74, 085308; **PRL** 61, 2015; **PRL** 92, 126603 |
| AHE (NCA) | intrinsic, side jump, adist, Gaussian skew, V₁³ skew, total, σ_xx | Sinitsyn et al., **PRB** 75, 045315 (2007) |
| **AHE (crossed)** | X diagram, Ψ = 0, full σ_xy = −8εm³/(ε²+3m²)² | Ado et al., EPL 111, 37004 (2015) |
| SHE | vertex cancellation in the Rashba 2DEG | Inoue, Bauer, Molenkamp, **PRB** 70, 041303 (2004) |

## What changed in v0.2 (the v0.1 limitations)

| v0.1 limitation | status in v0.2 |
|---|---|
| crossed X/Ψ diagrams missing | **done** — `transport.crossed` (real-space FFT form, any continuum model); X within 1% of Ado et al., Ψ → 0 |
| DMFT only Bethe; Padé only | **done** — `LatticeDMFT` (square, cubic, hypercubic, any ε_k or DOS) and `core.maxent` |
| no bold-line / polaron DiagMC | **done** — `qmc.polaron` (Prokof'ev–Svistunov), `diagrams.bold` (skeleton series, self-consistent bold scheme) |
| IPT overestimates U_c2 | **resolved (v0.3)** — T = 0 Lanczos ED-DMFT with 5–9 bath sites: U_c2 = 2.955 ± 0.002 D (odd baths) vs consensus 2.92–3.0 D; U_c1 = 2.375 ± 0.025 D (even baths) vs 2.39 ± 0.02 D. See Tutorial 8 for the bath-parity effect |
| stand-alone Gaussian skew rate | **partially** — see below |

## Remaining limitations (honest list)

1. **Stand-alone Gaussian skew (semiclassical w⁽⁴ᵃ⁾).** The interband self-energy pairings give exactly 1/3 of
   Sinitsyn's Eq. 71 and the |T₂|² interference exactly 1/6; the arc pairing, the crossed |T₂|² pairing and
   on-shell×on-shell terms do not have the sin(φ−φ′) form and do not supply the rest. The Gaussian skew term
   *is* computed (to 1e-5) by the Kubo route; only its semiclassical rate-by-rate decomposition is incomplete.
2. **ED-DMFT bath parity**: odd baths must be used for the metal (U_c2) and even baths for the insulator
   (U_c1); the wrong parity gives artifacts (Tutorial 8). U_c2 from a finite bath depends on the estimator
   (linear Z → 0 extrapolation vs. collapse point); the linear window 0.04 ≤ Z ≤ 0.15 is fixed in advance.
3. **Polaron DiagMC precision** is ~0.2% (8 chains × 1.5·10⁸ updates); the reference (Hahn et al.) is ~0.05%.
   Improve with longer runs or the exact energy estimator of Mishchenko et al.
4. **Bold scheme** evaluates skeleton diagrams up to order 4 by direct MC (no sampling over orders); the
   Kozik et al. misleading convergence is reproduced qualitatively at these orders.
5. ~~Single-site DMFT only; X/Ψ for continuum models only~~ → **resolved in v1.0** (DCA 2×2, lattice X/Ψ via BZ FFT).

## What changed in v1.0 (all 22 ROADMAP items)

| roadmap item | v1.0 status | benchmark (`benchmarks/RESULTS_v1.md`, 54/54 pass) |
|---|---|---|
| A1 CT-HYB | **done** (segment picture, density–density) | n, ⟨n↑n↓⟩, Σ(iω₀) vs ED |
| A2 Kanamori multi-orbital | **done** (atom, Hubbard-I, sparse ED ≤ 8 sites) | U−3J, U+J, 3U−9J |
| A3 DCA | **done** (2D, N_c = 1, 4; second-order cluster solver) | N_c = 1 ≡ DMFT |
| A4 DFT+DMFT | **done** (hr.dat, FLL/AMF/occupation dc; Hubbard-I, ED) | U = 0 limit, Mott |
| A5 continuation | **done** (Nevanlinna, stochastic) | interpolation, sum rules |
| A6 two-particle | **partial** (bubble, RPA; no vertex ladders / DΓA) | compressibility sum rule |
| A7 NRG | **done** (thermodynamics, Kondo) | resonant level, ln 2, Haldane slope |
| B8 CDet | **done** | atom exact, Anderson vs ED |
| B9 bold / misleading convergence | **done** (diagnostic) | Kozik et al. |
| B10 resummation | **done** | Domb–Sykes, Padé, conformal |
| B11 polaron Z, m* | **done** (Fröhlich; Holstein/SSH not yet) | weak coupling |
| C12 3D Wannier + modern OAM | **done** | Chern stacks; Xiao–Yao–Niu |
| C13 matrix disorder + T-matrix | **done** (finite-range disorder not yet) | Born + skew |
| C14 lattice crossed diagrams | **done** | k·p limit |
| C15 BCD, Edelstein, SOT | **done** | Sodemann–Fu; Edelstein |
| C16 finite-T Kubo–Bastin | **done** | Sommerfeld/Bastin consistency |
| D17 GPU backend | **done** (CuPy optional; MPI not yet) | — |
| D18 HDF5 checkpoints | **done** | restart test |
| D19 Γ-method, jackknife, Anderson/Broyden | **done** | AR(1) τ_int; same fixed point |
| D20 CI, CITATION, Zenodo, MkDocs | **done** | weekly benchmarks |
| D21 CLI + TOML | **done** | `feynlab dmft|hall|polaron|nrg|bench` |
| D22 units | **done** | (e²/h)/Å → S/cm |

## What changed in v1.1 (every v1.0 limitation addressed)

| v1.0 limitation | v1.1 status | benchmark (`benchmarks/RESULTS_v11.md`) |
|---|---|---|
| CT-HYB density-density only | **resolved**: `qmc.cthyb_matrix.MatrixCTHYB` (full local trace, any H_loc, spin flip + pair hopping) | n, docc, ⟨S_z²⟩, G vs 2-orbital ED within 1σ; sign = 1 |
| DFT+DMFT without charge self-consistency | **resolved at the Wannier level**: `solve_csc` (H(k) = H₀ + V[n], dc updated), `density_matrix_k`, `write_density_matrix` for external DFT codes; CT-HYB solver added; exact-tail electron count | CSC(U=0) = Hartree mean field to 6e-6; N(k) exact to 3e-7 |
| DCA cluster solver second order | **resolved**: cluster CT-INT (`ClusterCTINT`, `DCA2D.solve(solver='ctint')`) | dimer+bath vs ED within 1σ; N_c = 1 = single-site CT-INT |
| Nevanlinna without Hardy optimisation, n ≤ 40 | **resolved**: `hardy_optimize`, `hardy_auto`; multiprecision Pick test, dps = 3n | L1 0.04–0.06 (raw ~1) |
| stochastic continuation worse than MaxEnt | **resolved**: SAC with continuous positions, adaptive moves, Shao–Sandvik Θ criterion | L1 0.075 vs MaxEnt 0.25 |
| polaron m* few-% spread | **resolved**: momentum-fluctuation estimator | α = 0.2: 1.0340 ± 0.0002 (ref 1.0343) |
| semiclassical Gaussian skew partial | **resolved** (the v1.0 reference had a factor-2 transcription error) | Sinitsyn to 1e-6; Kubo residual to 1e-4 |
| two-particle = bubble/RPA only | **resolved**: ED vertex (Lehmann, divided differences), F/Γ in d/m channels, Schwinger–Dyson, DMFT χ(q), ladder DΓA, TPSC | Wick 1e-15; EOM vs ED 0.05%; χ_m(q=0) vs dm/dh 0.04%; DΓA → 2nd order as U → 0 |
| finite-range disorder | **resolved** (semiclassical: Born, side jump, skew, Gaussian skew) | σ_xx vs Bessel formula 1e-10 |
| orbital crossed diagrams for 3D Wannier | **resolved**: `LatticeCrossed` (any Bravais lattice, 2D/3D, operator currents) | oblique cell and 3D stacking reproduce the square-lattice result |
| MPI (chains only) | **resolved**: `utils.parallel.map_chains` (MPI / process / serial) in all QMC solvers | identical results across backends |
| Holstein / SSH polarons | **resolved**: variational ED (`qmc.lattice_polaron`) | Holstein E₀ = −2.46968472393287 (published 21-digit value) |
| adaptive 3D refinement not cost-optimal | **improved**: velocity-based Fermi-surface criterion; hot spots opt-in | 0.1% at 263k points vs 0.2% for a uniform 216k mesh |

## v1.2: first-principles charge self-consistency (Quantum ESPRESSO)

| v1.1 limitation | v1.2 status | check |
|---|---|---|
| CSC only at the Wannier level (model functional) | **resolved**: `feynlab.dft.qe_csc.QEDMFT` drives pw.x (QE ≥ 7.0 DMFT hook) + Wannier90; ΔN(k) archive written by feynlab, density updated by QE | U = 0: ΔN ≤ 1.5e-7, E_QE unchanged to 1e-8 Ry; SrVO₃ U = 5, J = 0.75 converged (`benchmarks/qe_csc/RESULTS.md`) |
| ED-DMFT for t₂g limited to 1 bath site/orbital | **resolved**: `qmc.multiorbital_lanczos` (sectors, continued fraction, finite-T Boltzmann average), 2 bath sites per orbital = 9 sites | full ED 7e-12 (T→0), 3e-6 (β = 20) |
| **bug fix**: `DFTDMFT` ED solver shifted bath levels by −μ twice | fixed; results of the `ed` solver before v1.2 are wrong whenever μ ≠ 0 | ED vs Lanczos vs CT-HYB at n = 0.8: Z = 0.5024 / 0.5025 / 0.497 |

## v1.3: closing the v1.2 list one by one

| v1.2 limitation | status | check |
|---|---|---|
| ladder DΓA without λ correction | **resolved**: `vertex.lambda_dga` (Moriya λ_m, γ-vertex form, bubble tail, box or exact target) | λ = 0 equals the plain ladder to 4e-15; sum rule exact; below T_N^DMFT χ_m(Q) = 29 > 0 where the plain ladder gives −11 (`benchmarks/lambda_dga`) |
| matrix CT-HYB practical only for 2–3 orbitals | **resolved**: `qmc.cthyb_sector.SectorCTHYB` — (N↑,N↓) blocks, bound screening, outer-state matrix–vector propagation, lazy trace, optional Haule truncation | trace = dense to 1e-10; 5 orbitals vs exact ED (2²⁰ states): ⟨S_z²⟩ 1.936 ± 0.008 vs 1.933; 15× faster at 3 orbitals (`benchmarks/cthyb_5orb`) |
| exact vertex only for small ED baths | **resolved**: two-particle G2 from segment CT-HYB (removal estimator, `solve(..., g2=(n_nu, n_om))`), any continuous Δ | vs ED: deviations at the statistical level (16 s vs 292 s); Schwinger–Dyson 1–2%; λ-DΓA box 32 (`benchmarks/vertex_cthyb`) |
| finite-range disorder only semiclassical | **resolved**: `transport.finite_range.FiniteRangeHall` — k-dependent SCBA and ladder vertex by FFT convolution + GMRES | d = 0 equals DisorderedHall to 1e-9; σ_xx vs Bessel Boltzmann ≤ 0.4%; extrinsic AHE vs semiclassical ≤ 0.6% for k_F d ≤ 2.6 (`benchmarks/finite_range_kubo`) |
| cluster CT-INT sign problem | **managed** (not removable, NP-hard): pooled sign-reweighted estimator with jackknife errors, pilot-run δ optimisation, exponential cost model | 2×2 plaquette U = 6t: ⟨s⟩(β) measured, Δf = 0.0076; optimal δ = 0.51 (δ = 0.3 → ⟨s⟩ = 0.10); jackknife ≈ spread of independent runs (`benchmarks/ctint_sign`) |
| ED vertex slow (matrix exponential per frequency) | partly: closed-form divided differences with a near-degenerate fallback, 3.5–4.5× faster | 5e-11 relative vs the matrix exponential |

## Remaining limitations after v1.2 (honest list)

* **QE charge self-consistency** is restricted to nspin = 1 (a QE limitation of the DMFT hook). The DMFT total
  energy (Galitskii–Migdal + double-counting correction) is not assembled; QE's printed energy is the DFT
  functional at the interacting density.
* **Vertex**: the ED (Lehmann) vertex is practical up to 2 bath sites. The CT-HYB vertex (v1.3) handles any bath,
  but it is single-orbital and paramagnetic, and uses the removal estimator (no worm sampling or improved
  estimators). The multi-orbital `SectorCTHYB` has no G2 measurement yet.
* **λ-DΓA** (v1.3, `vertex.lambda_dga`): only the Moriya λ_m is applied, with no λ_dm and no self-consistency.
  With the small ED box the result depends on the sum-rule target by up to a factor 2 in Im Σ
  (`benchmarks/lambda_dga/RESULTS.md`).
* **Cluster CT-INT** away from half filling keeps its sign problem: it is NP-hard in general. v1.3 measures it
  (jackknife), optimizes δ and predicts the cost, but large clusters at low T stay exponentially expensive. There
  are no submatrix updates.
* **Sector CT-HYB** (v1.3) handles 5 orbitals at about 1 ms per step. The trace is still recomputed per move
  (no tree or skip-list caching), H_loc must conserve (N↑, N↓), and Δ must be diagonal in flavour. The
  M-matrix G estimator is noisy at low expansion order with a strongly asymmetric Δ(τ); worm sampling is not
  implemented.
* **Finite-range Kubo–Středa** (v1.3) covers Gaussian (NCA) disorder only. The V1³ skew vertex and the crossed
  X/Ψ diagrams remain δ-correlated, and the model must be 2D.
