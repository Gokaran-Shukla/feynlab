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
5. Single-site DMFT only (no cluster DMFT); X/Ψ engine is for continuum models (lattice FFT version is a
   straightforward extension).

"Best in the world" is not something a benchmark suite can establish. What the suite does establish is that
every number it reports agrees with an exact or published result within the stated tolerance.
