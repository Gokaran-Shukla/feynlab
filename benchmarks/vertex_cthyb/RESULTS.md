# Two-particle vertex for a continuous bath (segment CT-HYB, removal estimator)

`CTHYBSolver.solve(..., g2=(n_nu, n_om))` measures

    G2_{σσ'}(ν,ν',ω) = (1/β) ⟨ M̂_σ(ν+ω,ν) M̂_σ'(ν',ν'+ω) − δ_σσ' M̂_σ(ν+ω,ν'+ω) M̂_σ(ν',ν) ⟩
    M̂(a,b) = Σ_ij M_ji e^{iν_a τ^e_i} e^{−iν_b τ^s_j}

(Gull et al., RMP 83, 349 (2011); Hafermann et al., PRB 85, 205106 (2012)). It returns a generalized-susceptibility
dict that `vertex.vertices` and `vertex.lambda_dga` take directly. Cost per measurement: O((2N)² k) with
N = n_nu + n_om. The bath can be any Δ(iω), with no discretization.

## 1. Exactness: discrete bath, CT-HYB vs Lehmann ED

Setup: 2 bath levels at ±1 with V = 0.6, U = 2, half filling, β = 5. Box: 8 fermionic × 7 bosonic frequencies.
Script: `ed_vs_cthyb_g2.py`.

| | ED (64 states) | CT-HYB (4 × 400k steps) |
|---|---|---|
| time | 292 s | **16 s** |
| χ_↑↑: median \|Δ\|/σ, 95th percentile | — | 0.67, 2.55 |
| χ_↑↓: median \|Δ\|/σ, 95th percentile | — | 0.98, 2.56 |
| χ_m(ω = 0) (box sum) | 1.0904 | 1.0850 |
| ⟨n↑n↓⟩ | 0.10547 | 0.10591 |

The deviations follow the statistical errors. A third bath level (256 states) already exhausts memory in the
Lehmann ED, while CT-HYB does not care about bath size. The U = 0 Wick test is a unit test
(`test_cthyb_g2_wick_at_U0`).

## 2. Continuous bath: square-lattice DMFT vertex and λ-DΓA

Setup: 2D Hubbard model, t = 0.25, U = 4t, T = 0.4t, half filling. DMFT with CT-HYB (no bath fit), then G2 with
2 × 2M steps. Scripts: `dmft_cthyb_vertex.py` and `lambda_dga_from_cthyb.py`.

**Schwinger–Dyson check.** Σ from the measured vertex (equation of motion) equals Σ from the Dyson equation
within 1–2% at box 16: at ν₀, −0.1274i vs −0.1297i. At box 32 it is within 2–4%, because the high-frequency
vertex is noisier.

| vertex source | box | sum target | χ_m(Q,0) plain | λ_m | χ_m(Q,0) λ | Im Σ(π,0) / (π/2,π/2) | with exact target: χ_m(Q,0), Im Σ(π,0) |
|---|---|---|---|---|---|---|---|
| ED, 2 bath sites | 16 × 15 | 0.2977 | 6.41 | 0.044 | 5.01 | −0.146 / −0.128 | 10.8, −0.271 |
| CT-HYB, continuous bath | 16 × 15 | 0.2946 | 8.00 | 0.072 | 5.07 | −0.133 / −0.114 | 11.3, −0.265 |
| CT-HYB, continuous bath | **32 × 31** | 0.3256 | 8.58 | 0.060 | 5.66 | −0.140 / −0.121 | 8.0, −0.200 |
| exact sum (n − 2d)/2 | ∞ | 0.3543 | | | | | |
| TPSC (reference) | | | | | 6.51 | | |

What the table shows:
* At the same box, the CT-HYB and ED vertices give the same λ-DΓA within about 10% for χ_m(Q) and Im Σ. The
  plain DMFT χ_m(Q) differs (6.4 vs 8.0); this is the effect of discretizing the bath with 2 sites.
* The larger box, affordable only with CT-HYB, moves the box target towards the exact sum rule
  (0.295 → 0.326 → 0.354). The spread between the two λ variants shrinks from a factor 2.2 to 1.4 in χ_m(Q),
  and from 2.0 to 1.4 in Im Σ. The box ambiguity noted in `benchmarks/lambda_dga` is therefore a finite-box
  effect, and it closes as the box grows.
