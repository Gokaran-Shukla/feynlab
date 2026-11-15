# 5-orbital full-Kanamori CT-HYB (`qmc.cthyb_sector.SectorCTHYB`)

**Model.**
- 5 orbitals with full Kanamori interactions (spin flip and pair hopping): U = 4, J = 0.6, U′ = U − 2J; β = 10.
- μ = 16.3, so the atomic ground state has N = 6 and is a Hund's S = 2 state.
- One bath level per orbital, Δ(iω) = V²/(iω − e_b). Lanczos ED of the same impurity (10 sites, 2²⁰ states, `sz_max=10`) is therefore exact.

**CT-HYB run.** 4 chains × 300k measured steps; the error bars are from the chain-to-chain spread.

Reproduce: `python benchmarks/cthyb_5orb/run_cthyb_5orb.py` (about 1.5 h on 2 cores).

## Accuracy against exact ED

| case | quantity | ED | CT-HYB | CT-HYB, e_cut = 10 (663 of 1024 states) |
|---|---|---|---|---|
| e_b = −2, V = 0.5 (⟨k⟩ = 0.30 per flavour) | n per flavour | 0.6090 | 0.6085 | — |
| | ⟨n↑n↓⟩ | 0.2179 | 0.2169 ± 0.0012 | — |
| | ⟨S_z²⟩ | 1.933 | 1.936 ± 0.008 | — |
| e_b = −1, V = 1.0 (⟨k⟩ = 1.97) | n per flavour | 0.6036 | 0.6019 | 0.6035 |
| | ⟨n↑n↓⟩ | 0.2075 | 0.204 ± 0.019 | 0.2075 ± 0.009 |
| | ⟨S_z²⟩ | 2.013 | 2.04 ± 0.05 | 1.99 ± 0.18 |
| | time (4 chains, serial) | — | 1291 s | 971 s |

**Exact trace.** On random operator strings, the sector trace equals the dense 64-state trace to 1e-10, and the lazy bound is never violated (`test_sector_trace_equals_dense_trace`).

**Green function (honest note).** In the e_b = −1 case, G(iω_n) agrees with ED within 1–2 chain error bars. In the e_b = −2 case, ⟨k⟩ is only 0.3 and |Δ(τ)| varies by e^{20} across [0, β]. There the standard estimator −⟨M e^{iω(τ_e − τ_s)}⟩/β has heavy tails, because M = F⁻¹ is huge where F is tiny. G is then biased by noise at the 0.1–0.3 level, while the local observables above remain accurate. This is a property of the estimator, which the dense solver shares, not of the trace. Improved estimators (worm sampling, or symmetric improved estimators) are the cure.

## Speed

| | dense `MatrixCTHYB` | `SectorCTHYB` |
|---|---|---|
| 3 orbitals, 4 × 305k steps | 349 s | 23 s (15× faster); n₀ = 0.487 ± 0.003 vs 0.494 ± 0.003, ⟨S_z²⟩ = 0.763 ± 0.016 vs 0.756 ± 0.007 |
| 5 orbitals | impractical: 1024² dense products per operator | about 1 ms per Monte Carlo step at ⟨k⟩ ≈ 2 per flavour (40 operators) |

Where the speed comes from (all rigorous, with relative error ≤ n_states × tol and tol = 1e-12):
1. **(N↑, N↓) blocks.** Each operator maps one block to exactly one block (largest block 100 states). Strings that annihilate a block are dropped before any arithmetic.
2. **Per-block norm bounds.** Blocks whose bound is negligible are skipped.
3. **Outer-state screening.** Only outer eigenstates with e^{−ℓ(E_n − E_min)} above tolerance are propagated, where ℓ is the wrap-around interval. This replaces matrix–matrix products with matrix–vector products.
4. **Lazy Metropolis.** A move is rejected with the trace bound alone when even the bound cannot pass (Sémon et al., PRB 90, 075149 (2014)).
5. **One propagation for all observables.** All local observables are measured from a single propagation.
6. **Optional Haule energy truncation** (`e_cut`). Checked above: e_cut = 2.5U gives ED values within errors, 25% faster.
