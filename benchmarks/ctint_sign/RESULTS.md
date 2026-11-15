# Cluster CT-INT sign problem: measured, optimised, priced

**Setup.**
- 2×2 plaquette with t = 1, U = 6 and μ = U/2 − 1 (hole doped).
- Weiss field from DCA coarse-graining of the non-interacting lattice.
- Each run: 2 chains × 200k steps.

Reproduce: `python benchmarks/ctint_sign/run_sign.py` (about 15 min).

The sign problem cannot be removed in general: it is NP-hard (Troyer & Wiese, PRL 94, 170201 (2005)). v1.3 adds three tools instead:
- correct statistics under the sign;
- a pilot-run choice of the auxiliary shift δ;
- a cost model that tells in advance whether a run is affordable.

## 1. Sign vs temperature, and the cost model

| β | ⟨s⟩ | ⟨k⟩ | Σ(π,0; iω₀) | jackknife error |
|---|---|---|---|---|
| 2 | 0.9998 ± 0.0002 | 20 | 1.97 − 1.33i | 0.03 / 0.04 |
| 4 | 0.967 ± 0.005 | 38 | 1.78 − 0.96i | 0.04 / 0.04 |
| 6 | 0.883 ± 0.018 | 57 | 1.66 − 0.74i | 0.05 / 0.05 |
| 8 | 0.759 ± 0.034 | 76 | 1.70 − 0.58i | 0.06 / 0.08 |
| 10 | 0.700 ± 0.029 | 95 | 1.81 − 0.62i | 0.07 / 0.08 |

The fit ⟨s⟩ = e^{−β N_c Δf} gives Δf = 0.0076. At fixed error, the predicted CPU cost relative to a sign-free run
is ×2.1 at β = 12, ×2.6 at β = 16 and ×3.4 at β = 20. On top of that comes the growth of ⟨k⟩ ∝ β, which costs
⟨k⟩² to ⟨k⟩³. At this doping the 2×2 plaquette is therefore usable down to T ≈ t/20. Larger clusters scale with
N_c in the exponent.

## 2. Auxiliary shift δ (α = ½ ± δ), pilot runs at β = 8

| δ | ⟨s⟩ | ⟨k⟩ | time (s) | efficiency ⟨s⟩²/t |
|---|---|---|---|---|
| 0.30 | 0.10 ± 0.01 | 42 | 2.4 | 0.004 |
| **0.51** | 0.74 ± 0.05 | 76 | 5.3 | **0.102** |
| 0.75 | 0.77 ± 0.08 | 136 | 14.9 | 0.040 |
| 1.00 | 0.85 ± 0.05 | 216 | 38.4 | 0.019 |

* A larger δ raises the sign, but ⟨k⟩ ∝ δ makes each step costlier. The net figure of merit peaks at the standard δ ≈ 0.5.
* δ < 0.5 is a trap: the sign collapses to 0.10.
* `optimize_delta` makes this choice automatically for any cluster, filling and temperature.

## 3. Error bars under the sign

The old estimator divided per chain (⟨sS⟩_c/⟨s⟩_c, then averaged), which is biased when ⟨s⟩ is small. `solve` now pools ⟨sS⟩ and ⟨s⟩ over all chains and propagates both, with their correlation, by a jackknife over n_bins × n_chains bins.

Check at β = 8 with 8 independent runs: the spread of Im Σ(π,0) is 0.131 and the mean jackknife error is 0.100 (Re: 0.119 vs 0.093). That is within the ±25% uncertainty of a spread estimated from 8 samples, with a mild underestimate from autocorrelation across bins. Use more bins or longer runs when ⟨s⟩ < 0.5.
