# Tutorial 22 — Living with the CT-INT sign problem

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


## Why it cannot be removed

Away from half filling, a fermionic cluster has configurations with negative weight. QMC then samples |w|, and every
observable becomes a ratio:

    ⟨O⟩ = ⟨s O⟩_{|w|} / ⟨s⟩_{|w|}

⟨s⟩ = Z/Z_{|w|} is a ratio of partition functions, so it is exponentially small in β × N_c. For a fixed error the
cost grows as 1/⟨s⟩². A generic cure would solve NP-hard problems (Troyer & Wiese, PRL 94, 170201 (2005)).
What a code can do is measure the sign honestly, minimize its effect, and tell you the price in advance.

## Three tools (`feynlab.qmc.ctint_cluster`)

1. **`ClusterCTINT.solve(..., n_bins=10)`.** It pools ⟨sS⟩ and ⟨s⟩ over all chains, so no per-chain ratio bias
   appears at small ⟨s⟩. It returns `sign`, `sign_err` and jackknife errors of Σ(K, iω) that include the sign
   fluctuations.
2. **`optimize_delta(factory, G0_K, deltas)`.** It runs pilot runs over the auxiliary shift δ and ranks them by
   efficiency ⟨s⟩²/t. A larger δ raises the sign but also the expansion order, so the optimum must be measured.
   δ < ½ is dangerous.
3. **`sign_cost_model(betas, signs, N_c)`.** It fits ⟨s⟩ = e^{−βN_cΔf} and predicts the sign and the relative
   cost 1/⟨s⟩² at lower T before you commit CPU time.

```python
from feynlab.qmc.ctint_cluster import ClusterCTINT, optimize_delta, sign_cost_model
best, table = optimize_delta(lambda d: ClusterCTINT(g, U, sites, L, K, delta=d), G0_K)
r = ClusterCTINT(g, U, sites, L, K, delta=best).solve(G0_K, n_bins=20)
df, predict = sign_cost_model([2, 4, 6, 8], [0.9998, 0.967, 0.883, 0.759], 4); predict(16.0)
```

## Results (`benchmarks/ctint_sign/RESULTS.md`, 2×2 plaquette, U = 6t, hole doped)

* **Sign vs temperature.** ⟨s⟩ falls from 0.9998 at β = 2 to 0.70 at β = 10, with Δf = 0.0076. The predicted
  extra cost is ×3.4 at β = 20.
* **Choice of δ.** δ = 0.51 is optimal (efficiency 0.10). δ = 1 has a better sign (0.85) but is 5× less
  efficient. δ = 0.3 collapses the sign to 0.10.
* **Error bars.** The jackknife errors reproduce the spread of independent runs within the uncertainty of the
  test (0.10 vs 0.13).

## Limits

* The sign problem itself remains. Large clusters (N_c ≥ 8) at low T and finite doping stay exponentially
  expensive.
* Submatrix updates and basis-optimized formulations are not implemented.
