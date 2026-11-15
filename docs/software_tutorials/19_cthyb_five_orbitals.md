# Tutorial 19 — Full-Kanamori CT-HYB for 5 orbitals (`SectorCTHYB`)

`MatrixCTHYB` multiplies dense 4^{n_orb} × 4^{n_orb} matrices for every operator in the local trace. That is fine
for 2 orbitals (16 states), slow for 3 (64), and impossible for a full d shell (1024 states).
`qmc.cthyb_sector.SectorCTHYB` samples the same configuration space with the same weights and estimators, but
evaluates the trace the way production codes do.

## 1. The trace

    Tr[e^{-(β-τ_1)H} O_1 … O_K e^{-τ_K H}] = Σ_{s0} tr_{s0}[ … ]

* **Blocks.** H_loc (Kanamori) conserves N↑ and N↓, so it is diagonalized inside (N↑, N↓) blocks. For 5 orbitals
  there are 36 blocks with at most 100 states each. Each c_f or c†_f maps one block to exactly one block.
* **Viability.** For a given operator string, the starting block s0 fixes every intermediate block. Strings that
  annihilate s0 are discarded with integer look-ups only.
* **Bounds.** A rigorous bound per block,
  `|tr_{s0}| ≤ d_{s0} Π ‖B‖ e^{−Δτ E_min}`, gives three savings:
  * negligible blocks are skipped;
  * in each kept block, only outer eigenstates with e^{−ℓ(E_n − E_min)} above tolerance are propagated
    (ℓ is the wrap-around interval), so the products are matrix–vector instead of matrix–matrix;
  * a move whose bound cannot pass the Metropolis test is rejected without computing the trace (lazy trace).
* **Truncation (optional).** `e_cut` drops atomic eigenstates with E − E_min > e_cut (Haule 2007). Check
  convergence in e_cut.

## 2. Usage

```python
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.qmc.cthyb_sector import SectorCTHYB
g = MatsubaraGrid(10.0, 512)
s = SectorCTHYB(g, n_orb=5, U=4.0, J=0.6, mu=16.3, e_cut=10.0)
r = s.solve(delta_iw, n_meas=300_000, n_chains=4)      # delta_iw: (10, n_iw) or (n_iw,)
r["G"], r["n_0"], r["docc_0"], r["Sz^2"], r["lazy_rejected"], r["n_states_kept"]
```

## 3. Validation (`benchmarks/cthyb_5orb/RESULTS.md`)

* **Trace.** Equal to the dense trace to 1e-10 on random strings; the bound always holds.
* **5 orbitals vs exact ED.** The test impurity has one bath level per orbital (2²⁰ states).
  * ⟨S_z²⟩: 1.936 ± 0.008 vs 1.933 exact.
  * ⟨n↑n↓⟩: 0.2169 ± 0.0012 vs 0.2179 exact.
  * n: within 0.3%.
  * e_cut = 10 keeps 663 of the 1024 states, is 25% faster, and stays within errors.
* **Speed.** 15× faster than the dense solver at 3 orbitals; about 1 ms per step at 5 orbitals.

## 4. Limits

* **Green function estimator.** G(iω) uses the standard M-matrix estimator. It is very noisy when the expansion
  order is low and Δ(τ) varies strongly across [0, β]. Use the local observables, or worm sampling (not
  implemented yet).
* **Hybridization.** Δ must be diagonal in flavour (no off-diagonal hybridization).
* **Interaction.** H_loc must conserve (N↑, N↓). Spin–orbit coupling needs the dense solver.
