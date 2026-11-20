# Tutorial 10 — CT-HYB, Kanamori multi-orbital physics and DFT+DMFT

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


Roadmap items A1, A2 and A4.

## 1. CT-HYB in the segment picture (`feynlab.qmc.cthyb`)

CT-INT expands in U. CT-HYB instead expands the partition function in the **hybridization** Δ(τ):

  Z = Z_bath Σ_k ∫ dτ₁…dτ'_k  det[Δ]↑ det[Δ]↓ · Tr_loc[T_τ c(τ'_k)c†(τ_k)…]

For a density–density interaction, each spin's operators form **segments** (occupied time intervals).
The local trace is then exp(−U·overlap + μ·length), with no matrix algebra
(Werner, Comanac, de' Medici, Troyer, Millis PRL 97, 076405 (2006)). The perturbation order is ~β·(kinetic energy),
so it **decreases** with U. This makes CT-HYB the method of choice at strong coupling, where CT-INT orders explode.

Moves implemented: insert/remove a segment and an anti-segment, both with detailed balance. Δ(τ) enters through
F(τ) = Δ(β−τ). The Green function is measured by the estimator G(τ) = −(1/β)⟨Σ M_ji δ(τ − (τ_i − τ'_j))⟩.

```python
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.qmc.cthyb import CTHYBSolver
from feynlab.qmc.ed import AndersonED
g = MatsubaraGrid(10.0, 1024)
ed = AndersonED(-2.0, 4.0, [-1.0, 0.0, 1.2], [0.5, 0.4, 0.6])
r = CTHYBSolver(g, U=4.0, n_iw_meas=6).solve(ed.hybridization(g.iw), mu=2.0, n_meas=1_000_000, n_chains=4)
print(r["docc"], ed.double_occupancy(10.0))      # 0.0780 vs 0.0782
```

Benchmarks (`RESULTS_v1.md`): occupation, double occupancy and Im Σ(iω₀) agree with ED within 3σ at U = 3 and 4.
In DMFT, `BetheDMFT(20, 2.0).solve('cthyb', ...)` gives Z = 0.235 vs 0.229 from finite-T ED-DMFT.

## 2. Kanamori interaction (`feynlab.qmc.multiorbital`)

  H_int = U Σ_a n_a↑n_a↓ + U′ Σ_{a≠b} n_a↑n_b↓ + (U′−J) Σ_{a<b,σ} n_aσ n_bσ − J (spin flip) + J (pair hopping),
  U′ = U − 2J (rotational invariance).

`KanamoriAtom(n_orb, U, J)` diagonalises the atom sector by sector. Its multiplets are the textbook result
(Georges, de' Medici, Mravlje, Annu. Rev. CMP 4, 137 (2013)):

| filling | multiplet energies |
|---|---|
| 2 orbitals, N = 2 | U − 3J (³T, ×3), U − J (×2), U + J (×1) |
| 3 orbitals (t₂g), N = 3 | ground state ⁴S at 3U − 9J, 4-fold degenerate |

`MultiOrbitalED(n_orb, U, J, eps, bath_eps, bath_V)` is the impurity solver. It uses sparse Jordan–Wigner
operators, blocks by N, and handles at most 8 sites. In its one-orbital limit it reproduces `AndersonED`
to 1e-14. `hubbard_one_sigma` is the multi-orbital Hubbard-I self-energy built from the atomic Green function.

**Physics to explore: Hund's metals.** At fixed U, increasing J lowers Z away from half filling
(the "Janus" effect). Try `DFTDMFT(..., n_elec=2, J=0..0.25U)` with the three-orbital model.

## 3. DFT+DMFT from Wannier90 (`feynlab.dmft.dftdmft`)

```python
from feynlab.dmft.dftdmft import DFTDMFT
m = DFTDMFT("wannier90_hr.dat", corr=[0, 1, 2], U=4.0, J=0.6, n_elec=2.0, beta=40, nk=(24, 24, 24), dc="FLL")
r = m.solve("hubbard1", n_iter=80, mix=0.3)    # or "ed"
print(r["Z"], r["n_total"])
```

* H(k) = Σ_R e^{ik·R} H_R / deg(R), read from `hr.dat` (`read_hr` / `write_hr`).
* Lattice Dyson equation G_loc = Σ_k [iω + μ − H(k) − Σ + Σ_dc]⁻¹, with μ fixed to `n_elec` by bisection.
* Double counting:
  * `FLL` = U(n − ½) − J(n_σ − ½)
  * `AMF`
  * `occupation`, which fixes Σ_dc so that Hubbard-I reproduces the DFT occupation, at the **centre** of the
    occupation plateau. In a Mott gap any shift inside the gap gives the same n, so bisection must not stop at
    the plateau edge.
* Validation: U = 0 gives Z = 1 and n_total = 2 exactly. U = 4, J = 0.5 at half filling gives a Mott insulator
  with Z ≈ 0.002.

**Limitations.** No charge self-consistency, and the hybridization-expansion solver is restricted to
density–density terms.

Quiz: **T11**.
