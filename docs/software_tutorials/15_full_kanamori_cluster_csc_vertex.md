# Tutorial 15 — Full-Kanamori CT-HYB, cluster CT-INT, charge self-consistency and two-particle vertices

Version 1.1 closes the solver-side limitations of v1.0. Every statement below is backed by a row of
`benchmarks/RESULTS_v11.md`.

## 1. Matrix CT-HYB: any local Hamiltonian (`qmc.cthyb_matrix.MatrixCTHYB`)

The segment solver (Tutorial 10) handles only density–density interactions. The matrix solver evaluates the
full local trace in the eigenbasis of H_loc:

  w(C) = Π_f det F_f · s_P Tr[e^{-(β-t₁)H} O₁ e^{-(t₁-t₂)H} O₂ … O_K e^{-t_K H}]

* F(τ) = −Δ(β−τ) = Δ(−τ) makes all single-orbital weights positive.
* s_P is the parity of the permutation that sorts the operators by decreasing time.
* Spin-flip and pair-hopping terms of the Kanamori interaction are therefore included exactly.

```python
from feynlab.qmc.cthyb_matrix import MatrixCTHYB
s = MatrixCTHYB(grid, n_orb=2, U=2.5, J=0.4, eps=[e0, e1], n_iw_meas=32)
r = s.solve(Delta_flavours, n_meas=1_500_000, n_chains=8, parallel="process")
r["G"], r["n_0"], r["docc_0"], r["Sz^2"]      # observables are measured by operator insertion at tau = 0
```

**Validation.** For two orbitals with spin flip and pair hopping, compared with multi-orbital ED:

* the average sign is 1;
* n, ⟨n↑n↓⟩, ⟨S_z²⟩ and G(iω₀) agree within one standard error.

In DMFT on a one-band lattice it gives Z = 0.250, against 0.247 from the segment solver.

**Two bugs worth learning from** (both caught by the ED comparison):

* *Skipped measurements.* A `continue` statement on impossible moves (removal at k = 0) skipped the
  measurement. That under-samples low orders and inflates G by 10%. Every Monte Carlo step must reach the
  measurement.
* *Stale ordering.* After a rejected removal the swapped operators must be restored. Otherwise the stored
  determinant and trace refer to different orderings, and the measured sign drops to 0.77.

**Inside DFT+DMFT** use `DFTDMFT(...).solve('cthyb', cthyb_kwargs=...)`. QMC noise in Σ grows like δG/|G|²,
so the solver keeps only the frequencies with a small propagated error. Beyond those it uses the exact Kanamori
Hartree–Fock tail Σ∞ = U n_{a,−σ} + Σ_b [U′ n_{b,−σ} + (U′−J) n_{b,σ}].

## 2. Cluster CT-INT and exact DCA (`qmc.ctint_cluster`, `DCA2D.solve(solver='ctint')`)

Vertices carry a cluster site. The determinant matrices use the cluster-translation-invariant Weiss field
G̃₀(r_p − r_q, τ_p − τ_q). G(K) is measured directly in momentum space.

* Validation against ED of a dimer in a bath (256 states): Σ(K, iω₀) agrees within one standard error for
  K = 0 and π.
* For N_c = 1 the DCA loop reproduces single-site CT-INT DMFT (Z = 0.4721 vs 0.4709).
* N_c = 4 now uses a numerically exact cluster solver instead of second-order perturbation theory.

## 3. Charge self-consistent DFT+DMFT (`DFTDMFT.solve_csc`)

```python
from feynlab.dmft.dftdmft import DFTDMFT, hartree_functional
m = DFTDMFT("wannier90_hr.dat", corr=[0, 1, 2], U=4, J=0.6, n_elec=2, beta=40, nk=(16, 16, 16))
r = m.solve_csc(hartree_functional(U_H), n_outer=40, mix_rho=0.5, solver="hubbard1")
m.write_density_matrix("N_k.h5")     # N(k) for an external DFT code (full ab-initio cycle)
```

* `density_matrix_k` returns N(k)_ab = T Σ_n G_ab(k, iω_n) e^{iω_n0⁺} with the analytic 1/(iω)² tail.
  At U = 0 it matches the Fermi occupations to 3·10⁻⁷.
* The outer loop updates H(k) = H₀(k) + V[n] and re-evaluates the double counting.
* With U_DMFT = 0 and a Hartree functional the loop must reproduce static mean field. It does: occupations agree
  to 6·10⁻⁶ with an independent SCF code.
* **Fixed along the way:** the electron count in `DFTDMFT` now uses the exact tail coefficient
  c₂ = ⟨H⟩ − μ + Re Σ(i∞). Without it, μ was set from a density that was 1.2% off.

## 4. Two-particle vertex from ED and DMFT susceptibilities (`feynlab.vertex`)

`ImpurityED(...).g2(beta, n_nu, n_om)` evaluates

  G2_{σσ′}(ν, ν′, ω) = ∫dτ₁dτ₂dτ₃ e^{−iντ₁} e^{i(ν+ω)τ₂} e^{−i(ν′+ω)τ₃} ⟨T c†_σ c_σ c†_σ′ c_σ′(0)⟩

exactly in the Lehmann representation:

* For each time ordering, the integral over the time simplex is the divided difference of e^{βz} at four nodes.
* It is evaluated as the corner element of exp(βJ) with J bidiagonal, which stays stable for degenerate energies.

Checks:

| check | result |
|---|---|
| U = 0 (Wick) | χ_↑↑ = −βG(ν)G(ν+ω)δ_νν′ and χ_↑↓ = 0 to 1e-15 |
| Schwinger–Dyson equation with F_↑↓ = (F_d − F_m)/2 | reproduces the ED self-energy to 0.05% at iν₀ |
| small U | F_d → +U, F_m → −U |
| DMFT uniform susceptibility (Bethe–Salpeter, 1/N_ν extrapolation) | χ_m(q=0) = 2.9089 vs dm/dh = 2.9102: the thermodynamic consistency of DMFT |

```python
from feynlab.vertex import ed_dmft, vertices, lattice_susceptibility, ladder_dga_sigma
ed  = ed_dmft(beta, U, n_bath=2, eps_k=eps_k)            # converged DMFT impurity (ED bath)
gs  = ed.generalized_susceptibilities(beta, n_nu=12, n_om=8)
v   = vertices(gs, beta)                                 # F_d, F_m, Gamma_d, Gamma_m
chi = lattice_susceptibility(v["Gamma_m"], chi0_q, beta) # DMFT chi_m(q)
S   = ladder_dga_sigma(eps_k, ed.mu, lambda x: ed.sigma(1j*x, beta), v, gs, beta, U, 1.0)
```

**Ladder DΓA.** The formula is
Σ(k,ν) = Σ_DMFT(ν) − (U/β²) Σ [(½F_d(q) − 3/2 F_m(q)) − (½F_d − 3/2 F_m)_loc] GGG.

* It combines the ph and transverse ph-bar ladders and removes the local double counting
  (Toschi, Katanin, Held 2007).
* At weak coupling it reduces to the non-local second-order self-energy; see `RESULTS_v11.md`.

## 5. TPSC (`vertex.tpsc`, `tpsc_sigma`)

The two-particle self-consistent approach of Vilk and Tremblay:

* U_sp and U_ch are fixed by the local spin and charge sum rules.
* The ansatz U_sp = U⟨n↑n↓⟩/(⟨n↑⟩⟨n↓⟩) closes the equations.
* The bosonic 1/Ω² tail is added analytically, so U → 0 gives ⟨n↑n↓⟩ = 1/4 exactly.

Exact internal identity: (T/N) Σ_k Σ(k) G⁰(k) = U⟨n↑n↓⟩. The residual halves when the frequency cutoff doubles
(1/N truncation) and extrapolates to zero.

## 6. Parallel Markov chains (`feynlab.utils.parallel`)

Every QMC solver (`CTINTSolver`, `CTHYBSolver`, `MatrixCTHYB`, `ClusterCTINT`) takes `parallel='serial' | 'process' | 'mpi' | 'auto'`.
Results are identical across backends (same seeds). Under `mpirun -n N`, chains are distributed over ranks
with mpi4py.

Quiz: **T16**.
