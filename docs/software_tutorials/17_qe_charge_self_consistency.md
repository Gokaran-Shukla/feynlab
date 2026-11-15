# Tutorial 17 — Fully charge self-consistent DFT+DMFT with Quantum ESPRESSO and Wannier90

This closes the last DFT-side limitation: feynlab now drives a real first-principles cycle. Quantum ESPRESSO
supplies the Kohn–Sham bands and the density, Wannier90 the projectors, and feynlab the DMFT.

## 1. The interface QE provides (QE ≥ 7.0, built with HDF5)

`pw.x` with `dmft = .true.` reads `<dmft_prefix>.h5` (source: `PW/src/add_dmft_occ.f90`, S. Beck):

| dataset | h5py shape | meaning |
|---|---|---|
| `dft_misc_input/band_window` | (1, n_k, 2) | first and last band (1-based) of the correlated window |
| `dft_update/delta_N` | (n_k, n_b, n_b, 2) | ΔN_ij(k) = N^DMFT_ij(k) − δ_ij f^DFT_i(k), real and imaginary parts, per spin |

What QE does with it:

* adds ΔN to its Fermi weights;
* diagonalises f + ΔN = v λ v†;
* rotates the Kohn–Sham states ψ′_a = Σ_m v_ma ψ_m and occupies them with λ_a.

The density is therefore Σ_mn ψ_m ψ_n* N_mn, with N_mn = ⟨c†_n c_m⟩. This is exactly feynlab's convention
(G_mn(τ) = −⟨T c_m c†_n⟩). Because QE transposes the matrix on reading, the dataset stores D[k, m, n] = ΔN_mn.
The test `test_qe_dmft_h5_layout` replays QE's Fortran index order to check this.

## 2. The cycle (`feynlab.dft.qe_csc.QEDMFT`)

```
scf  ──►  nscf (dmft=.true., nosym, explicit k)  ──►  wannier90 -pp / pw2wannier90 / wannier90  ──►  U(k)
  ▲                                                                                                │
  │       pw.x scf, restart, electron_maxstep = 1   ◄──  <seed>_dmft.h5 (ΔN)  ◄──  feynlab DMFT ◄──┘
```

1. H_W(k) = V(k)† diag(ε_k) V(k), with V = U_dis U (window → Wannier). Disentanglement is supported.
2. DMFT on H_W (`DFTDMFT`, solvers `hubbard1`, `ed`, `lanczos` or `cthyb`) gives Σ − Σ_dc. The DMFT self-energy
   is warm-started from the previous CSC cycle, so later cycles need only a few DMFT iterations.
3. The window density matrix is
   N(k) = T Σ_n [iω_n + μ − ε_k − V(Σ − Σ_dc)V†]⁻¹ e^{iω_n 0⁺}, with the exact 1/(iω)² tail.
   μ is fixed so that the window charge equals QE's.
4. ΔN = N − diag(f^DFT) is written to the archive; then QE performs one density update.

**Temperature consistency.** Use `smearing = 'fermi-dirac'` with `degauss` = k_BT (in Ry) of the DMFT. Then with
Σ = 0 the DMFT occupations are QE's own occupations, and ΔN vanishes identically.

## 3. Validation on SrVO₃ (PBE ultrasoft pseudopotentials generated with `ld1.x` from PSlibrary)

| check | result |
|---|---|
| Wannier90 t₂g functions | centred on V, spread 1.885 Å² each (3 isolated bands, 21–23) |
| U = 0: ΔN | max 1.5·10⁻⁷ (Matsubara truncation); μ_DMFT = E_F^QE |
| U = 0: QE reads the archive and updates | total energy −317.45621426 Ry vs DFT −317.45621425 Ry; scf accuracy 1.4·10⁻⁹ Ry |
| U = 5 eV, J = 0.75 eV, β = 40 eV⁻¹ (Lanczos ED, 2 bath sites per orbital) | see §3a and `benchmarks/qe_csc/RESULTS.md` |

### 3a. SrVO₃, U = 5 eV, J = 0.75 eV, β = 40 eV⁻¹ (Lanczos, 2 bath sites per orbital)

| cycle | Z (t2g) | E_F (eV) | QE scf accuracy (Ry) | E (Ry) | ΔE (Ry) | window charge |
|---|---|---|---|---|---|---|
| 0 | 0.428 | 12.7323 | 1.12e-04 | -317.427970 | — | 1.00000 |
| 1 | 0.430 | 12.7368 | 2.93e-05 | -317.428114 | -1.4e-04 | 1.00000 |
| 2 | 0.430 | 12.7455 | 1.19e-05 | -317.428157 | -4.3e-05 | 1.00000 |
| 3 | 0.431 | 12.7413 | 4.86e-06 | -317.428181 | -2.4e-05 | 1.00000 |
| 4 | 0.431 | 12.7448 | 2.40e-06 | -317.428186 | -5.0e-06 | 1.00000 |
| 5 | 0.431 | 12.7420 | 1.10e-06 | -317.428192 | -6.0e-06 | 1.00000 |
| 6 | 0.431 | 12.7440 | 6.00e-07 | -317.428191 | +1.0e-06 | 1.00000 |


* QE's scf accuracy falls geometrically, by a factor ≈ 2.2 per cycle, to 6·10⁻⁷ Ry after 7 cycles.
* E_F shifts by +12 meV.
* Z goes from 0.428 (one-shot) to 0.431. CSC effects are small for this isolated, symmetry-pinned t₂g manifold.
* Figure and details: `benchmarks/qe_csc/RESULTS.md`.

**Lesson from building it.** A fitted Δ(iω) = Σ V²/(iω − e) already measures the bath levels e from μ. Subtracting μ
again (the pre-v1.2 `DFTDMFT` ED solver did this) silently empties the bath as μ moves. With one bath site per
orbital the result looked like a Mott transition.

## 4. Running it

```bash
micromamba create -p /opt/qe -c conda-forge qe wannier90     # QE 7.6 (HDF5 + DMFT hook), Wannier90 4
feynlab example qe-csc > csc.toml                            # edit structure, window, U, J, beta
feynlab qe-csc csc.toml
```

Or from Python:

```python
from feynlab.dft.qe_csc import QEDMFT
q = QEDMFT("svo", open("scf.in").read(), (8, 8, 8), (21, 23), ["V:dxz;dyz;dxy"], U=5.0, J=0.75, beta=40.0,
           nbnd=30, solver="lanczos", solver_kwargs=dict(n_bath=2),
           pw_cmd="mpirun -np 4 pw.x -nk 4", pw2wan_cmd="mpirun -np 4 pw2wannier90.x")
history = q.cycle(n_csc=14, lattice_rows_A=a * np.eye(3), atoms_frac=atoms)
```

**Convergence criteria.** The loop stops on QE's estimated SCF accuracy, the change of the total energy and the
change of the correlated-window charge. Do not monitor the per-orbital t₂g occupations alone: in cubic SrVO₃
symmetry pins them at n/3.

## 5. Limits

* nspin = 1 only (QE: "DMFT update not implemented with nspin > 1"). Magnetic CSC needs a QE extension.
* The DMFT total energy (Galitskii–Migdal plus double-counting correction) is not yet assembled. QE's printed
  energy is the DFT functional evaluated at the interacting density, which is useful for monitoring convergence.
* The impurity solver sets the physics. ED with one bath site per orbital is too crude for t₂g: it drifts
  towards a spurious Mott state (Z ≈ 0.01). Use `lanczos` with ≥ 2 bath sites per orbital or `cthyb`.
* Each CSC cycle costs one NSCF + Wannier90 + a few DMFT iterations + one QE density step
  (≈ 20 min for SrVO₃ on 2 cores).

Reference: S. Beck, A. Hampel, O. Parcollet, C. Ederer, A. Georges, J. Phys.: Condens. Matter 34, 235601 (2022).

Quiz: **T18**.
