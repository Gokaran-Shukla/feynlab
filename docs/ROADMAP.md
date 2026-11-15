# feynlab roadmap — what would make it powerful and robust

> **v1.1 status:** all 22 items below are implemented, including the parts that were partial in v1.0
> (vertex ladders / DΓA / TPSC, full-Kanamori matrix CT-HYB, cluster CT-INT, charge self-consistency, Holstein/SSH
> polarons, finite-range disorder, crossed diagrams for any 2D/3D Wannier lattice, MPI chains).
> See `docs/software_tutorials/07_benchmarks_and_limitations.md`, `benchmarks/RESULTS_v1.md`, `benchmarks/RESULTS_v11.md`,
> tutorials 09–16 and quizzes T10–T17. **Next:** λ-corrected DΓA and dual fermions, worm-sampling vertices from CT-HYB,
> sliding-window/Krylov matrix CT-HYB for 5 orbitals, an interface to a DFT code (VASP/QE) for ab-initio CSC,
> finite-range disorder in the Kubo–Středa route.

Prioritised by (scientific impact) × (fit with the existing code). Items marked ★ connect directly to
orbitronics / correlated-materials research (Wannier-based OHC/SHC, kagome metals, altermagnets).

## A. Correlated-electron solvers (accuracy)
1. **CT-HYB impurity solver** (hybridization expansion, Werner et al. PRL 97, 076405 (2006)) — the workhorse for
   strong coupling and multi-orbital DMFT; complements CT-INT (weak/intermediate U).
2. **Multi-orbital interaction (Slater–Kanamori, U, U′, J)** in ED, CT-INT and CT-HYB → Hund's metals, ★ kagome
   CsV₃Sb₅/FeGe, t₂g oxides.
3. **Cluster DMFT / DCA** (2×2 plaquette) → d-wave, pseudogap; benchmark vs Maier et al. RMP 77, 1027 (2005).
4. **DFT+DMFT interface** to Wannier90 (`hr.dat`, projections) with double-counting options and charge
   self-consistency hooks.
5. **Analytic continuation**: stochastic (Mishchenko/Sandvik), Nevanlinna (Fei et al. PRL 126, 056402 (2021)).
6. **Two-particle quantities**: generalized susceptibilities, vertex Γ, DΓA / dual fermions → momentum-dependent
   self-energy beyond DMFT; ★ magnetic fluctuations in kagome/altermagnets.
7. **NRG impurity solver** (real-frequency, T = 0) — removes the bath-size issue of ED entirely.

## B. Diagrammatic Monte Carlo (the "Feynman diagram" core)
8. **Connected-determinant DiagMC (CDet)** (Rossi PRL 119, 045701 (2017)) — sums all connected diagrams of the
   2D Hubbard model order by order; benchmark vs Šimkovic et al. Science 2022-level data.
9. **Bold DiagMC with sampling over orders** (not just ≤4) + the Van Houcke et al. SciPost Phys. 16, 133 (2024)
   criterion for detecting misleading convergence.
10. **Resummation tools**: conformal maps / Padé of the perturbation series, radius-of-convergence analysis.
11. **Polaron extensions**: effective mass, Z-factor, lattice (Holstein, SSH) polarons, exact energy estimator.

## C. Transport (★ orbitronics)
12. ★ **Wannier90 + full 3D BZ** for OHC / SHC / quadrupole Hall with adaptive k-mesh refinement and the
    modern-theory orbital moment (beyond ACA) — merge with the CUDA OHC calculators.
13. ★ **Disorder models beyond δ-scalar**: spin-dependent / orbital-dependent impurities, finite-range (correlated)
    disorder (Ado et al. PRB 96, 235148 (2017)), T-matrix for strong scatterers (resonant skew scattering).
14. ★ **Crossed diagrams for lattice models** (FFT on the BZ) and for spin/orbital currents (orbital side jump /
    skew scattering — an open research topic).
15. **Nonlinear Hall / Berry-curvature dipole**, Edelstein effect, spin–orbit torque (Kubo–Bastin with vertex
    corrections) — ★ CIPS, WSe₂, MoO₂/Ni, NbP/Ni projects.
16. **Kubo–Bastin at finite temperature** and Fermi-sea terms for general models (not only σᴵᴵ = intrinsic − σᴵ).

## D. Robustness & engineering
17. **GPU back-end** (CuPy) for k-sums, FFTs and Lanczos; MPI / multiprocessing for Monte Carlo chains.
18. **HDF5 I/O + checkpointing** for long DMFT and QMC runs; restart files.
19. **Error analysis**: binning / autocorrelation (Γ-method), jackknife for derived quantities, automatic
    convergence diagnostics in DMFT loops (Anderson / Broyden mixing instead of linear mixing).
20. **Continuous integration with the full benchmark suite** on a schedule, coverage, type hints, docs site
    (Sphinx/MkDocs) with the tutorials rendered, `pip install feynlab` release, Zenodo DOI.
21. **Command-line interface & input files** (`feynlab dmft input.toml`) for reproducible runs.
22. **Units / conventions layer** (eV, Å, SI conductivities per layer or per volume) for comparison with experiment.

## Suggested next three
1. CT-HYB + multi-orbital Kanamori (opens real-materials DMFT).
2. Wannier90 3D OHC/SHC with vertex corrections and crossed diagrams (unique — nobody ships orbital skew
   scattering / orbital side jump today).
3. Connected-determinant DiagMC for the 2D Hubbard model (the strongest "Feynman-diagram" capability).
