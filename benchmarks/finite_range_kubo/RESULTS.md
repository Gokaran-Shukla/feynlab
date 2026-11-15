# Finite-range disorder in the Kubo–Středa route (`transport.finite_range.FiniteRangeHall`)

**Model.** Massive Dirac model, v = 1, Δ = 0.5, μ = 1 (k_F = 0.866).

**Disorder.** Gaussian impurities u(q) = e^{−q²d²/4}, with n_i V0² = 0.01, the same form factor as
`SemiclassicalHall(range_d=d)`.

**Method.** Non-crossing approximation: self-consistent Born Σ(k) and a k-dependent ladder Γ(k), solved with GMRES.
Grid: 1024² Cartesian points on [−2, 2]², with convolutions by zero-padded FFTs. About 1 min per d.

| d (k_F d) | σ_xx Kubo | σ_xx Boltzmann (Bessel, analytic) | extrinsic AHE Kubo (NCA) | semiclassical sj + adist + Gaussian skew | Γ |
|---|---|---|---|---|---|
| 0 | 85.87 | 85.71 (0.18%) | −0.5636 | −0.5663 (0.5%) | 0.0037 |
| 0.5 (0.43) | 108.84 | 108.48 (0.34%) | −0.5443 | −0.5460 (0.3%) | 0.0034 |
| 1.0 (0.87) | 206.06 | 205.18 (0.43%) | −0.4805 | −0.4801 (0.07%) | 0.0026 |
| 2.0 (1.73) | 1040.3 | 1037.4 (0.28%) | −0.2884 | −0.2867 (0.6%) | 0.0010 |
| 3.0 (2.60) | 3402.8 | 3390.1 (0.38%) | −0.2239 | −0.2229 (0.5%) | 0.0006 |

All values are in e²/h.

* **Delta limit.** With d = 0, the code reproduces `DisorderedHall` on the same grid to 1e-9. This is a unit test.
* **Diagrammatic vs semiclassical.** The diagrammatic and semiclassical routes now agree for finite-range
  disorder within 0.6%, which is the expected O(Γ/μ) level.
* **Physics.** Smooth disorder suppresses the Gaussian skew strongly (−0.138 → −0.017 from d = 0 to d = 3). The
  side jump and the anomalous distribution decrease more slowly (−0.214 → −0.103). The total extrinsic AHE falls
  by a factor 2.5, so the intrinsic part (−0.25) dominates in clean samples with smooth disorder.

**Not covered.** The third-order (V1³) skew vertex and the crossed X/Ψ diagrams are still δ-correlated only.
