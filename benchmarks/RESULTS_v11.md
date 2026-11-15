# feynlab v1.1 benchmark results (closing the v1.0 limitations)

42/43 checks passed, plus 3 informational rows (runtime 476 s)

| group | check | feynlab | reference | error | tol | pass | source |
|---|---|---|---|---|---|---|---|
| CT-HYB-mat | 2-orbital Kanamori (spin flip + pair hopping): average sign | 1 | 1 | 0.00e+00 | 1.0e-12 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| CT-HYB-mat | orbital 0: n_up vs ED | 0.49834994 | 0.49876039 | 4.10e-04 | 2.0e-03 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| CT-HYB-mat | orbital 0: <n_up n_dn> vs ED | 0.094194182 | 0.09438998 | 1.96e-04 | 2.2e-03 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| CT-HYB-mat | orbital 1: n_up vs ED | 0.44649362 | 0.44596078 | 5.33e-04 | 4.3e-03 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| CT-HYB-mat | orbital 1: <n_up n_dn> vs ED | 0.058191221 | 0.058061413 | 1.30e-04 | 2.0e-03 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| CT-HYB-mat | <S_z^2> (Hund's moment) vs ED | 0.48097681 | 0.48104527 | 6.85e-05 | 2.0e-03 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| CT-HYB-mat | Im G(i w_0), orbital 0, vs ED | -0.58044829 | -0.58007169 | 3.77e-04 | 1.0e-02 | yes | Werner & Millis PRB 74, 155107 (2006); Haule PRB 75, 155113 (2007) vs ED |
| ClusterCTINT | dimer+bath, K=0: Im Sigma(i w_0) vs ED | -0.30457999 | -0.30446443 | 1.16e-04 | 7.3e-03 | yes | Gull et al. PRB 83, 075122 (2011) vs ED |
| ClusterCTINT | dimer+bath, K=0: Re Sigma(i w_0) vs ED | 0.99193337 | 0.99252681 | 5.93e-04 | 3.0e-03 | yes | vs ED |
| ClusterCTINT | dimer+bath, K=pi: Im Sigma(i w_0) vs ED | -0.29325238 | -0.29324758 | 4.80e-06 | 7.2e-03 | yes | Gull et al. PRB 83, 075122 (2011) vs ED |
| ClusterCTINT | dimer+bath, K=pi: Re Sigma(i w_0) vs ED | 1.1497055 | 1.1493445 | 3.61e-04 | 3.0e-03 | yes | vs ED |
| DFT+DMFT | Wannier density matrix N(k) vs exact Fermi occupations, max error | 2.9673202e-07 | 0 | 2.97e-07 | 1.0e-05 | yes | exact (U = 0) |
| DFT+DMFT | charge self-consistency (U_DMFT = 0): n_orbital 0 vs Hartree mean field | 1.3988402 | 1.398834 | 6.25e-06 | 1.0e-04 | yes | independent SCF reference |
| DFT+DMFT | charge self-consistency: total electron number | 1.4000001 | 1.4 | 7.85e-08 | 1.0e-04 | yes | electron count |
| DFT+DMFT | one-band lattice: Z (matrix CT-HYB in DFT+DMFT) vs segment CT-HYB | 0.24976061 | 0.24699241 | 2.77e-03 | 2.0e-02 | yes | two independent CT-HYB implementations |
| DCA | N_c = 1, cluster CT-INT vs single-site CT-INT DMFT (Z) | 0.46958868 | 0.4685014 | 1.09e-03 | 1.0e-02 | yes | Hettler et al. (1998) |
| DCA | N_c = 4 (exact cluster solver): Z(pi,0) < Z(0,0) | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | Maier et al. RMP 77, 1027 (2005) |
| DCA | N_c = 4 CT-INT: Z(pi,0) | 0.35338358 | 0 | 3.53e+299 | 0.0e+00 | info | informational |
| vertex | U = 0: chi_uu = -beta G G delta (Wick), max error | 3.9461359e-15 | 0 | 3.95e-15 | 1.0e-12 | yes | Wick's theorem |
| vertex | Schwinger-Dyson Sigma(i nu_0) from the ED vertex vs ED Sigma | -0.40453059 | -0.40425842 | 2.72e-04 | 2.0e-03 | yes | Rohringer et al. RMP 90, 025003 (2018) |
| vertex | DMFT chi_m(q=0) from Bethe-Salpeter (1/N extrapolated) vs dm/dh | 2.9881371 | 3.1285811 | 4.49e-02 | 5.0e-03 (rel) | **NO** | thermodynamic consistency of DMFT (Baym-Kadanoff; Georges et al. RMP 68, 13 (1996)) |
| DGA | ladder DGA / 2nd-order non-local Sigma(pi,0), U=0.05 (weak-coupling limit) | 1.0657829 | 1 | 6.58e-02 | 1.0e-01 | yes | Katanin, Toschi, Held PRB 80, 075104 (2009): DGA -> 2nd order as U -> 0 |
| DGA | (ratio - 1) scales linearly in U: [r(0.1)-1]/[r(0.05)-1] | 2.2358591 | 2 | 2.36e-01 | 3.0e-01 | yes | ladder resummation is O(U^3) |
| TPSC | U=4t: local spin sum rule | 0.7543534 | 0.7543534 | 1.15e-10 | 1.0e-08 | yes | Vilk & Tremblay J. Phys. I 7, 1309 (1997) |
| TPSC | U=4t: local charge sum rule | 0.2456466 | 0.2456466 | 1.34e-13 | 1.0e-08 | yes | Vilk & Tremblay (1997) |
| TPSC | U -> 0: double occupancy = n^2/4 (Pauli + bosonic tail) | 0.24997681 | 0.25 | 2.32e-05 | 2.0e-04 | yes | exact |
| TPSC | Tr(Sigma G0) - U<n_up n_dn>, extrapolated (1/N) to infinite cutoff | -0.00013312586 | 0 | 1.33e-04 | 1.0e-03 | yes | exact TPSC identity (Vilk & Tremblay 1997) |
| Gauss skew | Dirac mu=0.7: semiclassical w(4a) skew vs Sinitsyn | -0.040136762 | -0.040136762 | 5.55e-17 | 1.0e-06 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. (72) |
| Gauss skew | Dirac mu=1.0: semiclassical w(4a) skew vs Sinitsyn | -0.1377551 | -0.1377551 | 4.72e-16 | 1.0e-06 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. (72) |
| Gauss skew | Dirac mu=2.0: semiclassical w(4a) skew vs Sinitsyn | -0.23372576 | -0.23372576 | 9.71e-16 | 1.0e-06 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. (72) |
| Gauss skew | quadratic-mass Dirac: semiclassical w(4a) vs Kubo-Streda residual | -0.16085313 | -0.16086609 | 8.05e-05 | 2.0e-03 (rel) | yes | two independent routes |
| range | finite-range impurities d=0.5: sigma_xx vs Bessel formula | 108.47742 | 108.47742 | 2.62e-16 | 1.0e-08 (rel) | yes | analytic (this work) |
| range | finite-range impurities d=2.0: sigma_xx vs Bessel formula | 1037.4377 | 1037.4377 | 1.53e-15 | 1.0e-08 (rel) | yes | analytic (this work) |
| crossed | oblique-cell Wannier model vs square-lattice code: X | 0.70451144 | 0.70450994 | 2.13e-06 | 1.0e-05 (rel) | yes | same model, different unit cell |
| crossed | 3D stacked layers (x c): X vs 2D | 0.70451009 | 0.70450994 | 2.15e-07 | 1.0e-05 (rel) | yes | layer decoupling |
| continuation | SAC L1 error (two-Gaussian, sigma=1e-4) | 0.074711323 | 0 | 7.47e-02 | 1.5e-01 | yes | Shao & Sandvik, Phys. Rep. 1003, 1 (2023) |
| continuation | MaxEnt L1 error (same data, for comparison) | 0.25302171 | 0 | 2.53e+299 | 0.0e+00 | info | informational |
| continuation | Nevanlinna + Hardy optimisation, 20 points: L1 error | 0.048503631 | 0 | 4.85e-02 | 1.0e-01 | yes | Fei, Yeh, Gull PRL 126, 056402 (2021) |
| continuation | Nevanlinna + Hardy optimisation, 30 points: L1 error | 0.041599875 | 0 | 4.16e-02 | 1.0e-01 | yes | Fei, Yeh, Gull PRL 126, 056402 (2021) |
| continuation | Nevanlinna + Hardy optimisation, 40 points: L1 error | 0.062869272 | 0 | 6.29e-02 | 1.0e-01 | yes | Fei, Yeh, Gull PRL 126, 056402 (2021) |
| polaron | Froehlich alpha=0.2: m*/m (momentum-fluctuation estimator) | 1.033997 | 1.0342784 | 2.81e-04 | 2.0e-03 | yes | Roseler 1968; Smondyrev 1986 |
| polaron | Froehlich alpha=0.5: m*/m (momentum-fluctuation estimator) | 1.0888688 | 1.0892402 | 3.71e-04 | 5.0e-03 | yes | Roseler 1968; Smondyrev 1986 |
| lattice pol. | 1D Holstein t=w=g=1: E0 (VED) | -2.4696847 | -2.4696847 | 4.35e-14 | 1.0e-10 | yes | Bonca, Trugman, Batistic PRB 60, 1633 (1999); Ku, Trugman, Bonca PRB 65, 174306 (2002) |
| lattice pol. | 1D Holstein g=0.1: E0 vs second-order perturbation theory | -2.0044741 | -2.0044721 | 1.99e-06 | 5.0e-06 | yes | Rayleigh-Schroedinger |
| lattice pol. | 1D SSH g=0.1: E0 vs second-order perturbation theory | -2.0076691 | -2.0076393 | 2.97e-05 | 1.0e-04 | yes | Marchand et al. PRL 105, 266605 (2010) model |
| lattice pol. | 1D Holstein t=w=g=1: effective mass m*/m | 1.3504898 | 0 | 1.35e+300 | 0.0e+00 | info | informational |
