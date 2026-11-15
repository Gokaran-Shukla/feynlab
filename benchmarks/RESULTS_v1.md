# feynlab v1.0 benchmark results (ROADMAP modules)

54/54 checks passed, plus 1 informational rows (runtime 649 s)

| group | check | feynlab | reference | error | tol | pass | source |
|---|---|---|---|---|---|---|---|
| errors | Gamma method tau_int, AR(1) rho=0.9 | 9.02284 | 9.5 | 5.02e-02 | 8.0e-02 (rel) | yes | Wolff, CPC 156, 143 (2004); exact (1+rho)/(2(1-rho)) |
| errors | jackknife <a>/<b> (a~1, b~2) | 0.500343 | 0.5 | 3.43e-04 | 1.7e-03 | yes | exact ratio of means |
| mixing | IPT U=2.6 beta=100: Z(Anderson) - Z(linear) | -2.15531e-10 | 0 | 2.16e-10 | 1.0e-06 | yes | same fixed point (Anderson, J. ACM 12, 547 (1965)) |
| mixing | IPT U=2.6: iterations linear / Anderson | 6 | 1 | 5.00e+00 | 0.0e+00 | info | speed-up factor (informational) |
| CT-HYB | U=3.0: occupation per spin vs ED | 0.47125 | 0.471163 | 8.72e-05 | 1.0e-02 | yes | Werner et al. PRL 97, 076405 (2006) vs exact diagonalization |
| CT-HYB | U=3.0: double occupancy vs ED | 0.0827708 | 0.0823542 | 4.17e-04 | 6.0e-03 | yes | Werner et al. PRL 97, 076405 (2006) vs exact diagonalization |
| CT-HYB | U=3.0: Im Sigma(i w_0) vs ED | -0.383161 | -0.386697 | 3.54e-03 | 3.0e-02 | yes | exact diagonalization |
| CT-HYB | U=4.0: occupation per spin vs ED | 0.500663 | 0.500483 | 1.80e-04 | 1.0e-02 | yes | Werner et al. PRL 97, 076405 (2006) vs exact diagonalization |
| CT-HYB | U=4.0: double occupancy vs ED | 0.0780475 | 0.0781527 | 1.05e-04 | 6.0e-03 | yes | Werner et al. PRL 97, 076405 (2006) vs exact diagonalization |
| CT-HYB | U=4.0: Im Sigma(i w_0) vs ED | -0.767218 | -0.76621 | 1.01e-03 | 4.2e-02 | yes | exact diagonalization |
| Kanamori | 2 orbitals, N=2: lowest multiplet (3-fold) = U - 3J | 2.2 | 2.2 | 8.88e-16 | 1.0e-10 | yes | Georges, de' Medici, Mravlje, Annu. Rev. CMP 4, 137 (2013) |
| Kanamori | 2 orbitals, N=2: highest multiplet = U + J | 4.6 | 4.6 | 0.00e+00 | 1.0e-10 | yes | Georges et al. (2013) |
| Kanamori | t2g^3 ground state 4S energy = 3U - 9J | 6.6 | 6.6 | 3.55e-15 | 1.0e-10 | yes | Georges et al. (2013) |
| Kanamori | t2g^3 ground-state degeneracy (S=3/2) | 4 | 4 | 0.00e+00 | 0.0e+00 | yes | Hund's rule |
| Kanamori | MultiOrbitalED one-orbital limit, max|G - G_ED| | 1.99437e-14 | 0 | 1.99e-14 | 1.0e-10 | yes | single-orbital AndersonED |
| DFT+DMFT | hr.dat write/read round trip, sum|H_R| | 1.7 | 1.7 | 0.00e+00 | 1.0e-05 | yes | Wannier90 format |
| DFT+DMFT | U=0: total occupation (target 2) | 2 | 2 | 3.10e-07 | 2.0e-03 | yes | electron count |
| DFT+DMFT | U=0: Z of orbital 1 | 1 | 1 | 0.00e+00 | 1.0e-03 | yes | non-interacting |
| DFT+DMFT | U=4 J=0.5 half filling, Hubbard-I: max Z (Mott) | 0.00170937 | 0 | 1.71e-03 | 2.0e-02 | yes | Mott insulator (Lichtenstein & Katsnelson PRB 57, 6884 (1998)) |
| DCA | N_c = 1 DCA = single-site DMFT (Z) | 0.486644 | 0.486638 | 5.91e-06 | 1.0e-04 | yes | Hettler et al. PRB 61, 12739 (2000) |
| DCA | N_c = 4: Z(pi,0) < Z(0,0) (antinodal suppression) | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | Maier et al. RMP 77, 1027 (2005) |
| two-part. | chi0(q=0, Omega=0) = dn/dmu (compressibility sum rule) | 0.208784 | 0.208784 | 1.07e-08 | 1.0e-06 | yes | exact identity |
| two-part. | RPA chi_s(Q) = chi0/(1 - U chi0) at Q=(pi,pi) | 0.710319 | 0.710319 | 1.11e-16 | 1.0e-12 | yes | RPA (Berk-Schrieffer) |
| NRG | U=0 resonant level: S_imp(T ~ Gamma) / exact | 1.01761 | 1 | 1.76e-02 | 6.0e-02 | yes | exact resonant-level entropy |
| NRG | U=0 resonant level: S_imp(T << Gamma) / exact | 1.03647 | 1 | 3.65e-02 | 6.0e-02 | yes | exact resonant-level entropy |
| NRG | Kondo U=10 Gamma: S_imp(T -> 0) | 1.19829e-08 | 0 | 1.20e-08 | 3.0e-02 | yes | Wilson, RMP 47, 773 (1975); Krishna-murthy et al. PRB 21, 1003 |
| NRG | Kondo U=10 Gamma: local-moment plateau S_imp = ln 2 | 0.687711 | 0.693147 | 5.44e-03 | 6.0e-02 | yes | Krishna-murthy, Wilkins, Wilson (1980) |
| NRG | Haldane slope d ln T_K / d[pi U/8G - pi G/2U] | -0.992686 | -1 | 7.31e-03 | 5.0e-02 | yes | Haldane, PRL 40, 416 (1978) |
| continuation | Nevanlinna interpolant reproduces input, max|dG| | 1.11022e-16 | 0 | 1.11e-16 | 1.0e-10 | yes | Fei, Yeh, Gull PRL 126, 056402 (2021) |
| continuation | Nevanlinna G(2i) vs exact (analytic continuation off axis) | 5.67181e-17 | 0 | 5.67e-17 | 1.0e-04 | yes | Fei, Yeh, Gull (2021) |
| continuation | Nevanlinna spectral normalisation | 0.995504 | 1 | 4.50e-03 | 3.0e-02 | yes | sum rule |
| continuation | stochastic continuation: normalisation | 1 | 1 | 4.21e-12 | 2.0e-02 | yes | Mishchenko et al. PRB 62, 6317 (2000); Sandvik PRB 57, 10287 (1998) |
| CDet | Hubbard atom a_1 (double occupancy series) vs exact | -0.3125 | -0.3125 | 0.00e+00 | 1.0e-06 | yes | Rossi PRL 119, 045701 (2017); exact -(1/4) tanh(beta U/4) |
| CDet | Hubbard atom a_3 vs exact | 0.16276 | 0.16276 | 0.00e+00 | 1.0e-06 | yes | exact |
| CDet | Anderson model a_1 vs ED Taylor fit | -0.0681258 | -0.067919 | 2.07e-04 | 1.0e-03 | yes | exact diagonalization |
| resum | Domb-Sykes radius of (1+x)^(-1/2) | 1 | 1 | 6.66e-16 | 1.0e-06 | yes | Domb & Sykes Proc. R. Soc. A 240, 214 (1957) |
| resum | Pade[8/8] beyond radius, x = 3 | 0.5 | 0.5 | 7.74e-09 | 1.0e-06 | yes | Baker & Graves-Morris (1996) |
| resum | conformal-map resummation, x = 3 | 0.5 | 0.5 | 4.30e-10 | 1.0e-06 | yes | Rossi, Ohgoe, Van Houcke, Werner PRL 121, 130405 (2018) |
| polaron | Froehlich alpha=0.2: Z vs 1 - alpha/2 | 0.899986 | 0.9 | 1.37e-05 | 3.5e-02 | yes | Mishchenko, Prokof'ev, Sakamoto, Svistunov PRB 62, 6317 (2000) |
| polaron | Froehlich alpha=0.2: m*/m vs 1 + alpha/6 + 0.0236 alpha^2 | 1.08143 | 1.03428 | 4.72e-02 | 6.0e-02 | yes | Roseler (1968); Smondyrev, Theor. Math. Phys. 68, 653 (1986) |
| T-matrix | T-matrix sigma_I = Born + third-order skew (n_i V0^2 fixed) | -1.72649 | -1.72625 | 1.40e-04 | 1.0e-02 (rel) | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. (72) |
| T-matrix | strong scatterers V0=+-0.5: even part ~ Born | -0.804424 | -0.813628 | 1.13e-02 | 3.0e-02 (rel) | yes | Born limit |
| T-matrix | identity impurity matrix reproduces scalar disorder | -0.813628 | -0.813628 | 0.00e+00 | 1.0e-12 | yes | identity |
| crossed | lattice (BZ-FFT) QWZ X vs continuum k.p limit | 0.794698 | 0.790464 | 5.36e-03 | 2.0e-02 (rel) | yes | Ado et al. EPL 111, 37004 (2015), k.p limit of QWZ |
| crossed | lattice QWZ Psi vs continuum k.p limit | -0.499375 | -0.524061 | 4.71e-02 | 8.0e-02 (rel) | yes | Ado et al. (2015); Psi != 0 from the k^2 mass |
| OAM | Dirac valley orbital moment m(k), max deviation | 2.22045e-16 | 0 | 2.22e-16 | 1.0e-10 | yes | Xiao, Yao, Niu PRL 99, 236809 (2007) |
| BCD | tilted Dirac: D_x Fermi-surface = Fermi-sea | -0.034916 | -0.034916 | 2.10e-08 | 1.0e-06 (rel) | yes | Sodemann & Fu PRL 115, 216806 (2015) |
| BCD | untilted Dirac: D_x = 0 (C3/mirror) | 2.35185e-18 | 0 | 2.35e-18 | 1.0e-10 | yes | Sodemann & Fu (2015) |
| Edelstein | Rashba |delta S_y| per eE tau = m alpha / 4 pi | 0.0238741 | 0.0238732 | 3.69e-05 | 5.0e-03 (rel) | yes | Edelstein, Solid State Commun. 73, 233 (1990) |
| Edelstein | Kubo intraband (odd) / tau = Boltzmann | -0.0238182 | -0.0238741 | 2.34e-03 | 1.0e-02 (rel) | yes | Freimuth, Bluegel, Mokrousov PRB 90, 174423 (2014) |
| Edelstein | Rashba, no magnetisation: intrinsic (even) S_y = 0 | -1.61146e-19 | 0 | 1.61e-19 | 1.0e-10 | yes | time-reversal symmetry |
| finite-T | Dirac sigma_int(mu=1, T=0.05): Fermi-sea vs Bastin convolution | -0.25213 | -0.252131 | 5.53e-07 | 2.0e-05 | yes | Sommerfeld/Bastin; Sinitsyn et al. Eq. (58) at T=0 |
| finite-T | Dirac sigma_int(mu=1, T=0.1): Fermi-sea vs Bastin convolution | -0.258741 | -0.258742 | 9.75e-07 | 2.0e-05 | yes | Sommerfeld/Bastin; Sinitsyn et al. Eq. (58) at T=0 |
| Wannier3D | stacked Chern layers: sigma_xy = C e^2/(h c) | 0.333333 | 0.333333 | 8.02e-08 | 1.0e-04 | yes | Halperin, Jpn. J. Appl. Phys. 26, 1913 (1987) |
| Wannier3D | metal mu=-0.9: adaptive 40^3 vs converged 200x200x160 | 0.210061 | 0.210397 | 1.60e-03 | 5.0e-03 (rel) | yes | Wang, Yates, Souza, Vanderbilt PRB 74, 195118 (2006) |

Figures: figures/nrg_kondo.png, figures/continuation.png
