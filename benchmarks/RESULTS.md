# feynlab benchmark results

102/102 checks passed, plus 2 informational rows  (runtime 1221 s, full mode)

IPT coexistence (beta=100): U_c1 ~ 2.6D, U_c2 ~ 3.4D; Lanczos ED-DMFT (6/7 bath sites, T=0): U_c1 ~ 2.38D, U_c2 ~ 2.96D (consensus 2.39D / 2.92-3.0D).

| group | check | feynlab | reference | error | tol | pass | source |
|---|---|---|---|---|---|---|---|
| matsubara | tau->iw transform, max |error| (free level) | 5.42811e-06 | 0 | 5.43e-06 | 1.0e-04 | yes | exact 1/(iw - e) |
| matsubara | iw->tau transform, max |error| (free level) | 8.74496e-09 | 0 | 8.74e-09 | 1.0e-06 | yes | exact -e^{-e tau}/(1+e^{-beta e}) |
| diagrams | # connected G diagrams, order 1 | 2 | 2 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022; Molinari PRB 71, 113102 (2005) |
| diagrams | # proper self-energy diagrams, order 1 | 2 | 2 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # improper diagrams, order 1 | 0 | 0 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # connected G diagrams, order 2 | 10 | 10 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022; Molinari PRB 71, 113102 (2005) |
| diagrams | # proper self-energy diagrams, order 2 | 6 | 6 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # improper diagrams, order 2 | 4 | 4 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # connected G diagrams, order 3 | 74 | 74 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022; Molinari PRB 71, 113102 (2005) |
| diagrams | # proper self-energy diagrams, order 3 | 42 | 42 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # improper diagrams, order 3 | 32 | 32 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # Sigma diagrams (Hartree G), order 2, 0 loops | 2 | 2 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 2, 1 loops | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 3, 0 loops | 10 | 10 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 3, 1 loops | 9 | 9 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 3, 2 loops | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 4, 0 loops | 74 | 74 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 4, 1 loops | 91 | 91 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 4, 2 loops | 23 | 23 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # Sigma diagrams (Hartree G), order 4, 3 loops | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | Molinari PRB 71, 113102 (2005), Eq. (11) |
| diagrams | # connected G diagrams, order 4 | 706 | 706 | 0.00e+00 | 0.0e+00 | yes | Mohr, physics/0512022 |
| diagrams | # Hubbard Sigma diagrams (HF-shifted G0), order 2 | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | standard (e.g. Yamada 1975: 1, 2, 12) |
| diagrams | # Hubbard Sigma diagrams (HF-shifted G0), order 3 | 2 | 2 | 0.00e+00 | 0.0e+00 | yes | standard (e.g. Yamada 1975: 1, 2, 12) |
| diagrams | # Hubbard Sigma diagrams (HF-shifted G0), order 4 | 12 | 12 | 0.00e+00 | 0.0e+00 | yes | standard (e.g. Yamada 1975: 1, 2, 12) |
| diagMC | Hubbard atom: sum of 3rd-order Sigma diagrams (Im, iw0) | -0.00296869 | 0 | 2.97e-03 | 1.5e-02 | yes | exact atom: Sigma = U/2 + U^2/(4iw) (orders>2 cancel) |
| diagMC | Hubbard atom: sum of 4th-order Sigma diagrams (Im, iw0) | 0.0274777 | 0 | 2.75e-02 | 1.1e-01 | yes | exact atom |
| diagMC | Hubbard atom: 2nd-order Sigma(iw0) Im | -0.795837 | -0.795775 | 6.23e-05 | 4.3e-03 | yes | exact U^2/(4 i w0) |
| diagMC | Anderson U=1: Sigma series to 4th order, Im Sigma(iw0) | -0.0173716 | -0.0173105 | 6.11e-05 | 4.3e-03 | yes | exact diagonalization (this package) |
| diagMC | Anderson U=1: Sigma series to 4th order, Im Sigma(iw1) | -0.0407375 | -0.0405785 | 1.59e-04 | 4.3e-03 | yes | exact diagonalization (this package) |
| diagMC | Anderson U=1: Sigma series to 4th order, Im Sigma(iw2) | -0.0506636 | -0.0504655 | 1.98e-04 | 4.3e-03 | yes | exact diagonalization (this package) |
| CT-INT | Anderson U=3 (no p-h symm.): Im Sigma(iw0) vs ED | -0.385203 | -0.386697 | 1.49e-03 | 1.8e-02 | yes | exact diagonalization |
| CT-INT | Anderson U=3 (no p-h symm.): Im Sigma(iw1) vs ED | -0.604451 | -0.604596 | 1.45e-04 | 1.8e-02 | yes | exact diagonalization |
| CT-INT | Anderson U=3 (no p-h symm.): Im Sigma(iw2) vs ED | -0.621926 | -0.620938 | 9.88e-04 | 2.0e-02 | yes | exact diagonalization |
| CT-INT | Anderson U=3 (no p-h symm.): Im Sigma(iw3) vs ED | -0.587355 | -0.584054 | 3.30e-03 | 2.2e-02 | yes | exact diagonalization |
| CT-INT | average sign (repulsive U, delta=0.51) | 1 | 1 | 0.00e+00 | 1.0e-09 | yes | Gull et al. RMP 83, 349 (2011) |
| weak | Yamada-Yosida: 1/Z - 1 from 2nd-order Sigma, U/Gamma=0.5 | 0.0136594 | 0.0134909 | 1.25e-02 | 3.0e-02 (rel) | yes | Yamada, Prog. Theor. Phys. 53, 970 (1975) |
| weak | Yamada-Yosida: 1/Z - 1 from 2nd-order Sigma, U/Gamma=1.0 | 0.0546377 | 0.0539636 | 1.25e-02 | 3.0e-02 (rel) | yes | Yamada, Prog. Theor. Phys. 53, 970 (1975) |
| weak | 2D Hubbard 2nd-order Sigma tail coefficient U^2 n(1-n) | 0.250026 | 0.25 | 2.62e-05 | 1.0e-03 | yes | exact moment sum rule |
| weak | spin-fluctuation (RPA/FLEX) -> 2nd order as U->0 (rel. diff, U=0.05) | 0.0139133 | 0 | 1.39e-02 | 2.0e-02 | yes | diagram expansion O(U^3) |
| DMFT | IPT Uc2 (metal disappears), beta=100 [IPT approximation; not a pass/fail check] | 3.4 | 2.94 | 1.56e-01 | 0.0e+00 | info | exact 2.92-3.0D; IPT is known to overestimate U_c2 -> use ED/CT-INT |
| DMFT | IPT Uc1 (insulator disappears), beta=100 [informational] | 2.6 | 2.39 | 8.79e-02 | 0.0e+00 | info | exact 2.39(2)D |
| DMFT | U=D, beta=10: Z(CT-INT DMFT) vs Z(IPT) | 0.742933 | 0.750106 | 7.17e-03 | 3.0e-02 | yes | IPT exact to O(U^2); agreement expected at weak coupling |
| DMFT | U=6D Hubbard-I: Sigma(iw0) -> atomic U^2/(4 iw0) | -143.239 | -143.239 | 2.33e-10 | 1.0e-02 (rel) | yes | atomic limit (exact at t=0) |
| e-gas | RPA plasmon dispersion coefficient alpha in w^2=wp^2+alpha vF^2 q^2 | 0.601676 | 0.6 | 1.68e-03 | 1.0e-02 | yes | Mattuck Ex. 13.6 (3/5) |
| e-gas | RPA correlation energy slope d e_c / d ln rs (rs->0) | 0.0310808 | 0.0310907 | 3.17e-04 | 5.0e-04 (rel) | yes | Gell-Mann & Brueckner, Phys. Rev. 106, 364 (1957) |
| e-gas | RPA correlation energy constant (Ha) | -0.0711254 | -0.0711 | 2.54e-05 | 5.0e-04 | yes | RPA high-density constant -0.0711 Ha |
| BCS | Delta(0)/Tc, lambda=0.2 | 1.76396 | 1.76388 | 4.54e-05 | 2.0e-03 (rel) | yes | BCS weak coupling pi e^-gamma = 1.7639 |
| BCS | Delta(0)/Tc, lambda=0.25 | 1.76447 | 1.76388 | 3.36e-04 | 2.0e-03 (rel) | yes | BCS weak coupling pi e^-gamma = 1.7639 |
| BCS | Tc / (w_D e^{-1/lambda}), lambda=0.2 | 1.13387 | 1.13387 | 1.82e-14 | 1.0e-03 (rel) | yes | BCS 2e^gamma/pi = 1.134 |
| Berry | QWZ m=-2.5: |sigma_xy| (e^2/h) at mu=0 | 4.08735e-17 | 0 | 4.09e-17 | 1.0e-03 | yes | Qi, Wu, Zhang PRB 74, 085308 (2006): |C|=1 for 0<|m|<2, 0 otherwise |
| Berry | QWZ m=-2.5: Kubo sigma_xy = -C(lower band, FHS lattice method) | -4.08735e-17 | 8.83487e-18 | 4.97e-17 | 1.0e-03 | yes | TKNN; Fukui-Hatsugai-Suzuki JPSJ 74, 1674 (2005) |
| Berry | QWZ m=-1.0: |sigma_xy| (e^2/h) at mu=0 | 1 | 1 | 3.33e-16 | 1.0e-03 | yes | Qi, Wu, Zhang PRB 74, 085308 (2006): |C|=1 for 0<|m|<2, 0 otherwise |
| Berry | QWZ m=-1.0: Kubo sigma_xy = -C(lower band, FHS lattice method) | -1 | -1 | 3.33e-16 | 1.0e-03 | yes | TKNN; Fukui-Hatsugai-Suzuki JPSJ 74, 1674 (2005) |
| Berry | QWZ m=1.0: |sigma_xy| (e^2/h) at mu=0 | 1 | 1 | 1.11e-16 | 1.0e-03 | yes | Qi, Wu, Zhang PRB 74, 085308 (2006): |C|=1 for 0<|m|<2, 0 otherwise |
| Berry | QWZ m=1.0: Kubo sigma_xy = -C(lower band, FHS lattice method) | 1 | 1 | 1.11e-16 | 1.0e-03 | yes | TKNN; Fukui-Hatsugai-Suzuki JPSJ 74, 1674 (2005) |
| Berry | QWZ m=2.5: |sigma_xy| (e^2/h) at mu=0 | 4.08735e-17 | 0 | 4.09e-17 | 1.0e-03 | yes | Qi, Wu, Zhang PRB 74, 085308 (2006): |C|=1 for 0<|m|<2, 0 otherwise |
| Berry | QWZ m=2.5: Kubo sigma_xy = -C(lower band, FHS lattice method) | -4.08735e-17 | -5.96354e-17 | 1.88e-17 | 1.0e-03 | yes | TKNN; Fukui-Hatsugai-Suzuki JPSJ 74, 1674 (2005) |
| Berry | Haldane (t2=0.1,phi=pi/2,M=0): |sigma_xy| (e^2/h) | 1 | 1 | 1.11e-16 | 1.0e-03 | yes | Haldane PRL 61, 2015 (1988) |
| Berry | Haldane trivial (M=1): lowest band Chern number | 7.0679e-17 | 0 | 7.07e-17 | 1.0e-06 | yes | Haldane PRL 61, 2015 (1988) |
| Berry | Rashba intrinsic SHC (e/2pi) mu=0.5 | 0.250001 | 0.25 | 1.15e-06 | 2.0e-03 | yes | Sinova et al. PRL 92, 126603 (2004): e/8pi |
| Berry | Rashba intrinsic SHC (e/2pi) mu=1.0 | 0.249998 | 0.25 | 1.78e-06 | 2.0e-03 | yes | Sinova et al. PRL 92, 126603 (2004): e/8pi |
| AHE | Dirac mu=0.7: intrinsic (e^2/h) | -0.357143 | -0.357143 | 4.51e-09 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 58 |
| AHE | Dirac mu=0.7: side jump | -0.138249 | -0.138249 | 2.78e-17 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 65 |
| AHE | Dirac mu=0.7: anomalous distribution | -0.138249 | -0.138249 | 0.00e+00 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 67 |
| AHE | Dirac mu=0.7: Gaussian (intrinsic) skew | -0.0401403 | -0.0401368 | 3.57e-06 | 3.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 72, 1st term |
| AHE | Dirac mu=0.7: conventional skew (V1/V0=0.3) | -0.00715201 | -0.00715201 | 2.60e-18 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 72, 2nd term |
| AHE | Dirac mu=0.7: total, Gaussian disorder (NCA) | -0.673781 | -0.673777 | 3.57e-06 | 3.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. B8 |
| AHE | Dirac mu=0.7: sigma_xx Boltzmann vs Drude (e^2/h) | 77.4194 | 77.4194 | 3.67e-16 | 1.0e-03 (rel) | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 54 (tau_tr) |
| AHE | Dirac mu=1.0: intrinsic (e^2/h) | -0.25 | -0.25 | 1.31e-08 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 58 |
| AHE | Dirac mu=1.0: side jump | -0.214286 | -0.214286 | 2.78e-17 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 65 |
| AHE | Dirac mu=1.0: anomalous distribution | -0.214286 | -0.214286 | 2.78e-17 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 67 |
| AHE | Dirac mu=1.0: Gaussian (intrinsic) skew | -0.137754 | -0.137755 | 1.48e-06 | 3.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 72, 1st term |
| AHE | Dirac mu=1.0: conventional skew (V1/V0=0.3) | -0.0350667 | -0.0350667 | 0.00e+00 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 72, 2nd term |
| AHE | Dirac mu=1.0: total, Gaussian disorder (NCA) | -0.816325 | -0.816327 | 1.49e-06 | 3.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. B8 |
| AHE | Dirac mu=1.0: sigma_xx Boltzmann vs Drude (e^2/h) | 171.429 | 171.429 | 3.32e-16 | 1.0e-03 (rel) | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 54 (tau_tr) |
| AHE | Dirac mu=2.0: intrinsic (e^2/h) | -0.125 | -0.125 | 2.81e-09 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 58 |
| AHE | Dirac mu=2.0: side jump | -0.197368 | -0.197368 | 8.33e-17 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 65 |
| AHE | Dirac mu=2.0: anomalous distribution | -0.197368 | -0.197368 | 8.33e-17 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 67 |
| AHE | Dirac mu=2.0: Gaussian (intrinsic) skew | -0.233716 | -0.233726 | 9.83e-06 | 3.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 72, 1st term |
| AHE | Dirac mu=2.0: conventional skew (V1/V0=0.3) | -0.118994 | -0.118994 | 5.55e-17 | 2.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 72, 2nd term |
| AHE | Dirac mu=2.0: total, Gaussian disorder (NCA) | -0.753453 | -0.753463 | 9.83e-06 | 3.0e-03 | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. B8 |
| AHE | Dirac mu=2.0: sigma_xx Boltzmann vs Drude (e^2/h) | 315.789 | 315.789 | 0.00e+00 | 1.0e-03 (rel) | yes | Sinitsyn et al. PRB 75, 045315 (2007) Eq. 54 (tau_tr) |
| SHE | Rashba SHC with ladder vertex, total (e/2pi), mu=0.5 | -0.000114346 | 0 | 1.14e-04 | 1.5e-02 | yes | Inoue, Bauer, Molenkamp PRB 70, 041303(R) (2004): exact cancellation |
| SHE | Rashba SHC with ladder vertex, total (e/2pi), mu=1.0 | -0.0096813 | 0 | 9.68e-03 | 1.5e-02 | yes | Inoue, Bauer, Molenkamp PRB 70, 041303(R) (2004): exact cancellation |
| AHE-X | Dirac eps=1: X diagram (e^2/h) | 0.485129 | 0.489796 | 9.53e-03 | 1.5e-02 (rel) | yes | Ado, Dmitriev, Ostrovsky, Titov, EPL 111, 37004 (2015) Eq. 17 |
| AHE-X | Dirac eps=1: Psi diagram (e^2/h) | 0.00319609 | 0 | 3.20e-03 | 1.0e-02 | yes | Ado, Dmitriev, Ostrovsky, Titov, EPL 111, 37004 (2015) Eq. 17 (Psi = 0) |
| AHE-X | Dirac eps=1: NCA total vs Ado Eq. 3 | -0.816325 | -0.816327 | 1.49e-06 | 3.0e-03 | yes | Ado, Dmitriev, Ostrovsky, Titov, EPL 111, 37004 (2015) Eq. 3 |
| AHE-X | Dirac eps=1: full sigma_xy (NCA + X + Psi) | -0.328 | -0.326531 | 1.47e-03 | 1.0e-02 | yes | Ado, Dmitriev, Ostrovsky, Titov, EPL 111, 37004 (2015) Eq. 4 |
| polaron | Feynman variational E0, alpha=1 | -1.01303 | -1.01303 | 3.54e-08 | 2.0e-06 | yes | Feynman, Phys. Rev. 97, 660 (1955); values from Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I |
| polaron | Feynman variational E0, alpha=3 | -3.13333 | -3.13333 | 4.43e-08 | 2.0e-06 | yes | Feynman, Phys. Rev. 97, 660 (1955); values from Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I |
| polaron | Feynman variational E0, alpha=5 | -5.44014 | -5.44014 | 5.79e-10 | 2.0e-06 | yes | Feynman, Phys. Rev. 97, 660 (1955); values from Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I |
| polaron | DiagMC E0 (bare series, all orders), alpha=0.5  [+-0.0007] | -0.504738 | -0.504081 | 6.58e-04 | 2.2e-03 | yes | exact weak-coupling series (Smondyrev 1986) |
| polaron | DiagMC E0 (bare series, all orders), alpha=1.0  [+-0.0023] | -1.01742 | -1.01662 | 8.01e-04 | 6.9e-03 | yes | Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I |
| polaron | DiagMC E0 (bare series, all orders), alpha=2.0  [+-0.0042] | -2.05926 | -2.06957 | 1.03e-02 | 1.3e-02 | yes | Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I |
| polaron | DiagMC E0 (bare series, all orders), alpha=3.0  [+-0.0081] | -3.17241 | -3.16829 | 4.12e-03 | 2.5e-02 | yes | Hahn, Klimin, Tempere, Devreese, Franchini, PRB 97, 134305 (2018), Table I |
| bold | atom betaU=0.5: skeleton series(G_exact)/Sigma_exact | 0.999314 | 1 | 6.86e-04 | 5.0e-03 | yes | Kozik, Ferrero, Georges, PRL 114, 156402 (2015): physical branch at weak coupling |
| bold | atom betaU=6.0: skeleton series(G_exact) misses Sigma_exact (1-ratio > 0.5) | 1 | 1 | 0.00e+00 | 0.0e+00 | yes | Kozik, Ferrero, Georges, PRL 114, 156402 (2015): convergence to unphysical branch at strong coupling |
| ED | sector-resolved ED vs dense ED, max|dG| | 8.1251e-13 | 0 | 8.13e-13 | 1.0e-10 | yes | exact |
| DMFT-ED | Lanczos ED-DMFT (7 bath = odd, T=0) U_c2, linear Z->0 (window 0.04<=Z<=0.15) | 2.95785 | 2.96 | 2.15e-03 | 4.0e-02 | yes | consensus U_c2 = 2.92-3.0D, U_c1 = 2.39(2)D: Karski, Raas, Uhrig PRB 72, 113110 (2005) and refs; Bulla PRL 83, 136 (1999) |
| DMFT-ED | Lanczos ED-DMFT (6 bath = even, T=0) U_c1 (insulator lost) | 2.375 | 2.39 | 1.50e-02 | 5.0e-02 | yes | consensus U_c2 = 2.92-3.0D, U_c1 = 2.39(2)D: Karski, Raas, Uhrig PRB 72, 113110 (2005) and refs; Bulla PRL 83, 136 (1999) |
| DMFT-ED | T=0 Lanczos ED vs full sector ED (4 bath), max|dG| | 1.11163e-10 | 0 | 1.11e-10 | 1.0e-08 | yes | exact |
| DMFT | LatticeDMFT(semicircle DOS) Z == BetheDMFT Z (U=2.5) | 0.189316 | 0.18932 | 3.86e-06 | 1.0e-04 | yes | identity |
| MaxEnt | two-peak spectrum: chi^2 per point | 0.960672 | 1 | 3.93e-02 | 2.0e-01 | yes | Jarrell & Gubernatis Phys. Rep. 269, 133 (1996) |
| MaxEnt | two-peak spectrum: L1 error of A(w) | 0.210052 | 0 | 2.10e-01 | 3.0e-01 | yes | synthetic test (noise 1e-4) |

Figures: figures/ahe_dirac_mechanisms.png, figures/ipt_mott_hysteresis.png, figures/ctint_vs_ed.png