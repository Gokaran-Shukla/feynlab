# Appendices A–N (intermediate part)
*Companion to R. D. Mattuck, A Guide to Feynman Diagrams in the Many-Body Problem (2nd ed.)*

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


| appendix | content (paraphrased) | where it lives in feynlab |
|---|---|---|
| canonical transformation | finding fictitious particles by changing variables | — (concept) |
| A Dirac formalism | bras, kets, representations | — |
| B time-development operator U(t) | interaction picture, Dyson series | `diagrams.evaluate` (imaginary-time version) |
| C ground-state energy from vacuum amplitude | ΔE = lim i d ln R/dt | Ch. 5 tutorial |
| D, E expansions of U and G | perturbation series of the propagator | `diagrams.enumerate` |
| F Wick's theorem | pairings → determinants | `qmc.ctint` (determinants = all pairings) |
| G graphical expansion | derivation of diagram rules | `Diagram.sign`, Feynman rules in `diagrams.evaluate` |
| H spectral density | Lehmann representation, sum rules | `qmc.ed.AndersonED.green` |
| I the iδ factor | causality, retarded vs time-ordered | `core.lattice.pade_continue(eta=...)` |
| J electron propagator in electron–phonon system | Migdal-type self-energy | — |
| K spin wave functions | spin sums, factor 2 per loop | `diagrams` (generic family) |
| L propagators & analytic properties | retarded/advanced/Matsubara relations | `core.matsubara`, Padé |
| M equation-of-motion decoupling as partial sums | Hubbard-I as a diagram subset | `selfenergy.local.hubbard_one_sigma` |
| N reduced graph expansion | resummed series bookkeeping | `diagrams` skeleton tools |

## Self-check
* Show that ⟨T c c†⟩ for H₀ = Σε c†c reproduces G₀(τ) = −e^{−ετ}(1 − f(ε)) for τ>0 (compare `core.matsubara.g0_tau_level`).
