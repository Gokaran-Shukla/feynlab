# Quizzes — Mattuck chapters and feynlab tutorials

## Ch 0 — The many-body problem for everybody

1. A quasi particle is best described as…
   - A. a bare electron with infinite lifetime
   - B. a real particle dressed by a cloud of disturbance it creates
   - C. a coherent oscillation of all particles
   - D. a hole below the Fermi surface
2. Which is a collective excitation rather than a quasi particle?
   - A. a dressed electron in a metal
   - B. a plasmon
   - C. a polaron
   - D. a quasi hole
3. The quasi-particle concept is useful only if…
   - A. its lifetime is long compared with ħ/ε
   - B. Z = 1
   - C. the interaction is weak everywhere
   - D. T = 0
4. Z < 1 in G ≈ Z/(ω−ε′+i/2τ) means…
   - A. the particle is heavier
   - B. only part of the spectral weight is in the quasi-particle pole
   - C. the lifetime is infinite
   - D. the system is insulating
5. Near the Fermi surface of a normal Fermi liquid 1/τ scales as…
   - A. (ε−ε_F)
   - B. (ε−ε_F)²
   - C. constant
   - D. 1/(ε−ε_F)

<details><summary>Answer key</summary>

1. **B** — Quasi particle = bare particle + polarisation cloud; collective modes are the whole-system oscillations.
2. **B** — Plasmons are collective density oscillations of the whole electron gas.
3. **A** — If it decays before completing an oscillation it is not a well-defined excitation.
4. **B** — The rest (1−Z) is the incoherent background.
5. **B** — Pauli blocking limits phase space quadratically.

</details>

## Ch 1 — Feynman diagrams

1. The single-particle propagator is…
   - A. a probability amplitude to go from (r₁,t₁) to (r₂,t₂)
   - B. an energy
   - C. the density
   - D. a scattering cross section
2. Which propagator gives the ground-state energy?
   - A. two-particle
   - B. particle–hole
   - C. vacuum amplitude (no-particle propagator)
   - D. single-particle
3. Number of topologically distinct connected G-diagrams at 2nd order (two-body interaction):
   - A. 4
   - B. 6
   - C. 10
   - D. 74
4. Collective modes are most naturally found from…
   - A. the vacuum amplitude
   - B. the particle–hole (polarisation) propagator
   - C. the single-particle propagator only
   - D. the density of states
5. Why is the factorial growth of diagrams not fatal?
   - A. diagrams cancel exactly
   - B. infinite subsets can be partially summed in closed form
   - C. high orders are always zero
   - D. computers are fast

<details><summary>Answer key</summary>

1. **A** — Its square gives the probability.
2. **C** — ΔE follows from the time dependence of the vacuum amplitude.
3. **C** — 2, 10, 74, 706 … (feynlab count_diagrams).
4. **B** — Poles of the polarisation propagator = collective modes.
5. **B** — Partial summation (Dyson, ladders, RPA) is the central technique.

</details>

## Ch 2 — Classical pinball propagator

1. Classical path probabilities combine by…
   - A. adding probabilities
   - B. adding amplitudes then squaring
   - C. multiplying probabilities of all paths
   - D. taking the maximum
2. Summing repeated scatterings with weight p at one obstacle gives…
   - A. P₀(1+p)
   - B. P₀/(1−p)
   - C. P₀ eᵖ
   - D. P₀ p
3. The geometric partial sum fails when…
   - A. p < 0
   - B. |p| ≥ 1
   - C. p = 0
   - D. P₀ > 1
4. In Ex. 2.1 with two strong scatterers, the number of order-n diagrams is…
   - A. n
   - B. 2n
   - C. 2ⁿ
   - D. n!
5. The classical partial sum foreshadows…
   - A. Wick's theorem
   - B. Dyson's equation
   - C. the linked-cluster theorem
   - D. BCS

<details><summary>Answer key</summary>

1. **A** — No interference classically.
2. **B** — Geometric series.
3. **B** — Radius of convergence.
4. **C** — Each step picks one of two scatterers.
5. **B** — Renormalised (dressed) propagator.

</details>

## Ch 3 — Quantum pinball propagator

1. The quantum propagator differs from the classical one because it is…
   - A. always real
   - B. an amplitude, so paths interfere
   - C. always smaller
   - D. time independent
2. Summing G₀ + G₀VG₀ + … for one level gives a pole at…
   - A. ε₀
   - B. ε₀ + V
   - C. ε₀ − V
   - D. V
3. Individual terms of the series diverge at ω = ε₀ because…
   - A. V is large
   - B. each G₀ has a pole at the unperturbed energy
   - C. of the iδ
   - D. the series is classical
4. In the random-impurity metal, the first-order (averaged) term gives…
   - A. a lifetime
   - B. a uniform energy shift n W₀
   - C. a mass enhancement only
   - D. zero
5. The second-order impurity term gives a lifetime 1/τ = …
   - A. n W₀
   - B. 2π n W₀² ρ(ε)
   - C. W₀²/ε
   - D. ρ(ε)/n

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — G = 1/(ω − ε₀ − V).
3. **B** — The divergence disappears only after infinite summation.
4. **B** — 
5. **B** — Golden rule / Born approximation; used for Γ in feynlab's Hall engines.

</details>

## Ch 4 — Quasi particles in Fermi systems

1. In diagrams a hole is…
   - A. a dashed line
   - B. a particle line running backward in time
   - C. an interaction
   - D. absent
2. The Hartree term represents…
   - A. exchange
   - B. the average field of the other particles
   - C. a lifetime
   - D. pairing
3. The Fock (exchange) term between opposite spins for a Hubbard contact interaction is…
   - A. U/2
   - B. U n
   - C. zero
   - D. −U n
4. A closed fermion loop contributes a factor…
   - A. +1
   - B. −1 (times spin sum)
   - C. i
   - D. 1/2
5. At high density the dominant diagrams in the electron gas are…
   - A. ladders
   - B. bubble chains (RPA)
   - C. crossed diagrams
   - D. tadpoles

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **C** — Exchange requires equal spins; the contact interaction acts only between opposite spins.
4. **B** — 
5. **B** — 

</details>

## Ch 5 — Vacuum amplitude & ground-state energy

1. The linked-cluster theorem says ΔE is given by…
   - A. all vacuum diagrams
   - B. only connected (linked) vacuum diagrams
   - C. only first-order diagrams
   - D. only ladders
2. Why must disconnected diagrams drop out of ΔE?
   - A. they are zero
   - B. ΔE must be extensive (∝ N)
   - C. they are imaginary
   - D. they violate momentum conservation
3. Topologically equivalent diagrams are…
   - A. counted separately
   - B. counted once
   - C. cancelled
   - D. multiplied by 2
4. ΔE is obtained from the vacuum amplitude R(t) via…
   - A. R(0)
   - B. lim i d ln R/dt
   - C. |R|²
   - D. ∫R dt
5. With momentum conservation, a line beginning and ending on the same interaction carries…
   - A. any momentum
   - B. zero momentum transfer
   - C. 2k_F
   - D. k_F

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **B** — feynlab: Diagram.canonical_key().
4. **B** — 
5. **B** — 

</details>

## Ch 6 — Bird's-eye view

1. Which tool answers 'what is T_c of a superconductor'?
   - A. vacuum amplitude
   - B. Nambu / anomalous propagators
   - C. Hartree–Fock
   - D. pinball sum
2. Finite temperature requires…
   - A. real-time only
   - B. imaginary time and Matsubara frequencies
   - C. no diagrams
   - D. bosons only
3. Self-consistent renormalisation is essential for…
   - A. the free electron gas
   - B. Kondo and critical phenomena
   - C. a single level
   - D. pinball
4. Collective modes come from…
   - A. polarisation propagator poles
   - B. Hartree term
   - C. vacuum amplitude
   - D. Fermi function
5. The kindergarten part of the book ends with chapter…
   - A. 3
   - B. 6
   - C. 9
   - D. 16

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **B** — 
4. **A** — 
5. **B** — 

</details>

## Ch 7 — Second quantization

1. {c_i, c_j†} equals…
   - A. 0
   - B. δ_ij
   - C. 1
   - D. c_i c_j
2. c†_i c†_i acting on any state gives…
   - A. the state
   - B. zero
   - C. twice the state
   - D. −the state
3. The sign of c_j|…⟩ is determined by…
   - A. the number of occupied states to the left of j
   - B. the energy
   - C. spin
   - D. nothing
4. A two-body operator in second quantization has the form…
   - A. Σ V c†c
   - B. ½ Σ V c†c†cc
   - C. Σ V cc
   - D. V N
5. Spin operator S = ½ Σ c†_α σ_{αβ} c_β uses…
   - A. Pauli matrices
   - B. Gell-Mann matrices
   - C. Dirac matrices
   - D. unit matrices

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — Pauli principle.
3. **A** — 
4. **B** — 
5. **A** — 

</details>

## Ch 8 — More about quasi particles

1. The general propagator contains a quasi-particle pole plus…
   - A. nothing
   - B. an incoherent background
   - C. a delta function at zero
   - D. a plasmon always
2. Quasi-particle lifetimes diverge at the Fermi surface because…
   - A. Pauli blocking of decay channels
   - B. interactions vanish
   - C. Z = 0
   - D. of impurities
3. In DMFT/IPT, Z decreases with U and vanishes at…
   - A. U = 0
   - B. the Mott transition U_c2
   - C. U = W/10
   - D. never
4. Fermi-liquid behaviour in Matsubara space appears as…
   - A. Im Σ(iωₙ) ∝ ωₙ at low frequency
   - B. Im Σ = const
   - C. Re Σ = 0
   - D. Σ ∝ 1/iωₙ
5. The pure Hartree model is useful because…
   - A. it is exactly soluble
   - B. it has superconductivity
   - C. it has no Fermi surface
   - D. it is classical

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — 
3. **B** — 
4. **A** — 
5. **A** — 

</details>

## Ch 9 — Single-particle propagator revisited

1. Wick's theorem reduces a time-ordered product of free operators to…
   - A. a single contraction
   - B. the sum over all pairings
   - C. zero
   - D. a determinant only for bosons
2. The spectral function satisfies ∫A(k,ω)dω = …
   - A. 0
   - B. 1
   - C. N
   - D. Z
3. The pair bubble is non-zero only when…
   - A. both lines are particles
   - B. one line is a particle and the other a hole
   - C. both are holes
   - D. q = 0
4. Second-order Hubbard self-energy in imaginary time:
   - A. U G₀(τ)
   - B. −U² G₀↑(τ)G₀↓(τ)G₀↓(−τ)
   - C. U² G₀(τ)³
   - D. U n
5. Number of closed loops F and sign of the 2nd-order Hubbard self-energy diagram:
   - A. F=0, +
   - B. F=1, −
   - C. F=2, +
   - D. F=1, +

<details><summary>Answer key</summary>

1. **B** — For fermions the sum over pairings is a determinant — the basis of CT-INT.
2. **B** — 
3. **B** — 
4. **B** — Implemented in feynlab.selfenergy.local.second_order_sigma.
5. **B** — One bubble → (−1)¹.

</details>

## Ch 10 — Dyson, RPA, ladders

1. A proper (1PI) self-energy part…
   - A. can be cut in two by cutting one particle line
   - B. cannot be cut in two by cutting one particle line
   - C. has no loops
   - D. is always first order
2. Dyson's equation is…
   - A. G = G₀ + G₀ΣG
   - B. G = G₀Σ
   - C. Σ = G₀ G
   - D. G = G₀ + Σ
3. Low-density Fermi systems are dominated by…
   - A. bubble chains
   - B. ladder diagrams (K-matrix)
   - C. crossed diagrams
   - D. tadpoles
4. The RPA screened interaction is…
   - A. v(1 − vΠ)
   - B. v/(1 − vΠ) (Π with sign convention χ)
   - C. v²Π
   - D. Π/v
5. Friedel oscillations come from the singularity of the Lindhard function at…
   - A. q = 0
   - B. q = k_F
   - C. q = 2k_F
   - D. q = ∞

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — 
3. **B** — 
4. **B** — Bubble chain sum.
5. **C** — δn(r) ∝ cos(2k_F r)/r³.

</details>

## Ch 11 — Self-consistent renormalization

1. A skeleton diagram contains…
   - A. no self-energy insertions on internal lines
   - B. only ladders
   - C. only bubbles
   - D. no loops
2. How many of the 12 fourth-order Hubbard Σ diagrams are skeletons?
   - A. 12
   - B. 9
   - C. 6
   - D. 3
3. If Σ = Aω near μ (A<0) the jump of n_k at k_F is…
   - A. 1
   - B. 1/(1−A)
   - C. A
   - D. 0
4. Using non-skeleton diagrams with dressed G would…
   - A. be exact
   - B. double count
   - C. violate causality
   - D. remove the Fermi surface
5. GF2 is…
   - A. self-consistent second-order perturbation theory
   - B. Hartree only
   - C. CT-INT
   - D. BCS

<details><summary>Answer key</summary>

1. **A** — 
2. **B** — feynlab: count_diagrams(4,'hubbard','skeleton',True) = 9.
3. **B** — 
4. **B** — 
5. **A** — 

</details>

## Ch 12 — Ground-state energy of electron gas

1. Individual 2nd-order ring diagrams of jellium…
   - A. converge
   - B. diverge logarithmically
   - C. vanish
   - D. are imaginary
2. Gell-Mann–Brueckner correlation energy leading term (Ry):
   - A. 0.0622 ln r_s
   - B. −0.916/r_s
   - C. 2.21/r_s²
   - D. 0.0311 r_s
3. RPA high-density constant (Ha) reproduced by feynlab:
   - A. −0.0711
   - B. −0.048
   - C. −0.096
   - D. 0
4. Brueckner theory for nuclear matter sums…
   - A. ring diagrams
   - B. ladder diagrams with hard cores
   - C. crossed diagrams
   - D. tadpoles
5. The 2.21/r_s² term is…
   - A. exchange
   - B. kinetic energy
   - C. correlation
   - D. Madelung

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — 0.0311 ln r_s in Hartree.
3. **A** — 
4. **B** — 
5. **B** — 

</details>

## Ch 13 — Collective excitations

1. Plasmon dispersion in RPA: ω² = ω_p² + α v_F² q² with α =
   - A. 1/3
   - B. 3/5
   - C. 1
   - D. 9/5
2. Collective modes are zeros of…
   - A. G
   - B. the dielectric function ε(q,ω)
   - C. Σ
   - D. the Fermi function
3. Inside the particle–hole continuum a plasmon…
   - A. is undamped
   - B. decays (Landau damping)
   - C. becomes a phonon
   - D. disappears at q=0
4. The polarisation propagator is a special case of…
   - A. the vacuum amplitude
   - B. the two-particle Green's function
   - C. the self-energy
   - D. the T-matrix only
5. The retarded polarisation propagator gives…
   - A. the ground-state energy
   - B. linear density response
   - C. the lifetime of holes
   - D. T_c

<details><summary>Answer key</summary>

1. **B** — Mattuck Ex. 13.6; feynlab benchmark gives 0.6017.
2. **B** — 
3. **B** — 
4. **B** — 
5. **B** — 

</details>

## Ch 14 — Finite temperature

1. Fermionic Matsubara frequencies are…
   - A. 2nπ/β
   - B. (2n+1)π/β
   - C. nπ/β
   - D. n/β
2. G(τ) for fermions is…
   - A. periodic in β
   - B. antiperiodic in β
   - C. constant
   - D. zero at τ=β/2
3. (1/β)Σₙ e^{iωₙ0⁺}/(iωₙ − ε) equals…
   - A. f(ε)
   - B. 1 − f(ε)
   - C. ε
   - D. β
4. The average particle number from the grand potential:
   - A. N = −∂Ω/∂μ
   - B. N = ∂Ω/∂T
   - C. N = Ω/μ
   - D. N = βΩ
5. Why does feynlab subtract 1/iω tails before Fourier transforms?
   - A. for speed only
   - B. because G(iω) ~ 1/iω decays too slowly for a converged sum
   - C. to remove interactions
   - D. to enforce causality

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **A** — 
4. **A** — 
5. **B** — 

</details>

## Ch 15 — Superconductivity

1. Weak-coupling BCS ratio Δ(0)/k_BT_c:
   - A. 1.134
   - B. 1.764
   - C. 3.52
   - D. π
2. Perturbation theory in the normal state breaks down at T_c because…
   - A. the particle–particle ladder diverges
   - B. Σ vanishes
   - C. phonons vanish
   - D. μ changes
3. Nambu formalism uses…
   - A. 4×4 spin matrices
   - B. a 2×2 particle–hole spinor with anomalous propagators
   - C. only holes
   - D. Matsubara bosons
4. Bogoliubov quasi-particle energy:
   - A. ξ_k
   - B. √(ξ_k² + Δ²)
   - C. Δ
   - D. ξ_k + Δ
5. T_c = 1.134 ω_D e^{−1/λ}: the prefactor is…
   - A. 2e^γ/π
   - B. π/e^γ
   - C. 1
   - D. √2

<details><summary>Answer key</summary>

1. **B** — 2Δ/k_BT_c = 3.52.
2. **A** — 
3. **B** — 
4. **B** — 
5. **A** — 

</details>

## Ch 16 — Phonons from a many-body viewpoint

1. Bare ions oscillate at…
   - A. zero frequency
   - B. the ion plasma frequency
   - C. ω_D
   - D. the electron plasma frequency
2. Electronic screening turns the ion plasma mode into…
   - A. an optical mode
   - B. an acoustic phonon
   - C. a plasmon
   - D. a magnon
3. Bohm–Staver sound velocity scales as…
   - A. √(m/M) v_F
   - B. v_F
   - C. M/m
   - D. k_F²
4. The phonon-mediated electron interaction is attractive…
   - A. at all frequencies
   - B. below ~ω_D
   - C. only at q=0
   - D. never
5. This chapter connects most directly to…
   - A. Ch. 15 superconductivity
   - B. Ch. 2 pinball
   - C. Ch. 7 only
   - D. Ch. 19

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **A** — 
4. **B** — 
5. **A** — 

</details>

## Ch 17 — Phase transitions

1. A superconductor breaks…
   - A. spin rotation
   - B. global gauge symmetry
   - C. time translation
   - D. nothing
2. The Stoner criterion for ferromagnetism:
   - A. Uχ₀ = 1
   - B. U = 0
   - C. χ₀ = 0
   - D. U = W
3. In the ordered phase one needs…
   - A. anomalous propagators
   - B. only bare propagators
   - C. classical mechanics
   - D. no self-consistency
4. On the half-filled square lattice RPA first diverges at…
   - A. q = 0
   - B. q = (π, π)
   - C. q = (π, 0)
   - D. never
5. At the transition the two-particle propagator…
   - A. vanishes
   - B. diverges
   - C. is constant
   - D. becomes real

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — 
3. **A** — 
4. **B** — Nesting → antiferromagnetism (feynlab raises RuntimeError).
5. **B** — 

</details>

## Ch 18 — Kondo problem

1. The Kondo divergence comes from a term proportional to…
   - A. T²
   - B. ln(D/T)
   - C. 1/T
   - D. e^{−T}
2. T_K =
   - A. D e^{−1/(2|J|ρ)}
   - B. Jρ D
   - C. D/J
   - D. √(JD)
3. Why can't finite-order perturbation theory give T_K?
   - A. T_K is non-analytic in J
   - B. J is large
   - C. D is infinite
   - D. it can
4. The parquet approximation sums…
   - A. leading logarithms
   - B. only ladders
   - C. only bubbles
   - D. tadpoles
5. Which feynlab solver can follow an Anderson impurity into the Kondo regime numerically exactly?
   - A. Hubbard-I
   - B. IPT
   - C. CT-INT
   - D. Hartree–Fock

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — 
3. **A** — 
4. **A** — 
5. **C** — 

</details>

## Ch 19 — Renormalization group

1. The RG transformations here leave invariant…
   - A. the form of the Dyson equation
   - B. the energy
   - C. the density
   - D. nothing
2. Lie (flow) equations are…
   - A. differential equations for transformed quantities
   - B. integral equations
   - C. algebraic
   - D. stochastic
3. Solution of dg/dl = −2ρg² with g₀<0:
   - A. g₀ e^{−l}
   - B. g₀/(1 + 2ρg₀l)
   - C. g₀ + l
   - D. constant
4. The divergence scale l* = 1/(2ρ|g₀|) corresponds to…
   - A. T_K = D e^{−l*}
   - B. ω_D
   - C. E_F
   - D. T_c of BCS
5. Mattuck illustrates the RG on…
   - A. the high-density electron gas
   - B. the Ising model
   - C. QCD
   - D. the Kondo lattice

<details><summary>Answer key</summary>

1. **A** — 
2. **A** — 
3. **B** — 
4. **A** — 
5. **A** — 

</details>

## Appendices

1. Wick's theorem for fermions at fixed vertex times produces…
   - A. a permanent
   - B. a determinant
   - C. a trace
   - D. a sum of exponentials
2. The iδ in G₀ = 1/(ω − ε + iδ) encodes…
   - A. causality / boundary conditions in time
   - B. interaction
   - C. spin
   - D. temperature
3. Spectral density is related to the retarded G by…
   - A. A = −Im G^R/π
   - B. A = Re G
   - C. A = |G|²
   - D. A = G^A
4. Equation-of-motion decoupling (App. M) at lowest order gives…
   - A. Hubbard-I
   - B. RPA
   - C. BCS
   - D. GW
5. The time-development operator in the interaction picture is expanded as…
   - A. a Dyson (time-ordered exponential) series
   - B. a Taylor series in t only
   - C. a Fourier series
   - D. a continued fraction

<details><summary>Answer key</summary>

1. **B** — This is why CT-INT samples determinants.
2. **A** — 
3. **A** — 
4. **A** — 
5. **A** — 

</details>

## T1 — Getting started

1. feynlab's Green's-function convention is…
   - A. G(τ) = +⟨T c c†⟩
   - B. G(τ) = −⟨T c(τ)c†(0)⟩
   - C. G(τ) = ⟨c†c⟩
   - D. retarded only
2. Dyson equation in feynlab:
   - A. G = G₀ + G₀ΣG
   - B. G = G₀ − Σ
   - C. Σ = G − G₀
   - D. G⁻¹ = G₀⁻¹ + Σ
3. Units of spin/orbital Hall conductivity in feynlab:
   - A. e²/h
   - B. e/2π (i.e. (ħ/e)e²/h)
   - C. S/m
   - D. dimensionless
4. Why are tail coefficients passed to iw_to_tau?
   - A. to fix the chemical potential
   - B. to treat the slow 1/iω decay analytically
   - C. to speed up FFT
   - D. they are ignored
5. The τ→iω transform integrates…
   - A. the piecewise-linear interpolant exactly
   - B. with a rectangle rule
   - C. with Monte Carlo
   - D. only at τ=0

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — 
3. **B** — 
4. **B** — 
5. **A** — 

</details>

## T2 — Diagrams

1. A generic order-n diagram in feynlab is stored as…
   - A. a bijection between line emitters and receivers
   - B. a matrix of energies
   - C. a string
   - D. a graph image
2. Why can feynlab count distinct diagrams by dividing by n!·2ⁿ?
   - A. connected diagrams with fixed external legs have no non-trivial automorphisms
   - B. all diagrams are equal
   - C. by convention
   - D. it cannot
3. Sum of the 12 fourth-order Σ diagrams for the half-filled Hubbard atom is…
   - A. U⁴/16
   - B. zero
   - C. U²/4
   - D. divergent
4. Molinari's loop-resolved count at order 3 (Hartree propagator):
   - A. 10 + 9ℓ + ℓ²
   - B. 2 + ℓ
   - C. 74 + 91ℓ
   - D. 42
5. `has_tadpole()` detects…
   - A. Hartree-type insertions attached by one interaction line
   - B. crossed lines
   - C. skeletons
   - D. vertex corrections

<details><summary>Answer key</summary>

1. **A** — 
2. **A** — 
3. **B** — The exact atomic Σ is purely second order.
4. **A** — 
5. **A** — 

</details>

## T3 — Self-energy weak→strong

1. For the symmetric Anderson model, 1/Z − 1 to O(U²) is…
   - A. (3 − π²/4)(U/πΓ)²
   - B. U/Γ
   - C. (U/Γ)²
   - D. π²/4
2. The RPA/FLEX kernel for U→0 reduces to…
   - A. U
   - B. U²χ₀ (2nd order)
   - C. zero
   - D. U³χ₀²
3. IPT coexistence region found by feynlab (β=100):
   - A. U_c1≈2.5D, U_c2≈3.3D
   - B. no coexistence
   - C. U_c≈1D
   - D. U_c≈10D
4. For U ≫ W at half filling, Σ(iω) →
   - A. 0
   - B. U/2 + U²/(4iω)
   - C. U
   - D. iω
5. A RuntimeError from sigma('rpa') signals…
   - A. a bug
   - B. a magnetic (Stoner) instability Uχ₀ ≥ 1
   - C. wrong units
   - D. non-convergence of FFT

<details><summary>Answer key</summary>

1. **A** — Yamada–Yosida.
2. **B** — 
3. **A** — Exact NRG U_c2 = 2.94D; IPT overestimates.
4. **B** — 
5. **B** — 

</details>

## T4 — CT-INT QMC

1. Each CT-INT configuration represents…
   - A. one Feynman diagram
   - B. the sum of all diagrams with the given vertices (a determinant)
   - C. a single path
   - D. a classical spin state
2. The auxiliary shift δ ≈ 0.51 is used to…
   - A. speed up updates
   - B. avoid the fermionic sign problem for repulsive U
   - C. fix the density
   - D. remove tails
3. CT-INT cost scales roughly as…
   - A. ⟨k⟩³ with ⟨k⟩ ~ βU
   - B. linear in β
   - C. exponential in N
   - D. independent of U
4. feynlab's CT-INT error bars come from…
   - A. independent Markov chains
   - B. jackknife over frequencies
   - C. they are not provided
   - D. the sign
5. CT-INT was validated in feynlab against…
   - A. exact diagonalization of an Anderson impurity
   - B. experiment
   - C. Hubbard-I
   - D. nothing

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **A** — 
4. **A** — 
5. **A** — 

</details>

## T5 — Hall transport

1. Which mechanism scales as 1/n_i (∝ σ_xx)?
   - A. intrinsic
   - B. side jump
   - C. conventional (V₁³) skew scattering
   - D. anomalous distribution
2. Side jump and anomalous distribution for the massive Dirac model are…
   - A. equal
   - B. opposite
   - C. zero
   - D. unrelated
3. In the Rashba 2DEG with δ-impurities the total SHC is…
   - A. e/8π
   - B. zero (vertex corrections cancel)
   - C. e²/h
   - D. divergent
4. The orbital-quadrupole Hall response of a nonmagnetic crystal is…
   - A. always large
   - B. zero, because {Q,v} is T-odd
   - C. equal to the OHC
   - D. equal to the SHC
5. Which diagrams are NOT included in feynlab's Hall engines?
   - A. Born self-energy
   - B. ladder vertex
   - C. crossed X and Ψ diagrams
   - D. third-order skew vertex

<details><summary>Answer key</summary>

1. **C** — 
2. **A** — Sinitsyn et al. 2007, reproduced to 1e-16.
3. **B** — 
4. **B** — 
5. **C** — Ado et al. EPL 2015 — stated limitation.

</details>

## T6 — Electron gas, BCS, companions

1. rpa_correlation_energy(1.0) ≈
   - A. −0.0788 Ha
   - B. −0.5 Ha
   - C. 0
   - D. +0.0311 Ha
2. plasmon_dispersion returns NaN when…
   - A. q is too small
   - B. the plasmon enters the particle–hole continuum
   - C. r_s > 1
   - D. always
3. friedel_oscillation(r) decays asymptotically as…
   - A. e^{−r}
   - B. cos(2k_F r)/r³
   - C. 1/r
   - D. r²
4. kondo_log returns the factor…
   - A. 1 + 2Jρ ln(D/T)
   - B. e^{−T/D}
   - C. J²
   - D. ln T
5. bcs.gap_T(T) for T > T_c returns…
   - A. 0
   - B. Δ(0)
   - C. NaN
   - D. T

<details><summary>Answer key</summary>

1. **A** — 
2. **B** — 
3. **B** — 
4. **A** — 
5. **A** — 

</details>

## T8 — v0.2: crossed diagrams, ED-DMFT, polaron, bold

1. The X and Ψ diagrams describe…
   - A. ladder vertex corrections
   - B. skew scattering on pairs of close impurities
   - C. Born self-energy
   - D. phonon scattering
2. For the massive Dirac model, adding X and Ψ to the NCA result gives σ_xy =
   - A. −4εm(ε²+m²)/(ε²+3m²)²
   - B. −8εm³/(ε²+3m²)²
   - C. 0
   - D. −m/2ε
3. Why does ED-DMFT enforce Σ V_l² = ∫ε²ρ(ε) in the bath fit?
   - A. speed
   - B. an underweighted bath shrinks the effective bandwidth and shifts U_c down
   - C. it is required by Pauli
   - D. to fix the sign
4. Kozik–Ferrero–Georges showed that at strong coupling the skeleton series for Σ[G]…
   - A. diverges
   - B. converges to an unphysical branch
   - C. is exact at order 2
   - D. vanishes
5. The polaron energy in DiagMC is extracted from…
   - A. the average order
   - B. the exponential decay of G(p=0,τ)
   - C. the Z factor only
   - D. Feynman's variational bound

<details><summary>Answer key</summary>

1. **B** — Ado et al. EPL 2015.
2. **B** — Ψ = 0 for this model; X nearly cancels the NCA.
3. **B** — 
4. **B** — Reproduced in feynlab for the Hubbard atom at βU = 6.
5. **B** — 

</details>

## T7 — Benchmarks & limitations

1. Which v0.1 limitation is only PARTIALLY resolved in v0.2?
   - A. X/Ψ diagrams
   - B. MaxEnt
   - C. the semiclassical rate-by-rate Gaussian skew decomposition
   - D. polaron DiagMC
2. Which analytic continuation is implemented?
   - A. MaxEnt
   - B. Padé
   - C. stochastic optimisation
   - D. none
3. DMFT in feynlab is…
   - A. cluster DMFT
   - B. single-site, Bethe lattice
   - C. GW+DMFT
   - D. lattice QMC
4. Practical limit of brute-force generic diagram enumeration:
   - A. order 4
   - B. order 20
   - C. order 1
   - D. unlimited
5. What does the test suite guarantee?
   - A. that feynlab is the best code in the world
   - B. that every benchmarked number agrees with an exact or published result
   - C. nothing
   - D. experimental agreement

<details><summary>Answer key</summary>

1. **C** — The Kubo route gives it exactly; the w(4a) decomposition reproduces 1/2.
2. **B** — 
3. **B** — 
4. **A** — 
5. **B** — 

</details>
