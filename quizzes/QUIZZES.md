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
5. Which diagrams were added to the NCA Hall engine in v0.2 (continuum) and v1.0 (lattice)?
   - A. Born self-energy
   - B. ladder vertex
   - C. crossed X and Ψ diagrams
   - D. third-order skew vertex

<details><summary>Answer key</summary>

1. **C** — 
2. **A** — Sinitsyn et al. 2007, reproduced to 1e-16.
3. **B** — 
4. **B** — 
5. **C** — Ado et al. EPL 2015; lattice version via BZ FFT.

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

## T9 — ED-DMFT and the Mott transition (T = 0)

1. Why can a T = 0 ED-DMFT solver use Lanczos instead of full diagonalization?
   - A. Lanczos is more accurate
   - B. only the ground state and its ±1-electron excitations enter G(iω) at T = 0
   - C. the Hamiltonian is diagonal
   - D. for Monte Carlo sampling
2. With an EVEN number of bath sites at half filling, the cluster ground state is…
   - A. a singlet
   - B. a spin doublet (odd number of sites)
   - C. a triplet
   - D. superconducting
3. Which bath parity should be used to locate U_c2 (where the metal disappears)?
   - A. even
   - B. odd
   - C. either
   - D. neither
4. Why does an odd-bath ED insulator disappear too late (spuriously)?
   - A. the fit diverges
   - B. the zero-energy bath level forms a Kondo resonance with the local moment, resolved at T = 0
   - C. the Lanczos method fails
   - D. the sum rule is violated
5. Why is the bath sum rule Σ V_l² = ∫ε²ρ(ε)dε enforced?
   - A. for speed
   - B. an under-weighted bath shrinks the effective bandwidth and pulls U_c down
   - C. it is required by Wick's theorem
   - D. to remove the sign problem
6. The consensus values for the Bethe-lattice Mott transition at T = 0 are…
   - A. U_c1 ≈ 2.39D, U_c2 ≈ 2.92–3.0D
   - B. U_c1 ≈ 3.3D, U_c2 ≈ 2.6D
   - C. U_c ≈ 1D
   - D. U_c1 = U_c2 = 4D

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — Odd number of electrons at half filling → doublet.
3. **B** — Odd baths include a zero-energy level and a singlet ground state that carries the Kondo quasiparticle.
4. **B** — 
5. **B** — 
6. **A** — 

</details>

## T7 — Benchmarks & limitations

1. Which v0.1 limitation is only PARTIALLY resolved in v0.2?
   - A. X/Ψ diagrams
   - B. MaxEnt
   - C. the semiclassical rate-by-rate Gaussian skew decomposition
   - D. polaron DiagMC
2. Which analytic-continuation methods does feynlab v1.0 provide?
   - A. MaxEnt only
   - B. MaxEnt, Nevanlinna and stochastic
   - C. Padé only
   - D. none
3. Which DMFT flavours does feynlab provide?
   - A. only single-site Bethe
   - B. single-site (Bethe/lattice), DCA clusters (exact cluster CT-INT), DFT+DMFT with charge self-consistency
   - C. GW+DMFT
   - D. lattice QMC only
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

## T10 — Infrastructure: errors, mixing, checkpoints, CLI

1. The naive error σ/√N of correlated Monte Carlo data must be multiplied by…
   - A. τ_int
   - B. √(2τ_int)
   - C. 2
   - D. 1/√τ_int
2. For an AR(1) process x_{t+1} = φx_t + noise, τ_int equals…
   - A. φ
   - B. 1/(1−φ)
   - C. (1+φ)/(2(1−φ))
   - D. φ²
3. The jackknife is needed for…
   - A. primary averages
   - B. errors of nonlinear functions of averages, e.g. ⟨a⟩/⟨b⟩
   - C. removing autocorrelation
   - D. speed
4. Why does linear mixing become slow near the Mott transition?
   - A. FFT errors
   - B. the largest Jacobian eigenvalue of the DMFT map approaches 1
   - C. the sign problem
   - D. memory
5. Anderson mixing is equivalent to…
   - A. Newton's method with the exact Jacobian
   - B. Pulay/DIIS extrapolation over previous residuals
   - C. simulated annealing
   - D. linear mixing with α = 1
6. Why does feynlab store input_json inside every HDF5 output?
   - A. compression
   - B. provenance: each result file carries its own input parameters
   - C. speed
   - D. GPU support

<details><summary>Answer key</summary>

1. **B** — Integrated autocorrelation time τ_int (Wolff's Γ-method).
2. **C** — feynlab's ar1_series test; φ = 0.9 gives 9.5.
3. **B** — 
4. **B** — Anderson/Broyden: 78 → 13 iterations at U = 2.6.
5. **B** — 
6. **B** — 

</details>

## T11 — CT-HYB, Kanamori, DFT+DMFT

1. CT-HYB expands the partition function in powers of…
   - A. U
   - B. the hybridization Δ(τ)
   - C. temperature
   - D. the hopping t of the lattice
2. Why is CT-HYB preferred at strong coupling?
   - A. its average order DECREASES with U
   - B. it has no sign problem ever
   - C. it works at T = 0
   - D. it needs no bath
3. In the segment picture the local trace for density–density U is…
   - A. a determinant
   - B. exp(−U·overlap + μ·length)
   - C. always 1
   - D. a Pfaffian
4. Two-orbital Kanamori atom with N = 2: the ground multiplet energy is…
   - A. U
   - B. U − 3J (spin triplet)
   - C. U + J
   - D. U − J
5. The t₂g³ (three electrons, three orbitals) ground state is…
   - A. a singlet at 3U − 3J
   - B. the 4-fold ⁴S state at 3U − 9J
   - C. a doublet at 3U
   - D. non-degenerate
6. In Hubbard-I DFT+DMFT, why must the 'occupation' double-counting sit at the plateau CENTRE?
   - A. for speed
   - B. inside a Mott gap every shift gives the same n; the edge is unstable
   - C. required by FLL
   - D. to reduce J

<details><summary>Answer key</summary>

1. **B** — Werner et al. PRL 97, 076405 (2006).
2. **A** — 
3. **B** — 
4. **B** — Hund's first rule.
5. **B** — 
6. **B** — 

</details>

## T12 — DCA, two-particle response, NRG

1. DCA with N_c = 1 reduces to…
   - A. Hartree–Fock
   - B. single-site DMFT
   - C. RPA
   - D. the atomic limit
2. In 2×2 DCA the more strongly renormalised patch is…
   - A. (0,0)
   - B. (π,0) — the antinode
   - C. (π,π)
   - D. none
3. The compressibility sum rule states…
   - A. χ₀(q→0, Ω=0) = ∂n/∂μ
   - B. χ₀ = 0
   - C. χ₀(π,π) = n
   - D. χ = U
4. The Wilson-chain hopping t_n decays as…
   - A. 1/n
   - B. Λ^{−n/2}
   - C. e^{−n²}
   - D. it is constant
5. For the symmetric Anderson model, S_imp on the local-moment plateau equals…
   - A. 0
   - B. ln 2
   - C. ln 4
   - D. 1
6. Haldane's T_K contains the exponential factor…
   - A. exp(−U/Γ)
   - B. exp(−πU/8Γ + πΓ/2U)
   - C. exp(−Γ/U)
   - D. exp(−1/Jρ) with no U

<details><summary>Answer key</summary>

1. **B** — feynlab: Z 0.486644 vs 0.486638.
2. **B** — Coupling to (π,π) spin fluctuations; seed of the pseudogap.
3. **A** — 
4. **B** — 
5. **B** — Free orbital: ln 4; Fermi liquid: 0.
6. **B** — feynlab NRG reproduces the slope −1.00.

</details>

## T13 — Analytic continuation

1. Analytic continuation is hard because…
   - A. the kernel 1/(iω−ω) is singular
   - B. the inverse problem is exponentially ill-conditioned
   - C. of the fermion sign
   - D. of the FFT
2. The Nevanlinna method exploits that…
   - A. G is a Nevanlinna (Herglotz) function, mapped to a Schur function
   - B. A is Gaussian
   - C. the noise is white
   - D. G is real
3. The Pick criterion tells you…
   - A. the best α in MaxEnt
   - B. whether a causal interpolant of the data exists
   - C. the number of poles
   - D. the temperature
4. Why does feynlab's Nevanlinna use mpmath with 64 digits?
   - A. style
   - B. the Schur recursion loses ~1 digit per point
   - C. speed
   - D. GPU
5. Stochastic continuation picks the fictitious temperature Θ where…
   - A. χ² is minimal
   - B. χ²/N ≈ 1
   - C. the entropy is maximal
   - D. A is smoothest
6. Recommended practice in feynlab:
   - A. trust any single method
   - B. MaxEnt with two default models + stochastic cross-check (Nevanlinna for near-exact data)
   - C. always Padé
   - D. fit Lorentzians by eye

<details><summary>Answer key</summary>

1. **B** — 
2. **A** — Fei, Yeh, Gull PRL 126, 056402 (2021).
3. **B** — 
4. **B** — 
5. **B** — 
6. **B** — 

</details>

## T14 — CDet, resummation, polaron Z and m*

1. Rossi's CDet recursion obtains the sum of connected diagrams at cost…
   - A. n!
   - B. 3ⁿ
   - C. 2^(n²)
   - D. n³
2. For the Hubbard atom, X(U) = ⟨n↑n↓⟩ − ¼ is…
   - A. −¼ tanh(βU/4)
   - B. U/4
   - C. e^{−βU}
   - D. 0
3. The Domb–Sykes plot of a_n/a_{n−1} vs 1/n gives…
   - A. the sum of the series
   - B. 1/R at 1/n → 0
   - C. the sign
   - D. Z
4. A conformal map is used to…
   - A. speed up DiagMC
   - B. sum a series beyond its radius of convergence
   - C. remove tails
   - D. fix μ
5. The polaron Z factor in feynlab is normalised using…
   - A. Feynman's bound
   - B. the exactly known order-0 diagram
   - C. experiment
   - D. sum rule on m*
6. At weak coupling the Fröhlich effective mass is…
   - A. 1 − α
   - B. 1 + α/6 + 0.0236α²
   - C. e^{α}
   - D. 1 + α

<details><summary>Answer key</summary>

1. **B** — Inclusion–exclusion over subsets.
2. **A** — CDet matches every Taylor coefficient.
3. **B** — 
4. **B** — 
5. **B** — 
6. **B** — Rösler 1968; Smondyrev 1986.

</details>

## T15 — Advanced transport

1. Expanding the T-matrix to third order in V produces…
   - A. side jump
   - B. conventional skew scattering ∝ n_i V³
   - C. intrinsic AHE
   - D. nothing new
2. Why is the crossed Ψ diagram nonzero for QWZ near Γ but zero for Ado's Dirac model?
   - A. numerical error
   - B. the quadratic mass Δ − k²/2 of the lattice model
   - C. temperature
   - D. spin–orbit
3. The gapped Dirac valley orbital moment is m(k) =
   - A. Δv²/[2(Δ² + v²k²)], equal in both bands
   - B. zero
   - C. opposite in the two bands
   - D. v k
4. The Berry-curvature dipole Fermi-sea and Fermi-surface forms agree because of…
   - A. integration by parts
   - B. time reversal
   - C. particle–hole symmetry
   - D. Wick's theorem
5. For the Rashba 2DEG with both bands occupied, |δS_y| per eEτ is…
   - A. mα/(4π), independent of μ
   - B. ∝ μ
   - C. zero
   - D. ∝ α²
6. Adaptive k refinement in WannierTB3D refines points where…
   - A. k is large
   - B. |Ω| is large or a band lies within a window of μ
   - C. randomly
   - D. at Γ only

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — 
3. **A** — Xiao, Yao, Niu PRL 99, 236809 (2007).
4. **A** — 
5. **A** — Edelstein 1990.
6. **B** — Wang, Yates, Souza, Vanderbilt PRB 74, 195118 (2006).

</details>

## T16 — v1.1 solvers: matrix CT-HYB, cluster CT-INT, CSC, vertices

1. Why can the matrix CT-HYB treat spin-flip and pair-hopping terms while the segment solver cannot?
   - A. it uses a larger bath
   - B. it evaluates the full local trace in the eigenbasis of H_loc instead of segment overlaps
   - C. it works at T = 0
   - D. it neglects them
2. A 'continue' statement on impossible moves (removal at k = 0) biased G by 10%. Why?
   - A. numba bug
   - B. those steps skipped the measurement, under-sampling low expansion orders
   - C. the determinant was wrong
   - D. the trace overflowed
3. In cluster CT-INT for DCA, the Weiss field in the determinant is…
   - A. local only
   - B. G~0(r_p - r_q, tau_p - tau_q), cluster-translation invariant
   - C. the lattice G(k)
   - D. the bare G0(k)
4. charge self-consistency at U_DMFT = 0 with a Hartree functional must reproduce…
   - A. the atomic limit
   - B. static (Hartree) mean field
   - C. DMFT with IPT
   - D. nothing
5. The imaginary-time integrals in the Lehmann representation of G2 are computed as…
   - A. Monte Carlo averages
   - B. divided differences of exp(beta z), via exp of a bidiagonal matrix
   - C. FFTs
   - D. Pade approximants
6. The DMFT Bethe-Salpeter chi_m(q = 0) equals dm/dh of the DMFT magnetisation because…
   - A. of the sum rule
   - B. DMFT is thermodynamically consistent (Phi-derivable with a local irreducible vertex)
   - C. coincidence
   - D. of particle-hole symmetry
7. In ladder DGA the kernel K = F_d/2 - 3F_m/2 - F_d^loc/2 + F_m^loc/2 reduces at lowest order to…
   - A. 2U
   - B. U
   - C. 0
   - D. 3U/2
8. TPSC fixes U_sp and U_ch by…
   - A. fitting QMC
   - B. the local spin and charge sum rules plus U_sp = U<n_up n_dn>/(<n_up><n_dn>)
   - C. RPA
   - D. Hubbard-I

<details><summary>Answer key</summary>

1. **B** — 
2. **B** — Every Monte Carlo step must reach the measurement.
3. **B** — 
4. **B** — feynlab: agreement 6e-6.
5. **B** — Stable for degenerate energies.
6. **B** — 2.9089 vs 2.9102 in feynlab.
7. **B** — So DGA -> second-order PT as U -> 0.
8. **B** — Vilk & Tremblay 1997.

</details>

## T17 — v1.1 transport, continuation, polarons

1. Why did v1.0 conclude that the semiclassical Gaussian skew was 'half missing'?
   - A. an approximation in the rate
   - B. the reference formula had been transcribed with an extra factor 2
   - C. numerical noise
   - D. missing crossed diagrams
2. For Gaussian-shaped (finite-range) impurities the Born rate is multiplied by…
   - A. exp(-q^2 d^2/2)
   - B. q^2
   - C. 1/q^2
   - D. nothing
3. Smooth disorder reduces the side-jump contribution because…
   - A. the DOS changes
   - B. large-momentum transfer, which carries interband coherence, is suppressed
   - C. time-reversal symmetry
   - D. the Fermi velocity grows
4. The generic-lattice crossed-diagram code was validated by…
   - A. experiment
   - B. writing QWZ in an oblique unit cell and in 3D stacked layers and recovering the square-lattice result
   - C. comparing with Boltzmann
   - D. it was not
5. SAC (Shao-Sandvik) chooses the sampling temperature Theta where…
   - A. chi^2 is minimal
   - B. <chi^2> = chi^2_min + a sqrt(2 chi^2_min), a ~ 0.5
   - C. Theta = 1
   - D. entropy is maximal
6. The momentum-fluctuation estimator gives m*/m from…
   - A. finite p differences
   - B. m/m* = 1 - lim <X^2>/(3 tau), X = sum q (t2 - t1) over phonon arcs
   - C. Z
   - D. the average order
7. Variational ED for the Holstein polaron builds its basis by…
   - A. random sampling
   - B. repeated application of H to |electron, no phonons> in a translation-invariant frame
   - C. Lanczos only
   - D. DMRG

<details><summary>Answer key</summary>

1. **B** — Validated by feeding the transcribed rate into the Boltzmann solver.
2. **A** — 
3. **B** — 
4. **B** — 
5. **B** — 
6. **B** — alpha = 0.2: 1.0340 +- 0.0002 vs 1.0343.
7. **B** — Bonca-Trugman-Batistic; E0 = -2.46968472393287 at t = w = g = 1.

</details>

## T18 — Charge self-consistent DFT+DMFT with Quantum ESPRESSO

1. Which QE input flag makes pw.x stop after the nscf step and later read a DMFT density correction?
   - A. lda_plus_u
   - B. dmft = .true.
   - C. tqr
   - D. noncolin
2. What does feynlab write into dft_update/delta_N?
   - A. the self-energy
   - B. Delta N_ij(k) = N^DMFT_ij(k) - delta_ij f^DFT_i(k) in the Bloch band window, per spin
   - C. the Wannier Hamiltonian
   - D. the total energy
3. How does QE turn f + Delta N into a new density?
   - A. it ignores Delta N
   - B. it diagonalises f + Delta N, rotates the Kohn-Sham states with the eigenvectors and occupies them with the eigenvalues
   - C. it mixes potentials
   - D. it reruns Wannier90
4. Why should QE use Fermi-Dirac smearing with degauss = k_B T of the DMFT?
   - A. speed
   - B. so that with Sigma = 0 the DMFT occupations equal QE's and Delta N vanishes identically
   - C. it is required by Wannier90
   - D. to avoid negative occupations
5. Why is the per-orbital t2g occupation a bad convergence criterion for SrVO3?
   - A. it is noisy
   - B. cubic symmetry pins it at n/3, so it cannot change
   - C. QE does not print it
   - D. it diverges
6. In the QE DMFT step, which pw.x settings are used?
   - A. from_scratch, many electron steps
   - B. restart_mode='restart', startingwfc/pot='file', electron_maxstep=1, small mixing_beta
   - C. nscf
   - D. bands
7. A fitted hybridisation Delta(iw) = sum_l V_l^2/(iw - e_l) is used in an ED impurity with chemical potential mu. In which frame are the e_l?
   - A. absolute energies, so subtract mu once more
   - B. already measured from mu: subtracting mu again empties or fills the bath as mu drifts
   - C. measured from the Fermi level of QE
   - D. it does not matter
8. Why was ED with one bath site per t2g orbital replaced by Lanczos with two?
   - A. speed
   - B. one site per orbital is too coarse to fit Delta(iw) and drives a spurious Mott state (Z ~ 0.01)
   - C. QE requires it
   - D. to reduce the sign problem

<details><summary>Answer key</summary>

1. **B** — QE >= 7.0, with dmft_prefix naming the HDF5 archive.
2. **B** — 
3. **B** — PW/src/add_dmft_occ.f90 and sum_band.f90.
4. **B** — feynlab: Delta N = 1.5e-7, energy unchanged to 1e-8 Ry.
5. **B** — Use the scf accuracy, the total energy and the window charge.
6. **B** — 
7. **B** — This was the v1.2 bug fix in DFTDMFT: the bath emptied as mu drifted to about -1 eV.
8. **B** — Two bath sites per orbital give 9 sites, handled by the sector Lanczos solver.

</details>

## T19 — lambda-corrected ladder DGA

1. Which physical defect of the plain ladder DGA does the Moriya lambda correction remove?
   - A. the Luttinger theorem violation
   - B. the divergence of chi_m(Q) at the DMFT Neel temperature (no Mermin-Wagner) and the violated local spin sum rule
   - C. the sign problem
   - D. the high-frequency tail of Sigma
2. How is chi_m modified?
   - A. chi -> chi + lambda
   - B. chi^lambda = 1/(1/chi + lambda)
   - C. chi -> lambda chi
   - D. only chi(Q) is changed
3. Which condition fixes lambda_m?
   - A. Luttinger count
   - B. (1/beta) sum_{q,w} chi_m^lambda(q,w) equals the local DMFT spin sum (1/beta) sum_w chi_m^loc(w)
   - C. chi(Q) = 1
   - D. Sigma(k) local part = Sigma_DMFT
4. Below the DMFT Neel temperature the plain ladder gives chi_m(Q,0) < 0. What must lambda_m satisfy?
   - A. lambda < 0
   - B. lambda > -1/chi_m(Q,0) > 0
   - C. lambda = 0
   - D. lambda = U
5. Why does feynlab write the ladder self-energy with triangular vertices gamma_r?
   - A. for speed only
   - B. lambda enters only through the physical chi_r(q,w) inside [1 - gamma_r(1 -+ U chi_r)]; at lambda = 0 it equals the F(q) form to machine precision
   - C. to avoid the local subtraction
   - D. because F(q) is unknown
6. Why is the default sum-rule target the DMFT local chi_m summed over the same frequency box (plus the bubble tail) rather than (n - 2d)/2?
   - A. it is cheaper
   - B. truncation errors of the vertex then cancel between lattice and local sides; the exact target applied to a truncated lattice chi is inconsistent
   - C. (n - 2d)/2 is wrong
   - D. TPSC requires it

<details><summary>Answer key</summary>

1. **B** — Katanin, Toschi, Held PRB 80, 075104 (2009).
2. **B** — A mass term in the Ornstein-Zernike form, same for all q and w.
3. **B** — 
4. **B** — At U = 8t, T = 0.31t: lambda_m = 0.124 > lambda_min = 0.089; chi_m(Q) = 29 instead of -11.
5. **B** — Sum_nu' chi_r(q)_{nu nu'} = chi0_nu gamma_r (1 - s_r U chi_r).
6. **B** — The spread between the two targets measures the box error.

</details>

## T20 — 5-orbital CT-HYB with a sector-blocked trace

1. Why does the dense matrix CT-HYB fail for a full d shell?
   - A. sign problem
   - B. every operator in the trace is a dense 1024 x 1024 product
   - C. the determinant is too large
   - D. the hybridization is complex
2. Which quantum numbers block the Kanamori H_loc exactly (spin flip and pair hopping included)?
   - A. N only
   - B. N_up and N_dn
   - C. each orbital occupation
   - D. S^2 and L^2
3. For a given operator string, what determines the sequence of blocks?
   - A. nothing; all blocks mix
   - B. the starting block at tau = 0: each c_f / c+_f maps one block to exactly one other block
   - C. the random number
   - D. the bath
4. What is the lazy trace evaluation?
   - A. skipping measurements
   - B. rejecting a move when even an upper bound of |Tr| cannot pass the Metropolis test, before the exact trace is computed
   - C. computing the trace in single precision
   - D. truncating the expansion order
5. Why can only a few outer eigenstates of a block be propagated?
   - A. they are degenerate
   - B. their contribution carries e^{-l E_n} with l the wrap-around interval; a rigorous bound tells which are negligible, turning matrix-matrix into matrix-vector products
   - C. QMC ignores them
   - D. they have zero norm
6. In the 5-orbital ED benchmark, why are n, docc and Sz^2 accurate while G(iw) is noisy at e_b = -2, V = 0.5?
   - A. a bug in the trace
   - B. the M = F^-1 estimator has heavy tails when the order is low and |Delta(tau)| varies by e^20; local observables come from the trace and are unaffected
   - C. ED is wrong
   - D. the bath is too large

<details><summary>Answer key</summary>

1. **B** — 4^5 = 1024 atomic states.
2. **B** — Spin flip and pair hopping move electrons between orbitals but keep N_up and N_dn.
3. **B** — Strings that annihilate a block are discarded with integer look-ups.
4. **B** — Semon et al. PRB 90, 075149 (2014); exact, not an approximation.
5. **B** — 
6. **B** — Worm sampling or improved estimators are the cure.

</details>

## T21 — Two-particle vertex from CT-HYB

1. Why is the Lehmann-ED vertex limited to small baths?
   - A. sign problem
   - B. the number of 4-state chains grows like (number of states)^4 before pruning; 256 states already exhaust memory
   - C. it needs a Fourier transform
   - D. ED has no frequencies
2. What is the CT-HYB removal estimator for G2_{ss'}(nu, nu', w)?
   - A. <M>^2
   - B. (1/beta)<Mh_s(nu+w,nu) Mh_s'(nu',nu'+w) - delta_ss' Mh_s(nu+w,nu'+w) Mh_s(nu',nu)>
   - C. the derivative of G with respect to Delta
   - D. a ladder sum
3. Check of the estimator at U = 0?
   - A. chi_ud = chi0
   - B. connected part vanishes: chi_ud = 0, chi_uu = chi0 delta_{nu nu'}
   - C. G2 = 0
   - D. chi_m diverges
4. How can one check a measured vertex without any reference solver?
   - A. compare with TPSC
   - B. the Schwinger-Dyson equation: Sigma from F_ud must equal Sigma from Dyson's equation
   - C. the sign
   - D. the order histogram
5. Why did a larger frequency box (only affordable with CT-HYB) matter for lambda-DGA?
   - A. it changed U
   - B. the box local spin sum approaches the exact (n - 2d)/2, and the two sum-rule variants of lambda-DGA converge (spread factor 2.2 -> 1.4)
   - C. it removes the need for lambda
   - D. for speed
6. What changed in the lattice DMFT susceptibility when the 2-site ED bath was replaced by the continuous CT-HYB bath at the same box?
   - A. nothing
   - B. the plain DMFT chi_m(Q) moved from 6.4 to 8.0: the bath discretization error
   - C. lambda became negative
   - D. Sigma changed sign

<details><summary>Answer key</summary>

1. **B** — Two bath sites take minutes; three run out of memory.
2. **B** — Removing two rows and two columns of the hybridization matrix gives a 2x2 minor of M.
3. **B** — This is a unit test.
4. **B** — 1-2% at U = 4t in the benchmark.
5. **B** — 
6. **B** — 

</details>
