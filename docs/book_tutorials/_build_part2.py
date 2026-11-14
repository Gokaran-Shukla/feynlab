"""Writes chapter tutorials 07-19 and appendices."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
CH = {}

CH["07"] = ("Occupation Number Formalism (Second Quantization)", r"""
## Learning goals
* Write many-fermion states as occupation lists |n₁n₂…⟩ and act on them with c†, c.
* Derive anticommutation rules and the sign rule (count occupied states to the left).
* Express one- and two-body operators: H = Σ ε c†c + ½ Σ V c†c†cc.
* Use particle–hole operators and the field operators ψ(r).

## Key ideas
1. Antisymmetry is built into {c_i, c_j†} = δ_ij — no Slater determinants needed.
2. Every operator is a sum of products of c†, c with the single-particle matrix elements as coefficients.
3. Bosons: same structure with commutators.

## Guide to exercises
* 7.1–7.3: practise signs; ⟨Ψ|N|Ψ⟩ = |A|²·1 + |B|²·3 for the given superposition.
* 7.6: contact interaction λδ(r₁−r₂) → (λ/Ω) Σ c†_{k+q} c†_{l−q} c_l c_k.
* 7.11: build S = ½ Σ c†_{α} σ_{αβ} c_β.

## feynlab hands-on
`feynlab.qmc.ed` builds Jordan–Wigner matrices for c and c† — exactly the formalism of this chapter:
```python
from feynlab.qmc.ed import _build_ops
c = _build_ops(3)
import numpy as np
print(np.allclose(c[0] @ c[1].T + c[1].T @ c[0], 0), np.allclose(c[0] @ c[0].T + c[0].T @ c[0], np.eye(8)))
```

## Self-check
* Show that c†_i c†_i = 0 follows from the anticommutator.
""")

CH["08"] = ("More about Quasi Particles", r"""
## Learning goals
* Solve an exactly soluble fermion model (the pure Hartree model) and see quasi particles appear exactly.
* Estimate quasi-particle lifetimes from phase space: 1/τ ∝ (ε − ε_F)².
* Write the general quasi-particle propagator: a pole with weight Z plus an incoherent background.

## Key ideas
1. Near the Fermi surface Pauli blocking restricts decay channels → lifetime diverges as (ε−ε_F)⁻².
2. G(k,ω) = Z_k/(ω − ε′_k + i/2τ_k) + (smooth background): the quasi-particle pole is only part of the spectral weight.

## Guide to exercises
* 8.3: compare the RPA damping with the excitation energy to find how far from k_F quasi particles survive.
* 8.5: the exact energy of the Hartree model equals the sum of quasi-particle energies plus diagonal interactions.

## feynlab hands-on
```python
from feynlab.dmft.bethe import BetheDMFT
d = BetheDMFT(beta=200, U=2.0); r = d.solve('ipt')
print("Z =", r['Z'])   # weight of the quasi-particle pole
import numpy as np
w = np.linspace(-3, 3, 600); A = d.spectral_function(w)   # coherent peak + Hubbard bands
```
Fermi-liquid check: −Im Σ(iωₙ) is linear in ωₙ at low frequency (→ Z), the analytic counterpart of 1/τ ∝ ε².

## Self-check
* Where did the missing spectral weight 1 − Z go in the DMFT spectrum?
""")

CH["09"] = ("The Single-Particle Propagator Re-visited", r"""
## Learning goals
* Define G(k,t) = −i⟨Ψ₀|T c_k(t) c_k†(0)|Ψ₀⟩ and its Lehmann (spectral) representation.
* Derive the diagram expansion with Wick's theorem; learn the **diagram rules** (sign (−1)^F, factors i, energy conservation).
* Topology of diagrams; the chemical-potential formalism; the single pair-bubble (2nd-order) approximation.

## Key ideas
1. Spectral density A(k,ω) ≥ 0, ∫A dω = 1 (sum rule). Quasi particles = peaks in A.
2. Wick's theorem turns the time-ordered product into all pairings → diagrams; closed loops give (−1).
3. The pair bubble is the lowest diagram giving a lifetime; it is the kernel of IPT and of the 2nd-order lattice self-energy.

## Key formula (Matsubara version used in feynlab)
Σ⁽²⁾(τ) = −U² G₀↑(τ) G₀↓(τ) G₀↓(−τ) for the Hubbard interaction.

## Guide to exercises
* 9.1: time-order with signs for each fermion swap.
* 9.5: decide topological equivalence by labelled connectivity — compare `canonical_key()`s.
* 9.9: the bubble vanishes unless one line is a particle and the other a hole.

## feynlab hands-on
```python
from feynlab.diagrams.enumerate import unique_diagrams
d = unique_diagrams(2, "hubbard")[0]
print(d.to_string(), d.n_loops, d.sign)       # the pair bubble: one loop, sign -1
from feynlab.qmc.ed import AndersonED          # Lehmann representation in code:
ed = AndersonED(-1.0, 2.0, [-1, 1], [0.5, 0.5]); print(ed.green([1j], beta=10.0))
```

## Self-check
* Prove ∫A(k,ω)dω = 1 from the anticommutator.
""")

CH["10"] = ("Dyson's Equation, Renormalization, RPA and Ladder Approximations", r"""
## Learning goals
* Distinguish proper (1PI) and improper self-energy parts; write Dyson's equation G = G₀ + G₀ΣG.
* Sum ladder diagrams (low density: K-matrix / T-matrix) and bubble chains (high density: RPA).
* Define the dressed (effective) interaction W = v/(1 − vΠ) and evaluate the pair bubble (Lindhard function), Friedel oscillations.

## Key ideas
1. Dyson's equation sums **all** improper diagrams once Σ (proper parts) is known.
2. Low density → ladders dominate (repeated scattering of two particles).
3. High density / long-range Coulomb → bubble chains dominate; they screen v_q = 4π/q² into a short-range interaction.
4. The static Lindhard function has a log singularity at q = 2k_F → Friedel oscillations cos(2k_F r)/r³.

## Guide to exercises
* 10.1/10.5: proper = cannot be cut into two by cutting one particle (or one interaction) line.
* 10.4: for a contact interaction the K-matrix depends only on total momentum and frequency → algebraic equation.
* 10.7: factorisable potential → closed-form K.

## feynlab hands-on
```python
from feynlab.diagrams.enumerate import count_diagrams
print(count_diagrams(2,'generic','proper'), count_diagrams(2,'generic','improper'))   # 6, 4
from feynlab.book.companions import friedel_oscillation
import numpy as np; r = np.linspace(2, 8, 40); print(friedel_oscillation(r, rs=2.0)[:5])
from feynlab.selfenergy.lattice import HubbardLattice2D
L = HubbardLattice2D(nk=32, U=1.0, beta=5.0); s, info = L.sigma('rpa')   # bubble chains on a lattice
```

## Self-check
* Why does each term of the bare Coulomb series diverge while the RPA sum is finite?
""")

CH["11"] = ("Self-Consistent Renormalization and the Existence of the Fermi Surface", r"""
## Learning goals
* Replace bare lines by dressed ones ("clothed skeletons") and avoid double counting.
* Understand why quasi particles exist when perturbation theory converges, and why the Fermi surface survives interactions (Luttinger: a jump Z in n_k).
* Dressed vertices.

## Key ideas
1. Skeleton expansion: Σ[G] uses only diagrams without self-energy insertions, evaluated with the full G.
2. Self-consistency (e.g. self-consistent Hartree–Fock, GF2, FLEX, DMFT) is a partial sum to infinite order.
3. If Im Σ ∝ (ω−μ)², n_k has a discontinuity Z at k_F → a sharp Fermi surface.

## Guide to exercises
* 11.1: diagrams with self-energy insertions on internal lines must be omitted from the skeleton series.
* 11.2: with Σ = Aω + i sgn(μ−ω)B(ω−μ)², the jump is Z = 1/(1 − A).
* 11.3: self-consistent HF reproduces itself when ε_k^{HF} contains V_{kk}-type terms summed over occupied states.

## feynlab hands-on
```python
from feynlab.diagrams.enumerate import count_diagrams
print([count_diagrams(n,'hubbard','proper',True) for n in (2,3,4)],
      [count_diagrams(n,'hubbard','skeleton',True) for n in (2,3,4)])   # [1,2,12] vs [1,2,9]
from feynlab.selfenergy.lattice import HubbardLattice2D
L = HubbardLattice2D(nk=16, U=1.0, beta=5.0); L.sigma('2nd', self_consistent=True)   # GF2
```

## Self-check
* Which 3 of the 12 fourth-order Hubbard diagrams are not skeletons, and why would including them double count in GF2?
""")

CH["12"] = ("Ground State Energy of the Electron Gas and Nuclear Matter", r"""
## Learning goals
* Write the energy diagrams of jellium; see that 2nd and higher orders diverge.
* Reproduce the Gell-Mann–Brueckner idea: sum the most divergent (ring) diagrams to get a finite correlation energy ∝ ln r_s.
* Brueckner theory for nuclear matter: ladder sums for hard cores.

## Key formulas
E/N = 2.21/r_s² − 0.916/r_s + 0.0622 ln r_s − 0.096 + … (Rydberg; high density).
In RPA alone the constant is −0.142 Ry (−0.0711 Ha); the remaining part comes from second-order exchange.

## feynlab hands-on
```python
from feynlab.egas.rpa import rpa_correlation_energy
import numpy as np
for rs in (1e-3, 1e-2, 1.0, 2.0): print(rs, rpa_correlation_energy(rs))
e3, e4 = rpa_correlation_energy(1e-3), rpa_correlation_energy(1e-4)
print("slope", (e3-e4)/np.log(10), "exact", (1-np.log(2))/np.pi**2)
```

## Self-check
* Why is the ring series dominant at small r_s? (Each bubble brings a factor r_s but the Coulomb line at small q brings 1/q².)
""")

CH["13"] = ("Collective Excitations and the Two-Particle Propagator", r"""
## Learning goals
* Define the two-particle Green's function and the polarisation (density-fluctuation) propagator.
* Connect the retarded polarisation propagator to linear response.
* Find collective modes as poles: the plasmon in RPA, its dispersion ω² = ω_p² + (3/5)v_F²q².
* Express the two-particle propagator through the scattering amplitude.

## Key ideas
1. Collective modes = zeros of the dielectric function ε(q,ω) = 1 + v_qΠ(q,ω) outside the particle–hole continuum.
2. Inside the continuum the mode decays into particle–hole pairs (Landau damping).

## Guide to exercises
* 13.2: sum the bubble chain → v/(1 − vΠ₀).
* 13.6: expand Π₀ for ω ≫ q v_F to order q²/ω⁴ → the 3/5 coefficient.

## feynlab hands-on
```python
import numpy as np
from feynlab.egas import rpa
rs = 2.0; kF = rpa.kF_of_rs(rs)
q = np.linspace(0.02, 0.6, 15)*kF
w = rpa.plasmon_dispersion(rs, q)     # NaN once the plasmon enters the continuum
print(np.c_[q/kF, w/rpa.plasma_frequency(rs)])
```

## Self-check
* Why is there no plasmon in a neutral Fermi liquid with short-range interactions (and what replaces it)?
""")

CH["14"] = ("Fermi Systems at Finite Temperature", r"""
## Learning goals
* Replace ground-state averages by thermal (grand-canonical) ones; imaginary time τ ∈ [0, β).
* Matsubara frequencies ωₙ = (2n+1)π/β (fermions), νₘ = 2mπ/β (bosons); frequency sums via contour integration.
* The finite-T pair bubble and vacuum amplitude.

## Key ideas
1. G(τ) is antiperiodic: G(τ−β) = −G(τ). All diagram rules carry over with ∫dτ over [0,β) and (1/β)Σₙ.
2. Frequency sums of rational functions are done by residues with the Fermi function as the weight.

## Guide to exercises
* 14.1–14.2: N̄ = −∂Ω/∂μ, and N = Σ_k f(ε_k − μ) fixes μ(T).
* 14.6: residues of (1/β)Σₙ 1/((iωₙ−a)(iωₙ+iν−b)) give (f(a)−f(b))/(a−b+iν).

## feynlab hands-on
```python
from feynlab.core.matsubara import MatsubaraGrid, fermi
g = MatsubaraGrid(beta=10.0, n_iw=1024)
eps = 0.3
print(g.density(1/(g.iw-eps), c2=eps), fermi(eps, 10.0))   # Matsubara sum == Fermi function
```
Everything in feynlab (IPT, CT-INT, lattice Σ) is formulated this way.

## Self-check
* Why are bosonic frequencies even multiples of π/β?
""")

CH["15"] = ("Diagram Methods in Superconductivity", r"""
## Learning goals
* Review BCS: Cooper pairing, the gap, Bogoliubov quasi particles ("bogolons").
* See why ordinary perturbation theory breaks down below T_c.
* Use the Nambu 2×2 formalism; include retardation; compute T_c.

## Key formulas
* Gap equation 1 = λ∫₀^{ω_D} dξ tanh(E/2T)/E, E = √(ξ²+Δ²).
* T_c = 1.134 ω_D e^{−1/λ}, Δ(0)/k_BT_c = π e^{−γ} ≈ 1.764.

## Guide to exercises
* 15.2: bogolon operators are canonical — check {α, α†} = u² + v² = 1.
* 15.3: each bogolon adds E_k to the ground-state energy.

## feynlab hands-on
```python
from feynlab.egas import bcs
for lam in (0.2, 0.3): print(lam, bcs.gap_T0(lam)/bcs.Tc(lam))
```

## Self-check
* Which ladder (particle–particle or particle–hole) diverges at T_c, and why does that signal the breakdown of perturbation theory?
""")

CH["16"] = ("Phonons from a Many-Body Viewpoint", r"""
## Learning goals
* Treat ions + electrons as a coupled many-body system; derive phonons as collective modes of the ions screened by electrons.
* See the Bohm–Staver sound velocity emerge from RPA screening of ion plasma oscillations.
* Recognise the effective electron–electron interaction mediated by phonons (input to Ch. 15).

## Key ideas
1. Bare ions oscillate at the ion plasma frequency Ω_p; electronic screening ε(q) converts this into an acoustic mode ω = Ω_p/√ε(q,0) ≈ c_s q.
2. The same screened interaction, retarded, is attractive below ω_D → Cooper pairing.

## feynlab hands-on
```python
import numpy as np
from feynlab.egas import rpa
rs, M = 2.0, 1836.0*23   # e.g. sodium-like ion mass (a.u.)
kF = rpa.kF_of_rs(rs); n = kF**3/(3*np.pi**2); Omega_p = np.sqrt(4*np.pi*n/M)
q = 0.05*kF; eps = 1 + 4*np.pi/q**2*rpa.lindhard_iw(q, 1e-12, kF)
print("c_s =", Omega_p/np.sqrt(eps)/q, " Bohm-Staver:", np.sqrt(1/(3*M))*kF)
```

## Self-check
* Why is c_s ∝ √(m/M) v_F?
""")

CH["17"] = ("Quantum Field Theory of Phase Transitions in Fermi Systems", r"""
## Learning goals
* Understand spontaneous symmetry breaking (gauge symmetry for superconductors, spin rotation for ferromagnets).
* Anomalous propagators and why the normal perturbation series breaks down in the ordered phase.
* Generalised (matrix) propagators; the Stoner ferromagnet with a δ-interaction; divergence of the two-particle propagator at T_c.

## Key ideas
1. At the transition the scattering amplitude (ladder or bubble sum) diverges: 1 − Uχ₀ = 0 (Stoner criterion) or 1 − λΠ_pp = 0 (Cooper).
2. Below T_c one must expand around the broken-symmetry state with anomalous propagators.

## feynlab hands-on
```python
from feynlab.selfenergy.lattice import HubbardLattice2D
for U in (1.0, 1.5, 2.0):
    try:
        L = HubbardLattice2D(nk=32, U=U, beta=5.0); _, info = L.sigma('rpa'); print(U, info['stoner_max'])
    except RuntimeError as e:
        print(U, "instability:", e)     # U chi0(pi,pi) >= 1 -> antiferromagnetic transition
```

## Self-check
* Why does the instability appear first at q = (π, π) on the half-filled square lattice?
""")

CH["18"] = ("Feynman Diagrams in the Kondo Problem", r"""
## Learning goals
* Second-order (Born) scattering off a magnetic impurity contains ln(D/T) — the Kondo divergence.
* Parquet (leading-log) summation and its breakdown at T_K.
* Self-consistent renormalisation with pseudofermions; the strong-coupling picture.

## Key formulas
ρ(T) ∝ J²[1 + 2Jρ₀ ln(D/T) + …]; leading-log sum diverges at T_K = D e^{−1/(2|J|ρ₀)}.

## feynlab hands-on
```python
import numpy as np
from feynlab.book.companions import kondo_log
f, TK = kondo_log(np.logspace(-6, -1, 6), J_rho=0.1); print(f, "T_K =", TK)
from feynlab.selfenergy.local import yamada_yosida_Z
print(yamada_yosida_Z(U=2.0, Delta=1.0))   # Anderson model: weak-coupling Fermi liquid behind the Kondo resonance
```
The CT-INT solver (Tutorial 4) can follow the Anderson impurity into the Kondo regime (large U/Γ, low T): watch Z collapse exponentially.

## Self-check
* Why does no finite order of perturbation theory capture T_K?
""")

CH["19"] = ("The Renormalization Group", r"""
## Learning goals
* RG as multiplicative transformations of the Dyson equation that leave its form invariant.
* Lie (flow) equations and their solution for the high-density electron gas.

## Key ideas
1. Rescaling the propagator and interaction by factors that compensate each other defines a group.
2. The flow equation for the effective coupling resums leading logarithms — the same structure as the Kondo parquet sum.

## feynlab hands-on
Numerically integrate a one-loop flow dg/dl = −2ρ g² (Kondo-like) and compare with the closed form g(l) = g₀/(1 + 2ρg₀l):
```python
import numpy as np
from scipy.integrate import solve_ivp
rho, g0 = 1.0, -0.1
sol = solve_ivp(lambda l, g: -2*rho*g**2, (0, 4.9), [g0], dense_output=True)
l = np.linspace(0, 4.9, 6); print(sol.sol(l)[0], g0/(1 + 2*rho*g0*l))
```
The divergence at l* = 1/(2ρ|g₀|) is T_K = D e^{−l*} of Chapter 18.

## Self-check
* Why does the RG reproduce the leading-log parquet sum but not the strong-coupling fixed point?
""")

CH["appendices"] = ("Appendices A–N (intermediate part)", r"""
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
""")

for k, (title, body) in CH.items():
    head = f"# Chapter {int(k)} — {title}" if k.isdigit() else f"# {title}"
    with open(os.path.join(HERE, f"ch{k}.md" if k.isdigit() else "appendices.md"), "w") as f:
        f.write(f"{head}\n*Companion to R. D. Mattuck, A Guide to Feynman Diagrams in the Many-Body Problem (2nd ed.)*\n{body}")
print("wrote", len(CH))
