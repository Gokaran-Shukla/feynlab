"""Writes chapter tutorials 00-06 (kindergarten part of Mattuck)."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
CH = {}

CH["00"] = ("The Many-Body Problem for Everybody", r"""
## Learning goals
* Say in one sentence what makes a many-body problem hard (the interaction couples every particle to every other).
* Distinguish **quasi particles** (a real particle dressed by a cloud of disturbance) from **collective excitations** (a coherent motion of the whole system, e.g. phonons, plasmons).
* Read the energy and lifetime of a quasi particle off a propagator.

## Key ideas (in our words)
1. *Fictitious bodies.* Strongly interacting real particles can often be replaced by weakly interacting fictitious ones — quasi particles and collective modes — with renormalised mass and finite lifetime. This is the whole strategy of the book.
2. *Dressing.* A quasi particle = bare particle + polarisation cloud. Its effective mass m* and lifetime τ follow from how strongly the cloud drags.
3. *Lifetime matters.* The concept is useful only if τ is long compared with the period ħ/ε of the excitation; near the Fermi surface this is guaranteed in a normal Fermi liquid (Im Σ ∝ (ε−μ)²).
4. *Collective modes* are not dressed particles at all; they carry quantum numbers of density or spin oscillations.

## Key formula
A quasi particle propagating with energy ε′ and lifetime τ gives
G(k,t) ≈ −i Z exp(−i ε′_k t − t/2τ_k)  (t > 0),
so |G|² decays as e^{−t/τ}; Z < 1 is the weight of the bare particle in the dressed state.

## feynlab hands-on
```python
import numpy as np
t = np.linspace(0, 40, 400); eps, tau, Z = 1.0, 10.0, 0.8
G = -1j*Z*np.exp(-1j*eps*t - t/(2*tau))
print("probability after one lifetime:", abs(G[np.searchsorted(t, tau)])**2 / Z**2)  # ~ e^-1
```
Later (Tutorial 3) you will compute Z directly: `BetheDMFT(beta=100,U=2).solve('ipt')['Z']`.

## Self-check
* Why is a quasi particle with τ·ε′ ≪ ħ physically meaningless?
* Give one quasi particle and one collective excitation in (a) a metal, (b) liquid ³He.
""")

CH["01"] = ("Feynman Diagrams, or How to Solve the Many-Body Problem by Means of Pictures", r"""
## Learning goals
* Define the single-particle propagator as a probability amplitude to go from (r₁,t₁) to (r₂,t₂).
* Build a propagator as a sum over "histories" drawn as diagrams (the drunken-man picture).
* Name the four propagators of the book: single-particle, two-particle, particle–hole, and the no-particle (vacuum) amplitude.

## Key ideas
1. *Propagators are the heroes.* Every measurable we need (quasi-particle energy, lifetime, collective modes, ground-state energy) can be extracted from a small set of propagators.
2. *Sum over paths.* The total amplitude is the sum of the amplitudes of all ways of getting from start to finish; each way is a picture (diagram) with a translation rule into a product of factors.
3. *Free propagation × interaction × free propagation …* is the universal grammar of every diagram.
4. *Particle–hole propagator* describes adding a particle and a hole simultaneously → collective modes (Ch. 13).
5. *Vacuum amplitude* → ground-state energy (Ch. 5, 12).

## Diagram grammar
line = free propagation (G₀), dot/wiggle = interaction, joining rule = multiply and integrate over intermediate positions/times (or sum over momenta).

## feynlab hands-on
```python
from feynlab.diagrams.enumerate import count_diagrams
for n in (1, 2, 3):
    print(n, count_diagrams(n, "generic", "G"))   # 2, 10, 74 histories (topologies) at order n
```

## Self-check
* Why does the number of diagrams explode factorially with order, and why is that not fatal? (Answer: partial summation, Ch. 2–4, 10.)
* Which propagator would you use to find a plasmon? To find the ground-state energy?
""")

CH["02"] = ("Classical Quasi Particles and the Pinball Propagator", r"""
## Learning goals
* Write the classical probability propagator P(r₂,r₁) of a ball in a pinball machine as a diagram series.
* Perform the first **partial summation**: summing an infinite subset of diagrams into a closed form.
* Understand why classical probabilities add (no interference) — the key difference from Ch. 3.

## Key ideas
1. *Free propagator × scattering probability × free propagator …* — every multi-bounce path is a diagram.
2. *Geometric series.* If each extra scattering at a given obstacle multiplies the path by the same factor p, the sum over 0,1,2,… repeats is P₀/(1−p): an infinite number of diagrams in one line.
3. *Renormalisation.* Repeated scattering by the same object turns the "bare" propagator into a "dressed" one — the classical seed of Dyson's equation.

## Key formula
P = P₀ + P₀ p P₀ + P₀ p P₀ p P₀ + … ⇒ P = P₀/(1 − p P₀)  (schematically; Mattuck's monkey/lion series).

## Guide to the book's exercises (2.1–2.4)
* 2.1: draw all paths with ≤ 2 strong scatterings; count at order n (each step chooses one of two strong scatterers → 2ⁿ).
* 2.3: with all free propagators equal to c, the series is geometric; sum it.
* 2.4: allow every animal → matrix geometric series; answer has the form (1 − cA)⁻¹.

## feynlab hands-on
```python
from feynlab.book.companions import pinball_propagator
exact, partial = pinball_propagator(P0=0.2, p_scatter=0.3, n_max=12)
print(exact, partial[:6])   # watch the partial sums converge to the resummed value
```

## Self-check
* When does the geometric partial sum fail to converge, and what does that mean physically?
""")

CH["03"] = ("Quantum Quasi Particles and the Quantum Pinball Propagator", r"""
## Learning goals
* Replace probabilities by amplitudes: G = ⟨final|e^{−iHt}|initial⟩, and |G|² the probability.
* Derive the quantum diagram series in (k, ω) space and see how energy denominators appear.
* Watch divergences of individual terms disappear after partial summation.
* Compute the energy shift and lifetime of an electron moving among random impurities (Sec. 3.5).

## Key ideas
1. *Fourier transform* turns time convolutions into products: G(ω) = G₀ + G₀VG₀ + G₀VG₀VG₀ + … .
2. *Disagreeable divergences.* Each term has poles at the unperturbed energy (G₀ = 1/(ω − ε₀ + iδ)) — higher powers diverge worse. Summing the series moves the pole: G = 1/(ω − ε₀ − V). A meaningful answer appears only after infinite summation.
3. *Random impurities.* Averaging over impurity positions restores translational invariance; the first-order term gives the average potential (energy shift N W₀/Ω), the second-order term has an imaginary part → finite lifetime (Fermi golden rule).

## Key formulas
* Dyson-type resummation: G(k,ω) = 1/(ω − ε_k − Σ(k,ω)).
* Impurity Born approximation: Re Σ = n W₀, 1/τ = 2π n W₀² ρ(ε) (contact impurities).

## Guide to exercises
* 3.1–3.4: particle in a box with a weak linear potential: build G₀, matrix elements V_{mn}, sum the diagonal geometric series → shifted levels; with only diagonal elements there is no lifetime.
* 3.7: the energy and lifetime formulas follow from first and second order in W with a step-function form factor.

## feynlab hands-on
```python
import numpy as np
from feynlab.book.companions import quantum_pinball, impurity_self_energy
w = np.linspace(-1, 1, 5)
print(quantum_pinball(w, eps0=0.1, V=0.2, n_order=3))   # truncated series
print(quantum_pinball(w, eps0=0.1, V=0.2))              # resummed: pole moved to 0.3
print(impurity_self_energy(k=1.0, n_imp=0.01, W0=0.5))  # (energy, 1/tau)
```
In `feynlab.transport.kubo_streda` the same Born self-energy (now a matrix in band space) sets the
broadening Γ that enters every Hall calculation.

## Self-check
* Why does the quantum series contain interference terms that the pinball series lacks?
* What happens to the lifetime if all impurities sit at identical positions in every sample (no averaging)?
""")

CH["04"] = ("Quasi Particles in Fermi Systems", r"""
## Learning goals
* Describe the non-interacting Fermi sea and particle–hole excitations.
* Learn the minimal occupation-number (second-quantisation) language.
* Extend the propagator to many fermions: particle lines (k > k_F) and hole lines (k < k_F), backward-in-time propagation.
* Understand Hartree and Hartree–Fock quasi particles and the first glimpse of the RPA.

## Key ideas
1. *Holes.* Removing an electron below k_F creates a hole that behaves like a particle with opposite momentum and energy; in diagrams a hole is a line running backwards in time.
2. *Virtual processes.* Intermediate states in diagrams need not conserve energy; momentum is conserved at every vertex.
3. *Hartree* (tadpole/bubble attached by one interaction line) = average field of all other particles; *Fock* (exchange) = Pauli correction. Both are first order, both are self-consistent fields.
4. *Electron gas.* At high density the dominant diagrams are chains of bubbles (RPA); they screen the Coulomb interaction (Sec. 4.9 previews Ch. 10).

## Key formulas
* Hartree–Fock energy: ε_k^{HF} = ε_k + Σ_{l<k_F}(V₀ − V_{k−l}) (direct − exchange).
* Particle and hole free propagators: G₀(k,ω) = θ(k−k_F)/(ω−ε_k+iδ) + θ(k_F−k)/(ω−ε_k−iδ).

## Guide to exercises
* 4.2: practise creation/annihilation signs on |11111000…⟩ — count occupied states to the left.
* 4.4: diagrams with a line beginning and ending on the same interaction without momentum transfer are forbidden for momentum-conserving interactions — check momentum at each vertex.
* 4.11: factorisable potential V_{kl} = λ u_k u_l → geometric series in λ Σ u² G₀ (closed form).

## feynlab hands-on
```python
from feynlab.diagrams.enumerate import unique_diagrams
for d in unique_diagrams(1, "generic", "proper", exclude_tadpoles=False):
    print(d.to_string(), "loops", d.n_loops)      # Hartree (1 loop) and Fock (0 loops)
from feynlab.selfenergy.local import hartree_sigma
print(hartree_sigma(U=2.0, n_other=0.5))           # Hubbard Hartree shift U<n>
```

## Self-check
* Why does a closed fermion loop carry a factor (−1)?
* Why is there no Fock term for the Hubbard contact interaction between opposite spins?
""")

CH["05"] = ("Ground State Energy and the Vacuum Amplitude", r"""
## Learning goals
* Interpret the vacuum amplitude R(t) = ⟨0|e^{−iHt}|0⟩ and extract the ground-state energy from its long-time behaviour.
* State and use the **linked-cluster theorem**: ln R is the sum of connected (linked) diagrams only.
* Classify diagrams as topologically distinct or equivalent.

## Key ideas
1. *E₀ from the phase.* For long times R(t) ~ e^{−iE₀t}, so ΔE = lim i (d/dt) ln R(t).
2. *Linked clusters.* Disconnected vacuum diagrams exponentiate; only connected ones contribute to ΔE (and the energy is extensive).
3. *Topological equivalence.* Diagrams that can be deformed into each other without cutting lines are counted once — the origin of symmetry factors.

## Key formula
ΔE = Σ (all connected, topologically distinct vacuum diagrams, with energy denominators).

## Guide to exercises
* 5.2: compare diagrams by their connectivity, not their drawing — feynlab's `canonical_key()` does exactly this.
* 5.5: with momentum conservation, a single line closing on one interaction carries zero momentum transfer; check which diagrams survive.

## feynlab hands-on
```python
from feynlab.diagrams.enumerate import _generic_labelled
n = 2
lab = [d for d in _generic_labelled(n)]
conn = [d for d in lab if d.is_connected()]
print(len(lab), len(conn))     # connected vs all labelled 2nd-order G-diagrams
```
(The ratio is the content of the linked-cluster theorem for the propagator: disconnected pieces cancel against the vacuum normalisation.)

## Self-check
* Why must the energy shift be proportional to N for an extensive system, and how does that rule out disconnected diagrams?
""")

CH["06"] = ("Bird's-Eye View of Diagram Methods", r"""
## Learning goals
* Map the field-theory ingredients (occupation numbers, propagators, diagrams, partial sums, self-consistency, finite T) onto the many-body questions they answer.
* Know which chapter treats which ingredient.

## Summary table (our paraphrase)
| ingredient | answers | chapters |
|---|---|---|
| occupation-number formalism | bookkeeping of many-fermion states | 4.3, 7 |
| single-particle propagator | quasi-particle energy & lifetime | 3, 4, 8–11 |
| vacuum amplitude | ground-state energy | 5, 12 |
| two-particle / polarisation propagator | collective modes, response | 13 |
| partial summation (Dyson, ladder, RPA) | strong/long-range interactions | 10 |
| self-consistent renormalisation | Fermi surface, Kondo, criticality | 11, 18 |
| finite temperature (Matsubara) | thermodynamics | 14 |
| anomalous propagators (Nambu) | superconductivity, broken symmetry | 15, 17 |

## feynlab map
Every row above corresponds to a module: `diagrams` (topology), `selfenergy`/`dmft` (single-particle propagator), `egas.rpa` (polarisation, ground-state energy), `egas.bcs` (Nambu), `core.matsubara` (finite T).

## Self-check
* Which propagator and which partial sum would you use for: (a) the plasmon, (b) the correlation energy of jellium, (c) T_c of a superconductor?
""")

for k, (title, body) in CH.items():
    with open(os.path.join(HERE, f"ch{k}.md"), "w") as f:
        f.write(f"# Chapter {int(k)} — {title}\n*Companion to R. D. Mattuck, A Guide to Feynman Diagrams in the Many-Body Problem (2nd ed.)*\n{body}")
print("wrote", len(CH))
