# Tutorial 3 — Self-energy from weak to very strong correlation

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


| regime (U vs bandwidth W) | method in feynlab | diagrams summed | module |
|---|---|---|---|
| weak, U ≪ W | Hartree–Fock | tadpole + exchange | `selfenergy.local.hartree_sigma` |
| weak | 2nd order (2OPT) | single pair bubble | `selfenergy.local.second_order_sigma`, `selfenergy.lattice` (`'2nd'`) |
| weak–moderate, near magnetism | RPA / one-shot FLEX | bubble chains + ph ladders | `selfenergy.lattice` (`'rpa'`), optional self-consistency |
| moderate, U ~ W | DMFT + IPT | 2nd order in a self-consistent bath | `dmft.BetheDMFT(...).solve('ipt')` |
| strong, U ≫ W | Hubbard-I / atomic limit | atomic self-energy | `solve('hubbard1')`, `local.atomic_sigma` |
| any (numerically exact local) | DMFT + CT-INT | **all** diagrams, sampled | `solve('ctint')` |

## Weak coupling: 2nd order and its exact checks

```python
from feynlab.core.matsubara import MatsubaraGrid
from feynlab.selfenergy import local
g = MatsubaraGrid(400.0, 4096)
G0 = 1/(g.iw + 1j*1.0)                 # symmetric Anderson model, flat band, Gamma = 1
s2 = local.second_order_sigma(g, G0, U=1.0)
Z = local.z_from_matsubara_extrap(g.iw, s2, 4)
print(1/Z - 1, local.yamada_yosida_Z(1.0, 1.0)**-1 - 1)     # 1/Z-1 = (3 - pi^2/4)(U/pi Gamma)^2
```

## Lattice self-energies (2D Hubbard)

```python
from feynlab.selfenergy.lattice import HubbardLattice2D
L = HubbardLattice2D(nk=32, U=1.0, beta=5.0, n_iw=256)   # half filling
sig, info = L.sigma('2nd')            # Sigma(k, iw_n) on the full k grid
sig_rpa, info = L.sigma('rpa')        # raises if U chi0 >= 1 (Stoner/AF instability)
print(info['stoner_max'])
```
All convolutions are FFTs in (r, τ): Σ(r,τ) = W(r,τ)·G(r,τ), with χ₀(r,τ) = G(r,τ)G(−r,β−τ).
The RPA/FLEX kernel W = U²[3/2 χ₀/(1−Uχ₀) + 1/2 χ₀/(1+Uχ₀) − χ₀] reduces to U²χ₀ (2nd order) for U→0 — checked in the benchmarks.
`self_consistent=True` turns '2nd' into GF2 and 'rpa' into FLEX.

## Moderate coupling: DMFT + IPT and the Mott transition

```python
from feynlab.dmft.bethe import BetheDMFT
for U in (2.0, 2.6, 3.0, 3.4):
    met = BetheDMFT(beta=100, U=U).solve('ipt')
    ins = BetheDMFT(beta=100, U=U).solve('ipt', insulating_start=True)
    print(U, met['Z'], ins['Z'])
```
You will find a coexistence window (≈2.4–2.6 < U < 3.2–3.4 in units of the half bandwidth D) — the
first-order Mott transition. Exact NRG values are U_c1 ≈ 2.5D and U_c2 ≈ 2.94D; IPT is known to
overestimate U_c2. See `benchmarks/figures/ipt_mott_hysteresis.png`.

Spectral function (Padé continuation): `d.spectral_function(np.linspace(-4,4,800))`.

## Strong coupling: Hubbard-I and the atomic limit

For U ≫ W, Σ(iω) → U/2 + U²/(4iω) (half filling): two Hubbard bands at ±U/2 and no quasiparticle.
`BetheDMFT(beta=50, U=6).solve('hubbard1')` reproduces it; compare to `local.atomic_sigma`.

## Numerically exact: DMFT + CT-INT

```python
d = BetheDMFT(beta=10, U=2.0, n_iw=1024)
r = d.solve('ctint', n_iter=10, mix=0.7,
            ctint_kwargs=dict(n_meas=400000, n_chains=4, n_iw_meas=64))
print(r['Z'], d.ctint_last['k_mean'], d.ctint_last['sign'])
```
CT-INT measures Σ on the lowest `n_iw_meas` frequencies; above them feynlab attaches the exact
high-frequency tail U/2 + U²n(1−n)/iω. Error bars come from independent Markov chains (`n_chains`).

## Choosing a method (rule of thumb)
* U/W < 0.3: 2OPT/RPA on the lattice — cheap, k-dependent.
* 0.3 < U/W < 1.5: DMFT; IPT for exploration, CT-INT for numbers you publish.
* U/W > 1.5: Mott physics; Hubbard-I for qualitative band positions, CT-INT for quantitative.
* Near a magnetic instability (Uχ₀ → 1): RPA/FLEX diverges — that is physics (Stoner), not a bug.

## v0.2 additions

### Exact-diagonalization DMFT (Caffarel–Krauth) — the accurate Mott transition
```python
from feynlab.dmft.bethe import BetheDMFT
d = BetheDMFT(beta=200, U=2.6)                  # beta = fictitious grid; the cluster itself is at T=0
r = d.solve('ed', n_iter=100, mix=0.5, ctint_kwargs=dict(n_bath=5))
print(r['Z'])                                    # ~0.08
```
Sector-resolved ED (N↑, N↓ blocks), a multi-start bath fit that enforces the exact sum rule Σ V_l² = ∫ε²ρ(ε),
and bounded bath levels. U_c2 ≈ 3.0D at T = 0 (consensus 2.92–3.0D), much closer than IPT (≈3.3D).

### Any lattice
```python
from feynlab.dmft.bethe import LatticeDMFT
LatticeDMFT.square(beta=20, U=2.0).solve('ipt')          # 2D square lattice (t = 1/4 -> W = 2)
LatticeDMFT.cubic(beta=20, U=2.0, nk=32).solve('ipt')    # simple cubic
LatticeDMFT.hypercubic(beta=20, U=2.0).solve('ipt')      # Gaussian DOS
LatticeDMFT.from_eps_k(beta=20, U=2.0, eps_k=my_dispersion)
```

### Maximum-entropy continuation
```python
import numpy as np
from feynlab.core.maxent import maxent
d = BetheDMFT(beta=50, U=2.0); r = d.solve('ipt')
w = np.linspace(-4, 4, 301)
A = maxent(d.grid.iw[:150], r['G'][:150], 1e-4, w)['A']        # spectral function A(w)
```
