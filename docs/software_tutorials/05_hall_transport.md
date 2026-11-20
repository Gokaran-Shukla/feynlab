# Tutorial 5 — Hall transport: intrinsic, side jump and skew scattering

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.

### for anomalous (AHC), spin (SHC), orbital (OHC) and orbital-quadrupole Hall conductivities

## 1. Models

```python
from feynlab.transport.models import MassiveDirac, Rashba2DEG, QWZ, Haldane, POrbitalSquare, WannierTB
d   = MassiveDirac(v=1.0, Delta=0.5)                  # 2D massive Dirac, Sinitsyn et al. 2007
r   = Rashba2DEG(m=1.0, alpha=0.3, Delta_ex=0.0)      # spin Hall (Delta_ex=0) or AHE (Delta_ex!=0)
p   = POrbitalSquare(lam=0.2, gamma=0.3, m_ex=0.0)    # px,py,pz x spin: OHC, SHC, quadrupole
# tb = WannierTB("seed_hr.dat", lat, operators={"Lz": Lz_matrix})   # your own Wannier90 model (kz=0)
```
Every model exposes `H(kx,ky)`, `dH(kx,ky)` and a dict `operators` (`'sz/2'`, `'Lz'`, `'Sz'`, `'Jz'`,
`'Qxy'`, `'Qx2-y2'`, `'Q3z2-r2'`, …). Write your own by subclassing `BandModel`.

## 2. Intrinsic responses (clean limit)

σ^O_xy = −e Σₙ∫ fₙ Ωₙ^O,  Ωₙ^O = −2 Im Σ_{m≠n} ⟨n|J^O_x|m⟩⟨m|v_y|n⟩/(εₙ−εₘ)², J^O_x = {O, v_x}/2.

```python
from feynlab.transport.intrinsic import intrinsic_hall, chern_numbers
print(intrinsic_hall(QWZ(-1.0), 0.0))                 # -1  (Chern insulator; = -C, see chern_numbers)
o = p.operators
for name in ("Lz", "Sz", "Qx2-y2"):
    print(name, intrinsic_hall(p, 0.5, op=o[name], nk=200))
```
Physics you should reproduce:
* OHC ≠ 0 even at λ = 0; SHC ∝ λ (orbital Hall is the "parent" of spin Hall).
* With time reversal (m_ex = 0) the charge AHE **and** the orbital-quadrupole Hall response vanish
  identically — {Q, v} is T-odd. Switch on `m_ex` to obtain them.
* Rashba 2DEG: intrinsic SHC = e/8π (0.25 in units e/2π) when both bands are occupied.

## 3. Disorder: the Kubo–Středa engine (any operator)

```python
import numpy as np
from feynlab.transport.kubo_streda import DisorderedHall
from feynlab.transport.intrinsic import k_grid_polar
mu = 1.0; kF = d.kF(mu)
grid = k_grid_polar(8.0, nr=300, nth=48, dense=[(kF, 0.4, 1500), (kF, 0.05, 1500)])
dh = DisorderedHall(d, mu, n_i=1.0, V0=np.sqrt(0.005), V1=0.03, grid=grid)
res = dh.decompose(intrinsic_value=d.ref_intrinsic(mu))
```
Diagrams included (non-crossing approximation): Born self-energy (matrix, anti-Hermitian part),
ladder vertex for the RA bubble (side jump + Gaussian skew), and the third-order vertex
n_i V₁³ [g^R M + M g^A] (conventional skew). `decompose` returns
`intrinsic, sigma_II, extrinsic_gaussian, skew_V1, total`, with
total = σᴵᴵ + σᴵ(ladder) and σᴵᴵ = intrinsic − σᴵ(clean). Pass `op=` for spin/orbital/quadrupole currents.

**Rashba check** (Inoue, Bauer & Molenkamp 2004): `dh.decompose(op=r.operators['sz/2'])` gives
intrinsic +0.25, extrinsic −0.25, total ≈ 0 — the famous vertex-correction cancellation.

## 4. Mechanism-resolved AHE: the Boltzmann engine

```python
from feynlab.transport.mechanisms import anomalous_hall_mechanisms
out = anomalous_hall_mechanisms(d, mu=1.0, band=1, V1_over_V0=0.3, intrinsic_value=d.ref_intrinsic(1.0))
for k in ("intrinsic", "side_jump", "anomalous_distribution", "skew_gaussian", "skew_V1", "total"):
    print(k, out[k])
```
* **side jump**: coordinate shift δr_{l′l} = A_{l′} − A_l − (∂_k+∂_{k′}) arg⟨u_{l′}|u_l⟩, computed
  gauge-invariantly from ∂H/∂k (no finite differences, no gauge fixing).
* **anomalous distribution**: the Boltzmann response to the side-jump energy shift.
* **skew_V1**: third-order asymmetric rate ∝ V₁³ — scales as 1/n_i (dominates in clean metals).
* **skew_gaussian**: Kubo NCA total − (intrinsic + side jump + adist), extrapolated to Γ→0.

All six numbers reproduce Sinitsyn et al., PRB 75, 045315 (2007), Eqs. 58, 65, 67, 72, B8 to ≲1e-5
(see `benchmarks/figures/ahe_dirac_mechanisms.png`).

## 5. Spin / orbital / quadrupole Hall with disorder

```python
from feynlab.transport.mechanisms import operator_hall_mechanisms
from feynlab.transport.intrinsic import k_grid_lattice
pm = POrbitalSquare(lam=0.2, gamma=0.3)
res = operator_hall_mechanisms(pm, 0.5, pm.operators["Lz"], nV2=0.02, grid=k_grid_lattice(pm, 160))
```
For lattice models use a fine k mesh: the RA bubble needs Δk ≲ Γ/v_F.

## 6. Scaling analysis (what experiments measure)
σ_xy^skew ∝ σ_xx (∝ 1/n_i), while intrinsic + side jump are σ_xx-independent. Run
`anomalous_hall_mechanisms` for several V₀ and plot σ_xy vs σ_xx to reproduce the
"skew → intrinsic/side-jump" crossover.

## 7. Known limitations (be honest in your papers)
* Crossed X and Ψ diagrams (Ado et al., EPL 111, 37004 (2015)) are **not** included; they are of the
  same order as side jump for Gaussian disorder and change the Dirac-model answer.
* Disorder is δ-correlated and orbital-independent; real impurities need the T-matrix of the actual potential.
* The Boltzmann engine needs a single star-shaped Fermi contour per call (one band); sum bands yourself.
* Real parts of the Born self-energy (band renormalisation) are dropped, as in the cited literature.

## 8. Beyond the non-crossing approximation: X and Ψ diagrams (v0.2)

Ado, Dmitriev, Ostrovsky & Titov (EPL 111, 37004 (2015); PRL 117, 046601 (2016)) showed that diagrams with
two *crossed* impurity lines — skew scattering on pairs of close impurities — are of the same order as side
jump. feynlab evaluates them in the real-space form (r = vector between the two impurities):

σ^X = (nV₀²)² ∫d²r Tr[J^x_r G^R_{−r} J^y_r G^A_{−r}],  σ^Ψ = (nV₀²)² ∫d²r Tr[J^x_r G^R_{−r} G^R_r J^y_{−r} + h.c.]

with ladder-dressed vertices J_r = ∫ G^A_p Γ G^R_p e^{ipr}, all by 2D FFTs, then extrapolated to Γ → 0.

```python
from feynlab.transport.models import MassiveDirac
from feynlab.transport.crossed import crossed_hall
d = MassiveDirac(1.0, 0.5)
cr = crossed_hall(d, 1.0, runs=((0.03, 500.0, 1024), (0.015, 900.0, 1536)))   # faster settings
print(cr["X"], cr["Psi"])     # -> 0.49 and 0 (Ado et al. Eq. 17), slightly less converged here
```
Or in one call: `anomalous_hall_mechanisms(d, 1.0, band=1, include_crossed=True)` → adds `X`, `Psi`,
`total_with_crossed` (→ −8εm³/(ε²+3m²)² for the Dirac model: the crossed diagrams nearly cancel the NCA
result in the metallic regime).
