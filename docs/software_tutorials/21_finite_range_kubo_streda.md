# Tutorial 21 — Finite-range disorder in the Kubo–Středa formalism

`DisorderedHall` assumes δ-correlated impurities. With that assumption the Born self-energy and the ladder vertex
do not depend on k, and the Bethe–Salpeter equation is an nb² × nb² linear system. Real disorder has a range d:
remote dopants, charged impurities and ripples. The impurity line then carries a form factor:

    W(q) = n_i V0² |u(q)|²,   Gaussian: u(q) = e^{−q²d²/4}

    Σ^R(k)  = ∫_k' W(k−k') M G^R(k') M†                      (self-consistent Born, anti-Hermitian part)
    Γ_y(k)  = v_y(k) + ∫_k' W(k−k') M G^R(k') Γ_y(k') G^A(k') M†

Both are convolutions in k. `transport.finite_range.FiniteRangeHall` evaluates them on a Cartesian grid by
zero-padded FFTs (exact linear convolution on the truncated domain) and solves the vertex equation with GMRES.
The rest follows `DisorderedHall`:

    σ^I = ½ ∫ Tr[2 J_x G^R Γ_y G^A − J_x G^R v_y G^R − J_x G^A v_y G^A]
    σ^II = intrinsic − σ^I(clean)

```python
from feynlab.transport.models import MassiveDirac
from feynlab.transport.finite_range import FiniteRangeHall
fr = FiniteRangeHall(MassiveDirac(1.0, 0.5), mu=1.0, n_i=1.0, V0=0.1, d=1.0, K=2.0, N=1024)
fr.sigma_xx(), fr.decompose()["extrinsic_gaussian"]
```

**Grid.** The spacing 2K/N should resolve the Fermi surface, roughly h ≲ Γ/v_F to a few Γ/v_F. Check convergence
in N.

**Validation** (`benchmarks/finite_range_kubo/RESULTS.md`, massive Dirac model):
* σ_xx agrees with the analytic Boltzmann (Bessel-function) result within 0.2–0.4% for k_F d = 0 to 2.6.
* The extrinsic AHE agrees with the semiclassical side jump + anomalous distribution + Gaussian skew within 0.6%.
* d = 0 reproduces `DisorderedHall` to 1e-9.

**Physics.** Smooth disorder kills the Gaussian (intrinsic) skew much faster than the side jump.

**Limits.** The third-order (V1³) skew vertex and the crossed X/Ψ diagrams are not implemented for finite range.
The model must be 2D.
