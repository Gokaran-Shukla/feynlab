# Tutorial 6 — Electron gas, superconductivity and the book companions

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


## RPA electron gas (Mattuck Chs. 10.4, 12.3, 13.6)

```python
import numpy as np
from feynlab.egas import rpa
rs = 2.0; kF = rpa.kF_of_rs(rs)
q = np.array([0.02, 0.04, 0.06]) * kF
w = rpa.plasmon_dispersion(rs, q)
print(np.polyfit(q**2, w**2 - rpa.plasma_frequency(rs)**2, 1)[0] / kF**2)   # -> 3/5 (book Ex. 13.6)
print(rpa.rpa_correlation_energy(1.0))     # -0.0788 Ha  (ring diagrams summed)
```
The correlation energy is the sum of ring (bubble) diagrams of Ch. 12: each order of the bare series
diverges logarithmically, the RPA sum is finite and gives e_c → 0.0311 ln r_s − 0.0711 Ha
(Gell-Mann & Brueckner 1957). The benchmarks check the slope (1−ln2)/π² and the constant.

`rpa.lindhard_iw(q, nu, kF)` and `rpa.lindhard_real(q, w, kF)` are the pair-bubble of Ch. 9.8/10.7.

## BCS (Ch. 15)

```python
import numpy as np
from feynlab.egas import bcs
lam = 0.2
print(bcs.gap_T0(lam) / bcs.Tc(lam))       # 1.764 = pi e^{-gamma}
print([bcs.gap_T(T, lam) for T in np.linspace(0.001, bcs.Tc(lam), 5)])
```
The gap equation is the self-consistent Nambu "bubble"; T_c is where the Cooper-channel ladder
(particle–particle K-matrix of Ch. 10) diverges — the breakdown of perturbation theory of Sec. 15.4.

## Book companions

```python
import numpy as np
from feynlab.book import companions as bk
bk.pinball_propagator(P0=0.2, p_scatter=0.3)        # Ch. 2: geometric partial sum
bk.quantum_pinball(np.linspace(-1,1,5), eps0=0.1, V=0.2, n_order=3)   # Ch. 3: truncated vs summed
bk.impurity_self_energy(k=1.0, n_imp=0.01, W0=0.5)  # Ch. 3.5: energy shift and lifetime
bk.geometric_dyson(G0=1/(0.5+0.1j), Sigma=0.2-0.05j, n_terms=6)   # Ch. 10.2 Dyson
bk.friedel_oscillation(np.linspace(2, 8, 50), rs=2.0)   # Ch. 10.7 cos(2kF r)/r^3
bk.kondo_log(np.logspace(-4, -1, 5), J_rho=0.1)      # Ch. 18.2: 1 + 2 J rho ln(D/T), T_K
```
Use them while reading: each reproduces a figure or worked equation of the corresponding chapter
(see the chapter tutorials in `docs/book_tutorials/`).
