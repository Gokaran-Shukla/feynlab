# Tutorial 18 — λ-corrected ladder DΓA

The plain ladder DΓA (`vertex.ladder_dga_sigma`) builds the lattice susceptibility from the local DMFT vertex. It
inherits two defects of DMFT:

* its spin susceptibility violates the local sum rule;
* it diverges at the DMFT Néel temperature, even in 2D, where the Mermin–Wagner theorem forbids order at T > 0.

Below that temperature χ_m(Q) turns negative and the self-energy becomes acausal. Moriya's λ correction fixes both
(Katanin, Toschi, Held, PRB 80, 075104 (2009); Rohringer et al., RMP 90, 025003 (2018), Sec. III.B).

## 1. Equations

Triangular vertices (γ_r = 1 in RPA), with s_d = +1 and s_m = −1:

    Σ_ν' χ_r(q)_{νν'} = χ0(q)_ν γ_r(ν; q) [1 − s_r U χ_r(q, ω)]

Self-energy (identical to the plain ladder at λ = 0):

    Σ(k,ν) = Σ_DMFT(ν) − Σ_EOM^loc(ν) + (U/β) Σ_{qω}/N_q G(k+q, ν+ω)
             × { ½[1 − γ_d(1 − U χ_d)] − 3/2 [1 − γ_m(1 + U χ_m^λ)] } + (local-vertex double-counting term)

λ correction and sum rule:

    χ_m^λ(q,ω) = [χ_m(q,ω)^{-1} + λ_m]^{-1}
    (1/β) Σ_{qω} χ_m^λ(q,ω) = (1/β) Σ_ω χ_m^loc(ω)

The sum rule fixes λ_m, with λ_m > λ_min = −1/χ_m(Q, 0).

* **Finite box.** With `n_tail > 0`, the bubble outside the frequency box is added to every physical susceptibility,
  and λ acts on the corrected χ.
* **Exact target.** `docc=⟨n↑n↓⟩` replaces the box target with (n − 2⟨d⟩)/2. See the note in
  `benchmarks/lambda_dga/RESULTS.md`: with a small box this is inconsistent.

## 2. Usage

```python
from feynlab.vertex import ed_dmft, vertices, lambda_dga
from feynlab.core.lattice import square_lattice_eps
eps = square_lattice_eps(16, 0.25); beta, U = 10.0, 2.0
ed = ed_dmft(beta, U, n_bath=2, eps_k=eps)
gs = ed.generalized_susceptibilities(beta, n_nu=8, n_om=7)
v = vertices(gs, beta)
sig = lambda x: ed.sigma(1j * x, beta, 0)
gs["G_of"] = lambda x: (1 / (1j * x[None, None] + U / 2 - eps[..., None] - sig(x))).mean((0, 1))
r = lambda_dga(eps, U / 2, sig, v, gs, beta, U, 1.0, n_tail=256)
r["lam_m"], r["Sigma"], r["chi_m_lambda"]
```

## 3. Results (2D Hubbard, half filling, t = 0.25)

| | U = 4t, T = 0.4t | U = 8t, T = 0.4t | U = 8t, T = 0.31t |
|---|---|---|---|
| plain ladder χ_m(Q,0) | 6.4 | 538 | −11.2 (acausal Σ) |
| λ-DΓA χ_m(Q,0) | 5.0 | 12.3 | 29.1 |
| λ_m | 0.044 | 0.080 | 0.124 |

* Below the DMFT T_N, λ-DΓA keeps χ_m(Q) finite and growing as T falls: Mermin–Wagner.
* At U = 4t it agrees with TPSC within 25%.
* Full table and figure: `benchmarks/lambda_dga/RESULTS.md`.

## 4. Limits

* Only λ_m is applied. λ_dm (Rohringer & Toschi, PRB 94, 125144 (2016)) and self-consistent λ are not implemented.
* The result depends on the box (see the target note). Larger boxes are the cure.
