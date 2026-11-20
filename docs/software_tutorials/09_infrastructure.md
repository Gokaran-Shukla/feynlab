# Tutorial 9 — Robust infrastructure: errors, mixing, checkpoints, CLI, GPU, units

<!-- gs-watermark -->
> **© 2026 Gokaran Shukla** · part of [feynlab](https://github.com/Gokaran-Shukla/feynlab) · MIT License — any reuse must keep this copyright notice and credit the author.


Roadmap items D17–D22. These are the parts that make a long QMC or DMFT run trustworthy and repeatable.

## 1. Error bars for correlated Monte Carlo data (`feynlab.utils.errors`)

Successive Monte Carlo measurements are autocorrelated. The naive error σ/√N underestimates the true error by
√(2 τ_int), where τ_int is the integrated autocorrelation time.

| tool | what it does | reference |
|---|---|---|
| `binning_analysis(x)` | error vs bin size; the plateau is the true error | Flyvbjerg & Petersen, JCP 91, 461 (1989) |
| `gamma_method(x, S=1.5)` | τ_int with Wolff's automatic window W | Wolff, CPC 156, 143 (2004) |
| `jackknife(f, *series)` | error of a *derived* quantity f(⟨a⟩, ⟨b⟩, …) | Efron (1982) |
| `ar1_series(n, φ)` | test process with exact τ_int = (1+φ)/(2(1−φ)) | — |

```python
from feynlab.utils.errors import ar1_series, gamma_method, jackknife
x = ar1_series(400000, 0.9)
mean, err, tau, W = gamma_method(x)        # tau ~ 9.0-9.5 (exact 9.5), W = summation window
ratio, dratio = jackknife(lambda a, b: a / b, series_a, series_b)
```

**Rule of thumb.** Bin so that each bin is ≥ 10 τ_int long. Report τ_int with every QMC result.

## 2. Self-consistency mixers (`feynlab.utils.mixing`)

A DMFT loop is a fixed-point problem x = F(x). Linear mixing x ← (1−α)x + αF(x) converges slowly near a Mott
transition, where the largest eigenvalue of the Jacobian approaches 1.

* `AndersonMixer(m=5, alpha=0.5)` solves a small least-squares problem over the last m residuals
  (Anderson, J. ACM 12, 547 (1965)). This is the same as DIIS/Pulay.
* `BroydenMixer` uses a quasi-Newton secant update of the inverse Jacobian (Johnson, PRB 38, 12807 (1988)).

```python
from feynlab.dmft.bethe import BetheDMFT
for m in ("linear", "anderson", "broyden"):
    r = BetheDMFT(100, 2.6).solve("ipt", mixer=m, mix=0.5, tol=1e-9, n_iter=500)
    print(m, r["iterations"], r["Z"])       # 78 / 13 / 12 iterations; identical Z to 1e-9
```

## 3. HDF5 checkpoints (`feynlab.utils.io`)

```python
r = BetheDMFT(50, 2.0).solve("ipt", n_iter=200, checkpoint="run.h5")  # restarts from the last iteration if run.h5 exists
from feynlab.utils.io import Checkpoint, save, load
ck = Checkpoint("qmc.h5", {"U": 2.0}); ck.save(it, sigma=sigma); it, data = ck.latest()
```
Every iteration's Σ(iω), the input parameters and the iteration index are stored. Killed jobs therefore lose
at most one iteration.

## 4. Command-line interface with TOML inputs

```bash
feynlab example dmft > dmft.toml     # also: hall, polaron, nrg
feynlab dmft dmft.toml               # results -> HDF5 file named in [output]
feynlab hall --model dirac           # quick run with the example input
feynlab nrg                          # Kondo thermodynamics with the example input
feynlab bench --suite v1             # core | v1 | all
```
The full input is saved in the HDF5 file as `input_json`, so each result file carries its own provenance.

## 5. GPU back-end and units

`feynlab.utils.backend` exposes `xp` (NumPy, or CuPy after `set_backend("cupy")`) and `to_numpy`.
Install it with `pip install feynlab[gpu]`. `feynlab.utils.units` converts between dimensionless units and
eV, Å, S/cm, and e²/h per layer. For example, `E2H_PER_A_TO_S_PER_CM = 3874.0` converts a 3D Hall conductivity
from (e²/h)/Å to S/cm.

## 6. Continuous integration and citation

`.github/workflows/benchmarks.yml` runs the tests on every push and the full benchmark suite every Monday
(cron `17 3 * * 1`). `CITATION.cff` and `.zenodo.json` make every GitHub release citable with a DOI.
`mkdocs.yml` renders these tutorials as a website (`mkdocs serve`).

Quiz: **T10** in `quizzes/quiz.html`.
