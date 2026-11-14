"""Array back-end selection: NumPy (default) or CuPy (GPU) with identical API.

    from feynlab.utils.backend import xp, to_numpy, set_backend
    set_backend("cupy")   # if available; falls back to numpy with a warning
Used by k-space kernels (intrinsic Hall, FFT convolutions) via `xp`.
"""
from __future__ import annotations

import warnings

import numpy as _np

xp = _np
NAME = "numpy"


def set_backend(name: str = "numpy"):
    global xp, NAME
    if name == "cupy":
        try:
            import cupy as cp  # noqa: F401
            cp.zeros(1)
            xp, NAME = cp, "cupy"
            return NAME
        except Exception as e:  # pragma: no cover
            warnings.warn(f"CuPy unavailable ({e}); using NumPy")
    xp, NAME = _np, "numpy"
    return NAME


def to_numpy(a):
    if NAME == "cupy":  # pragma: no cover
        import cupy as cp
        return cp.asnumpy(a)
    return _np.asarray(a)
