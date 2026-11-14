"""HDF5 storage and checkpoint/restart (h5py), with an .npz fallback."""
from __future__ import annotations

import json
import os
import time

import numpy as np

try:
    import h5py
except ImportError:  # pragma: no cover
    h5py = None


def _write(g, key, val):
    if isinstance(val, dict):
        sub = g.require_group(key)
        for k, v in val.items():
            _write(sub, k, v)
    elif isinstance(val, (str, bytes)):
        g.attrs[key] = val
    elif val is None:
        g.attrs[key] = "None"
    else:
        arr = np.asarray(val)
        if key in g:
            del g[key]
        g.create_dataset(key, data=arr)


def _read(g):
    out = {k: (v if not isinstance(v, bytes) else v.decode()) for k, v in g.attrs.items()}
    for k, v in g.items():
        out[k] = _read(v) if isinstance(v, h5py.Group) else v[()]
    return out


def save(path, data: dict, group: str = "/"):
    if h5py is None or not path.endswith((".h5", ".hdf5")):
        np.savez(path, **{k: np.asarray(v) for k, v in data.items()})
        return
    with h5py.File(path, "a") as f:
        g = f.require_group(group)
        for k, v in data.items():
            _write(g, k, v)


def load(path, group: str = "/"):
    if h5py is None or not path.endswith((".h5", ".hdf5")):
        return dict(np.load(path, allow_pickle=False))
    with h5py.File(path, "r") as f:
        return _read(f[group])


class Checkpoint:
    """Iteration checkpoints: ck.save(it, sigma=..., mu=...); ck.latest() -> (it, dict) or None."""

    def __init__(self, path, meta: dict | None = None):
        self.path = path
        if meta and not os.path.exists(path):
            save(path, {"meta_json": json.dumps(meta, default=str)})

    def save(self, iteration: int, **arrays):
        save(self.path, {f"it_{iteration:05d}": dict(arrays, time=time.time())})

    def latest(self):
        if not os.path.exists(self.path):
            return None
        d = load(self.path)
        its = sorted(k for k in d if k.startswith("it_"))
        if not its:
            return None
        return int(its[-1][3:]), d[its[-1]]
