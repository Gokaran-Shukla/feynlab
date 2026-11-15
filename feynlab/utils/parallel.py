"""Parallel Monte Carlo chains: MPI (mpi4py) when launched under mpirun, otherwise a process pool.

    from feynlab.utils.parallel import map_chains
    results = map_chains(worker, [args_chain0, args_chain1, ...], backend="auto")

* backend="mpi"     : chains are distributed round-robin over MPI ranks; results are gathered on every rank
                      (mpirun -n 8 python script.py).  Requires mpi4py.
* backend="process" : concurrent.futures.ProcessPoolExecutor (fork start method on Linux, so numba kernels
                      compiled in the parent are reused; cache=True kernels are reused everywhere).
* backend="serial"  : plain loop (debugging, deterministic order).
* backend="auto"    : "mpi" if running under MPI with more than one rank, else "process" if more than one
                      CPU and more than one chain, else "serial".
Independent chains are the natural parallelism of Markov-chain QMC (Gull et al. RMP 83, 349, Sec. VI); the
statistical error bars of every solver are computed across chains, so results are identical in distribution
for all backends.
"""
from __future__ import annotations

import os


def _mpi_comm():
    try:
        from mpi4py import MPI  # type: ignore
    except Exception:  # noqa: BLE001
        return None
    comm = MPI.COMM_WORLD
    return comm if comm.Get_size() > 1 else None


def available_backends():
    b = ["serial", "process"]
    if _mpi_comm() is not None:
        b.append("mpi")
    return b


def map_chains(worker, arg_list, backend: str = "auto", n_workers: int | None = None):
    """Evaluate worker(*args) for each args in arg_list; returns the list of results in input order."""
    arg_list = list(arg_list)
    if backend == "auto":
        if _mpi_comm() is not None:
            backend = "mpi"
        elif (os.cpu_count() or 1) > 1 and len(arg_list) > 1:
            backend = "process"
        else:
            backend = "serial"
    if backend == "serial":
        return [worker(*a) for a in arg_list]
    if backend == "mpi":
        comm = _mpi_comm()
        if comm is None:
            raise RuntimeError("backend='mpi' requires mpi4py and mpirun with more than one rank")
        rank, size = comm.Get_rank(), comm.Get_size()
        mine = {i: worker(*a) for i, a in enumerate(arg_list) if i % size == rank}
        allres = comm.allgather(mine)
        out = {}
        for d in allres:
            out.update(d)
        return [out[i] for i in range(len(arg_list))]
    if backend == "process":
        import multiprocessing as mp
        from concurrent.futures import ProcessPoolExecutor
        ctx = mp.get_context("fork") if "fork" in mp.get_all_start_methods() else None
        n = n_workers or min(len(arg_list), os.cpu_count() or 1)
        with ProcessPoolExecutor(max_workers=n, mp_context=ctx) as ex:
            futs = [ex.submit(worker, *a) for a in arg_list]
            return [f.result() for f in futs]
    raise ValueError(backend)
