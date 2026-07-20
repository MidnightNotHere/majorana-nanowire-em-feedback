"""
Shared process-parallel map for the R7/R8 sweep drivers.
=======================================================
The heavy work (scipy.linalg.eigh + matmuls) is BLAS-threaded, but for the
matrix sizes here (4N = 320-960) a single eigh only saturates ~1-2 cores, and
the OUTER sweeps (independent grid points / eta / N) run sequentially — so a
10-core machine sits at ~1.4 cores busy.

This module parallelizes those embarrassingly-parallel outer loops across
PROCESSES while pinning each worker's BLAS to a single thread, so N processes
do not oversubscribe the N physical cores. Net effect: near-linear speedup on
the independent sweeps.

Why processes + 1-thread BLAS (not threads): numpy releases the GIL inside
BLAS, but the point is to run *many independent eigh problems at once*; with
each worker single-threaded there is no thread contention, and 10 workers x 1
thread cleanly fills 10 cores.

Determinism note: single-threaded LAPACK is deterministic, so serial-vs-parallel
results are bit-identical AT THE SAME blas_threads. Capping BLAS to 1 thread can
shift results at the ~1e-13 ULP level versus an old multi-threaded run (a
different reduction order), never at any physically meaningful scale. The
equivalence check (`check_equivalence`) compares like-for-like (both 1-thread).
"""

from __future__ import annotations

import os


def _pin_blas(nthreads: int = 1) -> None:
    """Cap every common BLAS/OpenMP backend to `nthreads`. Set in the parent
    before spawning (children inherit os.environ on spawn) AND in each worker's
    initializer (which runs before the worker imports numpy), so the cap is in
    place no matter how OpenBLAS reads it."""
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[v] = str(int(nthreads))


def default_workers() -> int:
    """Evidence-based worker count for these heavy eigh sweeps. NOT the full
    physical-core count: measured on an i7-13620H (6 P + 4 E cores, 16 GB),
    ~4 single-thread-BLAS workers already give the full ~5.7x speedup on the
    R8-B sweeps — workers beyond the P-cores land on slow E-cores / HT siblings
    and add nothing, while 10 workers risk memory pressure (each worker re-imports
    the R4-R7 stack + holds Kubo matrices, ~0.4 GB). So cap at:
        min( physical_cores, 6, free_RAM / 0.5 GB ).
    Override explicitly via the module N_WORKERS / n_workers args when desired."""
    try:
        import psutil
        cores = psutil.cpu_count(logical=False) or os.cpu_count() or 1
        free_gb = psutil.virtual_memory().available / 1e9
    except Exception:
        cores = os.cpu_count() or 1
        free_gb = 4.0
    by_mem = max(1, int(free_gb / 0.5))           # ~0.5 GB budget per worker
    return max(1, min(int(cores), 6, by_mem))


def _apply(func_and_args):
    func, args = func_and_args
    return func(*args)


def parallel_map(func, args_list, n_workers: int | None = None,
                 blas_threads: int = 1, serial: bool = False) -> list:
    """
    Apply `func(*args)` for every args-tuple in `args_list`, across processes.

    Parameters
    ----------
    func         : top-level (picklable) function.
    args_list    : list of positional-arg TUPLES, one per call.
    n_workers    : process count (default = physical cores, capped at len).
    blas_threads : BLAS threads per worker (1 = no oversubscription).
    serial       : force the in-process serial path (used by the equivalence
                   check and as a fallback).

    Returns results in input order. Requires the calling program's entry point
    to be under `if __name__ == "__main__":` (Windows spawn).
    """
    args_list = list(args_list)
    if serial or n_workers == 1 or len(args_list) <= 1:
        return [func(*a) for a in args_list]

    n = min(int(n_workers) if n_workers else default_workers(), len(args_list))
    if n <= 1:
        return [func(*a) for a in args_list]
    # A FRESH (transient) pool per call. Measured to beat a single persistent
    # pool for this workload: a reused pool degrades over the ~40 heavy tasks a
    # full run submits (per-worker memory/BLAS-cache creep), whereas a fresh
    # small pool per big loop starts clean. Only worth it for loops with many
    # tasks — small loops (<= a few) should stay serial at the call site.
    return _pool_map(func, args_list, n, blas_threads)


def _pool_map(func, args_list, n_procs: int, blas_threads: int = 1) -> list:
    """Always run through a process pool of exactly `n_procs` workers (each
    with `blas_threads` BLAS threads). Used both by parallel_map and by the
    equivalence check (which compares 1-proc vs N-proc pools so BOTH sides use
    single-thread BLAS in fresh processes → truly bit-identical)."""
    _pin_blas(blas_threads)                       # parent env -> children inherit
    from concurrent.futures import ProcessPoolExecutor
    with ProcessPoolExecutor(max_workers=n_procs, initializer=_pin_blas,
                             initargs=(blas_threads,)) as ex:
        return list(ex.map(_apply, [(func, a) for a in args_list]))


# ════════════════════════════════════════════════════════════════════════
# Serial-vs-parallel equivalence check
# ════════════════════════════════════════════════════════════════════════
def _flatten_numbers(obj, out):
    """Collect all finite floats/ints from a nested dict/list result."""
    import math
    if isinstance(obj, bool):
        return
    if isinstance(obj, (int, float)):
        if isinstance(obj, float) and (math.isnan(obj) or math.isinf(obj)):
            return
        out.append(float(obj))
    elif isinstance(obj, dict):
        for v in obj.values():
            _flatten_numbers(v, out)
    elif isinstance(obj, (list, tuple)):
        for v in obj:
            _flatten_numbers(v, out)


def check_equivalence(func, args_list, n_workers=None, blas_threads=1,
                      atol=1e-12, rtol=1e-9) -> dict:
    """Compare a 1-worker pool against an N-worker pool — BOTH in fresh
    processes with single-thread BLAS — so the comparison isolates the
    parallelization itself (bit-identical LAPACK ⇒ max_abs_diff ~0), NOT the
    BLAS thread-count reduction-order difference. (Comparing against the
    already-loaded multi-thread parent would only agree to ~1e-9, which is that
    reduction-order floor, not a parallelization bug.)"""
    args_list = list(args_list)
    n = min(int(n_workers) if n_workers else default_workers(), len(args_list))
    serial = _pool_map(func, args_list, 1, blas_threads)
    par = _pool_map(func, args_list, max(2, n), blas_threads)
    a, b = [], []
    for rs, rp in zip(serial, par):
        _flatten_numbers(rs, a)
        _flatten_numbers(rp, b)
    if len(a) != len(b):
        return {"match": False, "reason": f"different result shapes {len(a)} vs {len(b)}"}
    max_abs = 0.0
    max_rel = 0.0
    for x, y in zip(a, b):
        d = abs(x - y)
        max_abs = max(max_abs, d)
        max_rel = max(max_rel, d / max(abs(x), 1e-300))
    ok = (max_abs <= atol) or (max_rel <= rtol)
    return {"match": bool(ok), "n_values": len(a),
            "max_abs_diff": max_abs, "max_rel_diff": max_rel}
