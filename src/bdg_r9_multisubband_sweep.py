"""
ISEF 2027 — R9 (Track 4): multi-subband loop-gain sweep — Tasks 1-3
====================================================================
Runs the two protocols locked, blind to s, in `bdg_r9_multisubband.py`'s
pre-registration commit (`r9_multisubband_methods.md` documents the full
spec). This is the SEPARATE, LATER commit the pre-registration promised: it is
the first place any s(mu;M) value is computed.

  TASK 1 — protocol (a): truncation convergence at fixed mu (cheap sanity
           check on the Hamiltonian construction; NOT the physics answer).
  TASK 2 — protocol (b): genuine added occupancy, s(mu;M) curves, reduced by
           the PRE-REGISTERED median-with-exclusion statistic (unmodified from
           the prior commit).
  TASK 3 — the answer, using ONLY that pre-registered statistic.

The self-consistency map is R5's `scalar_ring_map` (mean self-consistent
Peierls phase from the Biot-Savart self-field of the ring current), generalized
to the 4NM multi-subband Hamiltonian via `build_ring_lo_multisubband` /
`ring_current_density_ms`. The root-finding / warm-hint bisection / finite-
difference slope recipe is copied EXACTLY from R7's `_phi_star_a` (the same
tol, the same discontinuity-rejection rule |f|<1e-4, the same eps=1e-6 central
difference) — only the map itself is new.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import time

import numpy as np
from scipy.linalg import eigh

with contextlib.redirect_stdout(io.StringIO()):
    from bdg_r4_em_coupling import (                                  # noqa: E402
        compute_peierls_from_A, T_HOP_SI,
    )
    from bdg_r5_selfconsistent_ring import (                          # noqa: E402
        ring_biot_savart_self, _RING_KERNEL_CACHE,
    )
    from bdg_r7_inductance_sweep import A0, _phi_ext_bond              # noqa: E402
    from bdg_parallel import (                                        # noqa: E402
        parallel_map, check_equivalence, default_workers, _pin_blas,
    )
    import bdg_r9_multisubband as ms                                   # noqa: E402

# bdg_parallel does not export a BLAS-reset helper; R9-T1 defines its own local
# one (mirrored here) to restore the machine-default (unpinned) BLAS threading
# for the "serial_default_blas" baseline timing.


def _pin_blas_reset() -> None:
    for v in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS",
              "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ.pop(v, None)


N_WORKERS: int | None = None
_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "multisubband")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════════
# 0 — multi-subband scalar self-consistency map and slope (mirrors R5/R7 exactly)
# ════════════════════════════════════════════════════════════════════════
def scalar_ring_map_ms(N: int, mu: float, B: float, g: float,
                       phi_ext_bond: float, phi_uniform: float, M: int,
                       delta: float = 1.0, t_hop: float = T_HOP_SI,
                       lattice_const: float = A0, wire_radius=None,
                       W_nm: float = ms.WIRE_W_NM, couple: bool = True) -> float:
    """Multi-subband generalization of R5's `scalar_ring_map`: one application
    of the self-consistency map G to a uniform texture. Identical structure —
    build H at phi_tot -> eigh -> current -> Biot-Savart self-field -> Peierls
    phase — with `build_ring_lo_multisubband`/`ring_current_density_ms`
    replacing the single-subband builders. At M=1 this is `scalar_ring_map`."""
    phi_tot = np.full(N, float(phi_ext_bond) + float(phi_uniform))
    H = ms.build_ring_lo_multisubband(N, mu, B, phi_tot, M=M, W_nm=W_nm,
                                      delta=delta, couple=couple)
    evals, evecs = eigh(H)
    I_bonds = ms.ring_current_density_ms(evecs, evals, N, M, t_hop, phi_tot)
    A_t = ring_biot_savart_self(I_bonds, N, lattice_const, wire_radius)
    out = g * compute_peierls_from_A(A_t, lattice_const)
    return float(np.mean(out))


def _find_fixed_points_ms(N, mu, B, g, phi_ext, a, M, span=None, n_scan=49,
                          tol=1e-9, delta=1.0, W_nm=ms.WIRE_W_NM, couple=True) -> list:
    """Multi-subband port of R7's `_find_fixed_points_a` — identical bisection
    and discontinuity-rejection logic (|f|<1e-4), only the map argument M-aware."""
    if span is None:
        span = 1.5 * 2.0 * np.pi / N

    def f(p):
        return scalar_ring_map_ms(N, mu, B, g, phi_ext, p, M, delta=delta,
                                  lattice_const=a, W_nm=W_nm, couple=couple) - p

    grid = np.linspace(-span, span, n_scan)
    vals = np.array([f(p) for p in grid])
    roots = []
    for i in np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]:
        lo, hi, flo = float(grid[i]), float(grid[i + 1]), vals[i]
        while hi - lo > tol:
            m = 0.5 * (lo + hi)
            fm = f(m)
            if np.sign(fm) == np.sign(flo):
                lo, flo = m, fm
            else:
                hi = m
        m = 0.5 * (lo + hi)
        if abs(f(m)) < 1e-4:
            eps = 1e-6
            slope = (scalar_ring_map_ms(N, mu, B, g, phi_ext, m + eps, M, delta=delta,
                                        lattice_const=a, W_nm=W_nm, couple=couple)
                     - scalar_ring_map_ms(N, mu, B, g, phi_ext, m - eps, M, delta=delta,
                                         lattice_const=a, W_nm=W_nm, couple=couple)
                     ) / (2 * eps)
            roots.append({"phi": m, "slope": float(slope),
                          "picard_stable": bool(abs(slope) < 1.0)})
    roots.sort(key=lambda r: abs(r["phi"]))
    return roots


def _phi_star_ms(N, mu, B, g, phi_ext, a, M, hint=None, tol=1e-9, delta=1.0,
                 W_nm=ms.WIRE_W_NM, couple=True):
    """Multi-subband port of R7's `_phi_star_a`: principal (smallest-|phi|)
    self-consistent uniform self-phase, warm-hint bracket first, widen-to-3-
    periods fallback. Returns (phi, slope) or (nan, nan) — the detachment
    signature, exactly as at M=1."""
    period = 2.0 * np.pi / N

    def f(p):
        return scalar_ring_map_ms(N, mu, B, g, phi_ext, p, M, delta=delta,
                                  lattice_const=a, W_nm=W_nm, couple=couple) - p

    if hint is not None and np.isfinite(hint):
        w = 0.25 * period
        lo, hi = hint - w, hint + w
        flo, fhi = f(lo), f(hi)
        if np.sign(flo) != np.sign(fhi):
            while hi - lo > tol:
                m = 0.5 * (lo + hi)
                fm = f(m)
                if np.sign(fm) == np.sign(flo):
                    lo, flo = m, fm
                else:
                    hi = m
            m = 0.5 * (lo + hi)
            if abs(f(m)) < 1e-4:
                eps = 1e-6
                slope = (scalar_ring_map_ms(N, mu, B, g, phi_ext, m + eps, M, delta=delta,
                                            lattice_const=a, W_nm=W_nm, couple=couple)
                         - scalar_ring_map_ms(N, mu, B, g, phi_ext, m - eps, M, delta=delta,
                                              lattice_const=a, W_nm=W_nm, couple=couple)
                         ) / (2 * eps)
                return m, float(slope)
    roots = _find_fixed_points_ms(N, mu, B, g, phi_ext, a, M, tol=tol, delta=delta,
                                  W_nm=W_nm, couple=couple)
    if not roots:
        roots = _find_fixed_points_ms(N, mu, B, g, phi_ext, a, M, span=3.0 * period,
                                      n_scan=121, tol=tol, delta=delta, W_nm=W_nm,
                                      couple=couple)
    if not roots:
        return float("nan"), float("nan")
    return roots[0]["phi"], roots[0]["slope"]


def s_at_g_ms(N, mu, a, g, M, B, hint=None, delta=1.0, W_nm=ms.WIRE_W_NM,
             couple=True):
    """One principal-branch (phi*, s) evaluation at fixed (N, mu, B, g, M)."""
    _RING_KERNEL_CACHE.clear()
    phi_ext = _phi_ext_bond(N)
    p, s = _phi_star_ms(N, mu, B, g, phi_ext, a, M, hint=hint, delta=delta,
                        W_nm=W_nm, couple=couple)
    return p, s, (p if np.isfinite(p) else hint)


def bare_gap_ms(N, mu, B, M, delta=1.0, W_nm=ms.WIRE_W_NM, couple=True) -> float:
    """Clean (non-self-consistent, phi_self=0) multi-subband bulk gap at the
    external flux alone — used ONLY to flag candidate lobe-crossing points
    (anomalous gap dips) inside the protocol-(b) window, per the pre-registered
    exclusion rule's second clause. Costs one extra eigh per mu point."""
    phi = np.full(N, _phi_ext_bond(N))
    H = ms.build_ring_lo_multisubband(N, mu, B, phi, M=M, W_nm=W_nm, delta=delta,
                                      couple=couple)
    return float(np.min(np.abs(eigh(H, eigvals_only=True))))


# ════════════════════════════════════════════════════════════════════════
# TASK 1 — protocol (a): truncation convergence at fixed mu (sanity check)
# ════════════════════════════════════════════════════════════════════════
def run_protocol_a(N: int = 80, a: float = A0, g: float = ms.G_PROBE,
                   M_max: int = ms.M_MAX) -> dict:
    """Hold (mu, B) FIXED at the R6/R9-T1 single-channel operating point and
    raise ONLY the basis size M=1..M_max. Subbands 2..M stay unoccupied (mu is
    far below their thresholds) — they enter only virtually via second-order
    inter-subband admixture. Confirms the multi-subband CONSTRUCTION converges
    quickly/smoothly with M, validating it before protocol (b) is trusted. This
    is NOT the sec.6.1 physics answer (no new channel is occupied)."""
    mu, B = ms.protocol_a_operating_point()
    rows = []
    hint = None
    t0 = time.time()
    for M in range(1, M_max + 1):
        p, s, hint = s_at_g_ms(N, mu, a, g, M, B, hint=hint)
        rows.append({"M": M, "mu": mu, "B": B, "phi_star": p, "s": s,
                     "elapsed_s": time.time() - t0})
    s_vals = np.array([r["s"] for r in rows])
    finite = np.isfinite(s_vals)
    max_rel_change = float("nan")
    if finite.sum() >= 2:
        fv = s_vals[finite]
        max_rel_change = float(np.max(np.abs(np.diff(fv)) / np.abs(fv[:-1])))
    return {"N": N, "a_nm": a * 1e9, "g": g, "mu": mu, "B": B, "rows": rows,
            "all_finite": bool(finite.all()),
            "max_rel_change_between_successive_M": max_rel_change,
            "smooth_and_converged": bool(finite.all() and
                                        (np.isnan(max_rel_change) or max_rel_change < 0.5))}


# ════════════════════════════════════════════════════════════════════════
# TASK 2 — protocol (b): genuine added occupancy, s(mu;M)
# ════════════════════════════════════════════════════════════════════════
def _b_point_worker(N, mu, B, a, g, M, W_nm, delta, couple, hint=None):
    """One (mu, M) evaluation, with an OPTIONAL warm hint (bracket-and-bisect,
    ~30 eigh) — falls back to a from-scratch grid scan (~80 eigh) only if the
    hint bracket fails to bound a root. Also records the bare (non-self-
    consistent) gap for lobe-crossing detection. Returns (row_dict, hint_out)."""
    p, s, hint_out = s_at_g_ms(N, mu, a, g, M, B, hint=hint, delta=delta,
                               W_nm=W_nm, couple=couple)
    gap = bare_gap_ms(N, mu, B, M, delta=delta, W_nm=W_nm, couple=couple)
    row = {"M": int(M), "mu": float(mu), "B": float(B),
           "phi_star": float(p) if np.isfinite(p) else float("nan"),
           "s": float(s) if np.isfinite(s) else float("nan"),
           "bare_gap": gap}
    return row, hint_out


def _chain_worker_b(N, a, g, M, n_mu, W_nm, delta, couple) -> list:
    """One full protocol-(b) mu-sweep at fixed M, via WARM-HINT CONTINUATION
    along mu (top-level, picklable — one task per M, dispatched in parallel
    across the M_max independent chains). Mirrors R9 Track 1's `branch_curve`
    (continuation along g); here the continuation variable is mu. Only the
    FIRST point of each chain pays the expensive from-scratch grid-scan cost;
    every subsequent point reuses the previous point's phi* as a warm-hint
    bracket (~30 eigh instead of ~80), matching the scoping report's own
    per-point cost model. Detachments (no root even in the widened fallback)
    reset the hint to None so the NEXT point does not inherit a stale bracket
    from a dead branch."""
    mus = ms.protocol_b_window(M, n_mu, W_nm, ms.T_NORM)
    rows = []
    hint = None
    t0 = time.time()
    for i, mu in enumerate(mus):
        B = ms.protocol_b_field(float(mu), M, ms.OFFSET, W_nm, ms.T_NORM)
        row, hint = _b_point_worker(N, float(mu), B, a, g, M, W_nm, delta,
                                    couple, hint=hint)
        rows.append(row)
        print(f"[chain M={M} {i+1}/{len(mus)} t={time.time()-t0:6.1f}s] "
              f"mu={mu:.2f} s={row['s']:.4f}", flush=True)
    return rows


def _detect_lobe_crossings(mus: np.ndarray, gaps: np.ndarray) -> list:
    """Candidate lobe-crossing mu's: local minima of the bare gap curve falling
    below 30% of the window's median gap — a conservative, blind proxy (no
    full multi-subband phase diagram is built, per this task's own scope
    limits). Returns [] if none found (the expected outcome for protocol (b),
    since channel M is tracked near ITS OWN transition and the deep channels
    stay far from theirs)."""
    if len(gaps) < 3:
        return []
    med = float(np.median(gaps))
    thr = 0.3 * med
    out = []
    for i in range(1, len(gaps) - 1):
        if gaps[i] < thr and gaps[i] <= gaps[i - 1] and gaps[i] <= gaps[i + 1]:
            out.append(float(mus[i]))
    return out


# Per-M point count: measured cost is CUBIC in M (dim=4NM), so a UNIFORM n_mu
# would make the M=4 chain's wall-clock dominate everything else by ~8x. This
# asymmetric count is an EXPLICIT, REPORTED compute-driven reduction (same
# practice as R9-T3's narrowed N range) — it trades resolution at the most
# expensive M for a tractable total wall-clock, not a silent one.
DEFAULT_N_MU_BY_M = {1: 12, 2: 10, 3: 8, 4: 5}


def run_protocol_b(N: int = 80, a: float = A0, g: float = ms.G_PROBE,
                   n_mu: int | dict | None = None, M_max: int = ms.M_MAX,
                   n_workers: int | None = None) -> dict:
    """For each M=1..M_max: sweep the PRE-REGISTERED mu-window
    (`ms.protocol_b_window`) with field `ms.protocol_b_field(mu, M)`, computing
    s(mu) via warm-hint continuation ALONG mu within each M (`_chain_worker_b`
    — cheap after the first point), with the M_max independent chains
    dispatched in PARALLEL (each chain internally serial, since each point's
    hint depends on the previous one). `n_mu` may be a single int (same count
    for every M) or a dict {M: count} (default `DEFAULT_N_MU_BY_M`, an
    ASYMMETRIC count compensating for the M^3 cost growth — reported
    explicitly, not silently). Apply the PRE-REGISTERED exclusion (van Hove
    delta + detected lobe crossings) and the median statistic
    (`ms.preregistered_statistic`) — UNCHANGED from the pre-registration
    commit. Report both the excluded-aware and raw statistics."""
    if n_workers is None:
        n_workers = N_WORKERS
    if n_mu is None:
        n_mu_by_M = dict(DEFAULT_N_MU_BY_M)
    elif isinstance(n_mu, dict):
        n_mu_by_M = n_mu
    else:
        n_mu_by_M = {M: n_mu for M in range(1, M_max + 1)}
    jobs = [(N, a, g, M, n_mu_by_M[M], ms.WIRE_W_NM, 1.0, True)
            for M in range(1, M_max + 1)]
    chain_results = parallel_map(_chain_worker_b, jobs,
                                 n_workers=min(n_workers or default_workers(), M_max))

    by_M = {M: chain_results[M - 1] for M in range(1, M_max + 1)}

    results = {}
    for M in range(1, M_max + 1):
        rs = sorted(by_M[M], key=lambda r: r["mu"])
        mus = np.array([r["mu"] for r in rs])
        s_vals = np.array([r["s"] for r in rs])
        gaps = np.array([r["bare_gap"] for r in rs])
        lobes = _detect_lobe_crossings(mus, gaps)
        stat = ms.preregistered_statistic(mus, s_vals, M, extra_exclusions=lobes)
        results[f"M={M}"] = {"rows": rs, "lobe_crossings_detected": lobes,
                             "statistic": stat}
    return {"N": N, "a_nm": a * 1e9, "g": g, "n_mu_by_M": n_mu_by_M, "M_max": M_max,
            "n_workers_used": n_workers or default_workers(), "by_M": results}


# ════════════════════════════════════════════════════════════════════════
# PARALLELIZATION DISCIPLINE — serial vs parallel vs thread-capped, + gate
# ════════════════════════════════════════════════════════════════════════
def _b_point_worker_row_only(N, mu, B, a, g, M, W_nm, delta, couple) -> dict:
    """Independent (hint=None) single-point evaluation, top-level/picklable —
    the embarrassingly-parallel proxy workload used ONLY for the
    parallelization benchmark/equivalence gate (the real protocol-(b) sweep
    uses the cheaper warm-hint `_chain_worker_b` instead)."""
    row, _ = _b_point_worker(N, mu, B, a, g, M, W_nm, delta, couple, hint=None)
    return row


def parallel_benchmark(N: int = 80, a: float = A0, g: float = ms.G_PROBE,
                       M: int = 1, n_points: int = 4) -> dict:
    """Measure serial (default BLAS) vs serial (1-thread BLAS) vs parallel
    (1-thread BLAS) on an independent-point proxy of the protocol-(b) workload,
    before committing to a worker configuration — the same discipline as
    R9-T1's `parallel_benchmark` and R6/R9's measured-not-assumed cost-model
    culture. Deliberately run at the CHEAPEST M (M=1, dim=320): the question
    this benchmark answers is which process/BLAS-threading CONFIGURATION is
    fastest (a scheduling property of this hardware, already established by
    R7-R9 to be M-independent), not the per-point wall-clock at the production
    M — that would make the benchmark itself as expensive as the real sweep.
    Also runs the bit-identical serial-vs-parallel equivalence gate."""
    mus = ms.protocol_b_window(M, n_points)
    args = [(N, float(mu), ms.protocol_b_field(float(mu), M), a, g, M,
            ms.WIRE_W_NM, 1.0, True) for mu in mus]

    try:
        import psutil
        cores = psutil.cpu_count(logical=False) or os.cpu_count() or 1
        free_gb = psutil.virtual_memory().available / 1e9
    except Exception:
        cores = os.cpu_count() or 1
        free_gb = 4.0
    n_proc = max(2, min(int(cores), int(free_gb / 0.5), n_points))

    _pin_blas_reset()
    t0 = time.time(); _ = [_b_point_worker_row_only(*x) for x in args]; t_serial = time.time() - t0

    _pin_blas(1)
    t0 = time.time(); _ = [_b_point_worker_row_only(*x) for x in args]; t_serial_1thr = time.time() - t0
    _pin_blas_reset()

    t0 = time.time()
    _ = parallel_map(_b_point_worker_row_only, args, n_workers=n_proc, blas_threads=1)
    t_par_1thr = time.time() - t0

    eq = check_equivalence(_b_point_worker_row_only, args, n_workers=n_proc, blas_threads=1)
    _pin_blas_reset()

    times = {"serial_default_blas": t_serial, "serial_1thread_blas": t_serial_1thr,
             "parallel_1thread_blas": t_par_1thr}
    fastest = min(times, key=times.get)
    return {"M_tested": M, "n_points": n_points, "cores_physical": int(cores),
            "free_GB_at_bench": round(float(free_gb), 2), "n_proc_used": n_proc,
            "times_s": times, "fastest": fastest,
            "speedup_vs_serial_default": {k: (t_serial / v) for k, v in times.items()},
            "equivalence": eq}


# ════════════════════════════════════════════════════════════════════════
# DRIVER
# ════════════════════════════════════════════════════════════════════════
def generate_data(n_mu=None) -> dict:
    t0 = time.time()
    data = {}
    print(f"[{time.time()-t0:6.1f}s] parallelization benchmark...", flush=True)
    data["bench"] = parallel_benchmark()
    n_workers = data["bench"]["n_proc_used"] if data["bench"]["fastest"] != "serial_default_blas" else 1
    global N_WORKERS
    N_WORKERS = n_workers
    print(f"[{time.time()-t0:6.1f}s]   fastest={data['bench']['fastest']}, "
          f"equivalence match={data['bench']['equivalence']['match']}, using {n_workers} workers",
          flush=True)

    print(f"[{time.time()-t0:6.1f}s] Task 1 — protocol (a): truncation convergence...", flush=True)
    data["protocol_a"] = run_protocol_a()
    print(f"[{time.time()-t0:6.1f}s]   smooth_and_converged="
          f"{data['protocol_a']['smooth_and_converged']}", flush=True)

    print(f"[{time.time()-t0:6.1f}s] Task 2 — protocol (b): genuine occupancy s(mu;M)...", flush=True)
    data["protocol_b"] = run_protocol_b(n_mu=n_mu, n_workers=n_workers)
    for M in range(1, ms.M_MAX + 1):
        st = data["protocol_b"]["by_M"][f"M={M}"]["statistic"]
        print(f"[{time.time()-t0:6.1f}s]   M={M}: median_s_excluded="
              f"{st['median_s_excluded']:.4f}  median_s_raw={st['median_s_raw']:.4f}  "
              f"n_surviving={st['n_surviving']}/{st['n_raw_finite']}", flush=True)

    with open(_out("r9_multisubband_sweep_data.json"), "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"[{time.time()-t0:6.1f}s] done; wrote r9_multisubband_sweep_data.json", flush=True)
    return data


if __name__ == "__main__":
    import sys
    with contextlib.suppress(Exception):
        sys.stdout.reconfigure(encoding="utf-8")
    generate_data()
