"""
ISEF 2027 -- R9 (Track 4 extension): multi-subband loop-gain probe at M=5,6
============================================================================
Extends the pre-registered protocol of `bdg_r9_multisubband.py` /
`bdg_r9_multisubband_sweep.py` (M_MAX=4, "bounded-tractable ceiling, scoping
Sec 4.5") to two more occupied subbands. This is NOT a new protocol: every
formula used here -- `subband_offsets`, `van_hove_delta`, `protocol_b_window`,
`protocol_b_field`, `preregistered_statistic` -- is the identical closed-form,
locked-blind-to-s function already committed for M=1..4, evaluated at M=5,6.
No new judgment call is made; the window/exclusion/statistic protocol is a
pure function of M with no M-dependent tuning, so extending its domain is not
a post-hoc choice of what makes the answer come out favorably.

Discipline mirrored from the original commit sequence:
  (1) gates first -- PH symmetry, decoupled-subband limit, and the
      inter-subband selection rule, run at M=5 and M=6 specifically (not just
      M<=4), because a construction bug could in principle only appear at
      larger M (e.g. an index off-by-one that happens to cancel at small M).
  (2) the window/exclusion parameters for M=5,6 are computed and written out
      BEFORE any s(mu) is computed (this script prints and saves them in a
      block clearly separated from, and prior to, the sweep call).
  (3) the sweep itself, then the unmodified `ms.preregistered_statistic`
      reduction -- no new statistic, no new exclusion rule.

Only M=5,6 are computed here (M=1..4 already exist, gated and reported in
`r9_multisubband_sweep_data.json`); this avoids ~35 min of redundant re-
compute at M=4 while adding the two new points requested.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import time

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    from bdg_r9_multisubband import (                                   # noqa: E402
        build_ring_lo_multisubband, check_ph_symmetry_ms, subband_offsets,
        subband_epsilon1, subband_spacing, van_hove_delta, protocol_b_window,
        protocol_b_field, preregistered_statistic, validate_decoupled_limit,
        validate_intersubband_selection_rule, WIRE_W_NM, EPS_TOL, T_NORM,
        G_PROBE,
    )
    from bdg_r9_multisubband_sweep import (                             # noqa: E402
        _chain_worker_b, _detect_lobe_crossings,
    )
    from bdg_r7_inductance_sweep import A0                              # noqa: E402
    from bdg_parallel import parallel_map, default_workers              # noqa: E402

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "multisubband")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


M_NEW = [5, 6]
N = 80
N_MU_BY_M = {5: 12, 6: 10}   # asymmetric count, same compute-driven-reduction
                             # practice as DEFAULT_N_MU_BY_M for M<=4


# ════════════════════════════════════════════════════════════════════════
# STEP 1 -- construction gates at M=5,6, BEFORE any s(mu) is computed
# ════════════════════════════════════════════════════════════════════════
def run_gates(M_list=M_NEW) -> dict:
    out = {}
    worst_ph = 0.0
    rng = np.random.default_rng(5)
    for M in M_list:
        mu = float(rng.uniform(0.0, subband_offsets(M)[-1] + 2.0))
        B = float(rng.uniform(0.3, 1.5))
        phi = rng.uniform(-0.05, 0.05, 30)
        H = build_ring_lo_multisubband(30, mu, B, phi, M=M)
        eps = check_ph_symmetry_ms(H, 30, M)
        out[f"ph_M={M}"] = eps
        worst_ph = max(worst_ph, eps)

        dec = validate_decoupled_limit(N=30, M=M)
        out[f"decoupled_M={M}"] = dec

        sel = validate_intersubband_selection_rule(M=M)
        out[f"selection_rule_M={M}"] = sel

    out["worst_eps_PH"] = worst_ph
    out["all_pass"] = bool(
        worst_ph < 1e-12
        and all(out[f"decoupled_M={M}"]["pass"] for M in M_list)
        and all(out[f"selection_rule_M={M}"]["pass"] for M in M_list)
    )
    return out


# ════════════════════════════════════════════════════════════════════════
# STEP 2 -- window/exclusion parameters, computed and recorded BEFORE the
# sweep is run (pre-registration-extension record, not a post-hoc summary)
# ════════════════════════════════════════════════════════════════════════
def build_preregistration_extension(M_list=M_NEW, n_mu_by_M=N_MU_BY_M) -> dict:
    e1 = subband_epsilon1()
    delta = van_hove_delta()
    spec = {"e1": e1, "van_hove_delta": delta, "n_mu_by_M": n_mu_by_M,
            "windows": {}}
    for M in M_list:
        thr = subband_offsets(M)[-1]
        span = subband_spacing(M)
        mu_lo = thr - span / 3.0
        mu_hi = thr + span
        spec["windows"][f"M={M}"] = {
            "threshold_dEps_M": float(thr), "spacing_dEsb": float(span),
            "mu_min": float(mu_lo), "mu_max": float(mu_hi),
            "delta_over_window_width": float(2 * delta / (mu_hi - mu_lo)),
        }
    return spec


# ════════════════════════════════════════════════════════════════════════
# STEP 3 -- the sweep itself (identical machinery to run_protocol_b, M=5,6 only)
# ════════════════════════════════════════════════════════════════════════
def run_sweep(M_list=M_NEW, n_mu_by_M=N_MU_BY_M, N=N, a=A0, g=G_PROBE,
             n_workers=2) -> dict:
    jobs = [(N, a, g, M, n_mu_by_M[M], WIRE_W_NM, 1.0, True) for M in M_list]
    t0 = time.time()
    chain_results = parallel_map(_chain_worker_b, jobs,
                                 n_workers=min(n_workers, len(M_list)))
    elapsed = time.time() - t0

    results = {}
    for i, M in enumerate(M_list):
        rs = sorted(chain_results[i], key=lambda r: r["mu"])
        mus = np.array([r["mu"] for r in rs])
        s_vals = np.array([r["s"] for r in rs])
        gaps = np.array([r["bare_gap"] for r in rs])
        lobes = _detect_lobe_crossings(mus, gaps)
        stat = preregistered_statistic(mus, s_vals, M, extra_exclusions=lobes)
        results[f"M={M}"] = {"rows": rs, "lobe_crossings_detected": lobes,
                             "statistic": stat}
    return {"N": N, "a_nm": a * 1e9, "g": g, "n_mu_by_M": n_mu_by_M,
            "M_list": M_list, "n_workers_used": n_workers,
            "elapsed_s": elapsed, "by_M": results}


# ════════════════════════════════════════════════════════════════════════
# DRIVER
# ════════════════════════════════════════════════════════════════════════
def generate_data() -> dict:
    t0 = time.time()
    print(f"[{time.time()-t0:6.1f}s] STEP 1 -- construction gates at M=5,6...",
          flush=True)
    gates = run_gates()
    print(f"[{time.time()-t0:6.1f}s]   all_pass={gates['all_pass']}  "
          f"worst_eps_PH={gates['worst_eps_PH']:.3e}", flush=True)
    if not gates["all_pass"]:
        raise RuntimeError(
            "Construction gate FAILED at M=5 or M=6 -- refusing to compute "
            "s(mu;M) until this is fixed. See gates dict for detail."
        )

    print(f"[{time.time()-t0:6.1f}s] STEP 2 -- pre-registration extension "
          f"(window/exclusion, blind to s)...", flush=True)
    prereg = build_preregistration_extension()
    with open(_out("r9_multisubband_preregistration_M5_M6.json"), "w") as f:
        json.dump({"gates": gates, "prereg_extension": prereg}, f, indent=2,
                  default=str)
    for M in M_NEW:
        w = prereg["windows"][f"M={M}"]
        print(f"[{time.time()-t0:6.1f}s]   M={M}: window=[{w['mu_min']:.1f}, "
              f"{w['mu_max']:.1f}]  delta/width={w['delta_over_window_width']:.3f}",
              flush=True)

    print(f"[{time.time()-t0:6.1f}s] STEP 3 -- sweep s(mu;M) at M=5,6 "
          f"(n_mu={N_MU_BY_M})...", flush=True)
    sweep = run_sweep()
    for M in M_NEW:
        st = sweep["by_M"][f"M={M}"]["statistic"]
        print(f"[{time.time()-t0:6.1f}s]   M={M}: median_s_excluded="
              f"{st['median_s_excluded']:.4f}  median_s_raw={st['median_s_raw']:.4f}  "
              f"n_surviving={st['n_surviving']}/{st['n_raw_finite']}  "
              f"lobe_crossings={sweep['by_M'][f'M={M}']['lobe_crossings_detected']}",
              flush=True)

    data = {"gates": gates, "prereg_extension": prereg, "sweep": sweep}
    with open(_out("r9_multisubband_sweep_M5_M6.json"), "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"[{time.time()-t0:6.1f}s] done; wrote r9_multisubband_sweep_M5_M6.json",
          flush=True)
    return data


if __name__ == "__main__":
    import sys
    with contextlib.suppress(Exception):
        sys.stdout.reconfigure(encoding="utf-8")
    generate_data()
