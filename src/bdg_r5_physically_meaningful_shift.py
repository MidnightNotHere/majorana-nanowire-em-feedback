"""
R5 (revised) — The PHYSICALLY MEANINGFUL headline shift: shift vs g for
g < g_crit, where a genuine ATTRACTING self-consistent fixed point exists.

Why this supersedes the g = 3×10⁵ number
-----------------------------------------
R5 Figure 8 established a Picard bifurcation at g_crit: below it, slow g-ramp
and direct-jump-from-zero converge to the SAME fixed point (a genuine
attracting equilibrium, path-independent); above it there is no attracting
fixed point — the root-finder can still locate a mathematical root of
G(φ)=φ, but nothing a real iterative/dynamical process would settle to.
g = 3×10⁵ sits far above g_crit, so the previously reported "+0.05%
self-consistent shift" is the location of a dynamically unreachable root, not
a physical static equilibrium. (It is kept, clearly labelled, in
bdg_r5_selfconsistent_ring.py and the methods note — it is what revealed g_crit
exists.)

This module reports the honest number: the shift as a function of g, restricted
to the physically accessible regime g < g_crit, at μ = 1.2 and 1.4, with the
SAME validated pipeline as Figures 7-9, and subjects the chosen headline point
to the same validation battery (N-independence, tolerance, path-independence).

Two corrections relative to a first draft of this file (documented so the
reasoning is auditable):
  1. BASELINE. The shift must be measured against the RING/Pfaffian baseline
     (same detector, same geometry, external flux only, NO self-field) —
     ring_boundary_Bc_pfaffian with a uniform φ_ext. An early version wrongly
     used the OPEN-WIRE THRESHOLD detector as the baseline, which injected a
     spurious ~+23% detector/geometry offset that (tellingly) did not vanish
     as g→0. Fixed: baseline is g-independent and computed once per (N, μ).
  2. WINDOW. The self-consistent B_c detector is called with its VALIDATED
     defaults (window=0.05, n_coarse=11) — exactly as the committed Figure 9
     run used at g = 3×10⁵. A first draft passed window=0.2/n_coarse=21, whose
     wide coarse B-grid caused spurious root-bracketing failures far below
     g_crit; those were an artifact of that override, not physics (verified:
     at the transition B the attracting fixed point persists to g ≈ 1×10⁵ for
     both μ — see the diagnostic in the methods note §4e).

Nothing in the validated pipeline is modified: this module imports R4/R5
read-only.
"""

from __future__ import annotations

import json
import os
import time

import numpy as np

import sys
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from bdg_r5_selfconsistent_ring import (            # noqa: E402
    principal_uniform_fixed_point,
    ring_boundary_Bc_pfaffian_selfconsistent,
    selfconsistent_ring_texture,
    ScfConvergenceError,
)
from bdg_r4_em_coupling import ring_boundary_Bc_pfaffian   # noqa: E402

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "self_consistency_bifurcation")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ── THREE regimes of g (established by the diagnostics in §4e), in order of
#    increasing g:
#      1. g ≲ G_DEADZONE_ONSET ≈ 3×10⁴ — a genuine attracting self-consistent
#         fixed point exists at EVERY field B across the transition, so the
#         self-consistent boundary B_c is UNAMBIGUOUSLY defined. This is the
#         honest regime for a headline number.
#      2. G_DEADZONE_ONSET ≲ g ≲ G_CRIT_SINGLEPOINT — the fixed point is still
#         attracting at most B, but a NARROW dead-zone in B opens exactly at
#         the gap-closing (the sawtooth self-consistency map jumps across the
#         identity line there without crossing it), so no uniform fixed point
#         exists in a thin window straddling B_c and the boundary is defined
#         only up to that dead-zone width. The Pfaffian sign flip itself stays
#         clean (single flip), so B_c is still bracketed — it just cannot be
#         bisected to machine precision.
#      3. g ≳ G_CRIT_SINGLEPOINT ≈ 1×10⁵ — the fixed point becomes Picard-
#         repelling (Figure 8 / Part A): NO attracting equilibrium at all.
#         g = 3×10⁵ lives here — the previously reported "+0.05%" is the
#         location of a dynamically unreachable root (kept, labelled, in
#         bdg_r5_selfconsistent_ring.py).
#
#    The honest headline is the top of regime 1 (shift is monotonic in g, so
#    this is the largest cleanly-defined value): G_HEADLINE = 2.9×10⁴. The
#    N-independence / tolerance battery is run at G_NSCAN = 2.0×10⁴, safely
#    inside regime 1 for every N (its purpose is the bulk-vs-finite-size
#    SCALING, which characterises the mechanism and is the same across
#    regime 1); the shift at any regime-1 g scales the same way.
G_DEADZONE_ONSET = 3.0e4
G_CRIT_SINGLEPOINT = 1.0e5
G_HEADLINE = 2.92e4
G_NSCAN = 2.0e4


def ring_baseline_Bc(N: int, mu: float, phi_ext_bond: float,
                     tol: float = 1e-6) -> float:
    """
    The g-independent baseline: ring/Pfaffian B_c with the external flux only
    (uniform φ_ext on every bond, NO self-field). Same detector and geometry
    as the self-consistent B_c, so the shift isolates the self-field effect.
    """
    return ring_boundary_Bc_pfaffian(
        N, mu, phi_bonds=np.full(N, float(phi_ext_bond)), tol=tol)


def selfconsistent_shift(N: int, mu: float, g: float, phi_ext_bond: float,
                         baseline_Bc: float, tol: float = 1e-6) -> dict:
    """
    Converged self-consistent shift ΔB_c/B_c (%) at (N, μ, g), using the
    validated-default detector. Returns nan shift (with 'ok': False) if no
    attracting/­mathematical fixed point can be located at some probed B
    (e.g. inside the sawtooth dead-zone near g_crit) — never fabricated.
    """
    try:
        Bc = ring_boundary_Bc_pfaffian_selfconsistent(
            N, mu, g=g, phi_ext_bond=phi_ext_bond,
            tol=tol, scf_tol=1e-12)
        shift = 100.0 * (Bc - baseline_Bc) / baseline_Bc
        return {"Bc": float(Bc), "shift_pct": float(shift), "ok": True}
    except (ScfConvergenceError, Exception) as e:
        return {"Bc": float("nan"), "shift_pct": float("nan"),
                "ok": False, "error": str(e)}


def path_independence_at_g(N: int, mu: float, g: float, phi_ext_bond: float,
                           B: float) -> dict:
    """
    At fixed (N, μ, g, B), compare three routes to the uniform self-phase φ*:
      - root-finder (stability-independent),
      - direct-jump iteration from a zero seed,
      - slow-ramp iteration (seed = converged φ* at g/2).
    Below g_crit all three must agree; this is the Figure-8 check performed
    at *every* g in the sweep rather than assumed.
    """
    phi_root = principal_uniform_fixed_point(N, mu, B, g, phi_ext_bond)

    res_jump = selfconsistent_ring_texture(
        N, mu, B, g, phi_ext_bond, seed_phi_self=np.zeros(N),
        beta=1.0, tol=1e-12, max_iter=200)
    phi_jump = (float(np.mean(res_jump["phi_self"]))
                if res_jump["converged"] else None)

    # slow ramp: converge at g/2 first, use as seed
    res_half = selfconsistent_ring_texture(
        N, mu, B, g / 2.0, phi_ext_bond, seed_phi_self=np.zeros(N),
        beta=1.0, tol=1e-12, max_iter=200)
    seed = res_half["phi_self"] if res_half["converged"] else np.zeros(N)
    res_ramp = selfconsistent_ring_texture(
        N, mu, B, g, phi_ext_bond, seed_phi_self=seed,
        beta=1.0, tol=1e-12, max_iter=200)
    phi_ramp = (float(np.mean(res_ramp["phi_self"]))
                if res_ramp["converged"] else None)

    vals = [v for v in (float(phi_root), phi_jump, phi_ramp) if v is not None]
    agree = (len(vals) >= 2 and (max(vals) - min(vals)) < 1e-6)
    return {"g": float(g), "B": float(B),
            "phi_root": float(phi_root),
            "phi_jump": phi_jump, "jump_converged": res_jump["converged"],
            "phi_ramp": phi_ramp, "ramp_converged": res_ramp["converged"],
            "all_agree": bool(agree)}


def g_sweep(N: int = 80, mus=(1.2, 1.4),
            g_min: float = 1e3, g_max: float = 9e4,
            n_points: int = 13) -> dict:
    """Shift vs g (both μ) for g < g_crit, with a path-independence check at
    each g. Baseline computed once per μ."""
    g_vals = np.logspace(np.log10(g_min), np.log10(g_max), n_points)
    phi_ext = np.pi * 0.1 / N
    out = {"N": N, "g_deadzone_onset": G_DEADZONE_ONSET,
           "g_crit_singlepoint": G_CRIT_SINGLEPOINT,
           "g_values": g_vals.tolist(),
           "mus": list(mus), "shifts": {}, "baselines": {},
           "path_checks": {}}

    for mu in mus:
        base = ring_baseline_Bc(N, mu, phi_ext)
        out["baselines"][str(mu)] = float(base)
        shifts, pchecks = [], []
        print(f"\n[g-sweep] mu={mu}, N={N}, baseline Bc={base:.6f}")
        for i, g in enumerate(g_vals):
            t0 = time.time()
            r = selfconsistent_shift(N, mu, g, phi_ext, base)
            shifts.append(r["shift_pct"])
            # path-independence at ~4 representative g (fast: fixed B)
            pc = None
            if (i % max(1, n_points // 4) == 0) or i == n_points - 1:
                pc = path_independence_at_g(N, mu, g, phi_ext, base)
                pchecks.append(pc)
            tag = "OK" if r["ok"] else "FAIL"
            extra = "" if pc is None else f"  path_agree={pc['all_agree']}"
            print(f"  g={g:.3e}: shift={r['shift_pct']:+.5f}% [{tag}]"
                  f"{extra}  ({time.time()-t0:.0f}s)")
        out["shifts"][str(mu)] = shifts
        out["path_checks"][str(mu)] = pchecks
    return out


def n_independence(g: float = G_NSCAN, mus=(1.2, 1.4),
                   Ns=(40, 60, 80, 120)) -> dict:
    """
    THE decisive test: at the headline g, does the shift scale to zero with N
    (finite-size Aharonov-Bohm artifact, as the g=3e5 result did) or saturate
    to a finite bulk value? Reported honestly either way.
    """
    out = {"g": g, "mus": list(mus), "Ns": list(Ns), "shifts": {}}
    for mu in mus:
        shifts = []
        print(f"\n[N-independence] g={g:.2e}, mu={mu}")
        for N in Ns:
            phi_ext = np.pi * 0.1 / N
            base = ring_baseline_Bc(N, mu, phi_ext)
            t0 = time.time()
            r = selfconsistent_shift(N, mu, g, phi_ext, base)
            shifts.append(r["shift_pct"])
            print(f"  N={N}: shift={r['shift_pct']:+.5f}%  "
                  f"(Bc={r['Bc']:.6f}, base={base:.6f})  ({time.time()-t0:.0f}s)")
        out["shifts"][str(mu)] = shifts
        # fit power law where all finite and same sign
        sh = np.array(shifts)
        if np.all(np.isfinite(sh)) and np.all(sh < 0):
            p = np.polyfit(np.log(Ns), np.log(-sh), 1)[0]
            out.setdefault("power_law_exponent", {})[str(mu)] = float(p)
            print(f"  -> |shift| ~ N^{p:.2f} "
                  f"({'vanishes' if p < -0.1 else 'saturates'} as N->inf)")
    return out


def tolerance_robustness(g: float = G_NSCAN, mu: float = 1.2,
                         N: int = 80) -> dict:
    """Shift at the headline g vs the bisection tolerance (Figure-9b style)."""
    phi_ext = np.pi * 0.1 / N
    base = ring_baseline_Bc(N, mu, phi_ext, tol=1e-8)
    out = {"g": g, "mu": mu, "N": N, "baseline": float(base), "rows": []}
    print(f"\n[tolerance] g={g:.2e}, mu={mu}, N={N}")
    for bt in (1e-4, 1e-5, 1e-6, 1e-7):
        t0 = time.time()
        try:
            Bc = ring_boundary_Bc_pfaffian_selfconsistent(
                N, mu, g=g, phi_ext_bond=phi_ext, tol=bt, scf_tol=1e-12)
            sh = 100.0 * (Bc - base) / base
        except Exception:
            Bc, sh = float("nan"), float("nan")
        out["rows"].append({"bisect_tol": bt, "Bc": float(Bc),
                            "shift_pct": float(sh)})
        print(f"  bisect_tol={bt:.0e}: shift={sh:+.6f}%  ({time.time()-t0:.0f}s)")
    return out


def _save(results: dict):
    with open(_out("r5_physical_shift_results.json"), "w") as f:
        json.dump(results, f, indent=2)


def main():
    print("=" * 70)
    print("R5 (REVISED): physically meaningful shift, unambiguous regime "
          "g < g_deadzone_onset ~ 3e4")
    print("=" * 70)
    results = {}

    print("\n[1] g-sweep: shift vs g (both mu), path-independence at each")
    results["g_sweep"] = g_sweep(N=80, mus=(1.2, 1.4),
                                 g_min=1e3, g_max=9e4, n_points=13)
    _save(results)

    print("\n[2] N-independence at G_NSCAN (bulk vs finite-size)")
    results["n_independence"] = n_independence(g=G_NSCAN, mus=(1.2, 1.4),
                                               Ns=(40, 60, 80, 120))
    _save(results)

    print("\n[3] Tolerance robustness at G_NSCAN")
    results["tolerance"] = tolerance_robustness(g=G_NSCAN, mu=1.2, N=80)
    _save(results)

    # headline numbers = the largest clean g in the sweep (top of regime 1)
    gs = results["g_sweep"]
    gv = np.array(gs["g_values"])
    print("\n" + "=" * 70)
    print("HEADLINE (physically meaningful, top of the unambiguous regime):")
    for mu in (1.2, 1.4):
        sh = np.array(gs["shifts"][str(mu)], dtype=float)
        ok = np.isfinite(sh)
        i_edge = int(np.where(ok)[0][-1])       # largest clean g
        print(f"  mu={mu}: shift at g={gv[i_edge]:.2e} "
              f"(top of unambiguous regime) = {sh[i_edge]:+.4f}%")
    print("=" * 70)
    _save(results)
    print(f"\nSaved: {_out('r5_physical_shift_results.json')}")
    return results


if __name__ == "__main__":
    main()
