"""
R9 (Track 4): instrumented bracketing scan across the two
persistent-detachment onsets found in §11 of r9_multisubband_methods.md.

Section 11's investigation ruled out two specific numerical explanations for
the M=1/M=4 persistent detachment tails (large-B; search-bracket-too-narrow)
using only the coarse fields already recorded by the original sweep (mu, B,
phi_star, bare_gap). It could not, from that data alone, distinguish two
remaining possibilities for what happens to the principal-branch fixed point
at the onset:

  (i)  a genuine absence of any sign change in f(phi) = G(phi) - phi across
       the scanned bracket (no root exists nearby at all), or
  (ii) a sign change IS found and bisected, but the converged crossing has
       |f(phi*)| >= 1e-4 and is REJECTED by the discontinuity-rejection rule
       already in `_find_fixed_points_ms` — i.e. the same "sawtooth map jumps
       across the identity line without crossing it" mechanism R4/R5/R7
       documented at a gap-closing (`r4_pfaffian_boundary_methods.md`), just
       without an accompanying spectral gap closure here.

This module is a DIAGNOSTIC-ONLY instrumented port of `_find_fixed_points_ms`
/ `_phi_star_ms` from `bdg_r9_multisubband_sweep.py` that also returns which
of (i)/(ii) occurred, plus a DENSE (~9x the original grid) mu-scan bracketing
each of the two identified onsets:

  M=1: mu in [21.40, 34.0]  (original: last finite 21.40, first persistent
       detach 30.91, original step 9.51)
  M=4: mu in [515.70, 545.0] (original: last finite 515.70, first persistent
       detach 538.12, original step 22.42)

Does NOT alter the pre-registered protocol, window, delta, or statistic, and
does NOT touch `r9_multisubband_sweep_data.json` or the reported s-vs-M
result. This is new, additive diagnostic data answering the open question
Section 11.4 flagged, using the exact same physics (`scalar_ring_map_ms`,
`protocol_b_field`, same tol/eps/discontinuity-rejection constant) as the
original sweep — only the instrumentation and the mu-grid density are new.
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
    from bdg_r7_inductance_sweep import A0, _phi_ext_bond           # noqa: E402
    from bdg_r5_selfconsistent_ring import _RING_KERNEL_CACHE       # noqa: E402
    import bdg_r9_multisubband as ms                                # noqa: E402
    from bdg_r9_multisubband_sweep import (                         # noqa: E402
        scalar_ring_map_ms, bare_gap_ms,
    )

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "multisubband")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════════
# Instrumented root search — identical bisection/discontinuity-rejection
# logic to `_find_fixed_points_ms`, plus diagnostics on WHY a bracket failed.
# ════════════════════════════════════════════════════════════════════════
def _find_fixed_points_diag(N, mu, B, g, phi_ext, a, M, span=None, n_scan=49,
                            tol=1e-9, delta=1.0, W_nm=ms.WIRE_W_NM,
                            couple=True) -> tuple[list, dict]:
    """Identical to `_find_fixed_points_ms` (same span default, same n_scan,
    same |f(m)|<1e-4 discontinuity-rejection constant), but also returns a
    diagnostic dict: how many sign changes were found in the scan grid, and
    for each one, whether it was accepted as a root or rejected by the
    discontinuity rule (with the rejected |f(m)| value recorded)."""
    if span is None:
        span = 1.5 * 2.0 * np.pi / N

    def f(p):
        return scalar_ring_map_ms(N, mu, B, g, phi_ext, p, M, delta=delta,
                                  lattice_const=a, W_nm=W_nm, couple=couple) - p

    grid = np.linspace(-span, span, n_scan)
    vals = np.array([f(p) for p in grid])
    roots = []
    diag = {"n_sign_changes": 0, "n_accepted": 0, "n_rejected_discontinuity": 0,
            "rejected_abs_f": []}
    for i in np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]:
        diag["n_sign_changes"] += 1
        lo, hi, flo = float(grid[i]), float(grid[i + 1]), vals[i]
        while hi - lo > tol:
            m = 0.5 * (lo + hi)
            fm = f(m)
            if np.sign(fm) == np.sign(flo):
                lo, flo = m, fm
            else:
                hi = m
        m = 0.5 * (lo + hi)
        fm_final = f(m)
        if abs(fm_final) < 1e-4:
            diag["n_accepted"] += 1
            eps = 1e-6
            slope = (scalar_ring_map_ms(N, mu, B, g, phi_ext, m + eps, M, delta=delta,
                                        lattice_const=a, W_nm=W_nm, couple=couple)
                     - scalar_ring_map_ms(N, mu, B, g, phi_ext, m - eps, M, delta=delta,
                                         lattice_const=a, W_nm=W_nm, couple=couple)
                     ) / (2 * eps)
            roots.append({"phi": m, "slope": float(slope),
                          "picard_stable": bool(abs(slope) < 1.0)})
        else:
            diag["n_rejected_discontinuity"] += 1
            diag["rejected_abs_f"].append(float(abs(fm_final)))
    roots.sort(key=lambda r: abs(r["phi"]))
    return roots, diag


def _phi_star_diag(N, mu, B, g, phi_ext, a, M, hint=None, tol=1e-9, delta=1.0,
                   W_nm=ms.WIRE_W_NM, couple=True) -> tuple:
    """Diagnostic port of `_phi_star_ms`: identical warm-hint-then-fresh-then-
    widened cascade, but returns (phi, slope, diag) where diag merges the
    fresh+widened scan diagnostics (the warm-hint attempt, when it succeeds,
    needs no diagnosis — it directly found and accepted a bracket)."""
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
            fm_final = f(m)
            if abs(fm_final) < 1e-4:
                eps = 1e-6
                slope = (scalar_ring_map_ms(N, mu, B, g, phi_ext, m + eps, M, delta=delta,
                                            lattice_const=a, W_nm=W_nm, couple=couple)
                         - scalar_ring_map_ms(N, mu, B, g, phi_ext, m - eps, M, delta=delta,
                                             lattice_const=a, W_nm=W_nm, couple=couple)
                         ) / (2 * eps)
                return m, float(slope), {"path": "warm_hint", "n_sign_changes": 1,
                                         "n_accepted": 1, "n_rejected_discontinuity": 0,
                                         "rejected_abs_f": []}
            # warm-hint bracket DID contain a sign change but it was rejected —
            # itself diagnostic; fall through to the fresh/widened scan below
            # but keep this fact visible.
            warm_hint_rejected = float(abs(fm_final))
        else:
            warm_hint_rejected = None
    else:
        warm_hint_rejected = None

    roots, diag1 = _find_fixed_points_diag(N, mu, B, g, phi_ext, a, M, tol=tol,
                                           delta=delta, W_nm=W_nm, couple=couple)
    if not roots:
        roots, diag2 = _find_fixed_points_diag(N, mu, B, g, phi_ext, a, M,
                                                span=3.0 * period, n_scan=121,
                                                tol=tol, delta=delta, W_nm=W_nm,
                                                couple=couple)
    else:
        diag2 = {"n_sign_changes": 0, "n_accepted": 0,
                 "n_rejected_discontinuity": 0, "rejected_abs_f": []}
    merged = {
        "path": "fresh_then_widened" if hint is None else "hint_failed_then_fresh",
        "warm_hint_rejected_abs_f": warm_hint_rejected,
        "n_sign_changes": diag1["n_sign_changes"] + diag2["n_sign_changes"],
        "n_accepted": diag1["n_accepted"] + diag2["n_accepted"],
        "n_rejected_discontinuity": (diag1["n_rejected_discontinuity"]
                                     + diag2["n_rejected_discontinuity"]),
        "rejected_abs_f": diag1["rejected_abs_f"] + diag2["rejected_abs_f"],
        "fresh_n_sign_changes": diag1["n_sign_changes"],
        "widened_n_sign_changes": diag2["n_sign_changes"],
    }
    if not roots:
        return float("nan"), float("nan"), merged
    return roots[0]["phi"], roots[0]["slope"], merged


# ════════════════════════════════════════════════════════════════════════
# Bracket scan driver
# ════════════════════════════════════════════════════════════════════════
def run_bracket_scan(M: int, mu_lo: float, mu_hi: float, n_points: int,
                     hint0: float, N: int = 80, a: float = A0,
                     g: float = ms.G_PROBE, W_nm: float = ms.WIRE_W_NM,
                     delta: float = 1.0, couple: bool = True) -> list:
    """Dense warm-hint-continued mu-scan across [mu_lo, mu_hi], seeded with
    `hint0` — the ALREADY-KNOWN phi_star from the last finite point in the
    original (coarser) sweep, so continuity with the original data is
    preserved rather than re-derived. Same field-tracking rule
    (`ms.protocol_b_field`), same self-consistency map, same tol/eps as the
    original sweep — only the mu-grid density and the diagnostics differ."""
    mus = np.linspace(mu_lo, mu_hi, n_points)
    rows = []
    hint = hint0
    t0 = time.time()
    phi_ext = _phi_ext_bond(N)
    for i, mu in enumerate(mus):
        B = ms.protocol_b_field(float(mu), M, ms.OFFSET, W_nm, ms.T_NORM)
        _RING_KERNEL_CACHE.clear()
        p, s, diag = _phi_star_diag(N, float(mu), B, g, phi_ext, a, M,
                                    hint=hint, delta=delta, W_nm=W_nm, couple=couple)
        gap = bare_gap_ms(N, float(mu), B, M, delta=delta, W_nm=W_nm, couple=couple)
        row = {"M": M, "mu": float(mu), "B": float(B),
              "phi_star": float(p) if np.isfinite(p) else float("nan"),
              "s": float(s) if np.isfinite(s) else float("nan"),
              "bare_gap": gap, "diag": diag}
        rows.append(row)
        hint = p if np.isfinite(p) else None
        det = "DETACHED" if not np.isfinite(s) else ""
        print(f"[M={M} bracket {i+1}/{n_points} t={time.time()-t0:6.1f}s] "
              f"mu={mu:9.3f} s={s if np.isfinite(s) else float('nan'):>8.4f} "
              f"nsignchg={diag['n_sign_changes']} nrej={diag['n_rejected_discontinuity']} "
              f"{det}", flush=True)
    return rows


def generate_bracket_data() -> dict:
    t0 = time.time()
    data = {}
    # M=1 bracket: last finite mu=21.40 (phi_star=0.020960...), first
    # persistent-tail detach mu=30.91. Scan a bit past it too.
    print(f"[{time.time()-t0:6.1f}s] M=1 bracket scan (mu 21.40 -> 34.0)...", flush=True)
    data["M1_bracket"] = {
        "mu_lo": 21.40, "mu_hi": 34.0, "n_points": 13, "hint0": 0.020960,
        "rows": run_bracket_scan(1, 21.40, 34.0, 13, hint0=0.020960),
    }
    with open(_out("r9_multisubband_bracket_scan.json"), "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"[{time.time()-t0:6.1f}s] wrote partial (M=1 done); starting M=4...", flush=True)

    # M=4 bracket: last finite mu=515.70 (phi_star=-0.011030...), first
    # persistent-tail detach mu=538.12. Scan a bit past it too.
    print(f"[{time.time()-t0:6.1f}s] M=4 bracket scan (mu 515.70 -> 545.0)...", flush=True)
    data["M4_bracket"] = {
        "mu_lo": 515.70, "mu_hi": 545.0, "n_points": 13, "hint0": -0.011030,
        "rows": run_bracket_scan(4, 515.70, 545.0, 13, hint0=-0.011030),
    }
    with open(_out("r9_multisubband_bracket_scan.json"), "w") as f:
        json.dump(data, f, indent=2, default=str)
    print(f"[{time.time()-t0:6.1f}s] done; wrote r9_multisubband_bracket_scan.json", flush=True)
    return data


if __name__ == "__main__":
    import sys
    with contextlib.suppress(Exception):
        sys.stdout.reconfigure(encoding="utf-8")
    generate_bracket_data()
