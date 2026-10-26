"""
R9 (Track 4): the FULL self-consistent boundary calculation
at M=2, not just the cheap loop-gain-slope proxy.

R9 Track 4's headline result (r9_multisubband_methods.md Sec 8-9) was
measured with the CHEAP proxy: the local slope s = dG/dphi of the self-
consistency map AT the pre-registered (mu,B) operating point (~30 eigh
calls). That proxy showed s crosses 1 (Picard-repelling) at M=2, mirroring
exactly what R5 found for the single-channel case above g_crit ~ 8.67e4.

But R5's own headline result was never "s>1" by itself -- it was the
measured SELF-CONSISTENT BOUNDARY SHIFT (Delta B_c/B_c), found via Pfaffian-
sign bisection with the self-consistent self-phase substituted in at every
probed B (`ring_boundary_Bc_pfaffian_selfconsistent`). This module builds
the M-subband generalization of THAT function -- not a new physics engine,
a direct, validated PORT of the exact R5/R7 algorithm onto the ALREADY-
VALIDATED M-subband Hamiltonian (`bdg_r9_multisubband.build_ring_lo_multisubband`,
already gated bit-for-bit at M=1, decoupled-limit, and PH-symmetry in the
Track 4 pre-registration commit) -- and the R9-Track-1-style Anderson-mixing
seed-independence check, to determine whether any self-consistent solution
found is genuinely DYNAMICALLY ATTRACTING or merely mathematically locatable
by root-finding (exactly the distinction R9-T1 drew for R5's own s>1 branches).

Does not modify any R1-R9 code path. Purely additive; every new function
here is validated against its EXACT R1-R9 M=1 counterpart before its M=2
output is trusted, per this project's standing validation discipline.
"""

from __future__ import annotations

import contextlib
import io
import time

import numpy as np
from scipy.linalg import eigh

with contextlib.redirect_stdout(io.StringIO()):
    from bdg_r4_em_coupling import (                                   # noqa: E402
        compute_peierls_from_A, T_HOP_SI, _ring_pfaffian_sign,
        _build_ring_lo_phased,
    )
    from bdg_r5_selfconsistent_ring import (                           # noqa: E402
        ring_biot_savart_self, _RING_KERNEL_CACHE,
        ring_boundary_Bc_pfaffian_selfconsistent, G0_DEFAULT,
        ScfConvergenceError,
    )
    from bdg_r7_inductance_sweep import A0, _phi_ext_bond, G_R7        # noqa: E402
    import bdg_r9_multisubband as ms                                   # noqa: E402
    from bdg_r9_multisubband_sweep import (                            # noqa: E402
        scalar_ring_map_ms, _find_fixed_points_ms, _phi_star_ms,
    )

_OUTPUT_DIR = __import__("os").path.join(
    __import__("os").path.dirname(__import__("os").path.abspath(__file__)),
    "..", "results", "multisubband")


def _out(name: str) -> str:
    import os
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════════
# 1 — M-subband ring Pfaffian sign (generalizes `_ring_pfaffian_sign`)
# ════════════════════════════════════════════════════════════════════════
def _ring_pfaffian_sign_ms(N: int, mu: float, B: float,
                           phi_bonds: np.ndarray, M: int,
                           delta: float = 1.0, W_nm: float = ms.WIRE_W_NM,
                           couple: bool = True) -> float:
    """Sign of the Pfaffian of the periodic M-subband ring LO BdG
    Hamiltonian (`ms.build_ring_lo_multisubband`) in the Majorana basis.
    Identical recipe to `_ring_pfaffian_sign` (R4) -- rescale by max|H|,
    apply the SAME Omega = (1/sqrt2)[[I,I],[iI,-iI]] transform on the full
    Nambu space (dimension 2NM here instead of 2N, since the electron block
    of the multi-subband Hamiltonian is 2NM x 2NM, not 2N x 2N) -- and
    assert real-antisymmetry (PH symmetry) before taking the Pfaffian. At
    M=1 this must reproduce `_ring_pfaffian_sign` bit-for-bit (validated
    below), since `build_ring_lo_multisubband` is already gated entrywise
    identical to `_build_ring_lo_phased` at M=1."""
    from pfapack import pfaffian as pf

    H = ms.build_ring_lo_multisubband(N, mu, B, phi_bonds, M=M, W_nm=W_nm,
                                      delta=delta, couple=couple)
    scale = float(np.max(np.abs(H)))
    nambu_half = 2 * N * M                      # electron-block dimension
    Ib = np.eye(nambu_half, dtype=complex)
    Omega = (1.0 / np.sqrt(2.0)) * np.block([[Ib, Ib], [1j * Ib, -1j * Ib]])
    A = 1j * Omega @ (H / scale) @ Omega.conj().T

    if np.max(np.abs(A.imag)) > 1e-8 or not np.allclose(A.real, -A.real.T, atol=1e-8):
        raise ValueError(
            "M-subband Majorana-basis matrix is not real antisymmetric "
            "(particle-hole symmetry broken?) -- Pfaffian sign undefined."
        )
    return float(np.sign(np.real(pf.pfaffian(A.real))))


def validate_pfaffian_ms_M1_gate(N: int = 40, n_points: int = 6) -> dict:
    """MANDATORY gate: `_ring_pfaffian_sign_ms(..., M=1, ...)` must reproduce
    `_ring_pfaffian_sign` bit-for-bit across a spread of (mu, B, phi_bonds)
    -- run BEFORE any M=2 Pfaffian sign is trusted."""
    rng = np.random.default_rng(11)
    max_diff = 0.0
    rows = []
    for _ in range(n_points):
        mu = float(rng.uniform(-2.0, 2.0))
        B = float(rng.uniform(0.0, 3.0))
        phi_bonds = rng.uniform(-0.05, 0.05, size=N)
        v_ref = _ring_pfaffian_sign(N, mu, B, phi_bonds, 1.0)
        v_ms = _ring_pfaffian_sign_ms(N, mu, B, phi_bonds, M=1, delta=1.0)
        max_diff = max(max_diff, abs(v_ref - v_ms))
        rows.append({"mu": mu, "B": B, "ref": v_ref, "ms": v_ms})
    ok = max_diff == 0.0
    if not ok:
        raise AssertionError(
            f"_ring_pfaffian_sign_ms M=1 gate FAILED: max_diff={max_diff}")
    return {"gate": "pfaffian_sign_ms_M1", "max_diff": max_diff, "pass": ok, "rows": rows}


# ════════════════════════════════════════════════════════════════════════
# 2 — M-subband self-consistent Bc bisection (generalizes R5/R7's function)
# ════════════════════════════════════════════════════════════════════════
def ring_boundary_Bc_pfaffian_selfconsistent_ms(
        N: int, mu: float, M: int,
        g: float = ms.G_PROBE,
        phi_ext_bond: float | None = None,
        delta: float = 1.0,
        ana_center: float | None = None,
        window: float = 0.05,
        n_coarse: int = 11,
        tol: float = 1e-6,
        root_tol: float = 1e-9,
        warm_hint: bool = True,
        W_nm: float = ms.WIRE_W_NM,
        couple: bool = True,
        a: float = A0,
        diagnostics: dict | None = None,
        _max_widenings: int = 4) -> float:
    """M-subband port of `ring_boundary_Bc_pfaffian_selfconsistent`: at every
    B the bisection probes, the self-consistent uniform self-phase phi*(B)
    is solved via `_phi_star_ms` (identical warm-hint-then-widened-scan
    recipe already validated throughout Track 4), then the M-subband
    Pfaffian invariant (`_ring_pfaffian_sign_ms`) is evaluated at
    phi_bonds = phi_ext + phi*(B). Coarse-scans and bisects for the first
    sign flip exactly as R5/R7 do, widening the window if none is found.

    ana_center: the coarse-scan center. Defaults to sqrt(mu^2+delta^2) (R5's
    OWN convention, needed for the M=1 bit-for-bit gate below). At M>1 the
    caller should pass the channel-of-interest's own naive single-channel
    condition, sqrt((mu-Delta_eps_n)^2+delta^2) -- NOT the M=1 default,
    which for M>1 chemical potentials deep in the window is a numerically
    unrelated scale and would send the search in the wrong place. This is a
    scan STARTING POINT only (with automatic widening); it does not presume
    the answer.

    Raises ScfConvergenceError (same class R5/R7 use) if `_phi_star_ms`
    finds no fixed point at some probed B -- a deliberate refusal to report
    a Bc built on a non-converged self-field, exactly as upstream."""
    if abs(g - ms.G_PROBE) > 1e-6 * ms.G_PROBE:
        raise ValueError(f"Track 4 uses g={ms.G_PROBE:.3e} everywhere; got g={g:.3e}")
    if phi_ext_bond is None:
        phi_ext_bond = _phi_ext_bond(N)
    ana = np.sqrt(mu ** 2 + delta ** 2) if ana_center is None else float(ana_center)
    mem: dict = {"hint": None}

    def phi_star(B: float) -> float:
        hint = mem["hint"] if warm_hint else None
        p, _slope = _phi_star_ms(N, mu, B, g, phi_ext_bond, a, M, hint=hint,
                                 tol=root_tol, delta=delta, W_nm=W_nm, couple=couple)
        if not np.isfinite(p):
            raise ScfConvergenceError(
                f"no self-consistent fixed point found at N={N}, mu={mu}, "
                f"M={M}, B={B:.6f}, g={g:.3g}")
        mem["hint"] = p
        return p

    def s(B: float) -> tuple:
        p = phi_star(B)
        phi_bonds = np.full(N, phi_ext_bond + p)
        val = _ring_pfaffian_sign_ms(N, mu, B, phi_bonds, M, delta=delta,
                                     W_nm=W_nm, couple=couple)
        if val == 0.0:
            val = _ring_pfaffian_sign_ms(N, mu, B + 1e-9, phi_bonds, M,
                                         delta=delta, W_nm=W_nm, couple=couple)
        return val, p

    h = window
    for _ in range(_max_widenings + 1):
        Bg = np.linspace(max(0.02, ana - h), ana + h, n_coarse)
        out = [s(B) for B in Bg]
        signs = np.array([o[0] for o in out])
        flip = np.where(signs[:-1] != signs[1:])[0]
        if flip.size == 0:
            h *= 2.0
            mem["hint"] = None
            continue
        i = int(flip[0])
        B_lo, B_hi = Bg[i], Bg[i + 1]
        s_lo = signs[i]
        mem["hint"] = out[i][1]
        p_at = out[i][1]
        while B_hi - B_lo > tol:
            Bm = 0.5 * (B_lo + B_hi)
            val, p_at = s(Bm)
            if val == s_lo:
                B_lo = Bm
            else:
                B_hi = Bm
        Bc = 0.5 * (B_lo + B_hi)
        if diagnostics is not None:
            diagnostics["phi_star_at_Bc"] = p_at
            diagnostics["coarse_signs"] = signs.tolist()
            diagnostics["n_flips_coarse"] = int(flip.size)
            diagnostics["ana_center_used"] = ana
            diagnostics["final_window"] = h
        return Bc
    if diagnostics is not None:
        diagnostics["ana_center_used"] = ana
        diagnostics["final_window"] = h
        diagnostics["widenings_exhausted"] = True
    return float("nan")


def ring_boundary_Bc_pfaffian_selfconsistent_ms_robust(
        N: int, mu: float, M: int, g: float = ms.G_PROBE,
        ana_center: float | None = None, window: float = 0.05,
        n_coarse: int = 11, tol: float = 1e-6, root_tol: float = 1e-9,
        delta: float = 1.0, W_nm: float = ms.WIRE_W_NM, couple: bool = True,
        a: float = A0, _max_widenings: int = 4,
        n_jitter_attempts: int = 5, diagnostics: dict | None = None) -> tuple:
    """R7's `ring_Bc_selfconsistent_a` retry discipline, ported: a SINGLE
    `ScfConvergenceError` from one coarse-grid path does not by itself prove
    no self-consistent boundary exists -- it may mean that particular path
    unluckily probed a B where the fixed point genuinely doesn't exist
    (a real, narrow dead-zone, R4/R5's own established phenomenon) while a
    slightly different grid would step over it. Retries with jittered
    (window, n_coarse) pairs (same jitter schedule R7 used); if EVERY
    attempt raises, that is the robust, multi-path-confirmed answer: a
    genuine dead-zone over the probed range, not a single unlucky path.

    Returns (Bc_or_nan, outcome_str, diagnostics_dict) where outcome_str is
    one of: 'found' (a genuine Bc located), 'genuine_dead_zone' (every
    jittered path raised ScfConvergenceError), 'no_flip_in_window' (every
    path completed without error but found no sign flip, even after
    _max_widenings widenings within each path)."""
    attempts = [(window, n_coarse),
                (window * 1.11, n_coarse),
                (window * 0.9, n_coarse + 2),
                (window * 1.23, n_coarse + 4),
                (window * 0.78, n_coarse + 6)][:max(1, n_jitter_attempts)]
    errors = []
    no_flip_count = 0
    for k, (w, nc) in enumerate(attempts):
        diag_k: dict = {}
        try:
            Bc = ring_boundary_Bc_pfaffian_selfconsistent_ms(
                N, mu, M, g=g, ana_center=ana_center, window=w, n_coarse=nc,
                tol=tol, root_tol=root_tol, warm_hint=True, delta=delta,
                W_nm=W_nm, couple=couple, a=a, diagnostics=diag_k,
                _max_widenings=_max_widenings)
            if np.isfinite(Bc):
                if diagnostics is not None:
                    diagnostics.update(diag_k)
                    diagnostics["attempt_index"] = k
                    diagnostics["attempts_tried"] = k + 1
                return Bc, "found", diag_k
            no_flip_count += 1
            if diagnostics is not None:
                diagnostics[f"attempt_{k}_no_flip"] = diag_k
        except ScfConvergenceError as e:
            errors.append(str(e))
            if diagnostics is not None:
                diagnostics[f"attempt_{k}_error"] = str(e)
    if diagnostics is not None:
        diagnostics["attempts_tried"] = len(attempts)
        diagnostics["n_errors"] = len(errors)
        diagnostics["n_no_flip"] = no_flip_count
        diagnostics["all_errors"] = errors
    if len(errors) == len(attempts):
        return float("nan"), "genuine_dead_zone", diagnostics or {}
    return float("nan"), "no_flip_in_window", diagnostics or {}


def validate_Bc_ms_M1_gate(N: int = 40, mu: float = 1.2,
                           g: float = G_R7) -> dict:
    """MANDATORY gate: `ring_boundary_Bc_pfaffian_selfconsistent_ms(..., M=1,
    ...)` must reproduce R5's `ring_boundary_Bc_pfaffian_selfconsistent`
    bit-for-bit at the SAME (N, mu, g, window, n_coarse) -- run BEFORE the
    M=2 boundary result is trusted. Uses N=40 (not 80) purely for gate
    speed; the algorithm has no N-dependence that would make N=40 pass and
    N=80 fail (both call the exact same code paths)."""
    _RING_KERNEL_CACHE.clear()
    diag_ref: dict = {}
    Bc_ref = ring_boundary_Bc_pfaffian_selfconsistent(
        N, mu, g=g, window=0.05, n_coarse=11, tol=1e-6, root_tol=1e-10,
        diagnostics=diag_ref)
    _RING_KERNEL_CACHE.clear()
    diag_ms: dict = {}
    Bc_ms = ring_boundary_Bc_pfaffian_selfconsistent_ms(
        N, mu, M=1, g=g, window=0.05, n_coarse=11, tol=1e-6, root_tol=1e-10,
        diagnostics=diag_ms)
    diff = abs(Bc_ref - Bc_ms)
    ok = diff == 0.0
    if not ok:
        raise AssertionError(
            f"ring_boundary_Bc_pfaffian_selfconsistent_ms M=1 gate FAILED: "
            f"Bc_ref={Bc_ref!r} Bc_ms={Bc_ms!r} diff={diff}")
    return {"gate": "Bc_selfconsistent_ms_M1", "Bc_ref": Bc_ref, "Bc_ms": Bc_ms,
            "diff": diff, "pass": ok}


# ════════════════════════════════════════════════════════════════════════
# 3 — Scalar Anderson mixing + plain Picard on the M-subband map, for the
#     dynamical-attractiveness / seed-independence test (R9-T1's method)
# ════════════════════════════════════════════════════════════════════════
def anderson_scalar_ms(N: int, mu: float, B: float, g: float, phi_ext_bond: float,
                       M: int, seed_phi: float = 0.0, m_history: int = 4,
                       beta: float = 1.0, reg: float = 1e-10,
                       tol_relA: float = 1e-6, max_iter: int = 200,
                       delta: float = 1.0, W_nm: float = ms.WIRE_W_NM,
                       couple: bool = True, a: float = A0,
                       fp_resid_check: float = 1e-3,
                       record_history: bool = False) -> dict:
    """Scalar (uniform-manifold) Type-II Anderson acceleration of the
    M-subband self-consistency map `scalar_ring_map_ms`, mirroring R9
    Track 1's `anderson_selfconsistent_ring_texture` exactly (same
    algorithm, same convergence criterion ||dA||/||A||<tol_relA, same
    fp_resid_check for flagging a genuine fixed point) but on the scalar
    uniform-manifold reduction already used throughout Track 4 (valid here
    for the identical reason R5 established it for M=1: a uniform input
    texture produces a uniform output on this translation-invariant ring,
    verified separately -- see `validate_uniform_reduction_ms`)."""
    x = float(seed_phi)
    xs: list = []
    fs: list = []
    hist: dict = {"rel_A": [], "fp_residual": []}
    converged = diverged = False
    rel_A = np.inf
    fp_res = np.inf
    it = 0
    x_new = x
    for it in range(1, max_iter + 1):
        gx = scalar_ring_map_ms(N, mu, B, g, phi_ext_bond, x, M, delta=delta,
                                t_hop=T_HOP_SI, lattice_const=a, W_nm=W_nm,
                                couple=couple)
        f = gx - x
        fp_res = abs(f)
        xs.append(x)
        fs.append(f)
        mk = min(m_history, len(fs) - 1)
        if mk == 0:
            x_new = x + beta * f
        else:
            dF = np.array([fs[-i] - fs[-i - 1] for i in range(1, mk + 1)])
            dX = np.array([xs[-i] - xs[-i - 1] for i in range(1, mk + 1)])
            # scalar Tikhonov-regularized least squares (1-D dF)
            AtA = float(dF @ dF)
            lam = reg * (AtA if AtA > 0 else 1.0)
            gamma = (dF * f) / (AtA + lam) if (AtA + lam) != 0 else np.zeros_like(dF)
            x_new = x + beta * f - float((dX + beta * dF) @ gamma)
        rel_A = abs(x_new - x) / max(abs(x_new), 1e-30)
        if record_history:
            hist["rel_A"].append(rel_A)
            hist["fp_residual"].append(fp_res)
        if rel_A < tol_relA:
            converged = True
            x = x_new
            break
        if not np.isfinite(rel_A) or abs(x_new) > 1e3:
            diverged = True
            break
        x = x_new
    return {"phi": x, "converged": converged, "diverged": diverged, "n_iter": it,
           "rel_A": rel_A, "fp_residual": fp_res,
           "fixed_point": bool(converged and fp_res < fp_resid_check),
           "history": hist if record_history else None}


def simple_mixing_scalar_ms(N: int, mu: float, B: float, g: float,
                            phi_ext_bond: float, M: int, seed_phi: float = 0.0,
                            beta: float = 1.0, tol_relA: float = 1e-6,
                            max_iter: int = 200, delta: float = 1.0,
                            W_nm: float = ms.WIRE_W_NM, couple: bool = True,
                            a: float = A0, record_history: bool = False) -> dict:
    """Plain (optionally damped) Picard iteration, scalar M-subband map --
    the SAME ||dA||/||A|| stop rule as `anderson_scalar_ms` so iteration
    counts are directly comparable (R9-T1's Task 3 comparison)."""
    x = float(seed_phi)
    hist: dict = {"rel_A": [], "fp_residual": []}
    converged = diverged = False
    rel_A = np.inf
    fp_res = np.inf
    it = 0
    for it in range(1, max_iter + 1):
        gx = scalar_ring_map_ms(N, mu, B, g, phi_ext_bond, x, M, delta=delta,
                                t_hop=T_HOP_SI, lattice_const=a, W_nm=W_nm,
                                couple=couple)
        f = gx - x
        fp_res = abs(f)
        x_new = x + beta * f
        rel_A = abs(x_new - x) / max(abs(x_new), 1e-30)
        if record_history:
            hist["rel_A"].append(rel_A)
            hist["fp_residual"].append(fp_res)
        if rel_A < tol_relA:
            converged = True
            x = x_new
            break
        if not np.isfinite(rel_A) or abs(x_new) > 1e3:
            diverged = True
            break
        x = x_new
    return {"phi": x, "converged": converged, "diverged": diverged, "n_iter": it,
           "rel_A": rel_A, "fp_residual": fp_res,
           "fixed_point": bool(converged and fp_res < 1e-3),
           "history": hist if record_history else None}


def path_independence_check_ms(N: int, mu: float, B: float, g: float, M: int,
                               a: float = A0, delta: float = 1.0,
                               W_nm: float = ms.WIRE_W_NM, couple: bool = True,
                               m_history: int = 4, max_iter: int = 200,
                               rng_seed: int = 2027) -> dict:
    """R9-T1's seed-independence test (`task2_convergence_above_gcrit`),
    ported to the M-subband scalar map: run Anderson mixing from several
    seeds and compare to the root-finder's principal root. Records whether
    all converged seeds land on the SAME fixed point (genuinely attracting,
    seed-independent) or scatter / fail to converge (only mathematically
    locatable, not dynamically attracting -- the R5-above-g_crit outcome)."""
    _RING_KERNEL_CACHE.clear()
    phi_ext = _phi_ext_bond(N)
    period = 2.0 * np.pi / N
    root, root_slope = _phi_star_ms(N, mu, B, g, phi_ext, a, M, hint=None,
                                    delta=delta, W_nm=W_nm, couple=couple)
    rng = np.random.default_rng(rng_seed)
    seeds = {
        "zero": 0.0,
        "random": float(rng.normal(scale=5e-3)),
        "far_neg": -0.02,
        "far_pos": 0.02,
        "at_root": root if np.isfinite(root) else 0.0,
    }
    seed_results = {}
    phis = []
    for sname, s0 in seeds.items():
        res = anderson_scalar_ms(N, mu, B, g, phi_ext, M, seed_phi=s0,
                                 m_history=m_history, tol_relA=1e-6,
                                 max_iter=max_iter, delta=delta, W_nm=W_nm,
                                 couple=couple, a=a)
        seed_results[sname] = {
            "converged": res["converged"], "n_iter": res["n_iter"],
            "phi": res["phi"], "fp_residual": res["fp_residual"],
            "fixed_point": res["fixed_point"], "rel_A": res["rel_A"],
            "matches_root": bool(np.isfinite(root) and abs(res["phi"] - root) < 1e-6),
        }
        if res["fixed_point"]:
            phis.append(res["phi"])
    pic = simple_mixing_scalar_ms(N, mu, B, g, phi_ext, M, seed_phi=0.0,
                                  beta=1.0, tol_relA=1e-6, max_iter=max_iter,
                                  delta=delta, W_nm=W_nm, couple=couple, a=a)
    conv_phis = np.array(phis)
    seed_spread = (float(np.max(conv_phis) - np.min(conv_phis))
                  if conv_phis.size >= 2 else 0.0)
    n_step = int(sum(sr["converged"] for sr in seed_results.values()))
    n_fp = int(sum(sr["fixed_point"] for sr in seed_results.values()))
    return {
        "N": N, "mu": mu, "B": B, "g": g, "M": M,
        "root_bisection": float(root), "root_slope": float(root_slope),
        "root_exists": bool(np.isfinite(root)),
        "root_picard_stable": bool(np.isfinite(root_slope) and abs(root_slope) < 1.0),
        "anderson_seeds": seed_results,
        "n_seeds_step_criterion": n_step, "n_seeds_fixed_point": n_fp,
        "n_seeds": len(seeds), "seed_spread": seed_spread,
        "seed_independent": bool(n_fp >= 2 and seed_spread < 1e-6),
        "picard_converged": pic["converged"], "picard_n_iter": pic["n_iter"],
        "picard_fp_residual": pic["fp_residual"],
    }


# ════════════════════════════════════════════════════════════════════════
# 4 — Uniform-manifold reduction check at M>1 (R5 proved this exactly for
#     M=1; Track 4 assumed, but never separately verified, that it still
#     holds at M>1 when it built the scalar map)
# ════════════════════════════════════════════════════════════════════════
def validate_uniform_reduction_ms(N: int = 40, M: int = 2, mu: float = 60.0,
                                  B: float = 40.0, phi_ext_bond: float = 0.01,
                                  phi_uniform: float = 0.01,
                                  delta: float = 1.0, couple: bool = True) -> dict:
    """Confirm a UNIFORM input texture produces a UNIFORM output current (and
    hence a uniform G) for the M-subband ring, the same way R5 confirmed it
    at M=1 (std(G)/|G| ~ 1e-13). Builds the FULL non-uniform-capable
    current-density machinery directly (not the scalar shortcut) at a
    uniform phi to check this independently, rather than assuming it."""
    phi_tot = np.full(N, phi_ext_bond + phi_uniform)
    H = ms.build_ring_lo_multisubband(N, mu, B, phi_tot, M=M, delta=delta, couple=couple)
    evals, evecs = eigh(H)
    I_bonds = ms.ring_current_density_ms(evecs, evals, N, M, T_HOP_SI, phi_tot)
    std_over_mean = float(np.std(I_bonds) / (abs(np.mean(I_bonds)) + 1e-300))
    return {"N": N, "M": M, "I_bonds_std": float(np.std(I_bonds)),
           "I_bonds_mean": float(np.mean(I_bonds)),
           "std_over_mean": std_over_mean,
           "uniform_to_high_precision": bool(std_over_mean < 1e-9)}


# ════════════════════════════════════════════════════════════════════════
# DRIVER — the full Task 1-3 investigation this module exists to answer:
# does the loop-gain-slope proxy's s>1 finding at M=2 (mu=131.762, the
# near-peak point of the pre-registered protocol-b sweep) translate into a
# measurable shift of the FULL self-consistent boundary, at one
# representative point, with the mandatory validation/N-independence
# discipline applied throughout?
# ════════════════════════════════════════════════════════════════════════
MU_M2_PEAK = 131.762      # the pre-registered protocol-b peak s point (s=3.0444)
B_M2_PEAK = 53.2455       # the pre-registered protocol_b_field() value there


def run_all_gates() -> dict:
    """All four MANDATORY validation gates, in order — every one must pass
    before any M=2 output from this module is trusted."""
    g1 = validate_pfaffian_ms_M1_gate()
    g2 = validate_uniform_reduction_ms(N=40, M=2, mu=60.0, B=40.0)
    g3 = validate_Bc_ms_M1_gate(N=80, mu=1.2)
    if not (g1["pass"] and g2["uniform_to_high_precision"] and g3["pass"]):
        raise AssertionError(f"validation gate failure: {g1} {g2} {g3}")
    return {"pfaffian_M1": g1, "uniform_reduction_M2": g2, "Bc_M1": g3}


def run_path_independence(N: int = 80) -> dict:
    """Task 1.3: is the root the loop-gain proxy found at (MU_M2_PEAK,
    B_M2_PEAK) genuinely dynamically attracting, or only mathematically
    locatable? R9-T1's seed-independence method, ported."""
    return path_independence_check_ms(N, MU_M2_PEAK, B_M2_PEAK, G_R7, M=2)


def run_boundary_at_N(N: int) -> dict:
    """The full self-consistent Bc bisection at M=2, mu=MU_M2_PEAK, for one
    N — plus the clean-Hamiltonian same-lobe comparison, if a genuine
    boundary is found. This is the expensive part (tens of minutes to
    ~1.5h/point on this machine, per the M^3 cost scaling already
    established in the Track 4 scoping report)."""
    diag: dict = {}
    Bc, outcome, _ = ring_boundary_Bc_pfaffian_selfconsistent_ms_robust(
        N, MU_M2_PEAK, M=2, g=G_R7, ana_center=B_M2_PEAK, window=3.0,
        n_coarse=11, tol=1e-6, root_tol=1e-9, _max_widenings=6,
        n_jitter_attempts=5, diagnostics=diag)
    result = {"N": N, "mu": MU_M2_PEAK, "M": 2, "g": G_R7,
             "outcome": outcome, "Bc": Bc, "diagnostics": diag}
    if outcome == "found" and np.isfinite(Bc):
        phi_ext = _phi_ext_bond(N)

        def clean_sign(B):
            return _ring_pfaffian_sign_ms(N, MU_M2_PEAK, B, np.full(N, phi_ext), 2)

        span = 2.0
        Bg = np.linspace(Bc - span, Bc + span, 33)
        signs = [clean_sign(b) for b in Bg]
        flips = [i for i in range(len(signs) - 1) if signs[i] != signs[i + 1]]
        if flips:
            i = min(flips, key=lambda i: abs(0.5 * (Bg[i] + Bg[i + 1]) - Bc))
            lo, hi, s_lo = Bg[i], Bg[i + 1], signs[i]
            while hi - lo > 1e-6:
                m = 0.5 * (lo + hi)
                lo, hi = (m, hi) if clean_sign(m) == s_lo else (lo, m)
            Bc_clean = 0.5 * (lo + hi)
            result["Bc_clean_same_lobe"] = Bc_clean
            result["shift_abs"] = Bc - Bc_clean
            result["shift_rel"] = (Bc - Bc_clean) / Bc_clean
    return result


def generate_data(Ns: tuple = (80, 120, 60)) -> dict:
    """Full driver: gates -> path-independence at N=80 -> boundary bisection
    at every N in `Ns` (default order matches the order this investigation
    was actually run in: N=80 first, since it is the project's standard N;
    then N=120 and N=60 for the N-independence check once N=80 gave a
    result to check). Writes each N's result to its own JSON file as it
    completes (so a long run's partial progress is never lost), plus a
    combined summary file at the end."""
    out: dict = {"gates": run_all_gates(),
                "path_independence_N80": run_path_independence(80),
                "by_N": {}}
    for N in Ns:
        res = run_boundary_at_N(N)
        out["by_N"][str(N)] = res
        with open(_out(f"r9_multisubband_bc_M2_N{N}.json"), "w") as f:
            import json
            json.dump(res, f, indent=2, default=str)
    with open(_out("r9_multisubband_bc_M2_summary.json"), "w") as f:
        import json
        json.dump(out, f, indent=2, default=str)
    return out


if __name__ == "__main__":
    import json
    import sys
    with contextlib.suppress(Exception):
        sys.stdout.reconfigure(encoding="utf-8")
    data = generate_data()
    for N, res in data["by_N"].items():
        print(f"N={N}: outcome={res['outcome']} Bc={res.get('Bc')} "
              f"shift_rel={res.get('shift_rel')}")
