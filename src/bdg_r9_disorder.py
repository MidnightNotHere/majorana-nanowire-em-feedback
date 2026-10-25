"""
ISEF 2027 — R9 (Track 3): on-site disorder in the BdG / Pfaffian / self-
consistent-ring stack
===============================================================================
Every prior phase (R1-R8, R9 Tracks 1-2) used a perfectly translationally-
ordered lattice. Real InAs/Al nanowires are disordered (impurities, interface
roughness, charge inhomogeneity) -- a well-documented concern in the Majorana
nanowire literature, where disorder-induced trivial sub-gap states have been
proposed as an explanation for some contested zero-bias-peak claims (see
Task 3 below for the verified primary sources). This module adds on-site
disorder as a STRICTLY ADDITIVE layer on top of the existing clean
Hamiltonian: nothing in R1-R8 (or R9 Tracks 1-2) is modified; every clean
builder/detector is imported READ-ONLY.

────────────────────────────────────────────────────────────────────────────
THE DISORDER PRIMITIVE (the one place disorder is defined)
────────────────────────────────────────────────────────────────────────────
On-site disorder is a random potential eps_i, i=1..N, added to the chemical-
potential term at each site, i.e. mu -> mu - eps_i (equivalently the on-site
energy (2t-mu) -> (2t-mu) + eps_i in the electron block). Every clean builder
in this repo (build_bdg_lo, _build_ring_lo_phased) uses the SAME kron(spin,
site) Nambu ordering: first 2N = electrons (both spins), last 2N = holes,
with SITE the fast (inner) index. So the on-site electron term at site i,
both spins, sits at diagonal positions i and N+i; disorder is injected as

    V(eps) = block_diag( kron(I2, diag(eps)),  -kron(I2, diag(eps)) )   (4N x 4N, REAL)
    H_dis  = H_clean + V(eps)

PARTICLE-HOLE SYMMETRY PROOF (why the hole block needs -eps, not +eps):
check_ph_symmetry_complex (bdg_eigensolver_r1_r2_LOAdapted.py) requires
Sigma_x @ H.conj() @ Sigma_x == -H, where Sigma_x = [[0,I_2N],[I_2N,0]] is the
electron<->hole block swap. For a real V = [[A,0],[0,-A]] (A = kron(I2,diag(eps))),
Sigma_x @ V @ Sigma_x = [[-A,0],[0,A]] = -V exactly (a direct 2x2 block-matrix
computation), and since eps is real, V.conj() = V, so Sigma_x @ V.conj() @
Sigma_x = -V for ANY realization eps. Adding V to a clean H that already
satisfies the PH identity therefore preserves it exactly, for ANY eps -- this
is verified numerically (not just algebraically) in validate_W0() below via
check_ph_symmetry_complex on a real W>0 sample.

Because every downstream clean function (the Pfaffian invariant, the ring
boundary detector, the open-wire sub-gap detector, the self-consistent map)
builds its own H internally and then diagonalizes it, disorder-aware
"shells" are written here that call the clean builder, add V(eps), and run
the IDENTICAL downstream algorithm. At eps=0 (W=0) add_onsite_disorder is a
no-op, so every disordered shell reproduces its clean counterpart exactly
(not approximately) -- this is the mandatory W=0 validation gate
(validate_W0), which MUST pass before any W>0 result is trusted.

────────────────────────────────────────────────────────────────────────────
ENSEMBLE AVERAGING (mandatory -- a single realization is never a valid result)
────────────────────────────────────────────────────────────────────────────
Every reported disorder quantity is computed over an ENSEMBLE of independent
realizations at fixed (mu, B, W), each with its own on-site draw eps ~
Uniform[-W/2, W/2]. Seeding is via numpy SeedSequence.spawn (see
seed_sequences()) -- deterministic, reproducible, independent per realization,
with NO reliance on a global np.random state (required for correct behaviour
under multiprocessing / Windows spawn). Every reported number carries its
realization count and its spread (std / SEM for continuous quantities,
Bernoulli std for the topological-fraction indicator); convergence in
realization count is checked explicitly (convergence_check) rather than
assumed.

────────────────────────────────────────────────────────────────────────────
THE ONE APPROXIMATION IN THE SELF-CONSISTENT (Task 4.3) CHANNEL -- FLAGGED
────────────────────────────────────────────────────────────────────────────
R5's self-consistent solver reduces the ring self-field problem to a SCALAR
map on the uniform manifold, and that reduction is EXACT there because a
clean, translation-invariant ring guarantees uniform-in -> uniform-out
(verified by R5's `find_uniform_fixed_points`/`ring_mode_gains` machinery).
Disorder breaks translational invariance, so a uniform input texture no
longer produces a perfectly uniform current/self-field response in general.
`_disordered_scalar_ring_map` below evaluates the FULL per-bond response and
projects it back onto the uniform manifold via np.mean -- a controlled
approximation (mean-field over the disorder-broken response), NOT an exact
reduction. `validate_W0` confirms this projection is exact at W=0 (where the
response IS uniform); `_scshift_worker` additionally records the realized
per-bond NON-uniformity (std/mean of the raw per-bond output) at a sample
point so the size of this approximation is visible, not hidden. Solving the
full N-dimensional non-uniform self-consistency per disorder realization
(Newton-Krylov, as R5's nonuniform_fixed_point_search does for the clean
diagnostic) would multiply the cost of the already most-expensive R9 track
by another N-dimensional root solve per realization x grid point x W --
explicitly out of scope here; this is the stated limitation for Task 4.3.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import time

import numpy as np
from scipy.linalg import eigh, eigvalsh

# ── Read-only imports from the validated R1/R4/R5/R7 modules -------------
with contextlib.redirect_stdout(io.StringIO()):
    from bdg_eigensolver_r1_r2_LOAdapted import (        # noqa: E402
        build_bdg_lo, check_ph_symmetry_complex, LO_PARAMS,
    )
    from bdg_r4_em_coupling import (                     # noqa: E402
        _build_ring_lo_phased, _ring_pfaffian_sign, ring_boundary_Bc_pfaffian,
        apply_peierls_phases, compute_current_density, biot_savart_self,
        compute_peierls_from_A, T_HOP_SI, LATTICE_A,
    )
    from bdg_r5_selfconsistent_ring import (             # noqa: E402
        ring_current_density, ring_biot_savart_self,
        ring_boundary_Bc_pfaffian_selfconsistent, ScfConvergenceError,
        G0_DEFAULT, _RING_KERNEL_CACHE,
    )
    from bdg_r7_inductance_sweep import G_R7, A0, _phi_ext_bond  # noqa: E402
    from bdg_parallel import (                           # noqa: E402
        parallel_map, check_equivalence, default_workers, _pin_blas,
    )

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "disorder")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════════
# TASK 1.1/1.2 — the disorder primitive + reproducible seeding
# ════════════════════════════════════════════════════════════════════════
def disorder_epsilon(N: int, W: float, rng: np.random.Generator) -> np.ndarray:
    """One realization of on-site disorder: eps_i ~ Uniform[-W/2, W/2], i.i.d.
    over the N sites. W=0 returns exact zeros (the W=0 gate's fast path)."""
    if W == 0.0:
        return np.zeros(N, dtype=float)
    return rng.uniform(-W / 2.0, W / 2.0, size=N)


def add_onsite_disorder(H: np.ndarray, N: int, eps: np.ndarray) -> np.ndarray:
    """Add the on-site disorder potential V(eps) (see module docstring for the
    PH-symmetry proof) to a 4N x 4N BdG Hamiltonian H. eps=0 (or None) is an
    EXACT no-op (returns H unchanged, not a copy with zero added) -- this is
    what makes the W=0 validation gate bit-identical rather than merely
    close."""
    if eps is None or not np.any(eps):
        return H
    eps = np.asarray(eps, dtype=float)
    a = np.concatenate([eps, eps])          # both spins, kron(I2, diag(eps))
    v_diag = np.concatenate([a, -a])        # electron +eps, hole -eps (PH)
    return H + np.diag(v_diag).astype(H.dtype)


def seed_sequences(base_seed: int, n_real: int) -> list:
    """Independent, reproducible per-realization seeds via numpy's recommended
    multiprocessing-safe scheme: SeedSequence(base_seed).spawn(n_real). Each
    realization k gets seeds[k], reconstructed in the worker via
    np.random.default_rng(seeds[k]) -- NEVER a global np.random call, and
    NEVER an uncontrolled default (e.g. os.urandom-seeded) that could repeat
    across workers or be non-reproducible between runs."""
    return np.random.SeedSequence(int(base_seed)).spawn(int(n_real))


# ════════════════════════════════════════════════════════════════════════
# TASK 1.1 — disorder-aware shells (W=0 must reproduce the clean function
# EXACTLY -- see validate_W0)
# ════════════════════════════════════════════════════════════════════════
def _disordered_ring_pfaffian_sign(N: int, mu: float, B: float,
                                   eps: np.ndarray,
                                   phi_bonds: np.ndarray | None = None,
                                   delta: float = 1.0) -> float:
    """Disordered port of bdg_r4_em_coupling._ring_pfaffian_sign: identical
    Majorana-basis Pfaffian-sign algorithm, applied to H + V(eps)."""
    from pfapack import pfaffian as pf       # local import; heavy module

    H = _build_ring_lo_phased(N, mu, B, phi_bonds, delta)
    H = add_onsite_disorder(H, N, eps)
    scale = float(np.max(np.abs(H)))
    M = 2 * N
    Ib = np.eye(M, dtype=complex)
    Omega = (1.0 / np.sqrt(2.0)) * np.block([[Ib, Ib], [1j * Ib, -1j * Ib]])
    A = 1j * Omega @ (H / scale) @ Omega.conj().T

    if np.max(np.abs(A.imag)) > 1e-8 or not np.allclose(A.real, -A.real.T, atol=1e-8):
        raise ValueError(
            "Majorana-basis matrix is not real antisymmetric after disorder "
            "injection (PH symmetry broken?) -- Pfaffian sign undefined."
        )
    return float(np.sign(np.real(pf.pfaffian(A.real))))


def _disordered_ring_boundary_Bc(N: int, mu: float, eps: np.ndarray,
                                 delta: float = 1.0,
                                 phi_bonds: np.ndarray | None = None,
                                 window: float = 0.6, n_coarse: int = 61,
                                 tol: float = 1e-4,
                                 _max_widenings: int = 4) -> float:
    """Disordered port of ring_boundary_Bc_pfaffian: identical coarse-scan +
    bisection algorithm, invariant evaluated with disorder eps."""
    ana = np.sqrt(mu ** 2 + delta ** 2)

    def s(B):
        val = _disordered_ring_pfaffian_sign(N, mu, B, eps, phi_bonds, delta)
        if val == 0.0:
            val = _disordered_ring_pfaffian_sign(N, mu, B + 1e-6, eps, phi_bonds, delta)
        return val

    h = window
    for _ in range(_max_widenings + 1):
        Bg = np.linspace(max(0.02, ana - h), ana + h, n_coarse)
        signs = np.array([s(B) for B in Bg])
        flip = np.where(signs[:-1] != signs[1:])[0]
        if flip.size == 0:
            h *= 2.0
            continue
        i = int(flip[0])
        B_lo, B_hi = Bg[i], Bg[i + 1]
        s_lo = signs[i]
        while B_hi - B_lo > tol:
            Bm = 0.5 * (B_lo + B_hi)
            if s(Bm) == s_lo:
                B_lo = Bm
            else:
                B_hi = Bm
        return 0.5 * (B_lo + B_hi)
    return float("nan")


def disordered_openwire_subgap(N: int, mu: float, B: float, eps: np.ndarray,
                               delta: float = 1.0, phi_ext: float = 0.0,
                               threshold: float = 0.25) -> dict:
    """Disordered port of the openwire_boundary_Bc spectral detector (the
    engine of R4's finite-size sub-gap-state artifact / Task 4.2's spurious-
    state question): build the clean open wire, add on-site disorder, apply
    any external flux, and count states within `threshold` of zero energy.
    At W=0, phi_ext=0 this reduces exactly to a single point of
    openwire_boundary_Bc's internal `_scan` (bit-identical: add_onsite_disorder
    and apply_peierls_phases are both exact no-ops at eps=0/phi_ext=0)."""
    with contextlib.redirect_stdout(io.StringIO()):
        H = build_bdg_lo(N, mu=mu, B_zeeman=B, delta=delta)
    H = add_onsite_disorder(H, N, eps)
    H = apply_peierls_phases(H, N, phi_ext)
    evals = eigvalsh(H)
    abs_e = np.abs(evals)
    return {"min_abs_E": float(np.min(abs_e)),
            "n_subgap": int(np.sum(abs_e < threshold))}


def _disordered_scalar_ring_map(N: int, mu: float, B: float, g: float,
                                phi_ext_bond: float, phi_uniform: float,
                                eps: np.ndarray, delta: float = 1.0,
                                t_hop: float = T_HOP_SI,
                                lattice_const: float = LATTICE_A,
                                wire_radius: float | None = None,
                                diagnostics: dict | None = None) -> float:
    """Disordered port of R5 scalar_ring_map. See module docstring: disorder
    breaks the clean model's exact uniform-in/uniform-out reduction, so this
    projects the (generally non-uniform) per-bond response onto the uniform
    manifold via np.mean -- a controlled approximation, exact at W=0. If
    `diagnostics` is given, records the realized per-bond non-uniformity
    (std/|mean|) of the raw response so the approximation's size is visible."""
    phi_tot = np.full(N, float(phi_ext_bond) + float(phi_uniform))
    H = _build_ring_lo_phased(N, mu, B, phi_tot, delta)
    H = add_onsite_disorder(H, N, eps)
    evals, evecs = eigh(H)
    I_bonds = ring_current_density(evecs, evals, N, t_hop, phi_tot)
    A_t = ring_biot_savart_self(I_bonds, N, lattice_const, wire_radius)
    out = g * compute_peierls_from_A(A_t, lattice_const)
    if diagnostics is not None:
        scale = max(abs(float(np.mean(out))), 1e-300)
        diagnostics["texture_nonuniformity"] = float(np.std(out) / scale)
    return float(np.mean(out))


def _disordered_find_uniform_fixed_points(N, mu, B, g, eps,
                                          phi_ext_bond: float | None = None,
                                          span: float | None = None,
                                          n_scan: int = 49, tol: float = 1e-12,
                                          delta: float = 1.0) -> list:
    """Disordered port of R5 find_uniform_fixed_points (identical bisection /
    genuine-root-vs-jump discrimination), on the disordered scalar map."""
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    if span is None:
        span = 1.5 * 2.0 * np.pi / N

    def f(p):
        return _disordered_scalar_ring_map(N, mu, B, g, phi_ext_bond, p, eps, delta) - p

    grid = np.linspace(-span, span, n_scan)
    vals = np.array([f(p) for p in grid])
    roots = []
    for i in np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]:
        lo, hi = float(grid[i]), float(grid[i + 1])
        flo = vals[i]
        while hi - lo > tol:
            m = 0.5 * (lo + hi)
            fm = f(m)
            if np.sign(fm) == np.sign(flo):
                lo, flo = m, fm
            else:
                hi = m
        m = 0.5 * (lo + hi)
        if abs(f(m)) < 1e-4:
            fd = 1e-6
            slope = (_disordered_scalar_ring_map(N, mu, B, g, phi_ext_bond, m + fd, eps, delta)
                     - _disordered_scalar_ring_map(N, mu, B, g, phi_ext_bond, m - fd, eps, delta)
                     ) / (2 * fd)
            roots.append({"phi": m, "slope": float(slope),
                          "picard_stable": bool(abs(slope) < 1.0)})
    roots.sort(key=lambda r: abs(r["phi"]))
    return roots


def _disordered_principal_uniform_fixed_point(N, mu, B, g, eps,
                                              phi_ext_bond: float | None = None,
                                              hint: float | None = None,
                                              tol: float = 1e-12,
                                              delta: float = 1.0) -> float:
    """Disordered port of R5 principal_uniform_fixed_point (identical
    warm-hint bracket -> full-scan fallback -> widened-scan fallback)."""
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    if hint is not None and np.isfinite(hint):
        period = 2.0 * np.pi / N
        w = 0.25 * period

        def f(p):
            return _disordered_scalar_ring_map(N, mu, B, g, phi_ext_bond, p, eps, delta) - p

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
                return m
    roots = _disordered_find_uniform_fixed_points(N, mu, B, g, eps, phi_ext_bond,
                                                  tol=tol, delta=delta)
    if not roots:
        roots = _disordered_find_uniform_fixed_points(
            N, mu, B, g, eps, phi_ext_bond,
            span=3.0 * 2.0 * np.pi / N, n_scan=121, tol=tol, delta=delta)
    return roots[0]["phi"] if roots else float("nan")


def _disordered_Bc_selfconsistent(N: int, mu: float, eps: np.ndarray,
                                  g: float = G_R7,
                                  phi_ext_bond: float | None = None,
                                  delta: float = 1.0, window: float = 0.05,
                                  n_coarse: int = 11, tol: float = 1e-6,
                                  root_tol: float = 1e-10,
                                  warm_hint: bool = True,
                                  diagnostics: dict | None = None,
                                  _max_widenings: int = 3) -> float:
    """_max_widenings=3 (bumped from a post-hoc-review value of 2, still below
    the clean solver's 4): measured directly (the
    quickfull smoke run) that a disordered realization failing to converge
    costs MORE than one that succeeds, because each widening repeats a full
    coarse B-scan x root-search -- so a realization destined to fail (a
    genuine, reportable outcome, see _scshift_worker) should fail FAST, not
    burn 2x the wall-clock on retries unlikely to help. This trades a
    slightly higher measured convergence_failed rate for a bounded, honest
    per-task cost; the failure fraction itself remains a reported result."""
    """Disordered port of ring_boundary_Bc_pfaffian_selfconsistent (solver=
    'root' only -- root-finding needs no Picard stability, so it applies
    unchanged to the disordered map; see module docstring for the projection
    approximation this uses). Raises ScfConvergenceError (imported from R5,
    NOT redefined) if the disordered self-field cannot be solved at some
    probed B, exactly mirroring the clean function's refusal to report a
    B_c built on a non-converged texture.

    BUG FOUND + FIXED ON POST-HOC REVIEW: the coarse B-scan below used to be
    centered on the fixed CLEAN analytic value `ana = sqrt(mu^2+delta^2)`,
    inherited unchanged from the clean solver -- correct there ONLY because
    the clean rigid boundary coincides with `ana` exactly. Task 4.1
    independently established that disorder itself displaces the RIGID
    (non-self-consistent) boundary away from `ana` by a mean of +0.45 to
    +1.0 (up to +2.9) across the realistic W band -- far outside this
    solver's +-0.20 max search half-width -- so the original version was
    searching in the wrong place, not failing to find a self-consistent
    fixed point because none exists. A direct diagnostic (comparing this
    solver's own widest-searched window against the wide, Task-4.1-trusted
    rigid detector's boundary location for the SAME 20 seeds at each W)
    confirmed this is the dominant failure mode (15-18/20 boundaries fall
    outside the old window at W=10-46), not a genuine self-consistency-loop
    instability. Fix: center the scan on the PER-REALIZATION rigid boundary
    (cheap: one call to `_disordered_ring_boundary_Bc`, which found a
    boundary in 100/100 realizations across every W tested) instead of the
    fixed clean `ana`. At W=0 these coincide exactly, so the W=0 gate is
    unaffected. The self-consistent correction ON TOP OF the (already
    disorder-shifted) rigid boundary is still expected to stay small, per
    the R4-R5-R7 null lineage -- so the search window itself did not need
    to grow much, only its center did (a small `_max_widenings` bump 2->3
    was made after the fix, since a post-fix sample run still showed a
    residual few-percent-to-tens-of-percent rate of GENUINE
    ScfConvergenceError -- the self-consistency root find failing at a
    specific probed B even when correctly centered -- and a modest widening
    margin catches the cases where the self-consistent correction itself,
    not just the rigid displacement, is larger than the +-0.05 default. This
    residual rate is real and is reported as its own (smaller, secondary)
    finding, not eliminated by the centering fix -- see the methods note."""
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    if eps is None or not np.any(eps):
        # W=0 fast path: use the EXACT analytic value, not the rigid
        # detector's tol=1e-4 bisection -- this is what keeps the mandatory
        # W=0 gate BIT-IDENTICAL (not merely close) to the clean solver,
        # matching add_onsite_disorder's own eps=0 exact-no-op convention.
        ana = np.sqrt(mu ** 2 + delta ** 2)
    else:
        ana = _disordered_ring_boundary_Bc(N, mu, eps, delta=delta)
        if not np.isfinite(ana):
            raise ScfConvergenceError(
                f"no rigid (non-self-consistent) boundary found at N={N}, mu={mu}, "
                f"g={g:.3g} (eps realization has no Pfaffian sign flip within the "
                f"wide +-0.6 rigid-detector window) -- cannot search for a "
                f"self-consistent correction without a starting point.")
    mem = {"hint": None}

    def phi_star(B):
        hint = mem["hint"] if warm_hint else None
        p = _disordered_principal_uniform_fixed_point(
            N, mu, B, g, eps, phi_ext_bond, hint=hint, tol=root_tol, delta=delta)
        if not np.isfinite(p):
            raise ScfConvergenceError(
                f"no disordered self-consistent fixed point at N={N}, mu={mu}, "
                f"B={B:.6f}, g={g:.3g}")
        mem["hint"] = p
        return p

    def s(B):
        p = phi_star(B)
        phi_bonds = np.full(N, phi_ext_bond + p)
        val = _disordered_ring_pfaffian_sign(N, mu, B, eps, phi_bonds, delta)
        if val == 0.0:
            val = _disordered_ring_pfaffian_sign(N, mu, B + 1e-9, eps, phi_bonds, delta)
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
        return Bc
    return float("nan")


# ════════════════════════════════════════════════════════════════════════
# TASK 1.3 — W=0 exact-reproduction gate: MUST pass before any W>0 result
# ════════════════════════════════════════════════════════════════════════
def validate_W0(N: int = 80, mu: float = 1.2, verbose: bool = True) -> dict:
    """Confirm every disorder-aware shell reproduces its clean counterpart
    EXACTLY at W=0 (eps=0), plus a numerical (not just algebraic) PH-symmetry
    check at W>0. run_full() calls this first and refuses to proceed to any
    W>0 computation if it fails."""
    res = {}
    B = np.sqrt(mu ** 2 + 1.0) - 0.05
    eps0 = np.zeros(N)

    # (1) invariant: exact equality across a B-grid straddling B_c
    Bgrid = np.linspace(B - 0.1, B + 0.1, 9)
    diffs = [abs(_ring_pfaffian_sign(N, mu, Bt) - _disordered_ring_pfaffian_sign(N, mu, Bt, eps0))
             for Bt in Bgrid]
    res["invariant_max_abs_diff"] = float(max(diffs))

    # (2) ring boundary Bc: same bisection, zero perturbation -> bit-identical
    bc_clean = ring_boundary_Bc_pfaffian(N, mu)
    bc_dis = _disordered_ring_boundary_Bc(N, mu, eps0)
    res["boundary_absdiff"] = float(abs(bc_clean - bc_dis))

    # (3) self-consistent Bc shift: same root-finding map, zero perturbation
    scbc_clean = ring_boundary_Bc_pfaffian_selfconsistent(N, mu, g=G_R7)
    scbc_dis = _disordered_Bc_selfconsistent(N, mu, eps0, g=G_R7)
    res["selfconsistent_absdiff"] = float(abs(scbc_clean - scbc_dis))

    # (4) open-wire sub-gap detector: bit-identical spectrum at eps=0
    with contextlib.redirect_stdout(io.StringIO()):
        H_clean = build_bdg_lo(N, mu=mu, B_zeeman=B, delta=1.0)
    e_clean = float(np.min(np.abs(eigvalsh(H_clean))))
    sg = disordered_openwire_subgap(N, mu, B, eps0)
    res["subgap_absdiff"] = float(abs(e_clean - sg["min_abs_E"]))

    # (5) PH symmetry at a genuine W>0 sample (numerical, not just algebraic)
    rng = np.random.default_rng(12345)
    eps_sample = disorder_epsilon(N, 5.0, rng)
    H_test = add_onsite_disorder(H_clean.copy(), N, eps_sample)
    res["ph_symmetry_eps_at_W5"] = float(check_ph_symmetry_complex(H_test, N))

    res["ALL_PASS"] = bool(
        res["invariant_max_abs_diff"] == 0.0
        and res["boundary_absdiff"] < 1e-9
        and res["selfconsistent_absdiff"] < 1e-9
        and res["subgap_absdiff"] < 1e-12
        and res["ph_symmetry_eps_at_W5"] < 1e-10
    )
    if verbose:
        print(json.dumps(res, indent=2))
    return res


# ════════════════════════════════════════════════════════════════════════
# TASK 4.1 engine — topological classification relative to CLEAN references
# ════════════════════════════════════════════════════════════════════════
# The Pfaffian sign is only meaningful RELATIVE to a reference (the absolute
# +1/-1 label carries no physical meaning by itself). We fix that reference
# per (N, mu, delta) using two CLEAN (W=0) points: deep in the clean
# topological phase and deep in the clean trivial phase (B=0), at the SAME mu.
# A realization is then classified "topological" iff its sign matches the
# clean topological reference's sign -- sign-convention-agnostic and does not
# depend on any assumption about which absolute sign "means" topological.
_TOPO_REF_CACHE: dict = {}


def _topo_refs(N: int, mu: float, delta: float, B_deep_topo: float,
              B_ref_trivial: float) -> tuple:
    key = (N, mu, delta, B_deep_topo, B_ref_trivial)
    if key not in _TOPO_REF_CACHE:
        st = _disordered_ring_pfaffian_sign(N, mu, B_deep_topo, np.zeros(N), delta=delta)
        sv = _disordered_ring_pfaffian_sign(N, mu, B_ref_trivial, np.zeros(N), delta=delta)
        _TOPO_REF_CACHE[key] = (st, sv)
    return _TOPO_REF_CACHE[key]


def classify_topological(N: int, mu: float, B: float, eps: np.ndarray,
                         delta: float = 1.0, B_deep_topo: float | None = None,
                         B_ref_trivial: float = 0.0):
    """True/False if this realization's sign matches the clean-topological /
    clean-trivial reference; None if the two references coincide (degenerate,
    e.g. exactly mu=0 pathological choices of B_ref)."""
    if B_deep_topo is None:
        B_deep_topo = 3.0 * np.sqrt(mu ** 2 + delta ** 2)
    sign_topo, sign_triv = _topo_refs(N, mu, delta, B_deep_topo, B_ref_trivial)
    if sign_topo == sign_triv:
        return None
    sign_here = _disordered_ring_pfaffian_sign(N, mu, B, eps, delta=delta)
    return bool(sign_here == sign_topo)


# ════════════════════════════════════════════════════════════════════════
# TASK 2 — top-level (picklable) parallel workers
# ════════════════════════════════════════════════════════════════════════
def _topo_worker(N, mu, B, W, seed_seq, B_deep_topo=None, B_ref_trivial=0.0,
                 delta=1.0):
    """Task 4.1 unit: one disorder realization's topological classification."""
    rng = np.random.default_rng(seed_seq)
    eps = disorder_epsilon(N, W, rng)
    topo = classify_topological(N, mu, B, eps, delta, B_deep_topo, B_ref_trivial)
    return {"W": float(W), "topological": topo}


def _subgap_worker(N, mu, B, W, seed_seq, threshold=0.25, delta=1.0):
    """Task 4.2 unit: one disorder realization's open-wire sub-gap spectrum."""
    rng = np.random.default_rng(seed_seq)
    eps = disorder_epsilon(N, W, rng)
    sg = disordered_openwire_subgap(N, mu, B, eps, delta=delta, threshold=threshold)
    return {"W": float(W), "min_abs_E": sg["min_abs_E"], "n_subgap": sg["n_subgap"]}


def _boundary_worker(N, mu, W, seed_seq, delta=1.0):
    """Rigid (non-self-consistent) ring boundary Bc under disorder, vs the
    clean baseline -- cheap engine used for the benchmark and as a secondary
    (non-feedback) boundary-smearing cross-check."""
    rng = np.random.default_rng(seed_seq)
    eps = disorder_epsilon(N, W, rng)
    bc = _disordered_ring_boundary_Bc(N, mu, eps, delta=delta)
    base = ring_boundary_Bc_pfaffian(N, mu, delta)
    shift = 100.0 * (bc - base) / base if np.isfinite(bc) else float("nan")
    return {"W": float(W), "Bc": float(bc), "baseline_Bc": float(base),
            "shift_pct": float(shift)}


_SC_BASELINE_CACHE: dict = {}


def _sc_baseline(N, mu, g, delta=1.0):
    """Clean self-consistent Bc baseline is W-independent; each solve costs
    ~30s at N=80 (many eigh calls in the root search), so it is cached per
    (N,mu,g,delta) instead of being recomputed on every ensemble member --
    that redundancy alone was doubling the measured per-task cost of
    _scshift_worker before this cache was added."""
    key = (N, mu, g, delta)
    if key not in _SC_BASELINE_CACHE:
        _SC_BASELINE_CACHE[key] = ring_boundary_Bc_pfaffian_selfconsistent(
            N, mu, g=g, delta=delta)
    return _SC_BASELINE_CACHE[key]


def _scshift_worker(N, mu, W, seed_seq, g=G_R7, delta=1.0):
    """Task 4.3 unit: self-consistent (inductive-feedback) boundary shift
    under disorder vs the clean self-consistent baseline. Records
    convergence failures explicitly (does not silently drop them).

    NB (bug found + fixed in the first --full run): `_disordered_Bc_
    selfconsistent` has TWO distinct failure modes -- it raises
    ScfConvergenceError when the per-B self-consistent root search itself
    fails, but returns a silent `float("nan")` (no exception) when the
    coarse B-scan never brackets a Pfaffian sign flip at all even after
    every widening. The original code here only treated the FIRST mode as
    `convergence_failed=True`; the second mode returned nan while
    convergence_failed stayed False, which massively undercounted the true
    failure rate in the raw --full run (measured 16-19 silent failures out
    of 20 at realistic W, vs 1-2 explicitly flagged). Fixed by checking
    `np.isfinite(bc)` unconditionally -- both failure modes are now the
    same reported outcome."""
    rng = np.random.default_rng(seed_seq)
    eps = disorder_epsilon(N, W, rng)
    base = _sc_baseline(N, mu, g, delta)
    diag: dict = {}
    try:
        bc = _disordered_Bc_selfconsistent(N, mu, eps, g=g, diagnostics=diag)
    except ScfConvergenceError:
        bc = float("nan")
    if np.isfinite(bc):
        shift, failed = 100.0 * (bc - base) / base, False
    else:
        shift, failed = float("nan"), True
    return {"W": float(W), "sc_Bc": float(bc), "baseline_sc_Bc": float(base),
            "shift_pct": float(shift), "convergence_failed": bool(failed)}


# ════════════════════════════════════════════════════════════════════════
# TASK 1.2 / 2 — ensemble averaging: dispatch, summarize, check convergence
# ════════════════════════════════════════════════════════════════════════
def ensemble_run(worker, arg_builder, W: float, base_seed: int,
                 n_real: int, n_workers=None, serial: bool = False,
                 seeds: list | None = None) -> list:
    """Dispatch `worker(*arg_builder(W, seed))` for n_real independent
    realizations via bdg_parallel.parallel_map (process pool, 1-thread BLAS
    per worker -- see bdg_parallel.py). `arg_builder(W, seed) -> tuple` is
    supplied by the CALLER and must match that worker's own positional
    signature exactly -- this replaces an earlier generic
    `fixed_args + (N, W, seed)` concatenation that silently produced the
    WRONG positional order for 3 of the 4 workers (each worker places W and
    seed_seq at a different position, and one caller was missing `mu`
    entirely); explicit per-call builders make the mapping impossible to get
    wrong silently. `seeds` lets convergence_check reuse a single larger
    cumulative seed list."""
    if seeds is None:
        seeds = seed_sequences(base_seed, n_real)
    args_list = [arg_builder(W, sd) for sd in seeds[:n_real]]
    return parallel_map(worker, args_list, n_workers=n_workers, serial=serial)


def summarize_numeric(results: list, key: str) -> dict:
    """Mean/std/SEM across an ensemble for a continuous quantity, excluding
    non-finite (failed/dead-zone) entries but COUNTING them explicitly."""
    raw = [r.get(key, float("nan")) for r in results]
    vals = np.array([v for v in raw if np.isfinite(v)], dtype=float)
    n_failed = sum(1 for v in raw if not np.isfinite(v))
    if vals.size == 0:
        return {"n_real": 0, "n_failed": n_failed, "mean": float("nan"),
                "std": float("nan"), "sem": float("nan"), "distribution": []}
    mean = float(np.mean(vals))
    std = float(np.std(vals, ddof=1)) if vals.size > 1 else 0.0
    sem = std / np.sqrt(vals.size) if vals.size > 1 else float("nan")
    return {"n_real": int(vals.size), "n_failed": n_failed, "mean": mean,
            "std": std, "sem": sem, "distribution": vals.tolist()}


def summarize_fraction(results: list, key: str = "topological") -> dict:
    """Fraction-classified-topological across an ensemble, with the correct
    BERNOULLI spread sqrt(p(1-p)/n) -- NOT a Gaussian std, since the
    underlying per-realization quantity is a +-1 sign, not a continuous
    number."""
    vals = [r[key] for r in results if r.get(key) is not None]
    n = len(vals)
    if n == 0:
        return {"n_real": 0, "frac_topological": float("nan"), "bernoulli_std": float("nan")}
    p = float(sum(1 for v in vals if v) / n)
    std = float(np.sqrt(p * (1.0 - p) / n))
    return {"n_real": n, "frac_topological": p, "bernoulli_std": std}


def convergence_check(worker, arg_builder, W: float,
                      base_seed: int, schedule=(20, 30, 40, 50, 75, 100),
                      n_workers=None, key: str = "shift_pct",
                      is_fraction: bool = False) -> dict:
    """Task 1.2's mandatory convergence-in-realization-count check: draws the
    LARGEST count in `schedule` once via seed_sequences (so every smaller
    count is a strict prefix -- a true nested-ensemble comparison, not
    resampling), then reports mean+spread at each cumulative count and flags
    whether the running mean has stabilized (moved by less than its own
    spread since the previous count). `arg_builder(W, seed) -> tuple` must
    match the worker's positional signature (see ensemble_run)."""
    seeds_all = seed_sequences(base_seed, max(schedule))
    args_all = [arg_builder(W, sd) for sd in seeds_all]
    results_all = parallel_map(worker, args_all, n_workers=n_workers)
    rows = []
    prev_mean = None
    for n in schedule:
        subset = results_all[:n]
        if is_fraction:
            s = summarize_fraction(subset, key)
            m, spread = s["frac_topological"], s["bernoulli_std"]
        else:
            s = summarize_numeric(subset, key)
            m, spread = s["mean"], s["sem"]
        stabilized = (prev_mean is not None and np.isfinite(spread)
                     and abs(m - prev_mean) < spread)
        rows.append({"n_real": n, "mean": m, "spread": spread,
                    "stabilized_vs_prev": bool(stabilized)})
        prev_mean = m
    return {"schedule": rows}


# ════════════════════════════════════════════════════════════════════════
# TASK 3 — realistic disorder scale (verified against primary sources)
# ════════════════════════════════════════════════════════════════════════
# R6-verified induced gap (Chang et al. 2015, arXiv:1411.6255) -- same value
# used throughout R6-R9.
DELTA_EV = 0.19e-3

# Ahn, Pan, Woods, Stanescu, Das Sarma, "Estimating disorder and its adverse
# effects in semiconductor Majorana nanowires", Phys. Rev. Materials 5,
# 124602 (2021), arXiv:2109.00007 -- VERIFIED directly (two independent
# fetches of the ar5iv full text, both returning the identical verbatim
# sentence): "Onsite random potential drawn from an uncorrelated Gaussian
# distribution with zero mean and standard deviation sigma_mu = 2.5 meV",
# on a lattice with lattice constant a = 3.33 nm; separately, transport-based
# broadening estimates from three labs' InAs/InSb samples give Gamma =
# hbar/2tau ~ 2-3 meV. Headline conclusion (verbatim, abstract): "the current
# level of disorder in semiconductor Majorana nanowires is at least an order
# of magnitude higher than that necessary for the emergence of topological
# Majorana zero modes."
SIGMA_MU_MEV_LIT = 2.5           # Ahn et al. 2021, Gaussian on-site fit
GAMMA_BROADENING_MEV_LIT = (2.0, 3.0)   # Ahn et al. 2021, transport broadening

# Two honest conversions from that Gaussian sigma to OUR uniform W (Task 1
# specifies Uniform[-W/2, W/2], the literature source uses a Gaussian):
#   (a) "same energy scale" -- treat sigma_mu directly as our W's natural
#       scale (loose, order-of-magnitude convention, no distribution-shape
#       correction);
#   (b) "variance-matched" -- Var(Uniform[-W/2,W/2]) = W^2/12 = sigma_mu^2
#       => W = sigma_mu*sqrt(12).
# The two give a factor ~3.46x apart; BOTH are reported (never collapsed to
# a single number) because the conversion factor is itself model-dependent
# -- this is the softest number in the whole track (see methods note).
def realistic_W_band_delta_units() -> dict:
    """Realistic W band in units of Delta (R6-verified Delta=0.19 meV), from
    the VERIFIED Ahn et al. 2021 sigma_mu=2.5 meV Gaussian fit, both
    conversions, plus the transport-broadening cross-check."""
    same_scale = SIGMA_MU_MEV_LIT / (DELTA_EV * 1e3)
    variance_matched = (SIGMA_MU_MEV_LIT * np.sqrt(12.0)) / (DELTA_EV * 1e3)
    broadening_lo = GAMMA_BROADENING_MEV_LIT[0] / (DELTA_EV * 1e3)
    broadening_hi = GAMMA_BROADENING_MEV_LIT[1] / (DELTA_EV * 1e3)
    return {
        "same_scale_W_over_Delta": float(same_scale),
        "variance_matched_W_over_Delta": float(variance_matched),
        "broadening_band_W_over_Delta": [float(broadening_lo), float(broadening_hi)],
        "note": ("Ahn et al. 2021 (arXiv:2109.00007) find realistic InAs/Al "
                "disorder is an order of magnitude ABOVE what is needed for "
                "topological MZMs -- this band (~10-46 Delta) reflects that "
                "headline finding, not an arbitrary large-W choice."),
    }


# Sweep grid: weak (well below the gap) -> realistic (literature band above,
# ~11-46 Delta) -> exploratory (further beyond, to see whether the trend
# found in the realistic band saturates or keeps changing).
W_WEAK = (0.1, 0.5, 1.0, 2.0)
W_REALISTIC = (10.0, 15.0, 25.0, 46.0)
W_EXPLORATORY = (60.0, 90.0)


def print_task3():
    band = realistic_W_band_delta_units()
    print("\n=== TASK 3: literature-verified realistic disorder scale ===")
    print(f"Delta = {DELTA_EV*1e3:.3f} meV (Chang 2015, R6-verified)")
    print(f"Ahn, Pan, Woods, Stanescu, Das Sarma, PRMaterials 5, 124602 (2021), "
         f"arXiv:2109.00007: sigma_mu = {SIGMA_MU_MEV_LIT} meV Gaussian on-site fit "
         f"(a=3.33nm), transport broadening Gamma ~ {GAMMA_BROADENING_MEV_LIT} meV")
    print(f"  'same energy scale' convention:      W/Delta = {band['same_scale_W_over_Delta']:.2f}")
    print(f"  'variance-matched uniform' convention: W/Delta = {band['variance_matched_W_over_Delta']:.2f}")
    print(f"  broadening-based band:              W/Delta = {band['broadening_band_W_over_Delta']}")
    print(f"W_weak (exploratory-small)  = {W_WEAK}")
    print(f"W_realistic (literature band) = {W_REALISTIC}")
    print(f"W_exploratory (beyond realistic) = {W_EXPLORATORY}")
    return band


# ════════════════════════════════════════════════════════════════════════
# TASK 2 — parallelization benchmark + pre-full-run time estimate
# ════════════════════════════════════════════════════════════════════════
def benchmark_disorder(N: int = 80, mu: float = 1.2, n_tasks: int = 12,
                       W: float = 15.0, base_seed: int = 12345) -> dict:
    """Task 2.1/2.3: serial-vs-parallel equivalence + BLAS-oversubscription
    benchmark, on BOTH the cheap engine (subgap detector, one eigh) and the
    expensive engine (self-consistent Bc, many eigh per root search x several
    B points) -- since Task 4.3 uses the expensive engine, sizing the run
    from the cheap engine alone would be misleading."""
    out = {}
    B = np.sqrt(mu ** 2 + 1.0) - 0.05
    seeds = seed_sequences(base_seed, n_tasks)

    cheap_args = [(N, mu, B, W, sd, 0.25, 1.0) for sd in seeds]
    print(f"  [benchmark] cheap-engine equivalence check ({len(cheap_args)} tasks)...", flush=True)
    out["equivalence_cheap"] = check_equivalence(_subgap_worker, cheap_args)

    # Pre-warm the clean SC baseline cache ONCE (outside all timings below) so
    # every timed call measures only the disordered solve, not a redundant
    # ~30s clean re-solve -- see _sc_baseline().
    _sc_baseline(N, mu, G_R7, 1.0)
    n_exp = min(3, n_tasks)   # the expensive engine is ~100s/task; keep the
                              # benchmark itself tractable (measured, not guessed)
    expensive_args = [(N, mu, W, sd, G_R7, 1.0) for sd in seeds[:n_exp]]
    print(f"  [benchmark] expensive-engine equivalence check ({len(expensive_args)} tasks, "
         f"~100s/task -> ~{2*len(expensive_args)*100/60:.0f} min)...", flush=True)
    out["equivalence_expensive"] = check_equivalence(_scshift_worker, expensive_args)

    def _time(fn):
        t0 = time.perf_counter()
        fn()
        return time.perf_counter() - t0

    for label, worker, args in (("cheap_subgap", _subgap_worker, cheap_args),
                                ("expensive_scshift", _scshift_worker, expensive_args)):
        print(f"  [benchmark] timing {label} ({len(args)} tasks): serial 1-thread...", flush=True)
        _pin_blas(1)
        t_serial1 = _time(lambda: [worker(*a) for a in args])
        print(f"    serial_1thread_s = {t_serial1:.1f}", flush=True)
        _pin_blas(os.cpu_count() or 8)
        print(f"  [benchmark] timing {label}: serial default-BLAS...", flush=True)
        t_serialdef = _time(lambda: [worker(*a) for a in args])
        print(f"    serial_defaultBLAS_s = {t_serialdef:.1f}", flush=True)
        _pin_blas(1)
        nw = default_workers()
        print(f"  [benchmark] timing {label}: parallel ({nw} workers, 1-thread BLAS)...", flush=True)
        t_par = _time(lambda: parallel_map(worker, args, n_workers=nw, blas_threads=1))
        print(f"    parallel_capped_s = {t_par:.1f}", flush=True)
        out[label] = {"serial_1thread_s": t_serial1, "serial_defaultBLAS_s": t_serialdef,
                     "parallel_capped_s": t_par, "n_workers": nw, "n_tasks": len(args),
                     "per_task_s_serial1": t_serial1 / len(args),
                     "per_task_s_parallel": t_par / len(args)}
    return out


def estimate_full_run(per_task_s: float, n_real: int, n_grid_points: int,
                      n_workers: int) -> dict:
    """Task 2.4: honest pre-full-run wall-clock estimate, to be printed and
    inspected BEFORE the full Task 4 sweep is launched (not discovered
    mid-run)."""
    total_tasks = int(n_real) * int(n_grid_points)
    est_wall_s = total_tasks * per_task_s / max(1, n_workers)
    return {"total_tasks": total_tasks, "per_task_s": per_task_s,
            "n_workers": n_workers, "est_wall_s": est_wall_s,
            "est_wall_min": est_wall_s / 60.0}


# ════════════════════════════════════════════════════════════════════════
# TASK 4.1 — does the ensemble-averaged / fraction-topological invariant
# smear smoothly, sharply, or erratically as W grows?
# ════════════════════════════════════════════════════════════════════════
def boundary_smearing_vs_W(N: int = 80, mu: float = 1.2, n_real: int = 50,
                           Ws=None, base_seed: int = 2026, n_workers=None,
                           convergence_schedule=(20, 30, 40, 50, 75, 100)
                           ) -> dict:
    if Ws is None:
        Ws = list(W_WEAK) + list(W_REALISTIC) + list(W_EXPLORATORY)
    B = np.sqrt(mu ** 2 + 1.0) - 0.05     # near the clean boundary, off it
    # measured (see benchmark_disorder): this engine costs ~0.05s/task, and
    # process-pool dispatch (~0.2s/task overhead) is SLOWER than serial for
    # tasks this cheap -- so this driver always runs serial regardless of
    # n_workers, per the same "small/cheap loops stay serial" rule bdg_parallel
    # documents (parallel_sweeps memory).
    # _topo_worker(N, mu, B, W, seed_seq, B_deep_topo=None, B_ref_trivial=0.0, delta=1.0)
    build = lambda W_, sd: (N, mu, B, W_, sd)
    rows = []
    for W in Ws:
        results = ensemble_run(_topo_worker, build, W, base_seed, n_real,
                               serial=True)
        s = summarize_fraction(results, "topological")
        s["W"] = float(W)
        s["is_realistic"] = bool(W in W_REALISTIC)
        rows.append(s)
    conv = convergence_check(_topo_worker, build, W_REALISTIC[0], base_seed,
                            schedule=convergence_schedule, n_workers=1,
                            key="topological", is_fraction=True)
    return {"N": N, "mu": mu, "B": float(B), "rows": rows, "convergence": conv}


# ════════════════════════════════════════════════════════════════════════
# TASK 4.2 — spurious sub-gap states from disorder (R4 threshold-artifact
# connection)
# ════════════════════════════════════════════════════════════════════════
def subgap_states_vs_W(N: int = 80, mu: float = 1.2, n_real: int = 50,
                       Ws=None, threshold: float = 0.25,
                       base_seed: int = 2027, n_workers=None,
                       convergence_schedule=(20, 30, 40, 50, 75, 100)
                       ) -> dict:
    if Ws is None:
        Ws = list(W_WEAK) + list(W_REALISTIC) + list(W_EXPLORATORY)
    B = np.sqrt(mu ** 2 + 1.0) - 0.05
    # measured ~0.05s/task (benchmark_disorder, "cheap_subgap"): serial beats
    # pooled dispatch here too -- see the identical note in
    # boundary_smearing_vs_W above.
    # _subgap_worker(N, mu, B, W, seed_seq, threshold=0.25, delta=1.0)
    build = lambda W_, sd: (N, mu, B, W_, sd, threshold)
    rows = []
    for W in Ws:
        results = ensemble_run(_subgap_worker, build, W, base_seed, n_real,
                               serial=True)
        s_n = summarize_numeric(results, "n_subgap")
        s_e = summarize_numeric(results, "min_abs_E")
        rows.append({"W": float(W), "is_realistic": bool(W in W_REALISTIC),
                    "n_subgap_mean": s_n["mean"], "n_subgap_std": s_n["std"],
                    "n_subgap_n_real": s_n["n_real"],
                    "min_abs_E_mean": s_e["mean"], "min_abs_E_std": s_e["std"]})
    conv = convergence_check(_subgap_worker, build, W_REALISTIC[0], base_seed,
                            schedule=convergence_schedule, n_workers=1,
                            key="n_subgap")
    return {"N": N, "mu": mu, "B": float(B), "threshold": threshold,
            "rows": rows, "convergence": conv}


# ════════════════════════════════════════════════════════════════════════
# TASK 4.3 — does the R4/R5/R7 inductive-feedback null survive realistic W?
# ════════════════════════════════════════════════════════════════════════
def null_result_vs_W(N: int = 80, mu: float = 1.2, g: float = G_R7,
                     n_real: int = 30, Ws=None, base_seed: int = 2028,
                     n_workers=None,
                     convergence_schedule=(10, 15, 20, 25, 30)) -> dict:
    """Task 4.3: ensemble-averaged self-consistent boundary shift under
    disorder, at the SAME g=G_R7 (R7's validated subcritical feedback
    coupling) used throughout R4-R7-R9T2's null-result lineage. Reports the
    convergence-failure fraction explicitly (per module docstring/risk #2 --
    a non-converged realization is a reported outcome, not silently
    dropped).

    POST-HOC CORRECTION (superseding the original n_real=20 scope note): a
    subsequent review of the FIRST --full run's r9_disorder_data.json found
    that its Task 4.3 rows carried n_real=1-4 (out of 20 attempted) at
    realistic W -- i.e. the reported means there were NOT valid ensemble
    averages, contradicting this track's own mandatory-ensemble-averaging
    rule. Root-cause investigation (see r9_disorder_methods.md) found the
    dominant cause was a SEARCH-WINDOW BUG, not a physical self-consistency
    instability: _disordered_Bc_selfconsistent centered its coarse B-scan on
    the fixed clean `ana = sqrt(mu^2+delta^2)`, but Task 4.1 independently
    showed the RIGID boundary itself is displaced from `ana` by a mean of
    +0.45 to +1.0 (up to +2.9) at realistic W -- far outside the solver's old
    +-0.20 search ceiling. Fixed by centering on the per-realization rigid
    boundary instead (see _disordered_Bc_selfconsistent's docstring). This
    raised the realistic-W success rate from ~5-20% to ~35-100% in post-fix
    sampling; a smaller, genuine residual self-consistency-loop failure rate
    remains and is reported honestly, not eliminated.

    Given the corrected solver, per-task cost is now dominated by
    successful (not fast-failing) solves, so it did NOT drop after the fix
    -- see the re-run benchmark in the methods note for the honest updated
    pre-run estimate. n_real is raised from the original 20 to 30 (up from
    the Task-1 floor, though still short of Tasks 4.1/4.2's n=50 -- an
    explicitly stated, compute-driven limitation, not a silent one) and the
    convergence-check schedule is capped at 30 (vs 4.1/4.2's ceiling of 100)
    for the same reason. The W grid remains narrowed to one representative
    weak point (0, 1) plus the full realistic band (4.1/4.2 already cover
    the full weak grid cheaply)."""
    if Ws is None:
        Ws = (0.0, 1.0) + tuple(W_REALISTIC)
    # _scshift_worker(N, mu, W, seed_seq, g=G_R7, delta=1.0)
    build = lambda W_, sd: (N, mu, W_, sd, g)
    rows = []
    for W in Ws:
        n_this = n_real if W > 0 else 5    # W=0 is deterministic; 5 is a
                                            # sanity replicate, not an ensemble
        t0 = time.time()
        results = ensemble_run(_scshift_worker, build, W, base_seed, n_this,
                               n_workers=n_workers)
        print(f"    [4.3] W={W:.2f} done ({n_this} reps, {time.time()-t0:.0f}s)", flush=True)
        s = summarize_numeric(results, "shift_pct")
        n_fail = sum(1 for r in results if r.get("convergence_failed"))
        s["W"] = float(W)
        s["is_realistic"] = bool(W in W_REALISTIC)
        s["n_convergence_failed"] = n_fail
        rows.append(s)
    conv = convergence_check(_scshift_worker, build, W_REALISTIC[0], base_seed,
                            schedule=convergence_schedule, n_workers=n_workers,
                            key="shift_pct")
    return {"N": N, "mu": mu, "g": g, "rows": rows, "convergence": conv}


# ════════════════════════════════════════════════════════════════════════
# main / drivers
# ════════════════════════════════════════════════════════════════════════
def write_summary(data: dict) -> str:
    jpath = _out("r9_disorder_data.json")
    with open(jpath, "w") as fh:
        json.dump(data, fh, indent=2, default=lambda o: float(o)
                  if isinstance(o, (np.floating,)) else
                  (o.tolist() if isinstance(o, np.ndarray) else str(o)))
    return jpath


def make_figures(data: dict) -> list:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths = []

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))

    b = data["task4_1_boundary_smearing"]["rows"]
    Ws = [r["W"] for r in b]
    fr = [r["frac_topological"] for r in b]
    er = [r["bernoulli_std"] for r in b]
    colors = ["tab:blue" if r["is_realistic"] else "tab:gray" for r in b]
    axs[0].errorbar(Ws, fr, yerr=er, fmt="none", ecolor="gray", capsize=2)
    axs[0].scatter(Ws, fr, c=colors, zorder=3)
    axs[0].set_xscale("symlog")
    axs[0].set_xlabel(r"$W$ ($\Delta$)"); axs[0].set_ylabel("fraction topological")
    axs[0].set_title("Task 4.1: boundary smearing vs W")

    s = data["task4_2_subgap_states"]["rows"]
    Ws2 = [r["W"] for r in s]
    nmean = [r["n_subgap_mean"] for r in s]
    nstd = [r["n_subgap_std"] for r in s]
    axs[1].errorbar(Ws2, nmean, yerr=nstd, fmt="o-", ms=4)
    axs[1].set_xscale("symlog")
    axs[1].set_xlabel(r"$W$ ($\Delta$)"); axs[1].set_ylabel(r"$\langle$ n sub-gap states $\rangle$")
    axs[1].set_title("Task 4.2: spurious sub-gap states vs W")

    n = data["task4_3_null_result"]["rows"]
    Ws3 = [r["W"] for r in n]
    sh = [r["mean"] for r in n]
    se = [r["sem"] if np.isfinite(r["sem"]) else 0.0 for r in n]
    axs[2].errorbar(Ws3, sh, yerr=se, fmt="o-", ms=4)
    axs[2].axhline(0, color="k", lw=0.6)
    axs[2].set_xlabel(r"$W$ ($\Delta$)"); axs[2].set_ylabel(r"$\Delta B_c/B_c$ (%)")
    axs[2].set_title("Task 4.3: inductive-feedback null vs W")

    fig.tight_layout()
    p = _out("r9_disorder_summary.png")
    fig.savefig(p, dpi=130)
    plt.close(fig)
    paths.append(p)
    return paths


def run_full(N: int = 80, n_real: int = 50, n_real_sc: int = 20,
            n_workers=None, quick: bool = False) -> dict:
    t_all = time.time()
    if quick:
        # NB: keep N=80 even in --quickfull. The clean self-consistent solver
        # (R5/R7) has known N-dependent dead-zones off N=80 at this (mu,g)
        # (memory: "clean floor N~80") -- confirmed here directly: N=60
        # raises ScfConvergenceError on the CLEAN baseline itself. Shrinking
        # only the realization counts keeps the smoke test on solid ground.
        n_real, n_real_sc = 8, 4

    data = {}
    print("\n=== W=0 VALIDATION (MUST PASS before any W>0 result) ===")
    data["validation"] = validate_W0(N=max(N, 80), mu=1.2, verbose=True)
    if not data["validation"]["ALL_PASS"]:
        raise RuntimeError("W=0 validation FAILED -- refusing to run W>0 disorder "
                          "experiments until every disordered shell reproduces "
                          "its clean counterpart exactly.")

    print("\n=== TASK 3: literature-verified disorder scale ===")
    data["task3_literature"] = print_task3()

    print("\n=== TASK 2: parallelization benchmark + time estimate ===")
    bench = benchmark_disorder(N=N, n_tasks=12 if not quick else 6)
    data["benchmark"] = bench
    print(json.dumps(bench, indent=2, default=str))
    n_grid_cheap = len(W_WEAK) + len(W_REALISTIC) + len(W_EXPLORATORY)
    n_grid_sc = 2 + len(W_REALISTIC)     # null_result_vs_W's reduced default grid
    # Task 4.1/4.2 drivers dispatch SERIAL (measured: pool overhead > per-task
    # cost for this cheap engine -- see boundary_smearing_vs_W/subgap_states_vs_W),
    # so the honest estimate uses the SERIAL per-task cost, not the (slower,
    # pool-overhead-dominated) parallel number.
    est_cheap = estimate_full_run(bench["cheap_subgap"]["per_task_s_serial1"],
                                  n_real, n_grid_cheap * 2, 1)
    est_sc = estimate_full_run(bench["expensive_scshift"]["per_task_s_parallel"],
                               n_real_sc, n_grid_sc, bench["expensive_scshift"]["n_workers"])
    data["pre_run_estimate"] = {"cheap_tasks_4_1_and_4_2": est_cheap,
                                "expensive_task_4_3": est_sc,
                                "total_est_wall_min": est_cheap["est_wall_min"] + est_sc["est_wall_min"]}
    print(f"PRE-RUN ESTIMATE: Task4.1+4.2 ~{est_cheap['est_wall_min']:.1f} min, "
         f"Task4.3 ~{est_sc['est_wall_min']:.1f} min, "
         f"TOTAL ~{data['pre_run_estimate']['total_est_wall_min']:.1f} min")

    print("\n=== TASK 4.1: boundary smearing vs W ===")
    data["task4_1_boundary_smearing"] = boundary_smearing_vs_W(
        N=N, n_real=n_real, n_workers=n_workers,
        convergence_schedule=(20, 30, 40, 50) if quick else (20, 30, 40, 50, 75, 100))
    for r in data["task4_1_boundary_smearing"]["rows"]:
        print(f"  W={r['W']:.2f}: frac_topo={r['frac_topological']:.3f}"
             f"+-{r['bernoulli_std']:.3f} (n={r['n_real']})"
             f"{' [realistic]' if r['is_realistic'] else ''}")

    print("\n=== TASK 4.2: spurious sub-gap states vs W ===")
    data["task4_2_subgap_states"] = subgap_states_vs_W(
        N=N, n_real=n_real, n_workers=n_workers,
        convergence_schedule=(20, 30, 40, 50) if quick else (20, 30, 40, 50, 75, 100))
    for r in data["task4_2_subgap_states"]["rows"]:
        print(f"  W={r['W']:.2f}: <n_subgap>={r['n_subgap_mean']:.2f}"
             f"+-{r['n_subgap_std']:.2f} (n={r['n_subgap_n_real']})"
             f"{' [realistic]' if r['is_realistic'] else ''}")

    print("\n=== TASK 4.3: inductive-feedback null vs W ===")
    data["task4_3_null_result"] = null_result_vs_W(
        N=N, n_real=n_real_sc, n_workers=n_workers,
        convergence_schedule=(5, 10, 15, 20) if not quick else (5, 8))
    for r in data["task4_3_null_result"]["rows"]:
        print(f"  W={r['W']:.2f}: shift={r['mean']:+.4f}%+-{r['sem']:.4f}"
             f" (n={r['n_real']}, failed={r['n_convergence_failed']})"
             f"{' [realistic]' if r['is_realistic'] else ''}")

    data["wall_time_s"] = time.time() - t_all
    jpath = write_summary(data)
    figs = make_figures(data)
    print(f"\nWrote {jpath}\n      {figs}\nTotal wall: {data['wall_time_s']:.1f}s")
    return data


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true", help="W=0 gate only")
    ap.add_argument("--literature", action="store_true", help="Task 3 report only")
    ap.add_argument("--benchmark", action="store_true",
                    help="Task 2 parallel/BLAS benchmark only (surfaces the "
                        "pre-full-run time estimate WITHOUT running Task 4)")
    ap.add_argument("--full", action="store_true", help="full Task 3+4 run")
    ap.add_argument("--quickfull", action="store_true", help="full run, tiny grid/N")
    ap.add_argument("--n-real", type=int, default=50)
    ap.add_argument("--n-real-sc", type=int, default=30)
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    if args.literature or args.validate or args.full or args.quickfull:
        print_task3()
    if args.validate:
        validate_W0(N=80, mu=1.2)
    if args.benchmark:
        b = benchmark_disorder()
        print(json.dumps(b, indent=2, default=str))
        est = estimate_full_run(b["expensive_scshift"]["per_task_s_parallel"],
                               args.n_real_sc, 1 + len(W_WEAK) + len(W_REALISTIC),
                               b["expensive_scshift"]["n_workers"])
        print(f"Pre-run estimate for Task 4.3 alone: {est['est_wall_min']:.1f} min")
    if args.full or args.quickfull:
        run_full(n_real=args.n_real, n_real_sc=args.n_real_sc,
                n_workers=args.workers, quick=args.quickfull)
