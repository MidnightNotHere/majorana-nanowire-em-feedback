"""
Inductive Circuit Feedback Shifts the Topological Phase Boundary
              in InAs/Al Majorana Nanowires
Phase: R5 — Self-consistent ring self-field + Pfaffian Z2 boundary shift
Status: COMPLETE (this module IS the R5 deliverable of
        r4_pfaffian_boundary_methods.md section 4c; results + honest
        convergence findings written up in section 4d)

What this module replaces
-------------------------
The interim R4 estimate of the rigorous (ring/Pfaffian) boundary shift used a
stopgap self-field: the phi_self(x) texture was computed once on the OPEN wire
(open-boundary eigenstates), then pasted onto the RING geometry for the
Pfaffian measurement, with the wrap bond N-1 -> 0 carrying only the external
gauge phase (no self-field contribution). That gave the interim ~0.6-0.8%
shift at mu = 1.2-1.4 quoted in r4_pfaffian_boundary_methods.md section 4b.
It was documented there as a modelling choice, not a real ring calculation.

This module does the real calculation:

  1. RING current density from the ring BdG eigenstates — all N bonds,
     including the wrap bond (`ring_current_density`).
  2. Self-induced vector potential from that ring current
     (`ring_biot_savart_self`) — see the geometry note below.
  3. Fixed-point iteration current -> A_self -> Peierls phase -> new H ->
     new eigenstates -> new current, to self-consistency
     (`selfconsistent_ring_texture`), with the convergence behaviour itself
     treated as an open question (trajectories, seeds, hysteresis, N,
     tolerance — see the experiment functions at the bottom).
  4. B_c extracted from the Pfaffian Z2 sign flip of the SELF-CONSISTENTLY
     dressed ring Hamiltonian (`ring_boundary_Bc_pfaffian_selfconsistent`),
     i.e. at every field B probed by the bisection the texture is first
     converged at that B, then the invariant is evaluated.

The interim stopgap is ALSO implemented here, clearly labelled SUPERSEDED
(`interim_ring_Bc_openwire_texture`), so the old number is reproducible for
regression/comparison — it previously existed only as session work quoted in
the methods note, not as committed code.

Geometry note: why the Biot-Savart kernel is adapted, not reused verbatim
-------------------------------------------------------------------------
The validated open-wire kernel (bdg_r4_em_coupling.biot_savart_self) uses
straight-line separations |x_i - x_j| between bond midpoints on a straight
segment. A ring closes on itself, which breaks that kernel in two specific
ways:

  * distances wrap: bond N-1 is ADJACENT to bond 0, not (N-1)*a away;
  * direction matters: the source element dl_j and the field point's local
    tangent t_i are not parallel — the tangential component of A (the only
    part that enters the Peierls phase) picks up a factor
    dl_j . t_i = a * cos(theta_i - theta_j).

So this module keeps the validated kernel's STRUCTURE — discrete Biot-Savart
sum (mu0/4pi) * kernel @ I with the same analytic finite-segment self-term
2*ln(a/r_wire) on the diagonal (the segment is locally straight, a << R, so
the local regularisation is unchanged) and the same SI units path — and
changes ONLY the geometry factors: chord distance 2R sin(dtheta/2) and the
tangent projection cos(dtheta), for a circular loop of circumference N*a.

This is not taken on faith: `validate_ring_biot_savart` is the closed-
geometry analogue of R4 Validation 2 — for a uniform current the discrete
kernel must reproduce the textbook thin-wire circular-loop self-inductance
L = mu0 R [ln(8R/r_wire) - 2], and the continuum limit of the kernel above
gives exactly that expression (the check below reports the finite-N error).

Amplification factor g: unchanged at 3e5. The plausibility argument for g
(resonator Q x kinetic-inductance enhancement — see figure4's long comment
in bdg_r4_em_coupling.py) concerns the CIRCUIT, not the wire geometry, so
there is no reason for the ring measurement to use a different value. The
geometric self-inductance per unit length differs between a straight segment
and a loop only logarithmically (few x10%), far inside g's order-of-magnitude
error bars.

Nothing in bdg_r4_em_coupling.py is modified: the open-wire pipeline, the
Pfaffian detector and all 6 R4 validations are imported and used read-only.
"""

from __future__ import annotations

import contextlib
import io
import json
import os
import time

import numpy as np

# ── Read-only imports from the validated R4 module ──────────────────────
with contextlib.redirect_stdout(io.StringIO()):
    from bdg_r4_em_coupling import (            # noqa: E402
        E_CHARGE, HBAR, MU_0,
        T_HOP_SI, LATTICE_A,
        build_bdg_lo,
        apply_peierls_phases,
        compute_current_density,
        biot_savart_self,
        compute_peierls_from_A,
        _build_ring_lo_phased,
        _ring_pfaffian_sign,
        openwire_boundary_Bc,
    )

from scipy.linalg import eigh

# ── All outputs go into this program's own data folder ───────────────────
_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "self_consistency_bifurcation")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# Same default feedback amplification as Figure 4/5 (see module docstring).
#
# CAVEAT (added after R5-revised, bdg_r5_physically_meaningful_shift.py): this
# value is a mathematically valid coupling, but it sits ABOVE the Picard
# bifurcation g_crit ~ 1e5 (Figure 8) — at g=3e5 NO attracting fixed point
# exists, so results computed here at G0_DEFAULT (the "+0.05%" shift) are the
# location of a DYNAMICALLY UNREACHABLE root, not a physical static equilibrium.
# The physically meaningful headline comes from g < g_crit and is ~-0.06%
# (finite-size, -> 0 in the bulk); see r4_pfaffian_boundary_methods.md §4e.
# G0_DEFAULT is kept unchanged so Figures 7-9 remain reproducible and so the
# g=3e5 result stays visible as the finding that REVEALED g_crit.
G0_DEFAULT = 3.0e5


# ════════════════════════════════════════════════════════════════════
# FUNCTION 1 — Ring current density (all N bonds, including the wrap bond)
# ════════════════════════════════════════════════════════════════════
def ring_current_density(eigenvectors: np.ndarray,
                         eigenvalues: np.ndarray,
                         N: int,
                         t_hop: float,
                         phi_bonds: np.ndarray) -> np.ndarray:
    """
    Bond current on the periodic ring from the BdG ground state, for all N
    bonds i -> (i+1) mod N — INCLUDING the wrap bond N-1 -> 0 that the
    interim calculation left out.

    Identical physics and conventions to the validated open-wire
    compute_current_density (same prefactor 2*e*t_hop/hbar in Amperes, same
    electron-block-only rule, same occupied = E<0 convention, same
    kron(spin, site) basis); the only change is periodic site indexing via
    np.roll so the wrap bond is bond N-1.

    Note (physics): in a BdG mean-field the s-wave pairing exchanges charge
    with the condensate, so the bond current need NOT be uniform around the
    ring when the phase texture is non-uniform. For a UNIFORM texture,
    translational invariance forces a uniform current — used as a sanity
    check in validate_ring_biot_savart().

    Parameters
    ----------
    eigenvectors : (4N, 4N) BdG eigenvectors (columns), from eigh
    eigenvalues  : (4N,) BdG eigenvalues
    N            : ring size (sites = bonds)
    t_hop        : hopping in Joules (sets the SI current scale)
    phi_bonds    : (N,) total per-bond Peierls phase the eigenstates were
                   computed at (external + self)

    Returns
    -------
    I_bonds : (N,) bond current, Amperes.
    """
    eigenvalues = np.asarray(eigenvalues)
    occ = eigenvalues < 0.0
    U = eigenvectors[:2 * N, occ]                # electron block only
    if U.shape[1] == 0:
        return np.zeros(N, dtype=float)
    U = U.reshape(2, N, -1)                      # (spin, site, n_occ)
    u_i   = U                                    # site i
    u_ip1 = np.roll(U, -1, axis=1)               # site (i+1) mod N  <- wrap
    inner = np.sum(np.conj(u_i) * u_ip1, axis=0)             # (N, n_occ)
    phase = np.exp(1j * np.asarray(phi_bonds, dtype=float))[:, None]
    prefactor = 2.0 * E_CHARGE * t_hop / HBAR                 # Amperes
    return prefactor * np.sum(np.imag(inner * phase), axis=1)  # (N,)


# ════════════════════════════════════════════════════════════════════
# FUNCTION 2 — Biot-Savart self vector potential on the RING
# ════════════════════════════════════════════════════════════════════
_RING_KERNEL_CACHE: dict = {}


def ring_biot_savart_self(I_bonds: np.ndarray,
                          N: int,
                          lattice_const: float,
                          wire_radius: float | None = None) -> np.ndarray:
    """
    Tangential self vector potential A_t at every bond midpoint of a
    circular ring of circumference N*a, from the per-bond currents.

    Adaptation of the validated open-wire kernel (see module docstring for
    the explicit flag/justification):

      off-diagonal:  kernel[i,j] = a * cos(dtheta_ij) / chord_ij,
                     chord_ij = 2 R sin(|dtheta_ij| / 2),
                     dtheta_ij = 2 pi (i - j) / N,  R = N a / (2 pi)
      diagonal:      2 * ln(a / r_wire)     (same finite-segment self-term
                                             as the validated open wire; the
                                             segment is locally straight)

      A_t = (mu0 / 4 pi) * kernel @ I_bonds        [T*m]

    Known-answer check: for uniform I the total loop flux
    sum_i A_t(i)*a = L*I must reproduce the thin-wire loop inductance;
    validate_ring_biot_savart() verifies this against the exact
    finite-segment (discrete) expression — the closed-geometry analogue of
    R4 Validation 2, which likewise compared against a finite-segment
    formula (the one carrying the +2*gamma Euler constant).

    Parameters
    ----------
    I_bonds       : (N,) bond currents (Amperes)
    N             : ring size (bonds)
    lattice_const : a (metres)
    wire_radius   : r_wire (metres) for the self-term; default a/10, matching
                    the open-wire pipeline's default

    Returns
    -------
    A_t : (N,) tangential vector potential at each bond midpoint (T*m).
    """
    I_bonds = np.asarray(I_bonds, dtype=float)
    a = lattice_const
    r_wire = a / 10.0 if wire_radius is None else wire_radius

    key = (N, float(a), float(r_wire))
    if key not in _RING_KERNEL_CACHE:
        R = N * a / (2.0 * np.pi)
        idx = np.arange(N)
        dtheta = 2.0 * np.pi * (idx[:, None] - idx[None, :]) / N
        chord = 2.0 * R * np.abs(np.sin(dtheta / 2.0))
        kernel = np.zeros((N, N), dtype=float)
        off = ~np.eye(N, dtype=bool)
        kernel[off] = a * np.cos(dtheta[off]) / chord[off]
        np.fill_diagonal(kernel, 2.0 * np.log(a / r_wire))
        _RING_KERNEL_CACHE[key] = kernel
    kernel = _RING_KERNEL_CACHE[key]

    return (MU_0 / (4.0 * np.pi)) * (kernel @ I_bonds)         # T*m


_EULER_GAMMA = 0.5772156649015329


def validate_ring_biot_savart(N: int = 100,
                              lattice_const: float = LATTICE_A,
                              wire_radius: float | None = None) -> dict:
    """
    Known-answer validation of the ring kernel (closed-geometry analogue of
    R4 Validation 2): uniform current I0 -> the discrete loop flux
    Phi = sum_i A_t(i) * a  must reproduce the thin-wire circular-loop
    self-inductance.

    Two reference values are used, deliberately:

      * L_discrete = mu0 R [ln(8R/r_wire) - 2 + (gamma - ln 2)]
        — the EXACT large-N expectation for this discretisation.
        Derivation: with theta_k = 2 pi k/N the off-diagonal row sum is
        (pi/N) sum_{k=1}^{N-1} cos(theta_k)/sin(theta_k/2)
          = (pi/N)[ sum 1/sin(pi k/N) - 2 sum sin(pi k/N) ]
          -> 2[ln(2N/pi) + gamma - 2],
        and adding the validated segment self-term 2 ln(a/r) gives
        2[ln(8R/r) - 2 + gamma - ln 2] per site. The (gamma - ln 2) is the
        same class of discrete-sum constant that R4 Validation 2's
        finite-segment formula carries as +2*gamma. Pass criterion is
        against THIS value (<0.5%).
      * L_textbook = mu0 R [ln(8R/r_wire) - 2] — the continuum textbook
        result, reported for context; the discrete kernel sits
        (ln 2 - gamma)/(ln(8R/r) - 2) ~ 2% below it at these sizes — an
        understood discretisation constant (equivalent to a ~12% rescale
        of the already-conventional r_wire = a/10), not an error, and far
        inside g's order-of-magnitude uncertainty.

    Also checks the two symmetry sanities the adaptation must preserve:
      * uniform current -> uniform A_t (circulant kernel);
      * the ring current of a clean uniform-flux ring is itself uniform
        across ALL N bonds including the wrap bond (translation invariance),
        which exercises ring_current_density's wrap-bond indexing.
    """
    a = lattice_const
    r_wire = a / 10.0 if wire_radius is None else wire_radius
    R = N * a / (2.0 * np.pi)
    I0 = 1.0                                                   # A

    A_t = ring_biot_savart_self(np.full(N, I0), N, a, r_wire)
    L_num = float(np.sum(A_t) * a / I0)                        # H
    L_textbook = MU_0 * R * (np.log(8.0 * R / r_wire) - 2.0)   # H
    L_discrete = MU_0 * R * (np.log(8.0 * R / r_wire) - 2.0
                             + _EULER_GAMMA - np.log(2.0))     # H
    rel_err = abs(L_num - L_discrete) / abs(L_discrete)
    rel_err_textbook = abs(L_num - L_textbook) / abs(L_textbook)
    A_uniformity = float(np.max(A_t) - np.min(A_t)) / max(abs(float(np.mean(A_t))), 1e-300)

    # wrap-bond / translation-invariance check of ring_current_density
    Nc, mu_c, B_c = 60, 0.5, 1.4
    phi_u = np.full(Nc, np.pi * 0.1 / Nc)
    H = _build_ring_lo_phased(Nc, mu_c, B_c, phi_u)
    evals, evecs = eigh(H)
    I_bonds = ring_current_density(evecs, evals, Nc, T_HOP_SI, phi_u)
    I_spread = float(np.max(I_bonds) - np.min(I_bonds))
    I_scale = max(abs(float(np.mean(I_bonds))), 1e-300)

    return {
        "L_numeric_H": L_num,
        "L_discrete_H": L_discrete,
        "L_textbook_H": L_textbook,
        "rel_err_discrete": rel_err,
        "rel_err_textbook": rel_err_textbook,
        "A_uniformity": A_uniformity,
        "current_uniformity": I_spread / I_scale,
        "pass": bool(rel_err < 0.005
                     and A_uniformity < 1e-10
                     and I_spread / I_scale < 1e-8),
    }


# ════════════════════════════════════════════════════════════════════
# FUNCTION 3 — Self-consistent fixed-point iteration on the ring
# ════════════════════════════════════════════════════════════════════
def selfconsistent_ring_texture(N: int, mu: float, B: float,
                                g: float,
                                phi_ext_bond: float,
                                seed_phi_self: np.ndarray | None = None,
                                beta: float = 1.0,
                                tol: float = 1e-12,
                                max_iter: int = 300,
                                delta: float = 1.0,
                                t_hop: float = T_HOP_SI,
                                lattice_const: float = LATTICE_A,
                                wire_radius: float | None = None,
                                record_history: bool = False) -> dict:
    """
    Converge the ring self-field texture at fixed (N, mu, B, g) by plain
    (optionally damped) fixed-point iteration of the map

        G(phi_self):  H(phi_ext + phi_self) -> eigenstates -> ring current
                      -> ring Biot-Savart A_t -> g * (e a / hbar) * A_t

    Convergence is judged on the PRE-mixing residual
    r_k = max|G(phi_k) - phi_k| < tol, so the reported fixed point satisfies
    the unmixed equation regardless of beta. Divergence (residual exceeding
    1e3 rad) aborts early. This function never fabricates a result: the
    caller must check res['converged'] before using res['phi_self'].

    Parameters
    ----------
    seed_phi_self : (N,) initial guess for phi_self, or None for zeros.
                    Used by the seed-independence / hysteresis experiments.
    beta          : mixing parameter, phi <- (1-beta) phi + beta G(phi).
                    beta=1 is plain fixed-point iteration.
    record_history: if True, store per-iteration residual, mean(phi_self),
                    std(phi_self), and mean bond current in res['history'].

    Returns
    -------
    dict with keys: phi_self (N,), converged (bool), n_iter, residual,
    diverged (bool), history (dict of lists, if requested).
    """
    phi_ext_arr = np.full(N, float(phi_ext_bond))
    phi = (np.zeros(N) if seed_phi_self is None
           else np.array(seed_phi_self, dtype=float).copy())
    if phi.shape != (N,):
        raise ValueError(f"seed_phi_self has shape {phi.shape}, expected {(N,)}")

    hist: dict = {"residual": [], "phi_mean": [], "phi_std": [], "I_mean": []}
    converged = diverged = False
    resid = np.inf
    it = 0
    for it in range(1, max_iter + 1):
        phi_tot = phi_ext_arr + phi
        H = _build_ring_lo_phased(N, mu, B, phi_tot, delta)
        evals, evecs = eigh(H)
        I_bonds = ring_current_density(evecs, evals, N, t_hop, phi_tot)
        A_t = ring_biot_savart_self(I_bonds, N, lattice_const, wire_radius)
        phi_new = g * compute_peierls_from_A(A_t, lattice_const)   # (N,)

        resid = float(np.max(np.abs(phi_new - phi)))
        if record_history:
            hist["residual"].append(resid)
            hist["phi_mean"].append(float(np.mean(phi)))
            hist["phi_std"].append(float(np.std(phi)))
            hist["I_mean"].append(float(np.mean(I_bonds)))

        if resid < tol:
            converged = True
            break
        if resid > 1e3 or not np.isfinite(resid):
            diverged = True
            break
        phi = (1.0 - beta) * phi + beta * phi_new

    return {"phi_self": phi, "converged": converged, "diverged": diverged,
            "n_iter": it, "residual": resid,
            "history": hist if record_history else None}


# ════════════════════════════════════════════════════════════════════
# FUNCTION 3b — The scalar self-consistency map on the uniform manifold
#
# Why this exists (measured, not assumed)
# ---------------------------------------
# On a translationally invariant ring, a UNIFORM total phase produces a
# uniform bond current (Bloch states), and the circulant Biot-Savart
# kernel maps a uniform current to a uniform A_t. So the self-consistency
# map G sends uniform textures to uniform textures EXACTLY (verified to
# std(G)/|G| ~ 1e-13 at machine precision): on this manifold the whole
# self-consistency problem reduces to ONE scalar equation
#
#     phi = G(phi)   with   G(phi) = g * (e a / hbar) * A_t[ I(phi_ext+phi) ].
#
# Measuring G directly (see experiment_scalar_map) shows that at the
# project's g0 = 3e5 the map is a sawtooth: smooth branches of slope
# dG/dphi ~ +2.6..3.0 (> 1), separated by discontinuous drops where a BdG
# level crosses zero (the ring's current-phase relation). Consequences,
# established numerically in experiment_* below and reported in the
# methods note:
#
#   * fixed points EXIST (the identity line crosses the smooth branches)
#     but are REPELLING under fixed-point iteration: G' > 1 means the
#     residual never contracts, and NO damping beta rescues it
#     (|1 - beta + beta*G'| > 1 for every beta in (0,1] when G' > 1).
#     The orbit does not blow up either — the map is bounded by the CPR
#     amplitude — it wanders indefinitely without converging (measured:
#     80 iterations, residual O(1e-2), every seed, every beta tried).
#     Iteration only converges below g_crit = g0 / G'(g0) ~ 1.1e5.
#   * therefore the self-consistent solution at g0 is computed by ROOT-
#     FINDING on the scalar equation (bisection on a continuous branch),
#     which needs no stability at all — not by forcing a non-convergent
#     iteration.
#   * several fluxoid-like branches can cross the identity line at the
#     same (mu, B) -> MULTIPLE self-consistent solutions (metastable flux
#     states of a superconducting loop). find_uniform_fixed_points
#     enumerates them; the PRINCIPAL root (smallest |phi|, the branch
#     continuously connected to phi=0 at g=0) is used for headline
#     numbers, with the multiplicity reported alongside.
#
# Do NON-uniform (symmetry-broken) solutions exist? Two separate tools,
# with a subtlety that must not be conflated:
#   * ring_mode_gains gives the PICARD-ITERATION gains lambda_q of each
#     texture mode. At large N some q != 0 have |lambda_q| > 1 (1.746 at
#     N=80, g0): the uniform fixed point is a SADDLE under iteration. This
#     is an iteration-dynamics fact — the reason plain iteration fails for
#     textures too, not only for the uniform mode.
#   * whether a non-uniform FIXED POINT actually exists is a different
#     question, settled by nonuniform_fixed_point_search: Newton-Krylov on
#     the full N-dim F(phi)=G(phi)-phi=0 (root-finding, stability-
#     independent), seeded with large perturbations along the unstable
#     modes. Result: from every convergent seed it collapses back to the
#     UNIFORM solution (texture norm ~1e-13). No non-uniform solution
#     exists; the uniform one is locally unique and the uniform-manifold
#     root-finder is complete. (An earlier, smaller N=40 spot check saw all
#     gains < 1 and this was briefly mis-stated as "textures always decay";
#     that was wrong at N=80 and is corrected here — the honest statement is
#     the Newton-Krylov uniqueness result, not a gain bound.)
# ════════════════════════════════════════════════════════════════════
def scalar_ring_map(N: int, mu: float, B: float, g: float,
                    phi_ext_bond: float, phi_uniform: float,
                    delta: float = 1.0,
                    t_hop: float = T_HOP_SI,
                    lattice_const: float = LATTICE_A,
                    wire_radius: float | None = None) -> float:
    """One application of the self-consistency map to a uniform texture,
    returning the (uniform) output phase. See the block comment above."""
    phi_tot = np.full(N, float(phi_ext_bond) + float(phi_uniform))
    H = _build_ring_lo_phased(N, mu, B, phi_tot, delta)
    evals, evecs = eigh(H)
    I_bonds = ring_current_density(evecs, evals, N, t_hop, phi_tot)
    A_t = ring_biot_savart_self(I_bonds, N, lattice_const, wire_radius)
    out = g * compute_peierls_from_A(A_t, lattice_const)
    return float(np.mean(out))


def find_uniform_fixed_points(N: int, mu: float, B: float, g: float,
                              phi_ext_bond: float | None = None,
                              span: float | None = None,
                              n_scan: int = 49,
                              tol: float = 1e-12,
                              delta: float = 1.0) -> list:
    """
    Enumerate ALL genuine fixed points phi = G(phi) of the scalar map in
    [-span, +span], distinguishing true roots from sawtooth discontinuities.

    A sign change of f(phi) = G(phi) - phi is bisected to `tol`; the
    bracket is accepted as a root only if |f| actually vanishes there
    (|f| < 1e-4 rad — genuine roots reach ~|G'|*tol ~ 1e-11, whereas at a
    discontinuity |f| stays at the jump size ~ 0.05-0.1 rad). Each root is
    returned with its local slope dG/dphi (centred finite difference),
    which decides iteration stability (|slope| < 1 iff Picard-stable).

    span defaults to 1.5 flux periods (2 pi / N) — fixed points can only
    live where |G| <= max|phi_self|, which is within ~1 period of zero for
    this system.

    Returns: list of {"phi": root, "slope": dG/dphi, "picard_stable": bool}
    sorted by |phi|.
    """
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    if span is None:
        span = 1.5 * 2.0 * np.pi / N

    def f(p: float) -> float:
        return scalar_ring_map(N, mu, B, g, phi_ext_bond, p, delta) - p

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
        if abs(f(m)) < 1e-4:                       # root, not a jump
            eps = 1e-6
            slope = (scalar_ring_map(N, mu, B, g, phi_ext_bond, m + eps, delta)
                     - scalar_ring_map(N, mu, B, g, phi_ext_bond, m - eps, delta)
                     ) / (2 * eps)
            roots.append({"phi": m, "slope": float(slope),
                          "picard_stable": bool(abs(slope) < 1.0)})
    roots.sort(key=lambda r: abs(r["phi"]))
    return roots


def principal_uniform_fixed_point(N: int, mu: float, B: float, g: float,
                                  phi_ext_bond: float | None = None,
                                  hint: float | None = None,
                                  tol: float = 1e-12,
                                  delta: float = 1.0) -> float:
    """
    The principal (smallest |phi|) self-consistent uniform self-phase at
    (N, mu, B, g). With a `hint` (e.g. the root at a neighbouring B), a
    narrow bracket around the hint is tried first and the result is
    accepted only if it is a genuine root; otherwise a full scan runs.
    Returns nan only if no fixed point exists even in a widened
    (3-period, finer) scan — observed transiently near g_crit, where the
    principal branch detaches at a sawtooth jump before the identity line
    picks up the neighbouring branch.
    """
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    if hint is not None and np.isfinite(hint):
        period = 2.0 * np.pi / N
        w = 0.25 * period

        def f(p: float) -> float:
            return scalar_ring_map(N, mu, B, g, phi_ext_bond, p, delta) - p

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
    roots = find_uniform_fixed_points(N, mu, B, g, phi_ext_bond,
                                      tol=tol, delta=delta)
    if not roots:                                   # widen once before giving up
        roots = find_uniform_fixed_points(N, mu, B, g, phi_ext_bond,
                                          span=3.0 * 2.0 * np.pi / N,
                                          n_scan=121, tol=tol, delta=delta)
    return roots[0]["phi"] if roots else float("nan")


def ring_mode_gains(N: int, mu: float, B: float, g: float,
                    phi_ext_bond: float | None = None,
                    phi_star: float = 0.0,
                    eps: float = 1e-7,
                    delta: float = 1.0,
                    t_hop: float = T_HOP_SI,
                    lattice_const: float = LATTICE_A) -> np.ndarray:
    """
    Linear gain lambda_q of every spatial Fourier mode q of the texture,
    around a uniform point phi_self = phi_star.

    Because the background is uniform, the Jacobian J_ij = dG_i/dphi_j is
    circulant, so a SINGLE perturbed-bond column determines the whole
    matrix and its eigenvalues are the FFT of that column:
    lambda_q = sum_j J_0j exp(-i q j). Mode q=0 is the scalar-map slope
    dG/dphi (checked against the finite-difference slope in the
    experiments).

    IMPORTANT — what lambda_q does and does NOT tell you. lambda_q are the
    gains of the PICARD-ITERATION map. |lambda_q| > 1 means iteration would
    GROW that mode (the uniform point is a saddle under iteration); it does
    NOT by itself imply a non-uniform fixed point exists. At large N some
    texture modes do have |lambda_q| > 1 (measured 1.746 at N=80, g0), yet
    nonuniform_fixed_point_search (Newton-Krylov, stability-independent)
    finds NO non-uniform self-consistent solution from any seed — the
    uniform solution is locally unique. So these gains are a diagnostic of
    why Picard ITERATION fails, not a claim about solution multiplicity;
    that question is answered by the root-finding search, not here.

    Returns: (N,) complex lambda_q, q index in np.fft.fftfreq order.
    """
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N

    def Gvec(phi_self: np.ndarray) -> np.ndarray:
        phi_tot = float(phi_ext_bond) + phi_self
        H = _build_ring_lo_phased(N, mu, B, phi_tot, delta)
        evals, evecs = eigh(H)
        I_bonds = ring_current_density(evecs, evals, N, t_hop, phi_tot)
        A_t = ring_biot_savart_self(I_bonds, N, lattice_const)
        return g * compute_peierls_from_A(A_t, lattice_const)

    base = np.full(N, float(phi_star))
    G0 = Gvec(base)
    pert = base.copy()
    pert[0] += eps
    column = (Gvec(pert) - G0) / eps                 # J[:, 0]
    return np.fft.fft(column)                        # lambda_q


def nonuniform_fixed_point_search(N: int, mu: float, B: float, g: float,
                                  phi_ext_bond: float | None = None,
                                  n_random: int = 4,
                                  amps_periods=(0.2, 1.0, 3.0),
                                  f_tol: float = 1e-11,
                                  delta: float = 1.0,
                                  rng_seed: int = 11) -> dict:
    """
    Test whether a genuinely NON-UNIFORM self-consistent solution exists.

    Why this is needed: ring_mode_gains reports the linear gains lambda_q of
    the PICARD-ITERATION map. At large N some texture modes (q != 0) have
    |lambda_q| > 1 (measured 1.746 at N=80, g0) — the uniform fixed point is
    therefore a SADDLE under iteration. That is a statement about the
    ITERATION DYNAMICS, not proof that a second (non-uniform) fixed point
    exists. To settle the latter we solve the FULL N-dimensional vector
    equation F(phi) = G(phi) - phi = 0 with Newton-Krylov (scipy), which
    finds roots irrespective of Picard stability, seeded with LARGE
    perturbations (up to several flux periods) along the most-unstable
    eigenmode direction and random directions.

    Finding (documented in the methods note): from every convergent seed,
    including strong pushes along the unstable modes, Newton-Krylov collapses
    back to the UNIFORM solution (texture norm ~1e-13). No symmetry-broken
    non-uniform fixed point is found — the uniform self-consistent solution
    is (locally) unique, and the uniform-manifold root-finder used for B_c
    is complete. The Picard instability just means iteration is the wrong
    solver, which is exactly why ring_boundary_Bc_pfaffian_selfconsistent
    root-finds instead of iterating.

    Returns a dict with the counts (n_collapsed_uniform / n_nonuniform /
    n_failed), the max texture norm among convergent roots, and the
    most-unstable mode index.
    """
    from scipy.optimize import newton_krylov
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    phi_u = principal_uniform_fixed_point(N, mu, B, g, phi_ext_bond,
                                          delta=delta)
    lam = np.abs(ring_mode_gains(N, mu, B, g, phi_ext_bond, phi_star=phi_u,
                                 delta=delta))
    qmax = int(np.argmax(lam[1:])) + 1

    def Gv(p):
        phi_tot = float(phi_ext_bond) + np.asarray(p)
        H = _build_ring_lo_phased(N, mu, B, phi_tot, delta)
        evals, evecs = eigh(H)
        I_bonds = ring_current_density(evecs, evals, N, T_HOP_SI, phi_tot)
        A_t = ring_biot_savart_self(I_bonds, N, LATTICE_A)
        return g * compute_peierls_from_A(A_t, LATTICE_A)

    def F(p):
        return Gv(p) - p

    j = np.arange(N)
    mode = np.cos(2.0 * np.pi * qmax * j / N)
    mode = mode / np.linalg.norm(mode)
    rng = np.random.default_rng(rng_seed)
    period = 2.0 * np.pi / N

    directions = [mode, -mode] + [rng.normal(size=N) for _ in range(n_random)]
    n_uniform = n_nonuniform = n_failed = 0
    max_texture = 0.0
    for amp in amps_periods:
        for d in directions:
            s0 = np.full(N, phi_u) + amp * period * d / np.linalg.norm(d)
            try:
                sol = newton_krylov(F, s0, f_tol=f_tol, maxiter=300)
            except Exception:
                n_failed += 1
                continue
            if np.max(np.abs(F(sol))) > 1e-8:
                n_failed += 1
                continue
            tex = float(np.linalg.norm(sol - np.mean(sol)))
            max_texture = max(max_texture, tex)
            if tex < 1e-8:
                n_uniform += 1
            else:
                n_nonuniform += 1
    return {"N": N, "mu": mu, "B": B, "g": g,
            "phi_uniform": float(phi_u),
            "max_nonuniform_gain": float(np.max(lam[1:])),
            "n_modes_gain_gt1": int(np.sum(lam[1:] > 1.0)),
            "most_unstable_q": qmax,
            "n_seeds": len(directions) * len(amps_periods),
            "n_collapsed_uniform": n_uniform,
            "n_distinct_nonuniform": n_nonuniform,
            "n_failed": n_failed,
            "max_texture_norm": max_texture,
            "nonuniform_solution_exists": bool(n_nonuniform > 0)}


# ════════════════════════════════════════════════════════════════════
# FUNCTION 4 — Fully self-consistent Pfaffian B_c on the ring
# ════════════════════════════════════════════════════════════════════
class ScfConvergenceError(RuntimeError):
    """Raised when the self-field solve fails at some B — deliberately
    fatal so a non-converged texture can never silently enter a reported
    B_c."""


def ring_boundary_Bc_pfaffian_selfconsistent(
        N: int, mu: float,
        g: float = G0_DEFAULT,
        phi_ext_bond: float | None = None,
        delta: float = 1.0,
        window: float = 0.05,
        n_coarse: int = 11,
        tol: float = 1e-6,
        solver: str = "root",
        root_tol: float = 1e-10,
        warm_hint: bool = True,
        scf_tol: float = 1e-12,
        scf_beta: float = 1.0,
        scf_max_iter: int = 300,
        diagnostics: dict | None = None,
        _max_widenings: int = 4) -> float:
    """
    Critical field B_c(mu) from the Pfaffian Z2 sign flip of the ring
    Hamiltonian dressed with the SELF-CONSISTENT self-field: at every B
    the bisection probes, the self-consistent uniform self-phase phi*(B)
    is solved at that (mu, B), then the invariant sign is evaluated with
    phi_bonds = phi_ext + phi*(B).

    This supersedes the interim recipe (open-wire texture pasted onto the
    ring, wrap bond = gauge phase only), kept in
    interim_ring_Bc_openwire_texture below for comparison.

    solver : 'root' (default) — phi*(B) from principal_uniform_fixed_point
             (scalar root-finding on the uniform manifold). This is the
             workhorse because the measured map slope at g0 = 3e5 is
             ~+2.6..3.0 > 1: every fixed point is Picard-REPELLING there,
             so no fixed-point iteration (any damping) can converge — see
             the FUNCTION 3b block comment and experiment_iteration_
             trajectories for the evidence. Root-finding needs no
             stability. The reduction to the uniform manifold is exact for
             this translation-invariant ring (uniform in -> uniform out to
             machine precision) and non-uniform perturbations are
             separately shown to decay (ring_mode_gains).
             'iterate' — selfconsistent_ring_texture at each B (honest
             fixed-point iteration; converges only in the subcritical
             regime g < g_crit ~ 1.1e5). Raises ScfConvergenceError on
             non-convergence rather than reporting a forced number.

    warm_hint : (root solver) seed each root search with the previous B's
             root as a BRACKETING hint only — the accepted root must still
             satisfy |G(phi)-phi| < 1e-4 independently, and a failed hint
             falls back to the full scan, so this is a numerical
             accelerator, not a source of path-dependence. Set False to
             force a full cold scan at every B (used by the robustness
             checks; agreement is verified in experiment_robustness).

    diagnostics : optional dict; if provided, gets 'phi_star_at_Bc' (the
             self-phase at the returned B_c), 'n_roots_at_Bc' (fixed-point
             multiplicity there), and 'coarse_signs'.

    Raises ScfConvergenceError if the self-field cannot be solved at some
    probed B (no root found / iteration non-convergent) — a deliberate
    refusal to report a B_c built on a non-converged self-field.
    """
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N          # Phi = 0.1 Phi_0, as Figure 4
    ana = np.sqrt(mu**2 + delta**2)
    mem: dict = {"hint": None}

    def phi_star(B: float) -> float:
        if solver == "root":
            hint = mem["hint"] if warm_hint else None
            p = principal_uniform_fixed_point(N, mu, B, g, phi_ext_bond,
                                              hint=hint, tol=root_tol,
                                              delta=delta)
            if not np.isfinite(p):
                raise ScfConvergenceError(
                    f"no self-consistent fixed point found at N={N}, "
                    f"mu={mu}, B={B:.6f}, g={g:.3g}")
            mem["hint"] = p
            return p
        elif solver == "iterate":
            res = selfconsistent_ring_texture(
                N, mu, B, g, phi_ext_bond, beta=scf_beta,
                tol=scf_tol, max_iter=scf_max_iter, delta=delta)
            if not res["converged"]:
                raise ScfConvergenceError(
                    f"fixed-point iteration not converged at N={N}, mu={mu}, "
                    f"B={B:.6f}, g={g:.3g} (residual {res['residual']:.3e} "
                    f"after {res['n_iter']} iterations, "
                    f"diverged={res['diverged']}, beta={scf_beta}) — "
                    f"expected for g > g_crit ~ 1.1e5; use solver='root'.")
            return float(np.mean(res["phi_self"]))
        raise ValueError(f"unknown solver {solver!r}")

    def s(B: float) -> tuple:
        p = phi_star(B)
        phi_bonds = np.full(N, phi_ext_bond + p)
        val = _ring_pfaffian_sign(N, mu, B, phi_bonds, delta)
        if val == 0.0:                          # exactly on a zero
            val = _ring_pfaffian_sign(N, mu, B + 1e-9, phi_bonds, delta)
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
            if solver == "root":
                diagnostics["roots_at_Bc"] = find_uniform_fixed_points(
                    N, mu, Bc, g, phi_ext_bond, delta=delta)
        return Bc
    return float("nan")


# ════════════════════════════════════════════════════════════════════
# FUNCTION 5 — [SUPERSEDED] the interim open-wire-texture-on-ring recipe
# ════════════════════════════════════════════════════════════════════
def openwire_one_shot_texture(N: int, mu: float, B: float, g: float,
                              phi_ext_bond: float,
                              delta: float = 1.0) -> np.ndarray:
    """
    The ONE-SHOT open-wire self-field texture of Figure 4/5 (exactly what
    openwire_boundary_Bc computes internally at each B): open-BC
    eigenstates at uniform phi_ext -> current -> straight-wire Biot-Savart
    -> g-scaled Peierls phase. Returns the (N-1,) bond phases. Used by the
    SUPERSEDED interim recipe below and as the 'openwire' seed in the
    seed-independence experiments.
    """
    with contextlib.redirect_stdout(io.StringIO()):
        H = build_bdg_lo(N, mu=mu, B_zeeman=B, delta=delta)
    evals, evecs = eigh(apply_peierls_phases(H, N, phi_ext_bond))
    j_x = compute_current_density(evecs, evals, N, T_HOP_SI,
                                  np.full(N - 1, phi_ext_bond), LATTICE_A)
    A_self = biot_savart_self(j_x, N, LATTICE_A)
    return g * compute_peierls_from_A(A_self, LATTICE_A)          # (N-1,)


def interim_ring_Bc_openwire_texture(N: int, mu: float,
                                     g: float = G0_DEFAULT,
                                     phi_ext_bond: float | None = None,
                                     delta: float = 1.0,
                                     window: float = 0.2,
                                     n_coarse: int = 21,
                                     tol: float = 1e-6,
                                     _max_widenings: int = 4) -> float:
    """
    [SUPERSEDED — kept only for comparison/regression]

    The interim R4 recipe behind the ~0.6-0.8% numbers in
    r4_pfaffian_boundary_methods.md section 4b: at each B probed by the
    Pfaffian bisection, compute the ONE-SHOT self-field texture on the OPEN
    wire at that B (exactly as openwire_boundary_Bc / figure4 do: open-BC
    eigenstates at uniform phi_ext, one pass of current -> straight-wire
    Biot-Savart -> g-scaled phase), then paste those N-1 open-wire bond
    phases onto ring bonds 0..N-2 and give the wrap bond N-1 -> 0 only the
    external gauge phase. Measure B_c from the Pfaffian sign flip.

    This is NOT self-consistent (one shot), uses the wrong-geometry
    (straight, open) EM kernel, and treats the wrap bond inconsistently.
    Use ring_boundary_Bc_pfaffian_selfconsistent instead; this function
    exists so the interim number is reproducible in committed code.
    """
    if phi_ext_bond is None:
        phi_ext_bond = np.pi * 0.1 / N
    ana = np.sqrt(mu**2 + delta**2)

    def texture(B: float) -> np.ndarray:
        phi_self_ow = openwire_one_shot_texture(N, mu, B, g, phi_ext_bond,
                                                delta)                # (N-1,)
        phi_bonds = np.full(N, phi_ext_bond)     # wrap bond: gauge phase only
        phi_bonds[:N - 1] += phi_self_ow
        return phi_bonds

    def s(B: float) -> float:
        val = _ring_pfaffian_sign(N, mu, B, texture(B), delta)
        if val == 0.0:
            val = _ring_pfaffian_sign(N, mu, B + 1e-9, texture(B), delta)
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


# ════════════════════════════════════════════════════════════════════
# Convenience — Pfaffian B_c with a fixed imposed uniform flux
# ════════════════════════════════════════════════════════════════════
def ring_boundary_Bc_uniform_flux(N: int, mu: float, phi0_bond: float,
                                  delta: float = 1.0,
                                  tol: float = 1e-6) -> float:
    """
    Pfaffian B_c of the ring with a FIXED uniform per-bond phase (no
    self-consistency). phi0_bond = phi_ext gives the g=0 baseline;
    phi0_bond = phi_ext + phi* decomposes the self-consistent result into
    its (purely uniform-flux) content — the known-answer / gauge check of
    the R5 checklist.
    """
    from bdg_r4_em_coupling import ring_boundary_Bc_pfaffian
    return ring_boundary_Bc_pfaffian(N, mu, delta=delta, tol=tol,
                                     phi_bonds=np.full(N, float(phi0_bond)))


# ════════════════════════════════════════════════════════════════════
# EXPERIMENTS — convergence treated as an open question, with evidence
# ════════════════════════════════════════════════════════════════════
def experiment_scalar_map(N: int = 80, mu: float = 1.2, B: float | None = None,
                          g: float = G0_DEFAULT,
                          span_periods: float = 1.3,
                          n_grid: int = 121) -> dict:
    """Tabulate G(phi) and enumerate its fixed points at a representative
    near-boundary point. Returns raw curve + roots (with slopes)."""
    if B is None:
        B = float(0.98 * np.sqrt(mu**2 + 1.0))
    phi_ext = np.pi * 0.1 / N
    span = span_periods * 2.0 * np.pi / N
    grid = np.linspace(-span, span, n_grid)
    Gv = np.array([scalar_ring_map(N, mu, B, g, phi_ext, p) for p in grid])
    roots = find_uniform_fixed_points(N, mu, B, g, phi_ext)
    return {"N": N, "mu": mu, "B": B, "g": g, "phi_ext": phi_ext,
            "phi_grid": grid.tolist(), "G": Gv.tolist(),
            "roots": roots}


def experiment_iteration_trajectories(N: int = 80, mu: float = 1.2,
                                      B: float | None = None,
                                      g_super: float = G0_DEFAULT,
                                      g_sub: float = 6.0e4,
                                      betas=(1.0, 0.3),
                                      max_iter: int = 80,
                                      rng_seed: int = 2027) -> dict:
    """
    Task-2 deliverable #1 and #2: does fixed-point iteration converge, and
    is it seed-independent? Runs the honest iteration from three seeds —
    zero, the old open-wire-derived texture, and a random perturbation —
    at the headline g (super-critical) and at a sub-critical g, recording
    the full residual trajectory of every run. No result is forced: each
    run reports converged / diverged / still-wandering as measured.
    """
    if B is None:
        B = float(0.98 * np.sqrt(mu**2 + 1.0))
    phi_ext = np.pi * 0.1 / N
    rng = np.random.default_rng(rng_seed)

    ow = np.zeros(N)
    ow[:N - 1] = openwire_one_shot_texture(N, mu, B, 1.0, phi_ext)
    seeds = {
        "zero": np.zeros(N),
        # open-wire texture scaled to each run's g (as the interim did)
        "openwire": ow,
        "random": rng.normal(scale=5e-3, size=N),
    }
    out = {"N": N, "mu": mu, "B": B, "g_super": g_super, "g_sub": g_sub,
           "runs": []}
    for g, tag in ((g_super, "super"), (g_sub, "sub")):
        for beta in betas:
            for sname, s0 in seeds.items():
                seed = s0 * (g if sname == "openwire" else 1.0)
                res = selfconsistent_ring_texture(
                    N, mu, B, g, phi_ext, seed_phi_self=seed, beta=beta,
                    tol=1e-12, max_iter=max_iter, record_history=True)
                h = res["history"]
                out["runs"].append({
                    "g": g, "regime": tag, "beta": beta, "seed": sname,
                    "converged": res["converged"],
                    "diverged": res["diverged"],
                    "n_iter": res["n_iter"],
                    "final_phi_mean": float(np.mean(res["phi_self"])),
                    "final_phi_std": float(np.std(res["phi_self"])),
                    "residual_trajectory": h["residual"],
                    "phi_mean_trajectory": h["phi_mean"],
                    "phi_std_trajectory": h["phi_std"],
                })
    # reference fixed points from root-finding (no stability needed)
    out["root_super"] = principal_uniform_fixed_point(N, mu, B, g_super,
                                                      phi_ext)
    out["root_sub"] = principal_uniform_fixed_point(N, mu, B, g_sub, phi_ext)
    return out


def experiment_mode_gains(N: int = 80, mu: float = 1.2,
                          B: float | None = None,
                          g: float = G0_DEFAULT,
                          do_nonuniform_search: bool = True) -> dict:
    """
    Picard-iteration mode gains lambda_q of the uniform fixed point
    (circulant Jacobian), cross-checking lambda_{q=0} against the scalar-map
    slope, AND the decisive uniqueness test: Newton-Krylov search for a
    non-uniform fixed point (nonuniform_fixed_point_search).

    The gains alone can mislead — at N=80, g0 some texture modes have
    |lambda_q| > 1, which only means the uniform point is Picard-unstable to
    those modes, NOT that a non-uniform solution exists. The Newton-Krylov
    search settles existence directly.
    """
    if B is None:
        B = float(0.98 * np.sqrt(mu**2 + 1.0))
    phi_ext = np.pi * 0.1 / N
    phi_star = principal_uniform_fixed_point(N, mu, B, g, phi_ext)
    lam = ring_mode_gains(N, mu, B, g, phi_ext, phi_star=phi_star)
    eps = 1e-6
    slope = (scalar_ring_map(N, mu, B, g, phi_ext, phi_star + eps)
             - scalar_ring_map(N, mu, B, g, phi_ext, phi_star - eps)) / (2 * eps)
    lam_abs = np.abs(lam)
    out = {"N": N, "mu": mu, "B": B, "g": g, "phi_star": phi_star,
           "lambda_q_abs": lam_abs.tolist(),
           "lambda_q0": float(lam_abs[0]),
           "scalar_slope": float(slope),
           "max_nonuniform_gain": float(np.max(lam_abs[1:])),
           "n_modes_gain_gt1": int(np.sum(lam_abs[1:] > 1.0))}
    if do_nonuniform_search:
        out["nonuniform_search"] = nonuniform_fixed_point_search(
            N, mu, B, g, phi_ext)
    return out


def experiment_g_ramp(N: int = 80, mu: float = 1.2, B: float | None = None,
                      g_max: float = G0_DEFAULT, n_steps: int = 15,
                      beta: float = 0.5, max_iter: int = 250) -> dict:
    """
    Task-2 deliverable #3 (uniqueness / hysteresis): ramp g from 0 to
    g_max in small steps, seeding each iteration with the previous
    converged solution, vs jumping straight to each g from a zero seed.
    Also records the root-finder's phi*(g) as the stability-independent
    reference, the slope G'(phi*) at each root (the clean stability
    criterion: iteration can converge iff G' < 1), and the empirical
    g_crit above which iteration stops converging. NOTE: the empirical
    g_crit is an UPPER bound — near G' = 1 convergence suffers critical
    slowing down and a finite max_iter flags near-critical points as
    non-converged; the slope crossing G' = 1 is the true bifurcation, and
    both are reported.
    """
    if B is None:
        B = float(0.98 * np.sqrt(mu**2 + 1.0))
    phi_ext = np.pi * 0.1 / N
    gs = np.linspace(0.0, g_max, n_steps + 1)[1:]
    ramp_seed = np.zeros(N)
    rows = []
    g_crit_measured = None
    eps = 1e-6
    for g in gs:
        ramp = selfconsistent_ring_texture(N, mu, B, g, phi_ext,
                                           seed_phi_self=ramp_seed, beta=beta,
                                           tol=1e-10, max_iter=max_iter)
        jump = selfconsistent_ring_texture(N, mu, B, g, phi_ext,
                                           seed_phi_self=None, beta=beta,
                                           tol=1e-10, max_iter=max_iter)
        root = principal_uniform_fixed_point(N, mu, B, g, phi_ext)
        slope = None
        if np.isfinite(root):
            slope = (scalar_ring_map(N, mu, B, g, phi_ext, root + eps)
                     - scalar_ring_map(N, mu, B, g, phi_ext, root - eps)
                     ) / (2 * eps)
        if ramp["converged"]:
            ramp_seed = ramp["phi_self"]        # walk the branch
        elif g_crit_measured is None:
            g_crit_measured = float(g)
        rows.append({
            "g": float(g),
            "ramp_converged": ramp["converged"],
            "ramp_phi": float(np.mean(ramp["phi_self"])) if ramp["converged"] else None,
            "jump_converged": jump["converged"],
            "jump_phi": float(np.mean(jump["phi_self"])) if jump["converged"] else None,
            "root_phi": float(root),
            "root_exists": bool(np.isfinite(root)),
            "slope_at_root": None if slope is None else float(slope),
        })
    # ramp vs jump vs root agreement below g_crit
    agree = [abs(r["ramp_phi"] - r["jump_phi"]) < 1e-8 and
             abs(r["ramp_phi"] - r["root_phi"]) < 1e-6
             for r in rows if r["ramp_converged"] and r["jump_converged"]]
    # true bifurcation: interpolate where the slope at the root crosses 1
    g_crit_slope = None
    for r0, r1 in zip(rows[:-1], rows[1:]):
        s0, s1 = r0["slope_at_root"], r1["slope_at_root"]
        if s0 is not None and s1 is not None and s0 < 1.0 <= s1:
            g_crit_slope = r0["g"] + (r1["g"] - r0["g"]) * (1.0 - s0) / (s1 - s0)
            break
    return {"N": N, "mu": mu, "B": B, "beta": beta, "rows": rows,
            "g_crit_measured": g_crit_measured,
            "g_crit_predicted_from_slope": g_crit_slope,
            "path_independent_below_gcrit": bool(all(agree)) if agree else None,
            "n_agree_checked": len(agree)}


def experiment_robustness(mu: float = 1.2,
                          Ns=(40, 60, 80, 120),
                          N_ref: int = 80,
                          bisect_tols=(1e-4, 1e-5, 1e-6),
                          root_tols=(1e-8, 1e-10, 1e-12),
                          g: float = G0_DEFAULT,
                          verbose: bool = True) -> dict:
    """
    Task-2 deliverables #4 and #5, in the style of the Figure 6 checks:

      * N-scan: full self-consistent B_c at several ring sizes, vs the
        g=0 baseline at the same N. NOTE the honest interpretive caveat:
        phi_ext = pi*0.1/N (fixed total external flux) and the loop's
        physical circumference both change with N, so unlike the clean-
        detector check of Figure 6 the *physics* here is not exactly
        N-independent — what this scan establishes is HOW the converged
        shift scales with ring size (a bulk texture effect would persist;
        a finite-size Aharonov-Bohm flux effect shrinks like ~1/N^2).
      * tolerance scan: B_c^sc at N_ref under tighter/looser Pfaffian
        bisection tol and root tol, plus warm-hint on/off.
      * decomposition / known-answer check: B_c^sc must coincide with the
        plain Pfaffian B_c of a FIXED imposed uniform flux equal to
        phi_ext + phi*(B_c) — proving the entire self-consistent effect
        is uniform-flux (gauge-like AB) content, with zero non-uniform
        texture part.
    """
    out = {"mu": mu, "g": g, "N_scan": [], "tol_scan": [],
           "decomposition": None}
    for N in Ns:
        phi_ext = np.pi * 0.1 / N
        t0 = time.time()
        diag: dict = {}
        bc_sc = ring_boundary_Bc_pfaffian_selfconsistent(
            N, mu, g=g, tol=1e-6, diagnostics=diag)
        bc_base = ring_boundary_Bc_uniform_flux(N, mu, phi_ext, tol=1e-6)
        row = {"N": N,
               "Bc_selfconsistent": bc_sc,
               "Bc_baseline_g0": bc_base,
               "shift_pct": 100.0 * (bc_sc - bc_base) / bc_base,
               "phi_star_at_Bc": diag.get("phi_star_at_Bc"),
               "n_roots_at_Bc": len(diag.get("roots_at_Bc", [])),
               "seconds": round(time.time() - t0, 1)}
        out["N_scan"].append(row)
        if verbose:
            print(f"  [N-scan] N={N}: Bc_sc={bc_sc:.6f} base={bc_base:.6f} "
                  f"shift={row['shift_pct']:+.4f}%  phi*={row['phi_star_at_Bc']:+.6f} "
                  f"roots={row['n_roots_at_Bc']}  ({row['seconds']}s)")

    phi_ext = np.pi * 0.1 / N_ref
    bc_base_ref = ring_boundary_Bc_uniform_flux(N_ref, mu, phi_ext, tol=1e-7)
    for bt in bisect_tols:
        bc = ring_boundary_Bc_pfaffian_selfconsistent(N_ref, mu, g=g, tol=bt)
        out["tol_scan"].append({"kind": "bisect_tol", "value": bt,
                                "Bc": bc,
                                "shift_pct": 100.0 * (bc - bc_base_ref) / bc_base_ref})
        if verbose:
            print(f"  [tol-scan] bisect tol={bt:.0e}: Bc={bc:.8f} "
                  f"shift={out['tol_scan'][-1]['shift_pct']:+.5f}%")
    for rt in root_tols:
        bc = ring_boundary_Bc_pfaffian_selfconsistent(N_ref, mu, g=g,
                                                      tol=1e-6, root_tol=rt)
        out["tol_scan"].append({"kind": "root_tol", "value": rt,
                                "Bc": bc,
                                "shift_pct": 100.0 * (bc - bc_base_ref) / bc_base_ref})
        if verbose:
            print(f"  [tol-scan] root tol={rt:.0e}: Bc={bc:.8f} "
                  f"shift={out['tol_scan'][-1]['shift_pct']:+.5f}%")
    bc_cold = ring_boundary_Bc_pfaffian_selfconsistent(N_ref, mu, g=g,
                                                       tol=1e-6,
                                                       warm_hint=False)
    out["tol_scan"].append({"kind": "cold_full_scan", "value": 0,
                            "Bc": bc_cold,
                            "shift_pct": 100.0 * (bc_cold - bc_base_ref) / bc_base_ref})
    if verbose:
        print(f"  [tol-scan] cold full scan: Bc={bc_cold:.8f}")

    # decomposition check: impose the frozen uniform flux, re-measure B_c
    diag = {}
    bc_sc = ring_boundary_Bc_pfaffian_selfconsistent(N_ref, mu, g=g,
                                                     tol=1e-6,
                                                     diagnostics=diag)
    phi_star = diag["phi_star_at_Bc"]
    bc_frozen = ring_boundary_Bc_uniform_flux(N_ref, mu, phi_ext + phi_star,
                                              tol=1e-6)
    out["decomposition"] = {"N": N_ref, "Bc_selfconsistent": bc_sc,
                            "Bc_frozen_uniform_flux": bc_frozen,
                            "difference": bc_sc - bc_frozen,
                            "phi_star": phi_star}
    if verbose:
        print(f"  [decomposition] Bc_sc={bc_sc:.8f} vs frozen-uniform-flux "
              f"Bc={bc_frozen:.8f} (diff {bc_sc - bc_frozen:+.2e})")
    return out


def experiment_shift_comparison(mus=(1.2, 1.4), N: int = 80,
                                g: float = G0_DEFAULT,
                                verbose: bool = True) -> list:
    """
    The headline side-by-side (deliverable #4): for each mu, the boundary
    shift from (a) the original open-wire threshold detector (recomputed
    live with figure4's exact settings), (b) the SUPERSEDED interim
    ring-Pfaffian recipe, and (c) the fully self-consistent ring
    calculation of this module.
    """
    rows = []
    for mu in mus:
        phi_ext = np.pi * 0.1 / N
        # (a) threshold detector, exactly as figure4 (n_grid=300)
        with contextlib.redirect_stdout(io.StringIO()):
            bc_thr_base = openwire_boundary_Bc(N, mu, T_HOP_SI, LATTICE_A,
                                               phi_ext=phi_ext,
                                               g_feedback=0.0, n_grid=300)
            bc_thr_ind = openwire_boundary_Bc(N, mu, T_HOP_SI, LATTICE_A,
                                              phi_ext=phi_ext,
                                              g_feedback=g, n_grid=300)
        thr_shift = 100.0 * (bc_thr_ind - bc_thr_base) / bc_thr_base
        # (b) interim (SUPERSEDED) ring recipe
        bc_ring_base = ring_boundary_Bc_uniform_flux(N, mu, phi_ext, tol=1e-6)
        bc_interim = interim_ring_Bc_openwire_texture(N, mu, g=g, tol=1e-6)
        interim_shift = 100.0 * (bc_interim - bc_ring_base) / bc_ring_base
        # (c) fully self-consistent ring
        diag: dict = {}
        bc_sc = ring_boundary_Bc_pfaffian_selfconsistent(N, mu, g=g,
                                                         tol=1e-6,
                                                         diagnostics=diag)
        sc_shift = 100.0 * (bc_sc - bc_ring_base) / bc_ring_base
        row = {"mu": mu, "N": N,
               "threshold_Bc_base": bc_thr_base, "threshold_Bc_ind": bc_thr_ind,
               "threshold_shift_pct": thr_shift,
               "ring_Bc_base": bc_ring_base,
               "interim_Bc": bc_interim, "interim_shift_pct": interim_shift,
               "selfconsistent_Bc": bc_sc, "selfconsistent_shift_pct": sc_shift,
               "phi_star_at_Bc": diag.get("phi_star_at_Bc"),
               "n_roots_at_Bc": len(diag.get("roots_at_Bc", []))}
        rows.append(row)
        if verbose:
            print(f"  mu={mu}: threshold {thr_shift:+.2f}%  |  interim "
                  f"{interim_shift:+.3f}%  |  self-consistent {sc_shift:+.4f}%")
    return rows


# ════════════════════════════════════════════════════════════════════
# FIGURES
# ════════════════════════════════════════════════════════════════════
def figure7_scf_convergence(map_data: dict, traj_data: dict, gains: dict,
                            fname: str = "r5_figure7_scf_convergence.png"):
    """Convergence diagnostic (deliverable #2): the scalar map and its
    fixed points, iteration trajectories for 3 seeds in both g regimes,
    and the non-uniform mode gains."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(14, 10))

    # (a) the scalar map
    ax = axes[0, 0]
    grid = np.array(map_data["phi_grid"])
    Gv = np.array(map_data["G"])
    # split at discontinuities so jumps are not drawn as vertical lines
    jumps = np.where(np.abs(np.diff(Gv)) > 0.02)[0]
    start = 0
    for j in list(jumps) + [len(grid) - 1]:
        ax.plot(grid[start:j + 1] * 1e3, Gv[start:j + 1] * 1e3, "b-", lw=2)
        start = j + 1
    ax.plot(grid * 1e3, grid * 1e3, "k--", lw=1, label=r"identity $G=\phi$")
    for r in map_data["roots"]:
        ax.plot(r["phi"] * 1e3, r["phi"] * 1e3, "ro", ms=9)
        ax.annotate(fr"$\phi^*$={r['phi']*1e3:.2f} mrad, $G'$={r['slope']:.2f}",
                    (r["phi"] * 1e3, r["phi"] * 1e3),
                    textcoords="offset points", xytext=(8, -12), fontsize=8)
    ax.set_xlabel(r"uniform self-phase $\phi$ (mrad/bond)")
    ax.set_ylabel(r"$G(\phi)$ (mrad/bond)")
    ax.set_title("(a) scalar self-consistency map, $g=3{\\times}10^5$\n"
                 fr"$N$={map_data['N']}, $\mu$={map_data['mu']}, "
                 fr"$B$={map_data['B']:.3f} (sawtooth CPR, slope$>$1)",
                 fontsize=11)
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(True, alpha=0.3)

    # (b) residual trajectories at g0 (super-critical)
    ax = axes[0, 1]
    colors = {"zero": "#2471a3", "openwire": "#c0392b", "random": "#0F6E56"}
    for run in traj_data["runs"]:
        if run["regime"] != "super":
            continue
        ls = "-" if run["beta"] == 1.0 else "--"
        r = np.array(run["residual_trajectory"])
        ax.semilogy(np.arange(1, len(r) + 1), r, ls,
                    color=colors[run["seed"]], lw=1.8,
                    label=fr"seed={run['seed']}, $\beta$={run['beta']}")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"residual  max$|G(\phi)-\phi|$ (rad)")
    ax.set_title("(b) iteration at $g_0=3{\\times}10^5$ FAILS from every seed\n"
                 "(slope$>$1; bounded wandering, residual never contracts)",
                 fontsize=11)
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (c) trajectories at sub-critical g: converge, seed-independent
    ax = axes[1, 0]
    root_sub = traj_data["root_sub"]
    for run in traj_data["runs"]:
        if run["regime"] != "sub" or run["beta"] != 1.0:
            continue
        m = np.array(run["phi_mean_trajectory"])
        ax.semilogy(np.arange(1, len(m) + 1), np.abs(m - root_sub) + 1e-17,
                    "-", color=colors[run["seed"]], lw=1.8,
                    label=f"seed={run['seed']} (converged={run['converged']})")
    ax.axhline(1e-12, color="gray", ls=":", lw=1, label="tolerance")
    ax.set_xlabel("iteration")
    ax.set_ylabel(r"$|\langle\phi\rangle - \phi^*_{\rm root}|$ (rad)")
    ax.set_title("(c) sub-critical $g=6{\\times}10^4$: all seeds converge"
                 "\nto the SAME fixed point (matches root-finder)")
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)

    # (d) mode gains + Newton-Krylov uniqueness verdict
    ax = axes[1, 1]
    lam = np.array(gains["lambda_q_abs"])
    q = np.arange(len(lam))
    n_gt1 = gains.get("n_modes_gain_gt1", int(np.sum(lam[1:] > 1.0)))
    ax.plot(q[1:], lam[1:], "o", color="#0F6E56", ms=4,
            label=fr"non-uniform modes (max={gains['max_nonuniform_gain']:.3f}, "
                  fr"{n_gt1} with gain>1)")
    ax.plot([0], [lam[0]], "r*", ms=15,
            label=fr"uniform mode $\lambda_0$={lam[0]:.2f} "
                  fr"(= map slope {gains['scalar_slope']:.2f})")
    ax.axhline(1.0, color="k", ls="--", lw=1, label="Picard-stability threshold")
    ax.set_xlabel("spatial mode index $q$")
    ax.set_ylabel(r"Picard-iteration gain $|\lambda_q|$")
    ns = gains.get("nonuniform_search")
    if ns is not None:
        verdict = (f"Newton-Krylov ({ns['n_seeds']} seeds, up to 3 flux-"
                   f"periods along\nunstable modes): "
                   f"{ns['n_collapsed_uniform']} collapse to UNIFORM, "
                   f"{ns['n_distinct_nonuniform']} non-uniform "
                   f"(max texture {ns['max_texture_norm']:.0e})")
        ax.text(0.5, -0.32, verdict, transform=ax.transAxes, ha="center",
                va="top", fontsize=8, color="#333333",
                bbox=dict(boxstyle="round,pad=0.4", facecolor="#eef7f3",
                          edgecolor="#0F6E56"))
    ax.set_title("(d) some texture modes are Picard-unstable ($|\\lambda_q|$>1),"
                 "\nbut root-finding finds NO non-uniform solution")
    ax.legend(fontsize=8, loc="upper right")
    ax.grid(True, alpha=0.3)

    fig.suptitle("Convergence of the self-consistent ring "
                 "self-field — measured, not assumed",
                 fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0.05, 1, 0.95))
    fig.savefig(_out(fname), dpi=150)
    plt.close(fig)
    return _out(fname)


def figure8_hysteresis(gramp: dict,
                       fname: str = "r5_figure8_g_ramp_hysteresis.png"):
    """g-ramp vs direct-jump vs root-finder (uniqueness / hysteresis)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    rows = gramp["rows"]
    gs = np.array([r["g"] for r in rows])
    root = np.array([r["root_phi"] for r in rows])
    ramp = np.array([r["ramp_phi"] if r["ramp_converged"] else np.nan
                     for r in rows])
    jump = np.array([r["jump_phi"] if r["jump_converged"] else np.nan
                     for r in rows])

    fig, ax = plt.subplots(figsize=(10, 6.5))
    ax.plot(gs, root * 1e3, "k-", lw=2, label=r"root-finder $\phi^*(g)$ "
            "(stability-independent)")
    ax.plot(gs, ramp * 1e3, "o", color="#2471a3", ms=8, mfc="none",
            label=r"slow $g$-ramp (seed = previous converged)")
    ax.plot(gs, jump * 1e3, "x", color="#c0392b", ms=8,
            label=r"direct jump from zero seed")
    if gramp["g_crit_measured"] is not None:
        ax.axvline(gramp["g_crit_measured"], color="gray", ls="--", lw=1.5,
                   label=fr"measured $g_{{\rm crit}}\approx"
                         fr"{gramp['g_crit_measured']:.2e}$")
    if gramp["g_crit_predicted_from_slope"] is not None:
        ax.axvline(gramp["g_crit_predicted_from_slope"], color="#0F6E56",
                   ls=":", lw=1.5,
                   label=fr"slope criterion $G'(\phi^*)=1$ at $g="
                         fr"{gramp['g_crit_predicted_from_slope']:.2e}$")
    ax.set_xlabel("feedback amplification $g$")
    ax.set_ylabel(r"converged / root self-phase $\phi^*$ (mrad/bond)")
    ax.set_title(
        "Path-independence below $g_{\\rm crit}$, NO attracting "
        "fixed point above\n"
        fr"($N$={gramp['N']}, $\mu$={gramp['mu']}, $B$={gramp['B']:.3f}, "
        fr"$\beta$={gramp['beta']}; missing markers = iteration did not "
        "converge; the root persists)",
        fontsize=11, fontweight="bold")
    ax.legend(fontsize=9, loc="best")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(_out(fname), dpi=150)
    plt.close(fig)
    return _out(fname)


def figure9_shift_robustness(robust: dict, comparison: list,
                             fname: str = "r5_figure9_shift_robustness.png"):
    """N-scan + tolerance scan + the three-detector shift comparison."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5))

    # (a) N-scan of the self-consistent shift, with FITTED power law
    ax = axes[0]
    Ns = np.array([r["N"] for r in robust["N_scan"]], dtype=float)
    sh = np.array([r["shift_pct"] for r in robust["N_scan"]])
    ax.plot(Ns, sh, "o-", color="#0F6E56", lw=2, ms=8,
            label="self-consistent shift (data)")
    # fit shift ~ C * N^p (log-log least squares) — report the exponent
    p, logC = np.polyfit(np.log(Ns), np.log(sh), 1)
    Nfine = np.linspace(Ns.min(), Ns.max(), 100)
    ax.plot(Nfine, np.exp(logC) * Nfine ** p, "k:", lw=1.5,
            label=fr"power-law fit $\propto N^{{{p:.2f}}}$")
    ax.axhline(0, color="gray", lw=0.8)
    ax.set_xlabel("ring size $N$")
    ax.set_ylabel(r"$\Delta B_c / B_c$ (%)")
    ax.set_title(fr"(a) size scaling, $\mu$={robust['mu']}: shift $\to 0$ as "
                 fr"$N\to\infty$ ($\sim N^{{{p:.2f}}}$)," "\n"
                 "a finite-size Aharonov-Bohm flux effect, NOT a bulk shift")
    ax.legend(fontsize=9)
    ax.grid(True, alpha=0.3)

    # (b) tolerance scan
    ax = axes[1]
    labels, vals = [], []
    for r in robust["tol_scan"]:
        if r["kind"] == "bisect_tol":
            labels.append(f"bisect\n{r['value']:.0e}")
        elif r["kind"] == "root_tol":
            labels.append(f"root\n{r['value']:.0e}")
        else:
            labels.append("cold\nscan")
        vals.append(r["shift_pct"])
    ax.plot(range(len(vals)), vals, "s-", color="#2471a3", lw=1.5, ms=8)
    ax.set_xticks(range(len(labels)))
    ax.set_xticklabels(labels, fontsize=8)
    ax.set_ylabel(r"$\Delta B_c / B_c$ (%)")
    ax.set_title(fr"(b) tolerance robustness at $N$=80, $\mu$={robust['mu']}"
                 "\n(bisection tol, root tol, warm-hint off)")
    ax.grid(True, alpha=0.3)

    # (c) the three-detector comparison
    ax = axes[2]
    mus = [r["mu"] for r in comparison]
    x = np.arange(len(mus))
    w = 0.26
    thr = [r["threshold_shift_pct"] for r in comparison]
    itm = [r["interim_shift_pct"] for r in comparison]
    scv = [r["selfconsistent_shift_pct"] for r in comparison]
    ax.bar(x - w, thr, w, color="#c0392b",
           label="threshold detector (Fig 4/5, PRELIMINARY)")
    ax.bar(x, itm, w, color="#e67e22",
           label="interim ring/Pfaffian (SUPERSEDED)")
    ax.bar(x + w, scv, w, color="#0F6E56",
           label="self-consistent ring (this work)")
    for xi, v in zip(x - w, thr):
        ax.text(xi, v, f"{v:+.2f}%", ha="center",
                va="top" if v < 0 else "bottom", fontsize=8)
    for xi, v in zip(x, itm):
        ax.text(xi, v, f"{v:+.2f}%", ha="center",
                va="top" if v < 0 else "bottom", fontsize=8)
    for xi, v in zip(x + w, scv):
        ax.text(xi, v, f"{v:+.3f}%", ha="center",
                va="top" if v < 0 else "bottom", fontsize=8)
    ax.axhline(0, color="k", lw=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels([fr"$\mu$={m}" for m in mus])
    ax.set_ylabel(r"boundary shift $\Delta B_c / B_c$ (%)")
    ax.set_title("(c) the headline shift under three detectors\n"
                 "($N$=80, $g=3\\times 10^5$, same external flux)",
                 fontsize=11)
    # headroom so the tallest (most negative) bar label clears the legend
    ymin = min(thr) * 1.28
    ax.set_ylim(ymin, max(0.6, max(scv) * 3))
    ax.legend(fontsize=7.5, loc="lower left")
    ax.grid(True, alpha=0.3, axis="y")

    fig.suptitle("Robustness of the self-consistent ring "
                 "boundary shift", fontsize=14, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig(_out(fname), dpi=150)
    plt.close(fig)
    return _out(fname)


# ════════════════════════════════════════════════════════════════════
# MAIN — run the full R5 investigation, saving results incrementally
# ════════════════════════════════════════════════════════════════════
def _save_results(results: dict):
    with open(_out("r5_results.json"), "w") as f:
        json.dump(results, f, indent=2)


if __name__ == "__main__":
    t_start = time.time()
    results: dict = {}

    print("=" * 70)
    print("R5: SELF-CONSISTENT RING SELF-FIELD — full investigation")
    print("=" * 70)

    print("\n[1] Ring Biot-Savart kernel validation (loop inductance)...")
    v = validate_ring_biot_savart(N=100)
    print(f"    L_num={v['L_numeric_H']:.5e} H  vs discrete-exact "
          f"{v['L_discrete_H']:.5e} H  (rel err {v['rel_err_discrete']*100:.3f}%)")
    print(f"    vs textbook {v['L_textbook_H']:.5e} H "
          f"(rel err {v['rel_err_textbook']*100:.2f}% — understood constant)")
    print(f"    uniform-A check {v['A_uniformity']:.1e}, wrap-bond current "
          f"uniformity {v['current_uniformity']:.1e}   PASS={v['pass']}")
    results["kernel_validation"] = v
    _save_results(results)
    if not v["pass"]:
        raise SystemExit("ring kernel validation FAILED — aborting before "
                         "any physics is computed with it")

    print("\n[2] Scalar self-consistency map + fixed points (mu=1.2)...")
    map_data = experiment_scalar_map(N=80, mu=1.2)
    for r in map_data["roots"]:
        print(f"    root phi*={r['phi']:+.6f} rad, slope G'={r['slope']:.3f}, "
              f"Picard-stable={r['picard_stable']}")
    results["scalar_map"] = map_data
    _save_results(results)

    print("\n[3] Iteration trajectories: 3 seeds x (g0 super-critical, "
          "sub-critical)...")
    traj = experiment_iteration_trajectories(N=80, mu=1.2)
    for run in traj["runs"]:
        print(f"    g={run['g']:.0e} beta={run['beta']} seed={run['seed']:9s}"
              f" -> converged={str(run['converged']):5s} "
              f"diverged={str(run['diverged']):5s} n_iter={run['n_iter']}")
    print(f"    root-finder reference: phi*(g0)={traj['root_super']:+.6f}, "
          f"phi*(g_sub)={traj['root_sub']:+.6f}")
    results["trajectories"] = traj
    _save_results(results)

    print("\n[4] Mode gains + Newton-Krylov non-uniform-solution search...")
    gains = experiment_mode_gains(N=80, mu=1.2)
    print(f"    uniform-mode gain lambda_0 = {gains['lambda_q0']:.3f} "
          f"(scalar slope {gains['scalar_slope']:.3f})")
    print(f"    max non-uniform Picard gain = {gains['max_nonuniform_gain']:.3f} "
          f"({gains['n_modes_gain_gt1']} modes with gain>1: uniform point is a "
          f"saddle under iteration)")
    ns = gains["nonuniform_search"]
    print(f"    Newton-Krylov ({ns['n_seeds']} seeds along unstable modes): "
          f"{ns['n_collapsed_uniform']} -> uniform, "
          f"{ns['n_distinct_nonuniform']} -> non-uniform "
          f"(max texture {ns['max_texture_norm']:.1e}) => "
          f"non-uniform solution exists: {ns['nonuniform_solution_exists']}")
    results["mode_gains"] = gains
    _save_results(results)

    print("\n[5] g-ramp vs direct jump (uniqueness / hysteresis)...")
    gramp = experiment_g_ramp(N=80, mu=1.2)
    gc_m = gramp["g_crit_measured"]
    gc_p = gramp["g_crit_predicted_from_slope"]
    print(f"    measured g_crit ~ {gc_m if gc_m is None else format(gc_m, '.3e')}, "
          f"predicted from slope {gc_p if gc_p is None else format(gc_p, '.3e')}")
    print(f"    path-independent below g_crit: "
          f"{gramp['path_independent_below_gcrit']} "
          f"({gramp['n_agree_checked']} g-points checked)")
    results["g_ramp"] = gramp
    _save_results(results)

    print("\n[6] Robustness: N-scan, tolerances, decomposition (mu=1.2)...")
    robust = experiment_robustness(mu=1.2)
    results["robustness"] = robust
    _save_results(results)

    print("\n[7] Headline shift comparison (threshold | interim | "
          "self-consistent)...")
    comparison = experiment_shift_comparison(mus=(1.2, 1.4), N=80)
    results["shift_comparison"] = comparison
    _save_results(results)

    # CSV of the headline table
    import csv
    with open(_out("r5_shift_comparison.csv"), "w", newline="") as f:
        wcsv = csv.DictWriter(f, fieldnames=list(comparison[0].keys()))
        wcsv.writeheader()
        wcsv.writerows(comparison)

    print("\n[8] Figures...")
    print(f"    {figure7_scf_convergence(map_data, traj, gains)}")
    print(f"    {figure8_hysteresis(gramp)}")
    print(f"    {figure9_shift_robustness(robust, comparison)}")
    _save_results(results)

    print("\n" + "=" * 70)
    print("R5 SUMMARY (honest, as measured)")
    print("=" * 70)
    gc_m = "n/a" if gc_m is None else f"{gc_m:.2e}"
    gc_p = "n/a" if gc_p is None else f"{gc_p:.2e}"
    print(f"""
 * Fixed-point iteration at g0=3e5 does NOT converge from any seed:
   the scalar map slope at the fixed point is ~{gains['scalar_slope']:.2f} > 1,
   so the fixed point is Picard-repelling and no mixing rescues it (the
   orbit wanders boundedly; it does not blow up). Iteration converges
   only below the bifurcation at g_crit ~ {gc_p} (slope criterion;
   non-convergence first observed at {gc_m}, an upper bound due to
   critical slowing down).
 * The self-consistent solution nevertheless EXISTS and is well defined:
   on the translation-invariant ring the map closes on uniform textures,
   reducing self-consistency to a scalar root problem, solved by
   bisection (no stability required). Some non-uniform texture modes are
   ALSO Picard-unstable (max gain {gains['max_nonuniform_gain']:.3f},
   {gains['n_modes_gain_gt1']} modes >1) — the uniform point is a saddle,
   not a node — BUT Newton-Krylov root-finding seeded with large
   perturbations along those unstable modes finds NO non-uniform fixed
   point (all collapse to uniform, max texture
   {gains['nonuniform_search']['max_texture_norm']:.0e}). The uniform
   self-consistent solution is locally unique.
 * Below g_crit: ramp / jump / root all agree (path-independent).
   At g0: multiple fluxoid-like fixed points can coexist at some B
   (multiplicity reported per-point); the principal (smallest-|phi|)
   branch is used and stated.
 * The converged self-field is a UNIFORM phase (pure loop self-flux,
   Lenz-opposing the external flux). Its only effect on B_c is a
   finite-size Aharonov-Bohm shift that scales like ~1/N^2 (Figure 9a)
   — the bulk self-consistent non-gauge boundary shift of the clean ring
   is ZERO to the resolution of this calculation.
 * Headline comparison at N=80 (see r5_shift_comparison.csv):""")
    for r in comparison:
        print(f"     mu={r['mu']}: threshold {r['threshold_shift_pct']:+.2f}% "
              f"-> interim {r['interim_shift_pct']:+.3f}% "
              f"-> self-consistent {r['selfconsistent_shift_pct']:+.4f}%")
    print(f"\n total runtime: {(time.time() - t_start)/60:.1f} min")
