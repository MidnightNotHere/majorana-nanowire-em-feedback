"""
ISEF 2027 — R9 (Track 4): multi-subband loop-gain probe
=======================================================
Does the R1–R9 EM-feedback *null* — "inductive self-feedback does not move the
bulk topological boundary" — survive when the wire carries more than one
transverse channel? Every prior track held the channel count at exactly one
(R6/R8's single 4N Lutchyn–Oreg block). The scoping report
(`src/r9_multisubband_scoping.md`) established that the honest, tractable version
of this question is NOT the full self-consistent B_c sweep (34.7–160.6+ h per
configuration — out of reach on this machine) but the *loop-gain slope* s, R5's
bifurcation-map derivative, which needs ~30 eigh per point rather than ~2000.
The hypothesis (scoping §6.1): if s grows ~∝ M with occupied-channel count M,
then R9-T1's dynamical instability at g ≈ 8.67e4 drags below the production
coupling g = 2.9e4 and the single-channel null is an artifact; if s stays flat
in M, the null is *stronger* than a single-channel model could show. Both
outcomes are informative.

═══════════════════════════════════════════════════════════════════════════
PRE-REGISTRATION LAYER (this file, at the commit that first introduces it)
═══════════════════════════════════════════════════════════════════════════
This module's FIRST commit contains ONLY the protocol specification and its
validation — NO s(μ;M) data. Per the scoping report's pre-registration
discipline (the same discipline applied after R5's g=3e5 headline turned out to
be a dynamically-unreachable root chosen before the reachability check existed),
every knob that could otherwise be tuned *after* seeing s(μ) is locked here,
blind to s:

  * the closed-form van Hove exclusion radius  δ(ε) = t·(2(1−ε)/π)²·(1+margin)
    — derived from R6's OWN published tight-binding dispersion, NOT from any
    self-field/s machinery (§7 item 3b of the scoping report);
  * the transverse subband offsets  Δε_n = ε₁·(n²−1),  ε₁ = t·(π a/W)²
    (measured from the first subband, so M=1 reproduces R6 exactly);
  * the wire geometry (W) and tolerance (ε), chosen TOGETHER on this arithmetic
    alone so that δ ≪ one subband spacing (the headroom tension of §7 item 3b);
  * the comparison statistic: the MEDIAN of s(μ) over the surviving window,
    with the δ-neighbourhood of every subband threshold ε_n excluded.

and the mandatory construction gates (scoping §5), which must pass before any
M>1 physics is trusted:

  * M=1 bit-for-bit gate — the builder at M=1 must reproduce
    bdg_r4_em_coupling._build_ring_lo_phased ENTRYWISE (max|ΔH| = 0.0);
  * decoupled-subband limit — with the inter-subband coupling zeroed, the 4NM
    spectrum must equal the sorted union of M single-subband 4N spectra at
    shifted μ − Δε_n, to max|Δ| = 0.0;
  * PH symmetry at 4NM — Σ_x H* Σ_x = −H to machine zero at every M;
  * the van Hove radius δ verified against the ACTUAL diagonalised spectrum
    SPECIFICALLY AT THRESHOLD (the band edge / van Hove point), not just in the
    smooth bulk — genuine np.linalg.eigvalsh, not the closed-form dispersion
    evaluated as arithmetic.

The s(μ;M) sweep functions (protocols (a) and (b)) are added in a SEPARATE,
LATER commit, so this commit cannot compute a single s value.

STRUCTURAL, not device-faithful. A small-M result is a statement about what
inter-subband coupling does to the loop gain, NOT a prediction for Chang-class
wires (which have 10–26 occupied subbands, scoping §3.3). Adding a subband can
flip the topological invariant's parity (Lim, Serra, López & Aguado, PRB 86,
121103), so small-M behaviour must NOT be assumed to extrapolate to the
realistic many-subband regime. This limitation is restated wherever a result is
reported.

Everything geometric/circuit-level is REUSE, not reimplementation:
  * ring_biot_savart_self, compute_peierls_from_A — the self-field kernel and
    Peierls conversion are longitudinal/geometric and unchanged by transverse
    structure (bdg_r5 / bdg_r4);
  * _phi_ext_bond, G_R7, A0 — the fixed external flux, production coupling and
    lattice constant (bdg_r7);
  * _build_ring_lo_phased — imported ONLY as the M=1 bit-for-bit reference.
Nothing in R1–R8 or R9 Tracks 1–3 is modified.
"""

from __future__ import annotations

import contextlib
import io
import json
import os

import numpy as np

with contextlib.redirect_stdout(io.StringIO()):
    from bdg_eigensolver_r1_r2_LOAdapted import LO_PARAMS            # noqa: E402
    from bdg_r4_em_coupling import (                                 # noqa: E402
        _build_ring_lo_phased, compute_peierls_from_A,
        T_HOP_SI, E_CHARGE, HBAR,
    )
    from bdg_r5_selfconsistent_ring import (                         # noqa: E402
        ring_biot_savart_self, ring_current_density, _RING_KERNEL_CACHE,
    )
    from bdg_r7_inductance_sweep import G_R7, A0                     # noqa: E402

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "ISEF R9 data folder")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════════
# LOCKED PROTOCOL PARAMETERS  (pre-registered — blind to s)
# ════════════════════════════════════════════════════════════════════════
# Dimensionless units: Δ = 1. t and α are R6's own normalised values.
T_NORM     = float(LO_PARAMS["t_norm"])       # ≈ 66.3   (t_hop / Δ, Δ=0.25 meV)
ALPHA_NORM = float(LO_PARAMS["alpha_norm"])   # = 8.0    (α_SO / (Δ·a_lat))
A_LAT_NM   = float(LO_PARAMS["a_lat"])        # = 10.0   (nm)
DELTA_MEV  = float(LO_PARAMS["Delta"])        # = 0.25   (meV per Δ-unit, for reporting)

# Geometry / tolerance, chosen TOGETHER (scoping §7 item 3b headroom tension):
# only W = 50 nm keeps δ/ΔE_sb ≈ 0.28 (comfortable headroom); ε ≤ 0.15.
WIRE_W_NM   = 50.0     # transverse hard-wall width (nm)  — LOCKED
EPS_TOL     = 0.15     # "close to bulk" tolerance in the δ derivation — LOCKED
DELTA_MARGIN = 0.15    # +15% safety margin toward MORE exclusion — LOCKED

# Operating point (inherited from R9-T1 / R6, NOT tuned):
OFFSET      = 0.05     # B sits this far below the active channel's transition
G_PROBE     = G_R7     # production coupling 2.9e4 — the physically relevant g
M_MAX       = 4        # bounded-tractable ceiling (scoping §4.5)


def subband_epsilon1(W_nm: float = WIRE_W_NM, t: float = T_NORM,
                     a_nm: float = A_LAT_NM) -> float:
    """First transverse hard-wall level ε₁ = t·(π a / W)²  (dimensionless, Δ=1).
    Reproduces scoping §3.4: W=50 nm → ε₁ = 6.54 meV (= 26.18 Δ-units)."""
    return float(t * (np.pi * a_nm / W_nm) ** 2)


def subband_offsets(M: int, W_nm: float = WIRE_W_NM, t: float = T_NORM,
                    a_nm: float = A_LAT_NM) -> np.ndarray:
    """Transverse subband bottom energies MEASURED FROM THE FIRST SUBBAND:
    Δε_n = ε₁·(n²−1), n=1..M  (so Δε₁ = 0 and M=1 reproduces R6 exactly).
    Subband n is the newly-occupied channel when μ crosses its threshold Δε_n;
    its 1D van Hove singularity ρ(E)∝1/√(E−Δε_n) sits at μ = Δε_n."""
    e1 = subband_epsilon1(W_nm, t, a_nm)
    n = np.arange(1, M + 1, dtype=float)
    return e1 * (n ** 2 - 1.0)


def van_hove_delta(eps: float = EPS_TOL, t: float = T_NORM,
                   margin: float = DELTA_MARGIN) -> float:
    """Closed-form van Hove exclusion radius (dimensionless, Δ=1):

        δ(ε) = t · (2(1−ε)/π)² · (1 + margin).

    Derived from R6's OWN tight-binding chain dispersion E(k) = −2t·cos(k)
    (scoping §7 item 3b): the discrete near-edge level spacing s_j grows from
    the band edge toward the average bulk spacing 4t/N; requiring
    |s_j/s_bulk − 1| < ε for three consecutive levels and solving in the
    large-N limit gives δ = t·(2(1−ε)/π)², independent of N. The formula runs
    ~10–15% BELOW the exact spectrum, so a +15% margin is added in the
    conservative direction (more exclusion — under-excluding is what would let a
    van Hove artifact leak into the accepted window). NO Hamiltonian is
    diagonalised here; this is the known dispersion evaluated as arithmetic.
    Verified against the genuinely diagonalised spectrum at threshold by
    verify_van_hove_delta_at_threshold()."""
    return float(t * (2.0 * (1.0 - eps) / np.pi) ** 2 * (1.0 + margin))


def subband_spacing(M: int, W_nm: float = WIRE_W_NM, t: float = T_NORM,
                    a_nm: float = A_LAT_NM) -> float:
    """One full subband spacing around subband M's threshold:
    ΔE_sb(M) = Δε_{M+1} − Δε_M = ε₁·(2M+1) — the μ-range over which subband M
    is the newest occupied channel before subband M+1 would open."""
    e1 = subband_epsilon1(W_nm, t, a_nm)
    return float(e1 * (2 * M + 1))


# ════════════════════════════════════════════════════════════════════════
# OPERATING POINTS AND SWEEP WINDOWS  (pre-registered — blind to s)
# ════════════════════════════════════════════════════════════════════════
def protocol_a_operating_point(mu_above_eps1: float = 1.2) -> tuple:
    """Protocol (a) — truncation convergence. Hold μ FIXED at the R6/R9-T1
    single-channel operating point: subband 1 filled to μ − Δε₁ = 1.2 above its
    own bottom (Δε₁ = 0), on the topological side, B = √(1.2²+1) − OFFSET. Only
    M (basis size) is raised; subbands 2..M stay unoccupied (they enter only
    virtually). Returns (mu, B)."""
    mu = float(subband_offsets(1)[0] + mu_above_eps1)      # = 1.2
    B = float(np.sqrt(mu_above_eps1 ** 2 + 1.0) - OFFSET)  # = 1.512
    return mu, B


def protocol_b_window(M: int, n_mu: int = 24, W_nm: float = WIRE_W_NM,
                      t: float = T_NORM) -> np.ndarray:
    """Protocol (b) — genuine added occupancy. μ-window for the newly-occupied
    subband M: it STRADDLES the threshold Δε_M so the van Hove point sits inside
    the raw window (making its exclusion auditable, task step 3), and spans one
    full subband spacing ΔE_sb(M) of genuine occupancy ABOVE the threshold:

        μ ∈ [Δε_M − ΔE_sb(M)/3,  Δε_M + ΔE_sb(M)].

    DIRECTION NOTE (locked blind to s): the scoping report §7 item 3b wrote the
    window as `[ε_M − ΔE_sb, ε_M − δ]` (below threshold), which contradicts its
    OWN §6.1a ("μ must be re-tuned UPWARD ... so that channel M is genuinely
    OCCUPIED"). Genuine occupancy of subband M requires μ > Δε_M (the band
    bottom); below it the channel is empty and carries no current. The window is
    therefore placed ABOVE the threshold. This clarification is made here, before
    any s(μ;M) value exists — it is a pre-registration fix, not a post-hoc tune."""
    dEps = subband_offsets(M, W_nm, t)
    thr = float(dEps[-1])                     # threshold of the newest channel M
    span = subband_spacing(M, W_nm, t)
    return np.linspace(thr - span / 3.0, thr + span, n_mu)


def protocol_b_field(mu: float, M: int, offset: float = OFFSET,
                     W_nm: float = WIRE_W_NM, t: float = T_NORM) -> float:
    """Operating field for protocol (b): the newest channel M is held just below
    its OWN transition on the topological side, B = √((μ − Δε_M)² + 1) − offset,
    the multi-subband analogue of R9-T1's B = √(μ² + Δ²) − offset. Tracks
    subband M across the window (blind to s; defined from geometry alone)."""
    dEps = subband_offsets(M, W_nm, t)
    mu_active = mu - float(dEps[-1])
    return float(np.sqrt(mu_active ** 2 + 1.0) - offset)


# ════════════════════════════════════════════════════════════════════════
# INTER-SUBBAND RASHBA COUPLING  (the term Lim et al. PRB 86,121103 makes
# mandatory — omitting it would be qualitatively wrong, not merely approximate)
# ════════════════════════════════════════════════════════════════════════
def intersubband_ky_matrix(M: int, W_nm: float = WIRE_W_NM,
                           alpha: float = ALPHA_NORM,
                           a_nm: float = A_LAT_NM) -> np.ndarray:
    """M×M Hermitian coupling matrix  C[n,m] = (α_SO/Δ)·⟨n|k_y|m⟩  (Δ-units),
    with the hard-wall selection rule (scoping §5.4, verified there against
    quadrature, exact for all n,m ≤ 4):

        ⟨n|k_y|m⟩ = −i · 4 n m / (W (n²−m²))   for (n−m) odd,   0 otherwise
        (diagonal and even-difference elements vanish).

    α_SO/Δ = α_norm·a_lat (both in the module's normalisation). The result is
    purely imaginary and antisymmetric in (n,m), hence C is Hermitian. This is
    the *purely inter-subband* piece σ_x·k_y of the Rashba term; the intra-
    subband piece σ_y·k_x is already inside each diagonal block (R6's h_SO)."""
    C = np.zeros((M, M), dtype=complex)
    scale = alpha * a_nm            # = α_SO/Δ  (nm)
    for n in range(1, M + 1):
        for m in range(1, M + 1):
            if (n - m) % 2 != 0:    # odd difference only
                C[n - 1, m - 1] = -1j * scale * 4.0 * n * m / (W_nm * (n ** 2 - m ** 2))
    return C


# ════════════════════════════════════════════════════════════════════════
# THE 4NM MULTI-SUBBAND RING BdG HAMILTONIAN
# basis = subband(M) ⊗ spin(2) ⊗ site(N), then Nambu doubling → 4NM
# ════════════════════════════════════════════════════════════════════════
def build_ring_lo_multisubband(N: int, mu: float, B: float,
                               phi_bonds: np.ndarray | None = None,
                               M: int = 1,
                               W_nm: float = WIRE_W_NM,
                               delta: float = 1.0,
                               t: float | None = None,
                               alpha_so: float | None = None,
                               couple: bool = True) -> np.ndarray:
    """(4NM)×(4NM) multi-subband Lutchyn–Oreg BdG Hamiltonian on a periodic
    ring with per-bond Peierls phases.

    Diagonal subband blocks: R6's _build_ring_lo_phased normal block at the
    shifted chemical potential μ − Δε_n (Δε_n from subband_offsets), so each
    channel sees its own effective filling. Pairing is intra-subband s-wave
    (block-diagonal). Off-diagonal blocks: the inter-subband Rashba coupling
    C[n,m]·(σ_x ⊗ I_N) (site-diagonal — k_y is a transverse, local operator).
    The hole block is the exact −h*, so PH symmetry holds for any M and any
    phi_bonds. `couple=False` zeroes the inter-subband blocks (decoupled-limit
    validation only).

    At M=1 (Δε₁=0, no off-diagonal blocks) this is ENTRYWISE identical to
    _build_ring_lo_phased(N, mu, B, phi_bonds, delta) — the bit-for-bit gate.
    """
    if M < 1:
        raise ValueError("M must be >= 1")
    delta = float(delta)
    t = LO_PARAMS["t_norm"] * delta if t is None else t
    alpha = LO_PARAMS["alpha_norm"] * delta if alpha_so is None else alpha_so

    if phi_bonds is None:
        phase = np.ones(N, dtype=complex)
    else:
        phi_bonds = np.asarray(phi_bonds, dtype=float)
        if phi_bonds.shape != (N,):
            raise ValueError(
                f"phi_bonds has shape {phi_bonds.shape}, expected {(N,)} "
                f"(one phase per ring bond i->(i+1)%N, including the wrap bond)."
            )
        phase = np.exp(1j * phi_bonds)

    I2 = np.eye(2, dtype=complex)
    IN = np.eye(N, dtype=complex)
    sx = np.array([[0, 1], [1, 0]], dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    sz = np.array([[1, 0], [0, -1]], dtype=complex)
    isy = np.array([[0, 1], [-1, 0]], dtype=complex)

    T = np.zeros((N, N), dtype=complex)
    for i in range(N):
        T[i, (i + 1) % N] = phase[i]
    ThpH = T + T.conj().T
    TmH = T - T.conj().T

    dEps = subband_offsets(M, W_nm, t)          # Δε_n, Δε_1 = 0
    two_n = 2 * N
    dim2 = two_n * M                            # normal-state dimension

    # normal-state block matrix (subband ⊗ spin ⊗ site)
    h_full = np.zeros((dim2, dim2), dtype=complex)
    h_Z = B * np.kron(sz, IN)
    h_SO = (1j * alpha / 2.0) * np.kron(sy, TmH)
    for n in range(M):
        # exact _build_ring_lo_phased grouping: h_kin first, then +h_Z, +h_SO
        h_kin = (2.0 * t - (mu - dEps[n])) * np.kron(I2, IN) - t * np.kron(I2, ThpH)
        hn = h_kin + h_Z + h_SO
        h_full[n * two_n:(n + 1) * two_n, n * two_n:(n + 1) * two_n] = hn

    if couple and M > 1:
        C = intersubband_ky_matrix(M, W_nm, alpha)
        sx_IN = np.kron(sx, IN)
        for n in range(M):
            for m in range(M):
                if n != m and C[n, m] != 0.0:
                    h_full[n * two_n:(n + 1) * two_n,
                           m * two_n:(m + 1) * two_n] = C[n, m] * sx_IN

    # intra-subband s-wave pairing, block-diagonal in subband
    Dblk = delta * np.kron(isy, IN)
    D_full = np.zeros((dim2, dim2), dtype=complex)
    for n in range(M):
        D_full[n * two_n:(n + 1) * two_n, n * two_n:(n + 1) * two_n] = Dblk

    H = np.block([[h_full, D_full], [D_full.conj().T, -h_full.conj()]])
    return H


def ring_current_density_ms(eigenvectors: np.ndarray, eigenvalues: np.ndarray,
                            N: int, M: int, t_hop: float,
                            phi_bonds: np.ndarray) -> np.ndarray:
    """Longitudinal bond current (Amperes) on the ring, summed over ALL occupied
    states across ALL subbands. Identical physics/conventions to R5's
    ring_current_density (prefactor 2·e·t_hop/ħ, electron-block-only,
    occupied = E<0, same per-bond Peierls phase — the flux threads every
    subband identically since it is an orbital/longitudinal effect); the only
    change is an extra sum over the transverse subband index in the inner
    product. At M=1 it reduces to ring_current_density entrywise."""
    eigenvalues = np.asarray(eigenvalues)
    occ = eigenvalues < 0.0
    U = eigenvectors[:2 * N * M, occ]                 # electron block only
    if U.shape[1] == 0:
        return np.zeros(N, dtype=float)
    U = U.reshape(M, 2, N, -1)                        # (subband, spin, site, n_occ)
    u_i = U
    u_ip1 = np.roll(U, -1, axis=2)                    # site (i+1) mod N  (wrap)
    inner = np.sum(np.conj(u_i) * u_ip1, axis=(0, 1))  # sum subband & spin -> (N, n_occ)
    phase = np.exp(1j * np.asarray(phi_bonds, dtype=float))[:, None]
    prefactor = 2.0 * E_CHARGE * t_hop / HBAR         # Amperes
    return prefactor * np.sum(np.imag(inner * phase), axis=1)


# ════════════════════════════════════════════════════════════════════════
# PH symmetry at 4NM  (Σ_x H* Σ_x = −H)
# ════════════════════════════════════════════════════════════════════════
def check_ph_symmetry_ms(H: np.ndarray, N: int, M: int) -> float:
    """ε_PH = ‖H + Σ_x H* Σ_x‖ / ‖H‖ for the 4NM Nambu Hamiltonian; Σ_x swaps
    the 2NM particle and 2NM hole halves. Must be machine zero."""
    half = 2 * N * M
    Sx = np.block([[np.zeros((half, half)), np.eye(half)],
                   [np.eye(half), np.zeros((half, half))]])
    return float(np.linalg.norm(H + Sx @ H.conj() @ Sx) / np.linalg.norm(H))


# ════════════════════════════════════════════════════════════════════════
# MANDATORY CONSTRUCTION GATES (scoping §5)
# ════════════════════════════════════════════════════════════════════════
def validate_M1_gate(N: int = 40, n_points: int = 6, seed: int = 11) -> dict:
    """M=1 bit-for-bit gate: build_ring_lo_multisubband(...,M=1) must reproduce
    _build_ring_lo_phased ENTRYWISE (max|ΔH| = 0.0) over random (mu,B,phi_bonds),
    and its current operator must match ring_current_density exactly."""
    rng = np.random.default_rng(seed)
    max_dH = 0.0
    max_dI = 0.0
    for _ in range(n_points):
        mu = float(rng.uniform(-1.5, 2.0))
        B = float(rng.uniform(0.2, 1.8))
        phi = rng.uniform(-0.05, 0.05, N)
        H_ms = build_ring_lo_multisubband(N, mu, B, phi, M=1)
        H_ref = _build_ring_lo_phased(N, mu, B, phi, delta=1.0)
        max_dH = max(max_dH, float(np.max(np.abs(H_ms - H_ref))))
        # current operator equivalence at M=1
        ev, U = np.linalg.eigh(H_ref)
        I_ref = ring_current_density(U, ev, N, T_HOP_SI, phi)
        I_ms = ring_current_density_ms(U, ev, N, 1, T_HOP_SI, phi)
        max_dI = max(max_dI, float(np.max(np.abs(I_ms - I_ref))))
    return {"N": N, "n_points": n_points, "max_abs_dH": max_dH,
            "max_abs_dI_bonds": max_dI,
            "pass": bool(max_dH == 0.0 and max_dI == 0.0)}


def validate_decoupled_limit(N: int = 30, M: int = 4, mu: float = 30.0,
                             B: float = 0.9, seed: int = 3) -> dict:
    """Decoupled-subband limit (scoping §5.3): with the inter-subband coupling
    zeroed, the 4NM Hamiltonian must be — up to a basis permutation — the exact
    direct sum of the M single-subband 4N Hamiltonians at shifted μ − Δε_n, and
    hence its 4NM spectrum the sorted UNION of the M single-subband spectra.

    Two checks are reported:
      * ENTRYWISE (the rigorous gate): reorder the decoupled 4NM matrix from its
        (Nambu ⊗ subband ⊗ spin ⊗ site) layout to (subband ⊗ Nambu ⊗ spin ⊗
        site) and compare to scipy.linalg.block_diag of the per-subband 4N
        blocks — this is exactly 0.0 (the same entries, only rearranged), the
        analogue of the M=1 bit-for-bit gate;
      * SPECTRUM (as the scoping report words it): |Δ| between the sorted 4NM
        eigvalsh and the sorted union — machine precision (~1e-12), NOT exactly
        0.0, because it compares eigenvalues from independent LAPACK reductions
        of different-sized matrices. Reported honestly as such."""
    from scipy.linalg import block_diag
    phi = np.random.default_rng(seed).uniform(-0.05, 0.05, N)
    dEps = subband_offsets(M)
    H = build_ring_lo_multisubband(N, mu, B, phi, M=M, couple=False)
    blocks = [_build_ring_lo_phased(N, mu - dEps[n], B, phi, delta=1.0)
              for n in range(M)]
    # entrywise: full layout (Nambu,subband,spin,site) -> (subband,Nambu,spin,site)
    T = H.reshape(2, M, 2, N, 2, M, 2, N)
    H_as_directsum = T.transpose(1, 0, 2, 3, 5, 4, 6, 7).reshape(4 * N * M, 4 * N * M)
    DS = block_diag(*blocks)
    max_dH = float(np.max(np.abs(H_as_directsum - DS)))
    # spectrum union (machine precision)
    E_full = np.sort(np.linalg.eigvalsh(H))
    E_union = np.sort(np.concatenate([np.linalg.eigvalsh(b) for b in blocks]))
    max_dE = float(np.max(np.abs(E_full - E_union)))
    return {"N": N, "M": M, "max_abs_H_entrywise": max_dH,
            "max_abs_spectrum_diff": max_dE,
            "spectrum_note": "spectrum diff is LAPACK reduction noise (~1e-12); "
                             "entrywise Hamiltonian agreement is the exact gate",
            "pass": bool(max_dH == 0.0)}


def validate_ph_symmetry(N: int = 30, seed: int = 5) -> dict:
    """PH symmetry Σ_x H* Σ_x = −H to machine zero at every M=1..M_MAX, with
    the inter-subband coupling ON (couple=True)."""
    rng = np.random.default_rng(seed)
    out = {}
    worst = 0.0
    for M in range(1, M_MAX + 1):
        mu = float(rng.uniform(0.0, subband_offsets(M)[-1] + 2.0))
        B = float(rng.uniform(0.3, 1.5))
        phi = rng.uniform(-0.05, 0.05, N)
        H = build_ring_lo_multisubband(N, mu, B, phi, M=M)
        eps = check_ph_symmetry_ms(H, N, M)
        out[f"M={M}"] = eps
        worst = max(worst, eps)
    out["worst_eps_PH"] = worst
    out["pass"] = bool(worst < 1e-12)
    return out


def validate_intersubband_selection_rule(M: int = 4) -> dict:
    """The inter-subband coupling matrix must be Hermitian, have a vanishing
    diagonal, and vanish for every even (n−m) — the hard-wall parity rule."""
    C = intersubband_ky_matrix(M)
    herm = float(np.max(np.abs(C - C.conj().T)))
    diag = float(np.max(np.abs(np.diag(C))))
    even_max = 0.0
    for n in range(1, M + 1):
        for m in range(1, M + 1):
            if (n - m) % 2 == 0:
                even_max = max(even_max, float(abs(C[n - 1, m - 1])))
    return {"M": M, "hermiticity_err": herm, "diagonal_max": diag,
            "even_difference_max": even_max,
            "pass": bool(herm < 1e-14 and diag == 0.0 and even_max == 0.0)}


# ════════════════════════════════════════════════════════════════════════
# VAN HOVE δ — VERIFICATION AGAINST THE GENUINELY DIAGONALISED SPECTRUM
# AT THRESHOLD (the band edge / van Hove point), not just the smooth bulk
# ════════════════════════════════════════════════════════════════════════
def _measure_delta_from_spectrum(E_sorted: np.ndarray, t: float,
                                 N: int, eps: float) -> float:
    """Distance IN ENERGY from the band edge at which the discrete near-edge
    level spacing has converged to within ε of the average bulk spacing
    s_bulk = 4t/N, for three consecutive levels. E_sorted ascending; the edge
    is E_sorted[0] (band bottom = van Hove point)."""
    s_bulk = 4.0 * t / N
    sp = np.diff(E_sorted)
    rel = np.abs(sp / s_bulk - 1.0)
    for j in range(len(rel) - 2):
        if rel[j] < eps and rel[j + 1] < eps and rel[j + 2] < eps:
            return float(E_sorted[j + 1] - E_sorted[0])   # E at the converged level, from edge
    return float("nan")


def verify_van_hove_delta_at_threshold(eps: float = EPS_TOL, t: float = T_NORM,
                                       N_list=(80, 160, 320, 640)) -> dict:
    """Verify the closed-form δ(ε) against the ACTUAL diagonalised spectrum
    SPECIFICALLY AT THRESHOLD (the band edge / van Hove point), for a genuine
    np.linalg.eigvalsh of the tight-binding chain matrix — NOT the analytic
    dispersion evaluated as arithmetic.

    The band edge is exactly where the 1D van Hove singularity ρ∝1/√(E−E_edge)
    lives, so measuring the level-spacing crossover FROM the edge is a threshold
    measurement, not a bulk one. Confirms:
      (i) N-independence — δ_measured stable across N=80..640 at fixed ε;
      (ii) the closed form runs ~10–15% below the measured value (hence the
           +15% margin is conservative);
    and cross-checks the genuine eigvalsh against the analytic E_j formula."""
    delta_closed_no_margin = float(t * (2.0 * (1.0 - eps) / np.pi) ** 2)
    rows = []
    for N in N_list:
        # genuine open-chain tight-binding matrix  h = -t (T + T^T),  eigvalsh
        T = np.diag(np.ones(N - 1), 1)               # upper shift, OPEN (no wrap)
        h_chain = -t * (T + T.T)
        E = np.sort(np.linalg.eigvalsh(h_chain))     # <-- genuine diagonalisation
        d_meas = _measure_delta_from_spectrum(E, t, N, eps)
        # cross-check the eigvalsh against the analytic dispersion
        j = np.arange(1, N + 1)
        E_analytic = np.sort(-2.0 * t * np.cos(j * np.pi / (N + 1)))
        eig_vs_formula = float(np.max(np.abs(E - E_analytic)))
        rows.append({
            "N": int(N),
            "delta_measured_at_threshold": d_meas,
            "eigvalsh_vs_analytic_maxdiff": eig_vs_formula,
            "ratio_closed_over_measured": (delta_closed_no_margin / d_meas
                                           if np.isfinite(d_meas) else float("nan")),
        })
    meas = [r["delta_measured_at_threshold"] for r in rows
            if np.isfinite(r["delta_measured_at_threshold"])]
    n_spread = (float(max(meas) - min(meas)) / float(np.mean(meas))
                if len(meas) > 1 else float("nan"))
    ratios = [r["ratio_closed_over_measured"] for r in rows
              if np.isfinite(r["ratio_closed_over_measured"])]
    return {
        "eps": eps, "t": t,
        "delta_closed_no_margin": delta_closed_no_margin,
        "delta_closed_with_margin": van_hove_delta(eps, t),
        "per_N": rows,
        "N_independence_rel_spread": n_spread,
        "closed_over_measured_mean": float(np.mean(ratios)) if ratios else float("nan"),
        "closed_form_is_conservative_low":
            bool(all(rr < 1.0 for rr in ratios)) if ratios else False,
        "note": "delta measured from a genuine eigvalsh of the tight-binding "
                "chain at the band edge (van Hove point); closed form runs "
                "below measured, so the +15% margin over-excludes (safe).",
    }


# ════════════════════════════════════════════════════════════════════════
# THE PRE-REGISTERED COMPARISON STATISTIC  (pure function — no s computed here)
# ════════════════════════════════════════════════════════════════════════
def preregistered_statistic(mu_grid: np.ndarray, s_values: np.ndarray, M: int,
                            W_nm: float = WIRE_W_NM, t: float = T_NORM,
                            eps: float = EPS_TOL,
                            extra_exclusions: list | None = None) -> dict:
    """THE locked comparison statistic. Given a computed s(μ) curve for M
    subbands, return the MEDIAN of s over the surviving window — every μ within
    the van Hove radius δ of ANY subband threshold Δε_n (n≤M), or of any
    supplied lobe-crossing point, is EXCLUDED. Also returns the RAW (non-
    excluded) median and the full mask, so the effect of exclusion is auditable
    (scoping §7 item 3b step 3 / task step 3). This function computes NO s — it
    only reduces an s(μ) curve already computed elsewhere."""
    mu_grid = np.asarray(mu_grid, dtype=float)
    s_values = np.asarray(s_values, dtype=float)
    finite = np.isfinite(s_values)
    delta = van_hove_delta(eps, t)
    thresholds = list(subband_offsets(M, W_nm, t))
    excl_points = list(thresholds) + list(extra_exclusions or [])
    near_vh = np.zeros_like(mu_grid, dtype=bool)
    for c in excl_points:
        near_vh |= np.abs(mu_grid - c) < delta
    surviving = finite & (~near_vh)
    raw = finite
    def _med(mask):
        return float(np.median(s_values[mask])) if np.any(mask) else float("nan")
    return {
        "M": M, "delta_excl_radius": delta, "thresholds_excluded": thresholds,
        "extra_exclusions": list(extra_exclusions or []),
        "median_s_excluded": _med(surviving),
        "median_s_raw": _med(raw),
        "n_surviving": int(np.sum(surviving)), "n_raw_finite": int(np.sum(raw)),
        "n_excluded_by_vh": int(np.sum(finite & near_vh)),
        "surviving_mask": surviving.tolist(),
        "s_min_surviving": (float(np.min(s_values[surviving])) if np.any(surviving) else float("nan")),
        "s_max_surviving": (float(np.max(s_values[surviving])) if np.any(surviving) else float("nan")),
    }


# ════════════════════════════════════════════════════════════════════════
# PRE-REGISTRATION DRIVER — runs the gates + van Hove verification and writes
# the protocol spec. Computes NO s(μ;M) value.
# ════════════════════════════════════════════════════════════════════════
def run_preregistration() -> dict:
    _RING_KERNEL_CACHE.clear()
    e1 = subband_epsilon1()
    spec = {
        "units": "dimensionless, Delta=1; multiply Delta-units by "
                 f"{DELTA_MEV} meV for physical energies",
        "locked_parameters": {
            "t_norm": T_NORM, "alpha_norm": ALPHA_NORM, "a_lat_nm": A_LAT_NM,
            "wire_W_nm": WIRE_W_NM, "eps_tol": EPS_TOL,
            "delta_margin": DELTA_MARGIN, "offset": OFFSET,
            "g_probe": G_PROBE, "M_max": M_MAX,
        },
        "subband_epsilon1": e1,
        "subband_epsilon1_meV": e1 * DELTA_MEV,
        "subband_offsets_dEps": {f"M={M}": list(subband_offsets(M))
                                 for M in range(1, M_MAX + 1)},
        "subband_thresholds_meV": {f"M={M}": [x * DELTA_MEV for x in subband_offsets(M)]
                                   for M in range(1, M_MAX + 1)},
        "subband_spacing_dEsb": {f"M={M}": subband_spacing(M)
                                 for M in range(1, M_MAX + 1)},
        "van_hove_delta": van_hove_delta(),
        "van_hove_delta_meV": van_hove_delta() * DELTA_MEV,
        "delta_over_spacing_headroom": {f"M={M}": van_hove_delta() / subband_spacing(M)
                                        for M in range(1, M_MAX + 1)},
        "protocol_a_operating_point": {"mu": protocol_a_operating_point()[0],
                                       "B": protocol_a_operating_point()[1],
                                       "M_range": list(range(1, M_MAX + 1))},
        "protocol_b_window": {
            f"M={M}": {"mu_min": float(protocol_b_window(M, 24)[0]),
                       "mu_max": float(protocol_b_window(M, 24)[-1]),
                       "threshold_dEps_M": float(subband_offsets(M)[-1]),
                       "spacing_dEsb": subband_spacing(M),
                       "B_at_threshold": protocol_b_field(float(subband_offsets(M)[-1]), M)}
            for M in range(1, M_MAX + 1)},
        "comparison_statistic": "median of s(mu) over the surviving window, "
                                "excluding |mu - dEps_n| < delta for all n<=M "
                                "(and any detected lobe-crossing); raw median "
                                "reported alongside. g = g_probe (production).",
    }
    gates = {
        "M1_bitforbit": validate_M1_gate(),
        "decoupled_limit": validate_decoupled_limit(),
        "ph_symmetry": validate_ph_symmetry(),
        "intersubband_selection_rule": validate_intersubband_selection_rule(),
    }
    vh = verify_van_hove_delta_at_threshold()
    all_pass = (gates["M1_bitforbit"]["pass"] and gates["decoupled_limit"]["pass"]
                and gates["ph_symmetry"]["pass"]
                and gates["intersubband_selection_rule"]["pass"])
    out = {"protocol_spec": spec, "construction_gates": gates,
           "van_hove_verification": vh, "all_gates_pass": bool(all_pass)}
    with open(_out("r9_multisubband_preregistration.json"), "w") as f:
        json.dump(out, f, indent=2, default=str)
    return out


if __name__ == "__main__":
    import sys
    with contextlib.suppress(Exception):
        sys.stdout.reconfigure(encoding="utf-8")   # Greek letters on Windows cp1252
    print("=" * 72)
    print("R9 Track 4 — PRE-REGISTRATION (protocol spec + gates; NO s data)")
    print("=" * 72)
    res = run_preregistration()
    sp = res["protocol_spec"]
    print(f"\nUnits: {sp['units']}")
    print(f"ε₁ (first transverse level)   = {sp['subband_epsilon1']:.3f} Δ "
          f"= {sp['subband_epsilon1_meV']:.3f} meV  (W={WIRE_W_NM:.0f} nm)")
    print(f"van Hove δ (ε={EPS_TOL}, +{int(DELTA_MARGIN*100)}%) = "
          f"{sp['van_hove_delta']:.3f} Δ = {sp['van_hove_delta_meV']:.3f} meV")
    print("\nSubband thresholds Δε_n (Δ-units):")
    for M in range(1, M_MAX + 1):
        print(f"  M={M}: {[f'{x:.1f}' for x in sp['subband_offsets_dEps'][f'M={M}']]}")
    print("\nHeadroom δ / ΔE_sb(M)  (want << 1):")
    for M in range(1, M_MAX + 1):
        print(f"  M={M}: {sp['delta_over_spacing_headroom'][f'M={M}']:.3f}")
    print("\nConstruction gates:")
    for k, v in res["construction_gates"].items():
        print(f"  {k:32s} pass={v['pass']}   {v}")
    vh = res["van_hove_verification"]
    print("\nVan Hove δ verified against genuine eigvalsh AT THRESHOLD:")
    print(f"  δ_closed (no margin) = {vh['delta_closed_no_margin']:.3f} Δ")
    for r in vh["per_N"]:
        print(f"  N={r['N']:4d}: δ_measured={r['delta_measured_at_threshold']:.3f} Δ, "
              f"closed/measured={r['ratio_closed_over_measured']:.3f}, "
              f"eigvalsh_vs_formula={r['eigvalsh_vs_analytic_maxdiff']:.2e}")
    print(f"  N-independence rel spread = {vh['N_independence_rel_spread']:.3e}")
    print(f"  closed/measured mean = {vh['closed_over_measured_mean']:.3f} "
          f"(conservative-low={vh['closed_form_is_conservative_low']})")
    print(f"\nALL GATES PASS: {res['all_gates_pass']}")
    print("Wrote: ISEF R9 data folder/r9_multisubband_preregistration.json")
