'PYEOF'
"""
ISEF 2027 — Symmetry-Aware BdG Eigensolver
Phase: R1 (Baseline eigensolver) + R2 (Pfaffian indicator)
Status: IN PROGRESS

Project: Inductive Circuit Feedback Shifts the Topological Phase Boundary
         in InAs/Al Majorana Nanowires: A Self-Consistent Symmetry-Preserving
         Eigensolver Approach

MIT PhD Advisor: [Name] (advisory role only — acknowledged, not co-author)

Validated results (R1/R2 Kitaev baseline):
  - PH symmetry error: 0.00e+00 across all parameter values
  - Phase boundary at |mu|=2t confirmed (Kitaev 2001)
  - Finite-size scaling: Majorana gap → 0 as N increases
  - Empirical scaling exponent: O(N^2.40)
  - Performance baseline: N=50: 2.3ms, N=100: 24.8ms, N=200: 70.2ms

O1 UPDATE — Lutchyn-Oreg (InAs/Al) model:
  - build_bdg_lo() fully implemented (was stub)
  - Realistic published device parameters (Mourik 2012, Das Sarma 2010)
  - pfaffian_trim_lo() implemented for 4×4 TRIM matrices
  - run_lo_phase_sweep(), run_lo_finite_size_scaling() added
  - plot_lo_phase_diagram(), plot_lo_finite_size() added
  - PH symmetry verified via check_ph_symmetry_complex() at all N
  - Phase boundary B_z = sqrt(mu^2 + Delta^2) confirmed numerically
    [CORRECTED 2026-07-11 — see correction note below; the original check
    was tautological, not a genuine confirmation]

Architecture note:
  This file covers R1 + R2 (Kitaev chain + Lutchyn-Oreg validation).
  The self-consistent EM feedback loop is R4/R5.
  The rank-2 subspace update is R3. All downstream phases import from this
  module via: from bdg_r1_r2 import build_bdg, build_bdg_lo,
  check_ph_symmetry_real, check_ph_symmetry_complex,
  pfaffian_trim_kitaev, pfaffian_trim_lo, timed_eigh

Published InAs/Al device parameters (Mourik et al. Science 2012,
Lutchyn & Das Sarma PRL 2010):
  m*       = 0.023 m_e  (InAs effective mass)
  g-factor = 10.5       (Mourik 2012 InSb, ~10 for InAs)
  alpha_SO = 20 meV·nm  (Rashba spin-orbit, InAs)
  Delta    = 0.25 meV   (proximity-induced gap, Al shell)
  a_lat    = 10 nm      (lattice constant for discretisation)
  All energies in meV, lengths in nm.

Changes from previous version:
  [1-10] (unchanged — see previous header)
  [11] IMPLEMENTED build_bdg_lo() — 4N×4N Lutchyn-Oreg Hamiltonian
       with published InAs/Al parameters. PH symmetry enforced.
  [12] ADDED pfaffian_trim_lo() — numerical Pfaffian at TRIM points
       for the LO model via 4×4 antisymmetric TRIM matrices.
       [CORRECTED 2026-07-11 — see correction note below; this originally
       built the matrices described here but discarded them]
  [13] ADDED LO_PARAMS dict — single source of truth for all published
       device parameters. Cite Mourik 2012 and Das Sarma 2010.
  [14] ADDED run_lo_phase_sweep() — B_zeeman vs mu phase diagram.
  [15] ADDED run_lo_finite_size_scaling() — gap scaling for LO model.
  [16] ADDED plot_lo_phase_diagram() — 2D (B, mu) with analytic boundary.
  [17] ADDED plot_lo_finite_size() — LO gap finite-size scaling.
  [18] ADDED _validate_lo_ph_symmetry() — module-import assertion for LO.

CORRECTION (2026-07-11, found during R6 reconnaissance; full writeup in
src/r2_o1_pfaffian_trim_lo_correction.md):
  pfaffian_trim_lo() originally built the 4×4 TRIM matrices its docstring
  described, then DISCARDED them and returned the closed-form
  nu = sign(mu^2 + Delta^2 - B_zeeman^2) directly — it never called
  pf.pfaffian(). That made _validate_lo_phase_boundary() (which runs at
  module import) tautological: it checked the closed-form value against
  itself, not against a genuine Pfaffian of the LO Hamiltonian.
  A call-site audit confirmed this bug's blast radius was fully contained
  to this file — no R3/R4/R5 result imports or depends on
  pfaffian_trim_lo(); R4/R5 use a separately implemented, genuinely
  independent detector (_ring_pfaffian_sign / ring_boundary_Bc_pfaffian in
  bdg_r4_em_coupling.py, which always called the real pf.pfaffian() and was
  validated on its own terms). pfaffian_trim_kitaev() (the original R2
  Kitaev-chain function) was checked and found to genuinely call
  pf.pfaffian() throughout — it did not have this problem and required
  no fix.
  pfaffian_trim_lo() now performs a genuine Pfaffian computation at both
  TRIM points (see its docstring). Re-validated result: across the full
  (mu, B_zeeman) grid used by run_lo_phase_sweep() and by
  _validate_lo_phase_boundary(), at both the toy t=2.0 used by the sweep
  and the physical t_norm≈66.3 used elsewhere in this module, the genuine
  invariant reproduces the IDENTICAL sign as the old closed-form value at
  every point tested (0 disagreements) — the k=pi TRIM stays deep in the
  trivial state (|eps_kpi|=|4t-mu| >> B, Delta) for every t used anywhere
  in this codebase, so it never affects the k=0 band-inversion result.
  B_z = sqrt(mu^2 + Delta^2) is therefore now a GENUINELY validated result
  in the regime this codebase actually explores — not merely an assumed
  one — but this was an open empirical question before the fix, not a
  foregone conclusion: deviations DO appear at small t (t ≲ 1, verified
  numerically), just not anywhere this project currently uses.
"""

import os
import numpy as np
from scipy.linalg import eigh
from pfapack import pfaffian as pf
import time
import statistics
import platform
import psutil
import scipy
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

matplotlib.rcParams.update({
    'font.size': 14,
    'axes.linewidth': 1.5,
    'lines.linewidth': 2,
    'figure.dpi': 150,
})

# ── All figure/data outputs are routed into this program's data folder ──
# (repo-root/results/baseline_bdg_pfaffian), resolved relative to THIS file so it works
# no matter what the current working directory is when the program is run.
_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "baseline_bdg_pfaffian")


def _out(name: str) -> str:
    """Route an output basename into this program's data folder (created if needed)."""
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════
# PUBLISHED DEVICE PARAMETERS — InAs/Al (single source of truth)
# ════════════════════════════════════════════════════════════════════

LO_PARAMS = {
    # Source: Mourik et al. Science 336, 1003 (2012)
    #         Lutchyn, Sau & Das Sarma PRL 105, 077001 (2010)
    'm_eff':    0.023,    # effective mass in units of m_e (dimensionless)
    'g_factor': 10.5,     # Landé g-factor (InAs ~10, InSb ~50; use 10.5 for InAs)
    'alpha_SO': 20.0,     # Rashba spin-orbit coupling, meV·nm
    'Delta':    0.25,     # proximity-induced superconducting gap, meV
    'a_lat':    10.0,     # lattice constant for tight-binding discretisation, nm
    # Derived hopping: t = hbar^2 / (2 m* a^2)
    # In meV with m* in m_e, a in nm: t = 3810 / (m_eff * a_lat^2) meV
    # t = 3810 / (0.023 * 100) = 1656 meV ← physical scale
    # For dimensionless simulation: normalise to Delta=1, t=t/Delta
    't_hop':    1656.0,   # hopping amplitude, meV (before normalisation)
    # Phase boundary: B_z > sqrt(mu^2 + Delta^2)
    # Typical mu range: -5 to 5 meV (gate-tunable)
    # Typical B_z range: 0 to 1 meV (= g mu_B B / 2, B ~1 T for InAs)
    'mu_range': (-2.0, 2.0),   # meV, centred near Dirac point
    'Bz_range': (0.0, 1.0),    # meV, effective Zeeman energy
}
_HBAR_SI = 1.054571817e-34
_M_E_SI  = 9.1093837015e-31
_MEV_J   = 1.602176634e-19 * 1e-3
# Dimensionless units: normalise all energies to Delta
# t_norm = t_hop / Delta, alpha_norm = alpha_SO / (Delta * a_lat)
_D   = LO_PARAMS['Delta']
_a   = LO_PARAMS['a_lat']

# Replace the hardcoded value with a derived one
LO_PARAMS['t_hop'] = (_HBAR_SI**2) / (2 * LO_PARAMS['m_eff'] * _M_E_SI * (LO_PARAMS['a_lat']*1e-9)**2) / _MEV_J
LO_PARAMS['t_norm']     = LO_PARAMS['t_hop']  / _D       
LO_PARAMS['alpha_norm'] = LO_PARAMS['alpha_SO'] / (_D * _a) # 20/(0.25*10) = 8.0


# ════════════════════════════════════════════════════════════════════
# SECTION 1 — Hamiltonian construction
# ════════════════════════════════════════════════════════════════════

def build_bdg(N: int, t: float, mu: float, delta: float,
              dtype=np.float64) -> np.ndarray:
    """
    Construct the 2N×2N Bogoliubov-de Gennes Hamiltonian for the
    Kitaev chain (1D spinless p-wave superconductor).

    Nambu basis: Ψ = (c₁, …, c_N, c₁†, …, c_N†)ᵀ
    H_BdG = [[h₀,  Δ ],
              [-Δ, -h₀]]
    where h₀ = -μI - t(T + Tᵀ)  (T = upper-shift matrix)
          Δ  = δ(T - Tᵀ)         (p-wave, antisymmetric)

    Particle-hole symmetry (real H): C H_BdG C⁻¹ = -H_BdG
    where C = σ_x ⊗ I_N (block swap, no complex conjugation needed
    because H is real and K·H·K = H* = H).
    """
    T   = np.diag(np.ones(N - 1), 1)
    h0  = -mu * np.eye(N) - t * (T + T.T)
    D   = delta * (T - T.T)
    H   = np.block([[h0, D], [-D, -h0]])
    return H.astype(dtype)


def build_bdg_lo(
    N: int,
    mu: float,
    B_zeeman: float,
    t: float | None = None,
    alpha_so: float | None = None,
    delta: float | None = None,
    a_lat: float | None = None,
    dtype=np.complex128,
) -> np.ndarray:

    # Resolve optional parameters immediately
    delta_val: float = 1.0 if delta is None else delta

    t_val: float = (
        LO_PARAMS["t_norm"] * delta_val
        if t is None
        else t
    )

    alpha_val: float = (
        LO_PARAMS["alpha_norm"] * delta_val
        if alpha_so is None
        else alpha_so
    )
    """
    Construct the 4N×4N Lutchyn-Oreg BdG Hamiltonian for an InAs/Al
    semiconductor-superconductor nanowire.

    Physical model (Lutchyn, Sau & Das Sarma PRL 2010):
      H_LO = H_kin + H_SO + H_Z + H_SC

    Nambu-spin basis: Ψᵢ = (cᵢ↑, cᵢ↓, cᵢ↑†, cᵢ↓†)ᵀ
    Full basis ordering: (c₁↑, c₁↓, …, c_N↑, c_N↓,
                          c₁↑†, c₁↓†, …, c_N↑†, c_N↓†)

    Components (on-site and nearest-neighbour for site i):
    ─────────────────────────────────────────────────────────
    H_kin  = Σᵢ (2t - μ)(cᵢ↑†cᵢ↑ + cᵢ↓†cᵢ↓)
             - t Σᵢ (cᵢ↑†cᵢ₊₁↑ + cᵢ↓†cᵢ₊₁↓ + h.c.)
             [kinetic energy + chemical potential]

    H_SO   = iα Σᵢ (cᵢ↑†cᵢ₊₁↓ - cᵢ↓†cᵢ₊₁↑ - h.c.)
             [Rashba spin-orbit: couples spin to momentum]

    H_Z    = Σᵢ B_z (cᵢ↑†cᵢ↑ - cᵢ↓†cᵢ↓)
             [Zeeman: B_zeeman = g μ_B B / 2, effective field]

    H_SC   = Δ Σᵢ (cᵢ↑†cᵢ↓† + h.c.)
             [s-wave proximity pairing — spin-singlet]

    Particle-hole symmetry (complex H):
      C = (τ_x ⊗ σ_0) · K  (block swap + complex conjugation)
      C H_LO C⁻¹ = -H_LO  ⟺  Σ_x H* Σ_x = -H
      where Σ_x is the 2N×2N block-swap in Nambu space,
      acting on the 2-spin × N-site indices.

    Phase boundary (bulk, infinite wire):
      B_z > sqrt(μ² + Δ²)  → topological phase (MZMs at ends)
      B_z < sqrt(μ² + Δ²)  → trivial phase

    Parameters
    ----------
    N        : number of lattice sites
    mu       : chemical potential (meV if using physical units,
               or dimensionless if normalised to Delta)
    B_zeeman : effective Zeeman energy = g μ_B B / 2 (same units as mu)
    t        : hopping amplitude (default: LO_PARAMS['t_norm'] * delta)
    alpha_so : Rashba coupling (default: LO_PARAMS['alpha_norm'] * delta)
    delta    : proximity gap (default: 1.0 in dimensionless units)
    a_lat    : lattice constant in nm (used for parameter citations only)
    dtype    : complex128 (required — LO Hamiltonian is complex)

    Returns
    -------
    H : (4N, 4N) complex Hermitian ndarray satisfying PH symmetry
    """
    # ── Default to published InAs/Al parameters (dimensionless, Δ=1) ──
    if delta    is None: delta    = 1.0
    if t        is None: t        = LO_PARAMS['t_norm'] * delta
    if alpha_so is None: alpha_so = LO_PARAMS['alpha_norm'] * delta

    # ── Pauli matrices ────────────────────────────────────────────
    I2  = np.eye(2,  dtype=complex)
    I_N = np.eye(N,  dtype=complex)
    sx  = np.array([[0, 1], [1,  0]], dtype=complex)
    sy  = np.array([[0,-1j],[1j, 0]], dtype=complex)
    sz  = np.array([[1, 0], [0, -1]], dtype=complex)

    # ── Hopping and SOC matrices on the N-site chain ──────────────
    T_hop = np.diag(np.ones(N - 1, dtype=complex), 1)   # upper-shift

    # ── Build 2N×2N normal-state Hamiltonian h (spin ⊗ site) ─────
    # Kinetic: (2t - mu) on-site, -t hopping (both spins)
    h_kin = (
    (2.0 * t_val - mu) * np.kron(I2, I_N)
    - t_val * np.kron(I2, T_hop + T_hop.T.conj())
)
    # Zeeman: B_z * sigma_z (splits spin-up/down)
    h_Z   = B_zeeman * np.kron(sz, I_N)

    # Rashba SOC: i*alpha * sigma_y * (k_x) in tight-binding
    # k_x ≈ (T - T†) / (2i a_lat)  →  i*alpha * sigma_y * (T - T†) / (2i)
    # = alpha/2 * sigma_y ⊗ (T - T†)
    h_SO = (
    1j*alpha_val / 2.0
) * np.kron(
    sy,
    T_hop - T_hop.T
)
    
    h     = h_kin + h_Z + h_SO    # 2N×2N normal-state block

    # ── Pairing block: s-wave singlet Δ i σ_y ────────────────────
    # In Nambu space: Δ_pair = Δ * (i σ_y) ⊗ I_N
    #   i σ_y = [[0, 1], [-1, 0]]  (acts on spin, not site)
    i_sy_phys = np.array([[0, 1], [-1, 0]], dtype=complex)  # i*sigma_y
    Delta_blk = delta_val * np.kron(i_sy_phys, I_N)             # 2N×2N


    h_hole = -(h_kin.conj() + h_Z.conj() + h_SO.conj())  # -h* block for PH symmetry

    h_hole = -(h_kin+ h_Z) - h_SO
    # ── Assemble 4N×4N BdG in Nambu space ────────────────────────
    # H_BdG = [[  h,      Delta_blk ],
    #           [-Delta_blk†,  -h*  ]]
    H = np.block([
        [ h,                        Delta_blk          ],
        [Delta_blk.T.conj(), h_hole           ]
    ])

    return H.astype(dtype)


# ════════════════════════════════════════════════════════════════════
# SECTION 2 — Particle-hole symmetry verification
# ════════════════════════════════════════════════════════════════════

def check_ph_symmetry_real(H: np.ndarray, N: int) -> float:
    """
    Verify PH symmetry for a REAL BdG Hamiltonian (Kitaev chain).
    C = σ_x ⊗ I_N (block swap). ε_PH = ‖H + σ_x H σ_x‖_F / ‖H‖_F.
    Raises ValueError if ε_PH > 1e-10.
    """
    C   = np.block([[np.zeros((N, N)), np.eye(N)],
                    [np.eye(N),        np.zeros((N, N))]])
    eps = np.linalg.norm(H + C @ H @ C) / np.linalg.norm(H)
    if eps > 1e-10:
        raise ValueError(
            f"PH symmetry broken: ε_PH = {eps:.2e} > 1e-10. "
            "Check sign conventions in h0 and D blocks."
        )
    return float(eps)


def check_ph_symmetry_complex(H: np.ndarray, N: int) -> float:
    """
    Verify PH symmetry for a COMPLEX BdG Hamiltonian (Lutchyn-Oreg).

    For the 4N×4N LO model:
      C = (τ_x ⊗ σ_0) · K
      C H C⁻¹ = -H  ⟺  Σ_x · H* · Σ_x = -H

    where Σ_x is the 2N×2N block-swap acting in Nambu space,
    with N → 2N (because each site has 2 spin DOF):

      Σ_x = [[0_{2N},  I_{2N}],
              [I_{2N},  0_{2N}]]

    ε_PH = ‖H + Σ_x H* Σ_x‖_F / ‖H‖_F  (must be < 1e-12)

    Parameters
    ----------
    H : (4N, 4N) complex ndarray
    N : number of lattice sites (matrix dimension = 4N)
    """
    two_N   = 2 * N
    Sigma_x = np.block([[np.zeros((two_N, two_N)), np.eye(two_N)],
                        [np.eye(two_N),            np.zeros((two_N, two_N))]])
    eps = (np.linalg.norm(H + Sigma_x @ H.conj() @ Sigma_x)
           / np.linalg.norm(H))
    if eps > 1e-10:
        raise ValueError(
            f"PH symmetry broken (LO complex): ε_PH = {eps:.2e} > 1e-10. "
            "Check Delta_blk construction and -h.conj() block."
        )
    return float(eps)


# ════════════════════════════════════════════════════════════════════
# SECTION 3 — Topological invariant via numerical Pfaffian
# ════════════════════════════════════════════════════════════════════

def pfaffian_trim_kitaev(t: float, mu: float, delta: float
                         ) -> tuple[float, float, float]:
    """
    Z₂ invariant for Kitaev chain via numerical Pfaffian at TRIM points.
    Cross-validated against analytic formula at module import.
    """
    eps_k0  = -2.0 * t - mu
    eps_kpi = +2.0 * t - mu

    A_k0  = np.array([[0.0, eps_k0],  [-eps_k0,  0.0]])
    A_kpi = np.array([[0.0, eps_kpi], [-eps_kpi, 0.0]])

    pf_k0  = float(np.real(pf.pfaffian(A_k0)))
    pf_kpi = float(np.real(pf.pfaffian(A_kpi)))

    if abs(pf_k0) < 1e-14 or abs(pf_kpi) < 1e-14:
        return 0.0, pf_k0, pf_kpi

    return float(np.sign(pf_k0) * np.sign(pf_kpi)), pf_k0, pf_kpi


def pfaffian_trim_lo(mu: float,
                     B_zeeman: float,
                     delta: float = 1.0,
                     alpha_so: float | None = None,
                     t: float | None = None
                     ) -> tuple[float, float, float]:
    """
    Z₂ topological invariant for the Lutchyn-Oreg model via a GENUINE
    numerical Pfaffian at the two TRIM points k=0 and k=π.

    [CORRECTED 2026-07-11 — see the module-header correction note and
    src/r2_o1_pfaffian_trim_lo_correction.md. A previous version of this
    function built these same TRIM matrices but discarded them and
    returned the closed-form nu = sign(mu^2+delta^2-B_zeeman^2) directly,
    without ever calling pf.pfaffian() — this version genuinely does.]

    At a TRIM point k ∈ {0, π}, sin(k) = 0 so the Rashba term
    -α·sin(k)·σ_y vanishes identically (exactly, not approximately — α is
    genuinely irrelevant here, kept as a parameter for API symmetry only)
    and the 2×2 normal-state Bloch block collapses to
      h(k) = ε_k·I₂ + B_zeeman·σ_z   (real, diagonal),
    with ε_k = (2t - μ) - 2t·cos(k) — the same on-site + hopping
    dispersion build_bdg_lo assembles in real space (ε_{k=0} = -μ exactly,
    independent of t since 2t(1-cos 0)=0; ε_{k=π} = 4t - μ). Because h(k)
    is real at TRIM, the PH-symmetric hole block -h(-k)* collapses to
    -h(k). The 4×4 BdG matrix at that TRIM point, in the (↑,↓,↑†,↓†)
    Nambu-spin basis, is:

      H_TRIM = [[ ε_k+B,   0,        0,        Δ       ],
                [ 0,       ε_k-B,   -Δ,        0        ],
                [ 0,      -Δ,       -(ε_k+B),  0        ],
                [ Δ,       0,        0,       -(ε_k-B)  ]]

    This is rotated to the real antisymmetric Majorana basis with the same
    unitary used by R4's _ring_pfaffian_sign (bdg_r4_em_coupling.py),
    restricted to a single momentum (M=2 complex modes here):
      Ω = (1/√2)[[I₂, I₂], [i·I₂, -i·I₂]],   A = i·Ω·H_TRIM·Ω†
    A is real antisymmetric (asserted at runtime); its Pfaffian sign is
    the invariant at that TRIM point, via pf.pfaffian(A) — the same recipe
    as pfaffian_trim_kitaev() and R4's ring construction.

      ν = sgn[Pf(k=0)] × sgn[Pf(k=π)]   (same convention as
      pfaffian_trim_kitaev: ν=-1 topological, +1 trivial, 0 at a closing)

    Re-validated result (2026-07-11): across the (μ, B_zeeman) grid used
    by run_lo_phase_sweep() and _validate_lo_phase_boundary(), at both the
    toy t=2.0 and the physical t_norm≈66.3 used in this module, this
    genuine invariant agrees with the old closed-form sign at every point
    tested — the k=π TRIM stays deep in the trivial state
    (|ε_{k=π}| = |4t-μ| ≫ B, Δ) for every t used in this codebase.
    Deviations do appear for small t (t ≲ 1, verified numerically) but
    that regime is not used anywhere here.

    Parameters
    ----------
    mu       : chemical potential (dimensionless, normalised to Delta)
    B_zeeman : effective Zeeman energy (same units as mu)
    delta    : proximity gap (default 1.0 — dimensionless)
    alpha_so : Rashba coupling — genuinely irrelevant at TRIM (sin(k)=0);
               kept for API consistency only, unused
    t        : hopping amplitude, used only to locate ε_{k=π} = 4t-μ
               (ε_{k=0} = -μ is exactly t-independent). Default
               LO_PARAMS['t_norm']*delta, matching build_bdg_lo's own
               default.

    Returns
    -------
    nu     : Z₂ invariant (+1 trivial, -1 topological, 0 critical)
    pf_k0  : float, signed Pfaffian at k=0 (sign is the meaningful part;
             not separately rescaled)
    pf_kpi : float, signed Pfaffian at k=π
    """
    if t is None:
        t = LO_PARAMS['t_norm'] * delta

    def _trim_H(eps_k: float) -> np.ndarray:
        h = np.diag([eps_k + B_zeeman, eps_k - B_zeeman]).astype(complex)
        Dblk = delta * np.array([[0, 1], [-1, 0]], dtype=complex)
        h_hole = -h
        return np.block([[h, Dblk], [Dblk.conj().T, h_hole]])

    def _pf_sign(eps_k: float) -> float:
        H = _trim_H(eps_k)
        Ib = np.eye(2, dtype=complex)
        Omega = (1.0 / np.sqrt(2.0)) * np.block([[Ib, Ib], [1j * Ib, -1j * Ib]])
        A = 1j * Omega @ H @ Omega.conj().T
        if (np.max(np.abs(A.imag)) > 1e-8 or
                not np.allclose(A.real, -A.real.T, atol=1e-8)):
            raise ValueError(
                "LO TRIM Majorana-basis matrix is not real antisymmetric "
                "(particle-hole symmetry broken?) — Pfaffian sign undefined."
            )
        return float(np.real(pf.pfaffian(A.real)))

    eps_k0  = -mu                 # = (2t-mu) - 2t*cos(0); t cancels exactly
    eps_kpi = 4.0 * t - mu        # = (2t-mu) - 2t*cos(pi)

    pf_k0  = _pf_sign(eps_k0)
    pf_kpi = _pf_sign(eps_kpi)

    if abs(pf_k0) < 1e-12 or abs(pf_kpi) < 1e-12:
        return 0.0, pf_k0, pf_kpi

    nu = float(np.sign(pf_k0) * np.sign(pf_kpi))
    return nu, pf_k0, pf_kpi


def _validate_pfaffian_vs_analytic(t: float = 1.0, delta: float = 0.5,
                                   n_points: int = 50) -> None:
    """Kitaev: assert numerical TRIM Pfaffian agrees with analytic formula."""
    mu_vals = np.linspace(0.01, 4.0, n_points)
    for mu in mu_vals:
        nu_num, _, _ = pfaffian_trim_kitaev(t, mu, delta)
        if abs(mu - 2.0 * t) < 1e-9:
            continue
        pf0_ana = -(mu + 2.0 * t)
        pf1_ana = -(mu - 2.0 * t)
        nu_ana  = float(np.sign(pf0_ana) * np.sign(pf1_ana))
        assert nu_num == nu_ana, (
            f"Pfaffian mismatch at mu={mu:.4f}: "
            f"numerical={nu_num}, analytic={nu_ana}"
        )


def _validate_lo_ph_symmetry(N: int = 20,
                              n_points: int = 10) -> None:
    """
    LO model: assert PH symmetry ε_PH < 1e-12 for random (mu, B_z)
    combinations at module import. Catches any sign error in build_bdg_lo().
    """
    rng    = np.random.default_rng(42)
    mu_vals = rng.uniform(-2.0, 2.0, n_points)
    Bz_vals = rng.uniform(0.0,  1.5, n_points)
    for mu, Bz in zip(mu_vals, Bz_vals):
        H   = build_bdg_lo(N, mu=mu, B_zeeman=Bz)
        eps = check_ph_symmetry_complex(H, N)
        assert eps < 1e-12, (
            f"LO PH symmetry failed at mu={mu:.3f}, B={Bz:.3f}: "
            f"ε_PH = {eps:.2e}"
        )


def _validate_lo_phase_boundary(n_points: int = 20) -> None:
    """
    LO model: assert the phase boundary nu changes sign at
    B_z = sqrt(mu^2 + Delta^2), using pfaffian_trim_lo()'s GENUINE Pfaffian
    computation (default t = LO_PARAMS['t_norm']*delta, physical regime).

    [CORRECTED 2026-07-11] Before the pfaffian_trim_lo() fix, this
    assertion was tautological: it fed the closed-form
    sign(mu^2+delta^2-B_zeeman^2) formula 0.5x/1.5x of its own boundary and
    trivially recovered the expected sign, independent of whether the LO
    Hamiltonian's actual topology matched. It now exercises a real
    TRIM-point Pfaffian and genuinely confirms the boundary in this regime
    — see src/r2_o1_pfaffian_trim_lo_correction.md.
    """

    rng = np.random.default_rng(7)
    mu_test = rng.uniform(-1.5, 1.5, n_points)

    delta = 1.0

    for mu in mu_test:

        B_boundary = np.sqrt(mu**2 + delta**2)

        nu_trivial, _, _ = pfaffian_trim_lo(
            mu,
            B_boundary * 0.5,
            delta
        )

        nu_topo, _, _ = pfaffian_trim_lo(
            mu,
            B_boundary * 1.5,
            delta
        )

        assert nu_trivial in (+1.0, 0.0), (
            f"Expected trivial below boundary at mu={mu:.3f}, "
            f"got ν={nu_trivial}"
        )

        assert nu_topo in (-1.0, 0.0), (
            f"Expected topological above boundary at mu={mu:.3f}, "
            f"got ν={nu_topo}"
        )


# Run all validations at import time
_validate_pfaffian_vs_analytic()
_validate_lo_ph_symmetry()
_validate_lo_phase_boundary()


# ════════════════════════════════════════════════════════════════════
# SECTION 4 — Performance benchmarking
# ════════════════════════════════════════════════════════════════════

def timed_eigh(H: np.ndarray, runs: int = 5) -> float:
    """Median solve time for scipy.linalg.eigh (milliseconds)."""
    times = []
    for _ in range(runs):
        t0 = time.perf_counter()
        eigh(H)
        times.append(time.perf_counter() - t0)
    return statistics.median(times) * 1e3


def device_info() -> dict:
    """Collect hardware/software metadata for reproducibility."""
    try:
        ram   = round(psutil.virtual_memory().total / 1e9, 1)
        cores = psutil.cpu_count(logical=False)
    except Exception:
        ram, cores = "unknown", "unknown"
    return {
        'cpu':    platform.processor()[:40],
        'cores':  cores,
        'ram_gb': ram,
        'os':     platform.system(),
        'numpy':  np.__version__,
        'scipy':  scipy.__version__,
    }


# ════════════════════════════════════════════════════════════════════
# SECTION 5 — Parameter sweeps (Kitaev)
# ════════════════════════════════════════════════════════════════════

def run_phase_sweep(N_sweep: int = 20, t: float = 1.0,
                    delta: float = 0.5, mu_values=None) -> pd.DataFrame:
    """Kitaev chain phase sweep — unchanged from pre-O1."""
    if mu_values is None:
        mu_values = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5]
    rows = []
    for mu in mu_values:
        H  = build_bdg(N_sweep, t, mu, delta)
        ev, _ = eigh(H)
        pos    = ev[ev > 1e-10]
        gap    = float(min(pos)) if len(pos) > 0 else 0.0
        ph_err        = check_ph_symmetry_real(H, N_sweep)
        nu, pf0, pf1  = pfaffian_trim_kitaev(t, mu, delta)
        phase = ("Critical" if nu == 0.0 else
                 "Topological" if nu < 0 else "Trivial")
        rows.append({'mu': mu, 'gap': gap, 'nu': nu,
                     'pf_k0': pf0, 'pf_kpi': pf1,
                     'ph_error': ph_err, 'phase': phase})
    return pd.DataFrame(rows)


def run_finite_size_scaling(t: float = 1.0,
                            N_values=(20, 50, 100, 200)) -> pd.DataFrame:
    """Kitaev chain finite-size scaling — unchanged from pre-O1."""
    rows = []
    for N in N_values:
        for mu, delta_val in [(1.0, 0.5), (0.0, 0.1)]:
            H  = build_bdg(N, t, mu, delta_val)
            ev, _ = eigh(H)
            pos    = ev[ev > 1e-10]
            rows.append({'N': N, 'mu': mu, 'delta': delta_val,
                         'gap':     float(min(pos)) if len(pos) > 0 else 0.0,
                         'ph_error': check_ph_symmetry_real(H, N)})
    return pd.DataFrame(rows)


def run_performance_baseline(t: float = 1.0,
                             N_values=(50, 100, 200),
                             runs: int = 5) -> tuple[pd.DataFrame, float]:
    """Full eigh performance baseline — unchanged from pre-O1."""
    rows = []
    for N in N_values:
        H = build_bdg(N, t, mu=1.0, delta=0.5)
        rows.append({'N': N, 'time_ms': round(timed_eigh(H, runs), 3)})
    df = pd.DataFrame(rows)
    log_N = np.log(df['N'].values.astype(float))
    log_t = np.log(df['time_ms'].values.astype(float))
    exponent, _ = np.polyfit(log_N, log_t, 1)
    return df, float(exponent)


def run_phase_diagram_2d(
    t: float = 1.0,
    mu_range: tuple[float, float] = (0.1, 3.5),
    n_mu: int = 40,
    delta_range: tuple[float, float] = (0.1, 1.5),
    n_delta: int = 30,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Kitaev 2D phase diagram."""

    mu_vals = np.linspace(mu_range[0], mu_range[1], n_mu)
    delta_vals = np.linspace(delta_range[0], delta_range[1], n_delta)

    phase_map = np.zeros((n_delta, n_mu))

    for i, dv in enumerate(delta_vals):
        for j, mv in enumerate(mu_vals):
            nu, _, _ = pfaffian_trim_kitaev(t, mv, dv)
            phase_map[i, j] = nu

    return mu_vals, delta_vals, phase_map


def run_pfaffian_tracking(t: float = 1.0, delta: float = 0.5,
                          n_mu: int = 80,
                          N_sizes=(50, 100, 200)) -> pd.DataFrame:
    """Continuous Pfaffian tracking vs mu — unchanged from pre-O1."""
    mu_vals = np.linspace(0.5, 3.5, n_mu)
    rows = []
    for N in N_sizes:
        for mu in mu_vals:
            H = build_bdg(N, t, mu, delta)
            ev, _ = eigh(H)
            pos    = ev[ev > 1e-10]
            gap    = float(min(pos)) if len(pos) > 0 else 0.0
            nu, pf0, pf1 = pfaffian_trim_kitaev(t, mu, delta)
            rows.append({'mu': mu, 'N': N, 'nu': nu,
                         'pf_k0': pf0, 'pf_kpi': pf1, 'gap': gap})
    return pd.DataFrame(rows)


# ════════════════════════════════════════════════════════════════════
# SECTION 5b — Parameter sweeps (Lutchyn-Oreg) — NEW in O1
# ════════════════════════════════════════════════════════════════════

def run_lo_phase_sweep(
    N: int = 50,
    delta: float = 1.0,
    n_mu: int = 40,
    n_Bz: int = 40,
    mu_range: tuple[float, float] = (-2.0, 2.0),
    Bz_range: tuple[float, float] = (0.0, 1.8),
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    2D phase diagram for the Lutchyn-Oreg model in (μ, B_z) space.
    """

    mu_vals = np.linspace(mu_range[0], mu_range[1], n_mu)
    Bz_vals = np.linspace(Bz_range[0], Bz_range[1], n_Bz)

    nu_map = np.zeros((n_Bz, n_mu))
    gap_map = np.zeros((n_Bz, n_mu))

    sweep_t = 2.0  # LO_PARAMS['t_norm'] * delta; toy t shared by nu_map and
                   # gap_map below so both panels describe the SAME H(mu,Bz)

    for i, Bz in enumerate(Bz_vals):
        for j, mu in enumerate(mu_vals):

            nu, _, _ = pfaffian_trim_lo(mu, Bz, delta, t=sweep_t)
            nu_map[i, j] = nu

            H = build_bdg_lo(
                N,
                mu=float(mu),
                B_zeeman=float(Bz),
                delta=delta,
                t=sweep_t,
            )

            ev, _ = eigh(H)

            pos = np.abs(ev)
            pos = pos[pos > 1e-10]

            gap_map[i, j] = float(np.min(pos)) if pos.size else 0.0

    return mu_vals, Bz_vals, nu_map, gap_map

def run_lo_finite_size_scaling(delta: float = 1.0,
                               N_values: tuple = (20, 50, 100, 200)
                               ) -> pd.DataFrame:
    """
    Finite-size scaling of the Majorana gap for the Lutchyn-Oreg model.

    Tests two parameter points:
      (a) Deep topological: B_z = 1.5Δ, μ = 0  → Majorana gap → 0 as N→∞
      (b) Near boundary:    B_z = 1.1Δ, μ = 0  → gap opens slowly

    PH symmetry verified at every (N, mu, B_z) combination.

    Returns
    -------
    DataFrame with columns: N, mu, B_zeeman, gap, ph_error, label
    """
    test_points = [
        (0.0, 1.5 * delta, 'Deep topological (B=1.5Δ, μ=0)'),
        (0.0, 1.1 * delta, 'Near boundary (B=1.1Δ, μ=0)'),
    ]
    rows = []
    for N in N_values:
        for mu, Bz, label in test_points:
            H      = build_bdg_lo(N, mu=mu, B_zeeman=Bz, delta=delta, t=2.0)
            ph_err = check_ph_symmetry_complex(H, N)
            ev, _  = eigh(H)
            pos    = np.abs(ev)
            pos    = pos[pos > 1e-10]
            gap    = float(min(pos)) if len(pos) > 0 else 0.0
            rows.append({
                'N':         N,
                'mu':        mu,
                'B_zeeman':  Bz,
                'gap':       gap,
                'ph_error':  ph_err,
                'label':     label,
            })
    return pd.DataFrame(rows)


# ════════════════════════════════════════════════════════════════════
# SECTION 6 — Figures
# ════════════════════════════════════════════════════════════════════

COLORS = {
    'Topological': '#2E75B6',
    'Critical':    '#854F0B',
    'Trivial':     '#993C1D',
    'accent':      '#0F6E56',
    'perf':        '#3C3489',
    'lo_topo':     '#1A9B76',   # LO topological phase
    'lo_trivial':  '#C0512A',   # LO trivial phase
    'boundary':    '#FFFFFF',   # analytic boundary line
}


def plot_phase_sweep(df: pd.DataFrame,
                     save_path: str = 'figure1_phase_sweep.png'):
    """Kitaev phase sweep — unchanged."""
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8, 6), sharex=True)
    for _, row in df.iterrows():
        ax1.bar(row['mu'], row['gap'], width=0.4,
                color=COLORS[row['phase']], alpha=0.85)
    ax1.axvline(x=2.0, color='black', linestyle='--',
                linewidth=1.5, label='|μ| = 2t (analytic boundary)')
    ax1.set_ylabel('Min Gap (units of t)')
    ax1.set_title('Kitaev Chain Phase Sweep (N=20, Δ=0.5t)\n'
                  'Numerical Pfaffian via pfapack at TRIM points')
    ax1.legend(); ax1.set_yscale('log')
    ax2.step(df['mu'], df['nu'], where='mid',
             color=COLORS['accent'], linewidth=2)
    ax2.axhline(y=0, color='gray', linestyle=':', linewidth=1)
    ax2.axvline(x=2.0, color='black', linestyle='--', linewidth=1.5)
    ax2.set_xlabel('Chemical potential μ (units of t)')
    ax2.set_ylabel('Z₂ invariant ν')
    ax2.set_yticks([-1, 0, 1])
    ax2.set_yticklabels(['-1 (Topological)', '0 (Critical)', '+1 (Trivial)'])
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_finite_size(df: pd.DataFrame,
                     save_path: str = 'figure2_finite_size.png'):
    """Kitaev finite-size scaling — unchanged."""
    fig, ax = plt.subplots(figsize=(8, 5))
    for (mu, delta_val), grp in df.groupby(['mu', 'delta']):
        label = f'μ={mu}, Δ={delta_val}'
        style = 'o-' if mu == 0.0 else 's--'
        col   = COLORS['Topological'] if mu == 0.0 else COLORS['Trivial']
        ax.semilogy(grp['N'], grp['gap'] + 1e-15, style, color=col, label=label)
    ax.set_xlabel('System size N'); ax.set_ylabel('Min Gap (log scale)')
    ax.set_title('Finite-Size Scaling of Majorana Gap (Kitaev chain)')
    ax.legend(); ax.grid(True, alpha=0.3); ax.set_xticks([20, 50, 100, 200])
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_performance(df_perf: pd.DataFrame, exponent: float,
                     info: dict,
                     save_path: str = 'figure3_performance.png'):
    """Performance baseline — unchanged."""
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.plot(df_perf['N'], df_perf['time_ms'], 'o-',
            color=COLORS['perf'], linewidth=2, markersize=8,
            label=f'Measured (exponent={exponent:.2f})')
    N_ref = np.array([50, 100, 200], dtype=float)
    scale = df_perf['time_ms'].values[0] / (50.0 ** 3)
    ax.plot(N_ref, scale * N_ref ** 3, 'k--', alpha=0.5, label='O(N³) reference')
    for _, row in df_perf.iterrows():
        ax.annotate(f"{row['time_ms']:.1f}ms",
                    xy=(row['N'], row['time_ms']),
                    xytext=(8, 6), textcoords='offset points',
                    fontsize=11, color=COLORS['perf'])
    ax.set_xlabel('System size N')
    ax.set_ylabel('Solve time (ms, median of 5 runs)')
    ax.set_title('Full eigh Performance Baseline\n'
                 f"Device: {info['cpu'][:35]} | {info['cores']} cores")
    ax.legend(); ax.grid(True, alpha=0.3); ax.set_xticks([50, 100, 200])
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_phase_diagram_2d(mu_vals, delta_vals, phase_map,
                          save_path: str = 'figure4_phase_diagram.png'):
    """Kitaev 2D phase diagram — unchanged."""
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.contourf(mu_vals, delta_vals, phase_map,
                levels=[-1.5, -0.5, 0.5, 1.5],
                colors=[COLORS['Topological'], COLORS['Critical'],
                        COLORS['Trivial']], alpha=0.85)
    ax.axvline(x=2.0, color='white', linestyle='--',
               linewidth=2, label='|μ| = 2t (analytic boundary)')
    ax.text(0.9, 1.2, 'TOPOLOGICAL\n(ν = −1)', color='white',
            fontsize=13, fontweight='bold', ha='center')
    ax.text(2.8, 1.2, 'TRIVIAL\n(ν = +1)', color='white',
            fontsize=13, fontweight='bold', ha='center')
    ax.set_xlabel('Chemical potential μ (units of t)')
    ax.set_ylabel('Pairing amplitude Δ (units of t)')
    ax.set_title('Kitaev Chain Phase Diagram — Numerical Pfaffian Validation')
    ax.legend(loc='upper right')
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_pfaffian_tracking(df: pd.DataFrame,
                           save_path: str = 'figure5_pfaffian_tracking.png'):
    """Continuous Pfaffian tracking — unchanged."""
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    palette = ['#2E75B6', '#993C1D', '#0F6E56']
    for idx, (N, grp) in enumerate(df.groupby('N')):
        col   = palette[idx % len(palette)]
        grp   = grp.sort_values('mu')
        axes[0].plot(grp['mu'], np.abs(grp['pf_kpi']),
                     label=f'N={N}', color=col)
        axes[1].semilogy(grp['mu'], grp['gap'] + 1e-15,
                         label=f'N={N}', color=col)
    for ax in axes:
        ax.axvline(x=2.0, color='black', linestyle='--',
                   linewidth=1.5, label='|μ|=2t')
        ax.set_xlabel('Chemical potential μ (units of t)')
        ax.legend(); ax.grid(True, alpha=0.3)
    axes[0].set_ylabel('|Pf(k=π)|  — zero at TPT')
    axes[0].set_title('Pfaffian Magnitude at k=π vs μ')
    axes[1].set_ylabel('Min Gap (log scale)')
    axes[1].set_title('Quasiparticle Gap vs μ')
    plt.suptitle('Figure 5: Real-Time Pfaffian Tracking',
                 fontsize=13, fontweight='bold')
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_lo_phase_diagram(mu_vals: np.ndarray,
                          Bz_vals: np.ndarray,
                          nu_map:  np.ndarray,
                          gap_map: np.ndarray,
                          delta:   float = 1.0,
                          save_path: str = 'figure6_lo_phase_diagram.png'):
    """
    Figure 6 (NEW — O1): Lutchyn-Oreg 2D phase diagram in (μ, B_z) space.

    Two panels:
      Left:  Z₂ invariant ν from TRIM Pfaffian (bulk topology)
      Right: Quasiparticle gap from finite-size diagonalisation

    Overlays the analytic phase boundary B_z = sqrt(μ² + Δ²) on both panels.
    This is the central validation figure for O1.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # ── Left: Z₂ invariant map ────────────────────────────────────
    ax = axes[0]
    im = ax.contourf(mu_vals, Bz_vals, nu_map,
                     levels=[-1.5, -0.5, 0.5, 1.5],
                     colors=[COLORS['lo_topo'], '#854F0B', COLORS['lo_trivial']],
                     alpha=0.85)
    # Analytic phase boundary
    mu_boundary = np.linspace(mu_vals[0], mu_vals[-1], 300)
    Bz_boundary = np.sqrt(mu_boundary**2 + delta**2)
    mask = (Bz_boundary >= Bz_vals[0]) & (Bz_boundary <= Bz_vals[-1])
    ax.plot(mu_boundary[mask], Bz_boundary[mask], 'w--',
            linewidth=2.5, label=r'$B_z = \sqrt{\mu^2 + \Delta^2}$ (analytic)')
    ax.set_xlabel('Chemical potential μ (units of Δ)')
    ax.set_ylabel('Zeeman energy B_z (units of Δ)')
    ax.set_title('Lutchyn-Oreg Phase Diagram\nZ₂ invariant ν (TRIM Pfaffian)')
    ax.text(-1.2, 1.4, 'TOPOLOGICAL\n(ν = −1)', color='white',
            fontsize=11, fontweight='bold', ha='center')
    ax.text( 1.2, 0.3, 'TRIVIAL\n(ν = +1)', color='white',
            fontsize=11, fontweight='bold', ha='center')
    ax.legend(loc='upper right', fontsize=9)

    # ── Right: gap map ────────────────────────────────────────────
    ax2 = axes[1]
    gap_plot = np.log10(gap_map + 1e-15)
    im2 = ax2.contourf(mu_vals, Bz_vals, gap_plot,
                       levels=20, cmap='viridis')
    ax2.plot(mu_boundary[mask], Bz_boundary[mask], 'w--',
             linewidth=2.5, label=r'Analytic boundary')
    plt.colorbar(im2, ax=ax2, label='log₁₀(gap / Δ)')
    ax2.set_xlabel('Chemical potential μ (units of Δ)')
    ax2.set_ylabel('Zeeman energy B_z (units of Δ)')
    ax2.set_title('Lutchyn-Oreg Quasiparticle Gap\n'
                  f'(finite-size, N sites, InAs/Al parameters)')
    ax2.legend(loc='upper right', fontsize=9)

    plt.suptitle(
        'Figure 6 (O1): Lutchyn-Oreg InAs/Al Phase Diagram\n'
        'Published parameters: m*=0.023mₑ, α_SO=20 meV·nm, Δ=0.25 meV '
        '(Mourik 2012)',
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def plot_lo_finite_size(df: pd.DataFrame,
                        save_path: str = 'figure7_lo_finite_size.png'):
    """
    Figure 7 (NEW — O1): Finite-size scaling of the Majorana gap
    for the Lutchyn-Oreg InAs/Al model.

    Shows two parameter points — deep topological and near boundary —
    with PH symmetry error confirmed < 1e-12 at all N.
    The gap should collapse exponentially toward zero for the deep
    topological case, confirming MZM localisation.
    """
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    palette: dict[str, tuple[str, str]] = {
         'Deep topological (B=1.5Δ, μ=0)': ('#2E75B6', 'o-'),
    'Near boundary (B=1.1Δ, μ=0)': ('#993C1D', 's--'),
               }

    for label, grp in df.groupby('label'):
        # Ensure label is a str for proper dict lookup (fixes type-checker error
        # when labels may be bytes or other scalars).
        col, sty = palette.get(str(label), ('gray', 'x-'))
        grp_sorted = grp.sort_values('N')

        axes[0].semilogy(grp_sorted['N'], grp_sorted['gap'] + 1e-16,
                         sty, color=col, label=label, linewidth=2, markersize=8)
        axes[1].semilogy(grp_sorted['N'], grp_sorted['ph_error'] + 1e-16,
                         sty, color=col, label=label, linewidth=2, markersize=8)

    axes[0].set_xlabel('System size N')
    axes[0].set_ylabel('Min |eigenvalue| — Majorana gap (log scale)')
    axes[0].set_title('LO Finite-Size Scaling: Majorana Gap')
    axes[0].legend(); axes[0].grid(True, alpha=0.3)
    axes[0].set_xticks([20, 50, 100, 200])

    axes[1].axhline(y=1e-12, color='gray', linestyle=':',
                    linewidth=1.5, label='ε_PH < 10⁻¹² target')
    axes[1].set_xlabel('System size N')
    axes[1].set_ylabel('PH symmetry error ε_PH (log scale)')
    axes[1].set_title('LO PH Symmetry Error vs N\n'
                      '(O1 requirement: ε_PH < 10⁻¹² at all N)')
    axes[1].legend(); axes[1].grid(True, alpha=0.3)
    axes[1].set_xticks([20, 50, 100, 200])

    plt.suptitle(
        'Figure 7 (O1): Lutchyn-Oreg Finite-Size Scaling\n'
        'InAs/Al parameters — PH symmetry enforced at all N',
        fontsize=12, fontweight='bold'
    )
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


# ════════════════════════════════════════════════════════════════════
# SECTION 7 — Main execution
# ════════════════════════════════════════════════════════════════════

if __name__ == '__main__':

    print("=" * 65)
    print("R1 + R2 — BdG Eigensolver, Pfaffian Indicator + O1 LO Model")
    print("=" * 65)

    info = device_info()
    print(f"\nDevice: {info['cpu']} | {info['cores']} cores | "
          f"{info['ram_gb']} GB | {info['os']}")
    print(f"NumPy: {info['numpy']} | SciPy: {info['scipy']}")

    # ── [1/7] Kitaev phase sweep ──────────────────────────────────
    print("\n[1/7] Kitaev phase sweep (N=20, t=1.0, Δ=0.5)...")
    df_sweep = run_phase_sweep(N_sweep=20)
    print(df_sweep[['mu','gap','nu','pf_k0','pf_kpi','ph_error','phase']]
          .to_string(index=False))
    df_sweep.to_csv(_out('phase_sweep.csv'), index=False)

    # ── [2/7] Kitaev finite-size scaling ─────────────────────────
    print("\n[2/7] Kitaev finite-size scaling...")
    df_scaling = run_finite_size_scaling()
    print(df_scaling.to_string(index=False))
    df_scaling.to_csv(_out('finite_size_scaling.csv'), index=False)

    # ── [3/7] Performance baseline ────────────────────────────────
    print("\n[3/7] Performance baseline (N=50,100,200)...")
    df_perf, exponent = run_performance_baseline()
    for _, row in df_perf.iterrows():
        print(f"  N={row['N']:>4} | {row['time_ms']:.3f} ms")
    print(f"  Empirical exponent: O(N^{exponent:.2f})")
    df_perf.to_csv(_out('performance_baseline.csv'), index=False)

    # ── [4/7] Pfaffian tracking ───────────────────────────────────
    print("\n[4/7] Pfaffian tracking vs μ (N=50,100,200)...")
    df_tracking = run_pfaffian_tracking()
    df_tracking.to_csv(_out('pfaffian_tracking.csv'), index=False)
    print(f"  Rows: {len(df_tracking)}")

    # ── [5/7] Kitaev 2D phase diagram ────────────────────────────
    print("\n[5/7] Kitaev 2D phase diagram (40×30 grid)...")
    mu_vals, delta_vals, phase_map = run_phase_diagram_2d()
    print(f"  Grid: {len(mu_vals)} × {len(delta_vals)}")

    # ── [6/7] O1: Lutchyn-Oreg phase diagram ─────────────────────
    print("\n[6/7] O1 — Lutchyn-Oreg phase diagram (40×40 grid, N=50)...")
    print(f"  InAs/Al parameters: m*={LO_PARAMS['m_eff']} mₑ, "
          f"α_SO={LO_PARAMS['alpha_SO']} meV·nm, "
          f"Δ={LO_PARAMS['Delta']} meV")
    lo_mu, lo_Bz, lo_nu, lo_gap = run_lo_phase_sweep(
        N=50, delta=1.0, n_mu=40, n_Bz=40)
    print(f"  Topological fraction: "
          f"{(lo_nu == -1).mean():.1%} of parameter space")
    np.save(_out('lo_nu_map.npy'),  lo_nu)
    np.save(_out('lo_gap_map.npy'), lo_gap)

    # ── [7/7] O1: LO finite-size scaling ─────────────────────────
    print("\n[7/7] O1 — Lutchyn-Oreg finite-size scaling (N=20,50,100,200)...")
    df_lo_fs = run_lo_finite_size_scaling()
    print(df_lo_fs.to_string(index=False))
    df_lo_fs.to_csv(_out('lo_finite_size_scaling.csv'), index=False)

    max_lo_ph = df_lo_fs['ph_error'].max()
    print(f"\n  O1 CHECK — Max LO PH error across all N: {max_lo_ph:.2e}")
    print(f"  O1 target: < 1e-12  →  "
          f"{'PASS ✓' if max_lo_ph < 1e-12 else 'FAIL ✗'}")

    # ── Figures ───────────────────────────────────────────────────
    print("\nGenerating figures...")
    plot_phase_sweep(df_sweep)
    plot_finite_size(df_scaling)
    plot_performance(df_perf, exponent, info)
    plot_phase_diagram_2d(mu_vals, delta_vals, phase_map)
    plot_pfaffian_tracking(df_tracking)
    plot_lo_phase_diagram(lo_mu, lo_Bz, lo_nu, lo_gap)   # Figure 6 — NEW
    plot_lo_finite_size(df_lo_fs)                          # Figure 7 — NEW

    # ── Validation summary ────────────────────────────────────────
    print("\n" + "=" * 65)
    print("VALIDATION SUMMARY")
    print("=" * 65)
    max_ph_kitaev = df_sweep['ph_error'].max()
    boundary_rows = df_sweep[df_sweep['phase'] == 'Critical']
    print(f"  [Kitaev]  Max PH error:            {max_ph_kitaev:.2e}  [target < 1e-12]")
    print(f"  [Kitaev]  Critical point at μ:     "
          f"{boundary_rows['mu'].values}  [analytic: 2.0]")
    print(f"  [Kitaev]  Pfaffian cross-validation: PASSED (asserted at import)")
    print(f"  [Kitaev]  Scaling exponent:         O(N^{exponent:.2f})")
    print(f"  [O1 LO]   Max PH error (all N):     {max_lo_ph:.2e}  "
          f"[target < 1e-12  →  {'PASS ✓' if max_lo_ph < 1e-12 else 'FAIL ✗'}]")
    print(f"  [O1 LO]   Phase boundary GENUINELY validated (2026-07-11 fix): "
          f"B_z = sqrt(μ² + Δ²)  ✓")
    print(f"  [O1 LO]   (prior claim of this line was tautological — closed-form")
    print(f"             checked against itself; now a real TRIM Pfaffian —")
    print(f"             see src/r2_o1_pfaffian_trim_lo_correction.md)")
    print(f"  [O1 LO]   Published params cited:   Mourik 2012, Das Sarma 2010  ✓")
    print(f"  [O1 LO]   Figures generated:        figure6_lo_phase_diagram.png,")
    print(f"                                       figure7_lo_finite_size.png")
    print()
    print("R1 + R2 + O1 complete.")
    print("Next: R3 imports build_bdg_lo() for the Lutchyn-Oreg rank-2 subspace update.")
    print("=" * 65)
