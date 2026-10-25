"""
ISEF 2027 -- Symmetry-Aware BdG Eigensolver
Phase: R6 -- Full Lutchyn-Oreg phase diagram, realistic InAs/Al parameters

Project: Inductive Circuit Feedback Shifts the Topological Phase Boundary
         in InAs/Al Majorana Nanowires: A Self-Consistent Symmetry-Preserving
         Eigensolver Approach

R6 Phase 1 reconnaissance found:
  - build_bdg_lo() (R1/R2) already implements a genuine 4N x 4N spin (x)
    Nambu Lutchyn-Oreg Hamiltonian, including a real Rashba SOC term. R6 is
    therefore NOT a 2N->4N structural extension -- it is a realistic-
    parameters + genuine-invariant + target-figure exercise.
  - Neither Lutchyn-Sau-Das Sarma (arXiv:1002.4033, PRL 105, 077001) nor
    Oreg-Refael-von Oppen (arXiv:1003.1145, PRL 105, 177002) contains a 2D
    (mu, B) phase diagram -- that is later literature. ORV Figure 1 (BdG
    spectrum vs Zeeman field, gap closing at B=sqrt(mu^2+Delta^2)) was
    selected as the tractable, structurally-compatible target: it maps
    directly onto build_bdg_lo with no new geometry required.
  - pfaffian_trim_lo() (R1/R2) was found to discard its own TRIM matrices
    and return a hardcoded closed form; fixed in a standalone commit before
    this module was written (see docs/methods/r2_o1_pfaffian_trim_lo_correction.md).
    This module uses the corrected, genuine version throughout.
  - LO_PARAMS' Delta=0.25 meV traces to Mourik et al. 2012 (InSb, not
    InAs/Al). Corrected here to Delta=0.19 meV, the full-shell value
    directly measured in Chang et al. 2015 (arXiv:1411.6255), the actual
    InAs/Al epitaxial reference.
  - A verified, unresolved tension: the Zeeman field needed to reach the
    topological transition at these parameters (~0.625 T at mu=0) exceeds
    the verified axial critical field of a comparable thin epitaxial-Al
    shell (~0.075 T, Chang et al. 2015) by ~8x. See Section 1 below and
    docs/methods/r6_lo_phase_diagram_methods.md for the full disclaimer -- every
    figure in this module that sweeps B_zeeman into that regime is an
    explicitly idealized/theoretical parameter regime, not a claim of
    achievable-in-one-device realism.

This module is purely additive: it imports build_bdg_lo, pfaffian_trim_lo,
check_ph_symmetry_complex read-only from R1/R2, and build_ring_lo read-only
from R4 (already validated there; independently cross-checked against
build_bdg_lo's own hole-block convention to machine precision as part of
this module's development). No R1-R5 file is modified.
"""

import os
import numpy as np
from scipy.linalg import eigh
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "results", "lutchyn_oreg_phase_diagram")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


from bdg_eigensolver_r1_r2_LOAdapted import (   # noqa: E402
    build_bdg_lo,
    check_ph_symmetry_complex,
    pfaffian_trim_lo,
    COLORS,
)
from bdg_r4_em_coupling import build_ring_lo    # noqa: E402


# ════════════════════════════════════════════════════════════════════
# SECTION 1 -- R6 realistic InAs/Al parameters (corrected Delta; NOT a
# modification of the shared LO_PARAMS dict used by R1-R5)
# ════════════════════════════════════════════════════════════════════

_HBAR_SI = 1.054571817e-34
_M_E_SI  = 9.1093837015e-31
_MEV_J   = 1.602176634e-19 * 1e-3
_MU_B_MEV_PER_T = 5.7883818060e-2   # Bohr magneton, meV/T (CODATA)

R6_PARAMS = {
    # m*: bulk InAs conduction-band effective mass. [VERIFIED] Ioffe NSM
    #     Semiconductors Archive (ioffe.ru/SVA/NSM/Semicond/InAs/bandstr.html),
    #     "m_Gamma = 0.023 m_0"; consistent with Vurgaftman, Meyer, Ram-Mohan,
    #     J. Appl. Phys. 89, 5815 (2001) (not independently re-opened this
    #     session due to PDF-rendering limits, but the Ioffe value matches).
    'm_eff':    0.023,
    # g-factor: chosen within the verified InAs range, NOT sourced to
    #     Mourik et al. 2012 (that experiment used InSb, not InAs -- the
    #     old LO_PARAMS citation was a material mismatch, found in R6
    #     Phase 1). [VERIFIED] range: ORV (arXiv:1003.1145) gives InAs
    #     g~8; LSDS (arXiv:1002.4033) gives InAs g<~35.
    'g_factor': 10.5,
    # alpha_SO: Rashba coupling, meV*nm. [VERIFIED] order-of-magnitude
    #     match to LSDS's alpha~0.1 eV*Angstrom = 10 meV*nm (this value is
    #     ~2x that estimate; kept as in the original LO_PARAMS, not
    #     re-derived, since Phase 1 did not flag this as needing
    #     correction beyond noting the factor of ~2).
    'alpha_SO': 20.0,
    # Delta: induced gap, meV. [VERIFIED, CORRECTED] Chang et al. 2015
    #     (arXiv:1411.6255), the actual InAs/Al epitaxial reference:
    #     full-shell devices, 190 ueV = 0.19 meV. Replaces the old
    #     LO_PARAMS value of 0.25 meV, which traced to Mourik 2012 (InSb).
    'Delta':    0.19,
    'a_lat':    10.0,
    # Axial critical field of a comparable thin epitaxial-Al shell.
    # [VERIFIED] Chang et al. 2015: "External magnetic field, B, was
    # applied along the nanowire axis... Bc~75 mT for epitaxial [devices]."
    'Bc_axial_measured_T': 0.075,
}

_m, _a, _D = R6_PARAMS['m_eff'], R6_PARAMS['a_lat'], R6_PARAMS['Delta']
R6_PARAMS['t_hop']  = (_HBAR_SI**2) / (2*_m*_M_E_SI*(_a*1e-9)**2) / _MEV_J
# alpha_so passed to build_bdg_lo/build_ring_lo is the per-bond Rashba
# ENERGY (meV) = alpha_SO(meV*nm) / a_lat(nm) -- independent of Delta.
R6_PARAMS['alpha_phys'] = R6_PARAMS['alpha_SO'] / _a
# Zeeman energy at the mu=0 topological boundary, and its Tesla equivalent
# via B_zeeman = g*mu_B*B/2.
R6_PARAMS['Bc_meV_mu0'] = _D
R6_PARAMS['Bc_tesla_mu0'] = _D / (R6_PARAMS['g_factor'] * _MU_B_MEV_PER_T / 2)
R6_PARAMS['idealized_regime_ratio'] = (
    R6_PARAMS['Bc_tesla_mu0'] / R6_PARAMS['Bc_axial_measured_T']
)

IDEALIZED_REGIME_CAPTION = (
    "Idealized/theoretical regime, not a reproduction of an achievable "
    "single-device operating point: the model holds Delta fixed (no "
    "field-dependent pair-breaking Delta(B)), and the plotted Zeeman range "
    "exceeds the verified axial critical field (~0.075 T, Chang et al. 2015) "
    "of a comparable thin epitaxial-Al shell by "
    f"~{R6_PARAMS['idealized_regime_ratio']:.1f}x at mu=0. See "
    "docs/methods/r6_lo_phase_diagram_methods.md."
)

# Short form for placing the Delta(B) caveat directly on the physics panel of
# Figure 1 (the gap-closing figure), where the topological-transition claim is
# actually made -- same "state the caveat where the claim is" discipline used
# for the g=3e5 dynamical-reachability note in R5.
DELTA_B_CAVEAT_SHORT = (
    "Caveat: Delta held fixed\n(no Delta(B) pair-breaking).\n"
    "Idealized transition; a real\nthin Al shell loses its gap\n"
    "near ~0.075 T (Chang 2015),\nwell below B_c here."
)


def zeeman_meV_to_tesla(B_meV: float) -> float:
    """Convert a Zeeman energy (meV) to the applied field (T) via
    B_zeeman = g*mu_B*B/2, at R6_PARAMS' g-factor."""
    return B_meV / (R6_PARAMS['g_factor'] * _MU_B_MEV_PER_T / 2)


def _bloch_H_lo(k: float, mu: float, B: float, delta: float,
                alpha: float, t: float) -> np.ndarray:
    """
    2-band (4x4) Bloch BdG Hamiltonian for the clean Lutchyn-Oreg ring, in
    the (up, dn, up-hole, dn-hole) basis at momentum k:

        h(k)   = [[eps_k+B, -i*alpha*sin k], [ i*alpha*sin k, eps_k-B]],
                 eps_k = (2t-mu) - 2t*cos k
        H(k)   = [[ h(k),        delta*i*sigma_y ],
                  [ (delta*i*sigma_y)^dag, -h(k) ]]

    The hole block is -h(k), which equals the proper BdG -h*(-k) because the
    real-space ring h is real (verified). This construction was derived by
    directly Fourier-transforming build_ring_lo's blocks and validated to
    reproduce build_ring_lo's full spectrum to <1e-12 across 40 random
    (N, mu, B, delta, alpha, t) draws at arbitrary alpha (an earlier
    hand-derived version had two sign errors -- SOC off-diagonal, and
    -h(-k) vs -h(k) -- that happened to partially cancel; this version is
    the FT-verified one). Used only for the continuum-k validation below;
    the figures use build_ring_lo directly.
    """
    eps = (2.0 * t - mu) - 2.0 * t * np.cos(k)
    h = np.array([[eps + B, -1j * alpha * np.sin(k)],
                  [1j * alpha * np.sin(k), eps - B]], dtype=complex)
    D = delta * np.array([[0, 1], [-1, 0]], dtype=complex)
    return np.block([[h, D], [D.conj().T, -h]])


# ════════════════════════════════════════════════════════════════════
# SECTION 2 -- Momentum-resolved dispersion (ORV Fig. 1 style)
# ════════════════════════════════════════════════════════════════════

def ring_dispersion(N: int, mu: float, B: float, delta: float,
                    alpha: float, t: float) -> tuple:
    """
    Momentum-resolved BdG dispersion on a clean periodic ring.

    Diagonalises build_ring_lo(N, mu, B, ...) (R4, already validated;
    cross-checked against build_bdg_lo's own hole-block convention to
    machine precision during this module's development) and assigns each
    eigenvector a quasi-momentum k by projecting its PARTICLE-sector
    (first 2N rows: spin-up sites, spin-down sites) component onto plane
    waves exp(-i k n) for the N allowed ring momenta k = 2*pi*j/N. This
    avoids re-deriving a Bloch Hamiltonian by hand (an earlier attempt at
    that had a sign error caught by cross-validation against the ring
    eigenvalues; this approach instead lets the already-validated ring
    Hamiltonian do the momentum-space physics).

    Returns
    -------
    ks     : (N,) allowed momenta
    ev     : (4N,) eigenvalues (unsorted, in eigh order)
    k_of_mode : (4N,) the momentum assigned to each eigenvalue
    """
    H = build_ring_lo(N, mu, B, phi_ext=0.0, delta=delta, t=t, alpha_so=alpha)
    ev, evec = eigh(H)
    ks = 2 * np.pi * np.arange(N) / N
    phases = np.exp(-1j * np.outer(ks, np.arange(N)))   # (N, N)

    v_up = evec[:N, :]     # (N, 4N)
    v_dn = evec[N:2*N, :]  # (N, 4N)
    proj = np.abs(phases @ v_up)**2 + np.abs(phases @ v_dn)**2   # (N, 4N)
    k_of_mode = ks[np.argmax(proj, axis=0)]
    return ks, ev, k_of_mode


# ════════════════════════════════════════════════════════════════════
# SECTION 3 -- Validation against known limits (run before trusting any
# figure below, per this project's established discipline)
# ════════════════════════════════════════════════════════════════════

def validate_zero_zeeman_kramers_degeneracy(N: int = 50, mu: float = 0.3,
                                            delta: float = 1.0,
                                            alpha: float = 8.0,
                                            t: float = 66.26,
                                            B_nonzero: float = 0.4) -> bool:
    """
    Known limit: at B_zeeman=0, time-reversal symmetry is restored, so the
    BdG spectrum must be exactly doubly (Kramers) degenerate. Verified by
    counting near-zero adjacent gaps in the SORTED spectrum (not by
    blindly pairing (0,1),(2,3),...: this Hamiltonian also has an
    unrelated accidental degeneracy from ring k<->-k reflection symmetry,
    which naive fixed-index pairing conflates with the Kramers pairing --
    caught by an earlier version of this check reporting a false FAIL).
    At B=0, essentially every eigenvalue should be part of a near-zero-gap
    pair (full double degeneracy); at B!=0, that count must strictly
    decrease (Kramers-protected pairs split; the unrelated
    reflection-symmetry degeneracy may persist and is not expected to
    vanish).
    """
    H0 = build_ring_lo(N, mu, 0.0, phi_ext=0.0, delta=delta, t=t, alpha_so=alpha)
    ev0 = np.sort(eigh(H0, eigvals_only=True))
    n_degen0 = int(np.sum(np.diff(ev0) < 1e-8))

    Hb = build_ring_lo(N, mu, B_nonzero, phi_ext=0.0, delta=delta, t=t, alpha_so=alpha)
    evb = np.sort(eigh(Hb, eigvals_only=True))
    n_degenb = int(np.sum(np.diff(evb) < 1e-8))

    full_degeneracy_at_zero = n_degen0 == 2 * N   # 4N eigenvalues -> 2N pairs
    breaks_at_nonzero = n_degenb < n_degen0

    ok = full_degeneracy_at_zero and breaks_at_nonzero
    print(f"  [R6 validation] Zero-Zeeman Kramers degeneracy: "
          f"near-zero adjacent gaps at B=0: {n_degen0}/{4*N-1} "
          f"(want exactly {2*N}, full double degeneracy); "
          f"at B={B_nonzero}: {n_degenb} (want fewer -- Kramers pairs split) "
          f"-> {'PASS' if ok else 'FAIL'}")
    return ok


def validate_zero_soc_gap_closes(delta: float = 1.0, t: float = 66.26,
                                 test_points=((0.0, 1.5), (0.0, 1.2),
                                              (0.3, 1.3), (-0.5, 1.5),
                                              (0.7, 1.5), (0.0, 2.0)),
                                 n_k: int = 20000,
                                 gap_tol: float = 5e-3) -> bool:
    """
    Known limit: without spin-orbit coupling, s-wave singlet pairing plus a
    pure Zeeman splitting cannot open a bulk gap away from the TRIM points
    -- SOC is essential for a genuinely non-degenerate (protected)
    topological phase in a 1D BdG chain. At alpha=0, the naive TRIM-only
    invariant (pfaffian_trim_lo, SOC-independent by construction since
    sin(k)=0 at TRIM) would still flip sign at B > sqrt(mu^2+delta^2), yet
    the actual bulk gap must be ZERO in the continuum everywhere in that
    'topological' region.

    This is checked as a GENERAL property, at MULTIPLE (mu, B) points that
    are all topological by the TRIM criterion, via a continuous-k scan of
    the FT-verified Bloch Hamiltonian _bloch_H_lo (NOT via a fixed ring
    N-grid: a finite ring can miss the exact band-crossing momentum k*,
    leaving an apparent finite gap that is purely a discretisation
    artifact -- an earlier single-point ring N-scaling version of this
    check gave a spurious 'not collapsing' at (mu=0, B=2.0) for exactly
    that reason, where the true continuum gap is ~1e-3). The residual gap
    on a fine continuous grid is O(dk) near the linear band crossing, hence
    the gap_tol threshold rather than exactly zero.
    """
    ks = np.linspace(0.0, np.pi, n_k)
    all_ok = True
    reports = []
    for mu, B in test_points:
        Bc = np.sqrt(mu**2 + delta**2)
        if B <= Bc:
            raise ValueError(f"test point (mu={mu}, B={B}) is not topological "
                             f"by the TRIM criterion (B must exceed Bc={Bc:.3f})")
        gap = min(float(np.min(np.abs(eigh(_bloch_H_lo(k, mu, B, delta, 0.0, t),
                                           eigvals_only=True)))) for k in ks)
        pt_ok = gap < gap_tol
        all_ok = all_ok and pt_ok
        reports.append(f"(mu={mu:+.1f},B={B:.1f}):{gap:.1e}")

    # Secondary illustration: on the ACTUAL ring, the gap at a representative
    # point also shrinks toward zero with N (once N samples near k*).
    ring_gaps = [float(np.min(np.abs(eigh(
        build_ring_lo(N, 0.0, 1.5, phi_ext=0.0, delta=delta, t=t, alpha_so=0.0),
        eigvals_only=True)))) for N in (80, 320, 640)]

    print(f"  [R6 validation] Zero-SOC bulk gap == 0 (continuum k-scan, "
          f"{len(test_points)} topological points, all should be < {gap_tol}): "
          + ", ".join(reports)
          + f" -> {'PASS' if all_ok else 'FAIL'}"
          + f"  [ring cross-check, gap vs N at (0,1.5): "
          + ", ".join(f"N={n}:{g:.3f}" for n, g in zip((80, 320, 640), ring_gaps)) + "]")
    return all_ok


# ════════════════════════════════════════════════════════════════════
# SECTION 4 -- Figures
# ════════════════════════════════════════════════════════════════════

def figure_orv_dispersion(mu: float = 0.0, save_path: str = 'r6_figure1_orv_dispersion.png'):
    """
    Reproduces the qualitative content of Oreg-Refael-von Oppen Figure 1
    (gap closing and reopening across the topological transition) using
    THIS codebase's own Hamiltonian and realistic InAs/Al parameters --
    not a copy of the published figure. Three panels: trivial
    (B=0.5*Bc), critical (B=Bc), topological (B=1.5*Bc), all at mu=0,
    momentum-resolved via ring_dispersion().

    At these realistic parameters t (~16.5 meV) is far larger than
    Delta/B (~0.1-0.4 meV), so the relevant low-energy physics occupies a
    narrow window of the lattice Brillouin zone near k=0 (band minimum) --
    the ORV papers work directly in that continuum limit. N=300 (rather
    than the ~100 used elsewhere in this module) is used here specifically
    so enough ring momenta fall inside that narrow window to resolve the
    dispersion; the k-axis is re-centred to (-k, k) around the band
    minimum and zoomed accordingly, rather than showing the full
    (mostly-irrelevant-at-this-energy-scale) [0, 2*pi) lattice BZ.
    """
    delta, alpha, t = R6_PARAMS['Delta'], R6_PARAMS['alpha_phys'], R6_PARAMS['t_hop']
    Bc = np.sqrt(mu**2 + delta**2)
    N = 300
    k_window = 4.0 * np.sqrt(2.0 * Bc / t)   # a few times the natural k-scale near k=0

    fig, axes = plt.subplots(1, 3, figsize=(15, 5), sharey=True)
    panels = [(0.5*Bc, 'Trivial  (B=0.5 Bc)'), (Bc, 'Critical  (B=Bc)'),
              (1.5*Bc, 'Topological  (B=1.5 Bc)')]
    for ax, (B, title) in zip(axes, panels):
        ks, ev, kmode = ring_dispersion(N, mu, B, delta, alpha, t)
        k_centered = np.where(kmode > np.pi, kmode - 2*np.pi, kmode)   # (-pi, pi]
        ax.scatter(k_centered, ev, s=4, color=COLORS['Topological'], alpha=0.6)
        ax.axhline(0, color='gray', linestyle=':', linewidth=1)
        ax.set_xlabel('k (near band minimum)')
        ax.set_title(f'{title}\nB={B:.4f} meV')
        ax.set_xlim(-k_window, k_window)
        ax.set_ylim(-3*delta, 3*delta)
    axes[0].set_ylabel('E (meV)')
    # State the Delta(B) caveat directly on the middle (critical) panel, where
    # the gap-closing/reopening topological claim is actually made -- not only
    # in the footnote several lines below.
    axes[1].text(0.98, 0.97, DELTA_B_CAVEAT_SHORT, transform=axes[1].transAxes,
                 ha='right', va='top', fontsize=7, style='italic',
                 bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.85))
    plt.suptitle(
        'R6 Figure 1 -- Lutchyn-Oreg BdG spectrum vs Zeeman field '
        '(mu=0, realistic InAs/Al parameters)\n'
        f'InAs m*={R6_PARAMS["m_eff"]}me, alpha_SO={R6_PARAMS["alpha_SO"]} meV*nm, '
        f'Delta={delta} meV (Chang et al. 2015)  |  '
        'Qualitative reproduction of Oreg-Refael-von Oppen (arXiv:1003.1145) Fig. 1 '
        '-- own regenerated plot, not the published figure',
        fontsize=11, fontweight='bold'
    )
    fig.text(0.5, -0.02, IDEALIZED_REGIME_CAPTION, ha='center', fontsize=8,
             wrap=True, style='italic')
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def run_gap_vs_Bz_sweep(mu_values=(0.0, 0.3), N: int = 100,
                        n_Bz: int = 60) -> pd.DataFrame:
    """
    Sweeps B_zeeman at fixed mu, tracking the ring bulk gap AND the
    genuine TRIM Pfaffian invariant, to show spectral gap-closing and the
    topological invariant flip coincide (as they must at N=100, k=0 is the
    protected momentum -- see r2_o1_pfaffian_trim_lo_correction.md).
    """
    delta, alpha, t = R6_PARAMS['Delta'], R6_PARAMS['alpha_phys'], R6_PARAMS['t_hop']
    rows = []
    for mu in mu_values:
        Bc = np.sqrt(mu**2 + delta**2)
        Bz_vals = np.linspace(0.0, 2.0*Bc, n_Bz)
        for Bz in Bz_vals:
            H = build_ring_lo(N, mu, Bz, phi_ext=0.0, delta=delta, t=t, alpha_so=alpha)
            ev = eigh(H, eigvals_only=True)
            gap = float(np.min(np.abs(ev)))
            nu, _, _ = pfaffian_trim_lo(mu, Bz, delta, t=t)
            rows.append({'mu': mu, 'B_zeeman': Bz, 'Bc_analytic': Bc,
                        'gap': gap, 'nu': nu})
    return pd.DataFrame(rows)


def figure_gap_vs_Bz(df: pd.DataFrame, save_path: str = 'r6_figure2_gap_vs_Bz.png'):
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=False)
    ax1, ax2 = axes
    for mu, grp in df.groupby('mu'):
        grp = grp.sort_values('B_zeeman')
        ax1.plot(grp['B_zeeman'], grp['gap'], 'o-', markersize=3,
                 label=f'mu={mu} meV')
        ax1.axvline(grp['Bc_analytic'].iloc[0], color='gray', linestyle=':', linewidth=1)
    ax1.set_xlabel('B_zeeman (meV)')
    ax1.set_ylabel('Bulk gap (meV)')
    ax1.set_title('Ring bulk gap vs Zeeman field (dotted: analytic Bc=sqrt(mu^2+Delta^2))')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    for mu, grp in df.groupby('mu'):
        grp = grp.sort_values('B_zeeman')
        ax2.step(grp['B_zeeman'], grp['nu'], where='mid', label=f'mu={mu} meV')
        ax2.axvline(grp['Bc_analytic'].iloc[0], color='gray', linestyle=':', linewidth=1)
    ax2.set_xlabel('B_zeeman (meV)')
    ax2.set_ylabel('nu (genuine TRIM Pfaffian)')
    ax2.set_yticks([-1, 0, 1])
    ax2.set_title('Genuine topological invariant vs Zeeman field')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.suptitle('R6 Figure 2 -- Gap closing and genuine invariant flip coincide',
                 fontsize=12, fontweight='bold')
    fig.text(0.5, -0.01, IDEALIZED_REGIME_CAPTION, ha='center', fontsize=8,
             wrap=True, style='italic')
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


def run_phase_diagram_2d_r6(n_mu: int = 40, n_Bz: int = 40,
                            mu_range: tuple = (-0.6, 0.6),
                            Bz_range: tuple = (0.0, 0.4)) -> tuple:
    """
    2D (mu, B_zeeman) phase diagram at realistic InAs/Al parameters, using
    the GENUINE (post-fix) pfaffian_trim_lo for the invariant and a finite
    ring for the gap map -- the R6-scoped analogue of run_lo_phase_sweep,
    kept separate since the parameter range differs (physical meV, not
    the toy dimensionless t=2.0 sweep in R1/R2).
    """
    delta, alpha, t = R6_PARAMS['Delta'], R6_PARAMS['alpha_phys'], R6_PARAMS['t_hop']
    mu_vals = np.linspace(mu_range[0], mu_range[1], n_mu)
    Bz_vals = np.linspace(Bz_range[0], Bz_range[1], n_Bz)
    nu_map = np.zeros((n_Bz, n_mu))
    gap_map = np.zeros((n_Bz, n_mu))
    N = 60
    for i, Bz in enumerate(Bz_vals):
        for j, mu in enumerate(mu_vals):
            nu, _, _ = pfaffian_trim_lo(float(mu), float(Bz), delta, t=t)
            nu_map[i, j] = nu
            H = build_ring_lo(N, float(mu), float(Bz), phi_ext=0.0,
                              delta=delta, t=t, alpha_so=alpha)
            ev = eigh(H, eigvals_only=True)
            gap_map[i, j] = float(np.min(np.abs(ev)))
    return mu_vals, Bz_vals, nu_map, gap_map


def figure_phase_diagram_2d_r6(mu_vals, Bz_vals, nu_map, gap_map,
                               save_path: str = 'r6_figure3_phase_diagram.png'):
    delta = R6_PARAMS['Delta']
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    ax = axes[0]
    ax.contourf(mu_vals, Bz_vals, nu_map, levels=[-1.5, -0.5, 0.5, 1.5],
               colors=[COLORS['lo_topo'], '#854F0B', COLORS['lo_trivial']], alpha=0.85)
    mu_b = np.linspace(mu_vals[0], mu_vals[-1], 300)
    Bz_b = np.sqrt(mu_b**2 + delta**2)
    mask = (Bz_b >= Bz_vals[0]) & (Bz_b <= Bz_vals[-1])
    ax.plot(mu_b[mask], Bz_b[mask], 'w--', linewidth=2.5,
           label='B=sqrt(mu^2+Delta^2) (analytic; genuinely re-derived, not assumed -- '
                 'see r2_o1_pfaffian_trim_lo_correction.md)')
    ax.set_xlabel('mu (meV)'); ax.set_ylabel('B_zeeman (meV)')
    ax.set_title('R6: genuine TRIM Pfaffian invariant (InAs/Al parameters)')
    ax.legend(loc='upper right', fontsize=7)

    ax2 = axes[1]
    gap_plot = np.log10(gap_map + 1e-15)
    im2 = ax2.contourf(mu_vals, Bz_vals, gap_plot, levels=20, cmap='viridis')
    ax2.plot(mu_b[mask], Bz_b[mask], 'w--', linewidth=2.5)
    plt.colorbar(im2, ax=ax2, label='log10(ring bulk gap / meV)')
    ax2.set_xlabel('mu (meV)'); ax2.set_ylabel('B_zeeman (meV)')
    ax2.set_title('R6: ring bulk gap (N=60), same parameters')

    plt.suptitle(
        'R6 Figure 3 -- Lutchyn-Oreg InAs/Al phase diagram, realistic parameters '
        f'(m*={R6_PARAMS["m_eff"]}me, Delta={delta} meV, Chang et al. 2015)',
        fontsize=12, fontweight='bold'
    )
    fig.text(0.5, -0.02, IDEALIZED_REGIME_CAPTION, ha='center', fontsize=8,
             wrap=True, style='italic')
    plt.tight_layout()
    plt.savefig(_out(save_path), dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {save_path}")


# ════════════════════════════════════════════════════════════════════
# SECTION 5 -- Comparison / discrepancy report
# ════════════════════════════════════════════════════════════════════

def comparison_report(df_sweep: pd.DataFrame) -> dict:
    """
    Numerically compares the observed gap-closing field (from the spectral
    sweep) and the genuine-invariant flip field against the analytic
    B=sqrt(mu^2+Delta^2), reporting any discrepancy plainly rather than
    assuming agreement.
    """
    results = {}
    for mu, grp in df_sweep.groupby('mu'):
        grp = grp.sort_values('B_zeeman').reset_index(drop=True)
        Bc_analytic = float(grp['Bc_analytic'].iloc[0])
        i_min = int(grp['gap'].idxmin())
        Bc_spectral = float(grp['B_zeeman'].iloc[i_min])
        signs = np.sign(grp['nu'].values)
        flips = np.where(np.diff(signs) != 0)[0]
        Bc_invariant = float(grp['B_zeeman'].iloc[flips[0]]) if len(flips) else float('nan')
        grid_step_meV = float(grp['B_zeeman'].iloc[1] - grp['B_zeeman'].iloc[0])
        grid_step_T = zeeman_meV_to_tesla(grid_step_meV)
        results[mu] = {
            'Bc_analytic': Bc_analytic,
            'Bc_spectral_gridres': Bc_spectral,
            'Bc_invariant_gridres': Bc_invariant,
            'spectral_discrepancy': Bc_spectral - Bc_analytic,
            'invariant_discrepancy': Bc_invariant - Bc_analytic,
            'grid_step_meV': grid_step_meV,
            'grid_step_T': grid_step_T,
            'spectral_discrepancy_in_gridsteps': (Bc_spectral - Bc_analytic) / grid_step_meV,
        }
        print(f"  mu={mu:+.3f} meV: analytic Bc={Bc_analytic:.4f} meV | "
              f"spectral gap-min at B={Bc_spectral:.4f} meV | "
              f"invariant flips at B={Bc_invariant:.4f} meV | "
              f"discrepancy {Bc_spectral-Bc_analytic:+.4f} meV "
              f"({(Bc_spectral-Bc_analytic)/grid_step_meV:+.2f} grid steps; "
              f"1 step = {grid_step_meV:.4f} meV = {grid_step_T:.4f} T)")
    return results


# ════════════════════════════════════════════════════════════════════
# SECTION 6 -- Main
# ════════════════════════════════════════════════════════════════════

if __name__ == '__main__':
    print("=" * 65)
    print("R6 -- Full Lutchyn-Oreg phase diagram, realistic InAs/Al parameters")
    print("=" * 65)
    print(f"  m*={R6_PARAMS['m_eff']}me  g={R6_PARAMS['g_factor']}  "
          f"alpha_SO={R6_PARAMS['alpha_SO']} meV*nm  Delta={R6_PARAMS['Delta']} meV "
          f"(Chang et al. 2015)")
    print(f"  t_hop={R6_PARAMS['t_hop']:.3f} meV  alpha_phys={R6_PARAMS['alpha_phys']:.3f} meV")
    print(f"  Bc at mu=0: {R6_PARAMS['Bc_meV_mu0']:.3f} meV = "
          f"{R6_PARAMS['Bc_tesla_mu0']:.3f} T  (measured axial shell Bc: "
          f"{R6_PARAMS['Bc_axial_measured_T']} T -- idealized-regime ratio "
          f"{R6_PARAMS['idealized_regime_ratio']:.1f}x)")

    print("\n[1/5] Validating known limits...")
    ok1 = validate_zero_zeeman_kramers_degeneracy()
    ok2 = validate_zero_soc_gap_closes()
    if not (ok1 and ok2):
        raise RuntimeError("R6 known-limit validation FAILED -- refusing to "
                          "trust downstream figures. See printed diagnostics above.")

    print("\n[2/5] Figure 1 -- ORV-style dispersion (trivial/critical/topological)...")
    figure_orv_dispersion()

    print("\n[3/5] Figure 2 -- gap vs B_zeeman sweep + genuine invariant...")
    df_sweep = run_gap_vs_Bz_sweep()
    df_sweep.to_csv(_out('r6_gap_vs_Bz.csv'), index=False)
    figure_gap_vs_Bz(df_sweep)

    print("\n[4/5] Figure 3 -- 2D (mu, B_zeeman) phase diagram...")
    mu_vals, Bz_vals, nu_map, gap_map = run_phase_diagram_2d_r6()
    np.save(_out('r6_nu_map.npy'), nu_map)
    np.save(_out('r6_gap_map.npy'), gap_map)
    figure_phase_diagram_2d_r6(mu_vals, Bz_vals, nu_map, gap_map)

    print("\n[5/5] Comparison report (spectral / invariant vs analytic Bc)...")
    results = comparison_report(df_sweep)

    print("\n" + "=" * 65)
    print("R6 VALIDATION SUMMARY")
    print("=" * 65)
    print(f"  Zero-Zeeman Kramers degeneracy:  {'PASS' if ok1 else 'FAIL'}")
    print(f"  Zero-SOC bulk gap closes (N-scaling): {'PASS' if ok2 else 'FAIL'}")
    for mu, r in results.items():
        print(f"  mu={mu:+.3f}: spectral discrepancy = {r['spectral_discrepancy']:+.4f} meV, "
              f"invariant discrepancy = {r['invariant_discrepancy']:+.4f} meV")
    print(f"  Idealized-regime disclaimer: needed Bc={R6_PARAMS['Bc_tesla_mu0']:.3f} T "
          f"vs measured axial shell Bc={R6_PARAMS['Bc_axial_measured_T']} T "
          f"({R6_PARAMS['idealized_regime_ratio']:.1f}x)")
    print()
    print("R6 complete.")
    print("=" * 65)
