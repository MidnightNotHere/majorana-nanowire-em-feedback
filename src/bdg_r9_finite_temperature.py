"""
ISEF 2027 — R9 (Track 2): finite temperature in the BdG / Pfaffian / Kubo stack
===============================================================================
Every prior phase (R1-R8, and R9 Track 1) used ZERO-temperature BdG ground
states (sharp Fermi sea, occupation = theta(-E)) and a T=0 Kubo formula. Real
dilution-refrigerator operation is at tens of mK. This module adds finite
temperature as a STRICTLY ADDITIVE layer: it introduces Fermi-Dirac occupation
of the BdG quasiparticle spectrum and re-tests whether the project's null
results survive realistic thermal broadening. Nothing in R1-R8 is modified;
R4/R5/R7/R8A/R8B are imported READ-ONLY.

────────────────────────────────────────────────────────────────────────────
UNITS / THE DIMENSIONLESS TEMPERATURE  theta = k_B T / Delta
────────────────────────────────────────────────────────────────────────────
Throughout R1-R8 all energies are in units of the induced pairing Delta
(delta = 1.0). The R6-verified induced gap is Delta = 0.19 meV (Chang et al.
2015, arXiv:1411.6255, epitaxial-Al full-shell; base temperature 20 mK). So the
natural dimensionless temperature is

        theta  ==  k_B T / Delta          (energies already in units of Delta)

and the Fermi function is f(E) = 1/(exp(E/theta)+1), with f -> theta(-E) as
theta -> 0. `temperature_scale_report()` (Task 1) computes Delta/k_B and the
operating T/Tc explicitly; those numbers drive the theta range swept here.

────────────────────────────────────────────────────────────────────────────
WHERE TEMPERATURE ENTERS (and where it CANNOT)
────────────────────────────────────────────────────────────────────────────
Finite T enters through the OCCUPATION of the BdG quasiparticle spectrum:
  * the self-consistent RING CURRENT (R5/R7 self-field loop) — thermal
    occupation changes the current-carrying state, hence the self-field;
  * the KUBO current-current correlator (R8-B Lehmann sum) — Fermi-Dirac
    weights (f_n - f_m) replace the T=0 step, opening thermally-allowed
    sub-gap channels.

The TOPOLOGICAL INVARIANT (Pfaffian Z2 sign) is a property of the BdG
*Hamiltonian's* ground state (fermion parity), NOT of the occupation — it does
not depend on temperature at fixed Hamiltonian. So at finite T the sharp Z2
transition is unchanged *as a function of the Hamiltonian*; temperature can move
the detected boundary ONLY through the self-consistent self-field texture
phi*(T) at which we evaluate the invariant. We therefore define the finite-T
boundary two complementary ways, both implemented and reported here:
  (A) HARD boundary  — the Pfaffian sign flip evaluated at the finite-T
      self-consistent self-field phi*(T). T enters purely via phi*(T). This is
      the physically meaningful channel for THIS project (does the inductive
      self-field move the bulk boundary — now with a thermal self-field).
  (B) SOFT crossover  — a thermal parity indicator
          P(B) = sgn_Pf(B) * tanh( E_gap(B) / (2 theta) )
      that smears the sharp Z2 into a crossover of width ~theta around the
      gap-closing, honestly acknowledging that at T>0 the free energy is
      analytic in B and the transition is a crossover, not a jump. Its width is
      a NEW finite-T observable with no T=0 analogue.

────────────────────────────────────────────────────────────────────────────
VALIDATION DISCIPLINE
────────────────────────────────────────────────────────────────────────────
Every finite-T quantity has a T=0 counterpart already validated in R5/R7/R8.
`validate_T0_limit()` confirms theta -> 0 reproduces the existing T=0 result
for the current, the self-consistent operating point, the Kubo chi, and the
boundary — BEFORE any nonzero-T number is trusted. Any nonzero-T finding is
held to the same N-independence / path-independence standard as R5's headline.
"""

from __future__ import annotations

import argparse
import contextlib
import io
import json
import os
import time

import numpy as np

# ── Read-only imports from the validated R4/R5/R7/R8 modules ────────────────
with contextlib.redirect_stdout(io.StringIO()):
    from bdg_eigensolver_r1_r2_LOAdapted import LO_PARAMS            # noqa: E402
    from bdg_r5_selfconsistent_ring import (                        # noqa: E402
        _build_ring_lo_phased, ring_current_density,
        ring_biot_savart_self, compute_peierls_from_A,
        eigh, T_HOP_SI, E_CHARGE, HBAR, _RING_KERNEL_CACHE,
    )
    from bdg_r4_em_coupling import _ring_pfaffian_sign               # noqa: E402
    from bdg_r7_inductance_sweep import (                           # noqa: E402
        G_R7, A0, _phi_ext_bond,
    )
    import bdg_r8_reactive_impedance as r8a                         # noqa: E402
    import bdg_r8b_kubo_impedance as r8b                            # noqa: E402
    from bdg_parallel import parallel_map, check_equivalence, default_workers

# ── Physical constants / gap provenance (Task 1) ────────────────────────────
KB_EV = 8.617333262e-5          # Boltzmann constant, eV/K
DELTA_EV = 0.19e-3              # induced gap, eV (Chang et al. 2015, R6-verified)
DELTA_OVER_KB = DELTA_EV / KB_EV      # ~2.205 K, the k_B*Tc ~ Delta scale
T_BASE_CHANG_MK = 20.0         # Chang et al. 2015 dilution-fridge base temp (mK)

G_CRIT_LO = getattr(r8a, "G_CRIT_LO", 1.1e5)   # R5 g_crit band (for pole tests)

_OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                           "..", "ISEF R9 data folder")


def _out(name: str) -> str:
    os.makedirs(_OUTPUT_DIR, exist_ok=True)
    return os.path.join(_OUTPUT_DIR, os.path.basename(str(name)))


# ════════════════════════════════════════════════════════════════════════
# TASK 1 — temperature scale and the realistic operating ratio
# ════════════════════════════════════════════════════════════════════════
def theta_of_T(T_kelvin: float) -> float:
    """Dimensionless temperature theta = k_B T / Delta (energies in units of Delta)."""
    return KB_EV * float(T_kelvin) / DELTA_EV


def T_of_theta(theta: float) -> float:
    """Inverse: physical temperature (K) for a dimensionless theta."""
    return float(theta) * DELTA_EV / KB_EV


def temperature_scale_report() -> dict:
    """Task 1: Delta/k_B, the BCS Tc, and T/Tc + Fermi weight at a range of
    realistic dilution-fridge temperatures. `Tc` here uses the project's stated
    convention k_B*Tc ~ Delta (so Tc = Delta/k_B); the BCS relation Delta =
    1.764 k_B Tc is also reported for reference. The Fermi weight f(E=Delta) =
    1/(exp(1/theta)+1) is the thermal population of a full-gap quasiparticle —
    the single number that decides whether thermal effects can matter."""
    Tc_conv = DELTA_OVER_KB                       # k_B Tc = Delta convention
    Tc_bcs = DELTA_EV / (1.764 * KB_EV)           # BCS: Delta = 1.764 k_B Tc
    rows = []
    for T_mK in (10, 20, 30, 50, 100, 200, 300, 500):
        T = T_mK / 1000.0
        th = theta_of_T(T)
        f_gap = 1.0 / (np.exp(1.0 / th) + 1.0) if th > 0 else 0.0
        rows.append({"T_mK": T_mK, "theta": th,
                     "T_over_Tc_conv": T / Tc_conv,
                     "fermi_at_gap": float(f_gap)})
    return {"Delta_eV": DELTA_EV, "Delta_over_kB_K": Tc_conv, "Tc_BCS_K": Tc_bcs,
            "theta_at_base_20mK": theta_of_T(0.020), "rows": rows}


# ════════════════════════════════════════════════════════════════════════
# TASK 2.0 — Fermi-Dirac occupation
# ════════════════════════════════════════════════════════════════════════
def fermi_dirac(E: np.ndarray, theta: float) -> np.ndarray:
    """Fermi-Dirac occupation f(E) = 1/(exp(E/theta)+1), numerically stable
    (logistic split so neither exp overflows). theta = k_B T / Delta.

    theta <= 0 returns the T=0 step theta(-E) == (E < 0), which is EXACTLY the
    occupation convention used everywhere in R1-R8 (occupied == E_n < 0), so the
    theta -> 0 limit reproduces the T=0 code paths bit-for-bit."""
    E = np.asarray(E, dtype=float)
    if theta <= 0.0:
        return (E < 0.0).astype(float)
    x = E / float(theta)
    out = np.empty_like(x)
    pos = x >= 0.0
    ex_neg = np.exp(-x[pos])
    out[pos] = ex_neg / (1.0 + ex_neg)
    ex_pos = np.exp(x[~pos])
    out[~pos] = 1.0 / (1.0 + ex_pos)
    return out


# ════════════════════════════════════════════════════════════════════════
# TASK 2.1a — finite-T ring current (self-field channel), in R5's Ampere units
# ════════════════════════════════════════════════════════════════════════
def ring_current_density_T(eigenvectors: np.ndarray, eigenvalues: np.ndarray,
                           N: int, t_hop: float, phi_bonds: np.ndarray,
                           theta: float) -> np.ndarray:
    """Finite-T bond current, the thermal generalization of R5
    `ring_current_density`. The T=0 function sums the electron-block bond
    current over occupied (E<0) quasiparticle states; here each state is
    weighted by its Fermi-Dirac occupation f(E_n), summed over the FULL BdG
    spectrum:

        I_bond = (2 e t/hbar) * Σ_n f(E_n) Im[ conj(u_i,n) u_{i+1,n} e^{i phi} ]

    At theta -> 0, f(E_n) = (E_n < 0), so this reduces to R5's sum over occupied
    states BIT-FOR-BIT (same prefactor, same electron-block rule, same wrap
    bond). Same conventions/units (Amperes) as the validated `ring_current_density`,
    so it drops straight into the R5/R7 Biot-Savart self-field loop.

    (`current_expectation_T` gives the operator-form thermal current used by the
    Kubo side; `validate_T0_limit` cross-checks the two are proportional by the
    exact T=0 unit constant at all theta, so the two channels stay consistent.)"""
    evals = np.asarray(eigenvalues)
    f = fermi_dirac(evals, theta)                       # (4N,)
    U = eigenvectors[:2 * N, :]                         # electron block, ALL states
    U = U.reshape(2, N, -1)                             # (spin, site, 4N)
    u_i = U
    u_ip1 = np.roll(U, -1, axis=1)                      # wrap bond N-1 -> 0
    inner = np.sum(np.conj(u_i) * u_ip1, axis=0)        # (N, 4N)
    phase = np.exp(1j * np.asarray(phi_bonds, dtype=float))[:, None]
    prefactor = 2.0 * E_CHARGE * t_hop / HBAR           # Amperes
    return prefactor * np.sum(f[None, :] * np.imag(inner * phase), axis=1)  # (N,)


def full_current_operators(N: int, phi_bonds: np.ndarray, delta: float = 1.0):
    """The FULL current operator Ĵ = ∂H/∂φ (kinetic + spin-orbit) and its
    diamagnetic partner Ĵ_dia = ∂²H/∂φ², built analytically for a uniform bond
    twist and verified against numerical dH/dφ, d²H/dφ² of `_build_ring_lo_phased`.

    Why the finite-T Kubo needs the FULL current (not R8-B's kinetic-only Ĵ):
    R8-B deliberately used the kinetic current so its T=0 DC gate would bug-catch
    the kinetic self-field response. But at finite T the thermal-REOCCUPATION
    term in d<Ĵ>_T/dφ samples dE_n/dφ = <n|∂H/∂φ|n> — the FULL current — because
    the levels shift under the full Hamiltonian. A kinetic-only Ĵ therefore
    CANNOT satisfy the finite-T omega->0 gate (verified: the kinetic gate fails
    ~17% at theta=0.3, the full-current gate holds to ~1e-6 at all theta). So the
    finite-T Lehmann sum uses Ĵ = ∂H/∂φ; the residual kinetic-vs-full difference
    is exactly R8-B's quantified SO band. (The SELF-FIELD channel — boundary,
    pole, slope s — keeps R5/R7's kinetic `ring_current_density_T` unchanged; its
    s is finite-differenced from the actual thermal kinetic current, so it is
    exact independently of this operator choice.)"""
    t = LO_PARAMS['t_norm'] * float(delta)
    alpha = LO_PARAMS['alpha_norm'] * float(delta)
    I2 = np.eye(2, dtype=complex)
    sy = np.array([[0, -1j], [1j, 0]], dtype=complex)
    T = r8b._T_matrix(N, phi_bonds)
    J = r8b.current_operator_full(N, phi_bonds, delta)          # Ĵ_kin + Ĵ_SO
    # Ĵ_dia = ∂²H/∂φ²: kinetic t(I2⊗(T+T†)) + SO -(α/2)(σ_y⊗ i(T−T†))
    Jd_e = t * np.kron(I2, T + T.conj().T) \
        - (alpha / 2.0) * np.kron(sy, 1j * (T - T.conj().T))
    Z = np.zeros_like(Jd_e)
    Jdia = np.block([[Jd_e, Z], [Z, -Jd_e.conj()]])
    return J, Jdia


def current_expectation_T(E: np.ndarray, U: np.ndarray, N: int,
                          phi_bonds: np.ndarray, theta: float,
                          delta: float = 1.0, full: bool = True) -> float:
    """Operator-form thermal current <Ĵ>_T = Σ_n f(E_n) <n|Ĵ|n>. `full=True`
    uses Ĵ = ∂H/∂φ (kinetic+SO), the operator the finite-T DC gate is defined
    against; `full=False` uses R8-B's kinetic Ĵ (for the T=0 reduction check,
    where at theta->0 it equals Σ_{E<0}<n|Ĵ|n> == R8-B's `_current_expectation`).
    Using PH antisymmetry <n̄|Ĵ|n̄> = −<n|Ĵ|n>, this equals the standard BdG
    thermal current −Σ_{ε>0} tanh(ε/2θ) <ε|Ĵ|ε>."""
    if full:
        J, _ = full_current_operators(N, phi_bonds, delta)
    else:
        J, _ = r8b.current_operators(N, phi_bonds, delta)
    diagJ = np.einsum('ij,ji->i', U.conj().T @ J, U).real
    f = fermi_dirac(E, theta)
    return float(np.sum(f * diagJ))


# ════════════════════════════════════════════════════════════════════════
# TASK 2.1b — finite-T self-consistency map and operating point
# ════════════════════════════════════════════════════════════════════════
def scalar_ring_map_T(N: int, mu: float, B: float, g: float,
                      phi_ext_bond: float, phi_uniform: float, theta: float,
                      delta: float = 1.0, t_hop: float = T_HOP_SI,
                      lattice_const: float = A0,
                      wire_radius: float | None = None) -> float:
    """Finite-T one-shot self-consistency map (thermal generalization of R5
    `scalar_ring_map`). Identical except the ring current is `ring_current_density_T`
    at temperature theta. At theta -> 0 reproduces `scalar_ring_map` bit-for-bit,
    so the whole R5/R7 fixed-point machinery below reduces to the T=0 solver."""
    phi_tot = np.full(N, float(phi_ext_bond) + float(phi_uniform))
    H = _build_ring_lo_phased(N, mu, B, phi_tot, delta)
    evals, evecs = eigh(H)
    I_bonds = ring_current_density_T(evecs, evals, N, t_hop, phi_tot, theta)
    A_t = ring_biot_savart_self(I_bonds, N, lattice_const, wire_radius)
    out = g * compute_peierls_from_A(A_t, lattice_const)
    return float(np.mean(out))


def _find_fixed_points_T(N, mu, B, g, phi_ext, a, theta, span=None,
                         n_scan=49, tol=1e-10, delta=1.0):
    """Finite-T port of R7 `_find_fixed_points_a`: all roots of G_T(phi)-phi on
    the uniform manifold with their local slope dG_T/dphi. (theta enters only
    through `scalar_ring_map_T`.)"""
    if span is None:
        span = 1.5 * 2.0 * np.pi / N

    def f(p):
        return scalar_ring_map_T(N, mu, B, g, phi_ext, p, theta,
                                 lattice_const=a) - p

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
            slope = (scalar_ring_map_T(N, mu, B, g, phi_ext, m + eps, theta,
                                       lattice_const=a)
                     - scalar_ring_map_T(N, mu, B, g, phi_ext, m - eps, theta,
                                         lattice_const=a)) / (2 * eps)
            roots.append({"phi": m, "slope": float(slope),
                          "picard_stable": bool(abs(slope) < 1.0)})
    roots.sort(key=lambda r: abs(r["phi"]))
    return roots


def _phi_star_T(N, mu, B, g, phi_ext, a, theta, hint=None, tol=1e-9, delta=1.0):
    """Finite-T port of R7 `_phi_star_a`: principal (smallest-|phi|) self-
    consistent self-phase and the loop slope s at temperature theta, with the
    same warm-hint bracket + widen-to-3-periods fallback. Returns (phi, slope)
    or (nan, nan) in the dead-zone. At theta -> 0 reproduces `_phi_star_a`."""
    period = 2.0 * np.pi / N

    def f(p):
        return scalar_ring_map_T(N, mu, B, g, phi_ext, p, theta,
                                 lattice_const=a) - p

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
                slope = (scalar_ring_map_T(N, mu, B, g, phi_ext, m + eps, theta,
                                           lattice_const=a)
                         - scalar_ring_map_T(N, mu, B, g, phi_ext, m - eps, theta,
                                             lattice_const=a)) / (2 * eps)
                return m, float(slope)
    roots = _find_fixed_points_T(N, mu, B, g, phi_ext, a, theta, tol=tol, delta=delta)
    if not roots:
        roots = _find_fixed_points_T(N, mu, B, g, phi_ext, a, theta,
                                     span=3.0 * period, n_scan=121,
                                     tol=tol, delta=delta)
    if not roots:
        return float("nan"), float("nan")
    return roots[0]["phi"], roots[0]["slope"]


def operating_point_T(N, mu, B, a, theta, g=G_R7, hint=None, delta=1.0):
    """Finite-T principal self-consistent (phi*, s) — the thermal analogue of
    R8-A `operating_point`. (nan, nan) in the dead-zone."""
    phi_ext = _phi_ext_bond(N)
    p, s = _phi_star_T(N, mu, B, g, phi_ext, a, theta, hint=hint)
    return p, s, (p if np.isfinite(p) else hint)


# ════════════════════════════════════════════════════════════════════════
# TASK 2.1c — finite-T Kubo current-current correlator
# ════════════════════════════════════════════════════════════════════════
def _fermi_deriv(E: np.ndarray, theta: float) -> np.ndarray:
    """df/dE = -(1/theta) f (1-f)  (<= 0). Zero at theta=0 (no thermal Drude)."""
    if theta <= 0.0:
        return np.zeros_like(np.asarray(E, float))
    f = fermi_dirac(E, theta)
    return -(1.0 / float(theta)) * f * (1.0 - f)


def kubo_chi_T(E: np.ndarray, U: np.ndarray, J: np.ndarray, Jdia: np.ndarray,
               w_grid: np.ndarray, eta: float, theta: float,
               tol_deg: float = 1e-9) -> dict:
    """Finite-T retarded current-current correlator (thermal generalization of
    R8-B `kubo_chi`):

        chi(w) = chi_dia
               + Σ_{m≠n} (f_n - f_m) |<m|Ĵ|n>|² / (E_n - E_m + w + i eta)   [inter/intra-band]
               + D · i eta / (w + i eta)                                    [thermal Drude]

    with THERMAL occupations f_n = f(E_n) (Fermi-Dirac), thermal diamagnetic
    term chi_dia = Σ_n f_n <n|Ĵ_dia|n>, and the INTRABAND DRUDE weight

        D = Σ_{|E_n-E_m|<tol_deg} (df/dE)_{E_n} |<m|Ĵ|n>|²   (<= 0),

    the sum over (near-)degenerate states — the m=n diagonal plus any exactly-
    degenerate multiplet. D is the piece of d<Ĵ>_T/dphi that comes from thermal
    RE-OCCUPATION as the levels shift (Feynman-Hellmann: dE_n/dphi = J_nn), which
    the bare Lehmann sum omits; adding it is what makes the finite-T DC gate
    (`check_dc_gate_T`) pass at all theta. It reduces to the standard current-
    current correlator at theta=0.

    Reduction to T=0 (R8-B): at theta -> 0, f_n = (E_n<0) so (f_n-f_m) is ±1 only
    for occ<->unocc pairs (R8-B's `|df|>0.5` mask), chi_dia -> Σ_{E<0}<n|Ĵ_dia|n>,
    and D = 0 (df/dE = 0 for every gapped state). So kubo_chi_T reduces to
    `kubo_chi` bit-for-bit at theta=0. At theta>0 the continuous (f_n-f_m) weights
    open NEW channels (sub-gap absorption from thermally-excited quasiparticles)
    and D turns on a Drude peak of width eta at w=0.

    Drude lineshape caveat (same as R8-B Task 3): the model has no intrinsic
    scattering, so the Drude width is set by the numerical eta, not a physical
    linewidth — its w=0 WEIGHT D is physical, its width is not."""
    f = fermi_dirac(E, theta)                            # (4N,)
    fp = _fermi_deriv(E, theta)                          # df/dE (4N,)
    M = U.conj().T @ J @ U
    W = np.abs(M) ** 2
    df = f[None, :] - f[:, None]                         # f_n - f_m  ([m,n])
    dE = E[None, :] - E[:, None]                         # E_n - E_m  ([m,n])
    Mdia = np.einsum('ij,ji->i', U.conj().T @ Jdia, U).real
    chi_dia = float(np.sum(f * Mdia))                    # thermal diamagnetic

    deg = np.abs(dE) <= tol_deg                          # (near-)degenerate block
    off = ~deg                                           # genuine E_n != E_m
    wv = (df * W)[off]
    dev = dE[off]
    # intraband Drude weight over the degenerate block (m=n and any multiplet)
    fp_avg = 0.5 * (fp[None, :] + fp[:, None])
    D = float(np.sum((fp_avg * W)[deg]))                 # <= 0

    chi = np.empty(w_grid.shape, dtype=complex)
    for k, w in enumerate(w_grid):
        drude = D * (1j * eta) / (w + 1j * eta)
        chi[k] = chi_dia + np.sum(wv / (dev + w + 1j * eta)) + drude
    chi0 = chi_dia + float(np.sum(wv / dev)) + D         # w->0: drude -> D
    return {"w": w_grid, "chi": chi, "chi_dia": chi_dia, "chi0": chi0,
            "drude_weight": D}


def impedance_spectrum_T(N, mu, B, a, w_grid, eta, theta, g=G_R7,
                         omega_unit=None, delta=1.0) -> dict:
    """Finite-T complex AC impedance Z(w) (thermal analogue of R8-B
    `impedance_spectrum`): s(w) = s_cached(theta)·chi_T(w)/chi_T(0),
    L_eff(w)=L_geo/(1-s(w)), Z(w)=i w L_eff. The static self-response slope
    s_cached is taken at the SAME temperature (finite-T operating point), so the
    whole spectrum is thermally consistent. Re Z = dissipation."""
    _RING_KERNEL_CACHE.clear()
    p, s_cached, _ = operating_point_T(N, mu, B, a, theta, g=g, delta=delta)
    phi_ext = _phi_ext_bond(N)
    phi_bonds = np.full(N, phi_ext + (0.0 if not np.isfinite(p) else p))
    H = _build_ring_lo_phased(N, mu, B, phi_bonds, delta)
    E, Uvec = eigh(H)
    L_geo = r8a.geometric_inductance(N, a)
    J, Jdia = full_current_operators(N, phi_bonds, delta)   # Ĵ=∂H/∂φ (see full_current_operators)
    ku = kubo_chi_T(E, Uvec, J, Jdia, w_grid, eta, theta)
    chi, chi0 = ku["chi"], ku["chi0"]
    s_w = s_cached * chi / chi0
    L_eff = L_geo / (1.0 - s_w)
    wphys = w_grid * (omega_unit if omega_unit else 1.0)
    Z = 1j * wphys * L_eff
    return {"w": w_grid, "chi": chi, "chi0": chi0, "chi_dia": ku["chi_dia"],
            "drude_weight": ku["drude_weight"],
            "s_w": s_w, "s_cached": float(s_cached) if np.isfinite(s_cached) else float("nan"),
            "L_geo_H": float(L_geo), "L_eff": L_eff, "Z": Z, "eta": float(eta),
            "theta": float(theta), "N": N, "mu": mu, "B": float(B), "a_nm": a * 1e9}


def check_dc_gate_T(N, mu, B, a, theta, g=G_R7, eps=1e-4, eta=1e-3,
                    delta=1.0) -> dict:
    """Finite-T DC-limit gate (thermal analogue of R8-B `check_dc_gate`): the
    static Kubo value chi_T(0) must equal the central finite difference of the
    THERMAL current <Ĵ>_T w.r.t. a uniform bond-phase twist, at the SAME theta.
    A finite-T Kubo bug breaks this. Verified at multiple theta (Task 2.1c)."""
    _RING_KERNEL_CACHE.clear()
    p, _, _ = operating_point_T(N, mu, B, a, theta, g=g, delta=delta)
    phi_ext = _phi_ext_bond(N)
    phi0 = np.full(N, phi_ext + (0.0 if not np.isfinite(p) else p))

    def Jexp_at(dphi):
        pb = phi0 + dphi
        H = _build_ring_lo_phased(N, mu, B, pb, delta)
        E, Uv = eigh(H)
        return current_expectation_T(E, Uv, N, pb, theta, delta, full=True)

    dJ_findiff = (Jexp_at(+eps) - Jexp_at(-eps)) / (2 * eps)
    H = _build_ring_lo_phased(N, mu, B, phi0, delta)
    E, Uv = eigh(H)
    J, Jdia = full_current_operators(N, phi0, delta)        # Ĵ=∂H/∂φ (full)
    ku = kubo_chi_T(E, Uv, J, Jdia, np.array([0.0]), eta, theta)
    chi0 = ku["chi0"]
    rel = abs(chi0 - dJ_findiff) / max(abs(dJ_findiff), 1e-30)
    return {"theta": float(theta), "dJ_dphi_findiff": dJ_findiff,
            "chi0_kubo": chi0, "rel_mismatch": rel, "passes": bool(rel < 1e-3)}


# ════════════════════════════════════════════════════════════════════════
# TASK 2.1 (boundary) — finite-T self-consistent Pfaffian B_c + soft crossover
# ════════════════════════════════════════════════════════════════════════
def ring_Bc_selfconsistent_T(N, mu, g, phi_ext_bond, a, theta,
                             delta=1.0, window=0.05, n_coarse=11,
                             tol=1e-6, root_tol=1e-10, diagnostics=None):
    """Finite-T self-consistent Pfaffian B_c: the EXACT R7/R5 warm-hint
    bisection (`_ring_Bc_sc_a_singlepath`) but with the finite-T self-phase
    `_phi_star_T`. At each probed B the principal self-phase phi*(B, theta) is
    continued from the previous B and the T=0 Pfaffian invariant is evaluated at
    phi_ext + phi*(B, theta). Temperature enters ONLY through phi*(theta) — the
    Z2 sign is a Hamiltonian property (see module header). At theta -> 0 this
    reproduces R7 `ring_Bc_selfconsistent_a` (hence R5) bit-for-bit.

    Raises ScfConvergenceError-style ValueError only if the principal branch
    detaches on every jittered path (a genuine dead-zone), same as R7."""
    ana = np.sqrt(mu ** 2 + delta ** 2)
    attempts = [(window, n_coarse), (window * 1.11, n_coarse),
                (window * 0.9, n_coarse + 2), (window * 1.23, n_coarse + 4),
                (window * 0.78, n_coarse + 6)]
    last = None
    for k, (w, nc) in enumerate(attempts):
        try:
            bc = _ring_Bc_sc_T_singlepath(N, mu, g, phi_ext_bond, a, theta,
                                          delta, w, nc, tol, root_tol, diagnostics)
            if np.isfinite(bc):
                if diagnostics is not None:
                    diagnostics["solve_retries"] = k
                return bc
        except _DeadZone as e:
            last = e
    raise _DeadZone(f"genuine dead-zone at N={N}, mu={mu}, a={a*1e9:.3f}nm, "
                    f"theta={theta:.4g}, g={g:.3g} ({last})")


class _DeadZone(RuntimeError):
    """Local dead-zone signal (mirrors R5 ScfConvergenceError; kept local so we
    never import a name that could couple to an R1-R8 code path)."""


def _ring_Bc_sc_T_singlepath(N, mu, g, phi_ext_bond, a, theta, delta,
                             window, n_coarse, tol, root_tol, diagnostics):
    """One warm-hint finite-T self-consistent B_c bisection (port of R7
    `_ring_Bc_sc_a_singlepath` with `_phi_star_T`)."""
    ana = np.sqrt(mu ** 2 + delta ** 2)
    mem = {"hint": None}
    slopes_offtrans = []

    def s(B):
        p, slope = _phi_star_T(N, mu, B, g, phi_ext_bond, a, theta,
                               hint=mem["hint"], tol=root_tol)
        if not np.isfinite(p):
            raise _DeadZone(f"no fixed point at B={B:.6f}")
        mem["hint"] = p
        phi_bonds = np.full(N, phi_ext_bond + p)
        val = _ring_pfaffian_sign(N, mu, B, phi_bonds, delta)
        if val == 0.0:
            val = _ring_pfaffian_sign(N, mu, B + 1e-9, phi_bonds, delta)
        return val, p, slope

    h = window
    for _ in range(5):
        Bg = np.linspace(max(0.02, ana - h), ana + h, n_coarse)
        out = [s(B) for B in Bg]
        signs = np.array([o[0] for o in out])
        slopes_offtrans = [out[0][2], out[-1][2]]
        flip = np.where(signs[:-1] != signs[1:])[0]
        if flip.size == 0:
            h *= 2.0
            mem["hint"] = None
            continue
        i = int(flip[0])
        B_lo, B_hi, s_lo = Bg[i], Bg[i + 1], signs[i]
        mem["hint"] = out[i][1]
        p_at = out[i][1]
        while B_hi - B_lo > tol:
            Bm = 0.5 * (B_lo + B_hi)
            val, p_at, _ = s(Bm)
            if val == s_lo:
                B_lo = Bm
            else:
                B_hi = Bm
        Bc = 0.5 * (B_lo + B_hi)
        if diagnostics is not None:
            diagnostics.update(phi_star_at_Bc=float(p_at),
                               slope_offtransition=[float(x) for x in slopes_offtrans])
        return Bc
    return float("nan")


def clean_ring_gap(N, mu, B, phi_bonds, delta=1.0) -> float:
    """Bulk gap min|E| of the ring at given bond phases (one eigh)."""
    H = _build_ring_lo_phased(int(N), float(mu), float(B),
                              np.asarray(phi_bonds, float), float(delta))
    return float(np.min(np.abs(eigh(H, eigvals_only=True))))


def thermal_parity_indicator(N, mu, Bgrid, theta, phi_bonds_fn=None,
                             delta=1.0) -> dict:
    """SOFT finite-T boundary (definition B in the header). For each B, form

        P(B) = sgn_Pf(B) * tanh( E_gap(B) / (2 theta) ),

    the thermal expectation of the ground-state fermion parity given a two-level
    (±E_gap) gap at temperature theta. |P| -> 1 deep in either phase and -> 0 at
    the gap-closing over a B-window of width ~theta: the SHARP T=0 Z2 becomes a
    finite-T crossover. Returns P(B), the zero-crossing (soft boundary), and the
    crossover width (B-range where |P| < 0.5). phi_bonds_fn(B) supplies the
    (optionally self-consistent) texture; default = clean ring."""
    Bgrid = np.asarray(Bgrid, float)
    P = np.empty_like(Bgrid)
    gap = np.empty_like(Bgrid)
    for i, B in enumerate(Bgrid):
        pb = (phi_bonds_fn(B) if phi_bonds_fn is not None
              else np.zeros(N))
        sgn = _ring_pfaffian_sign(N, mu, B, pb, delta)
        if sgn == 0.0:
            sgn = _ring_pfaffian_sign(N, mu, B + 1e-9, pb, delta)
        g_ = clean_ring_gap(N, mu, B, pb, delta)
        gap[i] = g_
        P[i] = sgn * (np.tanh(g_ / (2.0 * theta)) if theta > 0 else 1.0)
    # soft boundary = sign change of P (== Pfaffian flip); width = |P|<0.5 band
    soft = float("nan")
    fl = np.where(np.sign(P[:-1]) != np.sign(P[1:]))[0]
    if fl.size:
        i = int(fl[0])
        soft = float(0.5 * (Bgrid[i] + Bgrid[i + 1]))
    narrow = Bgrid[np.abs(P) < 0.5]
    width = float(narrow.max() - narrow.min()) if narrow.size >= 2 else 0.0
    return {"B": Bgrid, "P": P, "gap": gap, "soft_boundary": soft,
            "crossover_width": width, "theta": float(theta)}


# ════════════════════════════════════════════════════════════════════════
# TASK 2.2 — validation against the T -> 0 limit (do BEFORE trusting T>0)
# ════════════════════════════════════════════════════════════════════════
def validate_T0_limit(N=80, mu=1.2, a=A0, g=G_R7, verbose=True) -> dict:
    """Confirm every finite-T quantity reduces to the existing T=0 R5/R7/R8
    result as theta -> 0, AND that the two current channels stay proportional at
    finite theta. Returns a dict of max abs/rel differences; all must be tiny."""
    _RING_KERNEL_CACHE.clear()
    phi_ext = _phi_ext_bond(N)
    B = np.sqrt(mu ** 2 + 1.0) - 0.05          # near the transition, off the closing
    res = {}

    # (1) ring current: theta=0 vs R5 ring_current_density (bit-for-bit)
    phi_tot = np.full(N, phi_ext + 0.001)
    H = _build_ring_lo_phased(N, mu, B, phi_tot, 1.0)
    E, U = eigh(H)
    I_T0 = ring_current_density_T(U, E, N, T_HOP_SI, phi_tot, 0.0)
    I_R5 = ring_current_density(U, E, N, T_HOP_SI, phi_tot)
    res["current_T0_vs_R5_maxabs"] = float(np.max(np.abs(I_T0 - I_R5)))

    # (1b) operator current == C * ring current at all theta (unit constant)
    #      C fixed at theta=0; check ratio is theta-independent.
    consts = []
    for th in (0.0, 0.05, 0.2):
        Jexp = current_expectation_T(E, U, N, phi_tot, th, full=False)  # kinetic vs kinetic
        Iring = float(np.sum(ring_current_density_T(U, E, N, T_HOP_SI, phi_tot, th)))
        consts.append(Jexp / Iring if Iring != 0 else float("nan"))
    res["opcurrent_ratio"] = float(consts[0])
    res["opcurrent_ratio_theta_relspread"] = float(
        (np.nanmax(consts) - np.nanmin(consts)) / abs(consts[0]))

    # (2) scalar map: theta=0 vs R5 scalar_ring_map (bit-for-bit)
    from bdg_r5_selfconsistent_ring import scalar_ring_map
    m_T0 = scalar_ring_map_T(N, mu, B, g, phi_ext, 0.002, 0.0, lattice_const=a)
    m_R5 = scalar_ring_map(N, mu, B, g, phi_ext, 0.002, lattice_const=a)
    res["scalarmap_T0_vs_R5_absdiff"] = float(abs(m_T0 - m_R5))

    # (3) operating point: theta=0 vs R8-A operating_point (phi*, s)
    p_T, s_T, _ = operating_point_T(N, mu, B, a, 0.0, g=g)
    p_A, s_A, _ = r8a.operating_point(N, mu, B, a, g=g)
    res["phistar_T0_vs_R8A_absdiff"] = float(abs(p_T - p_A))
    res["slope_T0_vs_R8A_absdiff"] = float(abs(s_T - s_A))

    # (4) Kubo chi: theta=0 vs R8-B kubo_chi over a w-grid
    J, Jdia = r8b.current_operators(N, phi_tot, 1.0)
    w = np.linspace(0.0, 3.0, 25)
    kT = kubo_chi_T(E, U, J, Jdia, w, 0.02, 0.0)
    kB = r8b.kubo_chi(E, U, J, Jdia, w, 0.02)
    res["kubo_chi_T0_vs_R8B_maxabs"] = float(np.max(np.abs(kT["chi"] - kB["chi"])))
    res["kubo_chi0_T0_vs_R8B_absdiff"] = float(abs(kT["chi0"] - kB["chi0"]))

    # (5) DC gate at multiple theta (Task 2.1c: not just theta=0). Uses the full
    # current operator; holds across the whole physically-meaningful range
    # (theta <= ~0.15; above that the fixed-Delta BdG model is itself unphysical
    # as T -> Tc/2, since the pairing would be thermally suppressed).
    gates = {th: check_dc_gate_T(N, mu, B, a, th, g=g)["rel_mismatch"]
             for th in (0.0, 0.009, 0.05, 0.1, 0.15)}
    res["dc_gate_rel_by_theta"] = {f"{k:.3g}": float(v) for k, v in gates.items()}
    res["dc_gate_max_rel"] = float(max(gates.values()))

    ok = (res["current_T0_vs_R5_maxabs"] < 1e-18
          and res["scalarmap_T0_vs_R5_absdiff"] < 1e-15
          and res["phistar_T0_vs_R8A_absdiff"] < 1e-9
          and res["kubo_chi_T0_vs_R8B_maxabs"] < 1e-9
          and res["opcurrent_ratio_theta_relspread"] < 1e-10
          and res["dc_gate_max_rel"] < 1e-3)
    res["ALL_PASS"] = bool(ok)
    if verbose:
        print(json.dumps(res, indent=2, default=str))
    return res


# ════════════════════════════════════════════════════════════════════════
# TASK 3 — parallel workers (top-level & picklable for Windows spawn)
# ════════════════════════════════════════════════════════════════════════
def _bc_worker(N, mu, a, theta, g):
    """One finite-T self-consistent B_c (+ clean baseline + shift %). Returns a
    JSON-able dict; dead-zone -> shift nan. Used as the parallel grid unit."""
    _RING_KERNEL_CACHE.clear()
    phi_ext = _phi_ext_bond(N)
    ana = np.sqrt(mu ** 2 + 1.0)
    from bdg_r4_em_coupling import ring_boundary_Bc_pfaffian
    base = ring_boundary_Bc_pfaffian(N, mu)          # clean-ring baseline (T-independent)
    diag = {}
    try:
        bc = ring_Bc_selfconsistent_T(N, mu, g, phi_ext, a, theta, diagnostics=diag)
        shift = 100.0 * (bc - base) / base
    except _DeadZone:
        bc, shift = float("nan"), float("nan")
    return {"N": N, "mu": mu, "a_nm": a * 1e9, "theta": float(theta), "g": g,
            "baseline_Bc": float(base), "sc_Bc": float(bc), "shift_pct": float(shift),
            "phi_star_at_Bc": float(diag.get("phi_star_at_Bc", float("nan")))}


def _pole_worker(N, mu, B, a, theta, g):
    """Enhancement 1/(1-s) at one (g, theta) on the principal branch (Track-1
    style). Returns s and enhancement (or nan in the dead-zone)."""
    _RING_KERNEL_CACHE.clear()
    phi_ext = _phi_ext_bond(N)
    p, s = _phi_star_T(N, mu, B, g, phi_ext, a, theta, hint=None)
    if not np.isfinite(s):
        return {"g": g, "theta": float(theta), "s": float("nan"),
                "enh": float("nan")}
    enh = 1.0 / (1.0 - s) if s < 1.0 else float("inf")
    return {"g": g, "theta": float(theta), "s": float(s), "enh": float(enh)}


def _timing_worker(N, mu, B, a, theta, g):
    """Wall-time of a single finite-T operating-point solve, for the load-balance
    check (does per-point cost depend on theta?)."""
    _RING_KERNEL_CACHE.clear()
    phi_ext = _phi_ext_bond(N)
    t0 = time.perf_counter()
    operating_point_T(N, mu, B, a, theta, g=g)
    return {"theta": float(theta), "wall_s": time.perf_counter() - t0}


# ════════════════════════════════════════════════════════════════════════
# TASK 3 — parallel benchmark + BLAS oversubscription + load-balance check
# ════════════════════════════════════════════════════════════════════════
def benchmark_parallel(mu=1.2, a=A0, g=G_R7, thetas=(0.0, 0.05, 0.1, 0.15),
                       N=80) -> dict:
    """Task 3: (1) serial-vs-parallel equivalence on a small (mu,B,theta) subset
    (bit-identical at fixed 1-thread BLAS), (2) BLAS-oversubscription benchmark
    (serial-1-thread vs serial-default-BLAS vs parallel-capped) on a
    representative subset, and (3) per-point cost vs theta (does a wider
    Fermi-Dirac tail make higher-T points cost more, requiring dynamic
    scheduling?)."""
    out = {}
    B = np.sqrt(mu ** 2 + 1.0) - 0.05
    # (1) equivalence: pole-worker over a small (g, theta) subset
    eq_args = [(N, mu, B, a, th, gg)
               for th in (0.0, 0.1) for gg in (2.9e4, 6.0e4, 8.0e4)]
    out["equivalence"] = check_equivalence(_pole_worker, eq_args)

    # (2) BLAS benchmark on the same representative subset (operating-point solves)
    def _time(fn):
        t0 = time.perf_counter(); fn(); return time.perf_counter() - t0
    from bdg_parallel import _pin_blas
    # serial, BLAS pinned to 1 thread
    _pin_blas(1)
    t_serial1 = _time(lambda: [_pole_worker(*ar) for ar in eq_args])
    # serial, default multithreaded BLAS
    _pin_blas(os.cpu_count() or 8)
    t_serialdef = _time(lambda: [_pole_worker(*ar) for ar in eq_args])
    _pin_blas(1)
    # parallel, 1-thread workers
    nw = default_workers()
    t_par = _time(lambda: parallel_map(_pole_worker, eq_args, n_workers=nw, blas_threads=1))
    out["blas"] = {"serial_1thread_s": t_serial1, "serial_defaultBLAS_s": t_serialdef,
                   "parallel_capped_s": t_par, "n_workers": nw,
                   "n_tasks": len(eq_args)}

    # (3) per-point cost vs theta (load balance): time one operating-point solve
    #     at each theta, several reps
    timing = {}
    for th in thetas:
        reps = [_timing_worker(N, mu, B, a, th, g)["wall_s"] for _ in range(3)]
        timing[f"{th:.3g}"] = float(np.median(reps))
    vals = list(timing.values())
    out["per_theta_cost_s"] = timing
    out["cost_spread_ratio"] = float(max(vals) / max(min(vals), 1e-9))
    out["scheduling_note"] = (
        "static equal-chunking OK (cost ~theta-independent)"
        if out["cost_spread_ratio"] < 1.5 else
        "cost varies with theta -> prefer dynamic scheduling (chunksize=1)")
    return out


# ════════════════════════════════════════════════════════════════════════
# TASK 4.1 — does the bulk boundary shift measurably at finite T?
# ════════════════════════════════════════════════════════════════════════
def boundary_shift_vs_T(mu=1.2, a=A0, g=G_R7, N=80,
                        thetas=(0.0, 0.009, 0.05, 0.1, 0.15),
                        n_scan=(80, 100, 120), n_workers=None) -> dict:
    """Self-consistent Pfaffian boundary shift vs temperature. Compares the shift
    at each theta to the T=0 null (R4-R7), then applies the project-standard
    N-independence check at the operating theta and at the largest theta. T
    enters only through the finite-T self-field phi*(theta)."""
    B = None
    # (a) shift vs theta at fixed N
    args = [(N, mu, a, th, g) for th in thetas]
    rows = parallel_map(_bc_worker, args, n_workers=n_workers)
    # (b) N-independence at operating theta and largest theta
    nscan = {}
    for th in (thetas[1] if len(thetas) > 1 else thetas[0], thetas[-1]):
        nargs = [(Nn, mu, a, th, g) for Nn in n_scan]
        nscan[f"{th:.3g}"] = parallel_map(_bc_worker, nargs, n_workers=n_workers)
    return {"mu": mu, "a_nm": a * 1e9, "g": g, "N": N,
            "shift_vs_theta": rows, "nscan": nscan}


# ════════════════════════════════════════════════════════════════════════
# TASK 4.2 — does the g_crit pole survive / sharpen / wash out at finite T?
# ════════════════════════════════════════════════════════════════════════
def pole_vs_T(mu=1.2, a=A0, N=80, thetas=(0.0, 0.05, 0.1, 0.15),
              gs=None, B_offset=0.05, n_workers=None) -> dict:
    """Enhancement 1/(1-s) on the principal branch vs feedback g, at several
    theta, near R5's g_crit band (reuses R9 Track-1's leading-edge protocol).
    Thermal broadening reduces the current response -> reduces s -> should push
    the detachment/pole to HIGHER g (wash out) or leave it (survive)."""
    if gs is None:
        gs = list(np.linspace(1.0e4, 1.2e5, 16))
    B = np.sqrt(mu ** 2 + 1.0) - B_offset
    series = {}
    for th in thetas:
        args = [(N, mu, B, a, th, gg) for gg in gs]
        res = parallel_map(_pole_worker, args, n_workers=n_workers)
        # peak enhancement before the principal branch detaches (s finite, <1)
        enh = [r["enh"] for r in res]
        finite = [e for e in enh if np.isfinite(e)]
        series[f"{th:.3g}"] = {"gs": [float(x) for x in gs], "res": res,
                               "max_enh": float(max(finite)) if finite else float("nan")}
    return {"mu": mu, "a_nm": a * 1e9, "N": N, "B": float(B), "series": series}


# ════════════════════════════════════════════════════════════════════════
# TASK 4.3 — does the dissipation spectrum change qualitatively at finite T?
# ════════════════════════════════════════════════════════════════════════
def dissipation_vs_T(mu=1.2, a=A0, g=G_R7, N=80, B_offset=0.05,
                     thetas=(0.0, 0.05, 0.1, 0.15), eta=0.02,
                     wmax=4.0, nw=400) -> dict:
    """Kubo dissipation Re Z(w) / Re chi(w) across a frequency grid at several
    theta at one clean operating point. Looks for (i) NEW sub-gap absorption
    from thermally-excited quasiparticles (absent at T=0), and (ii) thermal
    smearing of the interband onset feature (R8-B T2). Uses the full-current
    Kubo (finite-T DC-gate-validated)."""
    B = np.sqrt(mu ** 2 + 1.0) - B_offset
    w = np.linspace(1e-3, wmax, nw)
    spectra = {}
    for th in thetas:
        sp = impedance_spectrum_T(N, mu, B, a, w, eta, th, g=g,
                                  omega_unit=r8b.OMEGA_UNIT)
        rez = np.real(sp["Z"])
        # sub-gap dissipation weight (w below the ~1.5*Delta interband feature)
        subgap = w < 1.0
        spectra[f"{th:.3g}"] = {
            "w": [float(x) for x in w],
            "reZ": [float(x) for x in rez],
            "re_chi": [float(x) for x in np.real(sp["chi"])],
            "im_chi": [float(x) for x in np.imag(sp["chi"])],
            "drude_present": bool(sp["drude_weight"] != 0.0),
            "drude_weight": float(sp["drude_weight"]),
            "subgap_reZ_peak": float(np.max(np.abs(rez[subgap]))) if subgap.any() else 0.0,
            "s_cached": sp["s_cached"]}
    return {"mu": mu, "a_nm": a * 1e9, "g": g, "N": N, "B": float(B),
            "eta": eta, "spectra": spectra}


# ════════════════════════════════════════════════════════════════════════
# main / drivers  (Task-4 experiments live in the driver functions above)
# ════════════════════════════════════════════════════════════════════════
def print_task1():
    rep = temperature_scale_report()
    print("\n=== TASK 1: temperature scale ===")
    print(f"Delta = {rep['Delta_eV']*1e3:.3f} meV (Chang 2015, R6-verified)")
    print(f"Delta/k_B = {rep['Delta_over_kB_K']:.4f} K   (k_B*Tc ~ Delta convention)")
    print(f"BCS Tc = Delta/(1.764 k_B) = {rep['Tc_BCS_K']:.4f} K")
    print(f"theta at 20 mK base = {rep['theta_at_base_20mK']:.5f}")
    print(f"{'T(mK)':>7} {'theta':>10} {'T/Tc':>10} {'f(E=Delta)':>14}")
    for r in rep["rows"]:
        print(f"{r['T_mK']:>7} {r['theta']:>10.5f} {r['T_over_Tc_conv']:>10.5f} "
              f"{r['fermi_at_gap']:>14.3e}")
    return rep


def make_figures(data: dict) -> list:
    """Three-panel finite-T summary: boundary shift vs T, pole vs g at several T,
    dissipation Re Z(w) vs w at several T."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    paths = []

    fig, axs = plt.subplots(1, 3, figsize=(16, 4.6))
    # panel 1: boundary shift vs theta
    bs = data["boundary"]["shift_vs_theta"]
    th = [r["theta"] for r in bs]
    sh = [r["shift_pct"] for r in bs]
    axs[0].plot(th, sh, "o-")
    axs[0].axhline(0, color="k", lw=0.6)
    axs[0].axvline(theta_of_T(0.020), color="tab:red", ls="--", lw=1,
                   label="20 mK (operating)")
    axs[0].set_xlabel(r"$\theta = k_B T/\Delta$"); axs[0].set_ylabel(r"$\Delta B_c/B_c$ (%)")
    axs[0].set_title("Task 4.1: bulk boundary shift vs T"); axs[0].legend(fontsize=8)

    # panel 2: pole enhancement vs g at several theta
    for k, s in data["pole"]["series"].items():
        gs = s["gs"]; enh = [r["enh"] for r in s["res"]]
        axs[1].plot(gs, enh, "o-", ms=3, label=fr"$\theta$={k}")
    axs[1].set_xlabel("g (feedback)"); axs[1].set_ylabel(r"$1/(1-s)$")
    axs[1].set_title("Task 4.2: g_crit pole vs T"); axs[1].legend(fontsize=8)

    # panel 3: dissipation |Re Z(w)| on log scale so the sub-gap (w<1) thermal
    # absorption — 4 orders below the interband peaks — is visible.
    for k, s in data["dissipation"]["spectra"].items():
        rez = np.abs(np.asarray(s["reZ"])) + 1e-12
        axs[2].semilogy(s["w"], rez, lw=1.2, label=fr"$\theta$={k}")
    axs[2].axvline(2.0, color="gray", ls=":", lw=1, label=r"pair-breaking $2\Delta$")
    axs[2].axvspan(0, 1.0, color="gray", alpha=0.08)
    axs[2].set_xlabel(r"$\hbar\omega/\Delta$"); axs[2].set_ylabel(r"|Re $Z$| ($\Omega$)")
    axs[2].set_title("Task 4.3: dissipation vs T (shaded = sub-gap)")
    axs[2].legend(fontsize=8)

    fig.tight_layout()
    p = _out("r9_finiteT_summary.png"); fig.savefig(p, dpi=130); plt.close(fig)
    paths.append(p)
    return paths


def write_summary(data: dict) -> str:
    """Write the machine-readable JSON + a short human summary of the run."""
    jpath = _out("r9_finiteT_data.json")
    with open(jpath, "w") as fh:
        json.dump(data, fh, indent=2, default=lambda o: float(o)
                  if isinstance(o, (np.floating,)) else str(o))
    return jpath


def run_full(N=80, n_workers=None, quick=False) -> dict:
    """Full R9-Track-2 run: Task 1 scale, T->0 validation, Task 3 benchmark,
    Task 4.1/4.2/4.3 experiments. Parallelised across the (mu,B,theta,g) grids."""
    t_all = time.time()
    thetas_bulk = (0.0, 0.009, 0.05, 0.1, 0.15)
    thetas_ac = (0.0, 0.05, 0.1, 0.15)
    if quick:
        N = 60
        thetas_bulk = (0.0, 0.1)
        thetas_ac = (0.0, 0.1)

    data = {"task1_scale": temperature_scale_report()}

    print("\n=== T->0 VALIDATION ===")
    data["validation"] = validate_T0_limit(N=max(N, 80), mu=1.2, verbose=True)

    print("\n=== TASK 3: parallel benchmark ===")
    data["benchmark"] = benchmark_parallel(N=N, thetas=thetas_ac)
    print(json.dumps(data["benchmark"], indent=2, default=str))

    print("\n=== TASK 4.1: boundary shift vs T ===")
    data["boundary"] = boundary_shift_vs_T(
        N=N, thetas=thetas_bulk,
        n_scan=(80, 100, 120) if not quick else (60, 80),
        n_workers=n_workers)
    for r in data["boundary"]["shift_vs_theta"]:
        print(f"  theta={r['theta']:.4g}: sc_Bc={r['sc_Bc']:.5f} "
              f"shift={r['shift_pct']:+.5f}%")

    print("\n=== TASK 4.2: pole vs T ===")
    data["pole"] = pole_vs_T(N=N, thetas=thetas_ac,
                             gs=(list(np.linspace(1e4, 1.2e5, 16)) if not quick
                                 else list(np.linspace(1e4, 9e4, 5))),
                             n_workers=n_workers)
    for k, s in data["pole"]["series"].items():
        print(f"  theta={k}: max_enh={s['max_enh']:.3f}")

    print("\n=== TASK 4.3: dissipation vs T ===")
    data["dissipation"] = dissipation_vs_T(N=N, thetas=thetas_ac,
                                           nw=400 if not quick else 120)
    for k, s in data["dissipation"]["spectra"].items():
        print(f"  theta={k}: subgap_reZ_peak={s['subgap_reZ_peak']:.3e} "
              f"drude={s['drude_present']}")

    data["wall_time_s"] = time.time() - t_all
    jpath = write_summary(data)
    figs = make_figures(data)
    print(f"\nWrote {jpath}\n      {figs}\nTotal wall: {data['wall_time_s']:.1f}s")
    return data


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true",
                    help="fast smoke test (Task 1 + validation, small N)")
    ap.add_argument("--validate", action="store_true", help="only T->0 validation")
    ap.add_argument("--full", action="store_true", help="full Task 3+4 run")
    ap.add_argument("--quickfull", action="store_true", help="full run, tiny grid")
    ap.add_argument("--workers", type=int, default=None)
    args = ap.parse_args()

    print_task1()
    if args.validate or args.quick:
        print("\n=== T->0 VALIDATION (N=80) ===")
        validate_T0_limit(N=80, mu=1.2)
    if args.full or args.quickfull:
        run_full(n_workers=args.workers, quick=args.quickfull)
