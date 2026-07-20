# ISEF-2027-bdg-eigensolver
Symmetry-aware reduced-order BdG eigensolver for topological phase transitions in superconducting materials

# ISEF 2027 — Symmetry-Aware BdG Eigensolver

**A Symmetry-Aware Reduced-Order Eigensolver for Tracking Topological 
Phase Transitions in Superconducting Nanowires Under Inductive 
Electromagnetic Feedback**

## Project Overview
This project develops three contributions toward fault-tolerant 
quantum computing:
1. A PH-symmetric rank-2 subspace eigensolver for BdG Hamiltonians
2. A real-time Pfaffian topological phase indicator
3. A self-consistent quasi-static EM feedback loop (Biot-Savart)

## Current Status
- ✅ R1: Baseline BdG eigensolver (validated vs Kitaev 2001)
- ✅ R2: Pfaffian topological indicator (finite-size scaling confirmed)
- ✅ R3: PH-symmetric rank-2 update
- ✅ R4: Biot-Savart EM coupling layer — pipeline complete and validated
  (6/6 checks: uniform-A Peierls phases, Biot-Savart vs. analytic, PH
  symmetry, zero-flux invariance, current-density units/sign, phase-boundary
  shift under flux). R4's own open-wire spectral-threshold detector
  over-estimates the boundary-shift magnitude ~10×, so a rigorous,
  threshold-free, N-independent ring/Pfaffian Z₂ detector was added within
  R4 itself (`ring_boundary_Bc_pfaffian`, independently validated); with the
  *same* self-field it gives ~0.6–0.8 % at μ = 1.2–1.4 (vs ~6–9 % from the
  threshold detector). Deriving the final, physically-grounded headline
  number was deliberately handed off to R5 — and that number has since been
  superseded three more times (R5 → R5-revised → R7). See
  `src/r4_pfaffian_boundary_methods.md` §1–4b.
- ✅ R5: Self-consistent ring self-field (`src/bdg_r5_selfconsistent_ring.py`;
  R4 imported read-only, all 6 R4 checks still pass). Replaced *both* interim
  placeholders — the open-wire-seeded texture and the gauge-only wrap bond —
  with a real ring calculation: periodic ring current (incl. the wrap bond),
  a Biot–Savart kernel adapted to the closed loop and validated against the
  analytic loop inductance, and a fully self-consistent ring/Pfaffian B_c.
  Found: plain iteration does **not** converge at g=3×10⁵ (map slope ≈3>1,
  Picard-repelling) — solved by root-finding; the fixed point is locally unique
  (Newton–Krylov finds no non-uniform solution) but the map is multi-valued
  (fluxoid branches; ramp = direct-jump, path-independent). **Important:** the
  g=3×10⁵ number reported here (+0.05 %) later proved to be **past a dynamical
  bifurcation** — see R5-revised.
- ✅ R5-revised: Physically meaningful headline
  (`src/bdg_r5_physically_meaningful_shift.py`). g = 3×10⁵ sits **above g_crit**,
  where **no attracting fixed point exists** (Figure 8) — so its "+0.05 %" is a
  dynamically unreachable root, not a static equilibrium. Mapping the shift vs g
  **restricted to g < g_crit** (Figure 10; correct ring baseline) gives a small,
  **negative**, monotonic, path-independent dose-response, → 0 as g → 0. Three
  regimes: **regime 1** (g ≲ 3×10⁴) unambiguous boundary; **regime 2**
  (3×10⁴–1×10⁵) a dead-zone opens at the transition; **regime 3** (g ≳ 1×10⁵) no
  attracting fixed point (where g=3×10⁵ lives). Headline (top of regime 1,
  g ≈ 2.9×10⁴, N=80): **ΔB_c/B_c ≈ −0.064 % (μ=1.2), −0.066 % (μ=1.4)**.
  **The decisive result (Figure 11):** even this **vanishes as N → ∞
  (∼N⁻¹·⁸)** — a finite-size Aharonov–Bohm artifact, not a bulk shift.
  **Honest bottom line — a null result:** across the *entire* physically
  accessible coupling range, the inductive self-field does **not** move the
  *bulk* topological boundary (tiny, finite-size, → 0 in the bulk). Full lineage:
  threshold −5.86/−9.36 % → interim −0.77/−0.61 % → g=3×10⁵ +0.05 %
  (unreachable) → **physical −0.064/−0.066 % (finite-size → 0)**. The g-literature
  plausibility range straddles all three regimes and never pinned g below g_crit;
  what happens at g > g_crit (device doesn't reach it, or needs a *dynamical*
  treatment) is an open R6 question. See `src/r4_pfaffian_boundary_methods.md` §4e.
- ✅ R6: Phase diagram under EM feedback (`src/bdg_r6_lo_phase_diagram.py`).
  Lutchyn–Oreg InAs/Al 2D (B, μ) phase diagram with a genuine TRIM Pfaffian
  invariant, reproducing the ORV (arXiv:1003.1145) gap-closing structure. See
  `src/r6_lo_phase_diagram_methods.md`.
- ✅ R7: Circuit inductance / wire length **L = N·a** as a third axis (B, μ, L)
  (`src/bdg_r7_inductance_sweep.py`; imports R4/R5 read-only; all R4 checks and
  the R5 detector still pass — solver validated bit-for-bit vs R5, |ΔB_c| = 0).
  Hypothesis-neutral re-opening at the one validated coupling g = 2.9×10⁴
  (top of R5 regime 1). Two orthogonal scans disentangle L from finite-size:
  a fixed-a=10 nm **L-axis** (varies L and N together) and a **fixed-L mesh
  refinement** (holds the physical length fixed, a = L/N). **Finding: L does
  NOT move the bulk topological boundary.** Where a static boundary exists
  (L ≳ 0.8 μm) the shift is tiny (−0.06 % at L=0.8 μm, μ=1.2/1.4) and **vanishes
  under fixed-L mesh refinement (∼N⁻³·¹¹)** — a discretisation/finite-size
  artifact, a continuation of the R4/R5 null result along the new L axis. The
  **one** genuinely L-dependent effect is not a bulk shift but the
  **dynamical-stability (dead-zone) boundary**: at this fixed g, short wires
  (L ≤ 0.6 μm) lose the clean static self-consistent boundary (regime-2). See
  `src/r4_pfaffian_boundary_methods.md` §4g.
- ✅ R8: **Reactive impedance Z_r as a second observable**, alongside the bulk
  gap, split into two tracks (`src/r8_reactive_impedance_methods.md`).
  - **Track A — static feedback-renormalized reactance**
    (`src/bdg_r8_reactive_impedance.py`): `Z_r = iωL_eff`, `L_eff =
    L_geo/(1−s)`, where `s` is confirmed (bit-exact) to be R5's bifurcation
    slope and `L_geo` the R5-validated ring inductance — reuses the cached
    converged state, **zero new eigendecompositions**. **Verdict:** at
    production g = 2.9×10⁴, Z_r is largely *redundant* with the gap + geometry
    (nearly flat across the transition; its 1/(1−s) enhancement is a finite-size
    effect → 1 under fixed-L mesh refinement). Its one orthogonal feature is
    *dynamical, not topological*: `L_eff` diverges as s→1 (g ≈ 8–9×10⁴, R5's
    g_crit band), which the g-independent gap cannot see.
  - **Track B / Phase 8B — finite-ω Kubo current–current correlator**
    (`src/bdg_r8b_kubo_impedance.py`): `χ_JJ(ω)` (Lehmann sum over the full BdG
    spectrum) → complex `Z(ω) = iωL_eff(ω)` with a resistive part from interband
    absorption (new physics). **DC-limit cross-check PASSES** (`χ_JJ(0)`
    reproduces Track A's static response to ~1×10⁻⁴ at all 6 clean points, incl.
    the near-B_c region) — **this lifts the PRELIMINARY label off Track A**.
    Adversarial follow-ups: the dissipation onset is a bulk-pairing feature that
    does **not** track the topological minigap; the η→0 non-convergence is
    **finite-size** (roughness ∝ N⁻⁰·⁵), not an intrinsic linewidth; the omitted
    spin-orbit current is propagated as an explicit **−2 to −4 %** band; gap
    provenance corrected to the R6-verified **Δ = 0.19 meV** (Chang et al. 2015).
    Net: Z(ω) confirms Track A in the read-out band and adds no trustworthy
    topological information there.
  - **Track C — self-consistent spin-orbit current**
    (`src/bdg_r8c_so_selfconsistent.py`): sources the self-field from the
    **total** current (kinetic + SO) and re-solves φ*, promoting Task 4's linear
    SO band to a self-consistent treatment. Result: L_eff drops **−2.3 to
    −4.0 %** (agreeing with the linear estimate to ~0.01 %, since φ*~1e-3 makes
    the feedback correction second order), and the SO shift is **itself a
    finite-size effect** (shrinks −3.2 %→−1.3 % under mesh refinement, →0 in the
    continuum). Confirms — does not overturn — the R8 verdict.
  - Independent sweeps are process-parallelised (`src/bdg_parallel.py`;
    serial-vs-parallel **bit-identical**).

## Requirements

pip install numpy scipy matplotlib pandas pfapack psutil

## Running

    python src/bdg_eigensolver_r1_r2_LOAdapted.py   # R1+R2
    python src/bdg_eigensolver_r3_lo_adapted.py     # R3
    python src/bdg_r4_em_coupling.py                # R4
    python src/bdg_r5_selfconsistent_ring.py        # R5 (~35 min full run)
    python src/bdg_r5_physically_meaningful_shift.py # R5-revised (g<g_crit headline; ~40 min)
    python src/bdg_r6_lo_phase_diagram.py           # R6 (LO 2D phase diagram)
    python src/bdg_r7_inductance_sweep.py           # R7 (B,mu,L sweep; resumable, eigh-bound)
    python src/bdg_r8_reactive_impedance.py         # R8 Track A (Z_r vs gap; add --regen to recompute)
    python src/bdg_r8b_kubo_impedance.py            # R8-B (finite-w Kubo impedance + tasks 1-5)
    python src/bdg_r8c_so_selfconsistent.py         # R8-C (self-consistent spin-orbit current)

Each program writes its figures and data outputs directly into its own
data folder (`R1_R2_data_files/`, `R3_data_files/`, `ISEF R4 data folder/`,
`ISEF R5 data folder/`, `R6_data_files/`, `ISEF R7 data folder/`,
`ISEF R8 data folder/`).
R7 is **resumable**: it caches every (μ, N, a) point to
`ISEF R7 data folder/r7_points_cache.json` and skips already-computed points,
so a long, eigh-bound run can be stopped and restarted freely.

## Results
See the per-program data folders above for all generated plots and
CSV/NPY outputs.

## Science Fair
Denver Metro Regional Science and Engineering Fair (DMRSEF)  
February 19, 2027
