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
  gap, in three tracks (`src/r8_reactive_impedance_methods.md`).
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
- ✅ R9 (Track 1): **characterizing the reactive-impedance "pole" near g_crit**
  (`src/bdg_r9_pole_characterization.py`; appended to
  `src/r8_reactive_impedance_methods.md`). Closes the one place R8 said `Z_r`
  carries orthogonal information. **Finding (revises R8's reading):** on the
  physical principal branch the enhancement `1/(1−s)` does **not** diverge — it
  rises to a **bounded ~5.7×** (s plateaus at ~0.825) and the static fixed
  point then **detaches** (fluxoid-sawtooth dead-zone) before s reaches 1, so
  the reactance grows by a bounded factor and becomes undefined rather than
  blowing up. Detachment (`g≈8.67×10⁴`) is the *same* uniform-mode bifurcation
  as R5's g_crit (its leading edge); R5's own slope-crossing-1 rule is in fact
  **unsatisfiable** in this band, so R5's published 1.1–1.2×10⁵ band likely
  comes from its coarser empirical Picard-failure criterion instead. It is
  **purely reactive** (DC gate intact near the pole, no Im χ counterpart), the
  SO correction to the loop gain stays small (its `L_eff` effect grows only
  through the pole's own sensitivity), and the enhancement cap **grows** under
  mesh refinement (5.6→9.8× over N=80→200) — so a genuine continuum divergence
  is neither demonstrated nor excluded. **Net: `Z_r` is even more redundant
  with the gap than R8 concluded** — its lone orthogonal feature is a soft,
  bounded, resolution-sensitive proximity marker to a dynamical instability,
  not a sharp new observable. Parallelised via `bdg_parallel` (measured:
  capping BLAS to 1 thread/worker is the dominant win — machine-default
  multithreaded BLAS is *slower* than 1 thread even serially for these small
  matrices; serial-vs-parallel bit-identical).
- ✅ R9 (Track 2): **finite temperature** in the BdG / Pfaffian / Kubo stack
  (`src/bdg_r9_finite_temperature.py`; `src/r9_finite_temperature_methods.md`).
  R1–R8 were all implicitly **T=0** (sharp Fermi sea, T=0 Kubo). This track adds
  Fermi–Dirac occupation of the BdG quasiparticle spectrum (self-field current,
  boundary, and Kubo Lehmann sum) as a **strictly additive** layer and asks
  whether the null results survive realistic thermal broadening. **Task 1 sets
  the scale first:** the R6-verified Δ = 0.19 meV (Chang et al. 2015; base
  temperature **20 mK**, verified against the primary source) gives
  Δ/k_B ≈ 2.20 K, so at operating conditions **T/T_c ≈ θ = k_BT/Δ ≈ 0.009** and
  the thermal quasiparticle population is **~10⁻⁴⁸** — a tiny perturbation. Every
  finite-T quantity is validated to reproduce its R5/R7/R8 T=0 counterpart as
  θ→0 (**bit-for-bit** for the current, self-consistency map, operating point,
  and B_c solver; ~10⁻¹⁶ for the Kubo χ), and the finite-T ω→0 DC gate holds to
  ~10⁻⁶ at multiple θ. **A finding en route:** the finite-T Kubo must use the
  **full** current Ĵ = ∂H/∂φ, not R8-B's kinetic-only Ĵ — the thermal
  reoccupation term samples the full current, so kinetic-only cannot satisfy the
  finite-T gate. **Finding (a null, consistent with R4–R9):** at realistic
  operating temperature **finite temperature does not matter** — the bulk
  boundary shift is unchanged (−0.131 % vs −0.132 % at T=0, still decaying
  ~N⁻¹·⁹ → 0 with size at all T), the g_crit "pole" keeps its bounded ~6×
  enhancement, and the Kubo dissipation is unchanged. Thermal effects turn on
  only **far above** operating temperature (θ ≳ 0.05, ≈110 mK, approaching the
  fixed-Δ model's validity ceiling): there the boundary shift weakens *further*
  (more null), the pole **washes out** (enh 5.9→2.3 by θ=0.15), and a new
  intraband Drude peak + sub-gap absorption appear in the dissipation. The T-axis
  is process-parallelised via `bdg_parallel` (serial-vs-parallel bit-identical;
  per-point cost θ-independent → static scheduling; 1-thread BLAS 2.2× faster
  than default here). Confirms — does not overturn — the project's null lineage.
- ✅ R9 (Track 3): **on-site disorder** in the BdG / Pfaffian / self-consistent-
  ring stack (`src/bdg_r9_disorder.py`; `src/r9_disorder_methods.md`). R1–R8
  and R9 Tracks 1–2 were all perfectly translationally ordered; this track adds
  on-site ε_i ~ Uniform[−W/2,W/2] as a **strictly additive** layer (PH-symmetry
  preserved exactly, proven and verified), with **mandatory ensemble
  averaging** (SeedSequence-based, reproducible seeding; every reported number
  carries a realization count and spread) over independent disorder draws. The
  W=0 gate reproduces every clean R4/R5 result **bit-identically** before any
  W>0 result is trusted. **Realistic W is literature-grounded** (re-verified
  primary source: Ahn, Pan, Woods, Stanescu, Das Sarma, PRMaterials 5, 124602
  (2021), arXiv:2109.00007 — σ_μ=2.5 meV Gaussian on-site fit → **W≈10–46Δ**,
  an order of magnitude above the R6-verified gap, matching that paper's own
  headline conclusion). **Finding (a genuine non-null — the one R9 track where
  the project's null lineage does NOT uniformly survive):** at realistic W the
  ensemble-averaged topological classification **smoothly smears** (up to ~30%
  of realizations near the clean boundary misclassified at the top of the
  realistic band, saturating there, not diverging); the mean spurious sub-gap
  count does **not** grow but its realization-to-realization **spread** does
  (a physical mechanism for the *sporadic* — not universal — disorder-induced
  sub-gap structure the literature already links to contested zero-bias-peak
  claims); and most strikingly, the R4/R5/R7 self-consistent inductive-feedback
  solver's **convergence rate collapses from 100% (W=0) to 5–20% (realistic
  W)** — where it does converge, the shift is ~100× the clean-model null, though
  that magnitude rests on only 1–4 successful samples per W and is reported
  with that explicit caveat. A real bug was caught and fixed in the
  convergence-failure bookkeeping itself (a silent-NaN return path bypassed
  the exception-based failure flag, undercounting the true failure rate
  ~15–18×) — the underlying shift statistics were unaffected. Serial-vs-
  parallel is bit-identical (both a cheap classification engine and the
  expensive self-consistent engine); the cheap engine is dispatched **serially**
  (measured: pool overhead exceeds its ~0.05s/task cost) while the expensive
  engine is parallelised (~2.3–4×  speedup); the pre-run time estimate (~21 min)
  undershot the real wall-clock (~77 min) because disorder-driven convergence
  failure — invisible to an isolated per-task benchmark — costs *more* than a
  successful solve, a limitation flagged explicitly rather than discovered
  silently mid-run.
- 📋 R9 (Track 4): **multi-subband extension — SCOPING ONLY, no implementation**
  (`src/r9_multisubband_scoping.md`). No simulation code was written or run; no
  file in `src/` was modified. Asks whether extending R6's single-subband 4N×4N
  Lutchyn–Oreg Hamiltonian to multiple transverse subbands is tractable and
  worth a track. **Three findings.** (1) *Literature* (verified primary
  sources): the odd-occupied-subband rule and the resulting multi-lobe phase
  diagram are established (Potter–Lee PRL 105, 227003; Lutchyn–Stanescu–Das
  Sarma PRL 106, 127001 / PRB 84, 144522), but inter-subband Rashba coupling
  **hybridizes Majorana pairs and closes the gap** (Lim et al., PRB 86,
  121103) — so the naive "M subbands → count the parity" picture is wrong, and
  omitting that term would be qualitatively wrong rather than merely
  approximate. (2) *A caveat against R6/R8, independent of whether this track
  runs*: the single-subband assumption is **not realistic** for Chang-class
  InAs/Al devices — Woods/Das Sarma/Stanescu (PRB 101, 045405) find realistic
  surface charge makes a few-subband regime *"difficult or impossible"* (≥10
  subbands occupied), and Chang et al. 2015's own measured n = 5×10¹⁸ cm⁻³
  implies many tens. This is an idealisation of the device, to be carried
  alongside the existing fixed-Δ caveat. (3) *Cost* (**measured**, not assumed):
  the correct generalization has dimension **4NM** (not 4NM²), so memory ∝ M²
  and dense-`eigh` cost ∝ **M³** — confirmed by a 1-thread benchmark on this
  machine (fitted exponent **2.98**; ratios 8/28/66/130/977× vs M³ =
  8/27/64/125/1000). The workload is *the same number of solves, each M³
  larger*, so **`bdg_parallel` gives zero additional relief**: 6 cores buy a
  flat ~2× on one large solve (33% efficiency) and 6 workers × 1 thread already
  beats 1 worker × 6 threads by ~3× at every dim ≥ 1280. **Honest verdict:**
  the multi-subband *phase diagram* is tractable to M ≈ 4–6, but multi-subband
  *self-consistent EM feedback* is **not** — M=3 turns R9 Track 3's measured
  77-min run into ~35 h, M=5 into ~6.7 days, before R7's ~5× throttling.
  **Recommendation: conditional yes at a bounded scope** — the one question no
  other track can ask (does the EM-feedback null survive more channels? the loop
  gain s plausibly grows ∝ M, which would drag R9-T1's g≈8.67×10⁴ instability
  *below* the production g=2.9×10⁴) is answerable in **minutes**, because the
  loop-gain slope needs ~30 `eigh`, not the full B_c bisection's ~2000. The full
  self-consistent sweep is explicitly excluded as out of reach.

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
    python src/bdg_r9_pole_characterization.py      # R9 Track 1 (Z_r pole near g_crit; add --regen to recompute, --quick for a fast smoke test)
    python src/bdg_r9_finite_temperature.py --full  # R9 Track 2 (finite temperature; --validate for T->0 checks only)
    python src/bdg_r9_disorder.py --full            # R9 Track 3 (on-site disorder; --validate for W=0 gate only, --benchmark for the parallel/time estimate)

Each program writes its figures and data outputs directly into its own
data folder (`R1_R2_data_files/`, `R3_data_files/`, `ISEF R4 data folder/`,
`ISEF R5 data folder/`, `R6_data_files/`, `ISEF R7 data folder/`,
`ISEF R8 data folder/`, `ISEF R9 data folder/`).
R7 is **resumable**: it caches every (μ, N, a) point to
`ISEF R7 data folder/r7_points_cache.json` and skips already-computed points,
so a long, eigh-bound run can be stopped and restarted freely.

## Results
See the per-program data folders above for all generated plots and
CSV/NPY outputs.


