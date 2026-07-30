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
  headline conclusion). **Finding — a nuanced, W-dependent verdict, reached
  after TWO rounds of post-hoc self-correction:** the ensemble-averaged
  topological classification **smoothly smears** with W (up to ~30% of
  realizations near the clean boundary misclassified at the top of the
  realistic band, saturating there); the mean spurious sub-gap count does
  **not** grow but its realization-to-realization **spread** does (a physical
  mechanism for the *sporadic*, not universal, disorder-induced sub-gap
  structure the literature already links to contested zero-bias-peak
  claims) — both unaffected by anything below. The R4/R5/R7 self-consistent
  inductive-feedback boundary shift is where two corrections happened: **(1)**
  a genuine bug (the disordered solver centered its search on the fixed
  *clean* boundary instead of each realization's own already-displaced rigid
  boundary, so at realistic W it searched in the wrong place) was caught and
  fixed, raising a first-pass 5–20%-success/n=1–4 "ensemble" to a properly-
  converged N=80 result showing the shift growing smoothly from +0.2% (W=1Δ)
  to +29%→+45%→+50%→+57% across the realistic band; **(2)** a follow-up
  N-independence check (the same finite-size scrutiny already applied to
  R4/R5/R7/R8-C) found that N=80 result itself substantially **inflated** the
  bulk effect: re-run at N=120/160 (the nearest verified-clean system sizes;
  N=60/90/100 hit a genuine or lucky-path dead-zone at this g), the W=46Δ
  shift collapses from +56.5% to +9.6–10.8% and its significance from 3.74σ
  to only 1.2–1.5σ from zero — **a confident non-null claim at high W does
  NOT survive** — while the W=10Δ shift stays real and highly significant
  (5.9–7.8σ) at every N, settling around +20–23%. The residual
  self-consistency convergence-failure rate is *also* a finite-size effect
  (reaches 100% success by N=160), and a direct check found **no**
  correlation between failure and boundary-displacement magnitude (ruling
  out a specific survivorship-bias concern). **Both superseded rounds are
  preserved, not deleted**, in `r9_disorder_data.SUPERSEDED_task4_3.json`
  and `r9_disorder_N_independence.json`/the methods note. Serial-vs-parallel
  is bit-identical (all engines); per-task cost did **not** drop after either
  fix (successful solves cost as much as the fast failures they replaced, and
  larger N costs ~N³ more), so pre-run estimates undershot real wall-clock at
  every stage — flagged explicitly each time rather than discovered silently
  mid-run.
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
  gain s plausibly grows with occupancy, which would drag R9-T1's g≈8.67×10⁴
  instability *below* the production g=2.9×10⁴) is answerable in **hours, not
  days**, because the loop-gain slope needs ~30 `eigh` per point, not the full
  B_c bisection's ~2000. **A post-hoc catch during scoping (own subsection,
  §6.1a–b):** this project's own R9-T1 measurement shows s is peaked, not
  smooth, near gap closings, and multi-subband adds a van Hove density-of-states
  singularity at every subband threshold — so a naive 3-point read of s vs M at
  a fixed μ would risk sampling an artifact rather than the real occupancy
  trend. The probe is scoped instead as a pre-registered μ-window / exclusion /
  median-statistic protocol (committed before any s value is computed), costing
  ~2 hours rather than minutes. The full self-consistent sweep is explicitly
  excluded as out of reach.
- 🔴 R9 (Track 4) — **NON-NULL: the loop gain s grows dramatically with
  occupied-subband count M, crossing the Picard-stability boundary (s=1) at the
  PRODUCTION coupling itself once M≥2 — though the trend is NOT strictly
  monotonic once M=4 is properly resolved (see below)**
  (`src/bdg_r9_multisubband.py`, `src/bdg_r9_multisubband_sweep.py`,
  `src/r9_multisubband_methods.md §8–9`). Pre-registration first (own commit, no
  s data): builds the 4NM multi-subband ring Hamiltonian
  (subband(M)⊗spin(2)⊗site(N), Nambu-doubled — dimension 4NM per the scoping
  report, *not* 4NM²) with the mandatory inter-subband Rashba coupling (Lim et
  al.'s gap-closing term); locks the van Hove exclusion radius δ, subband
  energetics, protocol-(a)/(b) windows, and the median-with-exclusion statistic
  **before any loop-gain value is computed**. All four construction gates pass
  exactly (M=1 bit-for-bit vs `_build_ring_lo_phased`; decoupled-limit entrywise
  0.0; PH symmetry 0.0 at M=1–4; selection rule exact), and δ is verified
  against a **genuinely diagonalized spectrum at threshold** (not just the
  smooth bulk): the closed form runs 10–21% below the measured near-edge
  crossover at N=80–640, confirming the +15% margin is conservative.
  **Task 1 (protocol a, sanity check):** at fixed μ=1.2 (N=80), s converges
  smoothly with M (0.2897→0.3118→0.3119→0.3127, a 7.6% total change, biggest
  step M=1→2) and the M=1 value independently reproduces R9-T1's own published
  s≈0.29 at g=2.9e4 to 3 significant figures — cross-validating the new
  machinery. **Task 2 (protocol b, the real test — all four chains ran to
  completion):** median-s (pre-registered, van-Hove-excluded) at production
  g=2.9×10⁴: **M=1: 0.0004 → M=2: 1.99 → M=3: 2.90 → M=4: 2.06**. At M=1, s
  stays safely sub-critical everywhere probed (max 0.84); at M≥2, s **exceeds 1
  at every single finite probed μ, individually** — a qualitative regime
  change, not a percentage shift, reached by adding a single additional
  transverse channel with **no change to g**. M=4's number was revised DOWN
  from an initial sparse estimate of 3.40 (1 surviving point of 2, out of 5
  pre-registered) after densifying that chain to 15 points within the
  *identical* pre-registered window: the fuller data shows the same
  rise-then-decline shape already visible in M=2/M=3 (peaking ~3.4 near
  threshold, declining to ~1.55 before detaching), and 5 of 7 finite points
  survive exclusion — a much better-supported number, reported instead of the
  original one, not alongside it. A second, complementary finding at M=4: 8 of
  15 densified points **detach entirely** — no self-consistent solution at all
  — including a persistent detached region for every μ≥538, placing a leading-
  edge detachment boundary at μ_detach(M=4)≈516–538, the multi-subband
  analogue of R9-T1's own g≈8.67×10⁴ boundary but reached here at *fixed*
  production g by raising μ instead. One pre-registered protocol weakness is
  surfaced honestly rather than patched: M=1's narrow, threshold-centered
  window leaves only 1/5 points outside the exclusion radius, so its raw
  median (0.0039) is reported alongside — the M=1→M=2 jump is overwhelming
  under either number. Parallelization discipline measured (not assumed):
  parallel-1-thread-BLAS fastest (13.2s vs 44–46s serial, ~3.5×),
  serial-vs-parallel equivalence exact (`max_abs_diff=0.0`).
  **Task 3 verdict:** s rises sharply M=1→M=2→M=3, then **declines** M=3→M=4
  (2.90→2.06, ~29%) — a genuine deviation from the weak sum-rule's monotone-
  increase envelope, not merely an unusual shape consistent with it (two
  readings offered in the methods note, neither privileged: the envelope's
  premise may not hold once inter-subband coupling is strong, or M=4's
  surviving points may be a non-generic slice of a detachment-dominated
  window). What is NOT in doubt: **every finite point at every M≥2, individually,
  sits above s=1** — the qualitative sub-critical→super-critical transition at
  production coupling is robust even though the specific M=3-vs-M=4 magnitude
  ordering is not. **This is the second non-null result in the R9 series**
  (after Track 3's disorder finding) and structurally the more severe one: a
  qualitative regime change, not merely a quantitative boundary shift. Read
  with the same care R9-T1 applied to its own s>1 branches: a bisection-found
  root with s>1 is a genuine static solution but Picard-repelling (dynamically
  unreachable as drawn), not proof the physical wire is unstable — the full
  self-consistent multi-subband question remains out of reach (scoping §4.5)
  and is not overclaimed here. **Explicitly structural, not device-faithful:**
  M≤4 is far below real devices' 10–26 occupied subbands, and Lim et al. show a
  subband can flip the invariant's parity, so this trend — including its
  non-monotonic M=3→M=4 step, itself a small-scale demonstration of that same
  non-smoothness warning — must not be extrapolated smoothly to the
  many-subband regime.
  **Post-hoc investigation (`r9_multisubband_methods.md §10`, no re-run,
  no protocol change):** separated "no root at all" (genuine detachment) from
  "root found but van-Hove-excluded" per M — detachment rate at coarse
  (pre-registered-grid) resolution is 58%(M=1) → 0%(M=2) → 12.5%(M=3) →
  53%(M=4), itself non-monotonic, not a clean increasing-with-M trend. **The
  M=2 and M=3 figures here were later found to be undercounts — see the
  correction below.** Two distinct detachment mechanisms found via the
  bare-gap signature: a few isolated points show a local gap dip (weakly
  consistent with, but below the pre-registered lobe-detector's own 30%-of-
  median threshold for, an inter-subband gap closing), while most detachment
  — including M=1's and M=4's persistent high-μ tails — shows no gap feature
  at all, more consistent with the self-consistency map's search range being
  exceeded (an R9-T1-style existence boundary, reached via μ/B here instead of
  g). M=1's n=1 survivor is now explained by two COMPOUNDING causes, not one:
  a structural window/δ-ratio effect (δ/span = 28.4% at M=1, the highest of
  any M) plus M=1's own high detachment rate independently removing most of
  the window before exclusion even applies — an explicit statistical-weight
  caveat for the M=1 point is now recorded in the methods note. The M=3→M=4
  dip is NOT primarily a van-Hove-exclusion artifact (raw and excluded medians
  are within 5% of each other at every M) but the verdict is still genuinely
  mixed: consistent with the sum-rule's explicit non-guarantee across
  different M-specific (μ,B) windows, compounded by a real coverage asymmetry
  (M=4's high detachment rate leaves its surviving sample never reaching the
  upper ~18% of its own window, unlike M=3's near-complete coverage) —
  reported as an open, disclosed question, not forced to a clean answer.
  **A third follow-up (`bdg_r9_multisubband_bracket_scan.py`,
  `r9_multisubband_methods.md §12`) resolved the detachment MECHANISM
  directly:** an instrumented, ~9×-denser bracketing scan across both
  onsets (M=1: μ 21.40→34.0; M=4: μ 515.70→545.0) shows that **every single
  detached point finds a genuine sign change in the self-consistency map
  and bisects it — none return "no crossing at all."** Every one is rejected
  because the converged crossing's |f(φ*)| lands above the pre-existing
  1e-4 acceptance threshold — the same "sawtooth jumps across the identity
  line without crossing it" mechanism already established as genuine
  (non-numerical-bug) physics in R4/R5/R7, just located on the φ-sweep axis
  rather than the μ-axis `bare_gap` had checked (§11.2's ruling-out is
  corrected on this one point; its other two rulings-out stand). M=1's
  transition shows real intermediate structure, including one isolated
  accepted-but-wildly-unstable root (s=−580) sandwiched between rejected
  points; M=4's transition, at the resolution tested, is a single sharp
  step. **This s=−580 outlier has NO relationship to the pre-registered M=1
  statistic** — it only exists in this diagnostic dense scan (not among the
  12 pre-registered μ values), so it contributes nothing to M=1's reported
  median; the actual surviving point behind that number remains μ=−26.16,
  s=0.000435, small and well-behaved.
  **A fourth follow-up (`r9_multisubband_methods.md §13`) answered the
  remaining open question — is the sawtooth mechanism unique to M=1/M=4?
  No.** The identical bracket scan run on M=2 (previously 0% detached) and
  M=3 found the same multi-crossing structure everywhere — M=3's densest
  point shows 7 sign changes and 3 simultaneously-accepted roots, *more*
  structure than any M=1 or M=4 point — and **both windows have hidden
  detachments the coarse grid stepped clean over**: M=2 at μ=183.43 (in the
  exact gap between its coarse grid's neighbors), and M=3's known "isolated"
  failure at μ=183.11 turns out to be two consecutive points (183.11 and
  188.92). **Corrected detachment rates, same methodology, dense-bracket
  resolution: M=2 0%→7.7%, M=3 12.5%→15.4%** (propagated back to the
  original §10.1 table with a visible correction box, project convention —
  not silently overwritten). The real signature distinguishing M=1/M=4 from
  M=2/M=3 turns out not to be *how much* sawtooth structure exists (M=3 has
  the most) but the **rejection rate per candidate crossing found**: M=1/M=4
  reject ~90% of what they find; M=2/M=3 reject only 55–65%, and usually have
  a 2nd or 3rd candidate to fall back on when one is rejected. Why M=1/M=4's
  map specifically rejects more of its candidates — a question about the
  shape of the current-phase relation itself — remains open. None of this
  changes any pre-registered number in §8–9; diagnosis only throughout.
  **Beyond the proxy (`bdg_r9_multisubband_bc.py`, `r9_multisubband_methods.md
  §14`): does the FULL self-consistent boundary calculation — not just the
  loop-gain slope — actually shift at M=2?** Ports R5/R7's exact
  Pfaffian-bisection algorithm (`ring_boundary_Bc_pfaffian_selfconsistent`)
  and R9-T1's Anderson-mixing seed-independence test onto the M-subband
  Hamiltonian, gated bit-for-bit at M=1 against both R4's Pfaffian sign and
  R5's own boundary result (`diff=0.0` on both), plus a fresh check that the
  uniform-manifold reduction Track 4's own proxy relies on still holds
  exactly at M=2 (confirmed to 1.9×10⁻¹³). At one representative point
  (μ=131.762, the s(μ;M=2) peak, N=80): the proxy's root is confirmed
  genuine and seed-independent via Anderson mixing (matching R9-T1's own
  method) but Picard-repelling, exactly as the proxy said. The full boundary
  bisection **does locate a genuine Bc=51.2044** — a **−0.82% shift** versus
  the clean same-lobe boundary, **~13× R4/R5/R7's own single-channel shift**
  at the identical production coupling — inside a genuine 4-lobe Pfaffian
  oscillation confirmed present in the CLEAN Hamiltonian too (the
  literature's predicted multi-subband phase structure, not a
  self-consistency artifact). **But this shift is NOT N-independent, and not
  in the way the single-channel case was:** N=60 and N=120 (the mandatory
  N-check) both return **`genuine_dead_zone`** — no self-consistent solution
  found by ANY of 5 independently-jittered search paths, robustly, across the
  same B-range where N=80 succeeded — the OPPOSITE of R4/R5/R7's own
  N-dependence (there, dead-zones shrink and vanish as N grows; here N=80's
  solution is bracketed by dead-zones on both sides, fail→succeed→fail, not
  a converging trend). **Verdict: outcome (c), genuinely more complicated
  than a clean confirm or null** — the −0.82% number cannot be reported as a
  continuum-limit physical quantity, but the underlying finding is not
  weaker for it: the very *existence* of a static self-consistent boundary
  at M=2 is itself fragile to finite-size variation, a qualitatively more
  severe failure of the static picture than any percentage shift would be,
  and one the cheap proxy alone could never have shown (it only ever asks
  about a root at one fixed B, never whether a *boundary* exists at all).
  This deepens, rather than merely restates, the proxy's own s>1 warning.

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


