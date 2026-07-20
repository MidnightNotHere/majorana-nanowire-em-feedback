# R8 — Reactive impedance Z_r as a second observable (methods)

Two explicit tracks, run and reported **separately**. Track A ships as R8;
Track B is scoped as its own phase (R8-B / foldable into R9) with a required
Track A↔B DC-limit cross-check.

---

## Track A — static feedback-renormalized reactance (ships as R8)

### Definition

    Z_r(ω) = i ω L_eff ,      L_eff = L_geo / (1 − s)

* **`s`** — the converged self-consistency-loop slope `s = dG/dφ|_{φ*}` of
  R5's scalar ring map `G(φ) = g (e/ℏ) a · A_self[j(φ)]` at its fixed point
  `φ*`.

  **Confirmation that `s` is literally R5's bifurcation slope** (required
  before code; not assumed). R7's cached `slope_offtransition` is produced by
  `_phi_star_a`, which central-differences `scalar_ring_map` at `φ*` — the
  identical construction to R5's `find_uniform_fixed_points` slope and to the
  quantity the README calls the "map slope ≈ 3" that is Picard-repelling at
  g = 3×10⁵. Verified numerically through the exact R5 code path:

  | g | s (this code) | cross-reference |
  |---|---|---|
  | 2.9×10⁴ | **0.2897** | == cached `slope_offtransition[0]`, N=80/μ=1.2 (bit-exact) |
  | 3.0×10⁵ | **3.0023** | == README "map slope ≈ 3, Picard-repelling" |

  So the `s` in `L_eff = L_geo/(1−s)` is the same object whose `s → 1` crossing
  *is* R5's dynamical bifurcation. It is already cached per grid point, so
  Track A reuses the converged state with **zero new eigendecompositions**.

* **`L_geo`** — geometric self-inductance of the ring loop (circumference
  `N·a`), read directly off the R5 Biot–Savart ring kernel that R5 already
  validated against the analytic thin-wire loop inductance
  (`validate_ring_biot_savart`, rel_err < 0.5 %). Geometry only, no BdG.
  Cross-checked here against `μ₀R[ln(8R/r)−2+γ−ln2]` to < 5×10⁻⁵.

### Why this is the right linearization (and not a new one)

Linearising `φ_self = G(φ_ext + φ_self)` about `φ*` gives the closed-loop
response `d(φ_total)/d(φ_ext) = 1/(1−s)`. Bond phase ∝ flux ∝ inductance, so
the self-consistent feedback renormalizes the bare geometric inductance by
that factor: `L_eff = L_geo/(1−s)`. The linearization already exists in the
codebase — it *is* `s` — so no new linearization step and no re-diagonalization
are introduced. Consequence: `L_eff` (hence `Z_r`) **diverges as `s → 1`**, at
the uniform-mode instability (R5's `g_crit` family). The clean bulk gap has no
such pole and is `g`-independent, so `Z_r` is a genuinely distinct observable.

### Honest limitation (quasi-static, reactive-only) — R6-style

Magnetostatic treatment of a gapped ground state ⇒ no dissipative channel ⇒
`Z_r` is **purely imaginary (reactive)**. The resistive (real) part is **not
defined here**; it requires a finite-ω quasiparticle (Kubo) response — that is
Track B. All ω-independent physics lives in `L_eff` (equivalently `1/(1−s)`).
The absolute Ω scale is the **model ring's** self-inductance (~pH, X_r ~ mΩ at
1 GHz), *not* an experimental mm-scale loop; the trustworthy content is the
renormalization factor, not the absolute value.

### g_crit and dead-zone discipline

* Any point with `s ≥ 1` or `g ≥ g_crit` (R5 band 1.1–1.2×10⁵) is **flagged**,
  not reported as a number (`pole_s_ge_1` / `flag_g_ge_gcrit`).
* Dead-zone points (regime `deadzone`, no static fixed point) get
  `Z_r = UNDEFINED` — there is no static equilibrium to define a static
  reactive response, mirroring "no boundary shift without a real static fixed
  point." At g = G_R7 the N=40, N=60 (μ=1.2) points are dead-zone; their gap
  shift is `nan` and their `Z_r` is `nan` — the two observables agree there.

---

## Track B — finite-ω Kubo current–current correlator (its own phase)

### Definition and the formalism caveat that governs the cross-check

The R8 prompt's literal form is `Z(ω) = iωL − iω²M²/Z_wire(ω)` with
`Z_wire(ω)` from a finite-frequency current–current (Kubo) correlator. **A
caveat that must be resolved before this can be cross-checked against Track A:**
that literal form is a *mutual-inductance reflected-impedance* circuit model
(loop ↔ wire coupling `M`), whereas Track A linearizes the wire's **self-field**
feedback (ring **self**-inductance `L_geo`). R4's own `compute_mutual_inductance`
notes `M = 0` for the perpendicular threading and that "the genuine inductive
coupling that shifts the phase boundary is the self-field feedback." So the two
are **different circuit topologies**; taken literally they compare different
things and the DC-limit cross-check would be ill-posed.

**Resolution adopted for Track B:** formulate `Z_wire(ω)` as the finite-ω
generalization of the *same* self-field response Track A linearizes. Concretely,
`s` is the ω→0 limit of a frequency-dependent gain `s(ω)` built from the wire's
current–current correlator, and

    L_eff(ω) = L_geo / (1 − s(ω)) ,   with  s(ω→0) = s  (Track A).

`s(ω)` is the dynamical version of `dj/dφ`: the Kubo current–current
susceptibility `χ_JJ(ω)` of the converged BdG Hamiltonian, folded through the
same Biot–Savart/Peierls geometric factor `g(e/ℏ)a·(kernel)` that defines `G`.
The DC-limit cross-check is then well-posed: `L_eff(ω→0)` **must** equal Track
A's `L_eff` at every shared grid point. (The literal M-form is retained only as
a separate reporting convenience once the self-consistency is established.)

### 1 — Frequency range/resolution and compute cost (given R7's throttling)

Spectral (Lehmann) representation of the current–current correlator on the
converged Hamiltonian:

    χ_JJ(ω) = (1/𝒩) Σ_{m,n} (f_m − f_n) |⟨m|Ĵ|n⟩|² / (ω + E_m − E_n + i0⁺)

* **Frequency grid.** Sweep ω from ~0 (DC, for the Track A cross-check) up to a
  few × the largest interband energy `2·max|E_m|` (so Kramers–Kronig closure
  captures all spectral weight). Proposed: `N_ω ≈ 400` log-spaced in
  [10⁻³, 10] × Δ for the physics band, plus a linear DC neighbourhood for the
  ω→0 extrapolation. Broadening `η` swept (η = 2×, 1×, 0.5× the mean level
  spacing) as the resolution-convergence knob.

* **Cost — the key scoping fact.** In the Lehmann form the ω-dependence is only
  the energy-denominator sum: given the eigenpairs `{E_m, |m⟩}` and the matrix
  elements `|⟨m|Ĵ|n⟩|²`, each ω costs `O((4N)²)` (cheap). So **a fine frequency
  grid does NOT multiply the expensive eigendecomposition** — you diagonalize
  once per grid point and the whole ω-sweep is a cheap post-process. The real
  cost per grid point is: (i) one self-consistent fixed-point solve to get `φ*`
  (seconds — the scalar map, **not** the full R7 Bc bisection); (ii) **one full
  eigendecomposition** (all 4N eigenvectors, not just the ground state); (iii)
  the current-operator rotation `Ĵ → ⟨m|Ĵ|n⟩`, `O((4N)³)`. At N ≤ 120 that is
  ~1–2 s/point. On the existing **6 clean grid points**, Track B is a
  **few-minutes** job, not a throttled multi-hour sweep — *provided* it is
  scoped to the existing grid and does not re-diagonalize per frequency. R7's
  ~5× throttling concern applies to (wrong) per-frequency re-diagonalization or
  to opening a new (B,μ,L) sweep; neither is required here. If a broader grid or
  a full continuum-limit ω-study is later wanted, that is where it becomes
  compute-limited — reportable as such rather than overstated.

### 2 — New diagonalizations required (full recompute, not Track A reuse)

Track A reuses **only** the scalar slope `s` (ground-state current + one map
derivative). Track B needs the **entire BdG spectrum and eigenvectors** at each
grid point (for `χ_JJ`), which Track A never computed and which the R7 cache
does **not** store. This is a full re-computation per grid point — understood
and stated, not a reuse of the cached state. The self-consistent `φ*` can be
warm-started from the cache, but the eigendecomposition and current
matrix-elements are new.

### Validation requirements (no Track B results without these)

* **Causality (Kramers–Kronig).** Re and Im of `Z_wire(ω)` (equivalently
  `χ_JJ`) must satisfy the KK dispersion relation across the computed band.
  Reported as a relative KK residual; a large residual ⇒ the ω-range is too
  narrow (missing spectral weight) — fix the range, don't report.
* **Frequency-resolution convergence.** Same discipline as R5/R7 N-scaling but
  for the ω-grid and broadening `η`: features must survive `N_ω` refinement and
  `η → 0`, else they are resolution artifacts.
* **DC-limit cross-check (required, gating).** `Z(ω→0)` from Track B must reduce
  to Track A's `Z_r` within tolerance at every shared grid point. Disagreement
  is treated as a **bug to find**, not a physics finding — exactly as the R6
  8.3× critical-field tension was investigated as a possible artifact before
  being reported as a genuine limitation. Passing this cross-check is what lifts
  the **PRELIMINARY** label from the Track A figure.

### Track B honest-outcome clause

If, at the frequency/broadening resolution actually achieved, the KK residual
or the ω-convergence is inconclusive, Track B is reported as
**"compute-limited, inconclusive at the grid resolution achieved"** rather than
overstating a thin sweep.

## Track B — RESULTS (BUILT & staged: `src/bdg_r8b_kubo_impedance.py`)

The scope above was executed. Outputs in `ISEF R8 data folder/`
(`R8B_KUBO_SUMMARY.md`, `R8B_FOLLOWUPS.md`, `r8b_figure_kubo_impedance.png`,
`r8b_figure_followups.png`, `r8b_*.json`). R4/R5/R7/R8-A imported read-only.

**Current operator convention (governs the gate).** `Ĵ` is the KINETIC hopping
current only — matching `ring_current_density`, which sources Track A's map and
omits the spin-orbit current. This keeps the DC cross-check a genuine
bug-catcher rather than failing by a known ~9 % SO difference (quantified below
as Task 4).

**DC-limit cross-check (the gate) — PASSES.** The Kubo static value `χ_JJ(0)`,
computed independently, reproduces the finite-difference `d⟨Ĵ⟩/dφ` (Track A's
electronic response) — so `s(ω→0) = s` and `L_eff(ω→0)` equals Track A's
`L_eff` — to **rel. mismatch ~1×10⁻⁴ at all 6 clean grid points**, and it holds
through the near-degenerate region approaching B_c down to gap = 0.004 (mismatch
actually improves). Numerically anchored: `χ₀ = 1797.14` vs finite-diff
`1797.03` at (N=80, μ=1.2). **This is the check that lifts PRELIMINARY off Track
A.**

**Complex impedance / new physics.** `Z(ω) = iωL_eff(ω)` is genuinely complex —
`Re Z ≠ 0` from interband quasiparticle absorption, absent in Track A's
static-reactive limit. But at the 1 GHz probe (ℏω/Δ ≈ 0.022) it is essentially
reactive (|Re Z / Im Z| ~ 10⁻³); dissipation only turns on at ℏω/Δ ~ 1. So
across the whole realistic read-out band Track B **confirms** Track A.

**Causality / resolution.** `χ₀` is fully converged (resolution-independent).
The finite-ω *dissipative* spectrum is broadening-limited: KK median residual
~1 % at moderate η, but the residual cannot discriminate η (broadening preserves
causality, so over-smoothing still passes KK — the naive η_c-minimisation
saturated and was discarded).

### Adversarial follow-up findings (tasks 1–5; `R8B_FOLLOWUPS.md`)

1. **Full-grid gate** — passes at 6/6 clean points; dead-zone points N/A; robust
   through the near-B_c near-degenerate region. Not a quiet failure.
2. **Onset vs local gap — does NOT track it.** Sweeping B through B_c (2·gap
   closes and reopens), the dominant Im χ peak drifts smoothly (ℏω/Δ ≈
   1.49–1.62, a Zeeman shift of the BULK pairing edge) and the topological
   pair-breaking edge carries < 0.1 % of the absorption weight. The finite-ω
   dissipation is a bulk-band feature, **not** a topological-transition marker.
   Reported as a falsifiable negative.
3. **η non-convergence is FINITE-SIZE, not intrinsic.** Decisive physics: this
   closed gapped model has no bath / no lifetime term, so η is pure
   regularisation. Numerical confirmation: Im χ roughness at fixed η **decreases
   ∝ N^−0.54** (N = 60→240) — the discrete spectrum fills toward a smooth
   continuum (cf. R5's Aharonov–Bohm artifact). No intrinsic linewidth.
4. **Spin-orbit current — propagated as an explicit error band** (disposition,
   not a footnote): including `Ĵ_SO` lowers L_eff by **−2.3 % to −4.0 %** across
   the grid. A fully self-consistent SO treatment (SO also sourcing the map) is
   scoped as **future Track C**, not silently folded in.
5. **Frequency-scale provenance.** Gap corrected to the R6-verified **Δ = 0.19
   meV** (Chang et al. 2015, arXiv:1411.6255) — the earlier 0.25 meV was the
   InSb-sourced value R6 had already corrected. Interband scale ℏω/Δ ~ 1 ≈ 46
   GHz; the realistic read-out band (0.1–10 GHz, dispersive resonators 4–8 GHz)
   sits far below it, so the reactive-only conclusion is probe-independent.
   Inherits R6's idealised-regime caveat (fixed Δ, 8.3× critical-field tension).

### Combined note — do A and B agree in the DC limit? **Yes.**

The gate passes at every shared grid point (rel. mismatch ~1×10⁻⁴). Track A's
`L_eff = L_geo/(1−s)` is confirmed as the ω→0 limit of the independent finite-ω
Kubo impedance, so **PRELIMINARY is lifted** from the Track A result. Track B
adds no trustworthy *topological* information in the read-out band; its one
qualitatively new feature (resistive part) is a bulk-band, finite-size-limited
effect. Net verdict unchanged: in the production regime `Z_r` is largely
redundant with the gap + geometry, distinct only as a proximity-to-g_crit marker.

### Performance note (compute infrastructure)

The independent follow-up sweeps are process-parallelised via
`src/bdg_parallel.py` (workers pin BLAS to 1 thread; serial-vs-parallel results
verified **bit-identical**, `max_abs_diff = 0.0`). Sweet spot ≈ 4–5 workers on
the i7-13620H (6 P + 4 E cores); small loops stay serial. See the
`parallel-sweeps` note.

## Track C — self-consistent spin-orbit current (BUILT & staged: `src/bdg_r8c_so_selfconsistent.py`)

R8-B Task 4 quantified the omitted spin-orbit (SO) current as a LINEAR error
band on `L_eff`, evaluated at the kinetic-only fixed point. Track C closes that
loop: it sources the Biot–Savart self-field from the **total** physical current
`I = I_kin + I_SO` and re-solves the self-consistent fixed point φ*_C, giving a
genuinely self-consistent slope `s_C` and reactance `L_eff^C = L_geo/(1−s_C)`.

**Why I_SO is a real source current.** The physical charge current conjugate to
the Peierls phase is `I = (e/ℏ)⟨∂H/∂φ⟩` — the full derivative, including the
Rashba term `h_SO = (iα/2)(σ_y⊗(T−T†))`, which carries the same phase as the
kinetic hopping. `ring_current_density` (kinetic only) is therefore an
approximation; Track C uses the complete source current.

**Implementation (validated-operator ratio).** The map operates on a UNIFORM
texture (uniform current), so the total per-bond current is the kinetic per-bond
current times the scalar ratio `⟨Ĵ_full⟩/⟨Ĵ_kin⟩`, reusing R8-B's already-
validated kinetic and kinetic+SO current operators (exact for the uniform mode;
the constant Ampere-scale cancels, R8-B Check 1). `so=False` reproduces Track
A bit-exact; the SO fraction reproduces R8-B's `so_current_fraction` (−8.9 %).
(A first analytic per-bond SO derivation was discarded after it collapsed to the
(T−T†) structure ≈ 0 instead of the correct (T+T†) — the ratio route sidesteps
that.)

**Results.** Across the clean grid the self-consistent SO lowers `L_eff` by
**−2.3 % to −4.0 %** — agreeing with Task 4's linear estimate to ~0.01 %, because
φ* ~ 1e-3 makes the fixed-point feedback on the ~8 % source change second-order.
**Finite-size discipline:** under fixed-L mesh refinement the SO shift **shrinks
monotonically** (−3.24 % → −1.26 %, N=80→240) — it does NOT converge to a bulk
correction but trends to 0 in the continuum, the same finite-size character as
Track A's enhancement (both `s_A` and `s_C` scale ∝ a, so `s_C − s_A ∝ a → 0`).

**Verdict.** Making the SO current self-consistent confirms — does not overturn
— Task 4: the SO correction is small, sub-dominant, and itself a finite-size
effect that vanishes in the continuum, not a bulk term the kinetic-only Tracks
A/B missed. The self-field is only rescaled ~8 %, so the (already null,
finite-size) bulk topological boundary shift is unchanged. The overall R8
verdict stands: `Z_r` is largely redundant with the gap + geometry. Outputs:
`ISEF R8 data folder/R8C_SO_SUMMARY.md`, `r8c_figure_so_selfconsistent.png`,
`r8c_so_data.json`.

---

## R9 Track 1 — characterizing the reactive-impedance "pole" near g_crit (BUILT & staged: `src/bdg_r9_pole_characterization.py`)

R8 flagged the ONE place `Z_r = iωL_geo/(1−s)` carries information the
`g`-independent bulk gap cannot: sweeping the feedback strength `g`, the
enhancement `1/(1−s)` grows steeply and the static fixed point "detaches" near
`g ~ 9×10⁴`, at the leading edge of R5's `g_crit` band — read as an `s→1`
pole. R9 Track 1 characterizes that feature quantitatively. Everything REUSES
the validated modules (Track A `bdg_r8_reactive_impedance`, Track B
`bdg_r8b_kubo_impedance`, Track C `bdg_r8c_so_selfconsistent`, and
`bdg_parallel`) as read-only imports; nothing in R1–R8 is modified. Outputs in
`ISEF R8 data folder/` (`R9_TRACK1_SUMMARY.md`, `r9_figure1_bounded_pole.png`,
`r9_figure2_kubo_trackC.png`, `r9_pole_data.json`).

**Headline (honest, and it revises R8's reading).** At every lattice actually
simulated the reactance does **not** diverge. Tracked with Track A's own
validated bisection path `_phi_star_a` (a root-finder, **not** Picard, so it
finds the fixed point regardless of stability) under warm-hint continuation,
the physical *principal* branch's enhancement `1/(1−s)` climbs smoothly from
~1.1 (`g=1×10⁴`) to a bounded maximum of **~5.71×** (N=80, μ=1.2, a=10 nm,
B=ana−0.05; `s` plateaus at **~0.825**), and the static fixed point then
**detaches** at `g≈8.67×10⁴` via the fluxoid **sawtooth** structure of the
self-consistency map (R7's "razor-thin no-root windows"), **before `s` ever
reaches 1**. The reactance grows by a bounded factor and becomes UNDEFINED (no
static equilibrium); the branch that reappears at higher `g` is already
Picard-repelling (`s>1`). So R8's "`L_eff` diverges as `s→1`" is, on the
accessible static branch, a **bounded rise cut off by a detachment/dead-zone**
— the outcome the R9 brief explicitly allowed for ("a pole that doesn't
survive scrutiny would itself be an important finding").

**Task 1 — scaling, fit, and resolution.** The power-law fit
`1/(1−s) ~ (g_c−g)^(−p)` is honestly ill-posed with a free `g_c` (`p` and `g_c`
trade off): best free fit `p≈1.7` at `g_c≈1.26×10⁵` (R²=0.9998, landing near
R5's band), but a fit anchored at the detachment itself gives a much poorer
`p≈0.34` (R²=0.92). The rising trend, extrapolated, IS consistent with a pole
at `g_crit`, but the physical branch detaches before ever reaching it — a
clean simple pole is **not realized** on the accessible branch. Root
detection near the boundary depends on the scan **window/span** used (a
narrow warm-hint bracket can miss a root a wide multi-period scan still
finds), though a controlled density-only comparison (same ±2-period span,
n=61/121/241) was fully converged (spread=0) — so every detachment number here
is accurate to its scan window, not absolute.

**Task 1.3 — detachment vs R5's g_crit (sharper than "leading edge vs upper
edge").** Both are the same loss of an attracting static uniform fixed point
(R5's uniform-mode bifurcation), but R5's own `g_crit_predicted_from_slope`
rule (`bdg_r5_selfconsistent_ring.py:1137`, interpolates the smallest-|φ| root
slope crossing 1) is **literally unsatisfiable** at fine `g`-resolution here:
the root has `s=0.84` at `g=1.0×10⁵` (still <1), vanishes again over
`1.03–1.17×10⁵`, and reappears only at `g=1.2×10⁵` already at `s=1.17` (>1) —
no adjacent scanned pair straddles `s=1`. So R5's *published* "g_crit ≈
1.1–1.2×10⁵" almost certainly comes from `g_crit_measured` (R5's coarser
**empirical** Picard-iteration-failure criterion, an upper bound by
construction), not a literal slope-crossing root — the same dead-zone this
work maps prevents pinning either definition to one sharp number. Verdict: the
pole marks the bifurcation's onset (leading edge, `g≈8.67×10⁴`); R5's g_crit a
point ~30% higher, past the same dead-zone.

**Task 1 — finite-size discipline (a genuine correction to an earlier
draft).** Two distinct statements must be kept apart. At **fixed production
g**, Track A already showed the enhancement is a finite-size effect
(`s∝a→0`, enhancement→1) — unchanged. But the **max enhancement reached before
detachment** (at the *moving* `g_pole∝1/a`) **grows** under fixed-L mesh
refinement: 5.6→6.6→7.7→8.9→9.8 over N=80→200 (L=0.8 μm, per-N g-grids
centered on each expected detachment for a fair comparison). The pole
**sharpens** as the mesh refines — the ~6× cap is **not a universal, N-
independent bound**. This does not establish a genuine continuum divergence
either: at every N studied (60–200) the branch still detaches at a *finite*
enhancement before `s=1`, and no continuum extrapolation is claimed (it would
need much larger N and carries R7/Track A's fixed-hopping caveat). The honest,
robust claim: at every lattice simulated the reactance is bounded and detaches
before `s=1`; the continuum fate of the cap is left open.

**Task 2 — is the pole in the dissipative (Kubo) channel? No.** Reusing Track
B's validated machinery (`_converged_state`, `current_operators`, `kubo_chi`)
at the `g`-dependent converged state from production up toward detachment: the
ω→0 DC gate (Track B ≡ Track A) **holds throughout** (worst relative mismatch
9.6×10⁻⁴ ≪ 10⁻³ — does not degrade or break down as `s` approaches its
plateau), and Im χ develops **no matching feature** — absorption weight at the
topological edge stays ≤7.2×10⁻⁴ with the dominant peak fixed at the
bulk-pairing scale. The pole lives entirely in Re χ / `L_eff`; it is purely
reactive.

**Task 3 — does Track C's SO self-consistency change near the pole?** Cleanly
separates two things. The SO correction to the **loop slope `s`** is
essentially constant at **−7.1% to −8.2%** across the whole `g`-range — the
same value Track C found at production coupling; the SO *physics* does not
grow or change character near `g_crit`. But its effect on **`L_eff`** grows
from −3.2% (production) to −23.5% (near detachment) — not because the SO
physics amplifies, but because `L_eff=L_geo/(1−s)` is a steep function of `s`
near the plateau, so a constant `s`-shift maps to a growing `L_eff` shift.
This is exactly the effect anticipated ("the loop itself becomes more
sensitive there, independent of the SO physics"): yes, but the amplification
is **geometric** (a property of the bounded pole), not new SO physics.

**Task 4 — parallelization (measured on this sweep).** BLAS is
scipy-openblas 0.3.31 (DYNAMIC_ARCH, Haswell, MAX_THREADS=24),
multithreaded-capable with **no thread-cap env var set by default**; 10
physical cores, ~2 GB free at bench time ⇒ workers (each re-importing the
R4–R8 stack, ~0.4 GB) are RAM-limited to `n_proc≈3`, matching
`bdg_parallel.default_workers()`. Benchmark (12 independent `g`-points via
`s_worker`): serial/default-BLAS 84.8 s; serial/1-thread-BLAS 85.9 s (≈1.0×);
parallel/default-BLAS 36.9 s (2.3×); **parallel/1-thread-BLAS 13.7 s
(6.2×, fastest)**. Two separable, measured effects: for these small 4N=320
matrices, machine-default multithreaded BLAS is *slower* than 1 thread even
serially (thread-sync overhead ≫ the tiny eigendecomposition), and
process-parallelism adds a further ~n_proc× — so "more BLAS threads" is
actively harmful here. Equivalence gate (serial vs parallel, both 1-thread
BLAS): `max_abs_diff=0` over 72 values, bit-identical, the same standard R8
held `bdg_parallel` to. Config used for the R9 sweeps: `default_workers()` +
1-thread BLAS per worker.

**Task 5 — is the pole a measurable, realistic device quantity? Marginal,
leaning negative.** Using R4/R5's plausibility decomposition
`g = Q·(L_k/L_m)` and the R6-verified InAs/Al scales (Δ=0.19 meV ⇒ interband
~46 GHz), the leading-edge pole `g≈8.67×10⁴` is ~3.0× the production coupling
— nominally reachable by pushing resonator `Q` or the kinetic-inductance
ratio (e.g. `L_k/L_m=30 ⇒ Q≈2889`, within the achievable Q~10³–10⁶ envelope).
Three reasons it is unlikely to be a practically useful diagnostic: (1) the
feature is **bounded** (a few-fold enhancement, not a divergence) — no sharp
resonance to lock onto; (2) it sits on a **detachment/dead-zone** whose
location is scan-window-dependent even numerically — a soft, hard-to-place
crossover in a real device; (3) reaching it means operating right at the
dynamical-instability threshold, where by construction no stable static
equilibrium exists — the wire would self-oscillate rather than sit at a
measurable operating point. Its honest value is as a *qualitative*
proximity-to-instability indicator, not a sharp new observable.

**Net effect on the R8 verdict.** R9 Track 1 **weakens** R8's one exception.
The apparent pole is not a realized divergence but a bounded (at every finite
lattice) enhancement cut off by the fluxoid-sawtooth detachment — the reactive
face of R5's loss of a static equilibrium. It is purely reactive (no
dissipative Kubo counterpart, DC gate intact), the SO correction to the loop
gain stays small (its `L_eff` effect grows only through the pole's own
geometric sensitivity), and the cap climbs with mesh refinement so a genuine
continuum divergence is neither demonstrated nor excluded. `Z_r`'s single
orthogonal feature is a soft, bounded, resolution-sensitive proximity marker
to a dynamical instability — informative physics about *where the static
self-consistent description breaks down*, but not a sharp new device
observable, and **even more redundant with the gap** than R8's cautious
verdict allowed.
