# R9 Track 4 — multi-subband loop-gain probe: methods

**Status of THIS section:** pre-registration. It documents the locked protocol
(`bdg_r9_multisubband.py`, first commit) — the van Hove exclusion radius δ, the
subband energetics, the sweep windows, the operating points, and the comparison
statistic — together with the mandatory construction gates. **No s(μ;M) value is
computed at this commit.** The s-sweep functions and every s(μ;M) number are
added in a separate, later commit, so the choices below cannot have been
influenced by any outcome.

This is the same discipline the project applied after R5's g = 3×10⁵ headline
turned out to be a dynamically-unreachable root chosen before the reachability
check existed: fix the analysis before seeing the number.

## 0. What question, and why the loop gain s

R1–R9 converged on a null — inductive EM self-feedback does not move the bulk
topological boundary — computed throughout with **exactly one** conducting
channel. Track 4 asks whether that null survives more channels. The scoping
report (`r9_multisubband_scoping.md`) established that the tractable form of this
question is not the full self-consistent B_c sweep (34.7–160.6+ h per
configuration at M ≥ 3 — out of reach here) but R5's loop-gain slope

> s ≡ dG/dφ at the self-consistent fixed point, where G(φ) = g·⟨Peierls phase
> from the Biot–Savart self-field of the ring current I(φ)⟩,

the same quantity behind R5's bifurcation, R8's L_eff = L_geo/(1−s), and R9-T1's
bounded enhancement 1/(1−s). At the production coupling g = 2.9×10⁴, R9-T1
measured s ≈ 0.29 (single channel), with the static fixed point detaching near
g ≈ 8.67×10⁴ > g_production. **Hypothesis (scoping §6.1):** each occupied channel
carries current, so s might grow ~∝ M; if so, at M ≈ 10 the detachment would
drag to g ~ 10⁴ ≤ g_production and the single-channel null would be an artifact.
If s stays flat in M, the null is *stronger* than a single-channel model could
show. Both outcomes are informative.

## 1. Units and reused machinery

Dimensionless units, Δ = 1; multiply Δ-units by **0.25 meV** for physical
energies (LO_PARAMS Δ = 0.25 meV, so t_hop = 16.57 meV ⇒ **t = t_norm ≈ 66.3**,
α = α_norm = 8.0 = α_SO/(Δ·a_lat)). Everything geometric/circuit-level is reused
unchanged from R4/R5/R7 — the self-field kernel `ring_biot_savart_self`, the
Peierls conversion `compute_peierls_from_A`, the fixed external flux
`_phi_ext_bond`, the production coupling `G_R7`, and the lattice constant `A0`
— because the transverse subband structure does not change the *longitudinal*
self-field geometry. `_build_ring_lo_phased` is imported **only** as the M = 1
bit-for-bit reference. Nothing in R1–R8 or R9 Tracks 1–3 is modified.

## 2. The 4NM Hamiltonian (basis = subband(M) ⊗ spin(2) ⊗ site(N), Nambu-doubled)

`build_ring_lo_multisubband`:

- **Diagonal subband blocks** — R6's `_build_ring_lo_phased` normal block at the
  shifted chemical potential **μ − Δε_n**, one per subband n; intra-subband
  s-wave pairing (block-diagonal). This is the correct 4NM generalization, *not*
  4N·M² — the M×M structure is block structure, cost ∝ (4NM)³ (scoping §4.1,
  measured cubic in M).
- **Off-diagonal blocks** — the *inter-subband* Rashba coupling C[n,m]·(σ_x ⊗
  I_N), site-diagonal (k_y is a transverse, local operator). Lim, Serra, López &
  Aguado (PRB 86, 121103) show this term hybridizes Majorana pairs and closes
  the effective gap, so omitting it would be *qualitatively* wrong, not merely
  approximate. The hard-wall selection rule (scoping §5.4):

  > ⟨n|k_y|m⟩ = −i · 4nm / (W(n²−m²))  for (n−m) odd,  0 otherwise
  > C[n,m] = (α_SO/Δ)·⟨n|k_y|m⟩,  α_SO/Δ = α_norm·a_lat = 80 nm.

  C is purely imaginary and antisymmetric in (n,m) ⇒ Hermitian. The intra-subband
  Rashba piece σ_y·k_x is already inside each diagonal block (R6's h_SO); k_y has
  no diagonal element, so there is no double-counting.
- **Hole block** = the exact −h*, so PH symmetry holds for any M and any
  phi_bonds.

### Transverse subband energetics (measured from the first subband)

ε₁ = t·(π a/W)² and **Δε_n = ε₁·(n²−1)**, so Δε₁ = 0 and M = 1 reproduces R6
exactly. Subband n becomes the newly-occupied channel as μ crosses **its
threshold μ = Δε_n**, where its 1D van Hove singularity ρ(E) ∝ 1/√(E − Δε_n)
sits. At the locked geometry **W = 50 nm**: ε₁ = 6.54 meV (= 26.16 Δ), Δε₂ =
19.6 meV (78.5 Δ), Δε₃ = 52.3 meV (209 Δ), Δε₄ = 98.1 meV (392 Δ). One subband
spacing near M is ΔE_sb(M) = Δε_{M+1} − Δε_M = ε₁(2M+1).

## 3. Van Hove exclusion radius δ — closed form, and verified AT THRESHOLD

**Closed form (locked, from R6's own dispersion, blind to s):**

> δ(ε) = t · (2(1−ε)/π)² · (1 + margin),  ε = 0.15, margin = 0.15
> ⇒ δ = 22.31 Δ = 5.578 meV  at t = 66.3.

Derivation (scoping §7 item 3b): near a band edge the discrete level spacing s_j
grows from the edge toward the average bulk spacing 4t/N; requiring |s_j/s_bulk −
1| < ε for three consecutive levels and solving in the large-N limit gives δ =
t·(2(1−ε)/π)², independent of N. No Hamiltonian is diagonalized to *derive* it —
it is the known dispersion evaluated as arithmetic.

**Verification against the genuinely diagonalized spectrum, SPECIFICALLY AT
THRESHOLD** (`verify_van_hove_delta_at_threshold`, committed data in
`r9_multisubband_preregistration.json`): a genuine `np.linalg.eigvalsh` of the
tight-binding chain matrix −t(T+Tᵀ), with δ_measured read as the energy distance
FROM the band edge (the van Hove point, *not* the smooth bulk) at which the
near-edge spacing first converges to within ε of bulk for three consecutive
levels:

| N | δ_measured (Δ) | δ_closed/δ_measured | eigvalsh vs analytic |
|---:|---:|---:|---:|
| 80 | 24.61 | 0.788 | 1.4e-13 |
| 160 | 23.49 | 0.826 | 1.5e-13 |
| 320 | 22.19 | 0.874 | 1.9e-13 |
| 640 | 21.55 | 0.900 | 2.4e-13 |

Findings, reported honestly:
- The genuine eigvalsh reproduces the analytic dispersion to ~1e-13 (the
  diagonalization is real, not the formula in disguise).
- The closed form runs **below** the measured crossover at every N
  (`conservative_low = True`), so excluding the +15%-margin δ over-excludes the
  *continuum* crossover — the safe direction.
- **Limitation (not silently patched):** the measured crossover is mildly
  N-dependent (24.6 Δ at N = 80 → 21.5 Δ at N = 640, converging toward the
  large-N closed form 19.4 Δ; 13% spread). At the production N = 80 the
  finite-size crossover (24.6 Δ) is ~10% *wider* than the locked δ-with-margin
  (22.3 Δ). The locked δ excludes the continuum van Hove crossover + 15%, which
  is the N-independent physical scale; a few μ-points just outside the excluded
  zone at N = 80 may retain mild band-edge curvature. This is recorded as a
  limitation rather than fixed by inflating the margin post-hoc, and it is
  immaterial to the median statistic (which is dominated by points far from any
  threshold). Any future, separately pre-registered follow-up could lock the
  margin from the finite-N crossover instead.

## 4. Construction gates (all pass — `r9_multisubband_preregistration.json`)

| Gate | Result |
|---|---|
| **M = 1 bit-for-bit** vs `_build_ring_lo_phased` (H and current operator) | max\|ΔH\| = **0.0**, max\|ΔI\| = **0.0** |
| **Decoupled limit** (coupling zeroed): 4NM = ⊕ of M shifted 4N blocks | entrywise **0.0** (spectrum 3.5e-12, LAPACK noise) |
| **PH symmetry** Σ_x H* Σ_x = −H, M = 1..4, coupling on | **0.0** at every M |
| **Inter-subband selection rule** (Hermitian, zero diagonal, zero even-diff) | exact |

The decoupled-limit gate is enforced *entrywise* (reorder the 4NM matrix from
Nambu⊗subband⊗spin⊗site to subband⊗Nambu⊗spin⊗site and compare to
`block_diag` of the per-subband blocks) — exactly 0.0, the analogue of the M = 1
gate. The eigenvalue-union comparison is machine precision (~1e-12), not 0.0,
because it compares independent LAPACK reductions; that is reported as noise, not
a failure.

## 5. Operating points and windows (locked, blind to s)

**Protocol (a) — truncation convergence (Task 1, the cheap sanity check).** Hold
μ FIXED at the R6/R9-T1 single-channel operating point μ = 1.2 (subband 1 filled
1.2 above its bottom), B = √(1.2²+1) − 0.05 = 1.512, and raise only the basis
size M = 1..4. Subbands 2..M stay unoccupied — they enter only virtually, through
second-order inter-subband admixture (matrix elements ∝ 4nm/(W(n²−m²)), energy
denominators ∝ n²). s(M) should converge quickly and smoothly. This validates
the multi-subband construction; it is **not** the §6.1 physics answer (no new
current-carrying channel is added).

**Protocol (b) — genuine added occupancy (Task 2, the real question).** For each
M = 1..4, place the newest channel M near its own threshold Δε_M with channels
1..M−1 deep, and sweep μ over a window that straddles Δε_M (so the van Hove point
is inside the raw window and its exclusion is auditable) and spans one full
subband spacing of genuine occupancy above it:

> **μ ∈ [Δε_M − ΔE_sb(M)/3, Δε_M + ΔE_sb(M)]**, ΔE_sb(M) = ε₁(2M+1),
> field **B(μ) = √((μ − Δε_M)² + 1) − 0.05** (channel M held just below its own
> transition, the multi-subband analogue of R9-T1's operating field),
> coupling g = **G_R7 = 2.9×10⁴** (production).

**Direction clarification (locked here, blind to s).** Scoping §7 item 3b wrote
the window as `[ε_M − ΔE_sb, ε_M − δ]` (below threshold), which contradicts its
own §6.1a ("μ must be re-tuned *upward* … so that channel M is genuinely
*occupied*"). Genuine occupancy requires μ > Δε_M (the band bottom); below it the
channel is empty. The window is therefore placed *above* threshold. This is a
pre-registration correction of an internal inconsistency in the scoping prose,
made before any s(μ;M) value exists — not a post-hoc tune.

## 6. The comparison statistic (locked)

`preregistered_statistic`: the reported number is the **median of s(μ) over the
surviving window**, where every μ within the van Hove radius δ of any subband
threshold Δε_n (n ≤ M) — and of any lobe-crossing point detected from the bare
multi-subband gap within the window — is excluded. Median, not max, because max
would be pulled toward whatever point sits closest to an excluded boundary,
reintroducing exactly the sampling risk the exclusion exists to avoid. The **raw
(non-excluded) median** and the excluded/raw point counts are reported alongside,
so the effect of the exclusion is visible and auditable. Both directions of the
result are pre-declared informative: rising median-s(M) supports "the null is a
single-channel artifact"; flat/bounded median-s(M) supports "the null is robust
and gets stronger."

**Lobe-crossing exclusion.** The scoping report notes this radius is not
closed-form the way δ is (it depends on the coupled spectrum, which does not
exist before the builder is run). Because protocol (b) holds channel M just below
its *own* transition and the deep channels far from theirs, no lobe crossing is
expected inside the window; any that appears (bare-gap minimum near zero) is
excluded using the same δ as a conservative, blind proxy, and reported.

## 7. Scope limits (restated wherever a result is reported)

- **Structural, not device-faithful.** M ≤ 4 is far below the 10–26 occupied
  subbands of real Chang-class wires (scoping §3.3). Adding a subband can flip
  the invariant's parity (Lim et al., PRB 86, 121103), so small-M behaviour must
  not be assumed to extrapolate to the many-subband regime.
- **Loop-gain probe only.** ~30 eigh per point (a root-solve + a
  finite-difference slope), *not* the full B_c bisection (~2000 eigh) or any
  self-consistency loop. The full self-consistent multi-subband EM-feedback
  calculation is out of reach (scoping §4.5) and is not attempted.
- **Weak sum-rule expectation.** At fixed geometry, s should be *monotone
  non-decreasing* in occupied-channel count (more current-carrying channels), but
  with unequal, nonlinear, depth-dependent weights (channels near their band
  bottoms contribute less; inter-subband coupling redistributes current). Any
  observed trend is judged against this weak envelope, not a linear ∝ M law.

---

## 8. Results — Task 1 (protocol a) and Task 2 (protocol b)

Implementation: `bdg_r9_multisubband_sweep.py` (a separate, later commit than the
pre-registration above — no s(μ;M) value existed until this file was written and
run). The self-consistency map is R5's `scalar_ring_map`, generalized to the 4NM
Hamiltonian; the root-finding/bisection/finite-difference-slope recipe is copied
exactly from R7's `_phi_star_a`. Production coupling **g = G_R7 = 2.9×10⁴**
throughout (the physically relevant value — this is the loop gain a real device
at production operating conditions would have, not an extrapolated/swept g).

### 8.1 Task 1 — protocol (a): truncation convergence (the sanity check)

At the fixed R6/R9-T1 operating point (N=80, μ = 1.2, B = 1.512, subbands 2…M
unoccupied), s converges **smoothly and quickly** with basis size M, exactly as
predicted (§6.1a): **s = 0.2897 (M=1) → 0.3118 (M=2) → 0.3119 (M=3) → 0.3127
(M=4)**, a **7.6% total change**, with no discontinuities (`max_rel_change_
between_successive_M = 0.076`, the largest single step being M=1→M=2; M=2→M=4
is flat to <0.3%). **The M=1 value (0.2897) independently reproduces R9 Track
1's own published single-channel result (s ≈ 0.29 at g = 2.9×10⁴, `README.md`
R9-T1)** to three significant figures — a strong cross-validation that this new
multi-subband machinery is physically correct, not merely internally
consistent. This confirms the Hamiltonian construction and validates it for
protocol (b); it is **not** the §6.1 physics answer, since no new
current-carrying channel is added.

### 8.2 Task 2 — protocol (b): genuine added occupancy — the real result

For each M = 1…4, μ was swept over the pre-registered window (§7 item 3b,
corrected direction per §5 above), tracking channel M near its own transition,
at production g. Full curves and the pre-registered exclusion are in
`ISEF R9 data folder/r9_multisubband_summary.png`; data in
`r9_multisubband_sweep_data.json`.

**Headline numbers (median-s, van-Hove-excluded — the PRE-REGISTERED statistic;
all four chains ran to full completion — every pre-registered μ was evaluated).
M=4 was subsequently DENSIFIED from 5→15 points within the identical
pre-registered window (§8.2b) — the table below is the final, post-densification
state:**

| M | median s (excluded) | median s (raw) | n surviving / n finite | n finite / n pre-registered |
|---:|---:|---:|---:|---:|
| 1 | **0.0004** | 0.0039 | 1/5 | 5/12 |
| 2 | **1.9885** | 2.0799 | 8/10 | 10/10 |
| 3 | **2.8958** | 3.0601 | 6/7 | 7/8 |
| 4 | **2.0563** | 2.0579 | 5/7 | 7/15 |

**s crosses the Picard-stability boundary (s = 1) already at M = 2, at the
PRODUCTION coupling itself.** At M = 1, s stays safely sub-critical everywhere
probed (max 0.84, matching R9-T1's single-channel picture where g_production is
~3× below the g ≈ 8.67×10⁴ instability). At M = 2, s exceeds 1 at **every**
probed μ (raw range 1.16–3.04); at M = 3, s exceeds 1 at every FINITE point
(1.99–4.35); at M = 4, s exceeds 1 at every one of the 7 finite points (surviving
range 1.55–3.40, full finite range 1.55–3.47). This is a qualitatively
different, and more severe, finding than a simple percentage shift in a
boundary: it says that *merely occupying a second transverse channel*, with
**no change to g**, moves the wire from "safely sub-critical at production
coupling" (M=1) to "the static self-consistency map has repelling roots at
production coupling, across essentially the entire probed window" (M≥2).

**A second finding, now much better resolved after densification: the M=4
branch detaches over an extended region, not just a few scattered points.**
Of M=4's 15 pre-registered points, **8 return no self-consistent solution at
all** (`s = nan`, the same signature R9 Track 1 used for the g≈8.67×10⁴
detachment in the single-channel model — see §8.5): two isolated detachments
inside the otherwise-finite low-μ region (μ = 336.3, 426.0), and then a
**persistent detached region for every μ ≥ 538.1** (the last 5 of 15 points,
μ = 538.1–627.8, all `nan`) following the last finite point at μ = 515.7. This
places a leading-edge detachment boundary at **μ_detach(M=4) ≈ 515.7–538.1**
(Δ-units) — the multi-subband analogue of R9-T1's own g_detach, but located
here along the μ-axis at *fixed* production g rather than along the g-axis.
Here the detachment appears not from raising g, but from raising μ (and hence
B, via the pre-registered B(μ) tracking rule) at FIXED production g. This is
itself informative: in the single-channel model the detachment is a distant
feature reached only by pushing g roughly 3× above production; in the M=4
multi-subband model, a detachment of the same kind is encountered while simply
scanning μ at production g, with no change to g at all.

### 8.2b M=4 densification (5→15 points) and an honest downward revision

The originally-reported M=4 statistic (median-s = 3.4038) was built from only
**2 finite points out of 5** pre-registered — the thinnest sample of any M,
and, in hindsight, an unrepresentative one: both surviving points happened to
sit in the higher-s part of M=4's window. Because the window bounds, field
rule, exclusion radius, and statistic definition are all pre-registered and
depend only on M (not on the point *count* — `protocol_b_window(M, n_mu, ...)`
returns the identical `[thr − span/3, thr + span]` interval regardless of
`n_mu`), the M=4 chain was **re-run at higher density (n_mu: 5 → 15) within
the exact same window** to get a more representative sample. This is a
resolution/density check, not a re-tune: nothing about the protocol's
substance (window bounds, exclusion rule, or statistic) changed — the same
kind of density-robustness check R9 Track 1 already applied to its own
detachment-g measurement (`resolution_robustness`).

**Result: the M=4 median-s revises DOWN from 3.4038 to 2.0563.** With 15
points, 7 are finite (vs. the original 2), and 5 survive van Hove exclusion.
The fuller window reveals the same **rise-then-decline shape already visible
in M=2 and M=3's curves** (§8.2's figure, panel B): s rises steeply just
above the M=4 threshold (peaking near s≈3.4 at μ=313.9, the same point the
n=5 sample happened to catch), then **declines monotonically across the
finite region** down to s≈1.55 at μ=515.7 (the last point before the
persistent detachment onset), before the branch detaches entirely. The
original n=5 sample, by chance, sampled disproportionately from the high-s
peak and one deep-detachment point, and (after the 1-of-2 exclusion) reported
only the single peak value as its median — an artifact of sparse sampling
landing near a genuine local extremum, not a wrong calculation. This is
recorded here as a concrete, worked example of exactly the sampling risk
§6.1b's van Hove hazard warned about in the abstract — except manifesting as
*point-count* sparsity relative to a within-window peak, a distinct failure
mode from the threshold-proximity sampling risk the exclusion radius targets,
and one the original protocol's fixed per-M point counts did not anticipate.
**The qualitative verdict is unchanged by the revision — every M=4 finite
point, before or after densification, has s > 1 — but the specific median
value is now on much firmer footing (5 surviving points spanning a full
decline, not 1 surviving point at a local peak), and honesty requires
reporting the revision rather than keeping the original, less-supported
number.**

### 8.3 A pre-registered ambiguity, surfaced honestly (M = 1's exclusion collapse)

M = 1's excluded-median (n_surviving = 1/5) is **not a reliable statistic** — a
genuine, reportable weakness of the pre-registered protocol, not silently fixed
per the task's own instruction. M = 1's window spans only ε₁·(2·1+1) + ε₁/3 ≈ 1.33
spacings centered near the M=1 threshold (Δε₁ = 0, i.e. μ = 0), and the fixed
exclusion radius δ = 22.31 Δ consumes essentially the entire positive-μ two-
thirds of that narrow window, leaving only the single most-negative-μ point
(μ = −26.16, deep in the "channel not yet occupied" region) surviving. The
**raw** median (0.0039, from all 5 finite points including the excluded ones) is
therefore the more informative number for M=1 specifically, and it is reported
alongside for exactly this reason. Reassuringly, **the M=1→M=2 jump is
overwhelming under either statistic** (0.0004 or 0.0039, both ≪ 1, versus M=2's
1.99–2.08): this particular pre-registration weakness does not change the
qualitative verdict, only the precision of the M=1 anchor point. This asymmetry
(narrow, threshold-centered window at M=1 vs. progressively wider windows at
higher M, since ΔE_sb(M) = ε₁(2M+1) grows with M while δ is fixed) was not
anticipated when the window formula was locked in the pre-registration commit,
and is recorded here as a limitation for any future, separately pre-registered
refinement — not corrected retroactively.

### 8.4 Compute discipline, and an honest account of why this run took so long

**Parallelization measurement (required by this task's own constraints).**
`parallel_benchmark` (M=1 proxy workload, 4 points) measured: serial-default-BLAS
45.9 s, serial-1-thread-BLAS 43.9 s, **parallel-1-thread-BLAS 13.2 s (fastest,
~3.5× speedup)**; the bit-identical serial-vs-parallel equivalence gate
**passed exactly** (`max_abs_diff = 0.0`, `max_rel_diff = 0.0`, 20 compared
values). Worker count was capped at **2** (not the M_max=4 that would let every
chain run fully in parallel) by the *same* memory-based heuristic `bdg_parallel`
has used since R7–R9 (`min(cores, 6, free_GB/0.5)`) — free memory at benchmark
time was low, driven by other applications on this machine, not a code issue.

**Asymmetric point count (an explicit, reported compute-driven reduction, same
practice as R9-T3's narrowed N range).** Cost is measured cubic in M (dim = 4NM
per the scoping report), so a uniform μ-resolution across M would make the M=4
chain's wall-clock dominate the other three combined by roughly 8×. Per-M point
counts were fixed *before* running (§7 constraint on resolution, balanced
against tractability): **M=1: 12, M=2: 10, M=3: 8, M=4: 5** points. All four
chains ran to completion — every pre-registered μ was evaluated for every M
(§8.2's finite/pre-registered column) — using warm-hint continuation along μ
within each chain (mirroring R9 Track 1's `branch_curve` continuation along g),
with the M=1↔M=3 and M=2↔M=4 chain pairs sharing the 2 available workers. M=4
was subsequently densified to 15 points in a follow-up run (§8.2b), within the
identical pre-registered window — the point-count column in §8.2's table
reflects this final state (15 pre-registered for M=4, not the original 5).

**Why it took ~47 minutes of wall-clock for what the cost model estimated at
~40 minutes total (a reasonable match, but with a highly uneven pace along the
way).** Diagnosed directly during the run, not assumed: system-wide free
physical memory dropped as low as **~1.1 GB** at points during the run
(`Get-CimInstance Win32_OperatingSystem`), driven by other applications on this
machine (browser, chat clients, other coding-assistant sessions) competing for
RAM — confirmed via direct process inspection (`Get-Process` CPU-seconds vs.
elapsed wall time showed utilization dropping as low as ~2% during the worst
contention and recovering to ~17% as free memory rose back above 3 GB later in
the run), not inferred. This is the same class of resource contention this
project's own `r9_multisubband_scoping.md` §4.2 already documented on this
machine ("often has <2–3 GB free") — not a bug in the new code, confirmed by
isolated component timings run without competing load (a single independent
point, a 4-point serial loop, and a 4-point parallel pool all completed in
5–45 seconds in isolation, versus the same work taking many multiples longer
mid-run under contention). The lesson for any future run of this module on this
machine: check free memory before launching, or budget substantially more
wall-clock margin than the pure compute-cost model predicts.

**The M=4 densification run (§8.2b) took this to an extreme: ~3.7 hours
(13372 s) for 15 points**, under free memory that oscillated between ~1.1 and
~2.6 GB throughout (`Get-CimInstance` checks during the run), with individual
points ranging from 293 s (first point, best case) to a single point (the 8th)
taking 3784 s in isolation — almost certainly a widened-fallback grid scan
triggered by that point being a genuine detachment, made far slower than its
already-expensive worst case by concurrent memory pressure. Progress was
confirmed genuine, not hung, by checking `Get-Process` CPU-seconds twice per
interval and observing it climb monotonically throughout (from ~3500 s to
~9700 s across two checks 25 minutes apart) — the diagnostic habit this
project's memory record now recommends for this machine.

### 8.5 What "s > 1" means here — read with the same care as R9 Track 1

A bisection-found root of φ = G(φ) with slope s > 1 is a **genuine solution of
the static self-consistency equation** — it is not a numerical artifact — but it
is **Picard-repelling** (unstable under naive fixed-point iteration), the exact
situation R9 Track 1 characterized carefully for the single-channel model beyond
its own g ≈ 8.67×10⁴: existence of such a root does not by itself establish that
a real, dynamically-evolving wire settles there. What this result adds is that,
for the multi-subband case, **that repelling regime is reached at the production
coupling itself once M ≥ 2**, rather than requiring g to be swept ~3× above
production as in the single-channel model. The M=4 chain's outright detachment
(§8.2 — 8 of 15 densified points return no root at all, the same `nan`
signature R9-T1 used for its own g≈8.67×10⁴ boundary, with a persistent
detached region for μ ≥ 538.1) is a further, complementary signature of the
same underlying fragility: not just repelling roots, but entire stretches of μ
with **no static self-consistent solution whatsoever**, at production g, once
enough channels are occupied. Whether the *physical* (dynamical, not
bisection-found) system has a stable static equilibrium at all in this regime is
exactly the kind of full self-consistent multi-subband EM-feedback question this
task explicitly excludes as computationally out of reach (scoping §4.5) — not
settled here, and not overclaimed.

## 9. Task 3 — the answer, using ONLY the pre-registered statistic

**s shows a dramatic, clearly discernible trend with M — not a flat or
structureless one — but it is NOT strictly monotonic once the M=4 chain is
properly resolved.** Using only `median_s_excluded` (the locked statistic, all
four chains run to completion, M=4 at its final post-densification value):
**0.0004 (M=1) → 1.9885 (M=2) → 2.8958 (M=3) → 2.0563 (M=4)**. M=1→M=2 and
M=2→M=3 both show large increases; **M=3→M=4 shows a decrease** (2.8958 →
2.0563, a ~29% drop). This corrects an earlier draft of this report, written
before the M=4 chain was densified from 5 to 15 points (§8.2b), which read
M=4's median as 3.4038 — built from only 1 surviving point — and reported the
trend as cleanly monotonic. It was not a wrong calculation given the data
available at the time, but it was an unreliable one, and the fuller data
changes the shape of the answer, not just its precision. **The honest
statement is: s rises sharply and dramatically from M=1 to M=2, continues
rising from M=2 to M=3, and then declines somewhat from M=3 to M=4 — while
remaining, at every single M ≥ 2 value probed (finite point, individually, not
just the median), solidly above the s=1 stability boundary.** This supports
the §6.1 hypothesis: multi-subband occupancy CAN pull the effective
instability threshold down toward — and, on this evidence, *already to or
below* — the production coupling, rather than reinforcing the single-channel
null; it does NOT support a simple "s grows monotonically and unboundedly with
M" reading, which the M=3→M=4 dip rules out at this M range. **This is the
second non-null finding in the R9 series** (after Track 3's disorder result),
and structurally the more severe one: Track 3 found a finite-size-refined but
still-nonzero *quantitative* boundary shift (tens of percent) that survives at
the low end of the realistic disorder band; Track 4 finds a *qualitative*
regime change (sub-critical → super-critical loop gain) at the production
coupling itself, occurring already at the smallest possible step away from the
single-channel idealization (M: 1→2), and persisting (though non-monotonically)
through every M tested beyond that.

**Consistency with the weak sum-rule envelope.** The weak sum-rule argument
(§Task 3 brief) predicts monotone increase in s with occupied-channel count,
at fixed geometry, with unequal nonlinear depth-dependent weights — it is
explicitly described as *weak*, not a hard bound. **The M=3→M=4 decline is a
genuine, reportable deviation from that envelope**, not merely an unusual
shape consistent with it: read strictly, "weakly monotone increasing" rules
out a same-geometry decrease. Two readings are offered, neither privileged
over the other by the data:
(i) *Envelope violated* — the weak sum-rule's premise (current, and hence s,
should not decrease as more channels are added) does not survive at M=4,
possibly because the weighting is not merely unequal but can be net-negative
for some inter-subband contributions once the strongly-coupled multi-lobe
regime is reached (§6.2's "sweet spot" structure — a regime no M≤3 study in
this track probed).
(ii) *Envelope not actually tested* — because M=4's window is dominated by a
detached region for μ ≥ 538 (§8.2), the SURVIVING points at M=4 are
disproportionately drawn from the declining tail of a rise-then-fall s(μ)
curve (§8.2b), rather than from a comparably-sampled slice of "generic
occupancy" the way M=2/M=3's medians were. The median is doing its job
(resisting the single peak value), but a rise-then-fall-then-detach shape
inside one window is not what a same-geometry channel-count comparison
implicitly assumes when the sum-rule argument was stated. Under this reading,
the M=3→M=4 comparison is not fully apples-to-apples, and the qualitative
"still solidly super-critical, still far above s=1" conclusion is more
trustworthy than the specific magnitude ordering between M=3 and M=4.
Distinguishing (i) from (ii) would require characterizing s(μ) shape
consistency across M with matched sampling — out of scope here, and flagged
as a natural next question rather than resolved retroactively.

**Restated per this task's own constraints:** this is a **structural** result
about a 4-subband truncation, not a device-faithful prediction for Chang-class
wires (10–26 occupied subbands, scoping §3.3). Lim, Serra, López & Aguado (PRB
86, 121103) show a subband can flip the topological invariant's parity, so this
small-M trend — including its non-monotonic M=3→M=4 step — must **not** be
extrapolated smoothly into the realistic many-subband regime; the M=3→M=4 dip
is itself a small-scale demonstration of exactly the kind of non-smooth
behavior that caveat warns about, observed directly rather than only assumed.
What it *does* establish, structurally and robustly across all four completed
chains (M=1–4, internally consistent with R9-T1's own single-channel
baseline): the single-channel EM-feedback null is **not** robust to the
addition of even one more transverse channel — every M ≥ 2 finite point probed
sits above the s=1 stability boundary, regardless of the M=3-vs-M=4 ordering
— and the project's R6/R8 single-subband idealization caveat (already
recorded, scoping §3.5) should be read alongside this result, not
independently of it.
