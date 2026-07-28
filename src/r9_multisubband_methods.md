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

*(A deeper, quantitative investigation of exactly why n=1 — separating the
window-geometry effect above from a second, independent cause — is in §10.3.)*

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

## 10. Post-hoc investigation: detachment structure, M=1's reliability, and
##     the M=3→M=4 dip

**Status of this section:** a post-hoc *investigation* of the already-computed,
already-reported result in §8–9. It reads fields the statistic function
(`preregistered_statistic`) already computes and stores in
`r9_multisubband_sweep_data.json` (`n_raw_finite`, `n_excluded_by_vh`,
`surviving_mask`) and the already-recorded `bare_gap` column — it does **not**
re-run any sweep, re-derive δ, change any window, or alter any median value
already reported in §8–9. Every number below is read directly from the
existing data.

### 10.1 Task 1 — detachment (no root at all) vs. van Hove exclusion (a root
###      exists but is thrown out), separated per M

These are two different things that §8.2's "n surviving / n finite" framing
does not visually separate. **Detached** means `_phi_star_ms` found no sign
change in `f(φ) = G(φ) − φ` even after the widened 3-period/121-point fallback
scan (the same detachment signature R9 Track 1 used, identical code path) —
no self-consistent solution exists there at all, independent of the van Hove
rule. **Excluded by van Hove** means a root *was* found (finite s) but the
point sits within δ of a subband threshold, so the pre-registered protocol
throws it out. Reading `n_raw_finite`, `n_excluded_by_vh`, and
`n_pre_registered − n_raw_finite` directly from the stored statistic:

| M | n pre-registered | n finite (root found) | **n detached (no root)** | **detachment rate** | n excluded by van Hove | n surviving |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 12 | 5 | **7** | **58.3%** | 4 | 1 |
| 2 | 10 | 10 | **0** | **0%** | 2 | 8 |
| 3 | 8 | 7 | **1** | **12.5%** | 1 | 6 |
| 4 | 15 | 7 | **8** | **53.3%** | 2 | 5 |

**The detachment rate does NOT increase monotonically with M** — it is
58.3% (M=1) → 0% (M=2) → 12.5% (M=3) → 53.3% (M=4), the same non-monotonic
shape already found for median-s itself (§9), not a clean "more channels ⇒
more detachment" trend. M=2 is the cleanest chain of the four (every single
pre-registered point returned a genuine root). This is reported as found —
the hypothesis in the task brief that detachment rate might increase steadily
with M is **not supported** by this M range.

### 10.2 Two distinct detachment mechanisms, distinguished by the bare-gap signature

`bare_gap_ms` records the non-self-consistent bulk gap (φ_self = 0) at *every*
point, independent of whether a self-consistent root exists there — so it can
be checked at detached points too, to ask whether each detachment coincides
with a genuine local dip in the underlying spectrum (weakly consistent with a
Lim-et-al.-type inter-subband gap closing) or not (more consistent with the
self-consistency map's search range being exceeded — an "existence boundary"
effect, the same character as R9-T1's own g-driven detachment, but reached
here by raising μ/B at fixed g instead).

Checking every detached point against its immediate neighbors' bare gap:

| M | μ (detached) | bare_gap | local minimum vs. neighbors? | ≥30% of window-median gap (the pre-registered lobe-detector's own threshold)? |
|---:|---:|---:|---:|---:|
| 1 | 11.89 | 0.0489 | No (flat, ≈ neighbors 0.0489/0.0490) | No |
| 1 | 30.91–78.48 (6 pts) | 0.0490 (flat) | No | No |
| 3 | 183.11 | 0.0666 | **Yes** (neighbors 0.1136, 0.3450) | No (30%-of-median = 0.0344; 0.0666 is above it) |
| 4 | 336.33 | 0.1199 | **Yes** (neighbors 0.7546, 0.2137) | No (30%-of-median = 0.0683; 0.1199 is above it) |
| 4 | 426.01 | 0.2276 | **Yes** (neighbors 0.5398, 0.3642) | No (30%-of-median = 0.0683; 0.2276 is above it) |
| 4 | 470.86 | 0.3101 | No (monotonic decline 0.3642→0.3101→0.2686) | No |
| 4 | 538.12–627.81 (5 pts) | 0.2089→0.1385 (smooth decline) | No | No |

**Two qualitatively different patterns are visible, though neither reaches the
strength of a confirmed lobe crossing:**

- **(A) "Local-dip-associated" detachment** — μ=183.11 (M=3), μ=336.33 and
  426.01 (M=4): the bare gap dips *relative to its immediate neighbors* at
  exactly the detached point, which is weakly suggestive of proximity to an
  inter-subband gap-closing feature (the Lim et al. mechanism this track's
  scoping named as mandatory to include). **This should not be overstated:**
  the same-code, already-pre-registered `_detect_lobe_crossings` function
  (used operationally to add extra exclusions in §8.2) requires a gap below
  **30% of the window's median gap**, and *none* of these three points meet
  that bar (`lobe_crossings_detected: []` for every M, already in the stored
  data) — consistent with the honest reading that these are real local
  minima, but modest ones, not the sharp, deep gap closing a true lobe
  crossing would produce. "Weakly consistent with" is the right strength of
  claim here, not "confirms."
- **(B) "Existence-boundary" detachment** — μ=11.89 and the entire persistent
  tail μ≥30.91 (M=1); μ=470.86 and the entire persistent tail μ≥538.12 (M=4):
  no local gap dip at all — the bare gap is either flat (M=1) or continues a
  smooth, already-established monotonic trend straight through the
  detachment onset (M=4), with no distinguishing local feature. This pattern
  is more consistent with the self-consistency map's search bracket (even
  the widened fallback) failing to contain a root once B grows large enough
  as μ climbs — the same *character* of detachment R9 Track 1 found by
  raising g, here produced instead by raising μ (and hence B, via the
  pre-registered tracking rule) at fixed production g.

**Resolution caveat, stated plainly:** the existing μ-grid spacing (9.5–35 Δ
depending on M) cannot resolve whether any of these transitions are sharp or
gradual at finer resolution than the grid itself provides — this
characterization is a correlational read of where detachment coincides with a
gap feature, not a located phase boundary. Per this task's explicit
instruction not to re-run the protocol, no additional points were computed to
sharpen it; a genuine sharp-vs-soft characterization (mirroring R9-T1's own
g_detach characterization) is flagged as a natural follow-up, not attempted
here.

### 10.3 Task 2 — why M=1's statistic rests on n=1, and why that specific
###      number is expected structurally, not a code error

M=1's window is `[thr − span/3, thr + span] = [−26.16, 78.48]` with
`thr = Δε₁ = 0` and `span = 78.48`. The van Hove radius is fixed at
`δ = 22.31` regardless of M. **δ/span is the highest of any M and shrinks
monotonically with M** (M=1: 28.4%, M=2: 17.1%, M=3: 12.2%, M=4: 9.5%,
computed directly from `subband_spacing(M)` — span grows with M while δ is
fixed by construction). Because M=1's threshold sits at the absolute
reference μ=0 and the window's lower edge is only span/3 = 26.16 below it,
almost the entire window from μ=0 up to μ=+22.31 (most of the window's
positive-μ two-thirds) falls inside the exclusion zone, and part of the
negative-μ third does too — leaving only a narrow sliver, μ∈[−26.16,−22.31],
one grid step wide, as geometrically capable of surviving exclusion *at all*,
independent of whether a root exists there. This is exactly the
"narrow-window-relative-to-δ" structural explanation the task brief
anticipated, and it is confirmed quantitatively, not just plausible.

**But that alone does not explain "only 1 of 12" — it compounds with §10.1's
detachment finding.** The full 12-point classification:

| μ | s | classification |
|---:|---:|---|
| −26.16 | 0.000435 | **SURVIVES** |
| −16.65 | 0.000938 | excluded (van Hove) |
| −7.13 | 0.003872 | excluded (van Hove) |
| 2.38 | 0.354924 | excluded (van Hove) |
| 11.89 | nan | **detached** |
| 21.40 | 0.841974 | excluded (van Hove) |
| 30.91–78.48 (6 pts) | nan | **detached** |

Of the 12 pre-registered points: **1 survives, 4 are excluded by van Hove,
and 7 are detached** — roughly half the window's fragility comes from the
window-geometry effect (§10.3's δ/span argument, which alone would leave few
survivors even if every point were finite) and the other half from M=1's
own high detachment rate (§10.1: 58.3%, the highest of any M), which
independently removes the window's entire high-μ half from contention before
exclusion is even applied. Both causes are real; neither alone explains the
n=1 result.

**Raw (non-excluded) values, shown explicitly (already summarized as the
`median_s_raw` figure in §8.2, shown here as the individual points that
number is built from):** the 5 finite M=1 values are 0.000435, 0.000938,
0.003872, 0.354924, 0.841974 — median (all 5) = **0.003872**, versus the
pre-registered excluded statistic (n=1) = **0.000435**. Both numbers are
**≪ 1** and **≪** every M≥2 value reported (M=2's lowest surviving value is
1.159) — the qualitative M=1→M=2 jump is not sensitive to which of these two
M=1 numbers is used, but the specific M=1 magnitude should not be read as a
precise measurement.

**Methods-note caveat (added here, as the task requires, not merely noted in
this investigation's prose):** *the M=1 point in the pre-registered statistic
carries much less statistical weight than the M=2/M=3/M=4 points — it rests
on a single surviving sample out of 12 pre-registered points, for the
structural window-geometry and detachment-rate reasons given above. Its
near-zero value should be read as **consistent with** R9 Track 1's much more
thoroughly characterized single-channel baseline (s≈0.29 at the R6/R9-T1
operating point, §8.1's protocol-a cross-check), not as an independent,
precise re-confirmation of it from this track's own protocol-b window.*

### 10.4 Task 3 — is the M=3→M=4 dip real physics, a sampling artifact, or
###      genuinely unclear?

**Checking the van Hove exclusion step first, since it's the most obvious
candidate mechanism.** Comparing raw (pre-exclusion) medians to the
pre-registered (excluded) medians at every M shows the exclusion step moves
each number only modestly:

| M | median s (raw) | median s (excluded) | change from exclusion |
|---:|---:|---:|---:|
| 2 | 2.0799 | 1.9885 | −4.4% |
| 3 | 3.0601 | 2.8958 | −5.4% |
| 4 | 2.0579 | 2.0563 | −0.1% |

**The van Hove exclusion step is therefore NOT the primary driver of the
M=3-vs-M=4 gap** — both M=3 and M=4's raw and excluded medians are close to
each other, and the M=3>M=4 ordering is already present in the *raw* medians,
before any exclusion is applied at all. This rules out a simple "different
exclusion counts explain the dip" reading.

**What the raw point-count and coverage DO show:** M=3's finite (non-detached)
points span its *entire* window continuously except for the single detached
point at μ=183.11 — i.e. essentially full coverage from μ=148 (near threshold)
to μ=392 (window's far edge). M=4's finite points, by contrast, are the
*survivors of a 53.3% detachment rate* (§10.1) — they span only μ=313.90 to
515.70, **never reaching the upper ~18% of M=4's own window** (μ=515.70 to
627.81 is entirely detached, the persistent tail of §10.2). This is a genuine
coverage asymmetry between the two M's samples, conceptually similar to the
survivorship-bias check R9 Track 3 already ran on its own disorder ensemble:
M=4's median characterizes "s where a solution exists" over a **narrower
fraction of its nominal window** than M=3's median does, because a large part
of M=4's window has no self-consistent solution to measure at all. Whether s
would have continued declining, plateaued, or risen again in that unmeasured
region is genuinely unknown — no root exists there to evaluate, so this is
not a resolvable question with the current data, only a disclosable
disanalogy between the two comparisons.

**Consistency with the weak sum-rule argument (scoping §6.1).** That argument
is explicitly *weak*: total current sums over occupied channels with
*unequal, nonlinear, depth-dependent* weights, giving only an **envelope**
expectation of increase, not a point-wise guarantee — and critically, it was
stated for adding a channel at a **fixed** (μ, B) operating point, not for
comparing two *different* M-specific windows centered on two different
thresholds (M=3's window sits at μ≈148–392, M=4's at μ≈314–628 — different
(μ, B) neighborhoods, not the same point with one more channel added). Read
this way, the M=3→M=4 decline in the window-median statistic is **not a
violation of a hard prediction** — the sum-rule argument does not, in its own
stated form, forbid it.

**Verdict — genuinely mixed, as the task anticipated it might be, not forced
to a clean answer:** the M=3>M=4 ordering is not clean confirmed physics, and
it is not a simple counting/exclusion artifact either (the exclusion step
barely moves either number). The most defensible reading is that it reflects
**(i)** the weak sum-rule's explicit non-guarantee of point-wise monotonicity
across genuinely different M-specific (μ, B) windows, compounded by **(ii)**
a real coverage asymmetry — M=4's much higher detachment rate systematically
narrows the μ-range over which its median could even be computed, relative to
M=3 — that is disclosed here rather than resolved. Distinguishing how much of
the gap is attributable to (i) versus (ii) would require a matched-coverage
resampling analysis, which is out of this investigation's scope (it would
mean recomputing statistics over a restricted, common μ-range across M,
arguably a new derived statistic, not a re-tune of the existing one, but still
new analysis beyond what was asked here) — flagged as a natural next question,
not attempted. **What is unaffected by any of this:** every individual finite
s value at every M≥2, in both M=3 and M=4's full data, sits above the s=1
stability boundary (§8.2, §9) — the qualitative sub-critical→super-critical
transition at production coupling does not depend on resolving the M=3-vs-M=4
ordering question.

## 11. Second post-hoc investigation: why is the detachment rate U-shaped
##     (M=1, M=4 high; M=2, M=3 low), not monotonically increasing?

**Status of this section:** a further post-hoc investigation of §10.1's
detachment-rate finding (58.3% / 0% / 12.5% / 53.3% for M=1–4). Reads only
already-computed fields (`mu`, `B`, `phi_star`, `bare_gap`, all already stored
per point in `r9_multisubband_sweep_data.json`) plus the pre-registered window
geometry (`protocol_b_window`, `subband_offsets`, `subband_spacing` —
unchanged, not re-derived). No sweep was re-run, no window/δ/statistic was
altered, and no number already reported in §8–10 is changed by anything
below.

### 11.1 Where, within each window, do M=1's and M=4's detached points fall?

Expressing each point's μ as a fraction of its own window,
`frac = (μ − μ_min) / (μ_max − μ_min)` (0 = window start, 1 = window end):

**M=1** (window [−26.16, 78.48]): finite at frac = 0, 0.09, 0.18, 0.27, 0.45;
**one isolated detachment at frac = 0.36** (μ=11.89, between two finite
points); then **detached continuously from frac = 0.55 through frac = 1.00**
(μ = 30.91 … 78.48, six consecutive points, the window's entire remaining
upper half) — a **persistent tail that runs all the way to the window's
high-μ edge and never recovers**.

**M=4** (window [313.90, 627.81], post-densification): finite at frac = 0,
0.14, 0.21, 0.29, 0.43, 0.57, 0.64; **three isolated detachments scattered
through the interior** (frac = 0.07, 0.36, 0.50 — μ = 336.33, 426.01, 470.86,
each with a finite point immediately before *and* after it); then **detached
continuously from frac = 0.71 through frac = 1.00** (μ = 538.12 … 627.81, five
consecutive points) — the same **persistent-tail-to-the-edge** shape as M=1.

**Direct comparison:** both M=1 and M=4 show the *same qualitative shape* — a
mix of scattered isolated single-point failures (each recovered by the very
next point) plus a **persistent tail that begins at some onset point and runs
uninterrupted to the window's high-μ edge, never recovering**. This shared
shape is not coincidental-looking (a "point of no return" that never reverses
across 5–6 consecutive points, twice, independently) and does look like a real
common phenomenon. **But the onset location, normalized by each channel's own
threshold and spacing, is NOT the same between the two:**

| M | last finite point, (μ−thr)/span | first point of the *persistent* tail, (μ−thr)/span |
|---:|---:|---:|
| 1 | 0.273 | **0.394** |
| 4 | 0.524 | **0.619** |
| 2 (reference) | 1.000 (window end) | — no persistent tail at all |
| 3 (reference) | 1.000 (window end) | — no persistent tail at all (only 1 isolated point) |

M=1's persistent tail begins at only 39% of its own channel spacing past
threshold; M=4's begins at 62% — genuinely different fractions, not a shared
numerical constant. M=2 and M=3, by contrast, track their principal branch
successfully all the way to **100%** of their own window (their own full
`thr + span` upper bound) without ever hitting this kind of persistent
failure. **So M=1 and M=4 share a qualitative failure shape (persistent,
non-recovering tail to the window edge) that M=2/M=3 simply never exhibit at
all — but the specific onset location is not quantitatively identical between
M=1 and M=4**, which is evidence against (though does not fully rule out) a
single shared numerical trigger value.

### 11.2 Ruling out the two most obvious numerical explanations, with direct evidence

**Hypothesis: "large B exceeds some numerical range."** Checking the actual
`B` values recorded at every point: M=1's detached points span B ≈ 12–78;
M=4's span B ≈ 34–235. But **M=2 reaches B up to 130.75 and M=3 reaches B up
to 183.06 — both larger than M=1's entire range and comparable to most of
M=4's — with M=2 showing ZERO detachment and M=3 only one isolated point.**
Large B by itself is **not sufficient** to trigger detachment (M=2/M=3 are a
direct counterexample at even larger B), and is therefore **ruled out** as
the explanation, not merely unsupported.

**Hypothesis: "the fixed-width search bracket (from `_find_fixed_points_ms`
/ `_phi_star_ms`, ±1.5 or ±3 periods = ±0.118 / ±0.236 in phase units) can't
reach the true root once B/μ has pushed it far away" — the mechanism this
task named by analogy to R9 Track 3's search-centered-on-the-wrong-value bug.**
Checking `phi_star` (the self-consistent phase actually found) at the last
finite point immediately *before* each onset: M=1's φ* = 0.021 at μ=21.40
(only **8.9%** of the ±3-period search width used); M=4's φ* = −0.011 at
μ=515.70 (only **4.7%** of that width). **In both cases φ* is small and
comfortably inside the search range right up until the very next point, where
it vanishes entirely** (no gradual approach to the search boundary is
visible anywhere in the data). If the search window were too narrow for a
root that had drifted far away, φ* would be expected to trend toward the edge
of the search range before failing — it does not. **This specific
window-construction/search-range hypothesis is therefore also directly
contradicted by the data, not merely unconfirmed.**

**Checking against the established R4/R5/R7 "razor-thin dead-zone at a
gap-closing" precedent, too.** R4's own methods note documents a genuine,
previously-validated phenomenon: "exactly at the gap-closing the principal
self-phase branch detaches ... the sawtooth self-consistency map jumps across
the identity line without crossing it" (`r4_pfaffian_boundary_methods.md`,
Solver section and Finding 1) — a real, physical (not numerical) dead-zone,
but one specifically tied to the spectral gap closing. Checking `bare_gap` at
every point in both persistent tails: **M=1's tail sits at bare_gap ≈ 0.049,
completely flat**; **M=4's tail declines smoothly from 0.209 to 0.139 with no
discontinuity or minimum** — **no spectral gap closure accompanies either
persistent tail.** This means the M=1/M=4 phenomenon is **not simply the same
mechanism as R4/R5/R7's own precedent**, which requires an actual gap closing
— it is evidently a distinct way for the principal-branch fixed point to stop
existing, one this track's own code has not previously characterized.

### 11.3 Task 2 — are M=2/M=3 the unremarkable baseline?

Yes, and this reframes the question correctly, as the task anticipated. R4's
own established result is that the self-consistency map is **attracting at
the generic transition field for every N tried**, with genuine detachment
confined to a narrow window specifically at a gap-closing — i.e., **low/zero
detachment away from a gap-closing is the documented, expected default**, not
something that itself needs explaining. M=2 (0/10 detached) and M=3 (1/8,
a single isolated point) are fully consistent with that established baseline.
**The real, correctly-framed question is not "why is the rate U-shaped" but
"why do M=1 and M=4 specifically deviate from the otherwise-expected clean
baseline that M=2 and M=3 both confirm."**

### 11.4 Task 3 — verdict

Evaluating the three options directly:

- **(a) A shared numerical/structural cause, unrelated to physics** — the
  *specific* mechanism this task named (window-construction/search-range
  fragility, by analogy to R9 Track 3's search-centering bug) is **directly
  ruled out** by §11.2: φ* stays small and well inside the search range at
  every point immediately preceding a detachment onset, with no trend toward
  the boundary. A large-B explanation is also ruled out by direct
  counterexample (M=2/M=3 reach comparable or larger B with no such failure).
  **(a), in the specific form proposed, is not supported by the evidence.**
- **(b) Two different, independently-explicable causes — one for M=1, one
  for M=4 (genuine physics for M=4, a distinct small-window edge effect for
  M=1)** — partially supported: M=1's situation *does* have a documented,
  M=1-unique structural feature (§10.3: μ<0 there means the single channel is
  literally unoccupied, and δ/span=28.4% is the highest of any M) that has no
  exact analogue at M=4. But the *persistent-tail* shape itself (§11.1) is
  shared in form between M=1 and M=4 in a way that looks like more than
  coincidence, and neither M=1's nor M=4's tail is explained by a confirmed
  mechanism (§11.2 rules out the two most obvious ones for both).
- **(c) Genuinely inconclusive given the available sample sizes** — **this is
  the best-supported reading, with real content, not a shrug:** the two most
  obvious candidate explanations (large-B, search-range) are **actively
  ruled out**, not merely unconfirmed, which is itself a genuine finding. What
  remains is a real, shared, non-coincidental-looking *shape* (a
  never-recovering tail to the window's high-μ edge, occurring only for M=1
  and M=4) whose underlying cause is not identified by the data or code paths
  available here, and whose onset location is quantitatively different
  between the two (39% vs 62% of the channel's own spacing past threshold —
  §11.1), arguing against (though not fully excluding) a single shared
  numerical constant as the trigger.

**Conclusion: (c), with two specific hypotheses actively excluded rather than
merely unresolved.** The detachment-RATE finding should be reported as an
**open, disclosed phenomenon distinct from and additional to** the
gap-dip-associated detachments already characterized in §10.2 — not
attributed to a search-window bug (ruled out), not attributed to the
established R4/R5/R7 gap-closing dead-zone (ruled out via the flat/smooth
bare_gap), and not confidently attributed to M-dependent physics either
(no confirmed mechanism, and the two onset fractions differ). **What additional
data would resolve this:** a denser μ-scan specifically bracketing each
persistent-tail onset (e.g. at M=1 between μ=21.40 and μ=30.91, and at M=4
between μ=515.70 and μ=538.12) would show whether the transition is sharp or
gradual and whether it coincides with any other computable quantity (e.g. a
higher-harmonic feature in the current-phase relation not captured by
`bare_gap`, which only reports the minimum |eigenvalue|, not the shape of
I(φ) itself) — flagged as a natural, scoped follow-up, not attempted here per
this task's explicit instruction not to re-run the protocol.

**Unaffected by any of this:** the pre-registered s-vs-M result (§8–9) is
built entirely from the points that DID converge, and does not depend on
explaining why the others did not — **the s-crosses-1-at-M≥2 finding remains
this track's real headline, unchanged.** Because this section's conclusion is
diagnosis, not remediation, no fix is proposed or applied here; whether it is
worth a future, separately-scoped investigation (with new instrumentation to
distinguish a true sign-non-change from a discontinuity-rejected sawtooth
jump, per §11.2's precedent check) is a decision for a future task, not this
one.
