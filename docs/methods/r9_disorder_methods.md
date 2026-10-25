# R9 — Track 3: On-site disorder in the BdG / Pfaffian / self-consistent-ring framework

**Module:** `src/bdg_r9_disorder.py` (additive; imports R1/R4/R5/R7 and
`bdg_parallel` **read-only**, modifies nothing in R1–R8 or R9 Tracks 1–2).
**Outputs:** `ISEF R9 data folder/` (`r9_disorder_data.json`,
`r9_disorder_summary.png`).
**Run:** `py -3.13 src/bdg_r9_disorder.py --full` (`--validate` for the W=0
gate only, `--benchmark` for the parallel/BLAS benchmark + pre-run time
estimate, `--quickfull` for a tiny-grid smoke test).

Every prior phase (R1–R8, R9 Tracks 1–2) used a **perfectly
translationally-ordered** lattice. Real InAs/Al nanowires are disordered
(impurities, interface roughness, charge inhomogeneity), and this is a
directly literature-motivated concern: disorder-induced trivial sub-gap
states have been proposed in the Majorana-nanowire literature as an
explanation for some contested zero-bias-peak claims (see Task 3 for the
verified primary source). This track adds on-site disorder as a strictly
additive layer and re-tests the project's key null results — the boundary-
shift null (R4/R5/R7), the SO-correction finding (R8-C), and Z_r redundancy
(R8) — for robustness against the one ingredient every real device has that
this model lacked until now.

---

## Task 1 — the disorder primitive, and why it preserves particle-hole symmetry

On-site disorder is a random potential `eps_i ~ Uniform[-W/2, W/2]`, i.i.d.
over the N sites, added to the chemical-potential term. Every clean builder
in this repo (`build_bdg_lo`, `_build_ring_lo_phased`) uses the same
`kron(spin, site)` Nambu ordering — electrons first 2N, holes last 2N, site
the fast index — so disorder is injected as

```
V(eps) = block_diag( kron(I2, diag(eps)),  -kron(I2, diag(eps)) )   (4N x 4N, REAL)
H_dis  = H_clean + V(eps)
```

**PH-symmetry proof.** `check_ph_symmetry_complex` requires
`Sigma_x @ H.conj() @ Sigma_x == -H`, where `Sigma_x = [[0,I_2N],[I_2N,0]]`
is the electron↔hole block swap. For `V = [[A,0],[0,-A]]`
(`A = kron(I2,diag(eps))`), a direct 2×2 block computation gives
`Sigma_x @ V @ Sigma_x = [[-A,0],[0,A]] = -V` exactly, and since `eps` is
real, `V.conj() = V`, so `Sigma_x @ V.conj() @ Sigma_x = -V` for **any**
realization. Adding `V` to a clean `H` that already satisfies the PH
identity therefore preserves it exactly, for any `eps` — verified
numerically (not just algebraically) via `check_ph_symmetry_complex` on a
genuine W=5Δ sample in `validate_W0()` (result: ε_PH = 0.0, machine zero).

**Downstream shells.** Because every clean detector (the Pfaffian invariant,
the ring boundary bisection, the open-wire sub-gap detector, the
self-consistent scalar map) builds its own H internally, disorder-aware
"shells" (`_disordered_ring_pfaffian_sign`, `_disordered_ring_boundary_Bc`,
`disordered_openwire_subgap`, `_disordered_scalar_ring_map` and its root-
finding wrappers, `_disordered_Bc_selfconsistent`) call the identical clean
builder, add `V(eps)`, and run the *identical* downstream algorithm. At
`eps=0` (W=0), `add_onsite_disorder` is an exact no-op, so every disordered
shell reproduces its clean counterpart **exactly**, not approximately.

### The one approximation, flagged explicitly

R5's self-consistent solver reduces the ring self-field problem to a
**scalar map on the uniform manifold**, and that reduction is *exact* there
because a clean, translation-invariant ring guarantees uniform-in →
uniform-out (verified by R5's `find_uniform_fixed_points`/`ring_mode_gains`
machinery). Disorder breaks translational invariance, so a uniform input
texture no longer produces a perfectly uniform current/self-field response
in general. `_disordered_scalar_ring_map` evaluates the full per-bond
response and projects it back onto the uniform manifold via `np.mean` — a
controlled mean-field-style approximation, **not** an exact reduction.
`validate_W0` confirms the projection is exact at W=0 (where the response
*is* uniform); `_disordered_scalar_ring_map`'s optional `diagnostics` records
the realized per-bond non-uniformity so the approximation's size is visible
rather than hidden. Solving the full N-dimensional non-uniform
self-consistency per realization (Newton-Krylov, as R5's
`nonuniform_fixed_point_search` does for the clean diagnostic) would
multiply the cost of the already most-expensive R9 track by another
N-dimensional root solve per realization × grid point × W — explicitly out
of scope, and the stated limitation for Task 4.3.

---

## Task 1.3 — W=0 exact-reproduction validation gate (passed before any W>0 result)

`validate_W0(N=80, mu=1.2)`:

| check | quantity | result |
|---|---|---|
| Pfaffian invariant | max\|Δsign\| over a 9-point B-grid straddling B_c | **0.0** (exact) |
| ring boundary B_c | `_disordered_ring_boundary_Bc(W=0)` vs `ring_boundary_Bc_pfaffian` | **0.0** (bit-identical) |
| self-consistent B_c | `_disordered_Bc_selfconsistent(W=0)` vs `ring_boundary_Bc_pfaffian_selfconsistent` | **0.0** (bit-identical) |
| open-wire sub-gap detector | min\|E\| at W=0 vs `build_bdg_lo` spectrum | **0.0** (bit-identical) |
| PH symmetry at genuine W=5Δ | `check_ph_symmetry_complex` | **0.0** (machine zero) |
| **ALL_PASS** | | **True** |

`run_full()` calls this first and raises before computing any W>0 quantity
if it fails — the mandatory gate per the task instructions.

---

## Task 2 — parallelization, seeding, and the pre-run time estimate

### Seeding scheme (explicit, reproducible)

Each disorder realization gets an independent seed via
`np.random.SeedSequence(base_seed).spawn(n_real)`, reconstructed in the
worker as `np.random.default_rng(seeds[k])`. **No realization ever touches
the global `np.random` state.** Convergence checks reuse the *same*
cumulative seed list (`seeds_all[:n]` for increasing `n`) so smaller
realization counts are strict prefixes of larger ones — a true
nested-ensemble comparison, not independent resampling. This is the
numpy-recommended multiprocessing-safe scheme: seeds are passed as data
(picklable `SeedSequence` objects), so worker placement under
`ProcessPoolExecutor` never changes which draw a given realization gets.

### Parallel infrastructure: reused, not rewritten

All dispatch goes through `bdg_parallel.parallel_map`/`check_equivalence`
(R8's infra), unchanged. Two **structurally different engines** live in
this track and were benchmarked separately, since they have very different
per-task costs:

- **cheap engine** (topological classification / open-wire sub-gap
  detector): one 4N×4N `eigh`/`eigvalsh` (+ a Pfaffian for the topo
  classifier) per realization.
- **expensive engine** (disordered self-consistent B_c, Task 4.3): a full
  R5-style root search (many `eigh` calls) at each of several B points
  during a Pfaffian bisection, repeated for the disordered map — this is
  the compute-dominant part of the whole track, exactly as flagged in the
  task instructions.

**Serial-vs-parallel equivalence** (bit-identical, both engines, both a
1-worker and an N-worker pool with 1-thread BLAS): `max_abs_diff = 0.0`,
`max_rel_diff = 0.0` for both the cheap and expensive engines.

**BLAS-oversubscription benchmark** (measured, N=80):

| engine | serial, 1-thread BLAS | serial, default BLAS | parallel (Nw workers, 1-thread) |
|---|---:|---:|---:|
| cheap (subgap) | ~0.05 s/task | ~0.05 s/task (no benefit — matrix too small) | **worse**: ~0.20–0.40 s/task (pool overhead dominates) |
| expensive (self-consistent) | ~144–161 s/task | ~138–160 s/task (no benefit) | ~35–65 s/task (**~2.3–4×** speedup at 2–5 workers) |

**Conclusions, and how the code acts on them:**
1. Default multithreaded BLAS gives **no benefit** on these matrix sizes
   (4N=320) for either engine — consistent with every prior R7–R9 track.
2. For the **cheap** engine, process-pool dispatch is *slower* than plain
   serial (spawn/IPC overhead exceeds the ~0.05s task cost) — so
   `boundary_smearing_vs_W` (Task 4.1) and `subgap_states_vs_W` (Task 4.2)
   **always dispatch serially**, regardless of `n_workers`, matching the
   "small/cheap loops stay serial" rule already documented for
   `bdg_parallel` on other tracks.
3. For the **expensive** engine, parallelization gives a genuine ~2.3–4×
   speedup (worker count fluctuates 2–5 depending on free RAM at benchmark
   time — this laptop often has <2–3 GB free), so `null_result_vs_W`
   (Task 4.3) dispatches through `bdg_parallel.parallel_map` with 1-thread
   BLAS per worker.

### A measured cost surprise, and the fix

The first benchmark run showed `_scshift_worker` costing ~130s/task —
**twice** the ~65s expected — because it recomputed the ~28-30s **clean**
self-consistent baseline on *every single ensemble member*, even though
that baseline is W-independent. Fixed with a per-`(N,mu,g,delta)` process-
local cache (`_sc_baseline`), roughly halving the effective cost. A second,
more important discovery came from the `--quickfull` smoke test (real
end-to-end run, not a synthetic benchmark): **a disorder realization that
fails to converge costs *more* than one that succeeds**, because the
disordered root-finder's window-widening fallback repeats a full
coarse-B-scan × root-search on every retry, and realistic-W disorder makes
non-convergence common. This is exactly the kind of cost the isolated
per-task benchmark could not see. In response, `_disordered_Bc_selfconsistent`'s
`_max_widenings` was reduced from the clean solver's 4 to **2** — a
realization that is going to fail should fail fast, not burn 2× the
wall-clock on retries unlikely to help; the failure itself remains a
reported, counted outcome (`convergence_failed`), not silently dropped.

### Honest pre-run time estimate (surfaced before the full Task 4 sweep)

Printed by `--full` immediately after the benchmark, BEFORE any Task 4
experiment ran: **Task 4.1+4.2 ~1.5 min, Task 4.3 ~19.1 min, TOTAL ~20.6
min.** The realized wall time was higher (**77.1 min total** for the full
run, including the ~19-min benchmark itself) — the gap between the estimate
and reality is itself a finding, and is exactly the convergence-failure
effect described above: the estimate assumed every task costs the
benchmark's measured *average* per-task time, but at realistic W most
Task 4.3 realizations fail, and failing realizations are not free even with
`_max_widenings=2`. This is reported here rather than smoothed over,
consistent with the project's standard of flagging when a measured cost
exceeds a stated pre-run estimate (cf. R7's thermal-throttle disclosure).

### Scope reduction (measured, reported explicitly — not silent)

Task 4.3 (the expensive engine) is, as flagged in the task instructions, the
compute-dominant experiment of this track. Sweeping the *full* W_weak ∪
W_realistic grid at `n_real=30` (Task 1's upper guideline) was estimated at
~2.7 hours from the raw per-task benchmark number alone, **and the
quickfull smoke test showed the real cost running higher still** once
convergence-failure retries are accounted for. The delivered Task 4.3 run
therefore narrows the W grid to **one representative weak point (W=1Δ) plus
the full realistic band** (W_weak's other three points are already covered
cheaply by Tasks 4.1/4.2). The first pass used `n_real=20` (Task 1's stated
floor); a **post-hoc review found this run's near-total convergence
failure was a search-window bug, not genuine physics** (see §4.3 below) —
after fixing it, `n_real` was raised to **30** (still short of Tasks
4.1/4.2's 50, an explicitly stated, compute-driven limitation, since the
corrected solver's per-task cost did not drop after the fix). This is a
deliberate, reported scope reduction driven by measured compute cost, in
the same spirit as R7's thermal-throttle compute ceiling — not a silent
narrowing.

---

## Task 3 — a literature-verified realistic disorder scale

**Primary source (independently re-fetched and quoted verbatim — this
project caught a fabricated citation before, so nothing here is asserted
from memory):** Seongjin Ahn, Haining Pan, Benjamin Woods, Tudor D.
Stanescu, Sankar Das Sarma, *"Estimating disorder and its adverse effects in
semiconductor Majorana nanowires,"* Phys. Rev. Materials **5**, 124602
(2021), arXiv:2109.00007.

Two independent fetches of the paper's full text returned the identical
verbatim sentence: *"Onsite random potential drawn from an uncorrelated
Gaussian distribution with zero mean and standard deviation σ_μ = 2.5 meV"*
on a tight-binding lattice with constant **a = 3.33 nm**. The paper
separately estimates transport-based broadening Γ = ℏ/2τ ≈ **2–3 meV**
from three labs' InAs/InSb mobility/impurity-density data (Eindhoven,
Copenhagen, Purdue), and states its headline conclusion in the abstract:
*"the current level of disorder in semiconductor Majorana nanowires is at
least an order of magnitude higher than that necessary for the emergence of
topological Majorana zero modes."*

### Converting σ_μ to our units — the softest number in this track

Our model uses **Δ = 0.19 meV** (Chang et al. 2015, R6-verified) and
**Uniform[-W/2,W/2]** disorder (per the task specification), not the
literature's Gaussian σ_μ. Both honest conversions are reported (never
collapsed to one number, since the conversion factor is itself
model-dependent — flagged as the track's softest number):

| convention | formula | W/Δ |
|---|---|---|
| "same energy scale" (loose) | `W = σ_μ` | **13.16** |
| "variance-matched uniform" | `W = σ_μ·√12` | **45.58** |
| broadening-based cross-check | `Γ = 2–3 meV` directly | **10.5 – 15.8** |

All three land in the same qualitative place: **W_realistic ≈ 10–46 Δ**,
i.e. disorder comparable to or an order of magnitude *larger* than the gap
itself — this is not an arbitrary large-W choice, it is a direct reflection
of Ahn et al.'s headline finding that real devices are already
over-disordered relative to what topological MZMs need.

**Sweep grid used (Task 3.2), realistic vs exploratory clearly separated:**

- **W_weak** (exploratory-small, well below the gap): 0.1, 0.5, 1.0, 2.0 Δ
- **W_realistic** (literature-grounded band above): 10, 15, 25, 46 Δ
- **W_exploratory** (beyond the literature band, to see the trend continue): 60, 90 Δ

(Lattice-constant caveat: Ahn et al.'s tight-binding discretization uses
a=3.33 nm vs this project's a=10 nm throughout R1–R9; the σ_μ energy scale
is compared directly without a further discretization correction, which is
part of why both conversions above are reported rather than one.)

---

## Task 4 — does the topological classification survive disorder?

Full run: N=80, μ=1.2, g=G_R7=2.9×10⁴ (the same subcritical feedback used
throughout R4–R7–R9-Track-2), B fixed at `sqrt(mu²+Δ²) - 0.05` — a point on
the **clean-trivial** side, just below the clean boundary (`ring_boundary_
Bc_pfaffian(N=80,mu=1.2) ≈ 1.560Δ`). n_real=50 for Tasks 4.1/4.2 (both
converged — see below); n_real=30 for Task 4.3 (compute-bound, corrected
after a post-hoc review — see the scope reduction above and §4.3 below).
Every number below carries its realization count and spread; raw
distributions are in `r9_disorder_data.json`.

### 4.1 Boundary smearing — **smooth, and real: disorder measurably erodes the boundary**

Fraction of realizations misclassified as topological (relative to
`classify_topological`'s clean-referenced sign, see Task 1) at this
clean-trivial point, n=50 per W:

| W (Δ) | frac. topological | Bernoulli std |
|---:|---:|---:|
| 0.1 | 0.000 | 0.000 |
| 0.5 | 0.000 | 0.000 |
| 1.0 | 0.040 | 0.028 |
| 2.0 | 0.080 | 0.038 |
| **10 (realistic)** | **0.080** | 0.038 |
| **15 (realistic)** | **0.060** | 0.034 |
| **25 (realistic)** | **0.180** | 0.054 |
| **46 (realistic)** | **0.320** | 0.066 |
| 60 (exploratory) | 0.300 | 0.065 |
| 90 (exploratory) | 0.300 | 0.065 |

**Convergence-in-realization-count** (at W=10, the lowest realistic point,
schedule 20→100): frac. = 0.050→0.067→0.075→**0.080**→0.080→0.070, spread
shrinking 0.049→0.026; stabilized (|Δmean| < spread) from n=30 onward —
the ensemble is well-converged at n=50.

**Reading:** the trend is **smooth and monotonic-with-noise**, not sharp
and not erratic — no evidence of a single critical `W_c` where the
boundary abruptly collapses; instead the misclassified fraction climbs
gradually through the realistic band and **plateaus around 0.30** in the
exploratory regime (46–90Δ), i.e. disorder does not keep eroding the
classification without bound — it saturates. But the saturated value is
**not small**: at the top of the literature-realistic band (W≈46Δ,
the variance-matched-uniform conversion of Ahn et al.'s σ_μ), **very
close to 1 in 3 realizations of a clean-trivial point are misclassified as
topological.** This is a genuine, reportable disorder-driven smearing of
the phase boundary, consistent with the literature's warning that real
InAs/Al disorder is large enough to matter for topological classification.

### 4.2 Spurious sub-gap states — **no mean increase, but growing realization-to-realization spread**

Open-wire spectral sub-gap count (threshold 0.25Δ), same (N,μ,B), n=50:

| W (Δ) | ⟨n_subgap⟩ | std |
|---:|---:|---:|
| 0.1 | 2.00 | 0.00 |
| 0.5 | 2.20 | 0.61 |
| 1.0 | 2.40 | 0.81 |
| 2.0 | 2.64 | 0.94 |
| **10 (realistic)** | 2.12 | 0.94 |
| **15 (realistic)** | 2.16 | 1.33 |
| **25 (realistic)** | 1.84 | 1.27 |
| **46 (realistic)** | 1.88 | 1.30 |
| 60 (exploratory) | 1.60 | 1.46 |
| 90 (exploratory) | 1.24 | 1.13 |

**Convergence** (W=10, schedule 20→100): mean 2.50→2.40→2.20→**2.12**→2.11→2.14,
spread shrinking 0.20→0.10; stabilized from n=30 onward.

**Reading — this is the R4 threshold-artifact connection, investigated
explicitly.** The ensemble-**mean** sub-gap count does **not** grow with
disorder — if anything it drifts mildly down at large W (2.6→1.2). So this
run does **not** find that disorder systematically ADDS more spurious
sub-gap states on average, which would have been the naive expectation.
What it *does* find is that the **spread** across realizations grows
substantially with W (std 0.00→1.30, comparable to the mean itself at
large W) — i.e. a fixed disorder strength produces a wide RANGE of sub-gap
spectra from one realization to the next: some disordered wires show many
more sub-gap states than the clean wire, others show none extra. This
matches the *qualitative* phenomenology Ahn et al. (2021) describe for real
devices — "most spectra show no significant features... consistent with
the sporadic nature of experimental observations" — a physical mechanism
for occasional, not universal, disorder-induced sub-gap structure that
could mimic R4's finite-size sub-gap artifact on a subset of realizations,
even though the ensemble mean stays flat.

### 4.3 The R4/R5/R7 inductive-feedback null at realistic W

**⚠ SUPERSEDED FINDING, CORRECTED BELOW.** The first version of this
section reported that the self-consistent solver "mostly fails to
converge" at realistic W (n_real=1-4 out of 20 attempted at every
realistic point) and read this as itself the headline finding. A **post-hoc
review caught that this was wrong**: the reported n=1-4 samples were not a
genuine ensemble average, and — more importantly — the near-total
convergence collapse was traced to a **bug in this track's own code**, not
a physical instability. The original (superseded) table and reasoning are
preserved verbatim in `r9_disorder_data.json`'s
`task4_3_null_result_SUPERSEDED` key, per this project's standard of
recording what changed and why (cf. the `pfaffian_trim_lo` fix, the
Anderson seed-count correction). The root-cause investigation and the
corrected result follow.

#### Root-cause investigation (Task 1 of the post-hoc review)

`_disordered_Bc_selfconsistent`'s coarse B-scan searched in a window
centered on the fixed **clean** analytic value `ana = sqrt(mu²+Δ²)`,
inherited unchanged from the clean R5 solver — correct there *only* because
the clean rigid boundary coincides with `ana` exactly. The disordered
version never re-examined that assumption. A direct diagnostic reproduced
the EXACT 20 seeds `null_result_vs_W` uses (`base_seed=2028`) and, for each
one, computed the wide, Task-4.1-trusted rigid (non-self-consistent)
boundary — the same one that reliably found a boundary in 100/100
realizations across every W tested — and checked whether it fell inside
the old solver's widest-ever search half-width (±0.20, from `window=0.05`
doubled twice):

| W (Δ) | rigid boundary INSIDE ±0.20 | mean displacement from `ana` | max displacement |
|---:|---:|---:|---:|
| 1 | 20/20 | +0.010 | +0.063 |
| 10 | 2/20 | +0.445 | +1.218 |
| 15 | 3/20 | +0.607 | +1.098 |
| 25 | 3/20 | +0.952 | +2.870 |
| 46 | 5/20 | +0.999 | +2.803 |

This is decisive: at realistic W, disorder itself (a finding **already
established independently by Task 4.1**) displaces the rigid boundary by an
amount that dwarfs the solver's search window — 15-18 of 20 boundaries per
W point fall completely outside where the solver was ever looking. The
solver wasn't finding "no self-consistent fixed point exists" — it was
looking in the wrong place. This is **Option A** (a fixable search bug),
not Option B (an intrinsic self-consistency-loop instability): the
self-consistency question (does a fixed point exist near *this
realization's actual* boundary?) was never actually being asked at
realistic W.

#### The fix

`_disordered_Bc_selfconsistent` now centers its search on the
**per-realization rigid boundary** (`_disordered_ring_boundary_Bc`, the
same cheap, 100%-reliable detector Task 4.1 already depends on) instead of
the fixed clean `ana`. At W=0 (`eps` all zero) the function fast-paths to
the exact analytic `ana` rather than the rigid detector's `tol=1e-4`
bisection value — this is what keeps the mandatory W=0 gate **bit-identical
to machine precision** (`selfconsistent_absdiff = 0.0`, re-verified after
the fix), matching `add_onsite_disorder`'s own eps=0 exact-no-op
convention, rather than merely close to it. `_max_widenings` was also
bumped from 2 to 3 (still below the clean solver's 4) after a small
post-fix sample showed a residual few-percent-scale rate of genuinely
larger self-consistent corrections benefiting from one more widening.

A small post-fix sample (n=6 per W, same seeds, before committing to the
full corrected re-run) confirmed the fix works and surfaced a second,
smaller, genuine effect: success rate rose from ~5-20% to **33-100%** across
W=1-46Δ, and among successes the shifts are large and sign-varying (e.g.
+14% to +144% at W=46Δ) — a **real, much stronger non-null signal** than
the original (buggy) run's small n=1-4 samples suggested. A **residual,
genuine `ScfConvergenceError` rate remains** even with correct centering
(the self-consistency root find itself sometimes fails at a specific probed
B) — this is reported honestly as a real, smaller, secondary effect, not
eliminated by the fix.

#### Corrected result

Ensemble self-consistent boundary shift vs the clean self-consistent
baseline (`ring_boundary_Bc_pfaffian_selfconsistent(N=80,μ=1.2,g=G_R7) =
1.559951Δ`), N=80, g=G_R7, **n_real=30** attempted per W (5 at W=0,
deterministic) — raised from the original n_real=20 given the corrected
solver's success rate now supports a real ensemble, though still short of
Tasks 4.1/4.2's n=50 (a stated, compute-driven limitation: the corrected
solver's per-task cost did **not** drop after the fix, since successful
solves cost more than the fast failures they replaced — re-benchmarked at
87.1s/task serial, 28.7s/task parallel with 5 workers, essentially
unchanged from the original benchmark). **Re-run wall time: 4023s (67.1
min)** for the full corrected sweep (6 W points + the convergence check),
close to the fresh pre-run estimate.

| W (Δ) | mean shift | std | SEM | n converged (of 30 attempted) |
|---:|---:|---:|---:|---:|
| 0.0 | +0.0000 % | 0.000 | 0.000 | 5/5 |
| 1.0 | +0.208 % | 1.091 | 0.265 | 17/30 |
| **10 (realistic)** | **+29.398 %** | 18.075 | 3.769 | **23/30** |
| **15 (realistic)** | **+44.549 %** | 25.050 | 5.747 | **19/30** |
| **25 (realistic)** | **+50.335 %** | 52.777 | 12.800 | **17/30** |
| **46 (realistic)** | **+56.514 %** | 69.164 | 15.093 | **21/30** |

**Convergence-in-realization-count** (W=10, schedule 10→30): mean
25.26→27.69→30.46→31.30→**29.40**, SEM 2.58→4.73→4.64→4.09→3.77; stabilized
(|Δmean| < spread) from n=15 onward — the corrected ensemble **is**
demonstrably converged at n=30, unlike the original (which could not even
report a defined mean before n=20).

**Reading — the corrected, properly-powered non-null finding.** With the
search-window bug fixed, realistic-W success rates are now **57-77%**
(17-23 of 30, up from the original 5-20%), giving genuine ensemble
averages, not single-sample artifacts. The corrected result is **smooth,
monotonic, and unambiguously non-null**: the mean self-consistent shift is
**positive and grows systematically with W** — from a small +0.2%
(W=1Δ, comparable in scale to the R4-R5-R7 clean-model null of
~0.06-0.13%) up to **+29% (W=10Δ) → +45% (W=15Δ) → +50% (W=25Δ) → +57%
(W=46Δ)** — i.e. two to three orders of magnitude larger than the clean
baseline throughout the literature-realistic band, with the rate of growth
itself slowing between W=25 and W=46 (+50%→+57%, a much smaller relative
step than the earlier ones), suggestive of the trend beginning to
saturate rather than diverging without bound. The spread also grows
substantially with W (std 1.1%→18%→25%→53%→69%) and a genuine minority of
realizations show a **reversed (negative) sign** even at high W (e.g. 3 of
17 at W=25Δ, several of 21 at W=46Δ) — so while the *ensemble mean* is a
robust, large, positive, growing effect, any *individual* disordered
device's shift is only loosely constrained by that mean.

**⚠ This N=80-only picture is revised below.** Everything in this
subsection is at a single system size, N=80 — the same finite-size
scrutiny already applied to every other quantitative boundary-shift result
in this project (R4, R5, R7, R8-C) had not yet been applied here. The
**N-independence check further down found that this N=80 result
substantially over-states the effect, especially at high W** (the growth
with W above is itself partly a finite-size artifact concentrated at
N=80) — read the "N-independence check" subsection below before treating
the +29%→+57% growth trend as a settled bulk result.

**A residual, genuine self-consistency-loop failure remains** (23-43% of
attempts across W, without a clear trend in W) — this is the smaller,
secondary effect the root-cause investigation predicted: even correctly
centered, the self-consistency root find (`phi_star`) itself sometimes
fails to converge at a specific probed B under disorder. **A plausible
(but not independently confirmed here) physical mechanism** ties this
residual failure to the *same* large, growing mean shift: R9 Track 1
established that this system has a genuine dynamical bifurcation (g_crit)
where the self-consistent map's local slope crosses 1 and no attracting
fixed point exists. Disorder plausibly concentrates current density
unevenly across bonds (breaking the clean model's uniform-response
guarantee — see Task 1's flagged projection approximation), which could
push the *local* effective feedback gain at some realizations past their
own effective g_crit even though the *nominal* g=G_R7 is unchanged and
subcritical in the clean system — consistent with both the large, growing
mean shift (realizations pushed closer to the bifurcation respond more
strongly) and the residual genuine convergence failures (realizations
pushed past it have no fixed point at all, exactly like R5's original
g_crit finding). This connects Track 3's finding directly to R9 Track 1's
pole characterization, but confirming it precisely would require
per-realization local mode-gain diagnostics, which is out of scope here
and flagged as an open question, not asserted as proven.

#### N-independence check — does the corrected shift survive finite-size scrutiny?

**A second post-hoc review asked the obvious next question**: every other
quantitative boundary-shift result in this project (R4's threshold
overestimate, R5's Aharonov–Bohm decay, R7's L-independence, R8-C's SO
correction) was put through explicit N-independence scrutiny and found to
shrink toward zero in the bulk (N→∞) limit. The corrected Task 4.3 result
above had not yet been checked this way — data JSON: `r9_disorder_N_
independence.json`.

**N-range feasibility (checked directly, not assumed).** At (μ=1.2,
g=G_R7), N=60 is a **genuine regime-2 dead-zone** — it fails even R7's own
5-path retried solver (`ring_Bc_selfconsistent_a`), matching R7's own
established finding that this dead-zone is N-dependent, not monotonic in
N. N=90 and N=100 fail the plain (unmodified) default-path solver too;
N=100 is recoverable with R7's retry logic (not a genuine dead-zone), but
modifying `_sc_baseline` to add retry logic was out of scope (the
already-fixed core solver was not to be touched), so it was skipped rather
than the pipeline being changed. N=110, 120, 130, 140, 160 all give a clean
baseline on the first (unmodified) attempt. Given measured ~N³ cost growth
(a single clean baseline at N=160 cost ~395s vs ~30s at N=80), the range
was narrowed to **N ∈ {120, 160}** in addition to the existing N=80 — two
additional points, not the suggested 2–3, an explicitly stated
compute-driven reduction.

**Result** (n_real=30 per (N,W), same base_seed=2028 seeding scheme, W=10
and W=46 bracketing the realistic band):

| N | W=10Δ shift | success (of 30) | σ from zero | W=46Δ shift | success (of 30) | σ from zero |
|---:|---:|---:|---:|---:|---:|---:|
| 80 | +29.40 % ± 3.77 | 23 (77%) | **7.80** | +56.51 % ± 15.09 | 21 (70%) | **3.74** |
| 120 | +19.84 % ± 3.36 | 29 (97%) | **5.90** | +9.56 % ± 8.24 | 27 (90%) | 1.16 |
| 160 | +22.55 % ± 3.67 | 30 (100%) | **6.15** | +10.85 % ± 7.21 | 30 (100%) | 1.50 |

(A convergence-in-realization-count schedule (10→30, an independent draw)
at the most extreme combination, N=160/W=46, confirms this specific point
is genuinely ensemble-converged, not a fluke: mean 10.20→5.86→12.76→
11.66→**10.85**, stabilizing from n=15 onward, consistent with the main
sweep's own +10.85%. A first attempt at this full check was lost when the
harness session was torn down mid-run — confirmed via zero residual
`python.exe` processes, not a crash — and the 3 already-completed points
were recovered verbatim from the surviving log rather than re-paying
1100–2900s of already-sunk compute per point; only the missing N=160/W=46
combination was actually recomputed.)

**Reading — a genuinely mixed, two-part verdict, not a clean yes/no:**

1. **Every reported number substantially over-stated the bulk effect at
   N=80**, especially at high W: the W=46Δ point drops from +56.5% (N=80)
   to +9.6% (N=120) — an **83% reduction** — then holds flat at +10.8%
   (N=160). The W=10Δ point drops from +29.4% to +19.8% (N=120) then rises
   slightly to +22.6% (N=160) — non-monotonic at the point-estimate level,
   but the three values overlap within ~1 SEM of each other, consistent
   with a real plateau around ~20-23% rather than a clean power law. **The
   residual convergence-failure rate is itself a finite-size effect**: it
   vanishes entirely by N=160 (100% success at both W points, up from
   70-77% at N=80) — meaning the "genuine self-consistency-loop
   instability" flagged as a secondary effect above is, at least in this N
   range, a **finite-size instability**, not an intrinsic feature of the
   bulk disordered system.
2. **Unlike R4/R5/R7/R8-C, the shift does *not* vanish toward zero** — it
   drops sharply from N=80 then plateaus at a smaller but still non-zero
   value. But **the statistical CONFIDENCE in that non-zero value depends
   strongly on W**, and this is the most important nuance the raw means
   alone would hide: at **W=10Δ the effect stays highly significant at
   every N tested (5.9–7.8σ from zero)** — a real, non-vanishing, robustly
   supported disorder effect, just smaller than the (inflated) N=80 number
   suggested. At **W=46Δ the significance collapses from 3.74σ (N=80) to
   only 1.16–1.50σ (N=120, N=160)** — statistically consistent with a null
   (zero) shift once the huge realization-to-realization spread (std
   36–70%, growing with W) is properly accounted for at a converged N. **A
   confident non-null claim at the high end of the realistic band does
   NOT survive this scrutiny**, even though the point estimate itself
   stays nominally positive.

**Honest summary for 4.3:** the qualitative claim "disorder breaks the
R4/R5/R7 null" **survives at W≈10Δ** (the lower end of the literature-
realistic band) with high confidence, at a corrected magnitude of ~20-23%
(not the inflated ~29% first reported at N=80) — a genuine, non-vanishing,
two-orders-of-magnitude-above-clean-null effect. At **W≈46Δ (the upper
end) the claim does not survive with confidence**: the N=80 number was
both quantitatively inflated (56.5%→~10%) and statistically overconfident
(3.74σ→~1.2-1.5σ), and a converged, larger-N ensemble cannot rule out a
null shift there. This is reported plainly, not softened: the original
corrected write-up's implicit suggestion that the effect *grows* smoothly
with W up to +57% does **not** hold up once N-independence is checked —
growth with W was itself partly a finite-size artifact concentrated at
N=80.

#### Survivorship-bias check — are the residual failures correlated with boundary displacement?

The corrected-result write-up flagged, but did not test, a hypothesis:
that the residual self-consistency-loop failures might preferentially
exclude the *most*-displaced (most disorder-affected) realizations,
biasing the reported mean toward smaller values. Tested directly at N=80,
W=46 (the point with the largest effect and most failures) by reproducing
the same 30 seeds and recording, for every realization, both the rigid
(non-self-consistent) boundary displacement and whether the actual
self-consistent solve converged:

- **|displacement| among CONVERGED realizations:** mean 1.013, std 0.925 (n=21)
- **|displacement| among FAILED realizations:** mean 0.988, std 0.836 (n=9)

**These are statistically indistinguishable** (difference of 0.025 is tiny
relative to either std of ~0.83-0.93). **The survivorship-bias hypothesis
is not supported**: whether a given disorder realization's self-consistent
solve converges does not appear to depend on how far that realization's
rigid boundary has moved. The residual failures behave more like a
roughly displacement-independent numerical fragility (consistent with the
"local g_crit" mechanism proposed above, which need not correlate simply
with the rigid boundary's own displacement) than a systematic exclusion of
the most-affected realizations — so the reported means are not obviously
biased in this particular way, though this does not rule out other,
untested forms of bias.

---

## Bottom line

**On-site disorder at literature-realistic strength meaningfully matters —
but the honest verdict is a nuanced, W-dependent one, not a uniform
"disorder breaks everything" story**, and this track went through two
rounds of self-correction to get there:

- **Boundary smearing (4.1, well-converged, n=50):** smooth (not sharp, not
  erratic) but real — up to ~30% of realizations near the clean boundary
  are misclassified at the top of the realistic band, saturating there
  rather than growing without bound. Unaffected by either Task 4.3
  correction below.
- **Spurious sub-gap states (4.2, well-converged, n=50):** the ensemble
  *mean* count does not grow with W, but the realization-to-realization
  *spread* grows substantially — a physical mechanism for the sporadic
  (not universal) disorder-induced sub-gap structure the literature already
  associates with the contested zero-bias-peak claims. Unaffected by
  either Task 4.3 correction.
- **Inductive-feedback shift (4.3) — corrected TWICE, by two separate
  post-hoc reviews:**
  1. *Round 1:* a **genuine coding bug** (a fixed, not per-realization,
     search center) was caught, root-caused, and fixed after the first pass
     reported an invalid n=1-4 "ensemble." This produced a properly-
     converged N=80 result showing the shift growing smoothly from +0.2%
     (W=1Δ) to +29%→+45%→+50%→+57% across the realistic band.
  2. *Round 2 (N-independence + survivorship check):* that N=80 result was
     **itself substantially inflated by a finite-size effect** — re-run at
     N=120 and N=160 (the nearest N above N=80 with a verified-clean SC
     baseline; N=60/90/100 hit a genuine or lucky-path dead-zone at this
     g), the shift at W=46Δ collapses from +56.5% (N=80) to +9.6-10.8%
     (N≥120), and its statistical significance drops from 3.74σ to only
     1.2-1.5σ from zero — **a confident non-null claim at high W does not
     survive.** At W=10Δ, by contrast, the effect stays real and highly
     significant (5.9-7.8σ) at every N tested, settling around +20-23%
     (smaller than N=80's +29.4%, but genuinely non-zero, two orders of
     magnitude above the clean-model null). The residual self-consistency
     convergence-failure rate is *also* a finite-size effect — it reaches
     100% success at N=160 for both W points tested — and a direct
     survivorship-bias check (does convergence failure correlate with
     boundary displacement?) found **no** correlation (|displacement| 1.01
     vs 0.99 for converged vs failed, statistically indistinguishable),
     so the reported means are not obviously biased in that particular way.

**Honest, final headline:** disorder measurably smears the topological
boundary (4.1) and adds sub-gap-state variance (4.2) at every realistic W
tested — those findings are unaffected by anything above. The
self-consistent inductive-feedback shift (4.3) is genuinely non-null and
robustly significant **at the lower end of the realistic band (W≈10Δ,
~20-23%)**, but a confident non-null claim **does not survive
N-independence scrutiny at the upper end (W≈46Δ)** — there, the original
finding was both quantitatively inflated (~5x) and statistically
overconfident once properly converged at larger N. This is reported
plainly because it is what a second, more careful look actually found —
not because it is the more exciting or the more comfortable answer.

**Every number above carries its realization count, spread, and (for 4.3)
its statistical significance from zero; the W=0 gate passes exactly (all
diffs 0.0) at every N tested; the parallel infrastructure is bit-identical
serial-vs-parallel; the realistic W scale is grounded in a directly
re-verified primary source (Ahn, Pan, Woods, Stanescu, Das Sarma,
PRMaterials 5, 124602 (2021), arXiv:2109.00007).** Stated compute-driven
limitations, not silent ones: Task 4.3's n_real=30 (vs 4.1/4.2's n=50); the
N-independence check's N-range narrowed to {80,120,160} (not the suggested
2-3 additional points beyond a wider range) given measured ~N³ cost growth;
and only one (N,W) combination (160,46) got a full convergence-in-
realization-count schedule, the others relying on N=80's own already-
demonstrated convergence at n=30 for the same disorder-distribution family.
**Two concrete process lessons from this track:** (1) a first version of
the Task 4.3 re-run driver lacked the `if __name__ == "__main__":` guard
required by Windows' `multiprocessing` spawn method — Python's own safety
check caught this cleanly (`BrokenProcessPool`, zero residual processes, no
data corruption) rather than a runaway process cascade; (2) the first
N-independence run was lost mid-computation when the harness session was
torn down (not a crash — zero residual processes) after completing 3 of 4
(N,W) combinations; those 3 results were recovered verbatim from the
surviving log rather than re-paying ~1100-2900s of sunk compute each, and
only the missing combination was re-run.
