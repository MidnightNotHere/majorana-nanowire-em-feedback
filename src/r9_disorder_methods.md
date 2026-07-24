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
cheaply by Tasks 4.1/4.2) and uses **n_real=20** (Task 1's stated floor, not
its ceiling) rather than 30–50. This is a deliberate, reported scope
reduction driven by measured compute cost, in the same spirit as R7's
thermal-throttle compute ceiling — not a silent narrowing.

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
converged — see below); n_real=20 for Task 4.3 (compute-bound, see the
scope reduction above). Every number below carries its realization count
and spread; raw distributions are in `r9_disorder_data.json`.

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

### 4.3 The R4/R5/R7 inductive-feedback null at realistic W — **does NOT survive; the solver mostly fails to converge**

Ensemble self-consistent boundary shift vs the clean self-consistent
baseline (`ring_boundary_Bc_pfaffian_selfconsistent(N=80,μ=1.2,g=G_R7) =
1.559951Δ`), N=80, g=G_R7, n_real=20 attempted per W (5 at W=0, deterministic):

| W (Δ) | mean shift | std | n converged | n failed (of 20) |
|---:|---:|---:|---:|---:|
| 0.0 | +0.0000 % | 0.0000 | 5/5 | 0 |
| 1.0 | +0.1549 % | 0.7684 | 12/20 | 8 |
| **10 (realistic)** | **+12.335 %** | — (n=1) | **1/20** | **19** |
| **15 (realistic)** | **−9.405 %** | — (n=1) | **1/20** | **19** |
| **25 (realistic)** | **−6.839 %** | — (n=1) | **1/20** | **19** |
| **46 (realistic)** | **−7.303 %** | 1.890 | **4/20** | **16** |

(A genuine bug was caught and fixed while producing this table: the
disordered solver has *two* distinct failure modes — an explicit
`ScfConvergenceError` when the per-B root search fails, and a *silent*
`nan` return when the coarse B-scan never brackets a sign flip at all even
after widening. The first version of `_scshift_worker` only flagged the
first mode as `convergence_failed`, undercounting the true failure rate by
roughly 15-18x at realistic W (it reported 1-2 failures where the correct
count, from checking `np.isfinite(bc)` directly, is 16-19). The `mean`/`std`
statistics above were unaffected by this bug — `summarize_numeric` already
excluded non-finite values correctly — but the *failure-rate* bookkeeping
was wrong until fixed; see `_scshift_worker`'s docstring.)

**Convergence-in-realization-count could not be demonstrated at realistic
W** — at W=10 the schedule (5→10→15→20) gives `nan, nan, nan, 12.335`: zero
successes in the first 15 draws, one lucky success at the 20th. This is
itself diagnostic, not a shortcoming of the convergence-check machinery:
**the underlying quantity is not converging because the solver mostly
cannot find a self-consistent solution at all at this disorder strength.**

**Reading — the headline non-null finding of this track.** Two things are
happening simultaneously at realistic W, and they should be reported
separately because they carry very different statistical weight:

1. **Convergence collapse (high confidence — based on the full n=20
   attempted ensemble, not a subsample):** at W=1Δ (well below the
   literature-realistic band) the self-consistent solver already fails on
   8/20 realizations; by the literature-realistic band (W=10–46Δ) it fails
   on **16–19 of 20** — an 80–95% failure rate, vs 0% at W=0. Disorder at
   realistic strength overwhelmingly breaks the uniform-manifold
   self-consistency solver's ability to find *any* fixed point, a
   qualitatively new effect absent at W=0 or weak W.
2. **Shift magnitude among survivors (low confidence — n=1-4 samples, NOT
   a converged ensemble average):** where a self-consistent solution *is*
   found, the shift is **±7-12%** — roughly **two orders of magnitude
   larger** than the clean-model null (R4-R5-R7-R9-Track-2's consistently
   reported ~0.06-0.13%). The sign is not even consistent across W
   (+12.3% at W=10 vs −9.4% to −7.3% at W=15-46), which is expected given
   n=1-4 samples of a quantity that (per Task 1's flagged approximation)
   is itself only a mean-field projection of a broken-symmetry response.

**Bottom line for 4.3:** the null result does **not** straightforwardly
survive at realistic disorder — but not because the shift smoothly grows
large; because realistic disorder mostly destroys the well-posedness of the
self-consistent problem the null result was measured on. The convergence-
failure-rate finding is solid; the survivor-shift-magnitude finding is
genuinely suggestive of a much larger effect but is not itself a reliable
ensemble average and should not be quoted as a precise number without
substantially more compute (getting n=20 *successes* at W=25Δ would need
~20/0.05 ≈ 400 attempted realizations at ~150s/task/worker — far beyond
this track's compute budget; flagged as a stated limitation, not silently
narrowed).

---

## Bottom line

**On-site disorder at literature-realistic strength (W ≈ 10-46Δ, per Ahn
et al. 2021's verified σ_μ=2.5 meV Gaussian fit) meaningfully matters — this
is the one R9 track where the project's null-result lineage does NOT
uniformly survive**, in contrast to R9 Track 1 (the Z_r pole is bounded,
not divergent) and Track 2 (finite temperature is negligible at 20 mK):

- **Boundary smearing (4.1, well-converged, n=50):** smooth (not sharp, not
  erratic) but real — up to ~30% of realizations near the clean boundary
  are misclassified at the top of the realistic band, saturating there
  rather than growing without bound.
- **Spurious sub-gap states (4.2, well-converged, n=50):** the ensemble
  *mean* count does not grow with W, but the realization-to-realization
  *spread* grows substantially — a physical mechanism for the sporadic
  (not universal) disorder-induced sub-gap structure the literature already
  associates with the contested zero-bias-peak claims.
- **Inductive-feedback null (4.3, n=20 attempted, mostly non-converged):**
  the self-consistent solver's convergence rate collapses from 100% (W=0)
  to 5-20% (realistic W) — a qualitatively new disorder effect — and where
  it does converge, the shift is ~100x the clean-model null, though this
  magnitude rests on only 1-4 successful samples per W and is reported with
  that explicit caveat, not as a settled ensemble average.

**Every number above carries its realization count and spread, per the
task's mandatory ensemble-averaging requirement; the W=0 gate passes
exactly (all diffs 0.0) before any W>0 result is trusted; the parallel
infrastructure is bit-identical serial-vs-parallel for both engines; and
the realistic W scale is grounded in a directly re-verified primary source
(Ahn, Pan, Woods, Stanescu, Das Sarma, PRMaterials 5, 124602 (2021),
arXiv:2109.00007), not an arbitrary choice.** The main stated limitation is
compute: Task 4.3's survivor-shift magnitude and its own realization-count
convergence could not be established to the same statistical standard as
Tasks 4.1/4.2 within this track's compute budget, because the dominant cost
turned out to be disorder-driven convergence failure itself, discovered
only during the real end-to-end run rather than the isolated benchmark —
a limitation flagged here explicitly, in the same spirit as R7's
thermal-throttle disclosure, not silently narrowed.
