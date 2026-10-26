# R9 Track 1 — characterizing the reactive-impedance pole near g_crit

## Headline (honest)

R8 flagged an apparent s→1 pole in Z_r=iωL_geo/(1−s) as g→g_crit. On close characterization, **at every lattice actually simulated the reactance does NOT diverge**: along the physical principal branch the enhancement 1/(1−s) rises to a maximum of **~5.7×** (N=80; s plateaus at ~0.825) and then the static fixed point **detaches** at g≈8.67e+04 via the fluxoid sawtooth structure of the map, **before s ever reaches 1**. So the reactance grows by a bounded factor and then becomes UNDEFINED (no static equilibrium) — R8's 'divergent pole' is, on the accessible static branch, a finite rise cut off by a dead-zone. The one caveat, made explicit below (finite-size section): the bound is **N-dependent** — the cap climbs as the mesh refines — so whether a genuine divergence emerges in the strict continuum limit is not settled by these N; what is robust is that no lattice studied realizes an actual pole.

## Task 1 — pole scaling & fit

* **Enhancement curve (N=80, μ=1.2, a=10 nm, B=ana−0.05):** 1/(1−s) climbs smoothly from ~1.1 (g=1e4) to ~5.7 before detachment.

* **Power-law fit** 1/(1−s) ~ (g_c−g)^(−p). A free-g_c threshold fit is ill-posed (p and g_c trade off), so both are reported: best free fit p=1.72 at g_c=1.26e+05 (R²=0.9998); anchored fits [g_c=8.75e+04→p=0.34 (R²=0.915)]. The best-fit g_c is ABOVE the last attainable g (8.58e+04) and lands near R5's g_crit band — i.e. the rising enhancement, extrapolated, IS consistent with a pole at g_crit, but the physical branch **detaches before ever reaching it** (max attained enhancement 5.7×). So a clean simple pole is not realized on the accessible branch; the free-fit exponent p≈1.7 is an extrapolation across the detachment dead-zone (and, being degenerate with g_c, only loosely determined), reported with that caveat rather than as an observed critical exponent.

* **Solver robustness (resolution, grid DENSITY at fixed ±2-period span):** coarse (n=61): 8.825e+04; medium (n=121): 8.825e+04; fine (n=241): 8.825e+04 (spread=0.0e+00). At this controlled comparison — same scan window, only the point density varied — the detachment g is **robust to density** (all three agree). This does NOT mean the delicacy is absent: a separate axis, the scan WINDOW/span itself (not tested by this density sweep), was observed during development to shift the apparent detachment — a narrow fallback bracket can miss the root the wide scan still finds, because the principal branch threads a sawtooth near-vertical segment. So density-convergence here is a genuine, positive robustness result; span-sensitivity remains the caveat, and every detachment number in this report should be read as accurate to the scan window used, not as absolute.

## Task 1.3 — detachment vs R5's g_crit: same physics, two thresholds

Both are the loss of an attracting static uniform fixed point (R5's uniform-mode bifurcation), but they are **not one sharp number**:
* **Leading edge (this work):** the principal (attracting-side) branch detaches at g≈8.67e+04 — where Z_r's enhancement maxes out. This is the g R8 called 'detachment ~9e4'.
* **R5's `g_crit_predicted_from_slope`** (smallest-|φ| root slope crossing 1, `bdg_r5_selfconsistent_ring.py:1137`) reproduced here at fine g-resolution: **literally unsatisfiable in this band.** The smallest-|φ| root has s=0.84 at g=1.00e5 (still <1), then the branch detaches again over 1.03–1.17e5, and reappears only at g=1.20e5 already at s=1.17 (>1) — there is NO adjacent pair of scanned g with the root straddling s=1, so R5's interpolation rule, applied at this resolution, never fires. This means R5's *published* 'g_crit ~1.1–1.2e5' most likely comes from `g_crit_measured` — the coarser, EMPIRICAL Picard-iteration-failure criterion R5 reports alongside the slope rule (an upper bound by construction) — not from a literal slope-crossing-1 root, which the same dead-zone this work maps also prevents from being pinned down. R5's published band is consistent with the UPPER edge of that same dead-zone (where the reappearing branch is cleanly repelling, s>1).
**Verdict:** they are the same bifurcation viewed two ways; the pole marks its onset (leading edge, g≈8.67e+04), R5's g_crit a point ~30% higher, sitting past the dead-zone this work resolves. Not identical numbers, same physics — and the dead-zone itself is why neither definition can be pinned to a single sharp g.

## Task 1 — finite-size discipline (is the ~6× cap bulk or artifact?)

| scan | N | max enh | detach g |
|---|---|---|---|
| L-axis (a=10nm) | 60 | 4.47 | 8.16e+04 |
| L-axis (a=10nm) | 80 | 5.57 | 8.69e+04 |
| L-axis (a=10nm) | 100 | 6.65 | 8.87e+04 |
| L-axis (a=10nm) | 120 | 7.75 | 8.69e+04 |
| L-axis (a=10nm) | 160 | 8.85 | 8.51e+04 |
| fixed-L=0.8μm | 80 | 5.57 | 8.69e+04 |
| fixed-L=0.8μm | 100 | 6.65 | 1.11e+05 |
| fixed-L=0.8μm | 120 | 7.75 | 1.30e+05 |
| fixed-L=0.8μm | 160 | 8.85 | 1.70e+05 |
| fixed-L=0.8μm | 200 | 9.84 | 2.08e+05 |

Two DISTINCT finite-size statements must be kept apart:
* **At FIXED production g**, Track A already showed the enhancement is a finite-size effect: s∝a→0 under fixed-L mesh refinement, so 1/(1−s)→1. That is unchanged and not what this table measures.
* **The MAX enhancement REACHED before detachment** (at the moving g_pole, which rises ∝~N as a=L/N shrinks) **GROWS** under fixed-L mesh refinement here (5.6→6.6→7.7→8.9→9.8 over N=80→200). So the pole SHARPENS as the mesh refines: the branch reaches a higher s (closer to 1) before detaching. This means the ~6× cap is **not a universal bound** — it is N-dependent. What it does NOT establish is a genuine continuum divergence: at every finite N the branch still detaches at a FINITE enhancement, before s=1, and the trend is not extrapolated to the continuum here (it would need much larger N and carries R7/Track A's fixed-t caveat). **Robust claim:** at all N simulated (60–200) the reactance is bounded and detaches before s=1; the continuum fate of the cap is left open, honestly, rather than asserted.

## Task 2 — is the pole in the dissipative (Kubo) channel? **No.**

Reusing Track B's validated Kubo machinery at g from production up to g≈8.50e+04: (1) **the ω→0 DC gate (Track B == Track A) HOLDS throughout** — worst rel. mismatch 9.6e-04 ≪ 1e-3, so it does NOT break down as s approaches its plateau. (2) **Im χ develops no matching feature**: the absorption weight at the topological edge stays ≤7.2e-04 and the dominant peak stays at the bulk-pairing scale — the enhancement is purely reactive, exactly as Track A's magnetostatic picture requires. The pole (such as it is) lives entirely in Re χ / L_eff, not in the dissipative channel.

## Task 3 — does Track C's SO self-consistency change near the pole? **The SO physics doesn't; its effect on L_eff inherits the pole's sensitivity.**

Rerunning Track C (total kinetic+SO current sources the self-field, re-solve φ*) from g=2.9e4 up toward detachment cleanly separates two things:
* **The SO correction to the loop slope s is essentially CONSTANT** at **-7.1% to -8.2%** across the whole g-range — the same ~−8% Track C found at the production coupling. The SO *physics* (its effect on the feedback gain) does **not** grow or change character near g_crit; it stays the ~8% source-current change, finite-size-limited exactly as before.
* **But its effect on L_eff GROWS**, from -3.2% (production) to -23.5% (near detachment). This is **not** the SO physics amplifying — it is the *pole's own sensitivity*: L_eff=L_geo/(1−s) is a steep function of s as s→its plateau, so a constant ~8% shift in s maps to a growing shift in L_eff. This is precisely the effect Task 3.2 anticipated ('the self-consistency loop itself becomes more sensitive there, independent of the SO physics'): the answer is yes, but the amplification is geometric (a property of the bounded pole), not a new SO effect. The underlying SO correction remains small and its finite-size character is unchanged.

## Task 4 — parallelization strategy (measured on THIS sweep)

* **BLAS backend:** scipy-openblas (DYNAMIC_ARCH, MAX_THREADS=24), multithreaded-capable, no thread-cap env var set by default ({'OMP_NUM_THREADS': None, 'MKL_NUM_THREADS': None, 'OPENBLAS_NUM_THREADS': None, 'NUMEXPR_NUM_THREADS': None, 'VECLIB_MAXIMUM_THREADS': None}). 10 physical cores; 1.95 GB free at bench ⇒ RAM-safe worker count n_proc=3 (each worker re-imports the R4–R8 stack ~0.4 GB, so memory, not cores, is the binding constraint here).

* **Benchmark (12 independent g-points, `s_worker`, n_proc=3):** serial/default-BLAS 84.8s; serial/1-thread-BLAS 85.9s (1.0×); parallel/default-BLAS 36.9s (2.3×); parallel/1-thread-BLAS 13.7s (6.2×). **Fastest: parallel_1thread_blas.** Two separable effects, both measured: (i) for these small 4N=320 matrices the machine-default multithreaded BLAS is *slower* than 1-thread even serially (thread-sync overhead ≫ the tiny eigendecomposition), and (ii) process-parallelism across g adds a further near-n_proc× on top — so capping BLAS to 1 thread/worker is the dominant win, and parallel/1-thread is fastest. 'More BLAS threads' is actively harmful here; measured, not assumed.

* **Equivalence gate:** serial vs parallel (both 1-thread BLAS) max_abs_diff=0.0e+00 over 72 values — bit-identical, the same standard R8 held bdg_parallel to.

* **Config used for the R9 sweeps:** `bdg_parallel` with `default_workers()` (≈min(cores,6,free_GB/0.5)) and 1-thread BLAS per worker — the measured-fastest, oversubscription-free choice, and the one whose results are bit-identical to serial.

## Task 5 — is this a measurable, realistic device quantity? (honest)

Using R4/R5's plausibility decomposition g = Q·(L_k/L_m) and the R6-verified InAs/Al scales (Δ=0.19 meV ⇒ interband ~46 GHz), the leading-edge pole g≈8.67e+04 is ~3.0× the production coupling — reachable in principle by pushing the resonator Q or kinetic-inductance ratio:

| L_k/L_m | Q for g_production | Q for g_pole |
|---|---|---|
| 10 | 2900 | 8668 |
| 30 | 967 | 2889 |
| 100 | 290 | 867 |

**Resolvability verdict — marginal, leaning negative.** Three strikes against using it as a device diagnostic: (1) the feature is **bounded** (~6× reactance enhancement, not a divergence), so there is no sharp resonance to lock onto; (2) it sits on a **detachment/dead-zone** whose g-location is resolution/lattice-dependent (the numerical detachment g shifts with scan resolution), so in a real device it is a soft, hard-to-place crossover, not a clean pole; (3) reaching it requires operating the feedback loop right at its dynamical-instability threshold, where by construction no stable static equilibrium exists — the wire would be driven into the self-oscillating regime rather than parked at a measurable operating point. So while g_pole is *nominally* within the achievable Q·(L_k/L_m) envelope, the pole is unlikely to be a practically useful, sharp diagnostic; its honest value is as a *qualitative* proximity-to-instability indicator, and even that is a bounded, soft feature.

## Net verdict for the paper's 'so what'

R9 Track 1 **weakens** R8's one exception. On the accessible static branch the reactive impedance's apparent pole near g_crit is **not a realized divergence** but a bounded ~6× enhancement (at N=80) cut off by the fluxoid-sawtooth detachment (the reactive face of R5's loss of a static equilibrium). It is purely reactive (no dissipative Kubo counterpart, DC gate intact), and the SO correction to the loop gain stays small there (its L_eff effect grows only because L_eff is sensitive to s near the pole, not because the SO physics changes). The cap does climb as the mesh refines (the pole sharpens with N), so a genuine continuum divergence is not excluded — but it is also not demonstrated, and at every N simulated the reactance is finite and detaches before s=1. So for the paper's 'so what': Z_r's single orthogonal feature is a soft, bounded (at any finite lattice), resolution-sensitive proximity marker to a dynamical instability — informative physics about *where the static self-consistent description breaks down*, but not a sharp new device observable, and **even more redundant with the gap** than R8's cautious verdict allowed.
