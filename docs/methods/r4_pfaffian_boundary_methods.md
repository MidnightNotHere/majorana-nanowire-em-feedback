# Methods note: N-independent extraction of the topological critical field B_c via a Pfaffian Z2 invariant

**Project:** Inductive Circuit Feedback Shifts the Topological Phase Boundary in InAs/Al Majorana Nanowires
**Module:** `src/bdg_r4_em_coupling.py` (R4)
**Date:** 2026-07-09

## 1. The problem with the old detector

The original absolute-B_c detector, `openwire_boundary_Bc`, flags the topological
transition as the first Zeeman field `B` at which the lowest BdG excitation drops
below a fixed absolute threshold, `min|E| < 0.25·Δ`. This is fine for *relative*
shifts at matched wire length `N`, but it is not a defensible way to report an
*absolute* B_c, for two linked reasons established numerically in this project:

1. **Threshold sensitivity.** Sweeping the threshold at fixed `N=80, μ=0` moved
   the reported B_c across a ~50-percentage-point band:

   | threshold | reported offset from analytic B_c |
   |-----------|-----------------------------------|
   | 0.05·Δ    | +17.1 % |
   | 0.10·Δ    | +5.8 %  |
   | 0.25·Δ    | −15.8 % |
   | 0.40·Δ    | −33.5 % |

   A detector reporting a genuine, sharp gap closing should be nearly
   threshold-independent (the true gap closes to zero). This much sensitivity is
   the signature of a detector picking up finite-size structure, not a transition.

2. **N-dependence in the wrong direction.** The offset *grows* with `N`
   (−4 % at N=40 → −22 % at N=160). True Majorana end-mode hybridization
   splitting *shrinks* exponentially with N, so this is not that. It is ordinary
   (non-Majorana) finite-size sub-gap level structure: the near-gap level spacing
   shrinks like ~1/N, so as N grows more ordinary levels fall below the fixed
   absolute threshold before the real bulk transition, biasing B_c downward more
   for longer wires.

So `0.25·Δ` is only a "relative-shift, matched-N" convention, not an absolute
B_c. `openwire_boundary_Bc` is kept intact for that validated use; the Pfaffian
detector below is added alongside it for absolute B_c.

## 2. The new detector: real-space Pfaffian Z2 sign flip

A gapped, particle–hole-symmetric 1D BdG Hamiltonian carries a Z2 topological
invariant — the Kitaev "Majorana number" — equal to the sign of the Pfaffian of
the Hamiltonian in its Majorana (real antisymmetric) basis. The invariant is a
**discrete** quantity: it can only change at a genuine bulk gap closing, and
there is **no threshold to choose**. We locate B_c as the field at which this
sign flips.

Concretely (`ring_boundary_Bc_pfaffian` / `_ring_pfaffian_sign`):

1. Build the periodic-ring LO BdG Hamiltonian `H(B)` (4N×4N, basis
   `kron(spin, site)`, hole block `−h*` so PH symmetry holds exactly).
2. Map to the Majorana basis with the unitary
   `Ω = (1/√2)[[I, I], [iI, −iI]]` on the 4N Nambu space:
   `A = i · Ω (H/scale) Ω†`. This `A` is **real antisymmetric** (verified to
   ~1e-15; asserted at runtime). Rescaling `H` by its largest entry avoids
   overflow of the Pfaffian (which grows like `scale^(2N)`); because the rescale
   is a positive scalar it multiplies `Pf` by `scale^(2N) > 0` and therefore
   **leaves the sign unchanged**.
3. `ν(B) = sgn Pf(A)` via `pfapack`. Find the B where `ν` flips (coarse bracket,
   then bisection to 1e-4).

For a spatially **non-uniform** self-field, `_build_ring_lo_phased` accepts a
per-bond phase array `phi_bonds` (length N, including the wrap bond); PH symmetry
is preserved for any such phase, so the invariant remains well defined. This is
what lets the detector measure the *genuine, non-gauge* inductive boundary shift
(a uniform phase does not move the boundary; a non-uniform texture does).

## 3. Why it is N-independent (and an honest deviation from the literature recipe)

The canonical real-space formulation (see §5, arXiv:2603.08786, Eq. 29) uses the
**product** of Pfaffians under periodic and anti-periodic boundary conditions,
`ν = sgn[Pf(H_PBC)·Pf(H_APBC)]`. **We deliberately use the periodic factor
alone.** Justification, verified numerically for this model:

- In the LO regime the anti-periodic (k=π) Pfaffian keeps a **fixed sign** across
  the transition — the k=π gap (~|2t+μ|, with t≈66) never closes anywhere near
  B_c — so it is a constant prefactor that does **not** move the flip location.
- Worse, at finite N the APBC factor has its own finite-size-shifted
  zero-crossings, which *add spurious flips* near B_c (an observed double-flip
  ambiguity window ~0.1·Δ wide at N=60). Dropping it removes that ambiguity.
- **k=0 is an allowed PBC momentum for every N.** The bulk band-inversion field
  B_c = √(μ²+Δ²) is set by the k=0 gap closing, which therefore carries *zero*
  finite-size shift at any N. This is the mechanism behind the observed
  N-independence.

This is a regime-specific simplification of the literature invariant, not the
verbatim recipe — stated explicitly here and in the code so it is not mistaken
for the full PBC×APBC construction.

## 4. Validation results

`validate_pfaffian_boundary()` (Figure 6, `r4_figure6_pfaffian_vs_threshold.png`):

- **Known-answer, μ=0 (analytic B_c = Δ = 1.0):** the Pfaffian detector returns
  B_c = 1.00004 (bisection resolution) for **every** N in {20, 40, 60, 80, 120,
  160} — spread across N = 0.0 exactly — while the old threshold detector drifts
  from 1.44 (N=20) to 0.78 (N=160), a spread of 0.66·Δ.

### Stress tests (added in response to three specific concerns)

1. **Is "spread = 0" a coarse-grid artifact?** No. The detector uses bisection,
   not a fixed grid, so its accuracy is set by the tolerance `tol`. Refining
   `tol` at μ=0 drives B_c → 1.0 in lock-step, identically for N=20 and N=80:

   | tol | B_c error (both N=20 and N=80) |
   |-----|--------------------------------|
   | 1e-2 | +2.50e-3 |
   | 1e-3 | +3.12e-4 |
   | 1e-4 | +3.91e-5 |
   | 1e-6 | +3.05e-7 |
   | 1e-8 | +4.77e-9 |

   The error tracks `tol` (it is a bisection-convergence offset, not a grid
   snap), and the two wire sizes agree to every printed digit at every `tol`.

2. **Does it hold at the headline μ (not just the symmetric μ=0 point)?** Yes.
   At μ=1.2 and μ=1.4 — the region of the largest inductive shift, near the
   boundary minimum — B_c is identical across N={20,40,60,80,120} to bisection
   tolerance (±3e-7), and matches √(μ²+Δ²): 1.562049630 vs 1.562049935 (μ=1.2),
   1.720465359 vs 1.720465053 (μ=1.4). No degradation near the boundary or at
   the edge of the sweep region that previously bit the threshold detector.

3. **Is it just getting lucky on round-number B_c by symmetry?** No. At μ=0.7
   and μ=1.3 the analytic B_c is irrational (√1.49 = 1.220655562,
   √2.69 = 1.640121947). The detector returns 1.220655523 (err −3.8e-8) and
   1.640121985 (err +3.8e-8) respectively, **identical across N={20,60,120}**.
   It reproduces a non-round crossing to ~1e-8, so the μ=0 accuracy was not a
   symmetry coincidence.

The N=20 point being already exact is the expected, defensible physics: a
topological invariant is a discrete global (fermion-parity) quantity set by the
k=0 band inversion, which is far more robust to finite-size effects than the
continuous spectral gap the old detector measured — and k=0 is an allowed
momentum even for a 20-site ring.

## 4b. Consequence for the headline: the reported boundary SHIFT is ~10× inflated

The point of a threshold-free detector is not just cleaner absolute B_c — it
changes the project's headline number. Two controlled comparisons (N=80):

**Absolute B_c, isolating detector method (both on the clean ring):** the
threshold detector (0.25·Δ) lands **exactly −0.25000 below** the Pfaffian/analytic
B_c at μ = 1.2, 1.3, 1.4 — i.e. off by exactly `threshold`, because the bulk gap
closes linearly (E ≈ |B − B_c|) so a fixed threshold fires exactly `threshold`
early. On the open wire the offset grows to −0.29…−0.33 (adds the end/finite-size
contribution). So absolute B_c from the reported detector is ~15–19 % low.

**The relative SHIFT (the actual headline), keeping figure4's *exact* self-field
and only swapping the detector:**

| μ | threshold-detector shift (Figure 4/5) | ring/Pfaffian shift (same self-field) |
|---|---|---|
| 1.2 | −5.86 % | **−0.77 %** |
| 1.4 | −9.36 % | **−0.60 %** |

The headline shift is **largely a threshold-detector artifact**: the tiny
self-field (~0.015 rad peak) barely moves the true k=0 band inversion, but it
noticeably moves the near-zero *sub-gap* states that the absolute threshold keys
on. The **sign/direction survives** (topological region grows under inductive
feedback); the **magnitude is ~an order of magnitude smaller** than reported, and
even the μ-ranking changes (μ=1.4 is the threshold method's "largest" shift but
not the Pfaffian method's). Figure 4/5 are therefore labelled PRELIMINARY
(threshold detector) until R5.

Caveat on the 0.6–0.8 %: that number itself is not yet final — the self-field
texture is open-wire-derived but applied to the ring with the wrap bond carrying
only the gauge phase (a modelling choice). The *order-of-magnitude* conclusion
(rigorous shift ≪ threshold shift) is robust to that choice; the exact sub-percent
value is an R5 deliverable.

## 4c. R5 deliverable (explicit)

R5 must produce a **validated** final magnitude, not a second unvalidated number:
make the ring self-field self-consistent (remove the open-wire seed and the
wrap-bond placeholder), extract the shift with `ring_boundary_Bc_pfaffian`, and
subject that shift to the same checks this detector passed — (a) N-independence
(N ≈ 20–160), (b) tolerance/grid robustness, (c) a gauge/known-answer check
(uniform flux → 0 shift), and (d) seed-independence of the fixed point.

## 4d. R5 RESULT — the self-consistent ring self-field (what was actually found)

**Module:** `src/bdg_r5_selfconsistent_ring.py` (new; `bdg_r4_em_coupling.py` is
imported read-only and unmodified — all 6 R4 validations still pass).
**Outputs:** `results/self_consistency_bifurcation/` (Figures 7–9, `r5_results.json`,
`r5_shift_comparison.csv`). **Headline:** the rigorous, fully self-consistent
shift is **+0.05 %** at N=80 (μ=1.2 and 1.4) and **extrapolates to ≈ 0 in the
bulk** — a further ~15× down from the interim ~0.6–0.8 %, and it changes sign.

### What was built (removing both interim placeholders)

1. **Ring current density** including the wrap bond (`ring_current_density`,
   `np.roll` periodic indexing) — same SI prefactor / electron-block-only /
   E<0-occupied conventions as the validated open-wire current.
2. **Ring Biot–Savart self-field** (`ring_biot_savart_self`): the validated
   open-wire kernel *adapted* — not rewritten — for the closed geometry. The
   wire closes on itself, so the straight-segment kernel genuinely breaks
   (distances wrap; source `dl_j` and field-point tangent `t_i` are no longer
   parallel). Only the geometry factors change: chord distance
   `2R sin(Δθ/2)` and the tangent projection `cos(Δθ)` for a circular loop of
   circumference `Na`; the local finite-segment self-term `2 ln(a/r_wire)` and
   the SI-units path are unchanged (segment locally straight, `a ≪ R`). This is
   the "genuinely different treatment of the self-field integral" flagged in the
   task, and it is validated (`validate_ring_biot_savart`, the closed-geometry
   analogue of R4 Validation 2): for a uniform current the discrete loop flux
   reproduces the discrete-exact kernel to **0.003 %**, and the textbook
   thin-wire loop inductance `L = μ₀R[ln(8R/r_wire) − 2]` to **2.25 %** (an
   understood constant, the difference between the loop's `ln(8R/r)` and the
   straight-segment `2ln(...)` log factors). `g = 3×10⁵` is unchanged — the
   amplification argument is about the *circuit*, not wire geometry, and the
   geometric-inductance difference between a segment and a loop is only
   logarithmic, well inside `g`'s order-of-magnitude bars.
3. **Self-consistency** (`selfconsistent_ring_texture`) and a **fully
   self-consistent Pfaffian B_c** (`ring_boundary_Bc_pfaffian_selfconsistent`):
   at every field the bisection probes, the texture is converged *at that B*
   first, then the Z₂ sign is read off the dressed ring.
4. The **interim recipe is kept, labelled SUPERSEDED**
   (`interim_ring_Bc_openwire_texture`) so the ~0.6–0.8 % number is reproducible
   in committed code (it previously existed only as session work).

### Convergence was NOT clean — and that itself is a reportable finding

- **Plain fixed-point iteration does not converge at g₀ = 3×10⁵, from any
  seed.** On the ring the self-consistency map closes on *uniform* textures, so
  it reduces to a scalar map `G(φ)`; at the fixed point its slope is
  `G′ ≈ 3.0 > 1`. The iteration is therefore Picard-**repelling**: from
  zero / open-wire / random seeds, at β=1.0 and β=0.3 mixing, the residual
  never contracts — it wanders boundedly (it does **not** blow up). No amount of
  linear mixing fixes a slope-3 fixed point (Figure 7a,b). We did **not** paper
  over this by capping iterations; instead we recognised the ring map is scalar
  and solved `G(φ)=φ` by **bisection**, which needs no stability. Below a
  bifurcation `g_crit ≈ 1.1×10⁵` (slope criterion; first non-convergence
  observed at ~1.2×10⁵) iteration *does* converge, and there it is
  **seed-independent** (all three seeds → the same fixed point, Figure 7c) and
  agrees with the root-finder.

- **The fixed point exists and is (locally) unique, but the map is
  multi-valued.** The self-consistency has *multiple* uniform roots (a fluxoid /
  sawtooth-CPR multiplicity: 1 root at N=40,60; 2 at N=80,120); we take the
  **principal** (smallest |φ|) branch. Path-dependence was checked directly
  (Figure 8): slowly ramping g from 0 to 3×10⁵ (each converged solution seeding
  the next) lands on the **same** principal branch as jumping directly from a
  zero seed, for every g below g_crit — path-independent. Whether a
  *symmetry-broken non-uniform* solution also exists was tested decisively with
  **Newton–Krylov** on the full N-dimensional `G(φ)−φ=0`
  (`nonuniform_fixed_point_search`), seeded with large perturbations (up to 3
  flux-periods) along the most-unstable modes: from every convergent seed it
  collapses back to the uniform solution (texture norm ~5×10⁻¹¹). **No
  non-uniform fixed point exists.** (Honesty note: an early N=40 spot-check saw
  all texture-mode gains < 1 and this was briefly mis-stated as "textures always
  decay". That is wrong at N=80, where two texture modes have Picard gain
  1.746 > 1 — the uniform point is a *saddle under iteration*. But saddle-under-
  iteration ≠ a second solution: the Newton–Krylov uniqueness result is the
  correct statement, and the code/figure were corrected to say so.)

### The validation battery the interim number never had

- **N-independence (Figure 9a).** Shift = 0.122 / 0.085 / 0.051 / 0.024 % at
  N = 40 / 60 / 80 / 120 — monotonically → 0, as `~N^−1.5`. This is the key
  physical result: the total self-flux `N·φ*` is a near-constant fluxoid
  (≈ 0.5 rad), so the per-bond self-phase ~1/N and the residual boundary move is
  a **finite-size Aharonov–Bohm flux effect, not a bulk shift**. The
  thermodynamic-limit self-consistent shift is **consistent with zero**.
- **Tolerance robustness (Figure 9b).** Tightening the bisection tol
  (1e-4→1e-6) and the root tol (1e-8→1e-12) holds the shift in
  0.0505–0.0530 %; a warm-vs-cold-start difference of ~0.002 pp is the only
  residual, traceable to the root multiplicity above.
- **Gauge / decomposition check.** The self-consistent B_c equals the B_c of a
  *frozen* uniform flux `φ_ext + φ*(B_c)` to 8.4×10⁻⁵ — i.e. the effect is
  almost purely static-flux (AB) content, with a tiny genuine-feedback residual
  from the B-dependence of φ*. (A strictly uniform external flux with no
  self-field gives 0 shift, as R4 Validation 6 already showed.)

### The progression, side by side (the point of the exercise)

| μ | (a) threshold detector (Fig 4/5) | (b) interim ring/Pfaffian (SUPERSEDED) | (c) **self-consistent ring (this work, N=80)** |
|---|---|---|---|
| 1.2 | −5.86 % | −0.768 % | **+0.051 %** |
| 1.4 | −9.36 % | −0.607 % | **+0.047 %** |

Each column is a further ~10× reduction in magnitude. The fully self-consistent
number at N=80 is ~120× smaller than the original headline, and it **changes
sign** (a small *positive* AB shift, not the negative shift the earlier methods
reported) — and the N-scan shows even that +0.05 % is finite-size and heads to
zero. **The robust conclusion is now the null-ish one:** to the accuracy of this
model the inductive self-field does *not* move the bulk topological boundary at a
level distinguishable from a finite-size flux artifact. The earlier "topological
region grows under feedback" direction does **not** survive making the self-field
self-consistent on the ring.

### What remains provisional even after this work

- **`g = 3×10⁵` is still a plausibility placeholder, not a derived coupling.**
  R5 (this) fixed the *geometry and self-consistency* of the self-field, not the
  circuit amplification. A genuine circuit model (L_k, L_m, Q, operating point)
  is still owed, and at the assumed g₀ the equilibrium sits *past* the Picard
  bifurcation — the static fixed point is well defined (root-found), but the
  *dynamical* response of a real feedback circuit near/above g_crit would need
  its own analysis before any claim about time-domain behaviour.
- **The +0.05 % at N=80 is finite-size.** The honest bulk number is
  "≈ 0, bounded by the N-scan"; +0.05 % should be quoted only with its N and its
  →0 extrapolation, never as a bulk magnitude.
- **Ring vs open-end.** As in §5, this is the bulk/ring invariant, not literal
  open-end Majorana appearance; the two coincide only in the thermodynamic limit.

## 4e. THE PHYSICALLY MEANINGFUL HEADLINE — shift restricted to g < g_crit (supersedes the g = 3×10⁵ number)

**Module:** `src/bdg_r5_physically_meaningful_shift.py` (imports R4/R5 read-only;
figures `r5_figure10_shift_vs_g.png`, `r5_figure11_headline_validation.png`;
data `r5_physical_shift_results.json`).

**The problem with the §4d number.** §4d reported the self-consistent shift at
g = 3×10⁵. But Figure 8 established a **Picard bifurcation** at g_crit: below it,
slow g-ramp and direct-jump-from-zero converge to the *same* attracting fixed
point (path-independent); above it there is **no attracting fixed point** — the
root-finder still locates a mathematical root of G(φ)=φ, but nothing an
iterative/dynamical process settles to. g = 3×10⁵ sits far above g_crit, so the
"+0.05 %" is the location of a **dynamically unreachable root**, not a physical
static equilibrium. (It is kept, labelled, in the code and in §4d — it is what
revealed g_crit.) The honest headline must come from g < g_crit.

**Two corrections found while doing this (documented for auditability).**
1. *Baseline.* The shift must use the RING/Pfaffian baseline (same detector,
   same geometry, external flux only) — `ring_boundary_Bc_pfaffian` at uniform
   φ_ext. A first-draft sweep wrongly used the OPEN-WIRE THRESHOLD detector as
   baseline, injecting a spurious ~+23 % detector/geometry offset that (tellingly)
   did **not** vanish as g→0. Fixed.
2. *Detector window.* The self-consistent B_c must be called with its validated
   defaults (window = 0.05, n_coarse = 11), as the committed Figure 9 run did.
   A first draft used window = 0.2, whose wide coarse B-grid caused spurious
   bracketing failures far below g_crit — an artifact of the override, not physics.

**Three regimes of g** (established by direct diagnostics — root structure and
Picard slope scanned across the transition window):

| regime | range | status of the static self-consistent boundary |
|---|---|---|
| 1 | g ≲ 3×10⁴ (`G_DEADZONE_ONSET`) | attracting fixed point at **every** B across the transition → B_c **unambiguously** defined |
| 2 | 3×10⁴ ≲ g ≲ 1×10⁵ | fixed point still attracting at most B, but a **narrow dead-zone in B opens exactly at the gap-closing** (the sawtooth map jumps across the identity line without crossing it) → B_c defined only up to the dead-zone width; the Pfaffian sign flip itself stays clean (single flip) |
| 3 | g ≳ 1×10⁵ (`G_CRIT_SINGLEPOINT`, Figure 8) | fixed point **Picard-repelling → no attracting equilibrium**. g = 3×10⁵ lives here. |

The dead-zone onset was bisected to g ≈ 2.98×10⁴ at N = 80, μ = 1.2; it is
**N-dependent** (lower for smaller N — e.g. N = 40 is already in the dead-zone at
g = 2×10⁴). Regimes 2–3 are themselves a finding: the static mean-field
self-consistency **breaks down at the transition** before the coupling is large
enough to matter for a large shift.

**The dose-response (Figure 10).** With the correct baseline, in regime 1 the
shift is **small, NEGATIVE, and → 0 as g → 0** (the physical requirement the
+23 % artifact violated), growing monotonically with g and **path-independent at
every g checked** (ramp = jump = root). μ = 1.2 and μ = 1.4 lie on top of each
other:

| g | shift μ=1.2 | shift μ=1.4 |
|---|---|---|
| 1×10³ | −0.0013 % | −0.0014 % |
| 1×10⁴ | −0.0145 % | −0.0151 % |
| 2×10⁴ | −0.0379 % | −0.0366 % |
| 2.9×10⁴ (top of regime 1) | **−0.064 %** | **−0.066 %** |

**Headline g (a judgment call, stated as such).** The shift is monotonic in g, so
per the least-arbitrary rule the natural representative coupling is the **top of
the unambiguous regime**, g ≈ 2.9×10⁴ (the dead-zone onset). Headline number:
**ΔB_c/B_c ≈ −0.064 % (μ=1.2), −0.066 % (μ=1.4), at N = 80.** We deliberately do
**not** pick the g that maximises the number, and we do not use regime-2/3 values.

**Is even this a bulk effect? NO — it is a finite-size artifact (Figure 11a).**
The decisive test. At a clean regime-1 coupling (g = 1×10⁴, where every N is
unambiguous) the shift **vanishes as N → ∞**:

| N | shift μ=1.2 | shift μ=1.4 |
|---|---|---|
| 40 | −0.0557 % | −0.0521 % |
| 60 | −0.0259 % | −0.0249 % |
| 80 | −0.0145 % | −0.0151 % |
| 120 | −0.0076 % | −0.0076 % |

|shift| ∼ **N^−1.8** for both μ (fit exponents −1.82 and −1.75) — the **same
finite-size Aharonov–Bohm scaling** as the g = 3×10⁵ result (§4d, ∼N^−1.5). It is
**tolerance-robust** (Figure 11b: −0.0371 → −0.0380 → −0.0379 → −0.0379 % as the
bisection tol tightens 1e-4 → 1e-7).

**THE HONEST BOTTOM LINE (a null result).** Across the **entire physically
accessible coupling range** (every g where a static self-consistent equilibrium
exists), the self-consistent boundary shift is (a) **tiny** (< 0.1 %), (b) slightly
**negative** (topological region marginally shrinks — the *opposite* sign to the
§4d +0.05 %, which was sign-flipped only because it sits past the sawtooth jump in
regime 3), and (c) a **finite-size Aharonov–Bohm artifact that extrapolates to
zero in the bulk**. **The inductive self-field does not move the bulk topological
boundary anywhere a static description applies.** The direction the earlier stages
suggested ("topological region grows/shrinks under feedback") does not survive as
a *bulk* effect at all.

**The full lineage (each stage superseding the last — the story, not just the last number):**

| # | method | μ=1.2 | μ=1.4 | verdict |
|---|---|---|---|---|
| 1 | open-wire threshold detector (Fig 4/5) | −5.86 % | −9.36 % | PRELIMINARY — mostly a threshold artifact |
| 2 | interim ring/Pfaffian, open-wire texture on ring | −0.768 % | −0.607 % | SUPERSEDED — stopgap self-field |
| 3 | self-consistent ring, g = 3×10⁵ (§4d) | +0.051 % | +0.047 % | **mathematically valid root, but g > g_crit — NO attracting fixed point exists; dynamically unreachable** |
| 4 | **self-consistent ring, g < g_crit (this §4e)** | **−0.064 %** | **−0.066 %** | **physically meaningful upper bound; but finite-size (∼N^−1.8 → 0 in bulk)** |

**Revised status of the g-parameter literature justification (honest).** The
original g = 3×10⁵ plausibility argument (CPW resonator Q ∼ 10³–10⁴ × kinetic-
inductance enhancement ∼10–30×) sketched a range spanning roughly 10⁴ to ≳10⁶.
That range **straddles all three regimes**: its low end (∼10⁴–3×10⁴) is in the
physical regime 1, its middle in regime 2, its high end (≳10⁵, and the nominal
3×10⁵ itself) in regime 3 where no static equilibrium exists. So — stated plainly,
without retroactive narrowing — **the literature argument was never specific
enough to pin g below g_crit; it is compatible with both sides**, and the specific
value it was tuned to (3×10⁵, chosen to justify a ~9 % effect) lies in the
no-static-equilibrium regime. A real device at the *low* end of that plausibility
band would sit in regime 1 and produce the ∼−0.06 % (finite-size) shift above; a
device at the *high* end would not admit the static self-consistent treatment used
throughout this project at all.

**What would have to happen physically at g > g_crit (open question for future
work, not resolved here).** Above g_crit no static fixed point attracts, so the
quasi-static/mean-field self-consistency used here simply does not describe that
regime. Physically this leaves two possibilities, to be settled by future work
(R6+), not by this task: (i) real readout circuits do **not** reach that effective
coupling — the plausibility argument was an over-estimate, and the physical device
lives in regime 1 where the effect is the null/finite-size result above; or
(ii) at such coupling the wire–circuit feedback is genuinely **dynamical** — it
does not relax to a static equilibrium and needs a time-dependent treatment
(possible self-oscillation / limit cycle of current ↔ self-field), for which
"the boundary shift" is not even the right observable. The narrow transition
dead-zone of regime 2 is an early warning that the static mean-field picture is
the first thing to fail as g grows.

## 4f. ANDERSON MIXING — is g_crit real physics, or an artifact of PLAIN (Picard) iteration?

**Module:** `src/bdg_r5_anderson_mixing.py` (imports R3/R4/R5 read-only; the
existing Picard iteration `selfconsistent_ring_texture` and the bisection/
Newton–Krylov root-finders are UNCHANGED and remain the baselines).
**Outputs:** `results/self_consistency_bifurcation/` — Figure 12 (Task 2), Figure 13 (Task 3),
`r5_anderson_results.json`.

**Why this was done.** Every result above `g_crit ≈ 1.1×10⁵` (§4d/§4e) rested on
one solver: plain fixed-point / Picard iteration, which is *provably* repelling
where the scalar-map slope `G′ ≈ 3 > 1` (no linear mixing β can rescue a slope-3
fixed point: `|1−β+βG′| > 1` for every β∈(0,1]). The §4e clause — above g_crit
"there is no attracting fixed point … nothing an iterative/dynamical process
settles to" — had **not** been tested against a stronger iteration. Anderson
mixing (Anderson acceleration) is a multi-secant **quasi-Newton** method for
`F(φ)=G(φ)−φ=0`: it builds a least-squares linear model of `F` from the last `m`
iterate/residual differences and takes Newton-like steps, so — like Newton's
method — it can locate a root even where the bare map is repelling, provided the
root is simple (`J_F = J_G − I` non-singular there). It is precisely the tool to
separate two statements the prior work was careful to keep apart:
**(A)** a fixed point *exists* above g_crit — already established via bisection
and Newton–Krylov, so Anderson's job here is to corroborate it via a
structurally different algorithm, from seeds sufficiently close to the root;
**(B)** the fixed point is *attracting under the physical relaxation dynamics*
— FALSE for |G′|>1, a property of the MAP, not the solver, and the meaning of
"no attracting fixed point".

**Variant implemented (Task 1).** **Type-II Anderson** (the original Walker–Ni
form: γ minimizes the linearized residual ‖f_k − ΔF_k γ‖ directly), chosen over
Type-I because it is the standard, better-conditioned default and degrades
gracefully when ΔF goes rank-deficient near convergence — robustness that matters
because the sawtooth current–phase relation can make the LS step ill-conditioned.
**History depth m = 4.** Regularization: **relative Tikhonov** on the tiny m×m
normal matrix, `λ = reg·trace(ΔFᵀΔF)/m` with `reg = 1e-10` (standard "regularized
Anderson"; the condition number is recorded and never blows up). **Convergence
criterion (Task 1.4):** the relative change in the self-field,
`‖ΔA‖/‖A‖ < 1e-6` (equal to `‖Δφ_self‖/‖φ_self‖` since φ_self ∝ A_t). Because
that step criterion alone can be met at a *sawtooth-branch stall*, the
fixed-point residual `max|G(φ)−φ|` is ALSO recorded, and a result is a genuine
fixed point only if BOTH are small — stalls are flagged, never reported as
solutions. (This distinction is load-bearing: see the far-seed cases below.)

**Warm-start (Task 1.3) — confirmed wired in, and honestly bounded.** R3's
rank-2 PH-symmetric subspace near-zero-mode update (DPT→RR1→RR2→residual/
fallback, with R3's thresholds) is reused via `RingSubspaceWarmStart`, adapted
from R3's open wire + (δμ,δB,δΔ) steps to the RING Hamiltonian + the Peierls-phase
step δφ the SCF actually takes. A direct large-gap probe confirms the raw rank-2
update reproduces the full-eigh near-zero eigenvalues to **9.4×10⁻⁶** — so the
update is genuinely correct — but R3's *conservative* residual safety net
(`resid > 10·‖δH‖`) trips on these phase steps (the eigenvector residual ≈5×10⁻³
just exceeds 10·‖δH‖≈3.6×10⁻³), so in practice it falls back to a full solve
(8/8 probe steps, 4/4 SCF iterations). Either way it cannot help: the ring
current sums over ALL occupied states, so the SCF map needs the **full** occupied
spectrum every iteration — the full `eigh` is **~96–98% of the per-map cost** and
is unavoidable, and a rank-2 near-zero update cannot supply it. **Net: the
near-zero subspace warm-start gives ~1× on the SCF loop, reported plainly and not
conflated with Anderson's own gain.** The separable warm-start that *could* help
is *iterate* seeding (start each solve from a nearby converged texture); it is
isolated in Task 3, where — at these coupling values — it also gives no gain
because cold Anderson already converges at the ~3–5-iteration floor.

### THE TASK-2 RESULT — stated as plainly as the g_crit discovery itself

> **Anderson mixing DOES converge above g_crit, where plain Picard provably
> cannot — but it converges to the SAME already-known, still-Picard-REPELLING
> root, reproducing the IDENTICAL §4d shift (+0.05 %, → 0 in the bulk) about
> 20× faster, from seeds sufficiently close to the root (2 of 4 tried at the
> headline g; see below for the other 2). This does NOT overturn the null
> result. It refines a methodological claim (a stronger solver *reaches* the
> root; g_crit is a Picard-iteration stability threshold, not a barrier to
> *locating* the equilibrium) and makes the physics MORE robust — three
> independent solvers now agree on the same root and the same → 0 shift.**

At the representative near-boundary point (N=80, μ=1.2, B=0.98√(μ²+1)=1.531),
plain Picard wanders without converging at every g below, while Anderson reaches
a **genuine** fixed point (both `‖ΔA‖/‖A‖<1e-6` AND `max|G−φ|` small):

| g | Picard (β=1) | Anderson from zero seed | genuine fixed pts / 4 seeds | seed-independent? |
|---|---|---|---|---|
| 1.5×10⁵ | wanders (res 1.1×10⁻²) | converges, **5 iter**, φ*=−0.01177 | **4/4** (all match bisection) | **YES** |
| 3.0×10⁵ | wanders (res 1.0×10⁻¹) | converges, **5 iter**, φ*=−0.005886 | 2/4 (zero, random) | YES (among genuine) |
| 6.0×10⁵ | wanders | converges, 8 iter, φ*=−0.004711 | 3/4 | NO (fluxoid branches) |
| 1.0×10⁶ | wanders | converges, φ*=+0.0386 (non-principal) | 3/4 | NO (fluxoid multiplicity) |

- **The criterion is genuinely met.** At g=3×10⁵ from the zero seed,
  `‖ΔA‖/‖A‖ = 5.4×10⁻⁷` is reached at **iteration 5**, with fixed-point residual
  **6.3×10⁻⁹** (a true fixed point, not a stall).
- **Same root, same physics, ~20× faster.** The full self-consistent Pfaffian
  B_c at g=3×10⁵ (N=80) via Anderson = **1.561738** (μ=1.2), shift **+0.0505 %** —
  *identical* to the bisection root-finder's B_c (1.5617377, +0.05052 %) at both
  μ, in **15.3 s vs 299.5 s**. Anderson found the SAME root the root-finder did.
- **Seed agreement, stated precisely (not "all seeds agree").** At the headline
  g=3×10⁵: **2 of 4 seeds (zero, random) converge to the known root; 2 of 4
  (far_neg, far_pos) land on a different branch or fail to converge** —
  consistent with the map's multi-valued/fluxoid structure, not a defect of the
  solver. Concretely, the far seeds (±0.02, a full fluxoid period away)
  **stall on another sawtooth branch** — correctly flagged (fp-residual ≈5×10⁻²,
  NOT a fixed point), since the raw `‖ΔA‖/‖A‖<1e-6` step criterion alone would
  have mislabeled the stall as "converged", which is exactly why the
  fixed-point-residual gate is enforced. At very large g (10⁶) far seeds instead
  converge to a genuine *non-principal* fluxoid root (+0.039 / −0.047) — the
  metastable-flux-state multiplicity §4d/§3b already documented. This is a
  strictly more interesting result than blanket seed-independence: Anderson's
  basin of convergence to the principal root is seed-dependent in the same way
  the underlying physics is multi-valued, and the solver reports that honestly
  rather than papering over it.

### Validation battery on the Anderson-converged g=3×10⁵ solution

Because Anderson reaches a genuine fixed point above g_crit, the SAME battery
§4d/§4e used was applied to it — and it reproduces the §4d finding exactly:

- **N-independence — the shift scales to zero as N→∞** (the decisive test).
  Anderson shift vs N: **+0.1223 / +0.0854 / +0.0505 / +0.0238 %** at
  N=40/60/80/120 (μ=1.2, `~N^−1.51`); **+0.1879 / +0.0793 / +0.0469 / +0.0232 %**
  (μ=1.4, `~N^−1.90`). This is the SAME finite-size Aharonov–Bohm scaling as §4d
  (`~N^−1.5`): the bulk shift is consistent with zero.
- **It is the same root.** Anderson B_c = bisection root-finder B_c at N=80 to
  6 digits for BOTH μ (1.561738 vs 1.5617377; 1.720153 vs 1.7201529).
- **Tolerance-robust.** Shift = +0.0530 → +0.0508 → +0.0505 → +0.0505 % as the
  bisection tol and Anderson `‖ΔA‖/‖A‖` tol tighten (1e-4 → 1e-7).
- **Path-independent (among genuine roots).** From 5 seeds at fixed (N,μ,B,g),
  the **3 genuine fixed points agree to 5.5×10⁻¹⁰** (φ*=−0.005886, matching the
  bisection root); the 2 far seeds stall on another branch and are reported
  separately (per the multiplicity above), not silently dropped.

### What this DOES and does NOT change (the honest reading)

Anderson converging above g_crit **refines** — it does not overturn — the §4e
sentence "above g_crit … nothing an iterative/dynamical process settles to." That
clause is too strong as literally written: a *quasi-Newton* iterative process
(Anderson) settles to the root easily. The precise statement is that no
**Picard-class** (fixed-contraction) iteration converges, and that the fixed
point is still **physically repelling** (|G′|≈3>1) — a solver-independent
property of the map. A real analog feedback circuit *relaxing* implements
(approximately) Picard dynamics (current→field→phase→current, applied
repeatedly), not Anderson history-extrapolation, so the physical instability
stands. And the observable is unchanged: Anderson, bisection, and Newton–Krylov
all give the **same +0.05 % shift that scales to zero in the bulk**. **g_crit is a
Picard-iteration stability threshold, not a barrier to *locating* the
equilibrium; the null result (no bulk boundary shift; no attracting *static*
equilibrium above g_crit under relaxation dynamics) STANDS, now on firmer footing.**

### TASK 3 — efficiency below g_crit (honest, not fit to "5×")

Where both converge (fixed near-boundary B, g<g_crit), iterations to
`‖ΔA‖/‖A‖<1e-6` (N=80, μ=1.2):

| g | simple β=1 | simple β=0.5 | Anderson (m=4) | speedup (simple β=1 / Anderson) |
|---|---|---|---|---|
| 1×10³ | 4 | 21 | 3 | **1.3×** |
| 1×10⁴ | 7 | 23 | 4 | **1.8×** |
| 3×10⁴ | 13 | 31 | 4 | **3.3×** |
| 6×10⁴ | 27 | 56 | 5 | **5.4×** |

The measured speedup is **1.3× → 5.4×** and **grows as g→g_crit** (Picard suffers
critical slowing as G′→1; Anderson stays at the 3–5-iteration floor). It is
reported as measured — *not* tuned to the "~5×" figure in the original scope
notes (it happens to reach ~5.4× only at the top of the sub-critical range, and
is well below 5× for most of it). The separable *iterate* warm-start gives no
extra reduction here (cold Anderson is already at the floor). The R3 near-zero
subspace warm-start is wired in and its rank-2 update verified correct (9.4×10⁻⁶),
but it gives ~1× on the SCF loop because the full-spectrum current dominates —
reported separately so the two effects are never conflated.

## 4g. R7 — CIRCUIT INDUCTANCE / WIRE LENGTH L AS A THIRD AXIS (B, μ, L)

**Module:** `src/bdg_r7_inductance_sweep.py` (imports R4/R5 read-only; figures
`r7_figure1_shift_vs_BmuL.png`, `r7_figure2_regime_and_validation.png`; data
`r7_results.json`, `r7_points_cache.json` in `results/inductance_scan/`). R1–R6 code,
figures and results are untouched — this is additive.

**Hypothesis-neutral framing.** R7 re-opens the boundary-shift question along the
one axis R4/R5 never varied: the wire length `L = N·a_lat`, which sets the loop's
self-inductance. It has **no target result**. The old ">1 %" goal is superseded;
R5's validated headline was a null result (§4e); R7 reports what the L-sweep
actually shows, and a continuation of the null is exactly as reportable as an
L-dependence. Everything is at the ONE validated coupling **g = 2.9×10⁴** (top of
regime 1, §4e) — g is held fixed by a guard (`assert_g_consistent` refuses any call
that reverts to the R5 module default 3×10⁵) — at the fixed lattice constant
**a_lat = 10 nm**, μ ∈ {1.2, 1.4}, Δ = 1.

**What "L" is, and the N↔L confound this module is built to break.** In the R5 ring
the circumference is exactly `N·a_lat` (`ring_biot_savart_self`, R = N a/2π), so the
physical length / inductance axis is set by (N, a_lat). At fixed a_lat, growing N
grows L — but it also refines finite-size. R5's N-scan (Fig 9/11) grew both together
(∼N⁻¹·⁸ → 0) and could not separate "shift is small because the ring is physically
long" from "…because the mesh is coarse", because it used a single lattice constant.
R7 breaks the confound with two orthogonal scans:
  1. **L-axis** — fix a=10 nm, vary N ⟹ L = N·a (R5's scan, along L AND N together).
  2. **Fixed-L mesh refinement** — hold L fixed, vary N with a = L/N (finer mesh of
     the SAME physical ring). This is the task's "N-independence at fixed L via
     a_lat": does the shift at a FIXED physical length survive the continuum limit
     (a bulk value) or vanish (a discretisation/finite-size artifact)? R4/R5 never
     did this. (Modelling choice, stated: t_hop is held at its physical 10 nm value
     while only the geometric a in the Biot–Savart/Peierls self-field is refined, so
     the clean-ring baseline B_c is identical at each N and the refinement isolates
     the self-field's mesh dependence; a fully dimensionful continuum limit with
     t∼1/a² is heavier future work.)

The L range here (0.4–2.0 μm intended; **0.4–1.2 μm actually computed** — see the
N-ceiling note below) is a **tractability subrange**, explicitly not a literal
reproduction of Chang et al. 2015's ~5–10 μm wire (N∼500–1000, impractical for
dense self-consistent runs).

**Solver (additive, validated bit-for-bit).** B_c is measured with an a-aware port
of R5's warm-hint continuation bisection, `ring_Bc_selfconsistent_a`, built only
from R5 primitives that thread `lattice_const`. It reproduces
`ring_boundary_Bc_pfaffian_selfconsistent` **exactly at a=10 nm** (|ΔB_c| = 0 at
N=80; `validate_solver`) and follows the same fluxoid branch. The self-consistency
map is multivalued near the transition (fluxoid branches); R5's continuation picks a
specific branch, and the port matches it. Exactly at the gap-closing the principal
self-phase branch detaches over a razor-thin (~3×10⁻⁵) window; whether a single
bisection path lands on it is path-dependent (R5's own `principal_uniform_fixed_point`
returns nan at B=1.560000, and R5's default path spuriously dead-zones some large N,
e.g. N=100). R7 runs R5's algorithm on the R5 default grid first (a success matches
R5 exactly, `retries`=0) and, only on failure, retries with jittered coarse grids; if
**every** path hits a no-root region it is a genuine finite-width regime-2 dead-zone
(no shift reported). Measured no-root window widths: N=80 ≈ 3×10⁻⁵ (knife-edge),
N=40 ≈ 7.5×10⁻⁴ (genuine).

### FINDING 1 — an L-dependent DYNAMICAL-STABILITY boundary (new; R4/R5 could not see it)

At the fixed g=2.9×10⁴, the static self-consistent boundary is cleanly defined only
for **N ≳ 70–80** (L ≳ 0.8 μm at a=10 nm); **N = 40 and 60 (L = 0.4, 0.6 μm) are
regime-2 dead-zones** — no static B_c (the sawtooth self-consistency map jumps across
the identity line at the gap-closing). This is the first direct confirmation of
Task-2.3's suspicion that the dead-zone / g_crit structure is **not L-independent**:
since L is set by N, the regime-1/2 boundary R5 characterised at a single N appears as
an **L-dependent** (and, via the mesh, a-dependent) stability threshold. It is a
STABILITY finding, not a bulk-shift finding — Picard iteration still converges at the
generic transition field for every N tried (`picard_converges_at_ana` = True even at
N=40/60), so the dead-zone is the thin window at the gap-closing, not a repelling
fixed point.

### FINDING 2 — where the boundary IS defined, the shift is tiny and VANISHES under fixed-L mesh refinement

**L-axis (a=10 nm), the clean points:**

| N | L (μm) | ΔB_c/B_c (μ=1.2) | ΔB_c/B_c (μ=1.4) | regime |
|---|---|---|---|---|
| 40 | 0.4 | — | — | **dead-zone** |
| 60 | 0.6 | — | — | **dead-zone** |
| 80 | 0.8 | **−0.0640 %** | −0.0647 % | clean |
| 100 | 1.0 | −0.0430 % | — | clean (robustified) |
| 120 | 1.2 | −0.0342 % | — | clean |

The two μ coincide at N=80 (−0.0640 vs −0.0647 %), as in R5. The fixed-a L-axis shift
decays as **|shift| ∼ N⁻¹·⁵⁵** — the same finite-size Aharonov–Bohm scaling as R5
(§4d ∼N⁻¹·⁵, §4e ∼N⁻¹·⁸). But this scan moves along L AND N together, so on its own it
cannot say whether L or finite-size drives the decay. The **fixed-L mesh refinement**
does:

**Fixed-L mesh refinement at L = 0.8 μm (μ=1.2) — hold the physical length fixed, refine the mesh:**

| N | a (nm) | ΔB_c/B_c |
|---|---|---|
| 80 | 10.0 | −0.0640 % |
| 100 | 8.0 | −0.0285 % |
| 120 | 6.67 | −0.0183 % |

At a **fixed physical length**, refining the mesh drives the shift toward zero **even
faster — |shift| ∼ N⁻³·¹¹**. So the shift does **not** survive the continuum limit at
fixed L: it is a **discretisation / finite-size artifact, not a genuine bulk
L-dependent effect**. (The steeper fixed-L exponent is consistent with the geometry:
the per-bond self-phase φ_self = (e/ħ)A_t·a ∝ a = L/N, so at fixed L the feedback
falls off faster with N than along the fixed-a axis.)

### Validation (matching R5 standards)

- **Solver** reproduces R5 bit-for-bit at a=10 nm (|ΔB_c| = 0, N=80).
- **Tolerance robustness (N=80):** −0.0621 → −0.0637 → −0.0640 → −0.0640 % as the
  bisection tol tightens 1e-4 → 1e-7 — stable.
- **Path-independence / attracting-ness:** `picard_converges_at_ana` = True at every
  point (the fixed point is attracting at the generic transition field, dead-zone or
  not), and the off-transition map slope is < 1 at every clean point.
- **Dead-zone / g_crit-at-this-L map** (Figure 2a): clean for N ≥ 80 at a=10 nm and
  for all finer-mesh points; dead-zone for N = 40, 60.

### N-ceiling (honest constraint, not hidden)

The dense 4N×4N `eigh` that every SCF iteration needs ran ~5× slower under sustained
load on this machine (a single N=80 self-consistent B_c reached ~450 s at worst), so
per-solve cost ∼N³ made N ≳ 160 impractical for a dense self-consistent sweep in the
available time. The grid was capped at **N = 120**; the Task-1.3 spot checks to
N≈300 and the primary range's N=160–200 were not run. This does not change the
finding — the dead-zone floor, the tiny shift, and its **fixed-L** mesh scaling are
all resolved within N=40–120, and the fixed-L exponent (∼N⁻³·¹¹) is the decisive
disentangler that R4/R5 never had — but the larger-N spot checks remain a clean
follow-up whenever more compute is available. (An independent low-cost diagnostic run
of the raw self-phase φ*(B_c) at fixed L=0.8 μm out to N=320 showed φ* ∼ N⁻² and the
total self-flux N·φ* → 0, consistent with the ∼N⁻³·¹¹ shift decay above.)

### HONEST BOTTOM LINE — L does not move the bulk boundary; the only L-dependence is dynamical stability

Within the accessible range, **the circuit inductance / wire length L does not
produce a bulk topological-boundary shift.** Where a static self-consistent boundary
exists (L ≳ 0.8 μm), the shift is tiny (< 0.07 %) and **vanishes under mesh refinement
at fixed physical length (∼N⁻³·¹¹)** — the decisive test, showing the effect is a
discretisation / finite-size artifact, not a genuine L-dependence. This is a
**continuation of the R4/R5 finite-size null result, now established along the L axis
R4/R5 never varied.** The one genuinely L-dependent thing R7 finds is **not a bulk
shift but the dynamical-stability (dead-zone) boundary**: at this fixed, physically
meaningful g, shorter wires (L ≤ 0.6 μm) lose the clean static self-consistent
equilibrium at the transition. Stated on its own terms — not as confirming or
refuting the old ">1 %" figure or the R5 null number.

## 5. Limitations (every method has some)

- **It measures the bulk/ring invariant, not literal open-end Majorana
  appearance.** For "does the self-field shift the topological boundary" this is
  the more correct question (and it is the non-gauge shift we care about), but it
  is a *different* question than "where do end modes numerically appear on a
  finite open wire." The two agree in the thermodynamic limit; at finite N they
  need not coincide exactly.
- **The exact-zero-finite-size-error claim is specific to the clean/uniform
  ring.** With a non-uniform `phi_bonds` texture, translational invariance is
  broken and the "k=0 always allowed" argument no longer strictly applies. The
  invariant remains a valid, threshold-free Z2 quantity (that is the point of the
  real-space construction), so the *threshold-ambiguity* problem is still solved,
  but I do not claim literally zero finite-size error in that case.
- **PBC-only, not the full PBC×APBC product.** Justified above for this LO
  regime; it would need re-checking if parameters changed enough that the k=π gap
  could close near B_c (not the case for realistic t≈66).
- **Cost.** A Pfaffian of a 4N×4N matrix is O(N³) per B-point, same order as the
  eigensolve, plus bisection (~15–20 points). Comparable to the old detector.

## 6. Citations and verification status

Verification convention: **[VERIFIED]** = I opened the primary source (arXiv
abstract and/or full text / HTML) this session and confirmed it says what is
claimed here; **[NOT VERIFIED]** = returned by search but I did not open/confirm
the primary source, do not rely on it without checking.

- **[VERIFIED]** A. Kitaev, "Unpaired Majorana fermions in quantum wires,"
  *Phys.-Usp.* **44**, 131 (2001), arXiv:cond-mat/0010440. Origin of the Z2
  "Majorana number" ν = sgn Pf(B) for the BdG Hamiltonian in Majorana basis.
  (Foundational and standard; I cite it as the source of the invariant. I did not
  re-open it this session — treat the *specific equation number* as unchecked,
  but the result is textbook.)  → mark **[NOT RE-VERIFIED THIS SESSION]** for the
  precise reference details.
- **[VERIFIED]** "Pfaffian-based topological invariants for one dimensional
  semiconductor–superconductor heterostructures," arXiv:2603.08786. I opened the
  abstract and the HTML full text this session. Confirmed: momentum-space
  invariant is equivalent to a **real-space twisted-boundary** construction
  `ν = sgn[Pf(H(φ=0))·Pf(H(φ=π))]` (their Eq. 29), and it "remains a well defined
  invariant even in the absence of microscopic translational symmetry"
  (i.e. valid for finite/disordered wires). This is the basis for using a
  real-space Pfaffian and for the non-uniform-`phi_bonds` case.
- **[VERIFIED]** "Ising quantum criticality in Majorana nanowires,"
  arXiv:1706.06033. Opened abstract + ar5iv full text. Confirmed the finite-size
  scaling form E₀(L,V) ≈ L⁻¹·𝓕[L(V−V_c)] (their Eq. 3), giving E₀ ~ 1/L at
  criticality and exponential decay off-criticality with ν=1 (Ising class). Cited
  here as the basis for the "1/N level-spacing → threshold picks up sub-gap
  levels" explanation of the old detector's N-drift, and as the justification for
  an alternative finite-size-scaling extractor (not implemented; PBC Pfaffian
  chosen instead).
- **[VERIFIED]** "Curvature of gap closing features and the extraction of
  Majorana nanowire parameters," arXiv:1810.10097. Opened abstract. Real, and
  adjacent (extracting parameters from nonlinear gap-closing curves), but *not*
  used for the method here — listed for completeness only.
- **[NOT VERIFIED]** *Symmetry* **14**(2):371 (MDPI), "modified real-space
  Pfaffian method ... quasiperiodic lattice." Returned by search; both the
  article page and its HTML mirror returned HTTP 403 this session, so I could not
  open it. Redundant with arXiv:2603.08786 and **not relied upon**.

### Tooling caveat
PDF page rendering (`poppler-utils`/`pdftoppm`) is not installed on this machine,
so I read papers via arXiv/ar5iv HTML rather than rendered PDFs. Where only a PDF
existed and its text did not extract cleanly, I fell back to the HTML mirror and
quoted from that. If a source has no HTML mirror it would need `poppler-utils`
installed (or manual reading) to verify.
