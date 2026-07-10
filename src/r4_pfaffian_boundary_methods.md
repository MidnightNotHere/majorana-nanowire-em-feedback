# Methods note: N-independent extraction of the topological critical field B_c via a Pfaffian Z2 invariant

**Project:** ISEF 2027 — Inductive Circuit Feedback Shifts the Topological Phase Boundary in InAs/Al Majorana Nanowires
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
**Outputs:** `ISEF R5 data folder/` (Figures 7–9, `r5_results.json`,
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
