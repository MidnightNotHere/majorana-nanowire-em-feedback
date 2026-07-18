# R7 logbook summary — circuit inductance / wire length L as a third axis (B, μ, L)

**Module:** `src/bdg_r7_inductance_sweep.py`  ·  **Methods:** `src/r4_pfaffian_boundary_methods.md` §4g
**Figures:** `r7_figure1_shift_vs_BmuL.png`, `r7_figure2_regime_and_validation.png`  ·  **Data:** `r7_results.json`, `r7_points_cache.json`

## What was asked
Hypothesis-neutral re-opening of the inductive boundary-shift question along the
one axis R4/R5 never varied — the wire length **L = N·a_lat** (which sets the loop
self-inductance) — at the single validated coupling **g = 2.9×10⁴** (top of R5's
regime 1) and fixed **a_lat = 10 nm**. No target result; report what the sweep shows.

## Design (how L is disentangled from finite-size N)
Because the ring circumference is exactly N·a_lat, L and the finite-size N are locked
in R5's original scan. R7 uses two orthogonal scans:
1. **L-axis** — fix a=10 nm, vary N ⟹ L = N·a (moves along L AND N together).
2. **Fixed-L mesh refinement** — hold L fixed, vary N with a = L/N (finer mesh of the
   SAME physical ring). This is the decisive test R4/R5 never ran: does the shift at a
   fixed physical length survive the continuum limit, or vanish?

Solver `ring_Bc_selfconsistent_a` reproduces R5's `ring_boundary_Bc_pfaffian_selfconsistent`
**bit-for-bit at a=10 nm** (|ΔB_c| = 0). Every point carries an N-scaling / regime check.

## What was found
**1. An L-dependent DYNAMICAL-STABILITY boundary (new).** At g=2.9×10⁴ the static
self-consistent boundary is cleanly defined only for L ≳ 0.8 μm (N ≳ 80); **L = 0.4
and 0.6 μm (N = 40, 60) are regime-2 dead-zones** — no static B_c. Since L is set by N,
the regime-1/2 boundary R5 saw at one N is really an L-dependent stability threshold
(Task 2.3). This is a stability effect, not a bulk shift: Picard iteration still
converges at the generic transition field everywhere (dead-zone = the thin gap-closing
window, not a repelling fixed point).

**2. Where the boundary exists, the shift is tiny and vanishes under fixed-L refinement.**

| L-axis (a=10 nm) | L (μm) | ΔB_c/B_c μ=1.2 | μ=1.4 |
|---|---|---|---|
| N=80 | 0.8 | −0.0640 % | −0.0647 % |
| N=100 | 1.0 | −0.0430 % | — |
| N=120 | 1.2 | −0.0342 % | — |

fixed-a decay: |shift| ∼ N⁻¹·⁵⁵ (same finite-size AB scaling as R5).

| Fixed-L refine, L=0.8 μm | a (nm) | ΔB_c/B_c |
|---|---|---|
| N=80 | 10.0 | −0.0640 % |
| N=100 | 8.0 | −0.0285 % |
| N=120 | 6.67 | −0.0183 % |

**fixed-L decay: |shift| ∼ N⁻³·¹¹** — at a fixed physical length the shift vanishes
under mesh refinement (even faster than along the fixed-a axis). ⟹ the shift is a
**discretisation / finite-size artifact, not a bulk L-dependence.**

**3. Validation.** Solver = R5 to |ΔB_c| = 0; tolerance-robust (−0.0621→−0.0640 % as
tol tightens 1e-4→1e-7); `picard_converges_at_ana` = True at every point.

## Bottom line
Within the accessible range, **the circuit inductance / wire length L does not move
the bulk topological boundary** — where a static boundary exists the shift is < 0.07 %
and vanishes under fixed-L mesh refinement (∼N⁻³·¹¹). This is a **continuation of the
R4/R5 finite-size null result, now along the L axis R4/R5 never varied.** The single
genuinely L-dependent effect is **not a bulk shift but the dynamical-stability
(dead-zone) boundary**: short wires (L ≤ 0.6 μm) lose the clean static equilibrium at
this fixed g.

## Caveat / follow-up
Per-solve cost is dense-`eigh`-bound (∼N³) and ran ~5× slower under load on this
machine, so the grid was **capped at N=120**; the Task-1.3 spot checks to N≈300 and
the primary range's N=160–200 remain a clean follow-up (a low-cost φ*(B_c) diagnostic
out to N=320 already showed φ* ∼ N⁻² and total self-flux → 0, consistent with the
∼N⁻³·¹¹ shift decay). L range computed: 0.4–1.2 μm (intended 0.4–2.0 μm); a
tractability subrange, not a literal device length.
