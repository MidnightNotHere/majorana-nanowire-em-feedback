# R8-C — self-consistent spin-orbit current: findings

## SO current cross-check (vs R8-B Task 4)

Track C sources the self-field from the TOTAL current I_kin·(⟨Ĵ_full⟩/⟨Ĵ_kin⟩), reusing R8-B's VALIDATED kinetic and kinetic+SO current operators (exact for the uniform texture the map uses; the constant Ampere-scale cancels, R8-B Check 1). The SO fraction I_SO/I_kin = -8.2% (so_fraction = -8.9%) reproduces R8-B's operator-side `so_current_fraction` (−8.9 %) by construction.

## Track A (kinetic) vs Track C (total current) across the grid

| point | s_A | s_C | L_eff A (pH) | L_eff C (pH) | shift |
|---|---|---|---|---|---|
| μ1.2_N80_a10.0nm | 0.2897 | 0.2659 | 1.084 | 1.049 | -3.24% |
| μ1.2_N100_a10.0nm | 0.2991 | 0.2725 | 1.436 | 1.384 | -3.67% |
| μ1.2_N100_a8.0nm | 0.2393 | 0.2180 | 1.059 | 1.030 | -2.73% |
| μ1.2_N120_a10.0nm | 0.3097 | 0.2812 | 1.814 | 1.741 | -3.98% |
| μ1.2_N120_a6.7nm | 0.2066 | 0.1875 | 1.052 | 1.027 | -2.34% |
| μ1.4_N80_a10.0nm | 0.3089 | 0.2831 | 1.114 | 1.074 | -3.60% |

**Self-consistent SO shifts L_eff by -4.0% to -2.3%** across the grid (SO current opposes the kinetic, lowering the feedback gain s and hence L_eff). This is consistent with R8-B Task 4's linear estimate (−2.3 to −4.0 %): making the SO current self-consistent (re-solving φ*) does **not** qualitatively change the magnitude — the SO correction is small and nearly linear here (φ* ~ 1e-3, so the fixed-point feedback on top of the ~8 % source change is second order).

## Finite-size discipline

* **Fixed-L mesh refinement (L=0.8 μm, a=L/N):** the SO L_eff shift **shrinks monotonically** -3.24→-2.73→-2.34→-1.82→-1.49→-1.26% (N=80→240) as the mesh refines. It does NOT converge to a finite bulk correction — it trends toward 0 in the continuum, the SAME finite-size character as Track A's enhancement and R7's boundary shift. (Mechanism: both s_A and s_C scale ∝ a, so their difference s_C−s_A ∝ a → 0.) So the SO current is not a hidden bulk contribution the kinetic-only model missed; its reactance effect vanishes in the continuum limit just like the kinetic self-field renormalization itself.

## Honest verdict

Making the spin-orbit current self-consistent (Track C) lowers the feedback-renormalized reactance by ~2–4 % across the L-axis grid, **confirming — not overturning — R8-B Task 4's linear error band** (the self-consistent and linear numbers agree to ~0.01 %, because φ* ~ 1e-3 makes the fixed-point feedback on the ~8 % source change second-order). Two clean conclusions: (1) the SO correction is **small and sub-dominant** — Z_r remains dominated by the geometric inductance; (2) it is **itself a finite-size effect** — it shrinks under fixed-L mesh refinement and → 0 in the continuum, so it is not a bulk term the kinetic-only Tracks A/B were missing. The self-field is only rescaled ~8 %, so the (already null, finite-size) bulk topological boundary shift is unchanged. Net: Track C promotes the SO current from a quantified caveat to a self-consistently included correction, and shows it changes nothing qualitative — the R8 verdict (Z_r largely redundant with gap + geometry) stands.
