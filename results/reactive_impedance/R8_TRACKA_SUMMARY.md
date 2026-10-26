# R8 Track A — reactive impedance Z_r vs bulk gap: findings

_Quasi-static reactive-only Z_r. L_eff cross-VALIDATED by the R8-B finite-ω Kubo impedance: its ω->0 limit reproduces this L_eff to <1e-4 at every clean grid point (DC-limit gate passed)._

## Definition (confirmed, not assumed)

`Z_r = iωL_eff`, `L_eff = L_geo/(1−s)`. **`s` is confirmed to be R5's bifurcation slope** (bit-exact: 0.2897 at g=2.9e4 == cached `slope_offtransition[0]`; 3.00 at g=3e5 == README "map slope ≈3"). `L_geo` from the R5-validated ring kernel (analytic cross-check < 5e-5). Both reused with **zero new eigendecompositions** for the grid; `Z_r` is purely reactive (quasi-static limitation).

## Side-by-side grid (reused from cache)

| point | regime | Δ-shift % (gap obs.) | s | L_geo (pH) | 1/(1−s) | L_eff (pH) | X_r (mΩ) | flag |
|---|---|---|---|---|---|---|---|---|
| mu1.200_N40_a10.000nm | deadzone | nan | nan | 0.329 | — | — | — | undefined_no_fixed_point |
| mu1.200_N60_a10.000nm | deadzone | nan | nan | 0.543 | — | — | — | undefined_no_fixed_point |
| mu1.200_N80_a10.000nm | clean | -0.064 | 0.274 | 0.770 | 1.377 | 1.060 | 6.66 | ok |
| mu1.200_N100_a10.000nm | clean_robustified | -0.043 | 0.285 | 1.007 | 1.399 | 1.408 | 8.85 | ok |
| mu1.200_N100_a8.000nm | clean | -0.029 | 0.228 | 0.805 | 1.295 | 1.043 | 6.55 | ok |
| mu1.200_N120_a10.000nm | clean | -0.034 | 0.297 | 1.252 | 1.422 | 1.780 | 11.19 | ok |
| mu1.200_N120_a6.667nm | clean | -0.018 | 0.198 | 0.835 | 1.247 | 1.041 | 6.54 | ok |
| mu1.400_N80_a10.000nm | clean | -0.065 | 0.292 | 0.770 | 1.413 | 1.088 | 6.83 | ok |

## Where they agree, where they diverge

* **Across the topological transition (B-sweep, N=120, μ=1.2, L=1.20 μm):** the gap closes cleanly to ~0 at B_c=1.5615, an unambiguous sharp marker. `X_r(B)` is **nearly flat**: 11.42 mΩ (trivial) → 10.98 mΩ (topological), a **-3.8% step** — a weak feature, not a sharp signature. So for *detecting the topological boundary at g=2.9e+04, Z_r is largely redundant with (and dominated by) the geometric inductance* L_geo.

* **Dead-zone agreement:** the N=40, N=60 (μ=1.2) points are regime `deadzone` — gap-shift `nan` AND `Z_r` undefined. The two observables agree that there is no static equilibrium there.

* **Where Z_r carries orthogonal information (structure the gap cannot see):** the g-sweep shows `1/(1−s)` growing steeply (1.11 at g=10⁴ → 4.56 at g=8×10⁴); the static uniform fixed point then **detaches at g≈9.0e+04** (last stable point g≈8.0e+04), i.e. `L_eff` **diverges as s→1** at the leading edge of R5's g_crit band (1.1–1.2×10⁵). The clean bulk gap is **g-independent** — it has no such pole. Z_r's distinctive content is proximity to the *dynamical* bifurcation, not the topological one.

## Finite-size discipline (N-scaling)

Same discipline as R5/R7, applied to the enhancement 1/(1−s):

* **L-axis (a=10 nm; L & N grow together):** enhancement RISES with N (1.379→1.506) — confounded (both L and mesh change), exactly like R5/R7's L-axis boundary shift.

* **Fixed-L mesh refinement (L=0.8 μm, a=L/N, N=80→240):** enhancement FALLS monotonically (1.408→1.315→1.260→1.193→1.155→1.130) as the mesh refines. The slope s is very nearly proportional to the lattice constant (s/a ≈ 0.0290→0.0346 /nm, roughly constant), so **s → 0 and the enhancement → 1 (no renormalization) in the continuum limit.** Like R7's boundary shift, the reactive renormalization is dominated by the lattice discretisation and does not survive mesh refinement — a finite-size effect, not a bulk one.

* **Modelling caveat (inherited from R7):** this refinement holds t_hop at its 10 nm value and refines only the *geometric* a in the self-field (R7's stated choice). A fully dimensionful continuum limit (t ∝ 1/a²) is a separate, heavier computation; the s→0 trend is under the same fixed-t modelling choice R7 used, and carries the same caveat.

## Honest verdict — does Z_r add information beyond the bulk gap?

In the validated production regime (g=2.9e+04), **largely no — it is mostly redundant with the gap plus geometry.** (1) Across the topological transition Z_r is nearly featureless: a ~4% step in X_r where the gap closes cleanly to zero, and its magnitude is set mainly by the geometric inductance L_geo — so as a phase-boundary detector it is weaker than, and redundant with, the gap. (2) Its feedback-renormalization enhancement 1/(1−s) is itself a finite-size effect: under fixed-L mesh refinement s ∝ a → 0 and the enhancement → 1, the same artifact-like vanishing R7 found for the boundary shift (same fixed-t caveat). The **one** genuinely orthogonal thing Z_r carries is *dynamical, not topological*: L_eff = L_geo/(1−s) diverges as s→1 (here at g≈8–9×10⁴, the leading edge of R5's g_crit band, where the static fixed point detaches), whereas the clean bulk gap is g-INDEPENDENT and structurally cannot see this. So Z_r is a redundant topological marker but a non-redundant *proximity-to-the-dynamical-bifurcation* marker — and that added value only appears when g is pushed toward g_crit, not at the production coupling. This verdict is Track-A (quasi-static, reactive-only); its L_eff is now **cross-validated** — the R8-B finite-ω Kubo impedance reduces to exactly this L_eff as ω→0 at every clean grid point (DC-limit gate passed, rel. mismatch ~1e-4), so the PRELIMINARY caveat is lifted.
