# R8-B follow-ups — adversarial-review tasks 1–5

_Frequency axis uses the R6-verified induced gap Δ = 0.19 meV (Chang et al. 2015; corrected from the superseded 0.25 meV InSb value). Only the GHz labelling depends on this; the dimensionless physics is unchanged._

## Task 1 — DC gate across the FULL grid (+ near-B_c stress)

The DC cross-check passes at **6/6 non-dead-zone grid points** (rel. mismatch ~1e-4). Dead-zone points (N=40, 60) are N/A — no static fixed point, gate undefined, consistent with Track A. **Near-B_c stress (N=80):** the gate holds all the way into the near-degenerate region as the gap closes — mismatch stays ~1e-4 and actually improves toward B_c:

  * B_c−0.100 (gap=0.0989): rel 8.0e-05, passes=True
  * B_c−0.050 (gap=0.0489): rel 6.0e-05, passes=True
  * B_c−0.020 (gap=0.0189): rel 4.9e-05, passes=True
  * B_c−0.010 (gap=0.0089): rel 4.5e-05, passes=True
  * B_c−0.005 (gap=0.0039): rel 4.3e-05, passes=True

So the gate does not quietly fail in the small-denominator region — the place it was most likely to.

## Task 2 — does the dissipation onset track the local gap? **No.**

Sweeping B across B_c=1.562 (2·gap closes and reopens, 0.24→0.004→0.24), the dominant Im χ peak **drifts monotonically and smoothly** (1.62→1.49→1.49 in ℏω/Δ) — a Zeeman shift of the BULK pairing edge, not the V-shaped topological minigap. The topological pair-breaking edge (2·minigap) carries only 0.0064–0.15% of the absorption weight. **Conclusion:** the finite-ω dissipation is a bulk-band property; it is NOT a topological-transition marker (if it tracked the minigap the peak would dip to ~0 at B_c and recover — it does not). This is a falsifiable negative, reported as such.

## Task 3 — is the η non-convergence finite-size or intrinsic? **Finite-size.**

Two lines of evidence:

* **Physics (decisive):** this is a closed, gapped BdG model with NO bath and NO quasiparticle-lifetime term anywhere in the Hamiltonian. The +iη is a pure numerical regularisation with no physical value, so the η→0 non-convergence CANNOT be intrinsic broadening — there is no mechanism for it.

* **Numerical confirmation (falsifiable):** the KK residual cannot discriminate η (broadening preserves causality, so over-smoothing still passes KK — this is why the naive η_c-minimisation saturated and was discarded). The right probe is the raggedness of Im χ at a FIXED small η: it **decreases monotonically as N grows** (0.229→0.110 over N=60→240, ∝ N^-0.54). The discrete spectrum fills in toward a smooth continuum as N→∞ — exactly the finite-size signature (cf. R5's Aharonov–Bohm artifact). No amount of N gives an intrinsic linewidth; larger N just needs smaller η. **So the dissipative spectrum is finite-size-limited, not a physical lifetime.**

## Task 4 — disposition of the omitted spin-orbit current

**Propagated as an explicit error band (not a footnote).** Including the SO current in the response lowers s and hence L_eff by **-4.0% to -2.3%** across the grid (SO current opposes the kinetic one). So every Track-A/8B reactance carries a **2–4% downward SO band**:

| point | s (kinetic) | s (kin+SO) | L_eff shift |
|---|---|---|---|
| mu1.200_N80_a10.000nm | 0.290 | 0.266 | -3.2% |
| mu1.200_N100_a10.000nm | 0.299 | 0.272 | -3.7% |
| mu1.200_N120_a10.000nm | 0.310 | 0.281 | -4.0% |
| mu1.400_N80_a10.000nm | 0.309 | 0.283 | -3.6% |
| mu1.200_N100_a8.000nm | 0.239 | 0.218 | -2.7% |
| mu1.200_N120_a6.667nm | 0.207 | 0.187 | -2.3% |

A fully self-consistent SO treatment (SO current also sourcing the Biot–Savart map, and re-solving φ*) is a separate model — scoped as **Track C**, the same way B was scoped out of A rather than silently folded in. The DC gate above is convention-matched to kinetic-only by construction, so this band does not affect the A↔B agreement.

## Task 5 — frequency-scale provenance

* **Induced gap:** corrected to the R6-verified Δ = 0.19 meV (Chang et al. 2015, arXiv:1411.6255, epitaxial-Al full shell). The earlier 0.25 meV was the InSb-sourced value R6 had already corrected; reusing R6's verified number rather than re-deriving. This sets the interband scale ℏω/Δ~1 ≈ 46 GHz.

* **Readout band:** the 1 GHz probe (ℏω/Δ=0.022) is a representative RF-reflectometry frequency; dispersive resonators for these platforms are typically 4–8 GHz. The whole realistic readout band (0.1–10 GHz ⇒ ℏω/Δ ≈ 0.002–0.2) sits far below the ~46 GHz interband scale, so the reactive-only result is robust for ANY realistic probe, not just 1 GHz.

* **Idealised-regime caveat (inherited from R6):** this Δ is combined with a Zeeman range that exceeds the verified ~75 mT axial critical field of a comparable thin Al shell (the 8.3× tension), and the model holds Δ fixed in B — the same documented limitation R6 stamps on every figure.

## Net effect on the R8-B verdict

The core result stands and is strengthened: the DC gate holds across the whole grid including the near-degenerate region (T1); the finite-ω dissipation is real Kubo physics but is (a) a bulk-band feature that does not track the topological transition (T2) and (b) finite-size-limited, not an intrinsic linewidth (T3); the reactance carries an explicit ~2–4% SO band (T4); and the frequency scales are now tied to R6-verified provenance (T5). **Z(ω) confirms Track A in the readout band and adds no trustworthy topological information there.**
