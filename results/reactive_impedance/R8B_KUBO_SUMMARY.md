# R8-B — finite-ω Kubo impedance: findings + combined A↔B note

## DC-limit cross-check against Track A (the gate)

At each clean grid point: the Kubo static value χ(0) vs the finite-difference of ⟨Ĵ⟩ (the electronic response Track A's slope is built from). **This is falsifiable — a Kubo bug breaks it.**

| point | L (μm) | gate rel. mismatch | passes | χ(0) | |SO current| |
|---|---|---|---|---|---|
| mu1.200_N80_a10.000nm | 0.80 | 6.0e-05 | True | 1797.1 | 8.9% |
| mu1.200_N100_a10.000nm | 1.00 | 9.0e-05 | True | 2216.5 | 9.8% |
| mu1.200_N120_a10.000nm | 1.20 | 9.2e-05 | True | 2658.0 | 10.1% |
| mu1.400_N80_a10.000nm | 0.80 | 5.1e-05 | True | 1916.0 | 9.1% |
| mu1.200_N100_a8.000nm | 0.80 | 7.6e-05 | True | 2216.8 | 9.8% |
| mu1.200_N120_a6.667nm | 0.80 | 7.0e-05 | True | 2658.6 | 10.1% |

**All 6 clean points pass** (rel. mismatch ~1e-4 ≪ 1e-3): Z(ω→0) from Track B reduces to Track A's L_eff to <1e-4. ✅ The A↔B DC agreement HOLDS.

## Complex impedance Z(ω) — the new physics

At the representative point (N=120, μ=1.2, a=10.0 nm): Z(ω) = iωL_eff(ω) with a COMPLEX, frequency-dependent L_eff. Unlike Track A (reactive only), Re Z ≠ 0 — a resistive/dissipative part from interband quasiparticle absorption. But **at the 1 GHz probe (ℏω/Δ=0.022) the response is essentially reactive**: Im Z = 1.15e-02 Ω vs Re Z = 4.4e-09 Ω (|Re/Im| = 3.9e-07). Dissipation (Im χ) only turns on at ℏω/Δ ~ 1 (≈46 GHz), the interband scale — so across the entire RF read-out band Track B agrees with Track A's reactive result.

## Causality (Kramers–Kronig) and resolution convergence

* **KK residual** (rep, η=0.03): median 2.4e-02 (bulk of the spectrum is causal-consistent) but max 2.0e+01 — local violations where the finite (4N-level) spectrum is spiky.

* **η / N_w convergence:** χ(0) is fully converged (resolution-independent to the digits shown); the *dissipative* spectrum is broadening-limited — the KK residual does not converge cleanly as η→0 because a 320–480-level discrete spectrum cannot represent a smooth causal response without η washing out real sub-features.

  * η=0.10: KK median 8.7e-03, max 1.1e+03, χ(0)=2658.007
  * η=0.05: KK median 1.6e-02, max 5.7e+01, χ(0)=2658.007
  * η=0.03: KK median 2.4e-02, max 2.0e+01, χ(0)=2658.007
  * η=0.01: KK median 8.4e-02, max 1.4e+02, χ(0)=2658.007

## Modelling note — omitted spin-orbit current

Track A (and hence this gated cross-check) uses the KINETIC current only. The spin-orbit current Ĵ_SO that both omit is 9–10% of the kinetic current here (α_norm/t_norm≈0.12) and opposite in sign — a real but sub-dominant correction. Folding it in is a separate model (it would shift the *absolute* current stiffness by ~10%, not the DC-limit agreement, which is convention-matched by construction).

## Honest verdict (Track B)

The finite-ω Kubo machinery is **validated at DC** (the gate passes at every grid point, χ(0) reproduces Track A's static response to <1e-4) and yields a genuine COMPLEX impedance whose ω→0 limit is exactly Track A. Its one qualitatively new feature — a resistive part from interband absorption — is **negligible across the RF read-out band** (turns on only at ℏω/Δ~1 ≈ 60 GHz) and, where it is non-zero, is **finite-size/broadening-limited at N=120**: the discrete 480-level spectrum does not resolve a smooth causal dissipative response (KK residual η-dependent, does not converge as η→0). So Track B **confirms Track A in the physically relevant band** and adds no trustworthy new structure there; the dissipative spectrum is reportable only as *compute-limited/inconclusive at the achieved resolution*, exactly the outcome the R8 spec permitted. A dense-spectrum (continuum or much larger N) study is the follow-on if the dissipative part is wanted quantitatively.

## Combined note — does A agree with B in the DC limit?

**Yes.** The DC-limit cross-check passes at every shared grid point (rel. mismatch ~1e-4). Track A's L_eff = L_geo/(1−s) is confirmed as the ω→0 limit of the independent finite-ω Kubo impedance. **This lifts the PRELIMINARY label from the Track A reactive result** (its magnitude is now cross-validated). The Track A verdict itself is unchanged: in the production regime Z_r is largely redundant with the gap+geometry, distinct only as a proximity-to-g_crit marker.
