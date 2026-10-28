# Correction note: `pfaffian_trim_lo()` discarded its own TRIM matrices

**Project:** Inductive Circuit Feedback Shifts the Topological Phase Boundary in InAs/Al Majorana Nanowires
**Module:** `src/bdg_eigensolver_r1_r2_LOAdapted.py` (R1/R2/O1)
**Found:** 2026-07-11, during R6 reconnaissance (pre-implementation code audit)
**Fixed:** 2026-07-11, standalone corrective task (separate from R6)

## 1. What was wrong

`pfaffian_trim_lo()` computes the Lutchyn-Oreg model's Z₂ topological
invariant at the two TRIM (time-reversal-invariant momentum) points k=0 and
k=π. Its docstring described building real 4×4 antisymmetric TRIM matrices
in the Majorana basis and taking their Pfaffian — the same recipe used by
the (genuine) `pfaffian_trim_kitaev()` and by R4's `_ring_pfaffian_sign()`.

The code did build those matrices (`_trim_matrix()`, called for both k=0
and k=π) — but then **discarded them** and returned a hardcoded closed
form instead:

```python
A_k0  = _trim_matrix(2.0 * t_trim - mu)     # built...
A_kpi = _trim_matrix(-2.0 * t_trim - mu)    # ...and built...

pf_k0  = mu**2 + delta**2 - B_zeeman**2     # ...then never used. This is
pf_kpi = pf_k0                              # the closed-form answer, not
                                             # a Pfaffian of anything.
nu = float(np.sign(pf_k0))
```

`pf.pfaffian()` (the actual numerical routine) was never called. The
matrices `A_k0`/`A_kpi` had no effect on the return value.

## 2. How it was found

During R6 reconnaissance (extending the codebase to a full Lutchyn-Oreg
model for a published-phase-diagram comparison), a line-by-line read of
`build_bdg_lo()` and its neighboring functions — not just their
docstrings — turned up the discard. A follow-up audit was then run
specifically to determine whether this bug had contaminated any
already-reported R3/R4/R5 result.

**Call-site audit (completed before this fix, not repeated here):**
`pfaffian_trim_lo()` is called from exactly two places, both inside this
same file — `_validate_lo_phase_boundary()` (an assertion that runs at
module import) and `run_lo_phase_sweep()` (populates the ν-panel of
`results/baseline_bdg_pfaffian/figure6_lo_phase_diagram.png`, a **different file** from
R4's `results/em_feedback_detector_comparison/r4_figure6_pfaffian_vs_threshold.png` — the two
"figure 6"s should not be confused). R3 imports the function but never
calls it. **No R4 or R5 code calls it, directly or indirectly.** R4/R5's
Pfaffian-based results (`ring_boundary_Bc_pfaffian`,
`ring_boundary_Bc_pfaffian_selfconsistent`, the Anderson-mixing
cross-validation, Figures 6–13 in the R4/R5 data folders) all trace to
`_ring_pfaffian_sign()` in `bdg_r4_em_coupling.py`, a separate
implementation that always called the real `pf.pfaffian()` on the actual
(dressed) ring Hamiltonian and was validated on its own terms (known-answer
check at μ=0, N-independence across N=20–160, non-round-μ checks,
tolerance/bisection robustness — see `r4_pfaffian_boundary_methods.md`
§4). **The R4/R5 boundary-shift lineage, the g_crit finding, and the
three-solver (bisection / Newton–Krylov / Anderson) agreement are
unaffected by this bug.**

**Consequence for this file specifically:** the `_validate_lo_phase_boundary()`
assertion, and the printed claim "`[O1 LO] Phase boundary validated:
B_z = sqrt(μ²+Δ²)`", were **tautological** — the closed-form value was
checked against 0.5×/1.5× of itself, which passes by construction
regardless of whether the actual LO Hamiltonian's topology matches that
formula. This is a real overstatement in the module header (line ~27),
the changelog (`[12]`, and the printed `VALIDATION SUMMARY`), and the
corresponding README status line — corrected in this pass (see §5).

**`pfaffian_trim_kitaev()` was checked separately and found NOT to have
this problem.** It genuinely calls `pf.pfaffian(A_k0)` / `pf.pfaffian(A_kpi)`
and returns their real parts — verified directly (not just read) by tracing
`pfaffian_trim_kitaev(t, mu, delta)`'s returned `pf_k0`/`pf_kpi` numerically
against the actual matrix entries `eps_k0 = -2t-mu` / `eps_kpi = 2t-mu` at
several `(t, mu)` combinations, and confirming `pf.pfaffian([[0,a],[-a,0]])
= a` independently. **Scope of this specific check, stated precisely:** for
a 2×2 antisymmetric matrix, `Pf([[0,a],[-a,0]]) = a` is a trivial identity
(the Pfaffian of a 2×2 antisymmetric matrix is just its one independent
entry), so this trace confirms **the `pf.pfaffian()` call is genuinely
wired into the return value** — i.e. rules out the exact discard-and-hardcode
failure mode found in the LO function — but it does **not** by itself
validate that the Kitaev TRIM matrix construction is *physically* correct.
That is a separate, non-trivial claim, and it was already covered — before
this fix, unaffected by it — by `_validate_pfaffian_vs_analytic()` (asserted
at module import), which checks `pfaffian_trim_kitaev`'s output `ν` against
the independently-derived analytic Kitaev boundary `|μ|=2t` across a
`mu` sweep. These are two different checks for two different things: the
trace above shows the *mechanism* is real (no discard bug); the pre-existing
`_validate_pfaffian_vs_analytic` shows the *physics* is right. Both hold,
for different reasons — neither is being conflated with the other here.
**No fix was needed for the Kitaev function.**

## 3. What was fixed

`pfaffian_trim_lo()` now performs a genuine Pfaffian computation at both
TRIM points, matching `build_bdg_lo()`'s own real-space conventions exactly
(rather than reusing the old function's disconnected `t_trim=1` toy
matrices):

- At a TRIM point k ∈ {0, π}, `sin(k)=0`, so the Rashba term
  `-α·sin(k)·σ_y` vanishes identically and the 2×2 normal-state Bloch
  block collapses to `h(k) = ε_k·I₂ + B_zeeman·σ_z` (real, diagonal), with
  `ε_k = (2t-μ) - 2t·cos(k)` — the same on-site+hopping dispersion
  `build_bdg_lo()` assembles in real space. `ε_{k=0} = -μ` **exactly**,
  independent of `t` (since `2t(1-cos 0)=0`); `ε_{k=π} = 4t-μ`.
- Because `h(k)` is real at TRIM, the PH-symmetric hole block `-h(-k)*`
  collapses to `-h(k)`, giving the 4×4 BdG matrix at that point (basis
  `↑,↓,↑†,↓†`):
  ```
  H_TRIM = [[ ε_k+B,   0,        0,        Δ       ],
            [ 0,       ε_k-B,   -Δ,        0        ],
            [ 0,      -Δ,       -(ε_k+B),  0        ],
            [ Δ,       0,        0,       -(ε_k-B)  ]]
  ```
- This is rotated to the real antisymmetric Majorana basis with the same
  unitary R4's `_ring_pfaffian_sign()` uses (restricted to a single
  momentum, M=2 complex modes): `Ω = (1/√2)[[I₂,I₂],[iI₂,-iI₂]]`,
  `A = i·Ω·H_TRIM·Ω†`. Antisymmetry/realness of `A` are asserted at
  runtime (as in R4's construction).
- `pf.pfaffian(A.real)` is called for both k=0 and k=π (previously: never
  called at all). `ν = sgn[Pf(k=0)] · sgn[Pf(k=π)]` — the same convention
  as `pfaffian_trim_kitaev`.
- `t` is now an explicit parameter (default `LO_PARAMS['t_norm']*delta`,
  matching `build_bdg_lo`'s own default) rather than a disconnected
  hardcoded `t_trim=1`. `run_lo_phase_sweep()` was updated to pass the same
  `t=2.0` it already uses for the diagonalized gap map, so the ν-panel and
  gap-panel of `figure6_lo_phase_diagram.png` describe the identical
  Hamiltonian (previously accidentally consistent only because the old
  formula ignored `t` entirely).
- The unconditional debug `print()` inside the old function (which would
  have fired 1,600 times during `run_lo_phase_sweep`'s default 40×40 grid)
  was removed, matching `pfaffian_trim_kitaev`'s style.

## 4. Re-validated result — was this an open question? Yes, and it resolved cleanly

Before the fix, whether the genuine LO Pfaffian invariant actually flips at
`B_z = sqrt(μ²+Δ²)` was **not known** — the closed form had never been
checked against anything real. This was tested empirically, not assumed:

- **Full grid used by `run_lo_phase_sweep()`** (μ∈[-2,2]×40, B_z∈[0,1.8]×40,
  1,600 points), at the sweep's own toy `t=2.0`: **0 disagreements** between
  the genuine invariant and the old closed-form sign, at any point.
- **Range used by `_validate_lo_phase_boundary()`** (μ~U(-1.5,1.5)×20,
  tested at 0.5×/1.5× the analytic boundary, up to B≈2.70), at the
  function's default physical `t_norm≈66.3`: **0 disagreements**.
- **Mechanism — an exact algebraic identity, not a numerical coincidence.**
  For this `H_TRIM`/`Ω` construction, the Pfaffian at a TRIM point reduces
  in closed form to `Pf = B² - ε_k² - Δ²` — confirmed to hold to floating-point
  precision (max relative error 1.5×10⁻¹⁴) across 400 random `(ε_k, B, Δ)`
  triples spanning `ε_k∈[-500,500]`, `B∈[-50,50]`, `Δ∈[-20,20]`, i.e. this is
  an exact identity of the construction, not a fit to the sampled points.
  Consequently `pf_k0 = B²-μ²-Δ²` and `pf_kπ = B²-(4t-μ)²-Δ²`. Whenever
  `pf_kπ<0` — i.e. `(4t-μ)² > B²-Δ²` — the invariant collapses **exactly**:
  `ν = sgn(pf_k0)·sgn(pf_kπ) = sgn(B²-μ²-Δ²)·(-1) = sgn(μ²+Δ²-B²)`,
  identical to the old closed form, term for term, not merely matching in
  sign on average. That condition holds throughout every range used in this
  codebase with margin to spare: at toy `t=2.0`, `μ∈[-2,2]` gives
  `|4t-μ|=|8-μ|≥6`, while `sqrt(max(B²-Δ²,0))≈1.50` at the sweep's max
  `B_z=1.8`; at physical `t≈66.3` the margin is far larger still. **This is
  why the two invariants agree — not an emergent property of the specific
  grid scanned, but a provable consequence of `|4t-μ|` staying well above
  `B` and `Δ` for every parameter choice this codebase makes.**
- **A real deviation regime exists, and was located** (not merely
  asserted): scanning `t` downward at fixed extended range (μ∈[-3,3],
  B∈[0,4]), the k=π TRIM only starts flipping for `t ≲ 1` (e.g. at `t=1.0`,
  814/4,941 grid points already show a flipped `Pf(k=π)`). **No value of
  `t` used anywhere in this codebase (`t=2.0` toy, `t_norm≈66.3` physical)
  is anywhere near that regime.**
- **Direct regeneration confirms this at the artifact level, not just in a
  side computation:** re-running `bdg_eigensolver_r1_r2_LOAdapted.py`
  end-to-end regenerated `results/baseline_bdg_pfaffian/figure6_lo_phase_diagram.png`,
  `figure7_lo_finite_size.png`, `lo_nu_map.npy`, and
  `lo_finite_size_scaling.csv` — `git diff` against the pre-fix committed
  versions of every one of these is **empty** (byte-identical). The ν-panel
  of Figure 6 is pixel-identical to the pre-fix version.

  **Why the gap panel (right half of Figure 6) is unaffected — by
  construction, not by luck:** this fix changed two things — (1) how `ν` is
  computed (formula → genuine Pfaffian), and (2) made `pfaffian_trim_lo`
  accept an explicit `t`, which `run_lo_phase_sweep` now passes explicitly
  (`sweep_t=2.0`) to keep it matched to the diagonalized Hamiltonian.
  Neither change touches the gap-panel computation itself: `build_bdg_lo`
  was already being called with the literal `t=2.0` before this fix; the
  edit only renamed that literal to a shared `sweep_t` variable passed to
  *both* calls, so `build_bdg_lo` receives the identical numeric argument,
  in the identical code path, before and after. The gap map's
  byte-identical output is therefore not something that needed checking —
  no line feeding it changed. The ν-panel's byte-identical output *is* the
  non-trivial confirmation, and it is explained by the exact identity
  above, not merely by both `t`s "happening" to match: the identity holds
  for a wide range of `t` (anything with `|4t-μ|` large relative to `B,Δ`),
  so even a `t` mismatch between the two panels would very likely not have
  changed the sign map in this parameter regime — matching `t` was done for
  physical-consistency correctness (both panels should describe one
  Hamiltonian), not because it was required to reproduce this particular
  empirical result.

**Honest conclusion:** `B_z = sqrt(μ²+Δ²)` is now a **genuinely validated**
result for this codebase's actual usage — it was checked against a real
TRIM Pfaffian, not merely assumed, and it holds with zero exceptions found.
It was not a foregone conclusion going in (the closed form's independence
from `t` was a coincidence of the bug, and a real Pfaffian calculation
could in principle have disagreed, as it does demonstrably at `t≲1`). No
parameter was tuned to force this agreement — the check was run once, at
the codebase's existing parameter choices, and it passed.

## 5. Documentation corrected

- Module header (`~line 27`) and changelog entry `[12]` (`~line 50`) in
  `src/bdg_eigensolver_r1_r2_LOAdapted.py`: marked `[CORRECTED 2026-07-11]`
  with a pointer to this note, rather than silently rewritten.
- `pfaffian_trim_lo()`'s docstring: rewritten to describe what the function
  now actually does (§3 above), including the re-validation result and the
  known deviation regime.
- `_validate_lo_phase_boundary()`'s docstring: notes the assertion was
  previously tautological and is now genuine.
- `VALIDATION SUMMARY` print block in `__main__`: updated to flag the
  correction explicitly rather than repeat the old unqualified claim.
- **README:** as of this fix, `README.md` does not contain a separate "O1"
  status bullet making this specific claim — its closest line (`✅ R2:
  Pfaffian topological indicator (finite-size scaling confirmed)`) refers
  to the genuinely-validated Kitaev-chain result and does not itself
  overstate the LO boundary, so it was left as-is. No README text required
  correction for this specific issue.

## 6. Scope

Per the task constraints: this is a standalone fix to already-committed
O1/LO code. **No R3, R4, or R5 code, figures, or results were touched** —
confirmed by `git status` showing only `src/bdg_eigensolver_r1_r2_LOAdapted.py`
changed (plus this new file) after reverting two incidentally-regenerated,
unrelated Kitaev performance-benchmark artifacts whose values are expected
to fluctuate run-to-run with machine load and carry no scientific content
(`results/baseline_bdg_pfaffian/figure3_performance.png`,
`results/baseline_bdg_pfaffian/performance_baseline.csv` — reverted to their pre-run
committed versions). No R6 Phase 2 implementation work was started.
