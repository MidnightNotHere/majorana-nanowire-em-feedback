# R9 Track 4 — multi-subband loop-gain probe: methods

**Status of THIS section:** pre-registration. It documents the locked protocol
(`bdg_r9_multisubband.py`, first commit) — the van Hove exclusion radius δ, the
subband energetics, the sweep windows, the operating points, and the comparison
statistic — together with the mandatory construction gates. **No s(μ;M) value is
computed at this commit.** The s-sweep functions and every s(μ;M) number are
added in a separate, later commit, so the choices below cannot have been
influenced by any outcome.

This is the same discipline the project applied after R5's g = 3×10⁵ headline
turned out to be a dynamically-unreachable root chosen before the reachability
check existed: fix the analysis before seeing the number.

## 0. What question, and why the loop gain s

R1–R9 converged on a null — inductive EM self-feedback does not move the bulk
topological boundary — computed throughout with **exactly one** conducting
channel. Track 4 asks whether that null survives more channels. The scoping
report (`r9_multisubband_scoping.md`) established that the tractable form of this
question is not the full self-consistent B_c sweep (34.7–160.6+ h per
configuration at M ≥ 3 — out of reach here) but R5's loop-gain slope

> s ≡ dG/dφ at the self-consistent fixed point, where G(φ) = g·⟨Peierls phase
> from the Biot–Savart self-field of the ring current I(φ)⟩,

the same quantity behind R5's bifurcation, R8's L_eff = L_geo/(1−s), and R9-T1's
bounded enhancement 1/(1−s). At the production coupling g = 2.9×10⁴, R9-T1
measured s ≈ 0.29 (single channel), with the static fixed point detaching near
g ≈ 8.67×10⁴ > g_production. **Hypothesis (scoping §6.1):** each occupied channel
carries current, so s might grow ~∝ M; if so, at M ≈ 10 the detachment would
drag to g ~ 10⁴ ≤ g_production and the single-channel null would be an artifact.
If s stays flat in M, the null is *stronger* than a single-channel model could
show. Both outcomes are informative.

## 1. Units and reused machinery

Dimensionless units, Δ = 1; multiply Δ-units by **0.25 meV** for physical
energies (LO_PARAMS Δ = 0.25 meV, so t_hop = 16.57 meV ⇒ **t = t_norm ≈ 66.3**,
α = α_norm = 8.0 = α_SO/(Δ·a_lat)). Everything geometric/circuit-level is reused
unchanged from R4/R5/R7 — the self-field kernel `ring_biot_savart_self`, the
Peierls conversion `compute_peierls_from_A`, the fixed external flux
`_phi_ext_bond`, the production coupling `G_R7`, and the lattice constant `A0`
— because the transverse subband structure does not change the *longitudinal*
self-field geometry. `_build_ring_lo_phased` is imported **only** as the M = 1
bit-for-bit reference. Nothing in R1–R8 or R9 Tracks 1–3 is modified.

## 2. The 4NM Hamiltonian (basis = subband(M) ⊗ spin(2) ⊗ site(N), Nambu-doubled)

`build_ring_lo_multisubband`:

- **Diagonal subband blocks** — R6's `_build_ring_lo_phased` normal block at the
  shifted chemical potential **μ − Δε_n**, one per subband n; intra-subband
  s-wave pairing (block-diagonal). This is the correct 4NM generalization, *not*
  4N·M² — the M×M structure is block structure, cost ∝ (4NM)³ (scoping §4.1,
  measured cubic in M).
- **Off-diagonal blocks** — the *inter-subband* Rashba coupling C[n,m]·(σ_x ⊗
  I_N), site-diagonal (k_y is a transverse, local operator). Lim, Serra, López &
  Aguado (PRB 86, 121103) show this term hybridizes Majorana pairs and closes
  the effective gap, so omitting it would be *qualitatively* wrong, not merely
  approximate. The hard-wall selection rule (scoping §5.4):

  > ⟨n|k_y|m⟩ = −i · 4nm / (W(n²−m²))  for (n−m) odd,  0 otherwise
  > C[n,m] = (α_SO/Δ)·⟨n|k_y|m⟩,  α_SO/Δ = α_norm·a_lat = 80 nm.

  C is purely imaginary and antisymmetric in (n,m) ⇒ Hermitian. The intra-subband
  Rashba piece σ_y·k_x is already inside each diagonal block (R6's h_SO); k_y has
  no diagonal element, so there is no double-counting.
- **Hole block** = the exact −h*, so PH symmetry holds for any M and any
  phi_bonds.

### Transverse subband energetics (measured from the first subband)

ε₁ = t·(π a/W)² and **Δε_n = ε₁·(n²−1)**, so Δε₁ = 0 and M = 1 reproduces R6
exactly. Subband n becomes the newly-occupied channel as μ crosses **its
threshold μ = Δε_n**, where its 1D van Hove singularity ρ(E) ∝ 1/√(E − Δε_n)
sits. At the locked geometry **W = 50 nm**: ε₁ = 6.54 meV (= 26.16 Δ), Δε₂ =
19.6 meV (78.5 Δ), Δε₃ = 52.3 meV (209 Δ), Δε₄ = 98.1 meV (392 Δ). One subband
spacing near M is ΔE_sb(M) = Δε_{M+1} − Δε_M = ε₁(2M+1).

## 3. Van Hove exclusion radius δ — closed form, and verified AT THRESHOLD

**Closed form (locked, from R6's own dispersion, blind to s):**

> δ(ε) = t · (2(1−ε)/π)² · (1 + margin),  ε = 0.15, margin = 0.15
> ⇒ δ = 22.31 Δ = 5.578 meV  at t = 66.3.

Derivation (scoping §7 item 3b): near a band edge the discrete level spacing s_j
grows from the edge toward the average bulk spacing 4t/N; requiring |s_j/s_bulk −
1| < ε for three consecutive levels and solving in the large-N limit gives δ =
t·(2(1−ε)/π)², independent of N. No Hamiltonian is diagonalized to *derive* it —
it is the known dispersion evaluated as arithmetic.

**Verification against the genuinely diagonalized spectrum, SPECIFICALLY AT
THRESHOLD** (`verify_van_hove_delta_at_threshold`, committed data in
`r9_multisubband_preregistration.json`): a genuine `np.linalg.eigvalsh` of the
tight-binding chain matrix −t(T+Tᵀ), with δ_measured read as the energy distance
FROM the band edge (the van Hove point, *not* the smooth bulk) at which the
near-edge spacing first converges to within ε of bulk for three consecutive
levels:

| N | δ_measured (Δ) | δ_closed/δ_measured | eigvalsh vs analytic |
|---:|---:|---:|---:|
| 80 | 24.61 | 0.788 | 1.4e-13 |
| 160 | 23.49 | 0.826 | 1.5e-13 |
| 320 | 22.19 | 0.874 | 1.9e-13 |
| 640 | 21.55 | 0.900 | 2.4e-13 |

Findings, reported honestly:
- The genuine eigvalsh reproduces the analytic dispersion to ~1e-13 (the
  diagonalization is real, not the formula in disguise).
- The closed form runs **below** the measured crossover at every N
  (`conservative_low = True`), so excluding the +15%-margin δ over-excludes the
  *continuum* crossover — the safe direction.
- **Limitation (not silently patched):** the measured crossover is mildly
  N-dependent (24.6 Δ at N = 80 → 21.5 Δ at N = 640, converging toward the
  large-N closed form 19.4 Δ; 13% spread). At the production N = 80 the
  finite-size crossover (24.6 Δ) is ~10% *wider* than the locked δ-with-margin
  (22.3 Δ). The locked δ excludes the continuum van Hove crossover + 15%, which
  is the N-independent physical scale; a few μ-points just outside the excluded
  zone at N = 80 may retain mild band-edge curvature. This is recorded as a
  limitation rather than fixed by inflating the margin post-hoc, and it is
  immaterial to the median statistic (which is dominated by points far from any
  threshold). Any future, separately pre-registered follow-up could lock the
  margin from the finite-N crossover instead.

## 4. Construction gates (all pass — `r9_multisubband_preregistration.json`)

| Gate | Result |
|---|---|
| **M = 1 bit-for-bit** vs `_build_ring_lo_phased` (H and current operator) | max\|ΔH\| = **0.0**, max\|ΔI\| = **0.0** |
| **Decoupled limit** (coupling zeroed): 4NM = ⊕ of M shifted 4N blocks | entrywise **0.0** (spectrum 3.5e-12, LAPACK noise) |
| **PH symmetry** Σ_x H* Σ_x = −H, M = 1..4, coupling on | **0.0** at every M |
| **Inter-subband selection rule** (Hermitian, zero diagonal, zero even-diff) | exact |

The decoupled-limit gate is enforced *entrywise* (reorder the 4NM matrix from
Nambu⊗subband⊗spin⊗site to subband⊗Nambu⊗spin⊗site and compare to
`block_diag` of the per-subband blocks) — exactly 0.0, the analogue of the M = 1
gate. The eigenvalue-union comparison is machine precision (~1e-12), not 0.0,
because it compares independent LAPACK reductions; that is reported as noise, not
a failure.

## 5. Operating points and windows (locked, blind to s)

**Protocol (a) — truncation convergence (Task 1, the cheap sanity check).** Hold
μ FIXED at the R6/R9-T1 single-channel operating point μ = 1.2 (subband 1 filled
1.2 above its bottom), B = √(1.2²+1) − 0.05 = 1.512, and raise only the basis
size M = 1..4. Subbands 2..M stay unoccupied — they enter only virtually, through
second-order inter-subband admixture (matrix elements ∝ 4nm/(W(n²−m²)), energy
denominators ∝ n²). s(M) should converge quickly and smoothly. This validates
the multi-subband construction; it is **not** the §6.1 physics answer (no new
current-carrying channel is added).

**Protocol (b) — genuine added occupancy (Task 2, the real question).** For each
M = 1..4, place the newest channel M near its own threshold Δε_M with channels
1..M−1 deep, and sweep μ over a window that straddles Δε_M (so the van Hove point
is inside the raw window and its exclusion is auditable) and spans one full
subband spacing of genuine occupancy above it:

> **μ ∈ [Δε_M − ΔE_sb(M)/3, Δε_M + ΔE_sb(M)]**, ΔE_sb(M) = ε₁(2M+1),
> field **B(μ) = √((μ − Δε_M)² + 1) − 0.05** (channel M held just below its own
> transition, the multi-subband analogue of R9-T1's operating field),
> coupling g = **G_R7 = 2.9×10⁴** (production).

**Direction clarification (locked here, blind to s).** Scoping §7 item 3b wrote
the window as `[ε_M − ΔE_sb, ε_M − δ]` (below threshold), which contradicts its
own §6.1a ("μ must be re-tuned *upward* … so that channel M is genuinely
*occupied*"). Genuine occupancy requires μ > Δε_M (the band bottom); below it the
channel is empty. The window is therefore placed *above* threshold. This is a
pre-registration correction of an internal inconsistency in the scoping prose,
made before any s(μ;M) value exists — not a post-hoc tune.

## 6. The comparison statistic (locked)

`preregistered_statistic`: the reported number is the **median of s(μ) over the
surviving window**, where every μ within the van Hove radius δ of any subband
threshold Δε_n (n ≤ M) — and of any lobe-crossing point detected from the bare
multi-subband gap within the window — is excluded. Median, not max, because max
would be pulled toward whatever point sits closest to an excluded boundary,
reintroducing exactly the sampling risk the exclusion exists to avoid. The **raw
(non-excluded) median** and the excluded/raw point counts are reported alongside,
so the effect of the exclusion is visible and auditable. Both directions of the
result are pre-declared informative: rising median-s(M) supports "the null is a
single-channel artifact"; flat/bounded median-s(M) supports "the null is robust
and gets stronger."

**Lobe-crossing exclusion.** The scoping report notes this radius is not
closed-form the way δ is (it depends on the coupled spectrum, which does not
exist before the builder is run). Because protocol (b) holds channel M just below
its *own* transition and the deep channels far from theirs, no lobe crossing is
expected inside the window; any that appears (bare-gap minimum near zero) is
excluded using the same δ as a conservative, blind proxy, and reported.

## 7. Scope limits (restated wherever a result is reported)

- **Structural, not device-faithful.** M ≤ 4 is far below the 10–26 occupied
  subbands of real Chang-class wires (scoping §3.3). Adding a subband can flip
  the invariant's parity (Lim et al., PRB 86, 121103), so small-M behaviour must
  not be assumed to extrapolate to the many-subband regime.
- **Loop-gain probe only.** ~30 eigh per point (a root-solve + a
  finite-difference slope), *not* the full B_c bisection (~2000 eigh) or any
  self-consistency loop. The full self-consistent multi-subband EM-feedback
  calculation is out of reach (scoping §4.5) and is not attempted.
- **Weak sum-rule expectation.** At fixed geometry, s should be *monotone
  non-decreasing* in occupied-channel count (more current-carrying channels), but
  with unequal, nonlinear, depth-dependent weights (channels near their band
  bottoms contribute less; inter-subband coupling redistributes current). Any
  observed trend is judged against this weak envelope, not a linear ∝ M law.
