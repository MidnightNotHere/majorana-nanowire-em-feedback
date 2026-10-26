# Methods note: R6 — full Lutchyn-Oreg phase diagram, realistic InAs/Al parameters

**Project:** Inductive Circuit Feedback Shifts the Topological Phase Boundary in InAs/Al Majorana Nanowires
**Module:** `src/bdg_r6_lo_phase_diagram.py`
**Phase 1 reconnaissance:** structural audit found `build_bdg_lo()` already implements a genuine 4N×4N spin⊗Nambu Lutchyn-Oreg Hamiltonian (not a 2N→4N extension); target figure selected: Oreg–Refael–von Oppen (arXiv:1003.1145) Figure 1 (BdG spectrum vs Zeeman field, gap closing at `B=√(μ²+Δ²)`) as the primary deliverable, with a secondary (μ,B) 2D phase diagram using the corrected `pfaffian_trim_lo()` (see `docs/methods/r2_o1_pfaffian_trim_lo_correction.md`).

## 1. Parameter provenance and the idealized-regime disclaimer (locked in before implementation)

Parameters are drawn from InAs bulk material properties (effective mass,
Rashba coupling — Lutchyn/Sau/Das Sarma 2010, Oreg/Refael/von Oppen 2010)
and from epitaxial-Al induced pairing where directly measured (Δ≈0.19 meV,
Chang et al. 2015, arXiv:1411.6255), combined with a Zeeman field range
chosen to access the theoretical topological regime, as in
Oreg–Refael–von Oppen and Lutchyn–Sau–Das Sarma.

**We do not claim this specific parameter combination is simultaneously
achievable in a single real device.** At the measured axial critical field
of a comparable thin epitaxial-Al shell (~75 mT, Chang et al. 2015), the
wire's own superconductivity would collapse well before the Zeeman field
needed to reach this model's topological transition (~0.4–0.5 T at these
parameters) — verified directly against the primary source (field applied
along the nanowire axis, the orbitally-favorable geometry for a thin shell,
so ~75 mT is close to that device's ceiling, not a suppressed worst case).
A search for a verified higher in-plane critical field for a comparable
InAs/Al device did not turn up a number that survived direct
primary-source verification (two candidate numbers from search snippets
were checked against their claimed sources; one was confirmed **absent**
from the cited paper, the other could not be confirmed or denied with the
PDF tooling available). This is an idealized/theoretical parameter regime,
in the same spirit as the original papers' figures (both present their
key figures in dimensionless/rescaled natural units, separately from the
real-material estimates quoted in their text — neither claims their figure
values are simultaneously realized in one measured device), not a
reproduction of an experimentally realized operating point.

**Sharpened point — this is a limitation of the model, not just of the
device.** The code itself treats `Δ` as a fixed constant, independent of
`B_zeeman` (`h_Z = B_zeeman·kron(σ_z, I_N)` is a bare Zeeman term; no
orbital-coupling or field-dependent pair-breaking term of any kind exists
anywhere in `build_bdg_lo()`). So the computed phase boundary is not merely
*drawn at* an experimentally inaccessible field — it is *computed by* a
Hamiltonian that assumes the induced gap survives undiminished at every
field strength plotted, however large. A more complete model would
suppress `Δ(B)` toward zero as `B` approaches the shell's real critical
field, which would likely round off or eliminate the sharp transition this
model shows at large `B_zeeman`. That physics is not implemented here and
is an explicit, known, out-of-scope limitation of this idealized
treatment — not an oversight to be silently corrected later, but a
documented modeling choice to carry into every figure caption that plots
`B_zeeman` beyond the verified ~75 mT scale.

**Figure caption boilerplate (for every R6 figure sweeping B_zeeman into
the idealized regime):**

> Parameters: InAs bulk material properties (Lutchyn/Sau/Das Sarma 2010;
> Oreg/Refael/von Oppen 2010) and epitaxial-Al induced pairing
> (Chang et al. 2015) where directly measured, combined with a Zeeman
> field range chosen to access the theoretical topological regime as in
> the target papers. Idealized/theoretical regime, not a reproduction of
> an achievable single-device operating point: the model holds Δ fixed
> and does not include field-dependent pair-breaking, and the plotted
> field range exceeds the verified axial critical field (~75 mT,
> Chang et al. 2015) of a comparable thin epitaxial-Al shell.

## 2. Implementation status

**Module:** `src/bdg_r6_lo_phase_diagram.py`. Purely additive: imports
`build_bdg_lo`, `pfaffian_trim_lo`, `check_ph_symmetry_complex` read-only
from R1/R2, and `build_ring_lo` read-only from R4. No R1-R5 file was
modified; `git status` after this work shows only new files
(`bdg_r6_lo_phase_diagram.py`, this note, `results/lutchyn_oreg_phase_diagram/`).

### 2.1 Dispersion method — a course correction during development

The original plan was to hand-derive a Bloch Hamiltonian H(k) for the
dispersion figure. A first attempt had a sign error in the BdG hole-block
convention (caught by cross-validating its eigenvalues against a
self-contained periodic ring built the same way `build_bdg_lo` assembles
its blocks — max discrepancy 1.2 across 20 random parameter draws,
clearly wrong). Rather than keep hand-deriving momentum-space physics,
the module instead reuses R4's already-validated `build_ring_lo` directly
and extracts momentum-resolved dispersion by projecting each eigenvector's
particle-sector amplitude onto plane waves at the ring's allowed momenta
(`ring_dispersion()`). This was cross-checked against the independently-
fixed TRIM Pfaffian: the minimum bulk gap is consistently located at k=0
and matches `|B_zeeman - sqrt(mu^2+Delta^2)|` in both the trivial and
topological regimes tested (mu=0, B=0.5/1.5 Bc) — a genuine, non-trivial
agreement between two independently-computed quantities (a diagonalized
spectrum and a Pfaffian sign), not assumed.

Along the way, a suspected inconsistency between R4's `build_ring_lo` and
R1/R2's `build_bdg_lo` hole-block conventions was raised and then
DISPROVEN by direct entrywise comparison (periodicizing `build_bdg_lo`'s
own active hole-block formula and diffing against `build_ring_lo`'s
output at N=8: max difference 0.0, exact match). The apparent
inconsistency was an arithmetic error in manual reasoning (incorrectly
treating the SOC block as purely imaginary; it is real once combined with
the `i*sigma_y` prefactor), not a real bug in either function. Recorded
here so the false alarm isn't rediscovered later.

### 2.2 Known-limit validation results

Both required before any figure was trusted (`bdg_r6_lo_phase_diagram.py`
raises if either fails):

- **Zero-Zeeman Kramers degeneracy** (`validate_zero_zeeman_kramers_degeneracy`):
  at B_zeeman=0, exactly 2N of the 4N ring eigenvalues pair up to
  <1e-8 (full double degeneracy, time-reversal restored); at B=0.4 meV
  that count strictly decreases (100 -> 96 at N=50). Note: a naive
  fixed-index pairing of the sorted spectrum initially reported a false
  FAIL, because this Hamiltonian has an UNRELATED accidental degeneracy
  from ring k<->-k reflection symmetry that persists at B!=0 and gets
  conflated with Kramers pairing under naive pairing — fixed by counting
  near-zero adjacent gaps instead of assuming a fixed pairing pattern.
  **PASS.**
- **Zero-SOC bulk gap closes away from TRIM** (`validate_zero_soc_gap_closes`):
  checked as a GENERAL property at SIX (mu, B) points, all topological by
  the TRIM-only criterion -- (0,1.5), (0,1.2), (0.3,1.3), (-0.5,1.5),
  (0.7,1.5), (0,2.0) (units of Delta) -- not a single representative
  point. At alpha=0 the TRIM Pfaffian (SOC-independent by construction)
  still reports "topological", yet the true continuum bulk gap is ~0 at
  every one of these points: a continuum-k scan (20000 points) gives
  1.0e-3, 2.3e-4, 2.8e-4, 6.9e-5, 8.3e-4, 1.0e-3 meV respectively (the
  small residuals are O(dk) near the linear band crossing, i.e.
  grid-limited, not a real gap). This confirms SOC is genuinely essential
  for a bulk-gapped, protected topological phase in a 1D BdG chain -- the
  established physics result. **PASS.**

  **Method correction made during this check (recorded so it isn't
  redone):** an earlier version used a ring N-scaling at the single point
  (0,1.5) and inferred the general claim from it. When broadened to
  multiple points, (mu=0, B=2.0) reported a spurious "not collapsing"
  under ring N-scaling -- its ring gap plateaus at 0.086 through N=640 --
  purely because the discrete ring momenta 2*pi*j/N never landed near that
  point's band-crossing momentum k* (which shifts with B). The continuum-k
  scan shows the true gap there is 1.0e-3 (gapless, like the others). The
  validation was therefore switched from ring N-scaling to a continuous-k
  scan of an FT-verified Bloch Hamiltonian (`_bloch_H_lo`), robust to this
  sampling issue; a ring N-scaling at (0,1.5) (0.414 -> 0.101 -> 0.030 at
  N=80/320/640) is retained only as a secondary illustration. This is
  exactly the "checked at one point, assumed general" gap the reviewer
  flagged -- it was found to matter (one of six points behaved differently
  under the naive method) and was fixed, not reworded.

  `_bloch_H_lo` provenance: derived by directly Fourier-transforming
  build_ring_lo's blocks and validated to reproduce its full spectrum to
  <1e-12 across 40 random parameter draws at arbitrary alpha. An initial
  hand-derived version had TWO sign errors (SOC off-diagonal sign, and
  -h(-k) instead of the correct -h(k) = -h*(-k) hole block) that partially
  cancelled and left a ~0.36 spectrum discrepancy; the FT-verified version
  is exact -- same "anchor momentum-space physics to the validated
  real-space construction" lesson as in 2.1.

### 2.3 R6 realistic InAs/Al parameters used

`R6_PARAMS` (see module Section 1): m*=0.023 me [VERIFIED, Ioffe NSM
archive], g=10.5 [VERIFIED range, NOT the old Mourik/InSb citation],
alpha_SO=20.0 meV*nm [VERIFIED order-of-magnitude vs LSDS], **Delta=0.19
meV [VERIFIED, CORRECTED from the old 0.25 meV — Chang et al. 2015
full-shell value, replacing the InSb-sourced number]**. Derived: t_hop=
16.565 meV, alpha_phys=2.000 meV. At mu=0, Bc=0.190 meV=0.625 T vs the
verified measured axial critical field of 0.075 T (Chang et al. 2015) —
**idealized-regime ratio 8.3x**, exactly the tension flagged in Section 1
and stamped on every figure via `IDEALIZED_REGIME_CAPTION`.

### 2.4 Figures produced (`results/lutchyn_oreg_phase_diagram/`)

- **`r6_figure1_orv_dispersion.png`** — three-panel BdG dispersion at
  mu=0 (trivial/critical/topological), qualitatively reproducing ORV
  Fig. 1's gap-closes-and-reopens story with this codebase's own
  Hamiltonian and InAs/Al parameters (not a copy of the published
  figure). At these realistic parameters t (~16.5 meV) is far larger
  than Delta/B (~0.1-0.4 meV), so the relevant physics occupies a narrow
  window of the lattice Brillouin zone near k=0; the figure's k-axis is
  zoomed to that window (N=300, re-centred to (-pi,pi]) rather than
  showing the full, mostly-empty-at-this-scale [0,2pi) range — an
  intermediate version with the full BZ shown was visually nearly blank
  and was corrected before finalizing.
- **`r6_figure2_gap_vs_Bz.png`** — ring bulk gap and the genuine TRIM
  Pfaffian invariant vs B_zeeman at mu=0 and mu=0.3 meV. The gap closes
  exactly where the invariant flips, at both mu values.
- **`r6_figure3_phase_diagram.png`** — 2D (mu, B_zeeman) phase diagram at
  realistic InAs/Al parameters, genuine invariant (left) and ring bulk
  gap (right, log scale), analytic boundary overlaid.

### 2.5 Comparison to the analytic boundary (honest, not assumed)

| mu (meV) | analytic Bc (meV) | spectral gap-min (meV) | invariant flip (meV) | discrepancy | B-sweep grid step | discrepancy / step |
|---|---|---|---|---|---|---|
| 0.000 | 0.1900 | 0.1868 | 0.1868 | -0.0032 meV | 0.0064 meV = 0.0212 T | -0.50 |
| 0.300 | 0.3551 | 0.3491 | 0.3491 | -0.0060 meV | 0.0120 meV = 0.0396 T | -0.50 |

The spectral-gap-minimum location and the genuine-invariant flip location
agree with EACH OTHER exactly at both mu values (two independently
computed quantities landing on the identical B_zeeman). The discrepancy
from the analytic boundary is **exactly -0.50 grid steps at both points**:
the analytic Bc falls precisely midway between two sampled grid points, so
the nearest sampled minimum necessarily sits half a step below it. This is
a textbook signature of finite sweep resolution locating a continuous
minimum, NOT a physical deviation from B=sqrt(mu^2+Delta^2).

The grid step is given in physical units so the precision is judgeable
directly rather than hidden behind "one grid step": at mu=0 the B-sweep
resolution is **0.0064 meV = 0.0212 T** (the -0.0032 meV discrepancy =
0.0106 T, half a step); at mu=0.3 it is **0.0120 meV = 0.0396 T** (the
-0.0060 meV discrepancy = 0.0198 T, half a step). Halving the sweep step
would halve both discrepancies -- the tell that this is resolution, not
physics. No parameter was adjusted to produce this agreement; the sweep
was run once at its stated resolution. (Tesla conversion via
B_zeeman = g*mu_B*B/2, g=10.5; helper `zeeman_meV_to_tesla()`.)

### 2.6 What this does and does not establish

This reproduces ORV's qualitative spectral picture and the standard
(mu, B) phase boundary using the codebase's own genuine 4N Hamiltonian and
now-genuine invariant, at literature-sourced InAs/Al parameters. It does
**not** reproduce a specific numbered published figure pixel-for-pixel
(neither target paper contains one in this form — see Phase 1
reconnaissance), and it explicitly does **not** claim the plotted
parameter combination is achievable in one real device (Section 1). The
open item flagged in Section 1 (no Delta(B) suppression physics in the
model) remains unaddressed — future work, not silently assumed away.
