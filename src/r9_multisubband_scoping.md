# R9 — Track 4 (SCOPING ONLY): multi-subband extension of the Lutchyn–Oreg model

**Status:** analysis and planning. **No simulation code was written or run for
this document.** Nothing in `src/` was modified; no new Hamiltonian exists.
The only code executed was (a) closed-form arithmetic for transverse subband
energetics and (b) a wall-clock benchmark of `scipy.linalg.eigh` on *random*
Hermitian matrices, used to calibrate the cost model in §4. Both ran in a
scratch directory outside the repo and are hardware/arithmetic measurements,
not physics runs.

**Question this document answers:** is extending R6's single-subband 4N×4N
Lutchyn–Oreg Hamiltonian to a multi-subband model a tractable, well-motivated
next step — and if so, at what scope?

**Bottom line up front.**
1. The multi-subband literature is real, specific, and directly relevant
   (§2). It is not a speculative extension.
2. **The single-subband assumption in R6/R8 is almost certainly *not*
   realistic for the Chang et al. 2015 InAs/Al device class** (§3). This is a
   caveat that should be recorded against R6/R8 *regardless of whether Track 4
   is ever implemented* — it costs zero compute and it is the highest-value
   finding in this document.
3. Cost scales as **M³ per solve** (M = number of retained subbands), measured,
   not assumed (§4). The bottleneck is **"the same number of solves, each M³
   times larger"** — the one case where `bdg_parallel` provides **no additional
   relief at all** (§4.4).
4. Consequently: extending the **self-consistent EM-feedback stack** (R4/R5/R7 —
   the project's actual headline) to multi-subband is **not tractable** on this
   machine (§4.5, §7). Extending the **phase diagram only** (R6) is tractable to
   M ≈ 4–6.
5. Recommendation (§8): **conditional yes, at a sharply bounded scope.** The
   minimum viable scope is a phase-diagram structural study plus one cheap,
   targeted probe of the loop gain vs M — the single measurement that could
   overturn the project's null lineage. Full multi-subband self-consistency is
   explicitly out of reach and should be named as such, not deferred vaguely.

---

## 1. What R6 actually assumes

`build_bdg_lo` (R1/R2, used by R6/R8) builds a **4N×4N** Hamiltonian on the
basis site(N) ⊗ spin(2) ⊗ Nambu(2). There is exactly **one transverse
channel**. The transverse direction does not appear in the Hamiltonian at all:
it is integrated out by assumption, and `mu` is implicitly measured from the
bottom of that single subband.

R6's verified parameters (`R6_PARAMS`, `src/r6_lo_phase_diagram_methods.md` §2.3):
m\* = 0.023 mₑ, g = 10.5, α_SO = 20 meV·nm, **Δ = 0.19 meV** (Chang et al.
2015), a_lat = 10 nm, t_hop = 16.565 meV. R6 sweeps μ ∈ [−0.6, +0.6] meV and
B_zeeman ∈ [0, 0.4] meV — i.e. **μ is confined to within ±3Δ of the assumed
first-subband bottom.**

That is the assumption under test here.

---

## 2. Task 1.1 — literature grounding (verified against primary sources)

Every claim below was checked against the arXiv/open-access record in this
session. Where I could only verify the abstract, I say so and do not attribute
finer detail to the source.

### 2.1 The core result: subband *parity* controls topology

**Potter & Lee, "Multichannel Generalization of Kitaev's Majorana End States
and a Practical Route to Realize Them in Thin Films", Phys. Rev. Lett. 105,
227003 (2010), [arXiv:1007.4569](https://arxiv.org/abs/1007.4569).**
Verified abstract, verbatim: Majorana end states *"survive beyond the strict 1D
single-channel limit so long as the sample width does not exceed the
superconducting coherence length, and exist when an odd number of transverse
quantization channels are occupied. Consequently we find that the system
undergoes a sequence of topological phase transitions driven by changing the
chemical potential."*

This is the key structural claim: the phase diagram acquires **a sequence of
transitions**, not one boundary. Two conditions, both relevant to us — an
**odd** occupied-channel count, and a width bound set by ξ.

**Potter & Lee, "Majorana End-States in Multi-band Microstructures with Rashba
Spin-Orbit Coupling", Phys. Rev. B 83, 094525 (2011),
[arXiv:1011.6371](https://arxiv.org/abs/1011.6371).** Extends the above from an
idealized spinless p+ip superconductor to the *Rashba* case — i.e. to the
actual Lutchyn–Oreg ingredient set R6 uses. Verified abstract confirms the
odd-subband condition and the width-vs-ξ condition carry over, and that the
Rashba case has *"a number of new features"*; verified from the fetched text:
*"finite size quantization allows Majorana end-states even when the chemical
potential is outside of the induced Zeeman gap"*, and *"the bulk energy gap is
substantially reduced if the induced magnetization is too large."*

### 2.2 The multi-subband Lutchyn–Oreg phase diagram specifically

**Lutchyn, Stanescu & Das Sarma, "Search for Majorana fermions in multiband
semiconducting nanowires", Phys. Rev. Lett. 106, 127001 (2011),
[arXiv:1008.0629](https://arxiv.org/abs/1008.0629).** Verified abstract,
verbatim: *"We show that when odd number of subbands are occupied the system
realizes non-trivial topological state supporting Majorana modes localized at
the ends. We study the topological quantum phase transition in this system and
analytically calculate the phase diagram as a function of the chemical
potential and magnetic field. Our key finding is that multiband occupancy not
only lifts the stringent constraint of one-dimensionality but also allows to
have higher carrier density... in the limit of strong interband mixing there is
an optimal regime in the phase diagram ("sweet spot") where the topological
state is to a large extent insensitive to the presence of disorder."*

**Stanescu, Lutchyn & Das Sarma, "Majorana Fermions in Semiconductor
Nanowires", Phys. Rev. B 84, 144522 (2011),
[arXiv:1106.3078](https://arxiv.org/abs/1106.3078).** The long companion paper,
and the closest thing to a direct multi-subband generalization of what R6
computes. Verified abstract, verbatim: *"The non-trivial topological state
corresponds to a superconducting phase supporting an odd number of pairs of
Majorana modes localized at the ends of the wire, whereas the non-topological
state corresponds to a superconducting phase with no Majoranas or with an even
number of pairs."*

From the full text (fetched via the ar5iv HTML rendering — flagged as a
secondary rendering of the primary source, not the published PDF): the
single-band criterion is `|V_x| > √(Δ² + μ²)`; the occupancy rule is stated as
*"an odd number corresponds to a topologically non-trivial SC, while an even
number results in a standard superconductor"*; the phase diagram is organised
around **"sweet spots" — described in-text as multicritical points** where
subbands become degenerate at zero energy. Modelled geometry: **L_z = 5 nm,
L_y = 130 nm**, m\* = 0.04 m₀, α_R = 0.1 eV·Å, with an **inter-subband spacing
ΔE_sb ≈ 1.6 meV** and Δ₀ = 1 meV. It also finds that inhomogeneous
semiconductor–superconductor coupling generates *inter-subband pairing* terms
reaching *"about 50% of the diagonal contributions for neighboring sub-bands"*,
and that near the sweet spots the gap *vanishes* without inter-subband pairing.

> Caveat on this source: the ar5iv rendering garbles some equation numbering.
> I have quoted only statements that are self-consistent and consistent with
> the verified abstract. The specific claim "inter-subband pairing reaches ~50%
> of diagonal terms" should be re-verified against the published PDF before it
> is used in any figure caption or write-up.

### 2.3 The counter-result — inter-subband coupling is not benign

**Lim, Serra, López & Aguado, "Magnetic-field instability of Majorana modes in
multiband quantum wires", Phys. Rev. B 86, 121103(R) (2012),
[arXiv:1202.5057](https://arxiv.org/abs/1202.5057).** Verified abstract,
verbatim: *"Our results demonstrate that interband mixing coming from Rashba
spin orbit term hybridizes Majorana pairs originating from different transverse
modes while simultaneously closing the effective gap. Consequently, multiple
Majorana modes do not coexist in general. On the contrary, Majorana physics is
robust provided that only one single transverse mode contributes with a
Majorana pair."*

**This is the most important source for scoping.** It says the naive picture —
"M subbands ⇒ M independent Majorana pairs, count the parity" — is wrong once
the *inter-subband Rashba* term is included: the coupling hybridizes the pairs
*and closes the gap*. So a multi-subband implementation that omits the
inter-subband Rashba term would not merely be less accurate, it would be
qualitatively wrong in exactly the regime it was built to study. §5.4 makes
this a mandatory validation target.

### 2.4 Summary of what the literature establishes

| Claim | Source | Status |
|---|---|---|
| Majoranas survive multi-channel wires; **odd** channel count required | Potter–Lee PRL 105, 227003 | verified (abstract) |
| Phase diagram becomes a **sequence** of transitions vs μ | Potter–Lee PRL 105, 227003 | verified (abstract) |
| Odd/even subband rule holds for the Rashba (Lutchyn–Oreg) case | Potter–Lee PRB 83, 094525 | verified (abstract) |
| Multiband analytic (μ, B) phase diagram; "sweet spots" = multicritical points | Lutchyn–Stanescu–Das Sarma PRL 106, 127001; PRB 84, 144522 | verified (abstract; detail via ar5iv) |
| **Inter-subband Rashba hybridizes Majorana pairs and closes the gap** | Lim–Serra–López–Aguado PRB 86, 121103 | verified (abstract) |
| Width must not greatly exceed ξ | Potter–Lee (both) | verified (abstract) |

---

## 3. Task 1.2 — is R6/R8's single-subband assumption realistic for these devices?

**Short answer: no, and this should be flagged against R6/R8 now.**

### 3.1 What the primary source actually states

**Chang, Albrecht, Jespersen, Kuemmeth, Krogstrup, Nygård & Marcus, "Hard gap
in epitaxial semiconductor–superconductor nanowires", Nature Nanotechnology 10,
232 (2015), [arXiv:1411.6255](https://arxiv.org/abs/1411.6255)** — the source
R6 already uses for Δ = 0.19 meV, T_base = 20 mK and B_c ≈ 75 mT.

I re-read it for geometry. Verified statements: induced gap *"Δ\* = 190 μeV"*;
*"base temperature of 20 mK"*; *"B_c (∼75 mT for epitaxial, ∼250 mT for
control)"*; nanowires *"grown in the wurtzite [0001] direction"* with a
*"hexagonal InAs core"*, reaching *"a length of 5–10 μm"*; and for the
half-shell InAs core: **"carrier density (at high V_BG), n = 5×10¹⁸ cm⁻³;
mobility, μ = 3300 cm²/Vs; elastic scattering length, l_e = 100 nm."**

**The paper does not state a nanowire diameter.** I searched the full text,
figure captions and methods; the only nanometre figures are the Ti/Al contact
metallization (5/130 nm) and the 100 nm substrate oxide, neither of which is
the wire. So the subband count cannot be read off directly and must be
estimated. I do that two independent ways below, and state the uncertainty
rather than inventing a diameter.

### 3.2 Estimate A — from the measured carrier density

Taking the *measured* n = 5×10¹⁸ cm⁻³ and a parabolic band at R6's own
m\* = 0.023 mₑ: k_F = (3π²n)^(1/3) = 0.529 nm⁻¹, giving **E_F ≈ 464 meV**
above the band bottom. For a hard-wall cylindrical wire (Bessel-zero levels
E_{m,n} = ħ²j²_{m,n}/2m\*R²) that puts roughly **78 orbital subbands below E_F
at R = 35 nm, ~160 at R = 50 nm.**

Caveats, stated plainly: this is at *high back-gate voltage* (an accumulation
condition, the top of the gate range, not a Majorana operating point); InAs
conduction-band **non-parabolicity is severe at several hundred meV**, so the
true E_F is substantially lower than 464 meV and the true count substantially
below 78; and surface accumulation makes the radial density non-uniform. The
number is therefore an **upper bound on that device's occupancy at that gate
setting**, not the occupancy at a tuned operating point. Even discounted by an
order of magnitude it is not 1.

For reference, the same calculation at lower densities: n = 10¹⁸ cm⁻³ → ~23
subbands (R = 35 nm); n = 10¹⁷ → ~3; n = 10¹⁶ → 0–1. So single-subband
operation requires depleting roughly **two to three orders of magnitude** below
the density this device was characterised at.

### 3.3 Estimate B — from the dedicated literature (the stronger argument)

**Woods, Das Sarma & Stanescu, "Subband occupation in
semiconductor–superconductor nanowires", Phys. Rev. B 101, 045405 (2020),
[arXiv:1910.04362](https://arxiv.org/abs/1910.04362).** This paper exists to
answer exactly the question in Task 1.2, self-consistently (Schrödinger–Poisson),
for **InAs/Al**. Verified abstract, verbatim: *"Realistic surface charge
densities, which are responsible for band bending, are shown to significantly
increase the number of occupied subbands, making it difficult or impossible to
reach a regime where only a few subbands are occupied."*

From the full text (ar5iv rendering; secondary rendering caveat as above):
modelled InAs radii **R = 35 nm and R = 50 nm**, m\* = 0.026 m₀; at zero surface
charge ~26 subbands occupied at one representative gate setting; at high surface
charge *"the fewest number of occupied subbands is ten"*; average inter-subband
separation *"(ΔE)_ave ≈ 15−30 meV"* in the desirable thin-wire/low-charge
regime versus *"(ΔE)_ave ≈ 2−7 meV"* in the thick-wire/high-charge regime; and
the conclusion that *"gating the system to reach the few-subbands regime is
usually not possible (except for rather narrow parameter windows) due to the
onset of holes."*

They further note that **"high subband occupancy by itself works as an
effective emergent disorder"** — which connects directly to R9 Track 3, the one
place this project's null lineage already broke.

### 3.4 Cross-check against R6's own parameters

Using R6's m\* = 0.023 mₑ, the first hard-wall transverse level and the
(2,1)–(1,1) spacing of a square cross-section are:

| W (nm) | E₁ (1D hard wall) | square-section (2,1)−(1,1) spacing |
|---|---|---|
| 50 | 6.54 meV | 19.6 meV |
| 70 | 3.34 meV | 10.0 meV |
| 75 | 2.91 meV | 8.72 meV |
| 100 | 1.64 meV | 4.91 meV |

For a *cylindrical* wire at R = 35 nm (Woods 2020's geometry), E₁ = 7.82 meV
and E₂ − E₁ = 12.0 meV.

So the subband spacing is **~5–20 meV, i.e. ~25–100 × Δ**, while R6 sweeps μ
over ±0.6 meV ≈ ±3Δ. **Read one way this vindicates R6:** if μ genuinely sits
at the bottom of the first subband, the second is tens of Δ away and cannot
mix. **Read the honest way it does not:** nothing pins μ there. Chang's measured
density and Woods' Schrödinger–Poisson result both say that in a real InAs/Al
wire μ sits *hundreds* of meV — tens of subbands — above that point, and that
gating down to the first subband is *"difficult or impossible"*.

### 3.5 The caveat that should be recorded against R6/R8

> R6/R8 model a single transverse channel with μ measured from its band
> bottom. The transverse subband spacing at R6's own parameters is ~5–20 meV
> (25–100 Δ), so *conditional on the wire being tuned to the lowest subband*
> the single-channel treatment is internally consistent. However, the primary
> literature for this exact material system (Woods, Das Sarma & Stanescu, PRB
> 101, 045405 (2020)) finds that realistic surface charge makes reaching a
> few-subband regime *"difficult or impossible"*, with ten or more subbands
> occupied in realistic conditions; and Chang et al. 2015 itself reports a
> measured carrier density (5×10¹⁸ cm⁻³ at high V_BG) corresponding to many
> tens of occupied subbands. The single-subband assumption is therefore an
> **idealisation of the device, not a description of it** — of the same
> character as, and to be carried alongside, the fixed-Δ / no-Δ(B)
> idealisation already documented in `r6_lo_phase_diagram_methods.md` §1.

This is a documentation change to R6/R8's caveat set, not a result, and it
holds whether or not Track 4 proceeds. It is the cheapest deliverable here.

---

## 4. Task 2.1 — computational cost: derivation, measurement, and honest verdict

### 4.1 The correct generalization and its dimension

The task brief floated "4N × n_subbands × n_subbands". **That is not the right
scaling, and the distinction matters by a factor of M³.** Deriving it properly:

The single-subband basis is site(N) ⊗ spin(2) ⊗ Nambu(2), dimension 4N.
Adding a transverse subband index n = 1…M inserts **one more tensor factor**:

> basis = site(N) ⊗ **subband(M)** ⊗ spin(2) ⊗ Nambu(2) ⟹ **dim = 4NM**

The Hamiltonian is therefore **(4NM) × (4NM)**. It *contains* an M×M block
structure in subband space (each block 4N×4N: diagonal blocks are the R6
Hamiltonian at shifted chemical potential μ − ε_n, off-diagonal blocks are the
inter-subband coupling), and that block picture is presumably what "n_sub ×
n_sub" was reaching for — but the **linear dimension grows as M, not M²**.
Consequences:

- **element count / memory ∝ (4NM)² = 16N²M²** — quadratic in M;
- **dense diagonalization cost ∝ (4NM)³ = 64 N³M³** — **cubic in M**.

An alternative that avoids mode truncation entirely — discretizing the
transverse direction on the *same* real-space lattice — is treated in §4.6 and
is far worse.

### 4.2 Measured cost model (this machine, not a textbook constant)

`scipy.linalg.eigh` (eigenvalues **and** eigenvectors, as the repo's solvers
use), complex128 random Hermitian matrices, BLAS pinned to 1 thread — the
production configuration `bdg_parallel` uses:

| dim D | 1-thread eigh | 6-thread eigh |
|---:|---:|---:|
| 320 (N=80, M=1) | 0.072 s | 0.036 s |
| 640 | 0.591 s | 0.156 s |
| 960 | 2.013 s | 0.382 s |
| 1280 | 4.780 s | 2.573 s |
| 1600 | 9.390 s | 4.788 s |
| 1920 | 14.621 s | 7.911 s |
| 2400 | 29.859 s | 14.688 s |
| 3200 | 70.431 s | 35.103 s |

**Fitted exponent at 1 thread: 2.98** — a clean cubic, exactly as derived. The
measured ratios relative to D = 320 are 1× / 8× / 28× / 66× / 130× / 203× /
414× / 977× against the predicted M³ = 1 / 8 / 27 / 64 / 125 / 216 / 422 /
1000. The M³ law is **measured, not assumed.**

Fitted constant: **c = 2.15×10⁻⁹ s per dim³**, i.e.

> **T₁(N, M) ≈ 1.38×10⁻⁷ · N³ · M³ seconds**  (1-thread eigh)

which reproduces the measured points to ≲7%.

*(An earlier benchmark at the machine's default thread count returned a fitted
exponent of 2.25 — sub-cubic, and misleading. That was BLAS multithreading
masking the true scaling at large D while being actively harmful at small D
(0.150 s at D = 320 versus 0.036 s at 6 threads, consistent with the
oversubscription-onto-E-cores effect the repo already documents for R7–R9).
The 1-thread column is the one that describes production, and it is the one
used throughout below.)*

### 4.3 Per-solve cost and memory

**Wall clock for one eigh, 1 thread (seconds):**

| N \ M | 1 | 2 | 3 | 4 | 5 | 6 | 8 | 10 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| 60 | 0.03 | 0.24 | 0.80 | 1.90 | 3.71 | 6.42 | 15.2 | 29.7 |
| 80 | 0.07 | 0.56 | 1.90 | 4.51 | 8.80 | 15.2 | 36.1 | 70.4 |
| 120 | 0.24 | 1.90 | 6.42 | 15.2 | 29.7 | 51.3 | 122 | 238 |
| 160 | 0.56 | 4.51 | 15.2 | 36.1 | 70.4 | 122 | 288 | 563 |

**Memory is *not* the binding constraint per solve.** H alone plus ~4× LAPACK
workspace: 0.16 GB at (N=80, M=10), 1.47 GB at (N=120, M=10), 2.62 GB at
(N=160, M=10). One solve fits.

Memory binds the **worker count**, not the solve. `bdg_parallel.default_workers()`
hardcodes a **0.5 GB-per-worker** budget, calibrated for 4N = 320–960. At
(N=120, M=10) each worker needs ~1.5 GB, so that heuristic would silently
over-subscribe — and the repo already notes this laptop *"often has <2–3 GB
free"* (measured during this session: 2.5 GB free of 16.9 GB total). **Any
multi-subband work must re-derive that per-worker budget**; leaving it at 0.5 GB
would produce swapping misread as thermal throttling.

### 4.4 The parallelism question, answered directly

The brief asks whether the bottleneck is "many small solves" (parallelism helps)
or "fewer but much larger solves" (it does not). **Neither, precisely — and the
precise answer is the pessimistic one:**

> Multi-subband keeps **the same number of solves** (the same (B, μ) grid, the
> same bisection and root-finding iteration counts) and makes **each one M³
> times more expensive.**

So the workload is "**the same many solves, each much larger**". The
consequence for `bdg_parallel`:

- Parallelising across grid points still works exactly as well as it does now,
  and yields the same **~2.3–4×** the repo has already measured. **That speedup
  is already banked in the current single-subband numbers.** It is a constant
  factor and it does not grow with M.
- It provides **zero additional relief** against the M³ growth. The M³ lives
  entirely *inside* one `eigh` call.

Could BLAS multithreading rescue a single large solve instead? Measured, no:

| dim | 1-thr | 6-thr | speedup | parallel eff. | throughput 6×1thr | throughput 1×6thr |
|---:|---:|---:|---:|---:|---:|---:|
| 1280 | 4.780 s | 2.573 s | 1.86× | 31% | 1.255 solves/s | 0.389 |
| 1600 | 9.390 s | 4.788 s | 1.96× | 33% | 0.639 | 0.209 |
| 1920 | 14.621 s | 7.911 s | 1.85× | 31% | 0.410 | 0.126 |
| 2400 | 29.859 s | 14.688 s | 2.03× | 34% | 0.201 | 0.068 |
| 3200 | 70.431 s | 35.103 s | 2.01× | 33% | 0.085 | 0.029 |

Six cores buy a **flat ~2×** on a single large `eigh` — **~33% parallel
efficiency**, and it does not improve with size (the exponent is unchanged:
3.53 fitted at 6 threads vs 2.94 at 1 thread over D ≥ 960, the difference being
that the 6-thread curve is a near-constant factor below, not a shallower
slope). Meanwhile **6 workers × 1 thread delivers ~3× the throughput of
1 worker × 6 threads at every dim ≥ 1280.**

> **Verdict:** the repo's existing strategy (many single-threaded workers across
> grid points) remains optimal at every multi-subband size considered, and
> is already saturated. There is no scheduling change, no thread-count
> retune, and no `bdg_parallel` improvement that recovers any part of the M³.
> **Per-solve cost dominates, and it is not parallelizable on this hardware.**

*(Two dims, 640 and 960, show anomalously good 6-thread scaling — 3.8× and
5.3×. That band is almost certainly a cache-residency artifact and I do not
read anything into it; the ≥1280 asymptote is the regime that matters, and it
is flat at ~2×.)*

### 4.5 What this means against the constraints R7 already hit

R7 capped its grid at **N = 120** because per-solve cost is dense-`eigh`-bound
and *"ran ~5× slower under load on this machine"* (thermal throttling,
`R7_SUMMARY.md`). R9 Track 3's full run took **77.1 minutes** measured, against
a 20.6-minute pre-run estimate. Those are the real constraints. Projecting the
two engines the project actually runs:

**(A) R6-style phase diagram** — 40×40 = 1600 grid points at N = 60, roughly
two dense operations per point (one ring `eigh`, one TRIM Pfaffian; both O(D³)):

| M | serial | @4 workers | vs M=1 |
|---:|---:|---:|---:|
| 1 | 0.03 h | 0.01 h | 1× |
| 2 | 0.21 h | 0.05 h | 8× |
| 3 | 0.71 h | 0.18 h | 27× |
| 4 | 1.69 h | 0.42 h | 64× |
| 5 | 3.30 h | 0.83 h | 125× |
| 6 | 5.70 h | 1.43 h | 216× |
| 8 | 13.5 h | 3.38 h | 512× |
| 10 | 26.4 h | 6.60 h | 1000× |

**Tractable to M ≈ 4–6** (well under a day at 4 workers), painful at M = 8,
out of reach at M = 10 once R7's ~5× throttling factor is applied.

Note a scheduling flip worth acting on: the repo currently dispatches the cheap
classification engine **serially**, because at ~0.05 s/task the process-pool
overhead exceeds the task cost. At M ≥ 3 the per-task cost rises to ~1.6 s, so
**that decision reverses** — the cheap engine should be parallelised in the
multi-subband case. It should be re-benchmarked, not assumed either way.

**(B) The self-consistent EM-feedback engine** (R4/R5/R7/R9-T3) — measured in
this repo at **~144–161 s/task at N = 80, M = 1** (`r9_disorder_methods.md`,
"BLAS-oversubscription benchmark"). That is ~2000 `eigh` calls per B_c point
(bisection × root search). Scaling the *measured* 150 s/task by M³:

| M | per B_c point | R9 Track 3's full 77.1-min run, rescaled |
|---:|---:|---:|
| 1 | 150 s | 1.3 h |
| 2 | 1200 s (20 min) | 10.3 h |
| 3 | 4050 s (1.1 h) | **34.7 h** |
| 4 | 9600 s (2.7 h) | **82.2 h** |
| 5 | 18750 s (5.2 h) | **160.6 h (6.7 days)** |
| 10 | 150000 s (41.7 h) | ~55 days |

**and these are pre-throttling numbers.** R7's measured ~5× slowdown under
sustained load applies to exactly this kind of run. Note also that R9 Track 3
found convergence *failures* cost **more** than successes — and inter-subband
coupling closes the gap (§2.3), which is precisely the condition that makes the
self-consistent root-finder struggle. The realistic multiplier is therefore
worse than M³, not better.

> **This is the honest conclusion the brief asked for: extending the
> self-consistent EM-feedback stack to multi-subband is NOT tractable on this
> machine.** Even M = 3 — the minimum that shows any multi-subband structure at
> all — turns a 77-minute track into a 35-hour one before throttling, and the
> project's compute constraints are already binding at M = 1. This should not
> be scoped as "expensive but doable"; it is out of reach, and §8 scopes around
> it rather than through it.

### 4.6 The alternative that looks natural and is not

Instead of truncating to M modes, one could discretize the transverse
direction on the same a = 10 nm lattice R6 already uses. This is appealing
(no truncation, no basis choice, reuses the existing real-space machinery) and
it is **hopeless**:

| wire width | planar (1 transverse dir) | 3D wire (2 transverse dirs) |
|---|---|---|
| 50 nm | 4 channels, D = 1280, 4.8 s/eigh | 16 channels, D = 5120, 0.08 h/eigh, 0.42 GB |
| 70 nm | 6 channels, D = 1920, 14.6 s/eigh | 36 channels, D = 11520, **0.91 h/eigh**, 2.12 GB |
| 100 nm | 9 channels, D = 2880, 51 s/eigh | 81 channels, D = 25920, **10.4 h/eigh**, **10.8 GB** |

At a realistic 70 nm width, **one single diagonalization** of the 3D wire takes
~55 minutes — comparable to R9 Track 3's *entire* run — and a phase diagram
needs 1600 of them. The a = 10 nm lattice does not let you *choose* a small
channel count; the geometry chooses ~36–81 for you.

**Therefore: mode-space truncation (Option A) is the only viable route**, and
its truncation parameter M becomes a validated quantity in its own right
(§5.7), not a free knob.

---

## 5. Task 2.2 — what a real implementation would have to validate

Following this project's established pattern: R7 validated bit-for-bit against
R5 at a = 10 nm (|ΔB_c| = 0); R9 Track 2 gated on θ → 0 reproducing T = 0
bit-for-bit; R9 Track 3 gated on W = 0 reproducing every clean result
bit-identically *before any W > 0 number was trusted*. Track 4's analogue:

### 5.1 The M = 1 gate — mandatory, bit-for-bit, before any M > 1 result

The multi-subband builder at M = 1 must reproduce `build_bdg_lo` /
`build_ring_lo` **entrywise, max |ΔH| = 0.0** — not "to 1e-12". Then every
downstream detector must reproduce its R6 counterpart exactly: the TRIM
Pfaffian sign across a B-grid straddling B_c, the ring bulk gap, R6's 40×40
phase diagram, and (only if the EM stack is touched) R5's self-consistent B_c.
A `run_full()` that raises before computing anything at M > 1 if this fails,
exactly as `validate_W0()` does in Track 3.

### 5.2 PH symmetry at 4NM

`check_ph_symmetry_complex` requires Σ_x H\* Σ_x = −H. Σ_x must be rebuilt for
the new tensor ordering. **Getting the `kron` ordering wrong is by far the most
likely bug in this extension** — R6 already lost time to exactly this class of
error twice (`r6_lo_phase_diagram_methods.md` §2.1: a sign error in the hole
block, and two sign errors in a hand-derived Bloch Hamiltonian that partially
cancelled). PH symmetry is the cheap check that catches it; it must be machine
zero at every M, not small.

### 5.3 Decoupled-subband limit — a genuinely new check

Zero the inter-subband coupling. The 4NM spectrum must then be **exactly** the
sorted union of M single-subband 4N spectra evaluated at shifted chemical
potentials μ − ε_n, to max discrepancy 0.0. This has no single-subband analogue
and it simultaneously validates (a) the block structure, (b) the placement of
the subband bottoms ε_n, and (c) the basis ordering. It is the single highest-
value new test.

### 5.4 Inter-subband coupling — correctness of the term itself

This is the check §2.3 makes mandatory. For a wire along x with hard-wall
confinement along y, the Rashba term α(σ_x k_y − σ_y k_x) splits into an
intra-subband piece (σ_y k_x — already in R6's `h_SO`) and a **purely
inter-subband** piece (σ_x k_y). For hard-wall modes ψ_n = √(2/W) sin(nπy/W):

> ⟨n|∂_y|m⟩ = **4nm / (W(n² − m²))** for n − m **odd**, and **0** for n − m
> even (including n = m).

I verified this closed form numerically against direct quadrature — exact
agreement at all 16 pairs for n, m ≤ 4, including the vanishing of every
even-difference element and of the diagonal. It is antisymmetric under n ↔ m,
as required for i∂_y to be Hermitian.

Required checks: the diagonal vanishes; even-difference elements vanish;
odd-difference elements match the closed form; the assembled coupling block is
Hermitian; and — the physics check — **turning the coupling on must reduce the
gap near a lobe crossing**, reproducing the Lim et al. hybridization
qualitatively. A version that got the parity selection rule wrong would still
pass PH symmetry and still look plausible, which is exactly why this needs its
own test.

### 5.5 The invariant reproduces the literature's odd/even rule

In the decoupled limit the total TRIM Pfaffian sign must equal the **product**
of the M per-subband signs, so that the topological phase is exactly the
odd-count condition of Potter–Lee and Lutchyn–Stanescu–Das Sarma (§2.1–2.2).
This validates the generalized invariant against an independently established
result rather than against itself.

### 5.6 Multi-lobe boundary vs analytics

Gap closings must land on B = √((μ − ε_n)² + Δ²) for each n, at the same
standard R6 already met — R6 demonstrated agreement to **exactly −0.50 grid
steps** at both μ values tested, and identified that as a finite-sweep-resolution
signature rather than a physical deviation
(`r6_lo_phase_diagram_methods.md` §2.5). Track 4 should meet the same bar and
report the grid step in physical units, per that precedent.

### 5.7 Truncation convergence in M — the hard one, with no analogue

Every prior convergence parameter in this project (N, mesh refinement, θ, tol)
is a *discretization* that converges smoothly. **M is not.** Adding one subband
can flip the total invariant's parity — the topological classification can
change discontinuously with the truncation parameter. So:

- results must be shown stable under M → M+1, M+2, **per quantity**, and
- **a converged-in-M result is not the same as a device-faithful result.** If
  the real device has 10–26 occupied subbands (§3.3), a study truncated at
  M = 4 has *chosen* the parity that determines the topology. This is a
  scientific limitation, not a compute limitation, and no amount of compute
  inside the tractable range removes it.

**Implication for scoping:** a small-M study must be framed as a **structural**
question ("what does inter-subband coupling do to the phase-diagram topology
and to the feedback loop gain?") and **not** as a device-faithful prediction
for Chang-class wires. Framing it the second way would be an overclaim the
compute budget cannot support. §8 scopes accordingly.

### 5.8 Kramers degeneracy at B = 0

Generalize R6's check: exactly 2NM of the 4NM ring eigenvalues must pair to
<1e-8 at B = 0. **R6's specific warning applies with more force here.** R6
found that naive fixed-index pairing gave a false FAIL because the ring has an
*accidental* k ↔ −k degeneracy that persists at B ≠ 0 and gets conflated with
Kramers pairing; it fixed this by counting near-zero adjacent gaps instead
(§2.2). With M subbands there are **more** accidental degeneracies (subband
crossings, and exact degeneracies between decoupled blocks), so the naive
method would fail worse. Reuse R6's counting fix; do not re-derive it.

### 5.9 Serial-vs-parallel equivalence

`bdg_parallel.check_equivalence` at 1-thread BLAS on both sides, expecting
`max_abs_diff = 0.0`, as every track since R8 has done.

---

## 6. Task 2.3 — what genuinely new question does this answer?

R1–R9 have converged on a null: **inductive EM self-feedback does not move the
bulk topological boundary.** R4→R5→R5-revised→R7 showed the shift is a
finite-size artifact vanishing as N⁻¹·⁵ to N⁻³·¹; R8 showed Z_r is largely
redundant with the gap; R9 Track 1 showed the g_crit "pole" is bounded, not
divergent; R9 Track 2 showed temperature does not matter at 20 mK. R9 Track 3
is the only non-null, and it is about disorder.

**Every one of those results was computed with exactly one conducting channel.**

### 6.1 The primary new question — does the null survive more channels?

This is the case for Track 4 being its own track rather than a footnote on an
existing one.

The entire null rests on the self-field being **too weak to matter**, because
it is sourced by the current in a single channel. The chain is: loop current
→ self-flux φ\* → the dimensionless loop gain **s** (R5's bifurcation slope,
R8's L_eff = L_geo/(1−s), R9 Track 1's bounded enhancement 1/(1−s)).

In a wire with M occupied subbands, each channel carries current, so **the loop
current at a given phase twist should grow roughly ∝ M**, and therefore
**s should grow roughly ∝ M.** If that holds even approximately:

- R9 Track 1 measured s plateauing at **~0.825**, with the static fixed point
  detaching at **g ≈ 8.67×10⁴** — *above* the production coupling g = 2.9×10⁴.
  That ordering is why the dynamical instability has been a curiosity rather
  than a problem.
- At M ≈ 10 (§3.3's realistic occupancy) the **same** instability would be
  reached at roughly **g ~ 10⁴** — i.e. **at or below the production coupling**,
  putting a real device *inside* the regime where no static self-consistent
  boundary exists.

**That is a concrete, quantitative, falsifiable route by which this project's
central null result could fail** — and no other track tests it, because Tracks
1–3 all hold the channel count at one. It is also *cheap to test* (§7), which
is what makes the track viable at all.

I want to be careful about the status of this: **it is a hypothesis with a
plausible mechanism, not a result.** The ∝ M scaling is a heuristic (channels
near their band bottoms contribute less than channels well above; inter-subband
coupling redistributes current). Measuring whether s actually scales with M —
and with what exponent — *is the experiment*, and it is well-posed precisely
because both outcomes are informative: linear scaling means the null is a
single-channel artifact; sub-linear or flat scaling means the null is robust and
the project's headline gets *stronger*, not weaker.

### 6.2 Secondary — perturbing a multicritical point, not a generic one

Multi-subband turns one smooth boundary into a set of overlapping lobes
(§2.1–2.2), which **cross**. At a crossing the minigap is anomalously small, so
susceptibility to *any* perturbation is maximal there. Every prior track applied
its perturbation (EM feedback, temperature, disorder) at a **generic point on a
single smooth boundary**. Multi-subband is the first opportunity to apply them
at a **multicritical point** — Lutchyn–Stanescu–Das Sarma's "sweet spots". This
is structurally new and cannot be folded into an existing track, because the
structure being probed does not exist at M = 1.

### 6.3 Tertiary — the link to Track 3

Woods, Das Sarma & Stanescu state that *"high subband occupancy by itself works
as an effective emergent disorder"* (§3.3). R9 Track 3 found the project's only
non-null: at realistic disorder the self-consistent solver's convergence rate
collapses from 100% to 5–20%. If high subband occupancy acts as emergent
disorder, these may be **the same effect reached by two routes** — a testable
connection, and a reason to run Track 4's occupancy axis against Track 3's W
axis. Flagged as a connection worth testing, **not** as a claim.

### 6.4 What would *not* justify a separate track

For completeness: simply re-plotting R6's phase diagram with lobes in it is
**not** sufficient justification. That result is already in the literature
(§2.2), and reproducing it is a validation step (§5.5–5.6), not a contribution.
The contribution is §6.1 and §6.2 — applying *this project's* feedback,
impedance and stability machinery to a multi-channel wire, which is what nobody
in the cited literature does.

---

## 7. Minimum viable scope

Designed to fit the compute reality of §4, in dependency order. Times are at
4 workers, from the §4 cost model, and **do not** include R7's ~5× throttling
factor — apply it mentally to anything over an hour.

| # | Item | N | M | Est. compute | Gate |
|---|---|---|---|---|---|
| 0 | Record the §3.5 caveat against R6/R8 | — | — | **zero** | — |
| 1 | Mode-space builder + validation suite §5.1–5.5, §5.8–5.9 | 40–60 | 1–3 | minutes | **M=1 bit-for-bit or stop** |
| 2 | Multi-subband phase diagram, 40×40, incl. lobe structure & analytic check §5.6 | 60 | 1→4 | ~0.7 h cumulative | truncation scan §5.7 |
| 3 | **Loop gain s vs M** (§6.1) — the headline probe | 80 | 1→4 | **~3.5 min/point** | s(M=1) reproduces R9-T1 |
| 4 | Gap & feedback at a lobe crossing vs a generic point (§6.2) | 60–80 | 3–4 | ~1–2 h | — |
| — | *Full self-consistent B_c sweep vs M* | 80 | ≥3 | **35–160 h — EXCLUDED** | — |

Item 3 is the crux and is **cheap**, which is the finding that makes this track
viable. The loop gain s is the *slope of the self-consistency map at a fixed
operating point* — it needs one root-solve plus a finite-difference derivative
(~30 `eigh`), **not** the full B_c bisection (~2000 `eigh`). That is a ~65×
saving and it is what separates item 3 from the excluded row:

| M | s-slope probe (N=80, ~30 eigh) | full B_c point (~2000 eigh) |
|---:|---:|---:|
| 1 | 2.1 s | 150 s |
| 2 | 16.9 s | 1200 s |
| 3 | 57.0 s | 4050 s |
| 4 | 135 s | 9600 s |
| 5 | 264 s | 18750 s |

So the single measurement most likely to overturn the project's null costs
**minutes**, while the exhaustive version costs **days**. Scope to the former.

Items 0–3 total well under a day of compute. **Item 4 is optional**; items 0–1
alone are worth doing even if the track is otherwise abandoned, since item 0 is
free and item 1 is the reusable infrastructure.

---

## 8. Recommendation

**Conditional yes — at the bounded scope of §7 items 0–3, and explicitly not
beyond it.**

**Why yes:**
- The literature is established, specific, and verified (§2), and one source
  (Lim et al.) says the naive picture is *wrong* in a way that makes the
  question non-trivial rather than confirmatory.
- The single-subband assumption is a genuine, documentable weakness in R6/R8
  (§3), and saying so costs nothing.
- There is one sharp, falsifiable, project-specific question — **does the
  EM-feedback null survive more channels?** (§6.1) — that no other track can
  ask, and it happens to be answerable in **minutes**, not days (§7 item 3),
  because the loop gain does not require the expensive B_c bisection.
- The tractable range (M ≤ 4–6 for phase diagrams) is enough to exhibit the
  lobe structure and inter-subband hybridization the literature describes.

**Why bounded — stated as plainly as the brief asked:**
- **Multi-subband self-consistent EM feedback is not tractable here.** M = 3
  turns R9 Track 3's 77-minute run into ~35 hours before throttling; M = 5 into
  ~6.7 days. The M³ growth is inside a single dense `eigh` and `bdg_parallel`
  cannot touch it (§4.4). This is out of reach, not merely expensive.
- **Parallelism offers no relief.** Six cores buy a flat ~2× on a single large
  solve (33% efficiency), and the existing across-grid-point strategy — already
  saturated at ~2.3–4× — remains strictly better. There is no configuration
  change that recovers any part of M³.
- **Truncation in M is not a convergence parameter in the usual sense** (§5.7).
  A small-M result is a *structural* statement, not a device-faithful one, and
  must be written up that way. With 10–26 subbands occupied in real devices,
  a study at M = 4 has chosen the parity that sets the topology.

**The single most valuable deliverable in this document costs zero compute:**
the §3.5 caveat. If Track 4 is deprioritised for any reason, §3.5 should still
be recorded against R6/R8, alongside the existing fixed-Δ idealisation note.

**Suggested decision rule:** do §7 items 0–1 (free + minutes). Then run item 3
at M = 1, 2, 3. If s grows roughly ∝ M, the multi-subband question is
*load-bearing* for this project's headline and items 2 and 4 are justified. If
s is flat in M, the null is more robust than a single-channel model could
show — which is a **publishable strengthening** of the existing result, obtained
for minutes of compute, and Track 4 can then be closed honestly rather than
abandoned silently.

---

## Appendix — provenance of every number in this document

- **Literature**: arXiv abstract pages for 1007.4569, 1011.6371, 1008.0629,
  1106.3078, 1202.5057, 1411.6255, 1910.04362, fetched this session. Full-text
  detail for 1106.3078 and 1910.04362 came from the **ar5iv HTML rendering**,
  flagged inline as a secondary rendering; abstracts are from the primary arXiv
  record. The published PDF of 1106.3078 did not render extractably.
- **Chang et al. 2015 device numbers** (Δ\* = 190 μeV, 20 mK, B_c ~75 mT,
  n = 5×10¹⁸ cm⁻³, l_e = 100 nm): quoted from the paper. **No wire diameter is
  stated in that paper** — §3.2 estimates from the measured density instead and
  says so.
- **Subband energetics** (§3.2, §3.4): closed-form hard-wall levels (1D
  square-well and cylindrical Bessel-zero) at m\* = 0.023 mₑ / 0.026 m₀.
  The lattice cross-check reproduced R6's t_hop = 16.565 meV exactly from
  ħ²/2m\*a² — confirming the arithmetic is on R6's own footing.
- **Inter-subband matrix element** (§5.4): closed form verified against direct
  numerical quadrature, exact at all 16 pairs with n, m ≤ 4.
- **Timings** (§4.2–4.6): `scipy.linalg.eigh` on random complex Hermitian
  matrices on this machine (10 physical cores = 6 P + 4 E, 16.9 GB, 2.5 GB free
  at benchmark time), at 1 and 6 BLAS threads, 3 repetitions for D ≤ 960 and 1
  above. These are **hardware measurements on random matrices, not physics
  runs.**
- **Repo-measured costs** used as anchors: ~144–161 s/task for the
  self-consistent engine at N = 80 and ~0.05 s/task for the cheap engine
  (`r9_disorder_methods.md`); R9 Track 3 full run 77.1 min (ibid.); R7's N = 120
  cap and ~5× throttling (`ISEF R7 data folder/R7_SUMMARY.md`); R9 Track 1's
  s ≈ 0.825 plateau and g ≈ 8.67×10⁴ detachment (`README.md`, R9 Track 1).
- **Not measured, explicitly a hypothesis**: the s ∝ M scaling in §6.1. That is
  the proposed experiment, not a result.
