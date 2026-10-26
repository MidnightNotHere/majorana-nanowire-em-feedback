# R9 — Track 2: Finite temperature in the BdG / Pfaffian / Kubo framework

**Module:** `src/bdg_r9_finite_temperature.py` (additive; imports R4/R5/R7/R8A/R8B
and `bdg_parallel` **read-only**, modifies nothing in R1–R8).
**Outputs:** `results/finite_temperature/` (`r9_finiteT_data.json`,
`r9_finiteT_summary.png`, `r9_full_run.log`).
**Run:** `py -3.13 src/bdg_r9_finite_temperature.py --full`
(`--validate` for the T→0 checks only; `--quickfull` for a tiny-grid smoke test).

Every prior phase (R1–R8, and R9 Track 1) used **zero-temperature** BdG ground
states (occupation = θ(−E)) and a **T=0** Kubo formula. This track introduces
Fermi–Dirac occupation of the BdG quasiparticle spectrum and asks whether the
project's null / near-null results survive realistic thermal broadening. It is
strictly additive: the T=0 code paths are untouched, and every finite-T quantity
is validated to reproduce them as T→0 before any nonzero-T number is trusted.

---

## Task 1 — the temperature scale comes first

The induced gap is the R6-verified **Δ = 0.19 meV** (Chang et al. 2015,
*Nature Nanotechnology* 10, 232; arXiv:1411.6255, epitaxial-Al full-shell —
"an induced gap of Δ\* = 190 μeV"). This fixes a real, checkable temperature
scale:

$$\frac{\Delta}{k_B} = \frac{0.19\ \text{meV}}{8.617\times10^{-5}\ \text{eV/K}}
= 2.205\ \text{K}.$$

(Using the project's stated convention `k_B T_c ~ Δ`, so `T_c ≈ 2.205 K`; the
BCS relation Δ = 1.764 k_B T_c would give T_c ≈ 1.25 K. Either way the scale is
~1–2 K.)

**Operating temperature (verified, not assumed).** Chang et al. 2015 state:
*"Measurements were carried out in a dilution refrigerator with a base
temperature of 20 mK."* (verified against the primary source, ar5iv/arXiv
1411.6255; no separate electron temperature is reported). So the realistic
operating point is **T = 20 mK**.

The dimensionless temperature used throughout (energies are in units of Δ) is

$$\theta \equiv \frac{k_B T}{\Delta},\qquad
\theta(20\ \text{mK}) = 0.0091.$$

| T | θ = k_BT/Δ | T/T_c | f(E=Δ) = thermal QP weight |
|---:|---:|---:|---:|
| 10 mK | 0.0045 | 0.0045 | 1.8×10⁻⁹⁶ |
| **20 mK (operating)** | **0.0091** | **0.0091** | **1.3×10⁻⁴⁸** |
| 50 mK | 0.0227 | 0.0227 | 7.1×10⁻²⁰ |
| 100 mK | 0.0454 | 0.0454 | 2.7×10⁻¹⁰ |
| 300 mK | 0.1361 | 0.1361 | 6.4×10⁻⁴ |
| 500 mK | 0.2268 | 0.2268 | 1.2×10⁻² |

**Ratio and verdict (drives the rest of the investigation).** At the operating
point **T/T_c ≈ 0.009** — thermal effects are a *tiny* perturbation. The single
number that decides it is the thermal population of a full-gap quasiparticle,
`f(E=Δ) = 1/(e^{1/θ}+1) ≈ e^{−Δ/k_BT}`: at 20 mK this is **~10⁻⁴⁸**, i.e.
utterly negligible. Thermal effects cannot matter at operating conditions unless
something is anomalously sensitive; the swept θ range is therefore chosen to
**bracket the operating point and climb to where thermal effects *do* turn on**
(θ up to ~0.15, i.e. up to ~330 mK, ≈16× the operating point) so the crossover
scale is located, not assumed.

**Model-validity ceiling.** These are *fixed-Δ* BdG calculations: the pairing is
held at its T=0 value. That is fine while T ≪ T_c but becomes unphysical as
T → T_c/2 (real pairing is thermally suppressed). We therefore trust finite-T
results only up to **θ ≈ 0.15** (≈330 mK) and flag anything above as outside the
model's domain. This ceiling sits ~16× above the operating point, so it does not
limit any operating-condition conclusion.

---

## Task 2 — finite-temperature occupation (implementation)

Temperature enters through the **occupation** of the BdG spectrum, via the
numerically-stable Fermi function `fermi_dirac(E, θ)` = 1/(e^{E/θ}+1), which
returns exactly the T=0 step θ(−E) = (E<0) at θ=0 (matching the R1–R8
"occupied == E<0" convention bit-for-bit).

### 2.1a Self-field / current channel (boundary, pole)
`ring_current_density_T` is the thermal generalization of R5
`ring_current_density`: each quasiparticle state is weighted by f(E_n) and the
sum runs over the **full** BdG spectrum,
`I = (2et/ℏ) Σ_n f(E_n) Im[…]`. At θ→0 this is R5's sum over occupied states,
**bit-for-bit**. It drops straight into the R5/R7 Biot–Savart self-field loop
via `scalar_ring_map_T` (thermal port of `scalar_ring_map`), whose principal
self-consistent fixed point `_phi_star_T` / `operating_point_T` are thermal
ports of R7 `_phi_star_a` / R8A `operating_point`. **The loop slope s (Track A's
observable, and R5's bifurcation slope) is finite-differenced from the actual
thermal current, so it is exact at every θ — no operator-convention subtlety
enters this channel.**

### 2.1b Boundary detection at finite T — the invariant is a T=0 concept
The topological invariant (Pfaffian Z₂ sign) is a property of the **BdG
Hamiltonian's ground state (fermion parity)**, *not* of the occupation — it does
not depend on temperature at fixed Hamiltonian. So at finite T the boundary is
defined two complementary ways (both implemented):

- **(A) Hard boundary** — `ring_Bc_selfconsistent_T`: the *same* R5/R7 warm-hint
  Pfaffian bisection, but the invariant is evaluated at the **finite-T
  self-consistent self-field** φ\*(T). **Temperature enters only through
  φ\*(T).** This is the physically meaningful channel for this project (does the
  now-thermal inductive self-field move the bulk boundary). At θ→0 it reproduces
  R7 `ring_Bc_selfconsistent_a` (hence R5) bit-for-bit.
- **(B) Soft crossover** — `thermal_parity_indicator`: the thermal expectation of
  the ground-state parity for a two-level (±E_gap) gap,
  `P(B) = sgn_Pf(B)·tanh(E_gap(B)/2θ)`. Since at T>0 the free energy is analytic
  in B, the sharp Z₂ jump becomes a **crossover of width ~θ** around the
  gap-closing; the crossover width is a genuinely new finite-T observable with no
  T=0 analogue. Its zero-crossing coincides with the Pfaffian flip.

### 2.1c Kubo current–current correlator at finite T
`kubo_chi_T` is the thermal Lehmann sum,
`χ(w) = χ_dia + Σ_{m≠n}(f_n−f_m)|⟨m|Ĵ|n⟩|²/(E_n−E_m+w+iη) + D·iη/(w+iη)`,
with thermal occupations f_n = f(E_n), thermal diamagnetic term
`χ_dia = Σ_n f_n⟨n|Ĵ_dia|n⟩`, and an **intraband thermal Drude weight**
`D = Σ_{|E_n−E_m|<tol} (df/dE)|⟨m|Ĵ|n⟩|²`. At θ→0 the (f_n−f_m) factor is ±1 only
for occupied↔empty pairs (R8-B's `|df|>0.5` mask), χ_dia → Σ_{E<0}, and **D = 0**
(df/dE=0 for every gapped state) — so `kubo_chi_T` reduces to R8-B `kubo_chi`
bit-for-bit. At θ>0 the continuous weights open **sub-gap channels from
thermally-excited quasiparticles**, and D turns on a Drude peak.

**A subtlety that is itself a finding — the finite-T Kubo needs the *full*
current operator.** The finite-T DC gate (χ(0) must equal `d⟨Ĵ⟩_T/dφ`) exposed
that R8-B's deliberately *kinetic-only* current cannot satisfy the gate at T>0:
the thermal-reoccupation term in `d⟨Ĵ⟩_T/dφ` samples `dE_n/dφ = ⟨n|∂H/∂φ|n⟩`,
the **full** (kinetic + spin-orbit) current. With kinetic-only Ĵ the gate fails
by ~17% at θ=0.3; with the full current Ĵ = ∂H/∂φ (`full_current_operators`,
verified against numerical dH/dφ, d²H/dφ²) it holds to **~10⁻⁶ at all θ**. The
finite-T Lehmann sum therefore uses the full current; the residual kinetic-vs-
full difference is exactly R8-B's already-quantified spin-orbit band. (The
self-field channel keeps R5/R7's kinetic convention unchanged, as above.)

---

## Task 2 validation — the T→0 limit reproduces R5/R7/R8 (do this first)

`validate_T0_limit()` (N=80, μ=1.2, B near the transition). **All pass:**

| quantity | finite-T vs existing T=0 | result |
|---|---|---|
| ring current | `ring_current_density_T(θ=0)` vs R5 `ring_current_density` | **0.0** (bit-exact) |
| operator/electron-block current proportionality | ratio θ-independence | **5.5×10⁻¹⁴** rel |
| self-consistency map | `scalar_ring_map_T(θ=0)` vs R5 `scalar_ring_map` | **0.0** (bit-exact) |
| operating point (φ\*, s) | `operating_point_T(θ=0)` vs R8-A `operating_point` | **0.0 / 0.0** |
| Kubo χ(w) | `kubo_chi_T(θ=0)` vs R8-B `kubo_chi` | **8.9×10⁻¹⁶** max |
| finite-T DC gate | χ(0) vs `d⟨Ĵ⟩_T/dφ`, θ=0/0.009/0.05/0.1/0.15 | **2.5e-7 / 3.9e-7 / 8.5e-7 / 2.2e-6 / 4.8e-6** (all pass <1e-3) |

The DC gate is verified at **multiple θ** (not just θ=0), per the task — it is a
falsifiable, per-θ check that the finite-T Lehmann implementation is correct.

---

## Task 3 — parallelization, BLAS, and load balance

The natural parallel axis is **temperature** (each (μ,B,θ,g) point is
independent). The existing `src/bdg_parallel.py` is reused unchanged
(process pool, one BLAS thread per worker); the finite-T grid workers
(`_bc_worker`, `_pole_worker`, `_timing_worker`) are top-level/picklable for
Windows spawn.

- **Serial-vs-parallel equivalence** (`check_equivalence` on a (g,θ) subset):
  **bit-identical, max_abs_diff = 0.0** — parallelization changes no number.
- **Process-parallel speedup:** with RAM available, 6 workers ran the 6-task
  subset in **6.3 s** vs **~16 s** serial — a **~2.5×** speedup (sub-linear
  because only 6 tasks over the 6 P-cores, plus spawn overhead; consistent with
  the R8/R9-T1 sweet spot of ~4–5 workers). The benchmark demonstrates the pool
  is available and bit-identical; the production sweep here was run **serially**
  (`--workers 1`) for determinism and memory-safety on this RAM-constrained
  laptop, with 1-thread BLAS (the fastest single-process config).
- **BLAS oversubscription (measured at N=80, "use whichever is fastest"):** in a
  clean fresh-process comparison, 1-thread BLAS = **9.4 s** vs default
  multithreaded BLAS = **20.4 s** for the same 3 operating-point solves —
  **default multithreaded BLAS is 2.2× SLOWER**. These matrices (4N = 320) are
  too small to amortize BLAS threading; the cross-thread reduction overhead
  dominates. **Config used: 1-thread BLAS per worker** (as `bdg_parallel`
  enforces). Consistent with the R9-Track-1 measurement. (The in-run benchmark's
  serial-1-thread vs serial-default came out ≈equal only because the run was
  launched under a `*_NUM_THREADS=1` environment, which OpenBLAS reads once at
  import — so both branches were effectively single-threaded; the fresh-process
  numbers above are the authoritative comparison.)
- **Load balance / per-point cost vs θ:** the finite-T sum includes **all**
  states at every θ (no adaptive Fermi-tail truncation), so per-point cost is
  **θ-independent** — measured spread ratio **1.15** across θ∈{0,…,0.15} (if
  anything slightly *cheaper* at higher θ). Static equal-chunking is therefore
  optimal; **no work-stealing / dynamic scheduling is needed** — the specific
  concern raised in the task (a wider Fermi tail making high-T points dominate)
  is checked and found not to apply, because no tail truncation is used.
- **Memory note:** `default_workers()` is memory-bound (each worker re-imports
  the ~0.4 GB R4–R7 stack); it caps to 1 worker when free RAM < ~0.5 GB. This
  laptop hovers at 1–3 GB free under load, so the effective worker count varies
  run-to-run (1 when tight, up to 6 when free); the equivalence/BLAS/load-balance
  conclusions are RAM-independent.

---

## Task 4 — do R8's key findings survive finite T?

Full run: N=80, μ=1.2, g=2.9×10⁴, θ∈{0, 0.009 (operating), 0.05, 0.1, 0.15}
(`r9_finiteT_data.json`, `r9_finiteT_summary.png`). The T=0 anchor of the
boundary solver is verified **bit-for-bit** against R7's shipped
`ring_Bc_selfconsistent_a` (|ΔB_c| = 0.00×10⁰; both give −0.13190% at N=80), so
the finite-T comparison is internally exact.

### 4.1 Bulk topological boundary shift vs T — **null survives**
| θ (T) | ΔB_c/B_c at N=80 |
|---|---|
| 0 | −0.1319 % |
| **0.009 (20 mK, operating)** | **−0.1313 %** |
| 0.05 | −0.1313 % |
| 0.1 | −0.1295 % |
| 0.15 (≈330 mK) | −0.1206 % |

At the operating point the shift is **indistinguishable from T=0** (changes by
0.0006 percentage points, ~0.4 % of an already-tiny value). Even at θ=0.15 (16×
operating) it only weakens ~9 % — thermal broadening reduces the current, hence
the self-field, hence the (finite-size) shift. **The project-standard
N-independence check holds at every T:** the shift decays toward zero with system
size at both operating and elevated T —

| N | θ=0.009 | θ=0.15 |
|---|---|---|
| 80 | −0.131 % | −0.121 % |
| 100 | −0.086 % | −0.080 % |
| 120 | −0.061 % | −0.057 % |

(~N⁻¹·⁹, the same finite-size Aharonov–Bohm decay R5 found at ~N⁻¹·⁸). So finite
temperature does **not** create a bulk boundary shift; it leaves the R4–R7 null
result intact (and, if anything, slightly more null). *(Boundary definition (B),
the soft thermal-parity crossover, likewise has its zero-crossing pinned at the
same B_c; its only new content is a crossover of width ~θ, which at operating θ
is far narrower than the finite-size B_c uncertainty.)*

### 4.2 The g_crit pole (R9 Track 1) vs T — **washes out**
Peak principal-branch enhancement 1/(1−s) near R5's g_crit band, vs T:

| θ | max 1/(1−s) |
|---|---|
| 0 | **5.94** |
| 0.05 | 3.79 |
| 0.1 | 2.78 |
| 0.15 | 2.35 |

The T=0 value (5.94) reproduces R9-Track-1's bounded ~5.7–5.9× enhancement.
Thermal broadening **reduces** the loop slope s (the thermal current response is
smaller), so the enhancement drops sharply and the detachment/pole is pushed to
higher g: the pole **washes out** with temperature, monotonically. At the
operating point (θ=0.009 ≪ 0.05) the enhancement is essentially the T=0 value —
the pole is intact at 20 mK but is **not sharpened** by temperature; it is a
soft, bounded, T=0-limit feature that thermal broadening only softens further.
This *reinforces* R9-Track-1's conclusion that the pole is not a sharp new
observable.

### 4.3 The dissipation spectrum (R8-B) vs T — **new sub-gap absorption at high T only**
Sub-gap Re Z peak (ℏω/Δ < 1) and the intraband Drude weight D vs T:

| θ | sub-gap peak Re Z | Drude weight D |
|---|---|---|
| 0 | 8.4×10⁻⁵ (η-tail only) | 0 |
| 0.05 | 1.0×10⁻⁴ | −5.3 |
| 0.1 | 1.6×10⁻⁴ | −46 |
| 0.15 | 4.3×10⁻⁴ | −216 |

At T=0 the sub-gap Re Z is only the η-broadened tail of the interband feature
(no true sub-gap absorption, consistent with R8-B). At finite T two new,
genuinely-T=0-absent features appear: (i) an **intraband Drude peak** (weight D,
zero at T=0, growing steeply with T) from thermally-excited quasiparticles, and
(ii) **sub-gap absorption** growing ~5× by θ=0.15, which also smears the sharp
interband onset. Both are **negligible at the operating point** (θ=0.009: the
thermal QP population is ~10⁻⁴⁸) and only switch on around θ≳0.05 (≈110 mK). This
is the one place finite T adds *qualitatively new* physics — but only well above
operating temperature.

---

## Bottom line

**Finite temperature does not matter at realistic operating conditions.** At the
Chang et al. 2015 base temperature (20 mK, θ = k_BT/Δ ≈ 0.009, T/T_c ≈ 0.009) the
thermal quasiparticle population is ~10⁻⁴⁸ and every R8-era result is unchanged
to well within its finite-size uncertainty: the bulk boundary shift stays at the
same finite-size-→0 null (−0.131 % vs −0.132 % at T=0, decaying ~N⁻¹·⁹ toward
zero at all T), the g_crit "pole" keeps its bounded ~6× enhancement (and is only
*softened*, never sharpened, by heat), and the Kubo dissipation is unchanged.

Thermal effects turn on only far above operating temperature (θ ≳ 0.05, ≈110 mK):
there the boundary shift weakens further (more null), the pole washes out, and a
new intraband Drude peak + sub-gap absorption appear in the dissipation — but this
regime approaches the model's fixed-Δ validity ceiling (θ ≈ 0.15) and is ~10–15×
hotter than any real device. **This is a null result, consistent with the entire
R4–R9 lineage: the finite-temperature machinery is implemented and validated to
the T=0 limit bit-for-bit, and it confirms — does not overturn — the project's
null / near-null findings.** The one honest new physics addition (thermal sub-gap
dissipation) exists only outside operating conditions.

*(Implementation note preserved as a finding: the finite-T Kubo requires the
FULL current operator Ĵ = ∂H/∂φ, not R8-B's kinetic-only Ĵ, because the thermal
reoccupation term samples the full current — a kinetic-only operator cannot
satisfy the finite-T ω→0 gate. See Task 2.1c.)*
