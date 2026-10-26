# majorana-nanowire-em-feedback

BdG eigensolver, self-consistency solvers, validation scripts, and data for self-consistent EM feedback in Majorana nanowires.

This repository accompanies the paper in [`paper/`](paper/) (build with `cd paper && latexmk -pdf main.tex`; a prebuilt
[`paper/main.pdf`](paper/main.pdf) is included). It holds the code behind the numbers and figures in the paper, the raw
numerical results, and the methods notes in which each calculation, its validation gates and its corrections are
recorded. Every figure of the paper can be redrawn from the saved results by a script in `src/` or `scripts/`; the
few files that no script writes are listed in [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md).

**Scientific question.** Standard treatments of the Majorana nanowire fix the topological phase boundary and compute the
induced supercurrent afterwards; the current is never fed back into the Hamiltonian that produced it. Here the loop is
closed: a particle-hole-symmetric Bogoliubov-de Gennes (BdG) eigensolver, a real-space Pfaffian Z<sub>2</sub> invariant and a
Biot-Savart ring self-field are combined so that the wire's own flux dresses the Peierls phase in its Hamiltonian, and the
self-consistent shift of the topological boundary is computed.

## Summary of findings

All statements below are made at the validated feedback strength unless noted, on a lattice model of N sites. The
methods notes in [`docs/methods/`](docs/methods/) give the full derivations, validation gates and the corrections made
along the way, including the superseded intermediate numbers.

- **Single channel: no bulk shift.** The headline shift is at the top of the dynamically stable regime (g ≈ 2.9 × 10⁴,
  N = 80): ΔB<sub>c</sub>/B<sub>c</sub> ≈ −0.064 % (μ = 1.2) and −0.066 % (μ = 1.4). It decays as ~N<sup>−1.8</sup> with system size,
  so it is a finite-size artefact (an Aharonov-Bohm effect on the finite ring) and vanishes in the bulk. The result is
  unchanged along the wire-length / inductance axis (R7) and at the operating temperature (R9, Track 2).
- **Detector choice matters.** A conventional spectral-gap threshold detector overstates the shift by about an order of
  magnitude relative to the threshold-free Pfaffian invariant (R4).
- **Existence versus stability.** The self-consistency loop has a dynamical bifurcation above which a fixed point still
  exists but is not Picard-attracting (R5). Existence and dynamical stability have to be checked independently.
- **Reactive impedance.** The feedback-renormalised reactance Z<sub>r</sub> is largely redundant with the bulk gap. Its
  apparent divergence near g<sub>crit</sub> is bounded (≈ 5.7× at N = 80) and ends when the static fixed point detaches; the cap
  grows under mesh refinement, so a continuum divergence is neither demonstrated nor excluded (R8, R9 Track 1).
- **Finite temperature.** At 20 mK the thermal quasiparticle population is ~10<sup>−48</sup>, and the null result survives. Thermal
  effects appear only far above the operating temperature (R9, Track 2).
- **Disorder.** On-site disorder smears the topological classification. A self-consistent shift is real at W ≈ 10 Δ but does
  **not** survive the N-independence check at W ≈ 46 Δ (R9, Track 3).
- **Multiple subbands: not a null.** For a second occupied transverse subband (M ≥ 2) a pre-registered loop-gain diagnostic
  crosses the Picard stability threshold at the production coupling. A full self-consistent boundary calculation at M = 2
  finds a shift at N = 80, but N = 60 and N = 120 return no self-consistent solution, so the *existence* of a static
  boundary is fragile to the discretisation; the M = 1 to 6 trend is non-monotonic (R9, Track 4).

**Scope.** Single-subband results are for an idealised regime: fixed gap, quasi-static feedback, a one-dimensional
lattice model, parameters of InAs/Al devices from the literature. See the paper's discussion and
[`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md) for the limitations.

## Repository layout

```
.
├── src/                    21 Python modules (flat; they import each other by module name)
├── results/                data and figures written by the modules, one folder per topic
│   ├── baseline_bdg_pfaffian/             BdG eigensolver and Pfaffian invariant (R1, R2)
│   ├── subspace_update/                   PH-symmetric rank-2 subspace update (R3)
│   ├── em_feedback_detector_comparison/   Biot-Savart coupling, threshold vs Pfaffian detector (R4)
│   ├── self_consistency_bifurcation/      self-consistent ring, Anderson mixing, shift vs coupling (R5)
│   ├── lutchyn_oreg_phase_diagram/        Lutchyn-Oreg (B, mu) phase diagram (R6)
│   ├── inductance_scan/                   wire length / inductance axis (R7)
│   ├── reactive_impedance/                Z_r: static, Kubo, self-consistent spin-orbit (R8)
│   │   └── pole_characterization/         the apparent Z_r pole near g_crit (R9, Track 1)
│   ├── finite_temperature/                (R9, Track 2)
│   ├── disorder/                          (R9, Track 3)
│   └── multisubband/                      (R9, Track 4)
├── docs/
│   ├── methods/            methods notes with derivations, gates and correction records
│   ├── paper_audit/        claim-by-claim verification log for the paper
│   ├── REPRODUCING.md      what to run, how long it takes, what it overwrites
│   ├── KNOWN_ISSUES.md     limitations, provenance gaps, caveats
│   └── RENAME_MAP.md       old-to-new path map and the R-code glossary
├── paper/                  main.tex, references.bib, main.bbl, figures/, main.pdf (self-contained)
├── scripts/                redraw_*.py: redraw figures from the saved results (seconds each)
├── requirements.txt        16 pinned packages
├── CITATION.cff
└── LICENSE
```

The research was carried out in numbered stages (R1 to R9); file names keep that prefix. The table of codes and the
directory each one maps to is in [`docs/RENAME_MAP.md`](docs/RENAME_MAP.md).

## Quickstart

Tested with **Python 3.13.3** on Windows 11. Python 3.10 or newer is expected to work (the source uses PEP 604 type
annotations, so Python 3.9 fails to import it); other Python versions and other operating systems have not been tested.

```bash
git clone https://github.com/MidnightNotHere/majorana-nanowire-em-feedback.git
cd majorana-nanowire-em-feedback
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
python -X utf8 src/bdg_eigensolver_r1_r2_LOAdapted.py
```

The last command runs the baseline BdG eigensolver and Pfaffian validation (about one minute) and writes seven data files
and seven figures to `results/baseline_bdg_pfaffian/`. `-X utf8` is needed on a Windows console whose code page is not UTF-8,
because the scripts print Greek letters.

Each script writes into `results/<topic>/` of the checkout it lives in, **overwriting the committed files of that topic**.
To keep the committed data intact, run in a copy of the repository.

Build the paper (needs a LaTeX distribution with `revtex4-2`, `latexmk` and TikZ):

```bash
cd paper
latexmk -pdf main.tex
```

## Reproducing the results

[`docs/REPRODUCING.md`](docs/REPRODUCING.md) lists every script with its arguments, its approximate run time, its outputs
and whether it is safe to run for a quick check. Most scripts finish in about a minute or less. The self-consistent
calculations are long: from about 15 minutes up to 4 hours for the multi-subband extension. Their results are stored in
`results/` so that the figures can be redrawn without rerunning them.

## Documentation

- [`docs/methods/`](docs/methods/): methods notes, one per major calculation. Corrections are recorded in place and the
  superseded numbers are kept.
- [`docs/paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md`](docs/paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md):
  where each claim in the paper comes from.
- [`docs/KNOWN_ISSUES.md`](docs/KNOWN_ISSUES.md): limitations and known gaps.

## Licence and citation

Released under the [MIT licence](LICENSE). The licence covers the code, the data and the results in this repository.

To cite this work use [`CITATION.cff`](CITATION.cff) (GitHub's "Cite this repository" button reads it). A DOI will be added
to it when the first release is archived.

Author: Arnav Rohilla.
