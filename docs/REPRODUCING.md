# Reproducing the results

This guide says what each script does, how long it takes, what it writes, and how the committed data and figures relate
to it. The run times were measured on one machine (Windows 11, 16 logical cores, Python 3.13.3, the packages pinned in
`requirements.txt`); other hardware will differ. Where a time comes from a record in the repository rather than from a
fresh measurement, the source is named.

## 1. Setup

```bash
python -m venv .venv
. .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -r requirements.txt
```

Only Python 3.13.3 has been tested. Python 3.10 or newer is expected to work, because the source uses PEP 604
annotations (`X | Y`); Python 3.9 fails to import the baseline module. Run the scripts with `python -X utf8` on a Windows
console that is not set to UTF-8: they print Greek letters, and the default code page raises `UnicodeEncodeError`.

The scripts are run from anywhere as `python src/<script>.py`. Each one resolves its output directory from its own
location (`<repo>/results/<topic>/`).

## 2. Read this before running anything

**Every script writes into `results/<topic>/` of the checkout it lives in, and overwrites the committed files there.**
(`src/bdg_r9_multisubband_plot_M6.py` also copies its figure into `paper/figures/`, and each `scripts/redraw_*.py` file
overwrites the committed image it redraws unless it is given an output path as its argument.)
To check a result without touching the committed data, run it in a copy of the repository:

```bash
cp -r majorana-nanowire-em-feedback /tmp/scratch-run     # Windows: robocopy /E
cd /tmp/scratch-run
python -X utf8 src/bdg_eigensolver_r1_r2_LOAdapted.py
```

or restore the committed files afterwards with `git restore results/`.

What to expect from a rerun:

- **Values.** The scripts of section 3 were rerun on the machine described above and reproduced their committed JSON,
  CSV and NPY values exactly (largest absolute difference 0), apart from the timing columns listed next. Other BLAS or
  LAPACK builds can differ at the level of 1e-13 (a different reduction order). The process-parallel sweeps pin BLAS to
  one thread per worker; serial and parallel runs agree bit for bit at equal thread count.
- **Timings.** Files and figures that record wall-clock time (`performance_baseline.csv`, `r3_lo_performance.csv`,
  `r3_lo_sweep_tracking.csv`, `figure3_performance.png`, `r3_lo_figure2_speedup.png`, `r3_lo_figure3_sweep_tracking.png`)
  change from machine to machine. The non-timing columns are identical.
- **Random numbers.** Disorder realisations use `numpy.random.SeedSequence` with fixed seeds, so ensembles are
  reproducible.

## 3. Scripts that finish in minutes

| Script | Arguments | Time | Writes (in `results/`) |
|---|---|--:|---|
| `src/bdg_eigensolver_r1_r2_LOAdapted.py` | | 63 s | `baseline_bdg_pfaffian/`: 5 CSV, 2 NPY, 7 PNG |
| `src/bdg_eigensolver_r3_lo_adapted.py` | | 51 s | `subspace_update/`: 3 CSV, 3 PNG |
| `src/bdg_r6_lo_phase_diagram.py` | | 53 s | `lutchyn_oreg_phase_diagram/`: gap map and phase map (NPY), `r6_gap_vs_Bz.csv`, 3 PNG |
| `src/bdg_r8_reactive_impedance.py` | | 3 s | redraws `r8_figure1_*`, `r8_figure2_*` and rewrites `R8_TRACKA_SUMMARY.md` from the saved `r8_trackA_data.json`; `--regen` recomputes the data (long) |
| `src/bdg_r9_pole_characterization.py` | | 3 s | redraws `r9_figure1_*`, `r9_figure2_*` and rewrites `R9_TRACK1_SUMMARY.md` from the saved `r9_pole_data.json`; `--regen` recomputes the data (long), `--quick` is a smoke test |
| `src/bdg_r9_multisubband.py` | | 2 s | `multisubband/r9_multisubband_preregistration.json`; runs the construction gates, all of which must pass |
| `src/bdg_r9_multisubband_plot_M6.py` | | 1 s | `multisubband/fig7_multisubband_M1to6_extension.png`, from the saved sweep data; also copies it to `paper/figures/` |
| `src/bdg_r9_finite_temperature.py` | `--validate` | 34 s | T to 0 validation gates only; writes no production file |
| `src/bdg_r9_disorder.py` | `--validate` | 111 s | W = 0 gate only; writes no production file |
| `scripts/redraw_paper_fig3.py` | `[output.png]` | 2 s | `paper/figures/fig3_null_validation.png`, from `self_consistency_bifurcation/r5_physical_shift_results.json` |
| `scripts/redraw_paper_fig4.py` | `[output.png]` | 3 s | `paper/figures/fig4_bounded_pole.png`, from the saved pole data |
| `scripts/redraw_paper_fig6.py` | `[output.png]` | 2 s | `paper/figures/fig6_multisubband_summary.png`, from `multisubband/r9_multisubband_sweep_data.json` |
| `scripts/redraw_r5_figure10.py` | `[output.png]` | 2 s | `self_consistency_bifurcation/r5_figure10_shift_vs_g.png`, from `r5_physical_shift_results.json` and `r5_results.json` |
| `scripts/redraw_r5_figure11.py` | `[output.png]` | 2 s | `self_consistency_bifurcation/r5_figure11_headline_validation.png`, from `r5_physical_shift_results.json` |
| `scripts/redraw_r9_multisubband_summary.py` | `[output.png]` | 2 s | `multisubband/r9_multisubband_summary.png`, from `r9_multisubband_sweep_data.json` |
| `scripts/redraw_r9_disorder_summary.py` | `[output.png]` | 3 s | `disorder/r9_disorder_summary.png`, from `r9_disorder_data.json` (calls `make_figures()` of `src/bdg_r9_disorder.py`) |

The times of the `scripts/` rows include starting Python and importing Matplotlib. They write images only, and each draws
the values stored in the named result file without recomputing anything.

## 4. Long-running scripts

These were not rerun during the repository restructure. The times are those recorded with the committed results
(source in the last column), or measured as a lower bound.

| Script | Arguments | Time | Source of the time | Writes |
|---|---|--:|---|---|
| `src/bdg_r4_em_coupling.py` | | > 10 min in total; figure 6 alone 676 s (an earlier measurement had put it above 15 min); figure 4 > 45 min, not finished; figure 5 not timed | measured | `em_feedback_detector_comparison/`: 6 PNG (figures 1 to 3 are the quick part; R4 saves no data file) |
| `src/bdg_r5_selfconsistent_ring.py` | | about 35 min | earlier README | `self_consistency_bifurcation/r5_results.json`, `r5_shift_comparison.csv` |
| `src/bdg_r5_physically_meaningful_shift.py` | | about 40 min | earlier README | `r5_physical_shift_results.json` (the data of `r5_figure10`, `r5_figure11` and paper Fig. 4; those are drawn by `scripts/redraw_*.py`, section 3) |
| `src/bdg_r5_anderson_mixing.py` | | not recorded | | `r5_anderson_results.json` |
| `src/bdg_r7_inductance_sweep.py` | | > 10 min | measured | `inductance_scan/r7_results.json`, 2 PNG; resumable through the `r7_points_cache.json` cache |
| `src/bdg_r8b_kubo_impedance.py` | | > 10 min | measured | `reactive_impedance/r8b_*`, `R8B_*.md` |
| `src/bdg_r8c_so_selfconsistent.py` | | > 10 min | measured | `reactive_impedance/r8c_*`, `R8C_SO_SUMMARY.md` |
| `src/bdg_r9_finite_temperature.py` | `--full` | 915.5 s | `r9_finiteT_data.json` (`wall_time_s`), `r9_full_run.log` | `finite_temperature/` |
| `src/bdg_r9_disorder.py` | `--full` | 4,627 s | `r9_disorder_data.json` (`wall_time_s`) | `disorder/r9_disorder_data.json`, `r9_disorder_summary.png` |
| `src/bdg_r9_multisubband_sweep.py` | | about 47 min | `docs/methods/r9_multisubband_methods.md` | `multisubband/r9_multisubband_sweep_data.json` (drawn by `scripts/redraw_r9_multisubband_summary.py` and `redraw_paper_fig6.py`) |
| `src/bdg_r9_multisubband_sweep_M56.py` | | 14,920 s (4.1 h) | `r9_multisubband_sweep_M5_M6.json` (`elapsed_s`) | `multisubband/r9_multisubband_sweep_M5_M6.json`, `..._preregistration_M5_M6.json` |
| `src/bdg_r9_multisubband_bc.py` | | 5,711 s (N = 80), 4,257 s (N = 120), 406 s (N = 60) | `docs/methods/r9_multisubband_methods.md`, section 13 table | `multisubband/r9_multisubband_bc_M2_*.json` |
| `src/bdg_r9_multisubband_bracket_scan.py` | | not recorded | | `multisubband/r9_multisubband_bracket_scan.json` |

`src/bdg_parallel.py` is a helper imported by the sweep drivers, not a script to run. Its worker count and memory budget
are fixed for the machine it was tuned on; see [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md).

## 5. Which files the paper uses

The figure files are numbered in the order they were made, not in the order they appear in the typeset paper (the paper's
Figs. 1 and 2 are TikZ diagrams drawn in `main.tex`).

| File in `paper/figures/` | Fig. in the paper | How it is produced |
|---|---|---|
| `fig1_phase_diagram.png` | 7 | byte-identical copy of `results/lutchyn_oreg_phase_diagram/r6_figure3_phase_diagram.png` (`src/bdg_r6_lo_phase_diagram.py`) |
| `fig2_pfaffian_vs_threshold.png` | 3 | byte-identical copy of `results/em_feedback_detector_comparison/r4_figure6_pfaffian_vs_threshold.png` (`src/bdg_r4_em_coupling.py`, `validate_pfaffian_boundary()`, 676 s; rerun on 2026-10-25 it reproduced the file byte for byte) |
| `fig3_null_validation.png` | 4 | `scripts/redraw_paper_fig3.py`, from `results/self_consistency_bifurcation/r5_physical_shift_results.json` (written by `src/bdg_r5_physically_meaningful_shift.py`, long); reproduces the committed image byte for byte |
| `fig4_bounded_pole.png` | 5 | `scripts/redraw_paper_fig4.py`, from `results/reactive_impedance/pole_characterization/r9_pole_data.json` |
| `fig5_disorder_summary.png` | 6 | byte-identical copy of `results/disorder/r9_disorder_summary.png` (`src/bdg_r9_disorder.py --full`, long; or `scripts/redraw_r9_disorder_summary.py` from the saved data) |
| `fig6_multisubband_summary.png` | 8 | `scripts/redraw_paper_fig6.py`, from `results/multisubband/r9_multisubband_sweep_data.json` (written by `src/bdg_r9_multisubband_sweep.py`, long); reproduces the committed image byte for byte |
| `fig7_multisubband_M1to6_extension.png` | 9 | byte-identical copy of `results/multisubband/fig7_multisubband_M1to6_extension.png` (`src/bdg_r9_multisubband_plot_M6.py`) |

The "Generated by" text in the captions of Figs. 4, 5 and 8 names the script that makes the data, not the one that draws
the figure; see [`KNOWN_ISSUES.md`](KNOWN_ISSUES.md), section 4. Three result figures have no `src/` module that draws
them and are redrawn by scripts from the saved data: `r5_figure10_shift_vs_g.png`, `r5_figure11_headline_validation.png`
and `r9_multisubband_summary.png` (section 3 gives the scripts). Their plotting code is a reconstruction; see
`KNOWN_ISSUES.md`, section 3.

[`paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md`](paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md) traces each number
in the text to the file it comes from.

## 6. Building the paper

```bash
cd paper
latexmk -pdf main.tex
```

This needs a LaTeX distribution with `revtex4-2` and `apsrev4-2`, `latexmk`, `bibtex` and TikZ. `paper/` is
self-contained: it needs only `main.tex`, `references.bib`, `main.bbl` and `figures/`. The title block uses
`\date{\today}`, so every rebuild differs from the committed PDF in the date line.
