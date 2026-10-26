# Known issues and limitations

Things a reader or a re-user of this repository should know. Scientific caveats of the physics model are in the paper
and in [`methods/`](methods/); this file covers the repository itself: environment, provenance, figures and documents.
Nothing here changes a result.

## 1. Environment and performance

1. **Tested environment.** Only Python 3.13.3 on Windows 11 with the packages pinned in `requirements.txt` has been
   tested. Python 3.10 or newer is expected to work; Python 3.9 cannot import `bdg_eigensolver_r1_r2_LOAdapted.py`
   (PEP 604 annotations). Linux and macOS, and installing `pfapack` there, have not been tested.
2. **Console encoding.** The scripts print Greek letters. On a Windows console that is not UTF-8, run them with
   `python -X utf8` or set `PYTHONUTF8=1`.
3. **`src/bdg_parallel.py` has a hard-coded worker and memory budget.** `default_workers()` returns
   `min(physical cores, 6, free RAM / 0.5 GB)`. The numbers were tuned on one machine (Intel i7-13620H, 6 performance
   plus 4 efficiency cores, 16 GB) and each worker's BLAS is pinned to one thread. On other hardware this can leave cores
   idle or, on a machine with little free memory, run fewer workers than expected. Pass `n_workers` explicitly (or set
   the module's `N_WORKERS`) to override it. The module docstrings are out of date: they speak of "R7/R8 sweep
   drivers", "a 10-core machine" and "10 workers x 1 thread", and say the `n_workers` default is the number of physical
   cores, whereas the code caps it at 6 and the R9 drivers use it too. This affects run time and memory only. Serial and
   parallel runs agree bit for bit at equal BLAS thread count; pinning BLAS to one thread can shift values at about 1e-13
   compared with an earlier multi-threaded run.
4. **Scripts overwrite the committed results.** Every module writes into `results/<topic>/` of its own checkout. See
   [`REPRODUCING.md`](REPRODUCING.md) for how to run without touching the committed data.
5. **Long runs.** Several scripts take from 15 minutes to over 4 hours (table in `REPRODUCING.md`). For
   `bdg_r5_anderson_mixing.py` and `bdg_r9_multisubband_bracket_scan.py` no run time is recorded anywhere in the
   repository. The run-time estimates written before some of the R9 runs undershot the real wall-clock time; the
   methods notes say so where it happened.
6. **Timing-dependent files differ between machines.** `performance_baseline.csv`, `r3_lo_performance.csv`,
   `r3_lo_sweep_tracking.csv` and the figures that plot them contain wall-clock measurements. All other columns reproduce
   exactly.
7. **`src/bdg_eigensolver_r1_r2_LOAdapted.py` starts with a stray line `'PYEOF'`.** It is a harmless string
   expression, but it means the module's `__doc__` is `'PYEOF'` and the descriptive docstring that follows is not
   attached to the module. It was left as it was because removing it is not a path, name or comment change.
8. **`results/finite_temperature/r9_full_run.log` is an archival run log** and contains the absolute Windows paths
   of the machine it ran on, including a user name and the earlier name of this repository. It is kept byte for byte, since
   editing it would falsify the record.

## 2. Provenance gaps

Static search of `src/` and `scripts/` finds no code that writes the following files. They are kept as the committed
record, and the methods notes describe how they were made.

- **Data:** `results/disorder/r9_disorder_N_independence.json`, `results/disorder/r9_disorder_survivorship_check.json`
  and `results/disorder/r9_disorder_data.SUPERSEDED_task4_3.json`. The last is the superseded first-round disorder
  result, deliberately preserved.
- **Figures:**
  - `results/self_consistency_bifurcation/r5_figure10_shift_vs_g.png` and `r5_figure11_headline_validation.png`;
  - `results/multisubband/r9_multisubband_summary.png`;
  - `results/disorder/r9_disorder_summary.SUPERSEDED_task4_3.png`;
  - `paper/figures/fig3_null_validation.png` (Fig. 4 of the typeset paper, a re-layout of `r5_figure11`) and
    `paper/figures/fig6_multisubband_summary.png` (Fig. 8, a re-layout of `r9_multisubband_summary.png`).
- **R4 saves no data file**, only figures. Figures 4 and 6 of that module take more than 15 minutes each to recompute.

The Data availability paragraph of the paper says the repository contains "every figure-generating script". That is not
yet true for the figures listed above; either the scripts that made them have to be added or the sentence softened.

`docs/paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md` has two inaccuracies and one statement that no longer holds.
They were left untouched, because the log is a historical record and only its path strings were updated:

- its table of paper figures numbers the files in the order they were made (`fig1_...` to `fig6_...`), not as the
  figures appear in the typeset paper, and it names `bdg_r5_physically_meaningful_shift.py` as the generator of
  `fig3_null_validation.png` and `bdg_r9_multisubband_sweep.py` as the generator of `fig6_multisubband_summary.png`;
  neither script writes those images;
- the paragraph headed "Fig. 3's clipped embedded title" is about `paper/figures/fig4_bounded_pole.png`; its "Fig. 3"
  is neither the file number nor the number in the current typeset paper, where that figure is Fig. 5. The ad hoc script
  it mentions, `regen_fig3.py`, was never part of the repository and could not be recovered.
  `scripts/redraw_paper_fig4.py` replaces it; see section 3;
- its list of "genuine placeholders" says the repository URL is not in the Data availability statement. It is now: the
  paragraph gives the repository name and URL (no DOI).

## 3. Figures

During the restructure, the phase, track and competition labels (for example "R5 Figure 12", "R9 Track 1") were removed
from text drawn inside figures, and the stale label "Woodbury" in four figures was replaced by "Subspace", which is what
the code already said. The affected figures were redrawn with the **same plotting code and the same saved data**.
Before each was accepted, a fingerprint of everything plotted (line coordinates, images, collections, patches) and the
axis limits and scales were compared with the original, and every changed pixel was checked to lie inside a changed
text or legend box. Details that remain:

- **Layout shifts.** Text length enters `tight_layout` and `bbox_inches="tight"`, so shorter or longer titles move
  things slightly without changing the data. In `r3_lo_figure1_accuracy.png` and `r3_lo_figure3_sweep_tracking.png`
  the canvas is wider than before (2414 to 2476 and 2923 to 3574 pixels) because the new titles are longer; in
  `r5_figure13_efficiency_benchmark.png` the text box is narrower; in `r6_figure3_phase_diagram.png` the left panel is
  about 7 % narrower because the legend now cites `docs/methods/...` instead of the old `src/...` path.
- **Antialiasing.** Rerunning the original code for `r6_figure2_gap_vs_Bz.png` gives 536 pixels that differ from the
  committed image by at most 2/255 in one channel, which is antialiasing noise.
- **Figures that still carry stage labels.** `r5_figure10_shift_vs_g.png`, `r5_figure11_headline_validation.png` and
  `results/multisubband/r9_multisubband_summary.png` have no producing script (section 2), so they could not be
  redrawn without recomputing or reconstructing code. Their titles still read "R5 Figure 10", "R5 Figure 11" and "R9
  Track 4 ...", and the last also names a methods-note file. The glossary in
  [`RENAME_MAP.md`](RENAME_MAP.md) explains the codes. `paper/figures/fig3_null_validation.png` and
  `fig6_multisubband_summary.png` are re-layouts of two of them and do not carry the labels.
- **`paper/figures/fig4_bounded_pole.png`** was redrawn by `scripts/redraw_paper_fig4.py` from `r9_pole_data.json`. With
  the original text strings the same code reproduces the previously committed image pixel for pixel; the new image
  differs from it only in the super-title ("R9 Track 1" removed) and in one legend entry ("R5 g_crit band" became
  "g_crit band").
- **`results/em_feedback_detector_comparison/r4_figure4_phase_shift_preview.png` was not redrawn.** The R4 module keeps
  no data file, so redrawing it means recomputing a 121-point grid of open-wire critical-field searches. A comparison run
  of the original and the restructured code had not finished after about 45 minutes and was abandoned. The committed
  image therefore still has the title "R4 Figure 4: ..." and, in its caveat line, the old path
  `src/r4_pfaffian_boundary_methods.md`, which is now `docs/methods/r4_pfaffian_boundary_methods.md`. The plotting code
  in `src/bdg_r4_em_coupling.py` has the new text, so rerunning that module produces the corrected figure.
- **`r4_figure6_pfaffian_vs_threshold.png` was recomputed** with the unchanged function, once with the original code and
  once with the restructured code, in scratch directories. The returned arrays, the plotted-data fingerprint and a render
  without text were identical, and the original code reproduced the committed image byte for byte. The new image differs
  only in its title ("R4 Figure 6: " removed). `paper/figures/fig2_pfaffian_vs_threshold.png` is a byte copy of it.
- **Redrawn from saved data, not rerun.** The other redrawn figures that need hours to recompute were drawn from the
  saved result files.

## 4. The paper

- `main.tex` still has `\affiliation{Author affiliation to be completed}`.
- `\date{\today}` means every build differs from the committed `paper/main.pdf` in the date line.
- The Data availability paragraph gives the repository name and URL but no DOI; add the DOI after the first release has
  been archived.
- The build reports 2 overfull boxes, 15 underfull boxes and 2 "Deferred float stuck during \clearpage processing"
  messages, and no undefined references. A build of the paper as it stood before the restructure reports the same
  overfull boxes and deferred floats and 11 underfull boxes; the four additional underfull boxes are line breaks in the
  Data availability paragraph, which now contains the repository URL. No layout change was made to remove them.
- The figures that carry no phase or track labels now were redrawn; the paper text itself was changed only in the Data
  availability paragraph and in the data-folder paths quoted in six figure captions.

## 5. Documents

- **Stage codes.** The methods notes, the result summaries and the paper audit log refer to the research stages as R1
  to R9 and to "Track" and "Task" numbers. They were not rewritten, because they record the order in which results were
  obtained and corrected. [`RENAME_MAP.md`](RENAME_MAP.md) maps each code to its directory.
- **Logbook.** `docs/logbook/` contains the research logbook as a Word document. It was moved without changes and
  contains notes from the project's earlier, competition-oriented phase and two mentions of a previous LaTeX hosting
  service. Whether it belongs in a public archive is the author's decision.
- **Superseded results are kept.** Files and sections marked `SUPERSEDED`, and the correction boxes in the methods notes,
  record earlier numbers that were later found to be wrong or inflated. They are intentionally not deleted.

## 6. Notes on the restructure

Besides moving files and updating paths, the restructure edited docstrings and comments in the modules to remove
competition framing. That included deleting the line "Status: IN PROGRESS" from the headers of
`bdg_eigensolver_r1_r2_LOAdapted.py` and `bdg_r4_em_coupling.py`, and an unfilled advisor placeholder from the header of
the first. No computation, parameter, tolerance, seed or numerical output was changed.
[`RENAME_MAP.md`](RENAME_MAP.md) lists every moved file.
