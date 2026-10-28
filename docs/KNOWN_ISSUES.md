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
   [`REPRODUCING.md`](REPRODUCING.md) for how to run without touching the committed data. The one script that also writes
   outside `results/` is `src/bdg_r9_multisubband_plot_M6.py`: it saves its figure in `results/multisubband/` and copies
   it to `paper/figures/fig7_multisubband_M1to6_extension.png`. Every script that saves a file creates the output folder
   first (`os.makedirs(..., exist_ok=True)` in each module's `_out()` helper, and before the `savefig` and copy calls of
   the plot script and of the `scripts/redraw_*.py` files), so none depends on the folder already existing.
5. **Long runs.** Several scripts take from 15 minutes to over 4 hours (table in `REPRODUCING.md`). For
   `bdg_r5_anderson_mixing.py` and `bdg_r9_multisubband_bracket_scan.py` no run time is recorded anywhere in the
   repository. The run-time estimates written before some of the R9 runs undershot the real wall-clock time; the
   methods notes say so where it happened.
6. **Timing-dependent files differ between machines.** `performance_baseline.csv`, `r3_lo_performance.csv`,
   `r3_lo_sweep_tracking.csv` and the figures that plot them contain wall-clock measurements. All other columns reproduce
   exactly.
7. **`src/bdg_eigensolver_r1_r2_LOAdapted.py` used to start with a stray line `'PYEOF'`** (a leftover heredoc
   terminator). The line was deleted on 2026-10-25. It was a string expression, so no computation depended on it, but it
   had made the module's `__doc__` equal to `'PYEOF'` and left the descriptive docstring after it unattached. The
   module's `__doc__` is now that descriptive docstring; nothing else about the module changed, and it imports with
   Python 3.13.
8. **`results/finite_temperature/r9_full_run.log` is an archival run log.** It recorded the absolute Windows path of the
   machine it ran on, with a user name, in two lines (82 and 83, the "Wrote ..." lines for `r9_finiteT_data.json` and
   `r9_finiteT_summary.png`). On 2026-10-25 the part of those paths before `src\..\` was replaced by the placeholder
   `<REPO_ROOT>`; no other line and no number of the log was changed. The component that followed, the pre-restructure
   name of the output folder of the run (it began with a competition name), was replaced by `<R9 data folder>` on
   2026-10-27; the folder is now `results/finite_temperature/` (see [`RENAME_MAP.md`](RENAME_MAP.md)).

## 2. Provenance gaps

Static search of `src/` and `scripts/` finds no code that writes the following files. They are kept as the committed
record, and the methods notes describe how they were made.

- **Data:** `results/disorder/r9_disorder_N_independence.json`, `results/disorder/r9_disorder_survivorship_check.json`
  and `results/disorder/r9_disorder_data.SUPERSEDED_task4_3.json`. The last is the superseded first-round disorder
  result, deliberately preserved.
- **Figures:** `results/disorder/r9_disorder_summary.SUPERSEDED_task4_3.png`, the figure of the superseded first-round
  disorder result. It was drawn from `r9_disorder_data.SUPERSEDED_task4_3.json` by an earlier version of `make_figures()`
  in `src/bdg_r9_disorder.py`; the current function draws a different image from that file (19,023 pixels differ), so the
  image is kept as the record and not redrawn. It is not in the paper. Its panel titles read "Task 4.1", "Task 4.2" and
  "Task 4.3" (the labels of the superseded first round); they are part of the record and were left as they are.
- **R4 saves no data file**, only figures. Recomputing figure 6 (`validate_pfaffian_boundary()`) took 676 s on
  2026-10-25; figure 4 had not finished after about 45 minutes (section 3); figure 5 was not timed.

Apart from that superseded image, every committed image has a script that draws it: an `src/` module, or one of the
`scripts/redraw_*.py` files for the figures whose plotting code is not in the repository. The Data availability paragraph
of the paper says the repository contains "every figure-generating script"; that holds for every figure of the paper, and
[`REPRODUCING.md`](REPRODUCING.md), section 5, gives the script and the data file of each. The plotting code of five of
those scripts is a reconstruction (section 3), not the code that first made the images.

`docs/paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md` had two inaccuracies and one statement that no longer holds.
Each is now followed in the log by a dated correction note (2026-10-25). The original text is not edited, because the
log is a historical record; only its path strings were updated earlier:

- its table of paper figures numbers the files in the order they were made (`fig1_...` to `fig6_...`), not as the
  figures appear in the typeset paper, and it names `bdg_r5_physically_meaningful_shift.py` as the generator of
  `fig3_null_validation.png` and `bdg_r9_multisubband_sweep.py` as the generator of `fig6_multisubband_summary.png`;
  neither script writes those images, and the log's remark that `r9_multisubband_summary.png` was "used as the paper's
  Fig. 6" is inexact for the same reason;
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
- **Two stage labels removed on 2026-10-27.** The last scan of all text drawn inside the figures (every visible text
  object of every figure that can be redrawn from the committed result files or by the `scripts/` files) found the words
  "TASK-2" at the start of the verdict box of `r5_figure12_anderson_above_gcrit.png` and the stage code "R8-B" in the
  footnote of `r8_figure1_Zr_vs_gap.png` and `r8_figure2_nscaling_grid.png`. They were removed from the two strings in
  `src/bdg_r5_anderson_mixing.py` and `src/bdg_r8_reactive_impedance.py`, and the three images were redrawn from the
  committed result files. Plotted data, axis limits and a render without text are identical to the previous images, the
  canvas size is the same, and every changed pixel lies inside the box of a changed text. `R8_TRACKA_SUMMARY.md` quotes
  the same footnote with the stage code; it is a saved result and was not edited.
- **Figures redrawn by reconstructed scripts.** The plotting code of `r5_figure10_shift_vs_g.png`,
  `r5_figure11_headline_validation.png`, `results/multisubband/r9_multisubband_summary.png` and of the paper's
  `fig3_null_validation.png` and `fig6_multisubband_summary.png` is in neither the repository nor its history.
  `scripts/redraw_r5_figure10.py`, `redraw_r5_figure11.py`, `redraw_r9_multisubband_summary.py`, `redraw_paper_fig3.py`
  and `redraw_paper_fig6.py` draw them from the committed result files only (each docstring names the file and the keys
  read). They were reconstructed from the committed images: sizes, colours, line widths, margins and text positions were
  fitted until a render with the committed text reproduced the committed image. The plotted arrays were compared with the
  values stored in the result files and are equal. Fidelity of the reconstruction, against the images as committed
  before the text change: `fig3_null_validation.png` and `fig6_multisubband_summary.png` byte for byte; `r5_figure11`
  0 differing pixels; `r5_figure10` 53 of 1.6 million pixels, each by 1/255 (the annotation leader line);
  `r9_multisubband_summary` 2,952 of 1,834,560 pixels, each by 1/255 (markers of panel C). The three result images were
  then regenerated with the text changes below; the two paper images are unchanged. The text that changed is the stage
  prefix of the titles ("R5 Figure 10:", "R5 Figure 11:", "R9 Track 4 --"), the stage code in one legend entry of the
  summary ("R9-T1 single-channel baseline" became "single-channel baseline") and its methods-note file name, now given as
  `docs/methods/r9_multisubband_methods.md`. A sixth script differs in kind: `r9_disorder_summary.png` is drawn by
  `make_figures()` in `src/bdg_r9_disorder.py`, which a full disorder run (4,627 s) calls at its end;
  `scripts/redraw_r9_disorder_summary.py` calls it on the saved `r9_disorder_data.json` and reproduces the committed
  image, and the paper's copy of it, byte for byte.
- **`paper/figures/fig4_bounded_pole.png`** was redrawn by `scripts/redraw_paper_fig4.py` from `r9_pole_data.json`. With
  the original text strings the same code reproduces the previously committed image pixel for pixel; the new image
  differs from it only in the super-title ("R9 Track 1" removed) and in one legend entry ("R5 g_crit band" became
  "g_crit band").
- **Regenerate before release: `results/em_feedback_detector_comparison/r4_figure4_phase_shift_preview.png` was not
  redrawn.** The command is in [`REPRODUCING.md`](REPRODUCING.md), section 4 (it calls
  `figure4_phase_shift_preview()` of `src/bdg_r4_em_coupling.py`; it was not run). The R4 module keeps
  no data file, so redrawing it means recomputing a 121-point grid of open-wire critical-field searches. A comparison run
  of the original and the restructured code had not finished after about 45 minutes and was abandoned. The committed
  image therefore still has the title "R4 Figure 4: ..." and, in its caveat line, the old path
  `src/r4_pfaffian_boundary_methods.md`, which is now `docs/methods/r4_pfaffian_boundary_methods.md`. The plotting code
  in `src/bdg_r4_em_coupling.py` has the new text, so rerunning that module produces the corrected figure.
- **`r4_figure6_pfaffian_vs_threshold.png` was recomputed** with the unchanged function, once with the original code and
  once with the restructured code, in scratch directories. The returned arrays, the plotted-data fingerprint and a render
  without text were identical, and the original code reproduced the committed image byte for byte. The new image differs
  only in its title ("R4 Figure 6: " removed). `paper/figures/fig2_pfaffian_vs_threshold.png` is a byte copy of it. A
  further run of the committed `validate_pfaffian_boundary()` on 2026-10-25 (676 s, output redirected to a scratch
  folder) reproduced both the committed `r4_figure6` image and the paper's copy byte for byte.
- **Redrawn from saved data, not rerun.** The other redrawn figures that need hours to recompute were drawn from the
  saved result files.

## 4. The paper

- The affiliation line of `main.tex` was set to `Independent Researcher, Colorado, USA` on 2026-10-27; it had been a
  placeholder. The two comment lines under it in the source (`% NOTE: affiliation left as a placeholder ...`) are now out
  of date; they were left as they were and are not in the PDF.
- `\date{\today}` means every build differs from the committed `paper/main.pdf` in the date line.
- The Data availability paragraph gives the repository name and URL but no DOI; add the DOI after the first release has
  been archived. It was checked against the repository on 2026-10-27: every figure of the paper has a script that draws it
  ([`FIGURE_PROVENANCE.md`](FIGURE_PROVENANCE.md)), so the words "every figure-generating script" were left unchanged. Two
  details are looser than the repository: the raw results include NPY files as well as the "JSON/CSV data files" the
  paragraph names, and R4 (paper Fig. 3) saves no data file.
- The build reports 2 overfull boxes, 17 underfull boxes and 2 "Deferred float stuck during \clearpage processing"
  messages, and no undefined references. A build of the paper as it stood before the restructure reports the same
  overfull boxes and deferred floats and 11 underfull boxes; four of the additional underfull boxes are line breaks in
  the Data availability paragraph, which now contains the repository URL, and two are loose lines in the longer captions
  of Figs. 4 and 5 (badness 1168 and 1762, added on 2026-10-26). No layout change was made to remove them.
- The figures that carry no phase or track labels now were redrawn; the paper text itself was changed only in the Data
  availability paragraph, in the data-folder paths quoted in six figure captions and, on 2026-10-26, in the "Generated
  by" strings of the captions of Figs. 3, 4, 5, 6 and 8 (next item).
- **Corrected on 2026-10-26: the "Generated by" strings of five captions.** Three captions ended in a string that named
  the script that makes the data, not the script that draws the figure: Fig. 4 (`bdg_r5_physically_meaningful_shift.py`
  writes `r5_physical_shift_results.json` and no image; the figure is drawn by `scripts/redraw_paper_fig3.py`), Fig. 5
  (the paper's image is a re-titled redraw, by `scripts/redraw_paper_fig4.py`, of the figure that
  `bdg_r9_pole_characterization.py` draws) and Fig. 8 (`bdg_r9_multisubband_sweep.py` writes
  `r9_multisubband_sweep_data.json` and no image; `scripts/redraw_paper_fig6.py` draws it). The caption of Fig. 3 said
  "data in `results/em_feedback_detector_comparison/`", but R4 saves no data file and that folder holds the figure only.
  The captions of Figs. 4, 5 and 8 now name the redraw script, the saved data file it reads and the module that writes
  that file; the caption of Fig. 3 says the module saves only the figure; the caption of Fig. 6 also names
  `scripts/redraw_r9_disorder_summary.py`. The captions of Figs. 7 and 9 were already correct. No other caption or body
  text changed. The rebuilt `paper/main.pdf` differs from the previous one in these strings, in the date line (`\today`)
  and in the line and page breaks around the three longer captions.
- **Rebuilt on 2026-10-27** after the affiliation was filled in: 17 pages, 2 overfull boxes, 17 underfull boxes, 2 deferred
  floats and no undefined references or citations (`latexmk -g -pdf main.tex`). A word-level comparison of the text of the
  new and the previous PDF finds two differences: the affiliation line and the date line (26 to 27 October).

## 5. Documents

- **Stage codes.** The methods notes, the result summaries and the paper audit log refer to the research stages as R1
  to R9 and to "Track" and "Task" numbers. They were not rewritten, because they record the order in which results were
  obtained and corrected. [`RENAME_MAP.md`](RENAME_MAP.md) maps each code to its directory.
- **Logbook.** The research logbook, a Word document, was removed from the working tree of this branch on 2026-10-25,
  because it contains notes from the project's earlier, competition-oriented phase and two mentions of a previous LaTeX
  hosting service. It is not deleted from history: it is at the root of the tag `pre-cleanup` under the name
  `Physics Research Logbook - Majorana Nanowire Inductive Feedback Project.docx`, and in the commits that moved it to
  `docs/logbook/`. To restore a copy: `git show "pre-cleanup:Physics Research Logbook - Majorana Nanowire Inductive
  Feedback Project.docx" > logbook.docx`. The commit hashes that the methods notes quote still refer to that history.
  **Logbook moved out of the archive tree; still present in git history.** On 2026-10-27 a copy was placed outside the
  repository, in the folder `../majorana-private-archive/` (next to the repository folder), under the same file name; it
  is byte-identical to the file at the root of `pre-cleanup` (990,642 bytes, SHA-256
  `9a06b4d5a5317bb56b1bcef31c70e6653878c5750953f4610c8e49d9c7d25a8d`). No history was rewritten. Pull request #1 had
  already merged the restructure branch into `main` before the logbook was removed from the branch, so the file was in
  the tree of `main` (as `docs/logbook/...docx`) from that merge until this branch is merged, and it stays in the public
  history of the repository.
- **Superseded results are kept.** Files and sections marked `SUPERSEDED`, and the correction boxes in the methods notes,
  record earlier numbers that were later found to be wrong or inflated. They are intentionally not deleted.

## 6. Notes on the restructure

Besides moving files and updating paths, the restructure edited docstrings and comments in the modules to remove
competition framing. That included deleting the line "Status: IN PROGRESS" from the headers of
`bdg_eigensolver_r1_r2_LOAdapted.py` and `bdg_r4_em_coupling.py`, and an unfilled advisor placeholder from the header of
the first. No computation, parameter, tolerance, seed or numerical output was changed.
[`RENAME_MAP.md`](RENAME_MAP.md) lists every moved file.

A second pass (2026-10-25) added the six `scripts/redraw_*.py` files described in section 3, regenerated three result
figures from them, removed the logbook, replaced the machine paths in `r9_full_run.log` (section 1, item 8), gave the
`Project:` line of three methods notes (R2-O1, R4, R6) the title of the paper in place of a competition name, deleted
the stray first line of `bdg_eigensolver_r1_r2_LOAdapted.py` (section 1, item 7), made the plot script and the redraw
scripts create their output folder, and annotated the paper verification log (section 2). The paper text and
`paper/main.pdf` were not changed in that pass. A third pass (2026-10-26) corrected the "Generated by" strings of the
captions of Figs. 3, 4, 5, 6 and 8 and rebuilt `paper/main.pdf` (section 4).
