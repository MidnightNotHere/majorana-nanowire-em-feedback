# Rename map

Every file tracked at the tag `pre-cleanup` (`1e0a2ab`), where it is now, and what, if anything, changed in its content. Moves were made with `git mv`, so `git log --follow <new path>` reaches the pre-restructure history. Commit hashes quoted in the methods notes (and in the logbook, which is no longer in the working tree) refer to that history and are still valid.

## Phase and track codes

The research was carried out in numbered stages (R1-R9). The directory names below describe the content instead. File names keep their `rN_` prefix, because the `src/` modules import one another by name and the methods notes and result summaries cite them by that name.

| Code | Content | Directory under `results/` |
|---|---|---|
| R1, R2 | Baseline BdG eigensolver (validated against Kitaev 2001) and the Pfaffian Z2 indicator | `baseline_bdg_pfaffian/` |
| R3 | PH-symmetric rank-2 subspace update of the BdG eigenbasis | `subspace_update/` |
| R4 | Biot-Savart EM coupling layer; spectral-threshold vs Pfaffian boundary detector | `em_feedback_detector_comparison/` |
| R5, R5-revised | Self-consistent ring self-field, Anderson mixing, physically meaningful shift, N-scaling null result | `self_consistency_bifurcation/` |
| R6 | Lutchyn-Oreg 2D (B, mu) phase diagram with a TRIM Pfaffian invariant | `lutchyn_oreg_phase_diagram/` |
| R7 | Wire length / circuit inductance L = N a as a third axis | `inductance_scan/` |
| R8 (A, B, C) | Reactive impedance Z_r: A analytic, B Kubo, C self-consistent spin-orbit | `reactive_impedance/` |
| R9 Track 1 | Characterisation of the apparent Z_r pole near g_crit | `reactive_impedance/pole_characterization/` |
| R9 Track 2 | Finite temperature | `finite_temperature/` |
| R9 Track 3 | Disorder | `disorder/` |
| R9 Track 4 | Multi-subband wire (M = 1 to 6) | `multisubband/` |

## File-by-file map (138 files)

| # | Old path (at `pre-cleanup`) | New path | Change |
|--:|---|---|---|
| 1 | `.gitignore` | `.gitignore` | rules rewritten |
| 2 | `ISEF R4 data folder/r4_figure1_peierls_validation.png` | `results/em_feedback_detector_comparison/r4_figure1_peierls_validation.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 3 | `ISEF R4 data folder/r4_figure2_biot_savart.png` | `results/em_feedback_detector_comparison/r4_figure2_biot_savart.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 4 | `ISEF R4 data folder/r4_figure3_current_density.png` | `results/em_feedback_detector_comparison/r4_figure3_current_density.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 5 | `ISEF R4 data folder/r4_figure4_phase_shift_preview.png` | `results/em_feedback_detector_comparison/r4_figure4_phase_shift_preview.png` | moved, byte-identical |
| 6 | `ISEF R4 data folder/r4_figure5_flux_density.png` | `results/em_feedback_detector_comparison/r4_figure5_flux_density.png` | moved, byte-identical |
| 7 | `ISEF R4 data folder/r4_figure6_pfaffian_vs_threshold.png` | `results/em_feedback_detector_comparison/r4_figure6_pfaffian_vs_threshold.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 8 | `ISEF R5 data folder/r5_anderson_results.json` | `results/self_consistency_bifurcation/r5_anderson_results.json` | moved, byte-identical |
| 9 | `ISEF R5 data folder/r5_figure10_shift_vs_g.png` | `results/self_consistency_bifurcation/r5_figure10_shift_vs_g.png` | moved; figure redrawn: in-figure text only, plotted data identical; drawn by `scripts/redraw_r5_figure10.py` (reconstructed plotting code, `KNOWN_ISSUES.md` section 3), stage codes removed from the in-figure text |
| 10 | `ISEF R5 data folder/r5_figure11_headline_validation.png` | `results/self_consistency_bifurcation/r5_figure11_headline_validation.png` | moved; figure redrawn: in-figure text only, plotted data identical; drawn by `scripts/redraw_r5_figure11.py` (reconstructed plotting code, `KNOWN_ISSUES.md` section 3), stage codes removed from the in-figure text |
| 11 | `ISEF R5 data folder/r5_figure12_anderson_above_gcrit.png` | `results/self_consistency_bifurcation/r5_figure12_anderson_above_gcrit.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 12 | `ISEF R5 data folder/r5_figure13_efficiency_benchmark.png` | `results/self_consistency_bifurcation/r5_figure13_efficiency_benchmark.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 13 | `ISEF R5 data folder/r5_figure7_scf_convergence.png` | `results/self_consistency_bifurcation/r5_figure7_scf_convergence.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 14 | `ISEF R5 data folder/r5_figure8_g_ramp_hysteresis.png` | `results/self_consistency_bifurcation/r5_figure8_g_ramp_hysteresis.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 15 | `ISEF R5 data folder/r5_figure9_shift_robustness.png` | `results/self_consistency_bifurcation/r5_figure9_shift_robustness.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 16 | `ISEF R5 data folder/r5_physical_shift_results.json` | `results/self_consistency_bifurcation/r5_physical_shift_results.json` | moved, byte-identical |
| 17 | `ISEF R5 data folder/r5_results.json` | `results/self_consistency_bifurcation/r5_results.json` | moved, byte-identical |
| 18 | `ISEF R5 data folder/r5_shift_comparison.csv` | `results/self_consistency_bifurcation/r5_shift_comparison.csv` | moved, byte-identical |
| 19 | `ISEF R7 data folder/R7_SUMMARY.md` | `results/inductance_scan/R7_SUMMARY.md` | moved; path references updated |
| 20 | `ISEF R7 data folder/r7_figure1_shift_vs_BmuL.png` | `results/inductance_scan/r7_figure1_shift_vs_BmuL.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 21 | `ISEF R7 data folder/r7_figure2_regime_and_validation.png` | `results/inductance_scan/r7_figure2_regime_and_validation.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 22 | `ISEF R7 data folder/r7_points_cache.json` | `results/inductance_scan/r7_points_cache.json` | moved, byte-identical |
| 23 | `ISEF R7 data folder/r7_results.json` | `results/inductance_scan/r7_results.json` | moved, byte-identical |
| 24 | `ISEF R8 data folder/R8B_FOLLOWUPS.md` | `results/reactive_impedance/R8B_FOLLOWUPS.md` | moved, byte-identical |
| 25 | `ISEF R8 data folder/R8B_KUBO_SUMMARY.md` | `results/reactive_impedance/R8B_KUBO_SUMMARY.md` | moved, byte-identical |
| 26 | `ISEF R8 data folder/R8C_SO_SUMMARY.md` | `results/reactive_impedance/R8C_SO_SUMMARY.md` | moved, byte-identical |
| 27 | `ISEF R8 data folder/R8_TRACKA_SUMMARY.md` | `results/reactive_impedance/R8_TRACKA_SUMMARY.md` | moved, byte-identical |
| 28 | `ISEF R8 data folder/R9_TRACK1_SUMMARY.md` | `results/reactive_impedance/pole_characterization/R9_TRACK1_SUMMARY.md` | moved, byte-identical |
| 29 | `ISEF R8 data folder/r8_figure1_Zr_vs_gap.png` | `results/reactive_impedance/r8_figure1_Zr_vs_gap.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 30 | `ISEF R8 data folder/r8_figure2_nscaling_grid.png` | `results/reactive_impedance/r8_figure2_nscaling_grid.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 31 | `ISEF R8 data folder/r8_trackA_data.json` | `results/reactive_impedance/r8_trackA_data.json` | moved, byte-identical |
| 32 | `ISEF R8 data folder/r8b_figure_followups.png` | `results/reactive_impedance/r8b_figure_followups.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 33 | `ISEF R8 data folder/r8b_figure_kubo_impedance.png` | `results/reactive_impedance/r8b_figure_kubo_impedance.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 34 | `ISEF R8 data folder/r8b_followups.json` | `results/reactive_impedance/r8b_followups.json` | moved, byte-identical |
| 35 | `ISEF R8 data folder/r8b_kubo_data.json` | `results/reactive_impedance/r8b_kubo_data.json` | moved, byte-identical |
| 36 | `ISEF R8 data folder/r8c_figure_so_selfconsistent.png` | `results/reactive_impedance/r8c_figure_so_selfconsistent.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 37 | `ISEF R8 data folder/r8c_so_data.json` | `results/reactive_impedance/r8c_so_data.json` | moved, byte-identical |
| 38 | `ISEF R8 data folder/r9_figure1_bounded_pole.png` | `results/reactive_impedance/pole_characterization/r9_figure1_bounded_pole.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 39 | `ISEF R8 data folder/r9_figure2_kubo_trackC.png` | `results/reactive_impedance/pole_characterization/r9_figure2_kubo_trackC.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 40 | `ISEF R8 data folder/r9_pole_data.json` | `results/reactive_impedance/pole_characterization/r9_pole_data.json` | moved, byte-identical |
| 41 | `ISEF R9 data folder/fig7_multisubband_M1to6_extension.png` | `results/multisubband/fig7_multisubband_M1to6_extension.png` | moved; pre-existing uncommitted work, snapshotted in 57c9107 (panel layout 1x3 changed to 3x1 by the author; not part of the restructure) |
| 42 | `ISEF R9 data folder/r9_disorder_N_independence.json` | `results/disorder/r9_disorder_N_independence.json` | moved, byte-identical |
| 43 | `ISEF R9 data folder/r9_disorder_data.SUPERSEDED_task4_3.json` | `results/disorder/r9_disorder_data.SUPERSEDED_task4_3.json` | moved, byte-identical |
| 44 | `ISEF R9 data folder/r9_disorder_data.json` | `results/disorder/r9_disorder_data.json` | moved, byte-identical |
| 45 | `ISEF R9 data folder/r9_disorder_summary.SUPERSEDED_task4_3.png` | `results/disorder/r9_disorder_summary.SUPERSEDED_task4_3.png` | moved, byte-identical |
| 46 | `ISEF R9 data folder/r9_disorder_summary.png` | `results/disorder/r9_disorder_summary.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 47 | `ISEF R9 data folder/r9_disorder_survivorship_check.json` | `results/disorder/r9_disorder_survivorship_check.json` | moved, byte-identical |
| 48 | `ISEF R9 data folder/r9_finiteT_data.json` | `results/finite_temperature/r9_finiteT_data.json` | moved, byte-identical |
| 49 | `ISEF R9 data folder/r9_finiteT_summary.png` | `results/finite_temperature/r9_finiteT_summary.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 50 | `ISEF R9 data folder/r9_full_run.log` | `results/finite_temperature/r9_full_run.log` | moved; edited; absolute machine paths in two lines replaced by `<REPO_ROOT>` placeholders |
| 51 | `ISEF R9 data folder/r9_multisubband_bc_M2_N120.json` | `results/multisubband/r9_multisubband_bc_M2_N120.json` | moved, byte-identical |
| 52 | `ISEF R9 data folder/r9_multisubband_bc_M2_N60.json` | `results/multisubband/r9_multisubband_bc_M2_N60.json` | moved, byte-identical |
| 53 | `ISEF R9 data folder/r9_multisubband_bc_M2_N80.json` | `results/multisubband/r9_multisubband_bc_M2_N80.json` | moved, byte-identical |
| 54 | `ISEF R9 data folder/r9_multisubband_bc_M2_summary.json` | `results/multisubband/r9_multisubband_bc_M2_summary.json` | moved, byte-identical |
| 55 | `ISEF R9 data folder/r9_multisubband_bracket_scan.json` | `results/multisubband/r9_multisubband_bracket_scan.json` | moved, byte-identical |
| 56 | `ISEF R9 data folder/r9_multisubband_preregistration.json` | `results/multisubband/r9_multisubband_preregistration.json` | moved, byte-identical |
| 57 | `ISEF R9 data folder/r9_multisubband_preregistration_M5_M6.json` | `results/multisubband/r9_multisubband_preregistration_M5_M6.json` | moved, byte-identical |
| 58 | `ISEF R9 data folder/r9_multisubband_protocol_a_M6.json` | `results/multisubband/r9_multisubband_protocol_a_M6.json` | moved, byte-identical |
| 59 | `ISEF R9 data folder/r9_multisubband_summary.png` | `results/multisubband/r9_multisubband_summary.png` | moved; figure redrawn: in-figure text only, plotted data identical; drawn by `scripts/redraw_r9_multisubband_summary.py` (reconstructed plotting code, `KNOWN_ISSUES.md` section 3), stage codes removed from the in-figure text |
| 60 | `ISEF R9 data folder/r9_multisubband_sweep_M5_M6.json` | `results/multisubband/r9_multisubband_sweep_M5_M6.json` | moved, byte-identical |
| 61 | `ISEF R9 data folder/r9_multisubband_sweep_data.json` | `results/multisubband/r9_multisubband_sweep_data.json` | moved, byte-identical |
| 62 | `Physics Research Logbook - Majorana Nanowire Inductive Feedback Project.docx` | (removed from the working tree) | removed from the working tree on 2026-10-25; still in history (tag `pre-cleanup`), see `KNOWN_ISSUES.md` section 5 |
| 63 | `R1_R2_data_files/figure1_phase_sweep.png` | `results/baseline_bdg_pfaffian/figure1_phase_sweep.png` | moved, byte-identical |
| 64 | `R1_R2_data_files/figure2_finite_size.png` | `results/baseline_bdg_pfaffian/figure2_finite_size.png` | moved, byte-identical |
| 65 | `R1_R2_data_files/figure3_performance.png` | `results/baseline_bdg_pfaffian/figure3_performance.png` | moved, byte-identical |
| 66 | `R1_R2_data_files/figure4_phase_diagram.png` | `results/baseline_bdg_pfaffian/figure4_phase_diagram.png` | moved, byte-identical |
| 67 | `R1_R2_data_files/figure5_pfaffian_tracking.png` | `results/baseline_bdg_pfaffian/figure5_pfaffian_tracking.png` | moved, byte-identical |
| 68 | `R1_R2_data_files/figure6_lo_phase_diagram.png` | `results/baseline_bdg_pfaffian/figure6_lo_phase_diagram.png` | moved, byte-identical |
| 69 | `R1_R2_data_files/figure7_lo_finite_size.png` | `results/baseline_bdg_pfaffian/figure7_lo_finite_size.png` | moved, byte-identical |
| 70 | `R1_R2_data_files/finite_size_scaling.csv` | `results/baseline_bdg_pfaffian/finite_size_scaling.csv` | moved, byte-identical |
| 71 | `R1_R2_data_files/lo_finite_size_scaling.csv` | `results/baseline_bdg_pfaffian/lo_finite_size_scaling.csv` | moved, byte-identical |
| 72 | `R1_R2_data_files/lo_gap_map.npy` | `results/baseline_bdg_pfaffian/lo_gap_map.npy` | moved, byte-identical |
| 73 | `R1_R2_data_files/lo_nu_map.npy` | `results/baseline_bdg_pfaffian/lo_nu_map.npy` | moved, byte-identical |
| 74 | `R1_R2_data_files/performance_baseline.csv` | `results/baseline_bdg_pfaffian/performance_baseline.csv` | moved, byte-identical |
| 75 | `R1_R2_data_files/pfaffian_tracking.csv` | `results/baseline_bdg_pfaffian/pfaffian_tracking.csv` | moved, byte-identical |
| 76 | `R1_R2_data_files/phase_sweep.csv` | `results/baseline_bdg_pfaffian/phase_sweep.csv` | moved, byte-identical |
| 77 | `R3_data_files/r3_lo_accuracy.csv` | `results/subspace_update/r3_lo_accuracy.csv` | moved, byte-identical |
| 78 | `R3_data_files/r3_lo_figure1_accuracy.png` | `results/subspace_update/r3_lo_figure1_accuracy.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 79 | `R3_data_files/r3_lo_figure2_speedup.png` | `results/subspace_update/r3_lo_figure2_speedup.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 80 | `R3_data_files/r3_lo_figure3_sweep_tracking.png` | `results/subspace_update/r3_lo_figure3_sweep_tracking.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 81 | `R3_data_files/r3_lo_performance.csv` | `results/subspace_update/r3_lo_performance.csv` | moved, byte-identical |
| 82 | `R3_data_files/r3_lo_sweep_tracking.csv` | `results/subspace_update/r3_lo_sweep_tracking.csv` | moved, byte-identical |
| 83 | `R6_data_files/r6_figure1_orv_dispersion.png` | `results/lutchyn_oreg_phase_diagram/r6_figure1_orv_dispersion.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 84 | `R6_data_files/r6_figure2_gap_vs_Bz.png` | `results/lutchyn_oreg_phase_diagram/r6_figure2_gap_vs_Bz.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 85 | `R6_data_files/r6_figure3_phase_diagram.png` | `results/lutchyn_oreg_phase_diagram/r6_figure3_phase_diagram.png` | moved; figure redrawn: in-figure text only, plotted data identical |
| 86 | `R6_data_files/r6_gap_map.npy` | `results/lutchyn_oreg_phase_diagram/r6_gap_map.npy` | moved, byte-identical |
| 87 | `R6_data_files/r6_gap_vs_Bz.csv` | `results/lutchyn_oreg_phase_diagram/r6_gap_vs_Bz.csv` | moved, byte-identical |
| 88 | `R6_data_files/r6_nu_map.npy` | `results/lutchyn_oreg_phase_diagram/r6_nu_map.npy` | moved, byte-identical |
| 89 | `README.md` | `README.md` | rewritten |
| 90 | `paper/PAPER_SUMMARY_AND_VERIFICATION_LOG.md` | `docs/paper_audit/PAPER_SUMMARY_AND_VERIFICATION_LOG.md` | moved; path references updated; three dated correction notes added (2026-10-25), original text unchanged |
| 91 | `paper/figures/fig1_phase_diagram.png` | `paper/figures/fig1_phase_diagram.png` | figure redrawn: byte copy of the regenerated result figure of the same content |
| 92 | `paper/figures/fig2_pfaffian_vs_threshold.png` | `paper/figures/fig2_pfaffian_vs_threshold.png` | figure redrawn: byte copy of the regenerated result figure of the same content |
| 93 | `paper/figures/fig3_null_validation.png` | `paper/figures/fig3_null_validation.png` | unchanged |
| 94 | `paper/figures/fig4_bounded_pole.png` | `paper/figures/fig4_bounded_pole.png` | figure redrawn by `scripts/redraw_paper_fig4.py` from the saved pole data: super-title and one legend label lost their stage prefix, plotted data identical |
| 95 | `paper/figures/fig5_disorder_summary.png` | `paper/figures/fig5_disorder_summary.png` | figure redrawn: byte copy of the regenerated result figure of the same content |
| 96 | `paper/figures/fig6_multisubband_summary.png` | `paper/figures/fig6_multisubband_summary.png` | unchanged |
| 97 | `paper/figures/fig7_multisubband_M1to6_extension.png` | `paper/figures/fig7_multisubband_M1to6_extension.png` | pre-existing uncommitted work, snapshotted in 57c9107 (panel layout 1x3 changed to 3x1 by the author; not part of the restructure) |
| 98 | `paper/main.aux` | (removed from index) | removed from the index: LaTeX build artefact, now git-ignored (recreated by `latexmk`) |
| 99 | `paper/main.bbl` | `paper/main.bbl` | unchanged |
| 100 | `paper/main.blg` | (removed from index) | removed from the index: LaTeX build artefact, now git-ignored (recreated by `latexmk`) |
| 101 | `paper/main.fdb_latexmk` | (removed from index) | removed from the index: LaTeX build artefact, now git-ignored (recreated by `latexmk`) |
| 102 | `paper/main.fls` | (removed from index) | removed from the index: LaTeX build artefact, now git-ignored (recreated by `latexmk`) |
| 103 | `paper/main.log` | (removed from index) | removed from the index: LaTeX build artefact, now git-ignored (recreated by `latexmk`) |
| 104 | `paper/main.pdf` | `paper/main.pdf` | pre-existing uncommitted work, snapshotted in 57c9107; rebuilt |
| 105 | `paper/main.tex` | `paper/main.tex` | pre-existing uncommitted work, snapshotted in 57c9107; figure-caption data paths and Data availability text updated |
| 106 | `paper/mainNotes.bib` | (removed from index) | removed from the index: LaTeX build artefact, now git-ignored (recreated by `latexmk`) |
| 107 | `paper/references.bib` | `paper/references.bib` | unchanged |
| 108 | `requirements.txt` | `requirements.txt` | re-encoded UTF-16 to UTF-8, same 16 pins |
| 109 | `src/bdg_eigensolver_r1_r2_LOAdapted.py` | `src/bdg_eigensolver_r1_r2_LOAdapted.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed); stray first line `'PYEOF'` removed |
| 110 | `src/bdg_eigensolver_r3_lo_adapted.py` | `src/bdg_eigensolver_r3_lo_adapted.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 111 | `src/bdg_parallel.py` | `src/bdg_parallel.py` | unchanged |
| 112 | `src/bdg_r4_em_coupling.py` | `src/bdg_r4_em_coupling.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 113 | `src/bdg_r5_anderson_mixing.py` | `src/bdg_r5_anderson_mixing.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 114 | `src/bdg_r5_physically_meaningful_shift.py` | `src/bdg_r5_physically_meaningful_shift.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 115 | `src/bdg_r5_selfconsistent_ring.py` | `src/bdg_r5_selfconsistent_ring.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 116 | `src/bdg_r6_lo_phase_diagram.py` | `src/bdg_r6_lo_phase_diagram.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 117 | `src/bdg_r7_inductance_sweep.py` | `src/bdg_r7_inductance_sweep.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 118 | `src/bdg_r8_reactive_impedance.py` | `src/bdg_r8_reactive_impedance.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 119 | `src/bdg_r8b_kubo_impedance.py` | `src/bdg_r8b_kubo_impedance.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 120 | `src/bdg_r8c_so_selfconsistent.py` | `src/bdg_r8c_so_selfconsistent.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 121 | `src/bdg_r9_disorder.py` | `src/bdg_r9_disorder.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 122 | `src/bdg_r9_finite_temperature.py` | `src/bdg_r9_finite_temperature.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 123 | `src/bdg_r9_multisubband.py` | `src/bdg_r9_multisubband.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 124 | `src/bdg_r9_multisubband_bc.py` | `src/bdg_r9_multisubband_bc.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 125 | `src/bdg_r9_multisubband_bracket_scan.py` | `src/bdg_r9_multisubband_bracket_scan.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 126 | `src/bdg_r9_multisubband_plot_M6.py` | `src/bdg_r9_multisubband_plot_M6.py` | pre-existing uncommitted work, snapshotted in 57c9107; edited: output paths, docstring/comment text, in-figure text literals (no computation changed); output folders created before saving |
| 127 | `src/bdg_r9_multisubband_sweep.py` | `src/bdg_r9_multisubband_sweep.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 128 | `src/bdg_r9_multisubband_sweep_M56.py` | `src/bdg_r9_multisubband_sweep_M56.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 129 | `src/bdg_r9_pole_characterization.py` | `src/bdg_r9_pole_characterization.py` | edited: output paths, docstring/comment text, in-figure text literals (no computation changed) |
| 130 | `src/r2_o1_pfaffian_trim_lo_correction.md` | `docs/methods/r2_o1_pfaffian_trim_lo_correction.md` | moved; path references updated; `Project:` line replaced by the paper title |
| 131 | `src/r4_pfaffian_boundary_methods.md` | `docs/methods/r4_pfaffian_boundary_methods.md` | moved; path references updated; `Project:` line replaced by the paper title |
| 132 | `src/r6_lo_phase_diagram_methods.md` | `docs/methods/r6_lo_phase_diagram_methods.md` | moved; path references updated; `Project:` line replaced by the paper title |
| 133 | `src/r8_reactive_impedance_methods.md` | `docs/methods/r8_reactive_impedance_methods.md` | moved; path references updated |
| 134 | `src/r9_disorder_methods.md` | `docs/methods/r9_disorder_methods.md` | moved; path references updated |
| 135 | `src/r9_finite_temperature_methods.md` | `docs/methods/r9_finite_temperature_methods.md` | moved; path references updated |
| 136 | `src/r9_multisubband_M56_run.log` | `results/multisubband/r9_multisubband_M56_run.log` | moved, byte-identical |
| 137 | `src/r9_multisubband_methods.md` | `docs/methods/r9_multisubband_methods.md` | moved; path references updated |
| 138 | `src/r9_multisubband_scoping.md` | `docs/methods/r9_multisubband_scoping.md` | moved; path references updated |

Not tracked and therefore not in the table: LaTeX build artefacts (`paper/main.aux`, `.log`, `.fls`, `.fdb_latexmk`, `.blg`, `.out`, `mainNotes.bib`), which are git-ignored, and the files added by the restructure itself (`LICENSE`, `CITATION.cff`, `.gitattributes`, `scripts/`, `docs/RENAME_MAP.md`, `docs/KNOWN_ISSUES.md`, `docs/REPRODUCING.md`).
