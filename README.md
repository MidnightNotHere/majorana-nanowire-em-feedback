# ISEF-2027-bdg-eigensolver
Symmetry-aware reduced-order BdG eigensolver for topological phase transitions in superconducting materials

# ISEF 2027 — Symmetry-Aware BdG Eigensolver

**A Symmetry-Aware Reduced-Order Eigensolver for Tracking Topological 
Phase Transitions in Superconducting Nanowires Under Inductive 
Electromagnetic Feedback**

## Project Overview
This project develops three contributions toward fault-tolerant 
quantum computing:
1. A PH-symmetric rank-2 Woodbury eigensolver for BdG Hamiltonians
2. A real-time Pfaffian topological phase indicator
3. A self-consistent quasi-static EM feedback loop (Biot-Savart)

## Current Status
- ✅ R1: Baseline BdG eigensolver (validated vs Kitaev 2001)
- ✅ R2: Pfaffian topological indicator (finite-size scaling confirmed)
- ✅ R3: PH-symmetric rank-2 update
- 🔄 R4: Biot-Savart EM coupling layer — pipeline validated (6/6 checks).
  **Headline boundary-shift magnitude is PRELIMINARY:** it comes from an
  open-wire spectral-threshold detector that over-estimates the shift ~10×.
  A rigorous, threshold-free, N-independent ring/Pfaffian Z₂ detector has
  been added (`ring_boundary_Bc_pfaffian`); with the *same* self-field it
  gives ~0.6–0.8 % at μ = 1.2–1.4 (vs ~6–9 % from the threshold detector).
  Only the **sign** (topological region grows under inductive feedback) is
  trustworthy so far. See `src/r4_pfaffian_boundary_methods.md`.
- ⬜ R5: Self-consistent fixed-point iteration — **explicit deliverable
  (not just a number):** make the ring self-field self-consistent (replace
  the one-shot open-wire-seeded texture and the placeholder wrap-bond
  choice), extract the boundary shift with the ring/Pfaffian detector, and
  subject that shift to the *same* validation the detector itself passed —
  (a) N-independence across N ≈ 20–160, (b) robustness to bisection/grid
  tolerance, (c) a known-answer / gauge check (uniform flux → 0 shift), and
  (d) seed-independence of the fixed point. The goal is a *validated* final
  magnitude, not a second unvalidated one.
- ⬜ R6: Phase diagram under EM feedback

## Requirements

pip install numpy scipy matplotlib pandas pfapack psutil

## Running

    python src/bdg_eigensolver_r1_r2_LOAdapted.py   # R1+R2
    python src/bdg_eigensolver_r3_lo_adapted.py     # R3
    python src/bdg_r4_em_coupling.py                # R4

Each program writes its figures and data outputs directly into its own
data folder (`R1_R2_data_files/`, `R3_data_files/`, `ISEF R4 data folder/`).

## Results
See the per-program data folders above for all generated plots and
CSV/NPY outputs.

## Science Fair
Denver Metro Regional Science and Engineering Fair (DMRSEF)  
February 19, 2027
