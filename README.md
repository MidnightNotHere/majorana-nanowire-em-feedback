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
  This ~0.6–0.8 % is itself now SUPERSEDED by R5. See
  `src/r4_pfaffian_boundary_methods.md`.
- ✅ R5: Self-consistent ring self-field (`src/bdg_r5_selfconsistent_ring.py`;
  R4 imported read-only, all 6 R4 checks still pass). Replaced *both* interim
  placeholders — the open-wire-seeded texture and the gauge-only wrap bond —
  with a real ring calculation: periodic ring current (incl. the wrap bond),
  a Biot–Savart kernel adapted to the closed loop and validated against the
  analytic loop inductance, and a fully self-consistent ring/Pfaffian B_c.
  **Findings (honest):** (i) plain fixed-point iteration does **not** converge
  at g₀=3×10⁵ from any seed (map slope ≈3>1, Picard-repelling) — solved by
  root-finding the scalar ring map instead; below g_crit≈1.1×10⁵ iteration
  converges and is seed-independent; (ii) the fixed point is locally unique
  (Newton–Krylov finds no non-uniform solution) though the map is multi-valued
  (fluxoid branches; principal branch used; ramp = direct-jump, path-
  independent); (iii) the shift **passes** the N-independence, tolerance and
  gauge/decomposition checks. **Result:** the self-consistent shift is
  **+0.05 % at N=80 (μ=1.2, 1.4) and → 0 in the bulk (~N⁻¹·⁵)** — a finite-size
  Aharonov–Bohm flux effect, ~15× below the interim number and of *opposite
  sign*. Progression: threshold −5.86/−9.36 % → interim −0.77/−0.61 % →
  self-consistent +0.051/+0.047 % (μ=1.2/1.4). Net: the earlier "topological
  region grows under feedback" claim does **not** survive self-consistency;
  the bulk boundary shift is consistent with zero. `g=3×10⁵` remains a
  circuit-amplification placeholder (not derived) — the one substantial
  caveat left.
- ⬜ R6: Phase diagram under EM feedback

## Requirements

pip install numpy scipy matplotlib pandas pfapack psutil

## Running

    python src/bdg_eigensolver_r1_r2_LOAdapted.py   # R1+R2
    python src/bdg_eigensolver_r3_lo_adapted.py     # R3
    python src/bdg_r4_em_coupling.py                # R4
    python src/bdg_r5_selfconsistent_ring.py        # R5 (~35 min full run)

Each program writes its figures and data outputs directly into its own
data folder (`R1_R2_data_files/`, `R3_data_files/`, `ISEF R4 data folder/`,
`ISEF R5 data folder/`).

## Results
See the per-program data folders above for all generated plots and
CSV/NPY outputs.

## Science Fair
Denver Metro Regional Science and Engineering Fair (DMRSEF)  
February 19, 2027
