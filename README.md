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
- ✅ R3: PH-symmetric rank-2 update (August 2026)
- 🔄 R4: Biot-Savart EM coupling layer
- ⬜ R5: Self-consistent fixed-point iteration
- ⬜ R6: Phase diagram under EM feedback

## Requirements

pip install numpy scipy matplotlib pandas pfapack psutil

## Running

python src/bdg_eigensolver_r1_r2.py

## Results
See `/figures` for all generated plots and `/data` for CSV outputs.

## Science Fair
Denver Metro Regional Science and Engineering Fair (DMRSEF)  
February 19, 2027
