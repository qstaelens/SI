## Combining matrix product states and mean-field theory to capture magnetic order in quasi-1D cuprates

## Overview

This repository contains the VASP input files used in:

"Combining matrix product states and mean-field theory to capture
magnetic order in quasi-1D cuprates"\
Q. Staelens, D. Verraes, D. Vrancken, T. Braeckevelt, J. Haegeman, V.
Van Speybroeck (2026) (arXiv:2602.21695v1)

The workflow consists of:

1.  DFT structural optimization
2.  DFT ground-state calculations
3.  Long-wave limit (LWL) calculations
4.  Band structure calculations
5.  Maximally localized Wannier functions (MLWFs)
6.  cRPA interaction calculations
7.  Extraction of hopping parameters

The resulting Hubbard parameters are used in tensor-network simulations
(performed separately using HubbardTN.jl, https://github.com/qstaelens/HubbardTN.jl).

------------------------------------------------------------------------

## Repository Structure

Structural Relaxation

Downfolding\
├── SCF\
├── LWL\
├── BANDS\
├── MLWFS\
├── cRPA\
└── Hopping

------------------------------------------------------------------------

## Computational Details

All calculations were performed using:

-   VASP v6.4.2
-   PBE exchange-correlation functional
-   PAW pseudopotentials
-   Wannier90 v3.1 via VASP2Wannier90 interface

Example module load command:

module load
VASP/6.4.2-gomkl-2023a-VASPsol-20210413-vtst-197-Wannier90-3.1.0

------------------------------------------------------------------------

## Workflow Description

### 1. Structural Relaxation

Location: Structural Relaxation

Purpose: Full relaxation of lattice parameters and internal atomic
coordinates.

Input files: 
- INCAR
- KPOINTS
- POSCAR
- POTCAR

Output: 
- CONTCAR (relaxed structure)

------------------------------------------------------------------------

### 2. SCF (Ground-State Calculation)

Folder: Downfolding/SCF/

Input files: 
- INCAR
- KPOINTS
- POSCAR
- POTCAR

Output: 
- WAVECAR
- CHGCAR

------------------------------------------------------------------------

### 3. LWL (Long-Wave Limit)

Folder: Downfolding/LWL/

Input files: 
- INCAR
- KPOINTS
- POSCAR
- POTCAR
- WAVECAR (from SCF)

Output: 
- WAVECAR
- WAVEDER

------------------------------------------------------------------------

### 4. BANDS (Band Structure)

Folder: Downfolding/BANDS/

Input files: 
- INCAR
- KPOINTS (high-symmetry path)
- POSCAR
- POTCAR
- WAVECAR (from LWL)
- CHGCAR (from LWL)

------------------------------------------------------------------------

### 5. MLWFS (Maximally Localized Wannier Functions)

Folder: Downfolding/MLWFS/

Input files: 
- INCAR
- KPOINTS
- POSCAR
- POTCAR
- WAVECAR (from LWL)
- WAVEDER (from LWL)

Additional scripts: 
- interpolating.py
- LMdecom.py

------------------------------------------------------------------------

### 6. cRPA

Folder: Downfolding/cRPA/

Input files: 
- INCAR
- KPOINTS
- POSCAR
- POTCAR
- WAVECAR (from LWL)
- WAVEDER (from LWL)

------------------------------------------------------------------------

### 7. Hopping

Folder: Downfolding/Hopping/

Purpose: Extraction of tight-binding hopping parameters from MLWFs.

Required files: 
 - Hopping.py
- toolbox.py

------------------------------------------------------------------------

## Contact

Quentin Staelens\
Ghent University\
Quentin.Staelens@UGent.be
