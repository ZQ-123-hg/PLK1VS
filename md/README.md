# Comprehensive Molecular Dynamics Simulation Workflow

A complete molecular dynamics (MD) simulation pipeline for protein-ligand binding studies with automated free energy calculations.

## Overview

This repository contains a comprehensive molecular dynamics simulation workflow designed for protein-ligand interaction studies. The pipeline implements a standardized protocol from initial structure preparation through binding free energy calculation using MM/GBSA methodology.

## Workflow Components

### 1. System Preparation

**Initial Structure**
- `2rku.pdb` - Starting protein structure from Protein Data Bank

**Prepared Systems**
- `protein.gro` - Processed protein structure
- `complex.gro` - Complete protein-ligand complex system

### 2. Ligand Library

Multiple ligand candidates are evaluated in this study:
- `BDE30671203_8.pdbqt` - Compound BDE30671203 (pose 8)
- `HIT104364049_11_out.pdbqt` - Compound HIT104364049 (pose 11)
- `Onvansertib_original_out.pdbqt` - Reference compound Onvansertib
- `PB4767006058_out.pdbqt` - Compound PB4767006058

### 3. Parameterization and Topology

**Force Field Parameters**
- `ligand.itp` - GROMACS-compatible ligand topology files
- `topol.top` - Complete system topology including protein and ligands

### 4. Simulation Protocol

The workflow follows a standard MD simulation protocol with the following stages:

#### Energy Minimization
- `em.mdp` - Initial energy minimization parameters
- `em2.mdp` - Secondary minimization stage

#### System Equilibration
- `ions.mdp` - Ion addition and neutralization
- `nvt.mdp` - Canonical (NVT) ensemble equilibration
- `npt.mdp` - Isothermal-isobaric (NPT) ensemble equilibration

#### Production Simulation
- `md.mdp` - Production molecular dynamics parameters

### 5. Automation and Analysis

**Workflow Management**
- `md_preprocess.sh` - Automated simulation execution script

**Free Energy Calculation**
- `mmpbsa.in` - MM/GBSA binding free energy calculation parameters

## Features

- **Standardized Protocol**: Implements best practices for MD simulation workflows
- **Multiple Ligand Support**: Simultaneous evaluation of multiple compound candidates  
- **Automated Processing**: Shell script automation for reproducible execution
- **Binding Affinity Prediction**: MM/GBSA methodology for quantitative binding assessment
- **Quality Control**: Multi-stage energy minimization and equilibration protocols

## Requirements

- GROMACS (recommended version 2020 or later)
- AmberTools (for MM/GBSA calculations)
- Standard computational chemistry utilities

## Usage

1. **System Setup**: Ensure all input files are properly prepared
2. **Execution**: Run the automated workflow using `bash md_preprocess.sh`
3. **Analysis**: Evaluate binding free energies from MM/GBSA output
4. **Results**: Compare ligand binding affinities for lead compound selection

## Output

The workflow generates comprehensive simulation trajectories and binding free energy estimates for each ligand candidate, enabling quantitative comparison of protein-ligand interactions and rational drug design decisions.

## License

[Specify appropriate license]

## Citation

[Add relevant citations for methodologies and software used]
