# Molecular Docking and Post-Processing Workflow

## Overview

This repository contains a comprehensive molecular docking and post-processing workflow designed for protein-ligand interaction analysis. The pipeline incorporates automated docking simulations followed by systematic filtering to identify optimal binding poses for PLK1 inhibition studies.

## Workflow Components

### Protein Receptors
- **2RKU.pdbqt**: Prepared protein receptor structure for docking simulations
- **2YAC.pdbqt**: Alternative prepared protein receptor structure for comparative docking analysis

### Configuration Files
- **config_2RKU.txt**: Docking parameters and binding site coordinates for 2RKU receptor
- **config_2YAC.txt**: Docking parameters and binding site coordinates for 2YAC receptor

### Distributed Processing Scripts
- **distributed_prepare_ligand.sh**: Enables parallel execution of ligand preparation procedures
- **distributed_unidock.sh**: Facilitates distributed docking calculations across multiple computing nodes

### Core Filtering System
- **liggrep_multi-threaded.py**: Multi-threaded filtering engine that processes docking results
- **liggrep_filters_PLK1.json**: Interaction criteria specification file for PLK1-specific binding analysis

## Process Flow

1. **Ligand Preparation**: Parallel processing of input ligands using the distributed preparation script
2. **Molecular Docking**: Automated docking simulations against prepared protein receptors (2RKU and 2YAC)
3. **Results Filtering**: Systematic identification and selection of docking poses based on predefined interaction criteria
4. **Contact Analysis**: Evaluation of specific contacts with key binding site amino acid residues crucial for PLK1 inhibition

## Key Features

- **Parallel Processing**: Distributed computing capabilities for enhanced performance
- **Multi-threaded Filtering**: Efficient post-processing of large-scale docking results
- **Specific Target Focus**: Optimized for PLK1 inhibitor identification
- **Systematic Selection**: Criteria-based pose selection ensuring biologically relevant interactions

## Technical Requirements

- Compatible docking software environment
- Multi-core processing capabilities for optimal performance
- Sufficient storage for intermediate and final results

## Usage

The workflow is designed to systematically process molecular docking results and identify compounds with optimal binding characteristics for PLK1 inhibition research.

---
