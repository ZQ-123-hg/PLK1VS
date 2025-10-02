# Computational Drug Discovery Pipeline

This repository provides a comprehensive computational pipeline for drug discovery and molecular analysis, integrating three complementary workflows for structure-based drug design and molecular property evaluation.

## Overview

### 1. Molecular Dynamics Simulation Workflow
A complete MD simulation pipeline for protein-ligand systems starting from initial protein structures. The workflow encompasses system preparation, ligand parameterization, topology generation, and standard simulation protocols including energy minimization, equilibration phases (NVT/NPT), and production runs. Binding free energy calculations using MM/GBSA methodology provide quantitative assessment of ligand-protein interactions for lead compound evaluation.

### 2. Molecular Docking and Post-Processing Pipeline  
A high-throughput molecular docking framework featuring automated ligand preparation and parallel docking execution against multiple protein targets. The workflow incorporates sophisticated post-processing filters based on specific protein-ligand interaction patterns, enabling systematic identification and selection of promising binding poses that satisfy predefined pharmacophore requirements for targeted protein inhibition.

### 3. Chemical Analysis and Filtering Tools
A collection of specialized Python utilities for comprehensive molecular property assessment and chemical space analysis. The toolkit includes drug-likeness evaluation (Lipinski's Rule of Five), ADMET property prediction, synthetic accessibility scoring, scaffold-based filtering for specific therapeutic areas, molecular similarity analysis, and chemical space clustering algorithms for compound library organization and diversity assessment.

## Applications
This integrated pipeline supports rational drug design through structure-based approaches, enabling researchers to efficiently screen compound libraries, optimize lead compounds, and evaluate their pharmacological profiles for drug development projects.
