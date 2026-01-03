## External dependency

This analysis makes use of a previously published Snakemake workflow
for DARLIN/CARLIN lineage tracing data processing developed by
the ShouWenWang Lab:

- **snakemake_DARLIN**
- GitHub repository: https://github.com/ShouWenWang-Lab/snakemake_DARLIN

Users should install and configure the DARLIN pipeline according to
the instructions provided in the original repository.
This project calls the DARLIN Snakemake workflow externally and does
not modify the original pipeline


# HSC_Aging_Lineage_Tracing_Project

Code for analysis of lineage tracing experiments.

## Overview
This repository contains R and shell scripts used for data processing,
analysis, and visualization in a lineage tracing project.

## Languages and tools
- R (>= 4.2)
- Bash / shell scripts

## Repository structure
```
scripts/
├── r/ # R scripts for analysis and visualization
├── shell/ # Shell scripts for running pipelines
data/
└── processed/ # Example or processed data (no raw data)
```
## Usage
Scripts are organized by analysis steps.
Detailed instructions will be added as the project is finalized.

## Data availability
Raw sequencing data are not included in this repository.
They will be made available via a public repository (e.g. GEO) upon publication.

## Contact
Maintained by Yuting Wang
