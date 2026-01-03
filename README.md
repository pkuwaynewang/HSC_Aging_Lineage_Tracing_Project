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

## External dependency

This analysis makes use of a previously published Snakemake workflow
for DARLIN/CARLIN lineage tracing data processing developed by
the ShouWenWang Lab:

- **snakemake_DARLIN**
- GitHub repository: https://github.com/ShouWenWang-Lab/snakemake_DARLIN

Users should install and configure the DARLIN pipeline according to
the instructions provided in the original repository.
This project calls the DARLIN Snakemake workflow externally and does
not modify the original pipeline.

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

## step 1: Call DARLIN allele
This project provides two levels of execution in bulk DARLIN analysis:

## Batch execution (recommended)
```
bash scripts/shell/Bulk_batch_run.sh

```
The batch script iterates over a list of samples and, for each sample:

Updates the configuration file (config/config.yaml) with the current sample.

Calls the single-run script (scripts/shell/run.sh) to execute the pipeline.

Single execution for both bulk and 10x scRNA-seq DARLIN analysis:
```
bash scripts/shell/Bulk_run.sh

bash scripts/shell/10x_scRNAseq_run.sh

```

## step 2: create middle files

After finishing above analysis, run the following codes in 
```
/path/to/CARLIN/results_cutoff_override_#

```
```
cp /scripts/shell/middle_file_1.sh ./
cp /scripts/shell/middle_file_2.sh ./

for dir in */; do
  # Check if it's a directory
  if [ -d "$dir" ]; then
    # Run code.sh with the directory path as an argument
    (cd "$dir" && bash ../middle_file_1.sh && bash ../middle_file_2.sh )
  fi
done
```
middle_file_1.sh with create  a text file, for each cell UMI followed by the barcode information, with the following format:

```
 Cell_UMI_1   Barcode_1
 Cell_UMI_2   Barcode_2
```
middle_file_2.sh with create  a file with each barcode followed by the number of cell (UMI) from this barcode:

```
Barcode_1     16
Barcode_2    12
```

## Step 3: Down stream analysis and ploting

```
Rscript scripts/r/coupling_score_downstream.R \
  --allele_bank_ca data/allele_banks/alleles_CA_bank.csv \
  --sample_list data/sample_lists/file.name.txt \
  --base_dir results/C_regions/CARLIN/results_cutoff_override_3 \
  --n_samples 14 \
  --min_clone_size 2 \
  --n_shuffles 10000 \
  --out_dir outputs/coupling_score \
  --prefix lifelong_C_regions
```
The code above will calculate and plot coupling score among cell types.



## Data availability
Raw sequencing data are not included in this repository.
They will be made available via a public repository (e.g. GEO) upon publication.

## Contact
Maintained by Yuting Wang
