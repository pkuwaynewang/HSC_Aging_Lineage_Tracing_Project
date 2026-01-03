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
├── Python/ # python scripts for running pipelines
config/
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


## Cell-type clustering and tree visualization

Assuming the coupling score matrix has already been generated from the
previous analysis, the following code performs hierarchical clustering
of cell types based on their pairwise distances and visualizes the result
as a circular tree.

```r
source("scripts/r/treeplot_celltype_clustering.R")

coupling_mat <- read.table(
  "outputs/coupling_score/coupling_score_coupling_score.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t",
  check.names = FALSE
)

obj <- build_celltype_treeplot(
  coupling_mat,
  out_path = "outputs/coupling_score/celltype_tree.png"
)

print(obj$plot)
```

## Cell fate bias calculation

Cell fate bias is quantified based on a barcode-by–cell-type matrix.
For each clone, a tab-separated input file (`cellular_fate_bias.txt`)
should be prepared with the following format:

```
allele clone_size fate_size
barcode_1  356  310
barcode_2  1342 15
...
```

- **allele**: unique barcode (clone identifier)
- **clone_size**: total number of cells in the clone
- **fate_size**: number of cells belonging to the target cell fate within the clone

Assuming that `cellular_fate_bias.txt` has been prepared, run the following command:

```bash
python scripts/python/cell_fate_bias.py data/example_cell_fate.txt

```

Output

The script generates an output file:

```
/results/example_cell_fate_results.txt
```
An example output is shown below:

```
allele	clone_size	fate_size	all_cell_number	fate_number	P1	P2	P_min	Q_value	bias	fate_bias
Barcode_1	8	0	6279	2302	0.025834420337706096	0.9741655796622939	0.025834420337706096	0.041335072540329754	1.3836812958773301	-1.3836812958773301
Barcode_2	8	0	6279	2302	0.025834420337706096	0.9741655796622939	0.025834420337706096	0.041335072540329754	1.3836812958773301	-1.3836812958773301
```
P1 / P2: probabilities for fate enrichment and depletion, respectively

P_min: minimum of P1 and P2

Q_value: Benjamini–Hochberg–adjusted p-value

bias: signed enrichment score

fate_bias: direction-specific fate bias score

Clones exhibiting significant cell fate bias can be identified by applying
user-defined thresholds on Q_value and fate_bias.


## Data availability
Raw sequencing data are not included in this repository.
They will be made available via a public repository (e.g. GEO) upon publication.

## Contact
Maintained by Yuting Wang
