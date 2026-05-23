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
- **P1 / P2**: probabilities for fate enrichment and depletion, respectively
- **P_min**: minimum of P1 and P2
- **Q_value**: Benjamini–Hochberg–adjusted p-value
- **bias**: signed enrichment score
- **fate_bias**: direction-specific fate bias score

Clones exhibiting significant cell fate bias can be identified by applying
user-defined thresholds on Q_value and fate_bias.


## Cell hash demultiplex 

### HashTag reference file

The HashTag reference file (`BioLegend_TotalA_Human_Barcode.csv`) must be
a comma-separated table containing at least the following columns:

- `ID`: HashTag identifier (e.g. A0301), matching the values provided to `--hashtag`
- `Barcode`: Feature barcode sequence corresponding to each HashTag

Additional annotation columns are allowed and will be ignored by the script.

```
Rscript scripts/r/hash_tag_demux.R \
  --read1 path/to/feature_barcode_R1.fastq.gz \
  --read2 path/to/feature_barcode_R2.fastq.gz \
  --barcode_ref /data/BioLegend_TotalA_Mouse_hashtag_Barcode.csv \
  --hashtag A0301,A0302,A0303, A0304, A0305, A0306 \
  --outdir outputs/hashtag_demux \
  --max_mismatch 1 \
  --min_reads_per_cb 4 \
  --ratio_cutoff 10 \
  --prefix test1
```
### Example output

The final demultiplexing result (`CB_hash_assign.csv`) contains one row per
10x cell barcode. Cells that do not meet the assignment criteria are labeled
as `unassigned`.

```
assignment,A0301,A0302,A0303,A0304,A0305,A0306,CB
A0301,99.99,0,0,0,0,0.0,AAACCTGAGTACGTTC
unassigned,46.4,53.6,0,0,0,0,AAACCTGAGTTCGATT
```


### Tubulin Polarization Score (Angular Bin Analysis)

This repository provides a lightweight Python tool to quantify **tubulin signal distribution (polarization)** in single-cell confocal images using an angular binning strategy.

The method is designed for **within-experiment comparisons under consistent imaging conditions** and outputs both total tubulin signal content and a polarization score describing spatial asymmetry.

---

## Concept and Rationale

Given a single, centered cell image:

1. The cell center is fixed at the **image center** (images are pre-cropped and centered using CellProfiler).
2. The image is divided into **angular bins** (θ bins, default = 180) spanning \([-π, π)\).
3. Tubulin signal is defined by an **absolute intensity threshold**.
4. For each angular bin, we compute:
   - **Total signal intensity** (sum of pixel intensities above threshold)
   - **Signal pixel count**

This produces a 1D circular profile describing how tubulin signal is distributed around the cell.

## Polarization Score Definition

Let `sumI[i]` be the total tubulin signal in angular bin *i*.

### Step 1: Relative distribution
We convert the per-bin sums into a probability distribution:

```
p_i = sumI[i] / sum(sumI)
```

This separates **distribution shape** from **total signal content**.

### Step 2: Concentration metrics
We compute complementary measures of spatial asymmetry:

- **Entropy concentration**
```
C_entropy = 1 − (Shannon entropy of p) / log(N)
```
- 0 → uniform distribution  
- 1 → highly concentrated

- **Top-k fraction**
Fraction of total signal contained in the top 5% of angular bins.

- **Maximum angular gap**
Largest contiguous angular region with no signal (normalized by total bins).

### Step 3: Polarization score
The final polarization score is a weighted sum:

```
Polarization Score =
w_entropy * C_entropy+
w_topk * TopK_fraction+
w_gap * Max_gap_fraction
```

**Interpretation:**

- **Higher score** → tubulin signal is more localized, asymmetric, or polarized  
- **Lower score** → tubulin signal is more evenly distributed and symmetric

## Output Metrics

Each image produces:

### Summary features
- `total_signal_sum` – total tubulin signal content
- `polarization_score` – spatial asymmetry metric
- `entropy_concentration`
- `topk_fraction`
- `max_bin_fraction`
- `coverage`
- `gap_count`
- `max_gap_norm`

### Per-bin data
- `sumI` – per-bin total signal
- `count` – per-bin signal pixel count

An optional QC plot visualizes:
1. Raw image
2. Signal mask
3. Angular signal distribution

---

## Command Line Usage

```bash
python tubulin_polarization.py image.tiff \
  --bins 180 \
  --thr 6000 \
  --w-entropy 0.6 \
  --w-topk 0.3 \
  --w-gap 0.1
```

### shRNA Small-Scale Screening in Aged DP HSCs

The analysis pipeline is adapted from the published method (PMID: 35243374) with the following modification:

Instead of counting sequencing reads for each shRNA, we quantify the actual number of cells associated with each shRNA by incorporating a DNA-specific unique molecular identifier (UMI, 10 bp). This approach enables direct cell number estimation rather than read-based abundance. See the Methods section for details.

## Command Line Usage
```
bash /scripts/shell/shRNA_screening.sh \
sample_1_R1.fastq.gz \ # R1 file
sample_1_R2.fastq.gz \ # R2 file
 sample_1 \ # File name
/data/shRNA_screening_references/screening_reference # Reference
```

## Clonal diversity analysis 

The analysis will help to calculate the clonal diveristy changes using Shannon Entropy, Gini Index and Simpson Index.
Before running the analysis, a barcode x cell type matrix Df_matrix should be prepared first, with value as the cell number.

## Commonda Line Usage

```
source('./scripts/r/Clonal_diversity_analysis.R')
results <- sapply(Df_matrix, calculate_diversity)

```

## Data availability
Raw sequencing data are not included in this repository.
They will be made available via a public repository (e.g. GEO) upon publication.

## Contact
Maintained by Yuting Wang
