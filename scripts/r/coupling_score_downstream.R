#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Downstream analysis: coupling score between mature lineages based on CARLIN/DARLIN
# allele sharing (publication-ready script)
#
# Inputs:
#   1) Allele bank CSV (e.g., alleles_CA_bank.csv)
#   2) Sample list (one sample per line)
#   3) Per-sample count files: <base_dir>/<sample>/Barcode_UMI_number.txt
#   4) Configurable thresholds (min clone size, rare allele threshold, etc.)
#
# Outputs:
#   - Coupling score matrix (cosine-normalized)
#   - BH-adjusted empirical p-value matrix (shuffle-based)
#   - Heatmap PDF/PNG
#
# Notes:
#   - This script does NOT include raw data. Use relative paths within the repo.
#   - For large datasets / many shuffles, consider running on HPC.
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(stringr)
  library(dplyr)
  library(ggplot2)
  library(reshape2)
  library(scales)
})

# ---------------------------
# CLI arguments
# ---------------------------
option_list <- list(
  make_option(c("--allele_bank_ca"), type = "character", default = "data/allele_banks/alleles_CA_bank.csv",
              help = "CSV: allele bank for CA (must include columns: allele, count). [default %default]"),
  make_option(c("--allele_bank_ta"), type = "character", default = NULL,
              help = "Optional CSV: allele bank for TA (not used in current steps). [default %default]"),
  make_option(c("--sample_list"), type = "character", default = "data/sample_lists/file.name.txt",
              help = "Text file: one sample ID per line. [default %default]"),
  make_option(c("--base_dir"), type = "character", default = "results/C_regions/CARLIN/results_cutoff_override_3",
              help = "Base directory containing <sample>/Barcode_UMI_number.txt. [default %default]"),
  make_option(c("--n_samples"), type = "integer", default = 14,
              help = "Number of samples/cell-types to use (first N lines of sample_list). [default %default]"),
  make_option(c("--celltype_names"), type = "character", default = NULL,
              help = "Comma-separated cell type names (length must equal n_samples). If NULL, uses sample IDs."),
  make_option(c("--min_clone_size"), type = "integer", default = 2,
              help = "Keep clones with rowSum(counts) > min_clone_size. [default %default]"),
  make_option(c("--rare_bank_count_lt"), type = "double", default = 1,
              help = "Define rare alleles as allele_bank_count < this threshold. [default %default]"),
  make_option(c("--min_events_gt"), type = "integer", default = 1,
              help = "Keep alleles with number of events (comma count) > this threshold. [default %default]"),
  make_option(c("--cap_for_plot"), type = "double", default = 0.3,
              help = "Cap coupling score for visualization. [default %default]"),
  make_option(c("--n_shuffles"), type = "integer", default = 10000,
              help = "Number of shuffles for empirical p-values. [default %default]"),
  make_option(c("--seed"), type = "integer", default = 100,
              help = "Random seed. [default %default]"),
  make_option(c("--out_dir"), type = "character", default = "outputs/coupling_score",
              help = "Output directory. [default %default]"),
  make_option(c("--prefix"), type = "character", default = "coupling_score",
              help = "Output file prefix. [default %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

# ---------------------------
# Helpers
# ---------------------------
stop_if_missing <- function(path, label = "file") {
  if (!file.exists(path)) stop(sprintf("Missing %s: %s", label, path), call. = FALSE)
}

dir.create(opt$out_dir, recursive = TRUE, showWarnings = FALSE)

# ---------------------------
# Core functions
# ---------------------------

read_allele_bank <- function(path) {
  stop_if_missing(path, "allele bank")
  df <- fread(path)
  # Expect columns at least: allele, count
  required <- c("allele", "count")
  missing_cols <- setdiff(required, colnames(df))
  if (length(missing_cols) > 0) {
    stop("Allele bank is missing required columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
  }
  df[, .(allele, count)]
}

read_sample_ids <- function(path, n = NULL) {
  stop_if_missing(path, "sample list")
  ids <- readLines(path, warn = FALSE)
  ids <- ids[nzchar(trimws(ids))]
  if (!is.null(n)) ids <- ids[seq_len(min(n, length(ids)))]
  if (length(ids) == 0) stop("No sample IDs found in sample list.", call. = FALSE)
  ids
}

read_counts_one_sample <- function(base_dir, sample_id) {
  f <- file.path(base_dir, sample_id, "Barcode_UMI_number.txt")
  stop_if_missing(f, sprintf("count file for sample %s", sample_id))
  df <- fread(f, header = FALSE, col.names = c("allele", sample_id))
  df
}

merge_counts <- function(list_of_tables) {
  merged <- Reduce(function(x, y) merge(x, y, by = "allele", all = TRUE), list_of_tables)
  merged[is.na(merged)] <- 0
  merged
}

filter_rare_alleles <- function(counts_df, allele_bank_df, rare_bank_count_lt, min_events_gt) {
  df <- merge(counts_df, allele_bank_df, by = "allele", all.x = TRUE)
  df[is.na(df)] <- 0
  df$event <- str_count(df$allele, ",")
  # rare: bank count < threshold AND events > threshold
  df_rare <- df %>%
    filter(.data$count < rare_bank_count_lt, .data$event > min_events_gt)
  df_rare
}

# Normalize and compute coupling score
#   1) column-normalize counts (by column sums)
#   2) treat each clone row as a composition (normalize by row sum)
#   3) compute cosine similarity between columns using dot products
compute_coupling_score <- function(count_mat, eps = 1e-10) {
  # count_mat: clones x celltypes
  stopifnot(is.matrix(count_mat) || is.data.frame(count_mat))
  X <- as.matrix(count_mat)

  # Step 1: normalize by column sums (like your apply(...,1, x/colSums))
  col_sums <- colSums(X)
  if (any(col_sums == 0)) stop("At least one column sum is 0; cannot normalize.", call. = FALSE)
  X1 <- sweep(X, 2, col_sums, "/")

  # Step 2: normalize each row to sum 1
  rs <- rowSums(X1) + eps
  X2 <- X1 / rs

  # Step 3: dot products between columns
  G <- crossprod(X2)  # = t(X2) %*% X2

  # Step 4: cosine-normalize
  denom <- sqrt(outer(diag(G), diag(G), "*")) + eps
  S <- G / denom

  S
}

plot_coupling_heatmap <- function(score_mat, cap, out_path, title = "Heatmap of coupling score") {
  score_cap <- pmin(score_mat, cap)
  df_long <- melt(score_cap)

  p <- ggplot(df_long, aes(Var1, Var2, fill = value)) +
    geom_tile() +
    geom_text(aes(label = sprintf("%.2f", value)), color = "black", size = 3) +
    scale_fill_gradientn(
      colors = c("white", "white", "darkred"),
      values = rescale(c(0, cap * 0.05, cap * 0.85)),
      limits = c(0, cap),
      guide = "colorbar"
    ) +
    labs(x = NULL, y = NULL, title = title) +
    theme_bw() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid = element_blank()
    )

  ggsave(out_path, plot = p, width = 7.5, height = 6.5, dpi = 300)
}

# Empirical p-values via shuffling within each column (preserves marginal distribution per column)

empirical_pvalues <- function(count_mat, real_score, n_shuffles, seed = 1, eps = 1e-10, verbose_every = 250) {
  set.seed(seed)
  X <- as.matrix(count_mat)

  # Precompute column sums for the step-1 normalization (same across shuffles if we only permute rows within columns)
  col_sums <- colSums(X)
  if (any(col_sums == 0)) stop("At least one column sum is 0; cannot normalize.", call. = FALSE)

  # Count exceedances: shuffled_score >= real_score
  exceed <- matrix(0L, nrow = ncol(X), ncol = ncol(X))
  for (i in seq_len(n_shuffles)) {
    Xs <- apply(X, 2, sample)  # shuffle each column independently
    S_shuf <- compute_coupling_score(Xs, eps = eps)
    exceed <- exceed + (S_shuf >= real_score)

    if (verbose_every > 0 && (i %% verbose_every == 0)) {
      message(sprintf("Completed shuffle %d / %d", i, n_shuffles))
    }
  }

  # +1 smoothing
  p <- (exceed + 1) / (n_shuffles + 1)
  p_adj <- matrix(p.adjust(as.vector(p), method = "BH"), nrow = nrow(p), byrow = FALSE)
  dimnames(p_adj) <- dimnames(p)

  list(p_raw = p, p_adj = p_adj)
}

# ---------------------------
# Main
# ---------------------------

message("Reading allele bank...")
allele_bank_ca <- read_allele_bank(opt$allele_bank_ca)

message("Reading sample list...")
sample_ids <- read_sample_ids(opt$sample_list, n = opt$n_samples)

if (!is.null(opt$celltype_names)) {
  ct <- strsplit(opt$celltype_names, ",")[[1]] |> trimws()
  if (length(ct) != length(sample_ids)) {
    stop("--celltype_names length must equal n_samples (", length(sample_ids), ").", call. = FALSE)
  }
  celltype_names <- ct
} else {
  # Default: use provided sample IDs
  celltype_names <- sample_ids
}

message("Loading per-sample counts...")
count_tables <- lapply(sample_ids, function(sid) read_counts_one_sample(opt$base_dir, sid))

message("Merging counts across samples...")
merged_counts <- merge_counts(count_tables)

message("Filtering rare alleles using allele bank + event threshold...")
df_rare_all <- filter_rare_alleles(
  counts_df = merged_counts,
  allele_bank_df = allele_bank_ca,
  rare_bank_count_lt = opt$rare_bank_count_lt,
  min_events_gt = opt$min_events_gt
)

# Prepare count matrix: clones x celltypes
# merged_counts columns: allele + sample_ids
df_rare_all <- as.data.frame(df_rare_all)
rownames(df_rare_all) <- df_rare_all$allele

# Keep only the count columns in the original order
count_df <- df_rare_all[, sample_ids, drop = FALSE]
colnames(count_df) <- celltype_names

# Filter by clone size threshold (row sum on raw counts)
count_df2 <- count_df[rowSums(count_df) > opt$min_clone_size, , drop = FALSE]
if (nrow(count_df2) == 0) stop("No clones remain after filtering (min_clone_size too high?).", call. = FALSE)

message(sprintf("Clones kept after filtering: %d", nrow(count_df2)))

# Compute coupling score
message("Computing coupling score matrix...")
score_mat <- compute_coupling_score(count_df2)

rownames(score_mat) <- colnames(score_mat) <- celltype_names

# Save coupling score
out_score_tsv <- file.path(opt$out_dir, paste0(opt$prefix, "_coupling_score.tsv"))
write.table(score_mat, out_score_tsv, sep = "\t", quote = FALSE, col.names = NA)
message("Saved coupling score: ", out_score_tsv)

# Plot heatmap
out_heatmap <- file.path(opt$out_dir, paste0(opt$prefix, "_coupling_heatmap.png"))
plot_coupling_heatmap(score_mat, cap = opt$cap_for_plot, out_path = out_heatmap)
message("Saved heatmap: ", out_heatmap)

# Empirical p-values
message(sprintf("Computing empirical p-values with %d shuffles (this may take a while)...", opt$n_shuffles))
pv <- empirical_pvalues(
  count_mat = count_df2,
  real_score = score_mat,
  n_shuffles = opt$n_shuffles,
  seed = opt$seed,
  verbose_every = 250
)

out_padj_tsv <- file.path(opt$out_dir, paste0(opt$prefix, "_padj.tsv"))
write.table(pv$p_adj, out_padj_tsv, sep = "\t", quote = FALSE, col.names = NA)
message("Saved BH-adjusted empirical p-values: ", out_padj_tsv)

message("Done.")

