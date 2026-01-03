#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# HashTag / Cell hashing demultiplexing from 10x Feature Barcode FASTQs
#
# This script:
#   1) Extracts 10x cell barcode (CB) from Read1 (default: first 16 bp)
#   2) Extracts hashtag/feature barcode (FB) from Read2 (default: first 15 bp)
#   3) Matches FB to reference hashtag barcodes (allowing mismatches)
#   4) Builds a CB x HashTag count matrix, CPM-normalizes it
#   5) Assigns each CB to a HashTag (or "unassigned") based on a ratio threshold
#
# Inputs:
#   --read1: R1 FASTQ (cell barcode)
#   --read2: R2 FASTQ (feature barcode / hashtag)
#   --barcode_ref: CSV with columns {ID, Barcode}
#   --hashtag: comma-separated hashtag IDs to consider (e.g., A0251,A0252)
#
# Outputs (in outdir):
#   <prefix>CB_hash_count.csv
#   <prefix>CB_hash_count_CPM.csv
#   <prefix><HashTagID>_distribution.pdf   (one per hashtag; min-max scaled)
#   <prefix>CB_hash_assign.csv
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ShortRead)
  library(optparse)
  library(edgeR)
  library(ggplot2)
  library(dplyr)
})

# ---------------------------
# CLI arguments
# ---------------------------
option_list <- list(
  make_option(c("-a", "--read1"), type = "character", default = NULL,
              help = "R1 FASTQ (contains 10x cell barcode)"),
  make_option(c("-b", "--read2"), type = "character", default = NULL,
              help = "R2 FASTQ (contains hashtag/feature barcode sequence)"),
  make_option(c("-o", "--outdir"), type = "character", default = NULL,
              help = "Output directory"),
  make_option(c("-f", "--hashtag"), type = "character", default = NULL,
              help = "Comma-separated HashTag IDs to process (e.g., A0251,A0252)"),
  make_option(c("-r", "--barcode_ref"), type = "character", default = NULL,
              help = "Reference CSV containing columns: ID, Barcode"),
  make_option(c("--cb_len"), type = "integer", default = 16,
              help = "Length of cell barcode extracted from R1 [default %default]"),
  make_option(c("--fb_len"), type = "integer", default = 15,
              help = "Length of feature/hashtag barcode extracted from R2 [default %default]"),
  make_option(c("--max_mismatch"), type = "integer", default = 1,
              help = "Maximum mismatches allowed for hashtag matching [default %default]"),
  make_option(c("--min_reads_per_cb"), type = "integer", default = 4,
              help = "Minimum matched reads per CB to be kept for a hashtag [default %default]"),
  make_option(c("--ratio_cutoff"), type = "double", default = 10,
              help = "Assignment cutoff: CPM(h) / (sum CPM(others)+eps) > ratio_cutoff [default %default]"),
  make_option(c("--eps"), type = "double", default = 1e-6,
              help = "Small epsilon to avoid division by zero [default %default]"),
  make_option(c("--prefix"), type = "character", default = "",
              help = "Output filename prefix (optional) [default '%default']")
)

opt_parser <- OptionParser(
  usage = "Usage: %prog [options]\nDescription: Extract hash tag information and 10x cell barcode for all cells.",
  option_list = option_list
)
opt <- parse_args(opt_parser)

# ---------------------------
# Basic validation
# ---------------------------
die <- function(msg) {
  print_help(opt_parser)
  stop(msg, call. = FALSE)
}

if (is.null(opt$read1)) die("Please provide --read1 (R1 FASTQ).")
if (is.null(opt$read2)) die("Please provide --read2 (R2 FASTQ).")
if (is.null(opt$outdir)) die("Please provide --outdir.")
if (is.null(opt$hashtag)) die("Please provide --hashtag (comma-separated IDs).")
if (is.null(opt$barcode_ref)) die("Please provide --barcode_ref (CSV with ID,Barcode).")

if (!file.exists(opt$read1)) stop("R1 FASTQ not found: ", opt$read1, call. = FALSE)
if (!file.exists(opt$read2)) stop("R2 FASTQ not found: ", opt$read2, call. = FALSE)
if (!file.exists(opt$barcode_ref)) stop("barcode_ref not found: ", opt$barcode_ref, call. = FALSE)

# Normalize outdir and create it
outdir <- opt$outdir
if (!grepl("/$", outdir)) outdir <- paste0(outdir, "/")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

prefix <- opt$prefix
if (nzchar(prefix) && !grepl("_$", prefix)) prefix <- paste0(prefix, "_")

hash_tag_id <- trimws(unlist(strsplit(opt$hashtag, ",")))
if (length(hash_tag_id) == 0) stop("No hashtag IDs parsed from --hashtag.", call. = FALSE)

message("Hashtag IDs: ", paste(hash_tag_id, collapse = ", "))
message("R1: ", opt$read1)
message("R2: ", opt$read2)
message("Reference: ", opt$barcode_ref)
message("Output dir: ", outdir)

# ---------------------------
# Helper functions
# ---------------------------
min_max_fun <- function(v) {
  v.min <- min(v)
  v.max <- max(v)
  if (v.max == v.min) return(rep(0, length(v)))
  (v - v.min) / (v.max - v.min)
}

read_hashtag_reference <- function(path) {
  ref <- read.csv(path, header = TRUE, check.names = FALSE, stringsAsFactors = FALSE)
  required <- c("ID", "Barcode")
  missing <- setdiff(required, colnames(ref))
  if (length(missing) > 0) {
    stop("barcode_ref is missing required columns: ", paste(missing, collapse = ", "), call. = FALSE)
  }
  ref
}

extract_prefix_seq <- function(fastq_path, len) {
  # Read sequences from FASTQ and take the first len bases
  seqs <- as.character(sread(readFastq(fastq_path)))
  substr(seqs, 1, len)
}

match_hashtag_counts <- function(CB, FB, ref_barcode, max_mismatch, min_reads_per_cb) {
  # Compute edit distance between the reference barcode and all FB reads
  # Keep reads with <= max_mismatch (implemented as < max_mismatch+1)
  d <- adist(ref_barcode, FB)
  keep <- d <= max_mismatch
  if (!any(keep)) {
    return(data.frame(CB = character(0), value = numeric(0)))
  }
  cb_kept <- CB[keep]
  tab <- table(cb_kept)
  tab <- tab[tab >= min_reads_per_cb]
  data.frame(CB = names(tab), value = as.numeric(tab), stringsAsFactors = FALSE)
}

assign_hashtag <- function(cpm_mat, ratio_cutoff, eps) {
  hash_ids <- colnames(cpm_mat)
  out <- data.frame(row.names = rownames(cpm_mat))
  out$assignment <- "unassigned"

  # ratio columns for transparency/debugging
  for (h in hash_ids) {
    others <- setdiff(hash_ids, h)
    out[[h]] <- cpm_mat[, h] / (rowSums(cpm_mat[, others, drop = FALSE]) + eps)
  }

  # Assign if any ratio > ratio_cutoff; pick the first passing hashtag
  out <- out %>%
    rowwise() %>%
    mutate(
      assignment = {
        passing <- hash_ids[c_across(all_of(hash_ids)) > ratio_cutoff]
        if (length(passing) > 0) passing[1] else assignment
      }
    ) %>%
    ungroup()

  out$CB <- rownames(out)
  out
}

# ---------------------------
# Main
# ---------------------------
ref <- read_hashtag_reference(opt$barcode_ref)

# Validate IDs exist in reference
missing_ids <- setdiff(hash_tag_id, ref$ID)
if (length(missing_ids) > 0) {
  stop("The following hashtag IDs are not found in barcode_ref$ID: ",
       paste(missing_ids, collapse = ", "), call. = FALSE)
}

message("Reading FASTQs and extracting barcodes...")
message("Process R1 (CB)...")
CB <- extract_prefix_seq(opt$read1, opt$cb_len)

message("Process R2 (FB)...")
FB <- extract_prefix_seq(opt$read2, opt$fb_len)

if (length(CB) != length(FB)) {
  stop("R1 and R2 have different numbers of reads after loading: ",
       length(CB), " vs ", length(FB), call. = FALSE)
}

message("Matching hashtags and counting reads per CB...")
cell_hash_count_list <- list()

for (b in hash_tag_id) {
  ref_barcode <- ref$Barcode[ref$ID == b][1]
  message(b, ": reference barcode = ", ref_barcode)

  df_count <- match_hashtag_counts(
    CB = CB,
    FB = FB,
    ref_barcode = ref_barcode,
    max_mismatch = opt$max_mismatch,
    min_reads_per_cb = opt$min_reads_per_cb
  )

  if (nrow(df_count) == 0) {
    # still create an empty column later via merge
    df_count <- data.frame(CB = character(0), value = numeric(0), stringsAsFactors = FALSE)
  }

  colnames(df_count)[2] <- b
  cell_hash_count_list[[b]] <- df_count
}

# Merge CB x hashtag counts
cell_hash_count_merge <- Reduce(function(x, y) merge(x, y, by = "CB", all = TRUE), cell_hash_count_list)
cell_hash_count_merge[is.na(cell_hash_count_merge)] <- 0
rownames(cell_hash_count_merge) <- cell_hash_count_merge$CB
cell_hash_count_merge$CB <- NULL

# Write raw counts
out_counts <- paste0(outdir, prefix, "CB_hash_count.csv")
write.table(cell_hash_count_merge, file = out_counts, sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)
message("Wrote raw count matrix: ", out_counts)

# CPM normalization
cpm_mat <- edgeR::cpm(cell_hash_count_merge)
out_cpm <- paste0(outdir, prefix, "CB_hash_count_CPM.csv")
write.table(cpm_mat, file = out_cpm, sep = ",",
            col.names = TRUE, row.names = TRUE, quote = FALSE)
message("Wrote CPM-normalized matrix: ", out_cpm)

# Min-max normalize per cell for visualization only
cpm_minmax <- t(apply(cpm_mat, 1, min_max_fun))
cpm_minmax <- as.data.frame(cpm_minmax, check.names = FALSE)

# Plot distributions per hashtag
message("Saving hashtag distribution plots...")
for (h in colnames(cpm_minmax)) {
  p <- ggplot(cpm_minmax, aes(x = .data[[h]])) +
    geom_histogram(bins = 100) +
    theme_bw() +
    labs(title = paste0(h, " HashTag distribution (min-max scaled)"),
         x = "Min-max scaled CPM", y = "Number of cells")

  ggsave(
    filename = paste0(outdir, prefix, h, "_distribution.pdf"),
    plot = p, width = 8, height = 6
  )
}

# Assign hashtag per CB using ratio rule
message("Assigning hashtags to cell barcodes...")
CB_hash_assign <- assign_hashtag(
  cpm_mat = as.data.frame(cpm_mat, check.names = FALSE),
  ratio_cutoff = opt$ratio_cutoff,
  eps = opt$eps
)

out_assign <- paste0(outdir, prefix, "CB_hash_assign.csv")
write.table(CB_hash_assign, file = out_assign, sep = ",",
            col.names = TRUE, row.names = FALSE, quote = FALSE)
message("Wrote assignment table: ", out_assign)

message("Done.")

