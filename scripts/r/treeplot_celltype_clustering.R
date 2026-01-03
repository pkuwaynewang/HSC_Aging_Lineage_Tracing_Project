#!/usr/bin/env Rscript

# ------------------------------------------------------------------------------
# Tree plot (circular) for clustering cell types based on coupling score matrix
#
# Input:
#   - coupling_mat: a square matrix/data.frame of coupling scores 
#     with rownames and colnames = cell type labels.		
#
# Method:
#   - Convert coupling score S to distance D = 1 - S
#   - Hierarchical clustering on D (default: complete linkage)
#   - Convert to phylo and plot using ggtree + ggtreeExtra
#
# Output:
#   - Circular dendrogram with colored tips and a small bar "fruit" track
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(ggplot2)
  library(ape)
  library(phylobase)
  library(treeio)
  library(ggtree)
  library(ggtreeExtra)
})

# ---------------------------
# Utilities
# ---------------------------

validate_coupling_matrix <- function(coupling_mat, allow_out_of_range = TRUE, tol = 1e-8) {
  if (is.data.frame(coupling_mat)) coupling_mat <- as.matrix(coupling_mat)

  if (!is.matrix(coupling_mat)) {
    stop("`coupling_mat` must be a matrix or data.frame.", call. = FALSE)
  }
  if (nrow(coupling_mat) != ncol(coupling_mat)) {
    stop("`coupling_mat` must be a square matrix (nrow == ncol).", call. = FALSE)
  }
  if (is.null(rownames(coupling_mat)) || is.null(colnames(coupling_mat))) {
    stop("`coupling_mat` must have rownames and colnames (cell type labels).", call. = FALSE)
  }
  if (!all(rownames(coupling_mat) == colnames(coupling_mat))) {
    stop("Row/column names must match and be in the same order.", call. = FALSE)
  }
  if (any(!is.finite(coupling_mat))) {
    stop("`coupling_mat` contains NA/Inf. Please clean the matrix first.", call. = FALSE)
  }

  # Optional range check
  if (!allow_out_of_range) {
    if (any(coupling_mat < -tol) || any(coupling_mat > 1 + tol)) {
      stop("`coupling_mat` has values outside [0, 1].", call. = FALSE)
    }
  }
  invisible(coupling_mat)
}

make_distance_from_coupling <- function(coupling_mat, cap_to_unit_interval = TRUE) {
  coupling_mat <- validate_coupling_matrix(coupling_mat)

  S <- coupling_mat

  # In case numerical noise produces tiny negatives / >1 values, optionally cap
  if (cap_to_unit_interval) {
    S <- pmin(pmax(S, 0), 1)
  }

  D <- 1 - S

  # Ensure distance matrix properties
  diag(D) <- 0
  D <- (D + t(D)) / 2  # enforce symmetry

  as.dist(D)
}

default_palette_14 <- function() {
  # Your original palette retained for reproducibility
  c(
    "#899499",
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728",
    "#9467bd", "#8c564b", "#e377c2",
    "#17becf", "#bcbd22",
    "#EE3377", "#009988", "#CC3311", "#33BBEE"
  )
}

cluster_to_treedata <- function(dist_obj, method = "complete", tip_value = 0.1) {
  hc <- hclust(dist_obj, method = method)

  tip_df <- data.frame(
    id    = hc$labels,
    value = tip_value,
    types = hc$labels,
    stringsAsFactors = FALSE
  )

  phy <- as.phylo(hc)
  phy4d_obj <- phylo4d(phy, tip_df)
  as.treedata(phy4d_obj)
}

plot_circular_tree <- function(tree_data,
                               palette = NULL,
                               layout = "circular",
                               show_branch_length = FALSE,
                               fruit_offset = 0.5,
                               title = NULL) {
  if (is.null(palette)) palette <- default_palette_14()

  p <- ggtree(
    tree_data,
    layout = layout,
    branch.length = if (show_branch_length) "branch.length" else "none"
  ) +
    geom_tiplab(aes(angle = angle, color = types), size = 3) +
    ggtreeExtra::geom_fruit(
      mapping = aes(x = value, color = types, fill = types),
      geom = ggplot2::geom_col,
      offset = fruit_offset,
      orientation = "y",
      stat = "identity",
      width = 0.7
    ) +
    scale_color_manual(values = palette) +
    scale_fill_manual(values = palette) +
    theme(legend.position = "none")

  if (!is.null(title)) p <- p + ggtitle(title)
  p
}

save_plot <- function(p, out_path, width = 7, height = 7, dpi = 300) {
  dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
  ggsave(out_path, plot = p, width = width, height = height, dpi = dpi)
}

# ---------------------------
# Main callable function
# ---------------------------

#' Build and plot circular clustering tree from a coupling score matrix.
#'
#' @param coupling_mat Square matrix of coupling scores (rownames/colnames = cell types)
#' @param hclust_method hclust linkage method (default "complete")
#' @param out_path Optional output path (.png/.pdf). If NULL, doesn't save.
#' @param palette Optional named/unnamed vector of colors (length >= number of tips)
#' @param ... Additional parameters passed to plotting helpers
#' @return A list containing: dist, hclust, treedata, plot
build_celltype_treeplot <- function(coupling_mat,
                                    hclust_method = "complete",
                                    out_path = NULL,
                                    palette = NULL,
                                    tip_value = 0.1,
                                    fruit_offset = 0.5,
                                    title = "Cell-type clustering based on coupling score") {

  dist_obj <- make_distance_from_coupling(coupling_mat)
  tree_data <- cluster_to_treedata(dist_obj, method = hclust_method, tip_value = tip_value)

  p <- plot_circular_tree(
    tree_data,
    palette = palette,
    fruit_offset = fruit_offset,
    title = title
  )

  if (!is.null(out_path)) {
    save_plot(p, out_path)
  }

  list(
    dist = dist_obj,
    treedata = tree_data,
    plot = p
  )
}

# ---------------------------
# Example usage (for script execution)
# ---------------------------
# If you want this file to run as a standalone script, you can source a matrix
# and call build_celltype_treeplot(). 

# Example:
# results2 <- read.table("outputs/coupling_score.tsv", header=TRUE, row.names=1, sep="\t", check.names=FALSE)
# obj <- build_celltype_treeplot(results2, out_path="outputs/celltype_tree.png")
# print(obj$plot)

