
calculate_diversity <- function(cell_counts) {
  # Remove zeros
  counts <- cell_counts[cell_counts > 0]
  
  if (length(counts) == 0) {
    return(list(
      gini_index = NA,
      clone_count = 0,
      total_cells = 0,
      clonal_density = NA,
      shannon_entropy = 0,
      simpson_index = 0
    ))
  }
  
  # Calculate total cells and proportions
  total_cells <- sum(counts)
  props <- counts / total_cells
  
  # Clone count
  clone_count <- length(counts)
  
  # Clonal density (clones per cell)
  clonal_density <- clone_count / total_cells
  
  # Gini index
  sorted_props <- sort(props)
  n <- length(sorted_props)
  gini_index <- (2 * sum((1:n) * sorted_props)) / (n * sum(sorted_props)) - (n + 1) / n
  
  # Shannon entropy
  shannon_entropy <- -sum(props * log(props))
  
  # Simpson index (1 - Simpson's D)
  simpson_index <- 1 - sum(props^2)
  
  return(list(
    gini_index = gini_index,
    clone_count = clone_count,
    total_cells = total_cells,
    clonal_density = clonal_density,
    shannon_entropy = shannon_entropy,
    simpson_index = simpson_index
  ))
}
