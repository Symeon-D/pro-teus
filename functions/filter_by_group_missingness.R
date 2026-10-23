#' Filter proteins by missingness within groups
#'
#' This function filters rows (proteins/genes) based on missing-value thresholds
#' applied separately within each experimental group. A protein is kept if it
#' has at least a user-defined proportion of non-missing values in **any** group.
#'
#' @param input_matrix A numeric matrix of proteomics intensities (rows = proteins, columns = samples).
#' @param sample_info A data.frame containing metadata for each sample. Must include the grouping column.
#' @param threshold Numeric value between 0 and 1. Proportion of samples within a group
#'   that must have non-missing values for a protein to be retained. Default is 0.7.
#' @param group_column Character string giving the column in `sample_info` containing the group labels.
#'
#' @return A filtered numeric matrix containing only proteins that meet the missingness requirement
#'   in at least one group.
#'
#' @examples
#' # filtered <- filter_by_group_missingness(matrix_data, sample_info, threshold = 0.7)
#'
#' @export
filter_by_group_missingness <- function(
    input_matrix,
    sample_info,
    threshold = 0.7,
    group_column = "Group"
) {
  # Extract groups
  groups <- sample_info[[group_column]]
  group_levels <- unique(groups)

  # Threshold counts per group
  group_counts <- table(groups)
  threshold_counts <- round(threshold * group_counts)

  # Logical matrix: rows x groups
  keep_matrix <- sapply(group_levels, function(g) {
    group_samples <- groups == g
    # Count non-NA values per row for the samples of this group
    rowSums(!is.na(input_matrix[, group_samples, drop = FALSE])) >= threshold_counts[g]
  })

  # Keep rows that pass at least one group
  keep_genes <- rowSums(keep_matrix) > 0

  # Filter matrix
  filtered_data <- input_matrix[keep_genes, , drop = FALSE]

  cat("✅ Kept", sum(keep_genes), "genes out of", nrow(input_matrix), "\n")

  return(filtered_data)
}
