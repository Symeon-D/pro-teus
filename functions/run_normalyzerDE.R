#' Run NormalyzerDE from an in-memory expression matrix and sample table
#'
#' Converts a numeric expression matrix (proteins x samples) into the TSV formats
#' required by `NormalyzerDE::normalyzer()`, cleans/filters the design table to
#' the required `sample` and `group` columns, writes the temporary (or persistent)
#' input files, and runs NormalyzerDE.
#'
#' @param expr_matrix A numeric matrix of expression values with proteins in rows
#'   and samples in columns. `rownames(expr_matrix)` must be the protein IDs.
#' @param sample_table A `data.frame` containing at least a sample identifier
#'   column and a grouping column, plus any additional metadata columns (allowed).
#' @param noLogTransform Don't log-transform the input.
#' @param sample_col Name of the sample column in `sample_table`. If `NULL`,
#'   the function will try to auto-detect (case-insensitive) one of:
#'   `sample`, `Sample`, `samples`, `Samples`.
#' @param group_col Name of the group column in `sample_table`. If `NULL`,
#'   the function will try to auto-detect (case-insensitive) one of:
#'   `group`, `Group`, `condition`, `Condition`.
#' @param jobName A label passed to `NormalyzerDE::normalyzer(jobName=...)`.
#' @param outputDir Directory passed to `NormalyzerDE::normalyzer(outputDir=...)`.
#'   Defaults to `getwd()`.
#' @param keep_input_files If `TRUE`, the cleaned `*.tsv` inputs are written into
#'   `outputDir` and kept. If `FALSE`, inputs are written to `tempdir()`.
#' @param file_prefix Prefix for the written input TSVs.
#'
#' @return A list with paths to the generated TSV files and the return value of
#'   `NormalyzerDE::normalyzer()`.
#'
#' @export
run_normalyzerDE <- function(
    expr_matrix,
    sample_table,
    noLogTransform = TRUE,
    sample_col = NULL,
    group_col  = NULL,
    jobName = "Proteomics_Normalization_Check",
    outputDir = getwd(),
    keep_input_files = TRUE,
    file_prefix = "normalyzer") {

  # ---- dependencies ----
  if (!requireNamespace("NormalyzerDE", quietly = TRUE)) {
    stop("Package `NormalyzerDE` is required but not installed.")
  }

  # ---- validate expr_matrix ----
  if (!(is.matrix(expr_matrix) || is.data.frame(expr_matrix))) {
    stop("`expr_matrix` must be a matrix or data.frame.")
  }
  expr_matrix <- as.matrix(expr_matrix)

  rn <- rownames(expr_matrix)
  cn <- colnames(expr_matrix)

  storage.mode(expr_matrix) <- "numeric"     # in-place coercion, keeps dimnames
  expr_matrix[is.nan(expr_matrix)] <- NA_real_

  rownames(expr_matrix) <- rn
  colnames(expr_matrix) <- cn

  protein_ids <- rownames(expr_matrix)

  # Build NormalyzerDE data table: first column MUST be ProteinID
  data_df <- data.frame(
    ProteinID = protein_ids,
    expr_matrix,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # ---- validate and clean sample_table ----
  if (!is.data.frame(sample_table)) {
    stop("`sample_table` must be a data.frame.")
  }

  st <- sample_table
  st_names_lower <- tolower(colnames(st))

  auto_pick_col <- function(targets_lower) {
    hit <- which(st_names_lower %in% targets_lower)
    if (length(hit) >= 1) colnames(st)[hit[1]] else NULL
  }

  if (is.null(sample_col)) {
    sample_col <- auto_pick_col(c("sample", "samples"))
  }
  if (is.null(group_col)) {
    group_col <- auto_pick_col(c("group", "condition"))
  }

  if (is.null(sample_col) || !sample_col %in% colnames(st)) {
    stop("Could not identify `sample_col`. Please provide `sample_col` explicitly.")
  }
  if (is.null(group_col) || !group_col %in% colnames(st)) {
    stop("Could not identify `group_col`. Please provide `group_col` explicitly.")
  }

  design_df <- st[, c(sample_col, group_col), drop = FALSE]
  colnames(design_df) <- c("sample", "group")
  design_df$sample <- as.character(design_df$sample)
  design_df$group  <- as.character(design_df$group)

  # Filter design to samples that exist in the data
  design_df <- design_df[design_df$sample %in% colnames(expr_matrix), , drop = FALSE]

  # Ensure every data column has a design row; reorder to match data columns
  m <- match(colnames(expr_matrix), design_df$sample)
  if (anyNA(m)) {
    missing_samples <- colnames(expr_matrix)[is.na(m)]
    stop(
      "The following samples exist in `expr_matrix` but not in `sample_table` after filtering: ",
      paste(missing_samples, collapse = ", ")
    )
  }
  design_df <- design_df[m, , drop = FALSE]

  # ---- write input files ----
  if (!dir.exists(outputDir)) dir.create(outputDir, recursive = TRUE, showWarnings = FALSE)

  write_dir <- if (keep_input_files) outputDir else tempdir()
  data_path   <- file.path(write_dir, paste0(file_prefix, "_data_clean.tsv"))
  design_path <- file.path(write_dir, paste0(file_prefix, "_design_clean.tsv"))

  utils::write.table(data_df,   data_path,   sep = "\t", row.names = FALSE, quote = FALSE)
  utils::write.table(design_df, design_path, sep = "\t", row.names = FALSE, quote = FALSE)

  # ---- run NormalyzerDE ----
  res <- NormalyzerDE::normalyzer(
    jobName    = jobName,
    sampleColName = "sample",
    groupColName = "group",
    noLogTransform = noLogTransform,
    designPath = design_path,
    dataPath   = data_path,
    outputDir  = outputDir
  )

  invisible(list(
    data_path   = data_path,
    design_path = design_path,
    outputDir   = outputDir,
    jobName     = jobName,
    result      = res
  ))
}
