#' QC: assess whether global normalization is warranted (median shift, spread shift, missingness coupling)
#'
#' Computes three quick diagnostics that are useful for deciding whether to normalize a
#' proteomics expression matrix (proteins x samples):
#'
#' 1) `median_range`: global location shift across samples (max median - min median)  
#' 2) `cor_median_missing`: coupling between sample medians and missingness fraction  
#' 3) `iqr_range`: global spread shift across samples (max IQR - min IQR)
#'
#' Interpretation logic (default thresholds, log-scale):
#' - Larger `median_range` suggests global sample shifts (loading/run effects).  
#' - Large `|cor_median_missing|` suggests sample-level bias is tied to data completeness (quality effect).  
#' - Larger `iqr_range` suggests differing distribution spreads (variance heterogeneity), often favoring
#'   variance-stabilizing normalization (e.g., VSN) over simple centering.
#'
#' Notes:
#' - These are heuristics, not proofs. If your biology plausibly induces global shifts (rare, but possible),
#'   interpret cautiously and corroborate with RLE/MA/correlation plots.
#'
#' @param x Numeric matrix/data.frame of expression values with proteins in rows and samples in columns.
#'   Missing values should be `NA`. Non-numeric entries are coerced to numeric (may become `NA`).
#' @param log_base Log base for interpreting location/spread units. Default `2` (log2 units).
#' @param cor_method Correlation method: `"pearson"` (default), `"spearman"`, or `"kendall"`.
#' @param thresholds List of thresholds used for qualitative interpretation.
#'   Defaults are reasonable for log2-scale proteomics:
#'   `median = c(low=0.2, moderate=0.5)`,
#'   `cor = c(moderate=0.4, strong=0.6)`,
#'   `iqr = c(low=0.1, moderate=0.25)`.
#' @param verbose If `TRUE`, prints an interpretation summary to the console.
#'
#' @return A list with:
#' - `sample_medians`, `missing_frac`, `sample_iqr` (named vectors)
#' - `median_range`, `iqr_range`, `cor_median_missing`
#' - `median_fold_change_equivalent` (computed as `log_base^median_range`)
#' - `decision` (a concise recommendation string)
#' - `interpretation` (a longer narrative string you can paste into methods/QC notes)
#'
#' @export
#'
#' @examples
#' res <- qc_normalization_need(filtered_data)
#' res$median_range
#' cat(res$interpretation)
qc_normalization_need <- function(
    x,
    log_base = 2,
    cor_method = c("pearson", "spearman", "kendall"),
    thresholds = list(
      median = c(low = 0.2, moderate = 0.5),
      cor    = c(moderate = 0.4, strong = 0.6),
      iqr    = c(low = 0.1, moderate = 0.25)
    ),
    verbose = TRUE
) {
  cor_method <- match.arg(cor_method)
  
  if (!(is.matrix(x) || is.data.frame(x))) {
    stop("`x` must be a matrix or data.frame (proteins x samples).")
  }
  
  x <- as.matrix(x)
  
  # Coerce to numeric while preserving dimnames
  rn <- rownames(x)
  cn <- colnames(x)
  x <- suppressWarnings(matrix(
    as.numeric(x),
    nrow = nrow(x),
    ncol = ncol(x),
    dimnames = list(rn, cn)
  ))
  
  if (is.null(colnames(x)) || any(colnames(x) == "")) {
    stop("`x` must have non-empty column names (sample names).")
  }
  
  # Core metrics
  sample_medians <- apply(x, 2, stats::median, na.rm = TRUE)
  median_range   <- diff(range(sample_medians, na.rm = TRUE))
  
  missing_frac <- colMeans(is.na(x))
  
  cor_median_missing <- suppressWarnings(stats::cor(
    sample_medians, missing_frac,
    use = "complete.obs",
    method = cor_method
  ))
  
  sample_iqr <- apply(x, 2, stats::IQR, na.rm = TRUE)
  iqr_range  <- diff(range(sample_iqr, na.rm = TRUE))
  
  # Fold-change equivalent for the median range (only meaningful if x is log-scale)
  median_fc_equiv <- log_base^median_range
  
  # Qualitative labels
  median_label <- if (median_range <= thresholds$median["low"]) {
    "small"
  } else if (median_range <= thresholds$median["moderate"]) {
    "moderate"
  } else {
    "large"
  }
  
  cor_abs <- abs(cor_median_missing)
  cor_label <- if (!is.finite(cor_abs)) {
    "undefined"
  } else if (cor_abs < thresholds$cor["moderate"]) {
    "weak"
  } else if (cor_abs < thresholds$cor["strong"]) {
    "moderate"
  } else {
    "strong"
  }
  
  iqr_label <- if (iqr_range <= thresholds$iqr["low"]) {
    "small"
  } else if (iqr_range <= thresholds$iqr["moderate"]) {
    "moderate"
  } else {
    "large"
  }
  
  # Decision logic (teachable and explicit)
  score <- 0L
  if (median_label == "large") score <- score + 2L
  if (median_label == "moderate") score <- score + 1L
  if (cor_label == "strong") score <- score + 2L
  if (cor_label == "moderate") score <- score + 1L
  if (iqr_label == "large") score <- score + 2L
  if (iqr_label == "moderate") score <- score + 1L
  
  decision <- if (score >= 4L) {
    "Normalization is strongly recommended (multiple indicators of structured technical differences)."
  } else if (score >= 2L) {
    "Normalization is recommended (at least one meaningful indicator of technical shift/heterogeneity)."
  } else {
    "Normalization may be optional (indicators are small/weak), but confirm with RLE/MA/correlation plots."
  }
  
  interpretation <- paste0(
    "QC normalization heuristics (proteins x samples):\n",
    "- `median_range` = ", signif(median_range, 4), " (", median_label, " global location shift). ",
    "On log-base-", log_base, " scale, this corresponds to ~", signif(median_fc_equiv, 4),
    "× difference between the highest- and lowest-median samples.\n",
    "- `cor(median, missingness)` = ", signif(cor_median_missing, 4), " (", cor_label,
    " coupling). A ", cor_label, " absolute correlation suggests sample-level intensity shifts are tied ",
    "to completeness/quality, implying structured technical bias rather than random scatter.\n",
    "- `iqr_range` = ", signif(iqr_range, 4), " (", iqr_label, " spread/scale differences). Larger IQR ",
    "differences indicate heterogeneous variability across samples; variance-stabilizing normalization ",
    "(e.g., VSN) is often more appropriate than simple median-centering.\n",
    "Decision: ", decision
  )
  
  if (isTRUE(verbose)) {
    cat(interpretation, "\n")
  }
  
  invisible(list(
    sample_medians = sample_medians,
    missing_frac = missing_frac,
    sample_iqr = sample_iqr,
    median_range = median_range,
    iqr_range = iqr_range,
    cor_median_missing = cor_median_missing,
    median_fold_change_equivalent = median_fc_equiv,
    decision = decision,
    interpretation = interpretation
  ))
}
