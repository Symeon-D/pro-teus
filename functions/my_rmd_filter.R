#' Robust Mahalanobis Distance (rMd) sample filter - pmartR wrapper, adapted for the Pro-teus app
#'
#' The original my_rmd_filter() saved PNG / PDF files and made a PCA. This version
#' returns everything as R objects instead (no files, no PDF, no PCA: the app has its own PCA).
#' It is split in two so the app can compute the rMd ONCE and then re-flag instantly when the
#' p-value cutoff changes:
#'
#'   rmd_compute()        heavy step (pmartR::rmd_filter); run once per data set / metric choice
#'   rmd_flagged()        samples flagged at a cutoff       (pmartR: summary(filter, pvalue_threshold)$filtered_samples)
#'   rmd_make_plot()      labelled rMd plot at a cutoff
#'   rmd_sample_plot()    per-sample boxplots of every metric, the sample marked with an X
#'                        (pmartR: plot(filter, sampleID = s))
#'   my_rmd_filter()      convenience wrapper: compute + flag + all plots at one cutoff
#'
#' Metrics: MAD, Kurtosis, Skewness, Correlation, Proportion_Missing (2-5 of them).
#' Proportion_Missing is uninformative once the data are imputed (all zero).

rmd_compute <- function(Data, sample_info,
                        group_col = "Group", sample_col = "Sample",
                        metrics = NULL) {
  if (!requireNamespace("pmartR", quietly = TRUE))
    stop("pmartR is required for the rMd filter (install with BiocManager / install.packages('pmartR')).")
  if (!group_col %in% colnames(sample_info)) stop("group_col '", group_col, "' not found in sample_info.")
  if (!sample_col %in% colnames(sample_info)) stop("sample_col '", sample_col, "' not found in sample_info.")
  miss <- setdiff(sample_info[[sample_col]], colnames(Data))
  if (length(miss)) stop("These samples are in sample_info but not in Data: ", paste(miss, collapse = ", "))

  valid <- c("MAD", "Kurtosis", "Skewness", "Correlation", "Proportion_Missing")
  if (!is.null(metrics)) {
    bad <- setdiff(metrics, valid)
    if (length(bad)) stop("Invalid metrics: ", paste(bad, collapse = ", "))
    if (length(metrics) < 2) stop("At least 2 metrics are required.")
  }

  e_data <- as.data.frame(Data[, sample_info[[sample_col]], drop = FALSE])
  e_data <- cbind(ProteinID = rownames(e_data), e_data)
  rownames(e_data) <- NULL

  f_data <- sample_info
  f_data$SampleID <- f_data[[sample_col]]
  grp <- "pmart_grp_internal_"
  f_data[[grp]] <- f_data[[group_col]]
  colnames(f_data)[colnames(f_data) == "Group"] <- "group_original_"   # avoid pmartR's reserved name

  pro_obj <- suppressWarnings(pmartR::as.proData(
    e_data = e_data, f_data = f_data, edata_cname = "ProteinID", fdata_cname = "SampleID",
    data_scale = "log2"))
  pro_obj <- pmartR::group_designation(omicsData = pro_obj, main_effects = grp)

  rmd_filt <- pmartR::rmd_filter(omicsData = pro_obj, ignore_singleton_groups = TRUE, metrics = metrics)

  df <- as.data.frame(rmd_filt)
  colnames(df)[1] <- "Sample"
  df[[grp]] <- NULL
  df$Group <- sample_info[[group_col]][match(df$Sample, sample_info[[sample_col]])]
  df <- df[, c("Sample", "Group", setdiff(colnames(df), c("Sample", "Group"))), drop = FALSE]
  list(rmd_df = df, rmd_filt = rmd_filt, grp_col = grp, metrics = metrics)
}

# samples flagged at p < p_cutoff (pmartR's own summary(), with a plain p-value fallback)
rmd_flagged <- function(res, p_cutoff = 0.05) {
  fl <- tryCatch({
    s <- NULL
    utils::capture.output(s <- suppressMessages(summary(res$rmd_filt, pvalue_threshold = p_cutoff)))
    as.character(s$filtered_samples)
  }, error = function(e) NULL)
  if (is.null(fl)) fl <- res$rmd_df$Sample[!is.na(res$rmd_df$pvalue) & res$rmd_df$pvalue < p_cutoff]
  intersect(res$rmd_df$Sample, fl)
}

# table of the flagged samples with all their metrics
rmd_flag_table <- function(res, p_cutoff = 0.05) {
  d <- res$rmd_df[res$rmd_df$Sample %in% rmd_flagged(res, p_cutoff), , drop = FALSE]
  d[order(d$pvalue), , drop = FALSE]
}

# rMd (log2) per sample, flagged samples labelled
rmd_make_plot <- function(res, p_cutoff = 0.05, group_colors = NULL) {
  g <- res$grp_col
  p <- suppressWarnings(plot(res$rmd_filt, order_by = g, color_by = g,
                             pvalue_threshold = p_cutoff, bw_theme = TRUE))
  if (!is.null(group_colors))
    p <- p + ggplot2::scale_color_manual(values = group_colors) +
             ggplot2::scale_fill_manual(values = group_colors)
  fl <- rmd_flagged(res, p_cutoff)
  if (length(fl) && requireNamespace("ggrepel", quietly = TRUE)) {
    ld <- p$data[p$data$SampleID %in% fl, , drop = FALSE]
    ld <- ld[!duplicated(ld$SampleID), , drop = FALSE]
    xv <- rlang::get_expr(p$mapping$x)
    p <- p + ggrepel::geom_label_repel(data = ld, ggplot2::aes(x = !!xv, y = Log2.md, label = SampleID),
                                       size = 3, box.padding = 0.4, max.overlaps = 20,
                                       show.legend = FALSE, inherit.aes = FALSE)
  }
  p + ggplot2::labs(colour = "Group", fill = "Group", shape = "Group",
                    title = "Robust Mahalanobis distance per sample",
                    subtitle = paste0("Samples with p < ", p_cutoff, " are labelled (", length(fl), ")"))
}

# boxplots of every metric over all samples, the chosen sample marked with an X
rmd_sample_plot <- function(res, sample, p_cutoff = 0.05) {
  p <- suppressWarnings(plot(res$rmd_filt, sampleID = sample))
  r <- res$rmd_df[res$rmd_df$Sample == sample, , drop = FALSE]
  p + ggplot2::labs(title = paste0("Sample: ", sample),
                    subtitle = paste0("Group: ", r$Group[1], "  |  rMd p-value: ", signif(r$pvalue[1], 3),
                                      "  |  flagged at p < ", p_cutoff)) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", size = 16),
                   plot.subtitle = ggplot2::element_text(size = 12, colour = "firebrick"),
                   axis.text = ggplot2::element_text(size = 12), axis.title = ggplot2::element_text(size = 13),
                   strip.text = ggplot2::element_text(size = 12), legend.text = ggplot2::element_text(size = 12),
                   legend.title = ggplot2::element_text(size = 12))
}

# Convenience wrapper (same name as the original): everything at one cutoff, nothing written to disk
my_rmd_filter <- function(Data, sample_info, group_col = "Group", sample_col = "Sample",
                          metrics = NULL, p_cutoff = 0.05, group_colors = NULL) {
  res <- rmd_compute(Data, sample_info, group_col, sample_col, metrics)
  fl <- rmd_flagged(res, p_cutoff)
  list(rmd_results = transform(res$rmd_df, outlier = res$rmd_df$Sample %in% fl),
       flagged_samples = fl,
       rmdFilt_object = res$rmd_filt,
       rmd_plot = tryCatch(rmd_make_plot(res, p_cutoff, group_colors), error = function(e) NULL),
       metric_plots = setNames(lapply(fl, function(s) tryCatch(rmd_sample_plot(res, s, p_cutoff), error = function(e) NULL)), fl),
       computed = res)
}
