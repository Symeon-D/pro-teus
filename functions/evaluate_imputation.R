suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(tidyr) })

#' Evaluate MNAR Imputation Methods for Proteomics Data
#'
#' @description
#' Evaluates and compares multiple missing value imputation methods suited for
#' proteomics data, with a focus on Missing Not At Random (MNAR) and global
#' structure approaches. Generates nine diagnostic plots to guide method selection:
#' missing value pattern, density distributions, PCA, NRMSE, coefficient of
#' variation (within-group), intensity shift, NAguideR classic metrics
#' (NRMSE / SOR / ACC / PCC), and DE concordance (Jaccard similarity).
#'
#' @param data A numeric matrix or data frame with proteins as rows and samples
#'   as columns. Missing values should be encoded as \code{NA}.
#' @param sample_info A data frame with at least two columns identifying samples
#'   and their group membership.
#' @param sample_col Character. Column name in \code{sample_info} that matches
#'   \code{colnames(data)}. Default: \code{"Sample"}.
#' @param group_col Character. Column name in \code{sample_info} that defines
#'   experimental groups. Default: \code{"Group"}.
#' @param group_colors Named character vector of colours for each group. If
#'   \code{NULL} (default), colours are assigned automatically.
#' @param methods Character vector of imputation methods to evaluate. Accepts
#'   any method supported by \code{my_impute()}, or one of the shorthand strings
#'   \code{"all"}, \code{"mnar"}, \code{"global"}, \code{"mar"}.
#'   Default: \code{c("mindet", "minprob", "qr", "perseus", "bpca", "svd", "knn", "lls", "rf")}.
#' @param mask_fraction Numeric scalar in (0, 1). Fraction of observed values
#'   to artificially mask for NRMSE / SOR / ACC calculation. Default: \code{0.10}.
#' @param log2_transform Logical. Apply log2 transformation before imputation.
#'   Set \code{FALSE} if data are already log2-transformed. Default: \code{TRUE}.
#' @param run_de Logical. If \code{TRUE}, run \code{limma}-based DE analysis on
#'   each imputed matrix and compute Jaccard similarity of significant hits
#'   (FDR < 0.05). Requires \code{limma}. Default: \code{TRUE}.
#' @param fdr_cutoff Numeric. FDR threshold for the DE concordance plot.
#'   Default: \code{0.05}.
#' @param seed Integer. Random seed for reproducibility. Default: \code{42}.
#' @param output_dir Character. Directory for PDF output. \code{NULL} to skip
#'   saving. Default: \code{"imputation_plots"}.
#' @param pdf_name Character. Name of the output PDF file.
#'   Default: \code{"imputation_evaluation.pdf"}.
#'
#' @return Invisibly returns a named list:
#'   \describe{
#'     \item{\code{imputed}}{Named list of imputed matrices, one per method.}
#'     \item{\code{plot_missing}}{Heatmap of missing value pattern.}
#'     \item{\code{plot_density}}{Density plot of observed vs imputed values.}
#'     \item{\code{plot_pca}}{PCA plots before and after imputation.}
#'     \item{\code{plot_nrmse}}{Global NRMSE bar chart (random mask).}
#'     \item{\code{plot_cv}}{Within-group CV violin plot.}
#'     \item{\code{plot_downshift}}{Downshift from observed mean bar chart.}
#'     \item{\code{plot_classic}}{NAguideR metrics: NRMSE / SOR / ACC / PCC.}
#'     \item{\code{plot_de_concordance}}{Jaccard heatmap of DE hits (if \code{run_de = TRUE}).}
#'   }
#'
#' @details
#' \subsection{Interpreting the plots}{
#'   \enumerate{
#'     \item \strong{Missing pattern}: Check whether missingness is concentrated
#'       in specific samples or conditions.
#'     \item \strong{Density}: MNAR methods should produce a left-shifted peak
#'       below the observed distribution.
#'     \item \strong{PCA}: Verify that group separation is preserved after imputation.
#'     \item \strong{NRMSE (random mask)}: Compare \emph{within} category only —
#'       MNAR methods are penalised by random masks by design.
#'     \item \strong{Within-group CV}: Flag any method that inflates variability
#'       relative to the complete-case baseline.
#'     \item \strong{Downshift}: Confirms MNAR methods impute at the low end.
#'     \item \strong{NAguideR metrics}: Per-feature NRMSE, SOR, ACC, and
#'       inter-sample PCC heatmap.
#'     \item \strong{DE concordance}: High Jaccard among MNAR methods = robust,
#'       stable differential expression results.
#'   }
#' }
#'
#' @references
#' Lazar C, et al. (2016). J Proteome Res, 15(4):1116-25.
#' Wang S, et al. (2020). NAguideR. Nucleic Acids Res, 48(14):e83.
#' Jin L, et al. (2021). Sci Rep, 11(1):1760.
#'
#' @export
evaluate_imputation <- function(data,
                                sample_info,
                                sample_col     = "Sample",
                                group_col      = "Group",
                                group_colors   = NULL,
                                methods        = c("mindet", "minprob", "qr",
                                                   "perseus", "bpca", "svd",
                                                   "knn", "lls", "rf"),
                                mask_fraction  = 0.10,
                                log2_transform = TRUE,
                                run_de         = TRUE,
                                fdr_cutoff     = 0.05,
                                seed           = 42,
                                output_dir     = "imputation_plots",
                                pdf_name       = "imputation_evaluation.pdf",
                                impute_args    = list()) {   # e.g. list(knn = list(k = 5))

  set.seed(seed)

  # ── resolve method shorthands ──────────────────────────────────────────────
  all_methods <- c(
    "mindet", "minprob", "qr", "perseus", "zero", "minimum",
    "bpca", "svd", "mle", "colmedian", "rowmedian",
    "knn", "seqknn", "trknn", "lls", "irmi", "grr", "gms",
    "mice_norm", "mice_cart", "rf", "impseq", "impseqrob"
  )
  if (length(methods) == 1 && methods == "all") {
    methods <- all_methods
    cat("Running ALL", length(methods), "imputation methods.\n")
  } else if (length(methods) == 1 && methods == "mnar") {
    methods <- c("mindet", "minprob", "qr", "perseus", "zero", "minimum")
    cat("Running MNAR methods only:", paste(methods, collapse = ", "), "\n")
  } else if (length(methods) == 1 && methods == "global") {
    methods <- c("bpca", "svd", "mle", "colmedian", "rowmedian")
    cat("Running Global Structure methods only:", paste(methods, collapse = ", "), "\n")
  } else if (length(methods) == 1 && methods == "mar") {
    methods <- c("knn", "seqknn", "trknn", "lls", "irmi", "grr", "gms",
                 "mice_norm", "mice_cart", "rf", "impseq", "impseqrob")
    cat("Running MAR/Local methods only:", paste(methods, collapse = ", "), "\n")
  }

  # ── colour palettes ────────────────────────────────────────────────────────
  mnar_methods   <- c("mindet", "minprob", "qr", "perseus", "zero", "minimum")
  global_methods <- c("bpca", "svd", "mle", "colmedian", "rowmedian")
  mar_methods    <- c("knn", "seqknn", "trknn", "lls", "irmi", "grr", "gms",
                      "mice_norm", "mice_cart", "rf", "impseq", "impseqrob")

  method_type <- function(m) {
    dplyr::case_when(
      m %in% mnar_methods   ~ "MNAR",
      m %in% global_methods ~ "Global Structure",
      m %in% mar_methods    ~ "MAR / Local",
      TRUE                  ~ "Other"
    )
  }

  type_colors <- c(
    "MNAR"             = "#E07B54",
    "Global Structure" = "#5B8DB8",
    "MAR / Local"      = "#6AAB6E",
    "Other"            = "grey60",
    "observed"         = "#3A3A3A"
  )

  nrmse_shift_colors <- c(
    "MNAR"             = "#9B6FA0",
    "Global Structure" = "#3A9E8F",
    "MAR / Local"      = "#B8A830",
    "Other"            = "grey60"
  )

  # ── force matrix + log2 ────────────────────────────────────────────────────
  data <- as.matrix(data)
  if (log2_transform) {
    data[data <= 0] <- NA # Prevent -Inf or NaN from log2 of zeros/negatives
    data <- log2(data)
    data[is.infinite(data)] <- NA
    cat("Applied log2 transform.\n")
  } else {
    cat("Skipping log2 transform (data assumed already log2-transformed).\n")
  }

  # ── output directory ───────────────────────────────────────────────────────
  if (!is.null(output_dir) && !dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
    cat("Created output directory:", output_dir, "\n")
  }

  # ── sample_info helpers ────────────────────────────────────────────────────
  groups <- sample_info[[group_col]][match(colnames(data),
                                           sample_info[[sample_col]])]
  unique_groups <- unique(na.omit(groups))
  if (is.null(group_colors)) {
    pal          <- scales::hue_pal()(length(unique_groups))
    group_colors <- setNames(pal, unique_groups)
  }

  # ── run_imputation wrapper ─────────────────────────────────────────────────
  run_imputation <- function(mat, method) {
    do.call(my_impute, c(list(input_matrix = mat, method = method, seed = seed, verbose = FALSE),
                         impute_args[[method]]))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 1 — Missing value pattern
  # ══════════════════════════════════════════════════════════════════════════
  plot_missing <- function(mat) {

    na_mat    <- is.na(mat)
    row_order <- order(rowSums(na_mat), decreasing = TRUE)
    melt_na   <- reshape2::melt(na_mat[row_order, ])
    colnames(melt_na) <- c("Protein", "Sample", "Missing")

    p_heat <- ggplot(melt_na, aes(x = Sample, y = Protein, fill = Missing)) +
      geom_raster() +
      scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "tomato"),
                        labels = c("Observed", "Missing")) +
      labs(title = "Missing Value Pattern", x = "", y = "Proteins", fill = "") +
      theme_minimal(base_size = 11) +
      theme(axis.text.y      = element_blank(),
            axis.ticks.y     = element_blank(),
            axis.text.x      = element_text(angle = 45, hjust = 1),
            legend.position  = "top",
            plot.margin      = margin(5, 10, 5, 5))

    na_pct        <- data.frame(Sample = colnames(mat),
                                NA_pct = colMeans(is.na(mat)) * 100)
    na_pct$Sample <- factor(na_pct$Sample, levels = na_pct$Sample)

    p_bar <- ggplot(na_pct, aes(x = Sample, y = NA_pct, fill = NA_pct)) +
      geom_col() +
      scale_fill_gradient(low = "steelblue", high = "tomato") +
      labs(title = "% Missing per Sample", x = "", y = "% Missing") +
      theme_bw(base_size = 11) +
      theme(axis.text.x      = element_text(angle = 45, hjust = 1),
            legend.position  = "none",
            plot.margin      = margin(5, 5, 5, 10)) +
      coord_flip()

    cowplot::plot_grid(p_heat, p_bar, ncol = 2, rel_widths = c(1, 1))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 2 — Density: faceted, one panel per method, observed overlaid
  # ══════════════════════════════════════════════════════════════════════════
  plot_densities <- function(mat, imputed_list) {

    obs_vals <- as.vector(mat[!is.na(mat)])

    df <- do.call(rbind, lapply(names(imputed_list), function(m) {
      rbind(
        data.frame(value  = obs_vals,
                   Series = "Observed",
                   Type   = "observed",
                   Method = m),
        data.frame(value  = as.vector(imputed_list[[m]][is.na(mat)]),
                   Series = m,
                   Type   = method_type(m),
                   Method = m)
      )
    }))

    series_levels <- c("Observed", names(imputed_list))
    col_vec <- c(
      "Observed" = type_colors["observed"][[1]],
      setNames(type_colors[method_type(names(imputed_list))], names(imputed_list))
    )
    lty_vec <- c(
      "Observed" = "solid",
      setNames(rep("dashed", length(imputed_list)), names(imputed_list))
    )
    alpha_vec <- c(
      "Observed" = 0.15,
      setNames(rep(0.30, length(imputed_list)), names(imputed_list))
    )

    df$Series <- factor(df$Series, levels = series_levels)

    ggplot(df, aes(x = value, fill = Series, colour = Series,
                   linetype = Series, alpha = Series)) +
      geom_density(linewidth = 0.9) +
      facet_wrap(~ Method, ncol = 3, scales = "free_y") +
      scale_fill_manual(values     = col_vec) +
      scale_colour_manual(values   = col_vec) +
      scale_linetype_manual(values = lty_vec) +
      scale_alpha_manual(values    = alpha_vec) +
      labs(title    = "Density: Observed vs Imputed Values per Method",
           subtitle = "Black solid = observed. Coloured dashed = imputed (by method type).\nMNAR methods should show a left-tail peak.",
           x = "log2 Intensity", y = "Density",
           fill = "", colour = "", linetype = "", alpha = "") +
      theme_bw(base_size = 11) +
      theme(legend.position = "bottom") +
      guides(fill     = guide_legend(nrow = 1),
             colour   = guide_legend(nrow = 1),
             linetype = guide_legend(nrow = 1),
             alpha    = "none")
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 3 — PCA: shared legend, one panel per method
  # ══════════════════════════════════════════════════════════════════════════
  plot_pca <- function(mat, imputed_list) {

    do_pca <- function(m, label) {
      m <- as.matrix(m)
      m <- m[complete.cases(m), , drop = FALSE] # Drop incomplete rows instead of min-filling
      if (nrow(m) < 3) return(NULL) # Skip if too few proteins remain

      pca_res <- prcomp(t(m), scale. = TRUE, center = TRUE)
      data.frame(
        PC1    = pca_res$x[, 1],
        PC2    = pca_res$x[, 2],
        Sample = colnames(m),
        Group  = groups,
        Panel  = label
      )
    }

    pca_dfs <- c(
      list(do_pca(mat, "Before imputation")),
      lapply(names(imputed_list), function(m) do_pca(imputed_list[[m]], m))
    )
    pca_dfs <- Filter(Negate(is.null), pca_dfs) # Filter out failed PCAs

    if (length(pca_dfs) == 0) return(NULL)

    make_panel <- function(df_sub, label) {
      ggplot(df_sub, aes(PC1, PC2, colour = Group)) +
        geom_point(size = 3) +
        ggrepel::geom_text_repel(aes(label = Sample), size = 2,
                                 show.legend = FALSE, max.overlaps = 15) +
        scale_color_manual(values = group_colors) +
        labs(title = label, x = "PC1", y = "PC2", colour = group_col) +
        theme_bw(base_size = 10) +
        theme(legend.position = "none")
    }

    panel_labels <- sapply(pca_dfs, function(df) df$Panel[1])
    plots        <- mapply(make_panel, pca_dfs, panel_labels, SIMPLIFY = FALSE)
    ncols        <- min(3, length(plots))

    legend_plot <- ggplot(pca_dfs[[1]], aes(PC1, PC2, colour = Group)) +
      geom_point(size = 4) +
      scale_color_manual(values = group_colors) +
      labs(colour = group_col) +
      theme_bw() +
      guides(colour = guide_legend(nrow = 1,
                                   override.aes = list(size = 4, shape = 16)))
    shared_legend <- cowplot::get_legend(legend_plot)

    cowplot::plot_grid(
      cowplot::ggdraw() +
        cowplot::draw_label("PCA per Imputation Method", fontface = "bold", size = 13),
      shared_legend,
      cowplot::plot_grid(plotlist = plots, ncol = ncols),
      ncol = 1,
      rel_heights = c(0.04, 0.06, 1)
    )
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 4 — NRMSE bar chart (global random mask)
  # ══════════════════════════════════════════════════════════════════════════
  plot_nrmse <- function(imp_list_masked, mask_idx, true_vals) {
    if (length(imp_list_masked) == 0) return(NULL)

    nrmse_df <- do.call(rbind, lapply(names(imp_list_masked), function(m) {
      pred  <- as.matrix(imp_list_masked[[m]])[mask_idx]
      nrmse <- sqrt(mean((pred - true_vals)^2, na.rm = TRUE)) /
        sd(true_vals, na.rm = TRUE)
      data.frame(Method = m, NRMSE = nrmse)
    }))

    nrmse_df        <- nrmse_df[!is.na(nrmse_df$NRMSE), ]
    nrmse_df        <- nrmse_df[order(nrmse_df$NRMSE), ]
    nrmse_df$Method <- factor(nrmse_df$Method, levels = nrmse_df$Method)
    nrmse_df$Type   <- method_type(as.character(nrmse_df$Method))

    ggplot(nrmse_df, aes(x = Method, y = NRMSE, fill = Type)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = round(NRMSE, 3)), vjust = -0.4, size = 3.5) +
      scale_fill_manual(values = nrmse_shift_colors) +
      labs(title    = paste0("NRMSE on artificially masked values (",
                             mask_fraction * 100, "%)"),
           subtitle = "Compare within category only — MNAR methods are penalised on random masks by design.",
           x = "", y = "NRMSE", fill = "Method type") +
      theme_bw(base_size = 11) +
      theme(axis.text.x      = element_text(angle = 30, hjust = 1),
            legend.position  = "top")
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 5 — Within-group CV
  # ══════════════════════════════════════════════════════════════════════════
  plot_cv <- function(mat, imputed_list) {

    get_within_cv <- function(m, grps) {
      unlist(lapply(unique(grps), function(g) {
        sub <- m[, grps == g, drop = FALSE]
        apply(sub, 1, function(x) {
          if (sum(!is.na(x)) < 2) return(NA_real_)
          mn <- mean(x, na.rm = TRUE)
          if (is.na(mn) || mn == 0) return(NA_real_)
          sd(x, na.rm = TRUE) / abs(mn) * 100
        })
      }))
    }

    cv_list       <- list()
    complete_rows <- mat[complete.cases(mat), , drop = FALSE]
    if (nrow(complete_rows) > 0)
      cv_list[["Complete\n(no imputation)"]] <- get_within_cv(complete_rows, groups)

    for (m in names(imputed_list))
      cv_list[[m]] <- get_within_cv(imputed_list[[m]], groups)

    df <- do.call(rbind, lapply(names(cv_list), function(nm)
      data.frame(CV = cv_list[[nm]], Method = nm, stringsAsFactors = FALSE)))
    df <- df[!is.na(df$CV) & df$CV > 0, ]

    method_order <- c("Complete\n(no imputation)", names(imputed_list))
    df$Method    <- factor(df$Method, levels = method_order)

    ggplot(df, aes(x = Method, y = CV, fill = Method)) +
      geom_violin(alpha = 0.55, draw_quantiles = c(0.25, 0.5, 0.75), trim = TRUE) +
      geom_boxplot(width = 0.08, outlier.shape = NA, fill = "white") +
      scale_y_log10(labels = scales::comma) +
      labs(title    = "Within-group CV per imputation method",
           subtitle = "Compare to 'Complete (no imputation)' baseline — good imputation should not inflate CV",
           x = "", y = "CV (%) — log10 scale") +
      theme_bw(base_size = 11) +
      theme(legend.position  = "none",
            axis.text.x      = element_text(angle = 30, hjust = 1))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 6 — Downshift bar
  # ══════════════════════════════════════════════════════════════════════════
  plot_downshift <- function(mat, imputed_list) {

    obs_vals <- as.vector(mat[!is.na(mat)])

    stats_df <- do.call(rbind, lapply(names(imputed_list), function(m) {
      imp_vals <- as.vector(imputed_list[[m]][is.na(mat)])
      data.frame(
        Method          = m,
        Shift_from_mean = mean(obs_vals) - mean(imp_vals),
        Type            = method_type(m)
      )
    }))

    ggplot(stats_df,
           aes(x = reorder(Method, Shift_from_mean),
               y = Shift_from_mean, fill = Type)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = round(Shift_from_mean, 2)),
                hjust = -0.2, size = 3.5) +
      coord_flip() +
      scale_fill_manual(values = nrmse_shift_colors) +
      labs(title    = "Downshift from observed mean by method",
           subtitle = "MNAR methods should show large downshift.\nGlobal structure / MAR should show near-zero shift.",
           x = "", y = "Mean(observed) - Mean(imputed)",
           fill = "Method type") +
      theme_bw(base_size = 12) +
      theme(legend.position = "top")
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 7 — NAguideR classic metrics: NRMSE / SOR / ACC / PCC heatmap
  # ══════════════════════════════════════════════════════════════════════════
  plot_classic_metrics <- function(mat, imputed_list, imp_list_masked, mask_idx, true_vals) {
    if (length(imp_list_masked) == 0) return(NULL)

    # per-feature NRMSE
    feature_nrmse <- do.call(cbind, lapply(names(imp_list_masked), function(m) {
      imp <- imp_list_masked[[m]]
      sapply(seq_len(nrow(mat)), function(i) {
        row_start <- (i - 1L) * ncol(mat)
        row_mask  <- mask_idx[mask_idx > row_start & mask_idx <= row_start + ncol(mat)]
        if (length(row_mask) == 0L) return(NA_real_)
        col_idx  <- row_mask - row_start
        true_v   <- mat[i, col_idx]
        pred_v   <- imp[i, col_idx]
        sd_v     <- sd(mat[i, ], na.rm = TRUE)
        if (is.na(sd_v) || sd_v == 0) return(NA_real_)
        sqrt(mean((pred_v - true_v)^2, na.rm = TRUE)) / sd_v
      })
    }))
    colnames(feature_nrmse) <- names(imp_list_masked)
    mean_nrmse <- colMeans(feature_nrmse, na.rm = TRUE)

    # SOR
    sor_vals <- sapply(names(imp_list_masked), function(m) {
      col_idx <- which(colnames(feature_nrmse) == m)
      ranks   <- apply(feature_nrmse, 1, function(r) {
        if (all(is.na(r))) return(NA_real_)
        rank(r, na.last = "keep")[col_idx]
      })
      sum(ranks, na.rm = TRUE)
    })

    # ACC
    acc_vals <- sapply(names(imp_list_masked), function(m) {
      pred <- as.matrix(imp_list_masked[[m]])[mask_idx]
      suppressWarnings(cor(true_vals, pred, use = "complete.obs",
                           method = "pearson"))
    })

    metrics_df <- data.frame(
      Method = names(imp_list_masked),
      NRMSE  = mean_nrmse,
      SOR    = sor_vals,
      ACC    = acc_vals,
      Type   = method_type(names(imp_list_masked))
    )
    metrics_df        <- metrics_df[order(metrics_df$NRMSE), ]
    metrics_df$Method <- factor(metrics_df$Method, levels = metrics_df$Method)

    p_nrmse <- ggplot(metrics_df, aes(x = Method, y = NRMSE, fill = Type)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = round(NRMSE, 3)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = nrmse_shift_colors) +
      labs(title = "a1) Mean per-feature NRMSE", subtitle = "Lower = better",
           x = "", y = "Mean NRMSE", fill = "") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "none")

    sor_df        <- metrics_df[order(metrics_df$SOR), ]
    sor_df$Method <- factor(sor_df$Method, levels = sor_df$Method)

    p_sor <- ggplot(sor_df, aes(x = Method, y = SOR, fill = Type)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = round(SOR, 0)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = nrmse_shift_colors) +
      labs(title = "a2) SOR (NRMSE-based Sum of Ranks)", subtitle = "Lower = better",
           x = "", y = "SOR", fill = "") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "none")

    acc_df        <- metrics_df[order(-metrics_df$ACC), ]
    acc_df$Method <- factor(acc_df$Method, levels = acc_df$Method)

    p_acc <- ggplot(acc_df, aes(x = Method, y = ACC, fill = Type)) +
      geom_col(width = 0.6) +
      geom_text(aes(label = round(ACC, 3)), vjust = -0.3, size = 3) +
      scale_fill_manual(values = nrmse_shift_colors) +
      coord_cartesian(ylim = c(0, 1)) +
      labs(title = "a3) ACC — Pearson correlation: original vs imputed",
           subtitle = "Higher = better",
           x = "", y = "Pearson r", fill = "") +
      theme_bw(base_size = 10) +
      theme(axis.text.x = element_text(angle = 30, hjust = 1),
            legend.position = "none")

    pcc_df <- do.call(rbind, lapply(names(imputed_list), function(m) {
      cc      <- suppressWarnings(cor(imputed_list[[m]], use = "complete.obs"))
      melt_cc <- reshape2::melt(cc)
      colnames(melt_cc) <- c("Sample1", "Sample2", "PCC")
      melt_cc$Method <- m
      melt_cc
    }))

    p_pcc <- ggplot(pcc_df, aes(x = Sample1, y = Sample2, fill = PCC)) +
      geom_tile() +
      facet_wrap(~ Method, ncol = min(3, length(imputed_list))) +
      scale_fill_gradient(low = "white", high = "red",
                          limits = c(min(pcc_df$PCC, na.rm = TRUE), 1), name = "Pearson r") +
      labs(title    = "a4) PCC heatmap — inter-sample correlation after imputation",
           subtitle = "High uniform correlation = stable imputation",
           x = "", y = "") +
      theme_bw(base_size = 9) +
      theme(axis.text.x      = element_text(angle = 45, hjust = 1, size = 6),
            axis.text.y      = element_text(size = 6),
            legend.position  = "right",
            strip.text       = element_text(size = 8))

    top_row <- cowplot::plot_grid(p_nrmse, p_sor, p_acc, nrow = 1)
    title   <- cowplot::ggdraw() +
      cowplot::draw_label("NAguideR Classic Evaluation Metrics",
                          fontface = "bold", size = 13)

    cowplot::plot_grid(title, top_row, p_pcc,
                       ncol = 1, rel_heights = c(0.04, 0.35, 0.61))
  }

  # ══════════════════════════════════════════════════════════════════════════
  # PLOT 8 — DE concordance: Jaccard similarity of FDR < fdr_cutoff hits
  # ══════════════════════════════════════════════════════════════════════════
  plot_de_concordance <- function(imputed_list) {

    if (!requireNamespace("limma", quietly = TRUE))
      stop("run_de = TRUE requires the 'limma' package.", call. = FALSE)

    grp_fac <- factor(groups)
    combos  <- combn(levels(grp_fac), 2, simplify = FALSE)

    run_limma <- function(mat) {
      design   <- model.matrix(~ 0 + grp_fac)
      colnames(design) <- levels(grp_fac)
      cont_str <- sapply(combos, function(ab) paste(ab[2], ab[1], sep = "-"))
      cont_mat <- limma::makeContrasts(contrasts = cont_str, levels = design)
      fit <- limma::lmFit(mat, design)
      fit <- limma::contrasts.fit(fit, cont_mat)
      fit <- limma::eBayes(fit)
      ps  <- sapply(cont_str, function(ct)
        limma::topTable(fit, coef = ct, number = Inf,
                        sort.by = "none")$adj.P.Val)
      if (is.matrix(ps)) apply(ps, 1, min) else as.numeric(ps)
    }

    cat("  [DE concordance] running limma on each imputed matrix...\n")
    pval_list <- lapply(names(imputed_list), function(m) {
      cat("   ", m, "...")
      res <- tryCatch(run_limma(imputed_list[[m]]),
                      error = function(e) { cat(" FAILED:", e$message); NULL })
      cat("\n"); res
    })
    names(pval_list) <- names(imputed_list)
    pval_list <- Filter(Negate(is.null), pval_list)

    sig_lists <- lapply(pval_list, function(p) which(p < fdr_cutoff))
    n_sig     <- sapply(sig_lists, length)

    nm  <- names(sig_lists)
    n   <- length(nm)
    jac <- matrix(0, n, n, dimnames = list(nm, nm))
    for (i in seq_len(n)) for (j in seq_len(n)) {
      inter      <- length(intersect(sig_lists[[i]], sig_lists[[j]]))
      union_n    <- length(union(sig_lists[[i]], sig_lists[[j]]))
      jac[i, j]  <- if (union_n == 0) 1 else inter / union_n
    }

    melt_j <- reshape2::melt(jac)
    colnames(melt_j) <- c("Method1", "Method2", "Jaccard")

    ggplot(melt_j, aes(Method1, Method2, fill = Jaccard)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = round(Jaccard, 2)), size = 3.2) +
      scale_fill_gradient(low = "white", high = "#C0392B",
                          limits = c(0, 1), name = "Jaccard") +
      labs(
        title    = paste0("DE Concordance — Jaccard similarity of FDR < ",
                          fdr_cutoff, " hits"),
        subtitle = paste("Sig. proteins per method:",
                         paste(nm, n_sig, sep = "=", collapse = "  |  ")),
        x = "", y = ""
      ) +
      theme_bw(base_size = 11) +
      theme(axis.text.x      = element_text(angle = 35, hjust = 1),
            legend.position  = "right")
  }

  # ══════════════════════════════════════════════════════════════════════════
  # RUN ALL METHODS
  # ══════════════════════════════════════════════════════════════════════════
  cat("Running imputation methods:", paste(methods, collapse = ", "), "\n")
  cat("  Method types:\n")
  for (m in methods) cat("   ", m, "->", method_type(m), "\n")
  cat("\n")

  imputed_list <- list()
  for (m in methods) {
    cat("  ->", m, "... ")
    imputed_list[[m]] <- tryCatch(
      run_imputation(data, m),
      error = function(e) { cat("FAILED:", e$message, "\n"); NULL }
    )
    if (!is.null(imputed_list[[m]])) cat("done\n")
  }
  imputed_list <- Filter(Negate(is.null), imputed_list)
  cat("Successfully imputed:", paste(names(imputed_list), collapse = ", "), "\n\n")

  # ══════════════════════════════════════════════════════════════════════════
  # RUN MASKED IMPUTATION (ONCE)
  # ══════════════════════════════════════════════════════════════════════════
  cat("Running masked imputation for accuracy metrics...\n")
  set.seed(seed)
  observed_idx  <- which(!is.na(data))
  mask_idx      <- sample(observed_idx, round(length(observed_idx) * mask_fraction))
  true_vals     <- data[mask_idx]
  masked_data   <- data
  masked_data[mask_idx] <- NA

  imp_list_masked <- list()
  for (m in names(imputed_list)) {
    cat("  -> Masked", m, "... ")
    imp_list_masked[[m]] <- tryCatch(
      run_imputation(masked_data, m),
      error = function(e) { cat("FAILED\n"); NULL }
    )
    if (!is.null(imp_list_masked[[m]])) cat("done\n")
  }
  imp_list_masked <- Filter(Negate(is.null), imp_list_masked)
  cat("\n")

  # ══════════════════════════════════════════════════════════════════════════
  # GENERATE ALL PLOTS
  # ══════════════════════════════════════════════════════════════════════════
  cat("Generating plots...\n")

  n_methods   <- length(imputed_list)
  nrows_pca   <- ceiling((n_methods + 1) / min(3, n_methods + 1))
  height_pca  <- max(10, 4 * nrows_pca + 1)
  nrows_dens  <- ceiling(n_methods / 3)
  height_dens <- max(8, 4 * nrows_dens)
  nrows_pcc   <- ceiling(n_methods / 3)
  height_pcc  <- max(12, 3 * nrows_pcc + 6)

  cat("  Plot 1: missing value pattern\n")
  p1 <- plot_missing(data)

  cat("  Plot 2: density distributions\n")
  p2 <- plot_densities(data, imputed_list)

  cat("  Plot 3: PCA\n")
  p3 <- plot_pca(data, imputed_list)

  cat("  Plot 4: NRMSE (random mask)\n")
  p4 <- plot_nrmse(imp_list_masked, mask_idx, true_vals)

  cat("  Plot 5: within-group CV\n")
  p5 <- plot_cv(data, imputed_list)

  cat("  Plot 6: downshift\n")
  p6 <- plot_downshift(data, imputed_list)

  cat("  Plot 7: NAguideR classic metrics...\n")
  p7 <- plot_classic_metrics(data, imputed_list, imp_list_masked, mask_idx, true_vals)

  p8 <- NULL
  if (run_de) {
    cat("  Plot 8: DE concordance\n")
    p8 <- tryCatch(
      plot_de_concordance(imputed_list),
      error = function(e) { cat("  DE concordance FAILED:", e$message, "\n"); NULL }
    )
  }

  # ══════════════════════════════════════════════════════════════════════════
  # SAVE TO SINGLE PDF
  # ══════════════════════════════════════════════════════════════════════════
  if (!is.null(output_dir)) {
    pdf_path <- file.path(output_dir, pdf_name)

    pages <- list(
      list(plot = p1, w = 16, h = 10),
      list(plot = p2, w = 16, h = height_dens),
      list(plot = p3, w = 16, h = height_pca),
      list(plot = p4, w = 12, h = 7),
      list(plot = p5, w = 13, h = 7),
      list(plot = p6, w = 10, h = 7),
      list(plot = p7, w = 16, h = height_pcc)
    )
    if (!is.null(p8)) pages <- c(pages, list(list(plot = p8, w = 11, h = 8)))

    if (requireNamespace("pdftools", quietly = TRUE)) {
      tmp_files <- sapply(seq_along(pages), function(i) {
        if (is.null(pages[[i]]$plot)) return(NA)
        f <- tempfile(fileext = ".pdf")
        grDevices::pdf(f, width = pages[[i]]$w, height = pages[[i]]$h)
        print(pages[[i]]$plot)
        grDevices::dev.off()
        f
      })
      tmp_files <- tmp_files[!is.na(tmp_files)]
      pdftools::pdf_combine(tmp_files, output = pdf_path)
      file.remove(tmp_files)
      cat("Saved combined PDF (variable page sizes) to:", pdf_path, "\n")
    } else {
      grDevices::pdf(pdf_path, onefile = TRUE, width = 16, height = 10)
      for (pg in pages) if (!is.null(pg$plot)) print(pg$plot)
      grDevices::dev.off()
      cat("Saved PDF to:", pdf_path, "\n")
      cat("  Tip: install pdftools for variable page heights.\n")
    }
  }

  invisible(list(
    imputed             = imputed_list,
    plot_missing        = p1,
    plot_density        = p2,
    plot_pca            = p3,
    plot_nrmse          = p4,
    plot_cv             = p5,
    plot_downshift      = p6,
    plot_classic        = p7,
    plot_de_concordance = p8
  ))
}
