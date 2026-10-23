#' PCA Plot for Proteomics Expression Matrices
#'
#' Generates a publication-ready PCA plot from a proteomics (or general
#' expression) matrix. The function supports feature selection based on
#' variance, NIPALS/PPCA PCA when missing values are present, customizable color
#' and shape aesthetics, optional ellipses per group, custom sample ordering,
#' annotation settings, and automatic high-resolution PNG export.
#'
#' @param expression_matrix Numeric matrix of expression values
#'   (rows = proteins/genes, columns = samples).
#' @param sample_info Data frame containing sample-level metadata. Must include:
#'   \itemize{
#'     \item `Sample` — sample names matching column names of the matrix
#'     \item `Group` — biological group (default color variable)
#'   }
#' @param pcsToUse Integer vector of length 2 specifying which PCs to plot.
#'   Default: \code{1:2}.
#' @param pca_method Method for pca . Can be either classic or robust.
#' @param top_n Number of most variable features to use. `"all"` uses all
#'   rows. Default: `"all"`.
#' @param scale Logical.
#' @param show_names Logical. If TRUE, sample names are labeled. Default: TRUE.
#' @param clean_names Logical. If TRUE, removes `_1`, `_2` suffixes and labels only 1 replicate per pair. Default: FALSE.
#' @param add_ellipse Logical. If TRUE, draws 95% confidence ellipses per group.
#' @param fill_ellipse Logical. If TRUE, ellipses are filled with semi-transparent
#'   group colors.
#' @param grid_lines Logical. If FALSE, grid lines are removed. Default: FALSE.
#' @param width Width of the saved PCA PNG (in inches). Default: 18.
#' @param height Height of the saved PCA PNG (in inches). Default: 15.
#' @param output_path File path for saving the PCA image. Default:
#'   `"PCA_plot.png"`.
#' @param custom_order Optional character vector specifying custom ordering of
#'   groups (applies only when \code{color_by == "group"}).
#' @param group_colors Optional named vector assigning colors to groups.
#'   Example: \code{c("WT"="steelblue","KO"="firebrick")}.
#' @param color_by Column in `sample_info` used for point colors.
#'   Default: `"Group"`.
#' @param shape_by Optional column in `sample_info` used for point shapes.
#'   If NULL, all points share the same shape.
#' @param mv_handle Method for missing values. Options: "nipals", "ppca", "pcana".
#' @param nipals_scale Scaling method for pcaMethods. Default: "pareto".
#' @param pcana_method Method for rrcovNA::PcaNA. Default: "class".
#'
#' @return A list containing:
#' \itemize{
#'    \item **pca_plot** — ggplot object of the PCA plot
#'    \item **pca_data** — data frame with PCA coordinates and metadata
#'    \item **pca_res** — raw PCA result from `prcomp` or `pcaMethods::pca`
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point theme_bw
#' @importFrom pcaMethods pca scores
#' @importFrom ggrepel geom_label_repel
#' @importFrom grid unit
#' @import matrixStats
#' @export
my_plot_PCA_matrix <- function(expression_matrix,
                               sample_info,
                               pcsToUse = 1:2,
                               pca_method = "classic",
                               top_n ='all',
                               scale = TRUE,
                               show_names = TRUE,
                               clean_names = FALSE,
                               add_ellipse = FALSE,
                               fill_ellipse = FALSE,
                               grid_lines = FALSE,
                               width = 18,
                               height = 15,
                               output_path = "PCA_plot.png",
                               custom_order = NULL,
                               group_colors = NULL,
                               color_by = "Group",
                               shape_by = NULL,
                               mv_handle = "ppca",
                               nipals_scale = "pareto",
                               pcana_method = "class") {

  # 1. Remove proteins with all NA before variance calculation
  expression_matrix <- expression_matrix[rowSums(!is.na(expression_matrix)) > 0, , drop = FALSE]
  row_vars   <- matrixStats::rowVars(as.matrix(expression_matrix), na.rm = TRUE)

  if (top_n != 'all') {
    top_n_val    <- min(as.numeric(top_n), sum(!is.na(row_vars)))
    top_features <- names(sort(row_vars[!is.na(row_vars)], decreasing = TRUE))[1:top_n_val]
    expression_matrix <- expression_matrix[top_features, , drop = FALSE]
  }

  # 2. Βεβαιώσου ότι η σειρά ταιριάζει με τις στήλες του πίνακα
  sample_info <- sample_info[match(colnames(expression_matrix), sample_info$Sample), ]

  # 3. PCA (samples στις γραμμές) — auto-handles NA via selected method
  expression_matrix[!is.finite(as.matrix(expression_matrix))] <- NA  # Clean Inf/NaN

  mv_handle <- match.arg(arg = mv_handle, choices = c("nipals", "ppca", "pcana"), several.ok = FALSE)
  pca_method <- match.arg(arg = pca_method, choices = c("classic", "robust"), several.ok = FALSE)

  if (anyNA(expression_matrix)) {
    if (mv_handle %in% c("nipals", "ppca")) {
      #-------------------------- pcaMethods (NIPALS/PPCA) -------------------------
      message(paste0("⚠️ Missing values detected — using ", toupper(mv_handle), " PCA (pcaMethods)."))
      pca_res_obj <- pcaMethods::pca(t(expression_matrix),
                                     method = mv_handle,
                                     nPcs = max(pcsToUse),
                                     scale = nipals_scale)

      scores_mat <- pcaMethods::scores(pca_res_obj)
      pca_res <- list(x = scores_mat)
      r2_vals <- summary(pca_res_obj)["R2", ]
      percentVar <- round(100 * r2_vals[pcsToUse])
      #-----------------------------------------------------------------------------
    } else if (mv_handle == "pcana") {
      #------------------------------ PCA NA ---------------------------------------
      message("⚠️ Missing values detected — using Robust PCA (PcaNA from rrcovNA).")
      pca_res_robust <- rrcovNA::PcaNA(t(expression_matrix), k = 10, method = pcana_method, scale = scale)
      scores_mat <- rrcovNA::getScores(pca_res_robust)
      pca_res <- list(x = scores_mat)
      eigenvalues <- rrcovNA::getEigenvalues(pca_res_robust)
      total_var <- sum(eigenvalues)
      percentVar <- round(100 * (eigenvalues[pcsToUse] / total_var))
      #-----------------------------------------------------------------------------
    }
  } else {
    if (pca_method == "classic") {
      #------------------------------ STANDARD PCA ---------------------------------
      message("✅ No missing values — using standard PCA.")
      pca_res <- prcomp(t(expression_matrix), scale. = scale, rank. = 10)
      percentVar <- round(100 * (pca_res$sdev^2 / sum(pca_res$sdev^2)))[pcsToUse]
      #-----------------------------------------------------------------------------
    } else if (pca_method == "robust") {
      #------------------------------ ROBUST PCA -----------------------------------
      message("✅ No missing values — using Robust PCA (PcaHubert).")
      if (!requireNamespace("rrcov", quietly = TRUE)) stop("rrcov required for PcaHubert")
      pca_res_robust <- rrcov::PcaHubert(t(expression_matrix), k = 10, scale = scale)
      scores_mat <- rrcov::getScores(pca_res_robust)
      pca_res <- list(x = scores_mat)
      eigenvalues <- rrcov::getEigenvalues(pca_res_robust)
      total_var <- sum(eigenvalues)
      percentVar <- round(100 * (eigenvalues[pcsToUse] / total_var))
      #-----------------------------------------------------------------------------
    }
  }

  # 5. Data frame για ggplot
  pca_data <- data.frame(
    xaxis = pca_res$x[, pcsToUse[1]],
    yaxis = pca_res$x[, pcsToUse[2]],
    Group = sample_info$Group,
    name = sample_info$Sample
  )

  # Add color_by and shape_by columns
  if (!color_by %in% colnames(sample_info)) stop(paste("Missing column:", color_by))
  grp_vals <- as.character(sample_info[[color_by]])
  grp_lv <- if (!is.null(group_colors) && all(unique(grp_vals) %in% names(group_colors)))
    names(group_colors)[names(group_colors) %in% unique(grp_vals)] else unique(grp_vals)   # same order as bar / box plots
  pca_data[[color_by]] <- factor(grp_vals, levels = grp_lv)

  if (!is.null(shape_by)) {
    if (!shape_by %in% colnames(sample_info)) stop(paste("Missing column:", shape_by))
    pca_data[[shape_by]] <- as.factor(sample_info[[shape_by]])
  }

  # Apply custom order if set (only if color_by == "Group")
  if (!is.null(custom_order) && color_by == "Group") {
    pca_data[[color_by]] <- factor(pca_data[[color_by]], levels = custom_order)
  }

  # Clean sample labels (e.g., convert "33_1" -> "33" and leave duplicates blank)
  if (clean_names) {
    clean_label <- sub("_[0-9]+$", "", pca_data$name)
    pca_data$display_label <- ifelse(duplicated(clean_label), "", clean_label)
  } else {
    pca_data$display_label <- pca_data$name
  }

  # Aesthetics selection
  use_shape <- !is.null(shape_by) && color_by != shape_by

  if (use_shape) {
    aes_args <- aes(x = xaxis, y = yaxis, color = .data[[color_by]], fill = .data[[color_by]], shape = .data[[shape_by]])
  } else {
    aes_args <- aes(x = xaxis, y = yaxis, color = .data[[color_by]], fill = .data[[color_by]])
  }

  # Initialize plot
  pca_plot <- ggplot(pca_data, aes_args)

  # 6. Plot layers
  if (show_names) {
    pca_plot <- pca_plot +
      ggrepel::geom_text_repel(
        aes(label = display_label),
        fontface      = "bold",
        size          = 3.8,
        box.padding   = unit(0.35, "lines"),
        point.padding = unit(0.8, "lines"),
        segment.color = "grey60",
        show.legend   = FALSE
      )
  }

  # Add points
  if (use_shape) {
    pca_plot <- pca_plot + geom_point(size = 4, stroke = 1.0)

    # Custom shapes for Sex or standard filled shapes
    if (shape_by == "Sex") {
      # Custom gender shape symbols
      pca_plot <- pca_plot + scale_shape_manual(values = c("Female" = 17, "Male" = 15))
    } else {
      pca_plot <- pca_plot + scale_shape_manual(values = c(21, 22, 24, 25))
    }
  } else {
    pca_plot <- pca_plot + geom_point(size = 4, shape = 21, stroke = 0.8)
  }

  # Apply custom colors if provided
  if (!is.null(group_colors)) {
    # Ensure group_colors works dynamically with whatever column is passed to color_by
    pca_plot <- pca_plot +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors)
  }

  # Add ellipse
  if (add_ellipse) {
    pca_plot <- pca_plot +
      stat_ellipse(aes(group = .data[[color_by]],
                       fill = if (fill_ellipse) .data[[color_by]] else NA),
                   geom = "polygon",
                   type = "t",
                   level = 0.95,
                   linetype = "dashed",
                   linewidth = 1,
                   alpha = if (fill_ellipse) 0.2 else 0) +
      guides(fill = "none")
  }

  # Αφαίρεση grid lines αν ζητηθεί
  if (!grid_lines) {
    pca_plot <- pca_plot +
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank())
  }

  # Labels αξόνων
  pca_plot <- pca_plot +
    xlab(paste0("PC", pcsToUse[1], ": ", percentVar[1], "% variance")) +
    ylab(paste0("PC", pcsToUse[2], ": ", percentVar[2], "% variance")) +
    theme_bw()

  # Save
  ggsave(output_path, plot = pca_plot, width = width, height = height)

  return(list(pca_plot = pca_plot,
              pca_data = pca_data,
              pca_res = pca_res))
}
