#' Create a Volcano Plot for Proteomics Differential Expression Results
#'
#' This function generates a publication-ready volcano plot from proteomics
#' differential expression results (typically from limma). It classifies
#' proteins into up-regulated, down-regulated, low-regulation, or non-significant
#' groups, applies fold-change and significance thresholds, labels selected
#' proteins, and saves the resulting plot as a PNG file.
#'
#' @param DE_res A data frame containing differential expression results.
#'   Must include columns: `logFC`, a p-value column (e.g. `P.Value` or
#'   `adj.P.Val`), and a `Gene` column for labeling.
#' @param genes_of_interest Optional vector of genes to highlight (currently unused).
#' @param filter_choice Character string specifying the significance column
#'   (e.g. `"adj.P.Val"` or `"P.Value"`). Default: `"adj.P.Val"`.
#' @param y_axis_choice Column used for the y-axis transformation.
#'   Default: `"P.Value"`.
#' @param l2FC_thresh Numeric log2 fold-change threshold. Default: `1`.
#' @param p_thresh Numeric p-value/FDR significance threshold. Default: `0.05`.
#' @param plot_title Main plot title. Default: `"Volcano Plot"`.
#' @param show_names `"none"` (default), `"all"`, or a numeric value indicating
#'   how many top genes to label by significance.
#' @param which_genes Which set of genes to label: `"up"`, `"down"`, or `"all"`
#'   (default `"all"`).
#' @param grid_lines Logical. If `FALSE`, removes background grid lines.
#' @param save_prefix Optional prefix for the output PNG filename.
#' @param width Numeric width of the saved plot in inches. Default: `12`.
#' @param height Numeric height of the saved plot in inches. Default: `10`.
#' @param label_size Numeric. Default:`2`
#'
#' @return A `ggplot2` volcano plot object.
#'         A PNG file is saved as:
#'         `"Volcano_plot_<save_prefix>.png"`.
#'
#' @details
#' The function:
#' \enumerate{
#'   \item Categorizes proteins based on fold change and significance.
#'   \item Computes \eqn{-\log_{10}(p)} for the y-axis.
#'   \item Draws significance and fold-change threshold lines.
#'   \item Labels proteins according to `show_names` and `which_genes`.
#'   \item Saves a high-resolution PNG to disk.
#' }
#'
#' Labeling rules:
#' \itemize{
#'   \item `"none"` — no labels
#'   \item `"all"` — label *all* significant genes
#'   \item numeric — label *top N* by significance
#' }
#'
#' @examples
#' # Basic usage
#' # v <- my_volcano_plot_proteomics(DE_res)
#'
#' # Label top 10 upregulated genes
#' # v <- my_volcano_plot_proteomics(DE_res, show_names = 10, which_genes = "up")
#'
#' @import ggplot2
#'
#' @export
my_volcano_plot_proteomics <- function(
    DE_res,
    genes_of_interest = c(),
    filter_choice = "adj.P.Val",
    y_axis_choice = "P.Value",
    l2FC_thresh = 1,
    p_thresh = 0.05,
    plot_title = "Volcano Plot",
    show_names = "none",
    which_genes = "all",
    grid_lines = TRUE,
    save_prefix = NULL,
    width = 12,
    height = 10,
    label_size = 2
) {

  DE_res$y_value <- -log10(DE_res[[y_axis_choice]])
  # --- Format Thresholds (Max 2 Decimals for l2FC_thresh) ---
  l2fc_str <- format(round(l2FC_thresh, 2), nsmall = 0)

  # --- Calculate Horizontal Line Directly ---
  if (filter_choice == y_axis_choice) {
    y_intercept_val <- -log10(p_thresh)
  } else {
    sig_rows     <- which(DE_res[[filter_choice]] < p_thresh)
    non_sig_rows <- which(DE_res[[filter_choice]] >= p_thresh)

    if (length(sig_rows) > 0 && length(non_sig_rows) > 0) {
      # Transition boundary: Find the non-sig protein closest to p_thresh
      first_nonsig_idx <- non_sig_rows[which.min(DE_res[[filter_choice]][non_sig_rows])]
      y_intercept_val  <- -log10(DE_res[[y_axis_choice]][first_nonsig_idx])
    } else if (length(sig_rows) == 0) {
      # ZERO significant proteins: Place intercept above the top-ranked protein
      # or default to -log10(p_thresh) if p_value max is lower
      max_y_data      <- max(DE_res$y_value, na.rm = TRUE)
      y_intercept_val <- max(max_y_data + 0.1, -log10(p_thresh))
    } else {
      # ALL proteins are significant
      y_intercept_val <- -log10(p_thresh)
    }
  }

  # --- Identify Indices for Categories ---
  sigUp   <- which(DE_res$logFC > l2FC_thresh & DE_res[[filter_choice]] < p_thresh)
  sigDown <- which(DE_res$logFC < (-l2FC_thresh) & DE_res[[filter_choice]] < p_thresh)

  # --- Build Display Labels ---
  up_label    <- "up"
  down_label  <- "down"

  nonSig_count <- nrow(DE_res) - length(sigUp) - length(sigDown)
  nonsig_label <- "ns"

  # --- Assign Significance Status ---
  DE_res$Significance <- nonsig_label
  if (length(sigUp) > 0)   DE_res[sigUp, "Significance"]   <- up_label
  if (length(sigDown) > 0) DE_res[sigDown, "Significance"] <- down_label

  # --- Compute Y-axis ---
  DE_res$y_value <- -log10(DE_res[[y_axis_choice]])

  # --- Sort for labeling ---
  DE_res <- DE_res[order(DE_res[[filter_choice]]), ]

  # --- Select Labels ---
  label_data <- data.frame()

  if (show_names != "none") {
    if (length(genes_of_interest) > 0) {
      label_data <- DE_res[DE_res$Gene %in% genes_of_interest, ]
    } else {
      if (which_genes == "up") {
        label_data <- DE_res[sigUp, ]
      } else if (which_genes == "down") {
        label_data <- DE_res[sigDown, ]
      } else if (which_genes == "all") {
        label_data <- DE_res[c(sigUp, sigDown), ]
      } else {
        warning("Invalid value for which_genes; using 'all'")
        label_data <- DE_res[c(sigUp, sigDown), ]
      }
    }
    if (is.numeric(show_names)) {
      label_data <- head(label_data[order(-label_data$y_value), ], show_names)
    } else if (show_names != "all") {
      warning("Invalid value for show_names; using 'none'")
      label_data <- data.frame()
    }
  }

  # --- Color Mapping ---
  color_val_map <- setNames(
    c("#B31B21", "#1465AC", "grey70"),
    c(up_label, down_label, nonsig_label)
  )

  # --- Plot Construction ---
  v_plot <- ggplot(DE_res, aes(x = logFC, y = y_value, color = Significance)) +
    geom_point(size = 1.2) +
    geom_hline(yintercept = y_intercept_val, linetype = "dashed", color = "black") +
    geom_vline(xintercept = c(-l2FC_thresh, l2FC_thresh), linetype = "dashed", color = "black") +
    xlab(expression("Log"[2] * " Fold Change")) +
    ylab(paste0("-log10(", y_axis_choice, ")")) +
    scale_color_manual(
      breaks = c(up_label, down_label, nonsig_label),
      values = color_val_map
    ) +
    guides(color = guide_legend(override.aes = list(size = 4))) +
    ggtitle(plot_title) +
    theme_light() +
    (if (!grid_lines) theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) else theme()) +
    theme(
      text = element_text(size = 8),
      plot.title = element_text(size = 10, hjust = 0.5),
      legend.text = element_text(size = 7),
      plot.title.position = "panel"
    )

  if (nrow(label_data) > 0) {
    v_plot <- v_plot +
      ggrepel::geom_label_repel(
        data = label_data, aes(label = Gene, color = Significance),
        box.padding = unit(0.2, "lines"), point.padding = unit(0.1, "lines"),
        force = 20, max.overlaps = Inf,
        segment.color = "black", show.legend = FALSE,
        label.size = 0.15,
        size = label_size
      )
  }

  out_name <- if (!is.null(save_prefix)) paste0("Volcano_plot_", save_prefix, ".png") else "Volcano_plot.png"
  ggsave(out_name, plot = v_plot, width = width, height = height)

  return(v_plot)
}
