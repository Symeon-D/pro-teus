#' Custom GO Enrichment Dot Plot
#'
#' This function creates a dot plot for visualizing GO enrichment results.
#'
#' @param GO_clusters A data frame containing the GO cluster results.
#' @param ont The ontology (e.g., "BP", "MF", "CC") used for the analysis.
#' @param top_n The number of top results to display in the plot.
#' @param term_size_legend A character string indicating the legend type for term size. Can be either "numeric" or "percentage".
#' @param sub_cluster (Optional) A specific cluster name to filter the data by. Defaults to NULL.
#' @param lolipop Whether to plot the lines of the lilipop.
#' @param multi_cluster Whether to plot a dot plot with all the clusters together.
#'
#' @return A ggplot object representing the GO enrichment dot plot.
#'
#' @export
my_custom_gpro_dotplot <- function(GO_clusters, ont, top_n, term_size_legend, sub_cluster, lolipop=FALSE, multi_cluster,
                                 shape_map = NULL, x_angle = -45, lollipop_x = "padj",
                                 color_metric = "auto") {

  # Filter by sub-cluster if provided
  if (!is.null(sub_cluster)) {
    GO_clusters <- GO_clusters[GO_clusters$group_name == sub_cluster, ]
  }

  # Prepare term names and title
  # lollipop x axis: "padj" = -log10(padj) (colour = fold enrichment); "fold" = log2 fold enrichment (colour = -log10(padj))
  x_var <- if (!multi_cluster && lollipop_x == "fold") "log2(Fold_Enrichment)" else "-log10(padj)"
  x_title <- if (x_var == "log2(Fold_Enrichment)") "log2(Fold enrichment)" else "-log10(padj)"
  # colour bar: "auto" keeps the classic choice, "padj" / "fold" force the metric
  color_var <- if (color_metric == "padj") "-log10(padj)" else if (color_metric == "fold") "log2(Fold_Enrichment)" else
    if (multi_cluster) "-log10(padj)" else if (x_var == "-log10(padj)") "log2(Fold_Enrichment)" else "-log10(padj)"
  if (!color_var %in% colnames(GO_clusters)) color_var <- "-log10(padj)"
  color_title <- if (color_var == "-log10(padj)") "-log10(padj)" else "Fold Enrichment"
  GO_clusters$term_name <- factor(GO_clusters$term_name, levels = unique(GO_clusters$term_name[order(GO_clusters[[x_var]])]))
  title_text <- ont
  if (!is.null(sub_cluster)) {
    title_text <- paste(title_text, "| Cluster:", sub_cluster)
  }

  # Optionally compute percentage
  if (term_size_legend == "percentage") {
    GO_clusters$term_size_percentage <- (GO_clusters$intersection_size / GO_clusters$query_size) * 100
  }

  # Base aesthetics and plot
  if (multi_cluster) {
    plt <- ggplot(GO_clusters, aes(x = group_name, y = term_name))
  } else {
    plt <- ggplot(GO_clusters, aes(x = .data[[x_var]], y = term_name))
    if (lolipop) {
      plt <- plt + geom_segment(aes(x = 0, xend = .data[[x_var]], y = term_name, yend = term_name, color = round(.data[[color_var]], 1)), linewidth = 1, alpha = 0.7)
    }
  }

  # Point aesthetic mapping
  # "Gene count" = number of INPUT genes found in the pathway (intersection), not the pathway size
  size_var <- if (term_size_legend == "percentage") "term_size_percentage" else "intersection_size"

  if (multi_cluster && !is.null(shape_map)) {
    # shape_map: named vector group -> 21 (circle), 24 (up arrow), 25 (down arrow)
    plt <- plt + geom_point(aes(size = .data[[size_var]], color = round(.data[[color_var]], 1),
                                fill = round(.data[[color_var]], 1), shape = group_name), show.legend = TRUE) +
      scale_shape_manual(values = shape_map, guide = "none", drop = FALSE)
  } else {
    plt <- plt + geom_point(aes(size = .data[[size_var]], color = round(.data[[color_var]], 1)), show.legend = TRUE)
  }

  # Handle size scaling
  unique_sizes <- sort(unique(GO_clusters[[size_var]]))
  if (length(unique_sizes) > 5) {
    size_breaks <- unique_sizes[round(seq(1, length(unique_sizes), length.out = 5))]
  } else {
    size_breaks <- unique_sizes
  }

  size_label <- if (term_size_legend == "percentage") "Gene Ratio (%)" else "Gene Count"

  # legend markers follow the cluster markers: one shape everywhere -> that shape in the single legend;
  # a mix -> one size legend per shape, each with ITS OWN clusters' representative values
  size_legend_override <- NULL; mixed <- FALSE; size_limits <- range(GO_clusters[[size_var]], na.rm = TRUE)
  if (multi_cluster && !is.null(shape_map)) {
    shp_by_group <- shape_map[as.character(levels(GO_clusters$group_name))]
    shp_u <- unique(as.numeric(shp_by_group[!is.na(shp_by_group)]))
    if (length(shp_u) == 1) size_legend_override <- list(shape = shp_u, fill = "grey30", colour = "grey30")
    if (length(shp_u) > 1) mixed <- TRUE
  }
  plt <- plt + scale_size(
    name = size_label,
    limits = size_limits,
    breaks = size_breaks,
    labels = round(size_breaks, 1),
    guide = if (mixed) "none" else guide_legend(title.position = "top")
  )

  # Handle color scaling
  color_range <- range(GO_clusters[[color_var]])
  if (length(unique(GO_clusters[[color_var]])) > 1) {
    plt <- plt + scale_color_gradientn(
      colours = c("blue", "red"),
      limits = c(floor(color_range[1]), ceiling(color_range[2])),
      name = color_title
    )
  } else {
    plt <- plt + labs(color = color_title)
  }

  if (multi_cluster && !is.null(shape_map) && length(unique(GO_clusters[[color_var]])) > 1) {
    plt <- plt + scale_fill_gradientn(colours = c("blue", "red"),
                                      limits = c(floor(color_range[1]), ceiling(color_range[2])), guide = "none")
  }

  # Shared plot theming
  plt <- plt +
    ggtitle(title_text) +
    theme_bw(base_size = 15) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    guides(size = if (mixed) "none" else if (is.null(size_legend_override)) guide_legend(order = 1) else
             guide_legend(order = 1, override.aes = size_legend_override),
           color = guide_colorbar(order = 2)) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))

  if (multi_cluster) {
    plt <- plt + theme(axis.text.x = element_text(angle = x_angle,
                                                  vjust = if (abs(x_angle) == 90) 0.5 else if (x_angle == 0) 1 else 0,
                                                  hjust = if (x_angle > 0) 1 else if (x_angle < 0) 0 else 0.5)) +
      scale_x_discrete(drop = FALSE) + xlab("") + ylab("")
  } else {
    plt <- plt +
      scale_x_continuous(
        breaks = seq(min(0, floor(min(GO_clusters[[x_var]], na.rm = TRUE))), ceiling(max(GO_clusters[[x_var]], na.rm = TRUE)), by = 1),
        limits = c(min(0, floor(min(GO_clusters[[x_var]], na.rm = TRUE))), ceiling(max(GO_clusters[[x_var]], na.rm = TRUE)))
      ) +
      xlab(x_title) + ylab("")
  }

  plt <- plt + theme(legend.spacing.y = unit(0.05, "cm"), legend.justification = "top", legend.box.just = "left")

  if (mixed) {
    pick5 <- function(v) { u <- sort(unique(v)); if (length(u) > 5) u[round(seq(1, length(u), length.out = 5))] else u }
    # colour bar (the only legend left in the main plot)
    leg_col <- cowplot::get_legend(plt + theme(legend.position = "right"))
    plt0 <- plt + theme(legend.position = "none")
    legs <- lapply(shp_u[1:2], function(sh) {
      grp <- names(shp_by_group)[as.numeric(shp_by_group) == sh]
      v <- GO_clusters[[size_var]][as.character(GO_clusters$group_name) %in% grp]
      if (!length(v)) v <- GO_clusters[[size_var]]
      b <- pick5(v)
      g <- ggplot(data.frame(x = 1, y = 1, s = b), aes(x, y, size = s)) +
        geom_point(shape = sh, fill = "grey30", colour = "grey30") +
        scale_size(name = NULL, limits = size_limits, breaks = b, labels = round(b, 1)) +
        guides(size = guide_legend(override.aes = list(shape = sh, fill = "grey30", colour = "grey30"))) +
        theme_void(base_size = 15) + theme(legend.position = "right")
      cowplot::get_legend(g)
    })
    ttl <- cowplot::ggdraw() + cowplot::draw_label(size_label, fontface = "plain", size = 15, hjust = 0.5)
    topleft <- function(g) cowplot::ggdraw() + cowplot::draw_grob(g, halign = 0, valign = 1)   # no centring -> no big gaps
    lw <- function(g) tryCatch(grid::convertWidth(sum(g$widths), "in", valueOnly = TRUE), error = function(e) 1)
    size_row <- cowplot::plot_grid(topleft(legs[[1]]), topleft(legs[[2]]), nrow = 1,
                                   rel_widths = c(lw(legs[[1]]), lw(legs[[2]])))
    # stack with their NATURAL heights (title, size legends, colour bar) and a spacer below, so the
    # size legend sits right on top of the colour bar instead of floating at the top of the figure
    hin <- function(g) tryCatch(grid::convertHeight(sum(g$heights), "in", valueOnly = TRUE), error = function(e) 1.6)
    H_est <- max(3.8, 0.26 * length(unique(GO_clusters$term_name)) + 1.6)            # figure height in inches
    h_t <- 0.3; h_s <- max(hin(legs[[1]]), hin(legs[[2]])) + 0.05; h_c <- hin(leg_col) + 0.05
    right <- cowplot::plot_grid(ttl, size_row, topleft(leg_col), NULL, ncol = 1,
                                rel_heights = c(h_t, h_s, h_c, max(0.2, H_est - h_t - h_s - h_c)))
    plt <- cowplot::plot_grid(plt0, right, nrow = 1, rel_widths = c(1, 0.24))
  }

  return(plt)
}
