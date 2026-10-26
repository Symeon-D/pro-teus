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

  # legend markers follow the cluster markers: one shape everywhere -> one column of that shape; two different
  # shapes -> ONE legend with two columns (5 representative sizes per shape, taken from that shape's own clusters).
  # A single guide keeps the same position / title as the one-shape case and keeps the columns close together.
  size_limits <- range(GO_clusters[[size_var]], na.rm = TRUE)
  size_guide <- guide_legend(order = 1, title.position = "top")
  size_breaks_final <- size_breaks
  if (multi_cluster && !is.null(shape_map)) {
    shp_by_group <- shape_map[as.character(levels(GO_clusters$group_name))]
    shp_u <- unique(as.numeric(shp_by_group[!is.na(shp_by_group)]))
    if (length(shp_u) == 1) {
      size_guide <- guide_legend(order = 1, title.position = "top",
                                 override.aes = list(shape = shp_u, fill = "grey30", colour = "grey30"))
    } else if (length(shp_u) > 1) {
      pick5 <- function(v) {
        u <- sort(unique(v))
        if (length(u) >= 5) u[round(seq(1, length(u), length.out = 5))]
        else if (length(u) > 1) seq(min(u), max(u), length.out = 5)
        else seq(size_limits[1], size_limits[2], length.out = 5)
      }
      reps <- lapply(shp_u[1:2], function(sh) {
        grp <- names(shp_by_group)[as.numeric(shp_by_group) == sh]
        v <- GO_clusters[[size_var]][as.character(GO_clusters$group_name) %in% grp]
        pick5(if (length(v)) v else GO_clusters[[size_var]])
      })
      b1 <- reps[[1]]; b2 <- reps[[2]]
      # breaks must be distinct: nudge a value of group 2 that equals one of group 1 (inside the scale limits)
      b2 <- ifelse(b2 %in% b1, ifelse(b2 >= size_limits[2], b2 * (1 - 1e-6), b2 * (1 + 1e-6)), b2)
      size_breaks_final <- c(b1, b2)                   # NOT sorted: column 1 = shape 1, column 2 = shape 2
      size_guide <- guide_legend(order = 1, title.position = "top", ncol = 2, byrow = FALSE,
                                 override.aes = list(shape = c(rep(shp_u[1], length(b1)), rep(shp_u[2], length(b2))),
                                                     fill = "grey30", colour = "grey30"))
    }
  }
  plt <- plt + scale_size(
    name = size_label,
    limits = size_limits,
    breaks = size_breaks_final,
    labels = round(size_breaks_final, 1),
    guide = size_guide
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
    guides(size = size_guide, color = guide_colorbar(order = 2)) +
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
  if (utils::packageVersion("ggplot2") >= "3.5.0")      # keep the two size columns close together
    plt <- plt + theme(legend.key.spacing.x = unit(0.08, "cm"))

  return(plt)
}
