# ==========================================================
# Heatmap functions (ComplexHeatmap) - extended version of
# my_prepare_top_annotation_proteomics() / my_proteomics_heat_map()
# ==========================================================
library(ComplexHeatmap)
library(circlize)
library(grid)

# Largest font (pt) at which `labels` fit inside a box of width_mm x height_mm
fit_font <- function(labels, width_mm, height_mm = Inf, max_font = 20, min_font = 4) {
  labels <- as.character(labels)
  n <- max(nchar(labels), 1)
  by_w <- width_mm * 0.92 / (n * 0.3528 * 0.58)       # avg glyph width ~0.58 em
  by_h <- height_mm * 0.75 / 0.3528
  max(min_font, min(max_font, by_w, by_h))
}

my_prepare_top_annotation_proteomics <- function(sample_info,
                                                 group_colors,
                                                 custom_order = NULL,
                                                 col_fontsize = NA,
                                                 group_fontsize = NA,
                                                 merge = FALSE,
                                                 extra_vars = NULL,
                                                 extra_colors = NULL,
                                                 show_sample_names = TRUE,
                                                 show_group_names = TRUE,
                                                 col_width_mm = 8,
                                                 legend_fontsize = 14) {
  force(show_sample_names); force(show_group_names); force(col_fontsize); force(group_fontsize)

  sample_info$Group <- factor(sample_info$Group, levels = unique(sample_info$Group))
  block_h <- 8   # mm, height of the sample / extra-variable boxes

  if (merge) {
    col_split <- factor(levels(sample_info$Group), levels = levels(sample_info$Group))
    fill_colors <- group_colors[as.character(levels(sample_info$Group))]
    block_labels <- as.character(levels(sample_info$Group))
  } else {
    col_split <- factor(sample_info$Sample, levels = sample_info$Sample)
    fill_colors <- group_colors[as.character(sample_info$Group)]
    block_labels <- sample_info$Sample
  }
  fill_set_groups <- gpar(fill = unname(fill_colors))

  # Gaps between groups
  gap_vect <- c()
  if (merge) {
    gap_vect <- rep(1, length(levels(sample_info$Group)) - 1)
  } else {
    group_vec <- as.character(sample_info$Group)
    for (i in seq_len(length(group_vec) - 1))
      gap_vect <- c(gap_vect, if (group_vec[i] == group_vec[i + 1]) 0.5 else 1)
  }

  # One font for all boxes, small enough for the longest label (auto-fit)
  # font size: a number = used exactly; NA (blank in the UI) = auto-fit to the boxes
  box_font <- if (is.na(col_fontsize)) fit_font(block_labels, col_width_mm, block_h, max_font = 15) else col_fontsize

  annot_args <- list(empty = anno_empty(border = FALSE, height = unit(10, "mm")))
  if (merge) {
    annot_args$groups <- anno_block(gp = fill_set_groups, labels = NULL,
                                    height = unit(block_h, "mm"))
  } else {
    annot_args$samples <- anno_block(
      gp = fill_set_groups, height = unit(block_h, "mm"),
      labels = if (show_sample_names) block_labels else NULL,
      labels_gp = gpar(col = "black", fontsize = box_font))
  }

  legends <- list(Legend(title = "Group", at = names(group_colors),
                         legend_gp = gpar(fill = unname(group_colors)),
                         title_gp = gpar(fontsize = legend_fontsize, fontface = "bold"),
                         labels_gp = gpar(fontsize = legend_fontsize)))

  if (!is.null(extra_vars) && ncol(extra_vars) > 0) {
    for (v in colnames(extra_vars)) {
      vals <- as.character(extra_vars[[v]])
      annot_args[[v]] <- anno_block(
        gp = gpar(fill = unname(extra_colors[[v]][vals])),
        height = unit(block_h, "mm"),
        labels = if (show_sample_names) vals else NULL,
        labels_gp = gpar(col = "black",
                         fontsize = if (is.na(col_fontsize)) fit_font(vals, col_width_mm, block_h, max_font = 15) else col_fontsize))
      legends[[length(legends) + 1]] <- Legend(
        title = v, at = names(extra_colors[[v]]),
        legend_gp = gpar(fill = unname(extra_colors[[v]])),
        title_gp = gpar(fontsize = legend_fontsize, fontface = "bold"),
        labels_gp = gpar(fontsize = legend_fontsize))
    }
  }
  top_annot_labels <- do.call(HeatmapAnnotation, annot_args)

  group_vec <- as.character(sample_info$Group)
  group_labels <- names(table(sample_info$Group))
  group_indices <- split(seq_along(group_vec), group_vec)

  # Group names above the blocks; font is fitted to the width of each group span
  draw_group_annotations <- function() {
    if (!show_group_names) return(invisible())
    base_font <- 15      # same readable size as in the merged version (shrinks only if the span is too narrow)
    for (g in seq_along(group_labels)) {
      idx <- group_indices[[group_labels[g]]]
      tryCatch({
        seekViewport("global")
        seekViewport(paste0("annotation_empty_", idx[1]))
        loc1 <- deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
        seekViewport(paste0("annotation_empty_", idx[length(idx)]))
        loc2 <- deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
        seekViewport("global")
        w_mm <- convertWidth(loc2$x - loc1$x, "mm", valueOnly = TRUE)
        fs <- if (is.na(group_fontsize)) min(base_font, fit_font(group_labels[g], w_mm, 9, max_font = 15)) else group_fontsize
        # header box as wide as all colour boxes of the group, label centred in it
        grid.rect(x = (loc1$x + loc2$x) / 2, y = (loc1$y + loc2$y) / 2,
                  width = loc2$x - loc1$x, height = loc2$y - loc1$y, just = "center",
                  gp = gpar(fill = NA, col = "grey55", lwd = 0.8))
        grid.text(group_labels[g], x = (loc1$x + loc2$x) / 2, y = (loc1$y + loc2$y) / 2,
                  just = "center", gp = gpar(fontsize = fs, fontface = "bold", col = "black"))
      }, error = function(e) message("Group label '", group_labels[g], "' not drawn: ", conditionMessage(e)))
    }
  }

  list(top_annot_labels = top_annot_labels, col_split = col_split,
       gap_vect = gap_vect, draw_group_annotations = draw_group_annotations,
       legends = legends)
}

# Builds everything but does NOT draw (so Shiny can draw on screen or to file)
my_proteomics_heat_map_build <- function(input_object,
                                         sample_metadata,
                                         top_n = "all",
                                         merge = FALSE,
                                         cluster_rows = FALSE,
                                         show_row_dend = FALSE,
                                         dist_param = "euclidean",
                                         hclust_param = "average",
                                         scale_rows = TRUE,
                                         show_genes = FALSE,
                                         custom_order = NULL,
                                         group_colors,
                                         zscore_limits = c(-2, 2),
                                         color_low = "blue",
                                         color_mid = "white",
                                         color_high = "red",
                                         col_fontsize = NA,
                                         row_fontsize = 12,
                                         group_fontsize = NA,
                                         width = NULL,
                                         height = NULL,
                                         row_split = NULL,
                                         row_km = NULL,
                                         extra_cols = character(0),
                                         extra_colors = NULL,
                                         show_sample_names = TRUE,
                                         show_group_names = TRUE,
                                         leg_width = 10,
                                         leg_height = NULL,
                                         leg_font = 16) {

  main_matrix <- as.matrix(input_object)
  sample_metadata <- sample_metadata[colnames(main_matrix), , drop = FALSE]
  if (any(is.na(sample_metadata$Sample)))
    stop("Sample metadata does not fully match matrix column names.")

  if (merge) {
    main_matrix <- limma::avearrays(
      main_matrix,
      ID = factor(sample_metadata$Group, levels = unique(sample_metadata$Group)))
    sample_metadata <- sample_metadata[!duplicated(sample_metadata$Group), ]
    sample_metadata$Sample <- as.character(sample_metadata$Group)
    rownames(sample_metadata) <- sample_metadata$Sample
  }

  if (top_n == "all") top_n <- nrow(main_matrix)
  row_vars <- matrixStats::rowVars(main_matrix, na.rm = TRUE)
  keep_non_constant <- !is.na(row_vars) & row_vars > 0
  main_matrix <- main_matrix[keep_non_constant, , drop = FALSE]
  row_vars <- row_vars[keep_non_constant]
  top_n <- min(as.numeric(top_n), length(row_vars))
  keep_top <- rank(-row_vars, ties.method = "first") <= top_n
  main_matrix <- main_matrix[keep_top, , drop = FALSE]

  dist_mat <- as.matrix(dist(main_matrix, method = dist_param))
  na_rows <- rowSums(is.na(dist_mat)) > 0
  main_matrix <- main_matrix[!na_rows, , drop = FALSE]
  sample_metadata <- sample_metadata[colnames(main_matrix), , drop = FALSE]

  if (scale_rows) main_matrix <- t(scale(t(main_matrix), center = TRUE, scale = TRUE))
  main_matrix[!is.finite(main_matrix)] <- NA            # Inf / NaN -> NA

  # k-means cannot cope with NA / NaN / Inf: keep complete rows only (reported back to the app)
  dropped <- 0L
  if (!is.null(row_km)) {
    bad <- rowSums(!is.finite(main_matrix)) > 0
    dropped <- sum(bad)
    main_matrix <- main_matrix[!bad, , drop = FALSE]
    if (nrow(main_matrix) <= row_km)
      stop("K-means needs complete rows: only ", nrow(main_matrix), " proteins have no missing values. ",
           "Impute the data first or use fewer clusters.")
  }

  # ---- sizing (mm) ----
  longest_lab <- max(nchar(c(as.character(sample_metadata$Sample),
                             as.character(sample_metadata$Group))))
  .overhead_mm   <- 70
  .row_height_mm <- if (show_genes) 3.5 else 1.5
  # merged blocks must be wide enough to hold the group names
  .col_width_mm  <- if (merge) max(10, min(30, 3 + longest_lab * 1.6)) else 8
  extra_mm <- 8 * length(extra_cols)
  .height_mm <- if (!is.null(height)) height else
    max(120, min(.overhead_mm + extra_mm + nrow(main_matrix) * .row_height_mm, 500))
  gene_mm <- if (show_genes) max(nchar(rownames(main_matrix))) * row_fontsize * 0.3528 * 0.55 else 0
  .width_mm  <- if (!is.null(width)) width else
    max(150, min(60 + gene_mm + ncol(main_matrix) * .col_width_mm, 420))
  body_col_mm <- max(3, (.width_mm - 60 - gene_mm) / ncol(main_matrix))

  top_annot_results <- my_prepare_top_annotation_proteomics(
    sample_info = sample_metadata, group_colors = group_colors,
    custom_order = custom_order, col_fontsize = col_fontsize,
    group_fontsize = group_fontsize, merge = merge,
    extra_vars = if (length(extra_cols)) sample_metadata[, extra_cols, drop = FALSE] else NULL,
    extra_colors = extra_colors,
    show_sample_names = show_sample_names, show_group_names = show_group_names,
    col_width_mm = body_col_mm, legend_fontsize = leg_font)

  if (is.null(leg_height)) leg_height <- .height_mm * 0.12

  ht <- ComplexHeatmap::Heatmap(
    main_matrix,
    name = "Z-score",
    col = circlize::colorRamp2(breaks = c(zscore_limits[1], 0, zscore_limits[2]),
                               colors = c(color_low, color_mid, color_high)),
    na_col = "gray10",
    cluster_columns = FALSE,
    cluster_rows = cluster_rows,
    clustering_distance_rows = dist_param,
    clustering_method_rows = hclust_param,
    row_split = row_split,
    row_km = row_km,
    show_column_dend = FALSE,
    show_row_dend = show_row_dend,
    row_dend_width = if (show_row_dend) unit(15, "mm") else unit(0, "mm"),
    show_row_names = show_genes,
    row_names_side = "right",
    row_names_gp = gpar(fontsize = row_fontsize),
    top_annotation = top_annot_results$top_annot_labels,
    # no right-hand padding annotation: it pushed the gene labels away from the matrix
    column_split = top_annot_results$col_split,
    cluster_column_slices = TRUE,
    column_gap = unit(top_annot_results$gap_vect, "mm"),
    show_column_names = FALSE,
    border = FALSE,
    column_title = NULL,
    use_raster = FALSE,
    heatmap_legend_param = list(title_position = "topcenter",
                                legend_height = unit(leg_height, "mm"),
                                grid_width = unit(leg_width, "mm"),
                                title_gp = gpar(fontsize = leg_font, fontface = "bold"),
                                labels_gp = gpar(fontsize = leg_font)))

  list(ht = ht, main_matrix = main_matrix, annot = top_annot_results, dropped = dropped,
       width_mm = .width_mm, height_mm = .height_mm, top_n = top_n, merge = merge)
}

# Draw on the CURRENT device and return the drawn heatmap
draw_heatmap <- function(core) {
  ht_drawn <- ComplexHeatmap::draw(
    core$ht,
    annotation_legend_list = core$annot$legends,
    padding = unit(c(2, 2, 2, 6), "mm"))
  seekViewport("global")
  core$annot$draw_group_annotations()
  ht_drawn
}

extract_clusters <- function(ht_drawn, main_matrix) {
  extracted_order <- ComplexHeatmap::row_order(ht_drawn)
  cluster_assignment <- vector("integer", length = nrow(main_matrix))
  names(cluster_assignment) <- rownames(main_matrix)
  if (is.list(extracted_order)) {
    for (i in seq_along(extracted_order))
      cluster_assignment[rownames(main_matrix)[extracted_order[[i]]]] <- i
  } else {
    cluster_assignment[] <- 1
  }
  data.frame(Cluster = cluster_assignment, main_matrix,
             check.names = FALSE, stringsAsFactors = FALSE)
}
