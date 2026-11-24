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
                                                 col_width_mm = 8,       # size of ONE box along the heatmap (column width, or row height when sideways)
                                                 legend_fontsize = 14,
                                                 direction = "normal",
                                                 rot_names = FALSE) {        # TRUE: sample names written sideways (normal direction, fixed boxes)
  force(show_sample_names); force(show_group_names); force(col_fontsize); force(group_fontsize)
  side <- identical(direction, "sideways")

  sample_info$Group <- factor(sample_info$Group, levels = unique(sample_info$Group))

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

  # thickness of the boxes: 8 mm high (normal) or as wide as the longest label (sideways)
  rot <- side || isTRUE(rot_names)
  # without names the colour boxes are as thick as the group header (10 mm sideways, 8 mm on top)
  thick_of <- function(labs) if (!isTRUE(show_sample_names)) (if (side) 10 else 8)
                             else if (rot) max(12, min(45, 2 + 1.9 * max(nchar(labs), 1))) else 8
  block_h <- thick_of(block_labels)
  fit_box <- function(labs, max_font = 15) {
    th <- thick_of(labs)
    if (rot) fit_font(labs, th, col_width_mm, max_font = max_font)
    else     fit_font(labs, col_width_mm, th, max_font = max_font)
  }

  # Gaps between groups
  gap_vect <- c()
  if (merge) {
    gap_vect <- rep(1, length(levels(sample_info$Group)) - 1)
  } else {
    group_vec <- as.character(sample_info$Group)
    for (i in seq_len(length(group_vec) - 1))
      gap_vect <- c(gap_vect, if (group_vec[i] == group_vec[i + 1]) 0.5 else 1)
  }

  # font size: a number = used exactly; NA (blank in the UI) = auto-fit to the boxes
  box_font <- if (is.na(col_fontsize)) fit_box(block_labels) else col_fontsize

  # ---- group header font: one value for all groups; when the boxes are too narrow for horizontal text
  #      (small fixed cells) the names are written vertically in a taller header ----
  group_vec <- as.character(sample_info$Group)
  group_labels <- names(table(sample_info$Group))
  group_indices <- split(seq_along(group_vec), group_vec)
  spans <- if (merge) rep(1, length(group_labels)) else as.numeric(table(sample_info$Group)[group_labels])
  group_font_h <- min(15, vapply(seq_along(group_labels),
    function(g) fit_font(group_labels[g], spans[g] * col_width_mm, 9, max_font = 15), numeric(1)))
  # names stay horizontal (typed fonts always do); only when the auto font would be tiny are they written vertically,
  # in a header that never grows beyond 30 mm (the font shrinks instead)
  vertical <- !side && show_group_names && is.na(group_fontsize) && group_font_h < 4.5
  empty_mm <- 10
  group_font_auto <- group_font_h
  if (vertical) {
    longest_g <- max(nchar(group_labels), 1)
    group_font_auto <- min(12, (30 - 4) / (longest_g * 0.3528 * 0.6), vapply(seq_along(group_labels),
      function(g) fit_font(group_labels[g], Inf, spans[g] * col_width_mm, max_font = 12, min_font = 3), numeric(1)))
    group_font_auto <- max(3, group_font_auto)
    empty_mm <- max(10, min(30, longest_g * group_font_auto * 0.3528 * 0.6 + 4))
  }
  group_fs <- if (is.na(group_fontsize)) group_font_auto else group_fontsize

  blk <- function(gp, labels, fs, th) {
    if (side) anno_block(gp = gp, which = "row", width = unit(th, "mm"), labels = labels,
                         labels_gp = gpar(col = "black", fontsize = fs), labels_rot = 0)
    else      anno_block(gp = gp, which = "column", height = unit(th, "mm"), labels = labels,
                         labels_gp = gpar(col = "black", fontsize = fs), labels_rot = if (isTRUE(rot_names)) 90 else 0)
  }

  annot_args <- list(empty = if (side) anno_empty(which = "row", border = FALSE, width = unit(10, "mm"))
                              else     anno_empty(which = "column", border = FALSE, height = unit(empty_mm, "mm")))
  left_mm <- 10 + block_h
  if (merge) {
    annot_args$groups <- blk(fill_set_groups, NULL, box_font, block_h)
  } else {
    annot_args$samples <- blk(fill_set_groups, if (show_sample_names) block_labels else NULL, box_font, block_h)
  }

  legends <- list(Legend(title = "Group", at = names(group_colors),
                         legend_gp = gpar(fill = unname(group_colors)),
                         title_gp = gpar(fontsize = legend_fontsize, fontface = "bold"),
                         labels_gp = gpar(fontsize = legend_fontsize)))

  if (!is.null(extra_vars) && ncol(extra_vars) > 0) {
    for (v in colnames(extra_vars)) {
      vals <- as.character(extra_vars[[v]])
      th <- thick_of(vals)
      annot_args[[v]] <- blk(gpar(fill = unname(extra_colors[[v]][vals])),
                             if (show_sample_names) vals else NULL,
                             if (is.na(col_fontsize)) fit_box(vals) else col_fontsize, th)
      left_mm <- left_mm + th
      legends[[length(legends) + 1]] <- Legend(
        title = v, at = names(extra_colors[[v]]),
        legend_gp = gpar(fill = unname(extra_colors[[v]])),
        title_gp = gpar(fontsize = legend_fontsize, fontface = "bold"),
        labels_gp = gpar(fontsize = legend_fontsize))
    }
  }
  annot_obj <- if (side) do.call(rowAnnotation, annot_args) else do.call(HeatmapAnnotation, annot_args)

  # Group names above (or, sideways, left of) the blocks: one header box per group, drawn INSIDE the viewport of the
  # group's first "empty" annotation slice (ComplexHeatmap::decorate_annotation), so it always sits exactly on its samples
  draw_group_annotations <- function() {
    if (!show_group_names) return(invisible())
    for (g in seq_along(group_labels)) local({
      lab <- group_labels[g]
      idx <- group_indices[[lab]]
      inner <- if (length(idx) > 1) sum(gap_vect[idx[-length(idx)]]) else 0        # gaps between the group's own samples (mm)
      tryCatch(ComplexHeatmap::decorate_annotation("empty", slice = idx[1], {
        if (side) {
          gh <- length(idx) * convertHeight(unit(1, "npc"), "mm", valueOnly = TRUE) + inner
          grid.rect(x = unit(0.5, "npc"), y = unit(1, "npc"), width = unit(1, "npc"), height = unit(gh, "mm"),
                    just = "top", gp = gpar(fill = NA, col = "grey55", lwd = 0.8))
          grid.text(lab, x = unit(0.5, "npc"), y = unit(1, "npc") - unit(gh / 2, "mm"), rot = 90,
                    gp = gpar(fontsize = group_fs, fontface = "bold", col = "black"))
        } else {
          gw <- length(idx) * convertWidth(unit(1, "npc"), "mm", valueOnly = TRUE) + inner
          grid.rect(x = unit(0, "npc"), y = unit(0.5, "npc"), width = unit(gw, "mm"), height = unit(1, "npc"),
                    just = "left", gp = gpar(fill = NA, col = "grey55", lwd = 0.8))
          grid.text(lab, x = unit(gw / 2, "mm"), y = unit(0.5, "npc"), rot = if (vertical) 90 else 0,
                    gp = gpar(fontsize = group_fs, fontface = "bold", col = "black"))
        }
      }), error = function(e) message("Group label '", lab, "' not drawn: ", conditionMessage(e)))
    })
  }

  list(top_annot_labels = annot_obj, col_split = col_split,
       gap_vect = gap_vect, draw_group_annotations = draw_group_annotations,
       legends = legends, left_mm = left_mm, empty_mm = empty_mm, vertical = vertical,
       fonts = list(sample = box_font, group = group_font_auto))
}

# Builds everything but does NOT draw (so Shiny can draw on screen or to file)
#  direction: "normal" (proteins = rows) | "sideways" (transposed: samples = rows) | "circular" (circlize)
#  shape:     "rect" (fills the window) | "square" (fixed square boxes) | "circle" (fixed size, dots)
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
                                         row_fontsize = NA,
                                         group_fontsize = NA,
                                         width = NULL,               # typed canvas size (mm); NULL = fit the window
                                         height = NULL,
                                         win_w = NULL,               # the window (scroll area) the plot is fitted into (mm)
                                         win_h = NULL,
                                         row_split = NULL,
                                         row_km = NULL,
                                         extra_cols = character(0),
                                         extra_colors = NULL,
                                         show_sample_names = TRUE,
                                         show_group_names = TRUE,
                                         leg_width = 10,
                                         leg_height = NULL,
                                         leg_font = 16,
                                         direction = "normal",
                                         shape = "rect",
                                         cell_mm = NULL,             # NULL = auto-fit the box / circle size to the window
                                         circle_by_value = FALSE,
                                         circ_opts = list(),
                                         gene_side = NULL,            # normal: "right" | "left"; sideways: "down" | "up"
                                         show_group_legend = TRUE) {
  win_w <- win_w %||% 280; win_h <- win_h %||% 160

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

  col_fun <- circlize::colorRamp2(breaks = c(zscore_limits[1], 0, zscore_limits[2]),
                                  colors = c(color_low, color_mid, color_high))
  legend_title <- if (scale_rows) "Z-score" else "Value"

  # ---------------- circular: handled by draw_circular() ----------------
  if (identical(direction, "circular")) {
    nring <- ncol(main_matrix); n <- nrow(main_matrix)
    grp <- as.character(sample_metadata$Group)
    gap <- as.numeric(circ_opts$gap %||% 8); start <- as.numeric(circ_opts$start %||% 180)   # circlize insists on class "numeric"
    # canvas: the circle sits on the left, the Z-score bar on the right; everything is fitted into the window
    lg_mm <- leg_width + 8 + leg_font * 0.3528 * 3 + 12
    W0 <- width %||% win_w; H0 <- height %||% win_h
    side_mm <- max(40, min(H0, W0 - lg_mm))
    w <- if (!is.null(width)) width else side_mm + lg_mm
    h <- if (!is.null(height)) height else side_mm
    # radial budget: unit radius = 1; an empty centre, then one ring per sample (or per group when merged),
    # a small spacer between the samples of different groups
    hole <- 0.18
    nsp <- if (merge) 0L else max(0L, length(rle(grp)$lengths) - 1L)
    m1 <- min(0.004, 0.08 / (nring + nsp))
    avail <- 1 - hole - 0.03
    th <- (avail - (nring + nsp) * (m1 + 0.0001)) / (nring + 0.5 * nsp)
    th <- min(0.18, max(0.002, th))
    sp <- 0.5 * th
    # protein labels: auto font from the arc available per protein; the canvas is widened to hold them
    lab_chars <- max(nchar(rownames(main_matrix)), 1)
    ext <- if (show_genes) 1.2 else 1.04
    gene_font <- row_fontsize
    for (it in 1:3) {
      R_mm <- side_mm / (2 * ext)
      arc_mm <- 2 * pi * R_mm * (360 - gap) / 360 / n
      gene_font <- if (is.na(row_fontsize)) max(2, min(14, arc_mm * 0.75 / 0.3528)) else row_fontsize
      if (show_genes) {
        lab_mm <- min(0.22 * side_mm, lab_chars * gene_font * 0.3528 * 0.55)
        ext <- 1 / (1 - 2 * lab_mm / side_mm) + 0.03
      }
    }
    R_mm <- side_mm / (2 * ext)
    # group names live in the opening of the circle, level with the rings of their group: fit them there
    rl <- rle(if (merge) as.character(seq_len(nring)) else grp)
    ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1
    rpos <- numeric(nring); pos <- 1
    for (k in seq_len(nring)) {
      if (k > 1 && !merge && grp[k] != grp[k - 1]) pos <- pos - sp - (m1 + 0.0001)
      pos <- pos - (m1 + 0.0001) - th / 2; rpos[k] <- pos; pos <- pos - th / 2
    }
    g_names <- if (merge) colnames(main_matrix) else grp[starts]
    g_fit <- vapply(seq_along(starts), function(g) {
      rc <- mean(rpos[starts[g]:ends[g]])
      arc <- 2 * pi * rc * R_mm * gap / 360
      band <- (rpos[starts[g]] - rpos[ends[g]] + th) * R_mm
      fit_font(g_names[g], arc, band, max_font = 12, min_font = 3)
    }, numeric(1))
    ring_auto <- min(g_fit)
    ring_font <- if (is.na(group_fontsize)) ring_auto else group_fontsize
    if (is.null(leg_height)) leg_height <- side_mm * 0.45
    leg_height <- min(leg_height, h - leg_font * 0.3528 * 3 - 12)
    circ <- c(circ_opts[setdiff(names(circ_opts), c("gap", "start"))], list(gap = gap, start = start, ext = ext, th = th, sp = sp, m1 = m1, hole = hole,
                              gene_font = gene_font, ring_font = ring_font, show_genes = show_genes,
                              show_ring = isTRUE(show_group_names), merge = merge, groups = if (merge) colnames(main_matrix) else grp,
                              side_mm = side_mm, w_mm = w, h_mm = h))
    return(list(circular = TRUE, ht = NULL, main_matrix = main_matrix, col_fun = col_fun,
                circ = circ, annot = NULL, dropped = 0L,
                fonts = list(sample = NA_real_, group = ring_auto, gene = gene_font),
                width_mm = w, height_mm = h, body_h_mm = side_mm * 0.45, top_n = top_n, merge = merge,
                direction = "circular", legend_title = legend_title, fixed = FALSE, cell_mm = NULL,
                leg = list(width = leg_width, height = leg_height, font = leg_font)))
  }

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

  # ---- sizing (mm): EVERYTHING is fitted into the window (win_w x win_h) unless a size is typed ----
  side <- identical(direction, "sideways")
  fixed <- shape %in% c("square", "circle")
  nr <- nrow(main_matrix); nc <- ncol(main_matrix)
  n_extra <- length(extra_cols)
  dend_mm <- if (isTRUE(cluster_rows) && isTRUE(show_row_dend)) 15 else 0
  n_split <- if (!is.null(row_km)) row_km else if (!is.null(row_split)) row_split else 1
  split_gap <- max(0, n_split - 1)
  lab_chars <- max(nchar(c(names(group_colors), unlist(lapply(extra_colors, names)), extra_cols, "Group")), 1)
  leg_est <- list(z = leg_width + 8 + leg_font * 0.3528 * 3 + 4,
                  g = if (isTRUE(show_group_legend)) 8 + lab_chars * leg_font * 0.3528 * 0.55 + 6 else 0)
  legend_mm <- leg_est$z + leg_est$g + (if (isTRUE(show_group_legend)) 8 else 0) + 8      # width of the legend column(s)
  W0 <- width %||% win_w; H0 <- height %||% win_h
  extra_df <- if (length(extra_cols)) sample_metadata[, extra_cols, drop = FALSE] else NULL
  mk_annot <- function(cl) my_prepare_top_annotation_proteomics(
    sample_info = sample_metadata, group_colors = group_colors,
    custom_order = custom_order, col_fontsize = col_fontsize,
    group_fontsize = group_fontsize, merge = merge,
    extra_vars = extra_df, extra_colors = extra_colors,
    show_sample_names = show_sample_names, show_group_names = show_group_names,
    col_width_mm = cl, legend_fontsize = leg_font, direction = direction,
    rot_names = fixed && !side)

  cell <- if (!is.null(cell_mm)) cell_mm else max(0.3, min(W0 / nc, H0 / nr, 12))      # starting guess, refined below
  gene_font_est <- if (!is.na(row_fontsize)) row_fontsize else 8
  for (it in 1:5) {
    ta <- mk_annot(cell)
    gaps_mm <- sum(ta$gap_vect)
    gene_mm <- if (show_genes) max(nchar(rownames(main_matrix))) * gene_font_est * 0.3528 * 0.55 + 2 else 0
    if (!side) {
      hdr <- ta$empty_mm + (ta$left_mm - 10) + (1 + n_extra) + 6
      wav <- W0 - legend_mm - gene_mm - dend_mm - 6 - gaps_mm                  # room for the cells, along the samples
      hav <- H0 - hdr - split_gap                                              # ... along the proteins
      s_samp <- wav / nc; s_prot <- hav / nr
    } else {
      wav <- W0 - ta$left_mm - legend_mm - 12 - split_gap                      # proteins
      hav <- H0 - 10 - gene_mm - dend_mm - 4 - gaps_mm                         # samples
      s_samp <- hav / nc; s_prot <- wav / nr
    }
    s_samp <- max(0.2, s_samp); s_prot <- max(0.05, s_prot)
    cell_new <- if (!fixed) s_samp else if (!is.null(cell_mm)) cell_mm else max(0.05, min(s_samp, s_prot))
    gene_cell <- if (fixed) cell_new else s_prot
    gene_font_est <- if (is.na(row_fontsize)) fit_font(rownames(main_matrix), Inf, gene_cell, max_font = 14, min_font = 3) else row_fontsize
    done <- abs(cell_new - cell) < 0.02
    cell <- cell_new
    if (done) break
  }
  top_annot_results <- ta
  gene_font <- gene_font_est
  # safety net: an auto-fitted plot must fit the window completely (never scroll)
  needs <- function(cl) if (!side) c(w = nc * cl + gaps_mm + gene_mm + dend_mm + legend_mm + 6, h = nr * cl + split_gap + hdr)
                        else       c(w = nr * cl + split_gap + ta$left_mm + legend_mm + 12, h = nc * cl + gaps_mm + 10 + gene_mm + dend_mm + 4)
  if (fixed && is.null(cell_mm)) for (k in 1:6) {
    nd <- needs(cell); over <- c(nd[["w"]] - W0, nd[["h"]] - H0)
    if (all(over <= 0.5)) break
    cell <- max(0.05, cell - max(over / (if (!side) c(nc, nr) else c(nr, nc))) * 1.02)
  }

  if (!side) {
    body_w <- if (fixed) nc * cell + gaps_mm else wav + gaps_mm
    body_h_mm <- if (fixed) nr * cell + split_gap else hav + split_gap
    need_w <- body_w + gene_mm + dend_mm + legend_mm + 6
    need_h <- body_h_mm + hdr
  } else {
    body_w <- if (fixed) nr * cell + split_gap else wav + split_gap
    body_h_mm <- if (fixed) nc * cell + gaps_mm else hav + gaps_mm
    need_w <- body_w + ta$left_mm + legend_mm + 12
    need_h <- body_h_mm + 10 + gene_mm + dend_mm + 4
  }
  # canvas: fitted shapes wrap the plot tightly (never bigger than the window unless the user asked for a bigger box);
  # rectangles fill the window; a typed width / height always wins (a typed box size is never cropped)
  if (fixed) {
    .width_mm <- if (!is.null(width)) (if (is.null(cell_mm)) width else max(width, need_w)) else need_w
    .height_mm <- if (!is.null(height)) (if (is.null(cell_mm)) height else max(height, need_h)) else need_h
  } else {
    .width_mm <- W0; .height_mm <- H0
  }
  # the Z-score bar must stay inside the canvas (title + padding take ~ 3 lines of the legend font)
  leg_cap <- max(15, .height_mm - leg_font * 0.3528 * 3 - 14)
  leg_height <- min(if (is.null(leg_height)) body_h_mm else leg_height, leg_cap)

  # ---- legends: Z-score bar and the group / extra-variable legends side by side ----
  zleg <- Legend(title = legend_title, col_fun = col_fun,
                 legend_height = unit(leg_height, "mm"), grid_width = unit(leg_width, "mm"),
                 title_position = "topcenter",
                 title_gp = gpar(fontsize = leg_font, fontface = "bold"),
                 labels_gp = gpar(fontsize = leg_font))
  grp_pack <- if (isTRUE(show_group_legend))
    packLegend(list = top_annot_results$legends, direction = "vertical", gap = unit(4, "mm")) else NULL
  legend_obj <- list(z = zleg, grp = grp_pack)

  # ---- cell drawing ----
  mm_in <- if (side) t(main_matrix) else main_matrix
  zmax <- suppressWarnings(min(max(abs(zscore_limits)), max(abs(mm_in), na.rm = TRUE)))
  if (!is.finite(zmax) || zmax <= 0) zmax <- max(abs(zscore_limits))
  cell_fun <- NULL
  rect_gp <- gpar(col = NA)
  if (identical(shape, "square")) rect_gp <- gpar(col = "white", lwd = 0.6)
  if (identical(shape, "circle")) {
    rect_gp <- gpar(type = "none")
    cell_d <- cell                                       # the cell is cell x cell mm: the circle FILLS it (as in Morpheus)
    cell_fun <- function(j, i, x, y, w, h, fill) {
      v <- mm_in[i, j]
      # "reflects |value|": the diameter fraction runs linearly from 20 % to 100 % of the cell
      f <- if (isTRUE(circle_by_value) && is.finite(v)) 0.2 + 0.8 * min(1, abs(v) / zmax) else 1
      grid.circle(x, y, r = unit(cell_d / 2 * f * 0.96, "mm"), gp = gpar(fill = fill, col = NA))
    }
  }

  common_args <- list(
    matrix = mm_in, name = "Z-score", col = col_fun, na_col = "gray10",
    show_heatmap_legend = FALSE, border = FALSE, use_raster = FALSE,
    rect_gp = rect_gp, cell_fun = cell_fun)

  if (!side) {
    ht_args <- c(common_args, list(
      cluster_columns = FALSE, cluster_rows = cluster_rows,
      clustering_distance_rows = dist_param, clustering_method_rows = hclust_param,
      row_split = row_split, row_km = row_km,
      show_column_dend = FALSE, show_row_dend = show_row_dend,
      row_dend_width = if (show_row_dend) unit(15, "mm") else unit(0, "mm"),
      show_row_names = show_genes, row_names_side = if (identical(gene_side, "left")) "left" else "right",
      row_names_gp = gpar(fontsize = gene_font),
      top_annotation = top_annot_results$top_annot_labels,
      column_split = top_annot_results$col_split,
      cluster_column_slices = TRUE,
      column_gap = unit(top_annot_results$gap_vect, "mm"),
      show_column_names = FALSE, column_title = NULL))
    if (fixed) ht_args$width  <- unit(body_w, "mm")
    if (fixed) ht_args$height <- unit(body_h_mm, "mm")
  } else {
    ht_args <- c(common_args, list(
      cluster_rows = FALSE, cluster_columns = cluster_rows,
      clustering_distance_columns = dist_param, clustering_method_columns = hclust_param,
      column_split = row_split, column_km = if (is.null(row_km)) 1 else row_km,
      show_row_dend = FALSE, show_column_dend = show_row_dend,
      column_dend_height = if (show_row_dend) unit(15, "mm") else unit(0, "mm"),
      show_row_names = FALSE,
      show_column_names = show_genes, column_names_side = if (identical(gene_side, "up")) "top" else "bottom", column_names_rot = 90,
      column_names_gp = gpar(fontsize = gene_font),
      left_annotation = top_annot_results$top_annot_labels,
      row_split = top_annot_results$col_split, cluster_row_slices = FALSE,
      row_gap = unit(top_annot_results$gap_vect, "mm"),
      row_title = NULL, column_title = NULL))
    if (fixed) ht_args$width  <- unit(body_w, "mm")
    if (fixed) ht_args$height <- unit(body_h_mm, "mm")
  }
  ht <- do.call(ComplexHeatmap::Heatmap, ht_args)

  list(ht = ht, main_matrix = main_matrix, annot = top_annot_results, dropped = dropped,
       fonts = c(top_annot_results$fonts, list(gene = gene_font)), legend = legend_obj, leg_est = leg_est,
       width_mm = .width_mm, height_mm = .height_mm, body_h_mm = body_h_mm, leg_h_mm = leg_height + leg_font * 0.3528 * 3,
       fixed = fixed, cell_mm = if (fixed) cell else NULL, top_n = top_n, merge = merge, direction = direction, circular = FALSE)
}

# ---------------- circular heatmap (circlize, after "circular heatmap.R") ----------------
# one ring per sample (or per group when replicates are merged); sample names are not written - the GROUP names sit in the
# opening of the circle, level with the rings of their group; a small spacer separates the samples of different groups.
# If circlize reports "not enough space for cells" the rings are thinned and the drawing is retried.
draw_circular <- function(core) {
  m <- core$main_matrix
  o <- core$circ
  grp <- o$groups
  if (isTRUE(o$cluster)) {
    ok <- rowSums(!is.finite(m)) == 0
    if (sum(ok) > 2) {
      m <- m[ok, , drop = FALSE]
      hc <- hclust(dist(m, method = o$dist %||% "euclidean"), method = o$linkage %||% "average")
      m <- m[hc$order, , drop = FALSE]
    }
  }
  m <- m[nrow(m):1, , drop = FALSE]                 # as in the original script: first protein first, clockwise
  on.exit(circlize::circos.clear(), add = TRUE)
  gap <- o$gap; ring_cex <- o$ring_font / 12; gene_cex <- o$gene_font / 12
  nring <- ncol(m)
  # the circle gets a square figure region on the left of the device; the Z-score bar is drawn right next to it
  xlim <- ylim <- c(-o$ext, o$ext)
  par(fig = c(0, min(1, o$side_mm / o$w_mm), 0, 1))
  on.exit(par(fig = c(0, 1, 0, 1)), add = TRUE)

  attempt <- function(sc) {
    circlize::circos.clear()
    th <- o$th * sc; sp <- o$sp * sc
    circlize::circos.par(start.degree = o$start, gap.degree = gap,
                         track.margin = c(o$m1 * sc, 0.0001), points.overflow.warning = FALSE,
                         canvas.xlim = xlim, canvas.ylim = ylim)
    circlize::circos.heatmap.initialize(m, cluster = FALSE)
    tr <- integer(nring)
    for (k in seq_len(nring)) {
      if (k > 1 && !isTRUE(o$merge) && grp[k] != grp[k - 1])                       # spacer between groups
        circlize::circos.track(ylim = c(0, 1), track.height = sp, bg.border = NA)
      if (k == 1 && isTRUE(o$show_genes))
        circlize::circos.heatmap(m[, k, drop = FALSE], col = core$col_fun, na.col = "gray", track.height = th,
                                 rownames.side = "outside", rownames.cex = gene_cex)
      else
        circlize::circos.heatmap(m[, k, drop = FALSE], col = core$col_fun, na.col = "gray", track.height = th)
      tr[k] <- circlize::get.current.track.index()
    }
    if (isTRUE(o$show_ring)) {
      rl <- rle(as.character(grp)); ends <- cumsum(rl$lengths); starts <- ends - rl$lengths + 1
      for (g in seq_along(starts)) local({
        ng <- rl$lengths[g]; lab <- rl$values[g]
        mid <- starts[g] + (if (ng %% 2 == 1) (ng - 1) / 2 else ng / 2 - 1)        # the ring in the middle of the group
        at_edge <- ng %% 2 == 0                                                     # even: centre between two rings
        circlize::circos.track(track.index = tr[mid], panel.fun = function(x, y) {
          xl <- CELL_META$xlim; yl <- CELL_META$ylim
          x_gap <- (gap / 2) / ((360 - gap) / diff(xl))        # centre of the opening, in data units
          circlize::circos.text(xl[1] - x_gap, if (at_edge) yl[1] else mean(yl), lab, cex = ring_cex, facing = "reverse.clockwise")
        }, bg.border = NA)
      })
    }
    invisible(TRUE)
  }
  done <- FALSE
  for (sc in c(1, 0.85, 0.7, 0.55, 0.4, 0.25, 0.12)) {
    done <- tryCatch(attempt(sc), error = function(e) FALSE)
    if (isTRUE(done)) break
  }
  if (!isTRUE(done)) {
    circlize::circos.clear()
    plot.new(); text(0.5, 0.5, "The circular heatmap does not fit: use fewer samples,\nor increase the figure size.", cex = 1.1)
    return(invisible(NULL))
  }
  lg <- core$leg
  try({
    lgd <- ComplexHeatmap::Legend(title = core$legend_title, col_fun = core$col_fun,
                                  legend_height = unit(lg$height, "mm"), grid_width = unit(lg$width, "mm"),
                                  title_position = "topcenter",
                                  title_gp = gpar(fontsize = lg$font, fontface = "bold"),
                                  labels_gp = gpar(fontsize = lg$font))
    ComplexHeatmap::draw(lgd, x = unit(o$side_mm + 2, "mm"), just = "left")
  }, silent = TRUE)
  invisible(NULL)
}

# vertical centre (inches from the bottom of the device) of the heatmap body
body_centre_in <- function(ht_drawn) {
  ro <- ComplexHeatmap::row_order(ht_drawn)
  nrs <- if (is.list(ro)) length(ro) else 1
  yy <- function(vp, npc) {
    seekViewport("global"); seekViewport(vp)
    convertY(deviceLoc(x = unit(0, "npc"), y = unit(npc, "npc"))$y, "inches", valueOnly = TRUE)
  }
  top <- yy("heatmap_Z-score_1_1", 1)
  bot <- yy(paste0("heatmap_Z-score_", nrs, "_1"), 0)
  seekViewport("global")
  mean(c(top, bot))
}

# Draw on the CURRENT device and return the drawn heatmap (NULL for the circular plot).
# The Z-score bar and the group legend are drawn as two separate columns at the right edge,
# both centred vertically on the heatmap body.
draw_heatmap <- function(core) {
  if (isTRUE(core$circular)) return(draw_circular(core))
  zl <- core$legend$z; gl <- core$legend$grp
  # legend widths: ComplexHeatmap's width() generic is not exported -> take it from the namespace, estimate if that fails
  wmm <- function(l, est) tryCatch(
    convertWidth(utils::getFromNamespace("width", "ComplexHeatmap")(l), "mm", valueOnly = TRUE),
    error = function(e) est)
  wz <- wmm(zl, core$leg_est$z); wg <- if (is.null(gl)) 0 else wmm(gl, core$leg_est$g)
  gap_mm <- if (is.null(gl)) 0 else 8
  rpad <- wz + wg + gap_mm + 8
  ht_drawn <- ComplexHeatmap::draw(
    core$ht, show_heatmap_legend = FALSE, show_annotation_legend = FALSE,
    padding = unit(c(2, 2, 2, rpad), "mm"))
  yc_in <- tryCatch(body_centre_in(ht_drawn), error = function(e) NA_real_)
  seekViewport("global")
  dev_h <- convertHeight(unit(1, "npc"), "inches", valueOnly = TRUE)
  half <- (core$leg_h_mm %||% 0) / 25.4 / 2 + 0.2                  # half the Z-score legend (bar + title), in inches
  if (!is.finite(yc_in)) yc_in <- dev_h / 2
  yc <- unit(min(max(yc_in, half), max(half, dev_h - half)), "inches")
  if (!is.null(gl)) ComplexHeatmap::draw(gl, x = unit(1, "npc") - unit(3, "mm"), y = yc, just = c("right", "center"))
  ComplexHeatmap::draw(zl, x = unit(1, "npc") - unit(3 + wg + gap_mm, "mm"), y = yc, just = c("right", "center"))
  seekViewport("global")
  core$annot$draw_group_annotations()
  ht_drawn
}

extract_clusters <- function(ht_drawn, main_matrix, direction = "normal") {
  if (is.null(ht_drawn))
    return(data.frame(Cluster = 1L, main_matrix, check.names = FALSE, stringsAsFactors = FALSE,
                      row.names = rownames(main_matrix)))
  extracted_order <- if (identical(direction, "sideways")) ComplexHeatmap::column_order(ht_drawn)
                     else ComplexHeatmap::row_order(ht_drawn)
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
