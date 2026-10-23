#' Protein Detection Venn Diagram
#'
#' Draws a Venn diagram of detected proteins across groups (2–7 groups) using
#' \pkg{ggVennDiagram}'s low-level API with granular text sizing and subtractive blending.
#'
#' @param input_matrix Numeric matrix, features x samples (log2). \code{NA}
#'   indicates a missing/undetected value.
#' @param sample_info Data frame with columns \code{"Sample"} and \code{"Group"}.
#' @param threshold Minimum fraction of samples in a group that must have a
#'   non-NA value for a protein to count as detected in that group. Default: \code{0.7}.
#' @param group_colors Optional named character vector of hex colors, one per group.
#' @param group_labels Optional display labels. Underscores replaced with spaces if \code{NULL}.
#' @param title Plot title. Default: \code{"Protein detection overlap across groups"}.
#' @param subtitle Optional custom subtitle string. Default: \code{NULL}.
#' @param font_scale Numerical factor to scale all text elements globally if individual sizes are NULL. Default: \code{1.0}.
#' @param title_size Explicit font size for the plot title (in pt). Overrides \code{font_scale} if provided.
#' @param subtitle_size Explicit font size for the plot subtitle (in pt). Overrides \code{font_scale} if provided.
#' @param group_label_size Explicit font size for group display labels (in \code{geom_text} size units). Overrides \code{font_scale} if provided.
#' @param count_size Explicit font size for intersection count numbers (in \code{geom_text} size units). Overrides \code{font_scale} if provided.
#' @param count_fontface Font face for numbers (\code{"plain"}, \code{"bold"}, \code{"italic"}). Default: \code{"plain"}.
#' @param expand_limits Viewport expansion factor. Default: \code{1.0}.
#' @param fill_alpha Transparency/intensity factor for section fills (0.1 to 1.0). Default: \code{0.45}.
#' @param filename Output filename (PNG).
#' @param width,height Base plot dimensions in inches. Default: 8 x 7.
#' @param dpi Resolution. Default: 300.
#'
#' @return The ggplot object, invisibly. The plot is saved to \code{filename}.
#'
#' @export
venn_of_detected <- function(input_matrix,
                             sample_info,
                             threshold        = 0.7,
                             group_colors     = NULL,
                             group_labels     = NULL,
                             title            = "Protein detection overlap across groups",
                             subtitle         = NULL,
                             font_scale       = 1.0,
                             title_size       = NULL,
                             subtitle_size    = NULL,
                             group_label_size = NULL,
                             count_size       = NULL,
                             count_fontface   = "plain",
                             expand_limits    = 1.5,
                             fill_alpha       = 0.45,
                             filename         = "venn_diagram_of_detected_proteins.png",
                             width            = 8,
                             height           = 7,
                             dpi              = 300) {

  # ── Dependencies ────────────────────────────────────────────────────────────
  for (pkg in c("ggplot2", "ggVennDiagram")) {
    if (!requireNamespace(pkg, quietly = TRUE))
      stop(pkg, " is required: install.packages('", pkg, "')")
  }

  # ── Validation ───────────────────────────────────────────────────────────────
  if (!all(c("Sample", "Group") %in% colnames(sample_info)))
    stop("sample_info must contain columns 'Sample' and 'Group'.")
  missing_samp <- setdiff(sample_info$Sample, colnames(input_matrix))
  if (length(missing_samp) > 0)
    stop("Samples in sample_info but not in input_matrix: ",
         paste(missing_samp, collapse = ", "))

  groups <- unique(sample_info$Group)
  n_grps <- length(groups)
  if (n_grps < 2) stop("Need at least 2 groups.")
  if (n_grps > 7) stop("Maximum 7 groups supported.")
  if (!is.null(group_labels) && length(group_labels) != n_grps)
    stop("group_labels length (", length(group_labels),
         ") must equal number of groups (", n_grps, ").")

  # ── Colors ───────────────────────────────────────────────────────────────────
  default_palette <- c("#74c476", "#8884d8", "#bc80bd", "#fdae6b",
                       "#4393c3", "#d6604d", "#8c6d31")
  fill_colors <- setNames(default_palette[seq_len(n_grps)], groups)
  if (!is.null(group_colors)) {
    matched <- intersect(names(group_colors), groups)
    fill_colors[matched] <- group_colors[matched]
  }

  # ── Protein lists ─────────────────────────────────────────────────────────────
  group_samples  <- split(sample_info$Sample, sample_info$Group)
  group_proteins <- lapply(groups, function(grp) {
    samps <- group_samples[[grp]]
    sub   <- input_matrix[, samps, drop = FALSE]
    rownames(sub)[rowSums(!is.na(sub)) / length(samps) >= threshold]
  })

  display_names <- if (!is.null(group_labels)) group_labels else gsub("_", " ", groups)
  # internal unique names (display titles may be blank or duplicated)
  internal_names        <- paste0("G", seq_len(n_grps))
  names(group_proteins) <- internal_names
  names(fill_colors)    <- internal_names

  # ── Subtitle logic & Font Sizing ─────────────────────────────────────────────
  subtitle_text <- if (!is.null(subtitle)) {
    subtitle
  } else {
    paste0("Proteins detected in \u2265", threshold * 100, "% of samples per group")
  }

  # Calculate fallback default font sizing based on font_scale and expand_limits
  effective_font <- font_scale * max(1, expand_limits * 0.7)

  final_title_size       <- if (!is.null(title_size))       title_size       else 13 * effective_font
  final_subtitle_size    <- if (!is.null(subtitle_size))    subtitle_size    else 10 * effective_font
  final_group_label_size <- if (!is.null(group_label_size)) group_label_size else 5.2 * effective_font
  final_count_size       <- if (!is.null(count_size))       count_size       else 4.5 * effective_font

  # ── Shared Theme & Labels ────────────────────────────────────────────────────
  shared_labs <- ggplot2::labs(
    title    = title,
    subtitle = subtitle_text
  )

  shared_theme <- ggplot2::theme(
    plot.title       = ggplot2::element_text(hjust = 0.5, size = final_title_size, face = "bold"),
    plot.subtitle    = ggplot2::element_text(hjust = 0.5, color = "grey50", size = final_subtitle_size),
    plot.background  = ggplot2::element_rect(fill = "white", color = NA),
    panel.background = ggplot2::element_rect(fill = "white", color = NA)
  )

  # ── Process Geometry via ggVennDiagram ───────────────────────────────────────
  venn_obj  <- ggVennDiagram::Venn(group_proteins)
  venn_data <- ggVennDiagram::process_data(venn_obj)

  # Subtractive multiplicative color blending
  blend_colors_subtractive <- function(hex_vec, alpha = 0.45) {
    rgb_mat <- t(grDevices::col2rgb(hex_vec)) / 255
    mult_rgb <- apply(rgb_mat, 2, prod)
    final_rgb <- (1 - alpha) + alpha * mult_rgb
    grDevices::rgb(final_rgb[1], final_rgb[2], final_rgb[3])
  }

  region_edge_df  <- ggVennDiagram::venn_regionedge(venn_data)
  set_edge_df     <- ggVennDiagram::venn_setedge(venn_data)
  set_label_df    <- ggVennDiagram::venn_setlabel(venn_data)
  region_label_df <- ggVennDiagram::venn_regionlabel(venn_data)

  # two sets: put the circles side by side (the engine can return them stacked vertically)
  if (n_grps == 2) {
    cen <- stats::aggregate(cbind(X, Y) ~ id, data = set_edge_df, FUN = mean)
    if (abs(diff(cen$Y)) > abs(diff(cen$X))) {
      swapxy <- function(d) { t <- d$X; d$X <- d$Y; d$Y <- t; d }
      region_edge_df <- swapxy(region_edge_df); set_edge_df  <- swapxy(set_edge_df)
      set_label_df   <- swapxy(set_label_df);   region_label_df <- swapxy(region_label_df)
    }
  }

  # titles: custom text, placed ABOVE each circle
  set_label_df$name <- display_names[match(set_label_df$name, internal_names)]
  y_span <- diff(range(set_edge_df$Y))
  if (n_grps == 2) {      # 3+ circles keep the original symmetric title positions (top and bottom)
    for (i in seq_len(nrow(set_label_df))) {
      edge <- set_edge_df[set_edge_df$id == set_label_df$id[i], ]
      set_label_df$X[i] <- mean(range(edge$X))
      set_label_df$Y[i] <- max(edge$Y) + 0.07 * y_span
    }
  }

  unique_ids <- unique(region_edge_df$id)
  color_map  <- setNames(
    sapply(unique_ids, function(rid) {
      idx <- as.integer(strsplit(rid, "/")[[1]])
      grp_colors <- fill_colors[idx]
      blend_colors_subtractive(grp_colors, alpha = fill_alpha)
    }),
    unique_ids
  )

  region_edge_df$computed_fill <- color_map[region_edge_df$id]
  set_edge_df$set_color        <- fill_colors[as.integer(set_edge_df$id)]
  set_label_df$set_color       <- fill_colors[as.integer(set_label_df$id)]

  x_rng <- range(c(region_edge_df$X, set_edge_df$X, set_label_df$X))
  y_rng <- if (n_grps == 2) range(c(region_edge_df$Y, set_edge_df$Y, set_label_df$Y + 0.06 * y_span)) else range(c(region_edge_df$Y, set_edge_df$Y, set_label_df$Y))

  x_mid <- mean(x_rng)
  y_mid <- mean(y_rng)
  x_half <- (diff(x_rng) / 2) * expand_limits
  y_half <- (diff(y_rng) / 2) * expand_limits

  # ── Plot Construction ────────────────────────────────────────────────────────
  p <- ggplot2::ggplot() +
    ggplot2::geom_polygon(
      ggplot2::aes(X, Y, fill = computed_fill, group = id),
      data  = region_edge_df,
      color = NA
    ) +
    ggplot2::scale_fill_identity(guide = "none") +
    ggplot2::geom_path(
      ggplot2::aes(X, Y, color = set_color, group = id),
      data        = set_edge_df,
      linewidth   = 1.2 * max(1, expand_limits * 0.5),
      show.legend = FALSE
    ) +
    ggplot2::scale_color_identity(guide = "none") +
    # Group Labels
    ggplot2::geom_text(
      ggplot2::aes(X, Y, label = name, color = set_color),
      data        = set_label_df,
      vjust       = if (n_grps == 2) 0 else 0.5,
      size        = final_group_label_size,
      fontface    = "bold",
      show.legend = FALSE
    ) +
    # Intersection Count Numbers
    ggplot2::geom_text(
      ggplot2::aes(X, Y, label = count),
      data     = region_label_df,
      size     = final_count_size,
      color    = "black",
      fontface = count_fontface
    ) +
    ggplot2::coord_fixed(
      xlim = c(x_mid - x_half, x_mid + x_half),
      ylim = c(y_mid - y_half, y_mid + y_half),
      clip = "off"
    ) +
    ggplot2::theme_void() +
    shared_labs + shared_theme

  # ── Dynamic Canvas Cropping ─────────────────────────────────────────────────
  export_width  <- width / sqrt(expand_limits)
  export_height <- height / sqrt(expand_limits)

  ggplot2::ggsave(filename = filename, plot = p,
                  width = export_width, height = export_height, dpi = dpi, bg = "white")
  message("Venn diagram saved to: ", normalizePath(filename, mustWork = FALSE))
  invisible(p)
}
