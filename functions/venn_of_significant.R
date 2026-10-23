#' Venn Diagram of Significant DE Proteins Across Contrasts
#'
#' Draws a Venn diagram showing the overlap of differentially expressed
#' proteins across selected contrasts, using a \code{decideTests}-style
#' matrix as input. All text sizes and canvas dimensions are fully
#' controllable via parameters.
#'
#' @param decide_matrix A matrix of DE calls with proteins as rows and
#'   contrasts as columns. Values should be \code{1} (up), \code{-1}
#'   (down), or \code{0} (not significant), as produced by
#'   \code{limma::decideTests()}.
#' @param contrasts Character vector of contrast names to include. Must
#'   match column names of \code{decide_matrix}. Default: all columns.
#' @param direction One of \code{"either"}, \code{"up"}, or \code{"down"}.
#'   Controls which proteins are included in each set.
#' @param contrast_colors Character vector of colors, one per contrast
#'   (applied positionally). Must be at least as long as \code{contrasts}.
#' @param contrast_labels Optional display labels for the contrasts.
#'   Must have the same length as \code{contrasts}. If \code{NULL},
#'   underscores in contrast names are replaced with spaces.
#' @param text_size Numeric. Size of the count labels inside each Venn
#'   region. Default: \code{6.5}.
#' @param set_name_size Numeric. Size of the contrast name labels placed
#'   outside the circles. Default: \code{7.5}.
#' @param title_size Numeric. Font size of the plot title. Default: \code{18}.
#' @param subtitle_size Numeric. Font size of the plot subtitle.
#'   Default: \code{13}.
#' @param caption_size Numeric. Font size of the total-proteins caption
#'   in the bottom-right corner. Default: \code{6}.
#' @param label_line_height Numeric. Extra space in points added to the
#'   bottom margin of the subtitle per additional line in multi-line set
#'   labels (those containing \code{\\n}). This pushes the circles downward,
#'   giving the set labels room between the subtitle and the circle tops.
#'   Increase if labels are still clipped; decrease if there is too much
#'   gap. Default: \code{15}.
#' @param filename Output filename. Default: \code{"venn_DE_proteins.png"}.
#'
#' @return The ggplot object invisibly. The diagram is also saved to
#'   \code{filename}.
#'
#' @examples
#' \dontrun{
#' venn_of_significant(
#'   decide_matrix   = decideTests(fit2),
#'   contrasts       = c("KO_vs_WT", "KO_Stress_vs_KO"),
#'   direction       = "either",
#'   contrast_colors = c("#1f78b4", "#33a02c"),
#'   text_size       = 7,
#'   set_name_size   = 8,
#'   title_size      = 16,
#'   subtitle_size   = 11,
#'   width           = 9,
#'   height          = 8
#' )
#' }
#'
#' @export
venn_of_significant <- function(decide_matrix,
                                contrasts       = NULL,
                                direction       = c("either", "up", "down"),
                                contrast_colors = c("#1f78b4", "#33a02c"),
                                contrast_labels = NULL,
                                # ── Text sizes ──────────────────────────────
                                text_size       = 6.5,   # counts inside regions
                                set_name_size   = 7.5,   # contrast name labels
                                title_size      = 18,    # plot title
                                subtitle_size   = 13,    # plot subtitle
                                caption_size    = 6,     # total proteins note
                                label_line_height = 15,  # extra pts of subtitle bottom margin per line
                                filename        = "venn_DE_proteins.png") {

  if (!requireNamespace("ggvenn",   quietly = TRUE))
    stop("ggvenn is required: install.packages('ggvenn')")
  if (!requireNamespace("ggplot2",  quietly = TRUE))
    stop("ggplot2 is required: install.packages('ggplot2')")

  direction <- match.arg(direction)

  # Default to all contrasts if not specified
  if (is.null(contrasts)) contrasts <- colnames(decide_matrix)
  if (length(contrasts) > 3)
    warning("Venn diagrams with >3 sets are hard to read. ",
            "Consider subsetting contrasts.")

  if (length(contrast_colors) < length(contrasts))
    stop("contrast_colors must have at least as many colors as contrasts.")

  if (!is.null(contrast_labels) && length(contrast_labels) != length(contrasts))
    stop("contrast_labels must have the same length as contrasts.")

  # Subset matrix to requested contrasts
  mat <- as.matrix(decide_matrix)[, contrasts, drop = FALSE]

  # Build per-contrast protein lists based on direction
  protein_lists <- lapply(contrasts, function(ct) {
    vals <- mat[, ct]
    switch(direction,
           "either" = rownames(mat)[vals != 0],
           "up"     = rownames(mat)[vals == 1],
           "down"   = rownames(mat)[vals == -1])
  })

  total_de <- length(unique(unlist(protein_lists)))

  # Display labels
  display_names <- if (!is.null(contrast_labels)) contrast_labels else
    gsub("_", " ", contrasts)
  names(protein_lists) <- display_names

  colors    <- contrast_colors[seq_along(contrasts)]
  dir_label <- switch(direction,
                      "either" = "up- or down-regulated",
                      "up"     = "up-regulated only",
                      "down"   = "down-regulated only")

  # Extra canvas height for multi-line labels.
  # Each extra line beyond the first adds 0.4 inches to ggsave height.
  # The coord is left untouched so circle size is never affected.
  .max_lines     <- max(sapply(display_names,
                               function(x) length(strsplit(x, "\n")[[1]])))
  .subtitle_b_pt <- if (.max_lines > 1)
    5.5 + (.max_lines - 1) * label_line_height
  else
    5.5

  p <- ggvenn::ggvenn(
    protein_lists,
    show_stats     = "c",
    fill_color     = colors,
    fill_alpha     = 0.5,
    stroke_color   = colors,
    stroke_size    = 1.5,
    text_size      = text_size,
    set_name_size  = set_name_size,
    set_name_color = colors
  ) +
    { if (.max_lines > 1) ggplot2::coord_cartesian(clip = "off") else NULL } +
    ggplot2::labs(
      title    = "DE protein overlap across contrasts",
      subtitle = paste0("Proteins ", dir_label,
                        " at FDR < 0.05 and |log2FC| > log2(1.1)")
    ) +
    ggplot2::annotate(
      "text",
      x        = Inf,
      y        = -Inf,
      label    = paste0("Total unique DE proteins: ", total_de),
      hjust    = 1,
      vjust    = -1,
      size     = caption_size,
      color    = "black",
      fontface = "italic"
    ) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(hjust = 0.5, size = title_size),
      plot.subtitle = ggplot2::element_text(
        hjust  = 0.5,
        color  = "grey50",
        size   = subtitle_size,
        margin = ggplot2::margin(b = .subtitle_b_pt, unit = "pt")
      )
    )

  ggplot2::ggsave(filename = filename, plot = p,
                  width = 8, height = 7, dpi = 300,
                  bg = "white")

  message("Venn diagram saved to: ",
          normalizePath(filename, mustWork = FALSE))

  invisible(p)
}
