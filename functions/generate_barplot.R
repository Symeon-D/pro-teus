#' Generate High-Resolution Protein Count Barplot
#'
#' @description
#' Creates a barplot visualizing the number of quantified proteins
#' per sample.
#'
#' @param data A numeric matrix or data frame containing protein intensity data.
#'   Rows should represent proteins and columns should represent samples.
#' @param sample_info A data frame containing sample metadata.
#'   Must include a column named `Group` that matches the keys provided in `group_colors`.
#' @param group_colors A named vector of colors where names correspond to the `Group` levels
#'   found in `sample_info`. Example: `c("Vehicle" = "#a6cee3", "Treatment" = "#e31a1c")`.
#' @param ceiling_num Numeric. The rounding factor for the Y-axis upper limit.
#'   The maximum protein count is rounded up to the nearest multiple of this number.
#'   Default is `1000`.
#' @param plot_prefix Character string. A suffix to append to the output filename
#'   and plot title (e.g., "1st_QC" or "Post_Normalization"). Default is `"1st_QC"`.
#' @param width Width of the plot.
#' @param height height of the plot.
#' @return Generates a PNG file named `Number_of_Proteins_Barplot_[plot_prefix].png`
#'   in the current working directory. The function returns `NULL` invisibly.
#'
#' @details
#' The plot uses `layout()` to split the graphical device into two unequal parts (ratio 5:0.7).
#' The Y-axis limit is calculated dynamically: \code{ceiling(max_count / ceiling_num) * ceiling_num}.
#' This ensures the axis extends to a clean, round number just above the highest bar.
#'
#' @importFrom graphics barplot par layout axis title legend plot.new
#' @importFrom grDevices png dev.off
#' @export

# ==============================================================================
# Function: Generate Barplot
# ==============================================================================
generate_barplot <- function(data,
                             sample_info,
                             group_colors,
                             ymax = NULL,
                             ymin = 0,
                             ceiling_num = 100,
                             plot_prefix = "",
                             width = 6000,
                             height =3500 ){

  filename <- file.path(paste0("Barplot_", plot_prefix, ".png"))

  identified_counts <- colSums(!is.na(data))
  barpos <- barplot(identified_counts, plot = FALSE)
  ## Auto ymax ONLY if user does not provide it
  if (is.null(ymax)) {
    ymax <- ceiling(max(identified_counts, na.rm = TRUE) / 100) * 100   # always the next hundred
  }

  ## Shift counts for cropped plotting
  cropped_counts <- pmax(identified_counts, ymin) - ymin
  png(filename, width = width, height = height, res = 300)

  layout(matrix(c(1,2), ncol = 2), widths = c(5,0.7))
  mb <- max(5, 0.62 * max(nchar(colnames(data))) + 3)   # bottom margin grows with sample-name length
  par(mar = c(mb,5,5,0))
  par(mgp = c(3.5, 0.5, 0))

  barplot(cropped_counts,
          las = 2,
          main = "",
          xlab = "",
          ylab = "Number of quantified proteins",
          cex.main = 1.6,
          cex.lab = 1.4,
          col = group_colors[sample_info$Group],
          xaxs = "i",
          xlim = c(-0.5, max(barpos) + 0.5),
          ylim = c(0, ymax - ymin),
          yaxs = "i",
          yaxt = "n"
  )

  ## Correct axis: positions are shifted, labels are real values
  top <- ymax - ymin
  axis_ticks <- pretty(c(0, top), n = 5)
  axis_ticks <- axis_ticks[axis_ticks <= top]
  if (tail(axis_ticks, 1) < top) {                 # e.g. 4883 proteins: ticks ... 4000 and a last tick at 4900
    if (top - tail(axis_ticks, 1) < 0.06 * top) axis_ticks <- head(axis_ticks, -1)   # avoid two crowded labels
    axis_ticks <- c(axis_ticks, top)
  }
  title(xlab = "Samples", line = mb - 1.3, cex.lab = 1.4)
  axis(side = 2,
       at = axis_ticks,
       labels = axis_ticks + ymin,
       las = 1)

  title(main = "Proteins quantified per sample",
        line = 2.5,
        cex.main = 1.8)

  par(mar = c(0,0,0,0))
  plot.new()

  legend("left",
         inset = c(0, 0),
         legend = c("Group", names(group_colors)),
         fill = c(NA, group_colors),
         border = c(NA, "gray30"),
         text.font = c(2, rep(1, length(group_colors))),
         cex = c(1.3, rep(1.2, length(group_colors))),
         bty = "n",
         x.intersp = 0.3)

  dev.off()
}
