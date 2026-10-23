#' Generate High-Resolution Log2 Intensity Boxplot
#'
#' @description
#' Creates boxplot of Log2-transformed protein intensities.
#'
#' @param data Numeric matrix. The protein intensity data.
#'   The function automatically performs Log2 transformation internally.
#' @param sample_info Data frame. Must contain `Sample` and `Group` columns.
#' @param group_colors Named vector. Colors assigned to each `Group`.
#' @param plot_prefix Character. Suffix for the output filename. Default is "1st_QC".
#'
#' @return Generates a PNG file named `Log2_LFQ_box_plot_[prefix].png`. Returns NULL.
#'
#' @importFrom graphics boxplot layout par legend plot.new
#' @importFrom grDevices png dev.off
#'
#' @export

# ==============================================================================
# Function: Generate Boxplot
# ==============================================================================
generate_boxplot <- function(data,
                             sample_info,
                             group_colors,
                             plot_prefix = "") {

  {
    filename <- file.path(paste0("Boxplot_", plot_prefix, ".png"))

    png(filename, width = 6000, height = 3000, res = 300)
    box_colors <- group_colors[sample_info$Group]
    names(box_colors) <- sample_info$Sample
    layout(matrix(c(1,2), ncol = 2), widths = c(5,0.6))
    mb <- max(5, 0.62 * max(nchar(colnames(data))) + 3)
    par(mar = c(mb,6,3,0))
    par(mgp = c(3.5, 0.7, 0))
    boxplot(data,
            main = "Log2 LFQ intensity distributions per sample",
            ylab = "Log2(LFQ intensity)",
            xlab = "",
            cex.main = 1.5,
            cex.lab = 1.2,
            las = 2,
            outline = FALSE,
            col = box_colors)
    title(xlab = "Samples", line = mb - 1.3, cex.lab = 1.2)

    par(mar = c(0,0,0,0))
    plot.new()
    legend("left",inset = c(0, 0),
           legend = c("Group", names(group_colors)),
           fill = c(NA, group_colors),
           border = c(NA, "gray30"),
           text.font = c(2, rep(1, length(group_colors))),
           cex = c(1.3, rep(1.2, length(group_colors))),
           bty = "n",
           x.intersp = 0.3)

    dev.off()
  }
}
