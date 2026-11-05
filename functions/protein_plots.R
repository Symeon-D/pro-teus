# ==========================================================
# Plots for selected proteins: grid of subplots, one per protein
# ==========================================================
PROTEIN_PLOT_TYPES <- c("Boxplot" = "boxplot", "Violin" = "violin", "Bar (mean +/- SE)" = "bar",
                        "Dot / strip" = "dot", "Line (group means +/- SE)" = "line")

protein_long <- function(m, proteins, meta) {
  proteins <- intersect(proteins, rownames(m))
  grp <- meta$Group[match(colnames(m), meta$Sample)]
  df <- do.call(rbind, lapply(proteins, function(p)
    data.frame(Protein = p, Sample = colnames(m), Group = grp, Value = as.numeric(m[p, ]), stringsAsFactors = FALSE)))
  df$Protein <- factor(df$Protein, levels = proteins)               # keeps the order the user chose
  df$Group <- factor(df$Group, levels = unique(grp))
  df
}

# brackets: data.frame(Protein, g1, g2, label) -> bracket between two groups with stars (from limma);
# notes: named character vector (Protein -> text) for contrasts that are not a plain two-group comparison
protein_plot_grid <- function(m, proteins, meta, type = "boxplot", ncol = 3, free_y = TRUE, colors = NULL,
                              ylab = "log2 intensity", brackets = NULL, notes = NULL) {
  df <- protein_long(m, proteins, meta)
  if (!nrow(df)) return(NULL)
  p <- ggplot2::ggplot(df, ggplot2::aes(Group, Value, colour = Group, fill = Group))
  jit <- ggplot2::geom_jitter(width = 0.15, height = 0, size = 1.6, alpha = 0.9, show.legend = FALSE)
  p <- p + switch(type,
    boxplot = list(ggplot2::geom_boxplot(alpha = 0.35, outlier.shape = NA, show.legend = FALSE), jit),
    violin  = list(ggplot2::geom_violin(alpha = 0.35, show.legend = FALSE, scale = "width"), jit),
    bar     = list(ggplot2::stat_summary(fun = mean, geom = "col", alpha = 0.45, width = 0.7, show.legend = FALSE),
                   ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.25, colour = "grey20"), jit),
    dot     = list(jit, ggplot2::stat_summary(fun = mean, geom = "crossbar", width = 0.55, linewidth = 0.3, colour = "grey20", fill = NA)),
    line    = list(ggplot2::stat_summary(ggplot2::aes(group = 1), fun = mean, geom = "line", colour = "grey30", linewidth = 0.7),
                   ggplot2::stat_summary(fun.data = ggplot2::mean_se, geom = "errorbar", width = 0.2, colour = "grey30"),
                   ggplot2::stat_summary(fun = mean, geom = "point", size = 3, show.legend = FALSE), jit))
  if (!is.null(colors) && all(levels(df$Group) %in% names(colors)))
    p <- p + ggplot2::scale_colour_manual(values = colors) + ggplot2::scale_fill_manual(values = colors)
  if (!is.null(brackets) && nrow(brackets)) {
    st <- do.call(rbind, lapply(levels(df$Protein), function(pr) {
      v <- df$Value[df$Protein == pr]
      data.frame(Protein = pr, ymax = max(v, na.rm = TRUE), rng = max(diff(range(v, na.rm = TRUE)), 1e-6), stringsAsFactors = FALSE)
    }))
    b <- merge(brackets, st, by = "Protein")
    b$x1 <- match(b$g1, levels(df$Group)); b$x2 <- match(b$g2, levels(df$Group))
    b <- b[!is.na(b$x1) & !is.na(b$x2), , drop = FALSE]
    if (nrow(b)) {
      b$span <- abs(b$x2 - b$x1); b <- b[order(b$Protein, b$span), , drop = FALSE]
      b$k <- stats::ave(seq_len(nrow(b)), b$Protein, FUN = seq_along)
      b$y <- b$ymax + b$k * 0.14 * b$rng; b$tick <- 0.03 * b$rng
      b$Protein <- factor(b$Protein, levels = levels(df$Protein))
      p <- p +
        ggplot2::geom_segment(data = b, ggplot2::aes(x = x1, xend = x2, y = y, yend = y), inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
        ggplot2::geom_segment(data = b, ggplot2::aes(x = x1, xend = x1, y = y, yend = y - tick), inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
        ggplot2::geom_segment(data = b, ggplot2::aes(x = x2, xend = x2, y = y, yend = y - tick), inherit.aes = FALSE, colour = "grey20", linewidth = 0.4) +
        ggplot2::geom_text(data = b, ggplot2::aes(x = (x1 + x2) / 2, y = y, label = label), inherit.aes = FALSE, vjust = -0.1, size = 4) +
        ggplot2::geom_blank(data = b, ggplot2::aes(y = y + 0.12 * rng), inherit.aes = FALSE)
    }
  }
  lab <- if (length(notes)) ggplot2::as_labeller(function(x) ifelse(x %in% names(notes), paste0(x, "\n", notes[x]), x)) else "label_value"
  p + ggplot2::facet_wrap(~Protein, ncol = ncol, scales = if (free_y) "free_y" else "fixed", labeller = lab) +
    ggplot2::labs(x = NULL, y = ylab) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1), legend.position = "none",
                   strip.text = ggplot2::element_text(face = "bold"))
}

protein_plot_height <- function(n, ncol, n_brackets = 0, row_px = 270) max(300, ceiling(max(n, 1) / max(ncol, 1)) * (row_px + 24 * n_brackets) + 40)
