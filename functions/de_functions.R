# ==========================================================
# Differential analysis (limma) + subset helpers
# ==========================================================
library(limma)

# ---- design: ~ 0 + Group + covariates ---------------------
# meta must contain Group_de (syntactically valid factor) and the covariates
build_design <- function(meta, covariates = character(0)) {
  df <- data.frame(Group = meta$Group_de, row.names = rownames(meta))
  for (cv in covariates) {
    v <- meta[[cv]]
    df[[make.names(cv)]] <- if (is.numeric(v)) v else factor(v)
  }
  covs <- setdiff(colnames(df), "Group")
  f <- as.formula(paste("~ 0 + Group",
                        if (length(covs)) paste("+", paste(covs, collapse = " + ")) else ""))
  design <- model.matrix(f, data = df)
  rownames(design) <- rownames(df)
  colnames(design) <- sub("^Group", "", colnames(design))
  attr(design, "formula_text") <- paste(deparse(f), collapse = "")
  design
}

# ---- contrast safety check (expressions are eval()'d by makeContrasts) --
valid_contrast_expr <- function(expr, allowed_names) {
  if (!nzchar(trimws(expr))) return(FALSE)
  if (!grepl("^[A-Za-z0-9_.+*/() -]+$", expr)) return(FALSE)
  if (grepl("[A-Za-z0-9._)]\\s*\\(", expr)) return(FALSE)          # no function calls
  tokens <- regmatches(expr, gregexpr("[A-Za-z.][A-Za-z0-9._]*", expr))[[1]]
  all(tokens %in% allowed_names)
}

# plain matrix from a decideTests result (keeps column names exactly, NA -> 0)
dt_matrix <- function(dt) {
  m <- matrix(as.vector(dt), nrow(dt), ncol(dt), dimnames = list(rownames(dt), colnames(dt)))
  m[is.na(m)] <- 0L; m
}
# contrast names never contain underscores (they are shown with spaces everywhere)
clean_contrast_name <- function(x) trimws(gsub("_", " ", x))

# contrast_df: name, expr, type ("t" = ordinary contrast; "F" = joint test of several EXISTING contrasts,
# whose names are listed in expr separated by " | ")
make_contrast_matrix <- function(contrast_df, design) {
  if (!"type" %in% names(contrast_df)) contrast_df$type <- "t"
  t_rows <- contrast_df[contrast_df$type != "F", , drop = FALSE]
  f_rows <- contrast_df[contrast_df$type == "F", , drop = FALSE]
  if (!nrow(t_rows)) stop("Define at least one ordinary contrast (an F-test combines existing contrasts).")
  cm <- makeContrasts(contrasts = t_rows$expr, levels = design)
  colnames(cm) <- t_rows$name
  fsets <- setNames(lapply(f_rows$expr, function(e) trimws(strsplit(e, "|", fixed = TRUE)[[1]])), f_rows$name)
  for (nm in names(fsets)) {
    miss <- setdiff(fsets[[nm]], t_rows$name)
    if (length(miss)) stop("F-test '", nm, "' refers to contrasts that no longer exist: ", paste(miss, collapse = ", "))
    if (length(fsets[[nm]]) < 2) stop("F-test '", nm, "' needs at least two contrasts.")
  }
  attr(cm, "fsets") <- fsets
  attr(cm, "order_names") <- contrast_df$name
  cm
}

# counts per contrast (Down / NotSig / Up); F-test (global effect) columns only have NotSig / Up(=significant)
summarize_dt <- function(dt) {
  m <- dt_matrix(dt)
  out <- sapply(colnames(m), function(cn) c(Down = sum(m[, cn] == -1), NotSig = sum(m[, cn] == 0), Up = sum(m[, cn] == 1)))
  as.data.frame(out, check.names = FALSE)
}

# ---- full limma pipeline (mirrors the original script) ----
run_limma_pipeline <- function(y, design, contrast_matrix,
                               block = NULL, use_aw = FALSE,
                               test = c("treat", "eBayes"),
                               fc = 1.2, trend = TRUE, robust = TRUE, lfc = 0,
                               dt_method = "separate", adjust.method = "BH",
                               p.value = 0.05, progress = function(m) NULL) {
  test <- match.arg(test)

  aw <- NULL
  if (use_aw) {
    progress("arrayWeights")
    aw <- arrayWeights(y, design = design)
  }

  corr <- NULL
  if (!is.null(block)) {
    progress("duplicateCorrelation")
    corfit <- duplicateCorrelation(y, design = design, block = block, weights = aw)
    corr <- corfit$consensus
    progress("lmFit")
    fit <- lmFit(y, design, block = block, correlation = corr, weights = aw)
  } else {
    progress("lmFit")
    fit <- lmFit(y, design, weights = aw)
  }

  cfit <- contrasts.fit(fit, contrast_matrix)

  progress(test)
  if (test == "treat") {
    final_fit <- treat(cfit, fc = fc, trend = trend, robust = robust)
    tops <- lapply(colnames(cfit), function(cn)
      topTreat(final_fit, coef = cn, number = Inf, adjust.method = adjust.method))
    dt_lfc <- 0                      # threshold is already inside the treat p-values
  } else {
    final_fit <- eBayes(cfit, trend = trend, robust = robust)
    tops <- lapply(colnames(cfit), function(cn)
      topTable(final_fit, coef = cn, number = Inf, adjust.method = adjust.method))
    dt_lfc <- lfc
  }
  names(tops) <- colnames(cfit)
  tops <- lapply(tops, function(tt) {
    tt$Protein <- rownames(tt)
    tt[, c("Protein", setdiff(colnames(tt), "Protein")), drop = FALSE]
  })

  dt <- decideTests(final_fit, method = dt_method, adjust.method = adjust.method,
                    p.value = p.value, lfc = dt_lfc)

  # F-tests over several contrasts (e.g. Int_KO1 + Int_KO2, or KO1_vs_WT + KO2_vs_WT):
  #   topTable(fit, coef = c(...)) gives the joint F-test; decideTests(fit[, coef], method = "global")
  #   gives the direction of every component. F-tests need an eBayes fit (treat has no F statistic).
  fsets <- attr(contrast_matrix, "fsets"); order_names <- attr(contrast_matrix, "order_names")
  if (is.null(order_names)) order_names <- colnames(cfit)
  if (length(fsets)) {
    efit <- if (test == "eBayes") final_fit else eBayes(cfit, trend = trend, robust = robust)
    m <- matrix(0L, nrow(dt), length(order_names), dimnames = list(rownames(dt), order_names))
    for (nm in order_names) {
      if (nm %in% names(fsets)) {
        cols <- fsets[[nm]]
        tt0 <- topTable(efit, coef = cols, number = Inf, adjust.method = adjust.method, sort.by = "none")
        calls <- decideTests(efit[, cols], method = "global", adjust.method = adjust.method,
                             p.value = p.value, lfc = lfc)
        calls <- as.matrix(calls); calls[is.na(calls)] <- 0
        cl <- as.data.frame(calls[rownames(tt0), , drop = FALSE]); colnames(cl) <- paste0(cols, "_call")
        lf <- tt0[, seq_along(cols), drop = FALSE]; colnames(lf) <- cols      # one logFC column per contrast
        tt <- cbind(data.frame(Protein = rownames(tt0), F = tt0$F, P.Value = tt0$P.Value, adj.P.Val = tt0$adj.P.Val,
                               check.names = FALSE), lf, data.frame(AveExpr = tt0$AveExpr), cl)
        rownames(tt) <- tt$Protein
        tops[[nm]] <- tt
        # "significant" = joint F-test AND at least one component called (so the fold-change threshold applies)
        any_call <- rowSums(calls[rownames(tt), , drop = FALSE] != 0) > 0
        m[rownames(tt), nm] <- as.integer(!is.na(tt$adj.P.Val) & tt$adj.P.Val < p.value & any_call)
      } else { x <- dt[, nm]; x[is.na(x)] <- 0L; m[, nm] <- as.integer(x) }
    }
    dt <- m; tops <- tops[order_names]
  }

  list(fit = final_fit, cfit = cfit, dt = dt, tops = tops,
       p_value = p.value, fsets = fsets, single_names = setdiff(order_names, names(fsets)), ftest_names = names(fsets), correlation = corr, weights = aw,
       design = design, contrast_matrix = contrast_matrix, test = test)
}

# ---- predefined subsets (your function, with small guards) -
get_all_subsets_from_decideTests <- function(dt_results,
                                             method = c("venn", "contrast"),
                                             min_proteins = 5) {
  method <- match.arg(method)
  df <- as.data.frame(dt_matrix(dt_results), check.names = FALSE)
  df[is.na(df)] <- 0
  ids <- rownames(df); cols <- colnames(df); n_cols <- length(cols)
  subsets <- list()

  if (method == "contrast") {
    for (contrast in cols) {
      subsets[[paste0(contrast, "_up")]]   <- ids[df[[contrast]] == 1]
      subsets[[paste0(contrast, "_down")]] <- ids[df[[contrast]] == -1]
    }
  }

  if (method == "venn") {
    for (target in cols) {
      others <- setdiff(cols, target)
      others_all_zero <- if (length(others) > 0) rowSums(df[others] != 0) == 0 else TRUE
      subsets[[paste0("unique_", target, "_up")]]   <- ids[df[[target]] == 1  & others_all_zero]
      subsets[[paste0("unique_", target, "_down")]] <- ids[df[[target]] == -1 & others_all_zero]
    }
    if (n_cols >= 2) {
      pairs <- utils::combn(cols, 2)
      for (i in seq_len(ncol(pairs))) {
        c1 <- pairs[1, i]; c2 <- pairs[2, i]
        subsets[[paste0("common_", c1, "_&_", c2, "_both_up")]]   <- ids[df[[c1]] == 1  & df[[c2]] == 1]
        subsets[[paste0("common_", c1, "_&_", c2, "_both_down")]] <- ids[df[[c1]] == -1 & df[[c2]] == -1]
        subsets[[paste0("discordant_", c1, "_UP_", c2, "_DOWN")]] <- ids[df[[c1]] == 1  & df[[c2]] == -1]
        subsets[[paste0("discordant_", c2, "_UP_", c1, "_DOWN")]] <- ids[df[[c2]] == 1  & df[[c1]] == -1]
      }
    }
    if (n_cols > 2) {
      subsets[["common_all_contrasts_up"]]   <- ids[rowSums(df == 1)  == n_cols]
      subsets[["common_all_contrasts_down"]] <- ids[rowSums(df == -1) == n_cols]
      subsets[["common_all_contrasts_any"]]  <- ids[rowSums(df != 0)  == n_cols]
    }
  }
  if (!length(subsets)) return(list())
  subsets[vapply(subsets, length, integer(1)) >= min_proteins]
}

# ---- custom logic template --------------------------------
# rules: named character vector (names = contrasts) with values
#        "Ignore", "Up", "Down", "Any"
# mode : "AND" (all active rules) or "OR" (any active rule)
# only : if TRUE, contrasts set to "Ignore" must be non-significant
subset_from_rules <- function(dtm, rules, mode = c("AND", "OR"), only = FALSE) {
  mode <- match.arg(mode)
  cols <- names(rules)
  active <- rules != "Ignore"
  if (!any(active)) return(character(0))
  m <- vapply(cols[active], function(cn) {
    switch(rules[[cn]],
           Up = dtm[, cn] == 1, Down = dtm[, cn] == -1, Any = dtm[, cn] != 0)
  }, logical(nrow(dtm)))
  m <- matrix(m, nrow = nrow(dtm))
  hit <- if (mode == "AND") rowSums(m) == ncol(m) else rowSums(m) > 0
  if (only && any(!active))
    hit <- hit & rowSums(dtm[, cols[!active], drop = FALSE] != 0) == 0
  rownames(dtm)[hit]
}

# ==========================================================
# Volcano wrapper around your my_volcano_plot_proteomics()
# (it saves a PNG into the working dir; we run it in a temp dir and keep the ggplot)
# ==========================================================
volcano_from_top <- function(top, title, filter_choice = "adj.P.Val", y_axis_choice = "P.Value",
                             l2FC_thresh = 1, p_thresh = 0.05, show_names = "none",
                             which_genes = "all", grid_lines = FALSE, label_size = 2.5,
                             highlight = character(0), width = 6.5, height = 5.8) {
  top$Gene <- top$Protein
  p <- run_in_tmpdir(function()
    my_volcano_plot_proteomics(
      DE_res = top, filter_choice = filter_choice, y_axis_choice = y_axis_choice,
      l2FC_thresh = l2FC_thresh, p_thresh = p_thresh, plot_title = title,
      show_names = show_names, which_genes = which_genes, grid_lines = grid_lines,
      save_prefix = "tmp", width = width, height = height, label_size = label_size))$value
  hl <- top[top$Protein %in% highlight, , drop = FALSE]
  if (nrow(hl) > 0) {
    hl$y_value <- -log10(hl[[y_axis_choice]])
    p <- p +
      ggplot2::geom_point(data = hl, ggplot2::aes(x = logFC, y = y_value), inherit.aes = FALSE,
                          shape = 21, size = 3.2, stroke = 1.1, fill = NA, colour = "black") +
      ggrepel::geom_label_repel(data = hl, ggplot2::aes(x = logFC, y = y_value, label = Gene),
                                inherit.aes = FALSE, fontface = "bold", size = label_size + 0.5,
                                box.padding = 0.5, max.overlaps = Inf, segment.colour = "black",
                                fill = "white", label.size = 0.25)
  }
  p
}

# ==========================================================
# Diverging bar plot (your ggplot2 code; breaks/limits generalised,
# optional y-axis contrast labels because the original hides them)
# ==========================================================
deg_summary_from_dt <- function(dt, f_names = character()) {
  m <- dt_matrix(dt)
  dplyr::bind_rows(lapply(colnames(m), function(cn) {
    is_f <- cn %in% f_names
    data.frame(Raw_Contrast = cn, Is_F = is_f,
               Upregulated   = if (is_f) 0L else sum(m[, cn] == 1),
               Downregulated = if (is_f) 0L else sum(m[, cn] == -1),
               FTest         = if (is_f) sum(m[, cn] == 1) else 0L)       # F-test: no direction, just "changed"
  }))
}

diverging_bar_plot <- function(deg_summary, show_labels = TRUE, num_size_pt = 12) {
  library(dplyr); library(tidyr); library(ggplot2); library(scales)
  if (!"FTest" %in% names(deg_summary)) { deg_summary$FTest <- 0L; deg_summary$Is_F <- FALSE }
  deg_df <- deg_summary %>%
    dplyr::mutate(Clean_Contrast = factor(Raw_Contrast, levels = rev(unique(Raw_Contrast)))) %>%
    tidyr::pivot_longer(cols = c(Upregulated, Downregulated, FTest),
                        names_to = "Direction", values_to = "Count") %>%
    dplyr::filter(ifelse(Is_F, Direction == "FTest", Direction != "FTest")) %>%     # each row only its own kind of bar
    dplyr::mutate(Direction = ifelse(Direction == "FTest", "F-test (any change)", Direction),
                  Plot_Count = ifelse(Direction == "Downregulated", -abs(Count), abs(Count)))
  axis_max <- max(abs(deg_df$Plot_Count), 1)
  x_breaks <- pretty(c(-axis_max, axis_max), n = 6)
  axis_max <- max(abs(x_breaks))
  pad <- axis_max * 0.18          # room for the counts printed outside the bars

  ggplot(deg_df, aes(x = Plot_Count, y = Clean_Contrast, fill = Direction)) +
    geom_bar(stat = "identity", color = "grey20", linewidth = 0.3, width = 0.6) +
    scale_x_continuous(breaks = x_breaks, limits = c(-axis_max - pad, axis_max + pad),
                       labels = function(x) scales::comma(abs(x)), expand = c(0, 0)) +
    scale_fill_manual(values = c(Upregulated = "#d73027", Downregulated = "#4575b4", "F-test (any change)" = "#7f7f7f")) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.8) +
    geom_text(aes(label = scales::comma(abs(Count)), hjust = ifelse(Plot_Count >= 0, -0.25, 1.25)),
              size = num_size_pt / ggplot2::.pt, colour = "black") +
    labs(title = NULL, x = NULL, y = NULL) +
    theme_classic(base_size = 16) +
    theme(legend.position = "top", legend.title = element_blank(),
          axis.line.y = element_line(color = "black", linewidth = 0.8),
          axis.text.y = if (show_labels) element_text(size = 11, color = "black") else element_blank(),
          axis.ticks.y = element_line(color = "black", linewidth = 0.8),
          axis.ticks.length.y = unit(0.2, "cm"),
          axis.text.x = element_text(size = 12, color = "black"))
}

# ==========================================================
# Venn of significant proteins: same engine as venn_of_detected()
# (contrast lists are turned into a presence/absence matrix, so every option of
#  venn_of_detected - number size, circle-to-number ratio, ... - applies; no title/caption)
# ==========================================================
venn_sig_plot <- function(dt, contrasts, direction = "either", colors = NULL,
                          font_scale = 1, count_size = NULL, group_label_size = NULL,
                          expand_limits = 1.5, fill_alpha = 0.45, labels = NULL) {
  m <- dt_matrix(dt)[, contrasts, drop = FALSE]
  keep <- function(v) switch(direction, either = v != 0, up = v == 1, down = v == -1)
  mat <- vapply(contrasts, function(ct) ifelse(keep(m[, ct]), 1, NA_real_), numeric(nrow(m)))
  rownames(mat) <- rownames(m); colnames(mat) <- contrasts
  si <- data.frame(Sample = contrasts, Group = contrasts, stringsAsFactors = FALSE)
  if (is.null(colors)) colors <- setNames(c("#1f78b4", "#33a02c", "#e31a1c", "#ff7f00", "#6a3d9a", "#b15928", "#999999")[seq_along(contrasts)], contrasts)
  p <- run_in_tmpdir(function()
    venn_of_detected(mat, si, threshold = 0.5, group_colors = colors, group_labels = labels,
                     title = NULL, subtitle = "",
                     font_scale = font_scale, count_size = count_size,
                     group_label_size = group_label_size, expand_limits = expand_limits,
                     fill_alpha = fill_alpha, filename = "v.png"))$value
  p + ggplot2::theme(plot.title = ggplot2::element_blank(), plot.subtitle = ggplot2::element_blank())
}
