# ==========================================================
# Batch & SVA helpers (diagnostics, variance partitioning, surrogate variables, correction, PCA)
# All functions return data / ggplot objects and never write files.
# ==========================================================

pval_to_stars <- function(p) ifelse(is.na(p), "", ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", ""))))

# bias-corrected Cramer's V (Bergsma 2013): the plain V is biased upward with small tables / few samples
cramer_v <- function(tbl) {
  n <- sum(as.numeric(tbl)); r <- nrow(tbl); k <- ncol(tbl)
  if (n < 2 || r < 2 || k < 2) return(0)
  chi2 <- unname(suppressWarnings(stats::chisq.test(tbl, correct = FALSE))$statistic)
  phi2c <- max(0, chi2 / n - (k - 1) * (r - 1) / (n - 1))
  rc <- r - (r - 1)^2 / (n - 1); kc <- k - (k - 1)^2 / (n - 1)
  d <- min(kc - 1, rc - 1); if (d <= 0) return(0)
  sqrt(phi2c / d)
}

# binary / nominal / continuous / exclude, guessed from the column
guess_var_types <- function(meta, skip = character()) {
  n <- nrow(meta)
  vapply(names(meta), function(v) {
    x <- meta[[v]]; u <- length(unique(x[!is.na(x)]))
    if (v %in% skip || u < 2) return("exclude")
    if (is.numeric(x)) return(if (u == 2) "binary" else if (u <= 5) "nominal" else "continuous")
    if (u == 2) "binary" else if (u >= n / 2) "exclude" else "nominal"
  }, character(1))
}

as_num <- function(x) if (is.numeric(x)) x else as.numeric(factor(x))

# association between two variables -> list(stat, p, method)
pair_assoc <- function(a, b, ta, tb) {
  ok <- !is.na(a) & !is.na(b); a <- a[ok]; b <- b[ok]
  smooth <- c("continuous", "binary")
  if (ta %in% smooth && tb %in% smooth) {
    r <- suppressWarnings(stats::cor.test(as_num(a), as_num(b), method = "spearman", exact = FALSE))
    return(list(stat = unname(r$estimate), p = r$p.value, method = "spearman"))
  }
  if (ta == "continuous" || tb == "continuous") {                   # continuous vs nominal: ANOVA effect size
    y <- if (ta == "continuous") a else b; g <- factor(if (ta == "continuous") b else a)
    ss <- summary(stats::lm(as.numeric(y) ~ g))
    if (is.null(ss$fstatistic)) return(list(stat = NA, p = NA, method = "R2"))
    return(list(stat = sqrt(max(0, ss$adj.r.squared)),
                p = stats::pf(ss$fstatistic[1], ss$fstatistic[2], ss$fstatistic[3], lower.tail = FALSE), method = "R2"))
  }
  tbl <- table(a, b)
  if (min(dim(tbl)) < 2) return(list(stat = 0, p = NA, method = "cramers_v"))
  p <- tryCatch(suppressWarnings(stats::fisher.test(tbl, simulate.p.value = TRUE, B = 2000))$p.value, error = function(e) NA)
  list(stat = cramer_v(tbl), p = p, method = "cramers_v")
}

# ---------- 1. variable - variable relationships ----------
var_var_matrix <- function(meta, types) {
  vars <- names(types)[types != "exclude"]
  k <- length(vars)
  if (k < 2) return(NULL)
  cm <- pm <- matrix(NA_real_, k, k, dimnames = list(vars, vars)); diag(cm) <- 1; diag(pm) <- 0
  for (i in 1:(k - 1)) for (j in (i + 1):k) {
    r <- pair_assoc(meta[[vars[i]]], meta[[vars[j]]], types[[vars[i]]], types[[vars[j]]])
    cm[i, j] <- cm[j, i] <- r$stat; pm[i, j] <- pm[j, i] <- r$p
  }
  d <- reshape2::melt(cm); d$p <- reshape2::melt(pm)$value
  d$stars <- pval_to_stars(d$p)
  d$label <- ifelse(d$stars != "", paste0(round(d$value, 2), "\n", d$stars), as.character(round(d$value, 2)))
  d$Var1 <- factor(d$Var1, levels = vars); d$Var2 <- factor(d$Var2, levels = vars)
  plt <- ggplot2::ggplot(d, ggplot2::aes(Var2, Var1, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5, lineheight = 0.85) +
    ggplot2::scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3", midpoint = 0,
                                  limits = c(-1, 1), name = "r / V") +
    ggplot2::labs(title = "Variable-variable relationships", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, vjust = 1),
                   plot.title = ggplot2::element_text(face = "bold"))
  list(cor = cm, p = pm, plot = plt)
}

# ---------- 2. PCs vs covariates (complete rows only) ----------
pc_covariates <- function(m, meta, types, n_pc = 10) {
  vars <- names(types)[types != "exclude"]
  m <- m[stats::complete.cases(m), , drop = FALSE]; m <- m[apply(m, 1, stats::var) > 0, , drop = FALSE]
  n_pc <- min(n_pc, ncol(m) - 1, nrow(m))
  if (n_pc < 2 || !length(vars)) return(NULL)
  pca <- stats::prcomp(t(m), center = TRUE, scale. = FALSE, rank. = n_pc)
  PCs <- as.data.frame(pca$x); pvar <- round(summary(pca)$importance[2, seq_len(n_pc)] * 100, 2)
  pcn <- paste0("PC", seq_len(n_pc))
  cm <- pm <- matrix(NA_real_, n_pc, length(vars), dimnames = list(pcn, vars))
  for (i in seq_len(n_pc)) for (v in vars) {
    x <- meta[[v]]
    if (types[[v]] == "continuous") {
      r <- suppressWarnings(stats::cor.test(PCs[[i]], as.numeric(x), method = "spearman", exact = FALSE))
      cm[i, v] <- abs(r$estimate); pm[i, v] <- r$p.value
    } else {
      ss <- summary(stats::lm(PCs[[i]] ~ factor(x)))
      if (!is.null(ss$fstatistic)) {
        cm[i, v] <- sqrt(max(0, ss$adj.r.squared))
        pm[i, v] <- stats::pf(ss$fstatistic[1], ss$fstatistic[2], ss$fstatistic[3], lower.tail = FALSE)
      }
    }
  }
  imp <- matrix(0, n_pc, length(vars), dimnames = list(pcn, vars))
  for (i in seq_len(n_pc)) for (j in seq_along(vars))
    if (!is.na(pm[i, j]) && pm[i, j] < 0.05) imp[i, j] <- cm[i, j]^2 * pvar[i]
  cd <- reshape2::melt(cm); cd$stars <- pval_to_stars(reshape2::melt(pm)$value)
  cd$label <- ifelse(cd$stars == "", round(cd$value, 2), paste0(round(cd$value, 2), "\n", cd$stars))
  cd$Var1 <- factor(cd$Var1, levels = pcn); cd$Var2 <- factor(cd$Var2, levels = vars)
  p_cor <- ggplot2::ggplot(cd, ggplot2::aes(Var2, Var1, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5, lineheight = 0.85) +
    ggplot2::scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3", midpoint = 0,
                                  limits = c(0, 1), name = "|r| / sqrt(R2)") +
    ggplot2::labs(title = "PC ~ covariate correlations", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))
  id <- reshape2::melt(imp); id$label <- round(id$value, 2)
  id$Var1 <- factor(id$Var1, levels = pcn); id$Var2 <- factor(id$Var2, levels = vars)
  p_imp <- ggplot2::ggplot(id, ggplot2::aes(Var2, Var1, fill = value)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5) +
    ggplot2::scale_fill_gradient(low = "white", high = "firebrick3", name = "Impact") +
    ggplot2::labs(title = "PC covariate impact", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))
  list(plot_cor = p_cor, plot_impact = p_imp, n_rows = nrow(m))
}

# ---------- 3. variance partitioning ----------
# cat_vars: random effects (1|var); cont_vars: fixed numeric terms; sv: optional data.frame of SVs (fixed terms)
run_varpart <- function(m, meta, cat_vars, cont_vars = character(), sv = NULL) {
  m <- m[stats::complete.cases(m), , drop = FALSE]; m <- m[apply(m, 1, stats::var) > 0, , drop = FALSE]
  md <- meta[match(colnames(m), meta$Sample), , drop = FALSE]
  if (!is.null(sv)) md <- cbind(md, sv[colnames(m), , drop = FALSE])
  orig <- c(cat_vars, cont_vars, if (!is.null(sv)) colnames(sv)); safe <- make.names(orig, unique = TRUE)
  dd <- data.frame(row.names = colnames(m))
  for (i in seq_along(orig)) dd[[safe[i]]] <- if (orig[i] %in% cat_vars) factor(md[[orig[i]]]) else as.numeric(md[[orig[i]]])
  terms <- c(paste0("(1|", safe[seq_along(cat_vars)], ")"), safe[-seq_along(cat_vars)])
  form <- stats::as.formula(paste("~", paste(terms, collapse = " + ")))
  vp <- tryCatch(variancePartition::fitExtractVarPartModel(m, form, dd, quiet = TRUE, BPPARAM = BiocParallel::SerialParam()),
                 error = function(e) variancePartition::fitExtractVarPartModel(m, form, dd, quiet = TRUE))
  colnames(vp) <- ifelse(colnames(vp) %in% safe, orig[match(colnames(vp), safe)], colnames(vp))
  med <- apply(vp, 2, stats::median) * 100; mn <- apply(vp, 2, mean) * 100
  summ <- data.frame(Variable = names(med), Median_Variance_Pct = round(med, 2), Mean_Variance_Pct = round(mn, 2), row.names = NULL)
  summ <- summ[order(-summ$Median_Variance_Pct), ]
  list(vp = vp, summary = summ, n_rows = nrow(m))
}
plot_varpart <- function(vp, title = "Variance explained by covariates") {
  variancePartition::plotVarPart(variancePartition::sortCols(vp)) + ggplot2::ggtitle(title) +
    ggplot2::theme(plot.title = ggplot2::element_text(face = "bold", hjust = 0.5))
}

# ---------- 4. surrogate variables ----------
# protect: variables of interest (only in the full model); adjust: known technical covariates (full AND null model,
# as recommended: the null model holds the adjustment variables only)
run_sva <- function(m, md, protect, adjust = character(0), method = "leek", n_sv = NULL, seed = 2026) {
  if (anyNA(m)) stop("SVA needs a matrix without missing values (run imputation first).")
  m <- m[apply(m, 1, stats::var) > 0, , drop = FALSE]
  f_full <- stats::as.formula(paste("~", paste0("`", c(protect, adjust), "`", collapse = " + ")))
  mod_full <- stats::model.matrix(f_full, data = md)
  mod_null <- if (length(adjust)) stats::model.matrix(stats::as.formula(paste("~", paste0("`", adjust, "`", collapse = " + "))), data = md) else stats::model.matrix(~ 1, data = md)
  if (qr(mod_full)$rank < ncol(mod_full)) stop("The chosen variables are confounded (the model is not full rank): remove one of them.")
  set.seed(seed)
  n_est <- sva::num.sv(as.matrix(m), mod_full, method = method)
  n_use <- if (is.null(n_sv) || is.na(n_sv)) n_est else as.integer(n_sv)
  if (n_use < 1) return(list(n_est = n_est, n_used = 0L, sv = NULL))
  obj <- sva::sva(dat = as.matrix(m), mod = mod_full, mod0 = mod_null, n.sv = n_use)
  sv <- as.data.frame(obj$sv); colnames(sv) <- paste0("SV", seq_len(ncol(sv))); rownames(sv) <- colnames(m)
  list(n_est = n_est, n_used = n_use, sv = sv)
}

# SV x covariate association table + heatmap
sv_associations <- function(sv, meta, types) {
  vars <- names(types)[types != "exclude"]
  rows <- list()
  for (s in colnames(sv)) for (v in vars) {
    t_v <- types[[v]]
    r <- pair_assoc(sv[[s]], meta[[v]], "continuous", t_v)
    rows[[length(rows) + 1]] <- data.frame(SV = s, Variable = v, stat = round(r$stat, 3), p = signif(r$p, 3),
                                           sig = ifelse(is.na(r$p), "", ifelse(r$p < 0.001, "***", ifelse(r$p < 0.01, "**", ifelse(r$p < 0.05, "*", "ns")))),
                                           stringsAsFactors = FALSE)
  }
  tab <- do.call(rbind, rows)
  d <- tab; d$label <- ifelse(d$sig %in% c("", "ns"), round(d$stat, 2), paste0(round(d$stat, 2), "\n", d$sig))
  d$SV <- factor(d$SV, levels = colnames(sv)); d$Variable <- factor(d$Variable, levels = vars)
  plt <- ggplot2::ggplot(d, ggplot2::aes(Variable, SV, fill = stat)) +
    ggplot2::geom_tile(colour = "white", linewidth = 0.5) +
    ggplot2::geom_text(ggplot2::aes(label = label), size = 3.5, lineheight = 0.85) +
    ggplot2::scale_fill_gradient2(low = "steelblue4", mid = "white", high = "firebrick3", midpoint = 0, limits = c(-1, 1), name = "r / sqrt(R2)") +
    ggplot2::labs(title = "Surrogate variables ~ covariates", x = NULL, y = NULL) +
    ggplot2::theme_minimal(base_size = 12) + ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 35, hjust = 1, vjust = 1))
  list(table = tab, plot = plt)
}

# ---------- 5. correction for visualisation / protein plots ----------
# removes exactly the covariates in `covs` (columns of `md`), protecting the group design
correct_for_covariates <- function(m, md, covs, group) {
  covs <- covs[covs %in% colnames(md)]
  if (!length(covs)) return(m)
  d <- md[, covs, drop = FALSE]
  for (v in covs) if (!is.numeric(d[[v]])) d[[v]] <- factor(d[[v]])
  covmat <- stats::model.matrix(stats::as.formula(paste("~", paste0("`", covs, "`", collapse = " + "))), data = d)[, -1, drop = FALSE]
  design <- stats::model.matrix(~ 0 + factor(group))
  limma::removeBatchEffect(as.matrix(m), covariates = covmat, design = design)
}

# ---------- 6. PCA through the app's PCA function (colour by / shape by) ----------
pca_palette <- function(meta, color_by, group_colors) {
  if (identical(color_by, "Group") && !is.null(group_colors)) return(group_colors)
  lv <- unique(as.character(meta[[color_by]])); auto_colors(lv)
}

pca_one <- function(m, meta, color_by = "Group", shape_by = NULL, group_colors = NULL, title = NULL) {
  m <- m[stats::complete.cases(m), , drop = FALSE]; m <- m[apply(m, 1, stats::var) > 0, , drop = FALSE]
  if (ncol(m) < 3 || nrow(m) < 2) return(NULL)
  md <- meta[match(colnames(m), meta$Sample), , drop = FALSE]
  sh <- if (is.null(shape_by) || identical(shape_by, "none") || identical(shape_by, color_by)) NULL else shape_by
  res <- my_plot_PCA_matrix(expression_matrix = m, sample_info = md, pca_method = "classic", top_n = "all", scale = TRUE,
                            show_names = nrow(md) <= 40, add_ellipse = FALSE, fill_ellipse = FALSE, grid_lines = FALSE,
                            width = 8, height = 6, output_path = tempfile(fileext = ".png"),
                            group_colors = pca_palette(md, color_by, group_colors), color_by = color_by, shape_by = sh,
                            mv_handle = "ppca")
  if (!is.null(title)) res$pca_plot + ggplot2::ggtitle(title) else res$pca_plot
}

# two PCA plots (before / after) with one shared legend
pca_before_after <- function(m_pre, m_post, meta, color_by = "Group", shape_by = NULL, group_colors = NULL) {
  p1 <- pca_one(m_pre, meta, color_by, shape_by, group_colors, "Before correction")
  p2 <- pca_one(m_post, meta, color_by, shape_by, group_colors, "After correction")
  if (is.null(p1) || is.null(p2)) return(NULL)
  leg <- cowplot::get_legend(p1 + ggplot2::theme(legend.position = "bottom"))
  cowplot::plot_grid(cowplot::plot_grid(p1 + ggplot2::theme(legend.position = "none"),
                                        p2 + ggplot2::theme(legend.position = "none"), nrow = 1),
                     leg, ncol = 1, rel_heights = c(1, 0.12))
}
