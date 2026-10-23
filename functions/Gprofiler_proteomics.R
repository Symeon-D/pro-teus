#' Run g:Profiler enrichment for PART clusters or a gene vector
#'
#' Runs functional enrichment analysis via `gprofiler2::gost()` using either:
#' (i) cluster-specific gene sets stored in a PART results object, or
#' (ii) a user-provided gene vector.
#'
#' This function is **compute-only** and returns a single results table (`GO_df`).
#' Plotting is intentionally excluded and should be done with `gprofiler_wrapper()`.
#'
#' @param PART_results Optional list produced by `my_compute_PART_proteomics()`.
#'   Must contain `part_data` (with a `gene_cluster` column and rownames as gene IDs)
#'   and `cluster_map` (with `cluster` and `cluster_col` for color mapping).
#'   If provided, enrichment is run per cluster.
#' @param gene_vector Optional character vector of gene identifiers (symbols/Ensembl,
#'   depending on `organism` and g:Profiler expectations). Used when `PART_results` is `NULL`.
#' @param background Optional custom background (universe). Either:
#'   a character vector of gene IDs, or a matrix/data.frame whose rownames are gene IDs.
#'   If `NULL`, g:Profiler uses the annotated universe (`domain_scope = "annotated"`).
#' @param ordered_query Logical; passed to `gprofiler2::gost(ordered_query=...)`.
#'   Set to `TRUE` only when the ordering of `query` is meaningful (ranked list).
#' @param organism Character; g:Profiler organism code (e.g., `"mmusculus"`, `"hsapiens"`).
#' @param user_threshold Numeric; significance threshold for `gost()` (used only if `sig_only = TRUE`).
#' @param correction_method Character; multiple-testing correction method for g:Profiler.
#'   Common values include `"g_SCS"` and `"fdr"` (see `?gprofiler2::gost`).
#' @param all_as_one_cluster Logical; if `TRUE` and `PART_results` is used, also runs enrichment
#'   on all genes combined (cluster `"All"`).
#' @param verbose Logical; if `TRUE`, prints progress and error messages.
#' @param ont Character; source/database to query (as used in `gost(sources=...)`), e.g. `"GO:BP"`,
#'   `"GO:MF"`, `"GO:CC"`, `"REAC"`. Must match `gost()` output `source`.
#' @param top_n Optional integer; if provided, keeps only the top `top_n` terms per cluster
#'   after all filtering steps (based on current row order).
#' @param sig_only Logical; passed to `gost(significant=...)`. If `TRUE`, `gost()` returns only
#'   significant terms (according to `user_threshold` and `correction_method`).
#' @param min_term_size Optional integer; remove terms with `term_size < min_term_size`.
#' @param max_term_size Optional integer; remove terms with `term_size > max_term_size`.
#' @param min_gene_ratio Optional numeric; remove terms with `Gene_Ratio < min_gene_ratio`,
#'   where `Gene_Ratio = intersection_size / query_size`.
#' @param min_fold_enrichment Optional numeric; remove terms with `log2(Fold_Enrichment) < min_fold_enrichment`.
#' @param min_genes Optional integer; remove terms with `intersection_size < min_genes`.
#'
#' @return A data.frame (`GO_df`) with one row per enriched term per cluster. Columns include:
#'   `term_id`, `term_name`, `p_value`, `-log10(padj)`, `term_size`, `query_size`,
#'   `intersection_size`, `effective_domain_size`, `Fold_Enrichment`,
#'   `log2(Fold_Enrichment)`, and `group_name` (cluster label).
#'
#' @details
#' Internally, the function:
#' - Determines `domain_scope` (`"annotated"` vs `"custom_annotated"`) from `background`.
#' - Runs `gprofiler2::gost()` with `tryCatch()` so failures in a cluster do not crash the whole run.
#' - Extracts only results for `ont` from the `gost()` result table.
#' - Computes:
#'   - `-log10(padj)` using the returned `p_value` column (note: naming kept as in your pipeline),
#'   - `Gene_Ratio = intersection_size / query_size`,
#'   - `Fold_Enrichment = (intersection_size * effective_domain_size) / (query_size * term_size)`,
#'   - `log2(Fold_Enrichment)`.
#' - Applies optional post-filters (`max_term_size`, `min_gene_ratio`, `min_fold_enrichment`, `min_genes`)
#'   and then `top_n` per cluster.
#'
#' @examples
#' \dontrun{
#' # 1) Using PART clustering results
#' GO_df <- run_gprofiler(PART_results = PART_res, ont = "GO:BP", correction_method = "g_SCS")
#'
#' # 2) Using a gene vector (single cluster "All")
#' GO_df <- run_gprofiler(gene_vector = c("Mapt", "App", "Bace1"), organism = "mmusculus", ont = "REAC")
#' }
#'
#' @seealso `gprofiler_wrapper()` for plotting, `gprofiler2::gost()`.
#' @export
run_gprofiler <- function(PART_results = NULL,
                          gene_vector = NULL,
                          background = NULL,
                          ordered_query = FALSE,
                          organism = "mmusculus",
                          user_threshold = 0.05,
                          correction_method = "gSCS",
                          all_as_one_cluster = FALSE,
                          verbose = TRUE,
                          ont = "GO:BP",
                          top_n = NULL,
                          sig_only = TRUE,
                          min_term_size = 10,
                          max_term_size = 1000,
                          min_gene_ratio = NULL,
                          min_fold_enrichment = NULL,
                          min_genes = 3) {

  suppressPackageStartupMessages({
    library(gprofiler2)
    library(dplyr)
    library(openxlsx)
    library(ggplot2)
  })

  # background can be: NULL, character vector, or matrix/data.frame (uses rownames)
  if (is.null(background)) {
    bg <- NULL
  } else if (is.matrix(background) || is.data.frame(background)) {
    bg <- rownames(background)
    if (is.null(bg)) stop("`background` was a matrix/data.frame but has NULL rownames(). Provide a character vector instead.")
  } else if (is.character(background)) {
    bg <- background
  } else {
    stop("`background` must be NULL, a character vector, or a matrix/data.frame with rownames().")
  }

  # define domain_scope in the parent scope (your pasted code passes domain_scope into run_gost)
  domain_scope <- if (is.null(bg)) "annotated" else "custom_annotated"

  # Helper to safely run gost
  run_gost <- function(query,
                       custom_bg = NULL,
                       organism = "mmusculus",
                       significant = sig_only,
                       domain_scope = NULL,
                       user_threshold = 0.05,
                       sources = ont,
                       correction_method = "gSCS",
                       ordered_query = FALSE,
                       verbose = TRUE) {

    if (is.null(domain_scope)) {
      domain_scope <- if (is.null(custom_bg)) "annotated" else "custom_annotated"
    }

    tryCatch({
      gprofiler2::gost(
        query = query,
        organism = organism,
        significant = significant,
        domain_scope = domain_scope,
        user_threshold = user_threshold,
        correction_method = correction_method,
        custom_bg = custom_bg,
        sources = sources,
        ordered_query = ordered_query,
        evcodes = TRUE
      )
    }, error = function(e) {
      if (verbose) message("g:Profiler failed: ", conditionMessage(e))
      return(NULL)
    })
  }

  if (!is.null(PART_results)) {
    if (!is.list(PART_results) || !"part_data" %in% names(PART_results)) {
      stop("PART_results must be a list produced by my_compute_PART_proteomics() containing 'part_data'.")
    }

    part_data <- PART_results$part_data
    cluster_map <- PART_results$cluster_map
    clusters <- unique(part_data$gene_cluster)
    if (all_as_one_cluster) clusters <- c(clusters, "All")

    gpro_res_list <- list()

    for (clust in clusters) {
      gene_vect <- if (clust == "All") rownames(part_data) else rownames(part_data)[part_data$gene_cluster == clust]
      gostres <- run_gost(query = gene_vect,
                          organism = organism,
                          significant = sig_only,
                          domain_scope = domain_scope,
                          ordered_query = ordered_query,
                          user_threshold = user_threshold,
                          correction_method = correction_method,
                          custom_bg = bg,
                          sources = ont)
      if (!is.null(gostres)) gpro_res_list[[clust]] <- gostres
    }

  } else if (!is.null(gene_vector)) {

    gostres <- run_gost(query = gene_vector,
                        organism = organism,
                        significant = sig_only,
                        ordered_query = ordered_query,
                        user_threshold = user_threshold,
                        correction_method = correction_method,
                        custom_bg = bg,
                        sources = ont
    )

    gpro_res_list <- list(All = gostres)
    cluster_map <- data.frame(cluster = "All", cluster_col = "#1f77b4")

  } else {
    stop("You must provide either PART_results or gene_vector.")
  }

  # ----- Extract and process results -----
  GO_df <- data.frame(NULL)

  for (clust in names(gpro_res_list)) {
    gprof_dta <- gpro_res_list[[clust]]$result
    if (is.null(gprof_dta)) next

    gprof_dta$parents <- sapply(gprof_dta$parents, function(x) paste(unlist(x), collapse = "/"))

    if (ont %in% gprof_dta$source) {
      extracted <- gprof_dta[gprof_dta$source == ont, c("term_id", "term_name", "p_value",
                                                        "term_size", "effective_domain_size",
                                                        "query_size", "intersection_size","intersection")]

      extracted[["-log10(padj)"]] <- -log10(extracted$p_value)
      extracted$group_color <- rep(unique(cluster_map$cluster_col[cluster_map$cluster == clust]), nrow(extracted))
      extracted$group_name <- clust
      extracted$Gene_Ratio <- with(extracted, ifelse(query_size == 0, NA_real_,
                                                     intersection_size / query_size))
      extracted$Fold_Enrichment <- with(extracted, ifelse(query_size == 0 | term_size == 0, NA_real_,
                                                          (intersection_size * effective_domain_size) /
                                                            (query_size * term_size)))
      extracted[["log2(Fold_Enrichment)"]] <- log2(extracted$Fold_Enrichment)
      GO_df <- rbind(GO_df, extracted)
    }
  }

  if (nrow(GO_df) == 0) {
    if (verbose) message("No enriched terms found for ontology ", ont)
    return(GO_df)
  }

  if (!is.null(min_term_size)) GO_df <- GO_df[GO_df$term_size >= min_term_size, ]
  if (!is.null(max_term_size)) GO_df <- GO_df[GO_df$term_size <= max_term_size, ]
  if (!is.null(min_gene_ratio)) GO_df <- GO_df[GO_df$Gene_Ratio >= min_gene_ratio, ]
  if (!is.null(min_fold_enrichment)) GO_df <- GO_df[GO_df$`log2(Fold_Enrichment)` >= min_fold_enrichment, ]
  if (!is.null(min_genes)) GO_df <- GO_df[GO_df$intersection_size >= min_genes, ]

  # Optional top_n filtering
  if (!is.null(top_n)) {
    GO_df <- do.call(rbind, lapply(split(GO_df, GO_df$group_name), function(df) head(df, top_n)))
  }

  # Add enrichment metrics / final columns
  GO_df <- GO_df[, c("term_id", "term_name", "p_value", "-log10(padj)",
                     "term_size", "query_size", "intersection_size",
                     "effective_domain_size", "Fold_Enrichment",
                     "log2(Fold_Enrichment)", "group_name","intersection")]

  # IMPORTANT: No plotting call here (as requested)
  return(GO_df)
}

#' Plot g:Profiler enrichment results from a `GO_df` table
#'
#' Plotting wrapper that takes a processed enrichment results table (`GO_df`)
#' (typically returned by `run_gprofiler()`) and generates:
#' - a multi-cluster dot plot (optional), and
#' - per-cluster dot/lollipop plots.
#'
#' The plots are saved to disk using `ggplot2::ggsave()` and your internal
#' plotting function `my_custom_gpro_dotplot()`.
#'
#' @param GO_df A data.frame of enrichment results, usually the output of `run_gprofiler()`.
#'   Must contain (at minimum): `term_id`, `term_name`, `-log10(padj)`, `group_name`,
#'   `intersection_size`, `query_size`, `term_size`, and `log2(Fold_Enrichment)`.
#' @param ont Character; ontology/source label used for titles and filenames
#'   (e.g., `"GO:BP"`, `"GO:MF"`, `"GO:CC"`, `"REAC"`). Colons are replaced with underscores in filenames.
#' @param top_n Optional integer; forwarded to `my_custom_gpro_dotplot()` for labeling/titles.
#' @param save_dir Output directory for plot files (default `"Gprofiler_results"`).
#' @param outfix Saving outfix.
#' @param exclude_zero_log10_padj Logical; if `TRUE`, removes rows where `-log10(padj) <= 0`
#'   before plotting.
#' @param lolipop Logical; when `FALSE` (default) per-cluster plots are dotplots; when `TRUE`,
#'   per-cluster plots include lollipop segments (depends on `my_custom_gpro_dotplot()` implementation).
#' @param multi_cluster Logical; if `TRUE`, also saves a combined multi-cluster dotplot.
#' @param term_size_legend Character; passed to `my_custom_gpro_dotplot()`. Typical values:
#'   `"percentage"` or `"numeric"`.
#' @param custom_width Optional numeric; override width passed to `ggsave()`.
#' @param custom_height Optional numeric; override height passed to `ggsave()`.
#'
#' @return Invisibly returns `TRUE` after saving plots.
#'
#' @details
#' Files written (if data are available):
#' - `<save_dir>/<ont>_dotplot.png` (multi-cluster plot; only if `multi_cluster = TRUE`)
#' - `<save_dir>/<ont>_lolipop_<cluster>.png` (one per cluster; filename keeps your original naming)
#'
#' The function expects `my_custom_gpro_dotplot()` to exist in the calling environment.
#'
#' @examples
#' \dontrun{
#' GO_df <- run_gprofiler(PART_results = PART_res, ont = "GO:BP")
#' gprofiler_wrapper(GO_df, ont = "GO:BP", lolipop = TRUE, multi_cluster = TRUE)
#' }
#'
#' @seealso `run_gprofiler()`.
#' @export
GO_plot_wrapper <- function(GO_df,
                            ont,
                            top_n = NULL,
                            save_dir = "",
                            outfix = "",
                            lolipop = FALSE,
                            multi_cluster = FALSE,
                            term_size_legend = "percentage",
                            custom_width = NULL,
                            custom_height = NULL,
                            ball_size_boost = 0,
                            color_low = "blue",
                            color_high = "red",
                            lolipop_linewidth = 1) {

  suppressPackageStartupMessages({
    library(ggplot2)
  })

    #-----------------------------
  # 0️⃣ Directory setup
  #-----------------------------
  if (is.null(save_dir) || save_dir == "") {
    save_dir <- "Enrichment_plots"
  }

  # 3. Create the directory if it doesn't exist
  if (!dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
  }

  # ---------------- Plotting Section ----------------------------------------
  # (Copied from your pasted function with only the surrounding function wrapper added)

  # filter and order
  GO_top_cluster <- GO_df
  ont_clean <- gsub(":", "_", ont)

  if (nrow(GO_top_cluster) > 0) {
    # MULTI CLUSTER plot
    if (multi_cluster) {
      GO_top_cluster <- GO_top_cluster[order(GO_top_cluster$term_id, -GO_top_cluster$`-log10(padj)`), ]
      res_multi  <- my_custom_gpro_dotplot(GO_top_cluster,
                                          ont,
                                          top_n,
                                          term_size_legend,
                                          sub_cluster = NULL,
                                          lolipop = FALSE,
                                          multi_cluster = TRUE,
                                          ball_size_boost = ball_size_boost,
                                          color_low = color_low,
                                          color_high = color_high,
                                          lolipop_linewidth = lolipop_linewidth)
      plt_multi  <- res_multi$plot

      avg_term_length <- mean(nchar(as.character(unique(GO_top_cluster$term_id))))
      width <- ifelse(is.null(custom_width),
                      6 + (length(unique(GO_top_cluster$group_name)) * 0.3) + (avg_term_length * 0.02),
                      custom_width)
      #height <- ifelse(is.null(custom_height), 3 + nrow(GO_top_cluster) * 0.25, custom_height)
      height <- ifelse(is.null(custom_height), res_multi$height, custom_height)
      ggsave(file.path(save_dir, paste0(ont_clean,"_",outfix,"_dotplot.png")),
             plot = plt_multi, width = width, height = height, limitsize = FALSE)
    }

    # SUB-CLUSTER plots
    cluster_list <- unique(GO_top_cluster$group_name)
    for (cl in cluster_list) {
      cluster_data <- GO_top_cluster[GO_top_cluster$group_name == cl, ]
      if (nrow(cluster_data) == 0) next

      cluster_data <- cluster_data[order(cluster_data$term_id, -cluster_data$`-log10(padj)`), ]
      res_cl  <- my_custom_gpro_dotplot( GO_clusters = cluster_data,
                                         ont = ont,
                                         top_n = top_n,
                                         term_size_legend = term_size_legend,
                                         sub_cluster = cl,
                                         lolipop = lolipop,
                                         multi_cluster = FALSE,
                                         ball_size_boost = ball_size_boost,
                                         color_low = color_low,
                                         color_high = color_high,
                                         lolipop_linewidth = lolipop_linewidth)
      plt_cl  <- res_cl$plot

      #width <- ifelse(is.null(custom_width), 6 + length(unique(cluster_data$group_name)) * 0.6, custom_width)
      x_max <- max(cluster_data$`-log10(padj)`, na.rm = TRUE)
      max_text_length <- max(nchar(as.character(cluster_data$term_id)))
      width <- ifelse(is.null(custom_width), 8 + log(x_max + 1) * 5 + (max_text_length * 0.05), custom_width)
      #height <- ifelse(is.null(custom_height), 5 + nrow(cluster_data) * 0.1, custom_height)
      height <- ifelse(is.null(custom_height), res_cl$height, custom_height)
      ggsave(file.path(save_dir, paste0(ont_clean, "_lolipop_", cl,"_",outfix, ".png")),
             plot = plt_cl, width = width, height = height, limitsize = FALSE)
    }
  }

  invisible(TRUE)
}

#' Custom g:Profiler dotplot for single-cluster or multi-cluster enrichment results
#'
#' Build a `ggplot2` dotplot for a g:Profiler enrichment results table. The plot operates
#' in two modes:
#' - Single-cluster mode (`multi_cluster = FALSE`): x-axis is `-log10(padj)` and optional
#'   lollipop segments can be drawn from 0 to each point.
#' - Multi-cluster mode (`multi_cluster = TRUE`): x-axis is `group_name` (cluster label),
#'   y-axis is `term_name`, and points encode term size and significance across clusters.
#'
#' The function expects an input data frame that matches the structure produced by your
#' `run_gprofiler()` pipeline (e.g., includes `term_name`, `group_name`, `-log10(padj)`,
#' `term_size`, `intersection_size`, `query_size`, and `log2(Fold_Enrichment)`).
#'
#' @param GO_clusters A data.frame of enrichment results. Required columns:
#'   `term_name`, `group_name`, `-log10(padj)`, `term_size`, `intersection_size`,
#'   `query_size`, `log2(Fold_Enrichment)`. Additional columns are ignored.
#' @param ont Character. Ontology/source label used in the plot title (e.g., `"GO:BP"`, `"GO:MF"`, `"GO:CC"`, `"REAC"`).
#' @param top_n Integer or NULL. Displayed in the plot title as `top_n` (no internal filtering is performed here).
#' @param term_size_legend Character. Controls the point size encoding:
#'   - `"numeric"`: point size maps to `term_size`.
#'   - `"percentage"`: point size maps to `term_size_percentage`, computed as
#'     `(intersection_size / query_size) * 100`.
#' @param sub_cluster Character or NULL. If non-NULL, filters `GO_clusters` to rows
#'   where `group_name == sub_cluster` before plotting.
#' @param lolipop Logical. Only used when `multi_cluster = FALSE`.
#'   If TRUE, draws a `geom_segment()` from x = 0 to x = `-log10(padj)` for each term.
#' @param multi_cluster Logical. If TRUE, builds a multi-cluster dotplot with x = `group_name`.
#'   If FALSE, builds a single-cluster dotplot with x = `-log10(padj)`.
#'
#' @details
#' Aesthetics and encodings:
#' - y-axis: `term_name` (ordered by decreasing `-log10(padj)`).
#' - Point size:
#'   - `"percentage"`: `100 * intersection_size / query_size` (labelled "Gene Ratio (%)").
#'   - `"numeric"`: `term_size` (labelled "Gene Count").
#' - Point color:
#'   - `multi_cluster = TRUE`: rounds `-log10(padj)` to 1 decimal and uses it for color.
#'   - `multi_cluster = FALSE`: rounds `log2(Fold_Enrichment)` to 1 decimal and uses it for color
#'     (legend title "Fold Enrichment" as currently implemented).
#'
#' In single-cluster mode, the x-axis is forced to start at 0 and uses integer breaks up to
#' `ceiling(max(-log10(padj)))`. `coord_cartesian(clip = "off")` is used to avoid clipping.
#'
#' Dependencies:
#' - Uses `ggplot2`.
#' - Uses `stringr::str_wrap()` to wrap y-axis labels; ensure `stringr` is installed/available.
#'
#' @return A `ggplot` object.
#'
#' @examples
#' \dontrun{
#' # Single-cluster dotplot (optionally lollipop)
#' p1 <- my_custom_gpro_dotplot(GO_clusters = GO_df,
#'                              ont = "GO:BP",
#'                              top_n = 20,
#'                              term_size_legend = "percentage",
#'                              sub_cluster = "1",
#'                              lolipop = TRUE,
#'                              multi_cluster = FALSE)
#' print(p1)
#'
#' # Multi-cluster dotplot
#' p2 <- my_custom_gpro_dotplot(GO_clusters = GO_df,
#'                              ont = "GO:BP",
#'                              top_n = 20,
#'                              term_size_legend = "numeric",
#'                              sub_cluster = NULL,
#'                              lolipop = FALSE,
#'                              multi_cluster = TRUE)
#' print(p2)
#' }
#'
#' @export
my_custom_gpro_dotplot <- function(GO_clusters,
                                   ont,
                                   top_n,
                                   term_size_legend,
                                   sub_cluster,
                                   lolipop=FALSE,
                                   multi_cluster,
                                   ball_size_boost = 1,
                                   color_low = "blue",
                                   color_high = "red",
                                   lolipop_linewidth = 1) {

  # Filter by sub-cluster if provided
  if (!is.null(sub_cluster)) {
    GO_clusters <- GO_clusters[GO_clusters$group_name == sub_cluster, ]
  }

  # Prepare term names and title
  GO_clusters$term_name <- factor(GO_clusters$term_name, levels = unique(GO_clusters$term_name[order(GO_clusters$`-log10(padj)`)]))
  title_text <- paste0(ont)

  if (!is.null(top_n)) {
    actual_n <- min(top_n, nrow(GO_clusters))
    title_text <- paste0(title_text, " (top ", actual_n, ")")
    GO_clusters <- GO_clusters[1:actual_n, , drop = FALSE]
  }

  if (!is.null(sub_cluster) && sub_cluster != "All") {
    title_text <- paste(title_text, "| Cluster:", sub_cluster)
  }

  # Optionally compute percentage
  if (term_size_legend == "percentage") {
    GO_clusters$term_size_percentage <- (GO_clusters$intersection_size / GO_clusters$query_size) * 100
  }

  # Point aesthetic mapping
  size_var <- if (term_size_legend == "percentage") "term_size_percentage" else "intersection_size"
  color_var <- if (multi_cluster) "-log10(padj)" else "log2(Fold_Enrichment)"
  color_title <- if (multi_cluster) "-log10(padj)" else "log2(Fold Enrichment)"

  if (term_size_legend == "percentage") {
    has_sub1 <- any(GO_clusters[[size_var]] < 1, na.rm = TRUE)
    GO_clusters[[size_var]] <- pmax(GO_clusters[[size_var]], 1)
  }
  # ---------------------------------
  if (multi_cluster) {
    plt <- ggplot(GO_clusters, aes(x = group_name, y = term_name))
  } else {
    plt <- ggplot(GO_clusters, aes(x = `-log10(padj)`, y = term_name))
    if (lolipop) {
      plt <- plt + geom_segment(aes(x = 0, xend = `-log10(padj)`, y = term_name, yend = term_name, color = round(`log2(Fold_Enrichment)`, 1)), linewidth = lolipop_linewidth, alpha = 0)
    }
  }
  # NOW add geom_point (data is already clamped)
  plt <- plt + geom_point(aes(size = .data[[size_var]], color = round(.data[[color_var]], 1)), show.legend = TRUE)

  if (term_size_legend == "percentage") {
    # has_sub1 already set above
    unique_sizes <- sort(unique(GO_clusters[[size_var]]))  # already clamped
    if (length(unique_sizes) > 5) {
      size_breaks <- unique_sizes[round(seq(1, length(unique_sizes), length.out = 5))]
    } else {
      size_breaks <- unique_sizes
    }

    size_labels <- round(size_breaks, 1)
    if (has_sub1) size_labels[1] <- "<1"

    plt <- plt + scale_size(
      name   = "Gene Ratio (%)",
      limits = c(1, 100),
      breaks = size_breaks,
      labels = size_labels,
      range  = c(1, 10) * ball_size_boost,
      guide  = guide_legend(title.position = "top")
    )
  } else {
    unique_sizes <- sort(unique(GO_clusters[[size_var]]))
    if (length(unique_sizes) > 5) {
      size_breaks <- unique_sizes[round(seq(1, length(unique_sizes), length.out = 5))]
    } else {
      size_breaks <- unique_sizes
    }
    plt <- plt + scale_size(
      name   = "Gene Hits",
      limits = range(size_breaks, na.rm = TRUE),
      breaks = size_breaks,
      labels = round(size_breaks, 1),
      range  = c(1, 10) * ball_size_boost,
      guide  = guide_legend(title.position = "top")
    )
  }

  # Handle color scaling
  color_range <- range(GO_clusters[[color_var]])
  if (length(unique(GO_clusters[[color_var]])) > 1) {
    plt <- plt + scale_color_gradientn(
      colours = c(color_low, color_high),
      limits = c(floor(color_range[1]), ceiling(color_range[2])),
      name = color_title
    )
  } else {
    plt <- plt + labs(color = color_title)
  }

  # --- Compute wrapped labels and recommended figure height ---
  term_levels   <- levels(GO_clusters$term_name)
  wrapped_labels <- stringr::str_wrap(term_levels, width = 80)
  line_counts   <- stringr::str_count(wrapped_labels, "\n") + 1
  names(wrapped_labels) <- term_levels
  #max_lines     <- max(line_counts)
  #recommended_height <- length(term_levels) * 0.1 * max_lines + 2
  total_lines <- sum(line_counts)           # actual total lines across all labels
  recommended_height <- total_lines * 0.3 + 2   # 0.3 inches per line + 2 for title/legend

  # Shared plot theming
  plt <- plt +
    ggtitle(title_text) +
    theme_bw(base_size = 15) +
    theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)) +
    guides(size = guide_legend(order = 1), color = guide_colorbar(order = 2)) +
    #scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60))
    scale_y_discrete(labels = wrapped_labels) +
    theme(axis.text.y = element_text(lineheight = 0.8,size = 20))

  if (multi_cluster) {
    plt <- plt + theme(axis.text.x = element_text(angle = -45, vjust = 0)) + xlab("") + ylab("")
  } else {
    # --- CHANGED SECTION START ---
    xmax <- ceiling(max(GO_clusters$`-log10(padj)`, na.rm = TRUE))

    plt <- plt +
      scale_x_continuous(
        breaks = seq(0, xmax, by = 1),
        limits = c(0, xmax),
        expand = c(0, 0)  # remove empty space at start of x-axis
      ) +
      coord_cartesian(clip = "off") +  # ensure nothing gets clipped
      xlab("-log10(padj)") + ylab("")
  }

  return(list(plot = plt, height = recommended_height))
}
