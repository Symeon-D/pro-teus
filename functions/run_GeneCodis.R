#' Functional Enrichment via GeneCodis 4
#'
#' @param PART_results List containing 'part_data' (from PART analysis).
#' @param gene_vector Character vector of gene symbols (if not using PART).
#' @param background Optional custom background gene universe (character vector,
#'   or matrix/data.frame whose rownames are used as the universe).
#' @param organism Character string in GeneCodis format, e.g., "Mus musculus"
#'   or "Homo sapiens". See genecodis.genyo.es for all supported organisms.
#' @param ont Character vector of databases (e.g., "GO_BP", "KEGG_pathways").
#' @param sig_threshold Numeric p-value threshold (default 0.05).
#' @param min_genes Minimum genes in intersection (default 5).
#' @param max_term_size Optional integer; remove terms with `term_size > max_term_size`.
#' @param top_n Optional integer; keep only top n terms per group by padj.
#' @param verbose Logical, print progress messages.
#'
#' @return A processed data frame of enrichment results.
run_Gene_Codis <- function(PART_results = NULL,
                           gene_vector = NULL,
                           background = NULL,
                           organism = "Mus musculus",
                           ont = "GO_BP",
                           sig_threshold = 0.05,
                           min_genes = 5,
                           max_term_size = 1000,
                           top_n = NULL,
                           verbose = TRUE) {

  # -- Patch and source GeneCodis script --------------------------------------
  gc_lines <- readLines("https://raw.githubusercontent.com/GENyO-BioInformatics/GeneCodis/public_master/GC4libR.R")
  gc_lines <- gsub('^library\\(.*\\)', '', gc_lines)
  gc_lines <- gsub('content\\(request, *"text"\\)',
                   'httr::content(request, as="text", encoding="UTF-8")',
                   gc_lines)
  gc_lines <- gsub('content\\(request, *"parsed"\\)',
                   'httr::content(request, as="parsed")',
                   gc_lines)
  gc_lines <- gsub('(?<!:)\\bcontent\\(',
                   'httr::content(',
                   gc_lines, perl = TRUE)
  suppressPackageStartupMessages({
    library(jsonlite); library(httr); library(dplyr); library(rlang)
  })
  tmp <- tempfile(fileext = ".R")
  writeLines(gc_lines, tmp)
  source(tmp, local = FALSE)
  # ---------------------------------------------------------------------------

  suppressPackageStartupMessages({
    library(jsonlite)
    library(httr)
    library(dplyr)
    library(rlang)
  })

  ont <- gsub(":", "_", ont)
  if (verbose) message("Using Ontologies: ", paste(ont, collapse = ", "))

  # --- Validate and resolve background ---------------------------------------
  # background can be: NULL, character vector, or matrix/data.frame (uses rownames)
  # launchAnalysis expects `customUniverse` as a list(); empty list() = no custom universe
  if (is.null(background)) {
    bg <- list()
  } else if (is.matrix(background) || is.data.frame(background)) {
    bg <- rownames(background)
    if (is.null(bg)) stop("`background` was a matrix/data.frame but has NULL rownames(). Provide a character vector instead.")
    bg <- as.list(bg)
  } else if (is.character(background)) {
    bg <- as.list(background)
  } else {
    stop("`background` must be NULL, a character vector, or a matrix/data.frame with rownames().")
  }
  universe_scope <- "annotated"
  if (length(bg) > 0 && verbose) message("Using custom background universe of ", length(bg), " genes.")
  # ---------------------------------------------------------------------------

  # --- Internal helper -------------------------------------------------------
  # NOTE: no extra arguments -- bg and universe_scope are closed over from above
  run_gc_single <- function(query, species, dbs) {
    tryCatch({
      launchAnalysis(
        organism          = species,
        inputType         = "genes",
        inputQuery        = query,
        customUniverse    = bg,
        annotationsDBs    = dbs,
        enrichmentStat    = "hypergeom",
        universeScope     = universe_scope,  # "annotated" or "whole" (custom universe is passed separately via customUniverse)
        inputCoannotation = "no",
        ReportName        = "input1"
      )
    }, error = function(e) {
      if (verbose) message("GeneCodis API failed: ", conditionMessage(e))
      return(NULL)
    })
  }
  # ---------------------------------------------------------------------------

  # --- Input Handling --------------------------------------------------------
  gc_res_list <- list()

  if (!is.null(PART_results)) {
    part_data <- PART_results$part_data
    clusters  <- unique(part_data$gene_cluster)

    for (clust in clusters) {
      if (verbose) message("Running GeneCodis for cluster: ", clust)
      gene_vect <- rownames(part_data)[part_data$gene_cluster == clust]

      res <- run_gc_single(gene_vect, organism, ont)
      if (!is.null(res)) gc_res_list[[as.character(clust)]] <- res

      # Respect rate limit (10 requests/min)
      if (length(clusters) > 1) Sys.sleep(7)
    }
  } else if (!is.null(gene_vector)) {
    if (verbose) message("Running GeneCodis for gene_vector...")
    res <- run_gc_single(gene_vector, organism, ont)
    if (!is.null(res)) gc_res_list[["All"]] <- res
  } else {
    stop("Provide either PART_results or gene_vector.")
  }
  # ---------------------------------------------------------------------------

  # --- Data Processing -------------------------------------------------------
  final_df <- data.frame()

  for (group in names(gc_res_list)) {
    table_key <- paste0("input1-", ont)
    raw_table <- gc_res_list[[group]]$stats_tables[[table_key]]

    if (is.null(raw_table) || nrow(raw_table) == 0) next

    processed <- raw_table %>%
      rename(
        term_id           = "annotation_id",
        term_name         = "description",
        p_value           = "pval",
        p_val_adj         = "pval_adj",
        query_size        = "input_size",
        intersection_size = "genes_found",
        term_size         = "term_genes",
        intersection      = "genes"
      ) %>%
      mutate(
        group_name              = group,
        p_value                 = as.numeric(p_value),
        p_val_adj               = as.numeric(p_val_adj),
        intersection_size       = as.numeric(intersection_size),
        query_size              = as.numeric(query_size),
        term_size               = as.numeric(term_size),
        `-log10(padj)`          = -log10(p_val_adj),
        Fold_Enrichment         = as.numeric(relative_enrichment),
        `log2(Fold_Enrichment)` = log2(Fold_Enrichment)
      ) %>%
      filter(
        p_val_adj         < sig_threshold,
        intersection_size >= if (is.null(min_genes)) 0 else min_genes,
        term_size         <  if (is.null(max_term_size)) Inf else max_term_size
      )

    final_df <- rbind(final_df, processed)
  }
  # ---------------------------------------------------------------------------

  if (nrow(final_df) == 0) {
    if (verbose) message("No significant terms found (padj < ", sig_threshold, ").")
    return(final_df)
  }

  if (!is.null(top_n)) {
    final_df <- final_df %>%
      arrange(group_name, p_val_adj) %>%
      group_by(group_name) %>%
      slice_head(n = top_n) %>%
      ungroup()
  }

  final_df <- final_df[, c("term_id", "term_name", "p_value", "-log10(padj)",
                           "term_size", "intersection_size", "query_size",
                           "Fold_Enrichment", "log2(Fold_Enrichment)", "group_name","intersection")]

  return(final_df)
}
