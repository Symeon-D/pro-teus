# run_fgsea_contrasts.R
# Gene Set Enrichment Analysis for limma eBayes contrasts
# Requires: fgsea, msigdbr, limma, data.table, dplyr

# ─────────────────────────────────────────────────────────────────────────────

#' Run fgsea on one or more limma contrasts
#'
#' Extracts a ranked protein list (by moderated t-statistic) from a fitted
#' limma eBayes object and runs fast pre-ranked GSEA via \code{fgsea::fgsea()}.
#'
#' Ranking by t-statistic is preferred over logFC because it incorporates
#' variance information and is more robust to high-variance proteins with large
#' but unreliable fold changes (Subramanian et al. 2005, PNAS 102:15545;
#' Korotkevich et al. 2021, bioRxiv fgsea preprint).
#'
#' @section Contrast selection:
#' \itemize{
#'   \item \strong{All contrasts:} leave \code{coef = NULL} (default). Every
#'         column in \code{fit$coefficients} is run.
#'   \item \strong{Subset:} pass a character vector, e.g.
#'         \code{coef = c("nGRKO_Str_vs_Con", "mGRKO_Str_vs_Con")}.
#'   \item \strong{Single contrast:} pass one string, e.g.
#'         \code{coef = "nGRKO_Str_vs_Con"}.
#' }
#' In all cases the return value is a named list (one element per contrast),
#' so downstream code does not need to change based on how many contrasts
#' were run.
#'
#' @section Pathway databases for GO (mouse):
#' \itemize{
#'   \item \strong{GO:BP} \code{category="C5", subcategory="GO:BP"} —
#'         broadest coverage. Many terms share member genes so a significant
#'         term does not constitute independent evidence from an overlapping
#'         term (Lim et al. 2021, Bioinformatics 37:4593).
#'   \item \strong{GO:MF} \code{subcategory="GO:MF"} — molecular function;
#'         useful for identifying affected receptor or enzyme classes.
#'   \item \strong{GO:CC} \code{subcategory="GO:CC"} — cellular component;
#'         useful for confirming synaptic localisation of hits.
#' }
#'
#' @section Critical caveats:
#' \itemize{
#'   \item \strong{No correction across contrasts.} BH correction is applied
#'         within each contrast independently. Running 7 contrasts inflates the
#'         family-wise false discovery rate across the full analysis
#'         (Bourgon et al. 2010, PNAS 107:10734). A pathway significant in
#'         multiple contrasts is more credible, but this is not formally tested.
#'   \item \strong{GO redundancy.} GO terms are hierarchical and heavily
#'         overlapping. Reporting 30 significant GO:BP terms that all share
#'         the same 15 leading edge genes overstates the breadth of the finding.
#'         Use \code{fgsea::collapsePathways()} to reduce redundancy before
#'         reporting (Korotkevich et al. 2021).
#'   \item \strong{NES comparability.} NES is comparable across gene sets
#'         within one contrast. Cross-contrast NES comparison is only valid
#'         when ranked lists have similar length — true here since all
#'         contrasts use the same 3166-protein matrix.
#'   \item \strong{Leading edge size.} A set can reach significance with very
#'         few leading edge genes if they cluster at the ranking extreme.
#'         Always inspect \code{leadingEdge} length relative to \code{size}.
#' }
#'
#' @param fit         A \code{MArrayLM} object returned by
#'                    \code{limma::eBayes()}, with contrasts already applied
#'                    via \code{limma::contrasts.fit()}.
#' @param pathways    Named list of character vectors. Each element is a gene
#'                    set (MGI symbols matching \code{rownames(fit)}). Produced
#'                    via \code{msigdbr} + \code{split()}. See examples.
#' @param coef        Character vector of contrast names, or \code{NULL}
#'                    (default) to run all contrasts in \code{fit$coefficients}.
#' @param min_size    Integer. Minimum gene set size after intersection with
#'                    ranked list. Default 10. Lower to 5 for small contrasts
#'                    at the cost of testing more unstable small sets.
#' @param max_size    Integer. Maximum gene set size to test. Default 500.
#'                    Large sets tend to be biologically non-specific
#'                    (Lim et al. 2021, Bioinformatics 37:4593).
#' @param eps         Numeric. Lower p-value boundary. Set to \code{0} for
#'                    arbitrarily small p-values via Monte Carlo (slower).
#'                    Default \code{1e-10}.
#' @param seed        Integer. Random seed for reproducibility. Default 42.
#' @param nperm       Integer. Permutations for p-value estimation
#'                    (\code{nPermSimple} in fgsea). Default 10000.
#' @param save_dir    Character or \code{NULL}. If provided, each result is
#'                    written to \code{<save_dir>/<contrast>_fgsea.csv}.
#'                    \code{leadingEdge} is collapsed to semicolon-separated
#'                    string for CSV compatibility. Default \code{NULL}.
#' @param padj_filter Numeric or \code{NULL}. If set, each returned table is
#'                    filtered to \code{padj <= padj_filter}. Set to
#'                    \code{NULL} to keep all tested pathways. Default
#'                    \code{NULL}.
#'
#' @return Named list of \code{data.table} objects, one per contrast. Names
#'   match the \code{coef} values used. Each table contains:
#'   \describe{
#'     \item{pathway}{Gene set name}
#'     \item{pval}{Nominal p-value}
#'     \item{padj}{BH-adjusted p-value within this contrast}
#'     \item{log2err}{Expected error of the log2(p-value) estimate}
#'     \item{ES}{Raw enrichment score}
#'     \item{NES}{Normalised enrichment score}
#'     \item{size}{Genes from set present in ranked list}
#'     \item{leadingEdge}{List column: genes driving the enrichment}
#'     \item{contrast}{Contrast name, added for easy \code{rbind}}
#'   }
#'   Sorted by \code{pval} ascending within each table.
#'
#' @examples
#' \dontrun{
#' library(fgsea)
#' library(msigdbr)
#' library(dplyr)
#'
#' # --- Prepare GO pathway databases (mouse) --------------------------------
#' pathways_BP <- msigdbr(species = "Mus musculus",
#'                         category = "C5", subcategory = "GO:BP") %>%
#'   split(x = .$gene_symbol, f = .$gs_name)
#'
#' pathways_MF <- msigdbr(species = "Mus musculus",
#'                         category = "C5", subcategory = "GO:MF") %>%
#'   split(x = .$gene_symbol, f = .$gs_name)
#'
#' pathways_CC <- msigdbr(species = "Mus musculus",
#'                         category = "C5", subcategory = "GO:CC") %>%
#'   split(x = .$gene_symbol, f = .$gs_name)
#'
#' # --- All contrasts -------------------------------------------------------
#' fgsea_BP <- run_fgsea(
#'   fit      = eBayes_fit,
#'   pathways = pathways_BP,
#'   save_dir = "results/fgsea_GOBP"
#' )
#'
#' # --- Subset: stress contrasts only ---------------------------------------
#' fgsea_stress <- run_fgsea(
#'   fit      = eBayes_fit,
#'   pathways = pathways_BP,
#'   coef     = c("GRfl_fl_Str_vs_Con", "nGRKO_Str_vs_Con", "mGRKO_Str_vs_Con")
#' )
#'
#' # --- Single contrast -----------------------------------------------------
#' fgsea_nGRKO <- run_fgsea(
#'   fit      = eBayes_fit,
#'   pathways = pathways_BP,
#'   coef     = "nGRKO_Str_vs_Con"
#' )
#'
#' # --- Inspect results -----------------------------------------------------
#'
#' # Significant upregulated GO:BP for nGRKO stress
#' fgsea_stress[["nGRKO_Str_vs_Con"]][padj < 0.05 & NES > 0][order(-NES)]
#'
#' # Pathways shared between nGRKO and WT stress
#' sig_nGRKO <- fgsea_stress[["nGRKO_Str_vs_Con"]][padj < 0.05, pathway]
#' sig_WT    <- fgsea_stress[["GRfl_fl_Str_vs_Con"]][padj < 0.05, pathway]
#' intersect(sig_nGRKO, sig_WT)   # shared
#' setdiff(sig_nGRKO, sig_WT)     # nGRKO-specific
#'
#' # Reduce GO redundancy before reporting
#' res_nGRKO   <- fgsea_stress[["nGRKO_Str_vs_Con"]]
#' ranks_nGRKO <- sort(
#'   setNames(
#'     limma::topTable(eBayes_fit, coef = "nGRKO_Str_vs_Con",
#'                     number = Inf, sort.by = "t")$t,
#'     rownames(limma::topTable(eBayes_fit, coef = "nGRKO_Str_vs_Con",
#'                              number = Inf, sort.by = "t"))
#'   ), decreasing = TRUE
#' )
#' collapsed <- fgsea::collapsePathways(
#'   res_nGRKO[padj < 0.05], pathways_BP, ranks_nGRKO
#' )
#' res_nGRKO[pathway %in% collapsed$mainPathways][order(pval)]
#'
#' # Summary across all contrasts
#' summarise_fgsea(fgsea_BP)
#' }
#'
#' @seealso \code{\link{summarise_fgsea}}, \code{\link[fgsea]{fgsea}},
#'   \code{\link[fgsea]{collapsePathways}}
#' @export
run_fgsea <- function(fit,
                      pathways,
                      coef        = NULL,
                      min_size    = 10,
                      max_size    = 500,
                      eps         = 1e-10,
                      seed        = 42,
                      nperm       = 10000,
                      save_dir    = NULL,
                      padj_filter = NULL) {
  
  # --- Input checks -----------------------------------------------------------
  if (!inherits(fit, "MArrayLM")) {
    stop("`fit` must be a MArrayLM object from limma::eBayes().")
  }
  if (!is.list(pathways) || is.null(names(pathways))) {
    stop("`pathways` must be a named list of character vectors.")
  }
  
  all_coefs <- colnames(fit$coefficients)
  
  if (is.null(coef)) {
    coefs_to_run <- all_coefs
    message("No `coef` specified — running all ", length(coefs_to_run),
            " contrasts.")
  } else {
    missing_coefs <- setdiff(coef, all_coefs)
    if (length(missing_coefs) > 0) {
      stop(sprintf(
        "Contrast(s) not found in fit$coefficients: %s\nAvailable: %s",
        paste(missing_coefs, collapse = ", "),
        paste(all_coefs, collapse = ", ")
      ))
    }
    coefs_to_run <- coef
    message("Running fgsea for ", length(coefs_to_run), " contrast(s): ",
            paste(coefs_to_run, collapse = ", "))
  }
  
  if (!is.null(save_dir) && !dir.exists(save_dir)) {
    dir.create(save_dir, recursive = TRUE)
    message("Created output directory: ", save_dir)
  }
  
  # --- Main loop --------------------------------------------------------------
  results <- vector("list", length(coefs_to_run))
  names(results) <- coefs_to_run
  
  for (cf in coefs_to_run) {
    message("\n─── ", cf, " ───")
    
    # Build ranked vector from moderated t-statistics
    tt <- limma::topTable(fit, coef = cf, number = Inf,
                          adjust.method = "BH", p.value = 1, sort.by = "t")
    
    n_na <- sum(is.na(tt$t))
    if (n_na > 0) {
      warning(sprintf(
        "%d proteins with NA t-statistics in '%s' excluded from ranking.",
        n_na, cf
      ))
      tt <- tt[!is.na(tt$t), ]
    }
    
    ranks <- tt$t
    names(ranks) <- rownames(tt)
    
    if (anyDuplicated(names(ranks))) {
      warning("Duplicate rownames in '", cf, "' — keeping first occurrence.")
      ranks <- ranks[!duplicated(names(ranks))]
    }
    
    ranks <- sort(ranks, decreasing = TRUE)
    
    # Run fgsea
    set.seed(seed)
    res <- fgsea::fgsea(
      pathways    = pathways,
      stats       = ranks,
      minSize     = min_size,
      maxSize     = max_size,
      eps         = eps,
      nPermSimple = nperm
    )
    
    res[, contrast := cf]
    res <- res[order(pval)]
    
    # Console summary
    n_sig_05 <- sum(res$padj < 0.05, na.rm = TRUE)
    n_sig_01 <- sum(res$padj < 0.01, na.rm = TRUE)
    n_up     <- sum(res$padj < 0.05 & res$NES > 0, na.rm = TRUE)
    n_down   <- sum(res$padj < 0.05 & res$NES < 0, na.rm = TRUE)
    message(sprintf(
      "  Tested: %d | padj<0.05: %d (%d\u2191 / %d\u2193) | padj<0.01: %d",
      nrow(res), n_sig_05, n_up, n_down, n_sig_01
    ))
    
    if (!is.null(padj_filter)) {
      res <- res[padj <= padj_filter]
    }
    
    if (!is.null(save_dir)) {
      out_path <- file.path(save_dir, paste0(cf, "_fgsea.csv"))
      res_save <- data.table::copy(res)
      res_save[, leadingEdge := sapply(leadingEdge, paste, collapse = ";")]
      data.table::fwrite(res_save, out_path)
      message("  Saved: ", out_path)
    }
    
    results[[cf]] <- res
  }
  
  message("\n\u2713 Done. ", length(coefs_to_run), " contrast(s) complete.")
  return(results)
}


# ─────────────────────────────────────────────────────────────────────────────

#' Summarise fgsea results across contrasts
#'
#' Prints and returns a compact count table of significant pathways per
#' contrast. Useful as a first overview before inspecting individual results.
#'
#' @param fgsea_list  Named list of \code{data.table} objects as returned by
#'                    \code{\link{run_fgsea}}.
#' @param padj_thresh Numeric. Adjusted p-value threshold. Default 0.05.
#'
#' @return A \code{data.frame} with columns \code{contrast}, \code{n_tested},
#'   \code{n_sig}, \code{n_up} (NES > 0), \code{n_down} (NES < 0). Returned
#'   invisibly and printed to console.
#'
#' @examples
#' \dontrun{
#' summarise_fgsea(fgsea_BP)
#' summarise_fgsea(fgsea_BP, padj_thresh = 0.01)
#' }
#' @export
summarise_fgsea <- function(fgsea_list, padj_thresh = 0.05) {
  
  summary_df <- do.call(rbind, lapply(names(fgsea_list), function(nm) {
    res <- fgsea_list[[nm]]
    data.frame(
      contrast = nm,
      n_tested = nrow(res),
      n_sig    = sum(res$padj < padj_thresh, na.rm = TRUE),
      n_up     = sum(res$padj < padj_thresh & res$NES > 0, na.rm = TRUE),
      n_down   = sum(res$padj < padj_thresh & res$NES < 0, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  }))
  
  cat(sprintf("\n=== fgsea summary (padj < %.2f) ===\n", padj_thresh))
  print(summary_df, row.names = FALSE)
  cat("\n")
  
  invisible(summary_df)
}