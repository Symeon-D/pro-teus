#' Run GO ORA using clusterProfiler with gprofiler-compatible output
#'
#' This function performs Over-Representation Analysis using clusterProfiler
#' and returns results in the same format as run_gprofiler().
#'
#' @param gene_vector Character vector of gene/protein IDs to test.
#' @param background Character vector of background gene/protein IDs.
#' @param ont Character. GO ontology: "GO:BP", "GO:CC", or "GO:MF".
#' @param organism Character. Organism database. Default is "org.Mm.eg.db" (mouse).
#'   Use "org.Hs.eg.db" for human, "org.Rn.eg.db" for rat.
#' @param keyType Character. Type of gene ID. Default is "SYMBOL".
#'   Other options: "ENSEMBL", "ENTREZID", "UNIPROT".
#' @param pvalueCutoff Numeric. P-value cutoff. Default is 0.05.
#' @param qvalueCutoff Numeric. FDR cutoff. Default is 0.05.
#' @param minGSSize Integer. Minimum gene set size. Default is 10.
#' @param maxGSSize Integer. Maximum gene set size. Default is 500.
#' @param group_name Character. Label for this gene set. Default is "All".
#' @param top_n Integer or NULL. If provided, returns only the top N most significant
#'   terms (ordered by adjusted p-value). Default is NULL (return all significant terms).
#'
#' @return A data frame matching gprofiler output format, or NULL if no results.
#' @export
run_clusterprofiler <- function(gene_vector,
                                background = NULL,
                                ont = c("GO:BP", "GO:CC", "GO:MF"),
                                organism = "org.Mm.eg.db",
                                keyType = "SYMBOL",
                                pvalueCutoff = 0.05,
                                qvalueCutoff = 0.05,
                                minGSSize = 10,
                                maxGSSize = 500,
                                group_name = "All",
                                top_n = NULL,
                                pAdjustMethod = "BH") {

  require(clusterProfiler)
  require(organism, character.only = TRUE)

  ont <- match.arg(ont)

  gene_vector <- gene_vector[!is.null(gene_vector) & gene_vector != ""]
  background <- background[!is.null(background) & background != ""]

  if (length(gene_vector) == 0) {
    message("No valid genes provided.")
    return(NULL)
  }

  # Extract GO ontology type (BP, CC, MF)
  go_ont <- sub("GO:", "", ont)

  # Run enrichGO
  ora_result <- tryCatch({
    enrichGO(
      gene = gene_vector,
      universe = background,
      OrgDb = organism,
      keyType = keyType,
      ont = go_ont,
      pvalueCutoff = pvalueCutoff,
      qvalueCutoff = qvalueCutoff,
      minGSSize = minGSSize,
      maxGSSize = maxGSSize,
      pAdjustMethod = pAdjustMethod,
      readable = FALSE
    )
  }, error = function(e) {
    message("Error in enrichGO: ", e$message)
    return(NULL)
  })

  # Check for results
  if (is.null(ora_result) || nrow(as.data.frame(ora_result)) == 0) {
    message("No significant enrichment results for ", ont)
    return(NULL)
  }

  # Convert to gprofiler-compatible format
  result_df <- as.data.frame(ora_result)

  # Parse GeneRatio: "3/9" -> query_size = 9, intersection = 3
  gene_ratio_split <- strsplit(result_df$GeneRatio, "/")
  intersection_size <- as.numeric(sapply(gene_ratio_split, `[`, 1))
  query_size <- as.numeric(sapply(gene_ratio_split, `[`, 2))

  # Parse BgRatio: "94/28891" -> term_size = 94, effective_domain = 28891
  bg_ratio_split <- strsplit(result_df$BgRatio, "/")
  term_size <- as.numeric(sapply(bg_ratio_split, `[`, 1))
  effective_domain_size <- as.numeric(sapply(bg_ratio_split, `[`, 2))

  GO_df <- data.frame(
    term_id = result_df$ID,
    term_name = result_df$Description,
    p_value = result_df$pvalue,
    `-log10(padj)` = -log10(result_df$p.adjust),
    term_size = term_size,
    query_size = query_size,
    intersection_size = intersection_size,
    effective_domain_size = effective_domain_size,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  # Calculate Fold Enrichment
  GO_df$Fold_Enrichment <- (GO_df$intersection_size / GO_df$query_size) /
    (GO_df$term_size / GO_df$effective_domain_size)
  GO_df$`log2(Fold_Enrichment)` <- log2(GO_df$Fold_Enrichment)
  GO_df$group_name <- group_name

  # Reorder columns to match gprofiler output
  GO_df <- GO_df[, c("term_id", "term_name", "p_value", "-log10(padj)",
                     "term_size", "query_size", "intersection_size",
                     "effective_domain_size", "Fold_Enrichment",
                     "log2(Fold_Enrichment)", "group_name")]

  # Apply top_n filtering if requested
  if (!is.null(top_n) && nrow(GO_df) > top_n) {
    # Already ordered by p.adjust from enrichGO, so just take first top_n rows
    GO_df <- GO_df[1:top_n, ]
  }

  return(GO_df)
}
