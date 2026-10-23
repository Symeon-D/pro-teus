#' Simplify GO Enrichment Results using rrvgo (With Scatter Plot)
#'
#' This function takes a dataframe of GO enrichment results, calculates semantic
#' similarity, groups redundant terms, and saves diagnostic plots.
#'
#' @param GO_df Dataframe containing g:Profiler results.
#' @param ont Ontology string (e.g., "GO:BP", "GO:MF").
#' @param organism Organism string ("mmusculus" or "hsapiens").
#' @param threshold Similarity threshold (default 0.7).
#' @param method Semantic similarity method (default "Wang").
#' @param save_plots Logical. Whether to save Scatter/Treemap plots.
#' @param output_dir Directory path to save the plots.
#' @param outfix Outfix name.
#' @param verbose Logical. Print progress messages.
#' @param show_cluster_label Logical. If TRUE (default), the most significant
#'   term label is printed at the top-left of the scatter plot cluster.
#'   Set to FALSE to suppress these labels and rely solely on the in-cluster
#'   addLabel annotations.
#'
#' @return A simplified dataframe with redundant terms removed and renamed.
#' @export
simplify_with_rrvgo <- function(GO_df,
                                ont,
                                organism,
                                threshold = 0.7,
                                method = "Wang",
                                save_plots = TRUE,
                                output_dir = ".",
                                outfix = "",
                                verbose = TRUE,
                                show_cluster_label = TRUE) {

  if (!grepl("GO", ont)) {
    if(verbose) message("   [rrvgo] Skipping: Ontology is not GO (found ", ont, ").")
    return(list(simplified_df = GO_df, full_df = GO_df))
  }

  if (nrow(GO_df) == 0) return(list(simplified_df = GO_df, full_df = GO_df))

  if (verbose) message("🔄 [rrvgo] Starting simplification for ", ont, "...")

  org_db <- NULL
  if (organism == "mmusculus") {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) stop("Install org.Mm.eg.db")
    library(org.Mm.eg.db); org_db <- "org.Mm.eg.db"
  } else if (organism == "hsapiens") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE)) stop("Install org.Hs.eg.db")
    library(org.Hs.eg.db); org_db <- "org.Hs.eg.db"
  } else if (organism == "rnorvegicus") {
    if (!requireNamespace("org.Rn.eg.db", quietly = TRUE)) stop("Install org.Rn.eg.db")
    library(org.Rn.eg.db); org_db <- "org.Rn.eg.db"
  } else if (organism == "celegans") {
    if (!requireNamespace("org.Ce.eg.db", quietly = TRUE)) stop("Install org.Ce.eg.db")
    library(org.Ce.eg.db); org_db <- "org.Ce.eg.db"
  } else {
    warning("[rrvgo] Organism not supported. Skipping.")
    return(list(simplified_df = GO_df, full_df = GO_df))
  }

  # Initialise cluster column as "Ungrouped" for all rows up front
  GO_df$cluster <- "Ungrouped"

  simplified_list <- list()
  full_list        <- list()

  for(cl in unique(GO_df$group_name)) {
    sub_df <- GO_df[GO_df$group_name == cl, ]

    if(nrow(sub_df) > 1) {

      result <- tryCatch({
        if(verbose) message("   [rrvgo] Processing cluster: ", cl, " (", nrow(sub_df), " terms)")

        # A. Calculate Similarity Matrix
        simMatrix <- rrvgo::calculateSimMatrix(
          x = sub_df$term_id,
          orgdb = org_db,
          ont = gsub("GO:", "", ont),
          method = method
        )

        # B. Split found vs missing (missing stay "Ungrouped")
        found_ids     <- rownames(simMatrix)
        missing_count <- nrow(sub_df) - length(found_ids)

        if(missing_count > 0 && verbose) {
          message("      ! ", missing_count, " terms not found in OrgDb — labeled 'Ungrouped'.")
        }

        sub_df_found   <- sub_df[sub_df$term_id %in% found_ids, ]
        sub_df_missing <- sub_df[!sub_df$term_id %in% found_ids, ]

        if(is.null(simMatrix) || nrow(sub_df_found) < 2) {
          stop("Insufficient terms remaining after OrgDb filtering.")
        }

        # C. Group Redundant Terms
        scores       <- setNames(-log10(sub_df_found$p_value), sub_df_found$term_id)
        reducedTerms <- rrvgo::reduceSimMatrix(simMatrix, scores,
                                               threshold = threshold,
                                               orgdb = org_db)

        # D. Build cluster number map: ordered by cluster size desc (ties: alphabetical),
        #    matching the treemap top-left → bottom-right visual layout
        cluster_sizes      <- sort(table(reducedTerms$parentTerm), decreasing = TRUE)
        cluster_number_map <- setNames(seq_along(cluster_sizes), names(cluster_sizes))

        # E. Generate and Save Plots
        if(save_plots) {
          clean_ont <- gsub(":", "_", ont)
          clean_cl  <- gsub("[^[:alnum:]]", "_", cl)

          scatter_file <- file.path(output_dir,
                                    paste0(clean_ont, "_rrvgo_scatter_", clean_cl, "_", outfix, ".png"))
          p <- rrvgo::scatterPlot(simMatrix, reducedTerms, size = "score", addLabel = TRUE, labelSize = 3)

          if (!show_cluster_label) {
            layer_classes <- sapply(p$layers, function(l) class(l$geom)[1])
            repel_idx <- which(layer_classes %in% c("GeomTextRepel", "GeomLabelRepel"))
            if (length(repel_idx) > 0) p$layers <- p$layers[-repel_idx]
          }

          p <- p +
            ggplot2::ggtitle(paste0("Semantic Space: ", ont, " | Cluster: ", cl)) +
            ggplot2::theme(plot.title = ggplot2::element_text(size = 14, face = "bold"))
          ggplot2::ggsave(scatter_file, plot = p, width = 8, height = 7, dpi = 300)

          tree_file <- file.path(output_dir,
                                 paste0(clean_ont, "_rrvgo_treemap_", clean_cl, "_", outfix, ".png"))
          png(tree_file, width = 2400, height = 2000, res = 300)
          rrvgo::treemapPlot(reducedTerms)
          dev.off()
        }

        # F. Annotate all found terms with their cluster number
        sub_df_found$cluster <- cluster_number_map[
          reducedTerms$parentTerm[match(sub_df_found$term_id, reducedTerms$go)]
        ]
        sub_df_found$cluster[is.na(sub_df_found$cluster)] <- "Ungrouped"

        # Full table: all terms (found + missing), cluster annotated
        full_sub_df <- rbind(sub_df_found, sub_df_missing)

        # Simplified table: original logic — one row per parent term, most significant
        sub_df_found$Parent_Term <- reducedTerms$parentTerm[
          match(sub_df_found$term_id, reducedTerms$go)
        ]
        sub_df_found <- sub_df_found[order(sub_df_found$p_value), ]
        sub_df_found$term_name <- sub_df_found$Parent_Term
        simplified_sub_df <- sub_df_found[!duplicated(sub_df_found$term_name), ]

        list(simplified = simplified_sub_df, full = full_sub_df)

      }, error = function(e) {
        if(verbose) message("      ⚠️ Cluster '", cl, "' skipped simplification: ", e$message)
        list(simplified = sub_df, full = sub_df)
      })

      simplified_list[[cl]] <- result$simplified
      full_list[[cl]]        <- result$full

    } else {
      if(verbose) message("   [rrvgo] Cluster ", cl, ": Only 1 term. No reduction needed.")
      sub_df$cluster        <- 1L
      sub_df$Parent_Term    <- sub_df$term_name
      simplified_list[[cl]] <- sub_df
      full_list[[cl]]        <- sub_df
    }
  }

  simplified_df           <- do.call(rbind, simplified_list)
  rownames(simplified_df) <- NULL
  full_df                 <- do.call(rbind, full_list)
  rownames(full_df)       <- NULL

  return(list(simplified_df = simplified_df, full_df = full_df))
}
