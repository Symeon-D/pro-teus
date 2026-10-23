#' Run ORA Analysis Pipeline
#'
#' This function performs Over-Representation Analysis (ORA) using selected method(s),
#' generates visualizations, and optionally applies GO term simplification via rrvgo.
#'
#' @param subsets_list A named list where each element is a character vector of gene symbols/IDs.
#' @param subset_name Character. Name for the subset when `subsets_list` is a single gene vector.
#'   If NULL and a single vector is provided, will attempt to use the variable name.
#' @param background A character vector of background genes (gene universe). Defaults to NULL.
#' @param top_n Number of top terms to show.
#' @param output_dir Character string specifying the root directory for output files.
#' @param min_genes Numeric. Minimum number of genes per functional term. Default is 5.
#' @param min_term_size Numeric. Minimm allowed size of a GO term. Default is 10.
#' @param max_term_size Numeric. Maximum allowed size of a GO term. Default is 500.
#' @param correction_method Character vector. Multiple testing correction method(s).
#'   Options include "gSCS", "fdr", "bonferroni". Default is "gSCS".
#'   Can specify multiple methods (e.g., c("gSCS", "fdr")) to run both.
#'   Note: Only applies to gprofiler; other tools use their default corrections.
#' @param ontologies A character vector of ontologies to test. Default is "GO:BP".
#' @param ora_method Character vector. ORA method(s) to use.
#'   Options: "gprofiler", "genecodis", "clusterprofiler". Default is "gprofiler".
#'   Can specify multiple methods (e.g., c("gprofiler", "clusterprofiler")).
#' @param simplify Logical. Whether to run rrvgo simplification. Default is FALSE.
#' @param cp_simplify Logical. Whether to run clusterProfiler's native simplify(). Default is FALSE.
#' @param cp_cutoff Numeric. Similarity threshold for clusterProfiler's simplify(). Default is 0.7.
#' @param organism Character. Organism for clusterProfiler and rrvgo.
#'   Default is "org.Mm.eg.db" for clusterProfiler, "mmusculus" for rrvgo.
#' @param keyType Character. Gene ID type for clusterProfiler. Default is "SYMBOL".
#' @param rrvgo_threshold Numeric. Similarity threshold for rrvgo. Default is 0.7.
#' @param save_df Logical. If true it saves the enrichment result in a csv file.
#' @param show_cluster_label Logical. If TRUE (default), the most significant term label
#'   is printed at the top-left of each rrvgo scatter plot cluster. Set to FALSE to
#'   suppress these labels.
#'
#' @return Invisibly returns NULL. Outputs are saved to disk.
#'   Output folder structure:
#'     output_dir/
#'       <subset_name>/
#'         gprofiler/
#'           <gSCS>/   (one sub-folder per correction method)
#'         genecodis/
#'         clusterprofiler/
#'
#' @export
#'
#' @examples
#' \dontrun{
#' gene_list <- list(Cluster1 = c("Eno1", "Gapdh", "Pkm"), Cluster2 = c("Actb", "Myh9"))
#'
#' # Run only gprofiler without simplification
#' run_ora_analysis_pipeline(
#'   subsets_list = gene_list,
#'   output_dir = "ORA_Results",
#'   ora_method = "gprofiler",
#'   simplify = FALSE
#' )
#'
#' # Run all methods with multiple correction methods
#' run_ora_analysis_pipeline(
#'   subsets_list = gene_list,
#'   output_dir = "ORA_Results",
#'   ora_method = c("gprofiler", "genecodis", "clusterprofiler"),
#'   correction_method = c("gSCS", "fdr"),
#'   simplify = TRUE,
#'   show_cluster_label = FALSE
#' )
#' }
run_ora_analysis_pipeline <- function(subsets_list,
                                      subset_name = NULL,
                                      background = NULL,
                                      top_n = NULL,
                                      output_dir,
                                      min_genes = 5,
                                      min_term_size = 10,
                                      max_term_size = 500,
                                      correction_method = c("gSCS"),
                                      ontologies = c("GO:BP"),
                                      ora_method = c("gprofiler"),
                                      simplify = FALSE,
                                      cp_simplify = FALSE,
                                      cp_cutoff = 0.7,
                                      organism = "org.Mm.eg.db",
                                      keyType = "SYMBOL",
                                      rrvgo_threshold = 0.7,
                                      save_df = TRUE,
                                      show_cluster_label = TRUE,
                                      cp_correction_method = "BH") {

  # =========================================================================
  # Helper: create dir if missing
  # =========================================================================
  make_dir <- function(...) {
    p <- file.path(...)
    if (!dir.exists(p)) dir.create(p, recursive = TRUE)
    p
  }

  # =========================================================================
  # Helper: safe plot wrapper
  # Wraps GO_plot_wrapper in tryCatch so that rendering failures (e.g. the
  # ragg "dimensions exceed maximum" error when a term list is very large)
  # do not abort the pipeline. The enrichment CSV is always written first so
  # the data is never lost even when the plot fails.
  # =========================================================================
  safe_plot <- function(label, ...) {
    tryCatch(
      GO_plot_wrapper(...),
      error = function(e) {
        message(paste0("      [PLOT SKIPPED] ", label, ": ", conditionMessage(e)))
      }
    )
  }

  # =========================================================================
  # Handle input: convert single gene vector to named list
  # =========================================================================
  if (is.character(subsets_list) && is.null(dim(subsets_list))) {
    if (is.null(subset_name)) {
      subset_name <- deparse(substitute(subsets_list))
      if (grepl("\\$", subset_name)) subset_name <- sub(".*\\$", "", subset_name)
      if (grepl("\\[\\[", subset_name)) subset_name <- gsub('.*\\[\\["?|"?\\]\\].*', "", subset_name)
    }
    subsets_list <- setNames(list(subsets_list), subset_name)
    message(paste("Single gene vector converted to list with name:", subset_name))
  }

  # Validate ora_method
  valid_methods <- c("gprofiler", "genecodis", "clusterprofiler")
  ora_method <- match.arg(ora_method, valid_methods, several.ok = TRUE)

  # Validate correction_method
  valid_corrections <- c("gSCS", "fdr", "bonferroni")
  correction_method <- match.arg(correction_method, valid_corrections, several.ok = TRUE)
  cp_correction_method <- match.arg(cp_correction_method,
    c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY"), several.ok = TRUE)

  # Organism format conversion
  rrvgo_organism <- switch(organism,
                           "org.Mm.eg.db" = "mmusculus",
                           "org.Hs.eg.db" = "hsapiens",
                           "org.Rn.eg.db" = "rnorvegicus",
                           "org.Ce.eg.db" = "celegans",
                           organism)

  gc_organism <- switch(organism,
                        "org.Mm.eg.db" = "Mus musculus",
                        "org.Hs.eg.db" = "Homo sapiens",
                        "org.Rn.eg.db" = "Rattus norvegicus",
                        "org.Ce.eg.db" = "Caenorhabditis elegans",
                        organism)

  make_dir(output_dir)

  outfix_base <- paste0("min_", min_genes, "_TS_",
                        if (is.null(max_term_size)) "NULL" else max_term_size)

  for (current_subset_name in names(subsets_list)) {
    message(paste("\n=== Processing subset:", current_subset_name, "==="))

    genes <- subsets_list[[current_subset_name]]

    if (length(genes) < min_genes) {
      message(paste("  Skipping - only", length(genes), "genes (min:", min_genes, ")"))
      next
    }

    # Root folder for this subset
    subset_path <- make_dir(output_dir, current_subset_name)

    for (select_ontology in ontologies) {
      message(paste("  Ontology:", select_ontology))

      # Clean ontology string for use in filenames (e.g. "GO:BP" -> "GO_BP")
      clean_ont <- gsub(":", "_", select_ontology)

      # =====================================================================
      # GPROFILER — one sub-folder per ontology
      # =====================================================================
      if ("gprofiler" %in% ora_method) {
        for (corr_method in correction_method) {
          message(paste("    Running: gprofiler with", corr_method, "correction"))

          tool_path <- make_dir(subset_path, "gprofiler", clean_ont)
          outfix <- paste0(clean_ont, "_", outfix_base, "_", corr_method)

          GO_df <- run_gprofiler(
            organism = rrvgo_organism,   # PATCH: organism was never forwarded (defaulted to mmusculus)
            top_n = top_n,
            PART_results = NULL,
            gene_vector = genes,
            background = background,
            ont = select_ontology,
            correction_method = corr_method,
            min_term_size = min_term_size,
            max_term_size = max_term_size,
            min_genes = min_genes,
            min_fold_enrichment = NULL,
            all_as_one_cluster = FALSE
          )

          if (!is.null(GO_df) && nrow(GO_df) > 0) {


            # Save CSV first — data is never lost even if plot fails
            if (save_df)
              write.csv(GO_df,
                        file = file.path(tool_path, paste0(current_subset_name, "_GO_df_gprofiler_", outfix, ".csv")),
                        row.names = FALSE)

            safe_plot(
              label = paste0("gprofiler / ", corr_method, " / ", select_ontology),
              GO_df = GO_df, ont = select_ontology, save_dir = tool_path,
              outfix = outfix, top_n = top_n, lolipop = TRUE, multi_cluster = FALSE,
              term_size_legend = "percentage", custom_width = 18, ball_size_boost = 1,
              color_low = "blue", color_high = "red", lolipop_linewidth = 2
            )

            if (simplify) {
              message(paste("      Simplifying: gprofiler (", corr_method, ")"))
              rrvgo_result <- simplify_with_rrvgo(
                GO_df = GO_df, ont = select_ontology, organism = rrvgo_organism,
                threshold = rrvgo_threshold, save_plots = TRUE,
                output_dir = tool_path, outfix = outfix, verbose = TRUE,
                show_cluster_label = show_cluster_label
              )

              rrvgo_simplified_df <- rrvgo_result$simplified_df
              if (!is.null(rrvgo_simplified_df) && nrow(rrvgo_simplified_df) > 0) {
                safe_plot(
                  label = paste0("gprofiler rrvgo / ", corr_method, " / ", select_ontology),
                  GO_df = rrvgo_simplified_df, ont = select_ontology, save_dir = tool_path,
                  outfix = paste0("rrvgo_", outfix), top_n = top_n, lolipop = TRUE,
                  multi_cluster = FALSE, term_size_legend = "percentage", custom_width = 18,
                  ball_size_boost = 1, color_low = "blue", color_high = "red", lolipop_linewidth = 2
                )
              }

              rrvgo_full_df <- rrvgo_result$full_df
              if (save_df && !is.null(rrvgo_full_df) && nrow(rrvgo_full_df) > 0) {
                write.csv(rrvgo_full_df,
                          file = file.path(tool_path, paste0("GO_df_gprofiler_rrvgo_full_", outfix, ".csv")),
                          row.names = FALSE)
              }
            }
          } else {
            message(paste("      gprofiler (", corr_method, "): No significant terms"))
          }
        }
      }

      # =====================================================================
      # GENECODIS
      # =====================================================================
      if ("genecodis" %in% ora_method) {
        message("    Running: genecodis")

        tool_path <- make_dir(subset_path, "genecodis", clean_ont)
        outfix <- paste0(clean_ont, "_", outfix_base, "_GeneCodis")

        GO_df <- run_Gene_Codis(
          gene_vector = genes,
          top_n = top_n,
          background = background,
          organism = gc_organism,
          ont = select_ontology,
          min_genes = min_genes,
          max_term_size = max_term_size
        )

        if (!is.null(GO_df) && nrow(GO_df) > 0) {


          if (save_df)
            write.csv(GO_df,
                      file = file.path(tool_path, paste0("GO_df_genecodis_", outfix, ".csv")),
                      row.names = FALSE)

          safe_plot(
            label = paste0("genecodis / ", select_ontology),
            GO_df = GO_df, ont = select_ontology, save_dir = tool_path,
            outfix = outfix, top_n = top_n, lolipop = TRUE, multi_cluster = FALSE,
            term_size_legend = "percentage", custom_width = 18, ball_size_boost = 1,
            color_low = "blue", color_high = "red", lolipop_linewidth = 2
          )

          if (simplify) {
            message("      Simplifying: genecodis")
            rrvgo_result <- simplify_with_rrvgo(
              GO_df = GO_df, ont = select_ontology, organism = rrvgo_organism,
              threshold = rrvgo_threshold, save_plots = TRUE,
              output_dir = tool_path, outfix = outfix, verbose = TRUE,
              show_cluster_label = show_cluster_label
            )

            rrvgo_simplified_df <- rrvgo_result$simplified_df
            if (!is.null(rrvgo_simplified_df) && nrow(rrvgo_simplified_df) > 0) {
              safe_plot(
                label = paste0("genecodis rrvgo / ", select_ontology),
                GO_df = rrvgo_simplified_df, ont = select_ontology, save_dir = tool_path,
                outfix = paste0("rrvgo_", outfix), top_n = top_n, lolipop = TRUE,
                multi_cluster = FALSE, term_size_legend = "percentage", custom_width = 18,
                ball_size_boost = 1, color_low = "blue", color_high = "red", lolipop_linewidth = 2
              )
            }

            rrvgo_full_df <- rrvgo_result$full_df
            if (save_df && !is.null(rrvgo_full_df) && nrow(rrvgo_full_df) > 0) {
              write.csv(rrvgo_full_df,
                        file = file.path(tool_path, paste0("GO_df_genecodis_rrvgo_full_", outfix, ".csv")),
                        row.names = FALSE)
            }
          }
        } else {
          message("      genecodis: No significant terms")
        }
      }

      # =====================================================================
      # CLUSTERPROFILER
      # =====================================================================
      if ("clusterprofiler" %in% ora_method) for (cp_corr in cp_correction_method) {
        message("    Running: clusterprofiler with ", cp_corr, " correction")

        tool_path <- make_dir(subset_path, "clusterprofiler", clean_ont)
        outfix <- paste0(clean_ont, "_", outfix_base, "_ClusterProfiler_", cp_corr)

        GO_df <- run_clusterprofiler(
          gene_vector = genes,
          background = background,
          ont = select_ontology,
          organism = organism,
          keyType = keyType,
          minGSSize = if (is.null(min_term_size)) 1 else min_term_size,   # NULL = no limit
          maxGSSize = if (is.null(max_term_size)) 100000 else max_term_size,
          top_n = top_n,
          pAdjustMethod = cp_corr
        )

        if (!is.null(GO_df) && nrow(GO_df) > 0) {


          if (save_df)
            write.csv(GO_df,
                      file = file.path(tool_path, paste0("GO_df_clusterprofiler_", outfix, ".csv")),
                      row.names = FALSE)

          safe_plot(
            label = paste0("clusterprofiler / ", select_ontology),
            GO_df = GO_df, ont = select_ontology, save_dir = tool_path,
            outfix = outfix, top_n = top_n, lolipop = TRUE, multi_cluster = FALSE,
            term_size_legend = "percentage", custom_width = 18, ball_size_boost = 1,
            color_low = "blue", color_high = "red", lolipop_linewidth = 2
          )

          if (simplify) {
            message("      Simplifying: clusterprofiler")
            rrvgo_result <- simplify_with_rrvgo(
              GO_df = GO_df, ont = select_ontology, organism = rrvgo_organism,
              threshold = rrvgo_threshold, save_plots = TRUE,
              output_dir = tool_path, outfix = outfix, verbose = TRUE,
              show_cluster_label = show_cluster_label
            )

            rrvgo_simplified_df <- rrvgo_result$simplified_df
            if (!is.null(rrvgo_simplified_df) && nrow(rrvgo_simplified_df) > 0) {
              safe_plot(
                label = paste0("clusterprofiler rrvgo / ", select_ontology),
                GO_df = rrvgo_simplified_df, ont = select_ontology, save_dir = tool_path,
                outfix = paste0("rrvgo_", outfix), top_n = top_n, lolipop = TRUE,
                multi_cluster = FALSE, term_size_legend = "percentage", custom_width = 18,
                ball_size_boost = 1, color_low = "blue", color_high = "red", lolipop_linewidth = 2
              )
            }

            rrvgo_full_df <- rrvgo_result$full_df
            if (save_df && !is.null(rrvgo_full_df) && nrow(rrvgo_full_df) > 0) {
              write.csv(rrvgo_full_df,
                        file = file.path(tool_path, paste0("GO_df_clusterprofiler_rrvgo_full_", outfix, ".csv")),
                        row.names = FALSE)
            }
          }
        } else {
          message("      clusterprofiler: No significant terms")
        }
      }
    }
  }

  # clusterProfiler::simplify() is applied to the ORIGINAL result tables (never to rrvgo output);
  # the original and the simplified files are both kept
  if (cp_simplify) {
    message("\nApplying clusterProfiler::simplify() to the original result tables...")
    files <- list.files(output_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
    files <- files[!grepl("rrvgo|cpSimplified", basename(files))]
    for (f in files) {
      ont_f <- gsub("_", ":", basename(dirname(f)))            # GO_BP -> GO:BP
      df_f <- tryCatch(read.csv(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
      if (is.null(df_f) || nrow(df_f) < 2) next
      sres <- cp_simplify_table(df_f, ont_f, cp_cutoff, paste0(" (", basename(f), ")"), orgdb = organism)
      if (isTRUE(sres$ok) && nrow(sres$df) > 0)
        write.csv(sres$df, sub("\\.csv$", "_cpSimplified.csv", f), row.names = FALSE)
    }
  }

  message("\n=== Pipeline complete! Results saved in: ", output_dir, " ===")
  invisible(NULL)
}

#' Coerce an ORA data frame into a clusterProfiler enrichResult
#'
#' This function takes a flat data frame of enrichment results generated by tools
#' like \code{run_gprofiler} or \code{run_Gene_Codis} and coerces it into an S4
#' \code{enrichResult} object. This allows \code{clusterProfiler::simplify()}
#' to operate seamlessly on non-clusterProfiler results.
#'
#' @param df A data frame containing ORA results. The function automatically maps
#'   columns from \code{run_gprofiler} (e.g., \code{term_id}, \code{term_name},
#'   \code{p_value}) to clusterProfiler standards (\code{ID}, \code{Description},
#'   \code{p.adjust}, \code{Count}, \code{geneID}).
#' @param ont_string Character string representing the ontology (e.g., "GO:BP", "GO:CC").
#'   The "GO:" prefix is stripped internally.
#'
#' @return An S4 object of class \code{enrichResult}.
#'
#' @importFrom methods new
#' @export
df_to_enrichResult <- function(df, ont_string) {
  # Clean ontology string (e.g., "GO:BP" -> "BP")
  ont_clean <- gsub("GO:", "", ont_string)

  # 1. Map gprofiler column names to clusterProfiler standard names
  if ("term_id" %in% colnames(df) && !"ID" %in% colnames(df)) {
    df$ID <- df$term_id
  }
  if ("term_name" %in% colnames(df) && !"Description" %in% colnames(df)) {
    df$Description <- df$term_name
  }
  if ("p_value" %in% colnames(df) && !"p.adjust" %in% colnames(df)) {
    df$p.adjust <- df$p_value
    df$pvalue   <- df$p_value  # clusterProfiler sometimes calls both
    df$qvalue   <- df$p_value  # Added for safety during Wang semantic calculations
  }
  if ("intersection_size" %in% colnames(df) && !"Count" %in% colnames(df)) {
    df$Count <- as.numeric(df$intersection_size)
  }

  # 2. Format intersection genes (clusterProfiler uses '/' instead of ',')
  if ("intersection" %in% colnames(df) && !"geneID" %in% colnames(df)) {
    df$geneID <- gsub(",", "/", df$intersection)
  }

  # 3. Format Ratios (optional, but ensures plotting compatibility downstream)
  if ("intersection_size" %in% colnames(df) && "query_size" %in% colnames(df)) {
    df$GeneRatio <- paste0(df$intersection_size, "/", df$query_size)
  }
  if ("term_size" %in% colnames(df) && "effective_domain_size" %in% colnames(df)) {
    df$BgRatio <- paste0(df$term_size, "/", df$effective_domain_size)
  }

  # 4. Initialize the S4 enrichResult object
  suppressWarnings({
    mock_obj <- methods::new("enrichResult",
                             result = as.data.frame(df),
                             ontology = ont_clean,
                             pvalueCutoff = 1,
                             pAdjustMethod = "NONE",
                             organism = "UNKNOWN",
                             keytype = "UNKNOWN",
                             gene = character(0))
  })

  return(mock_obj)
}


#' clusterProfiler::simplify() for any ORA table (gprofiler / genecodis / clusterprofiler).
#' The table is coerced to an enrichResult, simplified, and the ORIGINAL rows of the retained
#' terms are returned, so the column structure used by the plotting functions is preserved.
cp_simplify_table <- function(GO_df, ont, cutoff, label = "", orgdb = NULL) {
  ont_short <- gsub("GO:", "", ont)
  id_col <- if ("term_id" %in% colnames(GO_df)) "term_id" else if ("ID" %in% colnames(GO_df)) "ID" else NULL
  p_col  <- if ("p_value" %in% colnames(GO_df)) "p_value" else if ("p.adjust" %in% colnames(GO_df)) "p.adjust" else NULL
  if (is.null(id_col) || is.null(p_col) || nrow(GO_df) < 2) {
    message("      [!] simplify skipped", label, ": needs term ids, p-values and at least 2 terms.")
    return(list(df = GO_df, ok = FALSE))
  }
  res <- tryCatch({
    # same idea as clusterProfiler::simplify(): Wang similarity between GO terms; of two terms with
    # similarity > cutoff only the more significant one is kept.
    semData <- GOSemSim::godata(OrgDb = orgdb %||% "org.Mm.eg.db", ont = ont_short, computeIC = FALSE)
    ids <- unique(GO_df[[id_col]])
    sim <- GOSemSim::mgoSim(ids, ids, semData = semData, measure = "Wang", combine = NULL)
    sim[is.na(sim)] <- 0
    pv <- GO_df[[p_col]][match(ids, GO_df[[id_col]])]
    ord <- order(pv); removed <- rep(FALSE, length(ids)); names(removed) <- ids
    for (k in seq_along(ord)) {
      i <- ord[k]; if (removed[i]) next
      later <- ord[-seq_len(k)]
      sim_i <- sim[ids[i], ids[later]]
      removed[later[sim_i > cutoff]] <- TRUE
    }
    keep <- names(removed)[!removed]
    list(df = GO_df[GO_df[[id_col]] %in% keep, , drop = FALSE], ok = TRUE)
  }, error = function(e) {
    message("      [!] simplify failed", label, ": ", conditionMessage(e))
    NULL
  })
  if (!is.null(res)) return(res)
  # fallback: clusterProfiler's own implementation
  tryCatch({
    fake_er <- df_to_enrichResult(GO_df, ont)
    simplified_er <- clusterProfiler::simplify(fake_er, cutoff = cutoff, by = "p.adjust", select_fun = min)
    list(df = GO_df[GO_df[[id_col]] %in% simplified_er@result$ID, , drop = FALSE], ok = TRUE)
  }, error = function(e) {
    message("      [!] clusterProfiler::simplify() also failed", label, ". Using unsimplified data. Error: ", conditionMessage(e))
    list(df = GO_df, ok = FALSE)
  })
}
