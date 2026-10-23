# ==========================================================
# GSEA helpers: GO gene sets from the organism database + wrapper around run_fgsea()
# ==========================================================
if (requireNamespace("data.table", quietly = TRUE)) suppressPackageStartupMessages(library(data.table))  # run_fgsea() uses data.table syntax

# Named list GO-term -> gene SYMBOLs (all ancestors included), for ont = "BP" | "MF" | "CC"
build_go_pathways <- function(orgdb = "org.Mm.eg.db", ont = "BP") {
  suppressPackageStartupMessages(library(orgdb, character.only = TRUE))
  db <- get(orgdb)
  m <- AnnotationDbi::select(db, keys = AnnotationDbi::keys(db, keytype = "SYMBOL"), keytype = "SYMBOL",
                             columns = c("GOALL", "ONTOLOGYALL"))
  m <- unique(m[!is.na(m$GOALL) & m$ONTOLOGYALL == ont, c("SYMBOL", "GOALL")])
  terms <- AnnotationDbi::select(GO.db::GO.db, keys = unique(m$GOALL), keytype = "GOID", columns = "TERM")
  nm <- setNames(terms$TERM, terms$GOID)
  sets <- split(m$SYMBOL, m$GOALL)
  names(sets) <- make.unique(unname(nm[names(sets)]))
  sets
}

# Runs run_fgsea() for several ontologies; returns list(ontology -> list(contrast -> data.table))
run_fgsea_app <- function(fit, orgdb, onts = "BP", coefs = NULL, min_size = 10, max_size = 500,
                          nperm = 10000, padj_filter = NULL, seed = 42, save_dir = NULL) {
  out <- list()
  for (o in onts) {
    message("=== GSEA: GO:", o, " ===")
    pw <- build_go_pathways(orgdb, o)
    out[[o]] <- run_fgsea(fit, pw, coef = coefs, min_size = min_size, max_size = max_size,
                          nperm = nperm, seed = seed, padj_filter = padj_filter,
                          save_dir = if (is.null(save_dir)) NULL else file.path(save_dir, paste0("GO_", o)))
  }
  out
}
