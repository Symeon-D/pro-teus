`%||%` <- function(a, b) if (is.null(a)) b else a
# make_manifest.R -- run from the app folder whenever dependencies change:
#   source("make_manifest.R")
# Creates manifest.json, the file Posit Connect Cloud uses to install packages.
stopifnot(file.exists("app.R"), dir.exists("functions"), dir.exists("parts"))
for (p in c("BiocManager", "rsconnect", "renv", "jsonlite"))
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

# CRITICAL: make Bioconductor visible, otherwise Bioconductor packages are recorded
# wrongly / left out (this is what broke the shinyapps.io deployment).
options(repos = BiocManager::repositories())
cat("R version:", as.character(getRversion()), " Bioconductor:", as.character(BiocManager::version()), "\n")

rsconnect::writeManifest(appDir = ".", appPrimaryDoc = "app.R", verbose = TRUE)

# ---- sanity checks on the result ----
m <- jsonlite::fromJSON("manifest.json", simplifyVector = FALSE)
src <- vapply(m$packages, function(x) as.character(x$Source %||% x$description$Source %||% "?"), "")
cat("\nPackages in manifest:", length(src), "\n"); print(table(src))

need <- c("ComplexHeatmap","limma","NormalyzerDE","pcaMethods","clusterProfiler","rrvgo","fgsea",
          "GO.db","GOSemSim","org.Hs.eg.db","org.Mm.eg.db","org.Rn.eg.db","org.Ce.eg.db",
          "callr","DT","colourpicker","shinyWidgets","pdftools","gprofiler2","circlize")
miss <- setdiff(need, names(m$packages))
if (length(miss)) stop("MISSING from manifest.json: ", paste(miss, collapse = ", "),
                       "\n  -> add library() lines to deps.R, install them locally, re-run.") else
  cat("\nOK: all key packages are in manifest.json\n")

bad <- names(src)[!src %in% c("CRAN", "Bioconductor", "github", "GitHub", "RSPM", "PPM", "Posit Package Manager")]
if (length(bad)) cat("\nWARNING: unusual source (may fail on Connect Cloud):\n", paste(bad, src[bad], sep = " = ", collapse = "\n "), "\n")
cat("\nNext: git add manifest.json && git commit && git push\n")
