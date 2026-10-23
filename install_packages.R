# Run once on your own computer before launching / deploying the app.
# (Connect Cloud does NOT run this file; it installs packages from manifest.json.)
options(repos = c(CRAN = "https://cloud.r-project.org"))

install.packages(c("shiny", "ggplot2", "colourpicker", "shinyWidgets", "matrixStats", "RColorBrewer",
                   "readxl", "dplyr", "tidyr", "scales", "ggrepel", "zip", "callr", "base64enc", "DT",
                   "data.table", "pdftools", "cowplot", "reshape2", "stringr", "httr", "openxlsx",
                   "gprofiler2", "ggVennDiagram", "ggvenn", "BiocManager",
                   "imputeLCMD", "missForest", "VIM", "mice", "norm", "rrcov", "rrcovNA", "MASS",
                   "circlize", "jsonlite", "rsconnect", "renv"))

BiocManager::install(c("ComplexHeatmap", "limma", "NormalyzerDE", "pcaMethods", "impute",
                       "clusterProfiler", "rrvgo", "fgsea", "GO.db", "GOSemSim",
                       "org.Mm.eg.db", "org.Hs.eg.db", "org.Rn.eg.db", "org.Ce.eg.db"),
                     ask = FALSE, update = FALSE)

# Optional extra imputers used by some methods: DreamAI, GMSimpute, SeqKnn / multiUS.
# Only add them here AND in deps.R if your code really needs them on the server.
