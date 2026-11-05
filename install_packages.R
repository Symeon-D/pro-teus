# Run once before launching / deploying the app
install.packages(c("shiny", "ggplot2", "colourpicker", "shinyWidgets", "matrixStats", "RColorBrewer",
                   "readxl", "dplyr", "tidyr", "scales", "ggrepel", "zip", "callr", "base64enc", "DT", "data.table", "pdftools", "cowplot", "reshape2",
                   "stringr", "httr", "openxlsx", "gprofiler2", "ggVennDiagram", "ggvenn", "BiocManager",
                   "imputeLCMD", "missForest", "VIM", "mice", "norm", "rrcov", "rrcovNA", "MASS"))
BiocManager::install(c("ComplexHeatmap", "limma", "NormalyzerDE", "pcaMethods", "impute",
                       "clusterProfiler", "rrvgo", "sva", "variancePartition", "fgsea", "GO.db", "GOSemSim", "org.Mm.eg.db", "org.Hs.eg.db"))
# optional extra imputers used by some methods: DreamAI, GMSimpute, SeqKnn / multiUS, pdftools
