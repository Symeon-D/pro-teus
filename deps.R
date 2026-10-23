# deps.R -- NEVER executed (guarded by if (FALSE)).
# Its only job: let rsconnect::writeManifest() see every package the app needs,
# including ones loaded indirectly (e.g. org.*.eg.db chosen by species at run time).
# Base packages (grid, grDevices, stats, utils, tools, parallel) are not listed.
if (FALSE) {
  # --- from app.R ---
  library(shiny); library(ggplot2); library(ComplexHeatmap); library(circlize)
  library(colourpicker); library(shinyWidgets); library(limma); library(matrixStats)
  library(RColorBrewer); library(readxl); library(dplyr); library(tidyr); library(scales)
  library(ggrepel); library(callr); library(base64enc); library(DT)

  # --- from install_packages.R ---
  library(zip); library(data.table); library(pdftools); library(cowplot); library(reshape2)
  library(stringr); library(httr); library(openxlsx); library(gprofiler2)
  library(ggVennDiagram); library(ggvenn); library(jsonlite)
  library(imputeLCMD); library(missForest); library(VIM); library(mice); library(norm)
  library(rrcov); library(rrcovNA); library(MASS)
  library(NormalyzerDE); library(pcaMethods); library(impute)
  library(clusterProfiler); library(rrvgo); library(fgsea); library(GO.db); library(GOSemSim)

  # --- species databases (loaded dynamically by organism choice) ---
  library(org.Hs.eg.db); library(org.Mm.eg.db); library(org.Rn.eg.db); library(org.Ce.eg.db)

  # --- ADD HERE anything scan_deps.R reports as missing, e.g.: ---
  library(AnnotationDbi); library(curl); library(rlang); library(GMSimpute); 
 library(multiUS); library(DreamAI)
}
