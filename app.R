# ==========================================================
# Proteomics Toolkit
#   1 Data & QC | 2 Normalization | 3 Imputation |
#   4 Differential analysis (+ subsets) | 5 Heatmap | 6 Pathway enrichment
# Folder layout: app.R, parts/ (UI + server pieces), functions/ (analysis functions)
# NOTE: the folder must NOT be called R/ - Shiny would auto-source it before input/output exist.
# ==========================================================
library(shiny)
library(ggplot2)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(grDevices)
library(colourpicker)
library(shinyWidgets)
library(limma)
library(matrixStats)
library(RColorBrewer)
library(readxl)
library(dplyr)
library(tidyr)
library(scales)
library(ggrepel)
library(callr)      
library(base64enc)  
library(DT)        

options(shiny.maxRequestSize = 200 * 1024^2)

# All analysis functions (yours + helpers) are loaded into the global environment
for (f in sort(list.files("functions", pattern = "\\.R$", full.names = TRUE))) source(f)
for (f in list.files("parts", pattern = "^ui_.*\\.R$", full.names = TRUE)) source(f)

# Tabs 2-6 are rendered through uiOutput(): bumping `app_gen` re-creates their UI, which puts every
# input back to its default (full reset when a new data set is loaded). Tab 1 holds the uploads, so it stays.
make_resettable <- function(tab, id) {
  body <- tab$children
  tab$children <- list(uiOutput(id))
  list(tab = tab, body = body)
}
RT <- list(body_norm    = make_resettable(tab_norm(),    "body_norm"),
           body_impute  = make_resettable(tab_impute(),  "body_impute"),
           body_de      = make_resettable(tab_de(),      "body_de"),
           body_heatmap = make_resettable(tab_heatmap(), "body_heatmap"),
           body_ora     = make_resettable(tab_ora(),     "body_ora"))

ui <- navbarPage(
  "Proteomics Toolkit", id = "nav",
  tab_data(), RT$body_norm$tab, RT$body_impute$tab, RT$body_de$tab, RT$body_heatmap$tab, RT$body_ora$tab,
  footer = div(
    style = paste("position:fixed; right:14px; bottom:14px; z-index:3000; background:#fff;",
                  "border:1px solid #ccc; border-radius:6px; padding:6px 10px;",
                  "box-shadow:0 2px 6px rgba(0,0,0,.2); max-width:420px;"),
    uiOutput("job_status", inline = TRUE),
    actionButton("stop_proc", "Stop process", icon = icon("stop"), class = "btn-danger btn-xs",
                 style = "margin-left:8px;"))
)

server <- function(input, output, session) {
  app_gen <- reactiveVal(0)
  for (id in names(RT)) local({
    i <- id
    output[[i]] <- renderUI({ app_gen(); RT[[i]]$body })
  })
  for (f in c("server_jobs.R", "server_data.R", "server_norm.R", "server_impute.R",
              "server_de.R", "server_heatmap.R", "server_ora.R", "server_gsea.R", "server_reset.R"))
    source(file.path("parts", f), local = TRUE)
}

shinyApp(ui, server)
