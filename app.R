# ==========================================================
# Pro-teus
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
library(callr)      # background jobs (Stop button)
library(base64enc)  # 'Open in ORA-cle' links
library(DT)         # editable data tables

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

# ---- look & feel: dark mode (button top right), GitHub link, public-server notice ----
PROTEUS_CSS <- "
  .dm-sun { display: none; } body.dark .dm-sun { display: inline; } body.dark .dm-moon { display: none; }
  #proteus-tools > li > a { padding-top: 15px; padding-bottom: 15px; font-size: 16px; }
  body.dark { background: #161a1f; color: #e4e6ea; }
  body.dark .navbar-default { background: #1f2430; border-color: #2b3140; }
  body.dark .navbar-default .navbar-brand, body.dark .navbar-default .navbar-nav > li > a { color: #d6d9df; }
  body.dark .navbar-default .navbar-nav > li > a:hover { background: #2b3140; color: #fff; }
  body.dark .navbar-default .navbar-nav > .active > a, body.dark .navbar-default .navbar-nav > .active > a:hover,
  body.dark .navbar-default .navbar-nav > .active > a:focus { background: #2b3140; color: #fff; }
  body.dark .well, body.dark .panel { background: #1f2430; border-color: #2b3140; color: #e4e6ea; }
  body.dark .form-control, body.dark .selectize-input, body.dark .selectize-dropdown, body.dark select,
  body.dark .selectize-input input { background: #0f1217; color: #e4e6ea; border-color: #39404f; }
  body.dark .selectize-dropdown .active { background: #2b3140; color: #fff; }
  body.dark .selectize-input.full { background: #0f1217; }
  body.dark .nav-tabs { border-color: #2b3140; }
  body.dark .nav-tabs > li > a:hover { background: #2b3140; border-color: #2b3140; }
  body.dark .nav-tabs > li.active > a, body.dark .nav-tabs > li.active > a:hover,
  body.dark .nav-tabs > li.active > a:focus { background: #1f2430; color: #fff; border-color: #2b3140 #2b3140 transparent; }
  body.dark .table { color: #e4e6ea; }
  body.dark .table-striped > tbody > tr:nth-of-type(odd) { background: #1b2029; }
  body.dark .table-hover > tbody > tr:hover { background: #2b3140; }
  body.dark .table > thead > tr > th, body.dark .table > tbody > tr > td { border-color: #2b3140; }
  body.dark .btn-default { background: #2b3140; color: #e4e6ea; border-color: #39404f; }
  body.dark .btn-link { color: #7db7ff; }
  body.dark .dropdown-menu, body.dark .modal-content { background: #1f2430; border-color: #39404f; color: #e4e6ea; }
  body.dark .dropdown-menu > li > a { color: #e4e6ea; }
  body.dark code, body.dark pre { background: #0f1217; color: #e4e6ea; border-color: #39404f; }
  body.dark .help-block, body.dark .text-muted, body.dark small { color: #a9afbb; }
  body.dark a { color: #7db7ff; } body.dark .navbar a, body.dark .btn { color: inherit; }
  body.dark hr { border-color: #2b3140; }
  body.dark .alert { background: #25303d; border-color: #39404f; color: #e4e6ea; }
  body.dark .shiny-notification { background: #1f2430; color: #e4e6ea; border-color: #39404f; }
  body.dark table.dataTable, body.dark .dataTables_wrapper, body.dark .dataTables_wrapper label { color: #e4e6ea; }
  body.dark table.dataTable tbody tr { background: #1f2430; }
  body.dark table.dataTable.stripe tbody tr.odd, body.dark table.dataTable.stripe tbody tr.odd > td { background: #1b2029; }
  body.dark table.dataTable thead th, body.dark table.dataTable thead td { border-color: #39404f; }
  body.dark .dataTables_wrapper .dataTables_paginate .paginate_button { color: #e4e6ea !important; }
  body.dark .shiny-plot-output img, body.dark .shiny-image-output img { background: #fff; border-radius: 4px; }
  body.dark .shiny-plot-output { background: #fff; border-radius: 4px; }
"
PROTEUS_JS <- sprintf("
  (function() {
    $(function() {
      if (localStorage.getItem('proteus_dark') === '1') $('body').addClass('dark');
      var tools = '<ul class=\"nav navbar-nav navbar-right\" id=\"proteus-tools\">' +
        '<li><a href=\"#\" id=\"dark_toggle\" title=\"Dark / light mode\"><span class=\"dm-moon\">' + %s +
        '</span><span class=\"dm-sun\">' + %s + '</span></a></li>' +
        '<li><a href=\"https://github.com/Symeon-D/pro-teus\" target=\"_blank\" rel=\"noopener\" title=\"Pro-teus on GitHub\">' + %s + '</a></li></ul>';
      $('.navbar-nav').first().after(tools);
      $(document).on('click', '#dark_toggle', function(e) {
        e.preventDefault();
        var on = $('body').toggleClass('dark').hasClass('dark');
        localStorage.setItem('proteus_dark', on ? '1' : '0');
      });
    });
  })();",
  jsonlite::toJSON(as.character(icon("moon")), auto_unbox = TRUE),
  jsonlite::toJSON(as.character(icon("sun")), auto_unbox = TRUE),
  jsonlite::toJSON(as.character(icon("github")), auto_unbox = TRUE))

ui <- navbarPage(
  "Pro-teus", id = "nav", windowTitle = "Pro-teus",
  header = tagList(
    tags$head(tags$style(HTML(PROTEUS_CSS)), tags$script(HTML(PROTEUS_JS))),
    div(style = "background:#fff7e0;border-bottom:1px solid #e6c200;padding:6px 14px;font-size:13px;color:#333;",
        "Public server: uploaded files are processed on third-party cloud infrastructure and deleted when your ",
        "session ends. Do not upload sensitive or unpublished data you are not comfortable sharing; use the demo ",
        "dataset or run the app locally instead.")),
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
