tab_data <- function() {
  tabPanel("1. Data & QC",
    tags$style(HTML("
      .data-scroll { height: calc(100vh - 300px); overflow: auto; scrollbar-width: thin; border: 1px solid #eee; }
      .data-scroll::-webkit-scrollbar { width: 8px; height: 8px; }
      .data-scroll::-webkit-scrollbar-thumb { background: #bbb; border-radius: 4px; }
      .expr-tbl th, .expr-tbl td { padding: 2px 5px !important; white-space: nowrap; }
      .qc-cell h5 { margin: 2px 0; font-weight: 600; }
      .qc-cell { margin-bottom: 4px; }")),
    uiOutput("pipeline_status_data"),
    fluidRow(
      # ---- left: uploads, group definition, filter, colours, QC options ----
      column(3, wellPanel(style = "max-height: 86vh; overflow-y: auto;",
        h4("Upload"),
        actionButton("load_demo", "Load demo data", icon = icon("flask"), class = "btn-info"),
        numericInput("global_seed", "Random seed (used by every module)", 42, width = "100%"),
        hr(),
        fileInput("expr_file", "Expression file", width = "100%",
                  accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
        uiOutput("gene_col_ui"),
        checkboxInput("log2", "log2 transform (only for raw, non-log intensities)", FALSE),
        fileInput("meta_file", "Sample info file", width = "100%",
                  accept = c(".csv", ".tsv", ".txt", ".xlsx", ".xls")),
        uiOutput("sample_col_ui"),
        h4("Group definition"),
        uiOutput("group_vars_ui"),
        h4("Filter by group missingness"),
        div(style = "display:flex; align-items:flex-end; gap:8px;",
            numericInput("miss_thr", "Min % valid", 70, min = 0, max = 100, step = 5, width = "90px"),
            div(style = "margin-bottom:15px;",
                actionButton("run_filter", "Filter data", class = "btn-primary btn-sm"),
                actionButton("reset_filter", "Reset", class = "btn-sm"))),
        textOutput("filter_msg"),
        uiOutput("group_colors_ui"),
        uiOutput("qc_opts_ui"))),
      # ---- right: 2 x 2 QC grid + data inspection ----
      column(9, tabsetPanel(id = "data_tabs", selected = "Expression data",
        tabPanel("Expression data", br(), uiOutput("data_status"),
                 helpText("Double-click a cell to edit it. Click a column header to rename the column (Enter to confirm)."),
                 DT::DTOutput("expr_dt")),
        tabPanel("Sample info (with Group)", br(),
                 helpText("Double-click a cell to fix a value; click a header to rename the column. ",
                          "Samples are matched to the expression columns by name."),
                 DT::DTOutput("meta_dt"), br(), h5("Groups"), tableOutput("group_counts")),
        tabPanel("QC plots", br(), uiOutput("qc_grid_ui"))))
    )
  )
}
