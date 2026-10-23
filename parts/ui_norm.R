NORM_METHODS <- c("median", "mean", "quantile", "cycloess", "VSN", "RLR", "GI")

tab_norm <- function() {
  tabPanel("2. Normalization", sidebarLayout(
    sidebarPanel(width = 3,
      checkboxInput("norm_is_log", "Data are log2-transformed", TRUE),
      helpText("They are converted back to the linear scale before NormalyzerDE (needed for VSN);",
               "all results are returned on the log2 scale."),
      hr(),
      h4("A. Evaluate methods"),
      # all methods are always evaluated: shown ticked and locked
      lapply(NORM_METHODS, function(m)
        div(class = "checkbox", tags$label(tags$input(type = "checkbox", checked = NA, disabled = NA), m))),
      actionButton("run_norm_bench", "Run benchmark", class = "btn-primary"),
      hr(),
      h4("B. Apply a method"),
      selectInput("norm_choice", "Normalization", NORM_METHODS),
      actionButton("apply_norm", "Normalize", class = "btn-success"),
      actionButton("skip_norm", "Skip normalization")
    ),
    mainPanel(width = 9,
      uiOutput("pipeline_status_norm"),
      tabsetPanel(id = "norm_tabs",
        tabPanel("Is normalization needed?", br(), uiOutput("norm_need_ui")),
        tabPanel("Benchmark", br(), uiOutput("norm_bench_ui")),
        tabPanel("Applied result", br(), textOutput("norm_applied_txt"),
                 plotOutput("norm_after_plot", height = "420px"))
      )
    )
  ))
}
