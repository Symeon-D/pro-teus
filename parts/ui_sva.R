tab_sva <- function() {
  tabPanel("4. Batch & SVA", sidebarLayout(
    sidebarPanel(width = 3, style = "max-height: 90vh; overflow-y: auto;",
      h4("Variables"),
      helpText("Types are guessed from the sample sheet; change them if needed. 'exclude' removes a variable from every analysis on this page."),
      uiOutput("sva_types_ui"),
      hr(),
      actionButton("sva_skip", "Skip to differential analysis", icon = icon("forward"), class = "btn-default")
    ),
    mainPanel(width = 9,
      uiOutput("pipeline_status_sva"),
      uiOutput("sva_notice"),
      tabsetPanel(id = "sva_tabs",
        tabPanel("Variable relationships", br(),
                 plotOutput("sva_varvar", height = "520px"), downloadButton("dl_sva_varvar", "PNG"),
                 helpText("Binary / continuous pairs: Spearman correlation (signed, -1 to 1). Nominal pairs: bias-corrected Cramer's V;",
                          "continuous vs nominal: square root of the adjusted R-squared of a one-way ANOVA (both unsigned, 0 to 1).",
                          "Stars use uncorrected p-values (* < 0.05, ** < 0.01, *** < 0.001). Technical replicates are treated as independent",
                          "samples, so the p-values are optimistic: use this as a diagnostic.")),
        tabPanel("PCs vs covariates", br(), uiOutput("sva_pc_note"),
                 fluidRow(column(6, plotOutput("sva_pc_cor", height = "520px"), downloadButton("dl_sva_pc_cor", "PNG")),
                          column(6, plotOutput("sva_pc_imp", height = "520px"), downloadButton("dl_sva_pc_imp", "PNG"))),
                 helpText("Each principal component is regressed on the covariate: square root of the adjusted R-squared with the F-test p-value",
                          "(continuous covariates: |Spearman|). Impact = r-squared x % variance explained by the PC, shown where p < 0.05.")),
        tabPanel("Variance partitioning", br(),
                 helpText("Median share of each protein's variance explained by every covariate (variancePartition; categorical = random effects, continuous = fixed)."),
                 actionButton("sva_run_vp", "Run variance partitioning", class = "btn-primary"), br(), br(),
                 uiOutput("sva_vp_ui")),
        tabPanel("Surrogate variables", br(),
                 fluidRow(
                   column(4, uiOutput("sva_protect_ui")),
                   column(4, radioButtons("sva_method", "Estimator for the number of SVs", c("leek", "be"), selected = "leek", inline = TRUE),
                             numericInput("sva_nsv", "Number of SVs (blank = use the estimate)", NA, min = 0, step = 1)),
                   column(4, br(), uiOutput("sva_run_ui"))),
                 uiOutput("sva_result_ui")),
        tabPanel("Post-correction check", br(),
                 helpText("Removes the carried SVs (and any known covariates ticked below) from the data and shows how the PCA and the variance partitioning change.",
                          "In the differential analysis the same covariates must be in the design; the heatmap and protein plots then use the corrected data automatically."),
                 fluidRow(column(5, uiOutput("sva_remove_ui")),
                          column(4, uiOutput("sva_pca_ui")),
                          column(3, br(), actionButton("sva_post_run", "Run post-correction check", class = "btn-primary"))),
                 uiOutput("sva_post_ui"))
      )
    )
  ))
}
