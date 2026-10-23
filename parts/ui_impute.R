IMP_CHOICES <- c("MinDet" = "mindet", "MinProb" = "minprob", "QRILC" = "qr", "Perseus" = "perseus",
                 "BPCA" = "bpca", "SVD" = "svd", "KNN" = "knn", "LLS" = "lls", "Random forest" = "rf",
                 "Zero" = "zero", "Minimum" = "minimum", "Column median" = "colmedian",
                 "Row median" = "rowmedian", "MLE" = "mle", "SeqKNN" = "seqknn", "IRMI" = "irmi",
                 "MICE (norm)" = "mice_norm", "MICE (cart)" = "mice_cart",
                 "ImpSeq" = "impseq", "ImpSeqRob" = "impseqrob")

tab_impute <- function() {
  tabPanel("3. Imputation", sidebarLayout(
    sidebarPanel(width = 4,
      h4("A. Evaluate methods"),
      imp_eval_methods_ui(),
      checkboxInput("imp_run_de", "Also compare differential-expression concordance (slower)", FALSE),
      actionButton("run_imp_eval", "Evaluate", class = "btn-primary"),
      downloadButton("dl_imp_pdf", "PDF"),
      hr(),
      h4("B. Apply a method"),
      selectInput("imp_final", "Imputation method", imp_grouped_choices(), selected = "perseus"),
      uiOutput("imp_apply_opts"),
      actionButton("apply_imp", "Impute", class = "btn-success"),
      actionButton("skip_imp", "Skip imputation")
    ),
    mainPanel(
      uiOutput("pipeline_status_imp"),
      textOutput("imp_na_txt"),
      tabsetPanel(id = "imp_tabs",
        tabPanel("Evaluation", br(),
          tabsetPanel(id = "imp_eval_tabs",
            # first sub-tab is always filled; the others stay empty until 'Evaluate' is pressed
            tabPanel("Missingness pattern", br(), plotOutput("imp_miss_plot", height = "950px")),
            tabPanel("Distributions (before / after)", br(), plotOutput("imp_density", height = "900px")),
            tabPanel("PCA", br(), plotOutput("imp_pca", height = "900px")),
            tabPanel("NRMSE", br(), plotOutput("imp_nrmse", height = "520px")),
            tabPanel("CV", br(), plotOutput("imp_cv", height = "520px")),
            tabPanel("Downshift", br(), plotOutput("imp_downshift", height = "520px")),
            tabPanel("Classic metrics", br(), plotOutput("imp_classic", height = "900px")),
            tabPanel("DE concordance", br(), plotOutput("imp_de_concordance", height = "600px")))),
        tabPanel("Applied result", br(), textOutput("imp_applied_txt"),
                 plotOutput("imp_after_plot", height = "420px"))
      )
    )
  ))
}

# parameters of my_impute() that matter per method: list(name = c(label, default, step))
IMP_PARAMS <- list(
  perseus = list(downshift = c("Downshift (sd)", 1.8, 0.1), new_width_cutoff = c("Width factor", 0.3, 0.05),
                 sigma_cutoff = c("Sigma cutoff", 6, 1)),
  mindet  = list(q = c("Quantile q", 0.01, 0.01)),
  minprob = list(q = c("Quantile q", 0.01, 0.01), tune.sigma = c("tune.sigma", 1, 0.1)),
  qr      = list(tune.sigma = c("tune.sigma", 1, 0.1)),
  knn     = list(k = c("k (neighbours)", 10, 1)), seqknn = list(k = c("k (neighbours)", 10, 1)),
  trknn   = list(k = c("k (neighbours)", 10, 1)), lls = list(k = c("k (neighbours)", 10, 1)),
  svd     = list(nPcs = c("Principal components", 3, 1)), bpca = list(nPcs = c("Principal components", 3, 1)))

imp_param_ui <- function(method, prefix) {
  sp <- IMP_PARAMS[[method]]
  if (is.null(sp)) return(NULL)
  tagList(lapply(names(sp), function(p)
    numericInput(paste0(prefix, method, "_", p), paste0(method, ": ", sp[[p]][1]),
                 as.numeric(sp[[p]][2]), step = as.numeric(sp[[p]][3]))))
}

IMP_CATEGORIES <- list(
  "Single value"     = c("zero", "minimum", "colmedian", "rowmedian", "mindet", "minprob", "qr", "perseus"),
  "Global structure" = c("svd", "bpca", "mle", "irmi", "mice_norm", "mice_cart", "impseq", "impseqrob"),
  "Local similarity" = c("knn", "seqknn", "lls", "rf"))

imp_method_label <- function(m) names(IMP_CHOICES)[match(m, IMP_CHOICES)]
imp_grouped_choices <- function()
  lapply(IMP_CATEGORIES, function(v) setNames(v, imp_method_label(v)))

# tick-box per method under bold category headers; each method's options appear right below it
imp_eval_methods_ui <- function(default = c("perseus", "svd", "knn")) {
  tagList(lapply(names(IMP_CATEGORIES), function(cat) tagList(
    tags$b(cat, style = "display:block; margin-top:10px;"),
    lapply(IMP_CATEGORIES[[cat]], function(m) tagList(
      checkboxInput(paste0("evm_", m), imp_method_label(m), value = m %in% default),
      conditionalPanel(sprintf("input.evm_%s", m), div(style = "margin-left:22px;", imp_param_ui(m, "ev_"))))))))
}
