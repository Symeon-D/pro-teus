# ======================= 4. BATCH & SVA =======================
sva_m <- reactive(final_expr())
sva_meta <- reactive({ cs <- colnames(sva_m()); sample_meta()[cs, , drop = FALSE] })
sva_vars <- reactive(setdiff(names(sva_meta()), c("Sample", input$sample_col)))
sva_guess <- reactive({
  v <- sva_vars(); g <- guess_var_types(sva_meta()[, v, drop = FALSE])
  if ("Group" %in% names(g)) g[["Group"]] <- "exclude"          # the combined group duplicates its variables
  g
})
output$sva_types_ui <- renderUI({
  vars <- sva_vars(); g <- sva_guess()
  tagList(lapply(seq_along(vars), function(i)
    selectInput(paste0("svt_", i), vars[i], c("binary", "nominal", "continuous", "exclude"), selected = g[[vars[i]]], width = "100%")))
})
sva_types <- reactive({
  vars <- sva_vars(); g <- sva_guess()
  setNames(vapply(seq_along(vars), function(i) input[[paste0("svt_", i)]] %||% g[[vars[i]]], character(1)), vars)
})

sva_complete <- reactive(!anyNA(sva_m()))
output$sva_notice <- renderUI({
  m <- sva_m()
  if (!anyNA(m)) return(NULL)
  div(class = "alert alert-warning",
      sprintf("SVA and the post-correction check need a matrix without missing values: run imputation first (tab 3). Everything else on this page already works; it uses the %d proteins (of %d) without missing values.",
              sum(stats::complete.cases(m)), nrow(m)))
})

dl_plot <- function(id, plot_fun, fname, w = 8, h = 6)
  output[[id]] <- downloadHandler(filename = function() fname,
    content = function(file) ggplot2::ggsave(file, plot_fun(), width = w, height = h, dpi = 300, bg = "white"))

# ---- variable relationships ----
sva_varvar_obj <- reactive({
  r <- var_var_matrix(sva_meta(), sva_types())
  validate(need(!is.null(r), "Keep at least two variables (a type other than 'exclude')."))
  r$plot
})
output$sva_varvar <- renderPlot(sva_varvar_obj(), res = 100)
dl_plot("dl_sva_varvar", function() sva_varvar_obj(), "variable_relationships.png", 7, 6)

# ---- PCs vs covariates ----
sva_pc <- reactive({
  m <- sva_m()
  validate(need(sum(stats::complete.cases(m)) >= 20, "Fewer than 20 proteins without missing values: impute first."))
  r <- pc_covariates(m, sva_meta(), sva_types())
  validate(need(!is.null(r), "Not enough samples / variables for the PCA."))
  r
})
output$sva_pc_note <- renderUI(helpText(sprintf("PCA on %d proteins without missing values.", sva_pc()$n_rows)))
output$sva_pc_cor <- renderPlot(sva_pc()$plot_cor, res = 100)
output$sva_pc_imp <- renderPlot(sva_pc()$plot_impact, res = 100)
dl_plot("dl_sva_pc_cor", function() sva_pc()$plot_cor, "PC_covariate_correlation.png", 7, 6)
dl_plot("dl_sva_pc_imp", function() sva_pc()$plot_impact, "PC_covariate_impact.png", 7, 6)

# ---- variance partitioning (background job) ----
vp_pre <- reactiveVal(NULL); vp_post <- reactiveVal(NULL); post_res <- reactiveVal(NULL); sva_res <- reactiveVal(NULL)
vp_vars <- function() { ty <- sva_types(); list(cat = names(ty)[ty %in% c("binary", "nominal")], cont = names(ty)[ty == "continuous"]) }

observeEvent(input$sva_run_vp, {
  v <- vp_vars()
  if (!length(c(v$cat, v$cont))) return(showNotification("Keep at least one variable.", type = "warning"))
  if (sum(stats::complete.cases(sva_m())) < 20) return(showNotification("Fewer than 20 proteins without missing values: impute first.", type = "warning"))
  start_job("Variance partitioning", "run_varpart", list(sva_m(), sva_meta(), v$cat, v$cont, NULL),
            on_done = function(res) vp_pre(res))
})
output$sva_vp_ui <- renderUI({
  r <- vp_pre(); validate(need(!is.null(r), "Press 'Run variance partitioning'."))
  tagList(helpText(sprintf("%d proteins without missing values.", r$n_rows)),
          plotOutput("sva_vp_plot", height = "480px"), downloadButton("dl_sva_vp", "PNG"),
          h5("Median / mean variance explained (%)"), tableOutput("sva_vp_tbl"))
})
output$sva_vp_plot <- renderPlot(plot_varpart(vp_pre()$vp), res = 100)
output$sva_vp_tbl <- renderTable(vp_pre()$summary, digits = 2)
dl_plot("dl_sva_vp", function() plot_varpart(vp_pre()$vp), "variance_partitioning.png", 8, 6)

# ---- surrogate variables ----
sva_all_vars <- reactive(c("Group", setdiff(sva_vars(), "Group")))      # every variable, the created Group first
output$sva_protect_ui <- renderUI({
  v <- sva_all_vars(); batchy <- grep("batch", v, ignore.case = TRUE, value = TRUE)
  tagList(
    checkboxGroupInput("sva_protect", "Variables of interest (protected)", choices = v, selected = "Group"),
    checkboxGroupInput("sva_adjust", "Known technical covariates (adjusted for)", choices = setdiff(v, "Group"), selected = batchy),
    helpText("Full model: interest + known covariates. Null model: known covariates only."))
})
output$sva_run_ui <- renderUI({
  if (sva_complete()) actionButton("sva_run", "Run SVA", icon = icon("play"), class = "btn-primary")
  else tagList(tags$button(class = "btn btn-default", disabled = NA, "Run SVA"),
               helpText("Not runnable: the matrix has missing values - run imputation first."))
})
observeEvent(input$sva_run, {
  if (!sva_complete()) return(showNotification("SVA needs a matrix without missing values: run imputation first.", type = "warning"))
  prot <- input$sva_protect %||% character(0); adj <- setdiff(input$sva_adjust %||% character(0), prot)
  if (!length(prot)) return(showNotification("Choose at least one variable of interest (normally the group).", type = "warning"))
  md <- sva_meta(); d <- data.frame(row.names = rownames(md))
  for (p in c(prot, adj)) { x <- md[[p]]; d[[p]] <- if (is.numeric(x) && length(unique(x)) > 5) as.numeric(x) else factor(x) }
  nsv <- input$sva_nsv; if (is.null(nsv) || is.na(nsv)) nsv <- NULL
  start_job("Surrogate variable analysis", "run_sva", list(sva_m(), d, prot, adj, input$sva_method %||% "leek", nsv),
            on_done = function(res) {
              sva_res(res)
              if (res$n_used == 0) showNotification(sprintf("No surrogate variables detected (estimate = %d, n = 0).", res$n_est), type = "message", duration = 10)
            })
})
sva_assoc <- reactive({ r <- sva_res(); req(r$sv); sv_associations(r$sv, sva_meta(), sva_types()) })
sva_assoc_h <- function() { r <- sva_res(); n <- if (is.null(r$sv)) 1 else ncol(r$sv); max(320, 62 * n + 170) }   # px, grows with the SVs
output$sva_result_ui <- renderUI({
  r <- sva_res(); validate(need(!is.null(r), "Press 'Run SVA'."))
  if (r$n_used == 0) return(div(class = "alert alert-info", sprintf("No surrogate variables detected (the %s estimate is %d), so there is nothing to carry forward.", input$sva_method %||% "leek", r$n_est)))
  a <- sva_assoc()$table; gv <- input$group_vars %||% character(0)
  bad <- unique(a$SV[a$Variable %in% gv & a$sig %in% c("*", "**", "***")])
  tagList(
    helpText(sprintf("Estimated %d SV(s); %d computed.", r$n_est, r$n_used)),
    if (length(bad)) div(class = "alert alert-warning", style = "padding:8px;",
      sprintf("%s is associated with a grouping variable: it may be absorbing real biology - check before including it.", paste(bad, collapse = ", "))),
    plotOutput("sva_assoc_plot", height = "auto"), downloadButton("dl_sva_assoc", "PNG"),
    helpText("SV vs covariate: Spearman correlation (binary / continuous) or square root of the adjusted R-squared of an ANOVA (nominal); stars use uncorrected p-values."),
    wellPanel(style = "margin-top:22px;",
      checkboxGroupInput("sva_pick", "SVs to carry into the differential analysis (none preselected)", choices = colnames(r$sv),
                         selected = character(0), inline = TRUE)),
    div(style = "margin-top:24px;", h5("Associations"),
        div(style = "max-height:260px; overflow:auto; margin-bottom:24px;", tableOutput("sva_assoc_tbl"))))
})
output$sva_assoc_plot <- renderPlot(sva_assoc()$plot, res = 100, height = function() sva_assoc_h())
output$sva_assoc_tbl <- renderTable(sva_assoc()$table, digits = 3)
dl_plot("dl_sva_assoc", function() sva_assoc()$plot, "SV_covariate_associations.png", 7, max(3.2, sva_assoc_h() / 100))

# SVs chosen for DE (data.frame, rows = samples) or NULL
sva_sv_df <- reactive({
  r <- sva_res(); p <- input$sva_pick
  if (is.null(r$sv) || !length(p)) NULL else r$sv[, intersect(p, colnames(r$sv)), drop = FALSE]
})

# ---- post-correction check ----
output$sva_remove_ui <- renderUI({
  ty <- sva_types(); vars <- setdiff(names(ty)[ty != "exclude"], c(input$group_vars, "Group"))
  checkboxGroupInput("sva_remove", "Known covariates to remove as well (optional)", choices = vars, selected = character(0), inline = TRUE)
})
output$sva_pca_ui <- renderUI({
  ch <- c("Group", setdiff(names(sva_meta()), c("Group", "Sample", input$sample_col))); req(length(ch) > 0)
  ctl_row(
    selectInput("sva_pca_color", "Colour by", ch, selected = isolate(input$sva_pca_color) %||% "Group", width = fit_w("Colour by", ch)),
    selectInput("sva_pca_shape", "Shape by", c("none", ch), selected = isolate(input$sva_pca_shape) %||% "none", width = fit_w("Shape by", c("none", ch))))
})
observeEvent(input$sva_post_run, {
  known <- input$sva_remove %||% character(0); sv <- sva_sv_df()
  if (is.null(sv) && !length(known)) return(showNotification("Carry at least one SV (previous tab) or tick a covariate to remove.", type = "warning"))
  m <- sva_m(); md <- sva_meta()
  md2 <- if (!is.null(sv)) cbind(md, sv[rownames(md), , drop = FALSE]) else md
  covs <- c(known, if (!is.null(sv)) colnames(sv))
  post_res(list(corrected = correct_for_covariates(m, md2, covs, md$Group), covs = covs)); vp_post(NULL)
  v <- vp_vars()
  if (sum(stats::complete.cases(m)) >= 20 && length(c(v$cat, v$cont))) {
    run_post <- function() start_job("Variance partitioning (after correction)", "run_varpart", list(m, md, v$cat, v$cont, sv),
                                     on_done = function(res) vp_post(res))
    if (is.null(vp_pre()))                     # the 'before' run was never done: do it first, then the 'after' run
      start_job("Variance partitioning (before correction)", "run_varpart", list(m, md, v$cat, v$cont, NULL),
                on_done = function(res) { vp_pre(res); run_post() })
    else run_post()
  }
})
sva_pca_obj <- reactive({
  r <- post_res(); req(r)
  p <- pca_before_after(sva_m(), r$corrected, sva_meta(), input$sva_pca_color %||% "Group", input$sva_pca_shape %||% "none", group_colors_all())
  validate(need(!is.null(p), "Not enough complete rows for the PCA."))
  p
})
sva_pc_before_obj <- reactive(sva_pc()$plot_cor + ggplot2::ggtitle("PC ~ covariate correlations (before correction)"))
sva_pc_after_obj <- reactive({
  r <- post_res(); req(r)
  x <- pc_covariates(r$corrected, sva_meta(), sva_types()); req(x)
  x$plot_cor + ggplot2::ggtitle("PC ~ covariate correlations (after correction)")
})
output$sva_post_ui <- renderUI({
  validate(need(!is.null(post_res()), "Press 'Run post-correction check'."))
  tagList(helpText(sprintf("Removed: %s.", paste(post_res()$covs, collapse = ", "))),
          plotOutput("sva_pca_plot", height = "480px"), downloadButton("dl_sva_pca", "PNG"),
          h5("PC ~ covariate correlations"),
          fluidRow(column(6, plotOutput("sva_pc_before", height = "480px"), downloadButton("dl_sva_pc_before", "PNG")),
                   column(6, plotOutput("sva_pc_after", height = "480px"), downloadButton("dl_sva_pc_after", "PNG"))),
          h5("Variance partitioning: before vs after"), uiOutput("sva_vp_cmp_ui"))
})
output$sva_pca_plot <- renderPlot(sva_pca_obj(), res = 100)
output$sva_pc_before <- renderPlot(sva_pc_before_obj(), res = 100)
output$sva_pc_after <- renderPlot(sva_pc_after_obj(), res = 100)
dl_plot("dl_sva_pca", function() sva_pca_obj(), "PCA_before_after_correction.png", 12, 6)
dl_plot("dl_sva_pc_before", function() sva_pc_before_obj(), "PC_covariate_correlation_before.png", 7, 6)
dl_plot("dl_sva_pc_after", function() sva_pc_after_obj(), "PC_covariate_correlation_after.png", 7, 6)
output$sva_vp_cmp_ui <- renderUI({
  if (is.null(vp_pre()) || is.null(vp_post()))
    return(helpText("Variance partitioning is running (before, then after) - or it needs at least 20 proteins without missing values."))
  tagList(fluidRow(column(6, plotOutput("sva_vp_before", height = "440px"), downloadButton("dl_sva_vp_before", "PNG")),
                   column(6, plotOutput("sva_vp_after", height = "440px"), downloadButton("dl_sva_vp_after", "PNG"))),
          h5("Median / mean variance explained (%)"), tableOutput("sva_cmp_tbl"))
})
output$sva_vp_before <- renderPlot(plot_varpart(vp_pre()$vp, "Before correction"), res = 100)
output$sva_vp_after <- renderPlot(plot_varpart(vp_post()$vp, "After correction (SVs as fixed terms)"), res = 100)
dl_plot("dl_sva_vp_before", function() plot_varpart(vp_pre()$vp, "Before correction"), "variance_partitioning_before.png", 8, 6)
dl_plot("dl_sva_vp_after", function() plot_varpart(vp_post()$vp, "After correction (SVs as fixed terms)"), "variance_partitioning_after.png", 8, 6)
output$sva_cmp_tbl <- renderTable({
  req(vp_pre(), vp_post())
  merge(vp_pre()$summary, vp_post()$summary, by = "Variable", all = TRUE, suffixes = c("_Before", "_After"))
}, digits = 2)

observeEvent(input$sva_skip, updateNavbarPage(session, "nav", selected = "5. Differential analysis"))
