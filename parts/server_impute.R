# ======================= 3. IMPUTATION =======================
imp_eval <- reactiveVal(NULL)

output$imp_na_txt <- renderText({
  m <- imp_input()
  sprintf("Input to imputation: %d proteins x %d samples, %d missing values (%.1f%%).",
          nrow(m), ncol(m), sum(is.na(m)), 100 * mean(is.na(m)))
})

# ---- missingness pattern: available from the start ----
output$imp_miss_plot <- renderPlot({
  m <- imp_input(); req(m)
  pm <- rowMeans(m, na.rm = TRUE); miss_f <- rowMeans(is.na(m))
  d_prot <- data.frame(mean_int = pm, miss = miss_f, has_na = ifelse(miss_f > 0, "With missing values", "Complete"))
  d_prot <- d_prot[is.finite(d_prot$mean_int), ]
  md <- sample_meta()
  d_samp <- data.frame(Sample = factor(colnames(m), levels = colnames(m)),
                       Missing = 100 * colMeans(is.na(m)), Group = md[colnames(m), "Group"])
  p1 <- ggplot(d_prot, aes(mean_int, fill = has_na, colour = has_na)) +
    geom_density(alpha = 0.35) +
    scale_fill_manual(values = c("Complete" = "grey50", "With missing values" = "#E07B54")) +
    scale_colour_manual(values = c("Complete" = "grey30", "With missing values" = "#B8532E")) +
    labs(title = "Intensity density: complete vs incomplete proteins",
         x = "Mean log2 intensity per protein (observed values)", y = "Density", fill = NULL, colour = NULL) +
    theme_bw(base_size = 12) + theme(legend.position = "top")
  p2 <- ggplot(d_prot[d_prot$miss > 0, ], aes(mean_int, 100 * miss)) +
    geom_point(alpha = 0.25, size = 0.8, colour = "#B8532E") +
    geom_smooth(method = "loess", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.8) +
    labs(title = "Missingness vs intensity (MNAR shows as a downward trend)",
         x = "Mean log2 intensity per protein", y = "% missing") + theme_bw(base_size = 12)
  # bottom row = the two graphs of the first evaluation page (missing-value heatmap + % missing per sample)
  cowplot::plot_grid(cowplot::plot_grid(p1, p2, nrow = 1), plot_missing_pattern(m, groups = sample_meta()[colnames(m), "Group"], group_colors = group_colors_all()),
                     ncol = 1, rel_heights = c(1, 1.5))
}, res = 100)

# ---- method-specific options (only what each method uses; k only for k-methods) ----
ev_methods <- reactive(Filter(function(m) isTRUE(input[[paste0("evm_", m)]]),
                              unlist(IMP_CATEGORIES, use.names = FALSE)))
output$imp_apply_opts <- renderUI(imp_param_ui(input$imp_final, "ap_"))

imp_args_for <- function(method, prefix) {
  sp <- IMP_PARAMS[[method]]
  if (is.null(sp)) return(list())
  setNames(lapply(names(sp), function(p) input[[paste0(prefix, method, "_", p)]] %||% as.numeric(sp[[p]][2])), names(sp))
}

# ---- evaluation (background job) ----
observeEvent(input$run_imp_eval, {
  req(length(ev_methods()) > 0)
  m <- imp_input(); md <- sample_meta()[colnames(m), , drop = FALSE]
  if (!anyNA(m)) return(showNotification("There are no missing values to impute.", type = "warning"))
  ms <- ev_methods()
  start_job("Imputation evaluation", "evaluate_imputation",
    list(data = m, sample_info = md, sample_col = "Sample", group_col = "Group",
         group_colors = group_colors_all(), methods = ms, mask_fraction = 0.15,
         log2_transform = FALSE, run_de = isTRUE(input$imp_run_de), seed = input$global_seed %||% 42,
         output_dir = NULL, impute_args = setNames(lapply(ms, imp_args_for, prefix = "ev_"), ms)),
    on_done = function(res) {
      imp_eval(res)
      updateTabsetPanel(session, "imp_tabs", selected = "Evaluation")
    })
})

for (nm in c("density", "pca", "nrmse", "cv", "downshift", "classic", "de_concordance")) local({
  n <- nm
  output[[paste0("imp_", n)]] <- renderPlot({
    p <- imp_eval()[[paste0("plot_", n)]]; req(p); print(p)
  }, res = 90)
})

output$dl_imp_pdf <- downloadHandler("imputation_evaluation.pdf", function(file) {
  res <- imp_eval(); req(res)
  pdf(file, width = 14, height = 9)
  for (n in c("missing", "density", "pca", "nrmse", "cv", "downshift", "classic", "de_concordance"))
    if (!is.null(res[[paste0("plot_", n)]])) print(res[[paste0("plot_", n)]])
  dev.off()
})

# ---- apply (background job) ----
observeEvent(input$apply_imp, {
  m <- imp_input(); meth <- input$imp_final
  if (!anyNA(m)) return(showNotification("No missing values - nothing to impute.", type = "warning"))
  start_job(paste("Imputation:", meth), "my_impute",
    c(list(input_matrix = m, method = meth, seed = (input$global_seed %||% 42), verbose = FALSE),
      imp_args_for(meth, "ap_")),
    on_done = function(res) {
      imp_rv(res); imp_label(meth)
      updateTabsetPanel(session, "imp_tabs", selected = "Applied result")
    })
})
observeEvent(input$skip_imp, {
  imp_rv(imp_input()); imp_label("none (missing values kept)")
  updateNavbarPage(session, "nav", selected = "4. Differential analysis")
})

output$imp_applied_txt <- renderText(if (is.null(imp_rv())) "No imputation applied yet." else paste("Applied:", imp_label()))
output$imp_after_plot <- renderPlot({
  m1 <- imp_rv(); req(m1); m0 <- imp_input()
  obs <- m0[!is.na(m0)]; imputed <- m1[is.na(m0)]; imputed <- imputed[is.finite(imputed)]
  req(length(imputed) > 0)
  df <- rbind(data.frame(Type = "Observed", Value = obs), data.frame(Type = "Imputed", Value = imputed))
  ggplot(df, aes(Value, fill = Type)) + geom_density(alpha = 0.5) +
    scale_fill_manual(values = c(Observed = "grey40", Imputed = "#E07B54")) +
    labs(title = paste("Observed vs imputed values:", imp_label()), x = "log2 intensity") + theme_bw(base_size = 12)
}, res = 100)
