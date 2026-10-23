# ======================= RESET WHEN NEW DATA ARRIVES =======================
# Loaded last: clears every downstream result when the data source changes.
observeEvent(list(input$expr_file, input$meta_file, input$load_demo),
             ignoreInit = TRUE, {
  if (!is.null(job$proc) && job$proc$is_alive()) { job$stopped <- TRUE; job$proc$kill_tree() }
  filt_rv(NULL); norm_rv(NULL); imp_rv(NULL); norm_label(NULL); imp_label(NULL); filter_msg("")
  qc$sel <- NULL; qc$bar <- qc$box <- qc$pca <- qc$venn <- NULL
  norm_cache(NULL); bench_sel(character(0)); report_pages(NULL); report_pdf(NULL)
  imp_eval(NULL)
  de_res(NULL); contrasts_rv(contrasts_rv()[0, ])
  pre_subsets(list()); custom_subsets(list())
  ora_dir(NULL); ora_sets_run(character(0))
  gsea_res(NULL); gsea_dir(NULL)
  expr_work(NULL); meta_work(NULL)
  updateNumericInput(session, "miss_thr", value = 70)
  app_gen(app_gen() + 1)      # re-creates tabs 2-6 and the data-page options with default values
})
