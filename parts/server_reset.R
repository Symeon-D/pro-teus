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
  sva_res(NULL); vp_pre(NULL); vp_post(NULL); post_res(NULL)
  expr_work(NULL); meta_work(NULL)
  updateNumericInput(session, "miss_thr", value = 70)
  app_gen(app_gen() + 1)      # re-creates tabs 2-6 and the data-page options with default values
})


# ======================= SAMPLES REMOVED / RESTORED =======================
# The matched sample list changes: pipeline stages restart from the raw (reduced) data and every result is cleared.
observeEvent(excluded(), ignoreInit = TRUE, {
  if (!is.null(job$proc) && job$proc$is_alive()) { job$stopped <- TRUE; job$proc$kill_tree() }
  filt_rv(NULL); norm_rv(NULL); imp_rv(NULL); norm_label(NULL); imp_label(NULL); filter_msg(""); norm_skip_msg(character(0))
  norm_cache(NULL); report_pages(NULL); report_pdf(NULL)
  imp_eval(NULL)
  de_res(NULL); pre_subsets(list()); custom_subsets(list())
  ora_dir(NULL); ora_sets_run(character(0)); gsea_res(NULL); gsea_dir(NULL)
  sva_res(NULL); vp_pre(NULL); vp_post(NULL); post_res(NULL)
  rmd_res(NULL)
  md <- tryCatch(sample_meta()[setdiff(all_matched(), excluded()), , drop = FALSE], error = function(e) NULL)
  if (!is.null(md)) {
    cnt <- table(md$Group); small <- names(cnt)[cnt < 2]
    if (length(small))
      showNotification(paste0("Warning: group(s) with fewer than 2 samples: ", paste(small, collapse = ", ")), type = "warning", duration = 10)
  }
})
observeEvent(list(input$expr_file, input$meta_file, input$load_demo), ignoreInit = TRUE, {
  excluded(character(0)); rmd_res(NULL); norm_skip_msg(character(0))
})
