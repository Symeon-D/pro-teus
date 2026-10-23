# ======================= 2. NORMALIZATION =======================
norm_cache   <- reactiveVal(NULL)
bench_sel    <- reactiveVal(character(0))
report_pages <- reactiveVal(NULL)   # PNG pages of the NormalyzerDE PDF report
report_pdf   <- reactiveVal(NULL)

# ---- "is normalization needed?" as compact cards ----
output$norm_need_ui <- renderUI({
  m <- norm_input()
  res <- tryCatch(qc_normalization_need(m, verbose = FALSE), error = function(e) NULL)
  validate(need(!is.null(res), "Could not compute the diagnostics for this matrix."))
  lvl <- function(x, lo, hi) if (x < lo) "low" else if (x < hi) "moderate" else "high"
  col <- c(low = "#5cb85c", moderate = "#f0ad4e", high = "#d9534f")
  card <- function(title, value, level, hint)
    div(style = "flex:1 1 220px; border:1px solid #ddd; border-radius:6px; padding:10px;",
        div(style = "font-size:12px; color:#777;", title),
        div(style = "font-size:24px; font-weight:600;", value),
        tags$span(class = "label", style = paste0("background:", col[[level]], ";"), toupper(level)),
        div(style = "font-size:12px; color:#555; margin-top:6px;", hint))
  cm <- res$cor_median_missing
  cl <- if (is.na(cm)) "low" else if (abs(cm) < 0.4) "low" else if (abs(cm) < 0.6) "moderate" else "high"
  verdict_level <- if (all(c(lvl(res$median_range, .2, .5), lvl(res$iqr_range, .1, .25), cl) == "low")) "low"
                   else if (any(c(lvl(res$median_range, .2, .5), lvl(res$iqr_range, .1, .25), cl) == "high")) "high" else "moderate"
  tagList(
    div(style = paste0("border-left:6px solid ", col[[verdict_level]], "; background:#f8f8f8; padding:10px 14px; margin-bottom:12px;"),
        strong("Verdict: "), res$decision),
    div(style = "display:flex; flex-wrap:wrap; gap:12px;",
      card("Spread of sample medians", sprintf("%.2f log2", res$median_range),
           lvl(res$median_range, .2, .5), sprintf("= %.2f-fold between the extreme samples", res$median_fold_change_equivalent)),
      card("Spread of sample IQRs", sprintf("%.2f", res$iqr_range), lvl(res$iqr_range, .1, .25),
           "Differences in distribution width"),
      card("Median vs missingness", if (is.na(cm)) "n/a" else sprintf("r = %.2f", cm), cl,
           "Samples with lower medians tend to have more missing values")),
    tags$button(type = "button", class = "btn btn-default btn-sm", style = "margin-top:12px;",
                `data-toggle` = "collapse", `data-target` = "#norm_interp", "Full interpretation"),
    div(id = "norm_interp", class = "collapse", style = "margin-top:8px;",
        tags$pre(style = "white-space:pre-wrap;", res$interpretation)))
})

# ---- NormalyzerDE in a background process (so it can be stopped) ----
run_norm_job <- function(after) {
  m <- norm_input(); md <- sample_meta()[colnames(m), , drop = FALSE]
  is_log <- isTRUE(input$norm_is_log)
  key <- paste(stage_key(m), is_log); cache <- norm_cache()
  if (!is.null(cache) && identical(cache$key, key)) { after(cache); return(invisible()) }
  out <- file.path(session_tmp, paste0("normalyzer_", format(Sys.time(), "%H%M%S")))
  lin <- if (is_log) 2^m else m                      # linear scale in; NormalyzerDE returns log2
  start_job("NormalyzerDE", "run_normalyzerDE",
    list(lin, md, noLogTransform = FALSE, sample_col = "Sample", group_col = "Group",
         jobName = "norm", outputDir = out, keep_input_files = TRUE),
    on_done = function(res) {
      files <- list.files(out, pattern = "-normalized\\.txt$", recursive = TRUE, full.names = TRUE)
      mats <- lapply(files, read_norm_file, ids = rownames(m), samples = colnames(m))
      names(mats) <- sub("-normalized\\.txt$", "", basename(files))
      mats <- mats[!vapply(mats, is.null, logical(1))]
      if (!length(mats)) stop("No normalised matrices could be read from: ", out)
      cache <- list(key = key, dir = out, mats = mats, input = if (is_log) m else safe_log2(m))
      norm_cache(cache)
      # NormalyzerDE report -> page images (fallback: embedded PDF)
      pdfs <- list.files(out, pattern = "\\.pdf$", recursive = TRUE, full.names = TRUE)
      report_pdf(if (length(pdfs)) pdfs[1] else NULL); report_pages(NULL)
      if (length(pdfs)) {
        addResourcePath("normrep", dirname(pdfs[1]))
        if (requireNamespace("pdftools", quietly = TRUE)) {
          pg <- tryCatch({
            n <- pdftools::pdf_info(pdfs[1])$pages
            pdftools::pdf_convert(pdfs[1], format = "png", dpi = 90, pages = seq_len(n),
                                  filenames = file.path(dirname(pdfs[1]), sprintf("report_page_%03d.png", seq_len(n))),
                                  verbose = FALSE)
          }, error = function(e) NULL)
          report_pages(pg)
        }
      }
      after(cache)
    })
}

observeEvent(input$run_norm_bench, {
  run_norm_job(function(cache) updateTabsetPanel(session, "norm_tabs", selected = "Benchmark"))
})

output$norm_bench_ui <- renderUI({
  validate(need(!is.null(norm_cache()), "Press 'Run benchmark' (all methods are evaluated)."))
  uiOutput("norm_report_ui")
})

output$norm_report_ui <- renderUI({
  pg <- report_pages(); pdf <- report_pdf()
  if (!is.null(pg) && length(pg))
    return(tagList(
      downloadButton("dl_norm_report", "Download PDF report"), br(), br(),
      div(style = "display:grid; grid-template-columns:1fr; gap:8px;",
          lapply(pg, function(f) tags$img(src = paste0("normrep/", basename(f)),
                                          style = "width:100%; border:1px solid #ddd;")))))
  if (!is.null(pdf))
    return(tagList(downloadButton("dl_norm_report", "Download PDF report"), br(), br(),
                   tags$iframe(src = paste0("normrep/", basename(pdf)), width = "100%", height = "800px",
                               style = "border:1px solid #ddd;")))
  helpText("NormalyzerDE did not produce a PDF report.")
})

output$dl_norm_report <- downloadHandler("NormalyzerDE_report.pdf", function(file) file.copy(report_pdf(), file))

# ---- apply one method ----
observeEvent(input$apply_norm, {
  choice <- input$norm_choice
  finish <- function() {
    imp_rv(NULL); imp_label(NULL)
    updateTabsetPanel(session, "norm_tabs", selected = "Applied result")
  }
  if (choice == "none") {
    norm_rv(norm_input()); norm_label("none"); finish()
  } else {
    run_norm_job(function(cache) {
      nm <- names(cache$mats)[tolower(names(cache$mats)) == tolower(choice)]
      if (!length(nm)) stop("Method '", choice, "' was not produced by NormalyzerDE.")
      norm_rv(cache$mats[[nm[1]]]); norm_label(nm[1]); finish()
    })
  }
})

output$norm_applied_txt <- renderText({
  if (is.null(norm_rv())) "No normalization applied yet." else paste("Applied:", norm_label())
})
output$norm_after_plot <- renderPlot({
  post <- norm_rv(); req(post); pre <- norm_input(); md <- sample_meta()
  long <- function(m, stage) data.frame(
    Stage = stage, Sample = factor(rep(colnames(m), each = nrow(m)), levels = colnames(post)),
    Group = rep(md[colnames(m), "Group"], each = nrow(m)), Value = as.vector(m))
  df <- rbind(long(pre[, colnames(post), drop = FALSE], "Before normalization"),
              long(post, paste("After:", norm_label())))
  df$Stage <- factor(df$Stage, levels = unique(df$Stage))
  # no normalization: before == after, so show only the un-normalized boxplot
  if (identical(norm_label(), "none")) df <- df[df$Stage == "Before normalization", ]
  ggplot(df, aes(Sample, Value, fill = Group)) + geom_boxplot(outlier.size = 0.3, linewidth = 0.3) +
    facet_wrap(~Stage, ncol = 2, scales = "free_y") +
    scale_fill_manual(values = group_colors_all()) + labs(y = "log2 intensity") +
    theme_bw(base_size = 12) + theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 7),
                                     legend.position = "top")
}, res = 100)

# ---- skip: keep the data as they are and go straight to the imputation tab ----
observeEvent(input$skip_norm, {
  norm_rv(norm_input()); norm_label("skipped")
  imp_rv(NULL); imp_label(NULL)
  updateNavbarPage(session, "nav", selected = "3. Imputation")
})
