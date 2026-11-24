# ======================= 1. DATA & QC =======================
demo_on <- reactiveVal(FALSE)
demo_cache <- reactive(make_demo_data())
observeEvent(input$load_demo, demo_on(TRUE))
observeEvent(input$expr_file, demo_on(FALSE), ignoreInit = TRUE)
observeEvent(input$meta_file, demo_on(FALSE), ignoreInit = TRUE)

det_fmt <- reactiveVal("")
id_label <- reactiveVal("Names")          # name of the gene / protein column in the user's file
# tables / downloads show the real column name instead of a generic one
as_names <- function(df) { if ("Protein" %in% names(df)) names(df)[names(df) == "Protein"] <- id_label(); df }

# gene-name column: auto-detected from the file, with a selector as fallback / override
expr_info <- reactive({
  req(input$expr_file)
  tryCatch(inspect_expression_file(input$expr_file$datapath, input$expr_file$name),
           error = function(e) read_fail("expression file", conditionMessage(e)))
})
output$gene_col_ui <- renderUI({
  if (demo_on()) return(NULL)
  req(input$expr_file); info <- expr_info()
  nm <- info$cols[info$candidates]
  ch <- setNames(as.character(info$candidates), ifelse(nzchar(nm), nm, "(unnamed first column)"))
  tagList(
    selectInput("gene_col", "Column holding gene names", choices = c("- select -" = "", ch),
                selected = if (is.na(info$auto)) "" else as.character(info$auto)),
    if (is.na(info$auto))
      div(class = "alert alert-warning", style = "padding:6px; font-size:12px;",
          "The gene-name column could not be identified automatically. Please indicate which column holds the gene names."))
})
read_fail <- function(what, msg) {
  showNotification(paste0("Could not read the ", what, ": ", msg), type = "error", duration = NULL)
  req(FALSE)       # stop quietly - the app stays responsive
}

expr_base <- reactive({
  if (demo_on()) {
    det_fmt("demo data"); id_label("Genes")
    m <- demo_cache()$expr
    if (isTRUE(input$log2)) m <- safe_log2(m)
    return(m)
  }
  req(input$expr_file)
  gc <- input$gene_col; req(!is.null(gc))
  if (!nzchar(gc)) {
    showNotification("Please indicate which column holds the gene names ('Column holding gene names').", type = "warning")
    req(FALSE)
  }
  r <- tryCatch(read_expression_file(input$expr_file$datapath, input$expr_file$name, id_col = as.integer(gc)),
                error = function(e) e)
  if (inherits(r, "error")) { det_fmt("")
    read_fail("expression file", paste0(conditionMessage(r),
      " Expected a protein matrix, MaxQuant proteinGroups.txt, DIA-NN pg_matrix or a Perseus export.")) }
  m <- r$m
  if (nrow(m) < 2 || ncol(m) < 2 || all(is.na(m))) { det_fmt("")
    read_fail("expression file", "it contains fewer than 2 proteins/samples or no numeric values.") }
  det_fmt(format_label(r$format))
  id_label(if (!is.null(r$id_name) && nzchar(r$id_name)) r$id_name else "Names")
  if (!r$log2_applied && isTRUE(input$log2)) m <- safe_log2(m)
  m
})

meta_base <- reactive({
  if (demo_on()) return(demo_cache()$meta)
  req(input$meta_file)
  d <- tryCatch(read_table_any(input$meta_file$datapath, input$meta_file$name), error = function(e) e)
  if (inherits(d, "error")) read_fail("sample info file", conditionMessage(d))
  if (ncol(d) < 2 || nrow(d) < 2) read_fail("sample info file", "it needs a sample-name column and at least one variable column.")
  d
})

output$sample_col_ui <- renderUI({
  req(meta_raw())
  cn <- colnames(meta_raw())
  keep <- isolate(input$sample_col)
  selectInput("sample_col", "Column holding sample names", cn,
              selected = if (!is.null(keep) && keep %in% cn) keep else if ("Sample" %in% cn) "Sample" else cn[1])
})

other_cols <- reactive({
  req(meta_raw(), input$sample_col)
  setdiff(colnames(meta_raw()), input$sample_col)
})

output$group_vars_ui <- renderUI({
  req(other_cols())
  oc <- other_cols()
  sel <- if (demo_on()) intersect(c("Genotype", "Condition"), oc) else oc[1]
  keep <- isolate(input$group_vars); if (length(keep) && all(keep %in% oc)) sel <- keep
  selectizeInput("group_vars", "Variable(s) defining the group (one or more)",
                 choices = oc, multiple = TRUE, selected = sel,
                 options = list(plugins = list("remove_button")))
})

sample_meta <- reactive({
  req(meta_raw(), input$sample_col, input$group_vars)
  md <- meta_raw()
  req(input$sample_col %in% colnames(md))
  gvs <- intersect(input$group_vars, colnames(md))      # a selection left over from the previous file is ignored
  req(length(gvs) > 0)
  md$Sample <- as.character(md[[input$sample_col]])
  gv <- md[, gvs, drop = FALSE]
  md$Group <- do.call(paste, c(lapply(gv, as.character), sep = " "))   # single space
  md <- md[!duplicated(md$Sample), ]
  rownames(md) <- md$Sample
  md
})

# samples removed from the analysis (Sample info tab, or the rMd outlier tab). Removing them here, where the
# expression columns are matched to the sample info, makes EVERY tab ignore them.
excluded <- reactiveVal(character(0))
all_matched <- reactive({
  req(expr_raw(), sample_meta())
  intersect(colnames(expr_raw()), rownames(sample_meta()))
})
common <- reactive({
  cs <- setdiff(all_matched(), excluded())
  validate(need(length(cs) > 1, if (length(all_matched()) > 1 && length(excluded()))
                  "Fewer than two samples are left: restore some samples in the 'Sample info' tab."
                else "No overlap between expression column names and sample names."))
  cs
})

output$excl_pick_ui <- renderUI({
  am <- all_matched(); ex <- excluded()
  selectizeInput("excl_pick", NULL, multiple = TRUE, width = "100%",
                 choices = setNames(am, ifelse(am %in% ex, paste0(am, "  (removed)"), am)),
                 options = list(placeholder = "choose samples...", plugins = list("remove_button")))
})
observeEvent(input$excl_add, {
  pk <- input$excl_pick; req(length(pk))
  excluded(union(excluded(), pk))
})
observeEvent(input$excl_restore, {
  pk <- input$excl_pick; req(length(pk))
  excluded(setdiff(excluded(), pk))
})
observeEvent(input$excl_restore_all, excluded(character(0)))
output$excl_status <- renderUI({
  ex <- excluded()
  am <- tryCatch(all_matched(), error = function(e) character(0)); req(length(am) > 0)
  md <- sample_meta()
  cnt_all <- table(md[am, "Group"]); cnt <- table(factor(md[setdiff(am, ex), "Group"], levels = names(cnt_all)))
  small <- names(cnt)[cnt < 2]
  tagList(
    if (length(ex)) tags$p(tags$b(sprintf("Removed from the analysis (%d): ", length(ex))), paste(ex, collapse = ", "))
    else tags$p(style = "color:#777;", "No sample is removed."),
    if (length(small))
      div(class = "alert alert-warning", style = "padding:6px; font-size:12px;",
          "Group(s) with fewer than 2 samples: ", paste0(small, " (", cnt[small], ")", collapse = ", "),
          ". Statistics that need replicates (limma contrasts, rMd, ...) cannot use them."))
})

output$data_status <- renderUI({
  req(expr_raw(), sample_meta())
  m <- expr_raw()
  tagList(
    if (nzchar(det_fmt())) tags$p(tags$b("Detected format: "), det_fmt()),
    tags$p(sprintf("%d proteins x %d columns | %d samples matched with sample info | %d groups | %.1f%% missing values.",
                   nrow(m), ncol(m), length(common()),
                   length(unique(sample_meta()[common(), "Group"])), 100 * mean(is.na(m)))),
    if (length(setdiff(colnames(m), rownames(sample_meta()))))
      tags$p(style = "color:#a94442;", "Not in sample info: ",
             paste(head(setdiff(colnames(m), rownames(sample_meta())), 6), collapse = ", ")))
})

# ---------------- group colours (set here, used by QC and heatmap) ----------------
group_levels_all <- reactive(unique(sample_meta()[common(), "Group"]))

output$group_colors_ui <- renderUI({
  lv <- group_levels_all(); d <- auto_colors(lv)
  tagList(h4("Group colours"),
          div(style = "display:flex; flex-wrap:wrap; gap:8px;",
              lapply(seq_along(lv), function(i)
                div(style = "width:48%;", colourInput(paste0("gcol_", i), lv[i], d[[i]])))))
})
group_colors_all <- reactive({
  lv <- group_levels_all(); d <- auto_colors(lv)
  setNames(vapply(seq_along(lv), function(i) input[[paste0("gcol_", i)]] %||% d[[i]], character(1)), lv)
})

# ---------------- pipeline state ----------------
raw_aligned <- reactive(expr_raw()[, common(), drop = FALSE])
filt_rv <- reactiveVal(NULL); norm_rv <- reactiveVal(NULL); imp_rv <- reactiveVal(NULL)
norm_label <- reactiveVal(NULL); imp_label <- reactiveVal(NULL); filter_msg <- reactiveVal("")
norm_skip_msg <- reactiveVal(character(0))       # "Normalization skipped: ..." shown in the status bar

observeEvent(raw_aligned(), {
  filt_rv(NULL); norm_rv(NULL); imp_rv(NULL); norm_skip_msg(character(0))
  norm_label(NULL); imp_label(NULL); filter_msg("")
})

norm_input <- reactive(filt_rv() %||% raw_aligned())
imp_input  <- reactive(norm_rv() %||% norm_input())
final_expr <- reactive(imp_rv() %||% norm_rv() %||% filt_rv() %||% raw_aligned())

status_ui <- renderUI({
  req(raw_aligned())
  dimtxt <- function(m) sprintf("%d x %d", nrow(m), ncol(m))
  sv_now <- tryCatch(sva_sv_df(), error = function(e) NULL)
  badge <- function(txt, on) tags$span(class = paste("label", if (on) "label-success" else "label-default"),
                                       style = "margin-right:4px; font-size:12px;", txt)
  tags$div(style = "margin-bottom:10px;",
    badge(paste("Raw", dimtxt(raw_aligned())), TRUE), "\u2192",
    badge(if (is.null(filt_rv())) "Filter: not applied" else paste("Filtered", dimtxt(filt_rv())), !is.null(filt_rv())), "\u2192",
    badge(if (is.null(norm_rv())) "Normalization: not applied" else paste("Normalized:", norm_label()), !is.null(norm_rv())), "\u2192",
    badge(if (is.null(imp_rv())) "Imputation: not applied" else paste("Imputed:", imp_label()), !is.null(imp_rv())), "\u2192",
    badge(if (is.null(sv_now)) "SVs: none" else paste("SVs:", paste(colnames(sv_now), collapse = ", ")), !is.null(sv_now)),
    if (length(excluded()))
      tags$span(class = "label label-warning", style = "margin-left:4px; font-size:12px;",
                sprintf("%d sample(s) removed", length(excluded()))),
    if (length(norm_skip_msg()))
      tags$span(class = "label label-danger", style = "margin-left:4px; font-size:12px;", norm_skip_msg()),
    tags$small(style = "color:#777; margin-left:8px;", "Downstream tabs use the latest stage."))
})
for (id in c("pipeline_status_data", "pipeline_status_norm", "pipeline_status_imp", "pipeline_status_sva", "pipeline_status_de"))
  output[[id]] <- status_ui

# ---------------- missingness filter ----------------
observeEvent(input$run_filter, {
  tryCatch({
    m <- raw_aligned(); md <- sample_meta()[colnames(m), , drop = FALSE]
    capture.output(res <- filter_by_group_missingness(m, md, threshold = input$miss_thr / 100,
                                                      group_column = "Group"))
    filt_rv(res); norm_rv(NULL); imp_rv(NULL); norm_label(NULL); imp_label(NULL)
    filter_msg(sprintf("Kept %d of %d proteins (>= %g%% valid values in at least one group).",
                       nrow(res), nrow(m), input$miss_thr))
  }, error = function(e) showNotification(conditionMessage(e), type = "error"))
})
observeEvent(input$reset_filter, {
  filt_rv(NULL); norm_rv(NULL); imp_rv(NULL); norm_label(NULL); imp_label(NULL)
  filter_msg("Filter removed.")
})
output$filter_msg <- renderText(filter_msg())

# ---------------- QC plots: created automatically, laid out as a 2 x 2 grid ----------------
data_ready <- reactive(tryCatch({ expr_raw(); sample_meta(); common(); TRUE }, error = function(e) FALSE))

output$qc_stage_ui <- renderUI({
  ch <- c("Raw" = "raw", if (!is.null(filt_rv())) c("Filtered" = "filtered"))
  radioButtons("qc_stage", "Data shown in the QC plots", ch,
               selected = if (!is.null(filt_rv())) "filtered" else "raw", inline = TRUE)
})
qc_matrix <- reactive(if (identical(input$qc_stage, "filtered") && !is.null(filt_rv())) filt_rv() else raw_aligned())

output$qc_opts_ui <- renderUI({
  app_gen()                      # re-created (defaults restored) whenever the app is reset
  req(data_ready())
  tagList(
    h4("QC plot options"),
    uiOutput("qc_stage_ui"),
    h5("PCA"),
    fluidRow(
      column(6, selectInput("pca_mv", "Missing values", c("ppca", "nipals", "pcana"))),
      column(6, checkboxInput("pca_names", "Sample names", TRUE),
                checkboxInput("pca_ellipse", "Ellipses", FALSE))),
    h5("Venn (detected proteins)"),
    uiOutput("venn_groups_ui"),
    uiOutput("venn_titles_ui"),
    numericInput("venn_thr", "Detected in >= % of samples", 70, 0, 100, 5),
    sliderInput("venn_alpha", "Fill opacity", 0.1, 0.9, 0.45, step = 0.05),
    sliderInput("venn_font", "Overall font scale", 0.5, 3, 1, step = 0.1),
    numericInput("venn_count", "Number size (blank = auto)", NA, min = 1, step = 0.5),
    numericInput("venn_label", "Title size (blank = auto)", NA, min = 1, step = 0.5),
    sliderInput("venn_ratio", "Circle size vs numbers (higher = smaller circles)", 1, 3, 1.5, step = 0.1))
})

output$venn_groups_ui <- renderUI({
  lv <- unique(sample_meta()[common(), "Group"])
  selectizeInput("venn_groups", "Groups to compare (2-7)", choices = lv, multiple = TRUE,
                 selected = lv[seq_len(min(2, length(lv)))], options = list(maxItems = 7))
})
# one title box per circle; blank is allowed
output$venn_titles_ui <- renderUI({
  g <- input$venn_groups; req(length(g) > 0)
  tagList(lapply(seq_along(g), function(i)
    textInput(paste0("vt_", i), paste("Title of circle", i), value = g[i])))
})

qc <- reactiveValues(sel = NULL, bar = NULL, box = NULL, pca = NULL, venn = NULL)   # kept for the reset logic

qc_group_md <- reactive({ m <- qc_matrix(); sample_meta()[colnames(m), , drop = FALSE] })

qc_bar_path <- reactive({
  req(data_ready()); m <- qc_matrix()
  run_in_tmpdir(function() generate_barplot(m, qc_group_md(), group_colors_all(), plot_prefix = "qc",
                                            width = 4200, height = 2400, ceiling_num = 100))$files[1]
})
qc_box_path <- reactive({
  req(data_ready()); m <- qc_matrix()
  run_in_tmpdir(function() generate_boxplot(m, qc_group_md(), group_colors_all(), plot_prefix = "qc"))$files[1]
})
qc_pca_obj <- reactive({
  req(data_ready()); m <- qc_matrix()
  res <- tryCatch(
    my_plot_PCA_matrix(m, qc_group_md(), pca_method = "classic", top_n = "all",
                       scale = TRUE, show_names = isTRUE(input$pca_names %||% TRUE),
                       add_ellipse = isTRUE(input$pca_ellipse), fill_ellipse = FALSE, grid_lines = FALSE,
                       width = 9, height = 7, output_path = tempfile(fileext = ".png"),
                       group_colors = group_colors_all(), color_by = "Group",
                       mv_handle = input$pca_mv %||% "ppca"),
    error = function(e) e)
  validate(need(!inherits(res, "error"), if (inherits(res, "error")) paste("PCA:", conditionMessage(res)) else ""))
  res$pca_plot
})
# Venn: reactive, redraws whenever any option changes
qc_venn_obj <- reactive({
  req(data_ready())
  m <- qc_matrix(); md <- qc_group_md()
  g <- input$venn_groups
  validate(need(length(g) >= 2, "Choose at least two groups for the Venn diagram."))
  titles <- vapply(seq_along(g), function(i) input[[paste0("vt_", i)]] %||% g[i], character(1))
  na_null <- function(x) if (is.null(x) || is.na(x)) NULL else x
  # group order follows the selection so titles match the circles
  mdv <- md[md$Group %in% g, , drop = FALSE]
  mdv$Group <- factor(mdv$Group, levels = g); mdv <- mdv[order(mdv$Group), ]; mdv$Group <- as.character(mdv$Group)
  res <- tryCatch(run_in_tmpdir(function() venn_of_detected(
    m, mdv, threshold = (input$venn_thr %||% 70) / 100, group_colors = group_colors_all(),
    group_labels = titles, filename = NULL,
    font_scale = input$venn_font %||% 1, count_size = na_null(input$venn_count),
    group_label_size = na_null(input$venn_label), expand_limits = input$venn_ratio %||% 1.5,
    fill_alpha = input$venn_alpha %||% 0.45))$value, error = function(e) e)
  validate(need(!inherits(res, "error"), if (inherits(res, "error")) paste("Venn:", conditionMessage(res)) else ""))
  res
})

output$qc_grid_ui <- renderUI({
  if (!data_ready())
    return(helpText("Load the expression file and the sample info (or use the demo data): the QC plots appear automatically."))
  dl <- function(id) downloadButton(id, "PNG", class = "btn-xs")
  cell <- function(title, ui, id) div(class = "qc-cell", h5(title), ui, dl(id))
  tagList(
    fluidRow(column(6, cell("Number of proteins per sample", imageOutput("qc_bar", height = "auto"), "dl_qc_bar")),
             column(6, cell("Box plot", imageOutput("qc_box", height = "auto"), "dl_qc_box"))),
    fluidRow(column(6, cell("PCA", plotOutput("qc_pca", height = "330px"), "dl_qc_pca")),
             column(6, cell("Venn diagram", plotOutput("qc_venn", height = "330px"), "dl_qc_venn"))),
    uiOutput("venn_regions_ui"))
})

# ---- Proteins in the diagram (detected-proteins Venn only): every region with its proteins ----
venn_vr <- reactive({
  req(data_ready()); g <- input$venn_groups
  validate(need(length(g) >= 2, ""))
  titles <- vapply(seq_along(g), function(i) input[[paste0("vt_", i)]] %||% g[i], character(1))
  md <- qc_group_md()
  venn_regions(detected_sets(qc_matrix(), md, g, (input$venn_thr %||% 70) / 100, titles))
})
output$venn_regions_ui <- renderUI({
  vr <- tryCatch(venn_vr(), error = function(e) NULL); req(vr)
  tagList(hr(), h5("Proteins in the diagram"),
    helpText("Every protein is listed in exactly one region (the exact combination of groups it was detected in). ",
             "Click a row to see its proteins, then copy them from the box."),
    fluidRow(column(6, DT::DTOutput("venn_reg_dt")),
             column(6, uiOutput("venn_reg_box"))),
    downloadButton("dl_venn_xlsx", "Excel workbook (3 sheets)", class = "btn-sm"))
})
output$venn_reg_dt <- DT::renderDT({
  vr <- venn_vr()
  DT::datatable(vr$summary[, c("Region", "Count")], rownames = FALSE, selection = "single",
                class = "compact stripe", options = list(pageLength = 15, dom = "tip", ordering = FALSE))
})
output$venn_reg_box <- renderUI({
  vr <- venn_vr(); i <- input$venn_reg_dt_rows_selected
  if (is.null(i)) return(helpText("Select a region on the left."))
  r <- vr$summary$Region[i]; p <- vr$members[[r]]
  tagList(tags$b(sprintf("%s: %d protein(s)", r, length(p))),
          tags$textarea(readonly = NA, rows = 12, style = "width:100%; font-family:monospace;",
                        onclick = "this.select();", paste(p, collapse = "\n")))
})
output$dl_venn_xlsx <- downloadHandler(
  filename = function() "Venn_detected_proteins.xlsx",
  content = function(file) venn_workbook(venn_vr(), file))

output$qc_bar <- renderImage({ list(src = qc_bar_path(), contentType = "image/png", width = "100%") }, deleteFile = FALSE)
output$qc_box <- renderImage({ list(src = qc_box_path(), contentType = "image/png", width = "100%") }, deleteFile = FALSE)
output$qc_pca  <- renderPlot({ print(qc_pca_obj()) }, res = 100)
output$qc_venn <- renderPlot({ print(qc_venn_obj()) }, res = 100)

output$dl_qc_bar  <- downloadHandler("QC_barplot.png", function(file) file.copy(qc_bar_path(), file))
output$dl_qc_box  <- downloadHandler("QC_boxplot.png", function(file) file.copy(qc_box_path(), file))
output$dl_qc_pca  <- downloadHandler("QC_PCA.png", function(file) ggplot2::ggsave(file, qc_pca_obj(), width = 9, height = 7, dpi = 300))
output$dl_qc_venn <- downloadHandler("QC_venn.png", function(file)
  ggplot2::ggsave(file, qc_venn_obj(), width = 8 / sqrt(input$venn_ratio %||% 1.5),
                  height = 7 / sqrt(input$venn_ratio %||% 1.5), dpi = 300, bg = "white"))

# ---- rMd outlier detection (pmartR): computed once, the slider only re-flags ----
rmd_res <- reactiveVal(NULL)
rmd_valid <- reactive({
  r <- rmd_res(); if (is.null(r)) return(NULL)
  if (identical(r$key, tryCatch(stage_key(qc_matrix()), error = function(e) NA))) r else NULL
})
output$rmd_msg <- renderUI({
  if (!requireNamespace("pmartR", quietly = TRUE))
    return(div(class = "alert alert-danger", style = "padding:6px; font-size:12px;",
               "The 'pmartR' package is not installed on this server."))
  NULL
})
observeEvent(input$rmd_run, {
  m <- qc_matrix(); md <- qc_group_md()
  mets <- input$rmd_metrics
  if (!anyNA(m) && "Proportion_Missing" %in% mets) {
    mets <- setdiff(mets, "Proportion_Missing")
    showNotification("Proportion_Missing was left out: there are no missing values (it would be all zero).", type = "message")
  }
  if (length(mets) < 2) return(showNotification("Choose at least two metrics.", type = "warning"))
  withProgress(message = "Computing robust Mahalanobis distances...", value = 0.5, {
    r <- tryCatch(rmd_compute(m, md, metrics = mets), error = function(e) e)
  })
  if (inherits(r, "error")) return(showNotification(paste("rMd filter:", conditionMessage(r)), type = "error", duration = NULL))
  r$key <- stage_key(m); rmd_res(r)
})
rmd_flag <- reactive({ r <- rmd_valid(); req(r); rmd_flagged(r, rmd_pv()) })
rmd_pv <- reactive({ p <- input$rmd_p; if (is.null(p) || is.na(p) || p <= 0 || p > 1) 0.05 else p })
output$rmd_plot <- renderPlot({
  r <- rmd_valid()
  validate(need(!is.null(r), "Press 'Compute rMd' (uses the data stage selected for the QC plots)."))
  rmd_make_plot(r, rmd_pv(), group_colors_all())
}, res = 100)
output$rmd_table <- renderTable({
  r <- rmd_valid(); req(r)
  d <- rmd_flag_table(r, rmd_pv())
  validate(need(nrow(d) > 0, "No sample is flagged at this cutoff."))
  d
}, digits = 4)
output$rmd_remove_ui <- renderUI({
  r <- rmd_valid(); req(r)
  n <- length(rmd_flag())
  if (n == 0) return(helpText("No outliers at this cutoff."))
  actionButton("rmd_remove", sprintf("Remove detected outliers (%d)", n), class = "btn-warning")
})
observeEvent(input$rmd_remove, {
  fl <- rmd_flag(); req(length(fl))
  excluded(union(excluded(), fl))
  showNotification(sprintf("%d sample(s) removed from the analysis: %s", length(fl), paste(fl, collapse = ", ")), type = "warning")
})
output$rmd_box_ui <- renderUI({
  fl <- rmd_flag(); if (!length(fl)) return(NULL)
  show <- head(fl, 24)
  tagList(h5("Why were they flagged? Boxplots of each metric (flagged sample = X)"),
          if (length(fl) > length(show)) helpText(sprintf("Showing the first %d of %d.", length(show), length(fl))),
          fluidRow(lapply(seq_along(show), function(i)
            column(12, plotOutput(paste0("rmd_box_", i), height = "430px")))))
})
for (i in 1:24) local({
  k <- i
  output[[paste0("rmd_box_", k)]] <- renderPlot({
    r <- rmd_valid(); fl <- rmd_flag(); req(r, k <= length(fl))
    rmd_sample_plot(r, fl[k], rmd_pv())
  }, res = 100)
})

# ---- automatic tab switching on the Data page ----
observeEvent(list(input$expr_file, input$meta_file, input$load_demo), ignoreInit = TRUE, {
  both <- demo_on() || (!is.null(input$expr_file) && !is.null(input$meta_file))
  sel <- if (both) "QC plots" else if (!is.null(input$meta_file)) "Sample info (with Group)" else "Expression data"
  updateTabsetPanel(session, "data_tabs", selected = sel)
})


# ======================= editable tables (DT) =======================
# edits are kept in *_work together with the key of the file they belong to: a new file / setting
# automatically discards them
expr_work <- reactiveVal(NULL); meta_work <- reactiveVal(NULL)
expr_hdr <- reactiveVal(0); meta_hdr <- reactiveVal(0)
df_key <- function(d) paste(dim(d), paste(names(d), collapse = "|"), sum(nchar(as.matrix(d)), na.rm = TRUE), collapse = "_")

expr_raw <- reactive({
  b <- expr_base(); w <- expr_work()
  if (!is.null(w) && identical(attr(w, "src"), stage_key(b))) w else b
})
meta_raw <- reactive({
  b <- meta_base(); w <- meta_work()
  if (!is.null(w) && identical(attr(w, "src"), df_key(b))) w else b
})

# column headers are edited in place: click a header, type the new name, press Enter (or click away)
header_rename_js <- function(id) sprintf(
  "var cont = $(table.table().container());
   function arm() { cont.find('thead th').attr('contenteditable', 'true').css({cursor: 'text'}); }
   arm(); table.on('draw.dt', arm);
   cont.on('focusin', 'thead th', function() { $(this).data('old', $(this).text().trim()); });
   cont.on('focusout', 'thead th', function() {
     var old = $(this).data('old'), nw = $(this).text().trim();
     if (old === undefined) return;
     if (nw === '') { $(this).text(old); return; }
     if (nw !== old) Shiny.setInputValue('%s_colrename', {col: $(this).index(), name: nw, t: Date.now()}, {priority: 'event'});
   });
   cont.on('keydown', 'thead th', function(e) { if (e.key === 'Enter') { e.preventDefault(); this.blur(); } });", id)

expr_df <- function(m) { d <- data.frame(rownames(m), m, check.names = FALSE, stringsAsFactors = FALSE); names(d)[1] <- id_label(); d }

output$expr_dt <- DT::renderDT({
  expr_hdr(); expr_base()
  m <- isolate(expr_raw()); req(m)
  d <- expr_df(m)
  DT::datatable(d, rownames = FALSE, selection = "none", editable = list(target = "cell"),
                class = "compact stripe nowrap", callback = DT::JS(header_rename_js("expr_dt")),
                options = list(scrollX = TRUE, scrollY = "52vh", pageLength = 50, ordering = FALSE,   # headers are editable, so no sorting
                               lengthMenu = c(25, 50, 100, 500), deferRender = TRUE)) %>%
    DT::formatRound(columns = seq_len(ncol(d))[-1], digits = 3)
})
proxy_expr <- DT::dataTableProxy("expr_dt")

observeEvent(input$expr_dt_cell_edit, {
  info <- input$expr_dt_cell_edit; m <- expr_raw()
  for (k in seq_len(nrow(info))) {
    r <- info$row[k]; j <- info$col[k] + 1; v <- trimws(as.character(info$value[k]))
    if (j == 1) {                                             # gene / protein name
      if (!nzchar(v) || v %in% rownames(m)[-r]) showNotification("Names must be non-empty and unique.", type = "warning")
      else rownames(m)[r] <- v
    } else {
      x <- if (!nzchar(v) || tolower(v) %in% c("na", "nan")) NA_real_ else suppressWarnings(as.numeric(v))
      if (is.na(x) && nzchar(v) && !tolower(v) %in% c("na", "nan")) showNotification("Please enter a number (or leave empty for NA).", type = "warning")
      else m[r, j - 1] <- x
    }
  }
  attr(m, "src") <- stage_key(expr_base()); expr_work(m)
  DT::replaceData(proxy_expr, expr_df(m), resetPaging = FALSE, rownames = FALSE)
})
observeEvent(input$expr_dt_colrename, {
  a <- input$expr_dt_colrename; m <- expr_raw(); j <- a$col + 1; nm <- a$name
  if (j == 1) { id_label(nm) }
  else {
    cn <- colnames(m)
    if (nm %in% cn[-(j - 1)]) return(showNotification("That column name is already used.", type = "warning"))
    cn[j - 1] <- nm; colnames(m) <- cn
    attr(m, "src") <- stage_key(expr_base()); expr_work(m)
  }
  expr_hdr(expr_hdr() + 1)
})

output$meta_dt <- DT::renderDT({
  meta_hdr(); meta_base()
  d <- isolate(meta_raw()); req(d)
  DT::datatable(d, rownames = FALSE, selection = "none", editable = list(target = "cell"),
                class = "compact stripe nowrap", callback = DT::JS(header_rename_js("meta_dt")),
                options = list(scrollX = TRUE, scrollY = "45vh", pageLength = 50, ordering = FALSE, deferRender = TRUE))
})
proxy_meta <- DT::dataTableProxy("meta_dt")

observeEvent(input$meta_dt_cell_edit, {
  info <- input$meta_dt_cell_edit; d <- meta_raw()
  for (k in seq_len(nrow(info)))
    d[info$row[k], info$col[k] + 1] <- DT::coerceValue(info$value[k], d[info$row[k], info$col[k] + 1])
  attr(d, "src") <- df_key(meta_base()); meta_work(d)
  DT::replaceData(proxy_meta, d, resetPaging = FALSE, rownames = FALSE)
})
observeEvent(input$meta_dt_colrename, {
  a <- input$meta_dt_colrename; d <- meta_raw(); j <- a$col + 1; nm <- a$name
  if (nm %in% names(d)[-j]) return(showNotification("That column name is already used.", type = "warning"))
  names(d)[j] <- nm
  attr(d, "src") <- df_key(meta_base()); meta_work(d); meta_hdr(meta_hdr() + 1)
})

output$group_counts <- renderTable({
  md <- sample_meta()[common(), , drop = FALSE]
  as.data.frame(table(Group = md$Group), responseName = "Samples")
})
