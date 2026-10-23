# ======================= GSEA (fgsea) =======================
gsea_res <- reactiveVal(NULL)    # list: ontology -> list(contrast -> data.table)
gsea_dir <- reactiveVal(NULL)
gsea_padj_used <- reactiveVal(NA_real_)     # NA = the run kept every pathway

gsea_missing <- reactive({
  need <- c("fgsea", "data.table", "GO.db", "AnnotationDbi", "callr", input$ora_org)
  need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
})

output$gsea_contrasts_ui <- renderUI({
  validate(need(!is.null(de_res()), "Run the differential analysis first (tab 4)."))
  cn <- de_res()$single_names
  tagList(
    if (length(gsea_missing()))
      div(class = "alert alert-warning", "Missing R packages: ", strong(paste(gsea_missing(), collapse = ", "))),
    checkboxGroupInput("gsea_contrasts", "Contrasts", choices = cn, selected = cn))
})

observeEvent(input$run_gsea, {
  if (is.null(de_res())) return(showNotification("Run the differential analysis first.", type = "error"))
  if (length(gsea_missing()))
    return(showNotification(paste("Missing packages:", paste(gsea_missing(), collapse = ", ")), type = "error"))
  if (!length(input$gsea_contrasts) || !length(input$gsea_ont))
    return(showNotification("Choose at least one contrast and one ontology.", type = "error"))
  nn <- function(x) if (is.null(x) || is.na(x)) NULL else x
  nn_or <- function(x, default) if (is.null(x) || is.na(x)) default else x    # blank = no limit
  # fgsea ranks by the moderated t: a treat fit has none, so the contrasts fit is re-moderated with eBayes
  fit <- if (de_res()$test == "eBayes") de_res()$fit else
    limma::eBayes(de_res()$cfit, trend = isTRUE(input$de_trend), robust = isTRUE(input$de_robust))
  out <- file.path(session_tmp, paste0("GSEA_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  padj_run <- nn(input$gsea_padj)
  start_job("GSEA", "run_fgsea_app",
    list(fit = fit, orgdb = input$ora_org, onts = input$gsea_ont, coefs = input$gsea_contrasts,
         min_size = nn_or(input$gsea_min, 1), max_size = nn_or(input$gsea_max, 1e6), nperm = nn_or(input$gsea_nperm, 10000),
         padj_filter = nn(input$gsea_padj), seed = input$global_seed %||% 42, save_dir = out),
    on_done = function(res) {
      gsea_res(res); gsea_dir(out); gsea_padj_used(if (is.null(padj_run)) NA_real_ else padj_run)
      if (dir.exists(out)) addResourcePath(paste0("gseares_", session$token), out)
      updateTabsetPanel(session, "gsea_tabs", selected = "Result files")
    })
})

# ---- helpers ----
gsea_long <- reactive({
  r <- gsea_res(); req(r)
  do.call(rbind, lapply(names(r), function(o) do.call(rbind, lapply(names(r[[o]]), function(cn) {
    d <- as.data.frame(r[[o]][[cn]]); if (!nrow(d)) return(NULL)
    d$leadingEdge <- vapply(d$leadingEdge, paste, "", collapse = ";")
    d$ontology <- paste0("GO:", o); d$contrast <- cn; d
  }))))
})

gsea_res_url <- function(rel) {
  comps <- strsplit(rel, "/", fixed = TRUE)[[1]]
  paste0("gseares_", session$token, "/", paste(vapply(comps, utils::URLencode, "", reserved = TRUE), collapse = "/"))
}

# same logic as the ORA result-file table: contrast | method | ontology | file | terms | download
output$gsea_files_ui <- renderUI({
  validate(need(!is.null(gsea_res()) && !is.null(gsea_dir()), "No GSEA results yet."))
  d <- gsea_dir()
  files <- list.files(d, pattern = "\\.csv$", recursive = TRUE)
  validate(need(length(files) > 0, "No result files were written."))
  trs <- lapply(files, function(f) {
    df <- safe_csv(file.path(d, f))
    total <- if (is.null(df)) 0L else nrow(df)
    n <- if (is.null(df) || !"padj" %in% colnames(df)) 0L else sum(df$padj < 0.05, na.rm = TRUE)
    tags$tr(tags$td(sub("_fgsea\\.csv$", "", basename(f))), tags$td("fgsea"),
            tags$td(sub("^GO_", "GO:", dirname(f))), tags$td(basename(f)), tags$td(total), tags$td(n),
            tags$td(tags$a(icon("download"), href = gsea_res_url(f), download = basename(f),
                           class = "btn btn-xs btn-default", title = "Download this CSV")))
  })
  tags$table(class = "table table-condensed table-hover",
    tags$thead(tags$tr(tags$th("Gene set (contrast)"), tags$th("Method"), tags$th("Ontology"),
                       tags$th("Files"), tags$th("Total terms"), tags$th("Terms (padj < 0.05)"), tags$th(""))),
    tags$tbody(trs))
})

output$dl_gsea_csvs <- downloadHandler("GSEA_all_CSVs.zip", function(file) {
  d <- gsea_dir(); req(d, dir.exists(d))
  f <- list.files(d, pattern = "\\.csv$", recursive = TRUE); req(length(f) > 0)
  zip::zip(file, files = f, root = d)
})

# ---- contrast / ontology pickers shared by the plots ----
gsea_pick_ui <- function(prefix) {
  r <- gsea_res(); validate(need(!is.null(r), "No GSEA results yet."))
  cn <- unique(unlist(lapply(r, names)))
  fluidRow(column(3, selectInput(paste0(prefix, "_ont"), "Ontology", setNames(names(r), paste0("GO:", names(r))))),
           column(5, selectInput(paste0(prefix, "_cn"), "Contrast", cn)),
           column(4, numericInput(paste0(prefix, "_padj"), "padj cut-off (blank = none)",
                                  gsea_padj_used(), min = 0, max = 1, step = 0.01)))
}
gsea_get <- function(prefix) {
  r <- gsea_res(); req(r, input[[paste0(prefix, "_ont")]], input[[paste0(prefix, "_cn")]])
  d <- r[[input[[paste0(prefix, "_ont")]]]][[input[[paste0(prefix, "_cn")]]]]
  validate(need(!is.null(d) && nrow(d) > 0, "No pathways for this selection."))
  d <- as.data.frame(d); d$leadingEdge <- vapply(d$leadingEdge, paste, "", collapse = ";")
  cut <- input[[paste0(prefix, "_padj")]]
  if (is.null(cut) || is.na(cut)) return(d)             # blank: every term is kept
  d[!is.na(d$padj) & d$padj <= cut, , drop = FALSE]
}

# ---- NES bar plot: top up / down pathways of one contrast ----
output$gsea_nes_ctl <- renderUI(tagList(gsea_pick_ui("gn"), numericInput("gn_top", "Pathways per direction", 10, min = 1)))
gsea_nes_obj <- reactive({
  d <- gsea_get("gn"); validate(need(nrow(d) > 0, "No significant pathways at this cut-off."))
  k <- input$gn_top %||% 10
  d <- rbind(head(d[d$NES > 0, ][order(-d$NES[d$NES > 0]), ], k), head(d[d$NES < 0, ][order(d$NES[d$NES < 0]), ], k))
  d$pathway <- factor(d$pathway, levels = d$pathway[order(d$NES)])
  ggplot(d, aes(NES, pathway, fill = -log10(padj))) + geom_col() +
    scale_fill_gradient(low = "grey80", high = "red", name = "-log10(padj)") +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 50)) +
    geom_vline(xintercept = 0) + labs(y = NULL, title = paste("GSEA:", input$gn_cn)) + theme_bw(base_size = 13)
})
# height from the wrapped label lines (as in the ORA lollipop): every pathway row is as tall as its longest label
lines_of <- function(x, width = 50) lengths(strsplit(stringr::str_wrap(as.character(x), width), "\n"))
gsea_nes_h <- function() {
  d <- tryCatch(gsea_nes_obj()$data, error = function(e) NULL)
  if (is.null(d)) return(400)
  slot <- 19 * max(lines_of(d$pathway)) + 10          # px per pathway row
  min(9000, max(380, slot * nrow(d) + 170))
}
output$gsea_nes_plot <- renderPlot(gsea_nes_obj(), res = 100, height = function() gsea_nes_h())
output$dl_gsea_nes <- downloadHandler(
  filename = function() paste0("GSEA_NES_", gsub("[^A-Za-z0-9]+", "_", input$gn_cn %||% "contrast"), "_GO_", input$gn_ont %||% "", ".png"),
  content = function(file) ggplot2::ggsave(file, gsea_nes_obj(), width = 9, height = gsea_nes_h() / 100, dpi = 300))

# ---- compare contrasts: NES dot plot of the top pathways of every contrast ----
output$gsea_cmp_ctl <- renderUI({
  r <- gsea_res(); validate(need(!is.null(r), "No GSEA results yet."))
  fluidRow(column(3, selectInput("gc_ont", "Ontology", setNames(names(r), paste0("GO:", names(r))))),
           column(3, numericInput("gc_top", "Top pathways per contrast", 5, min = 1)),
           column(3, numericInput("gc_padj", "padj cut-off (blank = none)", gsea_padj_used(), min = 0, max = 1, step = 0.01)))
})
gsea_cmp_obj <- reactive({
  r <- gsea_res(); req(r, input$gc_ont)
  lst <- r[[input$gc_ont]]; cut <- input$gc_padj
  parts <- lapply(names(lst), function(cn) {
    d <- as.data.frame(lst[[cn]]); if (!is.null(cut) && !is.na(cut)) d <- d[!is.na(d$padj) & d$padj <= cut, , drop = FALSE]
    if (!nrow(d)) return(NULL)
    d$contrast <- cn; head(d[order(d$padj), ], input$gc_top %||% 5)
  })
  top <- unique(unlist(lapply(Filter(Negate(is.null), parts), function(d) d$pathway)))
  validate(need(length(top) > 0, "No significant pathways at this cut-off."))
  all <- do.call(rbind, lapply(names(lst), function(cn) { d <- as.data.frame(lst[[cn]]); d$contrast <- cn; d }))
  d <- all[all$pathway %in% top, ]
  d$contrast <- factor(d$contrast, levels = names(lst))
  ggplot(d, aes(contrast, pathway, size = -log10(padj), colour = NES)) + geom_point() +
    scale_colour_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, 50)) +
    labs(x = NULL, y = NULL) + theme_bw(base_size = 13) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
})
gsea_cmp_h <- function() tryCatch({
  pw <- unique(as.character(gsea_cmp_obj()$data$pathway))
  min(9000, max(380, (19 * max(lines_of(pw)) + 10) * length(pw) + 170))
}, error = function(e) 400)
output$gsea_cmp_plot <- renderPlot(gsea_cmp_obj(), res = 100, height = function() gsea_cmp_h())
output$dl_gsea_cmp <- downloadHandler(
  filename = function() paste0("GSEA_compare_contrasts_GO_", input$gc_ont %||% "", ".png"),
  content = function(file) ggplot2::ggsave(file, gsea_cmp_obj(), width = 9, height = gsea_cmp_h() / 100, dpi = 300, limitsize = FALSE))
