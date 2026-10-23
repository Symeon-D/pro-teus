# ======================= 6. ORA =======================
# packages needed by the chosen methods
ora_missing <- reactive({
  need <- c("callr", "stringr", "httr")
  if ("gprofiler" %in% input$ora_method) need <- c(need, "gprofiler2", "openxlsx", "dplyr")
  if ("genecodis" %in% input$ora_method) need <- c(need, "jsonlite")
  if ("clusterprofiler" %in% input$ora_method) need <- c(need, "clusterProfiler", input$ora_org)
  if (isTRUE(input$ora_simplify)) need <- c(need, "rrvgo", input$ora_org)
  if (isTRUE(input$ora_cp_simplify)) need <- c(need, "clusterProfiler")
  need <- unique(need)
  need[!vapply(need, requireNamespace, logical(1), quietly = TRUE)]
})
output$ora_avail <- renderUI({
  miss <- ora_missing()
  if (!length(miss)) return(NULL)
  div(class = "alert alert-warning",
      "Missing R packages for the selected options: ", strong(paste(miss, collapse = ", ")),
      ". Install them (CRAN / Bioconductor) or untick the related options.")
})

# ---- gene sets: DE subsets + the clusters of the current heatmap (c1, c2, ...) ----
hm_clusters <- reactive({
  tt <- tryCatch(result_table(), error = function(e) NULL)
  if (is.null(tt)) return(list())
  s <- split(rownames(tt), tt$Cluster)
  if (length(s) < 2) return(list())            # a single cluster = all proteins, nothing new
  setNames(s, paste0("c", names(s)))
})
ora_all_sets <- reactive(c(all_subsets(), hm_clusters()))

output$ora_sets_ui <- renderUI({
  s <- ora_all_sets()
  selectizeInput("ora_sets", "Gene set(s)", multiple = TRUE,
                 selected = intersect(isolate(input$ora_sets), names(s)),
                 choices = setNames(names(s), paste0(names(s), " (", lengths(s), ")")))
})

# ---- run in a background process ----
ora_dir <- reactiveVal(NULL)
ora_sets_run <- reactiveVal(character(0))

observeEvent(input$run_ora, {
  s <- ora_all_sets()[input$ora_sets %||% character(0)]
  if (length(ora_missing()))
    return(showNotification(paste("Missing packages:", paste(ora_missing(), collapse = ", ")), type = "error"))
  if (!length(s)) return(showNotification("Select at least one gene set.", type = "error"))
  if (!length(input$ora_method) || !length(input$ora_ont))
    return(showNotification("Choose at least one method and one ontology.", type = "error"))
  nn <- function(x) if (is.null(x) || is.na(x)) NULL else x
  out <- file.path(session_tmp, paste0("ORA_", format(Sys.time(), "%Y%m%d_%H%M%S")))
  bg <- if (identical(input$ora_bg, "all")) rownames(final_expr()) else NULL
  run_names <- names(s)
  ok <- start_job("Enrichment analysis", "run_ora_analysis_pipeline",
    list(subsets_list = s, background = bg, top_n = nn(input$ora_top_n), output_dir = out,
         min_genes = input$ora_min_genes, min_term_size = nn(input$ora_min_ts),
         max_term_size = nn(input$ora_max_ts),
         correction_method = input$ora_corr %||% "gSCS", ontologies = input$ora_ont,
         ora_method = input$ora_method, simplify = isTRUE(input$ora_simplify),
         cp_simplify = isTRUE(input$ora_cp_simplify), cp_cutoff = input$ora_cp_cutoff %||% 0.7,
         organism = input$ora_org, keyType = "SYMBOL",
         rrvgo_threshold = input$ora_rrvgo %||% 0.7, save_df = TRUE,
         show_cluster_label = TRUE,
         cp_correction_method = input$ora_cp_corr %||% "BH"),
    on_done = function(res) {
      ora_sets_run(run_names)
      if (dir.exists(out)) addResourcePath(paste0("orares_", session$token), out)
      ora_dir(out)
      updateTabsetPanel(session, "ora_tabs", selected = "Result files")
    })
  if (isTRUE(ok)) ora_dir(NULL)
})

# ---- safe file helpers (never fail on a missing / unreadable file) ----
safe_csv <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  tryCatch(read.csv(path, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
}
res_url <- function(rel) {
  comps <- strsplit(rel, "/", fixed = TRUE)[[1]]
  paste0("orares_", session$token, "/", paste(vapply(comps, utils::URLencode, "", reserved = TRUE), collapse = "/"))
}
ont_of <- function(folder) sub("_", ":", folder)

ora_files <- reactive({
  d <- ora_dir(); req(d, dir.exists(d))
  list.files(d, pattern = "\\.csv$", recursive = TRUE)
})

# <gene set>/<tool>/<ontology>/<file>.csv
#   main = un-simplified-or-clusterProfiler-simplified enrichment table; kind = plain | cp | rrvgo (full table)
ora_index <- reactive({
  f <- ora_files(); parts <- strsplit(f, "/", fixed = TRUE)
  ok <- lengths(parts) >= 4
  f <- f[ok]; parts <- parts[ok]; b <- basename(f)
  data.frame(file = f, subset = vapply(parts, `[`, "", 1), tool = vapply(parts, `[`, "", 2),
             ont = vapply(parts, `[`, "", 3), main = !grepl("rrvgo|_full_|cpSimplified", b),
             kind = ifelse(grepl("rrvgo_full", b), "rrvgo", ifelse(grepl("cpSimplified", b), "cp", "plain")),
             stringsAsFactors = FALSE)
})

# number of terms in a result file: for simplified files this is the FILTERED number
n_terms_of <- function(df, kind) {
  if (is.null(df)) return(0L)
  if (kind == "rrvgo" && "cluster" %in% colnames(df)) {      # rrvgo: representatives (one per cluster + unclustered terms)
    cl <- as.character(df$cluster); un <- cl == "Ungrouped"
    return(as.integer(sum(un) + length(unique(paste(df$group_name[!un], cl[!un])))))
  }
  nrow(df)
}
FILTER_LABEL <- c(plain = "none", rrvgo = "rrvgo", cp = "cp:simplify")

# one row per CSV; gene sets without any result file still get a row (file "-", 0 terms)
ora_rows <- reactive({
  idx <- ora_index(); d <- ora_dir()
  do.call(rbind, lapply(ora_sets_run(), function(s) {
    ix <- idx[idx$subset == s, , drop = FALSE]
    if (!nrow(ix)) return(data.frame(set = s, tool = "-", ont = "-", filt = "-", rel = NA_character_, n = 0L, stringsAsFactors = FALSE))
    # same gene set: each unfiltered table is followed by its filtered versions (rrvgo, cp:simplify)
    base <- basename(ix$file)
    b2 <- ifelse(startsWith(base, paste0(s, "_")), substring(base, nchar(s) + 2), base)     # drop the gene-set prefix
    key <- gsub("_cpSimplified|rrvgo_full_", "", b2)
    ord <- order(ix$tool, ix$ont, key, match(ix$kind, c("plain", "rrvgo", "cp")))
    ix <- ix[ord, , drop = FALSE]
    data.frame(set = s, tool = ix$tool, ont = ont_of(ix$ont), filt = unname(FILTER_LABEL[ix$kind]), rel = ix$file,
               n = vapply(seq_len(nrow(ix)), function(k) n_terms_of(safe_csv(file.path(d, ix$file[k])), ix$kind[k]), integer(1)),
               stringsAsFactors = FALSE)
  }))
})

# ---- "Open in ORA-cle": anonymous Gists are no longer accepted by GitHub, so a token is needed ----
gist_token <- function() { t <- Sys.getenv("GITHUB_PAT"); if (!nzchar(t)) t <- Sys.getenv("GITHUB_TOKEN"); t }

create_gist <- function(txt, filename, description) {
  tk <- gist_token(); if (!nzchar(tk)) return(NULL)
  resp <- tryCatch(httr::POST("https://api.github.com/gists",
      httr::add_headers(Authorization = paste("Bearer", tk), Accept = "application/vnd.github+json"),
      body = list(description = description, public = FALSE,
                  files = setNames(list(list(content = txt)), filename)),
      encode = "json", httr::timeout(30)), error = function(e) NULL)
  if (is.null(resp) || httr::status_code(resp) != 201) return(NULL)
  httr::content(resp)$id
}

observeEvent(input$ora_open, {
  rows <- ora_rows(); r <- rows[input$ora_open$i, ]
  path <- if (is.na(r$rel)) NULL else file.path(ora_dir(), r$rel)
  if (is.null(path) || !file.exists(path)) {
    session$sendCustomMessage("ora_fail", "")
    return(showNotification("No result file for this gene set.", type = "warning"))
  }
  txt <- paste(readLines(path, warn = FALSE, encoding = "UTF-8"), collapse = "\n")
  id <- create_gist(txt, basename(path), paste("ORA-cle export:", r$set))
  url <- if (!is.null(id)) paste0(ORA_SITE, "#gist=", id) else {
    enc <- if (requireNamespace("curl", quietly = TRUE)) curl::curl_escape(enc2utf8(txt)) else utils::URLencode(enc2utf8(txt), reserved = TRUE)
    u <- paste0(ORA_SITE, "#data=", base64enc::base64encode(charToRaw(enc)))
    if (nchar(u) > 2e6) NULL else u
  }
  if (is.null(url)) {
    session$sendCustomMessage("ora_fail", "")
    showNotification("Could not create the link: set a GITHUB_PAT environment variable (token with 'gist' scope) for large tables.",
                     type = "error", duration = 10)
  } else {
    if (is.null(id) && !nzchar(gist_token()))
      showNotification("No GITHUB_PAT set: used a data URL (small tables only).", type = "message")
    session$sendCustomMessage("ora_go", url)
  }
})

output$ora_files_ui <- renderUI({
  validate(need(!is.null(ora_dir()), "No results yet."))
  rows <- ora_rows(); validate(need(nrow(rows) > 0, "No gene sets were analysed."))
  trs <- lapply(seq_len(nrow(rows)), function(i) {
    r <- rows[i, ]; has <- !is.na(r$rel)
    tags$tr(
      tags$td(r$set), tags$td(r$tool), tags$td(r$ont), tags$td(r$filt),
      tags$td(if (has) basename(r$rel) else "-"),
      tags$td(r$n),
      tags$td(if (has) tags$button(type = "button", class = "btn btn-xs btn-primary",
                                   onclick = sprintf("oraOpen(%d)", i), "Open in ORA-cle")
              else tags$button(type = "button", class = "btn btn-xs btn-default", disabled = NA, "Open in ORA-cle")),
      tags$td(if (has) tags$a(icon("download"), href = res_url(r$rel), download = basename(r$rel),
                              class = "btn btn-xs btn-default", title = "Download this CSV")))
  })
  tags$table(class = "table table-condensed table-hover",
    tags$thead(tags$tr(tags$th("Gene set"), tags$th("Method"), tags$th("Ontology"), tags$th("filtered by"), tags$th("Files"),
                       tags$th("Terms"), tags$th(""), tags$th(""))),
    tags$tbody(trs))
})

output$dl_ora_csvs <- downloadHandler("ORA_all_CSVs.zip", function(file) {
  f <- ora_files(); req(length(f) > 0)
  zip::zip(file, files = f, root = ora_dir())
})
output$dl_ora_zip <- downloadHandler("ORA_results.zip", function(file) {
  f <- list.files(ora_dir(), recursive = TRUE); req(length(f) > 0)
  zip::zip(file, files = f, root = ora_dir())
})

read_ora_csv <- function(rel, subset_name = NULL) {
  df <- safe_csv(file.path(ora_dir(), rel))
  validate(need(!is.null(df), "Result file not found."))
  validate(need(all(c("term_name", "-log10(padj)") %in% colnames(df)),
                "This result table lacks the columns needed for the plot (term_name, -log10(padj))."))
  if (!"group_name" %in% colnames(df)) df$group_name <- subset_name %||% "set"
  df
}

# ---- quick view: gene set + method + ontology, capped number of pathways ----
output$ora_quick_ui <- renderUI({
  validate(need(!is.null(ora_dir()), "No results yet."))
  sets <- ora_sets_run(); idx <- ora_index(); idx <- idx[idx$main, ]
  with_res <- intersect(sets, idx$subset)
  tagList(
    fluidRow(
      column(3, selectInput("qv_set", "Gene set", sets, selected = if (length(with_res)) with_res[1] else sets[1])),
      column(3, selectInput("qv_tool", "Method", unique(idx$tool))),
      column(3, selectInput("qv_ont", "Ontology", unique(idx$ont))),
      column(3, numericInput("ora_quick_n", "Pathways shown", 20, min = 1, step = 1))),
    fluidRow(
      column(4, radioButtons("ora_quick_x", "Lollipop x-axis", c("-log10(padj)" = "padj", "Fold enrichment" = "fold"), inline = TRUE)),
      column(4, radioButtons("ora_quick_col", "Colour bar", c("-log10(padj)" = "padj", "Fold enrichment" = "fold"), selected = "fold", inline = TRUE)),
      column(4, radioButtons("ora_quick_size", "Dot size", c("Gene ratio (%)" = "percentage", "Gene count" = "numeric"), inline = TRUE))))
})
quick_file <- reactive({
  req(input$qv_set, input$qv_tool, input$qv_ont)
  idx <- ora_index()
  idx <- idx[idx$main & idx$subset == input$qv_set & idx$tool == input$qv_tool & idx$ont == input$qv_ont, , drop = FALSE]
  if (nrow(idx)) idx[1, ] else NULL
})
# the maximum follows the number of terms in the selected result file
quick_avail <- reactive({
  r <- quick_file(); if (is.null(r)) return(0L)
  d <- safe_csv(file.path(ora_dir(), r$file)); if (is.null(d)) 0L else nrow(d)
})
observe({
  n <- quick_avail()
  if (n > 0) updateNumericInput(session, "ora_quick_n", max = n, value = min(isolate(input$ora_quick_n) %||% 20, n))
})
quick_n <- reactive({ n <- input$ora_quick_n %||% 20; if (is.na(n) || n < 1) 20 else round(n) })
quick_plot <- reactive({
  req(input$qv_set, input$qv_tool, input$qv_ont)
  idx <- ora_index()
  idx <- idx[idx$main & idx$subset == input$qv_set & idx$tool == input$qv_tool & idx$ont == input$qv_ont, , drop = FALSE]
  validate(need(nrow(idx) > 0, "No significant terms found"))
  row <- idx[1, ]
  df <- read_ora_csv(row$file, row$subset)
  validate(need(nrow(df) > 0, "No significant terms found"))
  xcol <- if (identical(input$ora_quick_x, "fold") && "log2(Fold_Enrichment)" %in% colnames(df)) "log2(Fold_Enrichment)" else "-log10(padj)"
  df <- head(df[order(-df[[xcol]]), ], quick_n())
  my_custom_gpro_dotplot(df, ont = ont_of(row$ont), top_n = quick_n(),
                         term_size_legend = input$ora_quick_size %||% "percentage",
                         sub_cluster = NULL, lolipop = TRUE, multi_cluster = FALSE,
                         lollipop_x = input$ora_quick_x %||% "padj", color_metric = input$ora_quick_col %||% "fold")
})
quick_h <- function() min(3000, max(360, 24 * quick_n() + 150))
output$ora_quick_plot <- renderPlot(quick_plot(), res = 100, height = function() quick_h())
output$dl_ora_quick <- downloadHandler(
  filename = function() paste0("ORA_lollipop_", tools::file_path_sans_ext(basename(quick_file()$file %||% "result")), ".png"),
  content = function(file) ggplot2::ggsave(file, quick_plot(), width = 9, height = quick_h() / 100, dpi = 300))

# ---- multi-cluster plot: any method, empty gene sets stay on the axis, per-cluster markers ----
output$ora_mc_ui <- renderUI({
  idx <- ora_index(); idx <- idx[idx$main, ]
  validate(need(!is.null(ora_dir()) && nrow(idx) > 0, "Run an enrichment first."))
  fluidRow(
    column(2, selectInput("mc_tool", "Method", unique(idx$tool))),
    column(2, selectInput("mc_ont", "Ontology", unique(idx$ont))),
    column(2, numericInput("mc_top", "Top terms per cluster", 10, min = 1)),
    column(3, radioButtons("mc_size", "Dot size", c("Gene ratio (%)" = "percentage", "Gene count" = "numeric"), inline = TRUE)),
    column(3, sliderInput("mc_angle", "x-axis label angle", -90, 90, 0, step = 15)))
  # second row: colour bar
})
output$mc_sets_ui <- renderUI({
  req(ora_dir())
  fluidRow(
    column(9, selectizeInput("mc_sets", "Clusters / gene sets to plot", multiple = TRUE, width = "100%",
                             choices = ora_sets_run(), selected = ora_sets_run())),
    column(3, radioButtons("mc_col", "Colour bar", c("-log10(padj)" = "padj", "Fold enrichment" = "fold"), inline = TRUE)))
})
output$mc_shapes_ui <- renderUI({
  sets <- input$mc_sets; req(length(sets) > 0)
  fluidRow(lapply(seq_along(sets), function(i)
    column(2, selectInput(paste0("mcs_", i), paste("Marker:", sets[i]),
                          c("Dot" = "21", "Up arrow" = "24", "Down arrow" = "25")))))
})
# max of "top terms" = largest table among the selected gene sets
observe({
  req(input$mc_sets, input$mc_tool, input$mc_ont)
  idx <- ora_index(); idx <- idx[idx$main & idx$tool == input$mc_tool & idx$ont == input$mc_ont &
                                   idx$subset %in% input$mc_sets, , drop = FALSE]
  n <- suppressWarnings(max(0L, vapply(idx$file, function(f) { d <- safe_csv(file.path(ora_dir(), f)); if (is.null(d)) 0L else nrow(d) }, integer(1))))
  if (n > 0) updateNumericInput(session, "mc_top", max = n, value = min(isolate(input$mc_top) %||% 10, n))
})
mc_data <- reactive({
  req(input$mc_sets, input$mc_tool, input$mc_ont, input$mc_top)
  idx <- ora_index(); idx <- idx[idx$main & idx$tool == input$mc_tool & idx$ont == input$mc_ont &
                                   idx$subset %in% input$mc_sets, , drop = FALSE]
  idx <- idx[!duplicated(idx$subset), ]
  parts <- lapply(seq_len(nrow(idx)), function(i) {
    d <- read_ora_csv(idx$file[i], idx$subset[i]); if (!nrow(d)) return(NULL)
    d$group_name <- idx$subset[i]
    head(d[order(-d$`-log10(padj)`), ], input$mc_top)
  })
  parts <- Filter(Negate(is.null), parts)
  validate(need(length(parts) > 0, "No significant terms in any of the selected gene sets."))
  df <- do.call(rbind, parts)
  df$group_name <- factor(df$group_name, levels = input$mc_sets)     # empty sets keep their place on the axis
  df
})
mc_plot <- reactive({
  df <- mc_data(); sets <- input$mc_sets
  shp <- vapply(seq_along(sets), function(i) as.numeric(input[[paste0("mcs_", i)]] %||% 21), numeric(1))
  my_custom_gpro_dotplot(df, ont = ont_of(input$mc_ont), top_n = input$mc_top,
                         term_size_legend = input$mc_size %||% "percentage", sub_cluster = NULL,
                         lolipop = FALSE, multi_cluster = TRUE,
                         shape_map = setNames(shp, sets), x_angle = input$mc_angle %||% 0,
                         color_metric = input$mc_col %||% "padj")
})
mc_height <- function() max(380, 26 * length(unique(mc_data()$term_name)) + 160)
output$ora_mc_plot <- renderPlot(mc_plot(), res = 100, height = function() tryCatch(mc_height(), error = function(e) 400))
output$dl_ora_mc <- downloadHandler(
  filename = function() paste0("ORA_multicluster_", input$mc_tool, "_", input$mc_ont, ".png"),
  content = function(file) {
    # same pixel size as the plot on screen (res = 100) -> identical proportions in the PNG
    w_px <- session$clientData$output_ora_mc_plot_width %||% 900
    ggplot2::ggsave(file, mc_plot(), width = w_px / 100, height = mc_height() / 100, dpi = 300, limitsize = FALSE)
  })

# ---- simplified results: rrvgo (lollipop of representatives + tree map + scatter) and clusterProfiler::simplify ----
ora_png_index <- reactive({
  d <- ora_dir(); req(d, dir.exists(d))
  f <- list.files(d, pattern = "\\.png$", recursive = TRUE)
  parts <- strsplit(f, "/", fixed = TRUE); ok <- lengths(parts) >= 4
  f <- f[ok]; parts <- parts[ok]; b <- basename(f)
  data.frame(file = f, subset = vapply(parts, `[`, "", 1), tool = vapply(parts, `[`, "", 2),
             ont = vapply(parts, `[`, "", 3),
             type = ifelse(grepl("rrvgo_treemap", b), "tree", ifelse(grepl("rrvgo_scatter", b), "scatter", "other")),
             stringsAsFactors = FALSE)
})

output$simp_ui <- renderUI({
  validate(need(!is.null(ora_dir()), "No results yet."))
  idx <- ora_index(); idx <- idx[idx$kind %in% c("rrvgo", "cp"), , drop = FALSE]
  validate(need(nrow(idx) > 0, "No simplified output: tick 'Simplify with rrvgo' and/or 'clusterProfiler::simplify()' and run the enrichment."))
  tagList(
    fluidRow(
      column(3, selectInput("sm_set", "Gene set", unique(idx$subset))),
      column(3, selectInput("sm_tool", "Method", unique(idx$tool))),
      column(3, selectInput("sm_ont", "Ontology", unique(idx$ont))),
      column(3, selectInput("sm_kind", "Simplification",
                            c("rrvgo" = "rrvgo", "clusterProfiler::simplify" = "cp")[c("rrvgo", "cp") %in% idx$kind]))),
    tabsetPanel(
      tabPanel("Lollipop plot", br(),
        fluidRow(
          column(3, numericInput("sm_n", "Pathways shown", 20, min = 1)),
          column(3, radioButtons("sm_x", "x-axis", c("-log10(padj)" = "padj", "Fold enrichment" = "fold"), inline = TRUE)),
          column(3, radioButtons("sm_col", "Colour bar", c("-log10(padj)" = "padj", "Fold enrichment" = "fold"), selected = "fold", inline = TRUE)),
          column(3, radioButtons("sm_size", "Dot size", c("Gene ratio (%)" = "percentage", "Gene count" = "numeric"), inline = TRUE))),
        plotOutput("sm_lollipop", height = "auto"), downloadButton("dl_sm_lollipop", "Download PNG")),
      tabPanel("Tree map", br(), uiOutput("sm_tree")),
      tabPanel("Scatter plot", br(), uiOutput("sm_scatter"))))
})
sm_df <- reactive({
  req(input$sm_set, input$sm_tool, input$sm_ont, input$sm_kind)
  idx <- ora_index()
  idx <- idx[idx$kind == input$sm_kind & idx$subset == input$sm_set & idx$tool == input$sm_tool & idx$ont == input$sm_ont, , drop = FALSE]
  validate(need(nrow(idx) > 0, "Not available for this selection."))
  df <- read_ora_csv(idx$file[1], idx$subset[1])
  if (input$sm_kind == "rrvgo") {
    # representatives: best term of every rrvgo cluster (unclustered terms are kept as they are)
    df <- df[order(df$p_value), ]
    if ("cluster" %in% colnames(df))
      df <- df[as.character(df$cluster) == "Ungrouped" | !duplicated(paste(df$group_name, df$cluster)), ]
  }
  validate(need(nrow(df) > 0, "No significant terms found"))
  attr(df, "ont") <- idx$ont[1]; attr(df, "file") <- idx$file[1]
  df
})
observe({
  n <- tryCatch(nrow(sm_df()), error = function(e) 0L)
  if (n > 0) updateNumericInput(session, "sm_n", max = n, value = min(isolate(input$sm_n) %||% 20, n))
})
sm_n <- reactive({ n <- input$sm_n %||% 20; if (is.na(n) || n < 1) 20 else round(n) })
sm_plot <- reactive({
  df <- sm_df()
  xcol <- if (identical(input$sm_x, "fold") && "log2(Fold_Enrichment)" %in% colnames(df)) "log2(Fold_Enrichment)" else "-log10(padj)"
  df <- head(df[order(-df[[xcol]]), ], sm_n())
  my_custom_gpro_dotplot(df, ont = ont_of(attr(sm_df(), "ont")), top_n = sm_n(), term_size_legend = input$sm_size %||% "percentage",
                         sub_cluster = NULL, lolipop = TRUE, multi_cluster = FALSE, lollipop_x = input$sm_x %||% "padj",
                         color_metric = input$sm_col %||% "fold")
})
sm_h <- function() min(3000, max(360, 24 * sm_n() + 150))
output$sm_lollipop <- renderPlot(sm_plot(), res = 100, height = function() sm_h())
output$dl_sm_lollipop <- downloadHandler(
  filename = function() paste0(if (identical(input$sm_kind, "cp")) "ORA_cp_simplify_lollipop_" else "ORA_rrvgo_lollipop_",
                               tools::file_path_sans_ext(basename(attr(sm_df(), "file") %||% "result")), ".png"),
  content = function(file) ggplot2::ggsave(file, sm_plot(), width = 9, height = sm_h() / 100, dpi = 300))

# tree map / scatter: PNGs drawn by the pipeline, shown to fit the window (no vertical scrollbar)
sm_images <- function(type) {
  idx <- ora_png_index(); req(input$sm_set, input$sm_tool, input$sm_ont)
  idx <- idx[idx$type == type & idx$subset == input$sm_set & idx$tool == input$sm_tool & idx$ont == input$sm_ont, , drop = FALSE]
  if (!nrow(idx) || !identical(input$sm_kind, "rrvgo")) return(helpText("Available for rrvgo simplification only."))
  tagList(lapply(seq_len(nrow(idx)), function(i) {
    u <- res_url(idx$file[i])
    div(style = "text-align:center;",
        tags$img(src = u, style = "max-width:100%; max-height:calc(100vh - 330px); object-fit:contain; border:1px solid #ddd;"),
        br(), tags$a(icon("download"), " Download PNG", href = u, download = basename(idx$file[i]), class = "btn btn-sm btn-default"))
  }))
}
output$sm_tree    <- renderUI(sm_images("tree"))
output$sm_scatter <- renderUI(sm_images("scatter"))

# ---- the main panel shows one tab per selected module (ORA and/or GSEA) ----
output$enrich_tabs <- renderUI({
  mods <- input$enrich_modules %||% "ora"
  panels <- list(ora = tabPanel("ORA", br(), ora_main_panel()), gsea = tabPanel("GSEA", br(), gsea_main_panel()))
  do.call(tabsetPanel, c(list(id = "enrich_main"), unname(panels[intersect(c("ora", "gsea"), mods)])))
})
