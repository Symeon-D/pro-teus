  # ======================= 2. DIFFERENTIAL =======================
  # built with its choices directly so that clicking the box always lists the variables
  output$de_covs_ui <- renderUI({
    app_gen(); req(other_cols())
    selectizeInput("de_covs", "Additional covariates (e.g. Batch)", choices = other_cols(),
                   multiple = TRUE, selected = isolate(input$de_covs))
  })

  de_meta <- reactive({
    cs <- common()
    md <- sample_meta()[cs, , drop = FALSE]
    md$Group_de <- factor(make.names(md$Group), levels = make.names(unique(md$Group)))
    md
  })
  de_levels <- reactive(levels(de_meta()$Group_de))

  de_design <- reactive({
    d <- tryCatch(build_design(de_meta(), input$de_covs %||% character(0)),
                  error = function(e) NULL)
    validate(need(!is.null(d), "Could not build the design matrix. Check the covariates."))
    d
  })

  output$design_txt <- renderPrint({
    d <- de_design()
    cat(attr(d, "formula_text"), "\n\n")
    cat("Design columns:", paste(colnames(d), collapse = ", "), "\n")
    ne <- nonEstimable(d)
    if (!is.null(ne)) cat("\nWARNING - non-estimable (confounded) coefficients:",
                          paste(ne, collapse = ", "), "\n")
  })

  # ======================= contrasts =======================
  contrasts_rv <- reactiveVal(data.frame(name = character(), expr = character(), type = character(),
                                         explain = character(), stringsAsFactors = FALSE))
  ct_refresh <- reactiveVal(0)

  # one row per cell of the design: values of the grouping variables + the design-column key
  group_info <- reactive({
    md <- sample_meta()[common(), , drop = FALSE]; gv <- input$group_vars
    cells <- unique(md[, gv, drop = FALSE]); cells[] <- lapply(cells, as.character)
    cells$key <- make.names(vapply(seq_len(nrow(cells)), function(i) paste(unlist(cells[i, gv]), collapse = " "), ""))
    cells
  })
  lab_of <- function(key) {                       # readable group name for a design key
    ci <- group_info(); r <- ci[ci$key == key, input$group_vars, drop = FALSE]
    if (!nrow(r)) key else paste(unlist(r), collapse = " ")
  }
  grp_choices <- reactive({ lv <- de_levels(); setNames(lv, vapply(lv, lab_of, "")) })

  explain_pair <- function(kb, ka) {              # plain-language description of B - A
    ci <- group_info(); gv <- input$group_vars
    rb <- ci[ci$key == kb, gv, drop = FALSE]; ra <- ci[ci$key == ka, gv, drop = FALSE]
    if (!nrow(rb) || !nrow(ra)) return(paste(kb, "vs", ka))
    vb <- unlist(rb); va <- unlist(ra); diff <- gv[vb != va]; same <- setdiff(gv, diff)
    if (length(diff) == 1) {
      txt <- sprintf("Difference of %s vs %s (%s)", vb[[diff]], va[[diff]], diff)
      if (length(same)) txt <- paste0(txt, " in ", paste(sprintf("%s (%s)", vb[same], same), collapse = " and "))
      txt
    } else sprintf("%s vs %s (they differ in %s)", paste(vb, collapse = " "), paste(va, collapse = " "), paste(diff, collapse = " and "))
  }
  auto_name_pair <- function(kb, ka)                # e.g. KO_Con_vs_WT_Con
    clean_contrast_name(paste(lab_of(kb), "vs", lab_of(ka)))

  # an already existing contrast (same name or same expression) is simply ignored
  add_contrast <- function(name, expr, type = "t", explain = "") {
    cur <- contrasts_rv()
    name <- clean_contrast_name(name)
    norm <- function(x) gsub("\\s+", "", x)
    if (name %in% cur$name || norm(expr) %in% norm(cur$expr)) return(invisible())
    comps <- trimws(strsplit(expr, "|", fixed = TRUE)[[1]])
    if (type == "F") {                                    # members must be existing ordinary contrasts
      if (length(comps) < 2 || !all(comps %in% cur$name[cur$type != "F"])) {
        showNotification("An F-test needs two or more existing contrasts.", type = "warning"); return(invisible())
      }
    } else if (!length(comps) || !all(vapply(comps, valid_contrast_expr, logical(1), allowed_names = colnames(de_design())))) {
      showNotification(paste0("Invalid expression: ", expr,
                              " (use group names, + - * / and parentheses only)"), type = "error")
      return(invisible())
    }
    contrasts_rv(rbind(cur, data.frame(name = name, expr = expr, type = type, explain = explain, stringsAsFactors = FALSE)))
  }

  # ---- 1. compare two groups (B - A) ----
  output$pair_builder <- renderUI({
    ch <- grp_choices()
    fluidRow(
      column(5, selectInput("c_b", "B", ch, selected = ch[min(2, length(ch))])),
      column(2, div(style = "text-align:center; font-size:28px; margin-top:22px;", HTML("&minus;"))),
      column(5, selectInput("c_a", "A (reference)", ch)))
  })
  observeEvent(input$add_c, {
    req(input$c_a, input$c_b)
    nm <- trimws(input$cname); if (!nzchar(nm)) nm <- auto_name_pair(input$c_b, input$c_a)
    add_contrast(nm, paste(input$c_b, "-", input$c_a), "t", explain_pair(input$c_b, input$c_a))
    updateTextInput(session, "cname", value = "")
  })
  observeEvent(input$add_all, {
    lv <- de_levels(); if (length(lv) < 2) return()
    for (p in utils::combn(lv, 2, simplify = FALSE))
      add_contrast(auto_name_pair(p[2], p[1]), paste(p[2], "-", p[1]), "t", explain_pair(p[2], p[1]))
  })

  # ---- advanced: manual interaction (D - C) - (B - A), F-test over existing contrasts, custom expression ----
  output$adv_builder <- renderUI({
    req(input$ctype); ch <- grp_choices()
    switch(input$ctype,
      inter = tagList(
        div(style = "border:1px solid #ddd; border-radius:4px; padding:6px 8px; margin-bottom:4px;",
            tags$b("(D - C)"),
            fluidRow(column(6, selectInput("ia_d", "D", ch, selected = ch[min(4, length(ch))])),
                     column(6, selectInput("ia_c", "C (reference)", ch, selected = ch[min(3, length(ch))])))),
        div(style = "text-align:center; font-size:24px; line-height:20px;", HTML("&minus;")),
        div(style = "border:1px solid #ddd; border-radius:4px; padding:6px 8px; margin-bottom:6px;",
            tags$b("(B - A)"),
            fluidRow(column(6, selectInput("ia_b", "B", ch, selected = ch[min(2, length(ch))])),
                     column(6, selectInput("ia_a", "A (reference)", ch)))),
        helpText("Interaction = (D - C) - (B - A): is the difference D - C different from B - A?")),
      fsel = {
        cur <- contrasts_rv(); tn <- cur$name[cur$type != "F"]
        if (length(tn) < 2) helpText("Define at least two contrasts first - the F-test combines existing contrasts.")
        else tagList(
          checkboxGroupInput("fsel_members", "Contrasts tested together (two or more)", choices = tn),
          helpText("Typical uses: several interactions (e.g. Int_KO1 + Int_KO2) or several baseline differences",
                   "(e.g. KO1_vs_WT + KO2_vs_WT). The joint F-test finds proteins that change in ANY of them;",
                   "the direction of each contrast comes from decideTests(method = 'global')."))
      },
      custom = tagList(
        textInput("c_custom", "Expression using group names", placeholder = "(KO_Tr - KO_Ct) - (WT_Tr - WT_Ct)"),
        helpText("Groups:", paste(names(ch), collapse = ", "), "- written with dots instead of spaces, e.g. WT.Con")))
  })
  observeEvent(input$add_adv, {
    nm <- trimws(input$cname_adv); n_now <- nrow(contrasts_rv()) + 1
    switch(input$ctype,
      inter = {
        req(input$ia_b, input$ia_a, input$ia_d, input$ia_c)
        add_contrast(if (nzchar(nm)) nm else paste0("Interaction ", n_now),
                     sprintf("(%s - %s) - (%s - %s)", input$ia_d, input$ia_c, input$ia_b, input$ia_a), "t",
                     sprintf("Difference of differences: (%s - %s) - (%s - %s)",
                             lab_of(input$ia_d), lab_of(input$ia_c), lab_of(input$ia_b), lab_of(input$ia_a)))
      },
      fsel = {
        mem <- input$fsel_members
        if (length(mem) < 2) return(showNotification("Select two or more contrasts for the F-test.", type = "warning"))
        add_contrast(if (nzchar(nm)) nm else paste0("F test ", sum(contrasts_rv()$type == "F") + 1),
                     paste(mem, collapse = " | "), "F",
                     sprintf("F-test over %s: does any of them differ from zero? Statistics only (F, adj. p); direction of each contrast from decideTests(method = 'global').",
                             paste(mem, collapse = ", ")))
      },
      custom = {
        ex <- trimws(input$c_custom); req(nzchar(ex))
        add_contrast(if (nzchar(nm)) nm else paste0("Contrast ", n_now), ex, "t", "Custom expression")
      })
    updateTextInput(session, "cname_adv", value = "")
  })

  # ---- contrasts table: expr | name (click to rename) | explanation | reorder / remove ----
  output$contrast_tbl <- renderUI({
    ct_refresh()
    cur <- contrasts_rv()
    if (!nrow(cur)) return(helpText("No contrasts defined yet."))
    n <- nrow(cur)
    btn <- function(i, act, label, disabled) tags$button(
      type = "button", class = "btn btn-xs btn-default", disabled = if (disabled) NA else NULL,
      onclick = sprintf("contrastAct(%d,'%s')", i, act), label)
    rows <- lapply(seq_len(n), function(i) tags$tr(
      tags$td(if (identical(cur$type[i], "F"))
                tags$span(tags$b("F-test: "), tags$code(cur$expr[i]))
              else tags$code(cur$expr[i])),
      tags$td(cur$name[i], contenteditable = "true", title = "Click to rename",
              style = "cursor:text; border-bottom:1px dashed #999; min-width:80px;",
              onblur = sprintf("contrastRename(this, %d)", i),
              onkeydown = "if (event.key === 'Enter') { event.preventDefault(); this.blur(); }"),
      tags$td(tags$small(cur$explain[i])),
      tags$td(style = "white-space:nowrap;",
              btn(i, "up", HTML("&#9650;"), i == 1), btn(i, "down", HTML("&#9660;"), i == n),
              tags$button(type = "button", class = "btn btn-xs btn-link", title = "Remove this contrast",
                          style = "color:#c0392b;", onclick = sprintf("contrastAct(%d,'del')", i), icon("trash")))))
    tagList(
      tags$table(class = "table table-condensed table-hover",
        tags$thead(tags$tr(tags$th("expr"), tags$th("name (editable)"), tags$th("explanation"), tags$th(""))), tags$tbody(rows)),
      actionButton("clr_all", "Clear all", icon = icon("trash"), class = "btn-xs"))
  })
  observeEvent(input$clr_all, contrasts_rv(contrasts_rv()[0, ]))
  observeEvent(input$c_rename, {
    a <- input$c_rename; cur <- contrasts_rv(); i <- a$i
    if (is.null(i) || i < 1 || i > nrow(cur)) return()
    nm <- clean_contrast_name(a$name)
    if (!nzchar(a$name) || nm %in% cur$name[-i]) {
      ct_refresh(ct_refresh() + 1)
      return(showNotification("The name is empty or already used by another contrast.", type = "warning"))
    }
    if (nm != cur$name[i]) {
      old <- cur$name[i]; cur$name[i] <- nm
      isf <- cur$type == "F"                        # F-tests refer to contrasts by name: keep them in sync
      cur$expr[isf] <- vapply(cur$expr[isf], function(e) {
        mem <- trimws(strsplit(e, "|", fixed = TRUE)[[1]]); paste(ifelse(mem == old, nm, mem), collapse = " | ") }, "")
      cur$explain[isf] <- gsub(old, nm, cur$explain[isf], fixed = TRUE)
      contrasts_rv(cur)
    } else ct_refresh(ct_refresh() + 1)
  })
  observeEvent(input$c_action, {
    a <- input$c_action; cur <- contrasts_rv(); n <- nrow(cur); i <- a$i
    if (is.null(i) || i < 1 || i > n) return()
    idx <- seq_len(n)
    if (a$a == "del") {
      idx <- idx[-i]
      gone <- cur$name[i]                           # drop F-tests that used the removed contrast
      uses <- cur$type == "F" & vapply(cur$expr, function(e) gone %in% trimws(strsplit(e, "|", fixed = TRUE)[[1]]), logical(1))
      if (any(uses)) { idx <- setdiff(idx, which(uses)); showNotification("F-test(s) using the removed contrast were removed too.", type = "message") }
    }
    if (a$a == "up" && i > 1) idx[c(i - 1, i)] <- idx[c(i, i - 1)]
    if (a$a == "down" && i < n) idx[c(i, i + 1)] <- idx[c(i + 1, i)]
    contrasts_rv(cur[idx, , drop = FALSE])
  })

  # ---- technical replicates (animals) ----
  output$block_col_ui <- renderUI({
    req(other_cols())
    selectInput("block_col", "Column identifying the biological unit / animal", other_cols())
  })

  de_block <- reactive({
    if (!isTRUE(input$de_block)) return(NULL)
    md <- de_meta()
    b <- if (identical(input$block_src, "col")) {
      req(input$block_col); as.character(md[[input$block_col]])
    } else {
      sub(input$block_regex, "", md$Sample)
    }
    b
  })

  output$block_preview <- renderPrint({
    b <- de_block(); req(b)
    cat(length(b), "samples ->", length(unique(b)), "animals\n")
    cat("Example:", paste(head(paste(de_meta()$Sample, "->", b), 3), collapse = " | "), "\n")
    if (length(unique(b)) == length(b))
      cat("Note: every sample is its own animal (no replicates detected).\n")
  })

  # percentage change corresponding to the fold change
  output$de_fc_pct <- renderText({
    fc <- input$de_fc %||% 1.2
    if (is.na(fc) || fc < 1) return("")
    sprintf("= +%.1f%% / -%.1f%%", (fc - 1) * 100, (1 - 1 / fc) * 100)
  })

  # ---- run ----
  de_res <- reactiveVal(NULL)

  observeEvent(input$run_de, {
    tryCatch({
      cdf <- contrasts_rv()
      validate(need(nrow(cdf) > 0, "Add at least one contrast first."))
      d <- de_design()
      ne <- nonEstimable(d)
      if (!is.null(ne)) stop("Design is not full rank; confounded: ", paste(ne, collapse = ", "))
      cm <- make_contrast_matrix(cdf, d)
      y <- final_expr()[, rownames(d), drop = FALSE]
      res <- withProgress(message = "Running limma", value = 0.5,
        run_limma_pipeline(y, d, cm, block = de_block(), use_aw = isTRUE(input$de_aw),
                           test = input$de_test, fc = input$de_fc,
                           trend = input$de_trend, robust = input$de_robust,
                           lfc = log2(input$de_fc), dt_method = input$dt_method,
                           adjust.method = input$dt_adjust, p.value = input$dt_p,
                           progress = function(m) setProgress(detail = m)))
      de_res(res)
      updateTabsetPanel(session, "de_tabs", selected = "Results")
      showNotification("Differential analysis finished.", type = "message")
    }, error = function(e) showNotification(conditionMessage(e), type = "error", duration = NULL))
  })

  fname <- function(x) gsub("[^A-Za-z0-9]+", "_", x)       # file names only; contrast names keep their spaces

  # ordinary contrasts: Up / Down / NS, with a download icon (full limma table) after the last column
  output$de_summary <- renderUI({
    r <- de_res(); req(r)
    m <- dt_matrix(r$dt)[, r$single_names, drop = FALSE]
    rows <- lapply(colnames(m), function(cn) tags$tr(
      tags$td(cn), tags$td(sum(m[, cn] == 1)), tags$td(sum(m[, cn] == -1)), tags$td(sum(m[, cn] == 0)),
      tags$td(tags$a(icon("download"), href = csv_data_uri(as_names(r$tops[[cn]])),
                     download = paste0("limma_", fname(cn), ".csv"), class = "btn btn-xs btn-default",
                     title = "Full limma result of this contrast"))))
    tags$table(class = "table table-condensed table-hover", style = "width:auto;",
      tags$thead(tags$tr(tags$th("Contrast"), tags$th("Up"), tags$th("Down"), tags$th("NS"), tags$th(""))),
      tags$tbody(rows))
  })

  # model parameters: directly under the contrasts matrix, above the F-test block
  output$de_meta_txt <- renderUI({
    r <- de_res(); req(r)
    tags$p(style = "color:#555;", sprintf("Test: %s | duplicateCorrelation consensus: %s | arrayWeights: %s",
           r$test, if (is.null(r$correlation)) "not used" else round(r$correlation, 3),
           if (is.null(r$weights)) "not used" else "used"))
  })

  # joint F-tests: same Up / Down / NS structure, one row per component contrast
  output$ftest_summary_ui <- renderUI({
    r <- de_res(); req(r); fn <- r$ftest_names
    if (!length(fn)) return(NULL)
    tagList(
      hr(), h5("Joint F-tests / multi-coefficient effects"),
      lapply(fn, function(nm) {
        tt <- r$tops[[nm]]
        tags$p(tags$b(nm), sprintf(" (joint F-test, adj. p < %g, plus a decideTests call) ", r$p_value),
               tags$a(icon("download"), nm, href = csv_data_uri(as_names(tt)),
                      download = paste0("F_test_", fname(nm), ".csv"), class = "btn btn-xs btn-default",
                      title = "Full statistics of this F-test"))
      }),
      tableOutput("f_summary_tbl"))
  })
  output$f_summary_tbl <- renderTable({
    r <- de_res(); req(r, length(r$ftest_names) > 0)
    do.call(rbind, lapply(r$ftest_names, function(nm) {
      tt <- r$tops[[nm]]; cols <- grep("_call$", names(tt), value = TRUE)
      data.frame(`F-test` = nm, Contrast = sub("_call$", "", cols),
                 Up = vapply(cols, function(cn) sum(tt[[cn]] == 1, na.rm = TRUE), numeric(1)),
                 Down = vapply(cols, function(cn) sum(tt[[cn]] == -1, na.rm = TRUE), numeric(1)),
                 NS = vapply(cols, function(cn) sum(tt[[cn]] == 0, na.rm = TRUE), numeric(1)),
                 check.names = FALSE, row.names = NULL)
    }))
  }, digits = 0)

  output$de_results_ui <- renderUI({
    req(de_res())
    cn_all <- names(de_res()$tops)                 # every result (incl. global F-test effects)
    cn <- de_res()$single_names                    # volcano plots need a log fold change
    lfc0 <- isolate(round(log2(input$de_fc), 2))
    tagList(
      downloadButton("dl_de_all", "All contrasts (CSV)"),
      hr(), h4("Volcano plots"),
      fluidRow(
        column(3, tags$label("Contrasts"), br(),
               shinyWidgets::dropdownButton(
                 label = "Select contrasts", circle = FALSE, status = "default", inline = TRUE,
                 checkboxGroupInput("vol_contrasts", NULL, choices = cn, selected = cn))),
        column(2, selectInput("vol_filter", "Significance on", c("adj.P.Val", "P.Value"))),
        column(2, selectInput("vol_y", "Y axis", c("P.Value", "adj.P.Val"))),
        column(2, numericInput("vol_lfc", "|log2FC| cutoff", lfc0, min = 0, step = 0.05)),
        column(3, numericInput("vol_p", "p cutoff", isolate(input$dt_p), min = 0, max = 1, step = 0.01))),
      fluidRow(
        column(2, numericInput("vol_nlab", "Number of labels", 10, min = 0, step = 1)),
        column(2, selectInput("vol_which", "Label which", c("all", "up", "down"))),
        column(2, numericInput("vol_label_size", "Label size", 2.5, min = 1, max = 8, step = 0.5)),
        column(2, br(), checkboxInput("vol_grid", "Grid lines", FALSE)),
        column(4, selectizeInput("vol_highlight", "Search / highlight proteins",
                                 choices = rownames(de_res()$dt), multiple = TRUE,
                                 options = list(maxOptions = 100, placeholder = "type a protein name...")))),
      uiOutput("volcano_cards"),
      hr(), h4("Differentially expressed proteins per contrast"),
      fluidRow(column(3, tags$label("Comparisons"), br(),
                      shinyWidgets::dropdownButton(label = "Select comparisons", circle = FALSE, status = "default", inline = TRUE,
                                                   checkboxGroupInput("bar_contrasts", NULL, choices = cn, selected = cn))),
               column(3, checkboxInput("bar_labels", "Show contrast names", TRUE)),
               column(3, numericInput("bar_num_size", "Count font size (pt)", 12, min = 4, max = 30)),
               column(3, br(), downloadButton("dl_bar", "Download PNG"))),
      plotOutput("deg_bar", height = "360px"),
      hr(), h4("Venn diagram of significant proteins"),
      fluidRow(
        column(4, tags$label("Contrasts"), br(),
               shinyWidgets::dropdownButton(label = "Select contrasts", circle = FALSE, status = "default", inline = TRUE,
                 checkboxGroupInput("venn_sig_contrasts", NULL, choices = cn_all,
                                    selected = cn[seq_len(min(2, length(cn)))]),
                 helpText("Pick 2-7 contrasts, or ONE F-test: it is shown as the contrasts it combines.",
                          "Choosing an F-test unticks the others and vice versa."))),
        column(3, radioButtons("venn_sig_dir", "Direction", c("either", "up", "down")))),
      fluidRow(
        column(3, sliderInput("vs_font", "Overall font scale", 0.5, 3, 1, step = 0.1)),
        column(3, numericInput("vs_count", "Number size (blank = auto)", NA, min = 1, step = 0.5)),
        column(3, numericInput("vs_label", "Contrast-name size (blank = auto)", NA, min = 1, step = 0.5)),
        column(3, sliderInput("vs_ratio", "Circle size vs numbers (higher = smaller circles)", 1, 3, 1.5, step = 0.1))),
      uiOutput("vs_titles_ui"),
      uiOutput("venn_note"),
      plotOutput("venn_sig_plot", height = "520px"),
      downloadButton("dl_venn_sig", "Download PNG"))
  })

  # ---- volcano cards: side by side, "Show data" toggle + PNG download under each ----
  output$volcano_cards <- renderUI({
    sel <- input$vol_contrasts; req(length(sel) > 0)
    cn <- names(de_res()$tops)
    div(style = "display:grid; grid-template-columns:repeat(2, minmax(0, 1fr)); gap:14px;",   # two equal columns
      lapply(sel, function(s) {
        i <- match(s, cn)
        div(style = "min-width:0;",
          plotOutput(paste0("vol_plot_", i), width = "100%", height = "520px"),
          div(actionButton(paste0("vol_tog_", i), "Show data", icon = icon("table"), class = "btn-xs"),
              downloadButton(paste0("vol_dl_", i), "Download PNG", class = "btn-xs")),
          conditionalPanel(sprintf("input.vol_tog_%d %% 2 == 1", i),
            div(style = "max-height:300px; overflow:auto; margin-top:6px;",
                tableOutput(paste0("vol_tbl_", i)))))
      }))
  })

  volcano_plot <- function(i) {
    req(de_res(), input$vol_filter, input$vol_y, input$vol_lfc, input$vol_p)
    nl <- input$vol_nlab %||% 0
    volcano_from_top(de_res()$tops[[i]], title = names(de_res()$tops)[i],
                     filter_choice = input$vol_filter, y_axis_choice = input$vol_y,
                     l2FC_thresh = input$vol_lfc, p_thresh = input$vol_p,
                     show_names = if (is.na(nl) || nl <= 0) "none" else nl,
                     highlight = input$vol_highlight %||% character(0),
                     which_genes = input$vol_which, grid_lines = isTRUE(input$vol_grid),
                     label_size = input$vol_label_size)
  }

  observeEvent(de_res(), {
    for (i in seq_along(de_res()$tops)) local({
      ii <- i
      output[[paste0("vol_plot_", ii)]] <- renderPlot(volcano_plot(ii), res = 100)
      output[[paste0("vol_tbl_", ii)]] <- renderTable(head(as_names(de_res()$tops[[ii]]), 300), digits = 4)
      output[[paste0("vol_dl_", ii)]] <- downloadHandler(
        filename = function() paste0("Volcano_", names(de_res()$tops)[ii], ".png"),
        content = function(file) ggplot2::ggsave(file, volcano_plot(ii), width = 6.5, height = 5.8, dpi = 300))
      observeEvent(input[[paste0("vol_tog_", ii)]], {
        updateActionButton(session, paste0("vol_tog_", ii),
                           label = if (input[[paste0("vol_tog_", ii)]] %% 2 == 1) "Hide data" else "Show data")
      }, ignoreInit = TRUE)
    })
  })

  bar_dt <- reactive({
    r <- de_res(); req(r)
    sel <- intersect(input$bar_contrasts %||% r$single_names, r$single_names)      # F-tests are never plotted here
    validate(need(length(sel) > 0, "Select at least one comparison."))
    dt_matrix(r$dt)[, sel, drop = FALSE]
  })
  output$deg_bar <- renderPlot({
    diverging_bar_plot(deg_summary_from_dt(bar_dt()), isTRUE(input$bar_labels), input$bar_num_size %||% 12)
  }, res = 100)
  output$dl_bar <- downloadHandler("DEP_summary_diverging_bar.png", function(file) {
    n <- ncol(bar_dt())
    ggplot2::ggsave(file, diverging_bar_plot(deg_summary_from_dt(bar_dt()), isTRUE(input$bar_labels), input$bar_num_size %||% 12),
                    width = if (isTRUE(input$bar_labels)) 7.5 else 5.5, height = max(3, 0.9 * n + 1.5), dpi = 300)
  })

  output$dl_de_all <- downloadHandler(
    filename = function() "DE_all_contrasts.csv",
    content = function(file) {
      dt <- de_res()$dt
      wide <- setNames(data.frame(rownames(dt), check.names = FALSE), id_label())
      for (cn in names(de_res()$tops)) {
        tt <- de_res()$tops[[cn]]; idx <- match(rownames(dt), tt$Protein)
        if ("F" %in% colnames(tt)) wide[[paste0(cn, "_F")]] <- tt$F[idx]
        wide[[paste0(cn, "_logFC")]]     <- tt$logFC[idx]
        wide[[paste0(cn, "_adj.P.Val")]] <- tt$adj.P.Val[idx]
        wide[[paste0(cn, "_status")]]    <- dt[, cn]
      }
      write.csv(wide, file, row.names = FALSE)
    })

  # ======================= 3. SUBSETS =======================
  dtm <- reactive({
    validate(need(!is.null(de_res()), "Run the differential analysis first (tab 2)."))
    dt_matrix(de_res()$dt)
  })

  pre_subsets <- reactiveVal(list())
  custom_subsets <- reactiveVal(list())
  make_pre_subsets <- function(methods, min_proteins) {
    unlist(lapply(methods, function(mt)
      get_all_subsets_from_decideTests(de_res()$dt, method = mt, min_proteins = min_proteins)),
      recursive = FALSE)
  }
  # all subsets (per-contrast AND venn) are generated automatically when the DE analysis finishes
  observeEvent(de_res(), {
    custom_subsets(list())
    pre_subsets(make_pre_subsets(c("contrast", "venn"), input$min_prot %||% 5) %||% list())
  })
  observeEvent(input$gen_pre, {
    s <- make_pre_subsets(input$pre_method %||% c("contrast", "venn"), input$min_prot %||% 5)
    if (!length(s)) showNotification("No subset reached the minimum size.", type = "warning")
    pre_subsets(s %||% list())
  })

  output$rules_ui <- renderUI({
    cols <- colnames(dtm())
    lapply(seq_along(cols), function(i)
      selectInput(paste0("rule_", i), cols[i],
                  c("Ignore" = "Ignore", "Up" = "Up", "Down" = "Down",
                    "Any significant" = "Any"), width = "100%"))
  })

  rule_hits <- reactive({
    m <- dtm(); cols <- colnames(m)
    rules <- vapply(seq_along(cols), function(i) input[[paste0("rule_", i)]] %||% "Ignore", character(1))
    names(rules) <- cols
    subset_from_rules(m, rules, input$rule_mode, isTRUE(input$rule_only))
  })
  output$rule_count <- renderText(paste(length(rule_hits()), "proteins match the current rules"))

  observeEvent(input$add_custom, {
    hits <- rule_hits()
    if (!length(hits)) return(showNotification("No proteins match these rules.", type = "warning"))
    nm <- trimws(input$rule_name)
    if (!nzchar(nm)) nm <- paste0("custom_", length(custom_subsets()) + 1)
    cs <- custom_subsets(); cs[[nm]] <- hits; custom_subsets(cs)
    updateTextInput(session, "rule_name", value = "")
  })

  all_subsets <- reactive(c(pre_subsets(), custom_subsets()))

  output$subset_tbl <- renderUI({
    s <- all_subsets()
    validate(need(length(s) > 0, "No subsets yet. Generate predefined ones or add a custom one."))
    trs <- lapply(names(s), function(nm) tags$tr(
      tags$td(nm), tags$td(length(s[[nm]])),
      tags$td(tags$a(icon("download"), href = csv_data_uri(setNames(data.frame(s[[nm]]), id_label())),
                     download = paste0(make.names(nm), ".csv"), class = "btn btn-xs btn-default",
                     title = "Download this subset"))))
    tags$table(class = "table table-condensed table-hover",
      tags$thead(tags$tr(tags$th("Subset"), tags$th("Proteins"), tags$th(""))), tags$tbody(trs))
  })

  output$dl_subsets <- downloadHandler(
    filename = function() "subsets.csv",
    content = function(file) {
      s <- all_subsets(); req(length(s) > 0)
      n <- max(lengths(s))
      write.csv(as.data.frame(lapply(s, function(x) { length(x) <- n; x }), check.names = FALSE),
                file, row.names = FALSE, na = "")
    })


  # ---- Venn of significant proteins: redraws automatically when any option changes ----
  # an F-test and ordinary contrasts exclude each other in the selection
  venn_prev <- reactiveVal(character())
  observeEvent(de_res(), venn_prev(character()))
  observeEvent(input$venn_sig_contrasts, ignoreNULL = FALSE, {
    sel <- input$venn_sig_contrasts %||% character(); fn <- de_res()$ftest_names %||% character()
    new <- setdiff(sel, venn_prev()); out <- sel
    if (length(new)) {
      nf <- intersect(new, fn)
      out <- if (length(nf)) nf[length(nf)] else setdiff(sel, fn)      # newest pick wins
    } else if (sum(sel %in% fn) > 1) out <- tail(intersect(sel, fn), 1)
    venn_prev(out)
    if (!setequal(out, sel)) updateCheckboxGroupInput(session, "venn_sig_contrasts", selected = out)
  })

  # standalone contrasts: the decideTests calls of the main analysis;
  # F-test: the calls of its components from decideTests(method = "global") inside that F-test
  venn_run <- reactive({
    r <- de_res(); req(r); sel <- input$venn_sig_contrasts %||% character(); fs <- r$fsets
    f <- intersect(sel, names(fs))
    if (length(f)) {
      nm <- f[1]; tt <- r$tops[[nm]]; mem <- fs[[nm]]
      calls <- as.matrix(tt[, paste0(mem, "_call"), drop = FALSE]); colnames(calls) <- mem; rownames(calls) <- tt$Protein
      list(members = mem, dt = calls, kind = "F", name = nm)
    } else list(members = sel, dt = r$dt, kind = "t", name = NA)
  })
  venn_members <- reactive(venn_run()$members)
  output$venn_note <- renderUI({
    v <- venn_run()
    tags$p(style = "color:#555; font-size:12px;",
           if (v$kind == "F") sprintf("F-test '%s': circles are its contrasts, with the calls of decideTests(method = 'global') computed inside this F-test.", v$name)
           else sprintf("Standalone contrasts: calls of the main analysis (decideTests, method = '%s').", input$dt_method %||% "separate"))
  })
  venn_sig_obj <- reactive({
    v <- venn_run(); sel <- v$members
    validate(need(length(sel) >= 2 && length(sel) <= 7, "Select 2-7 contrasts, or one F-test."))
    na_null <- function(x) if (is.null(x) || is.na(x)) NULL else x
    titles <- vapply(seq_along(sel), function(i) input[[paste0("vst_", i)]] %||% sel[i], character(1))
    venn_sig_plot(v$dt, sel, input$venn_sig_dir %||% "either", font_scale = input$vs_font %||% 1,
                  count_size = na_null(input$vs_count), group_label_size = na_null(input$vs_label),
                  expand_limits = input$vs_ratio %||% 1.5, labels = titles)
  })
  output$vs_titles_ui <- renderUI({
    sel <- venn_members(); req(length(sel) > 0)
    fluidRow(lapply(seq_along(sel), function(i)
      column(3, textInput(paste0("vst_", i), paste("Title of circle", i), value = gsub("_", " ", sel[i])))))
  })
  output$venn_sig_plot <- renderPlot({ print(venn_sig_obj()) }, res = 100)
  output$dl_venn_sig <- downloadHandler("Venn_significant.png", function(file)
    ggplot2::ggsave(file, venn_sig_obj(), width = 8 / sqrt(input$vs_ratio), height = 7 / sqrt(input$vs_ratio),
                    dpi = 300, bg = "white"))

  # ---- smart tab switching ----
  observeEvent(input$add_c,   updateTabsetPanel(session, "de_tabs", selected = "Design & contrasts"))
  observeEvent(input$add_adv, updateTabsetPanel(session, "de_tabs", selected = "Design & contrasts"))
  observeEvent(input$add_all, updateTabsetPanel(session, "de_tabs", selected = "Design & contrasts"))
  observeEvent(input$de_covs, updateTabsetPanel(session, "de_tabs", selected = "Design & contrasts"),
               ignoreInit = TRUE, ignoreNULL = FALSE)
