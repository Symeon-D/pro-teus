  # ======================= 2. DIFFERENTIAL =======================
  # built with its choices directly so that clicking the box always lists the variables
  output$de_covs_ui <- renderUI({
    app_gen(); req(other_cols())
    sv <- colnames(sva_sv_df() %||% data.frame())
    ch <- c(other_cols(), sv)
    selectizeInput("de_covs", "Additional covariates (e.g. Batch; SVs from the Batch & SVA tab are preselected)",
                   choices = ch, multiple = TRUE, selected = union(intersect(isolate(input$de_covs), ch), sv))
  })

  de_meta <- reactive({
    cs <- common()
    md <- sample_meta()[cs, , drop = FALSE]
    sv <- sva_sv_df()                                          # surrogate variables chosen on the Batch & SVA tab
    if (!is.null(sv)) md <- cbind(md, sv[rownames(md), , drop = FALSE])
    md$Group_de <- factor(make.names(md$Group), levels = make.names(unique(md$Group)))
    md
  })
  de_levels <- reactive(levels(de_meta()$Group_de))

  # data with exactly the covariates of the DE design removed (NULL when the design had none)
  corrected_expr <- reactive({
    r <- de_res(); req(r); covs <- r$covs %||% character(0)
    if (!length(covs)) return(NULL)
    md <- de_meta(); d <- r$design; sm <- rownames(d)
    correct_for_covariates(final_expr()[, sm, drop = FALSE], md[sm, , drop = FALSE], covs, md[sm, "Group_de"])
  })
  corrected_ok <- reactive(!is.null(tryCatch(corrected_expr(), error = function(e) NULL)))

  # ---- protein lists: significant proteins of the chosen contrasts, largest |t| first ----
  ordered_lists <- function(contrasts) {
    r <- de_res(); req(r)
    cns <- intersect(contrasts, r$single_names); allp <- rownames(r$dt)
    if (!length(cns)) return(list(sig = character(0), all = allp))
    dtm <- dt_matrix(r$dt)
    tm <- vapply(cns, function(cn) { tt <- r$tops[[cn]]; v <- abs(setNames(tt$t, tt$Protein)[allp]); v[is.na(v)] <- 0; v }, numeric(length(allp)))
    sg <- vapply(cns, function(cn) dtm[allp, cn] != 0, logical(length(allp)))
    tm <- matrix(tm, nrow = length(allp), dimnames = list(allp, cns)); sg <- matrix(sg, nrow = length(allp), dimnames = list(allp, cns))
    is_sig <- rowSums(sg) > 0
    sig <- allp[is_sig][order(-apply(tm[is_sig, , drop = FALSE] * sg[is_sig, , drop = FALSE], 1, max))]
    rest <- allp[!is_sig][order(-apply(tm[!is_sig, , drop = FALSE], 1, max))]
    list(sig = sig, all = c(sig, rest))
  }
  picker_inputs <- function(lists_fn) list(
    sig = reactive(if (is.null(de_res())) character(0) else lists_fn()$sig),
    all = reactive(if (is.null(de_res())) rownames(final_expr()) else lists_fn()$all))

  # ---- volcano search bar: N <-> bar ----
  vol_lists <- reactive(ordered_lists(input$vol_contrasts %||% de_res()$single_names))
  vi <- picker_inputs(vol_lists)
  vol_picker <- protein_picker_server("pp_vol", sig = vi$sig, all = vi$all)
  vol_skip <- reactiveVal(NULL)
  vol_topn <- function(n) {                                    # top n by |t| of every selected volcano
    r <- de_res(); cns <- intersect(input$vol_contrasts %||% character(0), r$single_names)
    unique(unlist(lapply(cns, function(cn) { tt <- r$tops[[cn]]; utils::head(tt$Protein[order(-abs(tt$t))], n) })))
  }
  observeEvent(list(input$vol_nlab, input$vol_contrasts, de_res(), vol_picker$ready()), {
    n <- input$vol_nlab; req(de_res())
    if (is.null(n) || is.na(n) || n < 0) return()
    if (!is.null(vol_skip()) && identical(as.numeric(n), as.numeric(vol_skip()))) { vol_skip(NULL); return() }   # N was set from the bar
    vol_picker$set(if (n == 0) character(0) else vol_topn(n))
  })
  observeEvent(vol_picker$selected(), {                         # bar -> N (exact with a single volcano selected)
    if (length(input$vol_contrasts %||% character(0)) != 1) return()
    len <- length(vol_picker$selected()); n <- input$vol_nlab
    if (!is.null(n) && !is.na(n) && as.numeric(n) != len) { vol_skip(len); updateNumericInput(session, "vol_nlab", value = len) }
  }, ignoreInit = TRUE)

  # ---- protein plots ----
  pp_lists <- reactive(ordered_lists(input$pp_contrasts %||% de_res()$single_names))
  pi_ <- picker_inputs(pp_lists)
  pp_picker <- protein_picker_server("pp_de", sig = pi_$sig, all = pi_$all)
  pp_default <- reactive({                                     # top 2 (largest |t|) of every selected contrast
    r <- de_res(); req(r); cns <- intersect(input$pp_contrasts %||% r$single_names, r$single_names)
    unique(unlist(lapply(cns, function(cn) { tt <- r$tops[[cn]]; utils::head(tt$Protein[order(-abs(tt$t))], 2) })))
  })
  observeEvent(list(input$pp_contrasts, de_res(), pp_picker$ready()), {    # default: top 2 of each contrast, shown in the bar
    req(de_res()); pp_picker$set(pp_default())
  })
  pp_proteins <- reactive(pp_picker$selected())
  pp_expr <- reactive(if (corrected_ok()) corrected_expr() else final_expr())      # corrected whenever the design had covariates

  # significance from the limma results (adjusted p of each contrast): plain two-group contrasts get a bracket,
  # other contrasts (interactions ...) a note under the protein name
  pp_stats <- reactive({
    r <- de_res(); req(r); cdf <- r$contrast_df; prot <- pp_proteins(); show_ns <- isTRUE(input$pp_ns)
    cns <- intersect(input$pp_contrasts %||% r$single_names, r$single_names)
    star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "ns")))
    br <- list(); nt <- list()
    for (cn in cns) {
      row <- cdf[cdf$name == cn, , drop = FALSE]; if (!nrow(row)) next
      tt <- r$tops[[cn]]; padj <- setNames(tt$adj.P.Val, tt$Protein)[prot]; lab <- star(padj)
      keep <- !is.na(padj) & (show_ns | lab != "ns"); if (!any(keep)) next
      mt <- regmatches(row$expr[1], regexec("^\\s*([A-Za-z0-9._]+)\\s*-\\s*([A-Za-z0-9._]+)\\s*$", row$expr[1]))[[1]]
      if (length(mt) == 3) br[[length(br) + 1]] <- data.frame(Protein = prot[keep], g1 = lab_of(mt[3]), g2 = lab_of(mt[2]), label = lab[keep])
      else nt[[length(nt) + 1]] <- data.frame(Protein = prot[keep], txt = paste0(cn, " ", lab[keep]))
    }
    brackets <- if (length(br)) do.call(rbind, br) else NULL
    notes <- if (length(nt)) { d <- do.call(rbind, nt); vapply(split(d$txt, d$Protein), paste, "", collapse = "; ") } else NULL
    list(brackets = brackets, notes = notes, n_brk = if (is.null(brackets)) 0 else max(table(brackets$Protein)))
  })
  pp_obj <- reactive({
    validate(need(length(pp_proteins()) > 0, "Select proteins in the search bar."))
    st <- pp_stats()
    p <- protein_plot_grid(pp_expr(), pp_proteins(), sample_meta(), input$pp_type %||% "boxplot",
                           input$pp_ncol %||% 4, isTRUE(input$pp_free %||% TRUE), group_colors_all(),
                           brackets = st$brackets, notes = st$notes)
    validate(need(length(pp_proteins()) > 0, "Select proteins in the search bar."), need(!is.null(p), "None of the selected proteins is in the data."))
    p
  })
  pp_h <- function() protein_plot_height(length(pp_proteins()), input$pp_ncol %||% 4, tryCatch(pp_stats()$n_brk, error = function(e) 0))
  output$pp_plot <- renderPlot(pp_obj(), res = 100, height = function() tryCatch(pp_h(), error = function(e) 400))
  output$dl_pp <- downloadHandler(filename = function() "protein_plots.png",
    content = function(file) {
      nc <- input$pp_ncol %||% 4; n <- length(pp_proteins())
      ggplot2::ggsave(file, pp_obj(), width = 3.3 * min(nc, max(n, 1)), height = pp_h() / 100 * 1.05,
                      dpi = 300, bg = "white", limitsize = FALSE)
    })

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
      res$covs <- input$de_covs %||% character(0)             # covariates of this design (for the corrected data)
      res$lfc <- log2(input$de_fc)                              # fold-change threshold of this run (volcano plots)
      res$contrast_df <- cdf                                   # contrasts of this run (for the protein-plot significance)
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
      ctl_row(
        ctl_dd("Contrasts", checkboxGroupInput("vol_contrasts", NULL, choices = cn, selected = cn)),
        selectInput("vol_y", "Y axis", c("adj.P.Val", "P.Value"), selected = "adj.P.Val", width = fit_w("Y axis", c("adj.P.Val", "P.Value"))),
        numericInput("vol_nlab", "Number of labels", 5, min = 0, step = 1, width = fit_w("Number of labels", "10")),
        numericInput("vol_label_size", "Label size", 2.5, min = 1, max = 8, step = 0.5, width = fit_w("Label size", "2.5")),
        checkboxInput("vol_grid", "Grid lines", FALSE)),
      protein_picker_ui("pp_vol", "Search / highlight proteins (the labelled proteins; the number above fills this bar)"),
      uiOutput("volcano_cards"),
      hr(), h4("Protein plots"),
      ctl_row(
        ctl_dd("Contrasts", checkboxGroupInput("pp_contrasts", NULL, choices = cn, selected = cn)),
        selectInput("pp_type", "Plot type", PROTEIN_PLOT_TYPES, width = fit_w("Plot type", names(PROTEIN_PLOT_TYPES))),
        numericInput("pp_ncol", "Subplots per row", 4, min = 1, max = 8, width = fit_w("Subplots per row", "4")),
        checkboxInput("pp_free", "Free y-axis", TRUE),
        checkboxInput("pp_ns", "show ns", FALSE)),
      protein_picker_ui("pp_de", "Proteins to plot (the top 2 of every selected contrast are selected by default)"),
      plotOutput("pp_plot", height = "auto"), downloadButton("dl_pp", "PNG"),
      hr(), h4("Differentially expressed proteins per contrast"),
      ctl_row(
        ctl_dd("Comparisons", checkboxGroupInput("bar_contrasts", NULL, choices = cn, selected = cn), btn = "Select comparisons"),
        checkboxInput("bar_labels", "Show contrast names", TRUE),
        numericInput("bar_num_size", "Count font size (pt)", 12, min = 4, max = 30, width = fit_w("Count font size (pt)", "12")),
        div(style = "margin-bottom:15px;", downloadButton("dl_bar", "Download PNG"))),
      plotOutput("deg_bar", height = "360px"),
      hr(), h4("Venn diagram of significant proteins"),
      ctl_row(
        ctl_dd("Contrasts", checkboxGroupInput("venn_sig_contrasts", NULL, choices = cn_all,
                                               selected = cn[seq_len(min(2, length(cn)))]),
               helpText("Pick 2-7 contrasts, or ONE F-test: it is shown as the contrasts it combines.",
                        "Choosing an F-test unticks the others and vice versa.")),
        radioButtons("venn_sig_dir", "Direction", c("either", "up", "down"), inline = TRUE)),
      ctl_row(
        sliderInput("vs_font", "Overall font scale", 0.5, 3, 1, step = 0.1, width = "190px"),
        numericInput("vs_count", "Number size (blank = auto)", NA, min = 1, step = 0.5, width = fit_w("Number size (blank = auto)", "5")),
        numericInput("vs_label", "Contrast-name size (blank = auto)", NA, min = 1, step = 0.5, width = fit_w("Contrast-name size (blank = auto)", "5")),
        sliderInput("vs_ratio", "Circle size vs numbers (higher = smaller circles)", 1, 3, 1.5, step = 0.1, width = "250px")),
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
    r <- de_res(); req(r, input$vol_y)
    volcano_from_top(r$tops[[i]], title = names(r$tops)[i],
                     filter_choice = "adj.P.Val", y_axis_choice = input$vol_y,
                     l2FC_thresh = r$lfc %||% log2(1.2), p_thresh = r$p_value %||% 0.05,   # the thresholds of the DE run
                     genes = vol_picker$selected(),            # the bar proteins are THE labels: boxed, coloured by side
                     grid_lines = isTRUE(input$vol_grid), label_size = input$vol_label_size)
  }

  observeEvent(de_res(), {
    for (i in seq_along(de_res()$tops)) local({
      ii <- i
      output[[paste0("vol_plot_", ii)]] <- renderPlot(volcano_plot(ii), res = 100)
      output[[paste0("vol_tbl_", ii)]] <- renderTable({ tt <- de_res()$tops[[ii]]; head(as_names(tt[order(-abs(tt$t)), ]), 300) }, digits = 4)
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
