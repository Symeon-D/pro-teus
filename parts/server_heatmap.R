  # ======================= 4. HEATMAP =======================
  output$inc_groups_ui <- renderUI({
    lv <- group_levels_all()
    selectizeInput("inc_groups", "Groups to include", choices = lv,
                   selected = lv, multiple = TRUE)
  })

  hm_samples <- reactive({
    req(input$inc_groups)
    cs <- common()
    grp <- sample_meta()[cs, "Group"]
    # columns follow the order of the groups in the "Groups to include" box
    cs <- unlist(lapply(input$inc_groups, function(g) cs[grp == g]), use.names = FALSE)
    validate(need(length(cs) > 1, "Select at least two samples."))
    cs
  })

  output$extra_vars_ui <- renderUI({
    req(other_cols())
    selectizeInput("extra_vars", "Extra variables shown as additional annotation boxes",
                   choices = other_cols(), multiple = TRUE)
  })

  group_levels <- reactive(unique(sample_meta()[hm_samples(), "Group"]))

  group_colors <- reactive(group_colors_all()[group_levels()])

  # extra variables: when replicates are merged, a variable that differs inside a group cannot be shown
  extra_dropped <- reactive({
    ev <- input$extra_vars
    if (!isTRUE(input$merge) || !length(ev)) return(character(0))
    md <- sample_meta()[hm_samples(), , drop = FALSE]
    ev[!vapply(ev, function(v) all(tapply(as.character(md[[v]]), md$Group,
                                          function(x) length(unique(x))) == 1), logical(1))]
  })
  extra_use <- reactive(setdiff(input$extra_vars %||% character(0), extra_dropped()))
  output$extra_msg <- renderUI({
    d <- extra_dropped(); if (!length(d)) return(NULL)
    div(class = "alert alert-info", style = "padding:8px; font-size:12px;",
        strong(paste(d, collapse = ", ")), " removed: with merged replicates each group becomes one column, ",
        "but this variable takes different values inside a group, so it cannot be displayed. ",
        "Untick 'Merge replicates' to show it per sample.")
  })

  extra_levels <- reactive({
    ev <- extra_use()
    md <- sample_meta()[hm_samples(), , drop = FALSE]
    setNames(lapply(ev, function(v) unique(as.character(md[[v]]))), ev)
  })

  output$extra_colors_ui <- renderUI({
    el <- extra_levels()
    if (!length(el)) return(NULL)
    pals <- c("Set1", "Dark2", "Paired", "Accent")
    tagList(lapply(seq_along(el), function(k) {
      lv <- el[[k]]; d <- auto_colors(lv, pals[(k - 1) %% 4 + 1])
      tagList(h5(paste("Colours:", names(el)[k])),
              lapply(seq_along(lv), function(i)
                colourInput(paste0("ecol_", k, "_", i), lv[i], d[[i]])))
    }))
  })

  extra_colors <- reactive({
    el <- extra_levels()
    pals <- c("Set1", "Dark2", "Paired", "Accent")
    out <- lapply(seq_along(el), function(k) {
      lv <- el[[k]]; d <- auto_colors(lv, pals[(k - 1) %% 4 + 1])
      setNames(vapply(seq_along(lv), function(i)
        input[[paste0("ecol_", k, "_", i)]] %||% d[[i]], character(1)), lv)
    })
    setNames(out, names(el))
  })

  # size-type inputs are applied once the user stops changing them (like the volcano labels), not at every step
  sz <- throttle(reactive(list(width_cm = input$width_cm, height_cm = input$height_cm, cell_mm = input$cell_mm,
                               leg_height = input$leg_height, leg_width = input$leg_width, leg_font = input$leg_font)), 500)

  core <- reactive({
    cs <- hm_samples()
    mat <- hm_expr()[, cs, drop = FALSE]
    if (identical(input$hm_source, "subset")) {
      sel <- intersect(hm_selected_proteins(), rownames(mat))
      validate(need(length(sel) > 1, "Select gene set(s) with at least 2 proteins."))
      mat <- mat[sel, , drop = FALSE]
    } else if (identical(input$hm_source, "genes")) {
      sel <- intersect(hm_genes(), rownames(mat))                # the user's order is kept
      validate(need(length(sel) > 1, "Select at least 2 proteins (type, or paste a list)."))
      mat <- mat[sel, , drop = FALSE]
    }
    md <- sample_meta()[cs, , drop = FALSE]
    tn <- trimws(input$top_n)
    sn <- input$split_n %||% 2; if (is.na(sn)) sn <- 0
    rs <- if (isTRUE(input$cluster_rows) && identical(input$split_by, "split") && sn > 0) sn else NULL
    rk <- if (isTRUE(input$cluster_rows) && identical(input$split_by, "km") && sn > 0) sn else NULL
    # fit the visible window: use the on-screen size of the plot unless a size is typed in
    dir <- input$hm_dir %||% "normal"
    shp <- if (identical(dir, "circular")) "rect" else (input$hm_shape %||% "rect")
    # the window (the scroll area around the plot, in px) is what everything is fitted into
    win <- input$hm_win
    win_w <- (win$w %||% 1100) / 96 * 25.4
    win_h <- (win$h %||% 620) / 96 * 25.4
    res <- my_proteomics_heat_map_build(
      input_object = mat, sample_metadata = md,
      top_n = {
        tnv <- suppressWarnings(as.numeric(tn))
        if (identical(input$hm_source, "genes") || is.na(tnv) || tnv < 1) "all" else tnv
      },
      merge = input$merge, cluster_rows = input$cluster_rows,
      show_row_dend = isTRUE(input$cluster_rows) && isTRUE(input$show_row_dend),
      dist_param = input$dist_param %||% "euclidean",
      hclust_param = input$hclust_param %||% "average",
      scale_rows = input$scale_rows, show_genes = isTRUE(input$show_genes) || identical(input$hm_source, "genes"),
      group_colors = group_colors(), zscore_limits = input$zlim,
      color_low = input$color_low, color_mid = input$color_mid,
      color_high = input$color_high,
      col_fontsize = if (font_auto$col) NA else (input$col_fontsize %||% NA),
      row_fontsize = if (font_auto$gene) NA else (input$row_fontsize %||% NA),
      group_fontsize = if (font_auto$grp) NA else (input$group_fontsize %||% NA),
      width = if (!size_auto$w && !is.na(sz()$width_cm %||% NA)) sz()$width_cm * 10 else NULL,
      height = if (!size_auto$h && !is.na(sz()$height_cm %||% NA)) sz()$height_cm * 10 else NULL,
      win_w = win_w, win_h = win_h,
      gene_side = input$gene_side, show_group_legend = isTRUE(input$show_group_legend %||% TRUE),
      row_split = rs, row_km = rk,
      extra_cols = extra_use(),
      extra_colors = extra_colors(),
      show_sample_names = isTRUE(input$show_sample_names),
      show_group_names = isTRUE(input$show_group_names),
      leg_width = sz()$leg_width %||% 10, leg_height = if (leg_auto$on) NULL else sz()$leg_height, leg_font = sz()$leg_font %||% 16,
      direction = dir, shape = shp,
      cell_mm = if (shp != "rect" && !cell_auto$on) sz()$cell_mm else NULL,
      circle_by_value = isTRUE(input$circ_by_val),
      circ_opts = list(cluster = isTRUE(input$cluster_rows), dist = input$dist_param, linkage = input$hclust_param,
                       gap = input$circ_gap %||% 8, start = input$circ_start %||% 180))
    if (res$dropped > 0)
      showNotification(sprintf("K-means: %d protein(s) with missing / non-finite values were left out.", res$dropped),
                       type = "warning", duration = 8)
    res
  })

  # gene sets: per-contrast (and custom) ones are pre-selected, Venn-created ones stay available
  is_venn_set <- function(nm) grepl("^(unique|common|discordant)_", nm)
  hm_subset_choices <- reactive(all_subsets())
  output$hm_sets_ui <- renderUI({
    s <- hm_subset_choices(); nm <- names(s)
    lab <- function(x) setNames(x, paste0(x, " (", lengths(s[x]), ")"))
    ch <- list()
    if (any(!is_venn_set(nm))) ch[["Per contrast"]] <- lab(nm[!is_venn_set(nm)])
    if (any(is_venn_set(nm)))  ch[["Venn"]]         <- lab(nm[is_venn_set(nm)])
    sel <- isolate(input$hm_sets)
    sel <- if (is.null(sel)) nm[!is_venn_set(nm)] else intersect(sel, nm)
    selectizeInput("hm_sets", "Gene set(s)", multiple = TRUE, choices = ch, selected = sel,
                   options = list(placeholder = "used when 'selected subset(s)' is chosen"))
  })
  hm_selected_proteins <- reactive(unique(unlist(hm_subset_choices()[input$hm_sets])))

  # draw_heatmap() returns the HeatmapList: renderPlot would print() it and redraw the whole
  # figure, wiping the group-name headers. invisible() prevents that.
  # The plot is always drawn at its canvas size (= the window when auto-fitted). A bigger canvas (typed width / height,
  # or a bigger box size) simply scrolls inside the window.
  canvas_px <- reactive({ cr <- core(); round(c(cr$width_mm, cr$height_mm) / 25.4 * 96) })

  # box / circle size: shows the auto-fitted size; changing it fixes the size (scrolls if bigger than the window)
  # values written by the app itself (auto-fit) are remembered, so their echo from the browser is never mistaken for a user edit
  ph <- new.env()
  prog_add <- function(k, v) ph[[k]] <- utils::tail(c(ph[[k]], v), 20)
  is_prog <- function(k, v, tol = 0.051) length(ph[[k]]) > 0 && any(abs(ph[[k]] - v) < tol)
  cell_auto <- reactiveValues(on = TRUE)
  observe({
    cr <- core(); v <- cr$cell_mm; req(is.finite(v))
    if (cell_auto$on) {
      v <- max(0.05, floor(v * 20) / 20); prog_add("cell", v)             # floor to the slider step: the plot never gets bigger than the window
      updateSliderInput(session, "cell_mm", value = v, min = 0.05, max = max(20, ceiling(v * 2)))
    }
  })
  observeEvent(input$cell_mm, {
    if (!is_prog("cell", input$cell_mm, 0.06)) cell_auto$on <- FALSE
  }, ignoreInit = TRUE)
  observeEvent(input$rst_cell_mm, cell_auto$on <- TRUE)
  output$hm_plot <- renderPlot({
    input$show_group_names; input$show_sample_names     # explicit dependencies: toggles redraw at once
    draw_heatmap(core()); invisible(NULL)
  }, res = 96, width = function() canvas_px()[1], height = function() canvas_px()[2])

  # circular: these options do not apply
  # (they are also unticked while greyed out, and get their old state back afterwards)
  circ_saved <- NULL
  observeEvent(input$hm_dir, {
    circ <- identical(input$hm_dir, "circular")
    ids <- c("show_sample_names", "show_group_legend", "show_row_dend")
    if (circ && is.null(circ_saved)) {
      circ_saved <<- lapply(ids, function(i) isTRUE(input[[i]])); names(circ_saved) <- ids
      for (i in ids) updateCheckboxInput(session, i, value = FALSE)
    } else if (!circ && !is.null(circ_saved)) {
      for (i in ids) updateCheckboxInput(session, i, value = circ_saved[[i]])
      circ_saved <<- NULL
    }
    session$sendCustomMessage("set_dis", list(ids = as.list(ids), dis = circ))
  }, ignoreNULL = FALSE)

  # Z-score bar: default = height of the heatmap body (until the user types / drags a value; "auto" goes back)
  leg_auto <- reactiveValues(on = TRUE)
  leg_prog <- reactiveVal(NA_real_)
  body_h_rv <- reactiveVal(NULL)                         # only changes when the body height really changes
  observe({
    h <- core()$body_h_mm
    if (is.finite(h) && !identical(round(h), isolate(body_h_rv()))) body_h_rv(round(h))
  })
  leg_bump <- reactiveVal(0)                             # a reset re-creates the slider at the default
  output$leg_height_ui <- renderUI({
    h <- body_h_rv(); leg_bump(); req(h)
    mx <- max(60, h)
    v <- if (isolate(leg_auto$on)) h else min(mx, isolate(input$leg_height) %||% h)
    prog_add("leg", v)
    sliderInput("leg_height", NULL, 15, mx, v, step = 1)
  })
  output$leg_rst_ui <- renderUI(rst_btn("leg_height", disabled = leg_auto$on))
  observeEvent(input$leg_height, {
    if (!is_prog("leg", input$leg_height, 0.6)) leg_auto$on <- FALSE
  }, ignoreInit = TRUE)
  observeEvent(input$rst_leg_height, { leg_auto$on <- TRUE; leg_bump(leg_bump() + 1) })

  # ---- figure size: the boxes show the auto-fitted size in cm; typing a value overrides, reset goes back ----
  size_auto <- reactiveValues(w = TRUE, h = TRUE)
  size_prog <- reactiveValues(w = NA_real_, h = NA_real_)
  observe({
    cr <- core()
    if (size_auto$w) { v <- round(cr$width_mm / 10, 1);  prog_add("w", v); updateNumericInput(session, "width_cm", value = v) }
    if (size_auto$h) { v <- round(cr$height_mm / 10, 1); prog_add("h", v); updateNumericInput(session, "height_cm", value = v) }
  })
  observeEvent(input$width_cm, {
    v <- input$width_cm; if (is.null(v) || is.na(v)) return()
    if (!is_prog("w", v)) size_auto$w <- FALSE
  }, ignoreInit = TRUE)
  observeEvent(input$height_cm, {
    v <- input$height_cm; if (is.null(v) || is.na(v)) return()
    if (!is_prog("h", v)) size_auto$h <- FALSE
  }, ignoreInit = TRUE)
  observeEvent(input$rst_width_cm, size_auto$w <- TRUE)
  observeEvent(input$rst_height_cm, size_auto$h <- TRUE)

  result_table <- reactive({
    cr <- core()
    pdf(NULL, width = cr$width_mm / 25.4, height = cr$height_mm / 25.4)
    on.exit(dev.off())
    d <- extract_clusters(draw_heatmap(cr), cr$main_matrix, cr$direction)
    if (!isTRUE(input$cluster_rows)) d$Cluster <- NULL          # a cluster column only exists when rows are clustered
    d
  })
  output$mat_tbl <- DT::renderDT({
    d <- result_table()
    out <- cbind(rownames(d), d); names(out)[1] <- id_label()
    DT::datatable(out, rownames = FALSE, selection = "none", class = "compact stripe nowrap",
                  options = list(scrollX = TRUE, scrollY = "52vh", pageLength = 50, deferRender = TRUE,
                                 lengthMenu = c(25, 50, 100, 500))) %>%
      DT::formatRound(columns = setdiff(seq_len(ncol(out)), c(1, which(names(out) == "Cluster"))), digits = 3)
  })

  file_stub <- reactive(paste0("Heatmap_top", core()$top_n, if (core()$merge) "_merged" else ""))

  output$dl_png <- downloadHandler(
    filename = function() paste0(file_stub(), ".png"),
    content = function(file) {
      cr <- core()
      png(file, width = cr$width_mm, height = cr$height_mm, units = "mm", res = 300)
      on.exit(dev.off())
      draw_heatmap(cr)
    })
  output$dl_csv <- downloadHandler(
    filename = function() paste0(file_stub(), "_matrix.csv"),
    content = function(file) {
      d <- result_table()
      out <- cbind(rownames(d), d); names(out)[1] <- id_label()
      write.csv(out, file, row.names = FALSE)
    })


  observeEvent(input$cluster_rows, {
    if (isTRUE(input$cluster_rows)) updateCheckboxInput(session, "show_row_dend", value = TRUE)
    if (!isTRUE(input$cluster_rows)) {
      updateSelectInput(session, "split_by", selected = "none")
      updateNumericInput(session, "split_n", value = 2)
    }
  }, ignoreInit = TRUE)

# ---- data source: covariate-corrected whenever the DE design had covariates (Batch and / or SVs) ----
hm_expr <- reactive({
  ce <- tryCatch(corrected_expr(), error = function(e) NULL)
  if (!is.null(ce)) ce else final_expr()
})
# specific proteins: significant proteins of ALL ordinary contrasts, largest |t| first
hm_lists <- reactive(ordered_lists(de_res()$single_names))
hi <- picker_inputs(hm_lists)
hm_picker <- protein_picker_server("pp_hm", sig = hi$sig, all = hi$all)
hm_genes <- hm_picker$selected

# ---- label fonts: the boxes show the auto-fitted size; typing a size switches to manual, "auto-fit" goes back ----
font_auto <- reactiveValues(col = TRUE, grp = TRUE, gene = TRUE)
font_prog <- reactiveValues(col = NA_real_, grp = NA_real_, gene = NA_real_)           # the value we wrote ourselves
observe({
  cr <- core()
  if (font_auto$gene && is.finite(cr$fonts$gene)) { v <- round(cr$fonts$gene, 1); prog_add("fg", v); updateNumericInput(session, "row_fontsize", value = v) }
  if (font_auto$col && is.finite(cr$fonts$sample)) { v <- round(cr$fonts$sample, 1); prog_add("fc", v); updateNumericInput(session, "col_fontsize", value = v) }
  if (font_auto$grp && is.finite(cr$fonts$group)) { v <- round(cr$fonts$group, 1);  prog_add("fgr", v); updateNumericInput(session, "group_fontsize", value = v) }
})
observeEvent(input$col_fontsize, {
  v <- input$col_fontsize; if (is.null(v) || is.na(v)) return()
  if (!is_prog("fc", v)) font_auto$col <- FALSE
}, ignoreInit = TRUE)
observeEvent(input$group_fontsize, {
  v <- input$group_fontsize; if (is.null(v) || is.na(v)) return()
  if (!is_prog("fgr", v)) font_auto$grp <- FALSE
}, ignoreInit = TRUE)
observeEvent(input$row_fontsize, {
  v <- input$row_fontsize; if (is.null(v) || is.na(v)) return()
  if (!is_prog("fg", v)) font_auto$gene <- FALSE
}, ignoreInit = TRUE)
observeEvent(input$rst_row_fontsize, font_auto$gene <- TRUE)
observeEvent(input$rst_col_fontsize, font_auto$col <- TRUE)
observeEvent(input$rst_group_fontsize, font_auto$grp <- TRUE)

# ---- reset buttons: enabled only when the value differs from its default ----
rst_send <- function(id, changed) session$sendCustomMessage("rst_state", list(id = id, dis = !changed))
observe({                       # auto-fitted boxes: changed = the user typed a value
  rst_send("col_fontsize", !font_auto$col); rst_send("row_fontsize", !font_auto$gene)
  rst_send("group_fontsize", !font_auto$grp)
  rst_send("width_cm", !size_auto$w); rst_send("height_cm", !size_auto$h)
  rst_send("cell_mm", !cell_auto$on)
})
rst_static <- list(                 # fixed defaults
  zlim = list("slider", c(-2, 2)), leg_width = list("slider", 10), leg_font = list("slider", 16),
  circ_gap = list("slider", 8), circ_start = list("slider", 180),
  top_n = list("text", ""), split_n = list("numeric", 2),
  color_low = list("colour", "blue"), color_mid = list("colour", "white"), color_high = list("colour", "red"))
same_val <- function(a, b, type) {
  if (is.null(a)) return(TRUE)
  if (type == "colour") return(tryCatch(identical(col2rgb(a), col2rgb(b)), error = function(e) FALSE))
  if (type == "text") return(identical(trimws(as.character(a)), b))
  isTRUE(all.equal(as.numeric(a), as.numeric(b)))
}
for (id in names(rst_static)) local({
  id <- id; type <- rst_static[[id]][[1]]; def <- rst_static[[id]][[2]]
  observeEvent(input[[paste0("rst_", id)]], {
    switch(type,
      slider  = updateSliderInput(session, id, value = def),
      numeric = updateNumericInput(session, id, value = def),
      text    = updateTextInput(session, id, value = def),
      colour  = colourpicker::updateColourInput(session, id, value = def))
  })
  observeEvent(input[[id]], rst_send(id, !same_val(input[[id]], def, type)), ignoreNULL = FALSE)
})

# gene-name position: sideways = below / above the heatmap, normal = right / left
output$gene_side_ui <- renderUI({
  req(isTRUE(input$show_genes))
  d <- input$hm_dir %||% "normal"
  if (d == "circular") return(NULL)
  if (d == "sideways") selectInput("gene_side", "Gene names position", c("down", "up"), selected = "down")
  else selectInput("gene_side", "Gene names position", c("right", "left"), selected = "right")
})

# ---- "Open in Morpheus": the processed matrix goes to a secret Gist (GITHUB_PAT needed, like ORA-cle),
#      then Morpheus is opened with ?dataset=<raw url>. GCT 1.3 carries the sample annotations (Group + extra variables).
gct13_text <- function(m, anno) {
  num <- function(x) ifelse(is.finite(x), sprintf("%.5g", x), "")
  rn <- make.unique(rownames(m))
  body <- vapply(seq_len(nrow(m)), function(i) paste(c(rn[i], rn[i], num(m[i, ])), collapse = "\t"), character(1))
  ann <- vapply(names(anno), function(a) paste(c(a, "", as.character(anno[[a]])), collapse = "\t"), character(1))
  paste(c("#1.3", paste(nrow(m), ncol(m), 1, ncol(anno), sep = "\t"),
          paste(c("id", "Description", colnames(m)), collapse = "\t"), ann, body), collapse = "\n")
}
create_gist_raw <- function(txt, filename, description) {
  tk <- gist_token(); if (!nzchar(tk)) return(NULL)
  resp <- tryCatch(httr::POST("https://api.github.com/gists",
      httr::add_headers(Authorization = paste("Bearer", tk), Accept = "application/vnd.github+json"),
      body = list(description = description, public = FALSE,
                  files = setNames(list(list(content = txt)), filename)),
      encode = "json", httr::timeout(30)), error = function(e) NULL)
  if (is.null(resp) || httr::status_code(resp) != 201) return(NULL)
  httr::content(resp)$files[[1]]$raw_url
}
observeEvent(input$morph_open, {
  cr <- tryCatch(core(), error = function(e) NULL)
  if (is.null(cr)) { session$sendCustomMessage("morph_fail", ""); return(showNotification("The heatmap is not ready yet.", type = "warning")) }
  m <- cr$main_matrix
  anno <- if (isTRUE(cr$merge)) data.frame(Group = colnames(m))
          else sample_meta()[colnames(m), unique(c("Group", extra_use())), drop = FALSE]
  anno[] <- lapply(anno, as.character)
  if (!nzchar(gist_token())) {
    session$sendCustomMessage("morph_fail", "")
    return(showNotification("Morpheus needs a public link to the data, which needs a GITHUB_PAT (token with 'gist' scope) - see the instructions. Meanwhile use 'Download GCT' and drag the file into Morpheus (File > Open).",
                            type = "error", duration = 15))
  }
  raw <- create_gist_raw(gct13_text(m, anno), "proteus_matrix.gct", "Pro-teus export for Morpheus")
  if (is.null(raw)) {
    session$sendCustomMessage("morph_fail", "")
    return(showNotification("Could not create the Gist (check the GITHUB_PAT token).", type = "error", duration = 10))
  }
  session$sendCustomMessage("morph_go", paste0("https://software.broadinstitute.org/morpheus/?json=",
    utils::URLencode(as.character(jsonlite::toJSON(list(dataset = raw), auto_unbox = TRUE)), reserved = TRUE)))
})

output$dl_gct <- downloadHandler(
  filename = function() "proteus_matrix.gct",
  content = function(file) {
    cr <- core(); m <- cr$main_matrix
    anno <- if (isTRUE(cr$merge)) data.frame(Group = colnames(m))
            else sample_meta()[colnames(m), unique(c("Group", extra_use())), drop = FALSE]
    anno[] <- lapply(anno, as.character)
    writeLines(gct13_text(m, anno), file)
  })
