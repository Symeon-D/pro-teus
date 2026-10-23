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

  core <- reactive({
    cs <- hm_samples()
    mat <- final_expr()[, cs, drop = FALSE]
    if (identical(input$hm_source, "subset")) {
      sel <- intersect(hm_selected_proteins(), rownames(mat))
      validate(need(length(sel) > 1, "Select gene set(s) with at least 2 proteins."))
      mat <- mat[sel, , drop = FALSE]
    }
    md <- sample_meta()[cs, , drop = FALSE]
    tn <- trimws(input$top_n)
    sn <- input$split_n %||% 2; if (is.na(sn)) sn <- 0
    rs <- if (isTRUE(input$cluster_rows) && identical(input$split_by, "split") && sn > 0) sn else NULL
    rk <- if (isTRUE(input$cluster_rows) && identical(input$split_by, "km") && sn > 0) sn else NULL
    # fit the visible window: use the on-screen size of the plot unless a size is typed in
    scr_w <- session$clientData$output_hm_plot_width
    scr_h <- session$clientData$output_hm_plot_height
    res <- my_proteomics_heat_map_build(
      input_object = mat, sample_metadata = md,
      top_n = if (tolower(tn) == "all" || tn == "") "all" else as.numeric(tn),
      merge = input$merge, cluster_rows = input$cluster_rows,
      show_row_dend = isTRUE(input$cluster_rows) && isTRUE(input$show_row_dend),
      dist_param = input$dist_param %||% "euclidean",
      hclust_param = input$hclust_param %||% "average",
      scale_rows = input$scale_rows, show_genes = input$show_genes,
      group_colors = group_colors(), zscore_limits = input$zlim,
      color_low = input$color_low, color_mid = input$color_mid,
      color_high = input$color_high,
      col_fontsize = input$col_fontsize %||% NA, row_fontsize = input$row_fontsize,
      group_fontsize = input$group_fontsize %||% NA,
      width = if (!is.na(input$width_mm)) input$width_mm else if (!is.null(scr_w)) scr_w / 96 * 25.4,
      height = if (!is.na(input$height_mm)) input$height_mm else if (!is.null(scr_h)) scr_h / 96 * 25.4,
      row_split = rs, row_km = rk,
      extra_cols = extra_use(),
      extra_colors = extra_colors(),
      show_sample_names = isTRUE(input$show_sample_names),
      show_group_names = isTRUE(input$show_group_names),
      leg_width = input$leg_width, leg_height = input$leg_height, leg_font = input$leg_font)
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
  output$hm_plot <- renderPlot({
    input$show_group_names; input$show_sample_names     # explicit dependencies: toggles redraw at once
    draw_heatmap(core()); invisible(NULL)
  }, res = 96)

  # Z-score bar height may grow up to the full plot height
  observe({
    h <- session$clientData$output_hm_plot_height; req(h)
    mm <- floor(h / 96 * 25.4)
    updateSliderInput(session, "leg_height", max = max(60, mm))
  })

  result_table <- reactive({
    cr <- core()
    pdf(NULL, width = cr$width_mm / 25.4, height = cr$height_mm / 25.4)
    on.exit(dev.off())
    extract_clusters(draw_heatmap(cr), cr$main_matrix)
  })
  output$mat_tbl <- renderTable({
    d <- head(result_table(), 200)
    out <- cbind(rownames(d), d); names(out)[1] <- id_label(); out
  }, rownames = FALSE)

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
