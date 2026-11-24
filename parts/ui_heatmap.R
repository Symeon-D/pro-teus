MORPH_JS <- "
function morphOpen() {
  window._morphWin = window.open('about:blank', '_blank');
  Shiny.setInputValue('morph_open', {t: Date.now()}, {priority: 'event'});
}
$(function() {
  var lastWin = '';
  setInterval(function() {                              // size of the heatmap window (rounded, so scrollbars do not cause jitter)
    var el = document.getElementById('hm_win');
    if (!el || el.clientWidth < 50 || !window.Shiny || !Shiny.setInputValue) return;
    var w = Math.round(el.clientWidth / 20) * 20, h = Math.round(el.clientHeight / 20) * 20, k = w + 'x' + h;
    if (k !== lastWin) { lastWin = k; Shiny.setInputValue('hm_win', {w: w - 24, h: h - 24}); }
  }, 400);
  Shiny.addCustomMessageHandler('morph_go', function(url) {
    if (window._morphWin && !window._morphWin.closed) window._morphWin.location.href = url;
    else window.open(url, '_blank');
  });
  Shiny.addCustomMessageHandler('set_dis', function(m) {
    m.ids.forEach(function(i) { $('#' + i).prop('disabled', m.dis).closest('.checkbox').css('opacity', m.dis ? 0.45 : 1); });
  });
  Shiny.addCustomMessageHandler('rst_state', function(m) {
    $('#rst_' + m.id).prop('disabled', m.dis);
  });
  Shiny.addCustomMessageHandler('morph_fail', function(msg) {
    if (window._morphWin && !window._morphWin.closed) window._morphWin.close();
  });
});"

tab_heatmap <- function() {
  tabPanel("6. Heatmap", tags$script(HTML(MORPH_JS)), sidebarLayout(
    sidebarPanel(width = 4, style = "max-height: 90vh; overflow-y: auto;",
      div(style = "display:flex; justify-content:space-between; align-items:center;",
          h4("Proteins & samples"),
          tags$button(type = "button", class = "btn btn-xs btn-default", onclick = "morphOpen()",
                      title = "Opens the processed matrix in the Broad Institute's Morpheus (new tab)",
                      "open in Morpheus")),
      radioButtons("hm_source", NULL,
                   c("All proteins in the expression file" = "all",
                     "Proteins from the selected subset(s)" = "subset",
                     "Specific proteins" = "genes"),
                   selected = "all"),
      uiOutput("hm_sets_ui"),
      conditionalPanel("input.hm_source == 'genes'", protein_picker_ui("pp_hm", "Proteins (order kept when rows are not clustered)")),
      uiOutput("inc_groups_ui"),
      uiOutput("extra_vars_ui"),
      uiOutput("extra_msg"),
      uiOutput("extra_colors_ui"),
      h4("Heatmap"),
      selectInput("hm_dir", "Direction",
                  c("Normal (proteins as rows)" = "normal",
                    "Sideways (samples as rows)" = "sideways",
                    "Circular (proteins around the circle)" = "circular")),
      conditionalPanel("input.hm_dir != 'circular'",
        selectInput("hm_shape", "Cell shape",
                    c("Rectangles (fill the window)" = "rect",
                      "Square boxes (fit the window)" = "square",
                      "Circles (fit the window)" = "circle")),
        conditionalPanel("input.hm_shape != 'rect'",
          lab_rst("Box / circle size (mm)", "cell_mm"),
          sliderInput("cell_mm", NULL, 0.05, 20, 6, step = 0.05),
          conditionalPanel("input.hm_shape == 'circle'",
            checkboxInput("circ_by_val", "Circle size reflects |value| (diameter 20-100% of the cell)", FALSE)))),
      conditionalPanel("input.hm_dir == 'circular'",
        lab_rst("Opening of the circle (degrees) - where the ring titles are written", "circ_gap"),
        sliderInput("circ_gap", NULL, 2, 90, 8, step = 1),
        lab_rst("Position of the opening (start angle, degrees)", "circ_start"),
        sliderInput("circ_start", NULL, 0, 360, 180, step = 5)),
      checkboxInput("merge", "Merge replicates (group means)", FALSE),
      lab_rst("Top N variable proteins (blank = all)", "top_n"),
      textInput("top_n", NULL, "", placeholder = "all"),
      checkboxInput("scale_rows", "Z-score rows", TRUE),
      checkboxInput("show_sample_names", "Show sample names", TRUE),
      checkboxInput("show_group_names", "Show group names", TRUE),
      checkboxInput("show_genes", "Show gene names", FALSE),
      uiOutput("gene_side_ui"),
      checkboxInput("show_group_legend", "Show group legend", TRUE),
      checkboxInput("cluster_rows", "Cluster rows", FALSE),
      conditionalPanel("input.cluster_rows",
        checkboxInput("show_row_dend", "Show row dendrogram", TRUE),
        selectInput("dist_param", "Distance",
                    c("euclidean", "manhattan", "maximum", "canberra", "minkowski")),
        selectInput("hclust_param", "Linkage",
                    c("average", "complete", "single", "ward.D", "ward.D2", "mcquitty")),
        selectInput("split_by", "Split rows by",
                    c("None" = "none", "Row split (slices)" = "split", "K-means" = "km")),
        conditionalPanel("input.split_by != 'none'",
          lab_rst("Number of clusters", "split_n"),
          numericInput("split_n", NULL, 2, min = 2, step = 1))),
      lab_rst("Z-score limits", "zlim"),
      sliderInput("zlim", NULL, -6, 6, c(-2, 2), step = 0.5),
      fluidRow(
        column(4, lab_rst("Low", "color_low"), colourInput("color_low", NULL, "blue")),
        column(4, lab_rst("Mid", "color_mid"), colourInput("color_mid", NULL, "white")),
        column(4, lab_rst("High", "color_high"), colourInput("color_high", NULL, "red"))),
      h4("Appearance"),
      lab_rst("Sample label font", "col_fontsize"),
      numericInput("col_fontsize", NULL, NA, min = 1, width = "130px"),
      lab_rst("Gene label font", "row_fontsize"),
      numericInput("row_fontsize", NULL, NA, min = 1, width = "130px"),
      lab_rst("Group label font", "group_fontsize"),
      numericInput("group_fontsize", NULL, NA, min = 1, width = "130px"),
      h5("Z-score scale bar"),
      lab_rst("Width (mm)", "leg_width"),
      sliderInput("leg_width", NULL, 3, 30, 10, step = 1),
      div(style = "display:flex; align-items:center; gap:8px; margin-bottom:-6px;",
          tags$label("Height (mm)", style = "margin:0;"), uiOutput("leg_rst_ui", inline = TRUE)),
      uiOutput("leg_height_ui"),
      lab_rst("Font size (pt)", "leg_font"),
      sliderInput("leg_font", NULL, 6, 30, 16, step = 1),
      h5("Figure size"),
      fluidRow(
        column(6, lab_rst("Width (cm)", "width_cm"), numericInput("width_cm", NULL, NA, min = 3, step = 0.5)),
        column(6, lab_rst("Height (cm)", "height_cm"), numericInput("height_cm", NULL, NA, min = 3, step = 0.5))),
      h4("Export"),
      downloadButton("dl_png", "PNG"),
      downloadButton("dl_gct", "GCT (for Morpheus)")
    ),
    mainPanel(
      tabsetPanel(
        tabPanel("Heatmap", br(),
                 # the scroll area: everything is fitted into it; it reports its size (px) to the server
                 div(id = "hm_win", style = "width:100%; height:calc(100vh - 180px); overflow:auto; border:1px solid #e3e3e3;",
                     tags$style("#hm_plot {display:block;}"), plotOutput("hm_plot", width = "auto", height = "auto"))),
        tabPanel("Processed matrix", br(),
                 downloadButton("dl_csv", "Download matrix + clusters (CSV)", class = "btn-sm"), br(), br(),
                 DT::DTOutput("mat_tbl")))
    )
  ))
}
