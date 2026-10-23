ORA_SITE <- "https://symeon-d.github.io/ora-cle/"

# opens the ORA-cle in a tab opened synchronously on click (avoids pop-up blockers), then
# the server sends the final URL once the Gist has been created
ORA_JS <- "
function oraOpen(i) {
  window._oraWin = window.open('about:blank', '_blank');
  Shiny.setInputValue('ora_open', {i: i, t: Date.now()}, {priority: 'event'});
}
$(function() {
  Shiny.addCustomMessageHandler('ora_go', function(url) {
    if (window._oraWin && !window._oraWin.closed) window._oraWin.location.href = url;
    else window.open(url, '_blank');
  });
  Shiny.addCustomMessageHandler('ora_fail', function(msg) {
    if (window._oraWin && !window._oraWin.closed) window._oraWin.close();
  });
});"

ora_main_panel <- function() {
  tabsetPanel(id = "ora_tabs",
    tabPanel("Result files", br(), uiOutput("ora_files_ui"),
             downloadButton("dl_ora_csvs", "Download All CSVs"),
             downloadButton("dl_ora_zip", "Download everything (CSVs + plots, ZIP)")),
    tabPanel("Quick view", br(), uiOutput("ora_quick_ui"),
             plotOutput("ora_quick_plot", height = "auto"),
             downloadButton("dl_ora_quick", "Download PNG")),
    tabPanel("Multi-cluster plot", br(), uiOutput("ora_mc_ui"), uiOutput("mc_sets_ui"),
             uiOutput("mc_shapes_ui"),
             plotOutput("ora_mc_plot", height = "auto"),
             downloadButton("dl_ora_mc", "Download PNG")),
    tabPanel("Simplified", br(), uiOutput("simp_ui")))
}

gsea_main_panel <- function() {
  tabsetPanel(id = "gsea_tabs",
    tabPanel("Result files", br(), uiOutput("gsea_files_ui"),
             downloadButton("dl_gsea_csvs", "Download All CSVs")),
    tabPanel("NES plot", br(), uiOutput("gsea_nes_ctl"),
             plotOutput("gsea_nes_plot", height = "auto"), downloadButton("dl_gsea_nes", "Download PNG")),
    tabPanel("Compare contrasts", br(), uiOutput("gsea_cmp_ctl"),
             plotOutput("gsea_cmp_plot", height = "auto"), downloadButton("dl_gsea_cmp", "Download PNG")))
}

tab_ora <- function() {
  has <- function(m) sprintf("input.enrich_modules && input.enrich_modules.indexOf('%s') > -1", m)
  tabPanel("6. Pathway enrichment", tags$script(HTML(ORA_JS)), sidebarLayout(
    sidebarPanel(width = 4, style = "max-height: 90vh; overflow-y: auto;",
      uiOutput("ora_avail"),
      checkboxGroupInput("enrich_modules", "Analysis", c("ORA" = "ora", "GSEA" = "gsea"),
                         selected = "ora", inline = TRUE),
      selectInput("ora_org", "Organism",
                  c("human (org.Hs.eg.db)" = "org.Hs.eg.db", "mouse (org.Mm.eg.db)" = "org.Mm.eg.db",
                    "rat (org.Rn.eg.db)" = "org.Rn.eg.db", "C. elegans (org.Ce.eg.db)" = "org.Ce.eg.db"),
                  selected = "org.Mm.eg.db"),
      # ---------------- ORA settings ----------------
      conditionalPanel(has("ora"),
        h4("ORA"),
        uiOutput("ora_sets_ui"),
        radioButtons("ora_bg", "Background",
                     c("None (tool default)" = "none", "All proteins in the expression file" = "all"),
                     selected = "none"),
        h5("Methods"),
        checkboxGroupInput("ora_method", NULL, c("gprofiler", "genecodis", "clusterprofiler"),
                           selected = "gprofiler", inline = TRUE),
        checkboxGroupInput("ora_ont", "Ontologies", c("GO:BP", "GO:MF", "GO:CC"),
                           selected = "GO:BP", inline = TRUE),
        conditionalPanel("input.ora_method && input.ora_method.indexOf('gprofiler') > -1",
          checkboxGroupInput("ora_corr", "Correction (g:Profiler)", c("gSCS", "fdr", "bonferroni"),
                             selected = "gSCS", inline = TRUE)),
        conditionalPanel("input.ora_method && input.ora_method.indexOf('clusterprofiler') > -1",
          checkboxGroupInput("ora_cp_corr", "Correction (clusterProfiler)",
                             c("BH", "bonferroni", "holm", "hochberg", "hommel", "BY"),
                             selected = "BH", inline = TRUE)),
        fluidRow(
          column(4, numericInput("ora_min_genes", "Min genes", 3, min = 1)),
          column(4, numericInput("ora_min_ts", "Min term size", 10, min = 1)),
          column(4, numericInput("ora_max_ts", "Max term size", NA, min = 1))),
        helpText("Blank = no limit."),
        numericInput("ora_top_n", "Top N terms (blank = all)", NA),
        checkboxInput("ora_simplify", "Simplify with rrvgo", FALSE),
        conditionalPanel("input.ora_simplify",
          numericInput("ora_rrvgo", "rrvgo threshold", 0.7, min = 0, max = 1, step = 0.05)),
        checkboxInput("ora_cp_simplify", "simplify with clusterProfiler::simplify", FALSE),
        conditionalPanel("input.ora_cp_simplify",
          numericInput("ora_cp_cutoff", "simplify cutoff", 0.7, min = 0, max = 1, step = 0.05)),
        actionButton("run_ora", "Run ORA", class = "btn-success")),
      # ---------------- GSEA settings ----------------
      conditionalPanel(has("gsea"),
        hr(), h4("GSEA"),
        helpText("Ranks proteins by the moderated t-statistic of each contrast (fgsea); GO gene sets come from the organism database."),
        uiOutput("gsea_contrasts_ui"),
        checkboxGroupInput("gsea_ont", "Ontologies", c("GO:BP" = "BP", "GO:MF" = "MF", "GO:CC" = "CC"),
                           selected = "BP", inline = TRUE),
        fluidRow(
          column(4, numericInput("gsea_min", "Min size", 10, min = 1)),
          column(4, numericInput("gsea_max", "Max size", 500, min = 1)),
          column(4, numericInput("gsea_nperm", "Permutations", 10000, min = 100, step = 1000))),
        numericInput("gsea_padj", "Keep padj <= (blank = keep all)", 0.05, min = 0, max = 1, step = 0.01),
        actionButton("run_gsea", "Run GSEA", class = "btn-success"))
    ),
    mainPanel(uiOutput("enrich_tabs"))
  ))
}
