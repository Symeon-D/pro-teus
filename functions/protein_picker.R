# ==========================================================
# Reusable protein selector (Shiny module)
#  - type to search, or PASTE a list (spaces / commas / semicolons / new lines / a column copied from a
#    limma table): matching names are selected automatically
#  - the list holds the significant proteins ordered by |t|; "include non significant" adds the others
#  - click a selected protein (it turns blue) and press Ctrl/Cmd+C to copy it; with nothing clicked,
#    Ctrl/Cmd+C copies all selected proteins (one per line)
#  - "Clear all" empties the selection
# Returns list(selected = reactive, set = function(values))
# ==========================================================
PP_JS <- "
function ppFallback(txt) {
  var ta = document.createElement('textarea'); ta.value = txt; ta.style.position = 'fixed'; ta.style.opacity = '0';
  document.body.appendChild(ta); ta.select(); try { document.execCommand('copy'); } catch (e) {} document.body.removeChild(ta);
}
function ppWrite(txt) {
  if (navigator.clipboard && window.isSecureContext) navigator.clipboard.writeText(txt).catch(function() { ppFallback(txt); });
  else ppFallback(txt);
}
function ppValues(self) {   // the highlighted (blue) items, or every selected item when none is highlighted
  var els = self.$control.children('.item.active'); if (!els.length) els = self.$control.children('.item');
  var v = []; els.each(function() { v.push($(this).attr('data-value')); }); return v;
}
if (!window._ppInit) {
  window._ppInit = true; window._ppActive = null;
  // Ctrl/Cmd+C copies the highlighted proteins of the selector that was clicked last
  document.addEventListener('keydown', function(e) {
    if (!((e.ctrlKey || e.metaKey) && (e.key === 'c' || e.key === 'C'))) return;
    var s = window._ppActive; if (!s) return;
    if (window.getSelection && window.getSelection().toString()) return;
    var el = document.activeElement; if (el && el.tagName === 'INPUT' && el.value) return;
    var v = ppValues(s); if (!v.length) return;
    e.preventDefault(); ppWrite(v.join('\\n'));
    Shiny.setInputValue(s._ppId + '_copied', {n: v.length, t: Date.now()}, {priority: 'event'});
  }, true);
  document.addEventListener('mousedown', function(e) {      // clicking elsewhere ends the copy target
    if (!(e.target.closest && e.target.closest('.pp-picker'))) window._ppActive = null;
  }, true);
}"

protein_picker_ui <- function(id, label = "Proteins") {
  ns <- shiny::NS(id); sel_id <- ns("sel")
  init <- sprintf("function() {
      var self = this, inputId = '%s', ctrl = self.$control[0], last = null;
      self._ppId = inputId;
      this.$control_input.on('paste', function(e) {
        var t = ((e.originalEvent || e).clipboardData || window.clipboardData).getData('text');
        if (t && /[\\s,;]/.test(t.trim())) {
          e.preventDefault();
          Shiny.setInputValue(inputId + '_paste', {text: t, t: Date.now()}, {priority: 'event'});
        }
      });
      // Ctrl/Cmd+click toggles one protein, Shift+click selects a range; a plain click highlights just one
      ctrl.addEventListener('mousedown', function(e) {
        window._ppActive = self;
        var it = e.target.closest ? e.target.closest('.item') : null;
        if (!it || !ctrl.contains(it) || e.target.closest('.remove')) return;
        var items = Array.prototype.slice.call(ctrl.querySelectorAll('.item'));
        if (e.shiftKey && last && items.indexOf(last) > -1) {
          e.preventDefault(); e.stopImmediatePropagation();
          var i1 = items.indexOf(last), i2 = items.indexOf(it), lo = Math.min(i1, i2), hi = Math.max(i1, i2);
          items.forEach(function(x, i) { x.classList.toggle('active', i >= lo && i <= hi); });
        } else if (e.ctrlKey || e.metaKey) {
          e.preventDefault(); e.stopImmediatePropagation();
          it.classList.toggle('active'); last = it;
        } else { last = it; return; }
        self.$activeItems = Array.prototype.slice.call(ctrl.querySelectorAll('.item.active'));
      }, true);
    }", sel_id)
  shiny::tagList(
    shiny::tags$script(shiny::HTML(PP_JS)),
    # tells the server that the selector was (re-)created, so the choices are sent again
    shiny::tags$script(shiny::HTML(sprintf("Shiny.setInputValue('%s', Date.now(), {priority: 'event'});", ns("ready")))),
    shiny::div(style = "display:flex; gap:10px; align-items:flex-end;",
      shiny::div(class = "pp-picker", style = "flex:1; min-width:0;",
        shiny::selectizeInput(sel_id, label, choices = NULL, multiple = TRUE, width = "100%",
          options = list(placeholder = "type a name, or paste a list (spaces, commas, new lines, a copied column)...",
                         plugins = list("remove_button"), maxOptions = 300, sortField = "$order", onInitialize = I(init)))),
      shiny::div(style = "margin-bottom:15px; white-space:nowrap;", shiny::checkboxInput(ns("show_all"), "include non significant", FALSE))),
    shiny::div(style = "margin:-8px 0 8px 0;",
      shiny::actionButton(ns("clear"), "Clear all", icon = shiny::icon("eraser"), class = "btn-xs"),
      shiny::tags$small(style = "margin-left:8px; color:#777;", shiny::textOutput(ns("note"), inline = TRUE))))
}

# sig / all: reactives returning character vectors, ALREADY ordered (largest |t| first);
# `all` holds every protein (significant ones first)
protein_picker_server <- function(id, sig, all) {
  sig_fn <- sig; all_fn <- all                      # (the argument `all` would otherwise mask base::all)
  shiny::moduleServer(id, function(input, output, session) {
    note <- shiny::reactiveVal(""); pending <- shiny::reactiveVal(NULL)
    base_choices <- shiny::reactive({
      s <- tryCatch(sig_fn(), error = function(e) character(0)); a <- tryCatch(all_fn(), error = function(e) character(0))
      if (isTRUE(input$show_all) || !length(s)) a else s
    })
    refresh <- function(selected) {
      ch <- union(selected, base_choices())
      shiny::updateSelectizeInput(session, "sel", choices = ch, selected = selected, server = TRUE)
    }
    shiny::observeEvent(list(base_choices(), input$show_all), {
      keep <- pending(); if (is.null(keep)) keep <- shiny::isolate(input$sel); pending(NULL)
      refresh(keep %||% character(0))
    })
    shiny::observeEvent(input$ready, refresh(shiny::isolate(input$sel) %||% character(0)))
    shiny::observeEvent(input$clear, { note(""); refresh(character(0)) })
    shiny::observeEvent(input$sel_copied, note(sprintf("copied %d protein(s)", input$sel_copied$n)))
    shiny::observeEvent(input$sel_paste, {
      tryCatch({
        a <- all_fn(); txt <- input$sel_paste$text
        toks <- unique(trimws(gsub("[\"']", "", unlist(strsplit(txt, "[[:space:],;]+")))))
        toks <- toks[nzchar(toks)]
        hit <- a[match(toks, a)]                                       # exact first
        miss <- toks[is.na(hit)]
        ci <- a[match(tolower(miss), tolower(a))]                       # then ignoring case
        matched <- unique(c(hit[!is.na(hit)], ci[!is.na(ci)]))
        gone <- miss[is.na(ci)]; gone <- gone[is.na(suppressWarnings(as.numeric(gone)))]   # numbers are not names
        new_sel <- unique(c(shiny::isolate(input$sel), matched))
        note(sprintf("pasted: %d matched%s", length(matched),
                     if (length(gone)) sprintf(", %d not found (%s%s)", length(gone), paste(utils::head(gone, 5), collapse = ", "), if (length(gone) > 5) ", ..." else "") else ""))
        if (!length(matched)) return()
        if (!isTRUE(input$show_all) && !base::all(matched %in% base_choices())) {   # anything may be selected, not only significant ones
          pending(new_sel); shiny::updateCheckboxInput(session, "show_all", value = TRUE)
        } else refresh(new_sel)
      }, error = function(e) note(paste("paste failed:", conditionMessage(e))))
    })
    output$note <- shiny::renderText(note())
    list(selected = shiny::reactive(input$sel %||% character(0)), set = function(values) refresh(values),
         ready = shiny::reactive(input$ready))
  })
}
