# ======================= BACKGROUND JOBS + STOP BUTTON =======================
# Heavy tasks (ORA, imputation, NormalyzerDE) run in a separate R process (callr)
# so that the "Stop process" button can really kill them.
# per-session scratch folder for analysis outputs; removed when the session ends
session_tmp <- file.path(tempdir(), paste0("tk_", session$token))
dir.create(session_tmp, recursive = TRUE, showWarnings = FALSE)
session$onSessionEnded(function() unlink(session_tmp, recursive = TRUE))

job <- reactiveValues(proc = NULL, label = NULL, log = NULL, out = NULL, stopped = FALSE)
job_cb <- NULL   # callback to run when the job finishes (kept outside reactiveValues)
job_fail_cb <- NULL   # optional: called instead of an error message when the job (or its callback) fails

start_job <- function(label, fun, args, on_done, on_fail = NULL) {
  if (!is.null(job$proc) && job$proc$is_alive()) {
    showNotification("Another process is still running. Stop it first.", type = "warning")
    return(invisible(FALSE))
  }
  logf <- tempfile(fileext = ".log"); outf <- tempfile(fileext = ".rds")
  src <- normalizePath(sort(list.files("functions", pattern = "\\.R$", full.names = TRUE)))
  proc <- callr::r_bg(
    function(fun, args, src, outf) {
      suppressPackageStartupMessages({
        library(ggplot2); library(dplyr); library(tidyr); library(scales); library(ggrepel); library(limma)
      })
      for (f in src) source(f)
      res <- do.call(fun, args)
      saveRDS(res, outf)
      invisible(TRUE)
    },
    args = list(fun = fun, args = args, src = src, outf = outf),
    stdout = logf, stderr = "2>&1", supervise = TRUE)
  job$proc <- proc; job$label <- label; job$log <- logf; job$out <- outf; job$stopped <- FALSE
  job_cb <<- on_done; job_fail_cb <<- on_fail
  invisible(TRUE)
}

job_log_tail <- function(n = 1) {
  if (is.null(job$log) || !file.exists(job$log)) return("")
  x <- tryCatch(readLines(job$log, warn = FALSE), error = function(e) character())
  x <- trimws(x[nzchar(trimws(x))])
  if (!length(x)) "" else paste(tail(x, n), collapse = " | ")
}

observe({
  req(job$proc)
  invalidateLater(1000)
  p <- job$proc
  if (!p$is_alive()) isolate({
    lbl <- job$label; outf <- job$out; cb <- job_cb; stopped <- job$stopped; fcb <- job_fail_cb
    job$proc <- NULL
    if (stopped) {
      showNotification(paste(lbl, "stopped."), type = "warning")
    } else if (identical(p$get_exit_status(), 0L) && file.exists(outf)) {
      tryCatch({ cb(readRDS(outf)); showNotification(paste(lbl, "finished."), type = "message") },
               error = function(e) {
                 if (is.function(fcb)) fcb(conditionMessage(e))
                 else showNotification(paste0(lbl, ": ", conditionMessage(e)), type = "error", duration = NULL)
               })
    } else {
      if (is.function(fcb)) fcb(job_log_tail(3))
      else showNotification(paste0(lbl, " failed: ", job_log_tail(3)), type = "error", duration = NULL)
    }
  })
})

observeEvent(input$stop_proc, {
  if (!is.null(job$proc) && job$proc$is_alive()) {
    job$stopped <- TRUE
    job$proc$kill_tree()
  } else showNotification("No process is running.", type = "message")
})

output$job_status <- renderUI({
  if (is.null(job$proc)) return(tags$span(style = "color:#777;", "Idle"))
  invalidateLater(1500)
  tags$span(icon("spinner", class = "fa-spin"), strong(paste0(" ", job$label, "... ")),
            tags$small(style = "color:#555;", substr(job_log_tail(1), 1, 70)))
})
