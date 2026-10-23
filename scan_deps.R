# scan_deps.R -- run from the app folder: source("scan_deps.R")
# Finds every package referenced in app.R, functions/ and parts/ and tells you
# which are (a) not installed locally, (b) not listed in deps.R.
files <- c("app.R",
           list.files("functions", pattern = "\\.R$", full.names = TRUE, recursive = TRUE),
           list.files("parts",     pattern = "\\.R$", full.names = TRUE, recursive = TRUE))
txt <- unlist(lapply(files, readLines, warn = FALSE))
txt <- sub("(^|[^\"'])#.*$", "\\1", txt)                 # drop comments (approximate)

grab <- function(rx, group = 2) {
  m <- regmatches(txt, gregexpr(rx, txt, perl = TRUE))
  unique(unlist(lapply(unlist(m), function(s) sub(rx, paste0("\\", group), s, perl = TRUE))))
}
found <- unique(c(
  grab("(library|require|requireNamespace|loadNamespace)\\(\\s*[\"']?([A-Za-z][A-Za-z0-9.]*)[\"']?"),
  grab("([A-Za-z][A-Za-z0-9.]*):::?[A-Za-z_.]", 1)
))
base_pkgs <- c("base","stats","utils","grDevices","graphics","grid","methods","tools","parallel",
               "datasets","splines","stats4","compiler","tcltk")
found <- setdiff(sort(found), base_pkgs)

deps_txt <- if (file.exists("deps.R")) paste(readLines("deps.R", warn = FALSE), collapse = "\n") else ""
in_deps  <- vapply(found, function(p) grepl(paste0("library\\(", gsub("\\.", "\\\\.", p), "\\)"), deps_txt), logical(1))
installed <- vapply(found, requireNamespace, logical(1), quietly = TRUE)

cat("\nPackages referenced in your code:", length(found), "\n")
cat("\n[1] NOT INSTALLED locally (install them, then re-run):\n"); print(found[!installed])
cat("\n[2] NOT LISTED in deps.R (add 'library(pkg)' lines inside the if (FALSE) block):\n"); print(found[!in_deps])

dyn <- grep("(paste0?|sprintf)\\([^)]*org\\.|loadNamespace\\(\\s*[a-z_]+\\s*\\)|library\\(\\s*[a-z_]+\\s*,\\s*character\\.only",
            txt, value = TRUE, perl = TRUE)
if (length(dyn)) { cat("\n[3] Dynamic package loading detected -- make sure every possible package is in deps.R:\n"); print(unique(trimws(dyn))) }
cat("\nHard-coded secrets check:\n")
sec <- grep("ghp_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,}", txt, value = TRUE)
if (length(sec)) { cat("!!! POSSIBLE TOKEN IN CODE - remove before pushing:\n"); print(sec) } else cat("none found\n")
