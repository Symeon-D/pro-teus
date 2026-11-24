# Regions of a Venn diagram: every protein belongs to exactly ONE region = the exact combination of sets it is in.
#   sets: named list of character vectors (names = circle titles; must be unique)
# returns list(summary = data.frame(Region, Sets, N_sets, Count), members = named list region -> proteins,
#              long = data.frame(Protein, Region, <one yes/no column per set>))
venn_regions <- function(sets) {
  nm <- names(sets); k <- length(sets)
  all_p <- unique(unlist(sets, use.names = FALSE))
  mem <- vapply(sets, function(s) all_p %in% s, logical(length(all_p)))
  if (!is.matrix(mem)) mem <- matrix(mem, nrow = length(all_p), dimnames = list(NULL, nm))
  rownames(mem) <- all_p
  combos <- unlist(lapply(seq_len(k), function(m) combn(seq_len(k), m, simplify = FALSE)), recursive = FALSE)
  label <- function(ix) if (length(ix) == 1) paste("Unique to", nm[ix]) else paste("Common to", paste(nm[ix], collapse = " + "))
  members <- list(); rows <- list()
  for (ix in combos) {
    inn <- rowSums(mem[, ix, drop = FALSE]) == length(ix)
    out <- if (length(ix) < k) rowSums(mem[, -ix, drop = FALSE]) == 0 else rep(TRUE, nrow(mem))
    p <- sort(all_p[inn & out])
    lab <- label(ix)
    members[[lab]] <- p
    rows[[length(rows) + 1]] <- data.frame(Region = lab, Sets = paste(nm[ix], collapse = " + "),
                                           N_sets = length(ix), Count = length(p), stringsAsFactors = FALSE)
  }
  summ <- do.call(rbind, rows)
  summ <- summ[order(summ$N_sets, -summ$Count), , drop = FALSE]; rownames(summ) <- NULL
  members <- members[summ$Region]
  long <- do.call(rbind, lapply(names(members), function(r) {
    p <- members[[r]]; if (!length(p)) return(NULL)
    d <- data.frame(Protein = p, Region = r, stringsAsFactors = FALSE)
    for (s in nm) d[[s]] <- ifelse(mem[p, s], "yes", "no")
    d
  }))
  if (is.null(long)) long <- data.frame(Protein = character(), Region = character())
  list(summary = summ, members = members, long = long)
}

# detected-protein sets as used by the QC Venn (a protein is detected in a group when it is present in >= thr of its samples)
detected_sets <- function(m, md, groups, thr, titles) {
  sets <- lapply(groups, function(g) {
    s <- md$Sample[md$Group == g]
    sub <- m[, s, drop = FALSE]
    rownames(sub)[rowSums(!is.na(sub)) / length(s) >= thr]
  })
  titles[is.na(titles) | !nzchar(titles)] <- groups[is.na(titles) | !nzchar(titles)]
  names(sets) <- make.unique(titles, sep = " #")
  sets
}

# Excel workbook with 3 sheets: counts per region / long table / one column per region
venn_workbook <- function(vr, file, thr = NULL) {
  wb <- openxlsx::createWorkbook()
  hs <- openxlsx::createStyle(textDecoration = "bold", fgFill = "#E8EEF7", border = "bottom")
  add <- function(name, d) {
    openxlsx::addWorksheet(wb, name); openxlsx::writeData(wb, name, d, headerStyle = hs)
    openxlsx::setColWidths(wb, name, cols = seq_len(ncol(d)), widths = "auto"); openxlsx::freezePane(wb, name, firstRow = TRUE)
  }
  add("Summary", vr$summary)
  add("Proteins by region", vr$long)
  n <- max(lengths(vr$members), 1)
  wide <- as.data.frame(lapply(vr$members, function(p) c(p, rep(NA_character_, n - length(p)))),
                        check.names = FALSE, stringsAsFactors = FALSE)
  add("One column per region", wide)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}
