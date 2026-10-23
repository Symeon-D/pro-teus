# ==========================================================
# Shared helpers (no dependency on the retired Proteomics.Functions package)
# ==========================================================

`%||%` <- function(a, b) if (is.null(a)) b else a

read_table_any <- function(path, name) {
  ext <- tolower(tools::file_ext(name))
  if (ext %in% c("xlsx", "xls")) {
    as.data.frame(readxl::read_excel(path))
  } else {
    sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
    read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "",
               check.names = FALSE, stringsAsFactors = FALSE, fill = TRUE)
  }
}

auto_colors <- function(levels_vec, pal = "Set2") {
  n <- length(levels_vec)
  base <- RColorBrewer::brewer.pal(max(3, min(n, 8)), pal)
  cols <- if (n <= 8) base[seq_len(n)] else grDevices::colorRampPalette(base)(n)
  setNames(cols, levels_vec)
}

# Several of your functions write PNGs into the working directory with fixed names.
# This runs such a function inside a throw-away directory and returns the result
# plus the PNG(s) it produced, so the app can show / offer them for download.
run_in_tmpdir <- function(fun) {
  td <- tempfile("run_"); dir.create(td)
  owd <- setwd(td); on.exit(setwd(owd), add = TRUE)
  val <- fun()
  list(value = val, files = list.files(td, pattern = "\\.png$", full.names = TRUE))
}

# ---------- demo dataset ------------------------------------
# 2 genotypes x 2 conditions x 3 animals x 2 technical replicates = 24 samples
make_demo_data <- function(seed = 1) {
  set.seed(seed)
  # protein names come from a real data set (data/demo_proteins.csv); every intensity is simulated here
  f <- "data/demo_proteins.csv"
  prot <- if (file.exists(f)) read.csv(f, stringsAsFactors = FALSE) else
    data.frame(Gene = paste0("Prot", 1:1000), Base = rnorm(1000, 15, 2))
  ids <- prot$Gene; mu <- prot$Base; n <- length(ids)

  # 2 x 2 factorial design: Genotype (WT, KO) x Condition (Con, Stress), 3 animals each, 2 technical replicates
  meta <- expand.grid(Replicate = 1:2, Animal_n = 1:3, Condition = c("Con", "Stress"),
                      Genotype = c("WT", "KO"), stringsAsFactors = FALSE)
  meta$Animal <- paste(meta$Genotype, meta$Condition, paste0("M", meta$Animal_n), sep = "_")
  meta$Sample <- paste(meta$Animal, meta$Replicate, sep = "_")
  meta$Batch  <- c("A", "B", "A")[meta$Animal_n]
  meta <- meta[order(factor(meta$Genotype, levels = c("WT", "KO")), factor(meta$Condition, levels = c("Con", "Stress")),
                     meta$Animal_n, meta$Replicate), ]
  meta <- meta[, c("Sample", "Genotype", "Condition", "Animal", "Batch")]
  rownames(meta) <- NULL

  # planted biology (gene families -> GO terms) so ORA / GSEA have something to find
  has <- function(p) grepl(p, ids)
  f_ribo <- has("^(Rpl|Rps)[0-9]"); f_syn <- has("^(Syn[0-9a-z]|Snap|Syt|Stx|Vamp|Sv2|Slc17a)")
  f_resp <- has("^(Ndufa|Ndufb|Ndufs|Ndufv|Cox[0-9]|Atp5|Uqcr|Sdh)")
  f_int  <- has("^(Hsp|Dnaj|Mrpl|Mrps|H2a|H2b|Hist|H3|H4)")
  rnd <- function(p, lo, hi) ifelse(runif(n) < p, sample(c(-1, 1), n, TRUE) * runif(n, lo, hi), 0)
  eff_stress <- ifelse(f_ribo, 1.0, ifelse(f_syn, -0.9, 0)) + rnd(0.05, 0.6, 1.3)     # stress, both genotypes
  eff_ko     <- ifelse(f_resp, -1.0, 0) + rnd(0.04, 0.6, 1.2)                          # KO baseline
  eff_int    <- ifelse(f_int, 1.3, 0) + rnd(0.02, 0.6, 1.2)                            # KO x Stress interaction

  m <- nrow(meta)
  animal_eff <- matrix(rnorm(n * length(unique(meta$Animal)), 0, 0.15), n, dimnames = list(NULL, unique(meta$Animal)))
  sample_off <- rnorm(m, 0, 0.3)
  x <- matrix(NA_real_, n, m, dimnames = list(ids, meta$Sample))
  for (j in seq_len(m)) {
    s <- as.numeric(meta$Condition[j] == "Stress"); k <- as.numeric(meta$Genotype[j] == "KO")
    x[, j] <- mu + eff_stress * s + eff_ko * k + eff_int * s * k +
      animal_eff[, meta$Animal[j]] + sample_off[j] + rnorm(n, 0, 0.12)
  }
  q <- quantile(mu, 0.15)
  p_miss <- plogis(-(x - q) * 1.0 - 1.5) * 0.9 + 0.03          # MNAR (low intensity) + MAR
  x[matrix(runif(n * m) < p_miss, n, m)] <- NA
  list(expr = x, meta = meta, id_name = "Genes")
}

# ---------- normalization benchmark utilities ---------------
read_norm_file <- function(f, ids, samples) {
  df <- tryCatch(read.delim(f, check.names = FALSE, stringsAsFactors = FALSE), error = function(e) NULL)
  if (is.null(df) || !all(samples %in% colnames(df))) return(NULL)
  m <- as.matrix(df[, samples, drop = FALSE]); storage.mode(m) <- "numeric"
  rownames(m) <- as.character(df[[1]])
  if (all(ids %in% rownames(m))) m <- m[ids, , drop = FALSE]
  m
}

pooled_cv <- function(m, groups) {
  lin <- 2^m
  unlist(lapply(unique(groups), function(g) {
    x <- lin[, groups == g, drop = FALSE]
    ok <- rowSums(!is.na(x)) >= 2
    if (!any(ok)) return(numeric(0))
    x <- x[ok, , drop = FALSE]
    apply(x, 1, sd, na.rm = TRUE) / rowMeans(x, na.rm = TRUE) * 100
  }))
}

stage_key <- function(m) paste(dim(m), sum(is.na(m)), round(sum(m, na.rm = TRUE), 3), collapse = "_")

# ---------- expression file readers: format and gene-name column are detected automatically ----------
first_token <- function(x) sub(";.*$", "", as.character(x))

# data-URI download links (no server round trip): used for the small per-row download icons
uri_escape <- function(txt) if (requireNamespace("curl", quietly = TRUE)) curl::curl_escape(txt) else utils::URLencode(txt, reserved = TRUE)
csv_data_uri <- function(df) {
  txt <- paste(utils::capture.output(utils::write.csv(df, row.names = FALSE)), collapse = "\n")
  paste0("data:text/csv;charset=utf-8,", uri_escape(enc2utf8(txt)))
}

# a column is "text-like" when most of its non-missing values are not numbers
is_text_col <- function(v) {
  v <- as.character(v); v <- v[!is.na(v) & nzchar(v) & !startsWith(v, "#!")]
  if (!length(v)) return(FALSE)
  mean(is.na(suppressWarnings(as.numeric(v)))) > 0.5
}

# Finds the column holding gene / protein names; returns its index or NA when it cannot be told
find_id_column <- function(df) {
  nm <- names(df); txt <- vapply(df, is_text_col, logical(1))
  pick <- function(idx) { idx <- idx[txt[idx]]; if (length(idx)) idx[1] else NA_integer_ }
  r <- pick(grep("^(pg\\.|PG\\.)?(genes?|gene[ ._]?names?|gene[ ._]?symbols?|symbols?)$", nm, ignore.case = TRUE))
  if (is.na(r)) r <- pick(grep("gene|symbol", nm, ignore.case = TRUE))
  if (is.na(r)) r <- pick(grep("^(protein[ ._]?(group|ids?|names?|accession)?|majority protein ids|accession|pg\\.protein)", nm, ignore.case = TRUE))
  if (is.na(r) && length(txt) && txt[1]) r <- 1L          # e.g. row names exported with an empty header
  r
}

# Looks at the file signature: Perseus ("#!{Type}" rows), MaxQuant (Majority protein IDs / LFQ intensity ...),
# DIA-NN (Protein.Group + Genes) or a plain matrix.
detect_expression_format <- function(path, name) {
  if (grepl("\\.(xlsx|xls)$", name, ignore.case = TRUE)) return(list(format = "table", quant = NULL))
  lines <- readLines(path, n = 8, warn = FALSE)
  if (!length(lines)) stop("The file is empty.")
  if (any(startsWith(lines, "#!") | grepl("\t#!\\{", lines))) return(list(format = "perseus", quant = NULL))
  sep <- if (grepl("\t", lines[1])) "\t" else ","
  cols <- gsub("^\"|\"$", "", strsplit(lines[1], sep, fixed = TRUE)[[1]])
  if ("Majority protein IDs" %in% cols || any(startsWith(cols, "LFQ intensity ")) ) {
    q <- if (any(startsWith(cols, "LFQ intensity "))) "LFQ intensity " else if (any(startsWith(cols, "Intensity "))) "Intensity " else "iBAQ "
    return(list(format = "maxquant", quant = q))
  }
  if ("Protein.Group" %in% cols && any(c("Genes", "Protein.Ids", "Protein.Names") %in% cols))
    return(list(format = "diann", quant = NULL))
  list(format = "table", quant = NULL)
}

format_label <- function(f) c(table = "generic matrix", perseus = "Perseus export", maxquant = "MaxQuant proteinGroups",
                              diann = "DIA-NN pg_matrix")[[f]]

# first rows only: format, column names, text-like candidate columns and the auto-detected gene column
inspect_expression_file <- function(path, name) {
  det <- detect_expression_format(path, name)
  if (grepl("\\.(xlsx|xls)$", name, ignore.case = TRUE)) {
    df <- read_table_any(path, name)
  } else {
    l1 <- readLines(path, n = 1, warn = FALSE)
    sep <- if (grepl("\t", l1)) "\t" else ","
    df <- utils::read.table(path, header = TRUE, sep = sep, quote = "\"", comment.char = "", check.names = FALSE,
                            stringsAsFactors = FALSE, fill = TRUE)
  }
  df <- df[!startsWith(as.character(df[[1]]), "#!"), , drop = FALSE]
  cols <- names(df)
  txt <- which(vapply(df, function(v) is_text_col(utils::head(v, 5000)), logical(1)))
  cand <- if (length(txt)) txt else seq_along(cols)     # every non-numeric column (all of them, not only the visible ones)
  list(format = det$format, cols = cols, candidates = cand, auto = find_id_column(df))
}

# returns list(m = numeric matrix proteins x samples, log2_applied, format);
# id_col = index of the column holding gene names (NULL = detect)
read_expression_file <- function(path, name, format = "auto", mq_quant = NULL, id_col = NULL) {
  if (format == "auto") { det <- detect_expression_format(path, name); format <- det$format; mq_quant <- det$quant }
  need_id <- function(df) {
    i <- if (!is.null(id_col) && !is.na(id_col) && id_col >= 1 && id_col <= ncol(df)) id_col else find_id_column(df)
    if (is.na(i)) stop("The column holding the gene names could not be identified - please select it.")
    i
  }
  ids_from <- function(df, i, alt = NULL) {
    ids <- first_token(df[[i]])
    if (!is.null(alt)) { bad <- is.na(ids) | ids == ""; ids[bad] <- alt[bad] }
    make.unique(as.character(ids))
  }

  if (format == "table") {
    df <- read_table_any(path, name)
    if (ncol(df) < 3) stop("Expected a protein-name column plus at least two sample columns.")
    i <- need_id(df)
    others <- setdiff(seq_len(ncol(df)), i)
    num <- vapply(others, function(j) !is_text_col(df[[j]]), logical(1))
    sc <- others[num]
    if (length(sc) < 2) stop("Could not find at least two numeric sample columns.")
    m <- matrix(suppressWarnings(as.numeric(as.matrix(df[, sc, drop = FALSE]))), nrow = nrow(df),
                dimnames = list(ids_from(df, i), names(df)[sc]))
    return(list(m = m, log2_applied = FALSE, format = "table", id_name = names(df)[i]))
  }

  df <- read.delim(path, check.names = FALSE, stringsAsFactors = FALSE, quote = "")
  if (format == "perseus") {
    first <- as.character(df[[1]])
    type_row <- which(startsWith(first, "#!{Type}"))
    tp <- NULL
    if (length(type_row)) {            # first cell is "#!{Type}" + the type of the first column
      tp <- as.character(unlist(df[type_row[1], ]))
      tp[1] <- sub("^#!\\{Type\\}", "", tp[1])
    }
    df <- df[!startsWith(first, "#!"), , drop = FALSE]
    i <- need_id(df)
    num_idx <- if (!is.null(tp)) which(tp == "E") else setdiff(seq_len(ncol(df)), i)
    num_idx <- setdiff(num_idx, i)
    if (length(num_idx) < 2) stop("Perseus file: fewer than two expression columns (type 'E').")
    m <- matrix(suppressWarnings(as.numeric(as.matrix(df[, num_idx, drop = FALSE]))), nrow = nrow(df),
                dimnames = list(ids_from(df, i), names(df)[num_idx]))
    m[!is.finite(m)] <- NaN
    return(list(m = m, log2_applied = TRUE, format = "perseus", id_name = names(df)[i]))   # Perseus data are already log2
  }
  if (format == "maxquant") {
    flag <- function(col) if (col %in% names(df)) df[[col]] == "+" else rep(FALSE, nrow(df))
    df <- df[!(flag("Reverse") | flag("Potential contaminant") | flag("Only identified by site")), , drop = FALSE]
    qcols <- names(df)[startsWith(names(df), mq_quant)]
    if (!length(qcols)) stop("No columns starting with '", mq_quant, "' found in this MaxQuant file.")
    i <- if (!is.null(id_col) && !is.na(id_col) && id_col <= ncol(df)) id_col else match("Gene names", names(df))
    if (is.na(i)) i <- need_id(df)
    alt <- if ("Majority protein IDs" %in% names(df)) first_token(df[["Majority protein IDs"]]) else seq_len(nrow(df))
    m <- as.matrix(df[, qcols, drop = FALSE]); storage.mode(m) <- "numeric"
    colnames(m) <- sub(mq_quant, "", qcols, fixed = TRUE)
    rownames(m) <- ids_from(df, i, alt)
  } else if (format == "diann") {
    ann <- intersect(c("Protein.Group", "Protein.Ids", "Protein.Names", "Genes",
                       "First.Protein.Description", "Protein.Description"), names(df))
    sc <- setdiff(names(df), ann)
    sc <- sc[vapply(df[sc], is.numeric, logical(1))]
    if (!length(sc)) stop("No numeric sample columns found (expected a DIA-NN pg_matrix file).")
    i <- if (!is.null(id_col) && !is.na(id_col) && id_col <= ncol(df)) id_col else match("Genes", names(df))
    if (is.na(i)) i <- need_id(df)
    alt <- if ("Protein.Group" %in% names(df)) first_token(df$Protein.Group) else seq_len(nrow(df))
    m <- as.matrix(df[, sc, drop = FALSE]); storage.mode(m) <- "numeric"
    colnames(m) <- sub("\\.(raw|mzML|d|dia|wiff|tdf)$", "", basename(gsub("\\\\", "/", sc)))
    rownames(m) <- ids_from(df, i, alt)
  } else stop("Unknown file format: ", format)
  m[!is.finite(m) | m <= 0] <- NA
  list(m = safe_log2(m), log2_applied = TRUE, format = format, id_name = names(df)[i])
}

# log2 with Inf / -Inf silently turned into NaN (so they count as missing values)
safe_log2 <- function(x) {
  y <- suppressWarnings(log2(x))
  y[!is.finite(y)] <- NaN
  y
}

# The two graphs of the first page of evaluate_imputation(): missing-value heatmap + % missing per sample
plot_missing_pattern <- function(mat, groups = NULL, group_colors = NULL) {
  na_mat    <- is.na(mat)
  row_order <- order(rowSums(na_mat), decreasing = TRUE)
  melt_na   <- reshape2::melt(na_mat[row_order, , drop = FALSE])
  colnames(melt_na) <- c("Protein", "Sample", "Missing")
  use_grp <- !is.null(groups) && !is.null(group_colors)
  if (use_grp) {
    gmap <- setNames(as.character(groups), colnames(mat))
    melt_na$State <- ifelse(melt_na$Missing, gmap[as.character(melt_na$Sample)], "Observed")
    cols <- c(group_colors, Observed = "grey92")
  }
  p_heat <- if (use_grp) {
    ggplot(melt_na, aes(x = Sample, y = Protein, fill = State)) + geom_raster() +
      scale_fill_manual(values = cols, breaks = names(group_colors)) +
      labs(title = "Missing Value Pattern (missing = group colour)", x = "", y = "Proteins", fill = "")
  } else {
    ggplot(melt_na, aes(x = Sample, y = Protein, fill = Missing)) + geom_raster() +
      scale_fill_manual(values = c("FALSE" = "steelblue", "TRUE" = "tomato"), labels = c("Observed", "Missing")) +
      labs(title = "Missing Value Pattern", x = "", y = "Proteins", fill = "")
  }
  p_heat <- p_heat + theme_minimal(base_size = 11) +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "top",
          plot.margin = margin(5, 10, 5, 5))
  na_pct <- data.frame(Sample = colnames(mat), NA_pct = colMeans(is.na(mat)) * 100)
  na_pct$Sample <- factor(na_pct$Sample, levels = na_pct$Sample)
  if (use_grp) na_pct$Group <- as.character(groups)
  p_bar <- (if (use_grp) ggplot(na_pct, aes(x = Sample, y = NA_pct, fill = Group)) + geom_col() + scale_fill_manual(values = group_colors)
            else ggplot(na_pct, aes(x = Sample, y = NA_pct, fill = NA_pct)) + geom_col() + scale_fill_gradient(low = "steelblue", high = "tomato")) +
    labs(title = "% Missing per Sample", x = "", y = "% Missing") +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none",
          plot.margin = margin(5, 5, 5, 10)) + coord_flip()
  cowplot::plot_grid(p_heat, p_bar, ncol = 2, rel_widths = c(1, 1))
}

# the gene / protein name column is called "Names" everywhere the user sees it
as_names <- function(df) { if ("Protein" %in% names(df)) names(df)[names(df) == "Protein"] <- "Names"; df }
