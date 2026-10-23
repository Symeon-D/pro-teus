#' Impute missing values using NAguideR-style imputation methods
#'
#' Impute missing values (NA/NaN/Inf) in an omics-style numeric matrix using a
#' broad set of methods aligned with those compiled in NAguideR (e.g., MinDet,
#' MinProb, QRILC, Perseus/PI, KNN, Seq-KNN, trKNN wrapper, SVD/BPCA/LLS,
#' EM/MLE, IRMI, GRR, GMS, MICE, and Random Forest).
#'
#' The function assumes the common omics orientation **features x samples**
#' (rows = proteins/genes, columns = samples). For methods that expect
#' **samples x features**, the matrix is transposed internally and transposed
#' back before returning.
#'
#' Non-finite values (`Inf`, `-Inf`, `NaN`) are coerced to `NA` prior to imputation.
#'
#' @param input_matrix Numeric matrix-like object (matrix/data.frame) with missing values.
#'   Rows are features and columns are samples.
#' @param method Imputation method to use. One of:
#'   \describe{
#'     \item{`perseus`}{Perseus/PI left-shifted normal imputation per column.}
#'     \item{`impseq`}{Sequential imputation via `rrcovNA::impSeq` (operates on samples x features internally).}
#'     \item{`impseqrob`}{Robust sequential imputation via `rrcovNA::impSeqRob`.}
#'     \item{`mindet`}{Minimal value imputation via `imputeLCMD::impute.MinDet`.}
#'     \item{`minprob`}{Minimal probabilistic imputation via `imputeLCMD::impute.MinProb`.}
#'     \item{`qr`}{QRILC imputation via `imputeLCMD::impute.QRILC`.}
#'     \item{`zero`}{Replace missing values by 0.}
#'     \item{`minimum`}{Replace missing values by the global minimum of the matrix.}
#'     \item{`colmedian`}{Replace missing values by column medians.}
#'     \item{`rowmedian`}{Replace missing values by row medians.}
#'     \item{`svd`}{SVD imputation via `pcaMethods` (`svdImpute`).}
#'     \item{`bpca`}{Bayesian PCA imputation via `pcaMethods` (`bpca`).}
#'     \item{`mle`}{EM/MLE imputation using `norm::prelim.norm`, `norm::em.norm`, `norm::imp.norm`.}
#'     \item{`knn`}{KNN imputation via `impute::impute.knn`.}
#'     \item{`seqknn`}{Sequential KNN via `SeqKnn::SeqKNN` (preferred) or `multiUS::seqKNNimp` fallback.}
#'     \item{`trknn`}{trKNN via `sim_trKNN_wrapper()` from NAguideR helper code (must be provided).}
#'     \item{`lls`}{Local Least Squares imputation via `pcaMethods::llsImpute`.}
#'     \item{`irmi`}{Iterative Robust Model-based Imputation via `VIM::irmi` (samples x features internally).}
#'     \item{`grr`}{Ridge regression-style imputation via `DreamAI::impute.RegImpute`.}
#'     \item{`gms`}{Group-wise Matrix Completion (Lasso/TS) via `GMSimpute::GMS.Lasso`.}
#'     \item{`mice_norm`}{Multiple imputation using `mice` with method `norm` (samples x features internally).}
#'     \item{`mice_cart`}{Multiple imputation using `mice` with method `cart` (samples x features internally).}
#'     \item{`rf`}{Random Forest imputation via `missForest::missForest` (samples x features internally).}
#'   }
#' @param sigma_cutoff Numeric. Perseus/PI: multiplier used to define the original width as `sigma_cutoff * sigma`.
#' @param new_width_cutoff Numeric. Perseus/PI: fraction of the original width used for the imputation distribution.
#' @param downshift Numeric. Perseus/PI: number of `sigma` units to left-shift the imputation mean from `mu`.
#' @param seed Integer. Random seed used by stochastic methods (`perseus`, `mice_*`, `rf`, etc.).
#' @param q Numeric in (0,1). Quantile parameter used by `mindet`/`minprob` (and conceptually aligned settings).
#' @param tune.sigma Numeric. Spread parameter used by `imputeLCMD` for `minprob` and `qr` (QRILC).
#' @param k Integer. Neighborhood size used by `knn`, `seqknn`, and `lls`.
#' @param nPcs Integer. Number of principal components used by `svd`/`bpca`.
#' @param center Logical. Centering option where supported (mainly PCA/LLS paths).
#' @param scale Logical. Scaling option where supported (mainly PCA/LLS paths).
#' @param trknn_source Character path to NAguideR helper script (e.g., `Imput_funcs.R`)
#'   that defines `sim_trKNN_wrapper()`. If `NULL`, the function must already exist in the session.
#' @param verbose Logical. Verbosity for methods that support it (e.g., `mice`, `missForest`).
#'
#' @return A numeric matrix of the same dimensions as `input_matrix` with missing values imputed.
#'   Row and column names are preserved when present.
#'
#' @examples
#' \dontrun{
#' # Basic:
#' x_imp <- my_impute(x, method = "mindet", q = 0.01)
#'
#' # Perseus/PI:
#' x_pi  <- my_impute(x, method = "perseus", downshift = 1.8, new_width_cutoff = 0.3)
#'
#' # KNN:
#' x_knn <- my_impute(x, method = "knn", k = 10)
#'
#' # BPCA:
#' x_bp  <- my_impute(x, method = "bpca", nPcs = 5)
#' }
#'
#' @export
my_impute <- function(
    input_matrix,
    method = c(
      "perseus", "impseq", "impseqrob",
      "mindet", "minprob", "qr",
      "zero", "minimum", "colmedian", "rowmedian",
      "svd", "bpca", "mle",
      "knn", "seqknn", "trknn", "lls",
      "irmi", "grr", "gms",
      "mice_norm", "mice_cart", "rf"
    ),
    # Perseus-style params (your existing ones)
    sigma_cutoff = 6,
    new_width_cutoff = 0.3,
    downshift = 1.8,
    seed = 123456,

    # Common params used by multiple NAguideR methods
    q = 0.01,              # for MinDet/MinProb/QRILC (NAguideR uses q-quantile idea) :contentReference[oaicite:2]{index=2}
    tune.sigma = 1,        # MinProb + QRILC in imputeLCMD commonly expose this
    k = 10,                # KNN / SeqKNN / LLS neighborhood/cluster size
    nPcs = 3,              # SVD/BPCA components
    center = TRUE,         # used where applicable
    scale = FALSE,         # used where applicable (pcaMethods handles centering/scaling internally)

    # trKNN support: NAguideR uses a helper called sim_trKNN_wrapper from "Imput_funcs.R" :contentReference[oaicite:3]{index=3}
    trknn_source = NULL,   # optional path to Imput_funcs.R if you have it

    verbose = TRUE
) {
  method <- match.arg(method)
  x <- as.matrix(input_matrix)

  # Preserve dimnames
  rn <- rownames(x); cn <- colnames(x)

  # Standardize non-finite to NA
  x[!is.finite(x)] <- NA_real_

  # Helpers ---------------------------------------------------------
  need_pkg <- function(pkg) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Method '%s' requires package '%s'. Install it first.",
        method, pkg
      ), call. = FALSE)
    }
  }

  # Many algorithms define rows=observations, cols=variables
  # For omics matrices (features x samples), "observations" are samples:
  to_obs_by_var <- function(m) t(m)    # samples x features
  back_to_feat_by_samp <- function(m) t(m)

  # Methods ---------------------------------------------------------
  if (method == "zero") {
    x[is.na(x)] <- 0
    rownames(x) <- rn; colnames(x) <- cn
    return(x)
  }

  if (method == "minimum") {
    mmin <- suppressWarnings(min(x, na.rm = TRUE))
    x[is.na(x)] <- mmin
    rownames(x) <- rn; colnames(x) <- cn
    return(x)
  }

  if (method == "colmedian") {
    med <- apply(x, 2, function(v) median(v, na.rm = TRUE))
    for (j in seq_len(ncol(x))) {
      idx <- which(is.na(x[, j]))
      if (length(idx)) x[idx, j] <- med[j]
    }
    rownames(x) <- rn; colnames(x) <- cn
    return(x)
  }

  if (method == "rowmedian") {
    med <- apply(x, 1, function(v) median(v, na.rm = TRUE))
    for (i in seq_len(nrow(x))) {
      idx <- which(is.na(x[i, ]))
      if (length(idx)) x[i, idx] <- med[i]
    }
    rownames(x) <- rn; colnames(x) <- cn
    return(x)
  }

  if (method == "mindet") {
    need_pkg("imputeLCMD")
    out <- imputeLCMD::impute.MinDet(x, q = q)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "minprob") {
    need_pkg("imputeLCMD")
    set.seed(seed)
    out <- imputeLCMD::impute.MinProb(x, q = q, tune.sigma = tune.sigma)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "qr") {
    need_pkg("imputeLCMD")
    set.seed(seed)
    out <- imputeLCMD::impute.QRILC(x, tune.sigma = tune.sigma)[[1]]
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "perseus") {
    set.seed(seed)
    need_pkg("MASS")

    out <- x
    # Fallback pool if a column is too sparse/constant
    pool <- as.numeric(x)
    pool <- pool[is.finite(pool)]
    pool_mu <- mean(pool)
    pool_sd <- stats::sd(pool)
    if (!is.finite(pool_sd) || pool_sd == 0) pool_sd <- 1e-6

    for (j in seq_len(ncol(x))) {
      col_vals <- x[, j]
      col_vals <- col_vals[is.finite(col_vals)]
      if (length(col_vals) < 2) {
        mu <- pool_mu; sigma <- pool_sd
      } else {
        fit <- tryCatch(MASS::fitdistr(col_vals, "normal"), error = function(e) NULL)
        if (is.null(fit)) {
          mu <- mean(col_vals); sigma <- stats::sd(col_vals)
        } else {
          mu <- unname(fit$estimate["mean"])
          sigma <- unname(fit$estimate["sd"])
        }
        if (!is.finite(sigma) || sigma == 0) sigma <- pool_sd
        if (!is.finite(mu)) mu <- pool_mu
      }

      width <- sigma_cutoff * sigma
      new_width <- width * new_width_cutoff
      new_sigma <- new_width / sigma_cutoff
      new_mean <- mu - downshift * sigma

      na_idx <- which(is.na(x[, j]))
      if (length(na_idx) > 0) {
        out[na_idx, j] <- stats::rnorm(length(na_idx), mean = new_mean, sd = new_sigma)
      }
    }

    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "impseq") {
    need_pkg("rrcovNA")
    X <- to_obs_by_var(x)
    out <- rrcovNA::impSeq(X)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "impseqrob") {
    need_pkg("rrcovNA")
    X <- to_obs_by_var(x)
    out <- rrcovNA::impSeqRob(X)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "knn") {
    # NAguideR uses impute::impute.knn :contentReference[oaicite:4]{index=4}
    need_pkg("impute")
    kk <- min(k, nrow(x) - 1L)
    if (kk < 1L) stop("KNN needs at least 2 rows (features).", call. = FALSE)
    tmp <- impute::impute.knn(x, k = kk)
    out <- tmp$data
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "seqknn") {
    # NAguideR originally used SeqKnn, but multiUS::seqKNNimp is the reliable, 
    # supported alternative that avoids modern package installation and deployment failures.
    if (requireNamespace("multiUS", quietly = TRUE)) {
      out <- as.matrix(multiUS::seqKNNimp(as.data.frame(x), k = k))
    } else {
      stop("seqknn imputation requires the 'multiUS' package.", call. = FALSE)
    }
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "trknn") {
    # NAguideR calls sim_trKNN_wrapper from a helper script Imput_funcs.R :contentReference[oaicite:5]{index=5}
    if (!exists("sim_trKNN_wrapper", mode = "function")) {
      if (!is.null(trknn_source)) {
        source(trknn_source)
      }
    }
    if (!exists("sim_trKNN_wrapper", mode = "function")) {
      stop(
        "trknn requires function 'sim_trKNN_wrapper' (NAguideR helper). ",
        "Provide 'trknn_source' pointing to Imput_funcs.R, or define sim_trKNN_wrapper() in your session.",
        call. = FALSE
      )
    }
    out <- sim_trKNN_wrapper(x)
    out <- as.matrix(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "lls") {
    need_pkg("pcaMethods")
    res <- pcaMethods::llsImpute(x, k = k, center = center, completeObs = FALSE)
    out <- pcaMethods::completeObs(res)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "svd") {
    need_pkg("pcaMethods")
    X <- to_obs_by_var(x)
    pc <- pcaMethods::pca(X, method = "svdImpute", nPcs = nPcs, center = center, scale = if (isTRUE(scale)) "pareto" else "none")
    out <- pcaMethods::completeObs(pc)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "bpca") {
    need_pkg("pcaMethods")
    X <- to_obs_by_var(x)
    pc <- pcaMethods::pca(X, method = "bpca", nPcs = nPcs, center = center, scale = if (isTRUE(scale)) "pareto" else "none")
    out <- pcaMethods::completeObs(pc)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "mle") {
    # norm::prelim.norm + em.norm + imp.norm :contentReference[oaicite:6]{index=6}
    need_pkg("norm")
    X <- to_obs_by_var(x)
    s <- norm::prelim.norm(X)
    theta <- norm::em.norm(s)
    out <- norm::imp.norm(s, theta, X)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "irmi") {
    # VIM::irmi (Iterative Robust Model-based Imputation) :contentReference[oaicite:7]{index=7}
    need_pkg("VIM")
    X <- as.data.frame(to_obs_by_var(x))
    out <- VIM::irmi(X)
    out <- as.matrix(out)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "grr") {
    # DreamAI::impute.RegImpute (glmnet ridge regression) :contentReference[oaicite:8]{index=8}
    need_pkg("DreamAI")
    out <- DreamAI::impute.RegImpute(x)
    out <- as.matrix(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "gms") {
    # GMSimpute::GMS.Lasso :contentReference[oaicite:9]{index=9}
    need_pkg("GMSimpute")
    out <- tryCatch(
      GMSimpute::GMS.Lasso(x, log.scale = TRUE, TS.Lasso = TRUE),
      error = function(e) NULL
    )
    if (is.null(out)) {
      # fallback orientation if user’s matrix is sample x feature
      out2 <- GMSimpute::GMS.Lasso(to_obs_by_var(x), log.scale = TRUE, TS.Lasso = TRUE)
      out <- back_to_feat_by_samp(out2)
    }
    out <- as.matrix(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

  if (method == "mice_norm" || method == "mice_cart") {
    need_pkg("mice")
    X <- as.data.frame(to_obs_by_var(x))
    mice_method <- if (method == "mice_norm") "norm" else "cart"
    set.seed(seed)
    imp <- mice::mice(X, m = 1, maxit = 5, method = mice_method, printFlag = verbose)
    out <- mice::complete(imp, 1)
    out <- as.matrix(out)
    out <- back_to_feat_by_samp(out)
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
  }

if (method == "rf") {
    need_pkg("missForest")
    set.seed(seed)
    tmp <- missForest::missForest(x, verbose = verbose)
    out <- tmp$ximp
    rownames(out) <- rn; colnames(out) <- cn
    return(out)
}

  stop("Unhandled method: ", method, call. = FALSE)
}
