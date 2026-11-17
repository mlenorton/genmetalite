#' Empirical covariance of a design matrix that may contain factors
#'
#' Expands any character/factor columns to treatment-coded dummy variables and
#' returns both the (optionally weighted) empirical covariance matrix and the
#' expanded numeric design matrix.
#'
#' @param X A data frame or matrix of covariates. Must not contain the outcome.
#' @param wt Optional numeric vector of observation weights, passed to
#'   [stats::cov.wt()]. `NULL` (default) means equal weights.
#' @param ... Further arguments passed to [stats::cov.wt()], e.g.
#'   `method = "ML"`.
#'
#' @details
#' Unused factor levels are dropped before expansion: an empty level produces an
#' all-zero dummy column, which makes the covariance matrix singular and causes
#' a confusing failure much later in the pipeline.
#'
#' Expansion uses [stats::model.matrix()], so a factor `edu` with levels
#' `1, 2, 3` becomes columns `edu2` and `edu3`. This is the same naming used for
#' the reference population elsewhere in the package, so covariance blocks match
#' by name across datasets.
#'
#' @return A list with components
#'   \item{cov_X_s}{the covariance matrix of the expanded design.}
#'   \item{X_dummy_mat}{the expanded numeric design matrix (no intercept).}
#'   \item{center}{the (weighted) column means.}
#'
#' @examples
#' d <- data.frame(age = rnorm(50), edu = factor(sample(1:3, 50, TRUE)))
#' out <- cov_factors(d)
#' colnames(out$X_dummy_mat)
#'
#' @export
cov_factors <- function(X, wt = NULL, ...) {
  X <- as.data.frame(X)
  X_dummy_mat <- dummy_matrix(X, add_intercept = FALSE)

  if (anyNA(X_dummy_mat)) {
    stop(
      "Covariates contain missing values. Covariance cannot be estimated ",
      "reliably from incomplete rows; drop or impute them before calling ",
      "cov_factors().",
      call. = FALSE
    )
  }

  emp <- if (is.null(wt)) {
    stats::cov.wt(X_dummy_mat, ...)
  } else {
    stats::cov.wt(X_dummy_mat, wt = wt, ...)
  }

  list(
    cov_X_s = emp$cov,
    X_dummy_mat = X_dummy_mat,
    center = emp$center
  )
}
