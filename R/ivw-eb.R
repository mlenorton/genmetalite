#' Inverse-variance weighted pooling of reduced models
#'
#' Pools each coefficient across studies by ordinary inverse-variance weighted
#' meta-analysis, ignoring the fact that the reduced models omit different
#' covariates. This is the naive comparator for [linear_estimate()]: it is
#' efficient but biased whenever the omitted covariates are correlated with the
#' retained ones, because each study estimates a different marginal effect.
#'
#' @param reduced_models A list of reduced models; see [run_regression()].
#'
#' @return A list with
#'   \item{beta}{named vector of pooled coefficients.}
#'   \item{var}{named vector of pooled variances (one per coefficient, since
#'     the coefficients are pooled one at a time).}
#'   \item{n_studies}{named vector counting how many studies contributed to each
#'     coefficient.}
#'
#' @examples
#' set.seed(5)
#' mk <- function(vars) {
#'   d <- data.frame(x1 = rnorm(400), x2 = rnorm(400), x3 = rnorm(400))
#'   d$y <- 1 + d$x1 + 0.5 * d$x2 - 0.2 * d$x3 + rnorm(400)
#'   run_regression(d$y, d[vars], dgp = "linear")
#' }
#' rms <- list(mk(c("x1", "x2")), mk(c("x1", "x3")))
#' round(IVW_estimate(rms)$beta, 3)
#'
#' @seealso [linear_estimate()], [EB_estimate()]
#' @export
IVW_estimate <- function(reduced_models) {
  n_cohorts <- length(reduced_models)
  if (n_cohorts == 0L) {
    stop("`reduced_models` is empty.", call. = FALSE)
  }

  all_var_names <- unique(unlist(lapply(reduced_models, function(r) {
    names(r$theta_hat_s)
  })))

  coeff <- matrix(
    NA_real_, n_cohorts, length(all_var_names),
    dimnames = list(NULL, all_var_names)
  )
  prec <- coeff

  for (s in seq_len(n_cohorts)) {
    th <- reduced_models[[s]]$theta_hat_s
    v <- diag(reduced_models[[s]]$cov_theta_hat_s)
    nms <- names(th)
    coeff[s, nms] <- as.numeric(th)
    ## Precision, not standard error: guard against zero/negative variances
    ## produced by a degenerate sandwich estimate.
    v <- v[nms]
    prec[s, nms] <- ifelse(is.finite(v) & v > 0, 1 / v, NA_real_)
  }

  ## A coefficient contributes only where BOTH estimate and precision exist.
  usable <- is.finite(coeff) & is.finite(prec)
  coeff[!usable] <- NA_real_
  prec[!usable] <- NA_real_

  sum_prec <- colSums(prec, na.rm = TRUE)
  if (any(sum_prec == 0)) {
    bad <- all_var_names[sum_prec == 0]
    stop(sprintf(
      "No usable variance estimates for coefficient(s): %s.",
      paste(bad, collapse = ", ")
    ), call. = FALSE)
  }

  beta <- colSums(coeff * prec, na.rm = TRUE) / sum_prec
  var <- 1 / sum_prec

  list(
    beta = stats::setNames(beta, all_var_names),
    var = stats::setNames(var, all_var_names),
    n_studies = stats::setNames(colSums(usable), all_var_names)
  )
}

#' Empirical-Bayes blend of the mapped and inverse-variance weighted estimates
#'
#' Shrinks the (approximately unbiased but higher-variance) mapped estimate from
#' [linear_estimate()] towards the (lower-variance but biased) inverse-variance
#' weighted estimate from [IVW_estimate()], coefficient by coefficient.
#'
#' @param linear_est Output of [linear_estimate()].
#' @param IVW_est Output of [IVW_estimate()].
#'
#' @return A list with
#'   \item{beta}{the blended coefficients.}
#'   \item{K}{the per-coefficient weight placed on `linear_est$beta`.}
#'   \item{bias}{the estimated bias of the IVW estimate.}
#'
#' @details
#' For each coefficient, with \eqn{b = \hat\beta_{\mathrm{IVW}} -
#' \hat\beta_{\mathrm{linear}}} and \eqn{v} the variance of the mapped estimate,
#' the weight is \eqn{K = b^2 / (v + b^2)} and the blend is
#' \eqn{K\hat\beta_{\mathrm{linear}} + (1 - K)\hat\beta_{\mathrm{IVW}}}. When the
#' apparent bias is large relative to the sampling error, \eqn{K \to 1} and the
#' mapped estimate is used; when the two estimates agree, \eqn{K \to 0} and the
#' more precise IVW estimate is used.
#'
#' No variance is returned. The shrinkage weight is itself estimated from the
#' data, so the naive plug-in variance would understate the uncertainty; use a
#' bootstrap if a standard error for the blended estimate is needed.
#'
#' @seealso [linear_estimate()], [IVW_estimate()]
#' @export
EB_estimate <- function(linear_est, IVW_est) {
  nms <- names(linear_est$beta)
  missing <- setdiff(nms, names(IVW_est$beta))
  if (length(missing) > 0L) {
    stop(sprintf(
      "The IVW estimate is missing coefficient(s): %s.",
      paste(missing, collapse = ", ")
    ), call. = FALSE)
  }

  ivw_beta <- IVW_est$beta[nms]
  bias <- ivw_beta - linear_est$beta
  v <- diag(linear_est$var)[nms]

  denom <- v + bias^2
  K <- ifelse(denom > 0, bias^2 / denom, 1)

  list(
    beta = stats::setNames(K * linear_est$beta + (1 - K) * ivw_beta, nms),
    K = stats::setNames(K, nms),
    bias = stats::setNames(bias, nms)
  )
}
