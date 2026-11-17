#' Convert an existing model fit into the reduced-model format
#'
#' Use this when a study's reduced model was fitted elsewhere (or by a
#' collaborator who cannot share individual-level data) and only the fitted
#' object is available.
#'
#' @param reg_out A fitted `lm` or `glm` object.
#' @param var_type Passed to [compute_var()]. Ignored when `vcov` is supplied.
#' @param vcov Optional coefficient covariance matrix to use instead of
#'   `compute_var()`, or a function of `reg_out` returning one. Use this when
#'   the study is itself a complex sample: the default sandwich assumes
#'   independent rows, so a clustered cohort needs a cluster-robust matrix such
#'   as `sandwich::vcovCL(fit, cluster = ~psu)`. See Details.
#' @param keep_intercept Logical. If `TRUE` (default) the intercept is retained
#'   in `theta_hat_s` and `cov_theta_hat_s`, matching [run_regression()]. Set to
#'   `FALSE` to drop it.
#' @param ... Unused; present for forward compatibility.
#'
#' @return A list of class `reduced_model`; see [run_regression()].
#'
#' @section When the study is itself a complex sample:
#' [linear_estimate()] takes each study's coefficient covariance
#' \eqn{V_s} as given and pools with \eqn{V_s^{-1}}, so \eqn{V_s} must reflect
#' however that study was sampled. The default Huber-White sandwich assumes
#' independent rows. If a cohort is clustered — repeated measures, households,
#' clinics, or a survey of its own — that assumption understates \eqn{V_s}, and
#' the pooled standard errors inherit the problem. Supply a cluster-robust
#' matrix instead:
#'
#' ```r
#' rm1 <- as_reduced_model(
#'   fit,
#'   vcov = function(f) sandwich::vcovCL(f, cluster = mydata$psu)
#' )
#' ```
#'
#' The reference population's design is handled separately and automatically;
#' see [ref_design()].
#'
#' @section Weights:
#' If the fit carried observation weights (`$weights` for `lm`,
#' `$prior.weights` for `glm`) they are used when computing `cov_X_s`, so the
#' covariance describes the same weighted population the coefficients do.
#' `glm` sets `prior.weights` to a vector of ones when no weights were supplied,
#' which is harmless.
#'
#' @examples
#' set.seed(7)
#' d <- data.frame(x1 = rnorm(300), x2 = rnorm(300))
#' d$y <- 1 + d$x1 + rnorm(300)
#' fit <- lm(y ~ x1 + x2, data = d)
#' rm1 <- as_reduced_model(fit)
#' names(rm1$theta_hat_s)
#'
#' @export
as_reduced_model <- function(reg_out,
                             var_type = c("robust", "standard", "best"),
                             keep_intercept = TRUE,
                             vcov = NULL,
                             ...) {
  var_type <- match.arg(var_type)
  if (!inherits(reg_out, c("lm", "glm"))) {
    stop(
      "`reg_out` must be an 'lm' or 'glm' fit. For other model classes, ",
      "assemble the reduced-model list yourself: it needs the elements ",
      "theta_hat_s, cov_theta_hat_s, cov_X_s and inds_theta_s.",
      call. = FALSE
    )
  }

  theta_hat_s <- stats::coef(reg_out)
  cov_theta_hat_s <- resolve_vcov(reg_out, vcov, var_type)
  X <- stats::model.matrix(reg_out)

  ## Standardise the intercept name before anything else, so that
  ## keep_intercept = FALSE removes exactly one, well-identified column.
  rename_intercept <- function(nm) {
    nm[is_intercept_name(nm)] <- "Intercept"
    nm
  }
  names(theta_hat_s) <- rename_intercept(names(theta_hat_s))
  dimnames(cov_theta_hat_s) <- list(
    rename_intercept(rownames(cov_theta_hat_s)),
    rename_intercept(colnames(cov_theta_hat_s))
  )
  colnames(X) <- rename_intercept(colnames(X))

  ## Weights, if any.
  wt <- NULL
  if (inherits(reg_out, "glm")) {
    if (!is.null(reg_out$prior.weights)) wt <- as.numeric(reg_out$prior.weights)
  } else if (!is.null(reg_out$weights)) {
    wt <- as.numeric(reg_out$weights)
  }
  if (!is.null(wt) && (all(wt == wt[1L]) || any(wt <= 0))) {
    wt <- NULL # constant or degenerate weights: use the unweighted covariance
  }

  ## cov_X_s never includes the intercept (a constant column has zero variance
  ## and would make the covariance matrix singular).
  X_slopes <- X[, !is_intercept_name(colnames(X)), drop = FALSE]
  cov_X <- if (is.null(wt)) {
    stats::cov(X_slopes)
  } else {
    stats::cov.wt(X_slopes, wt = wt)$cov
  }

  if (!keep_intercept) {
    keep <- !is_intercept_name(names(theta_hat_s))
    theta_hat_s <- theta_hat_s[keep]
    cov_theta_hat_s <- cov_theta_hat_s[
      names(theta_hat_s), names(theta_hat_s),
      drop = FALSE
    ]
  }

  structure(
    list(
      theta_hat_s = theta_hat_s,
      cov_theta_hat_s = cov_theta_hat_s,
      cov_X_s = cov_X,
      inds_theta_s = rep(1L, ncol(X_slopes)),
      converged = if (inherits(reg_out, "glm")) {
        isTRUE(reg_out$converged)
      } else {
        TRUE
      },
      n_obs = stats::nobs(reg_out),
      dgp = if (inherits(reg_out, "glm")) reg_out$family$family else "linear",
      fit = reg_out
    ),
    class = "reduced_model"
  )
}

#' @rdname as_reduced_model
#' @export
convert_reg_out_to_reduced_model_format <- function(reg_out,
                                                    var_type = "robust",
                                                    ...) {
  as_reduced_model(reg_out, var_type = var_type, ...)
}

#' Feed one genmetalite fit into another
#'
#' Converts the output of [calc_linear_beta()] into the reduced-model format, so
#' that a pooled estimate can itself act as one "study" in a second round of
#' pooling (for example, when combining consortia that have each already run
#' their own meta-analysis).
#'
#' @param genmetalite_out Output of [calc_linear_beta()].
#' @param cov_X Covariance matrix of the covariates in the pooled model. This
#'   must be supplied explicitly: the covariance stored inside
#'   `genmetalite_out` is study-specific and generally covers a different set of
#'   variables than the pooled coefficient vector.
#'
#' @return A list of class `reduced_model`.
#' @export
as_reduced_model_genmetalite <- function(genmetalite_out, cov_X) {
  if (is.null(genmetalite_out$linear_est)) {
    stop("`genmetalite_out` does not look like calc_linear_beta() output.",
      call. = FALSE
    )
  }
  theta_hat_s <- genmetalite_out$linear_est$beta
  cov_theta_hat_s <- genmetalite_out$linear_est$var

  if (missing(cov_X) || is.null(cov_X)) {
    stop(
      "`cov_X` is required: supply the covariance matrix of the covariates in ",
      "the pooled model (e.g. from cov_factors() on the reference data).",
      call. = FALSE
    )
  }
  nms <- drop_intercept_names(names(theta_hat_s))
  cov_X <- subset_by_name(cov_X, nms, "supplied cov_X")

  structure(
    list(
      theta_hat_s = theta_hat_s,
      cov_theta_hat_s = cov_theta_hat_s,
      cov_X_s = cov_X,
      inds_theta_s = rep(1L, length(nms)),
      converged = TRUE,
      n_obs = NA_integer_,
      dgp = genmetalite_out$dgp %||% "linear"
    ),
    class = "reduced_model"
  )
}

#' @rdname as_reduced_model_genmetalite
#' @export
convert_GENMETAL_out_to_reduced_model_format <- function(genmetalite_out,
                                                         cov_X) {
  as_reduced_model_genmetalite(genmetalite_out, cov_X)
}

#' Resolve a user-supplied vcov (matrix or function) against the default.
#' @noRd
resolve_vcov <- function(fit, vcov, var_type) {
  if (is.null(vcov)) {
    return(compute_var(fit, var_type))
  }
  v <- if (is.function(vcov)) vcov(fit) else vcov
  v <- as.matrix(v)
  k <- length(stats::coef(fit))
  if (nrow(v) != ncol(v)) {
    stop("`vcov` must be square.", call. = FALSE)
  }
  if (nrow(v) != k) {
    stop(sprintf(
      "`vcov` is %d x %d but the model has %d coefficient(s).",
      nrow(v), ncol(v), k
    ), call. = FALSE)
  }
  if (is.null(rownames(v))) {
    dimnames(v) <- list(names(stats::coef(fit)), names(stats::coef(fit)))
  }
  v
}
