#' Fit a reduced regression model
#'
#' Fits one study's reduced model and returns it in the "reduced model" list
#' format consumed by [linear_estimate()] and [calc_linear_beta()].
#'
#' @param y Numeric outcome vector.
#' @param X_s Data frame or matrix of the covariates measured in this study.
#'   Factors are expanded internally; do not pass an intercept column.
#' @param inds_theta_s Either `1` (the default, meaning "use every column of
#'   `X_s`") or a 0/1 vector of length `ncol(X_s)` selecting the columns to keep.
#' @param dgp Data generating process / link. One of `"linear"`, `"logistic"`,
#'   `"poisson"`, `"nb"` (negative binomial, requires the \pkg{MASS} package).
#' @param var_type Passed to [compute_var()]; `"robust"` by default. If the
#'   study is itself a clustered or complex sample, the default sandwich
#'   understates the coefficient covariance; fit the model yourself and use
#'   [as_reduced_model()] with its `vcov` argument instead.
#' @param maxit Maximum IRLS iterations for the `glm` fits.
#'
#' @return A list of class `reduced_model` with components
#'   \item{theta_hat_s}{named coefficient vector, **including** the intercept
#'     (named `Intercept`).}
#'   \item{cov_theta_hat_s}{coefficient covariance matrix with matching
#'     dimnames.}
#'   \item{cov_X_s}{empirical covariance of the expanded covariates
#'     (no intercept).}
#'   \item{inds_theta_s}{the column selector actually used.}
#'   \item{converged}{logical.}
#'   \item{n_obs}{number of observations used.}
#'   \item{dgp}{the `dgp` argument.}
#'
#' @section Intercept convention:
#' `run_regression()` **keeps** the intercept in `theta_hat_s` and
#' `cov_theta_hat_s`, but **excludes** it from `cov_X_s`. [linear_estimate()]
#' works in the slope-only basis and drops the intercept itself, so both
#' conventions are safe; what matters is that the same convention is used for
#' every study. [as_reduced_model()] follows the same convention.
#'
#' @examples
#' set.seed(42)
#' n <- 500
#' d <- data.frame(x1 = rnorm(n), x2 = rnorm(n))
#' y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(n)
#' rm1 <- run_regression(y, d, dgp = "linear")
#' round(rm1$theta_hat_s, 3)
#'
#' @seealso [as_reduced_model()], [calc_linear_beta()]
#' @export
run_regression <- function(y,
                           X_s,
                           inds_theta_s = 1,
                           dgp = c("linear", "logistic", "poisson", "nb"),
                           var_type = c("robust", "standard", "best"),
                           maxit = 300L) {
  dgp <- match.arg(dgp)
  var_type <- match.arg(var_type)

  X_s <- as.data.frame(X_s)
  n_obs <- length(y)
  if (n_obs != nrow(X_s)) {
    stop(sprintf(
      "length(y) is %d but X_s has %d rows.", n_obs, nrow(X_s)
    ), call. = FALSE)
  }

  ## Column selection ---------------------------------------------------------
  k_raw <- ncol(X_s)
  if (identical(inds_theta_s, 1) || identical(inds_theta_s, 1L)) {
    inds_theta_s <- rep(1L, k_raw)
  } else {
    if (length(inds_theta_s) != k_raw) {
      stop(sprintf(
        paste0(
          "`inds_theta_s` has length %d but X_s has %d columns. It must be ",
          "either the scalar 1 (use all columns) or a 0/1 vector with one ",
          "entry per column of X_s."
        ),
        length(inds_theta_s), k_raw
      ), call. = FALSE)
    }
    if (!all(inds_theta_s %in% c(0, 1))) {
      stop("`inds_theta_s` must contain only 0s and 1s.", call. = FALSE)
    }
    if (sum(inds_theta_s) == 0L) {
      stop("`inds_theta_s` selects no columns.", call. = FALSE)
    }
    if (sum(inds_theta_s) < k_raw) {
      X_s <- X_s[, inds_theta_s == 1L, drop = FALSE]
    }
  }

  ## Empirical covariance of the (expanded) covariates ------------------------
  cov_out <- cov_factors(X_s)
  emp_cov_X_s <- cov_out$cov_X_s
  X_mat <- cov_out$X_dummy_mat

  ## Fit ----------------------------------------------------------------------
  fit_data <- list(y = y, X_mat = X_mat)
  model_s <- switch(dgp,
    linear = stats::lm(y ~ X_mat, data = fit_data),
    logistic = stats::glm(y ~ X_mat,
      family = stats::binomial(link = "logit"),
      data = fit_data, control = list(maxit = maxit)
    ),
    poisson = stats::glm(y ~ X_mat,
      family = stats::poisson(link = "log"),
      data = fit_data, control = list(maxit = maxit)
    ),
    nb = {
      if (!requireNamespace("MASS", quietly = TRUE)) {
        stop(
          "dgp = 'nb' requires the 'MASS' package. Install it with ",
          "install.packages('MASS').",
          call. = FALSE
        )
      }
      MASS::glm.nb(y ~ X_mat, data = fit_data, maxit = maxit, init.theta = 1)
    }
  )

  theta_hat_s <- stats::coef(model_s)
  cov_theta_hat_s <- compute_var(model_s, var_type)

  ## Restore the covariate names. lm()/glm() prepend the name of the matrix
  ## argument to each column name -- except when the matrix has exactly ONE
  ## column, where the column name is dropped entirely and the coefficient is
  ## called just "X_mat". Stripping the prefix with a regular expression
  ## therefore yields an empty name for single-covariate studies, which then
  ## propagates as a nameless covariate through the whole pipeline. Mapping the
  ## generated names back explicitly avoids that.
  prefixed <- if (ncol(X_mat) == 1L) {
    "X_mat"
  } else {
    paste0("X_mat", colnames(X_mat))
  }
  lookup <- c(
    stats::setNames("Intercept", "(Intercept)"),
    stats::setNames(colnames(X_mat), prefixed)
  )
  clean <- function(nm) {
    out <- unname(lookup[nm])
    ifelse(is.na(out), nm, out)
  }
  names(theta_hat_s) <- clean(names(theta_hat_s))
  dimnames(cov_theta_hat_s) <- list(
    clean(rownames(cov_theta_hat_s)),
    clean(colnames(cov_theta_hat_s))
  )

  ## lm() objects have no $converged element; a completed lm() fit has, by
  ## definition, converged.
  converged <- if (inherits(model_s, "glm")) {
    isTRUE(model_s$converged)
  } else {
    TRUE
  }
  if (!converged) {
    warning(
      "The reduced model did not converge; its coefficients and variance ",
      "should not be trusted.",
      call. = FALSE
    )
  }

  structure(
    list(
      theta_hat_s = theta_hat_s,
      cov_theta_hat_s = cov_theta_hat_s,
      cov_X_s = emp_cov_X_s,
      inds_theta_s = inds_theta_s,
      converged = converged,
      n_obs = n_obs,
      dgp = dgp,
      fit = model_s
    ),
    class = "reduced_model"
  )
}

#' @export
print.reduced_model <- function(x, ...) {
  cat(sprintf(
    "<reduced_model>  dgp = %s,  n = %d,  %d coefficient(s)%s\n",
    x$dgp %||% "?", x$n_obs %||% NA_integer_, length(x$theta_hat_s),
    if (isFALSE(x$converged)) "  [DID NOT CONVERGE]" else ""
  ))
  print(round(x$theta_hat_s, 4L))
  invisible(x)
}
