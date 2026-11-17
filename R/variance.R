#' Variance-covariance matrix of fitted regression coefficients
#'
#' A thin wrapper that returns the coefficient covariance matrix of a fitted
#' `lm` or `glm` object under one of three conventions.
#'
#' @param reg_out A fitted model object, typically from [stats::lm()] or
#'   [stats::glm()].
#' @param var_type One of:
#'   \describe{
#'     \item{`"robust"`}{(default) the Huber-White sandwich estimator,
#'       [sandwich::sandwich()]. Consistent under model misspecification, which
#'       is the relevant setting here because the reduced models are deliberately
#'       misspecified relative to the full model.}
#'     \item{`"standard"`}{the model-based estimator, [stats::vcov()].}
#'     \item{`"best"`}{sandwich estimates for most coefficients, but model-based
#'       estimates substituted for indicator columns with a small effective
#'       sample size, where the sandwich estimator is known to be badly
#'       downward-biased. See Details.}
#'   }
#' @param min_obs For `var_type = "best"` with a linear model: indicator columns
#'   with fewer than `min_obs` non-zero entries use the model-based variance.
#' @param min_eff For `var_type = "best"` with a generalized linear model:
#'   indicator columns whose effective sample size
#'   \eqn{n_0 n_1 / (n_0 + n_1)} within outcome levels falls below `min_eff` use
#'   the model-based variance.
#'
#' @details
#' `var_type = "best"` is a heuristic, not an established estimator, and is
#' flagged as experimental. It exists because in the simulation settings that
#' motivated this package some studies contain rare factor levels for which the
#' sandwich variance is unusably small; substituting the model-based variance for
#' just those coefficients stabilises the generalized least squares pooling
#' downstream. The thresholds are deliberately exposed as arguments rather than
#' hard-coded. It only applies to columns that look like indicators (all entries
#' 0 or 1); continuous columns always keep the sandwich variance.
#'
#' @return A symmetric numeric matrix with dimnames taken from the coefficient
#'   names.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
#' d$y <- 1 + 0.5 * d$x1 - 0.3 * d$x2 + rnorm(200)
#' fit <- lm(y ~ x1 + x2, data = d)
#' round(compute_var(fit, "robust"), 4)
#' round(compute_var(fit, "standard"), 4)
#'
#' @seealso [run_regression()]
#' @export
compute_var <- function(reg_out,
                        var_type = c("robust", "standard", "best"),
                        min_obs = 100L,
                        min_eff = 10L) {
  var_type <- match.arg(var_type)

  if (var_type == "robust") {
    return(sandwich::sandwich(reg_out))
  }
  if (var_type == "standard") {
    return(stats::vcov(reg_out))
  }

  ## var_type == "best" ------------------------------------------------------
  cov_robust <- sandwich::sandwich(reg_out)
  cov_model <- stats::vcov(reg_out)

  ## vcov() can carry all-NA rows/columns for aliased coefficients; drop them.
  all_na <- apply(cov_model, 2L, function(x) all(is.na(x)))
  if (any(all_na)) {
    cov_model <- cov_model[!all_na, !all_na, drop = FALSE]
  }

  X <- stats::model.matrix(reg_out)
  indicator <- apply(X, 2L, function(col) all(col %in% c(0, 1)))
  indicator_names <- colnames(X)[indicator]

  low_obs_names <- character(0)

  if (inherits(reg_out, "glm")) {
    y <- stats::model.response(stats::model.frame(reg_out))
    y_levels <- unique(stats::na.omit(y))
    if (length(y_levels) == 2L && length(indicator_names) > 0L) {
      ## Effective sample size for an indicator column, within a binary outcome:
      ## n0 * n1 / (n0 + n1) where nk = # rows with indicator == 1 and y == k.
      counts <- vapply(indicator_names, function(nm) {
        col <- X[, nm]
        c(
          sum(col == 1 & y == y_levels[1L], na.rm = TRUE),
          sum(col == 1 & y == y_levels[2L], na.rm = TRUE)
        )
      }, numeric(2L))
      denom <- colSums(counts)
      eff <- ifelse(denom > 0, counts[1L, ] * counts[2L, ] / denom, 0)
      low_obs_names <- indicator_names[eff < min_eff]
    } else if (length(y_levels) != 2L) {
      warning(
        "var_type = 'best' expects a binary outcome for glm fits; ",
        "returning the robust sandwich variance unchanged.",
        call. = FALSE
      )
      return(cov_robust)
    }
  } else if (inherits(reg_out, "lm")) {
    if (length(indicator_names) > 0L) {
      n_obs_per_level <- colSums(X[, indicator_names, drop = FALSE])
      low_obs_names <- indicator_names[n_obs_per_level < min_obs]
    }
  } else {
    warning(
      "var_type = 'best' only knows how to handle 'lm' and 'glm' fits; ",
      "returning the robust sandwich variance unchanged.",
      call. = FALSE
    )
    return(cov_robust)
  }

  ## Only substitute names present in BOTH matrices (aliased coefficients may
  ## have been dropped from one of them).
  low_obs_names <- intersect(
    low_obs_names,
    intersect(rownames(cov_robust), rownames(cov_model))
  )

  if (length(low_obs_names) > 0L) {
    cov_robust[low_obs_names, colnames(cov_model)] <-
      cov_model[low_obs_names, , drop = FALSE]
    cov_robust[rownames(cov_model), low_obs_names] <-
      cov_model[, low_obs_names, drop = FALSE]
  }

  symmetrise(cov_robust)
}
