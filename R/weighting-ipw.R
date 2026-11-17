#' Reweight a reference population to a cohort by inverse probability weighting
#'
#' Fits a logistic "membership" model that discriminates cohort rows from
#' reference rows, and returns odds weights for the reference rows. Weighting
#' the reference population by these odds makes its covariate distribution
#' resemble the cohort's, which is the same target as [ebal_weights()] reached by
#' a different route.
#'
#' @param cohort_data Data frame of individual-level cohort covariates.
#' @param ref_data Data frame of individual-level reference covariates.
#' @param weighting.vars Character vector of variables to include in the
#'   membership model. `NULL` (default) uses every variable shared by the two
#'   datasets. Variables absent from either dataset are dropped with a warning
#'   rather than causing an error.
#' @param weighting.corr If `TRUE` (default), add all pairwise interactions
#'   among the weighting variables, so second-order covariate structure is
#'   balanced too.
#' @param weighting.max.moments Highest polynomial term of each continuous
#'   weighting variable to include (1 = linear only).
#' @param stabilize If `TRUE` (default), multiply the odds weights by the
#'   marginal odds \eqn{P(\mathrm{cohort}) / P(\mathrm{reference})} so that the
#'   weights average roughly one. Stabilization does not change the weighted
#'   covariance (a common scale factor cancels), but it makes the weights easier
#'   to inspect.
#' @param trim Optional two-element vector of quantiles, e.g. `c(0.01, 0.99)`,
#'   at which to winsorise the weights. `NULL` (default) means no trimming.
#' @param base.weight Optional base weights for the reference rows. **Survey
#'   weights go here.** They enter twice: as regression weights in the
#'   membership model, so the reference represents its population rather than
#'   its sample, and as a multiplier on the resulting odds weights, so the
#'   returned vector carries both the sampling design and the covariate shift.
#' @param maxit Maximum IRLS iterations for the membership model.
#'
#' @return A list with
#'   \item{weights}{numeric vector of length `nrow(ref_data)`.}
#'   \item{converged}{logical, from the membership `glm`.}
#'   \item{ps}{the fitted membership probabilities for the reference rows.}
#'   \item{model}{the fitted `glm` object.}
#'
#' @details
#' Writing \eqn{T = 1} for a cohort row and \eqn{T = 0} for a reference row, and
#' \eqn{e(x) = P(T = 1 \mid X = x)}, the density ratio between the cohort and
#' reference covariate distributions is proportional to \eqn{e(x)/(1 - e(x))}.
#' Weighting reference row \eqn{i} by that odds ratio therefore tilts the
#' reference distribution onto the cohort's, which is exactly what
#' [make_sigma_X_ls()] needs before computing a cohort-matched covariance.
#'
#' Extreme weights are the usual failure mode: if a few reference rows have
#' membership probabilities near one they dominate the weighted covariance.
#' A warning is issued when the effective sample size
#' \eqn{(\sum w)^2 / \sum w^2} falls below 10% of the reference sample size.
#'
#' @seealso [ebal_weights()]
#'
#' @examples
#' set.seed(3)
#' ref <- data.frame(age = rnorm(2000), bmi = rnorm(2000))
#' coh <- data.frame(age = rnorm(400, 0.6), bmi = rnorm(400, 0.3))
#' w <- ipw_weights(coh, ref, weighting.corr = FALSE)
#' # weighted reference mean moves toward the cohort mean
#' c(
#'   unweighted = mean(ref$age),
#'   weighted = weighted.mean(ref$age, w$weights),
#'   cohort = mean(coh$age)
#' )
#'
#' @export
ipw_weights <- function(cohort_data,
                        ref_data,
                        weighting.vars = NULL,
                        weighting.corr = TRUE,
                        weighting.max.moments = 1L,
                        stabilize = TRUE,
                        trim = NULL,
                        base.weight = NULL,
                        maxit = 100L) {
  cohort_data <- as.data.frame(cohort_data)
  ref_data <- as.data.frame(ref_data)

  shared <- intersect(names(cohort_data), names(ref_data))
  if (length(shared) == 0L) {
    stop("`cohort_data` and `ref_data` share no variable names.", call. = FALSE)
  }

  if (is.null(weighting.vars)) {
    weighting.vars <- shared
  } else {
    ## A variable the cohort did not measure cannot be balanced, and that is
    ## the normal situation here -- studies measure different subsets. Only a
    ## name absent from BOTH datasets indicates a typo worth warning about.
    unknown <- setdiff(
      weighting.vars, union(names(cohort_data), names(ref_data))
    )
    if (length(unknown) > 0L) {
      warning(
        "These weighting variables are in neither dataset and were dropped: ",
        paste(unknown, collapse = ", "), ".",
        call. = FALSE
      )
    }
    weighting.vars <- intersect(weighting.vars, shared)
    if (length(weighting.vars) == 0L) {
      stop("No usable weighting variables remain.", call. = FALSE)
    }
  }

  n_k <- nrow(cohort_data)
  n_ref <- nrow(ref_data)
  treatment <- c(rep(0L, n_ref), rep(1L, n_k))

  X <- rbind(
    ref_data[, weighting.vars, drop = FALSE],
    cohort_data[, weighting.vars, drop = FALSE]
  )
  X <- droplevels(X)

  ok <- stats::complete.cases(X)
  if (!all(ok)) {
    stop(sprintf(
      paste0(
        "%d of %d rows have missing values in the weighting variables. ",
        "Drop or impute them before calling ipw_weights()."
      ),
      sum(!ok), length(ok)
    ), call. = FALSE)
  }

  ## Build the membership-model formula --------------------------------------
  numeric_vars <- weighting.vars[
    vapply(X[weighting.vars], is.numeric, logical(1))
  ]
  terms_rhs <- weighting.vars
  if (weighting.max.moments > 1L && length(numeric_vars) > 0L) {
    for (m in 2:weighting.max.moments) {
      terms_rhs <- c(terms_rhs, sprintf("I(%s^%d)", numeric_vars, m))
    }
  }
  rhs <- paste(terms_rhs, collapse = " + ")
  if (isTRUE(weighting.corr) && length(weighting.vars) >= 2L) {
    pairs <- utils::combn(weighting.vars, 2L, simplify = FALSE)
    rhs <- paste(
      rhs, "+",
      paste(vapply(pairs, function(p) paste(p, collapse = ":"), character(1)),
        collapse = " + "
      )
    )
  }
  form <- stats::as.formula(paste("treatment ~", rhs))

  if (!is.null(base.weight)) {
    base.weight <- as.numeric(base.weight)
    if (length(base.weight) != n_ref) {
      stop(sprintf(
        "`base.weight` has length %d but `ref_data` has %d rows.",
        length(base.weight), n_ref
      ), call. = FALSE)
    }
    if (any(!is.finite(base.weight)) || any(base.weight < 0) ||
      sum(base.weight) <= 0) {
      stop("`base.weight` must be finite, non-negative and not all zero.",
        call. = FALSE
      )
    }
  }

  fit_df <- cbind(treatment = treatment, X)
  ## Reference rows carry their sampling weights so the membership model
  ## contrasts the cohort against the population the reference represents.
  ## quasibinomial rather than binomial: non-integer weights are legitimate
  ## here, and binomial() would emit a "non-integer #successes" warning on
  ## every call while fitting exactly the same coefficients.
  fit_family <- if (is.null(base.weight)) {
    stats::binomial(link = "logit")
  } else {
    stats::quasibinomial(link = "logit")
  }
  model_weights <- if (is.null(base.weight)) {
    rep(1, nrow(fit_df))
  } else {
    ## Scale the reference weights to the reference sample size, so the two
    ## groups keep their relative influence and the fit stays well conditioned.
    c(base.weight * (n_ref / sum(base.weight)), rep(1, n_k))
  }
  fit_df$.model_weight <- model_weights

  ps_model <- stats::glm(
    form,
    family = fit_family,
    data = fit_df,
    weights = fit_df$.model_weight,
    control = list(maxit = maxit)
  )

  if (!isTRUE(ps_model$converged)) {
    warning(
      "The membership model did not converge. This usually means the cohort ",
      "and reference populations are (nearly) perfectly separable on these ",
      "covariates, in which case reweighting cannot be trusted.",
      call. = FALSE
    )
  }

  e <- stats::predict(ps_model, type = "response")[seq_len(n_ref)]

  ## Guard against numerical 1s before dividing.
  eps <- .Machine$double.eps^0.5
  e <- pmin(pmax(e, eps), 1 - eps)
  w <- e / (1 - e)

  if (!is.null(base.weight)) {
    ## Odds weights tilt the reference POPULATION onto the cohort, so they
    ## multiply the sampling weights rather than replacing them.
    w <- w * base.weight
  }

  if (isTRUE(stabilize)) {
    ## Rescale to the total the reference started with: n_ref for an unweighted
    ## reference (so weights average one), or the population total for a survey
    ## reference (so that interpretation survives the reweighting). A common
    ## scale factor cancels out of every weighted covariance, so this changes
    ## nothing downstream -- it only makes the weights readable.
    target_total <- if (is.null(base.weight)) n_ref else sum(base.weight)
    w <- w * (target_total / sum(w))
  }

  if (!is.null(trim)) {
    if (length(trim) != 2L || any(trim < 0) || any(trim > 1) ||
      trim[1L] >= trim[2L]) {
      stop("`trim` must be two increasing probabilities in [0, 1].",
        call. = FALSE
      )
    }
    qs <- stats::quantile(w, probs = trim, names = FALSE)
    w <- pmin(pmax(w, qs[1L]), qs[2L])
  }

  n_eff <- sum(w)^2 / sum(w^2)
  if (n_eff < 0.1 * n_ref) {
    warning(sprintf(
      paste0(
        "Effective reference sample size is %.0f out of %d rows (%.1f%%). ",
        "A few rows dominate the weighted covariance; consider trimming, or ",
        "using ebal_weights(), which bounds weight dispersion by construction."
      ),
      n_eff, n_ref, 100 * n_eff / n_ref
    ), call. = FALSE)
  }

  list(
    weights = as.numeric(w),
    converged = isTRUE(ps_model$converged),
    ps = e,
    model = ps_model
  )
}
