#' Pool reduced regression models across studies
#'
#' The main entry point. Accepts either individual-level cohort data (in which
#' case the reduced models are fitted for you) or pre-fitted reduced models, plus
#' a reference population, and returns the mapped-and-pooled estimate together
#' with an inverse-variance weighted comparator and an empirical-Bayes blend.
#'
#' @param reduced_models Optional list of reduced models; see
#'   [run_regression()] and [as_reduced_model()]. If `NULL`, they are fitted
#'   from `cohort_data`.
#' @param cohort_data Optional list of individual-level cohort data frames. Each
#'   must contain the outcome in a column named `y`; every other column is
#'   treated as a covariate.
#' @param ref_data The reference population. Any of: a data frame; a list of
#'   `n_cohorts` data frames; a [ref_design()] object carrying weights, strata
#'   and clusters; or a `survey.design` object from the \pkg{survey} package
#'   (for example an NHANES design built with `svydesign()`). When a design is
#'   supplied its weights are used for the covariance and its strata and
#'   clusters for the variance, and `survey_strata`/`survey_psu` are not needed.
#' @param cov_ref_X Optional pre-computed reference covariance matrix, as an
#'   alternative to `ref_data`. Supplying this disables reweighting and drops
#'   the reference term from the variance.
#' @param cov_p_type `"partitioned"` (default) or `"reference"`; see
#'   [make_sigma_X_ls()].
#' @param dgp Data generating process for fitting the reduced models; see
#'   [run_regression()]. Ignored when `reduced_models` is supplied.
#' @param var_type Variance convention for the reduced models; see
#'   [compute_var()].
#' @param constr_vars Optional character vector of covariates to constrain to
#'   zero.
#' @param ref_weights Optional list of pre-computed reference weight vectors.
#' @param weighting.method `"ebal"` or `"IPW"`.
#' @param weighting.vars Variables to reweight the reference population on.
#'   Reweighting happens only if this is non-empty.
#' @param weighting.corr,weighting.max.moments,weighting.max.iterations
#'   Reweighting controls; see [ebal_weights()] and [ipw_weights()].
#' @param jacobian,check_jacobian,eps Passed to [linear_estimate()].
#' @param survey_strata,survey_psu Complex-survey design for the reference
#'   sample, as loose vectors. Only needed when `ref_data` is a plain data
#'   frame; a `ref_design()` or `survey.design` object carries this already, and
#'   these arguments are then ignored with a warning.
#' @param lonely_psu How to treat strata with a single PSU; see
#'   [survey_reference_cov()].
#'
#' @return A list of class `genmetalite` with
#'   \item{linear_est}{output of [linear_estimate()].}
#'   \item{IVW_est}{output of [IVW_estimate()], reordered to match.}
#'   \item{EB_est}{output of [EB_estimate()].}
#'   \item{reduced_models}{the reduced models used.}
#'   \item{empty_levels}{a data frame of factor levels a cohort measured but
#'     had no observations in; see [check_empty_levels()]. Zero rows is the
#'     usual case.}
#'   \item{sigma_X_ls, cov_cohort_X_s_ls, ref_weights}{intermediate quantities
#'     from [make_sigma_X_ls()]. With a survey reference, `ref_weights` are the
#'     combined sampling-and-balancing weights actually used.}
#'   \item{ref_strata, ref_psu}{the reference sampling design, or `NULL`.}
#'   \item{call, dgp, cov_p_type, constr_vars, weighting.vars}{the settings
#'     used.}
#'
#' @section Centering:
#' No centering of the input data is required or performed. The estimator works
#' with the covariance \eqn{\mathrm{Var}(X)}, which is invariant to location, and
#' with slope coefficients, which are likewise invariant. (Earlier versions of
#' this code mean-centered the numeric columns but not the dummy columns derived
#' from factors, which was inconsistent; it is now unnecessary either way.)
#'
#' @examples
#' set.seed(2026)
#' n_ref <- 4000
#' make_pop <- function(n) {
#'   x1 <- rnorm(n)
#'   x2 <- 0.6 * x1 + rnorm(n, sd = 0.8)
#'   x3 <- rnorm(n)
#'   data.frame(x1 = x1, x2 = x2, x3 = x3)
#' }
#' beta <- c(x1 = 1, x2 = -0.5, x3 = 0.3)
#' ref <- make_pop(n_ref)
#'
#' mk_cohort <- function(n, drop) {
#'   d <- make_pop(n)
#'   d$y <- 2 + as.matrix(d) %*% beta + rnorm(n)
#'   d[, setdiff(names(d), drop), drop = FALSE]
#' }
#' cohorts <- list(mk_cohort(800, "x3"), mk_cohort(800, "x2"))
#'
#' fit <- calc_linear_beta(cohort_data = cohorts, ref_data = ref)
#' round(fit$linear_est$beta, 3)
#' round(fit$IVW_est$beta, 3) # biased: each study omits a correlated covariate
#'
#' @section Survey references:
#' Passing a `survey.design` object makes the estimator design-aware end to end:
#' \eqn{\Sigma} becomes a survey-weighted covariance, and the reference term of
#' the variance uses the stratified ultimate-cluster linearization instead of the
#' independent-sampling formula. `fit$linear_est$design_df` then reports the
#' design degrees of freedom, which are often small enough that a \eqn{t}
#' reference distribution is more appropriate than a normal one.
#'
#' Reweighting composes with this. With `weighting.vars` set, the sampling
#' weights become the base weights of the entropy-balancing problem, so each
#' cohort gets one weight vector carrying both the sampling design and the
#' covariate shift, and the shared strata and clusters still apply. The
#' per-cohort influence contributions are summed across cohorts *before* the
#' design is applied, which keeps the correlation induced by their sharing the
#' same reference rows.
#'
#' @seealso [linear_estimate()], [make_sigma_X_ls()], [run_regression()],
#'   [ref_design()]
#' @export
calc_linear_beta <- function(reduced_models = NULL,
                             cohort_data = NULL,
                             ref_data = NULL,
                             cov_ref_X = NULL,
                             cov_p_type = c("partitioned", "reference"),
                             dgp = c("linear", "logistic", "poisson", "nb"),
                             var_type = c("robust", "standard", "best"),
                             constr_vars = NULL,
                             ref_weights = NULL,
                             weighting.method = c("ebal", "IPW"),
                             weighting.vars = NULL,
                             weighting.corr = TRUE,
                             weighting.max.moments = 1L,
                             weighting.max.iterations = 1000L,
                             jacobian = c("numerical", "analytic"),
                             check_jacobian = FALSE,
                             eps = 1e-6,
                             survey_strata = NULL,
                             survey_psu = NULL,
                             lonely_psu = c("remove", "adjust", "fail")) {
  cl <- match.call()
  cov_p_type <- match.arg(cov_p_type)
  dgp <- match.arg(dgp)
  var_type <- match.arg(var_type)
  weighting.method <- match.arg(weighting.method)
  jacobian <- match.arg(jacobian)
  lonely_psu <- match.arg(lonely_psu)

  if (is.null(reduced_models) && is.null(cohort_data)) {
    stop(
      "Supply either `cohort_data` (individual-level study data) or ",
      "`reduced_models` (pre-fitted reduced models).",
      call. = FALSE
    )
  }
  if (is.null(ref_data) && is.null(cov_ref_X)) {
    stop(
      "Supply either `ref_data` (individual-level reference data) or ",
      "`cov_ref_X` (a reference covariance matrix).",
      call. = FALSE
    )
  }

  if (!is.null(cohort_data) && is.data.frame(cohort_data)) {
    cohort_data <- list(cohort_data)
  }
  n_cohorts <- if (!is.null(reduced_models)) {
    length(reduced_models)
  } else {
    length(cohort_data)
  }
  if (n_cohorts == 0L) stop("No cohorts supplied.", call. = FALSE)

  ## ------------------------------------------------------------------
  ## Reweighting feasibility
  ## ------------------------------------------------------------------
  weighting <- length(weighting.vars) > 0L
  if (weighting && is.null(cohort_data)) {
    warning(
      "Reweighting needs individual-level cohort data, which were not ",
      "supplied, so `weighting.vars` is ignored. Either re-run with ",
      "`cohort_data`, or weight the reference data yourself and pass the ",
      "weights as `ref_weights` (or the weighted covariance as `cov_ref_X`).",
      call. = FALSE
    )
    weighting <- FALSE
  }
  if (weighting && !is.null(cov_ref_X)) {
    warning(
      "`cov_ref_X` was supplied, so there is nothing left to reweight and ",
      "`weighting.vars` is ignored. Reweight the reference data before ",
      "computing its covariance if that is what you want.",
      call. = FALSE
    )
    weighting <- FALSE
  }
  if (weighting && !is.null(ref_weights)) {
    warning(
      "Both `ref_weights` and `weighting.vars` were supplied; using the ",
      "supplied `ref_weights` and not recomputing them.",
      call. = FALSE
    )
    weighting <- FALSE
  }

  ## ------------------------------------------------------------------
  ## Harmonise factor levels and check the reference covers the studies
  ## ------------------------------------------------------------------
  ref <- normalise_reference(ref_data, n_cohorts)
  ref_list <- ref$data

  if (!is.null(ref) && (!is.null(ref$strata) || !is.null(ref$psu))) {
    if (!is.null(survey_strata) || !is.null(survey_psu)) {
      warning(
        "`ref_data` already carries a sampling design, so the separate ",
        "`survey_strata`/`survey_psu` arguments are ignored.",
        call. = FALSE
      )
    }
    survey_strata <- ref$strata
    survey_psu <- ref$psu
  }
  empty_levels <- data.frame(
    cohort = integer(0), variable = character(0), level = character(0),
    status = character(0), stringsAsFactors = FALSE
  )

  if (!is.null(cohort_data)) {
    cohort_data <- lapply(cohort_data, as.data.frame)
    for (k in seq_len(n_cohorts)) {
      if (!"y" %in% names(cohort_data[[k]])) {
        stop(sprintf(
          paste0(
            "cohort_data[[%d]] has no column named 'y'. The outcome must be ",
            "called 'y'; every other column is treated as a covariate."
          ), k
        ), call. = FALSE)
      }

      if (!is.null(ref_list)) {
        covs <- setdiff(names(cohort_data[[k]]), "y")
        check_vars(
          ref_list[[k]], cohort_data[[k]][, covs, drop = FALSE],
          action = "warn"
        )

        ## Align the level sets before anything expands a factor, so that a
        ## cohort which happens not to observe the reference's first level does
        ## not silently get a different baseline -- and hence dummy columns
        ## encoding different contrasts under the same names.
        cohort_data[[k]] <- harmonise_levels(cohort_data[[k]], ref_list[[k]])

        el <- check_empty_levels(
          ref_list[[k]], cohort_data[[k]][, covs, drop = FALSE],
          action = "warn", label = k
        )
        if (nrow(el) > 0L) {
          empty_levels <- rbind(empty_levels, cbind(cohort = k, el))
        }
      }

      cohort_data[[k]] <- droplevels(cohort_data[[k]])
    }
  }

  ## ------------------------------------------------------------------
  ## Fit reduced models if needed
  ## ------------------------------------------------------------------
  if (!is.null(reduced_models) && !is.null(cohort_data)) {
    message(
      "Both `reduced_models` and `cohort_data` were supplied. The supplied ",
      "models are used as-is; the raw data are used only for reweighting the ",
      "reference and for the study-block variance term."
    )
  } else if (is.null(reduced_models)) {
    reduced_models <- lapply(seq_len(n_cohorts), function(k) {
      d <- cohort_data[[k]]
      run_regression(
        y = d[["y"]],
        X_s = d[, setdiff(names(d), "y"), drop = FALSE],
        dgp = dgp,
        var_type = var_type
      )
    })
  }

  not_converged <- which(!vapply(
    reduced_models, function(r) isTRUE(r$converged), logical(1)
  ))
  if (length(not_converged) > 0L) {
    warning(
      "Reduced model(s) ", paste(not_converged, collapse = ", "),
      " did not converge; the pooled estimate inherits that problem.",
      call. = FALSE
    )
  }

  ## ------------------------------------------------------------------
  ## Covariance assembly and pooling
  ## ------------------------------------------------------------------
  sigma_out <- make_sigma_X_ls(
    n_cohorts = n_cohorts,
    cov_ref_X = cov_ref_X,
    reg_out = reduced_models,
    cov_p_type = cov_p_type,
    cohort_data = cohort_data,
    ref_data = if (is.null(ref)) NULL else ref_data,
    ref_weights = ref_weights,
    weighting = weighting,
    weighting.method = weighting.method,
    weighting.vars = weighting.vars,
    weighting.corr = weighting.corr,
    weighting.max.moments = weighting.max.moments,
    weighting.max.iterations = weighting.max.iterations
  )

  linear_est <- linear_estimate(
    n_cohorts = n_cohorts,
    reg_out = reduced_models,
    sigma_X_ls = sigma_out$sigma_X_ls,
    cov_cohort_X_s_ls = sigma_out$cov_cohort_X_s_ls,
    cohort_data = cohort_data,
    ref_pop_dummy = sigma_out$ref_pop_dummy,
    ref_weights = sigma_out$ref_weights,
    cov_p_type = cov_p_type,
    constr_vars = constr_vars,
    jacobian = jacobian,
    check_jacobian = check_jacobian,
    survey_strata = survey_strata,
    survey_psu = survey_psu,
    lonely_psu = lonely_psu,
    eps = eps
  )

  IVW_est <- IVW_estimate(reduced_models)
  keep <- names(linear_est$beta)
  IVW_est$beta <- IVW_est$beta[keep]
  IVW_est$var <- IVW_est$var[keep]
  IVW_est$n_studies <- IVW_est$n_studies[keep]

  EB_est <- EB_estimate(linear_est = linear_est, IVW_est = IVW_est)

  structure(
    list(
      linear_est = linear_est,
      IVW_est = IVW_est,
      EB_est = EB_est,
      reduced_models = reduced_models,
      empty_levels = empty_levels,
      sigma_X_ls = sigma_out$sigma_X_ls,
      cov_cohort_X_s_ls = sigma_out$cov_cohort_X_s_ls,
      ref_weights = sigma_out$ref_weights,
      ref_strata = sigma_out$ref_strata,
      ref_psu = sigma_out$ref_psu,
      weighting_converged = sigma_out$weighting_converged,
      call = cl,
      dgp = dgp,
      cov_p_type = cov_p_type,
      constr_vars = constr_vars,
      weighting.vars = weighting.vars,
      weighting.corr = weighting.corr,
      weighting.max.moments = weighting.max.moments
    ),
    class = "genmetalite"
  )
}

#' @export
print.genmetalite <- function(x, ...) {
  cat("<genmetalite>\n")
  cat(sprintf(
    "  %d studies,  dgp = %s,  cov_p_type = %s\n",
    length(x$reduced_models), x$dgp, x$cov_p_type
  ))
  cat("\nMapped-and-pooled estimate:\n")
  print(x$linear_est)
  cat("\nInverse-variance weighted comparator:\n")
  print(round(x$IVW_est$beta, 4L))
  invisible(x)
}

#' Coefficient table for a genmetalite fit
#'
#' @param object A `genmetalite` object from [calc_linear_beta()].
#' @param level Confidence level for the interval.
#' @param ... Unused.
#' @return A data frame with estimates, standard errors, z statistics, p-values
#'   and confidence limits.
#' @export
summary.genmetalite <- function(object, level = 0.95, ...) {
  b <- object$linear_est$beta
  se <- sqrt(pmax(diag(object$linear_est$var), 0))
  z <- ifelse(se > 0, b / se, NA_real_)
  crit <- stats::qnorm(1 - (1 - level) / 2)
  data.frame(
    term = names(b),
    estimate = as.numeric(b),
    std.error = as.numeric(se),
    statistic = as.numeric(z),
    p.value = 2 * stats::pnorm(-abs(z)),
    conf.low = as.numeric(b - crit * se),
    conf.high = as.numeric(b + crit * se),
    ivw.estimate = as.numeric(object$IVW_est$beta[names(b)]),
    row.names = NULL,
    stringsAsFactors = FALSE
  )
}
