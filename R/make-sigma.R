#' Assemble per-study covariance matrices for the full covariate set
#'
#' Each study measured only some of the covariates in the target full model, so
#' no study can estimate the whole covariate covariance matrix. This function
#' builds, for every study, a covariance matrix over the *full* covariate set by
#' combining what the study can estimate with what has to be borrowed from the
#' reference population.
#'
#' @param n_cohorts Number of studies.
#' @param cov_ref_X Optional pre-computed reference covariance. Either a single
#'   matrix (used for every study) or a list of `n_cohorts` matrices. If
#'   supplied, `ref_data` is not needed and no reweighting is possible.
#' @param reg_out List of reduced models (see [run_regression()]). Used for the
#'   study-side covariance blocks and for the covariate names.
#' @param cov_p_type How to fill the block among covariates the study *did*
#'   measure:
#'   \describe{
#'     \item{`"partitioned"`}{(default) use the study's own empirical
#'       covariance. Cross-blocks with unmeasured covariates still come from the
#'       reference.}
#'     \item{`"reference"`}{take every block from the reference population.}
#'   }
#' @param cohort_data Optional list of individual-level cohort data frames.
#'   Required only for reweighting.
#' @param ref_data Individual-level reference data. Either a single data frame
#'   (reused for every study), a list of `n_cohorts` data frames, a
#'   [ref_design()] object, or a `survey.design` object from the \pkg{survey}
#'   package. When a design is supplied, its sampling weights are used for the
#'   covariance and its strata and clusters are returned for the variance.
#' @param ref_weights Optional list of pre-computed reference weight vectors,
#'   one per study. Supplying these bypasses `weighting`.
#' @param weighting Logical; reweight the reference to match each cohort before
#'   computing its covariance.
#' @param weighting.method `"ebal"` (entropy balancing, see [ebal_weights()]) or
#'   `"IPW"` (see [ipw_weights()]).
#' @param weighting.vars,weighting.corr,weighting.max.moments,weighting.max.iterations
#'   Passed through to the weighting function.
#'
#' @return A list with
#'   \item{sigma_X_ls}{list of `n_cohorts` covariance matrices over the full
#'     covariate set, all in the same canonical variable order.}
#'   \item{cov_cohort_X_s_ls}{list of study-side covariance matrices.}
#'   \item{ref_weights}{list of reference weight vectors (`NULL` entries when
#'     unweighted).}
#'   \item{ref_pop_dummy}{list of dummy-expanded reference design matrices, or
#'     `NULL` if only `cov_ref_X` was supplied.}
#'   \item{ref_strata,ref_psu}{the reference sampling design, or `NULL`.}
#'   \item{design_weights}{the reference's own sampling weights before any
#'     reweighting, or `NULL`.}
#'   \item{weighting_converged}{logical vector, one entry per study.}
#'
#' @section Intercept:
#' All matrices are returned in the *slope-only* basis: any intercept column is
#' removed. A constant column has zero variance, so keeping it makes every
#' covariance matrix singular; [linear_estimate()] handles the intercept
#' separately.
#'
#' @section Survey references and reweighting:
#' When the reference carries sampling weights and `weighting = TRUE`, those
#' weights become the *base weights* of the entropy-balancing problem, so the
#' returned weights tilt the population the reference represents rather than the
#' raw sample. One weight vector per cohort therefore carries both the sampling
#' design and the covariate shift. The strata and clusters are unaffected by
#' reweighting and are passed through unchanged.
#'
#' @seealso [linear_estimate()], [calc_linear_beta()], [ref_design()]
#' @export
make_sigma_X_ls <- function(n_cohorts,
                            cov_ref_X = NULL,
                            reg_out = NULL,
                            cov_p_type = c("partitioned", "reference"),
                            cohort_data = NULL,
                            ref_data = NULL,
                            ref_weights = NULL,
                            weighting = FALSE,
                            weighting.method = c("ebal", "IPW"),
                            weighting.vars = NULL,
                            weighting.corr = TRUE,
                            weighting.max.moments = 1L,
                            weighting.max.iterations = 200L) {
  cov_p_type <- match.arg(cov_p_type)
  weighting.method <- match.arg(weighting.method)

  if (is.null(cov_ref_X) && is.null(ref_data)) {
    stop(
      "Supply either `ref_data` (individual-level reference data) or ",
      "`cov_ref_X` (a pre-computed reference covariance matrix).",
      call. = FALSE
    )
  }

  ref <- normalise_reference(ref_data, n_cohorts)
  ref_data <- ref$data
  design_weights <- ref$weights
  ref_strata <- ref$strata
  ref_psu <- ref$psu
  cohort_data <- as_cohort_list(cohort_data, n_cohorts, "cohort_data")

  ## `converged` defaults to TRUE: with no weighting there is nothing to
  ## converge, and downstream code should not see NA.
  weighting_converged <- rep(TRUE, n_cohorts)
  ref_pop_dummy <- NULL

  ## ---------------------------------------------------------------------
  ## 1. Reference covariance, optionally reweighted
  ## ---------------------------------------------------------------------
  if (!is.null(ref_data)) {
    if (!is.null(cov_ref_X)) {
      warning(
        "Both `ref_data` and `cov_ref_X` were supplied; using `cov_ref_X` and ",
        "ignoring `ref_data` for the covariance. `ref_data` is still used for ",
        "the influence-function variance term.",
        call. = FALSE
      )
    }

    ## Expand factors once per study and drop any intercept column.
    ref_pop_dummy <- lapply(ref_data, function(df) {
      m <- dummy_matrix(df, add_intercept = FALSE)
      m[, !is_intercept_name(colnames(m)), drop = FALSE]
    })

    ## Keep only the covariates the studies actually use. The full model is
    ## defined by the union of the reduced models' coefficients, so a reference
    ## column no study measured contributes nothing -- and a survey design's
    ## `variables` slot routinely carries stratum, PSU, weight and identifier
    ## columns that would otherwise be treated as covariates.
    if (!is.null(reg_out)) {
      needed <- unique(unlist(lapply(reg_out, function(r) {
        drop_intercept_names(names(r$theta_hat_s))
      })))
      ref_pop_dummy <- lapply(ref_pop_dummy, function(m) {
        keep <- intersect(colnames(m), needed)
        if (length(keep) == 0L) {
          stop(
            "None of the reduced models' covariates appear in the reference ",
            "data. Check that the reference uses the same variable names, and ",
            "the same factor coding, as the cohorts.",
            call. = FALSE
          )
        }
        m[, keep, drop = FALSE]
      })
    }

    ref_dims <- vapply(ref_pop_dummy, ncol, integer(1))
    if (length(unique(ref_dims)) > 1L) {
      stop(
        "The dummy-expanded reference data have different numbers of columns ",
        "across cohorts (", paste(ref_dims, collapse = ", "), "). This means a ",
        "factor level is absent from one of the reference subsets. Harmonise ",
        "the factor levels before calling make_sigma_X_ls().",
        call. = FALSE
      )
    }

    ## Weights -------------------------------------------------------------
    ## With a survey reference and no reweighting, the design weights ARE the
    ## reference weights.
    if (is.null(ref_weights) && !isTRUE(weighting)) {
      ref_weights <- design_weights
    }

    if (is.null(ref_weights) && isTRUE(weighting)) {
      if (is.null(cohort_data)) {
        stop(
          "`weighting = TRUE` needs individual-level `cohort_data` to weight ",
          "the reference towards.",
          call. = FALSE
        )
      }
      wfun <- if (weighting.method == "ebal") {
        function(k) {
          ebal_weights(
            cohort_data = cohort_data[[k]],
            ref_data = ref_data[[k]],
            weighting.vars = weighting.vars,
            weighting.corr = weighting.corr,
            weighting.max.moments = weighting.max.moments,
            ## NOTE: the caller's `weighting.max.iterations` maps onto
            ## ebal_weights()'s `max.iterations`. In the original script this
            ## argument fell into `...` and was silently ignored.
            max.iterations = weighting.max.iterations,
            ## Survey weights become the base weights of the balancing
            ## problem, so the result tilts the population the reference
            ## represents rather than the raw sample.
            base.weight = if (is.null(design_weights)) {
              NULL
            } else {
              design_weights[[k]]
            }
          )
        }
      } else {
        function(k) {
          ipw_weights(
            cohort_data = cohort_data[[k]],
            ref_data = ref_data[[k]],
            weighting.vars = if (is.list(weighting.vars)) {
              unique(unlist(weighting.vars))
            } else {
              weighting.vars
            },
            weighting.corr = weighting.corr,
            weighting.max.moments = weighting.max.moments,
            base.weight = if (is.null(design_weights)) {
              NULL
            } else {
              design_weights[[k]]
            }
          )
        }
      }
      wout <- lapply(seq_len(n_cohorts), wfun)
      ref_weights <- lapply(wout, `[[`, "weights")
      weighting_converged <- vapply(
        wout, function(z) isTRUE(z$converged), logical(1)
      )
    }

    if (is.null(cov_ref_X)) {
      cov_ref_X <- lapply(seq_len(n_cohorts), function(k) {
        w <- if (is.null(ref_weights)) NULL else ref_weights[[k]]
        Xk <- ref_pop_dummy[[k]]
        if (is.null(w)) {
          stats::cov(Xk)
        } else {
          if (length(w) != nrow(Xk)) {
            stop(sprintf(
              paste0(
                "ref_weights[[%d]] has length %d but the reference data for ",
                "cohort %d have %d rows."
              ),
              k, length(w), k, nrow(Xk)
            ), call. = FALSE)
          }
          ## cov.wt() on the whole list rather than one element was a bug in the
          ## original script; it silently produced a covariance of the wrong
          ## object whenever weighting was switched on.
          stats::cov.wt(Xk, wt = w, center = TRUE)$cov
        }
      })
    }
  }

  ## Normalise cov_ref_X into a list of length n_cohorts.
  if (is.matrix(cov_ref_X)) {
    cov_ref_X <- rep(list(cov_ref_X), n_cohorts)
  }
  if (length(cov_ref_X) != n_cohorts) {
    stop(sprintf(
      "`cov_ref_X` has %d element(s) but %d cohorts were supplied.",
      length(cov_ref_X), n_cohorts
    ), call. = FALSE)
  }
  cov_ref_X <- lapply(cov_ref_X, function(m) {
    keep <- !is_intercept_name(colnames(m))
    m[keep, keep, drop = FALSE]
  })

  if (is.null(ref_weights)) ref_weights <- vector("list", n_cohorts)

  ## ---------------------------------------------------------------------
  ## 2. Study-side covariance
  ## ---------------------------------------------------------------------
  cov_cohort_X_s <- NULL
  if (!is.null(reg_out)) {
    cov_cohort_X_s <- lapply(reg_out, function(r) {
      m <- r$cov_X_s
      keep <- !is_intercept_name(colnames(m))
      m[keep, keep, drop = FALSE]
    })
  }

  var_names_ref_pop <- unique(unlist(lapply(cov_ref_X, colnames)))
  if (is.null(var_names_ref_pop)) {
    stop("The reference covariance matrix has no dimnames.", call. = FALSE)
  }

  if (is.null(cov_cohort_X_s) || cov_p_type == "reference") {
    if (is.null(reg_out)) {
      stop(
        "`reg_out` is required: the covariate names of each reduced model ",
        "determine which covariance blocks that study contributes.",
        call. = FALSE
      )
    }
    if (is.null(cov_cohort_X_s)) {
      warning(
        "No study-level covariance available, so covariances within each ",
        "study are estimated from the reference population.",
        call. = FALSE
      )
    }
    cov_names_ls <- lapply(reg_out, function(r) {
      drop_intercept_names(names(r$theta_hat_s))
    })
    cov_cohort_X_s <- lapply(seq_len(n_cohorts), function(k) {
      subset_by_name(
        cov_ref_X[[k]], cov_names_ls[[k]],
        sprintf("reference covariance for cohort %d", k)
      )
    })
  }

  ## ---------------------------------------------------------------------
  ## 3. Assemble
  ## ---------------------------------------------------------------------
  sigma_X_ls <- lapply(seq_len(n_cohorts), function(k) {
    sig_ref <- cov_ref_X[[k]]
    vars_cohort <- intersect(var_names_ref_pop, colnames(cov_cohort_X_s[[k]]))

    unknown <- setdiff(colnames(cov_cohort_X_s[[k]]), var_names_ref_pop)
    if (length(unknown) > 0L) {
      stop(sprintf(
        paste0(
          "Cohort %d measures covariate(s) %s that do not appear in the ",
          "reference population. Add them to the reference data (or the ",
          "reference covariance matrix) and try again."
        ),
        k, paste(unknown, collapse = ", ")
      ), call. = FALSE)
    }

    if (cov_p_type == "reference") {
      return(sig_ref[var_names_ref_pop, var_names_ref_pop, drop = FALSE])
    }

    vars_ref_only <- setdiff(var_names_ref_pop, vars_cohort)

    s11 <- cov_cohort_X_s[[k]][vars_cohort, vars_cohort, drop = FALSE]
    if (length(vars_ref_only) == 0L) {
      out <- s11
    } else {
      s12 <- sig_ref[vars_cohort, vars_ref_only, drop = FALSE]
      s22 <- sig_ref[vars_ref_only, vars_ref_only, drop = FALSE]
      out <- rbind(
        cbind(s11, s12),
        cbind(t(s12), s22)
      )
    }
    out[var_names_ref_pop, var_names_ref_pop, drop = FALSE]
  })

  list(
    sigma_X_ls = sigma_X_ls,
    cov_cohort_X_s_ls = cov_cohort_X_s,
    ref_weights = ref_weights,
    ref_pop_dummy = ref_pop_dummy,
    ref_strata = ref_strata,
    ref_psu = ref_psu,
    design_weights = design_weights,
    weighting_converged = weighting_converged
  )
}
