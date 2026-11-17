#' Reweight a reference population to a cohort by entropy balancing
#'
#' Computes weights for the rows of `ref_data` so that the weighted reference
#' population matches `cohort_data` on the requested covariate moments, using
#' the entropy balancing algorithm of Hainmueller (2012).
#'
#' @param cohort_data Data frame of individual-level cohort covariates (the
#'   target the reference is matched *to*).
#' @param ref_data Data frame of individual-level reference covariates (the rows
#'   that get weights).
#' @param weighting.vars Either `NULL` (match on every variable shared by the two
#'   datasets), a character vector of variable names, or a list with one element
#'   per moment: `weighting.vars[[m]]` names the variables whose `m`-th moment is
#'   to be balanced. A `NULL` element means "all shared variables".
#' @param weighting.max.moments Highest moment to balance. Factor variables are
#'   dropped for moments above the first, since powers of a 0/1 dummy are the
#'   dummy itself and would be perfectly collinear.
#' @param weighting.corr If `TRUE`, also balance pairwise products
#'   (i.e. covariances) among `weighting.corr.vars` at the first moment.
#' @param weighting.corr.vars Variables to form pairwise products from. Defaults
#'   to the first-moment weighting variables.
#' @param max.iterations Maximum iterations for the first attempt. If the
#'   algorithm does not converge, one further attempt is made with
#'   `retry.iterations`.
#' @param retry.iterations Iteration cap for the retry.
#' @param base.weight Optional base weights for the reference rows. **Survey
#'   weights go here.** The balancing solution then tilts the population the
#'   reference represents rather than the raw sample, and the returned weights
#'   are a single vector carrying both the sampling design and the covariate
#'   shift. Defaults to `NULL`, meaning equal weights.
#' @param engine `"internal"` (default) uses the solver in [entropy_balance()];
#'   `"ebal"` calls [ebal::ebalance()] instead, which requires the \pkg{ebal}
#'   package. Both solve the same strictly convex problem and agree to
#'   tolerance; the internal solver is used by default because it handles
#'   `base.weight` transparently and needs no extra dependency.
#' @param tol Convergence tolerance on the largest standardised imbalance
#'   (`engine = "internal"` only).
#' @param trim If `TRUE`, post-process with [ebal::ebalance.trim()] to reduce
#'   the dispersion of the weights. Requires `engine = "ebal"`.
#' @param plot If `TRUE`, also return a \pkg{cobalt} love plot of covariate
#'   balance before and after weighting. Requires the \pkg{cobalt} package.
#' @param ... Ignored, with a warning, so that typos in argument names are
#'   visible rather than silent.
#'
#' @return A list with
#'   \item{weights}{numeric vector of length `nrow(ref_data)`, scaled to sum to
#'     `sum(base.weight)` when base weights are supplied, so that a survey
#'     weight's population-total meaning survives.}
#'   \item{converged}{logical.}
#'   \item{max_imbalance}{largest standardised difference between the achieved
#'     and target moments (`engine = "internal"` only).}
#'   \item{plot}{a `ggplot` object, or `NULL` if `plot = FALSE`.}
#'
#' @section Missing values:
#' Rows with missing values in any weighting variable are an error rather than a
#' silent deletion. [stats::model.matrix()] would otherwise drop them, returning
#' a design matrix shorter than the treatment vector and producing weights that
#' do not line up with the reference rows.
#'
#' @references
#' Hainmueller, J. (2012). Entropy Balancing for Causal Effects.
#' *Political Analysis*, 20(1), 25-46.
#'
#' @seealso [ipw_weights()], [entropy_balance()], [ref_design()]
#' @export
ebal_weights <- function(cohort_data,
                         ref_data,
                         weighting.vars = NULL,
                         weighting.max.moments = 1L,
                         weighting.corr = FALSE,
                         weighting.corr.vars = NULL,
                         max.iterations = 200L,
                         retry.iterations = 3000L,
                         base.weight = NULL,
                         trim = FALSE,
                         plot = FALSE,
                         engine = c("internal", "ebal"),
                         tol = 1e-8,
                         ...) {
  engine <- match.arg(engine)
  if (engine == "ebal" && !requireNamespace("ebal", quietly = TRUE)) {
    stop(
      "engine = 'ebal' requires the 'ebal' package. Install it with ",
      "install.packages('ebal'), or use engine = 'internal' (the default).",
      call. = FALSE
    )
  }
  if (isTRUE(trim) && engine != "ebal") {
    stop(
      "trim = TRUE is implemented by ebal::ebalance.trim(); set ",
      "engine = 'ebal' to use it.",
      call. = FALSE
    )
  }
  dots <- list(...)
  if (length(dots) > 0L) {
    warning(
      "Ignoring unused argument(s): ", paste(names(dots), collapse = ", "),
      ". Note that make_sigma_X_ls() passes `weighting.max.iterations`; ",
      "ebal_weights() calls the same setting `max.iterations`.",
      call. = FALSE
    )
  }

  cohort_data <- as.data.frame(cohort_data)
  ref_data <- as.data.frame(ref_data)

  n_k <- nrow(cohort_data)
  n_ref <- nrow(ref_data)
  treatment <- c(rep(0L, n_ref), rep(1L, n_k))

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

  shared.vars <- intersect(names(cohort_data), names(ref_data))
  if (length(shared.vars) == 0L) {
    stop("`cohort_data` and `ref_data` share no variable names.", call. = FALSE)
  }

  weighting.vars <- normalise_weighting_vars(
    weighting.vars, weighting.max.moments
  )

  X_with_fct <- rbind(
    ref_data[, shared.vars, drop = FALSE],
    cohort_data[, shared.vars, drop = FALSE]
  )
  X_with_fct <- droplevels(X_with_fct)
  factor.vars <- names(X_with_fct)[
    vapply(X_with_fct, function(x) is.factor(x) || is.character(x), logical(1))
  ]

  X <- NULL

  for (m in seq_len(weighting.max.moments)) {
    vars_m <- weighting.vars[[m]] %||% shared.vars
    vars_m <- intersect(vars_m, shared.vars)
    if (length(vars_m) == 0L) {
      stop(sprintf(
        "No shared variables left to balance at moment %d.", m
      ), call. = FALSE)
    }

    X_m <- X_with_fct[, vars_m, drop = FALSE]

    if (m > 1L) {
      ## Powers of a 0/1 dummy equal the dummy, so factors must go.
      X_m <- X_m[, setdiff(names(X_m), factor.vars), drop = FALSE]
      if (ncol(X_m) == 0L) {
        stop(
          "Only categorical variables were selected, so moments above the ",
          "first cannot be balanced (a dummy raised to any power is itself). ",
          "Set weighting.max.moments = 1 or include a continuous variable.",
          call. = FALSE
        )
      }
    }

    if (anyNA(X_m)) {
      stop(sprintf(
        paste0(
          "Missing values in weighting variable(s) %s. Entropy balancing needs ",
          "complete rows; drop or impute them first."
        ),
        paste(names(X_m)[vapply(X_m, anyNA, logical(1))], collapse = ", ")
      ), call. = FALSE)
    }

    if (m == 1L) {
      form_rhs <- paste(names(X_m), collapse = " + ")
      if (isTRUE(weighting.corr)) {
        corr_vars <- intersect(weighting.corr.vars %||% names(X_m), names(X_m))
        if (length(corr_vars) >= 2L) {
          pairs <- utils::combn(corr_vars, 2L, simplify = FALSE)
          form_rhs <- paste(
            form_rhs, "+",
            paste(vapply(
              pairs, function(p) paste(p[1L], ":", p[2L]), character(1)
            ), collapse = " + ")
          )
        }
      }
      form <- stats::as.formula(paste("~", form_rhs))
      X <- dummy_design(form, X_m)
    } else {
      X_dummy <- dummy_design(
        stats::as.formula(paste("~", paste(names(X_m), collapse = " + "))), X_m
      )
      X_pow <- X_dummy^m
      colnames(X_pow) <- paste0(colnames(X_dummy), "_pow", m)
      X <- cbind(X, X_pow)
    }
  }

  ## Drop constant columns: they carry no balance information and make the
  ## entropy-balancing Newton step singular.
  keep <- apply(X, 2L, function(col) stats::sd(col) > 0)
  if (any(!keep)) {
    warning(
      "Dropping constant column(s) from the balancing design: ",
      paste(colnames(X)[!keep], collapse = ", "), ".",
      call. = FALSE
    )
    X <- X[, keep, drop = FALSE]
  }

  ## Fit -----------------------------------------------------------------
  ref_rows <- seq_len(n_ref)
  Z_ref <- X[ref_rows, , drop = FALSE]
  target <- colMeans(X[-ref_rows, , drop = FALSE])
  max_imbalance <- NA_real_

  if (engine == "internal") {
    eb <- entropy_balance(
      Z = Z_ref, target = target, base_weight = base.weight,
      max_iterations = max(max.iterations, retry.iterations), tol = tol
    )
    weights_out <- eb$weights
    converged <- eb$converged
    max_imbalance <- eb$max_imbalance
    if (!converged) {
      warning(sprintf(
        paste0(
          "Entropy balancing did not converge (largest standardised ",
          "imbalance %.3g). This usually means the cohort's moments lie ",
          "outside the range the reference population can be reweighted to. ",
          "Consider balancing fewer moments or dropping nearly collinear ",
          "covariates."
        ), max_imbalance
      ), call. = FALSE)
    }
  } else {
    eb.out <- ebal::ebalance(
      Treatment = treatment, X = X,
      max.iterations = max.iterations,
      base.weight = base.weight
    )
    if (!isTRUE(eb.out$converged)) {
      eb.out <- ebal::ebalance(
        Treatment = treatment, X = X,
        max.iterations = retry.iterations,
        base.weight = base.weight
      )
    }
    if (isTRUE(trim)) {
      eb.out <- ebal::ebalance.trim(eb.out)
    }
    weights_out <- as.numeric(eb.out$w)
    converged <- isTRUE(eb.out$converged)
    if (!converged) {
      warning(
        "Entropy balancing did not converge after ", retry.iterations,
        " iterations. Consider balancing fewer moments, dropping nearly ",
        "collinear covariates, or increasing `retry.iterations`.",
        call. = FALSE
      )
    }
  }

  p_love <- NULL
  if (isTRUE(plot)) {
    if (!requireNamespace("cobalt", quietly = TRUE)) {
      warning("plot = TRUE requires the 'cobalt' package; skipping the plot.",
        call. = FALSE
      )
    } else {
      bal_df <- as.data.frame(cbind(treatment = treatment, X))
      p_love <- cobalt::love.plot(
        cobalt::bal.tab(
          stats::as.formula(paste(
            "treatment ~", paste(sprintf("`%s`", colnames(X)), collapse = " + ")
          )),
          data = bal_df,
          weights = c(weights_out, rep(1, n_k))
        )
      )
    }
  }

  list(
    weights = weights_out,
    converged = converged,
    max_imbalance = max_imbalance,
    plot = p_love
  )
}

#' Build a dummy design matrix from a formula, dropping the intercept
#' @noRd
dummy_design <- function(form, data) {
  old <- options(na.action = "na.pass")
  on.exit(options(old), add = TRUE)
  mm <- stats::model.matrix(form, data = droplevels(data))
  mm[, colnames(mm) != "(Intercept)", drop = FALSE]
}

#' Normalise the several accepted shapes of `weighting.vars` into a list of
#' length `n_moments`.
#' @noRd
normalise_weighting_vars <- function(weighting.vars, n_moments) {
  n_moments <- as.integer(n_moments)
  if (is.na(n_moments) || n_moments < 1L) {
    stop("`weighting.max.moments` must be a positive integer.", call. = FALSE)
  }
  if (is.null(weighting.vars) || length(weighting.vars) == 0L) {
    return(vector("list", n_moments))
  }
  if (is.character(weighting.vars)) {
    ## The same variables at every moment.
    return(rep(list(weighting.vars), n_moments))
  }
  if (!is.list(weighting.vars)) {
    stop(
      "`weighting.vars` must be NULL, a character vector, or a list of ",
      "character vectors (one per moment).",
      call. = FALSE
    )
  }
  out <- vector("list", n_moments)
  for (m in seq_len(n_moments)) {
    if (m <= length(weighting.vars)) out[[m]] <- weighting.vars[[m]]
  }
  out
}
