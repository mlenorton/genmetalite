#' Describe a reference population and how it was sampled
#'
#' Bundles the reference data with the sampling weights and design information
#' that go with them, so the design travels with the population instead of being
#' passed around as loose vectors. [calc_linear_beta()] and
#' [make_sigma_X_ls()] accept one of these wherever they accept `ref_data`.
#'
#' @param data A data frame of reference covariates. No outcome is needed.
#' @param weights Optional sampling weights, one per row.
#' @param strata Optional stratum identifiers, one per row.
#' @param psu Optional primary sampling unit identifiers, one per row.
#'
#' @return An object of class `genmetalite_ref_design`.
#'
#' @details
#' Three kinds of reference input are accepted throughout the package:
#' \describe{
#'   \item{a data frame}{an unweighted, independently sampled reference.}
#'   \item{a `survey.design` object}{from `survey::svydesign()`; weights,
#'     strata and clusters are read off it. See [as_ref_design()].}
#'   \item{a `ref_design` object}{this constructor, for when the design is known
#'     but the \pkg{survey} package is not being used.}
#' }
#'
#' Supplying weights here affects both halves of the estimate: the covariance
#' \eqn{\Sigma} becomes a weighted covariance, and the reference term of the
#' variance uses the design rather than the independent-sampling formula.
#'
#' @examples
#' set.seed(1)
#' d <- data.frame(x1 = rnorm(200), x2 = rnorm(200))
#' rd <- ref_design(
#'   d,
#'   weights = runif(200, 1, 5),
#'   strata = rep(1:5, each = 40),
#'   psu = rep(1:20, each = 10)
#' )
#' rd
#'
#' @seealso [as_ref_design()], [calc_linear_beta()]
#' @export
ref_design <- function(data, weights = NULL, strata = NULL, psu = NULL) {
  data <- as.data.frame(data)
  n <- nrow(data)
  if (n == 0L) {
    stop("`data` has no rows.", call. = FALSE)
  }

  check_len <- function(x, nm) {
    if (is.null(x)) {
      return(NULL)
    }
    if (length(x) != n) {
      stop(sprintf(
        "`%s` has length %d but `data` has %d rows.", nm, length(x), n
      ), call. = FALSE)
    }
    x
  }

  weights <- check_len(weights, "weights")
  if (!is.null(weights)) {
    weights <- as.numeric(weights)
    if (any(!is.finite(weights)) || any(weights < 0)) {
      stop("`weights` must be finite and non-negative.", call. = FALSE)
    }
    if (sum(weights) <= 0) {
      stop("`weights` sum to zero.", call. = FALSE)
    }
  }

  structure(
    list(
      data = data,
      weights = weights,
      strata = check_len(strata, "strata"),
      psu = check_len(psu, "psu")
    ),
    class = "genmetalite_ref_design"
  )
}

#' @export
print.genmetalite_ref_design <- function(x, ...) {
  cat(sprintf(
    "<reference design>  %d rows, %d variable(s)\n",
    nrow(x$data), ncol(x$data)
  ))
  cat(sprintf(
    "  weights: %s\n",
    if (is.null(x$weights)) {
      "none (equal)"
    } else {
      sprintf(
        "supplied, sum = %.4g, Kish n_eff = %.0f",
        sum(x$weights), sum(x$weights)^2 / sum(x$weights^2)
      )
    }
  ))
  cat(sprintf(
    "  design : %s\n",
    if (is.null(x$strata) && is.null(x$psu)) {
      "independent sampling"
    } else {
      sprintf(
        "%d stratum(a), %d PSU(s)",
        length(unique(x$strata %||% 1L)),
        nrow(unique(data.frame(
          s = x$strata %||% 1L, p = x$psu %||% seq_len(nrow(x$data))
        )))
      )
    }
  ))
  invisible(x)
}

#' Convert a reference population into a reference design
#'
#' Accepts a data frame, a `survey.design` object from the \pkg{survey} package,
#' or an object already made by [ref_design()], and returns a
#' `genmetalite_ref_design`.
#'
#' @param x A data frame, a `survey.design`/`survey.design2` object, or a
#'   `genmetalite_ref_design`.
#' @param ... Unused.
#'
#' @return A `genmetalite_ref_design` object.
#'
#' @section Reading a survey design:
#' Weights come from `stats::weights()` (falling back to `1/design$prob`), and
#' strata and clusters from the *first stage* of `design$strata` and
#' `design$cluster`. That is the "ultimate cluster" convention: for a multistage
#' design, only the first-stage units contribute to the variance, which is what
#' NCHS recommends for NHANES and what `survey::svydesign()` itself computes by
#' default.
#'
#' Two things are deliberately ignored, both in the conservative direction.
#' Finite population corrections are dropped, so the variance is the
#' with-replacement one. Post-stratification and calibration adjustments are not
#' propagated, so a calibrated design gets a warning and is treated as if only
#' its weights mattered; the reported variance is then somewhat too large.
#'
#' Replicate-weight designs (`svrepdesign()`) are rejected, because their
#' variance is computed by re-estimating over replicates rather than by
#' linearization, which is a different code path than this package implements.
#'
#' @examples
#' if (requireNamespace("survey", quietly = TRUE)) {
#'   data(api, package = "survey")
#'   des <- survey::svydesign(
#'     id = ~dnum, strata = ~stype, weights = ~pw,
#'     data = apistrat, nest = TRUE
#'   )
#'   as_ref_design(des)
#' }
#'
#' @seealso [ref_design()]
#' @export
as_ref_design <- function(x, ...) {
  if (inherits(x, "genmetalite_ref_design")) {
    return(x)
  }
  if (is.data.frame(x)) {
    return(ref_design(x))
  }
  if (inherits(x, c("svyrep.design", "svrepdesign"))) {
    stop(
      "Replicate-weight survey designs are not supported. Their variance is ",
      "computed by re-estimating over replicate weights rather than by ",
      "linearization. Either rebuild the design with survey::svydesign() ",
      "using strata and cluster identifiers, or supply the pieces directly ",
      "with ref_design(data, weights, strata, psu).",
      call. = FALSE
    )
  }
  if (inherits(x, c("survey.design", "survey.design2"))) {
    return(ref_design_from_survey(x))
  }
  stop(
    "`ref_data` must be a data frame, a list of data frames, a ",
    "survey.design object, or a ref_design() object; got an object of class ",
    paste(class(x), collapse = "/"), ".",
    call. = FALSE
  )
}

#' Pull data, weights, strata and first-stage clusters off a survey design.
#' @noRd
ref_design_from_survey <- function(des) {
  if (!requireNamespace("survey", quietly = TRUE)) {
    stop(
      "Reading a survey.design object requires the 'survey' package. ",
      "Install it with install.packages('survey'), or pass the pieces ",
      "directly with ref_design(data, weights, strata, psu).",
      call. = FALSE
    )
  }

  vars <- des$variables
  if (is.null(vars)) {
    stop("The survey design has no `variables` component.", call. = FALSE)
  }
  vars <- as.data.frame(vars)

  w <- tryCatch(as.numeric(stats::weights(des)), error = function(e) NULL)
  if (is.null(w) || length(w) != nrow(vars)) {
    if (is.null(des$prob)) {
      stop(
        "Could not read sampling weights off the survey design. Pass them ",
        "explicitly with ref_design(data, weights, ...).",
        call. = FALSE
      )
    }
    w <- 1 / as.numeric(des$prob)
  }

  ## First stage only (ultimate cluster).
  first_stage <- function(x) {
    if (is.null(x)) {
      return(NULL)
    }
    if (is.data.frame(x) || is.matrix(x)) {
      if (ncol(x) == 0L) {
        return(NULL)
      }
      return(as.character(x[[1L]]))
    }
    as.character(x)
  }

  st <- first_stage(des$strata)
  ps <- first_stage(des$cluster)

  if (isFALSE(des$has.strata)) st <- NULL
  ## svydesign() fills `cluster` with row numbers when id = ~0 or ~1, which is
  ## the no-clustering case and carries no information.
  if (!is.null(ps) && length(unique(ps)) == nrow(vars)) ps <- NULL

  if (!is.null(des$postStrata)) {
    warning(
      "The survey design is post-stratified or calibrated. Those adjustments ",
      "reduce variance, and this package does not propagate them, so the ",
      "reported reference variance will be conservative.",
      call. = FALSE
    )
  }
  if (!is.null(des$fpc) && !is.null(des$fpc$popsize)) {
    message(
      "Finite population corrections are ignored; the with-replacement ",
      "('ultimate cluster') variance is used, as recommended for NHANES."
    )
  }

  ref_design(vars, weights = w, strata = st, psu = ps)
}

#' Normalise any accepted reference input into the internal representation.
#'
#' Returns a list with per-cohort `data` and `weights`, plus `strata`/`psu`
#' shared across cohorts (a single reference sample can only have one design).
#' `shared` records whether every cohort sees the same reference rows, which
#' determines whether their influence contributions may be added together.
#' @noRd
normalise_reference <- function(ref_data, n_cohorts) {
  if (is.null(ref_data)) {
    return(NULL)
  }

  ## A plain list of data frames: genuinely separate reference samples, or the
  ## same one repeated. No design in that case.
  if (is.list(ref_data) && !is.data.frame(ref_data) &&
    !inherits(ref_data, "genmetalite_ref_design") &&
    !inherits(ref_data, c("survey.design", "survey.design2")) &&
    all(vapply(ref_data, is.data.frame, logical(1)))) {
    dfs <- as_cohort_list(ref_data, n_cohorts, "ref_data")
    same <- all(vapply(dfs, function(d) identical(d, dfs[[1L]]), logical(1)))
    return(list(
      data = lapply(dfs, function(d) droplevels(as.data.frame(d))),
      weights = NULL,
      strata = NULL,
      psu = NULL,
      shared = same
    ))
  }

  des <- as_ref_design(ref_data)
  d <- droplevels(des$data)

  list(
    data = rep(list(d), n_cohorts),
    weights = if (is.null(des$weights)) {
      NULL
    } else {
      rep(list(des$weights), n_cohorts)
    },
    strata = des$strata,
    psu = des$psu,
    shared = TRUE
  )
}
