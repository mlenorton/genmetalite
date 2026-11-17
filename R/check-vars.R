#' Check that one dataset's variables and factor levels appear in another
#'
#' Compares a study (or reference) dataset against a template dataset and
#' reports variables present in `.df2` but not `.df1`, and factor levels present
#' in `.df2` but not `.df1`. Levels missing from the reference population are a
#' common cause of downstream singular-covariance failures, because the
#' corresponding dummy column is constant.
#'
#' @param .df1 The template data frame (typically the reference population).
#' @param .df2 The data frame being checked (typically a study cohort).
#' @param action `"warn"` (default) to emit warnings, `"stop"` to raise an
#'   error at the first problem, or `"message"` to report quietly.
#'
#' @return Invisibly, a character vector of the problems found (empty if none).
#'   Called mainly for its side effect.
#'
#' @examples
#' ref <- data.frame(g = factor(c("a", "b")))
#' coh <- data.frame(g = factor(c("a", "c")))
#' check_vars(ref, coh)
#'
#' @export
check_vars <- function(.df1, .df2, action = c("warn", "stop", "message")) {
  action <- match.arg(action)
  report <- switch(action,
    warn = function(msg) warning(msg, call. = FALSE),
    stop = function(msg) stop(msg, call. = FALSE),
    message = function(msg) message(msg)
  )

  problems <- character(0)

  extra_vars <- setdiff(names(.df2), names(.df1))
  if (length(extra_vars) > 0L) {
    msg <- paste0(
      "These covariates are in the second dataset but not the first: ",
      paste(sort(extra_vars), collapse = ", "), "."
    )
    problems <- c(problems, msg)
    report(msg)
  }

  ## inherits() rather than class() == "factor": an ordered factor has class
  ## c("ordered", "factor"), which the equality test silently misses.
  factor_cols <- names(.df2)[vapply(.df2, is.factor, logical(1))]

  for (nm in intersect(factor_cols, names(.df1))) {
    missing_levels <- setdiff(
      as.character(unique(.df2[[nm]])),
      as.character(unique(.df1[[nm]]))
    )
    missing_levels <- missing_levels[!is.na(missing_levels)]
    if (length(missing_levels) > 0L) {
      msg <- sprintf(
        "%s: level%s %s present in the second dataset but not the first.",
        nm,
        if (length(missing_levels) == 1L) "" else "s",
        paste(sort(missing_levels), collapse = ", ")
      )
      problems <- c(problems, msg)
      report(msg)
    }
  }

  invisible(problems)
}

#' @rdname check_vars
#' @export
check_vars_stop <- function(.df1, .df2) {
  check_vars(.df1, .df2, action = "stop")
}

#' Find factor levels a cohort measured but has no observations in
#'
#' A cohort can measure a categorical covariate and still have an empty cell in
#' one of its levels — a small study with nobody in the highest education
#' category, say. The dummy column for that level is constant at zero, so it
#' cannot be kept, and the study ends up contributing no information about it.
#'
#' @param .ref The reference data frame.
#' @param .cohort The cohort data frame.
#' @param action `"warn"` (default), `"stop"`, `"message"` or `"none"`.
#' @param label Optional label for the cohort, used in messages.
#'
#' @return Invisibly, a data frame with columns `variable`, `level` and
#'   `status`. `status` is `"unobserved"` for an ordinary empty level and
#'   `"baseline_unobserved"` for the special case described below. Zero rows
#'   means nothing was found.
#'
#' @section What happens to an empty level:
#' The dummy column is dropped, and [linear_estimate()] then treats that level
#' as a covariate the study never measured, filling its covariance blocks from
#' the reference population. That is a real assumption: it says the study's
#' (nonexistent) members of that level would have behaved like the reference's.
#' It is usually the least-bad option, but it should never be silent, which is
#' why this check exists and why [calc_linear_beta()] records its result in the
#' `empty_levels` component of the returned object.
#'
#' @section The baseline case is different:
#' If the *reference's baseline* level is the one with no observations, there is
#' nothing valid to continue with, so this is an error rather than a warning.
#' Dropping it makes [base::droplevels()] promote a different level to baseline,
#' and the cohort's `edu3` column then means "3 versus 2" while the reference's
#' `edu3` means "3 versus 1". Since the pipeline matches covariance blocks by
#' name, those two different contrasts would be silently treated as the same
#' covariate. Fix it by releveling the factor so a level both datasets observe
#' is the baseline, or by collapsing the sparse levels.
#'
#' @examples
#' ref <- data.frame(edu = factor(c("1", "2", "3")))
#' coh <- data.frame(edu = factor(c("1", "2"), levels = c("1", "2", "3")))
#' check_empty_levels(ref, coh)
#'
#' @seealso [check_vars()]
#' @export
check_empty_levels <- function(.ref,
                               .cohort,
                               action = c("warn", "stop", "message", "none"),
                               label = NULL) {
  action <- match.arg(action)
  where <- if (is.null(label)) "This cohort" else sprintf("Cohort %s", label)

  found <- data.frame(
    variable = character(0), level = character(0), status = character(0),
    stringsAsFactors = FALSE
  )

  categorical <- function(x) is.factor(x) || is.character(x)
  shared <- intersect(names(.ref), names(.cohort))
  shared <- shared[
    vapply(.ref[shared], categorical, logical(1)) &
      vapply(.cohort[shared], categorical, logical(1))
  ]

  for (nm in shared) {
    ref_obs <- unique(stats::na.omit(as.character(.ref[[nm]])))
    coh_obs <- unique(stats::na.omit(as.character(.cohort[[nm]])))
    missing_levels <- setdiff(ref_obs, coh_obs)
    if (length(missing_levels) == 0L) next

    ## The baseline is the first level of the reference factor, or the first
    ## value in sort order for a character column, which is what
    ## model.matrix() would use.
    baseline <- if (is.factor(.ref[[nm]])) {
      levels(droplevels(.ref[[nm]]))[1L]
    } else {
      sort(ref_obs)[1L]
    }

    status <- ifelse(
      missing_levels == baseline, "baseline_unobserved", "unobserved"
    )
    found <- rbind(found, data.frame(
      variable = nm, level = missing_levels, status = status,
      stringsAsFactors = FALSE
    ))
  }

  if (nrow(found) == 0L || action == "none") {
    return(invisible(found))
  }

  ## The baseline case has no valid continuation, so it is always an error.
  bad <- found[found$status == "baseline_unobserved", , drop = FALSE]
  if (nrow(bad) > 0L) {
    stop(sprintf(
      paste0(
        "%s has no observations in the reference baseline level of %s ",
        "(level '%s'). Dropping it would promote a different level to ",
        "baseline, so this cohort's dummy columns would encode different ",
        "contrasts than the reference's while carrying the same names. ",
        "Relevel %s so that a level both datasets observe is first, or ",
        "collapse the sparse levels."
      ),
      where,
      paste(unique(bad$variable), collapse = ", "),
      paste(bad$level, collapse = "', '"),
      paste(unique(bad$variable), collapse = " and ")
    ), call. = FALSE)
  }

  msg <- sprintf(
    paste0(
      "%s has no observations in %s. Th%s dummy column%s cannot be kept, so ",
      "%s treated as covariate(s) this cohort did not measure and %s ",
      "covariance blocks are filled from the reference population. Check that ",
      "this is what you want."
    ),
    where,
    paste(sprintf("%s level '%s'", found$variable, found$level),
      collapse = ", "
    ),
    if (nrow(found) == 1L) "at" else "ose",
    if (nrow(found) == 1L) "" else "s",
    if (nrow(found) == 1L) "it is" else "they are",
    if (nrow(found) == 1L) "its" else "their"
  )

  switch(action,
    warn = warning(msg, call. = FALSE),
    stop = stop(msg, call. = FALSE),
    message = message(msg)
  )

  invisible(found)
}

#' Give a cohort's shared factors the reference's level set and ordering.
#'
#' Without this, a cohort that happens not to observe the reference's first
#' level gets a different baseline from `model.matrix()`, and its dummy columns
#' silently encode different contrasts under the same names.
#' @noRd
harmonise_levels <- function(cohort_df, ref_df) {
  shared <- intersect(names(cohort_df), names(ref_df))
  for (nm in shared) {
    rv <- ref_df[[nm]]
    cv <- cohort_df[[nm]]
    if (!(is.factor(rv) || is.character(rv))) next
    if (!(is.factor(cv) || is.character(cv))) next

    ref_levels <- if (is.factor(rv)) {
      levels(droplevels(rv))
    } else {
      sort(unique(stats::na.omit(rv)))
    }
    coh_values <- unique(stats::na.omit(as.character(cv)))
    ## Reference levels first, so the baseline is the reference's baseline;
    ## anything only the cohort has is appended rather than dropped.
    all_levels <- c(ref_levels, setdiff(coh_values, ref_levels))
    cohort_df[[nm]] <- factor(as.character(cv), levels = all_levels)
  }
  cohort_df
}
