## Internal helpers -----------------------------------------------------------
## Not exported. Kept small and dependency-free on purpose.

`%||%` <- function(x, y) if (is.null(x)) y else x

#' Is a name an intercept name?
#' @noRd
is_intercept_name <- function(nm) {
  grepl("^intercept$|^\\(intercept\\)$", tolower(nm))
}

#' Drop intercept-like entries from a character vector
#' @noRd
drop_intercept_names <- function(nm) nm[!is_intercept_name(nm)]

#' Does a data frame contain any categorical column?
#' @noRd
has_categorical <- function(df) {
  any(vapply(df, function(x) is.character(x) || is.factor(x), logical(1)))
}

#' Expand a data frame to a numeric dummy-coded matrix.
#'
#' Uses [stats::model.matrix()] with treatment contrasts (first level dropped)
#' and *always* removes the `(Intercept)` column. Unused factor levels are
#' dropped first, since they produce all-zero columns and hence singular
#' covariance matrices.
#'
#' `na.action = "na.pass"` is set deliberately: the default
#' (`getOption("na.action")`, normally `na.omit`) silently deletes rows with any
#' missing value, which makes the returned matrix shorter than the input and
#' misaligns it with the outcome vector or a weight vector.
#'
#' @param df A data frame.
#' @param add_intercept Logical; if `TRUE`, prepend a column of ones named
#'   `Intercept`.
#' @return A numeric matrix with `nrow(df)` rows.
#' @noRd
dummy_matrix <- function(df, add_intercept = FALSE) {
  df <- as.data.frame(df)
  if (ncol(df) == 0L) {
    stop("`df` has no columns to expand.", call. = FALSE)
  }
  df <- droplevels(df)
  old <- options(na.action = "na.pass")
  on.exit(options(old), add = TRUE)
  mm <- stats::model.matrix(~., data = df)
  mm <- mm[, colnames(mm) != "(Intercept)", drop = FALSE]
  if (add_intercept) {
    mm <- cbind(Intercept = 1, mm)
  }
  mm
}

#' Symmetrise a matrix (guards against accumulated floating-point asymmetry)
#' @noRd
symmetrise <- function(m) (m + t(m)) / 2

#' Solve with an informative error message
#' @noRd
safe_solve <- function(m, what = "matrix", b = NULL) {
  out <- tryCatch(
    if (is.null(b)) solve(m) else solve(m, b),
    error = function(e) e
  )
  if (inherits(out, "error")) {
    stop(
      sprintf(
        paste0(
          "Could not invert the %s (%d x %d). This usually means two ",
          "covariates are collinear, or a factor level is empty in one of ",
          "the datasets. Original error: %s"
        ),
        what, nrow(m), ncol(m), conditionMessage(out)
      ),
      call. = FALSE
    )
  }
  out
}

#' Coerce a "maybe list, maybe data frame" reference argument into a list of
#' length `n`, recycling a single data frame.
#' @noRd
as_cohort_list <- function(x, n, arg = deparse(substitute(x))) {
  if (is.null(x)) {
    return(NULL)
  }
  if (is.data.frame(x)) {
    return(rep(list(x), n))
  }
  if (!is.list(x)) {
    stop(sprintf("`%s` must be a data frame or a list of data frames.", arg),
      call. = FALSE
    )
  }
  if (length(x) == 1L && n > 1L) {
    return(rep(x, n))
  }
  if (length(x) != n) {
    stop(sprintf(
      "`%s` has length %d but %d cohorts were supplied.",
      arg, length(x), n
    ), call. = FALSE)
  }
  x
}

#' All-NULL list test
#' @noRd
all_null <- function(x) {
  is.null(x) || (is.list(x) && all(vapply(x, is.null, logical(1))))
}

#' Reorder/subset a square named matrix by name, with a clear error.
#' @noRd
subset_by_name <- function(m, nms, what = "matrix") {
  missing <- setdiff(nms, colnames(m))
  if (length(missing) > 0L) {
    stop(sprintf(
      "The %s is missing these variables: %s.",
      what, paste(missing, collapse = ", ")
    ), call. = FALSE)
  }
  m[nms, nms, drop = FALSE]
}
