#' Complex-survey variance of an influence-function matrix
#'
#' Taylor-linearization ("ultimate cluster", with-replacement) variance for a
#' statistic whose per-observation influence contributions are the rows of `IF`,
#' under a stratified, clustered design with sampling weights. This is the
#' estimator NCHS recommends for NHANES, and the one
#' `survey::svydesign(id = ~PSU, strata = ~STRATA, weights = ~W)` computes.
#'
#' @param IF Numeric matrix, one row per reference observation and one column
#'   per parameter. Centering is handled internally, so raw influence
#'   contributions are fine.
#' @param weights Numeric vector of sampling weights, one per row of `IF`.
#'   `NULL` means equal weights.
#' @param strata Stratum identifiers, one per row of `IF`. `NULL` means a single
#'   stratum.
#' @param psu Primary sampling unit identifiers, one per row of `IF`. `NULL`
#'   means each row is its own PSU (i.e. no clustering).
#' @param lonely_psu How to handle a stratum containing exactly one PSU:
#'   `"remove"` (the default; the stratum contributes nothing, with a warning),
#'   `"adjust"` (centre its PSU total at the overall mean PSU total rather than
#'   the undefined stratum mean, which is conservative), or `"fail"`.
#'
#' @return A list with
#'   \item{cov}{the variance-covariance matrix.}
#'   \item{df}{design degrees of freedom, number of PSUs minus number of
#'     strata. Use a \eqn{t} rather than a normal reference distribution when
#'     this is small, which is common for survey subgroups.}
#'   \item{n_psu,n_strata}{counts used for `df`.}
#'   \item{n_eff}{Kish effective sample size \eqn{(\sum w)^2 / \sum w^2},
#'     reflecting the weighting only, not the clustering.}
#'
#' @details
#' With linearized per-observation contributions to the weighted mean
#' \deqn{q_j = w_j (\mathrm{IF}_j - \overline{\mathrm{IF}}_w) / \textstyle\sum_j w_j,}
#' PSU totals \eqn{t_{hi} = \sum_{j \in hi} q_j} and stratum means
#' \eqn{\bar t_h}, the variance is
#' \deqn{\widehat V = \sum_h \frac{n_h}{n_h - 1}
#'       \sum_{i \in h} (t_{hi} - \bar t_h)(t_{hi} - \bar t_h)'.}
#'
#' With no `strata` or `psu` this reduces to the independent-sampling form
#' \eqn{\mathrm{crossprod}(w \cdot \mathrm{IF}_c) / (\sum_j w_j)^2}, times the
#' usual \eqn{n/(n-1)} factor. Note that the iid branch inside
#' [linear_estimate()] omits that finite-sample factor, so the two paths differ
#' by \eqn{n/(n-1)} — negligible at any realistic reference sample size.
#'
#' Survey weights are population totals, so \eqn{\sum w} is not a sample size
#' and the variance must never be divided by it: the variance of a weighted mean
#' scales like \eqn{\sum w^2 / (\sum w)^2}, which is why the squared-weight form
#' above is the right one.
#'
#' @examples
#' set.seed(11)
#' IF <- matrix(rnorm(200), ncol = 2)
#' out <- survey_reference_cov(
#'   IF,
#'   weights = rep(1, 100),
#'   strata = rep(1:5, each = 20),
#'   psu = rep(1:20, each = 5)
#' )
#' out$df
#'
#' @seealso [linear_estimate()]
#' @export
survey_reference_cov <- function(IF,
                                 weights = NULL,
                                 strata = NULL,
                                 psu = NULL,
                                 lonely_psu = c("remove", "adjust", "fail")) {
  lonely_psu <- match.arg(lonely_psu)
  IF <- as.matrix(IF)
  n <- nrow(IF)
  p <- ncol(IF)

  weights <- if (is.null(weights)) rep(1, n) else as.numeric(weights)
  strata <- strata %||% rep(1L, n)
  psu <- psu %||% seq_len(n)

  for (nm in c("weights", "strata", "psu")) {
    if (length(get(nm)) != n) {
      stop(sprintf(
        "`%s` has length %d but `IF` has %d rows.", nm, length(get(nm)), n
      ), call. = FALSE)
    }
  }
  if (any(!is.finite(weights)) || any(weights < 0)) {
    stop("`weights` must be finite and non-negative.", call. = FALSE)
  }
  total_w <- sum(weights)
  if (total_w <= 0) {
    stop("`weights` sum to zero.", call. = FALSE)
  }

  strata <- as.character(strata)
  psu <- as.character(psu)

  ## Linearized per-observation contributions to the weighted mean. Centering
  ## happens here rather than being assumed of the caller: a constant shift in
  ## IF adds a PSU-size-dependent amount to each PSU total, which subtracting
  ## the stratum mean does not remove unless every PSU carries the same total
  ## weight.
  IFbar <- colSums(weights * IF) / total_w
  q <- (weights * sweep(IF, 2L, IFbar, "-")) / total_w

  out <- psu_linearized_cov(q, strata, psu, lonely_psu)
  out$n_eff <- total_w^2 / sum(weights^2)
  out
}

#' Ultimate-cluster variance of already-linearized contributions.
#'
#' `q` holds one row per observation, already weight-scaled and centred, so that
#' the statistic of interest is `colSums(q)`. Keeping this separate from
#' [survey_reference_cov()] lets [linear_estimate()] combine contributions that
#' carry *different* weight vectors per cohort — which is exactly what happens
#' when one survey reference is entropy-balanced separately to each cohort —
#' before applying the design.
#' @noRd
psu_linearized_cov <- function(q,
                               strata = NULL,
                               psu = NULL,
                               lonely_psu = "remove") {
  q <- as.matrix(q)
  n <- nrow(q)
  p <- ncol(q)
  strata <- as.character(strata %||% rep(1L, n))
  psu <- as.character(psu %||% seq_len(n))

  ## PSU totals per stratum. Building the matrix with an explicit `nrow` keeps
  ## the result n_psu x p even when p == 1; `t(vapply(..., numeric(p)))`
  ## silently transposes to 1 x n_psu in that case, which then fails to conform.
  psu_totals <- function(idx) {
    ids <- unique(psu[idx])
    out <- matrix(0, nrow = length(ids), ncol = p)
    for (i in seq_along(ids)) {
      rows <- idx[psu[idx] == ids[i]]
      out[i, ] <- colSums(q[rows, , drop = FALSE])
    }
    out
  }

  strata_ids <- unique(strata)
  by_stratum <- lapply(strata_ids, function(h) psu_totals(which(strata == h)))
  names(by_stratum) <- strata_ids

  n_psu <- sum(vapply(by_stratum, nrow, integer(1)))
  n_strata <- length(strata_ids)
  grand_mean <- colMeans(do.call(rbind, by_stratum))

  V <- matrix(0, p, p)
  lonely <- 0L

  for (h in strata_ids) {
    th <- by_stratum[[h]]
    n_h <- nrow(th)

    if (n_h == 1L) {
      lonely <- lonely + 1L
      if (lonely_psu == "fail") {
        stop(sprintf(
          paste0(
            "Stratum '%s' has a single PSU, so its variance cannot be ",
            "estimated. Set lonely_psu = 'remove' to let it contribute 0, ",
            "'adjust' to centre it at the overall mean, or collapse it with a ",
            "neighbouring stratum."
          ), h
        ), call. = FALSE)
      }
      if (lonely_psu == "adjust") {
        dev <- sweep(th, 2L, grand_mean, "-")
        V <- V + crossprod(dev)
      }
      next
    }

    dev <- sweep(th, 2L, colMeans(th), "-")
    V <- V + (n_h / (n_h - 1)) * crossprod(dev)
  }

  if (lonely > 0L && lonely_psu == "remove") {
    warning(sprintf(
      paste0(
        "%d single-PSU stratum(a) contributed 0 to the reference variance ",
        "(lonely_psu = 'remove'), which understates it. Consider collapsing ",
        "those strata, or lonely_psu = 'adjust'."
      ), lonely
    ), call. = FALSE)
  }

  list(
    cov = V,
    df = n_psu - n_strata,
    n_psu = n_psu,
    n_strata = n_strata
  )
}
