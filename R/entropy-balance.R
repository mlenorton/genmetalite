#' Entropy balancing weights for a set of balance functions
#'
#' Finds weights on the rows of `Z` that reproduce `target` exactly while
#' staying as close as possible, in Kullback-Leibler divergence, to the base
#' weights. This is the estimator of Hainmueller (2012), solved through its
#' convex dual.
#'
#' @param Z Numeric matrix of balance functions, one row per reference
#'   observation and one column per moment to match.
#' @param target Numeric vector of length `ncol(Z)`: the weighted means the
#'   reference should be tilted to reproduce, normally the cohort's means.
#' @param base_weight Optional non-negative base weights, one per row of `Z`.
#'   Sampling weights go here: the solution then tilts the *population* the
#'   reference represents, not the raw sample. Defaults to equal weights.
#' @param max_iterations Maximum Newton iterations.
#' @param tol Convergence tolerance on the largest standardised imbalance.
#'
#' @return A list with
#'   \item{weights}{the balancing weights, scaled to sum to `sum(base_weight)`
#'     so that a survey weight's population-total interpretation survives.}
#'   \item{converged}{logical.}
#'   \item{iterations}{Newton iterations used.}
#'   \item{max_imbalance}{largest absolute standardised difference between the
#'     achieved and target moments.}
#'
#' @details
#' Writing \eqn{\tilde Z_j = Z_j - \mathrm{target}} and \eqn{q_j} for the base
#' weights, the weights have the form
#' \eqn{w_j \propto q_j \exp(-\tilde Z_j' \lambda)}, and \eqn{\lambda} minimises
#' the smooth convex dual \eqn{L(\lambda) = \log \sum_j q_j
#' \exp(-\tilde Z_j' \lambda)}. Its gradient is the weighted imbalance and its
#' Hessian is the weighted covariance of \eqn{\tilde Z}, so Newton's method with
#' a backtracking line search converges quickly and the solution is unique
#' whenever the target lies inside the convex hull of the rows of `Z`.
#'
#' If the target lies outside that hull — a cohort mean more extreme than any
#' reference observation — no set of positive weights can reproduce it, the
#' dual has no finite minimiser, and the function reports non-convergence rather
#' than returning near-degenerate weights.
#'
#' @references
#' Hainmueller, J. (2012). Entropy Balancing for Causal Effects.
#' *Political Analysis*, 20(1), 25-46.
#'
#' @examples
#' set.seed(4)
#' Z <- cbind(a = rnorm(1000), b = rnorm(1000))
#' out <- entropy_balance(Z, target = c(a = 0.5, b = -0.2))
#' c(
#'   achieved_a = weighted.mean(Z[, "a"], out$weights),
#'   achieved_b = weighted.mean(Z[, "b"], out$weights)
#' )
#'
#' @seealso [ebal_weights()]
#' @export
entropy_balance <- function(Z,
                            target,
                            base_weight = NULL,
                            max_iterations = 200L,
                            tol = 1e-8) {
  Z <- as.matrix(Z)
  n <- nrow(Z)
  k <- ncol(Z)

  if (length(target) != k) {
    stop(sprintf(
      "`target` has length %d but `Z` has %d column(s).", length(target), k
    ), call. = FALSE)
  }
  if (anyNA(Z) || anyNA(target)) {
    stop("`Z` and `target` must not contain missing values.", call. = FALSE)
  }

  q <- if (is.null(base_weight)) rep(1, n) else as.numeric(base_weight)
  if (length(q) != n) {
    stop(sprintf(
      "`base_weight` has length %d but `Z` has %d row(s).", length(q), n
    ), call. = FALSE)
  }
  if (any(!is.finite(q)) || any(q < 0) || sum(q) <= 0) {
    stop("`base_weight` must be finite, non-negative and not all zero.",
      call. = FALSE
    )
  }
  q_total <- sum(q)
  q <- q / q_total

  Zt <- sweep(Z, 2L, as.numeric(target), "-")

  ## Scale used to judge imbalance, so `tol` means the same thing whatever the
  ## covariates are measured in.
  scale_z <- apply(Zt, 2L, function(col) stats::sd(col))
  scale_z[!is.finite(scale_z) | scale_z <= 0] <- 1

  weights_at <- function(lambda) {
    e <- -as.numeric(Zt %*% lambda)
    e <- e - max(e) # stabilise exp()
    w <- q * exp(e)
    w / sum(w)
  }
  objective <- function(lambda) {
    e <- -as.numeric(Zt %*% lambda)
    m <- max(e)
    log(sum(q * exp(e - m))) + m
  }

  lambda <- rep(0, k)
  obj <- objective(lambda)
  converged <- FALSE
  iter <- 0L
  w <- weights_at(lambda)

  for (iter in seq_len(max_iterations)) {
    w <- weights_at(lambda)
    grad <- -as.numeric(crossprod(Zt, w)) # -E_w[Zt]
    imbalance <- max(abs(grad) / scale_z)
    if (imbalance < tol) {
      converged <- TRUE
      break
    }

    mu <- as.numeric(crossprod(Zt, w))
    H <- crossprod(Zt * sqrt(w)) - tcrossprod(mu)

    step <- tryCatch(solve(H, grad), error = function(e) NULL)
    if (is.null(step) || any(!is.finite(step))) {
      ## Ridge-regularised fallback for a singular Hessian, which happens when
      ## two balance functions are collinear.
      ridge <- 1e-8 * max(diag(H), 1)
      step <- tryCatch(solve(H + diag(ridge, k), grad),
        error = function(e) NULL
      )
      if (is.null(step)) {
        break
      }
    }

    ## Backtracking line search on the convex dual.
    alpha <- 1
    ok <- FALSE
    for (ls in 1:40) {
      cand <- lambda - alpha * step
      obj_new <- objective(cand)
      if (is.finite(obj_new) && obj_new <= obj) {
        lambda <- cand
        obj <- obj_new
        ok <- TRUE
        break
      }
      alpha <- alpha / 2
    }
    if (!ok) break
  }

  w <- weights_at(lambda)
  final_imbalance <- max(abs(as.numeric(crossprod(Zt, w))) / scale_z)
  if (final_imbalance < tol) converged <- TRUE

  list(
    weights = as.numeric(w * q_total),
    converged = converged,
    iterations = iter,
    max_imbalance = final_imbalance
  )
}
