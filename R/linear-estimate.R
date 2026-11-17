#' Map and pool reduced-model coefficients
#'
#' Maps each study's reduced-model coefficient vector onto the full-model
#' parameter using covariate second moments, pools the mapped vectors by
#' generalized least squares, and estimates the variance of the pooled estimate
#' with an influence-function sandwich.
#'
#' @param n_cohorts Number of studies.
#' @param reg_out List of reduced models; see [run_regression()].
#' @param sigma_X_ls List of per-study covariance matrices over the full
#'   covariate set, as produced by [make_sigma_X_ls()].
#' @param cov_cohort_X_s_ls List of study-side covariance matrices, also from
#'   [make_sigma_X_ls()].
#' @param cohort_data Optional list of individual-level cohort data frames, each
#'   containing the outcome in a column named `y`. Required for the study-block
#'   variance term under `cov_p_type = "partitioned"`.
#' @param ref_pop_dummy Optional list of dummy-expanded reference design
#'   matrices. Required for the reference variance term; without it the variance
#'   is the naive GLS variance and is too small unless the reference sample is
#'   very large.
#' @param ref_weights Optional list of reference weight vectors.
#' @param cov_p_type `"partitioned"` or `"reference"`; must match what was used
#'   to build `sigma_X_ls`.
#' @param constr_vars Optional character vector of covariate names to constrain
#'   to zero.
#' @param jacobian `"numerical"` (forward finite differences, the default) or
#'   `"analytic"` (closed form). They agree to numerical precision; the analytic
#'   form is much faster when `p` is large.
#' @param check_jacobian If `TRUE`, compute both and report the largest absolute
#'   discrepancy. Useful as a regression test, expensive otherwise.
#' @param survey_strata,survey_psu Optional stratum and primary-sampling-unit
#'   identifiers for the reference rows, for complex-survey reference samples
#'   such as NHANES. Both must be aligned to the rows of `ref_pop_dummy[[1]]`.
#' @param lonely_psu How to treat strata containing a single PSU: `"remove"`
#'   (contribute nothing), `"adjust"` (centre at the grand mean), or `"fail"`.
#' @param eps Step size for the numerical Jacobian.
#'
#' @return A list with
#'   \item{beta}{named vector of pooled slope coefficients.}
#'   \item{var}{estimated variance-covariance matrix of `beta`.}
#'   \item{M}{the naive (GLS) variance, i.e. `var` without the reference and
#'     study-block terms.}
#'   \item{A}{list of the per-study design maps \eqn{A_s}.}
#'   \item{slope_names}{covariate names, in the order used throughout.}
#'   \item{cov_p_type}{as supplied.}
#'   \item{design_df}{complex-survey degrees of freedom, or `NULL`.}
#'
#' @section Method:
#' Everything runs in the *slope-only* basis. After centering, the second-moment
#' matrix is block diagonal in the intercept and the slopes, so the intercept
#' decouples: dropping it keeps every matrix that has to be inverted
#' non-singular and means no intercept has to be supplied in the reference data.
#'
#' For study \eqn{s} with observed covariates \eqn{o} and unobserved covariates
#' \eqn{m}, the design map is
#' \deqn{A_s = [\, I_{k_s} \; | \; B_s \,], \qquad
#'       B_s = \Sigma_{oo}^{-1}\,\Sigma_{om},}
#' read from `sigma_X_ls[[s]]`. Under `cov_p_type = "partitioned"` the
#' \eqn{oo} block of that matrix is the study's own covariance and only
#' \eqn{\Sigma_{om}} is borrowed; under `"reference"` both come from the
#' reference. The pooled estimate is the GLS solution
#' \eqn{\hat\beta = M \sum_s A_s' V_s^{-1} \hat\theta_s} with
#' \eqn{M = (\sum_s A_s' V_s^{-1} A_s)^{-1}}.
#'
#' The variance has up to three terms:
#' \itemize{
#'   \item **(I) naive**, \eqn{M} itself, treating \eqn{A_s} as fixed;
#'   \item **(II) reference**, the extra uncertainty from estimating
#'     \eqn{\Sigma} in a finite reference sample, obtained by pushing the
#'     centered second-moment influence function
#'     \eqn{(x - \mu)(x - \mu)' - \Sigma} through the Jacobian
#'     \eqn{\partial\beta / \partial\,\mathrm{vec}\,\Sigma}. This term is
#'     \eqn{O(1/n_{\mathrm{ref}})} and is omitted, with a warning, when
#'     `ref_pop_dummy` is not supplied;
#'   \item **(I-b) study block**, present only under `"partitioned"`, because
#'     the \eqn{oo} block of \eqn{A_s} is then estimated on the same sample as
#'     \eqn{\hat\theta_s}. It is \eqn{O(1/n_s)} and is exactly zero when a study
#'     measures every covariate or when the coefficients of the covariates it
#'     misses are zero.
#' }
#'
#' @section Weights are treated as fixed:
#' If the reference population was reweighted by [ebal_weights()] or
#' [ipw_weights()], the estimation of those weights contributes a further
#' influence-function term that is **not** included here. Reported standard
#' errors are therefore mildly anti-conservative under reweighting; a joint
#' bootstrap over study and reference rows is the recommended check.
#'
#' @seealso [calc_linear_beta()], [make_sigma_X_ls()], [IVW_estimate()]
#' @export
linear_estimate <- function(n_cohorts,
                            reg_out,
                            sigma_X_ls,
                            cov_cohort_X_s_ls,
                            cohort_data = NULL,
                            ref_pop_dummy = NULL,
                            ref_weights = NULL,
                            cov_p_type = c("partitioned", "reference"),
                            constr_vars = NULL,
                            jacobian = c("numerical", "analytic"),
                            check_jacobian = FALSE,
                            survey_strata = NULL,
                            survey_psu = NULL,
                            lonely_psu = c("remove", "adjust", "fail"),
                            eps = 1e-6) {
  cov_p_type <- match.arg(cov_p_type)
  jacobian <- match.arg(jacobian)
  lonely_psu <- match.arg(lonely_psu)

  stopifnot(length(reg_out) == n_cohorts)

  ## ------------------------------------------------------------------
  ## 0. Slope-only basis
  ## ------------------------------------------------------------------
  cov_names <- unique(unlist(lapply(reg_out, function(r) {
    names(r$theta_hat_s)
  })))
  slope_names <- drop_intercept_names(cov_names)
  p <- length(slope_names)
  if (p == 0L) {
    stop("No non-intercept coefficients found in `reg_out`.", call. = FALSE)
  }

  Sig <- lapply(seq_len(n_cohorts), function(s) {
    subset_by_name(
      sigma_X_ls[[s]], slope_names,
      sprintf("assembled covariance for cohort %d", s)
    )
  })

  obs_names <- mis_names <- o_idx <- m_idx <- vector("list", n_cohorts)
  theta_s <- Vi <- Psi_s <- vector("list", n_cohorts)

  for (s in seq_len(n_cohorts)) {
    nm <- names(reg_out[[s]]$theta_hat_s)
    o <- intersect(slope_names, drop_intercept_names(nm))
    if (length(o) == 0L) {
      stop(sprintf("Cohort %d has no slope coefficients.", s), call. = FALSE)
    }
    obs_names[[s]] <- o
    mis_names[[s]] <- setdiff(slope_names, o)
    o_idx[[s]] <- match(o, slope_names)
    m_idx[[s]] <- match(mis_names[[s]], slope_names)
    theta_s[[s]] <- reg_out[[s]]$theta_hat_s[o]

    Vs <- reg_out[[s]]$cov_theta_hat_s
    Vi[[s]] <- safe_solve(
      Vs[o, o, drop = FALSE],
      sprintf("coefficient covariance of cohort %d", s)
    )

    ## Always subset by name with drop = FALSE. The original script special-cased
    ## a single observed covariate and used the whole matrix, which silently
    ## produced the wrong Psi_s whenever a study had more than one column in
    ## cov_cohort_X_s_ls but only one slope coefficient.
    Psi_s[[s]] <- subset_by_name(
      cov_cohort_X_s_ls[[s]], o,
      sprintf("study covariance for cohort %d", s)
    )
  }

  ## ------------------------------------------------------------------
  ## 1. Design map and closed-form estimate
  ## ------------------------------------------------------------------
  make_A <- function(Sig_s, s) {
    oi <- o_idx[[s]]
    mi <- m_idx[[s]]
    k <- length(oi)
    A <- matrix(0, k, p, dimnames = list(obs_names[[s]], slope_names))
    A[, oi] <- diag(k)
    if (length(mi) > 0L) {
      A[, mi] <- safe_solve(
        Sig_s[oi, oi, drop = FALSE],
        sprintf("Sigma_oo block for cohort %d", s),
        b = Sig_s[oi, mi, drop = FALSE]
      )
    }
    A
  }

  betaG_of_Sigma <- function(Sig_list) {
    A <- lapply(seq_len(n_cohorts), function(s) make_A(Sig_list[[s]], s))
    bread <- Reduce(`+`, lapply(seq_len(n_cohorts), function(s) {
      t(A[[s]]) %*% Vi[[s]] %*% A[[s]]
    }))
    M <- safe_solve(bread, "pooled information matrix")
    bvec <- Reduce(`+`, lapply(seq_len(n_cohorts), function(s) {
      t(A[[s]]) %*% Vi[[s]] %*% theta_s[[s]]
    }))
    list(beta = as.numeric(M %*% bvec), M = M, A = A)
  }

  base <- betaG_of_Sigma(Sig)
  betaG <- stats::setNames(base$beta, slope_names)
  M <- base$M
  dimnames(M) <- list(slope_names, slope_names)
  A_ls <- base$A

  V_IF <- M
  design_df <- NULL

  ## ------------------------------------------------------------------
  ## 2. Reference term (II)
  ## ------------------------------------------------------------------
  if (is.null(ref_pop_dummy)) {
    warning(
      "Individual-level reference data were not supplied, so the reference ",
      "term of the variance is omitted. Standard errors will be too small ",
      "unless the reference sample is very large relative to the studies.",
      call. = FALSE
    )
  } else {
    if (!is.list(ref_pop_dummy)) ref_pop_dummy <- list(ref_pop_dummy)
    if (length(ref_pop_dummy) == 1L && n_cohorts > 1L) {
      ref_pop_dummy <- rep(ref_pop_dummy, n_cohorts)
    }
    if (all_null(ref_weights)) {
      ref_weights <- lapply(ref_pop_dummy, function(X) rep(1, nrow(X)))
    }

    ## Which entries of Sigma are reference-sourced, per cohort.
    ref_used <- lapply(seq_len(n_cohorts), function(s) {
      msk <- matrix(FALSE, p, p)
      oi <- o_idx[[s]]
      mi <- setdiff(seq_len(p), oi)
      if (cov_p_type == "reference") {
        msk[oi, ] <- TRUE
      } else if (length(mi) > 0L) {
        msk[oi, mi] <- TRUE
        msk[mi, oi] <- TRUE
      }
      msk
    })

    D_list <- jacobian_by_cohort(
      method = jacobian, n_cohorts = n_cohorts, p = p, Sig = Sig,
      ref_used = ref_used, o_idx = o_idx, A_ls = A_ls, Vi = Vi,
      theta_s = theta_s, betaG = betaG, M = M,
      betaG_of_Sigma = betaG_of_Sigma, eps = eps
    )

    if (isTRUE(check_jacobian)) {
      other <- if (jacobian == "analytic") "numerical" else "analytic"
      D_other <- jacobian_by_cohort(
        method = other, n_cohorts = n_cohorts, p = p, Sig = Sig,
        ref_used = ref_used, o_idx = o_idx, A_ls = A_ls, Vi = Vi,
        theta_s = theta_s, betaG = betaG, M = M,
        betaG_of_Sigma = betaG_of_Sigma, eps = eps
      )
      diffs <- max(vapply(seq_len(n_cohorts), function(s) {
        max(abs(D_list[[s]] - D_other[[s]]))
      }, numeric(1)))
      message(sprintf(
        "Jacobian check (%s vs %s): max abs diff = %.3e (mean |D| = %.3e)",
        jacobian, other, diffs,
        mean(abs(do.call(cbind, D_list)))
      ))
    }

    ref_term <- reference_variance(
      ref_pop_dummy = ref_pop_dummy, ref_weights = ref_weights,
      slope_names = slope_names, Sig = Sig, D_list = D_list,
      n_cohorts = n_cohorts, p = p,
      survey_strata = survey_strata, survey_psu = survey_psu,
      lonely_psu = lonely_psu
    )
    V_IF <- V_IF + ref_term$cov
    design_df <- ref_term$df
  }

  ## ------------------------------------------------------------------
  ## 3. Study-block term (I-b), partitioned only
  ## ------------------------------------------------------------------
  if (cov_p_type == "partitioned") {
    if (is.null(cohort_data)) {
      warning(
        "cov_p_type = 'partitioned' but individual-level `cohort_data` were ",
        "not supplied, so the study-block variance term is omitted. The ",
        "variance will be too small when studies miss covariates whose ",
        "coefficients are non-zero.",
        call. = FALSE
      )
    } else {
      V_IF <- V_IF + study_block_variance(
        n_cohorts = n_cohorts, cohort_data = cohort_data, reg_out = reg_out,
        obs_names = obs_names, mis_names = mis_names, Sig = Sig,
        Psi_s = Psi_s, A_ls = A_ls, Vi = Vi, theta_s = theta_s,
        betaG = betaG, M = M, p = p
      )
    }
  }

  V_IF <- symmetrise(as.matrix(V_IF))
  dimnames(V_IF) <- list(slope_names, slope_names)

  ## ------------------------------------------------------------------
  ## 4. Optional zero constraints
  ## ------------------------------------------------------------------
  beta_out <- betaG
  if (!is.null(constr_vars)) {
    unknown <- setdiff(constr_vars, slope_names)
    if (length(unknown) > 0L) {
      stop(sprintf(
        "Constraint variable(s) %s do not match any covariate name.",
        paste(unknown, collapse = ", ")
      ), call. = FALSE)
    }
    C <- matrix(0, nrow = length(constr_vars), ncol = p)
    for (i in seq_along(constr_vars)) {
      C[i, match(constr_vars[i], slope_names)] <- 1
    }
    E <- safe_solve(C %*% M %*% t(C), "constraint matrix C M C'")
    P <- diag(p) - M %*% t(C) %*% E %*% C

    beta_out <- stats::setNames(as.numeric(P %*% betaG), slope_names)
    beta_out[abs(beta_out) < 1e-10] <- 0

    ## The constrained variance sandwiches the FULL influence-function variance
    ## by the projector. (The original script computed a constrained estimate
    ## and then returned the unconstrained one, so constraints had no effect.)
    V_IF <- P %*% V_IF %*% t(P)
    V_IF[abs(V_IF) < 1e-10] <- 0
    dimnames(V_IF) <- list(slope_names, slope_names)
  }

  structure(
    list(
      beta = beta_out,
      var = as.matrix(V_IF),
      M = M,
      A = A_ls,
      slope_names = slope_names,
      cov_p_type = cov_p_type,
      design_df = design_df,
      constr_vars = constr_vars
    ),
    class = "genmetalite_estimate"
  )
}

#' @rdname linear_estimate
#' @export
linear_estimate_centered <- linear_estimate

#' @export
print.genmetalite_estimate <- function(x, ...) {
  se <- sqrt(pmax(diag(x$var), 0))
  out <- cbind(
    Estimate = x$beta,
    `Std. Error` = se,
    `z value` = ifelse(se > 0, x$beta / se, NA_real_)
  )
  cat(sprintf(
    "<genmetalite_estimate>  cov_p_type = %s,  %d covariate(s)\n",
    x$cov_p_type, length(x$beta)
  ))
  print(round(out, 4L))
  invisible(x)
}


## Internal: Jacobians --------------------------------------------------------

#' Per-cohort Jacobian of beta with respect to vec(Sigma_s)
#'
#' Returns a list of `n_cohorts` matrices, each `p x p^2`. Splitting the
#' Jacobian by cohort (rather than perturbing every cohort's Sigma at once)
#' matters when the reference has been reweighted differently for each cohort,
#' so that the cohorts do not share one Sigma. When they do share it, the sum of
#' the per-cohort Jacobians equals the single joint Jacobian, so results are
#' unchanged.
#' @noRd
jacobian_by_cohort <- function(method, n_cohorts, p, Sig, ref_used, o_idx,
                               A_ls, Vi, theta_s, betaG, M, betaG_of_Sigma,
                               eps) {
  if (method == "numerical") {
    lapply(seq_len(n_cohorts), function(s) {
      D <- matrix(0, p, p * p)
      for (a in seq_len(p)) {
        for (b in seq_len(p)) {
          if (!ref_used[[s]][a, b]) next
          Sig_pert <- Sig
          Sig_pert[[s]][a, b] <- Sig_pert[[s]][a, b] + eps
          bp <- betaG_of_Sigma(Sig_pert)$beta
          D[, (b - 1L) * p + a] <- (bp - betaG) / eps
        }
      }
      D
    })
  } else {
    lapply(seq_len(n_cohorts), function(s) {
      D <- matrix(0, p, p * p)
      oi <- o_idx[[s]]
      r_s <- theta_s[[s]] - A_ls[[s]] %*% betaG
      Sig_oo <- Sig[[s]][oi, oi, drop = FALSE]
      ImDA <- diag(p)
      ImDA[oi, ] <- ImDA[oi, ] - A_ls[[s]]
      for (a in seq_len(p)) {
        if (!(a %in% oi)) next
        for (b in seq_len(p)) {
          if (!ref_used[[s]][a, b]) next
          dSig_rows <- matrix(0, length(oi), p)
          dSig_rows[match(a, oi), b] <- 1
          dA <- safe_solve(
            Sig_oo, "Sigma_oo block", b = dSig_rows %*% ImDA
          )
          dpsi <- t(dA) %*% Vi[[s]] %*% r_s -
            t(A_ls[[s]]) %*% Vi[[s]] %*% (dA %*% betaG)
          D[, (b - 1L) * p + a] <- as.numeric(M %*% dpsi)
        }
      }
      D
    })
  }
}


## Internal: variance terms ---------------------------------------------------

#' Reference-sample contribution to the variance
#'
#' Each cohort contributes a per-reference-row linearized influence
#' `q_s = (w_s / sum(w_s)) * (z_s - zbar_s) %*% t(D_s)`. When the cohorts share
#' one reference sample those rows are added, giving the total influence of each
#' reference observation on the pooled estimate, and the design (if any) is
#' applied to the sum. Adding first and applying the design second is what makes
#' per-cohort weights compatible with a shared survey design: the cohorts'
#' contributions are correlated because they come from the same rows, and
#' summing before squaring keeps that correlation.
#' @noRd
reference_variance <- function(ref_pop_dummy, ref_weights, slope_names, Sig,
                               D_list, n_cohorts, p,
                               survey_strata, survey_psu, lonely_psu) {
  Xs <- lapply(seq_len(n_cohorts), function(s) {
    X <- ref_pop_dummy[[s]]
    missing <- setdiff(slope_names, colnames(X))
    if (length(missing) > 0L) {
      stop(sprintf(
        "Reference data for cohort %d are missing column(s): %s.",
        s, paste(missing, collapse = ", ")
      ), call. = FALSE)
    }
    as.matrix(X[, slope_names, drop = FALSE])
  })

  nrows <- vapply(Xs, nrow, integer(1))
  same_rows <- length(unique(nrows)) == 1L
  use_survey <- !is.null(survey_strata) || !is.null(survey_psu)

  if (use_survey && !same_rows) {
    stop(
      "A survey design describes one reference sample, but the cohorts were ",
      "given reference datasets of different sizes. Supply a single ",
      "ref_design()/survey.design object, or drop the design information.",
      call. = FALSE
    )
  }

  contrib <- lapply(seq_len(n_cohorts), function(s) {
    X <- Xs[[s]]
    w <- ref_weights[[s]] %||% rep(1, nrow(X))
    if (length(w) != nrow(X)) {
      stop(sprintf(
        "ref_weights[[%d]] has length %d but the reference matrix has %d rows.",
        s, length(w), nrow(X)
      ), call. = FALSE)
    }
    mu <- colSums(w * X) / sum(w)
    Xc <- sweep(X, 2L, mu, "-")
    XX <- Xc[, rep(seq_len(p), each = p), drop = FALSE] *
      Xc[, rep(seq_len(p), times = p), drop = FALSE]
    Zc <- sweep(XX, 2L, as.numeric(Sig[[s]]), "-")
    IF <- Zc %*% t(D_list[[s]])
    IFbar <- colSums(w * IF) / sum(w)
    ## Weight-scaled and centred: the cohort's contribution to the estimate is
    ## colSums(q), so contributions from different cohorts are directly additive.
    (w * sweep(IF, 2L, IFbar, "-")) / sum(w)
  })

  if (!same_rows) {
    ## Genuinely different reference samples per cohort: independent, so add
    ## their variances rather than their influence rows.
    cov_ref <- Reduce(`+`, lapply(contrib, crossprod))
    return(list(cov = cov_ref, df = NULL))
  }

  U <- Reduce(`+`, contrib)

  if (!use_survey) {
    return(list(cov = crossprod(U), df = NULL))
  }

  n_ref <- nrow(U)
  for (nm in c("survey_strata", "survey_psu")) {
    v <- get(nm)
    if (!is.null(v) && length(v) != n_ref) {
      stop(sprintf(
        "`%s` has length %d but the reference sample has %d rows.",
        nm, length(v), n_ref
      ), call. = FALSE)
    }
  }

  sv <- psu_linearized_cov(U, survey_strata, survey_psu, lonely_psu)
  list(cov = sv$cov, df = sv$df)
}

#' Study-block contribution to the variance (partitioned mode only)
#' @noRd
study_block_variance <- function(n_cohorts, cohort_data, reg_out, obs_names,
                                 mis_names, Sig, Psi_s, A_ls, Vi, theta_s,
                                 betaG, M, p) {
  Ib <- matrix(0, p, p)
  for (s in seq_len(n_cohorts)) {
    o <- obs_names[[s]]
    m <- mis_names[[s]]
    if (length(m) == 0L) next # complete study: no correction
    beta_m <- betaG[m]
    if (all(abs(beta_m) < 1e-12)) next # nothing missing that matters

    dgp <- reg_out[[s]]$dgp %||% "linear"
    if (!identical(dgp, "linear") && !identical(dgp, "gaussian")) {
      warning(sprintf(
        paste0(
          "Cohort %d was fitted with dgp = '%s'. The study-block variance term ",
          "is derived for least-squares residuals and is skipped for ",
          "non-linear links."
        ),
        s, dgp
      ), call. = FALSE)
      next
    }

    df <- cohort_data[[s]]
    if (!"y" %in% names(df)) {
      warning(sprintf(
        "cohort_data[[%d]] has no column 'y'; skipping its study-block term.", s
      ), call. = FALSE)
      next
    }
    yv <- df[["y"]]
    Xall <- dummy_matrix(df[, setdiff(names(df), "y"), drop = FALSE])

    if (!all(o %in% colnames(Xall))) {
      warning(sprintf(
        paste0(
          "Cohort %d: could not align design column(s) %s between the raw data ",
          "and the reduced model; skipping its study-block term."
        ),
        s, paste(setdiff(o, colnames(Xall)), collapse = ", ")
      ), call. = FALSE)
      next
    }

    Xo <- scale(Xall[, o, drop = FALSE], center = TRUE, scale = FALSE)
    ns <- nrow(Xo)
    epsv <- as.numeric(yv - mean(yv) - Xo %*% theta_s[[s]])

    Phi_om <- Sig[[s]][o, m, drop = FALSE]
    Psi_inv <- safe_solve(Psi_s[[s]], sprintf("study covariance %d", s))
    Bs <- Psi_inv %*% Phi_om
    g <- as.numeric(Bs %*% beta_m)
    Pref <- t(A_ls[[s]]) %*% Vi[[s]] %*% Psi_inv

    C_rows <- (Xo * epsv) %*% t(Pref)
    DX <- Xo * as.numeric(Xo %*% g)
    DX <- sweep(DX, 2L, as.numeric(Psi_s[[s]] %*% g), "-")
    D_rows <- DX %*% t(Pref)

    cov_dd <- stats::cov(D_rows)
    cov_cd <- stats::cov(C_rows, D_rows)
    Ib <- Ib + (M %*% (cov_dd + cov_cd + t(cov_cd)) %*% M) / ns
  }
  Ib
}
