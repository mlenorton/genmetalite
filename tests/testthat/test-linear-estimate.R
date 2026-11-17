fit_pieces <- function(cohorts, ref, cov_p_type = "partitioned") {
  rms <- lapply(cohorts, function(d) {
    run_regression(d$y, d[, setdiff(names(d), "y"), drop = FALSE])
  })
  sig <- make_sigma_X_ls(
    n_cohorts = length(cohorts), reg_out = rms, ref_data = ref,
    cov_p_type = cov_p_type, cohort_data = cohorts
  )
  list(rms = rms, sig = sig)
}

test_that("linear_estimate recovers the full-model coefficients", {
  s <- sim_setup(n_s = 2000, n_ref = 10000, seed = 31)
  p <- fit_pieces(s$cohorts, s$ref)

  est <- linear_estimate(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = s$cohorts,
    ref_pop_dummy = p$sig$ref_pop_dummy
  )

  expect_s3_class(est, "genmetalite_estimate")
  expect_identical(names(est$beta), c("x1", "x2", "x3"))
  expect_equal(est$beta, true_beta, tolerance = 0.1)
})

test_that("the estimate is much less biased than IVW when covariates are omitted", {
  s <- sim_setup(n_s = 2000, n_ref = 10000, seed = 32)
  p <- fit_pieces(s$cohorts, s$ref)
  est <- linear_estimate(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = s$cohorts, ref_pop_dummy = p$sig$ref_pop_dummy
  )
  ivw <- IVW_estimate(p$rms)

  err_lin <- max(abs(est$beta - true_beta))
  err_ivw <- max(abs(ivw$beta[names(true_beta)] - true_beta))
  expect_lt(err_lin, err_ivw / 3)
})

test_that("the variance is symmetric and positive definite", {
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 33)
  p <- fit_pieces(s$cohorts, s$ref)
  est <- linear_estimate(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = s$cohorts, ref_pop_dummy = p$sig$ref_pop_dummy
  )
  expect_true(isSymmetric(unname(est$var)))
  expect_true(all(eigen(est$var, only.values = TRUE)$values > 0))
  ## The full variance must exceed the naive GLS variance.
  expect_true(all(diag(est$var) >= diag(est$M)))
})

test_that("analytic and numerical Jacobians agree", {
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 34)
  p <- fit_pieces(s$cohorts, s$ref)
  args <- list(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = s$cohorts, ref_pop_dummy = p$sig$ref_pop_dummy
  )
  num <- do.call(linear_estimate, c(args, list(jacobian = "numerical")))
  ana <- do.call(linear_estimate, c(args, list(jacobian = "analytic")))

  expect_equal(num$beta, ana$beta)
  expect_equal(num$var, ana$var, tolerance = 1e-5)
})

test_that("analytic and numerical Jacobians agree in reference mode too", {
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 35)
  p <- fit_pieces(s$cohorts, s$ref, cov_p_type = "reference")
  args <- list(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    ref_pop_dummy = p$sig$ref_pop_dummy, cov_p_type = "reference"
  )
  num <- do.call(linear_estimate, c(args, list(jacobian = "numerical")))
  ana <- do.call(linear_estimate, c(args, list(jacobian = "analytic")))
  expect_equal(num$var, ana$var, tolerance = 1e-5)
})

test_that("omitting the reference rows warns and drops the reference term", {
  s <- sim_setup(n_s = 600, n_ref = 3000, seed = 36)
  p <- fit_pieces(s$cohorts, s$ref)

  ## Regression test: the original code referred to `betaG` and `V_IF`, which
  ## only existed on the other branch, so this path failed outright with
  ## "object not found".
  expect_warning(
    est <- linear_estimate(
      n_cohorts = 3, reg_out = p$rms,
      sigma_X_ls = p$sig$sigma_X_ls,
      cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
      cohort_data = s$cohorts, ref_pop_dummy = NULL
    ),
    "reference term"
  )
  expect_true(is.matrix(est$var))
  expect_equal(est$beta, true_beta, tolerance = 0.12)
  ## The reference term is gone but the study-block term remains, so the
  ## variance still exceeds the naive one.
  expect_true(all(diag(est$var) >= diag(est$M)))
})

test_that("with neither reference rows nor cohort data, var equals M", {
  s <- sim_setup(n_s = 600, n_ref = 3000, seed = 361)
  p <- fit_pieces(s$cohorts, s$ref)
  est <- suppressWarnings(linear_estimate(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = NULL, ref_pop_dummy = NULL
  ))
  expect_equal(est$var, unname(est$M), ignore_attr = TRUE)
})

test_that("constraints are actually applied to the returned estimate", {
  ## Regression test: the original computed a constrained estimate and then
  ## overwrote it with the unconstrained one before returning.
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 37)
  p <- fit_pieces(s$cohorts, s$ref)
  est <- linear_estimate(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = s$cohorts, ref_pop_dummy = p$sig$ref_pop_dummy,
    constr_vars = "x3"
  )
  expect_equal(unname(est$beta[["x3"]]), 0)
  expect_equal(unname(est$var["x3", "x3"]), 0)
  expect_identical(est$constr_vars, "x3")
})

test_that("an unknown constraint variable is an error", {
  s <- sim_setup(n_s = 400, n_ref = 2000, seed = 38)
  p <- fit_pieces(s$cohorts, s$ref)
  expect_error(
    linear_estimate(
      n_cohorts = 3, reg_out = p$rms,
      sigma_X_ls = p$sig$sigma_X_ls,
      cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
      cohort_data = s$cohorts, ref_pop_dummy = p$sig$ref_pop_dummy,
      constr_vars = "nope"
    ),
    "do not match"
  )
})

test_that("complete studies give A_s equal to the identity", {
  set.seed(39)
  ref <- sim_pop(3000)
  cohorts <- list(sim_cohort(800), sim_cohort(800))
  p <- fit_pieces(cohorts, ref)
  est <- linear_estimate(
    n_cohorts = 2, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = cohorts, ref_pop_dummy = p$sig$ref_pop_dummy
  )
  for (A in est$A) expect_equal(unname(A), diag(3))
  ## With no missing covariates, the study-block term is exactly zero, so the
  ## variance should be the naive one plus a reference term that is itself zero
  ## (the Jacobian mask is empty).
  expect_equal(est$var, unname(est$M), ignore_attr = TRUE)
})

test_that("a study observing a single covariate is handled", {
  ## Regression test for the length(o) == 1 special case, which used the whole
  ## study covariance matrix instead of the 1 x 1 block.
  set.seed(40)
  ref <- sim_pop(4000)
  cohorts <- list(
    sim_cohort(1500, c("x2", "x3")),
    sim_cohort(1500, "x1"),
    sim_cohort(1500)
  )
  p <- fit_pieces(cohorts, ref)
  est <- linear_estimate(
    n_cohorts = 3, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = cohorts, ref_pop_dummy = p$sig$ref_pop_dummy
  )
  expect_equal(est$beta, true_beta, tolerance = 0.12)
})

test_that("factor covariates flow through end to end", {
  set.seed(41)
  n_ref <- 6000
  ref <- sim_pop_factor(n_ref)
  b <- c(age = 0.8, edu2 = -0.4, edu3 = 0.6)
  mk <- function(n, drop) {
    d <- sim_pop_factor(n)
    X <- model.matrix(~., d)[, -1, drop = FALSE]
    d$y <- as.numeric(1 + X %*% b + rnorm(n))
    d[, setdiff(names(d), drop), drop = FALSE]
  }
  cohorts <- list(mk(2000, "edu"), mk(2000, character(0)))
  p <- fit_pieces(cohorts, ref)
  est <- linear_estimate(
    n_cohorts = 2, reg_out = p$rms,
    sigma_X_ls = p$sig$sigma_X_ls,
    cov_cohort_X_s_ls = p$sig$cov_cohort_X_s_ls,
    cohort_data = cohorts, ref_pop_dummy = p$sig$ref_pop_dummy
  )
  expect_setequal(names(est$beta), c("age", "edu2", "edu3"))
  expect_equal(est$beta[names(b)], b, tolerance = 0.12)
})
