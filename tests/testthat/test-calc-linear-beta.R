test_that("calc_linear_beta recovers the truth end to end", {
  s <- sim_setup(n_s = 2000, n_ref = 10000, seed = 71)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)

  expect_s3_class(fit, "genmetalite")
  expect_equal(fit$linear_est$beta, true_beta, tolerance = 0.08)

  ## The IVW comparator is badly biased here, by construction.
  expect_gt(
    max(abs(fit$IVW_est$beta[names(true_beta)] - true_beta)),
    3 * max(abs(fit$linear_est$beta - true_beta))
  )
})

test_that("summary() gives a usable coefficient table", {
  s <- sim_setup(n_s = 1000, n_ref = 5000, seed = 72)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  tab <- summary(fit)

  expect_s3_class(tab, "data.frame")
  expect_identical(tab$term, c("x1", "x2", "x3"))
  expect_true(all(tab$std.error > 0))
  expect_true(all(tab$conf.low < tab$estimate))
  expect_true(all(tab$conf.high > tab$estimate))
  ## Nominal 95% intervals should cover the truth at these sample sizes.
  expect_true(all(tab$conf.low < true_beta & tab$conf.high > true_beta))
})

test_that("pre-fitted reduced models give the same answer as raw data", {
  s <- sim_setup(n_s = 1000, n_ref = 5000, seed = 73)
  rms <- lapply(s$cohorts, function(d) {
    run_regression(d$y, d[, setdiff(names(d), "y"), drop = FALSE])
  })

  from_raw <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  from_fits <- suppressWarnings(
    calc_linear_beta(reduced_models = rms, ref_data = s$ref)
  )
  ## Identical point estimates; the variance differs only by the study-block
  ## term, which needs raw data.
  expect_equal(from_raw$linear_est$beta, from_fits$linear_est$beta)
  expect_true(all(diag(from_raw$linear_est$var) >=
    diag(from_fits$linear_est$var) - 1e-12))
})

test_that("cov_p_type = 'reference' also recovers the truth", {
  s <- sim_setup(n_s = 2000, n_ref = 10000, seed = 74)
  fit <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = s$ref, cov_p_type = "reference"
  )
  expect_equal(fit$linear_est$beta, true_beta, tolerance = 0.08)
})

test_that("supplying a covariance matrix instead of reference rows works", {
  s <- sim_setup(n_s = 1000, n_ref = 8000, seed = 75)
  cov_ref <- cov(s$ref)
  fit <- suppressWarnings(calc_linear_beta(
    cohort_data = s$cohorts, cov_ref_X = cov_ref
  ))
  expect_equal(fit$linear_est$beta, true_beta, tolerance = 0.1)
  ## Without reference rows the reference variance term is unavailable, but the
  ## study-block term is still computed from the raw cohort data.
  expect_true(all(diag(fit$linear_est$var) >= diag(fit$linear_est$M)))

  ## With neither, the variance collapses to the naive GLS variance.
  bare <- suppressWarnings(calc_linear_beta(
    reduced_models = fit$reduced_models, cov_ref_X = cov_ref
  ))
  expect_equal(bare$linear_est$var, unname(bare$linear_est$M),
    ignore_attr = TRUE
  )
})

test_that("constraints propagate from calc_linear_beta", {
  s <- sim_setup(n_s = 1000, n_ref = 5000, seed = 76)
  fit <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = s$ref, constr_vars = "x3"
  )
  expect_equal(unname(fit$linear_est$beta[["x3"]]), 0)
})

test_that("results are invariant to shifting the covariates", {
  ## The estimator works with Var(X) and slope coefficients, both invariant to
  ## location. This pins down the removal of the old mean-centering step.
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 77)
  shift <- function(d) {
    for (v in intersect(names(d), c("x1", "x2", "x3"))) d[[v]] <- d[[v]] + 10
    d
  }
  a <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  b <- calc_linear_beta(
    cohort_data = lapply(s$cohorts, shift), ref_data = shift(s$ref)
  )
  expect_equal(a$linear_est$beta, b$linear_est$beta, tolerance = 1e-8)
  expect_equal(a$linear_est$var, b$linear_est$var, tolerance = 1e-8)
})

test_that("missing inputs raise informative errors", {
  s <- sim_setup(n_s = 200, n_ref = 800, seed = 78)
  expect_error(calc_linear_beta(ref_data = s$ref), "cohort_data")
  expect_error(calc_linear_beta(cohort_data = s$cohorts), "ref_data")
})

test_that("a cohort without a column named y is an error", {
  s <- sim_setup(n_s = 200, n_ref = 800, seed = 79)
  bad <- s$cohorts
  names(bad[[1]])[names(bad[[1]]) == "y"] <- "outcome"
  expect_error(
    calc_linear_beta(cohort_data = bad, ref_data = s$ref),
    "no column named 'y'"
  )
})

test_that("weighting without cohort data warns and is skipped", {
  s <- sim_setup(n_s = 400, n_ref = 2000, seed = 80)
  rms <- lapply(s$cohorts, function(d) {
    run_regression(d$y, d[, setdiff(names(d), "y"), drop = FALSE])
  })
  msgs <- character(0)
  withCallingHandlers(
    calc_linear_beta(
      reduced_models = rms, ref_data = s$ref, weighting.vars = "x1"
    ),
    warning = function(w) {
      msgs <<- c(msgs, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  expect_true(any(grepl(
    "Reweighting needs individual-level cohort data", msgs
  )))
})

test_that("the EB estimate is returned and sits between the two", {
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 81)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  lin <- fit$linear_est$beta
  ivw <- fit$IVW_est$beta[names(lin)]
  eb <- fit$EB_est$beta
  expect_true(all(eb >= pmin(lin, ivw) - 1e-12))
  expect_true(all(eb <= pmax(lin, ivw) + 1e-12))
})

test_that("print methods run without error", {
  s <- sim_setup(n_s = 300, n_ref = 1500, seed = 82)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  expect_output(print(fit), "genmetalite")
  expect_output(print(fit$linear_est), "cov_p_type")
  expect_output(print(fit$reduced_models[[1]]), "reduced_model")
})
