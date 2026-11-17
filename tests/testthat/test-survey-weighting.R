## Entropy balancing a survey reference: the case where a per-cohort weight
## vector has to coexist with a shared sampling design.

svy_setup <- function(seed = 401, n_str = 6, n_psu = 5, m = 60, n_s = 900) {
  set.seed(seed)
  N <- n_str * n_psu * m
  psu <- rep(seq_len(n_str * n_psu), each = m)
  strata <- rep(seq_len(n_str), each = n_psu * m)
  u <- rnorm(n_str * n_psu, sd = 0.9)[psu]
  x1 <- u + rnorm(N)
  x2 <- 0.6 * x1 + rnorm(N, sd = 0.8)
  x3 <- 0.4 * x2 + rnorm(N)
  list(
    ref = ref_design(
      data.frame(x1 = x1, x2 = x2, x3 = x3),
      weights = exp(rnorm(N, 0, 0.4)) * 100,
      strata = strata, psu = psu
    ),
    cohorts = list(
      sim_cohort(n_s, "x3"), sim_cohort(n_s, "x2"), sim_cohort(n_s, "x1")
    )
  )
}

test_that("entropy balancing composes with a survey design", {
  s <- svy_setup()
  fit <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = s$ref,
    weighting.vars = c("x1", "x2", "x3"), weighting.method = "ebal"
  )

  expect_true(all(fit$weighting_converged))
  expect_length(fit$ref_weights, 3L)
  ## One weight vector per cohort: the design is shared, the tilt is not.
  expect_false(isTRUE(all.equal(fit$ref_weights[[1]], fit$ref_weights[[2]])))
  ## The design still applies.
  expect_identical(fit$linear_est$design_df, 30L - 6L)
  expect_identical(fit$ref_strata, s$ref$strata)
  expect_identical(fit$ref_psu, s$ref$psu)
})

test_that("the combined weights balance each cohort and keep the survey total", {
  s <- svy_setup(seed = 402)
  fit <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = s$ref,
    weighting.vars = c("x1", "x2", "x3"), weighting.method = "ebal"
  )

  ref_df <- s$ref$data
  for (k in seq_along(s$cohorts)) {
    w <- fit$ref_weights[[k]]
    expect_equal(sum(w), sum(s$ref$weights), tolerance = 1e-6)
    ## Every covariate the cohort actually measured is matched exactly.
    for (v in intersect(names(ref_df), names(s$cohorts[[k]]))) {
      expect_equal(
        weighted.mean(ref_df[[v]], w), mean(s$cohorts[[k]][[v]]),
        tolerance = 1e-5
      )
    }
  }
})

test_that("per-cohort weights plus a design no longer error", {
  ## This combination used to be rejected outright: the survey path required a
  ## single weight vector shared by every cohort.
  s <- svy_setup(seed = 403)
  expect_no_error(
    calc_linear_beta(
      cohort_data = s$cohorts, ref_data = s$ref,
      weighting.vars = c("x1", "x2"), weighting.method = "ebal"
    )
  )
})

test_that("IPW also composes with a survey design", {
  s <- svy_setup(seed = 404)
  fit <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = s$ref,
    weighting.vars = c("x1", "x2", "x3"), weighting.method = "IPW"
  )
  expect_length(fit$ref_weights, 3L)
  expect_equal(sum(fit$ref_weights[[1]]), sum(s$ref$weights), tolerance = 1e-6)
  expect_identical(fit$linear_est$design_df, 30L - 6L)
})

test_that("balancing changes what the design contributes", {
  ## Without reweighting, clustering inflates the reference term (pinned in
  ## test-reference-design.R). After entropy balancing it need not: balancing
  ## matches the covariate means exactly, which removes precisely the
  ## between-PSU variation in those means that clustering would have inflated.
  ## So assert what is actually true -- the point estimate is untouched and both
  ## reference terms are valid covariances -- rather than a uniform inequality.
  s <- svy_setup(seed = 405)
  no_design <- ref_design(s$ref$data, weights = s$ref$weights)

  a <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = no_design,
    weighting.vars = c("x1", "x2", "x3")
  )
  b <- calc_linear_beta(
    cohort_data = s$cohorts, ref_data = s$ref,
    weighting.vars = c("x1", "x2", "x3")
  )

  ## Sigma depends on the weights alone, so the estimate is identical.
  expect_equal(a$linear_est$beta, b$linear_est$beta)
  expect_null(a$linear_est$design_df)
  expect_identical(b$linear_est$design_df, 30L - 6L)

  for (fit in list(a, b)) {
    ref_term <- fit$linear_est$var - fit$linear_est$M
    expect_true(all(diag(ref_term) > 0))
    expect_true(all(eigen(ref_term, only.values = TRUE)$values > -1e-10))
  }
})

test_that("without reweighting, clustering does inflate the reference term", {
  s <- svy_setup(seed = 409)
  a <- calc_linear_beta(
    cohort_data = s$cohorts,
    ref_data = ref_design(s$ref$data, weights = s$ref$weights)
  )
  b <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  expect_true(all(
    diag(b$linear_est$var - b$linear_est$M) >
      diag(a$linear_est$var - a$linear_est$M)
  ))
})

test_that("the summed influence keeps cross-cohort correlation", {
  ## Cohorts share the same reference rows, so their contributions are
  ## correlated. Summing the influence rows before applying the design captures
  ## that; adding the per-cohort variances separately would not.
  s <- svy_setup(seed = 406)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  ref_term <- fit$linear_est$var - fit$linear_est$M
  ## The reference term must still be a valid covariance matrix.
  expect_true(isSymmetric(unname(ref_term), tol = 1e-8))
  expect_true(all(eigen(ref_term, only.values = TRUE)$values > -1e-10))
})

test_that("a design cannot be combined with differently sized references", {
  s <- svy_setup(seed = 407)
  refs <- list(
    s$ref$data, s$ref$data[1:100, ], s$ref$data
  )
  expect_error(
    calc_linear_beta(
      cohort_data = s$cohorts, ref_data = refs,
      survey_strata = s$ref$strata, survey_psu = s$ref$psu
    ),
    "one reference sample"
  )
})

test_that("summary() reports design-aware standard errors", {
  s <- svy_setup(seed = 408)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  tab <- summary(fit)
  expect_true(all(tab$std.error > 0))
  expect_identical(tab$term, c("x1", "x2", "x3"))
})
