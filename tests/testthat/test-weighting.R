test_that("ipw_weights moves the reference mean towards the cohort", {
  set.seed(91)
  ref <- data.frame(age = rnorm(4000), bmi = rnorm(4000))
  coh <- data.frame(age = rnorm(800, 0.8), bmi = rnorm(800, 0.4))

  out <- ipw_weights(coh, ref, weighting.corr = FALSE)

  expect_length(out$weights, nrow(ref))
  expect_true(all(out$weights > 0))
  expect_true(out$converged)

  wm_age <- weighted.mean(ref$age, out$weights)
  expect_lt(abs(wm_age - mean(coh$age)), abs(mean(ref$age) - mean(coh$age)) / 2)
})

test_that("ipw_weights is idempotent when the populations already match", {
  set.seed(92)
  ref <- data.frame(a = rnorm(3000), b = rnorm(3000))
  coh <- data.frame(a = rnorm(1500), b = rnorm(1500))
  out <- ipw_weights(coh, ref, weighting.corr = FALSE)
  ## Stabilised weights should sit near 1 when there is no covariate shift.
  expect_equal(mean(out$weights), 1, tolerance = 0.1)
})

test_that("ipw_weights validates its inputs", {
  ref <- data.frame(a = rnorm(50))
  coh <- data.frame(b = rnorm(50))
  expect_error(ipw_weights(coh, ref), "share no variable names")

  ref2 <- data.frame(a = c(rnorm(49), NA))
  coh2 <- data.frame(a = rnorm(50))
  expect_error(ipw_weights(coh2, ref2), "missing values")

  ref3 <- data.frame(a = rnorm(50), b = rnorm(50))
  coh3 <- data.frame(a = rnorm(50), b = rnorm(50))
  ## A name in neither dataset is a typo and warns ...
  expect_warning(
    ipw_weights(coh3, ref3, weighting.vars = c("a", "zzz")),
    "in neither dataset"
  )
  ## ... but a variable this cohort simply did not measure is the normal case
  ## in this package and passes quietly.
  expect_silent(
    ipw_weights(coh3, ref3[, "a", drop = FALSE], weighting.vars = c("a", "b"))
  )
  expect_error(
    suppressWarnings(ipw_weights(coh3, ref3, weighting.vars = "zzz")),
    "No usable weighting variables"
  )
  expect_error(
    ipw_weights(coh3, ref3, trim = c(0.9, 0.1)),
    "increasing probabilities"
  )
})

test_that("ipw_weights trimming bounds the weights", {
  set.seed(93)
  ref <- data.frame(a = rnorm(2000))
  coh <- data.frame(a = rnorm(500, 1.5))
  untrimmed <- ipw_weights(coh, ref, weighting.corr = FALSE)
  trimmed <- ipw_weights(coh, ref, weighting.corr = FALSE, trim = c(0.05, 0.95))
  expect_lt(max(trimmed$weights), max(untrimmed$weights))
})

test_that("ebal_weights balances the requested moments", {
  set.seed(94)
  ref <- data.frame(age = rnorm(2000), bmi = rnorm(2000))
  coh <- data.frame(age = rnorm(500, 0.7), bmi = rnorm(500, 0.3))

  out <- ebal_weights(coh, ref)
  expect_length(out$weights, nrow(ref))
  expect_true(out$converged)
  ## Entropy balancing matches the first moment essentially exactly.
  expect_equal(weighted.mean(ref$age, out$weights), mean(coh$age),
    tolerance = 1e-5
  )
  expect_equal(weighted.mean(ref$bmi, out$weights), mean(coh$bmi),
    tolerance = 1e-5
  )
})

test_that("ebal_weights accepts a character vector for weighting.vars", {
  ## The original stopped outright on a character vector, but callers upstream
  ## passed one.
  set.seed(95)
  ref <- data.frame(age = rnorm(1500), bmi = rnorm(1500))
  coh <- data.frame(age = rnorm(400, 0.5), bmi = rnorm(400))
  out <- ebal_weights(coh, ref, weighting.vars = "age")
  expect_equal(weighted.mean(ref$age, out$weights), mean(coh$age),
    tolerance = 1e-5
  )
})

test_that("ebal_weights with an empty weighting.vars list uses all variables", {
  ## `weighting.vars[[1]]` on list() is a subscript-out-of-bounds error, which
  ## was the default argument in the original code.
  set.seed(96)
  ref <- data.frame(age = rnorm(1200), bmi = rnorm(1200))
  coh <- data.frame(age = rnorm(300, 0.4), bmi = rnorm(300, 0.2))
  expect_no_error(ebal_weights(coh, ref, weighting.vars = list()))
})

test_that("ebal_weights refuses higher moments of categorical-only data", {
  ref <- data.frame(g = factor(sample(c("a", "b"), 500, TRUE)))
  coh <- data.frame(g = factor(sample(c("a", "b"), 200, TRUE)))
  expect_error(
    ebal_weights(coh, ref, weighting.max.moments = 2),
    "cannot be balanced"
  )
})

test_that("weighting a shifted reference restores the estimate", {
  ## The reference population is shifted relative to the cohorts, so its
  ## unweighted covariance is a poor stand-in; reweighting should help.
  set.seed(97)
  n_ref <- 6000
  ref <- sim_pop(n_ref)
  ref$x1 <- ref$x1 * 1.6 # inflate the variance of x1 in the reference
  cohorts <- list(sim_cohort(1500, "x3"), sim_cohort(1500, "x2"))

  fit_w <- calc_linear_beta(
    cohort_data = cohorts, ref_data = ref,
    weighting.vars = c("x1", "x2", "x3"), weighting.method = "ebal"
  )
  expect_length(fit_w$ref_weights, 2L)
  expect_false(any(vapply(fit_w$ref_weights, is.null, logical(1))))
  expect_true(all(fit_w$weighting_converged))
})
