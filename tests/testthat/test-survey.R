test_that("survey_reference_cov reduces to the iid sandwich", {
  set.seed(61)
  n <- 200
  IF <- matrix(rnorm(n * 3), ncol = 3)
  IF <- sweep(IF, 2, colMeans(IF)) # centered, as the estimator supplies it
  w <- rep(1, n)

  ## One stratum, one PSU per row: the design formula equals the iid formula
  ## up to the finite-population factor n/(n-1).
  out <- survey_reference_cov(IF, weights = w)
  iid <- crossprod(w * IF) / sum(w)^2
  expect_equal(out$cov, iid * n / (n - 1))
  expect_equal(out$df, n - 1L)
  expect_equal(out$n_eff, n) # Kish n_eff is n under equal weights
})

test_that("a single parameter column is handled", {
  ## Regression test: aggregating PSU totals with t(vapply(., numeric(p)))
  ## transposes to 1 x n_psu when p == 1, and then fails to conform with the
  ## p x p accumulator.
  set.seed(611)
  n <- 400
  IF <- matrix(rnorm(n), ncol = 1)
  out <- survey_reference_cov(
    IF,
    weights = runif(n, 0.5, 5),
    strata = rep(1:5, length.out = n),
    psu = rep(1:20, length.out = n)
  )
  expect_identical(dim(out$cov), c(1L, 1L))
  expect_gt(out$cov[1, 1], 0)
  expect_equal(out$df, 15L)
})

test_that("the result is invariant to shifting the influence functions", {
  ## Centering is done internally. A constant shift adds a PSU-size-dependent
  ## amount to each PSU total, which subtracting the stratum mean does not
  ## remove unless every PSU carries the same total weight -- so a caller that
  ## forgot to centre would otherwise get a silently inflated variance.
  set.seed(612)
  n <- 300
  IF <- matrix(rnorm(n * 3), ncol = 3)
  args <- list(
    weights = runif(n, 0.5, 5),
    strata = rep(1:6, length.out = n),
    psu = rep(1:18, length.out = n)
  )
  a <- do.call(survey_reference_cov, c(list(IF), args))
  b <- do.call(survey_reference_cov, c(list(sweep(IF, 2, c(100, -50, 7), "+")), args))
  expect_equal(a$cov, b$cov)
})

test_that("weighting reduces the Kish effective sample size", {
  set.seed(613)
  IF <- matrix(rnorm(200), ncol = 2)
  w <- c(rep(1, 50), rep(20, 50))
  out <- survey_reference_cov(IF, weights = w)
  expect_lt(out$n_eff, 100)
  expect_equal(out$n_eff, sum(w)^2 / sum(w^2))
})

test_that("clustering inflates the variance when PSUs are homogeneous", {
  set.seed(62)
  n_psu <- 40
  m <- 10
  psu_effect <- rnorm(n_psu, sd = 2)
  IF <- matrix(
    rep(psu_effect, each = m) + rnorm(n_psu * m, sd = 0.2),
    ncol = 1
  )
  IF <- IF - mean(IF)
  psu <- rep(seq_len(n_psu), each = m)
  w <- rep(1, n_psu * m)

  clustered <- survey_reference_cov(IF, weights = w, psu = psu)$cov
  ignoring <- survey_reference_cov(IF, weights = w)$cov
  expect_gt(clustered[1, 1], 3 * ignoring[1, 1])
})

test_that("degrees of freedom are PSUs minus strata", {
  IF <- matrix(rnorm(100), ncol = 2)
  out <- survey_reference_cov(
    IF,
    strata = rep(1:5, each = 10), psu = rep(1:10, each = 5)
  )
  expect_equal(out$n_psu, 10L)
  expect_equal(out$n_strata, 5L)
  expect_equal(out$df, 5L)
})

test_that("lonely PSUs are handled per the requested policy", {
  IF <- matrix(rnorm(120), ncol = 2) # 60 rows
  strata <- c(rep(1, 40), rep(2, 20))
  psu <- c(rep(1:4, each = 10), rep(5, 20)) # stratum 2 has a single PSU

  expect_error(
    survey_reference_cov(IF, strata = strata, psu = psu, lonely_psu = "fail"),
    "single PSU"
  )
  expect_warning(
    v_rm <- survey_reference_cov(IF,
      strata = strata, psu = psu,
      lonely_psu = "remove"
    )$cov,
    "single-PSU"
  )
  v_adj <- survey_reference_cov(IF,
    strata = strata, psu = psu,
    lonely_psu = "adjust"
  )$cov
  ## "adjust" adds the lonely stratum's contribution, so it is never smaller.
  expect_true(all(diag(v_adj) >= diag(v_rm)))
})

test_that("input length mismatches are rejected", {
  IF <- matrix(rnorm(20), ncol = 2)
  expect_error(survey_reference_cov(IF, weights = rep(1, 5)), "length")
  expect_error(survey_reference_cov(IF, strata = rep(1, 3)), "length")
  expect_error(survey_reference_cov(IF, weights = rep(-1, 10)), "non-negative")
})

test_that("PSU ids repeated across strata are not merged", {
  ## PSU 1 in stratum A and PSU 1 in stratum B are different clusters.
  set.seed(63)
  IF <- matrix(rnorm(40), ncol = 1)
  strata <- rep(c("A", "B"), each = 20)
  psu <- rep(rep(1:2, each = 10), 2)
  out <- survey_reference_cov(IF, strata = strata, psu = psu)
  expect_equal(out$n_psu, 4L)
  expect_equal(out$n_strata, 2L)
})
