test_that("entropy_balance reproduces the target moments exactly", {
  set.seed(301)
  Z <- cbind(a = rnorm(2000), b = rnorm(2000), c = rnorm(2000))
  target <- c(a = 0.5, b = -0.2, c = 0.1)
  out <- entropy_balance(Z, target)

  expect_true(out$converged)
  expect_true(all(out$weights > 0))
  achieved <- colSums(out$weights * Z) / sum(out$weights)
  expect_equal(achieved, target, tolerance = 1e-7, ignore_attr = TRUE)
})

test_that("weights sum to the base-weight total", {
  set.seed(302)
  Z <- cbind(a = rnorm(1000))
  w0 <- runif(1000, 1, 50)
  out <- entropy_balance(Z, target = c(a = 0.3), base_weight = w0)
  expect_equal(sum(out$weights), sum(w0), tolerance = 1e-8)

  ## With no base weights, the total is the row count.
  out2 <- entropy_balance(Z, target = c(a = 0.3))
  expect_equal(sum(out2$weights), 1000, tolerance = 1e-8)
})

test_that("an already-matching target leaves the base weights alone", {
  set.seed(303)
  Z <- cbind(a = rnorm(1500), b = rnorm(1500))
  w0 <- runif(1500, 1, 10)
  target <- colSums(w0 * Z) / sum(w0)
  out <- entropy_balance(Z, target, base_weight = w0)
  expect_true(out$converged)
  expect_equal(out$weights, w0, tolerance = 1e-6)
})

test_that("base weights tilt the population, not the raw sample", {
  ## Half the rows carry 10x the sampling weight. The solution must match the
  ## target in the weighted (population) sense.
  set.seed(304)
  n <- 2000
  Z <- cbind(a = rnorm(n))
  w0 <- c(rep(1, n / 2), rep(10, n / 2))
  out <- entropy_balance(Z, target = c(a = 0.4), base_weight = w0)
  expect_equal(
    weighted.mean(Z[, "a"], out$weights), 0.4,
    tolerance = 1e-7
  )
  ## The unweighted mean of Z under these weights is NOT the target, which is
  ## the point: the base weights are respected.
  expect_false(isTRUE(all.equal(mean(Z[, "a"] * out$weights), 0.4)))
})

test_that("an unreachable target is reported rather than faked", {
  set.seed(305)
  Z <- cbind(a = rnorm(500))
  ## Far outside the convex hull of the data: no positive weights can get there.
  out <- suppressWarnings(entropy_balance(Z, target = c(a = 25)))
  expect_false(out$converged)
  expect_gt(out$max_imbalance, 1e-6)
})

test_that("entropy_balance validates its inputs", {
  Z <- cbind(a = rnorm(50))
  expect_error(entropy_balance(Z, target = c(1, 2)), "length 2")
  expect_error(entropy_balance(Z, target = 0, base_weight = rep(1, 3)), "length 3")
  expect_error(
    entropy_balance(Z, target = 0, base_weight = rep(-1, 50)), "non-negative"
  )
  Zna <- Z
  Zna[1] <- NA
  expect_error(entropy_balance(Zna, target = 0), "missing values")
})

test_that("ebal_weights matches the requested moments with the internal engine", {
  set.seed(306)
  ref <- data.frame(age = rnorm(3000), bmi = rnorm(3000))
  coh <- data.frame(age = rnorm(700, 0.7), bmi = rnorm(700, 0.3))

  out <- ebal_weights(coh, ref)
  expect_true(out$converged)
  expect_length(out$weights, nrow(ref))
  expect_equal(weighted.mean(ref$age, out$weights), mean(coh$age),
    tolerance = 1e-6
  )
  expect_equal(weighted.mean(ref$bmi, out$weights), mean(coh$bmi),
    tolerance = 1e-6
  )
})

test_that("ebal_weights balances factor levels", {
  set.seed(307)
  ref <- data.frame(g = factor(sample(c("a", "b", "c"), 3000, TRUE)))
  coh <- data.frame(
    g = factor(sample(c("a", "b", "c"), 800, TRUE, prob = c(.6, .3, .1)))
  )
  out <- ebal_weights(coh, ref)
  expect_true(out$converged)
  for (lv in c("b", "c")) {
    expect_equal(
      weighted.mean(ref$g == lv, out$weights), mean(coh$g == lv),
      tolerance = 1e-6
    )
  }
})

test_that("ebal_weights carries survey weights through as base weights", {
  set.seed(308)
  ref <- data.frame(age = rnorm(2500), bmi = rnorm(2500))
  coh <- data.frame(age = rnorm(600, 0.5), bmi = rnorm(600))
  sw <- exp(rnorm(2500, 0, 0.5)) * 100

  out <- ebal_weights(coh, ref, base.weight = sw)
  expect_true(out$converged)
  expect_equal(weighted.mean(ref$age, out$weights), mean(coh$age),
    tolerance = 1e-6
  )
  ## The population total the survey weights represent is preserved.
  expect_equal(sum(out$weights), sum(sw), tolerance = 1e-6)
  ## Rows with a big sampling weight keep a big weight.
  expect_gt(cor(sw, out$weights), 0.5)
})

test_that("ebal_weights rejects malformed base weights", {
  ref <- data.frame(a = rnorm(100))
  coh <- data.frame(a = rnorm(50))
  expect_error(ebal_weights(coh, ref, base.weight = rep(1, 7)), "length 7")
  expect_error(
    ebal_weights(coh, ref, base.weight = rep(-1, 100)), "non-negative"
  )
})

test_that("trim requires the ebal engine", {
  ref <- data.frame(a = rnorm(100))
  coh <- data.frame(a = rnorm(50))
  expect_error(ebal_weights(coh, ref, trim = TRUE), "engine = 'ebal'")
})

test_that("the internal and ebal engines agree when ebal is available", {
  skip_if_not_installed("ebal")
  set.seed(309)
  ref <- data.frame(age = rnorm(1500), bmi = rnorm(1500))
  coh <- data.frame(age = rnorm(400, 0.4), bmi = rnorm(400, 0.2))
  a <- ebal_weights(coh, ref, engine = "internal")
  b <- ebal_weights(coh, ref, engine = "ebal")
  ## Same strictly convex problem, so the same solution up to scale.
  expect_equal(
    a$weights / sum(a$weights), b$weights / sum(b$weights),
    tolerance = 1e-5
  )
})

test_that("ipw_weights accepts survey weights as base weights", {
  set.seed(310)
  ref <- data.frame(age = rnorm(3000), bmi = rnorm(3000))
  coh <- data.frame(age = rnorm(700, 0.8), bmi = rnorm(700))
  sw <- exp(rnorm(3000, 0, 0.4)) * 50

  out <- ipw_weights(coh, ref, base.weight = sw, weighting.corr = FALSE)
  expect_length(out$weights, nrow(ref))
  expect_equal(sum(out$weights), sum(sw), tolerance = 1e-6)
  ## Moves the weighted reference mean towards the cohort.
  expect_lt(
    abs(weighted.mean(ref$age, out$weights) - mean(coh$age)),
    abs(weighted.mean(ref$age, sw) - mean(coh$age)) / 2
  )
  expect_error(ipw_weights(coh, ref, base.weight = rep(1, 4)), "length 4")
})
