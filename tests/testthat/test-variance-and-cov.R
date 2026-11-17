test_that("compute_var returns the expected estimators", {
  d <- sim_cohort(500, seed = 11)
  fit <- lm(y ~ x1 + x2 + x3, data = d)

  v_std <- compute_var(fit, "standard")
  v_rob <- compute_var(fit, "robust")

  expect_equal(v_std, vcov(fit))
  expect_equal(v_rob, sandwich::sandwich(fit))
  expect_identical(dim(v_rob), dim(v_std))
  ## Under homoskedasticity the two should be close but not identical.
  expect_false(isTRUE(all.equal(v_rob, v_std)))
  expect_equal(max(abs(v_rob - v_std)) / max(abs(v_std)), 0, tolerance = 0.3)
})

test_that("compute_var('best') falls back gracefully and stays symmetric", {
  d <- sim_cohort(400, seed = 12)
  d$rare <- c(rep(1, 5), rep(0, nrow(d) - 5))
  fit <- lm(y ~ x1 + rare, data = d)

  v <- compute_var(fit, "best", min_obs = 100)
  expect_true(isSymmetric(unname(v)))
  ## The rare indicator has 5 observations, so it takes the model-based value.
  expect_equal(v["rare", "rare"], vcov(fit)["rare", "rare"])
  ## x1 is continuous, so it keeps the sandwich value.
  expect_equal(v["x1", "x1"], sandwich::sandwich(fit)["x1", "x1"])
})

test_that("compute_var('best') does not error on a non lm/glm object", {
  obj <- structure(list(), class = "something_else")
  expect_error(compute_var(obj, "best"))
})

test_that("cov_factors expands factors and drops the intercept", {
  d <- sim_pop_factor(500, seed = 13)
  out <- cov_factors(d)

  expect_setequal(colnames(out$X_dummy_mat), c("age", "edu2", "edu3"))
  expect_false("(Intercept)" %in% colnames(out$X_dummy_mat))
  expect_identical(dim(out$cov_X_s), c(3L, 3L))
  expect_equal(out$cov_X_s, cov(out$X_dummy_mat))
})

test_that("cov_factors drops unused factor levels", {
  d <- data.frame(
    a = rnorm(20),
    g = factor(rep(c("x", "y"), 10), levels = c("x", "y", "z"))
  )
  out <- cov_factors(d)
  ## Level "z" is empty; keeping it would add an all-zero, singular column.
  expect_setequal(colnames(out$X_dummy_mat), c("a", "gy"))
})

test_that("cov_factors refuses incomplete rows rather than silently dropping", {
  d <- data.frame(a = c(1, 2, NA, 4), b = c(1, 2, 3, 4))
  expect_error(cov_factors(d), "missing values")
})

test_that("cov_factors honours weights", {
  set.seed(14)
  d <- data.frame(a = rnorm(300), b = rnorm(300))
  w <- runif(300, 0.5, 2)
  out <- cov_factors(d, wt = w)
  expect_equal(out$cov_X_s, cov.wt(as.matrix(d), wt = w)$cov)
})
