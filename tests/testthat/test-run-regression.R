test_that("run_regression recovers the coefficients of a linear model", {
  d <- sim_cohort(4000, seed = 1)
  rm1 <- run_regression(d$y, d[, c("x1", "x2", "x3")], dgp = "linear")

  expect_s3_class(rm1, "reduced_model")
  expect_equal(
    unname(rm1$theta_hat_s[c("x1", "x2", "x3")]),
    unname(true_beta),
    tolerance = 0.06
  )
  expect_equal(unname(rm1$theta_hat_s[["Intercept"]]), true_intercept,
    tolerance = 0.06
  )
})

test_that("the intercept is named consistently and kept out of cov_X_s", {
  d <- sim_cohort(500, seed = 2)
  rm1 <- run_regression(d$y, d[, c("x1", "x2")], dgp = "linear")

  expect_identical(names(rm1$theta_hat_s)[1], "Intercept")
  expect_identical(rownames(rm1$cov_theta_hat_s)[1], "Intercept")
  expect_false(any(grepl("Intercept", colnames(rm1$cov_X_s))))
  expect_identical(colnames(rm1$cov_X_s), c("x1", "x2"))
})

test_that("coefficient names have no matrix-argument prefix", {
  d <- sim_cohort(300, seed = 3)
  rm1 <- run_regression(d$y, d[, c("x1", "x2")], dgp = "linear")
  expect_false(any(grepl("^X_mat", names(rm1$theta_hat_s))))
  expect_identical(names(rm1$theta_hat_s), rownames(rm1$cov_theta_hat_s))
})

test_that("factors are expanded with model.matrix naming", {
  d <- sim_pop_factor(600, seed = 4)
  d$y <- rnorm(600)
  rm1 <- run_regression(d$y, d[, c("age", "edu")], dgp = "linear")
  expect_setequal(colnames(rm1$cov_X_s), c("age", "edu2", "edu3"))
  expect_true(all(c("edu2", "edu3") %in% names(rm1$theta_hat_s)))
})

test_that("converged is TRUE for lm fits rather than NULL", {
  ## lm objects carry no $converged element; the original code returned NULL
  ## here, which then propagated as an unusable value.
  d <- sim_cohort(200, seed = 5)
  rm1 <- run_regression(d$y, d[, c("x1", "x2")], dgp = "linear")
  expect_true(isTRUE(rm1$converged))
})

test_that("logistic and poisson fits work", {
  set.seed(6)
  d <- sim_pop(3000)
  eta <- -0.5 + as.matrix(d) %*% true_beta
  yb <- rbinom(3000, 1, plogis(eta))
  rm_b <- run_regression(yb, d, dgp = "logistic")
  expect_equal(unname(rm_b$theta_hat_s[c("x1", "x2", "x3")]),
    unname(true_beta),
    tolerance = 0.2
  )

  yp <- rpois(3000, exp(0.2 + as.matrix(d) %*% (true_beta / 2)))
  rm_p <- run_regression(yp, d, dgp = "poisson")
  expect_equal(unname(rm_p$theta_hat_s[c("x1", "x2", "x3")]),
    unname(true_beta / 2),
    tolerance = 0.15
  )
})

test_that("dimension mismatches and bad selectors are rejected", {
  d <- sim_cohort(100, seed = 7)
  X <- d[, c("x1", "x2", "x3")]

  expect_error(run_regression(d$y[1:50], X), "length\\(y\\)")
  expect_error(run_regression(d$y, X, inds_theta_s = c(1, 0)), "length 2")
  expect_error(run_regression(d$y, X, inds_theta_s = c(1, 2, 3)), "0s and 1s")
  expect_error(run_regression(d$y, X, inds_theta_s = c(0, 0, 0)), "no columns")
  expect_error(run_regression(d$y, X, dgp = "nonsense"))
})

test_that("inds_theta_s selects columns", {
  d <- sim_cohort(400, seed = 8)
  X <- d[, c("x1", "x2", "x3")]
  rm1 <- run_regression(d$y, X, inds_theta_s = c(1, 0, 1))
  expect_setequal(names(rm1$theta_hat_s), c("Intercept", "x1", "x3"))
})

test_that("a single selected column does not collapse to a vector", {
  d <- sim_cohort(400, seed = 9)
  X <- d[, c("x1", "x2", "x3")]
  rm1 <- run_regression(d$y, X, inds_theta_s = c(1, 0, 0))
  expect_setequal(names(rm1$theta_hat_s), c("Intercept", "x1"))
  expect_equal(dim(rm1$cov_X_s), c(1L, 1L))
})
