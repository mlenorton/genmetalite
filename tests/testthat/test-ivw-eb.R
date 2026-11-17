test_that("IVW_estimate reproduces a hand-computed inverse-variance average", {
  named_diag <- function(v, nms) {
    m <- diag(v, nrow = length(v))
    dimnames(m) <- list(nms, nms)
    m
  }
  nms <- c("Intercept", "x1")
  rm1 <- list(
    theta_hat_s = c(Intercept = 0, x1 = 1.0),
    cov_theta_hat_s = named_diag(c(1, 0.01), nms)
  )
  rm2 <- list(
    theta_hat_s = c(Intercept = 0, x1 = 2.0),
    cov_theta_hat_s = named_diag(c(1, 0.04), nms)
  )
  out <- IVW_estimate(list(rm1, rm2))

  w <- c(1 / 0.01, 1 / 0.04)
  expect_equal(unname(out$beta[["x1"]]), sum(c(1, 2) * w) / sum(w))
  expect_equal(unname(out$var[["x1"]]), 1 / sum(w))
  expect_equal(unname(out$n_studies[["x1"]]), 2L)
})

test_that("IVW_estimate handles coefficients present in only some studies", {
  s <- sim_setup(n_s = 500, n_ref = 2000, seed = 51)
  rms <- lapply(s$cohorts, function(d) {
    run_regression(d$y, d[, setdiff(names(d), "y"), drop = FALSE])
  })
  out <- IVW_estimate(rms)
  expect_setequal(names(out$beta), c("Intercept", "x1", "x2", "x3"))
  ## Each of x1, x2, x3 is omitted by exactly one of the three studies.
  expect_equal(unname(out$n_studies[c("x1", "x2", "x3")]), c(2L, 2L, 2L))
  expect_equal(unname(out$n_studies[["Intercept"]]), 3L)
  expect_false(anyNA(out$beta))
})

test_that("IVW_estimate rejects an empty input", {
  expect_error(IVW_estimate(list()), "empty")
})

test_that("EB_estimate lies between the two inputs", {
  v <- diag(c(0.01, 0.01))
  dimnames(v) <- list(c("x1", "x2"), c("x1", "x2"))
  linear_est <- list(beta = c(x1 = 1.0, x2 = -0.5), var = v)
  ivw <- list(beta = c(x1 = 1.4, x2 = -0.5))
  eb <- EB_estimate(linear_est, ivw)

  expect_true(eb$beta[["x1"]] >= min(1.0, 1.4))
  expect_true(eb$beta[["x1"]] <= max(1.0, 1.4))
  ## Large apparent bias -> weight goes to the mapped estimate.
  expect_gt(eb$K[["x1"]], 0.9)
  ## No apparent bias -> weight goes to the IVW estimate.
  expect_equal(eb$K[["x2"]], 0, tolerance = 1e-12)
  expect_equal(eb$beta[["x2"]], -0.5)
})

test_that("EB_estimate errors when a coefficient is missing from IVW", {
  linear_est <- list(
    beta = c(x1 = 1, x2 = 2),
    var = diag(2, nrow = 2, ncol = 2)
  )
  dimnames(linear_est$var) <- list(c("x1", "x2"), c("x1", "x2"))
  expect_error(
    EB_estimate(linear_est, list(beta = c(x1 = 1))),
    "missing coefficient"
  )
})
