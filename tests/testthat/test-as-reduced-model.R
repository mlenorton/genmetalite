test_that("as_reduced_model matches run_regression on the same data", {
  d <- sim_cohort(800, seed = 101)
  fit <- lm(y ~ x1 + x2 + x3, data = d)

  a <- as_reduced_model(fit)
  b <- run_regression(d$y, d[, c("x1", "x2", "x3")], dgp = "linear")

  expect_equal(a$theta_hat_s, b$theta_hat_s, tolerance = 1e-10)
  expect_equal(unname(a$cov_theta_hat_s), unname(b$cov_theta_hat_s),
    tolerance = 1e-10
  )
  expect_equal(a$cov_X_s, b$cov_X_s, tolerance = 1e-10)
})

test_that("the intercept is renamed and optionally dropped", {
  d <- sim_cohort(300, seed = 102)
  fit <- lm(y ~ x1 + x2, data = d)

  with_int <- as_reduced_model(fit)
  expect_identical(names(with_int$theta_hat_s)[1], "Intercept")
  expect_false("Intercept" %in% colnames(with_int$cov_X_s))

  without <- as_reduced_model(fit, keep_intercept = FALSE)
  expect_false("Intercept" %in% names(without$theta_hat_s))
  expect_identical(dim(without$cov_theta_hat_s), c(2L, 2L))
})

test_that("glm fits are accepted", {
  set.seed(103)
  d <- sim_pop(600)
  d$y <- rbinom(600, 1, plogis(as.matrix(d) %*% true_beta))
  fit <- glm(y ~ x1 + x2 + x3, family = binomial, data = d)
  out <- as_reduced_model(fit)
  expect_equal(out$dgp, "binomial")
  expect_setequal(names(out$theta_hat_s), c("Intercept", "x1", "x2", "x3"))
})

test_that("non lm/glm objects are rejected with guidance", {
  expect_error(as_reduced_model(structure(list(), class = "nope")), "lm")
})

test_that("the legacy names still work", {
  d <- sim_cohort(200, seed = 104)
  fit <- lm(y ~ x1, data = d)
  expect_equal(
    convert_reg_out_to_reduced_model_format(fit)$theta_hat_s,
    as_reduced_model(fit)$theta_hat_s
  )
})

test_that("a genmetalite fit can be fed back in as a reduced model", {
  s <- sim_setup(n_s = 800, n_ref = 4000, seed = 105)
  fit <- calc_linear_beta(cohort_data = s$cohorts, ref_data = s$ref)
  cov_ref <- cov(s$ref)

  rm2 <- as_reduced_model_genmetalite(fit, cov_X = cov_ref)
  expect_s3_class(rm2, "reduced_model")
  expect_setequal(names(rm2$theta_hat_s), c("x1", "x2", "x3"))
  expect_identical(dim(rm2$cov_X_s), c(3L, 3L))

  expect_error(as_reduced_model_genmetalite(fit), "cov_X. is required")
  expect_error(
    as_reduced_model_genmetalite(fit, cov_X = cov_ref[1:2, 1:2]),
    "missing these variables"
  )
})

test_that("as_reduced_model output works inside calc_linear_beta", {
  s <- sim_setup(n_s = 1200, n_ref = 6000, seed = 106)
  rms <- lapply(s$cohorts, function(d) {
    vars <- setdiff(names(d), "y")
    as_reduced_model(lm(
      stats::as.formula(paste("y ~", paste(vars, collapse = " + "))),
      data = d
    ))
  })
  fit <- suppressWarnings(
    calc_linear_beta(reduced_models = rms, ref_data = s$ref)
  )
  expect_equal(fit$linear_est$beta, true_beta, tolerance = 0.1)
})
