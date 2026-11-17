make_reduced <- function(cohorts) {
  lapply(cohorts, function(d) {
    run_regression(d$y, d[, setdiff(names(d), "y"), drop = FALSE])
  })
}

test_that("partitioned mode uses the study block and the reference elsewhere", {
  s <- sim_setup(n_s = 600, n_ref = 3000, seed = 21)
  rms <- make_reduced(s$cohorts)

  out <- make_sigma_X_ls(
    n_cohorts = 3, reg_out = rms, ref_data = s$ref,
    cov_p_type = "partitioned"
  )

  expect_length(out$sigma_X_ls, 3L)
  expect_identical(colnames(out$sigma_X_ls[[1]]), c("x1", "x2", "x3"))

  ## Cohort 1 omits x3, so its x1/x2 block is its own covariance ...
  study_cov <- rms[[1]]$cov_X_s
  expect_equal(
    out$sigma_X_ls[[1]][c("x1", "x2"), c("x1", "x2")],
    study_cov[c("x1", "x2"), c("x1", "x2")]
  )
  ## ... and anything touching x3 comes from the reference.
  ref_cov <- cov(s$ref)
  expect_equal(out$sigma_X_ls[[1]]["x3", "x3"], ref_cov["x3", "x3"])
  expect_equal(out$sigma_X_ls[[1]]["x1", "x3"], ref_cov["x1", "x3"])
})

test_that("reference mode takes every block from the reference", {
  s <- sim_setup(n_s = 600, n_ref = 3000, seed = 22)
  rms <- make_reduced(s$cohorts)
  out <- make_sigma_X_ls(
    n_cohorts = 3, reg_out = rms, ref_data = s$ref,
    cov_p_type = "reference"
  )
  ref_cov <- cov(s$ref)[c("x1", "x2", "x3"), c("x1", "x2", "x3")]
  for (k in 1:3) expect_equal(out$sigma_X_ls[[k]], ref_cov)
})

test_that("all cohorts share one canonical variable ordering", {
  s <- sim_setup(n_s = 400, n_ref = 2000, seed = 23)
  rms <- make_reduced(s$cohorts)
  out <- make_sigma_X_ls(n_cohorts = 3, reg_out = rms, ref_data = s$ref)
  nms <- lapply(out$sigma_X_ls, colnames)
  expect_true(all(vapply(nms, identical, logical(1), nms[[1]])))
})

test_that("intercept columns are stripped from every covariance", {
  s <- sim_setup(n_s = 400, n_ref = 2000, seed = 24)
  ref <- s$ref
  ref$Intercept <- 1
  rms <- make_reduced(s$cohorts)
  out <- make_sigma_X_ls(n_cohorts = 3, reg_out = rms, ref_data = ref)
  ## A constant column has zero variance and would make every matrix singular.
  expect_false(any(grepl("Intercept", colnames(out$sigma_X_ls[[1]]))))
})

test_that("supplied ref_weights are used and not silently discarded", {
  ## The original code overwrote user-supplied ref_weights with NULL whenever
  ## `weighting` was FALSE, so the weighted covariance never reached the caller.
  s <- sim_setup(n_s = 400, n_ref = 2000, seed = 25)
  rms <- make_reduced(s$cohorts)
  w <- lapply(1:3, function(i) runif(nrow(s$ref), 0.5, 2))

  out <- make_sigma_X_ls(
    n_cohorts = 3, reg_out = rms, ref_data = s$ref, ref_weights = w,
    cov_p_type = "reference"
  )
  expect_false(any(vapply(out$ref_weights, is.null, logical(1))))
  expect_equal(out$ref_weights, w)

  X <- as.matrix(s$ref)
  expect_equal(
    out$sigma_X_ls[[1]],
    cov.wt(X, wt = w[[1]], center = TRUE)$cov[
      c("x1", "x2", "x3"), c("x1", "x2", "x3")
    ]
  )
})

test_that("weighted covariance uses the right reference matrix per cohort", {
  ## Regression test: the original passed the whole list of reference matrices
  ## to cov.wt() instead of the element for cohort k.
  s <- sim_setup(n_s = 400, n_ref = 1500, seed = 26)
  rms <- make_reduced(s$cohorts)
  w <- list(
    rep(1, nrow(s$ref)),
    c(rep(2, nrow(s$ref) / 2), rep(1, nrow(s$ref) / 2)),
    rep(3, nrow(s$ref))
  )
  out <- make_sigma_X_ls(
    n_cohorts = 3, reg_out = rms, ref_data = s$ref, ref_weights = w,
    cov_p_type = "reference"
  )
  ## Cohorts 1 and 3 have proportional weights, so identical covariances;
  ## cohort 2 must differ.
  expect_equal(out$sigma_X_ls[[1]], out$sigma_X_ls[[3]])
  expect_false(isTRUE(all.equal(out$sigma_X_ls[[1]], out$sigma_X_ls[[2]])))
})

test_that("a cohort covariate absent from the reference is an error", {
  s <- sim_setup(n_s = 300, n_ref = 1000, seed = 27)
  rms <- make_reduced(s$cohorts)
  ref_small <- s$ref[, c("x1", "x2")]
  expect_error(
    make_sigma_X_ls(n_cohorts = 3, reg_out = rms, ref_data = ref_small),
    "do not appear in the reference"
  )
})

test_that("missing reference input is an error", {
  expect_error(make_sigma_X_ls(n_cohorts = 2), "ref_data")
})

test_that("weighting without cohort data is an error", {
  s <- sim_setup(n_s = 300, n_ref = 1000, seed = 28)
  rms <- make_reduced(s$cohorts)
  expect_error(
    make_sigma_X_ls(
      n_cohorts = 3, reg_out = rms, ref_data = s$ref,
      weighting = TRUE, weighting.vars = "x1"
    ),
    "cohort_data"
  )
})

test_that("a single reference data frame is recycled across cohorts", {
  s <- sim_setup(n_s = 300, n_ref = 1000, seed = 29)
  rms <- make_reduced(s$cohorts)
  out <- make_sigma_X_ls(n_cohorts = 3, reg_out = rms, ref_data = s$ref)
  expect_length(out$ref_pop_dummy, 3L)
  expect_equal(out$ref_pop_dummy[[1]], out$ref_pop_dummy[[3]])
})
