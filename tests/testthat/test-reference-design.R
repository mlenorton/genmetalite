make_ref_pieces <- function(n_str = 6, n_psu = 5, m = 40, seed = 201) {
  set.seed(seed)
  N <- n_str * n_psu * m
  psu <- rep(seq_len(n_str * n_psu), each = m)
  strata <- rep(seq_len(n_str), each = n_psu * m)
  u <- rnorm(n_str * n_psu, sd = 0.9)[psu]
  x1 <- u + rnorm(N)
  x2 <- 0.6 * x1 + rnorm(N, sd = 0.8)
  x3 <- 0.4 * x2 + rnorm(N)
  list(
    df = data.frame(x1 = x1, x2 = x2, x3 = x3),
    psu = psu, strata = strata, w = exp(rnorm(N, 0, 0.4)) * 100
  )
}

svy_cohorts <- function(n = 900, seed = 202) {
  set.seed(seed)
  list(sim_cohort(n, "x3"), sim_cohort(n, "x2"), sim_cohort(n, "x1"))
}

test_that("ref_design validates its inputs", {
  d <- data.frame(a = rnorm(10))
  expect_s3_class(ref_design(d), "genmetalite_ref_design")
  expect_error(ref_design(d, weights = rep(1, 5)), "length 5")
  expect_error(ref_design(d, strata = rep(1, 3)), "length 3")
  expect_error(ref_design(d, weights = c(rep(1, 9), -1)), "non-negative")
  expect_error(ref_design(d, weights = rep(0, 10)), "sum to zero")
  expect_error(ref_design(d[0, , drop = FALSE]), "no rows")
})

test_that("as_ref_design passes through data frames and designs", {
  d <- data.frame(a = rnorm(10))
  expect_identical(as_ref_design(d)$data, d)
  expect_null(as_ref_design(d)$weights)

  rd <- ref_design(d, weights = runif(10, 1, 2))
  expect_identical(as_ref_design(rd), rd)

  expect_error(as_ref_design(1:10), "must be a data frame")
  expect_error(
    as_ref_design(structure(list(), class = "svyrep.design")),
    "Replicate-weight"
  )
})

test_that("a survey.design object is read correctly", {
  skip_if_not_installed("survey")
  r <- make_ref_pieces()
  des <- survey::svydesign(
    id = ~psu, strata = ~strata, weights = ~w,
    data = cbind(r$df, psu = r$psu, strata = r$strata, w = r$w),
    nest = TRUE
  )
  rd <- as_ref_design(des)

  expect_s3_class(rd, "genmetalite_ref_design")
  expect_equal(rd$weights, r$w, tolerance = 1e-8)
  expect_identical(length(unique(rd$strata)), 6L)
  expect_identical(length(unique(paste(rd$strata, rd$psu))), 30L)
  expect_output(print(rd), "reference design")
})

test_that("an unclustered survey design reports no PSUs", {
  skip_if_not_installed("survey")
  d <- data.frame(x1 = rnorm(200), w = runif(200, 1, 5))
  des <- survey::svydesign(id = ~1, weights = ~w, data = d)
  rd <- as_ref_design(des)
  ## id = ~1 means "no clustering"; svydesign fills cluster with row numbers,
  ## which carries no information and must not be mistaken for real PSUs.
  expect_null(rd$psu)
  expect_null(rd$strata)
  expect_equal(rd$weights, d$w, tolerance = 1e-8)
})

test_that("ref_design and the equivalent survey.design give identical results", {
  skip_if_not_installed("survey")
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()
  des <- survey::svydesign(
    id = ~psu, strata = ~strata, weights = ~w,
    data = cbind(r$df, psu = r$psu, strata = r$strata, w = r$w),
    nest = TRUE
  )
  a <- calc_linear_beta(cohort_data = cohorts, ref_data = des)
  b <- calc_linear_beta(
    cohort_data = cohorts,
    ref_data = ref_design(r$df, weights = r$w, strata = r$strata, psu = r$psu)
  )
  expect_equal(a$linear_est$beta, b$linear_est$beta)
  expect_equal(a$linear_est$var, b$linear_est$var)
  expect_identical(a$linear_est$design_df, b$linear_est$design_df)
})

test_that("bookkeeping columns in a design do not become covariates", {
  skip_if_not_installed("survey")
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()
  des <- survey::svydesign(
    id = ~psu, strata = ~strata, weights = ~w,
    data = cbind(r$df, psu = r$psu, strata = r$strata, w = r$w),
    nest = TRUE
  )
  ## The design's `variables` slot carries psu, strata and w alongside the
  ## covariates; only the covariates the studies use may enter the model.
  fit <- calc_linear_beta(cohort_data = cohorts, ref_data = des)
  expect_setequal(names(fit$linear_est$beta), c("x1", "x2", "x3"))
  expect_setequal(colnames(fit$sigma_X_ls[[1]]), c("x1", "x2", "x3"))
})

test_that("the design changes the variance but not the point estimate", {
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()

  weighted_only <- calc_linear_beta(
    cohort_data = cohorts, ref_data = ref_design(r$df, weights = r$w)
  )
  with_design <- calc_linear_beta(
    cohort_data = cohorts,
    ref_data = ref_design(r$df, weights = r$w, strata = r$strata, psu = r$psu)
  )

  ## Sigma depends only on the weights, so beta is identical ...
  expect_equal(weighted_only$linear_est$beta, with_design$linear_est$beta)
  ## ... but PSU-level clustering inflates the reference term.
  expect_true(all(
    diag(with_design$linear_est$var) > diag(weighted_only$linear_est$var)
  ))
  expect_null(weighted_only$linear_est$design_df)
  expect_identical(with_design$linear_est$design_df, 30L - 6L)
})

test_that("weights change the covariance and so the estimate", {
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()
  a <- calc_linear_beta(cohort_data = cohorts, ref_data = r$df)
  b <- calc_linear_beta(
    cohort_data = cohorts, ref_data = ref_design(r$df, weights = r$w)
  )
  expect_false(isTRUE(all.equal(a$linear_est$beta, b$linear_est$beta)))
})

test_that("loose survey_strata/survey_psu still work for a plain data frame", {
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()
  a <- calc_linear_beta(
    cohort_data = cohorts, ref_data = r$df,
    survey_strata = r$strata, survey_psu = r$psu
  )
  b <- calc_linear_beta(
    cohort_data = cohorts,
    ref_data = ref_design(r$df, strata = r$strata, psu = r$psu)
  )
  expect_equal(a$linear_est$var, b$linear_est$var)
})

test_that("a design object supersedes loose design arguments, with a warning", {
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()
  expect_warning(
    calc_linear_beta(
      cohort_data = cohorts,
      ref_data = ref_design(r$df, strata = r$strata, psu = r$psu),
      survey_strata = rep(1, nrow(r$df))
    ),
    "already carries a sampling design"
  )
})

test_that("design vectors of the wrong length are caught", {
  r <- make_ref_pieces()
  cohorts <- svy_cohorts()
  expect_error(
    calc_linear_beta(
      cohort_data = cohorts, ref_data = r$df,
      survey_psu = rep(1:5, length.out = 17)
    ),
    "has length 17 but the reference sample has"
  )
})
