## Cohorts that measured a factor but have an empty cell in one of its levels.

sim_factor_setup <- function(n_ref = 6000, n_s = 2000, seed = 121) {
  set.seed(seed)
  b <- c(age = 0.8, edu2 = -0.4, edu3 = 0.6)
  make <- function(n, probs) {
    d <- data.frame(
      age = rnorm(n),
      edu = factor(sample(1:3, n, replace = TRUE, prob = probs),
        levels = 1:3
      )
    )
    X <- model.matrix(~., d)[, -1, drop = FALSE]
    d$y <- as.numeric(1 + X %*% b[colnames(X)] + rnorm(n))
    d
  }
  list(
    beta = b,
    ref = make(n_ref, c(.5, .3, .2))[, c("age", "edu")],
    full = make(n_s, c(.5, .3, .2)),
    no_edu3 = make(n_s, c(.6, .4, 0))
  )
}

test_that("check_empty_levels finds a level with no observations", {
  ref <- data.frame(edu = factor(c("1", "2", "3")))
  coh <- data.frame(edu = factor(c("1", "2"), levels = c("1", "2", "3")))

  expect_warning(out <- check_empty_levels(ref, coh), "level '3'")
  expect_identical(nrow(out), 1L)
  expect_identical(out$variable, "edu")
  expect_identical(out$level, "3")
  expect_identical(out$status, "unobserved")
})

test_that("check_empty_levels is silent when every level is observed", {
  ref <- data.frame(edu = factor(c("1", "2", "3")), x = 1:3)
  coh <- data.frame(edu = factor(c("1", "2", "3")), x = 4:6)
  expect_silent(out <- check_empty_levels(ref, coh))
  expect_identical(nrow(out), 0L)
})

test_that("an unobserved reference baseline is an error", {
  ## Dropping it would promote another level to baseline, so the cohort's
  ## dummy columns would encode different contrasts under the same names.
  ref <- data.frame(edu = factor(c("1", "2", "3")))
  coh <- data.frame(edu = factor(c("2", "3"), levels = c("1", "2", "3")))
  expect_error(check_empty_levels(ref, coh), "baseline")
  expect_error(check_empty_levels(ref, coh), "Relevel")
})

test_that("action = 'none' collects without signalling", {
  ref <- data.frame(edu = factor(c("1", "2", "3")))
  coh <- data.frame(edu = factor(c("1", "2"), levels = c("1", "2", "3")))
  expect_silent(out <- check_empty_levels(ref, coh, action = "none"))
  expect_identical(nrow(out), 1L)
})

test_that("character columns are handled like factors", {
  ref <- data.frame(g = c("a", "b", "c"), stringsAsFactors = FALSE)
  coh <- data.frame(g = c("a", "b"), stringsAsFactors = FALSE)
  expect_warning(out <- check_empty_levels(ref, coh), "level 'c'")
  expect_identical(out$status, "unobserved")
})

test_that("calc_linear_beta warns and records the empty level", {
  s <- sim_factor_setup()
  cohorts <- list(s$no_edu3, s$full)

  expect_warning(
    fit <- calc_linear_beta(cohort_data = cohorts, ref_data = s$ref),
    "no observations in edu level '3'"
  )

  expect_s3_class(fit$empty_levels, "data.frame")
  expect_identical(nrow(fit$empty_levels), 1L)
  expect_identical(fit$empty_levels$cohort, 1L)
  expect_identical(fit$empty_levels$variable, "edu")
  expect_identical(fit$empty_levels$level, "3")

  ## Cohort 1 contributes no edu3 coefficient, but the pooled model still has
  ## one, borrowed from cohort 2 and the reference covariance.
  expect_false("edu3" %in% names(fit$reduced_models[[1]]$theta_hat_s))
  expect_setequal(names(fit$linear_est$beta), c("age", "edu2", "edu3"))
  expect_equal(fit$linear_est$beta[names(s$beta)], s$beta, tolerance = 0.12)
})

test_that("empty_levels has zero rows when nothing is wrong", {
  s <- sim_factor_setup(seed = 122)
  fit <- calc_linear_beta(
    cohort_data = list(s$full, s$full), ref_data = s$ref
  )
  expect_identical(nrow(fit$empty_levels), 0L)
})

test_that("baselines are harmonised so dummy names mean the same thing", {
  ## Cohort 1 observes only levels 2 and 3. Left alone, droplevels() would make
  ## "2" its baseline and its `edu3` column would mean "3 vs 2" while the
  ## reference's means "3 vs 1" -- matched by name, silently different.
  s <- sim_factor_setup(seed = 123)
  coh1 <- s$full[s$full$edu != "1", ]
  coh1$edu <- factor(as.character(coh1$edu), levels = c("2", "3"))

  expect_error(
    calc_linear_beta(cohort_data = list(coh1, s$full), ref_data = s$ref),
    "baseline"
  )
})

test_that("harmonise_levels keeps cohort-only levels rather than dropping them", {
  ref <- data.frame(g = factor(c("a", "b")))
  coh <- data.frame(g = factor(c("a", "c")))
  out <- genmetalite:::harmonise_levels(coh, ref)
  expect_identical(levels(out$g), c("a", "b", "c"))
  expect_identical(as.character(out$g), c("a", "c"))
})
