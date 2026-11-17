## Shared simulation helpers for the test suite.

## Correlated covariates, so that omitting one really does bias the others.
sim_pop <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  x1 <- stats::rnorm(n)
  x2 <- 0.6 * x1 + stats::rnorm(n, sd = 0.8)
  x3 <- 0.4 * x2 + stats::rnorm(n)
  data.frame(x1 = x1, x2 = x2, x3 = x3)
}

true_beta <- c(x1 = 1, x2 = -0.5, x3 = 0.3)
true_intercept <- 2

sim_cohort <- function(n, drop = character(0), seed = NULL, sd = 1) {
  d <- sim_pop(n, seed)
  d$y <- as.numeric(
    true_intercept + as.matrix(d) %*% true_beta + stats::rnorm(n, sd = sd)
  )
  d[, setdiff(names(d), drop), drop = FALSE]
}

sim_pop_factor <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  age <- stats::rnorm(n)
  edu <- factor(sample(1:3, n, replace = TRUE, prob = c(.5, .3, .2)))
  data.frame(age = age, edu = edu)
}

## A standard three-study setup where each study omits one covariate.
sim_setup <- function(n_s = 1500, n_ref = 8000, seed = 2026) {
  set.seed(seed)
  list(
    ref = sim_pop(n_ref),
    cohorts = list(
      sim_cohort(n_s, "x3"),
      sim_cohort(n_s, "x2"),
      sim_cohort(n_s, "x1")
    )
  )
}
