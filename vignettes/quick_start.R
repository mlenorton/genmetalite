install.packages("genmetalite", repos = NULL, type = "source")

# Case 1: You have raw data from cohorts
set.seed(2026)
n_ref <- 4000
make_pop <- function(n) {
  x1 <- rnorm(n)
  x2 <- 0.6 * x1 + rnorm(n, sd = 0.8)
  x3 <- rnorm(n)
  data.frame(x1 = x1, x2 = x2, x3 = x3)
}
beta <- c(x1 = 1, x2 = -0.5, x3 = 0.3)
ref <- make_pop(n_ref)

mk_cohort <- function(n, drop) {
  d <- make_pop(n)
  d$y <- 2 + as.matrix(d) %*% beta + rnorm(n)
  d[, setdiff(names(d), drop), drop = FALSE]
}
cohorts <- list(mk_cohort(800, "x3"), mk_cohort(800, "x2"))

fit1a <- calc_linear_beta(cohort_data = cohorts, ref_data = ref, dgp = "linear", cov_p_type = "reference") 
round(fit1a$linear_est$beta, 3)
round(fit1a$IVW_est$beta, 3) # biased: each study omits a correlated covariate

# Option 1: re-weight reference pop to match each study
fit1b <- calc_linear_beta(cohort_data = cohorts, ref_data = ref, 
                          dgp = "linear", cov_p_type = "reference",
                          weighting.vars = list(c("x1", "x2"), c("x1","x3"))) 
round(fit1b$linear_est$beta, 3)

# Option 2: re-weight reference pop to match each study and partition the variance matrix
fit1c <- calc_linear_beta(cohort_data = cohorts, ref_data = ref, 
                          dgp = "linear", cov_p_type = "partitioned",
                          weighting.vars = list(c("x1", "x2"), c("x1","x3"))) 
round(fit1c$linear_est$beta, 3)


# Case 2: You prefer to fit the models yourself ahead of time and feed in summary statistics
data1 <- mk_cohort(800, "x3")
data2 <- mk_cohort(800, "x2")
model1 <- lm(y ~ x1 + x2, data1)
model2 <- lm(y ~ x1 + x3, data2)

red_mod1 <- convert_reg_out_to_reduced_model_format(model1)
red_mod2 <- convert_reg_out_to_reduced_model_format(model2)
fit2 <- calc_linear_beta(reduced_models = list(red_mod1, red_mod2), ref_data = ref)
round(fit2$linear_est$beta, 3)

# Case 3: You have summary statistics only (perhaps in a .csv)
# ---- build one sublist per study --------------------------------------------
build_study <- function(df) {
  # Preserve CSV row order within the study
  theta_hat_s <- df$RR
  # Name the elements by coeff_name for traceability (like named vector)
  names(theta_hat_s) <- df$coeff_name
  
  se_s <- df$se
  # Diagonal covariance approximation: variances = se^2, off-diagonals = 0
  cov_theta_hat_s <- diag(se_s^2, nrow = length(se_s), ncol = length(se_s))
  dimnames(cov_theta_hat_s) <- list(df$coeff_name, df$coeff_name)
  
  list(
    theta_hat_s     = theta_hat_s,
    cov_theta_hat_s = cov_theta_hat_s,
    cov_X_s         = NULL,   # not available for reduced models
    inds_theta_s    = 1
  )
}

# 
dat <- data.frame(study = c(1,2,1,2),
                  coeff_name = c("x1","x1", "x2", "x3"),
                  RR = c(1.0429,0.7212, -0.5332, 0.3115),
                  se = c(0.04819, 0.03625, 0.04844, 0.03604))

study_order <- unique(dat$study)
reduced_models_list <- lapply(study_order, function(s) {
  build_study(dat[dat$study == s, , drop = FALSE])
})
names(reduced_models_list) <- study_order
fit3 <- calc_linear_beta(reduced_models = reduced_models_list, ref_data = ref, cov_p_type = "reference")
round(fit2$linear_est$beta, 3)
