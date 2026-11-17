# genmetalite

<!-- badges: start -->
<!-- badges: end -->

**Generalized meta-analysis of reduced regression models.**

Several studies each fit a regression model, but none of them measured every
covariate you care about. `genmetalite` pools their coefficients into a single
estimate of the *full* model, using covariate second moments from an external
reference population to correct for what each study left out.

The naive alternative — inverse-variance weighted meta-analysis of the
coefficients as reported — is biased whenever an omitted covariate is correlated
with a retained one, because each study is then estimating a different marginal
effect. `genmetalite` maps each study's coefficients onto a common target before
pooling, and propagates the extra uncertainty that mapping introduces.

## Installation

``` r
# install.packages("remotes")
remotes::install_github("mlenorton/genmetalite")
```

## Quick start

``` r
library(genmetalite)

set.seed(2026)

# A population with correlated covariates.
make_pop <- function(n) {
  x1 <- rnorm(n)
  x2 <- 0.6 * x1 + rnorm(n, sd = 0.8)
  x3 <- 0.4 * x2 + rnorm(n)
  data.frame(x1 = x1, x2 = x2, x3 = x3)
}
beta <- c(x1 = 1, x2 = -0.5, x3 = 0.3)

# A large reference sample with all three covariates, no outcome needed.
ref <- make_pop(20000)

# Three studies, each missing a different covariate.
make_study <- function(n, drop) {
  d <- make_pop(n)
  d$y <- as.numeric(2 + as.matrix(d) %*% beta + rnorm(n))
  d[, setdiff(names(d), drop), drop = FALSE]
}
cohorts <- list(
  make_study(1500, "x3"),
  make_study(1500, "x2"),
  make_study(1500, "x1")
)

fit <- calc_linear_beta(cohort_data = cohorts, ref_data = ref)
summary(fit)
#>   term   estimate  std.error statistic  conf.low  conf.high  ivw.estimate
#> 1   x1  1.0099     0.0289     34.9       0.9532     1.0666    0.8304
#> 2   x2 -0.5397     0.0336    -16.1      -0.6054    -0.4739   -0.1908
#> 3   x3  0.3643     0.0239     15.2       0.3174     0.4112    0.2959
```

The `estimate` column recovers the truth; `ivw.estimate` — ordinary
inverse-variance weighting — does not.

## How it works

For study *s* with observed covariates *o* and unobserved covariates *m*, the
reduced-model coefficients relate to the full-model parameter through

<!-- A_s = [ I | B_s ],  B_s = Sigma_oo^{-1} Sigma_om -->
```
theta_s  ~=  A_s beta,       A_s = [ I | B_s ],   B_s = Sigma_oo^-1 Sigma_om
```

so the pooled estimate is the generalized least squares solution

```
beta_hat = M sum_s A_s' V_s^-1 theta_hat_s,    M = ( sum_s A_s' V_s^-1 A_s )^-1
```

The pipeline is three steps, each of which you can run separately:

| Step | Function | What it does |
|---|---|---|
| 1 | `run_regression()` | Fit each study's reduced model (or bring your own with `as_reduced_model()`) |
| 2 | `make_sigma_X_ls()` | Assemble each study's covariance over the full covariate set |
| 3 | `linear_estimate()` | Map, pool, and compute the influence-function variance |

`calc_linear_beta()` runs all three and also returns `IVW_estimate()` and
`EB_estimate()` for comparison.

### Where the covariance blocks come from

`cov_p_type` controls the block among the covariates a study *did* measure:

- **`"partitioned"`** (default) — that block is the study's own empirical
  covariance; only the cross-blocks with unmeasured covariates are borrowed from
  the reference. More faithful to each study, but it makes `A_s` depend on the
  same sample as `theta_hat_s`, which adds a variance term.
- **`"reference"`** — every block comes from the reference population. Simpler
  variance, but it assumes the reference covariance describes each study.

### The variance

`linear_estimate()` reports up to three terms:

1. **Naive (I)** — `M`, treating the mapping `A_s` as known.
2. **Reference (II)** — the cost of estimating `Sigma` in a finite reference
   sample, `O(1/n_ref)`. Requires individual-level `ref_data`; if you supply only
   `cov_ref_X`, this term is dropped and you get a warning.
3. **Study block (I-b)** — only under `"partitioned"`, `O(1/n_s)`. Requires
   individual-level `cohort_data`. Exactly zero for a study that measured
   everything, or when the coefficients of what it missed are zero.

In a 300-replicate simulation (three studies of n = 1500, reference n = 20000,
each study missing one of three correlated covariates), the estimator is
unbiased and the reported standard errors match the empirical spread, with 95%
interval coverage of 0.953, 0.950 and 0.960.

With a clustered, weighted reference (60 PSUs in 10 strata) the design-aware
variance restores most of what ignoring the design loses: coverage of 0.932,
0.932 and 0.928, against 0.920, 0.864 and 0.920 when only the weights are used.
The residual shortfall below 0.95 is the usual first-order linearization
behaviour at modest PSU counts.

## The reference population and how it was sampled

`ref_data` accepts four things:

| Input | Use when |
|---|---|
| a data frame | the reference is an independent sample |
| `ref_design(data, weights, strata, psu)` | you know the design but don't use the survey package |
| a `survey.design` object | you built the design with `survey::svydesign()` |
| a list of data frames | genuinely separate reference samples, one per study |

Passing a design makes the estimator design-aware end to end: `Sigma` becomes a
survey-weighted covariance, and the reference variance term uses the stratified
ultimate-cluster linearization instead of the independent-sampling formula.

``` r
library(survey)
des <- svydesign(id = ~SDMVPSU, strata = ~SDMVSTRA, weights = ~WTMEC2YR,
                 data = nhanes, nest = TRUE)

fit <- calc_linear_beta(cohort_data = cohorts, ref_data = des)
fit$linear_est$design_df   # PSUs minus strata; often small
```

Use a *t* reference distribution on `design_df` rather than a normal when it is
small, which is common for survey subgroups. Finite population corrections are
ignored (with-replacement variance, as NCHS recommends), and calibrated or
post-stratified designs get a warning and a conservative variance.

## Reweighting the reference population

If the reference does not resemble a study, reweight it:

``` r
fit <- calc_linear_beta(
  cohort_data      = cohorts,
  ref_data         = des,                     # survey design, or a data frame
  weighting.vars   = c("age", "sex", "bmi"),
  weighting.method = "ebal"                   # or "IPW"
)
```

- `ebal_weights()` — entropy balancing. Matches the requested moments exactly
  while keeping weight dispersion low. Solved natively by `entropy_balance()`;
  `engine = "ebal"` calls `ebal::ebalance()` instead.
- `ipw_weights()` — odds weights from a logistic membership model. No extra
  dependency, but check the effective sample size it reports.

**Reweighting composes with the design.** The sampling weights become the *base
weights* of the balancing problem, so each cohort gets one weight vector
carrying both the sampling design and the covariate shift, with the population
total preserved. The strata and clusters still apply. Internally, the per-cohort
influence contributions are summed *before* the design is applied, which keeps
the correlation induced by every cohort reweighting the same reference rows.

**Caveat:** weights are treated as fixed. Estimating them contributes a further
influence-function term that is not included, so standard errors under
reweighting are mildly anti-conservative. A joint bootstrap over study and
reference rows is the recommended check.

## When a *study* is a complex sample

`linear_estimate()` takes each study's coefficient covariance as given and pools
with its inverse, so that matrix has to reflect however the study was sampled.
The default sandwich assumes independent rows. If a cohort is clustered, supply
a cluster-robust matrix:

``` r
rm1 <- as_reduced_model(
  fit_lm,
  vcov = function(f) sandwich::vcovCL(f, cluster = mydata$psu)
)
```

In simulation, clustered cohorts fitted with the plain sandwich dropped 95%
interval coverage to about 0.89 even with the reference design handled
correctly.

## Conventions worth knowing

- **The outcome column must be named `y`** in every `cohort_data` data frame.
  Every other column is treated as a covariate.
- **Factors are expanded with `model.matrix()`**, so a factor `edu` with levels
  `1, 2, 3` becomes columns `edu2` and `edu3`. Everything downstream is matched
  by *name*, never by position, and mismatches raise errors rather than
  silently misaligning.
- **No centering is needed.** The estimator uses `Var(X)` and slope
  coefficients, both invariant to location.
- **Everything runs in the slope-only basis.** Intercept columns are stripped
  from every covariance matrix, since a constant column has zero variance and
  would make them singular.

## Getting help

``` r
vignette("genmetalite")
?calc_linear_beta
```

## License

MIT © Emily Norton
