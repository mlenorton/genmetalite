# genmetalite 0.1.0

First packaged release, assembled from the `functions_for_simulations.R` and
`linear_estimate_centered.R` research scripts.

## New in this revision

* **The reference population can carry its sampling design.** `ref_data` now
  accepts a `survey.design` object from the \pkg{survey} package (an NHANES
  design built with `svydesign()`, say) or a `ref_design()` object built by
  hand, alongside the data frame and list-of-data-frames it always took. When a
  design is supplied, its weights are used for the covariance and its strata and
  first-stage clusters for the variance, so the estimator is design-aware end to
  end. `fit$linear_est$design_df` reports the design degrees of freedom.

  Only the covariates the studies actually use are kept from the reference: a
  survey design's `variables` slot routinely carries stratum, PSU, weight and
  identifier columns that would otherwise be treated as covariates.

  Finite population corrections are ignored (with-replacement variance, as NCHS
  recommends for NHANES) and calibration/post-stratification adjustments are not
  propagated, so a calibrated design gets a warning and a conservative variance.
  Replicate-weight designs are rejected with an explanation rather than silently
  mishandled.

* **Reweighting composes with the design.** Sampling weights become the *base
  weights* of the entropy-balancing problem, so each cohort gets one weight
  vector carrying both the sampling design and the covariate shift, with the
  population total preserved. `ipw_weights()` gained `base.weight` too, using it
  as regression weights in the membership model and as a multiplier on the odds
  weights.

  This combination — per-cohort weights plus a shared design — used to be
  rejected outright. It now works, because the per-cohort influence
  contributions are summed *before* the design is applied, which keeps the
  correlation induced by every cohort reweighting the same reference rows.
  Adding per-cohort variances instead would have understated it.

* **Entropy balancing has a native solver.** `entropy_balance()` solves the
  convex dual by Newton's method with a backtracking line search. It is now the
  default engine for `ebal_weights()`, so \pkg{ebal} is no longer needed;
  `engine = "ebal"` still calls `ebal::ebalance()`. Both solve the same strictly
  convex problem, so they agree to tolerance. The native solver handles base
  weights transparently, reports the largest standardised imbalance, and says so
  when a target lies outside the reference's convex hull instead of returning
  near-degenerate weights.

* **`as_reduced_model(vcov = )`.** Supply a coefficient covariance matrix, or a
  function returning one, for a study that is itself a complex sample. The
  default sandwich assumes independent rows; a clustered cohort needs something
  like `sandwich::vcovCL()`, and without it the pooled standard errors inherit
  the understatement. In simulation, clustered cohorts fitted with the plain
  sandwich dropped 95% interval coverage to about 0.89 even with the reference
  design handled correctly.

## Behaviour that changed

These are deliberate departures from the scripts, not refactors. Check them
against any saved simulation output before assuming results are comparable.

* **`linear_estimate()` is the centered implementation.** The estimator from
  `linear_estimate_centered.R` is now the single implementation, exported under
  both `linear_estimate()` and `linear_estimate_centered()`. The older
  uncentered `linear_estimate()`, which built the design map from the
  *uncentered* second-moment matrix `Phi` and returned two slightly different
  point estimates (`beta` and `my_beta`), is gone.

* **No mean-centering of inputs.** `calc_linear_beta()` no longer centers the
  numeric columns of the cohort and reference data. Centering was needed only by
  the `Phi` parametrization; the covariance parametrization is location
  invariant. The old centering was also applied to numeric columns but not to
  dummy columns derived from factors, which was inconsistent.

* **Everything runs in the slope-only basis.** Intercept columns are stripped
  from every covariance matrix. A constant column has zero variance, so keeping
  it made the matrices singular.

* **Factor expansion uses `model.matrix()` everywhere**, replacing
  `fastDummies::dummy_cols()` plus a regular-expression rename. The old rename
  (`_(\\d)$` -> `\\1`) only handled single-digit levels, so a factor with ten or
  more levels produced names that did not match the reference data.

* **The reference variance term is computed per cohort.** The Jacobian
  `d beta / d vec(Sigma)` is now formed separately for each cohort and the
  contributions summed, instead of using only `ref_pop_dummy[[1]]` and
  `ref_weights[[1]]`. When every cohort shares one reference and one weight
  vector the result is identical; when the reference is reweighted differently
  per cohort, the old code silently ignored all but the first.

* **`convert_reg_out_to_reduced_model_format()` keeps the intercept** by default,
  matching `run_regression()`. The two functions previously used opposite
  conventions. Pass `keep_intercept = FALSE` for the old behaviour. Both
  converters are also available under the shorter names `as_reduced_model()` and
  `as_reduced_model_genmetalite()`.

* **Empty factor levels are now reported, not silent.** A cohort that measured a
  categorical covariate but has no observations in one of its levels loses that
  dummy column, and the estimator then treats the level as a covariate the study
  never measured, filling its covariance blocks from the reference. That is a
  real assumption and it used to happen silently. `calc_linear_beta()` now warns
  and records what it found in the new `empty_levels` component of its return
  value; the check is also exported as `check_empty_levels()`.

* **Factor baselines are harmonised against the reference.** Cohort factors are
  releveled to the reference's level ordering before anything expands them. A
  cohort that happened not to observe the reference's first level previously got
  a different baseline from `model.matrix()`, so its `edu3` column meant
  "3 versus 2" while the reference's meant "3 versus 1" — and since the pipeline
  matches covariance blocks by name, the two were silently treated as the same
  covariate. If the reference baseline itself is unobserved in a cohort there is
  no valid way to continue, so that case is an error naming the fix.

* **`survey_reference_cov()` matches `survey_variance.R`.** Reconciled against
  the original script: the linearized contributions `q_j`, the PSU totals, the
  `n_h/(n_h-1)` factor, `df = #PSUs - #strata` and `n_eff` are all as written
  there, and the two agree to floating-point exactly on stratified, clustered
  and simple-random designs. Three additions: `lonely_psu = "adjust"` alongside
  `"remove"` and `"fail"`; validation of the weight vector; and the `p == 1` fix
  below.

* **`ipw_weights()` was rewritten.** The original could not run (see below); it
  now fits a logistic membership model and returns odds weights.

* **`compute_var(var_type = "best")` thresholds are arguments** (`min_obs`,
  `min_eff`) rather than hard-coded constants, and the heuristic is documented
  as experimental.

* **`calc_linear_beta()` no longer accepts `cov_cohort_X`.** The argument was
  never used. Passing it is now an error rather than silently ignored.

* **`cov_p_type` defaults to `"partitioned"`** and no longer accepts `NULL`.

## Bugs fixed

* `ipw_weights()` referenced an undefined object `X_test`, returned `eb.out$w`
  and `eb.out$converged` from a call that was commented out, and sourced
  `my_ipwpoint.R` by relative path. It could not run at all.

* `linear_estimate()` referenced `betaG` and `V_IF` outside the branch that
  defined them, so any call without individual-level reference data failed with
  "object not found" rather than returning the naive variance.

* `linear_estimate()` computed a constrained estimate `beta_hat_c` and then
  returned the *unconstrained* `betaG`, so `constr_vars` had no effect on the
  output. Constraints are now applied to both the estimate and its variance.

* `linear_estimate()` set `p <- ncol(ref_pop_dummy[[1]])`, overwriting the
  earlier `p <- length(cov_names)`. Any extra column in the reference data made
  the Jacobian the wrong shape.

* `linear_estimate()` special-cased a study with exactly one observed covariate
  by using the *whole* study covariance matrix instead of its 1x1 block.

* `make_sigma_X_ls()` called `cov.wt(ref_pop_dummy, ...)` on the whole list
  rather than `ref_pop_dummy[[k]]`, so every weighted covariance was computed
  from the wrong object.

* `make_sigma_X_ls()` tested `is.null(cov_cohort_X_s)` on a variable that was
  only assigned inside an `if (!is.null(reg_out))` block, erroring with "object
  not found" whenever `reg_out` was `NULL`.

* `make_sigma_X_ls()` reset `ref_weights` to a list of `NULL`s whenever
  `weighting` was `FALSE`, discarding weights the caller had supplied.

* `make_sigma_X_ls()` returned an undefined `ref_pop_dummy` when only
  `cov_ref_X` was supplied.

* `make_sigma_X_ls()` passed `weighting.max.iterations` to `ebal_weights()`,
  whose formal is `max.iterations`; the value fell into `...` and was ignored,
  so the iteration cap never took effect.

* `survey_reference_cov()` built its PSU totals with
  `t(vapply(upsu, ..., numeric(p)))`. When `p == 1` `vapply()` returns a plain
  vector and `t()` transposes it to `1 x n_psu`, so the stratum contribution
  failed to conform with the `p x p` accumulator and the function errored on any
  single-parameter model with more than one PSU per stratum.

* `run_regression()` returned `model_s$converged`, which is `NULL` for `lm`
  fits, so `converged` was unusable for linear models.

* `run_regression()` derived coefficient names by stripping the `X_s_dummy_mat`
  prefix. When a study had exactly one covariate, `lm()` names the coefficient
  just `X_s_dummy_mat` with no column name appended, so stripping produced an
  **empty** covariate name that propagated through the whole pipeline.

* `run_regression()` called `library(MASS)` inside the function for
  `dgp = "nb"`, attaching MASS to the global search path and masking
  `dplyr::select()` for the rest of the session. It now uses `MASS::glm.nb()`.

* `run_regression()` and `cov_factors()` called `stringr::str_replace()` without
  loading or importing stringr.

* `run_regression()` silently produced "object 'model_s' not found" for an
  unrecognised `dgp`; `dgp` and `var_type` are now validated with `match.arg()`.

* `cov_factors()` used `dplyr::relocate(any_of(ends_with("_bin")))`, which is
  not valid tidyselect (`any_of()` takes names, not a selection helper).

* `check_vars()` identified factor columns with `which(sapply(df, class) ==
  "factor")`. An ordered factor has class `c("ordered", "factor")`, so `sapply`
  returned a list and ordered factors were never checked.

* `ebal_weights()` errored with "subscript out of bounds" on its own default
  `weighting.vars = list()`, and rejected a character vector outright even
  though `make_sigma_X_ls()` passed one.

* `ebal_weights()` built its design with `model.matrix()` under the default
  `na.action`, which deletes incomplete rows. The resulting design matrix was
  shorter than the treatment vector, which is the cause of the "why are rows
  being dropped?!" note in the original. Missing values are now an error.

* `compute_var(var_type = "best")` read the outcome as `reg_out$data$y`, which
  is not the response for a `glm` fitted without a `data` argument; left
  `low_obs_names` undefined for non-`lm`/`glm` objects; and could index the
  sandwich matrix with names absent from `vcov()`.

* `calc_linear_beta()` added an `Intercept` column to the reference data only
  inside one branch of a name-comparison `if`, so the downstream
  `Xref[, cov_names]` lookup failed when the names happened to match.

* `calc_linear_beta()` passed `constr_vars` positionally into a call whose other
  arguments were named, relying on R's leftover-formal matching.

* Numerous functions returned `list(list = list(...))` — a stray `list =`
  argument to `return()`.

## Known limitations

* The study-block variance term is derived for least-squares residuals and is
  skipped, with a warning, for non-linear links.
* Reference weights are treated as fixed; estimating them adds an
  influence-function term that is not included.
* `compute_var(var_type = "best")` is a heuristic, not an established estimator.
