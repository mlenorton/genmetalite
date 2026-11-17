#' genmetalite: Generalized Meta-Analysis of Reduced Regression Models
#'
#' @description
#' `genmetalite` pools regression coefficients across studies that each fit a
#' *reduced* model containing only some of the covariates in a target *full*
#' model. Each study's reduced-model coefficient vector is mapped onto the
#' full-model parameter using covariate second moments, and the mapped estimates
#' are combined by generalized least squares.
#'
#' @section Workflow:
#' The usual entry point is [calc_linear_beta()], which accepts either
#' individual-level cohort data or pre-fitted reduced models:
#'
#' \enumerate{
#'   \item Reduced models are fitted with [run_regression()] (or supplied by the
#'     user through [as_reduced_model()]). Each reduced model contributes a
#'     coefficient vector, a robust variance matrix, and the empirical covariance
#'     of its own covariates.
#'   \item [make_sigma_X_ls()] assembles, for each study, the covariance matrix
#'     of the *full* covariate set. Blocks involving covariates the study did not
#'     measure must come from the reference population; the block among observed
#'     covariates can come either from the study itself
#'     (`cov_p_type = "partitioned"`) or from the reference
#'     (`cov_p_type = "reference"`).
#'   \item [linear_estimate()] maps and pools, and computes an
#'     influence-function variance.
#' }
#'
#' [IVW_estimate()] and [EB_estimate()] give an inverse-variance weighted
#' comparator and an empirical-Bayes style blend of the two.
#'
#' @section Reference reweighting:
#' When the reference population does not resemble a study, it can be reweighted
#' with [ebal_weights()] (entropy balancing) or [ipw_weights()] (inverse
#' probability of "membership" weighting) before its covariance is computed.
#'
#' @section Naming conventions:
#' Categorical covariates are expanded with [stats::model.matrix()] throughout,
#' so a factor `edu` with levels `1, 2, 3` yields columns `edu2` and `edu3`.
#' Coefficient names, covariance dimnames and reference-data column names must
#' agree; the package matches by name rather than position everywhere, and
#' errors rather than silently misaligning.
#'
#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom stats as.formula binomial coef complete.cases cov cov.wt
#'   glm lm model.frame model.matrix model.response na.omit nobs pnorm poisson
#'   predict qnorm quantile sd setNames vcov
#' @importFrom utils combn
## usethis namespace: end
NULL
