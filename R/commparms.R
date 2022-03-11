#' Common Parameters used by functions of \code{remiod}
#'
#' @param object object inheriting from class 'remoid'
#' @param trtvar the name of treatment variable. When necessary, its reference category,
#'               i.e. control arm, can be set in \code{refcats} argument.
#' @param method a method for obtaining multiple-imputed dataset. Options include
#'                \code{MAR}, \code{J2R}, \code{CR}, and \code{delta} adjustment.
#'                Default is MAR.
#' @param delta specific value used for Delta adjustment, applicable only
#'              for method="delta".
#' @param algorithm either algorithm \code{tang_seq} proposed by Tang (2018) or
#'                  \code{jags} the original method inherited in JAGS (Plummer 2003).
#' @param model_order optional. manually specify an order for imputation models.
#' @param ord_cov_dummy optional. specify whether ordinal variables should be treated as
#'                      categorical variables or continuous variables when they are
#'                      included as covariates in the sequential imputation models.
#'                      Default is `TRUE`, dummy variables will be created accordingly.
#' @param rinv a small number used to adjusting Fish information matrix
#' @param scheme scheme of distribution used for proposing coefficients of imputation models.
#'               scheme=1: beta ~ N( mean + inv(I)*score, inv(I));
#'               scheme=2: beta ~ N( mean , inv(I)).
#' @param subset subset of parameters/variables/nodes (columns in the MCMC
#'               sample). Follows the same principle as the argument
#'               \code{monitor_params} and \code{selected_parms}.
#' @param exclude_chains optional vector of the index numbers of chains that
#'                       should be excluded
#' @param start the first iteration of interest
#'              (see \code{\link[coda]{window.mcmc}})
#' @param end the last iteration of interest
#'            (see \code{\link[coda]{window.mcmc}})
#' @param progress.bar character string specifying the type of
#'                 progress bar. Possible values are "text" (default), "gui",
#'                 and "none" (see \code{\link[rjags]{update}}). Note: when
#'                 sampling is performed in parallel it is not possible to
#'                 display a progress bar.
#' @param seed optional; seed value (for reproducibility)
#' @param include logical, if TRUE, raw data will be included in imputed data sets
#'                with imputation ID = 0.
#' @param mi.setting a list of arguments related to multiple imputation, including
#'                   trtvar, algorithm, method, include, exclude_chains, thin, start,
#'                   end, and seed.
#' @name commParams
#' @keywords internal
NULL

