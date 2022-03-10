#' Visualizing the posterior sample
#' Creates a set of plots for visually evaluating convergence and mixing of the chains
#' from the MCMC sample of an object of class 'remiod'.
#'
#' @inheritParams JointAI::densplot
#' @inheritParams base::plot
#' @inheritDotParams graphics::matplot -x -y -type -xlab -ylab -pch -log -xlim -ylim
#' @param object an object inheriting from class 'remoid'
#' @param what select either trace or density plots from MCMC samples
#' @param mi.setting a list of arguments for extracting MI data set, which
#'                   will be used to update the one in `remoid` object.
#'                   Default is NULL, meaning no update to the mi.setting
#'                   in `remoid` object.
#' @param ... optional arguments.
#'
#' @import ggplot2
#' @importFrom JointAI traceplot densplot
#' @export
#'
mcmcplot <- function(object, what=c("trace","density"), subset = c(analysis_main=TRUE),
                     outcome = NULL, mi.setting=NULL, nrow=NULL, ncol=NULL,
                     use_ggplot=TRUE, mess=TRUE, warn=FALSE,...){
  if (!inherits(object, "remiod"))
    errormsg("Use only with objects of class remiod")

  old = object[['mi.setting']]
  if (!is.null(mi.setting)) miset = list_update(old, mi.setting)
  else miset = old

  ## MI setting args
  start = miset[['start']]
  end = miset[['end']]
  thin = miset[['thin']]
  seed = miset[['seed']]
  exclude_chains = miset[['exclude_chains']]

  what <- match.arg(what)
  plotfun <- switch(what,
                    "trace" = JointAI::traceplot,
                    "density" = JointAI::densplot)

  plotfun(object = object[['mc.mar']], subset=subset, outcome=outcome, use_ggplot=use_ggplot,
          exclude_chains=exclude_chains, start=start, end=end, thin=thin,
          nrow=nrow, ncol=ncol, mess=mess, warn=warn,...)

}
