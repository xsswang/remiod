#' Extract a specified number of multiple imputed datasets
#'
#' @inheritParams commParams
#' @param M number of imputed datasets
#' @param minspace minimum number of iterations between iterations to be chosen
#'                 as imputed values (to prevent strong correlation between
#'                 imputed datasets in the case of high autocorrelation of the
#'                 MCMC chains).
#' @param mess logical; should messages be given? Default is TRUE.
#'
#' @return A \code{data.frame} in which the imputed datasets are stacked onto
#'         each other.
#'
#' @examples
#' \donttest{
#' # data(schizow)
#'
#' test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
#'               trtvar = 'tx', algorithm = 'jags', method="MAR",
#'               ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
#'               n.iter = 100, thin = 2, warn = FALSE, seed = 1234)
#'
#' extdt = extract_MIdata(object=test, method="J2R",mi.setting=NULL, M=10, minspace=2)
#'
#' }
#' @export
extract_MIdata <- function(object, method=c("MAR","J2R","CR","delta"), delta=0,
                           mi.setting=NULL, M=10, minspace=2, mess=FALSE){

  if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
  method <- match.arg(method)

  old = object[["mi.setting"]]
  if (!is.null(mi.setting)) miset = list_update(old, mi.setting)
  else miset = old

  # apply specified adjustment onto MAR samples
  upmidt = updateMI(object=object, method=method, delta=delta, mi.setting=miset, mess=mess)

  # randomly draw which iterations should be used as imputation
  Nimd = max(upmidt$mi.data$Imputation_)
  if (Nimd / minspace < M)
    errormsg("The total number of iterations (%s) is too small to select %s
             iterations with spacing of >= %s.", Nimd, M, minspace)

  seed = miset[['seed']]
  if (!is.null(seed)) set.seed(seed)

  cand_iters <- seq(from = sample.int(minspace, size = 1), to = Nimd, by = minspace)
  imp_iters <- sort(sample(cand_iters, size = M))

  # reduce to the relevant numbers
  midt <- subset(upmidt$mi.data, Imputation_ %in%  imp_iters)
  # re-id imputation
  mida <- merge(data.frame(Imputation_=imp_iters, Imp_=1:length(imp_iters)), midt, by="Imputation_")
  mida$Imputation_ = NULL
  return(mida)
}
