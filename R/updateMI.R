#' Apply a MI method following initial run of \code{remoid} function
#'
#' @inheritParams commParams
#' @param object object inheriting from class \code{remoid}.
#' @param mi.setting a list of arguments for extracting MI data set, which
#'                   will be used to update the one in `remoid` object.
#'                   Default is NULL, meaning no update to the mi.setting
#'                   in `remoid` object.
#' @param mess logical; should messages be given? Default is TRUE.
#' @importFrom future multisession sequential plan
#' @import doFuture
#' @importFrom data.table rbindlist
#' @importFrom foreach foreach %dopar%
#'
#' @keywords internal
updateMI <- function(object, method=c("MAR","J2R","CR","delta"), delta=0,
                     mi.setting=NULL, mess=TRUE){
  if (!inherits(object, "remiod"))
    errormsg("Use only with 'remiod' objects.")

  if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
  method <- match.arg(method)

  old = object[["mi.setting"]]
  if (!is.null(mi.setting)) miset = list_update(old, mi.setting)
  else miset = old

  ## MI setting args
  trtvar = miset[['trtvar']]
  include = miset[['include']]
  start = miset[['start']]
  end = miset[['end']]
  seed = miset[['seed']]

  ord_cov_dummy = miset[['ord_cov_dummy']]
  # minspace = miset[['minspace']]

  algorithm = miset[['algorithm']]
  thin = miset[['thin']]
  exclude_chains = miset[['exclude_chains']]

  ## save MAR sampling
  mcsamp = object[["mc.mar"]]
  dtimp = object[["mi.data"]]

  if (method != "delta"){
    t0 <- Sys.time()
    if (algorithm == "tang_seq"){
      dimp = tang_MI_RB(object=mcsamp, dtimp=dtimp, treatment=trtvar, method=method,
                        ord_cov_dummy=ord_cov_dummy,
                        exclude_chains=exclude_chains, include=include)
    } else{
      dimp = get_MI_RB(object=mcsamp, delta=delta, treatment=trtvar, method=method,
                       exclude_chains=exclude_chains, thin=thin, include=include,
                       start=start, end=end, seed = seed)
    }
    t1 <- Sys.time()
  } else {
    future_info <- get_future_info()

    run_delta <- ifelse(future_info$parallel, delta_parallel, delta_seq)

    t0 <- Sys.time()
    dimp <- run_delta(object=mcsamp, dtimp=dtimp, algorithm=algorithm, delta=delta,
                      treatment=trtvar, method=method, n_workers = future_info$workers,
                      exclude_chains=exclude_chains, thin=thin, include=include,
                      start=start, end=end, seed = seed)
    t1 <- Sys.time()
  }

  if (mess)
    msg( paste0("Time to compute ", method," adjustment is %s min"), eval(round((t1 - t0)/60,2)) )

  out = structure(list(mc.mar=mcsamp, mi.setting=miset, mi.data=dimp), class="remiod")
  return(out)
}
