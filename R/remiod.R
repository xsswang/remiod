#' Reference-Based Multiple Imputation for Ordinal/Binary Response
#'
#' @inheritParams JointAI::model_imp
#' @inheritParams commParams
#'
#' @return A list includes (1) Information from JAGS modeling and MCMC samples
#'         and (2) A \code{data.frame} in which the original data (if
#'         \code{include = TRUE}) and the imputed datasets are stacked onto
#'         each other.\cr
#'         The variable \code{Imputation_} indexes the imputation, while
#'         \code{.rownr} links the rows to the rows of the original data.
#'         In cross-sectional datasets the
#'         variable \code{.id} is added as subject identifier.
#'
#' @rawNamespace import(coda, except=c(densplot,traceplot))
#' @import JointAI
#' @import utils
#' @import stats
#' @importFrom rjags coda.samples jags.model
#' @importFrom foreach foreach %dopar%
#' @import mathjaxr
#'
#' @examples
#'
#' \donttest{
#' data(schizow)
#'
#' test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
#'               trtvar = 'tx', algorithm = 'jags', method="MAR",
#'               ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
#'               n.iter = 10, thin = 2, warn = FALSE, seed = 1234)
#' }
#'
#' @export
remiod = function(formula, data, trtvar, refcats = NULL, family=NULL, method="MAR",
                  delta=0, algorithm=c("tang_seq","jags"), rinv=0.0001, scheme=2,
                  model_order=NULL, models=NULL, ord_cov_dummy = TRUE,
                  n.chains = 2, n.adapt = 100, n.iter = 1000, thin = 2,
                  start=NULL, end=NULL, seed=1234, exclude_chains=NULL,
                  subset=NULL, include=FALSE, mess=TRUE, warn=FALSE,
                  progress.bar=TRUE, ... ){

  if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
  #method <- match.arg(method)

  if(!missing(algorithm) & length(algorithm)>1) stop("Only one 'algorithm' allowed.")
  algorithm <- match.arg(algorithm)

  ## re-coding reference level of treatment variable to be 0
  if (!is.null(refcats[[trtvar]])) data[,trtvar] = ifelse(data[,trtvar]==refcats[[trtvar]], 0, 1)

  if (missing(formula)) errormsg("No model formula specified.")
  allvar = all.vars(formula)
  y = allvar[attr(terms(formula), "response")]
  # if (!is.ordered(data[,y]) | length(unique(data[,y])) != 2)
  #   errormsg("Response variable must be an binary or ordered variable.")

  if (algorithm == "jags"){
    ## control arm
    raw = subset(data, select=allvar)
    raw$id = 1:nrow(data)
    raw$rm = apply(raw, 1, function(x) sum(is.na(x)))
    rawm = subset(raw, rm>0 & raw[,trtvar]==1)

    varm = names(which(apply(raw,2,function(x) sum(is.na(x)))>0))

    if (!is.null(family)){
      ## list all linear predictors need to be monitored
      for (k in 1:length(varm)){
        mu_monit = paste0("mu_",varm[k],"[",rawm$id,"]")
        if (k==1) mu_monits = mu_monit
        else mu_monits = c(mu_monits,mu_monit)
      }

      mcsamp = glm_imp_custom(formula=formula, data=data, family=family, refcats=refcats,
                              model_order=model_order, models = models, ord_cov_dummy = ord_cov_dummy,
                              warn = warn, mess = mess, seed = seed, trtvar = trtvar,
                              n.chains=n.chains, n.adapt=n.adapt, n.iter=n.iter, thin=thin,
                              monitor_params = list(imps = TRUE, other_models = TRUE, other=c(mu_monits)),...)
    } else {
      ## list all linear predictors need to be monitored
      for (k in 1:length(varm)){
        eta_monit = paste0("eta_",varm[k],"[",rawm$id,"]")
        if (k==1) eta_monits = eta_monit
        else eta_monits = c(eta_monits,eta_monit)
      }

      if (is.null(models)){
        mcsamp = clm_imp_custom(formula = formula, data= data, model_order = model_order,
                                trtvar = trtvar, refcats=refcats, n.chains = n.chains,
                                n.adapt = n.adapt, n.iter = n.iter, thin = thin, warn = warn,
                                mess = mess, seed = seed, ord_cov_dummy = ord_cov_dummy,
                                monitor_params = list(imps = TRUE, other_models = TRUE, other=c(eta_monits)),...)
      } else {
        if (length(models)==1) {
          models = unlist(lapply(varm, function(x) models))
          names(models) = varm
        }

        if (models[1]=="opm"){
          ## list all cutoffs need to be monitored
          for (k in 1:length(varm)){
            c_monit = paste0("c_",varm[k])
            if (k==1) c_monits = c_monit
            else c_monits = c(c_monits,c_monit)
          }

          mcsamp = opm_imp_custom(formula = formula, data= data, model_order = model_order,
                                  trtvar = trtvar, refcats=refcats, n.chains = n.chains, models=models,
                                  n.adapt = n.adapt, n.iter = n.iter, thin = thin, warn = warn,
                                  mess = mess, seed = seed, ord_cov_dummy = ord_cov_dummy,
                                  monitor_params = list(imps = TRUE, other_models = TRUE,
                                                        other=c(eta_monits,c_monits)),...)
        }
       }
      }

    #if (n.iter < nimp) errormsg("Number of iterations is too small for requested number of imputations")
    if (n.iter==0) dimp = NULL
    else dimp = get_MI_RB(object=mcsamp, treatment=trtvar, method=method, delta = delta,
                          include = include, start=start, end=end, thin=thin,
                          seed = seed)

    mi.setting=list(trtvar=trtvar, ord_cov_dummy=ord_cov_dummy, method=method,
                    algorithm=algorithm, include=include, exclude_chains=exclude_chains,
                    thin=thin, start=start, end=end, seed = seed)

    mcsamp$comp_info[[3]] = packageVersion("remiod")
    names(mcsamp$comp_info)[[3]] = "remiod_version"
    out = structure(list(mc.mar=mcsamp, mi.setting=mi.setting, mi.data=dimp), class="remiod")
  } else {
    if (mess)
      msg("Prepare information for MCMC run with Tang's Sequential Modeling...")

    tsamp = remiod(formula = formula, data = data, algorithm="jags",  n.iter = 0,
                    family=family, refcats=refcats, model_order=model_order, models = models,
                    trtvar=trtvar, ord_cov_dummy=ord_cov_dummy, method=method, mess=FALSE)

    if (n.iter == 0L) out = tsamp
    #else if (n.iter < nimp) errormsg("Number of iterations is too small for requested number of imputations")
    else {
      init_seed <- unlist(lapply(get_rng(seed = seed, n_chains = n.chains), function(x) x[[2]]))

      object = tsamp$mc.mar
      inits = beta_ini(object=object, n.chains=n.chains, seed=seed)

      future_info <- get_future_info()

      run_tang <- ifelse(future_info$parallel, run_parallel_tang, run_sequential_tang)

      t0 <- Sys.time()
      tang_res <- run_tang(object = object, ord_cov_dummy = ord_cov_dummy, rinv = rinv,
                           method = method, trtvar = trtvar, scheme = scheme,
                           burnin = n.adapt, n.iter = n.iter, n.chains = n.chains,
                           inits = inits, thin = thin, seed = init_seed, mess = mess,
                           n_workers = future_info$workers, ...)

      tsamp$mc.mar$MCMC <- tang_res$mcmc

      # dtim = tang_res$dtimp
      # for (k in 2:n.chains){
      #   m1 = max(dtim[[k-1]]$Imputation_)
      #   dtim[[k]]$Imputation_ = dtim[[k]]$Imputation_ + m1
      # }

      tsamp$mi.data = tang_res$dtimp #data.table::rbindlist(dtim)
      t1 <- Sys.time()

      tsamp$mc.mar$comp_info[[1]] = t0
      tsamp$mc.mar$comp_info[[2]] = t1 - t0
      tsamp$mc.mar$comp_info[[4]] = future_info$call
      tsamp$mi.setting$algorithm = "tang_seq"
      }
    out = tsamp
    }
  return(out)
}
