

probfut = function(betax,betacut, ncat, xtemp,ntem){
  if (is.matrix(betax)) betx = betax
  else betx = matrix(unlist(betax),ncol=1)

  xbeta = as.matrix(xtemp) %*% betx

  if (ncat>2) {
    xbet = matrix(rep(xbeta,ncat-1),ncol=ncat-1,byrow = FALSE)
    xb = xbet + matrix(rep(c(0,unlist(betacut)),ntem), nrow=ntem, ncol=ncat-1,byrow = T)
    } else xb = xbeta
  probtemp =1/(1+exp(-xb));
  return(probtemp);
}

probfut_mnar = function(betax,betacut, ncat, xtemp,ntem, adj=0, deltafir=1){
  if (is.matrix(betax)) betx = matrix(betax,ncol=1)
  else betx = matrix(unlist(betax),ncol=1)

  xbeta = as.matrix(xtemp) %*% betx
  if (deltafir==1) xbeta = xbeta - adj

  if (ncat>2) {
    xbet = matrix(rep(xbeta,ncat-1),ncol=ncat-1,byrow = FALSE)
    xb = xbet + matrix(rep(c(0,unlist(betacut)),ntem), nrow=ntem, ncol=ncat-1,byrow = T)
  }
  probtemp =1/(1+exp(-xb));
  return(probtemp);
}


probi = function(object, beta, resp, ini, end, ord_cov_dummy){
  Mlist = get_Mlist(object)
  coefs = object$coef_list
  infolist = object$info_list
  Xmat = dsmat(Mlist=Mlist, ord_cov_dummy=ord_cov_dummy, data=resp)
  yv = rev(names(coefs))

  probtemp=1;
  for (k in ini : end){
    ncat = infolist[[yv[k]]]$ncat

    if (ncat>2) {
      sta = as.numeric(unlist(resp[yv[k]]))
      selvar = c("(Intercept)", coefs[[yv[k]]]$varname)
    } else {
      sta = unlist(resp[yv[k]])
      selvar = coefs[[yv[k]]]$varname
      }
    xtemp = subset(Xmat, select=selvar)

    betak = beta[['b']][[yv[k]]]
    cuts = beta[['cexps']][[yv[k]]]
    xbeta = sum(betak * xtemp)

    if (ncat>2) {
      if (sta==1) {pk = 1- 1/(1+exp(xbeta));}
      if (sta==2) {pk = 1/(1+exp(xbeta)) - 1/(1+exp(cuts[1]+ xbeta))}
      if (sta>2 & sta<ncat) {pk = 1/(1+exp(cuts[sta-2]+ xbeta)) - 1/(1+exp(cuts[sta-1]+ xbeta))}
      if (sta==ncat) {pk =1/(1+exp(cuts[sta-2]+ xbeta));}
    } else if (ncat==2) {
      #xbeta = matrix(betak,nrow=1) %*% matrix(xtemp,ncol=1)
      if (sta ==1) pk = 1/(1+exp(-xbeta))
      else pk =  1/(1+exp(xbeta))
      }
    probtemp= probtemp*pk;
  }
  if (is.list(probtemp)) probtemp = unlist(probtemp)
  return(probtemp);
}

mnorm = function(B, A, seed){
  Lx = chol(A);
  df=nrow(A);
  b = backsolve(t(U), B)  ## trisolv(2,Lx, B)
  mtemp = b + rnorm(df)
  varnorm = t(backsolve(t(U), mtemp))  # (trisolv(1,Lx,mtemp));
  return (varnorm);
}


#' calculate likelihood of each observation
#' @param prob the vector to be expanded
#' @noRd
likprob = function(prob){
  pcum = cumsum(prob)

  likj = numeric()
  for (j in 1:(length(ns)-1))
    likj[j] = prob[j]*log(pcum[j]/pcum[j+1]) + prob[j+1]*log(prob[j+1]/pcum[j+1])

  lik = sum(likj)
  return(lik)
}


#' tag the position of missing values
#' @param Mlist setting list of imputation models
#' @noRd
miss_tag = function(Mlist){

  data = Mlist$data
  if (!("pattern" %in% colnames(data))) {
    data$pattern = get_pattern(Mlist = Mlist)
  }

  yv = names(Mlist$models)
  nvisit = length(yv)
  dy = subset(data, select= rev(yv))
  totm = apply(dy,1,function(x) sum(is.na(x)))

  for (i in 1:nrow(dy)){
    locmis = which(is.na(dy[i,]))
    if (length(locmis)>0) locm = c(locmis, rep(NA,nvisit - length(locmis)))
    else locm = rep(NA, nvisit)
    if (i==1) locs = locm
    else locs= rbind(locs,locm)
  }

  colnames(locs) = paste0("tag",1:length(Mlist$models))
  out = data.frame(obs=1:nrow(dy), totm=totm, firstm=locs[,1], pattern=data$pattern, locs)
  out = subset(out, totm>0)
  return(out)
}


#' get initial values of beta and cutoffs
#' @noRd
beta_ini = function(object, n.chains=2, seed=123){
  init_seed <- get_rng(seed = seed, n_chains = n.chains)
  coefs = object[["coef_list"]]
  models = object[["models"]]

  binits = lapply(seq(n.chains), function(k){
            set.seed(init_seed[[k]][[2]][[1]])
            lapply(names(coefs), function(x) { if (models[[x]] == "clm") runif(nrow(coefs[[x]]) + 1)
                                        else  runif(nrow(coefs[[x]])) }) # +1 to add intercept
            })

  for (k in 1:n.chains){
    names(binits[[k]]) <- names(coefs)
  }

  ## add cutoffs for ordinal
  infolist = object$info_list
  cinits = lapply(seq(n.chains), function(k){
    set.seed(init_seed[[k]][[2]][[1]])
    lapply(infolist, function(x) { if (x$ncat>2) runif(x$ncat-2,-1,1)
                                   else NULL})
  })

  inits = list()
  for (i in 1:n.chains){
    cexps = lapply(cinits[[i]], function(x) { if (!is.null(x)) cumsum(exp(x))} )
    initi = list(b = binits[[i]], c = cinits[[i]], cexps = cexps)
    init = list(initi)
    inits = c(inits, init)
  }

  return(inits)
}

#' assign missing pattern
#' @noRd
get_pattern = function(Mlist){
  yv = names(Mlist$models)

  data = Mlist$data
  dy = subset(data, select = yv)
  pattern = rep(ncol(dy), nrow(dy))

  for (i in 1:nrow(dy)){
    for (h in 1:length(yv))
      if (sum(is.na(dy[i,1:h]))==h) pattern[i] = ncol(dy)-h
  }
  return(pattern)
}

###

#' create design matrix when \code{ord_cov_dummy = TRUE}
#' @noRd
dsmat = function(Mlist, ord_cov_dummy = TRUE, data=NULL){
  if (is.null(data)) data = Mlist$data
  fixed = Mlist$fixed
  auxvars = Mlist$auxvars
  refs = Mlist$refs
  terms_list = Mlist$terms_list
  models = Mlist$models
  if (!any(models == "clm")) ord_cov_dummy = TRUE

  # * design matrix -----------------------------
  if (ord_cov_dummy)
    Xmat <- model_matrix_combi(fmla = c(fixed, auxvars), data = data, refs = refs, terms_list =  terms_list)
  else {
    colnames(data)[grep("int",colnames(data),ignore.case=TRUE)] =  "(Intercept)"
    yv = rev(names(Mlist$models))
    covs = setdiff(all_vars(fixed), yv)
    Xmat <- subset(data,select = c("(Intercept)",covs, yv))
    Xmat[,yv] = prep_covoutcomes(Xmat[,yv])
  }
  return(Xmat)
}

#' naming monitored parameter
#' @noRd
naming_param = function(beta, coefs){
  yv = names(coefs)
  bvec = numeric()
  for (v in yv){
    if (grepl("Intercept",coefs[[v]]$varname[1])) names(beta[['b']][[v]]) <-  coefs[[v]]$coef
    else names(beta[['b']][[v]]) <- c(paste0('gamma_',v,'[1]'), coefs[[v]]$coef)

    if (is.null(beta[['c']][[v]])) bvec = c(bvec,beta[['b']][[v]])
    else {
      names(beta[['c']][[v]]) <- c(paste0('gamma_',v,'[',1+seq(length(beta[['c']][[v]])),']'))

      nc = length(beta[['cexps']][[v]])
      #beta[['cexps']][[v]] = beta[['cexps']][[v]] + unlist(rep(beta[['b']][[v]][1],nc))
      names(beta[['cexps']][[v]]) <- c(paste0('gamma_',v,'[',1+seq(length(beta[['c']][[v]])),']'))
      bvec = c(bvec,beta[['cexps']][[v]], beta[['b']][[v]])
    }
  }

  bm = matrix(bvec,nrow=1)
  colnames(bm) = names(bvec)
  return(bm)
}


beta_mat_to_list = function(beta, coefs, infolist){
  yv = rev(names(coefs))
  blist = list()
  clist = list()
  for (v in yv){
    nambeta <- c(paste0('gamma_',v,'[1]'), coefs[[v]]$coef)
    bvec = beta[,nambeta,drop=FALSE]
    bv = list(bvec)
    names(bv) = v
    blist = c(blist, bv)

    ncat = infolist[[v]]$ncat
    namcexps <- paste0('gamma_',v,'[',1+seq(ncat-2),']')
    cexpvec = beta[,namcexps,drop=FALSE]
    cexpv = list(cexpvec)
    names(cexpv) = v
    clist = c(clist, cexpv)
  }
  return(list(b=blist, cexps=clist))
}


#' @noRd
run_parallel_tang <- function(object, inits, seed = 1234, rinv=0.0001, scheme=2,
                              n.chains=2, n.iter=1000, burnin=100, thin=1,
                              ord_cov_dummy, n_workers, trtvar, method="MAR",
                              mess = TRUE, ...) {

  if (any(burnin > 0L, n.iter > 0L)) {

    if (mess)
      msg("Parallel sampling with %s workers started (%s).", eval(n_workers), Sys.time())

    res <- foreach::`%dopar%`(foreach::foreach(i = seq_along(inits)),
                              tang_seq_imp(object=object, beta.init=inits[[i]], seed=seed[i],
                                           rinv=rinv, scheme=scheme, ord_cov_dummy=ord_cov_dummy,
                                           method=method, trtvar=trtvar, n.iter=n.iter,
                                           burnin=burnin, thin=thin, mess=mess, ...)
    )

    mcmc_raw <- coda::as.mcmc.list(lapply(res, function(x) coda::mcmc(x$beta)))
    mcmcw <- window(mcmc_raw, start=burnin+1, thin=thin)

    dtimp = lapply(res, function(x) x$data.imp)
    # for (i in 1:length(dtimp)){
    #   ni = max(dtimp[[i]]$Imputation_)
    #   if (i>1) dtimp[[i]]$Imputation_ = dtimp[[i]]$Imputation_ + ni*(i-1)
    # }
    # dtim = data.table::rbindlist(dtimp)

    out = list(mcmc=mcmcw, dtimp = dtimp)
  } else out = list(mcmc=NULL, dtimp = NULL)
  return(out)
}

run_sequential_tang <- function(object, inits, seed = 1234, rinv=0.0001, scheme=2,
                                n.chains=2, n.iter=1000, burnin=100, thin=1, ord_cov_dummy,
                                trtvar, method="MAR", mess = TRUE, ...) {

  if (any(burnin > 0L, n.iter > 0L)) {

    if (mess)
      msg("Sequential sampling for %s chains started (%s).", eval(n.chains), Sys.time())

    for (i in 1:n.chains){
      mcs <-  tang_seq_imp(object=object, beta.init=inits[[i]], seed=seed[i],
                           rinv=rinv, scheme=scheme, ord_cov_dummy=ord_cov_dummy,
                           method=method, trtvar=trtvar, n.iter=n.iter,
                           burnin=burnin, thin=thin,mess=mess,...)

      dtim = list(mcs$data.imp)
      #ni = max(dtim$Imputation_)

      reslst= list(coda::mcmc(mcs$beta))
      if (i==1) {
        mcmc = reslst
        dtimp = dtim
      } else {
        mcmc = c(mcmc, reslst)
        #dtim$Imputation_ = dtim$Imputation_ + ni*(i-1)
        dtimp = c(dtimp, dtim)
        }
    }
    mcmc_raw <- coda::as.mcmc.list(mcmc)
    mcmcw <- window(mcmc_raw, start=burnin+1, thin=thin)

    out = list(mcmc=mcmcw, dtimp = dtimp)
  } else out = list(mcmc=NULL, dtimp = NULL)
  return(out)
}

