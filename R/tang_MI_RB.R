#' Implement the algorithm proposed by Tang
#'
#' @inheritParams commParams
#'
#' @param dtimp imputed complete data sets from \code{remiod} function.
#' @param treatment treatment variable.
#' @keywords internal

tang_MI_RB = function(object, dtimp, treatment, method="MAR", delta=0, ord_cov_dummy=FALSE,
                      exclude_chains=NULL, include=FALSE){
  Mlist = get_Mlist(object)
  infolist = object$info_list
  coefs = object$coef_list

  data = Mlist[["data"]]
  ntot = Mlist$N

  if (!("pattern" %in% colnames(data))) {
    data$pattern = get_pattern(Mlist = Mlist)
  }

  if (length(grep("Int", colnames(data),ignore.case = T))==0) {
    data = data.frame(Intercept=1, data)
  }

  dy = subset(data, select= names(Mlist$models))

  minpt = min(data$pattern)
  maxpt = max(data$pattern)
  yv = rev(names(infolist))
  pat = minpt:maxpt

  loccpatactive = list()
  locpatactive = list()
  for (k in 1:length(yv)){
    loc = which(data$pattern < k & data[,trtvar]==1)
    if (length(loc) > 0) loclist = list(loc)
    else loclist = list(NULL)
    names(loclist) = yv[k]
    loccpatactive = c(loccpatactive, loclist)

    lock = which(data$pattern == k-1 & data[,trtvar]==1)
    if (length(lock) > 0) klist = list(lock)
    else klist = list(NULL)
    names(klist) = yv[k]
    locpatactive = c(locpatactive, klist)
  }

  if (method == "MAR") {
    chains = seq_along(dtimp)
    if (!is.null(exclude_chains)) {
      chains <- chains[-exclude_chains]
      dtims = dtimp[chains]
    } else dtims = dtimp
  }
  else {
    bmcmc = object$MCMC
    chains = seq_along(bmcmc)
    if (!is.null(exclude_chains)) {
      chains <- chains[-exclude_chains]
      bmcmc = bmcmc[chains]
    }

    bmcmc = lapply(bmcmc, function(x) as.matrix(x))
    dtims = list()

    for (k in 1:length(bmcmc)){
      bmc = bmcmc[[k]]

      for(j in 1:nrow(bmc)){
        beta = beta_mat_to_list(beta=bmc[j,,drop=FALSE], coefs=coefs, infolist=infolist)

        datat = subset(dtimp[[k]], Imputation_ == j)
        Xmat = dsmat(Mlist=Mlist, ord_cov_dummy=ord_cov_dummy, data=datat)

        ### imputation of post-dropout data under MNAR;
        for(h in 1:length(yv)) {
          yh = yv[h]
          if (!is.null(loccpatactive[[yh]]) ){
            ntem = length(loccpatactive[[yh]])
            xtemp = subset(Xmat[loccpatactive[[yh]],,drop=FALSE], select=c("(Intercept)",coefs[[yh]]$varname))
            ncat = infolist[[yh]]$ncat

            betai = unlist(beta[['b']][[yh]])
            cuti = unlist(beta[['cexps']][[yh]])
            tcoef = coefs[[yh]]$coef[coefs[[yh]]$varname==trtvar]

            if (method=="CR"){
              betai[1,tcoef] = 0
              probb= probfut_mnar(betax=betai, betacut=cuti, ncat=ncat, xtemp=xtemp, ntem=ntem,
                                  adj=0, deltafir=0)
            } else if (method=="delta") {

              betai[1,tcoef] = betai[1,tcoef] - delta
              probb= probfut_mnar(betax=betai, betacut=cuti, ncat=ncat, xtemp=xtemp, ntem=ntem,
                                  adj=0, deltafir=0)
            } else if (method=="deltafir") {
              dadj = rep(0, length(loccpatactive[[yh]]))
              dloc = which(locpatactive[[yh]] %in% loccpatactive[[yh]])
              dadj[dloc] = delta
              probb= probfut_mnar(betax=betai, betacut=cuti, ncat=ncat, xtemp=xtemp, ntem=ntem,
                                  adj=dadj, deltafir=1)
            }

            pb = probb
            for (kk in 2:(ncat-1)) pb[,kk] =  probb[,kk]-probb[,kk-1]
            pb = cbind(pb, matrix(1 - apply(pb,1,sum),ncol=1))

            status= apply(pb,1, function(pbrow) which(rmultinom(1,1, pbrow)==1))
            datat[loccpatactive[[yh]], yh]=status;
            Xmat = dsmat(Mlist=Mlist, ord_cov_dummy=ord_cov_dummy, data=datat)
          }
        }

        if (j==1) dtip = datat
        else dtip = rbind(dtip, datat)
      }
      dtim = list(dtip)
      dtims = c(dtims, dtim)
    }
  }

  ## renumbering ids of imputed data sets
  if (length(dtims)>1){
  for (h in 2:length(dtims)){
    nh = max(dtims[[h-1]]$Imputation_)
    dtims[[h]]$Imputation_ = dtims[[h]]$Imputation_ + nh
  }}

  if (include){
    data = cbind(Imputation_=0, data)
    dtims = c(list(data), dtims)
  }

  if (method %in% c('delta','deltafir')){
    for (i in 1:length(dtims)) {
      dtims[[i]] <- cbind(dtims[[i]], "delta" = delta)
    }
  }

  ## stack imputed data sets from all chains
  dtall = data.table::rbindlist(dtims)

  return(dtall)
}
