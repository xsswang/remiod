#' Implement controlled multiple imputation algorithms proposed by Tang
#'
#' Internal function, creates multiple imputed datasets based on assigned
#' imputation method with the algorithm of Tang's sequential modeling.
#'
#' @inheritParams commParams
#'
#' @return multiple imputed datasets stacked onto each other (i.e., long format;
#' optionally including the original incomplete data).\cr
#' The variable \code{Imputation_} indexes the imputation, while
#'         \code{.rownr} links the rows to the rows of the original data.
#'         In cross-sectional datasets the
#'         variable \code{.id} is added as subject identifier.
#'
#' @param dtimp imputed complete data sets from \code{remiod} function.
#' @param treatment treatment variable.

tang_MI_RB = function(object, dtimp, treatment, method="MAR", delta=0, ord_cov_dummy=FALSE,
                      exclude_chains=NULL, include=FALSE, thin=1){
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
  colnames(data)[1] = "(Intercept)"

  yv = rev(names(infolist))

  if (method == "J2R"){
    ### keep raw data
    rawdt = data
    ### J2R assumes that all treatment benefits are gone immediately after dropout
    ### change rows with missing outcome to all missing
    respmis = which(data$pattern < max(data$pattern) & data[,treatment]==1)
    data[respmis, yv] = NA
    data$pattern[respmis] = 0
  }

  dy = subset(data, select= names(Mlist$models))

  minpt = min(data$pattern)
  maxpt = max(data$pattern)
  pat = minpt:maxpt

  loccpatactive = list()
  locpatactive = list()
  for (k in 1:length(yv)){
    loc = which(data$pattern < k & data[,treatment]==1)
    if (length(loc) > 0) loclist = list(loc)
    else loclist = list(NULL)
    names(loclist) = yv[k]
    loccpatactive = c(loccpatactive, loclist)

    lock = which(data$pattern == k-1 & data[,treatment]==1)
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
    thin.old = attr(bmcmc[[1]], "mcpar")[3]
    if (thin > thin.old) thin.new = thin %/% thin.old
    else thin.new = thin.old

    chains = seq_along(bmcmc)
    if (!is.null(exclude_chains)) {
      chains <- chains[-exclude_chains]
      bmcmc = bmcmc[chains]
      dtimc = dtimp[chains]
    }

    if (thin > 1) {
      bmcmt = window(bmcmc, thin=thin)

      dtimp =  lapply(dtimc, function(dx) {
        dk = seq(1, max(dx$Imputation_), by= thin.new)
        dxs = subset(dx, Imputation_ %in% dk)
        dxs$Imputation_ = 1 + (dxs$Imputation_ %/% thin.new)
        dxs
      })
    }

    bmcmc = lapply(bmcmt, function(x) as.matrix(x))
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
            tcoef = coefs[[yh]]$coef[coefs[[yh]]$varname==treatment]

            if (method=="CR"){
              betai[1,tcoef] = 0
              probb= probfut_mnar(betax=betai, betacut=cuti, ncat=ncat, xtemp=xtemp, ntem=ntem,
                                  adj=0, deltafir=0)
            } else if (method=="J2R"){
              betai[1,tcoef] = 0
              probb= probfut_mnar(betax=betai, betacut=cuti, ncat=ncat, xtemp=xtemp, ntem=ntem,
                                  adj=0, deltafir=0)
            }
            else if (method=="delta") {
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
        if (method=="J2R"){
          rawdti = data.frame(Imputation_ = datat$Imputation_, rawdt)
          rawdti[is.na(rawdti)] = datat[is.na(rawdti)]
          datat = rawdti
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
