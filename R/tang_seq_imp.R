#' Function to run MCMC iteration with the algorithm of Tang's sequential modeling
#'
#' @inheritParams JointAI::model_imp
#' @inheritParams commParams
#'
#' @import ordinal
#' @import progressr
#' @keywords internal
#' @noRd
tang_seq_imp = function(object, beta.init=NULL, ord_cov_dummy, seed = 1234, rinv=0.0001, scheme=2,
                        method="MAR", trtvar, n.chains=1,n.iter=10, burnin=1, thin=1,
                        progress.bar=FALSE, mess=TRUE,...){
  Mlist = get_Mlist(object)
  infolist = object$info_list
  refs = Mlist$refs
  coefs = object$coef_list

  nparm = 0
  for (k in 1:length(infolist)) nparm = nparm + infolist[[k]]$ncat - 1 + nrow(coefs[[k]])

  data = Mlist[["data"]]
  ntot = Mlist$N

  if (!("pattern" %in% colnames(data))) {
    data$pattern = get_pattern(Mlist = Mlist)
  }

  if (length(grep("Int", colnames(data),ignore.case = T))==0) {
    data = data.frame(`Intercept`=1, data)
  }
  colnames(data)[1] = "(Intecept)"

  dy = subset(data, select= names(Mlist$models))
  ncat = infolist[[1]]$ncat

  minpt = min(data$pattern)
  maxpt = max(data$pattern)
  nposttrt = numeric()

  if (minpt==0) {for (h in minpt:maxpt) nposttrt[h+1]= sum(data$pattern >= h)}
  else {for (h in minpt:maxpt) nposttrt[h]= sum(data$pattern >= h)}

  ## tag missing values
  misstag = miss_tag(Mlist = Mlist)

  obs = misstag$obs
  tmiss = misstag$totm
  indexall = misstag[, 3:4];
  indexmiss = misstag[, grep("tag", colnames(misstag))]

  ### intermittent missing data;
  itmisstag = subset(misstag, firstm < pattern)
  obs.itm = itmisstag$obs
  indexall.itm = itmisstag[, 3:4];
  indexmiss.itm = itmisstag[, grep("tag", colnames(misstag))]

  resplist = lapply(1:nrow(itmisstag), function(i) {
    respi = data[itmisstag$obs[i],,drop=FALSE]

    mitemn = which(is.na(respi))
    mitem = names(respi)[mitemn]
    out = lapply(mitem, function(x) {if (infolist[[x]]$ncat>2) 1:(infolist[[x]]$ncat)
      else unlist(unique(lapply(datat[,mitemn],levels)))})
    respiall = expand.grid(out)
    colnames(respiall) = mitem

    resps = rep.data.frame(respi, nrow(respiall))
    for (mit in mitem) {
      resps[, mit] = respiall[,mit]
      if (infolist[[mit]][['modeltype']]=="clm")
        resps[, mit] = factor(resps[, mit], ordered = attr(refs[[mit]], 'ordered'))
      else if (infolist[[mit]][['family']]=="binomial") resps[, mit] = factor(resps[, mit])
    }
    resps
  })

  ### get initial beta
  if (is.null(beta.init)) beta = beta_ini(object=object, n.chains=n.chains, seed=seed)[[1]]
  else beta = beta.init

  if (progress.bar) p <- progressr::progressor(steps = n.iter)
  ts = 0 # for sequencing imputed data set

  n.iter = n.iter + burnin
  for (t in 1:n.iter){
    set.seed(seed+t)

    datat = data

    ### update intermittent missing data;
    for (i in 1:nrow(itmisstag)){
      m_ini=indexall.itm[i,1]
      m_end=indexall.itm[i,2]

      temp_prob= probm(object=object, beta=beta, resp=resplist[[i]], ord_cov_dummy = ord_cov_dummy,
                       ini=m_ini, end=m_end)

      hh = which(rmultinom(1,size=1, prob=temp_prob/sum(temp_prob))[,1]==1)
      datat[obs.itm[i],] = resplist[[i]][hh,]
    }

    beta.upd = updatebeta(object, beta=beta, datat=datat, ord_cov_dummy = ord_cov_dummy,
                          rinv=rinv, scheme=scheme)

    beta = beta.upd[['beta']]

    # if (method=="CR"){
    #   if (datat$pattern < maxpt & datat[,trtvar]==1) datat[,trtvar] = 0
    # }

    Xmat = dsmat(Mlist=Mlist, ord_cov_dummy=ord_cov_dummy, data=datat)
    yv = rev(names(infolist))
    pat = minpt:maxpt

    ### imputation of post-dropout data under MAR;
    for(h in 1:length(nposttrt)) {
      ntem = ntot - nposttrt[h]
      s = pat[h]
      if (ntem >0 ){
        if (grepl("Intercept",coefs[[yv[s]]]$varname[1])) selvar = coefs[[yv[s]]]$varname
        else selvar = c("(Intercept)",coefs[[yv[s]]]$varname)
        xtemp = subset(Xmat[datat$pattern < s,], select=selvar)
        ncat = infolist[[yv[s]]]$ncat

        betai = unlist(beta[['b']][[yv[s]]])
        cuti = unlist(beta[['cexps']][[yv[s]]])

        probb= probfut(betax=betai, betacut=cuti, ncat=ncat, xtemp=xtemp, ntem=ntem)

        pb = probb
        if (ncat>2){
          for (kk in 2:(ncat-1)) pb[,kk] =  probb[,kk]-probb[,kk-1]
          pb = cbind(pb, matrix(1 - apply(pb,1,sum),ncol=1))
        } #else { pb = cbind(pb, 1-pb)}

        if (ncat>2) status= apply(pb,1, function(pbrow) which(rmultinom(1,1, pbrow)==1))
        else status = unlist(lapply(pb, function(x) rbinom(1,1,x) ))
        datat[datat$pattern < s, yv[s]]=status;
        #print(table(status))
        Xmat = dsmat(Mlist=Mlist, ord_cov_dummy=ord_cov_dummy, data=datat)
      }
    }

    ## Post-process of sampled parameters, to make it following the
    ## naming/parameterization rules of JAGS
    bvec = naming_param(beta = beta.upd[['beta']], coefs=coefs)
    if (t==1) bmats = bvec
    else bmats = rbind(bmats, bvec)

    if (t>burnin){
      tt = t - burnin

      if (tt %% thin == 0){
        ts = ts + 1
        dtimp = cbind(Imputation_ = ts, datat)
        if (ts==1) dtimps = dtimp
        else dtimps = rbind(dtimps, dtimp)
      }
    }
    rm(Xmat, datat)
    if (progress.bar) p(message = sprintf("MCMC iteration: %d", t))
  }

  return(list(beta = bmats, data.imp = dtimps))
}
