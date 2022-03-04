
get_likelihood = function(xtemp, ytemp, betai, cuti, ncat)
{
  ncov = ncol(xtemp);
  nitemp = nrow(xtemp);
  covdim = ncat-2+ncov;

  score=rep(0,covdim);
  Fish=matrix(0,covdim,covdim);
  loglik=0;

  xtemp = as.matrix(xtemp)
  if (ncat>2){
    cuti = unlist(cuti)
    cuts = cumsum(exp(cuti))
  }
  if (ncat==2 & min(ytemp)==0) ytemp = ytemp + 1

  bi = matrix(unlist(betai), ncol=1)
  xbeta = xtemp %*% bi

  picum=matrix(1,nitemp,ncat)
  picum[,1] = 1- 1/(1+exp(xbeta))

  if (ncat > 2){
    for( k in 2:(ncat-1) ) {
      for (j in 1:nitemp) {
        aeb = exp(cuts[k-1]+ xbeta[j,1])
        picum[j,k] = 1 - 1/(1+ aeb)
      }
    }}

  pitemp=picum
  for( k in 2:ncat ) {
    pitemp[,k]=picum[,k]-picum[,k-1];
  }

  partgamma = matrix(0, nitemp, covdim*ncat);
  xpand = as.matrix(xtemp)
  if ( ncat>2 ) xpand = cbind(matrix(0,nitemp,ncat-2), xpand)
  for (k in 1:(ncat-1)){
    if (ncat>2 & k>1) xpand[,k-1] <- unlist(rep(exp(cuti[k-1]),nitemp))
     #  xpand[,k-1] = j(nitemp,1,exp(betaitemp[k-1]));
    partgamma[, covdim*(k-1)+(1:covdim)]= xpand * ((picum[,k])*(1-picum[,k]));
  }

  partbeta=partgamma
  for (k in 1:ncat){
    if (k>1)
      partbeta[, covdim*(k-1)+(1:covdim)] = partgamma[, covdim*(k-1)+(1:covdim)] - partgamma[, covdim*(k-2)+(1:covdim)];
    Fish = Fish+ t(partbeta[, covdim*(k-1)+(1:covdim)] /pitemp[,k]) %*% (partbeta[, covdim*(k-1)+(1:covdim)]);
  }

  for (i in 1:nitemp){
    loglik= loglik + log(pitemp[i,ytemp[i]]);
    score=score+ t(partbeta[i, covdim*(ytemp[i]-1)+(1:covdim)]/pitemp[i,ytemp[i]]);
  }
  #for (i in 1:nitemp) {if (is.na(log(pitemp[i,ytemp[i]]))) print(i)}

  return(list(loglik=loglik, score=score, Fish=Fish))
}

