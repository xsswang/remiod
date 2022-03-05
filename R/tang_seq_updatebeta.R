
updatebeta = function(object, beta, datat, ord_cov_dummy, rinv=0.00001, scheme=2){
  Mlist = get_Mlist(object)
  coefs = object$coef_list
  infolist = object$info_list
  Xmat = dsmat(Mlist=Mlist, ord_cov_dummy=ord_cov_dummy, data=datat)

  yv = rev(names(coefs))
  if (Mlist$models[1]=="clm") ymat = dsmat(Mlist=Mlist, ord_cov_dummy=FALSE, data=datat)
  else ymat = prep_covoutcomes(datat[,yv])

  minpt = min(datat$pattern)
  maxpt = max(datat$pattern)
  if (minpt==0) pat = 1:maxpt
  else pat = minpt:maxpt

  accept=0

  for (i in 1:length(yv)){

    ncat = infolist[[yv[i]]]$ncat
    if (ncat>2) {
      selvar = c("(Intercept)", coefs[[yv[i]]]$varname)
      ytemp = ymat[datat$pattern >= pat[i], yv[i]]
    } else {
      selvar = coefs[[yv[i]]]$varname
      ytemp = ymat[datat$pattern >= pat[i], grep(yv[i],colnames(ymat))]
    }

    xtemp = subset(Xmat[datat$pattern >= pat[i],], select=selvar)

    covdim = ncat-2+ncol(xtemp);

    if (ncat>2) betatemp = unlist(c(beta[['c']][[yv[i]]], beta[['b']][[yv[i]]])) #betatemp;
    else betatemp = unlist(beta[['b']][[yv[i]]])
    #betatempnew = betatemp;

    liks = get_likelihood(xtemp, ytemp, betai=beta[['b']][[yv[i]]],
                          cuti=beta[['c']][[yv[i]]], ncat=ncat)

    loglik = liks$loglik
    fish = liks$Fish

    Lx = chol(fish + (rinv)*diag(covdim))
    # print("the following is Lx")
    # print(Lx)
    # cat("\n")

    if (scheme == 1) {
      m1 = backsolve(t(Lx), t(liks$score))
      mtemp = m1+ rnorm(covdim)
      mean  = backsolve(Lx, m1);
      move  = backsolve(Lx,mtemp);
    }else{
      mtemp = rnorm(covdim)
      move  = backsolve(Lx,mtemp);
    }
    # print("the following is move")
    # print(move)
    # cat("\n")
    betatempnew = betatemp  + (move)

    # print("the following is betatemp")
    # print(betatemp)
    # cat("\n")

    if (scheme == 1) {
      log_dem = -sum((betatemp)^2)*(rinv)/2+ loglik + sum(log(diag(Lx))) - t(mean-move) %*% (t(Lx) %*% Lx) %*% (mean-move)/2;
    }else {
      log_dem = -sum((betatemp)^2)*(rinv)/2+ loglik + sum(log(diag(Lx))) - t(move)%*% (t(Lx) %*% Lx) %*% move/2;
    }
    # print("the following is dem")
    # print(log_dem)
    # cat("\n")

    if (ncat>2) {
      cutnew= cumsum(exp( unlist(betatempnew[1:(ncat-2)])) )

      betainew = betatempnew[-c(1:(ncat-2))]
      cutinew  = betatempnew[c(1:(ncat-2))]
    } else {
      betainew = betatempnew
      cutinew = NULL
      cutnew = NULL
    }

    liks2 = get_likelihood(xtemp, ytemp, betai=betainew, cuti=cutinew, ncat=ncat)
    loglik2 = liks2$loglik
    fish2 = as.matrix(liks2$Fish)

    if (length(which(is.na(fish2)))>0) {
      next
      # beta[i,]=beta[i,];
      # betatemp[i,]=betatemp[i,];
    } else {
      egval = eigen(fish2, symmetric=TRUE)$values
      if (min(egval) > 0) Lx2 = chol(fish2 + (rinv)*diag(covdim))
      else {
        pd = Matrix::nearPD(fish2, corr=TRUE)
        Lx2 = chol(as.matrix(pd$mat) + (rinv)*diag(covdim))
      }

      if (scheme == 1) {
        m2=backsolve(Lx2, t(liks2$score));
        mean2= backsolve(Lx2,m2);
        log_num = -sum((betatempnew)^2)*(rinv)/2+ loglik2 + sum(log(diag(Lx))) - t(mean2-move) %*% (t(Lx) %*% Lx) %*% (mean2-move)/2;
      } else {
        log_num = -sum((betatempnew)^2)*(rinv)/2+ loglik2 + sum(log(diag(Lx))) - t(move)%*% (t(Lx) %*% Lx) %*% move/2;
      }

      t = runif(1)
      # print("the following is t")
      # print(t)
      # cat("\n")

      if (log(t)< log_num -log_dem) {
        #print("Accepted")
        accept = accept + 1;
        #betacexp = betacexpnew
        betatemp = betatempnew

        beta[['b']][[yv[i]]] = betainew  #betatempnew[-c(1:(ncat-2))]
        beta[['c']][[yv[i]]] = cutinew   #betatempnew[c(1:(ncat-2))]
        beta[['cexps']][[yv[i]]] = cutnew
      }
   }
  }
  return(list(beta=beta, accept=accept))
}
