#' Apply Delta adjustment to Update JAGS MCMC outputs under MAR for Cumulative Logistic Model
#'
#' Internal function to obtain delta-adjusted MCMC from an MAR object.
#' @inheritParams commParams
#' @param object an object of class remiod
#' @param treatment the variable name of treatment. Reference level of treatment should be coded as 0.
#' @param delta specific value used for Delta adjustment, applicable only
#'              for method="delta".
#' @param start first iteration to be used.
#' @param end last iteration to be used.
#' @param thin thinning to be applied.
#' @param subset subset of parameters (columns of the mcmc object) to be used.
#' @param exclude_chains optional vector of numbers, indexing MCMC chains to be excluded from the output.
#' @param mess logical, should messages be displayed?
#' @param seed optional seed value.
#' @param ... optional arguments pass from main function.
#'
#' @importFrom utils getFromNamespace
#'
#' @return A matrix of MCMC samples with all monitored parameters.A subset of
#'   the MCMC sample can be selected using \code{start}, \code{end} and
#'   \code{thin}.
#'

clm_MI_delta <- function(object, treatment, delta=0, start=NULL, end=NULL, exclude_chains=NULL,
                         thin=NULL, subset=FALSE, ord_cov_dummy=TRUE, seed=NULL, mess=FALSE,...)
  {
  if (!inherits(object, "remiod")) stop("Use only with 'remiod' objects.")

  Mlist = get_Mlist(object)

  info_list = object$info_list

  data = Mlist$data
  if (!("pattern" %in% colnames(data))) data$pattern = get_pattern(Mlist = Mlist)
  if(length(which(colnames(data)=='row'))==0) data$row = 1:nrow(data)

  fixed = Mlist$fixed
  auxvars = Mlist$auxvars
  refs = Mlist$refs

  coef_list = object$coef_list

  ## dummy vars -->> original varname
  vdvlist = lapply(names(refs), function(dv) data.frame(var=dv, varname=attr(refs[[dv]],"dummies")) )
  vdv = do.call(rbind.data.frame, vdvlist)

  scale_pars <- do.call(rbind, unname(Mlist$scale_pars))
  if (!is.null(scale_pars)) {
    scale_pars$center[is.na(scale_pars$center)] <- 0L
  }

  MCMC <- prep_MCMC(object, start = start, end = end, thin = thin,
                    subset = subset, exclude_chains = exclude_chains,
                    mess = mess)

  #mcparms = attr(MCMC,"dimnames")[[2]]

  rawdt = subset(data, select=names(info_list))
  misind1 = as.data.frame(which(is.na(rawdt), arr.ind = TRUE))

  idt = subset(data, select=c("row","pattern",treatment))
  misind2 = merge(misind1, idt, by="row", all.x=TRUE)
  misind2$colrev = length(info_list)+1 - misind2$col

  mvard = data.frame(mvar = names(info_list), col=1:length(names(info_list)) )
  misind3 = merge(misind2, mvard, by="col",all.x=TRUE)
  ### keep only treated subject, and ignore intermittent missingness
  misind = subset(misind3, misind3[,treatment]==1 & colrev >= pattern)

  misind$eta_nam = paste0("eta_",misind$mvar,"[",misind$row,"]")
  misind$mc_imp_col_nam = paste0("M_lvlone[",misind$row,",",misind$col,"]")

  ## J2R
  MCMC_delta = MCMC
  vimp = rev(names(info_list)[names(info_list) %in% unique(misind$mvar)])

  for (j in 1:length(vimp)){
    varname = vimp[j]
    coefs <- coef_list[[vimp[j]]]
    coefs <- merge(coefs, vdv, by="varname", all.x=TRUE)
    coefs$var[is.na(coefs$var)] = coefs$varname

    misind_j = subset(misind, mvar==vimp[j])

    for (h in 1:nrow(misind_j)){
      eta_nam_h = misind_j$eta_nam[h]
      eta = MCMC_delta[,eta_nam_h] + delta  # - MCMC[,coefs$coef[coefs$var==treatment]]

      ## sequential imputation: previously imputed values need to be adjusted in subsequent model
      ## when J=1, covariates are all observed
      ## following codes for j>1
      misind_h0 = subset(misind, row == misind_j$row[h])
      misind_h = subset(misind_h0, mvar %in% unique(coefs$var))

      cvar = unique(coefs$var)[which(unique(coefs$var) %in% misind_h$mvar)]

      if (length(cvar)>0) {
        if (ord_cov_dummy){
          dsmat0 = lapply(cvar, function(vk){
            ## in J2R, treatment effect need to removed for subjects in treated arm
            contr_vk = attr(refs[[vk]],"contr_matrix")
            contr_vk = data.frame(lev=rownames(contr_vk), contr_vk)
            mcc = MCMC[, misind_h$mc_imp_col_nam[misind_h$mvar==vk],drop=F]
            colnames(mcc) = "lev"
            left_join(mcc,contr_vk, by="lev")[,-1]
          })
          dsmat0 = as.matrix(do.call(cbind.data.frame, dsmat0))
          if (any(is.na(dsmat0))) dsmat0[is.na(dsmat0)] <- 0

          dsmat1 = lapply(cvar, function(vk){
            contr_vk = attr(refs[[vk]],"contr_matrix")
            contr_vk = data.frame(lev=rownames(contr_vk), contr_vk)
            mcc = MCMC_delta[, misind_h$mc_imp_col_nam[misind_h$mvar==vk],drop=F]
            colnames(mcc) = "lev"
            left_join(mcc,contr_vk, by="lev")[,-1]
          })
          dsmat1 = as.matrix(do.call(cbind.data.frame, dsmat1))
          if (any(is.na(dsmat1))) dsmat1[is.na(dsmat1)] <- 0
        } else {
          misind_ho = misind_h[match(cvar, misind_h$mvar), , drop = FALSE] # ordering monitored imp cols
          dsmat0 = MCMC[, misind_ho$mc_imp_col_nam,drop=F]
          dsmat1 = MCMC_delta[, misind_ho$mc_imp_col_nam,drop=F]
        }

        coef_t = coefs$coef[coefs$var %in% cvar]

        eta = eta - apply(MCMC[,coef_t,drop=FALSE] * dsmat0, 1, sum) +
          apply(MCMC[,coef_t,drop=FALSE] * dsmat1, 1, sum)
      }
      if(!is.matrix(eta)) eta = matrix(eta, ncol=1)

      gammas <- lapply( grep(paste0("gamma_", varname), colnames(MCMC), value = TRUE),
        function(k)
          matrix(nrow = nrow(eta), ncol = ncol(eta), data = rep(MCMC[, k], ncol(eta)),byrow = FALSE)
        )

      # add the category specific intercepts to the linear predictor
      lp <- lapply(seq_along(gammas), function(k) { gammas[[k]] + eta })

      mat1 <- matrix(nrow = nrow(eta), ncol = ncol(eta), data = 1L)
      mat0 <- mat1 * 0L

      if (info_list[[varname]]$rev) {
        names(lp) <- paste0("logOdds(", varname, "<=", seq_along(lp), ")")
        pred <- rev(c(lapply(rev(lp), plogis), list(mat0)))

        probs <- lapply(seq_along(pred)[-1L], function(k) {
          minmax_mat(pred[[k]] - pred[[k - 1L]])
        })

        probs <- c(probs,
                   list(
                     1L - minmax_mat(
                       apply(array(dim = c(dim(probs[[1L]]), length(probs)),
                                   unlist(probs)), c(1L, 2L), sum)
                     ))
        )
      } else {
        names(lp) <- paste0("logOdds(", varname, ">", seq_along(lp), ")")
        pred <- c(lapply(lp, plogis), list(mat0))

        probs <- lapply(seq_along(pred)[-1L], function(k) {
          minmax_mat(pred[[k - 1L]] - pred[[k]])
        })

        probs <- c(list(
          1L - minmax_mat( apply(array(dim = c(dim(probs[[1L]]), length(probs)),
                           unlist(probs)), c(1L, 2L), sum)
          )),
          probs)
      }
      names(probs) <- paste0("P_", varname, "_", levels(data[, varname]))

      probs = as.data.frame(probs)
      set.seed(seed)
      class = numeric()
      for (i in 1:nrow(probs)) class[i]= which( rmultinom(1,1,prob=probs[i,])==1 )
      clasm = matrix(class,ncol=1)

      ## corresponding MCMC col with imputed values
      mc_imp_col_nam = paste0("M_lvlone[",misind_j$row[h],",",misind_j$col[h],"]")
      MCMC_delta[,mc_imp_col_nam] <- clasm
      }
    }

  return(MCMC_delta)
}

