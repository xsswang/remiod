#' Apply Jump-to-Reference(J2R) Method to Update JAGS MCMC outputs under MAR for Generalized Linear Model
#'
#' Internal function to obtain Jump-to-Reference(J2R) MCMC from an MAR object.
#' @param object an object of class remiod
#' @param treatment the variable name of treatment. Reference level of treatment should be coded as 0.
#' @param start first iteration to be used
#' @param end last iteration to be used
#' @param thin thinning to be applied
#' @param subset subset of parameters (columns of the mcmc object) to be used
#' @param exclude_chains optional vector of numbers, indexing MCMC chains to be excluded from the output
#' @param mess logical, should messages be displayed?
#' @param seed optional seed value.
#' @param ... optional arguments pass from main function.
#'
#' @return A matrix of MCMC samples with all monitored parameters.A subset of
#'   the MCMC sample can be selected using \code{start}, \code{end} and
#'   \code{thin}.
#'

glm_MI_J2R <- function(object, treatment, start=NULL, end=NULL, thin=NULL,
                      exclude_chains=NULL, subset=FALSE, seed=5432, mess=FALSE,...)
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

  # X <- model_matrix_combi(fmla = c(fixed, auxvars), data = data, refs = refs,
  #                         terms_list = get_terms_list(fmla = c(fixed, auxvars), data = data))

  coef_list = object$coef_list

  ## dummy vars -->> original varname
  if (length(refs) > 0){
    vdvlist = lapply(names(refs), function(dv) data.frame(var=dv, varname=attr(refs[[dv]],"dummies")) )
    vdv = do.call(rbind.data.frame, vdvlist)
  }

  MCMC <- prep_MCMC(object, start = start, end = end, thin = thin,
                    subset = subset, exclude_chains = exclude_chains,
                    warn = warn, mess = mess)

  mcpm_raw = attr(MCMC,"dimnames")[[2]]

  rawdt = subset(data, select=names(info_list))

  ### change rows with missing outcome to all missing
  respmis = which(is.na(rawdt[,names(fixed)]) & data[,treatment]==1)
  rawdt[respmis,] = NA

  misind0 = as.data.frame(which(is.na(rawdt), arr.ind = TRUE))
  misind1 = subset(misind0, row %in% respmis)

  idt = subset(data, select=c("row","pattern",treatment))
  misind2 = merge(misind1, idt, by="row", all.x=TRUE)
  misind2$colrev = length(info_list)+1 - misind2$col

  mvard = data.frame(mvar = names(info_list), col=1:length(names(info_list)) )
  misind3 = merge(misind2, mvard, by="col",all.x=TRUE)
  ### keep only treated subject, and ignore intermittent missingness
  misind4 = subset(misind3, misind3[,treatment]==1)
  misind = subset(misind4, row %in% respmis | (!(row %in% respmis) & colrev >= pattern))

  misind$mu_nam = paste0("mu_",misind$mvar,"[",misind$row,"]")
  misind$mc_imp_col_nam = paste0("M_lvlone[",misind$row,",",misind$col,"]")

  ## J2R
  MCMC_J2R = MCMC
  vimp = rev(names(info_list)[names(info_list) %in% unique(misind$mvar)])

  for (j in 1:length(vimp)){
    varname = vimp[j]

    linkinv <- if (info_list[[varname]]$family %in%
                   c("gaussian", "binomial", "Gamma", "poisson")) {
                    get(info_list[[varname]]$family)(link = info_list[[varname]]$link)$linkinv
                  } else if (info_list[[varname]]$family %in% "lognorm") {
                    gaussian(link = "log")$linkinv
                  } else if (info_list[[varname]]$family %in% "beta") {
                    plogis
                  }

    coefs <- coef_list[[vimp[j]]]
    if (nrow(vdv)>0) {
      coefs <- merge(coefs, vdv, by="varname", all.x=TRUE)
      coefs$var[is.na(coefs$var)] = "inter"
    } else {
      coefs$var = coefs$varname
    }

    misind_j = subset(misind, mvar==vimp[j])

    for (h in 1:nrow(misind_j)){
      mu_nam_h = misind_j$mu_nam[h]
      mcparms = attr(MCMC,"dimnames")[[2]]
      if (mu_nam_h %in% mcparms) mu = MCMC[,mu_nam_h] - MCMC[,coefs$coef[coefs$varname==treatment]]
      else browser()

      ## for those not monitored M_lvlone
      mcimp = misind_j$mc_imp_col_nam[h] %in% mcparms
      if (!mcimp){

          if (info_list[[varname]]$family == "binomial"){
            mu_inv = get_bin(mu = MCMC[, mu_nam_h], linkinvf=linkinv, seed=seed)
          } else {
            mu_inv = linkinv(MCMC[,misind_j$mu_nam[h]])
          }
          MCMC = cbind(MCMC, mu_inv)
          colnames(MCMC)[ncol(MCMC)] = misind_j$mc_imp_col_nam[h]

          MCMC_J2R = cbind(MCMC_J2R, mu_inv)
          colnames(MCMC_J2R)[ncol(MCMC_J2R)] = misind_j$mc_imp_col_nam[h]
      }

      ## sequential imputation: previously imputed values need to be adjusted in subsequent model
      ## when J=1, covariates are all observed
      misind_h0 = subset(misind, row == misind_j$row[h])
      misind_h = subset(misind_h0, mvar %in% unique(coefs$var))

      ## the order of coefs$var must be maintained
      cvar = unique(coefs$var)[which(unique(coefs$var) %in% misind_h$mvar)]
      coef_t = coefs$coef[coefs$var %in% cvar]

      if (nrow(misind_h) > 0){
          misind_h = misind_h[match(cvar, misind_h$mvar),,drop=FALSE]

          mc_old = subset(MCMC, select=misind_h$mc_imp_col_nam)
          mu_old = apply(mc_old * MCMC[,coef_t,drop=FALSE],1,sum)

          mc_new = subset(MCMC_J2R, select=misind_h$mc_imp_col_nam)
          mu_new = apply(mc_new * MCMC[,coef_t,drop=FALSE],1,sum)
          mu = mu - matrix(mu_old, ncol=1) + matrix(mu_new, ncol=1)
      }

      if (info_list[[varname]]$family == "binomial")
        predc = get_bin(mu, linkinvf=linkinv, seed=seed)
      else predc = linkinv(mu)

      ## corresponding MCMC col with imputed values
      mc_imp_col_nam = misind_j$mc_imp_col_nam[h]
      if (!mcimp) {
        MCMC_J2R = cbind(MCMC_J2R, predc)
        colnames(MCMC_J2R)[ncol(MCMC_J2R)] = mc_imp_col_nam
      } else MCMC_J2R[,mc_imp_col_nam] <- predc
    }
  }

  MCMC_J2R_out = subset(MCMC_J2R, select=mcpm_raw)
  class(MCMC_J2R_out) <- c(class(MCMC),"glm")
  return(MCMC_J2R_out)
}

