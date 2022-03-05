#' Update MCMC outputs with delta adjustment.
#'
#' Obtains delta adjustment to MCMC from an object of class 'remiod'.
#' @param object an object of class remiod
#' @param treatment the variable name of treatment. Reference level of treatment should be coded as 0.
#' @param delta specific value used for Delta adjustment, applicable only
#'              for method="delta".
#' @param start first iteration to be used
#' @param end last iteration to be used
#' @param thin thinning to be applied
#' @param subset subset of parameters (columns of the mcmc object) to be used
#' @param exclude_chains optional vector of numbers, indexing MCMC chains to be excluded from the output
#' @param mess logical, should messages be displayed?
#' @param seed optional seed value.
#' @param ... optional arguments pass from main function.
#'
#' @keywords internal
#'
#' @return A matrix of MCMC samples with all monitored parameters.A subset of
#'   the MCMC sample can be selected using \code{start}, \code{end} and
#'   \code{thin}.
#'

glm_MI_delta <- function(object, treatment, delta=0, start=NULL, end=NULL, thin=NULL,
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
  vdvlist = lapply(names(refs), function(dv) data.frame(var=dv, varname=attr(refs[[dv]],"dummies")) )
  vdv = do.call(rbind.data.frame, vdvlist)

  MCMC <- prep_MCMC(object, start = start, end = end, thin = thin,
                    subset = subset, exclude_chains = exclude_chains,
                    warn = warn, mess = mess)

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

  misind$mu_nam = paste0("mu_",misind$mvar,"[",misind$row,"]")
  misind$mc_imp_col_nam = paste0("M_lvlone[",misind$row,",",misind$col,"]")

  ## delta adjustment
  MCMC_delta = MCMC
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
    coefs <- merge(coefs, vdv, by="varname", all.x=TRUE)
    coefs$var[is.na(coefs$var)] = "inter"

    misind_j = subset(misind, mvar==vimp[j])

    for (h in 1:nrow(misind_j)){
      mu_nam_h = misind_j$mu_nam[h]
      mu = MCMC_delta[,mu_nam_h] + matrix(rep(delta, nrow(MCMC_delta)),ncol=1)

      ## for those not monitored M_lvlone
      mcparms = attr(MCMC,"dimnames")[[2]]
      mcimp = misind_j$mc_imp_col_nam[h] %in% mcparms
      if (!mcimp){
        if (info_list[[varname]]$family == "binomial"){
          mu_inv = get_bin(mu = MCMC[, mu_nam_h], linkinvf=linkinv, seed=seed)
        } else {
          mu_inv = linkinv(MCMC[,misind_h$mu_nam[h]])
        }
        MCMC = cbind(MCMC, mu_inv)
        colnames(MCMC)[ncol(MCMC)] = misind_j$mc_imp_col_nam[h]

        MCMC_delta = cbind(MCMC_delta, mu_inv)
        colnames(MCMC_delta)[ncol(MCMC_delta)] = misind_j$mc_imp_col_nam[h]
      }

      ## sequential imputation: previously imputed values need to be adjusted in subsequent model
      ## when J=1, covariates are all observed
      ## following codes for j>1
      misind_h0 = subset(misind, row == misind_j$row[h])
      misind_h = subset(misind_h0, mvar %in% unique(coefs$var))

      ## the order of coefs$var must be maintained
      cvars = unique(coefs$var)
      cvar = cvars[which(cvars %in% misind_h$mvar)]

      coef_t = coefs$coef[coefs$var %in% cvar]

      if (nrow(misind_h) > 0){
        misind_h = misind_h[match(cvar, misind_h$mvar),,drop=FALSE]

        mc_old = subset(MCMC, select=misind_h$mc_imp_col_nam)
        mu_old = apply(mc_old * MCMC[,coef_t,drop=FALSE],1,sum)

        mc_new = subset(MCMC_delta, select=misind_h$mc_imp_col_nam)
        mu_new = apply(mc_new * MCMC[,coef_t,drop=FALSE],1,sum)
        mu = mu - matrix(mu_old, ncol=1) + matrix(mu_new, ncol=1)
      }

      if (info_list[[varname]]$family == "binomial")
        predc = get_bin(mu, linkinvf=linkinv, seed=seed)
      else predc = linkinv(mu)

      ## corresponding MCMC col with imputed values
      mc_imp_col_nam = misind_j$mc_imp_col_nam[h]
      if (!mcimp) {
        MCMC_delta = cbind(MCMC_delta, predc)
        colnames(MCMC_delta)[ncol(MCMC_delta)] = mc_imp_col_nam
        } else MCMC_delta[,mc_imp_col_nam] <- predc
      }
    }

  class(MCMC_delta) <- c(class(MCMC),"glm")
  return(MCMC_delta)
}

