#' Obtain ordinal results based on log-odds (eta) and cut-offs (gamma) from MCMC
#' Internal function to obtain ordinal results based on log-odds (eta) and cut-offs (gamma) from MCMC.
#' @param MCMC an matrix of MCMC samples.
#' @param impvar a name of imputation variable.
#' @param eta_name a name of eta in mcmc samples
#' @param rev logical. Reverse order or not.
#' @param seed optional. A seed value for randomness.
#' @param shift shift of the thresholds caused by centring covariates (see
#'   `threshold_shift()`); a vector with one value per row of `MCMC`, or 0.
#' @keywords internal
get_class = function(MCMC, impvar, eta_name, rev=FALSE, seed=NULL, shift=0){
  eta = MCMC[,eta_name,drop=FALSE]
  gammas <- lapply( grep(paste0("gamma_", impvar), colnames(MCMC), value = TRUE),
                    function(k)
                      matrix(nrow = nrow(eta), ncol = ncol(eta), data = rep(MCMC[, k], ncol(eta)),byrow = FALSE)
  )

  # add the category specific intercepts to the linear predictor;
  # thresholds are on the original scale, eta on the scale used in JAGS
  lp <- lapply(seq_along(gammas), function(k) { gammas[[k]] + shift + eta })

  mat1 <- matrix(nrow = nrow(eta), ncol = ncol(eta), data = 1L)
  mat0 <- mat1 * 0L

  if (rev) {
    names(lp) <- paste0("logOdds(", impvar, "<=", seq_along(lp), ")")
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
    names(lp) <- paste0("logOdds(", impvar, ">", seq_along(lp), ")")
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
  #names(probs) <- paste0("P_", impvar, "_", levels(data[, impvar]))

  probs = as.data.frame(probs)
  if (!is.null(seed)) set.seed(seed)
  class = numeric()
  for (i in 1:nrow(probs)) class[i]= which( rmultinom(1,1,prob=probs[i,])==1 )
  #class = apply(probs,1, function(x) which(rmultinom(1,1,prob=probs[x,])==1))
  clasm = matrix(class,ncol=1)
  return(clasm)
}


#' Shift of the thresholds caused by centring covariates
#'
#' Internal function. Covariates named in `scale_vars` are centred and scaled
#' before they are passed to JAGS, so the JAGS node `eta_<y>` is the linear
#' predictor of the scaled covariates. Since JointAI 1.1.0.9000, the
#' thresholds `gamma_<y>` in the MCMC sample are transformed back to the
#' original scale of the covariates:
#' `gamma = gamma_scaled - sum(beta * center)`, where `beta` are the
#' coefficients on the original scale. This function returns
#' `sum(beta * center)`, which has to be added to `gamma + eta`. With older
#' versions of JointAI, or when no covariate is scaled, it returns 0.
#'
#' Only proportional-odds coefficients are used (non-proportional effects are
#' not used in the imputation either).
#'
#' @param MCMC a matrix of MCMC samples (on the original scale).
#' @param coefs a data.frame of the regression coefficients of the ordinal
#'   model (an element of `object$coef_list`; columns `coef` and `varname`).
#' @param scale_pars a data.frame of scaling parameters with column `center`
#'   (centres, 0 or `NA` for columns that are not centred), with one row per
#'   design-matrix column, named by the column; `NULL` if nothing is scaled.
#' @return A vector with one value per row of `MCMC`.
#' @keywords internal
threshold_shift = function(MCMC, coefs, scale_pars) {
  shift = rep(0, nrow(MCMC))

  # up to JointAI 1.1.0, the thresholds were not transformed back
  if (is.null(scale_pars) || utils::packageVersion("JointAI") <= "1.1.0") {
    return(shift)
  }

  if (!is.null(coefs$outcat)) coefs = coefs[is.na(coefs$outcat), , drop = FALSE]

  center = scale_pars[coefs$varname, "center"]
  center[is.na(center)] = 0

  for (k in which(center != 0)) {
    shift = shift + MCMC[, coefs$coef[k]] * center[k]
  }
  shift
}
