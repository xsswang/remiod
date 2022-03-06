#' Obtain ordinal results based on log-odds (eta) and cut-offs (gamma) from MCMC
#' Internal function to obtain ordinal results based on log-odds (eta) and cut-offs (gamma) from MCMC.
#' @param MCMC an matrix of MCMC samples.
#' @param impvar a name of imputation variable.
#' @param eta_name a name of eta in mcmc samples
#' @param rev logical. Reverse order or not.
#' @param seed optional. A seed value for randomness.
#' @keywords internal
get_class = function(MCMC, impvar, eta_name, rev=FALSE, seed=NULL){
  eta = MCMC[,eta_name,drop=FALSE]
  gammas <- lapply( grep(paste0("gamma_", impvar), colnames(MCMC), value = TRUE),
                    function(k)
                      matrix(nrow = nrow(eta), ncol = ncol(eta), data = rep(MCMC[, k], ncol(eta)),byrow = FALSE)
  )

  # add the category specific intercepts to the linear predictor
  lp <- lapply(seq_along(gammas), function(k) { gammas[[k]] + eta })

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
