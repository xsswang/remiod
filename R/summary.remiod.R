#' Summarize the results from an object of class remiod
#'
#' Obtain and print the \code{summary}, (fixed effects) coefficients
#' (\code{coef}) and credible interval (\code{confint}).
#'
#' @param outcome specify outcome variable to select imputation model(s) to summarize.
#'                Default generates summaries for all models.
#' @param digits the minimum number of significant digits to be printed in values.
#' @param quantiles posterior quantiles
#' @inheritParams commParams
#' @inheritParams JointAI::model_imp
#'
#' @return summary information, including parameter posterior mean, posterior SD,
#'         quantiles, tail probability \code{tail-prob},  Gelman-Rubin criterion
#'         \code{GR-crit}, the ratio of the Monte Carlo error and posterior standard
#'         deviation) for specified parameters \code{MCE/SD}.
#'
#' @examples
#' \donttest{
#' # data(schizow)
#'
#' test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
#'               trtvar = 'tx', algorithm = 'jags', method="MAR",
#'               ord_cov_dummy = FALSE, n.adapt = 50, n.chains = 1,
#'               n.iter = 50, thin = 2, warn = FALSE, seed = 1234)
#'
#' summary(object = test, outcome = c("y6","y3"))
#' }
#'
#' @name summary
#' @export

summary <- function(object, ...) {
  UseMethod("summary", object)
}

#' @rdname summary
#' @export
summary.remiod <- function(object, start = NULL, end = NULL, thin = NULL,
                            quantiles = c(0.025, 0.975), outcome = NULL,
                            exclude_chains = NULL, warn = TRUE, mess = TRUE, ...) {
  object = object$mc.mar
  if (is.null(object$MCMC)) errormsg("There is no MCMC sample.")

  cl <- as.list(match.call())[-1]
  autoburnin <- if (is.null(cl$autoburnin)) FALSE else eval(cl$autoburnin)

  MCMC <- prep_MCMC(object, start = start, end = end, thin = thin,
                    subset = FALSE, exclude_chains = exclude_chains,
                    warn = warn, mess = mess)

  # create results matrices
  statnames <- c("Mean", "SD", paste0(quantiles * 100, "%"), "tail-prob.",
                 "GR-crit", "MCE/SD")

  vars <- if (is.null(outcome)) {
    names(object$coef_list)
  } else {
    #names(object$coef_list[[outcome]])
    outcome
  }
  coeflist = object$coef_list

  res_list <- sapply(vars, function(varname) {

    if (object$info_list[[varname]]$modeltype %in% c("clm"))
      cuts = grep(paste0("gamma_", varname), colnames(MCMC),value = TRUE)
    else cuts = NULL

    modelvars = intersect(colnames(MCMC), c(object$coef_list[[varname]]$coef,cuts))
    MCMCsub <- MCMC[, modelvars, drop = FALSE]

    if (ncol(MCMCsub) > 0) {

      grcrit <- if (length(object$MCMC) - length(exclude_chains) > 1) {
        GR_crit(object = object, start = start, end = end, thin = thin,
                warn = warn, mess = FALSE, multivariate = FALSE,
                exclude_chains = exclude_chains,
                subset = list(selected_parms = modelvars),
                autoburnin = autoburnin)[[1]][, "Upper C.I."]
      }

      mcerror <- if (length(object$MCMC) - length(exclude_chains) > 1) {
        try(MC_error(MCMC=MCMCsub, digits = 2, warn = FALSE, mess = FALSE))
      }

      stats <- matrix(nrow = length(colnames(MCMCsub)),
                      ncol = length(statnames),
                      dimnames = list(colnames(MCMCsub), statnames))

      stats[, "Mean"] <- apply(MCMCsub, 2, mean)
      stats[, "SD"] <- apply(MCMCsub, 2, sd)
      stats[, paste0(quantiles * 100, "%")] <- t(apply(MCMCsub, 2,
                                                       quantile, quantiles))
      stats[, "tail-prob."] <- apply(MCMCsub, 2, computeP)

      if (length(object$MCMC) - length(exclude_chains) > 1)
        stats[, "GR-crit"] <- grcrit[modelvars]

      if (length(object$MCMC) - length(exclude_chains) > 1) {
        if (!inherits(mcerror, "try-error"))
          stats[, "MCE/SD"] <- mcerror$data_scale[, "MCSE/SD"]
      }

      sigma <- if (object$info_list[[varname]]$family %in%
                   c("gaussian", "Gamma", "lognorm") &&
                   !is.null(object$info_list[[varname]]$family)) {
        sig <- grep(paste0("sigma_", varname), rownames(stats))

        if (length(sig) > 0) {
          stats[sig, -which(colnames(stats) == "tail-prob."), drop = FALSE]
        }
      }

      modelvy = paste0("beta_",coeflist[[varname]]$outcome,"_",coeflist[[varname]]$varnam_print)
      modelvars[which(modelvars %in% coeflist[[varname]]$coef)] = modelvy
      rownames(stats) = modelvars

      list(regcoef=stats, modeltype = object$info_list[[varname]]$modeltype,
           family = object$info_list[[varname]]$family, sigma = sigma)
    }
  }, simplify = FALSE)


  out <- list()
  out$call <- object$call
  out$start <- ifelse(is.null(start), start(object$MCMC),
                      max(start, start(object$MCMC)))
  out$end <- ifelse(is.null(end), end(object$MCMC), min(end, end(object$MCMC)))
  out$thin <- coda::thin(object$MCMC)
  out$nchain <- coda::nchain(object$MCMC) - sum(exclude_chains %in%
                                                  seq_along(object$MCMC))
  out$res <- res_list
  out$outcome <- outcome

  out$analysis_type <- object$analysis_type
  out$size <- object$Mlist$N

  class(out) <- "summary.remiod"
  return(out)
}

#' print summary outputs
#' @rdname summary
#' @param x an object of class \code{summary.remiod}
#' @export
print.summary.remiod <- function(x, digits = 3, ...) {

  if (!inherits(x, "summary.remiod"))
    errormsg("Use only with objects.", sQuote("summary.remiod"))

  cat("\n")

  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n", sep = "")

  for (k in seq_along(x$res)) {
    if (!is.null(x$res[[k]])) {
      cat("\n\n")
      if (sum(!sapply(x$res, is.null)) > 1 | !is.null(x$outcome))
        cat(paste0(
          "# ", paste0(c(rep("-", 69)), collapse = ""), " #\n",
          "  ", "Bayesian ",
          print_type(x$res[[k]]$modeltype, x$res[[k]]$family), " for ",
          dQuote(names(x$res)[k]), "\n",
          "# ", paste0(c(rep("-", 35)), collapse = " "), " #\n\n"
        ))


      if (!is.null(x$res[[k]]$regcoef)) {
        cat("Posterior summary:\n")
        print(x$res[[k]]$regcoef, digits = digits, na.print = "")
      }

      if (!is.null(x$res[[k]]$sigma))  {
        cat("\nPosterior summary of residual std. deviation:\n")
        print(x$res[[k]]$sigma, digits = digits, na.print = "")
      }
    }
  }

  cat("\n\n")
  if (sum(!sapply(x$res, is.null)) > 1)
    cat("#", paste0(c(rep("-", 59)), collapse = ""), "#\n\n")

  cat("MCMC settings:\n")
  cat("Iterations = ", x$start, ":", x$end, "\n", sep = "")
  cat("Sample size per chain =", (x$end - x$start) / x$thin +
        1, "\n")
  cat("Thinning interval =", x$thin, "\n")
  cat("Number of chains =", x$nchain, "\n")
  cat("\n")
  cat("Number of observations:", x$size["lvlone"], "\n")
  if (length(x$size) > 1) {
    i <- which(!names(x$size) %in% "lvlone")
    cat("Number of groups:\n",
        paste0("- ", names(x$size)[i], ": ", x$size[i], "\n")
    )
  }

  invisible(x)
}


#' @rdname summary
#' @export
coef.summary.remiod <- function(object, start = NULL, end = NULL, thin = NULL,
                                 subset = NULL, exclude_chains = NULL,
                                 warn = TRUE, mess = TRUE, ...) {

  if (!inherits(object, "summary.remiod"))
    errormsg("Use only with %s objects.", sQuote("summary.remiod"))

  Filter(Negate(is.null),
         lapply(object$res, "[[", "regcoef")
  )
}
