#' Extract specific parameters from MCMC samples
#'
#' @param object an object of class `MCMC`.
#' @param subset subset of parameters (columns of the mcmc object) to be used. See
#'               https://nerler.github.io/JointAI/articles/SelectingParameters.html
#'               for key-words of subseting parameters. Besides, `selected_parms`
#'               and `selected_vars` are new key-words for arbitrarily selecting
#'               parameters.
#' @param warn logical, should warning messages be displayed?
#' @param mess logical, should messages be displayed?
#'
#' @examples
#' \donttest{
#' data(schizow)
#'
#' test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
#'               trtvar = 'tx', algorithm = 'jags', method="MAR",
#'               ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
#'               n.iter = 0, thin = 1, warn = FALSE, seed = 1234)
#'
#' pms = c("beta[2]","alpha[2]","alpha[6]","alpha[9]")
#' mcsub = get_subset(object = test$mc.mar, subset=c(selected_parms = list(pms)))
#'
#' }
#' @export

get_subset = function (object, subset, warn = TRUE, mess = TRUE) {
  if (identical(subset, FALSE))
    return(object$MCMC)

  if (!is.list(subset)) subset <- as.list(subset)

  if ("selected_vars" %in% names(subset)){
    coef_list = object$coef_list
    coef_list = data.table::rbindlist(coef_list)
    coef_var = subset(coef_list, varname == subset[['selected_vars']])
    subset=c(subset, selected_parms = list(coef_var$coef))
  }

  if (is.null(subset$selected_parms)){
    if (length(subset) == 0 & !as.logical(as.list(object$monitor_params)$analysis_main))
      return(object$MCMC)
    if (length(subset) == 0 & as.logical(as.list(object$monitor_params)$analysis_main))
      subset$analysis_main <- TRUE
    if (!isFALSE(subset$analysis_main)) {
      subset$analysis_main <- TRUE
    }
    Mlist_new <- get_Mlist(object)
    Mlist_new$ppc <- as.list(subset)$ppc
    s <- do.call(get_params, c(list(Mlist = Mlist_new, info_list = object$info_list),
                               subset, mess = mess))
    if (is.null(s)) {
      errormsg("You have selected an empty subset of parameters.")
    }
    sub <- unique(unlist(c(sapply(paste0("^", s, "\\["),
                                  grep, colnames(object$MCMC[[1]]), value = TRUE),
                           colnames(object$MCMC[[1]])[na.omit(sapply(s, match, table = colnames(object$MCMC[[1]])))])))
  } else sub = subset$selected_parms
  if (length(sub) == 0)
    sub <- colnames(object$MCMC[[1]])
  return(object$MCMC[, sub, drop = FALSE])
}
