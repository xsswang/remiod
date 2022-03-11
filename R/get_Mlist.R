#' Prepare imputation-model-related information
#'
#' Internal function to extract information of imputation models.
#'
#' @param object object inheriting from class \code{remoid}.
#' @return a list include raw data, imputation models, model types, fixed effects,
#'         random effects if any, reference categories corresponding to categorical
#'         variables in models, and interaction terms.
#'

get_Mlist = function (object) {
  if (!(inherits(object, "remiod") | inherits(object,
                                              "remiod_errored")))
    errormsg("%s must be of class %s or %s.", dQuote("object"),
             dQuote("remiod"), dQuote("remiod_errored"))
  c(object[c("data", "models", "fixed", "random")],
    object$Mlist, list(M = object$data_list[paste0("M_",
                                                   names(object$Mlist$group_lvls))]))
}
