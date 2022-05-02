#' Listing the sequence of models used for imputation
#' @param object an object of class `remiod`
#' @param details logical. Default is FALSE, where listing all models in formula format. If TRUE,
#'        details of each models will be presented.
#' @param print logical. Default is TRUE to print all imputation models or detailed imputation models.
#' @return a list of formula of imputation models. If \code{details=TRUE}, information on the conditional
#'         distributions of the covariates in each imputation models. Note: the sequence of conditional
#'         models together specifies the joint distribution.
#'
#' @examples
#'
#' \donttest{
#' # data(schizow)
#'
#' test = remiod(formula = y6 ~ tx + y0 + y1 + y3, data = schizow,
#'               trtvar = 'tx', algorithm = 'jags', method="MAR",
#'               ord_cov_dummy = FALSE, n.adapt = 10, n.chains = 1,
#'               n.iter = 10, thin = 2, warn = FALSE, seed = 1234)
#'
#' list.models(test)
#' }
#' @export

list.models <- function (object, details = FALSE, print = TRUE) {
  if (!inherits(object, "remiod") )
    errormsg("Use only with 'remiod' objects.\n")

  if (details) {
    obj = object$mc.mar
    attr(obj, "class") <- "JointAI"
    if (print) JointAI::list_models(obj)
  } else {
    info_list = object$mc.mar$info_list
    coef_list = object$mc.mar$coef_list
    Mlist = get_Mlist(object$mc.mar)
    refs = Mlist$refs

    ## dummy vars -->> original varname
    if (length(refs)>0){
      vdvlist = lapply(names(refs), function(dv) data.frame(var=dv, varname=attr(refs[[dv]],"dummies")) )
      vdv = do.call(rbind.data.frame, vdvlist)
    }
    vimp = rev(names(info_list))

    for (j in 1:length(vimp)){
      varname = vimp[j]
      coefs <- coef_list[[vimp[j]]]

      if (length(refs)>0) coefs <- merge(coefs, vdv, by="varname", all.x=TRUE)
      else coefs$var = coefs$varname

      coefs$var[is.na(coefs$var)] = coefs$varname[is.na(coefs$var)]
      coefs = subset(coefs, !grepl("intercept", var,ignore.case = TRUE))
      dtj = data.frame(Order=j, model_formula = paste0(varname," ~ ", paste(unique(coefs$var), collapse = " + ")))
      if (j==1) modelseq = dtj
      else modelseq = rbind(modelseq, dtj)
      #cat( )
      #cat("\n")
    }
    if (print) print(modelseq, right=FALSE)
  }
}

