#' Create multiple imputed datasets based on assigned imputation method.
#'
#' Internal function, creates multiple imputed datasets based on assigned
#' imputation method returns multiple imputed datasets stacked
#' onto each other (i.e., long format; optionally including the original,
#' incomplete data).\cr
#'
#' @param object an object of class JointAI
#' @param treatment the variable name of treatment. Reference level of treatment should be coded as 0.
#' @param method a method for obtaining multiple-imputed dataset. Options
#'               include MAR, J2R, CR, and Delta adjustment.
#' @param delta specific value used for Delta adjustment, applicable only
#'              for method="delta".
#' @param include should the original, incomplete data be included? Default is
#'                \code{TRUE}.
#' @param start first iteration to be used.
#' @param end last iteration to be used.
#' @param thin thinning to be applied.
#' @param subset subset of parameters (columns of the mcmc object) to be used.
#' @param exclude_chains optional vector of numbers, indexing MCMC chains to be excluded from the output.
#' @param mess logical, should messages be displayed?
#' @param seed optional seed value.
#' @param ... optional arguments pass from main function.
#'
#' @return A \code{data.frame} in which the original data (if
#'         \code{include = TRUE}) and the imputed datasets are stacked onto
#'         each other.\cr
#'         The variable \code{Imputation_} indexes the imputation, while
#'         \code{.rownr} links the rows to the rows of the original data.
#'         In cross-sectional datasets the
#'         variable \code{.id} is added as subject identifier.
#'
#'
get_MI_RB <- function(object, treatment, method=c("MAR","J2R","CR","delta"), delta=0,
                      exclude_chains=NULL, start=NULL, end=NULL, seed=NULL,
                      thin=NULL, subset=FALSE, include=TRUE, mess=TRUE, ...)
  {
  if(!missing(method) & length(method)>1) stop("Only one 'method' allowed.")
  method <- match.arg(method)

  oldseed <- .Random.seed
  on.exit({
    .Random.seed <<- oldseed
  })

  # set seed value if provided
  if (!is.null(seed)) {
    set_seed(seed)
  }

  # extract original data and add
  # - column with row numbers (needed for plot_imp_distr())
  # - an id variable if there is none
  DF <- object$data
  DF$.rownr <- seq_len(nrow(DF))
  if (length(object$Mlist$groups) < 2) DF$.id <- seq_len(nrow(DF))

  # for keeping raw missingness
  DF0 = DF

  # extract variable levels
  Mlvls <- object$Mlist$Mlvls

  # names of variables that were imputed
  vars <- intersect(names(object$models), names(DF)[colSums(is.na(DF)) > 0])

  # get a summary of the relevant characteristics of the imputed variables
  varinfo <- lapply(object$info_list[vars], function(x) {
    data.frame(varname = x$varname,
               modeltype = x$modeltype,
               family = ifelse(!is.null(x$family), x$family, NA),
               stringsAsFactors = FALSE)
  })

  if (varinfo[[1]]$modeltype == 'clm'){
      mcUpdateFun = switch(method,
                           'MAR' = prep_MCMC,
                           'J2R' = clm_MI_J2R,
                           'CR' = clm_MI_CR,
                           'delta' = clm_MI_delta)
  } else {
    mcUpdateFun = switch(method,
                         'MAR' = prep_MCMC,
                         'J2R' = glm_MI_J2R,
                         'CR' = glm_MI_CR,
                         'delta' = glm_MI_delta)
    }

  MCMC = mcUpdateFun(object, treatment=treatment, delta=delta, seed = seed,
                     start = start, end = end, thin = thin,
                     subset = subset, exclude_chains = exclude_chains,
                     mess = mess)

  # prepare a list of copies of the original data
  df_list <- list()
  df_list[[1]] <- cbind("Imputation_" = 0, DF)

  # if (method == "J2R"){
  #   ### change rows with missing outcome to all missing
  #   respmis = which(is.na(DF[,names(object$fixed)]) & DF[,treatment]==1)
  #   DF[respmis,names(object$models)] = NA
  # }

  m = nrow(MCMC)
  for (i in 2:(m + 1)) {
    df_list[[i]] <- cbind("Imputation_" = i - 1, DF)
  }

  for (i in vars) {
    impval <- NULL

    # identify the names of the columns in MCMC corresponding to variable i
    pat <- paste0(Mlvls[i], "\\[[[:digit:]]*,",
                  match(i, colnames(object$data_list[[Mlvls[i]]])),
                  "\\]")

    if (!any(grepl(pat, colnames(MCMC))))
      errormsg("I cannot find imputed values for %s. Did you monitor them?",
               dQuote(i))

    impval <- MCMC[, grep(pat, colnames(MCMC), value = TRUE), drop = FALSE]

    if (length(impval) > 0) {
      rownrs <- gsub(",[[:digit:]]*\\]", "",
                     gsub("^[[:print:]]*\\[", "", colnames(impval)))

      for (j in (1:m) + 1) {
        iv <- impval[j - 1, na.omit(match(
          object$Mlist$groups[[gsub("M_", "", Mlvls[i])]],
          as.numeric(rownrs)
        ))]

        if (is.factor(df_list[[j]][, i])) {
          df_list[[j]][is.na(df_list[[j]][, i]), i] <-
            factor(iv, labels = levels(df_list[[j]][, i]),
                   levels = seq_along(levels(df_list[[j]][, i])) -
                     as.numeric(length(levels(df_list[[j]][, i])) == 2)
                   )
        } else {
          df_list[[j]][is.na(df_list[[j]][, i]), i] <- iv
        }
      }
    }
  }

  # if (method == 'J2R'){
  #   for (i in 1:length(df_list)) {
  #     DFi = DF0
  #     dfi = df_list[[i]]
  #     DFi[is.na(DFi)] = dfi[is.na(DFi)]
  #     df_list[[i]] <- DFi
  #   }
  # }

  if (!include) df_list <- df_list[-1]

  if (method == 'delta'){
    for (i in 1:length(df_list)) {
      df_list[[i]] <- cbind(df_list[[i]], "delta" = delta)
    }
  }

  # build dataset --------------------------------------------------------------
  imp_df <- data.table::rbindlist(df_list) #do.call(rbind, df_list)

  return(imp_df)
}
