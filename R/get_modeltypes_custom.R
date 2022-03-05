#' Specify the default (imputation) model types
#'
#' Modified from codes in JointAI to allow users customizing the orders of sequential models
#' Original code set the order of models based on number of missing values from each variable
#' Now, the order of models can be specified through \code{model_order}.
#'
#' @param fixed a two sided formula describing the fixed-effects part of the
#'              model (see \code{\link[stats]{formula}})
#' @param random only for multi-level models:
#'               a one-sided formula of the form \code{~x1 + ... + xn | g},
#'               where \code{x1 + ... + xn} specifies the model for the random
#'               effects and \code{g} the grouping variable
#' @param data a \code{data.frame} containing the original data
#'             (more details below)
#' @param auxvars optional; one-sided formula of variables that should be used
#'                as predictors in the imputation procedure (and will be imputed
#'                if necessary) but are not part of the analysis model(s).
#' @param timevar name of the variable indicating the time of the measurement
#'                of a time-varying covariate in a proportional hazards survival mode.
#' @param models optional; named vector specifying the types of models for
#'               (incomplete) covariates.
#'               This arguments replaces the argument \code{meth} used in
#'               earlier versions.
#'               If \code{NULL} (default) models will be determined
#'               automatically based on the class of the respective columns of
#'               \code{data}.
#' @param model_order optional; manually specify an order for imputation models.
#' @param no_model optional; vector of names of variables for which no model
#'                 should be specified. Note that this is only possible for
#'                 completely observed variables and implies the assumptions
#'                 of independence between the excluded variable and the
#'                 incomplete variables.
#' @param analysis_type type of analysis, e.g. \code{clm} for ordinal variable.
#' @param warn logical, should warnings be displayed?
#'
#' @note
#' the following scenario gives the sequential imputation models in an order of
#' y0, y1, y3, y2, y4,..., which is based on 'nmis' variable. The expected order
#' would be y0, y1, y2, y3, y4,.... To specify the order, what we need is just set
#' model_order = c('y0', 'y1', 'y2', 'y3', 'y4', 'y5')
#'
#'    rowID   out lvl nmis nlev ordered type L1
#' 8 rowID1 FALSE   1  428    3    TRUE  clm y5
#' 7 rowID1 FALSE   1  426    4    TRUE  clm y4
#' 5 rowID1 FALSE   1  423    4    TRUE  clm y2
#' 6 rowID1 FALSE   1   63    4    TRUE  clm y3
#' 4 rowID1 FALSE   1   11    4    TRUE  clm y1
#' 3 rowID1 FALSE   1    3    4    TRUE  clm y0
#'
#' @import survival


get_models_custom <- function(fixed, random = NULL, data, auxvars = NULL,
                       timevar = NULL, no_model = NULL, models = NULL, model_order = NULL,
                       analysis_type = NULL, warn = TRUE) {

  if (missing(fixed))
    errormsg("No formula specified.")

  if (missing(data))
    errormsg("No dataset given.")

  if (!is.null(auxvars) & class(auxvars) != 'formula')
    errormsg("The argument %s should be a formula.", dQuote("auxvars"))

  models_user <- models

  if (is.null(attr(fixed[[1]], 'type')))
    fixed <- extract_outcome_data(fixed, random = random, data = data,
                                  analysis_type = analysis_type,
                                  warn = FALSE)$fixed


  # if there is a time variable, add it to no_model
  if (!is.null(timevar)) {
    no_model <- c(no_model, timevar)
  }

  # check that all variables are found in the data
  allvars <- unique(c(all_vars(c(fixed, random, auxvars)), timevar))

  if (any(!names(models) %in% names(data))) {
    errormsg("Variable(s) %s were not found in the data." ,
             paste_and(dQuote(names(models)[!names(models) %in% names(data)])))
  }


  if (!is.null(no_model) &&
      any(colSums(is.na(data[, no_model, drop = FALSE])) > 0)) {
    errormsg("Variable(s) %s have missing values and imputation models are
             needed for these variables." ,
             paste(dQuote(no_model[colSums(is.na(data[, no_model,
                                                      drop = FALSE])) > 0]),
                   collapse = ", "))
  }

  if (any(!names(models_user) %in% allvars)) {
    errormsg("You have specified covariate model types for the variable(s) %s
             which are not part of the model.",
             paste_and(dQuote(names(models_user)[
               !names(models_user) %in% allvars])))
  }



  # extract the id variable from the random effects formula and get groups
  idvar <- extract_id(random, warn = warn)
  groups <- get_groups(idvar, data)

  random2 <- remove_grouping(random)


  # new version of allvars, without the grouping variable
  allvars <- unique(c(names(fixed),
                      all_vars(c(remove_lhs(fixed), random2, auxvars)),
                      names(models), timevar))

  group_lvls <- colSums(!identify_level_relations(groups))
  max_lvl <- max(group_lvls)

  if (length(allvars) > 0) {

    varinfo <- sapply(allvars, function(k) {
      x <- eval(parse(text = k), envir = data)
      out <- k %in% names(fixed)
      lvl <- group_lvls[
        check_varlevel(x, groups, group_lvls = identify_level_relations(groups))]
      nmis <- sum(is.na(x[match(unique(groups[[names(lvl)]]),
                                groups[[names(lvl)]])]))
      nlev <- length(levels(x))
      ordered <- is.ordered(x)
      data.frame(out = out, lvl = lvl, nmis = nmis, nlev = nlev,
                 ordered = ordered, type = NA)
    }, simplify = FALSE)

    varinfo <- melt_data.frame_list(varinfo, id.vars = colnames(varinfo[[1]]))



    varinfo$type[!varinfo$lvl %in% max_lvl & varinfo$nlev > 2 &
                   varinfo$ordered] <- "clmm"
    varinfo$type[varinfo$lvl %in% max_lvl & varinfo$nlev > 2 &
                   varinfo$ordered] <- "clm"
    varinfo$type[!varinfo$lvl %in% max_lvl & varinfo$nlev > 2 &
                   !varinfo$ordered] <- "mlogitmm"
    varinfo$type[varinfo$lvl %in% max_lvl & varinfo$nlev > 2 &
                   !varinfo$ordered] <- "mlogit"
    varinfo$type[!varinfo$lvl %in% max_lvl & varinfo$nlev == 2] <-
      "glmm_binomial_logit"
    varinfo$type[varinfo$lvl %in% max_lvl & varinfo$nlev == 2] <-
      "glm_binomial_logit"
    varinfo$type[!varinfo$lvl %in% max_lvl & varinfo$nlev == 0] <- "lmm"
    varinfo$type[varinfo$lvl %in% max_lvl & varinfo$nlev == 0] <- "lm"

    survmods <- sapply(fixed, 'attr', 'type') %in% c('coxph', 'survreg', 'JM')
    if (any(survmods)) {
      varinfo$type[varinfo$L1 %in% names(fixed[survmods])] <-
        sapply(fixed[survmods], 'attr', 'type')
    }

    if (!is.null(attr(fixed[[1]], 'type')))
      varinfo$type[varinfo$L1 %in% names(fixed)[1]] <- attr(fixed[[1]], 'type')

    varinfo <- varinfo[which(!varinfo$L1 %in% no_model), , drop = FALSE]


    types <- split(varinfo,
                   ifelse(varinfo$out, 'outcome',
                          ifelse(varinfo$nmis > 0,
                                 paste0('incomplete_lvl', varinfo$lvl),
                                 paste0('complete_lvl', varinfo$lvl)
                          )))


    types[which(names(types) != 'outcome')] <-
      lapply(types[which(names(types) != 'outcome')],
             function(x){
              if (is.null(model_order)) x[order(-x$lvl, x$nmis, decreasing = TRUE), , drop = FALSE]
              else x[match(rev(model_order), x$L1), , drop = FALSE]
      })



    NA_lvls <- unique(varinfo$lvl[varinfo$nmis > 0])


    models <- do.call(rbind,
                      c(types['outcome'],
                        if (any(!varinfo$out) & length(NA_lvls) > 0)
                        lapply(1:max(NA_lvls), function(k) {
                          set <- if (k == max(NA_lvls)) {
                            c('incomplete_lvl')
                          } else {
                            c('incomplete_lvl', 'complete_lvl')
                          }
                          do.call(rbind, types[paste0(set, k)])
                        })
                      ))

    models <- unlist(setNames(models$type, models$L1))


    models[names(models_user)] <- models_user

  } else {
    models <- NULL
  }
  models
}
