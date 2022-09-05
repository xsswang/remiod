#' Model set-up
#' Code is modified from JointAI to allow for additional control to imputation models
#'
#' @param formula a two sided model formula (see \code{\link[stats]{formula}})
#'                or a list of such formulas; (more details below).
#' @param fixed a two sided formula describing the fixed-effects part of the
#'              model (see \code{\link[stats]{formula}})
#' @param random only for multi-level models:
#'               a one-sided formula of the form \code{~x1 + ... + xn | g},
#'               where \code{x1 + ... + xn} specifies the model for the random
#'               effects and \code{g} the grouping variable
#' @param data a \code{data.frame} containing the original data
#'             (more details below)
#' @param family only for \code{glm_imp} and \code{glmm_imp}/\code{glmer_imp}:
#'               a description of the distribution and link function to
#'               be used in the model. This can be a character string naming a
#'               family function, a family function or the result of a call to
#'               a family function. (For more details see below and
#'               \code{\link[stats]{family}}.)
#' @param rd_vcov character string or list specifying the structure of the
#'                random effects variance covariance matrix, see details below.
#' @param monitor_params named list or vector specifying which parameters
#'                       should be monitored (more details below)
#' @param inits optional; specification of initial values in the form of a list
#'              or a function (see \code{\link[rjags]{jags.model}}).
#'              If omitted, starting values for the random number generator are
#'              created by \strong{JointAI}, initial values are then generated
#'              by JAGS.
#'              It is an error to supply an initial value for an observed node.
#' @param auxvars optional; one-sided formula of variables that should be used
#'                as predictors in the imputation procedure (and will be imputed
#'                if necessary) but are not part of the analysis model(s).
#'                For more details with regards to the behaviour with
#'                non-linear effects see the vignette.
#' @param models optional; named vector specifying the types of models for
#'               (incomplete) covariates.
#'               This arguments replaces the argument \code{meth} used in
#'               earlier versions.
#'               If \code{NULL} (default) models will be determined
#'               automatically based on the class of the respective columns of
#'               \code{data}.
#' @param refcats optional; either one of \code{"first"}, \code{"last"},
#'                \code{"largest"} (which sets the category for all categorical
#'                variables) or a named list specifying which category should
#'                be used as reference category per categorical variable.
#'                Options are the category label, the category number,
#'                or one of "first" (the first category),
#'                "last" (the last category) or "largest" (chooses the category
#'                with the most observations).
#'                Default is "first". If reference categories are specified for
#'                a subset of the categorical variables the default will be
#'                used for the remaining variables.
#'                (See also \code{\link{set_refcat}})
#' @param shrinkage optional; either a character string naming the shrinkage
#'                  method to be used for regression coefficients in all models
#'                  or a named vector specifying the type of shrinkage to be
#'                  used in the models given as names.
#' @param rev optional character vector; vector of ordinal outcome variable
#'            names for which the odds should be reversed, i.e.,
#'            \eqn{logit(y\le k)} instead of \eqn{logit(y > k)}.
#' @param nonprop optional named list of one-sided formulas specifying
#'                covariates that have non-proportional effects in cumulative
#'                logit models. These covariates should also be part of the
#'                regular model formula, and the names of the list should be
#'                the names of the ordinal response variables.
#' @param ... additional, optional arguments
#'            \describe{
#'            \item{`trunc`}{named list specifying limits of truncation for the
#'                 distribution of the named incomplete variables (see vignette)}
#'            \item{`hyperpars`}{list of hyper-parameters, as obtained by
#'                 \code{\link{default_hyperpars}()}}
#'            \item{`scale_vars`}{named vector of (continuous) variables that
#'                 will be centred and scaled (such that mean = 0 and sd = 1)
#'                 when they enter a linear predictor to improve
#'                 convergence of the MCMC sampling. Default is that all
#'                 numeric variables and integer variables with >20 different
#'                 values will be scaled.
#'                 If set to \code{FALSE} no scaling will be done.}
#'            \item{`custom`}{named list of JAGS model chunks (character strings)
#'                 that replace the model for the given variable.}
#'            \item{`append_data_list`}{list that will be appended to the list
#'                 containing the data that is passed to **rjags**
#'                 (`data_list`). This may be necessary if additional data /
#'                 variables are needed for custom (covariate) models.}
#'            \item{`progress.bar`}{character string specifying the type of
#'                 progress bar. Possible values are "text" (default), "gui",
#'                 and "none" (see \code{\link[rjags]{update}}). Note: when
#'                 sampling is performed in parallel it is not possible to
#'                 display a progress bar.}
#'            \item{`quiet`}{logical; if \code{TRUE} then messages generated by
#'                 \strong{rjags} during compilation as well as the progress bar
#'                 for the adaptive phase will be suppressed,
#'                 (see \code{\link[rjags]{jags.model}})}
#'            \item{`keep_scaled_mcmc`}{should the "original" MCMC sample (i.e.,
#'                 the scaled version returned by \code{coda.samples()}) be
#'                 kept? (The MCMC sample that is re-scaled to the scale of the
#'                 data is always kept.)}
#'            \item{`modelname`}{character string specifying the name of the
#'                  model file (including the ending, either .R or .txt). If
#'                  unspecified a random name will be generated.}
#'            \item{`modeldir`}{directory containing the model file or directory
#'                 in which the model file should be written. If unspecified a
#'                 temporary directory will be created.}
#'            \item{`overwrite`}{logical; whether an existing model file with
#'                 the specified \code{<modeldir>/<modelname>} should be
#'                 overwritten. If set to \code{FALSE} and a model already
#'                 exists, that model will be used. If unspecified (\code{NULL})
#'                 and a file exists, the user is asked for input on how to
#'                 proceed.}
#'            \item{`keep_model`}{logical; whether the created JAGS model file
#'                 should be saved or removed from (\code{FALSE}; default) when
#'                 the sampling has finished.}
#' }
#' @name model_imp_custom
#'
#' @details # The modification is to allow user-specified order on the sequence
#' of imputation models.
#'
#' default Sequence of models:
#' Models generated automatically (i.e., not mentioned in `formula` or `fixed`
#' are specified in a sequence based on the level of the outcome of the
#' respective model in the multi-level hierarchy and within each level
#' according to the number of missing values.
#'
#' After the modification, users can specify the order through \code{model_order}.
#'
#' the following scenario gives the sequential imputation models in an order of
#' y0, y1, y3, y2, y4,..., which is based on 'nmis' variable. The expected order
#' would be y0, y1, y2, y3, y4,.... To specify the order, what we need is just set
#' model_order = c('y0', 'y1', 'y2', 'y3', 'y4', 'y5')
#'
#'
#    rowID   out lvl nmis nlev ordered type L1
# 8 rowID1 FALSE   1  428    3    TRUE  clm y5
# 7 rowID1 FALSE   1  426    4    TRUE  clm y4
# 5 rowID1 FALSE   1  423    4    TRUE  clm y2
# 6 rowID1 FALSE   1   63    4    TRUE  clm y3
# 4 rowID1 FALSE   1   11    4    TRUE  clm y1
# 3 rowID1 FALSE   1    3    4    TRUE  clm y0
#'
#' @noRd


model_imp_custom <- function(formula = NULL, fixed = NULL, data, random = NULL,
                      trtvar = NULL, family = NULL, df_basehaz = NULL,
                      rd_vcov = "blockdiag",
                      n.chains = 3, n.adapt = 100, n.iter = 0, thin = 1,
                      monitor_params = c(analysis_main = TRUE), auxvars = NULL,
                      timevar = NULL, refcats = NULL,
                      models = NULL, no_model = NULL, trunc = NULL,
                      model_order = NULL, ord_cov_dummy = TRUE,
                      shrinkage = FALSE, custom = NULL,
                      nonprop = NULL, rev = NULL,
                      ppc = TRUE, seed = NULL, inits = NULL,
                      scale_vars = NULL, hyperpars = NULL,
                      modelname = NULL, modeldir = NULL,
                      keep_model = FALSE, overwrite = NULL,
                      quiet = TRUE, progress.bar = "text",
                      warn = TRUE, mess = TRUE,
                      keep_scaled_mcmc = FALSE,
                      analysis_type, assoc_type = NULL,
                      append_data_list = NULL, ...) {

  modimpcall <- as.list(match.call())[-1L]


  # checks & warnings -------------------------------------------------------
  if (!is.null(formula) & is.null(fixed) & is.null(random)) {
    formula <- check_formula_list(formula)
    fixed <- split_formula_list(formula)$fixed
    random <- split_formula_list(formula)$random
  }

  # check if the arguments meth, n.cores or parallel are provided
  # (no longer used)
  args <- as.list(match.call())
  if (!is.null(args$meth))
    errormsg("The argument %s has been replaced by the argument %s.",
              dQuote("meth"), dQuote("models"))

  if (!is.null(args$parallel) | !is.null(args$n.cores)) {
    errormsg("The arguments %s and %s are no longer used. To perform the
             computation in parallel, specify future::plan().
             For an example, see ?model_imp.",
             dQuote("parallel"), dQuote("n.cores"))
  }


  # data pre-processing --------------------------------------------------------
  data <- check_data(data, fixed, random, auxvars, timevar, mess)


  # * divide matrices ----------------------------------------------------------
  Mlist <- divide_matrices_custom(data, fixed, analysis_type = analysis_type,
                           random = random, models = models, auxvars = auxvars,
                           timevar = timevar, no_model = no_model, trtvar = trtvar,
                           model_order = model_order, ord_cov_dummy = ord_cov_dummy,
                           scale_vars = scale_vars, refcats = refcats,
                           nonprop = nonprop, rev = rev,
                           warn = warn, mess = mess, ppc = ppc,
                           shrinkage = shrinkage, df_basehaz = df_basehaz,
                           rd_vcov = rd_vcov)

  # * model dimensions ---------------------------------------------------------
  par_index_main <- get_model_dim(Mlist$lp_cols[names(Mlist$lp_cols) %in%
                                     names(Mlist$fixed)],
                     Mlist = Mlist)
  par_index_other <- get_model_dim(Mlist$lp_cols[!names(Mlist$lp_cols) %in%
                                         names(Mlist$fixed)],
                         Mlist = Mlist)

  # * model info ---------------------------------------------------------------
  info_list <- get_model_info(Mlist, par_index_main = par_index_main,
                              par_index_other = par_index_other,
                              trunc = trunc, assoc_type = assoc_type, custom =custom )

  if (!ord_cov_dummy){
    for (k in names(info_list)){
      if (info_list[[k]][["modeltype"]] == "clm" & info_list[[k]][["ncat"]] ==0)
        info_list[[k]][["ncat"]] = length(levels(Mlist$data[,k]))
    }
  }

  # * data list ----------------------------------------------------------------
  data_list <- get_data_list(Mlist, info_list, hyperpars, append_data_list)

  # write model ----------------------------------------------------------------
  modelfile <- make_filename(modeldir = modeldir, modelname = modelname,
                             keep_model = keep_model, overwrite = overwrite,
                             mess = mess)

  if (!file.exists(modelfile) || (file.exists(modelfile) &
                                  attr(modelfile, "overwrite") == TRUE)) {
      write_model(info_list = info_list, Mlist = Mlist, modelfile = modelfile)
  }

  # initial values -------------------------------------------------------------
  inits <- get_initial_values(inits = inits, seed = seed, n_chains = n.chains,
                              warn = warn)

  # parameters to monitor ------------------------------------------------------
  if (any(grepl("^beta\\[", unlist(monitor_params))) &
      any(!is.na(unlist(Mlist$scale_pars)))) {

    monitor_params <- c(lapply(monitor_params, function(x) {
      if (any(grepl("^beta[", x, fixed = TRUE))) {
        x[-grep("^beta[", x, fixed = TRUE)]
      } else {
        x
      }
    }),
    betas = TRUE)

    if (mess)
      msg("Note: %s was set in %s because re-scaling of the effects of
             the regression coefficients in the main model(s) requires all
             of them to be monitored.",
          dQuote("betas = TRUE"), dQuote("monitor_params"))
  }

  var_names <- do.call(get_params, c(list(Mlist = Mlist, info_list = info_list,
                                          mess = mess),
                                     monitor_params))

  # run JAGS -----------------------------------------------------------------
  # Message if no MCMC sample will be produced.
  if (n.iter == 0) {
    if (mess)
      msg("Note: No MCMC sample will be created when n.iter is set to 0.")
  }

  future_info <- get_future_info()

  if (n.iter == 0) run_jags = run_seq
  else run_jags <- ifelse(future_info$parallel, run_parallel, run_seq)

  t0 <- Sys.time()
  jags_res <- run_jags(n_adapt = n.adapt, n_iter = n.iter, n_chains = n.chains,
                       inits = inits, thin = thin,
                       n_workers = future_info$workers,
                       data_list = data_list, var_names = var_names,
                       modelfile = modelfile, quiet = quiet,
                       progress_bar = progress.bar, mess = mess, warn = warn)
  adapt <- jags_res$adapt
  mcmc <- jags_res$mcmc

  t1 <- Sys.time()

  if (n.iter > 0 & !inherits(mcmc, "mcmc.list"))
    warnmsg("There is no mcmc sample. Something went wrong.")

  # post processing ------------------------------------------------------------
  if (n.iter > 0 & !is.null(mcmc)) {
    MCMC <- mcmc

    if (any(!vapply(Mlist$scale_pars, is.null, FUN.VALUE = logical(1)),
            !is.na(unlist(Mlist$scale_pars)))) {
      coefs <- try(get_coef_names(info_list))

      for (k in seq_len(length(MCMC))) {
        MCMC[[k]] <- coda::as.mcmc(
          rescale(MCMC[[k]],
                  coefs = do.call(rbind, coefs),
                  scale_pars = do.call(rbind, unname(Mlist$scale_pars)),
                  info_list = info_list,
                  data_list = data_list,
                  groups = Mlist$groups))
        attr(MCMC[[k]], "mcpar") <- attr(mcmc[[k]], "mcpar")
      }
    }
  }


  # prepare output -------------------------------------------------------------
  mcmc_settings <- list(modelfile = modelfile,
                        n.chains = n.chains,
                        n.adapt = n.adapt,
                        n.iter = n.iter,
                        variable.names = if (exists("var_names")) var_names,
                        thin = thin,
                        inits = inits,
                        seed = seed)


  object <- structure(
    list(analysis_type = analysis_type,
         data = Mlist$data,
         models = Mlist$models,
         fixed = Mlist$fixed,
         random = Mlist$random,
         Mlist = Mlist[setdiff(names(Mlist), c("data", "models",
                                               "fixed", "random",
                                               "M"))],
         par_index_main = par_index_main,
         par_index_other = par_index_other,
         jagsmodel = structure(readChar(modelfile,
                                        file.info(modelfile)$size),
                               class = "modelstring"),
         mcmc_settings = mcmc_settings,
         monitor_params = c(monitor_params,
                            if (!"analysis_main" %in% names(monitor_params))
                              setNames(TRUE, "analysis_main")),
         data_list = data_list,
         hyperpars = if (is.null(hyperpars)) default_hyperpars() else hyperpars,
         info_list = info_list,
         coef_list = get_coef_names(info_list),
         model = if (n.adapt > 0) adapt,
         sample = if (n.iter > 0 & !is.null(mcmc) & keep_scaled_mcmc) mcmc,
         MCMC = if (n.iter > 0 & !is.null(mcmc)) coda::as.mcmc.list(MCMC),
         comp_info = list(start_time = t0,
                          duration = t1 - t0,
                          remiod_version = packageVersion("remiod"),
                          future = future_info$call),
         call = modimpcall$thecall
    ), class = "remiod")


  object$fitted.values <- try(fitted_values(object, mess = FALSE, warn = FALSE),
                              silent = TRUE)

  object$residuals <- try(residuals(object, type = "working", warn = FALSE),
                          silent = TRUE)

  if (inherits(object$fitted.values, "try-error"))
    object$fitted.values <- NULL
  if (inherits(object$residuals, "try-error"))
    object$residuals <- NULL

  if (inherits(adapt, "try-error"))
    class(object) <- "remiod_errored"

  if (!attr(modelfile, "keep_model")) {
    file.remove(modelfile)
  }

  return(object)
}


#'
#'
glm_imp_custom <- function(formula, family, data, trtvar = NULL,
                    n.chains = 3, n.adapt = 100, n.iter = 0, thin = 1,
                    monitor_params = c(analysis_main = TRUE), auxvars = NULL,
                    refcats = NULL,
                    models = NULL, no_model = NULL, model_order = NULL,
                    shrinkage = FALSE, ppc = TRUE, seed = NULL, inits = NULL,
                    warn = TRUE, mess = TRUE, ord_cov_dummy = TRUE,
                    ...) {

  if (missing(formula)) errormsg("No model formula specified.")
  if (missing(family))
    errormsg("The argument %s needs to be specified.", dQuote("family"))

  arglist <- prep_arglist(analysis_type = "glm",
                          family = family,
                          formals = formals(), call = match.call(),
                          sframe = sys.frame(sys.nframe()))

  do.call(model_imp_custom, arglist)
}


#'
#'
clm_imp_custom <- function(formula, data, trtvar = NULL,
                    n.chains = 3, n.adapt = 100, n.iter = 0, thin = 1,
                    monitor_params = c(analysis_main = TRUE), auxvars = NULL,
                    refcats = NULL, nonprop = NULL, rev = NULL,
                    models = NULL, no_model = NULL, model_order = NULL,
                    shrinkage = FALSE, ppc = TRUE, seed = NULL, inits = NULL,
                    warn = TRUE, mess = TRUE, ord_cov_dummy = TRUE, ...) {

  if (missing(formula)) errormsg("No model formula specified.")

  arglist <- prep_arglist(analysis_type = "clm",
                          formals = formals(), call = match.call(),
                          sframe = sys.frame(sys.nframe()))

  do.call(model_imp_custom, arglist)
}


opm_imp_custom <- function(formula, data, trtvar = NULL,
                           n.chains = 3, n.adapt = 100, n.iter = 0, thin = 1,
                           monitor_params = c(analysis_main = TRUE), auxvars = NULL,
                           refcats = NULL, nonprop = NULL, rev = NULL,
                           models = NULL, no_model = NULL, model_order = NULL,
                           shrinkage = FALSE, ppc = TRUE, seed = NULL, inits = NULL,
                           warn = TRUE, mess = TRUE, ord_cov_dummy = TRUE, ...) {

  if (missing(formula)) errormsg("No model formula specified.")

  arglist <- prep_arglist(analysis_type = "opm",
                          formals = formals(), call = match.call(),
                          sframe = sys.frame(sys.nframe()))

  do.call(model_imp_custom, arglist)
}



