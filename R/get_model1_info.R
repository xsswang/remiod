
# get model info for a single model
get_model1_info <- function(k, Mlist, par_index_main, par_index_other,
                            trunc = NULL, assoc_type = NULL,
                            custom = NULL, isgk = FALSE) {

  arglist <- as.list(match.call())[-1L]

  # model type, family & link -------------------------------------------------
  modeltype <- get_modeltype(Mlist$models[k])
  family <- get_family(Mlist$models[k])
  link <- get_link(Mlist$models[k])


  if (modeltype %in% c("glm", "clm", "opm", "mlogit") & !is.null(Mlist$random[[k]])) {
    errormsg("There is a random effects structure specified for the variable %s,
             but %s is being modeled with a model of type %s.",
             dQuote(k), dQuote(k), dQuote(modeltype))
  }

  # response matrix and column(s) --------------------------------------------
  resp_mat <- get_resp_mat(
    resp = k, Mlvls = Mlist$Mlvls,
    outnames = if (!is.null(Mlist$outcomes$outcomes[[k]])) {
      names(Mlist$outcomes$outcomes[[k]])
    } else {
      k
    })
  # resp_mat <- if (k %in% names(Mlist$Mlvls)) {
  #   # if the variable is a column of one of the design matrices, use the level
  #   # of that matrix
  #   Mlist$Mlvls[k]
  # } else if (attr(Mlist$fixed[[k]], "type") %in% c("survreg", "coxph", "JM")) {
  #   # if the model is a survival model (variable name is the survival expression
  #   # and not a single variable name) get the levels of the separate variables
  #   # involved in the survival expression
  #   if (all(names(Mlist$outcomes$outcomes[[k]]) %in% names(Mlist$Mlvls))) {
  #     Mlist$Mlvls[names(Mlist$outcomes$outcomes[[k]])]
  #   } else {
  #     errormsg("I have identified %s as a survival outcome, but I cannot find
  #              some of its elements in any of the matrices %s.",
  #              dQuote(k), dQuote("M"))
  #   }
  # } else {
  #   errormsg("I cannot find the variable %s in any of the matrices %s.",
  #            dQuote(k), dQuote("M"))
  # }

  resp_col <- if (k %in% names(Mlist$fixed) &&
                  attr(Mlist$fixed[[k]], "type") %in%
                  c("survreg", "coxph", "JM")) {
    vapply(names(Mlist$outcomes$outcomes[[k]]), function(x) {
      match(x, colnames(Mlist$M[[resp_mat[x]]]))
    }, FUN.VALUE = numeric(1L))
  } else {
    match(k, colnames(Mlist$M[[resp_mat[1L]]]))
  }

  # linear predictor columns -------------------------------------------------
  lp <- get_lp(k, Mlist)

  # check and error message for time-dependent cox with interaction between
  # time-varying covariate with incomplete baseline covariate
  if (modeltype == "coxph" &
      !is.null(Mlist$interactions) &
      !is.null(lp$M_lvlone)) {

    lapply(Mlist$interactions[intersect(names(Mlist$interactions),
                                        names(lp$M_lvlone))],
           function(x) {
             if (any(names(x$elmts) == "M_lvlone") &
                 any(names(x$elmts) != "M_lvlone") &
                 attr(x, "has_NAs")
             ) {
               errormsg("It seems that there is an interaction between a time-varying
                 covariate and a incomplete baseline covariate. This is not
                 yet implemented for proportional hazards models.")
             }
           })
  }


  # parameter elements ------------------------------------------------------
  parelmts <- get_parelmts(k, Mlist, par_index_main, par_index_other, lp)


  # scaling parameter matrices -----------------------------------------------
  scale_pars <- Mlist$scale_pars

  # dummy columns -------------------------------------------------------------
  dummy_cols <- if (k %in% names(Mlist$refs) &
                    (any(is.na(Mlist$M[[resp_mat[1L]]][, resp_col[1L]])) |
                     any(vapply(Mlist$fixed, "attr", "type",
                                FUN.VALUE = character(1L)
                     ) %in% "JM"))) {
    match(attr(Mlist$refs[[k]], "dummies"), colnames(Mlist$M[[resp_mat[[1L]]]]))
  }

  if (all(is.na(dummy_cols)))
    dummy_cols <- NULL


  # index name -----------------------------------------------------------------
  index <- get_indices(Mlist)


  # transformations ------------------------------------------------------------
  trafos <- paste_trafos(Mlist, varname = k,
                         index = index[gsub("M_", "", resp_mat[1L])],
                         isgk = isgk)

  # JM settings ----------------------------------------------------------------
  # * covariate names ----------------------------------------------------------
  covnames <- if (modeltype %in% "JM") {
    unique(unlist(
      lapply(lp, function(x) {
        vapply(names(x), replace_dummy, refs = Mlist$refs,
               FUN.VALUE = character(1L))
      })
    ))
  }

  # * time-varying covariates --------------------------------------------------
  tv_vars <- if (modeltype %in% "JM") {

    # find the (longitudinal) covariates involved in the lp of the survival part
    covars <- unlist(
      lapply(unlist(lapply(lp, names)),
             replace_trafo, Mlist$fcts_all)
    )
    covars <- vapply(covars, replace_dummy, refs = Mlist$refs,
                     FUN.VALUE = character(1L))
    covars <- covars[covars %in% unlist(lapply(Mlist$M, colnames))]


    rep_lvls <- names(which(Mlist$group_lvls < Mlist$group_lvls[
      gsub("M_", "", resp_mat[2L])]))

    tvars <- unique(unlist(c(lapply(lp[paste0("M_", rep_lvls)], names),
                             lapply(Mlist$lp_cols[covars], function(x)
                               names(unlist(unname(x[paste0("M_", rep_lvls)]))))
    )))


    # get the variables needed to re-fit the models for "covars" in the
    # Gauss-Kronrod quadrature
    tvars <- unique(unlist(lapply(tvars, replace_interaction, Mlist$interactions)))
    tvars <- unlist(lapply(tvars, replace_trafo, Mlist$fcts_all))

    tvars <- unique(vapply(tvars[!tvars %in% Mlist$timevar],
                           replace_dummy, refs = Mlist$refs,
                           FUN.VALUE = character(1L))
    )

    tvars <- tvars[Mlist$Mlvls[tvars] %in% paste0("M_", rep_lvls)]

    # get the model info for these variables
    setNames(lapply(tvars, function(i) {
      arglist_new <- arglist
      arglist_new$k <- replace_dummy(i, refs = Mlist$refs)
      arglist_new$isgk <- TRUE
      subinfo <- do.call(get_model1_info, arglist_new)
      subinfo$surv_lvl <- gsub("M_", "", resp_mat[2L])
      subinfo
    }), tvars)
  }



  # Hierarchical centering -----------------------------------------------------
  hc_list <- get_hc_info(varname = k,
                         resplvl = gsub("M_", "", resp_mat[length(resp_mat)]),
                         Mlist, parelmts, lp)

  nranef <- vapply(hc_list$hcvars, function(x) {
    as.integer(attr(x, "rd_intercept")) +
      if (!is.null(x$rd_slope_coefs)) {
        nrow(x$rd_slope_coefs)
        # if (any(!vapply(x$rd_slope_coefs, is.null, FUN.VALUE = logical(1L)))) {
        #   nrow(do.call(rbind, x$rd_slope_coefs))
      } else {
        0L
      }
  }, FUN.VALUE = integer(1L))


  # random effects variance covariance matrix
  rd_vcov <- nlapply(names(Mlist$rd_vcov), function(lvl) {
    x <- Mlist$rd_vcov[[lvl]]

    w <- which(lvapply(x, function(z) k %in% z))
    if (length(w) > 0L) {
      type <- names(w)

      attr(type, "ranef_index") <- attr(x[[w]], "ranef_index")
      attr(type, "name") <- attr(x[[w]], "name")
      type
    } else if (isTRUE(nranef[lvl] > 0L)) {
      "blockdiag"
    }
  })

  # shrinkage ------------------------------------------------------------------
  shrinkage <- if (k %in% names(Mlist$shrinkage)) {
    Mlist$shrinkage[k]
  } else if (is.null(names(Mlist$shrinkage))) {
    Mlist$shrinkage
  }


  assoc_type <- if (modeltype %in% "JM") {
    covrs <- unique(unlist(lapply(names(unlist(unname(lp))),
                                  replace_dummy, Mlist$refs)))
    get_assoc_type(intersect(tvars, covrs),
                   Mlist$models, assoc_type, Mlist$refs)
  } else if (modeltype %in% "coxph") {
    "obs.value"
  } else if (isTRUE(isgk)) {
    get_assoc_type(k, Mlist$models, assoc_type, Mlist$refs)
  }

  # collect all info ---------------------------------------------------------
  list(
    varname = if (modeltype %in% c("survreg", "coxph", "JM")) {
      clean_survname(k)
    } else {
      k
    },
    modeltype = modeltype,
    family = family,
    link = link,
    resp_mat = resp_mat,
    resp_col = resp_col,
    dummy_cols = dummy_cols,
    ncat = length(levels(Mlist$refs[[k]])),
    rev = k %in% Mlist$rev,
    lp = lp,
    parelmts = parelmts,
    scale_pars = scale_pars,
    index = index,
    parname = switch(as.character(k %in% names(Mlist$fixed)),
                     "TRUE" = "beta",
                     "FALSE" = "alpha"),
    hc_list = if (length(hc_list) > 0L) hc_list,
    nranef = nranef,
    rd_vcov = rd_vcov,
    group_lvls = Mlist$group_lvls,
    trafos = trafos,
    trunc = trunc[[k]],
    custom = custom[[k]],
    ppc = FALSE,
    shrinkage = shrinkage,
    refs = Mlist$refs[[k]],
    covnames = covnames,
    assoc_type  = assoc_type,
    tv_vars = tv_vars,
    N = Mlist$N,
    df_basehaz = Mlist$df_basehaz
  )
}



get_model_info = function (Mlist, par_index_main, par_index_other, trunc = NULL,
                           assoc_type = NULL, custom = NULL)
{
  args <- as.list(match.call())[-1L]
  setNames(lapply(names(Mlist$lp_cols), function(k) {
    do.call(get_model1_info, c(k = replace_dummy(k, refs = Mlist$refs),
                               args))
  }), names(Mlist$lp_cols))
}
