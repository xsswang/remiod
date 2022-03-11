#' Extract MCMC samples of monitored parameters from JointAI object. Code from JointAI.
#'
#' @param object an object of class JointAI
#' @param treatment the variable name of treatment. Reference level of treatment should be coded as 0.
#' @param delta specific value used for Delta adjustment, applicable only
#'              for method="delta".
#' @param start first iteration to be used
#' @param end last iteration to be used
#' @param thin thinning to be applied
#' @param subset subset of parameters (columns of the mcmc object) to be used
#' @param exclude_chains optional vector of numbers, indexing MCMC chains to be
#                   excluded from the output
#' @param warn logical, should warning messages be displayed?
#' @param mess logical, should messages be displayed?
#' @param ... optional arguments passed from \code{remiod}.
#'
#' @importFrom utils getFromNamespace globalVariables
#' @importFrom JointAI get_MIdat
#' @importFrom mcmcse mcse
#' @keywords internal
#'
prep_MCMC <- function(object, treatment=NULL, delta=0, start = NULL, end = NULL, thin = NULL,
                      subset = FALSE, exclude_chains = NULL, warn = TRUE, mess = TRUE, ...) {

  # set start, end and thin (or use values from MCMC sample)
  if (is.null(start)) {
    start <- start(object$MCMC)
  } else {
    start <- max(start, start(object$MCMC))
  }

  if (is.null(end)) {
    end <- end(object$MCMC)
  } else {
    end <- min(end, end(object$MCMC))
  }

  if (is.null(thin))
    thin <- coda::thin(object$MCMC)

  # obtain subset of parameters of the MCMC samples
  MCMC <- get_subset(object, subset=subset, warn = warn, mess = mess)

  # exclude chains, if set by user
  chains <- seq_along(MCMC)
  if (!is.null(exclude_chains)) {
    chains <- chains[-exclude_chains]
  }

  # restrict MCMC sample to selected iterations and chains
  MCMC <- do.call(rbind,
                  window(MCMC[chains],
                         start = start,
                         end = end,
                         thin = thin))
  return(MCMC)
}


## this function updates the elements of list 'object' to contain all of the elements
## of 'new', overwriting elements with the same name, and (optionally) copying unnamed
## elements.
list_update <- function(object, new, unnamed=FALSE)
{
  retval <- object

  for(name in names(new))
    retval[[name]] <- new[[name]]

  if(unnamed)
  {
    if(is.null(names(new)))
      names(new) <- rep("", length=length(new))
    for(i in (1:length(new))[names(new)==""] )
      retval <- append(retval, new[[i]])
  }

  retval
}



##
minmax_mat <- function(mat, minval = 1e-10, maxval = 1 - 1e-10) {
  apply(mat, 2, function(x) {
    pmax(minval, pmin(maxval, x))
  })
}

## transform mu to binary outcome

get_bin = function(mu, linkinvf=NULL, seed=NULL){
  pred = linkinvf(mu)
  if(!is.matrix(pred)) pmu = matrix(pred, ncol=1)
  else pmu = pred

  probs = cbind(pmu, 1-pmu)
  if (!is.null(seed)) set.seed(seed)
  class = numeric()
  for (i in 1:nrow(probs)) class[i] = rbinom(1,1,prob=probs[i,])
  predc = matrix(class,ncol=1)
  return(predc)
}

## sequential or parallel run MI for multiple delta

delta_seq = function(object, dtimp, treatment, algorithm, method="delta", delta = 0,
                     exclude_chains=NULL, start=NULL, end=NULL, thin=NULL,
                     subset=FALSE, seed=NULL, include=TRUE, mess=TRUE, ...)
  {
  for(i in 1:length(delta)) {
    if (algorithm == "tang_seq"){
      dti = tang_MI_RB(object=object, dtimp=dtimp, treatment=treatment, method=method,
                       delta=delta[i], exclude_chains=exclude_chains, include=include)
    } else {
      dti = get_MI_RB(object=object, treatment=treatment, method=method, delta=delta[i],
                      exclude_chains=exclude_chains, thin=thin, include=include,
                      start=start, end=end, seed=seed,...)
    }
    dimp = data.frame(delta=delta[i], dti)
  }
  if (inherits(dimp, "list")) dimp = as.data.frame(data.table::rbindlist(dimp))
  return(dimp)
}

delta_parallel = function(object, dtimp, treatment, algorithm, method="delta", delta = 0,
                          n_workers, exclude_chains=NULL, start=NULL, end=NULL,
                          thin=NULL, subset=FALSE, seed=NULL, include=TRUE, mess=TRUE, ...)
{

    if (mess)
      msg("Parallel computing delta adjustment with %s workers started (%s).",
          eval(n_workers), Sys.time())

  if (algorithm == "tang_seq"){
    res <- foreach::`%dopar%`(foreach::foreach(i = seq_along(delta)),
                              tang_MI_RB(object=object, dtimp=dtimp, treatment=treatment,
                                         method=method, delta=delta[i], include=include,
                                         exclude_chains=exclude_chains)
                              )
  } else {
    res <- foreach::`%dopar%`(foreach::foreach(i = seq_along(delta)),
                              get_MI_RB(object=object, treatment=treatment, method=method,
                                        delta=delta[i], exclude_chains=exclude_chains,
                                        start=start, end=end, thin=thin, include=include,
                                        seed=seed,...)
                              )
    }

   dimp = data.table::rbindlist(res)
   return(dimp)
}

##
GR_crit = function (object, confidence = 0.95, transform = FALSE, autoburnin = TRUE,
          multivariate = TRUE, subset = NULL, exclude_chains = NULL,
          start = NULL, end = NULL, thin = NULL, warn = TRUE, mess = TRUE,
          ...)
{
  # if (!inherits(object, "remiod"))
  #   errormsg("Object must be of class \"remiod\".")
  if (is.null(object$MCMC))
    errormsg("No MCMC sample.")
  if (is.null(start))
    start <- start(object$MCMC)
  if (is.null(end))
    end <- end(object$MCMC)
  if (is.null(thin))
    thin <- coda::thin(object$MCMC)
  MCMC <- get_subset(object, subset, warn = warn, mess = mess)
  chains <- seq_along(MCMC)
  if (!is.null(exclude_chains)) {
    chains <- chains[-exclude_chains]
  }
  MCMC <- window(MCMC[chains], start = start, end = end, thin = thin)
  plotnams <- get_plotmain(object, colnames(MCMC[[1]]), ylab = TRUE)
  #for (i in seq_len(length(MCMC))) colnames(MCMC[[i]]) <- plotnams
  coda::gelman.diag(x = MCMC, confidence = confidence, transform = transform,
                    autoburnin = autoburnin, multivariate = multivariate)
}


MC_error = function (MCMC, digits = 2, warn = TRUE, mess = TRUE,...)
{
  # plotnams <- get_plotmain(x, colnames(MCMCsub), ylab = TRUE)
  # colnames(MCMC) <- plotnams
  MCE1 <- t(apply(MCMC, 2, function(k) {
    mce <- try(mcmcse::mcse(k, ...), silent = TRUE)
    if (inherits(mce, "try-error")) {
      c(NA, NA)
    }
    else {
      unlist(mce)
    }
  }))
  colnames(MCE1) <- c("est", "MCSE")
  MCE1 <- cbind(MCE1, SD = apply(MCMC, 2, sd)[match(colnames(MCMC),
                                                    row.names(MCE1))])
  MCE1 <- cbind(MCE1, `MCSE/SD` = MCE1[, "MCSE"]/MCE1[, "SD"])

  out <- list(data_scale = MCE1, digits = digits)
  class(out) <- "MCElist"
  return(out)
}

##
get_subset = function (object, subset, warn = TRUE, mess = TRUE) {
  if (identical(subset, FALSE))
    return(object$MCMC)

  if (!is.list(subset)) subset <- as.list(subset)

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


##
check_data = function (data, fixed, random, auxvars, timevar, mess) {
  check_vars_in_data(names(data), fixed = fixed, random = random,
                     auxvars = auxvars, timevar = timevar)
  check_classes(data, fixed = fixed, random = random, auxvars = auxvars)
  data <- drop_levels(data = data, allvars = all_vars(c(fixed, random, auxvars)), mess = mess)

  # data <- convert_variables(data = data, allvars = all_vars(c(fixed, random, auxvars)), mess = mess)
  data
}

##
get_scale_pars = function (mat, groups, scale_vars, refs, fcts_all, interactions, data)
{
  if (is.null(mat) | (is.null(scale_vars)))
    return(NULL)
  vars <- find_scalevars(mat, refs, fcts_all, interactions, data)
  if (!is.null(scale_vars))
    vars <- intersect(vars, scale_vars)
  rows <- match(unique(groups), groups)
  do.call(rbind, sapply(colnames(mat), function(k) {
    if (k %in% vars) {
      scaled <- scale(mat[rows, k])
      data.frame(center = attr(scaled, "scaled:center"),
                 scale = attr(scaled, "scaled:scale"))
    }
    else {
      data.frame(center = NA, scale = NA)
    }
  }, simplify = FALSE))
}


# message functions ------------------------------------------------------------
errormsg <- function(x, ...) {
  stop(strwrap(gettextf(x, ...), prefix = "\n"), call. = FALSE)
}

msg <- function(x, ..., exdent = 0L) {
  message(strwrap(gettextf(x, ...), prefix = "\n", exdent = exdent))
}

warnmsg <- function(x, ..., exdent = 0L) {
  warning(strwrap(gettextf(x, ...), prefix = "\n", exdent = exdent),
          call. = FALSE, immediate. = TRUE)
}

paste_and <- function(x) {
  x1 <- paste0(x[-length(x)], collapse = ", ")

  if (length(x) > 1L) {
    paste(x1, x[length(x)], sep = " and ")
  } else {
    x
  }
}

nlapply <- function(x, fun, ...) {
  # a named version of lapply, intended to replace sapply(..., simplify = FALSE)

  l <- lapply(x, fun, ...)
  if (is.null(names(l)))
    if (!is.null(names(x))) {
      names(l) <- names(x)
    } else if (is.character(x)) {
      names(l) <- x
    }
  l
}

lvapply <- function(x, fun, ...) {
  vapply(x, fun, FUN.VALUE = logical(1L), ..., USE.NAMES = TRUE)
}

ivapply <- function(x, fun, ...) {
  vapply(x, fun, FUN.VALUE = integer(1L), ..., USE.NAMES = TRUE)
}

nvapply <- function(x, fun, ...) {
  vapply(x, fun, FUN.VALUE = numeric(1L), ..., USE.NAMES = TRUE)
}

cvapply <- function(x, fun, ...) {
  vapply(x, fun, FUN.VALUE = character(1L), ..., USE.NAMES = TRUE)
}

##

print_seq <- function(min, max) {

  m <- Map(function(min, max) {
    if (min == max) {
      max
    } else {
      paste0(min, ":", max)
    }
  }, min = min, max = max)

  unlist(m)
}

# used in divide_matrices, get_modeltypes, helpfunctions_checks,
# helpfunctions_formulas, plots, simulate_data
check_varlevel <- function(x, groups, group_lvls = NULL) {
  # identify the level of a variable
  # - x: a vector
  # - groups: a list of grouping information (obtained from get_groups())
  # - group_lvls: the grouping level matrix
  #               (obtained from identify_level_relations())

  # if there are no groups, make a list with group name "no_levels" so that the
  # syntax does not fail for single-level models
  if (!is.list(groups))
    groups <- list("no_levels" = groups)

  # check the clustering of the variable
  clus <- check_cluster(x, grouping = groups)

  # clus is a logical vector, which is TRUE if x varies in a given level and
  # FALSE when x is constant in the level


  if (sum(!clus) > 1L) {
    # if the variable is constant in more than one level, the exact level needs
    # to be determined using the level structure of the grouping
    if (is.null(group_lvls))
      group_lvls <- identify_level_relations(groups)

    names(which.max(colSums(!group_lvls[!clus, !clus, drop = FALSE])))
  } else if (sum(!clus) == 1L) {
    # if the variable is constant in exactly one level, that level is the
    # level of the variable
    names(clus)[!clus]
  } else {
    # if the variable varies in all levels, it is from level one
    "lvlone"
  }
}



# model_info helpers -----------------------------------------------------------
# used in get_model_info
replace_dummy <- function(nam, refs) {
  # check if a variable name is a dummy variable and replace it with the name
  # of the original variable
  # if the variable is a factor
  # - nam: one variable name
  # - refs: list of reference category information (part of Mlist)

  if (is.null(refs)) {
    return(nam)
  }

  dummies <- lapply(refs, "attr",  "dummies")

  if (any(lvapply(dummies, function(k) nam %in% k))) {
    names(dummies)[lvapply(dummies, function(k) nam %in% k)]
  } else {
    nam
  }
}


replace_interaction <- function(nam, interactions) {
  if (nam %in% names(interactions)) {
    attr(interactions[[nam]], "elements")
  } else {
    nam
  }
}

# merge data while keeping order

left_join <- function(x, y, ...) {
  merge_exec(x = x, y = y, all.x = TRUE, ...)
}
right_join <- function(x, y, ...) {
  merge_exec(x = x, y = y, all.y = TRUE, ...)
}
inner_join <- function(x, y, ...) {
  merge_exec(x = x, y = y, all = TRUE, ...)
}
full_join <- function(x, y, ...) {
  merge_exec(x = x, y = y, ...)
}

# workhorse:
merge_exec <- function(x, y, ...) {
  if (!is.data.frame(x)) x = as.data.frame(x)
  if (!is.data.frame(y)) y = as.data.frame(y)
  # set index
  x$join_id_ <- 1:nrow(x)
  # do the join
  joined <- merge(x = x, y = y, sort = FALSE, ...)
  # get suffices (yes, I prefer this over suffixes)
  if ("suffixes" %in% names(list(...))) {
    suffixes <- list(...)$suffixes
  } else {
    suffixes <- c("", "")
  }
  # get columns names in right order, so the 'by' column won't be forced first
  cols <- unique(c(colnames(x),
                   paste0(colnames(x), suffixes[1]),
                   colnames(y),
                   paste0(colnames(y), suffixes[2])))
  # get the original row and column index
  joined[order(joined$join_id),
         cols[cols %in% colnames(joined) & cols != "join_id_"]]
}

# model_imp helpers ------------------------------------------------------------

#' Replace a full with a block-diagonal variance covariance matrix
#' Check if a full random effects variance covariance matrix is specified
#' for a single variable. In that case, it is identical to a block-diagonal
#' matrix. Change the `rd_vcov` specification to `blockdiag` for clarity
#' (because then the variable name is used in the name of `b`, `D`, `invD`, ...)
#'
#' @param rd_vcov a valid random effects variance-covariance structure
#'                specification (i.e., checked using `expand_rd_vcov_full()`)
#' @return a valid random effects variance-covariance structure specification
#' @noRd

check_full_blockdiag <- function(rd_vcov) {

  if (!inherits(rd_vcov, "list") | any(!lvapply(rd_vcov, inherits, "list"))) {
    errormsg("%s should be a list (by grouping level) of lists
    (per covariance matrix).", dQuote("rd_vcov"))
  }

  nlapply(names(rd_vcov), function(lvl) {
    bd <- names(rd_vcov[[lvl]]) == "full" &
      ivapply(rd_vcov[[lvl]], length) == 1
    names(rd_vcov[[lvl]])[bd] <- "blockdiag"
    rd_vcov[[lvl]]
  })
}

#' Check / create the random effects variance-covariance matrix specification
#' @param rd_vcov variance covariance specification provided by the user
#' @param nranef list by level with named vectors of number of random effects
#'               per variable (obtained by `get_nranef()`)
#' @noRd

check_rd_vcov <- function(rd_vcov, nranef) {

  idvar <- names(nranef)

  rd_vcov <- expand_rd_vcov_full(rd_vcov,
                                 rd_outnam = nlapply(nranef, function(r) {
                                   names(r)[r > 0L]}))

  rd_vcov <- check_full_blockdiag(rd_vcov)


  if (any(unlist(lapply(rd_vcov, names)) == "full")) {
    for (lvl in idvar) {

      ## if a full vcov is used, determine the number of random effects
      for (k in which(names(rd_vcov[[lvl]]) == "full")) {

        nrd <- nranef[[lvl]][rd_vcov[[lvl]][[k]]]

        ranef_nr <- print_seq(
          min = cumsum(c(1, nrd))[-(length(nrd) + 1)],
          max = cumsum(nrd)
        )

        attr(rd_vcov[[lvl]][[k]], "ranef_index") <-
          setNames(ranef_nr, rd_vcov[[lvl]][[k]])
      }

      ## if there is more than one full vcov, number them
      if (sum(names(rd_vcov[[lvl]]) %in% "full") > 1) {
        rd_full <- which(names(rd_vcov[[lvl]]) %in% "full")
        for (k in seq_along(rd_full)) {
          attr(rd_vcov[[lvl]][[rd_full[k]]], "name") <- k
        }
      }
    }
  }
  rd_vcov
}


#' First validation for rd_vcov
#'
#' Check if rd_vcov is a list with elements for all grouping levels or does
#' not specify a grouping level. If valid, make sure it is a list per grouping
#' level by duplicating the contents if necessary.
#'
#' @param rd_vcov the random effects variance covariance structure provided by
#'                the user
#' @param idvar vector with the names of the grouping variables
#'              (without "lvlone")
#' @return A named list per grouping level where each elements contains
#'        information on how the random effects variance-covariance matrices on
#'        that level are structured. Per level it can be either a character
#'        string (e.g. `"full"`) or a list specifying structures per (groups) of
#'        variable(s) (e.g. `list(full = c("a", "b"), indep = "c")`)
#' @noRd
check_rd_vcov_list <- function(rd_vcov, idvar) {

  if (!inherits(rd_vcov, "list") | all(!idvar %in% names(rd_vcov))) {
    nlapply(idvar, function(x) rd_vcov)
  } else if (inherits(rd_vcov, "list") & any(!idvar %in% names(rd_vcov))) {
    errormsg("Please provide information on the variance-covariance structure
             of the random effects for all levels.")
  } else {
    rd_vcov
  }
}

#' Extract the number of random effects
#' @param idvar vector of the names of id variables
#' @param random a random effect formula or list of random effects formulas
#' @param data a `data.frame`
#' @return a list by grouping level (`idvar`) with a named vector of the number
#'         of random effects per variable (=names).
#' @noRd
get_nranef <- function(idvar, random, data) {
  nlapply(idvar, function(lvl) {
    if (inherits(random, "formula")) {
      rm_gr <- remove_grouping(random)
      if (lvl %in% names(rm_gr)) {
        ncol(model.matrix(remove_grouping(random)[[lvl]], data = data))
      } else 0L
    } else if (inherits(random, "list")) {
      if (length(random) == 1L) {
        rm_gr <- remove_grouping(random)
        nrd <- if (lvl %in% names(rm_gr)) {
          ncol(model.matrix(remove_grouping(random)[[lvl]], data = data))
        } else 0L
      } else {
        nrd <- ivapply(remove_grouping(random), function(x) {
          if (lvl %in% names(x)) {
            ncol(model.matrix(x[[lvl]], data = data))
          } else {0L}
        })
      }
      names(nrd) <- names(random)
      nrd
    } else {
      errormsg("I expected either a formula or list of formulas.")
    }
  })
}


#' Expand rd_vcov using variable names in case "full" is used
#'
#'
#' @param rd_vcov the random effects variance covariance structure provided by
#'                the user (`check_rd_vcov_list()` is called internally)
#' @param rd_outnam list by grouping level of the names of the outcome variables
#'                  that have random effects on this level
#' @noRd
#' @return A named list per grouping level where each elements contains
#'        information on how the random effects variance-covariance matrices on
#'        that level are structured. Per level there is a list of grouping structures
#'        containing the names of variables in each structure
#'        (e.g. `list(full = c("a", "b"), indep = "c")`)

expand_rd_vcov_full <- function(rd_vcov, rd_outnam) {
  idvar <- names(rd_outnam)

  rd_vcov <- check_rd_vcov_list(rd_vcov, idvar)

  nlapply(idvar, function(lvl) {
    if (is.character(rd_vcov[[lvl]]) & length(rd_vcov[[lvl]]) == 1) {

      setNames(list(rd_outnam[[lvl]]), rd_vcov[[lvl]])

    } else if (inherits(rd_vcov[[lvl]], "list")) {

      if (setequal(unlist(rd_vcov[[lvl]]), rd_outnam[[lvl]])) {
        rd_vcov[[lvl]]
      } else {
        errormsg("According to the random effects formulas, there are
                 random effects on the level %s for the models for %s but in
                 the structure specified for the random effects
                 variance-covariance matrices the variables %s have random
                 effects on this level.",
                 dQuote(lvl), paste_and(dQuote(rd_outnam[[lvl]])),
                 paste_and(dQuote(unlist(rd_vcov[[lvl]])))
        )
      }

    } else {
      errormsg("%s should be a character string or a list.",
               dQuote("rd_vcov[[lvl]]"))
    }
  })
}

# used in add_samples, get_params, model_imp
get_coef_names <- function(info_list) {
  # extract the names of the regression coefficients and the corresponding
  # variable names
  # - info_list: a model info list (obtained from get_model_info())

  nlapply(info_list, function(info) {

    # find all parameter elements with the same parameter name to find
    # out if this parameter needs to get indexed or not
    pars <- nlapply(info_list, function(k) {
      if (k$parname %in% info$parname)
        unlist(c(k$parelmts, lapply(k$parelmts, "attr", "nonprop")))
    })


    parelmts <- unlist(unname(info$parelmts), recursive = FALSE)

    if (!is.list(parelmts)) {
      parelmts <- list(parelmts)
      names(parelmts) <- NA
    }


    out <- if (any(ivapply(info$lp, length) > 0L)) {
      data.frame(outcome = unname(info$varname),
                 outcat = rep(names(parelmts), ivapply(parelmts, length)),
                 varname = names(unlist(unname(parelmts))),
                 coef = paste0(info$parname,
                               if (length(unlist(pars)) > 1L)
                                 paste0("[", unlist(parelmts), "]")
                 ),
                 stringsAsFactors = FALSE
      )
    }

    nonprop <- unlist(unname(lapply(info$parelmts, attr, "nonprop")),
                      recursive = FALSE)

    if (!is.null(unlist(nonprop))) {
      out <- rbind(out,
                   data.frame(outcome = unname(info$varname),
                              outcat = rep(names(nonprop),
                                           ivapply(nonprop, length)),
                              varname = unlist(lapply(nonprop, names)),
                              coef = paste0(info$parname,
                                            paste0("[", unlist(nonprop), "]")),
                              stringsAsFactors = FALSE)
      )
    }

    if (!is.null(out)) {
      out$varnam_print <- cvapply(seq_along(out$outcat), function(k) {
        switch(as.character(is.na(out$outcat[k])),
               "TRUE" = out$varname[k],
               "FALSE" = paste0(out$outcat[k], ": ", out$varname[k])
        )
      })
    }

    rownames(out) <- NULL
    out
  })
}



# data_list helpers --------------------------------------------------------

get_row <- function(dat, i) {
  row <- lapply(1:ncol(dat), function(j) {.subset2(dat, j)[i]})
  names(row) <- names(dat)
  attr(row, "class") <- "data.frame"
  attr(row, "row.names") <- 1L
  row
}


# seed value
set_seed <- function(seed) {
  if ((R.version$major > 3L |
       (R.version$major == 3L & R.version$minor >= 6.0)) &
      Sys.getenv("IS_CHECK") == "true") {
    suppressWarnings(set.seed(seed, sample.kind = "Rounding"))
  } else {
    set.seed(seed)
  }
}

## Import internal functions from JointAI package
##
utils::globalVariables(c("i","U","mvar","seed","mess","warn","linkinv","colrev","ord_cov_dummy ",
                         "pattern", "Imputation_","firstm", "M","trtvar","x","autoburnin"))

#get_Mlist <- getFromNamespace("get_Mlist","JointAI")
get_terms_list <- getFromNamespace("get_terms_list","JointAI")
prep_arglist  <- getFromNamespace("prep_arglist","JointAI")
get_model_info <- getFromNamespace("get_model_info","JointAI")
get_1model_dim <- getFromNamespace("get_1model_dim","JointAI")
get_model_dim <- getFromNamespace("get_model_dim","JointAI")

identify_level_relations <- getFromNamespace("identify_level_relations","JointAI")
check_formula_list <- getFromNamespace("check_formula_list","JointAI")
split_formula_list <- getFromNamespace("split_formula_list","JointAI")
get_data_list <- getFromNamespace("get_data_list","JointAI")
#check_data <- getFromNamespace("check_data","JointAI")
make_filename <- getFromNamespace("make_filename","JointAI")
#get_subset <- getFromNamespace("get_subset","JointAI")

model_matrix_combi <- getFromNamespace("model_matrix_combi","JointAI")
get_initial_values <- getFromNamespace("get_initial_values","JointAI")
get_params <- getFromNamespace("get_params","JointAI")
get_future_info <- getFromNamespace("get_future_info","JointAI")
write_model <- getFromNamespace("write_model","JointAI")

extract_outcome_data <- getFromNamespace("extract_outcome_data","JointAI")
extract_id <- getFromNamespace("extract_id","JointAI")
get_groups <- getFromNamespace("get_groups","JointAI")
remove_grouping <- getFromNamespace("remove_grouping","JointAI")
check_varlevel <- getFromNamespace("check_varlevel","JointAI")
melt_data.frame_list <- getFromNamespace("melt_data.frame_list","JointAI")

reformat_longsurvdata <- getFromNamespace("reformat_longsurvdata","JointAI")
outcomes_to_mat <- getFromNamespace("outcomes_to_mat","JointAI")
prep_covoutcomes <- getFromNamespace("prep_covoutcomes","JointAI")
get_refs <- getFromNamespace("get_refs","JointAI")

extract_fcts <- getFromNamespace("extract_fcts","JointAI")
match_interaction <- getFromNamespace("match_interaction","JointAI")
get_linpreds <- getFromNamespace("get_linpreds","JointAI")
get_nonprop_lp <- getFromNamespace("get_nonprop_lp","JointAI")
check_cluster <- getFromNamespace("check_cluster","JointAI")
all_vars <- getFromNamespace("all_vars","JointAI")
#get_scale_pars <- getFromNamespace("get_scale_pars","JointAI")
find_scalevars <- getFromNamespace("find_scalevars","JointAI")
expand_rd_vcov_full <- getFromNamespace("expand_rd_vcov_full","JointAI")

run_parallel <- getFromNamespace("run_parallel","JointAI")
run_seq <- getFromNamespace("run_seq","JointAI")
rescale <- getFromNamespace("rescale","JointAI")
fitted_values <- getFromNamespace("fitted_values","JointAI")
fill_locf <- getFromNamespace("fill_locf","JointAI")

plot_prep <- getFromNamespace("plot_prep","JointAI")
get_plotmain <- getFromNamespace("get_plotmain","JointAI")

check_vars_in_data <- getFromNamespace("check_vars_in_data","JointAI")
check_classes <- getFromNamespace("check_classes","JointAI")
drop_levels <- getFromNamespace("drop_levels","JointAI")
computeP <- getFromNamespace("computeP","JointAI")
print_type <- getFromNamespace("print_type","JointAI")

get_rng <- getFromNamespace("get_rng","JointAI")
get_future_info <- getFromNamespace("get_future_info","JointAI")

remove_lhs <- getFromNamespace("remove_lhs","JointAI")
prep_covoutcomes <- getFromNamespace("prep_covoutcomes","JointAI")
default_hyperpars <- getFromNamespace("default_hyperpars","JointAI")

## @keywords internal
bs <- splines::bs

## @keywords internal
ns <- splines::ns

## @keywords internal
Surv <- survival::Surv
