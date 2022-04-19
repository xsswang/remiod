
paste_phi = function (nr, env = parent.frame()) {
  vn <- env$.info$varname
  ind <- env$.index
  gk <- env$.isgk
  paste0("phi( c_",vn,"[", nr, "] - eta_",vn,"[",ind, "] ) ")
}

write_opm = function (info, index, isgk = FALSE, indent = 4L) {
  .info <- info
  .index <- index
  .isgk <- isgk
  probs <- cvapply(2L:(info$ncat - 1L), function(k, .info = info, .index = index, .isgk = isgk) {
    paste0(tab(indent), paste_p(k), " <- ", (if (isTRUE(.info$rev)) {
      paste0(paste_phi(k - 1L), " - ", paste_phi(k))
    }
    else {
      paste0(paste_phi(k), " - ", paste_phi(k - 1L))
    }))
  })

  paste0(tab(indent), paste_p(1L), " <- ", if (isTRUE(info$rev)) {
    1 - paste_phi(1L)
  } else {
    paste_phi(1L)
  }, "\n",
  paste(probs, collapse = "\n"), "\n",

  tab(indent), paste_p(info$ncat), " <- ", if (isTRUE(info$rev)) {
    paste_phi(info$ncat)
  } else {
    paste0("1 - ", paste_phi(info$ncat - 1L))
  })
}

write_priors_opm = function (info) {

  # sign <- if (isTRUE(info$rev)) {
  #   " + "
  # } else {
  #   " - "
  # }
  gammas <- cvapply(1L:(info$ncat - 1L), function(k) {
    if (k == 1L) {
      paste0(tab(), "gamma_", info$varname, "[",
             k, "] ~ dnorm(mu_delta_ordinal, tau_delta_ordinal)")
    }
    else {
      paste0(tab(), "gamma_", info$varname, "[",
             k, "] ~ dgamma (0.001, 0.001)")
    }
  })

  deltas <- cvapply(1L:(info$ncat - 1L), function(k) {
    paste0(tab(), "c_", info$varname, "[",
           k, "] <- sum( gamma_", info$varname, "[1:",k,"])")
  })

  paste0(c(gammas, "", deltas), collapse = "\n")
}




jagsmodel_opm = function (info) {
  if (info$ncat < 3L) {
    errormsg("A cumulative logit mixed model is supposed to be fitted for the\n
             variable %s but %s only has %s categories.",
             dQuote(info$varname), dQuote(info$varname), info$ncat)
  }
  if (!is.null(info$hc_list)) {
    errormsg("I found a random effects structure. Did you mean to use %s\n
             instead of %s?",
             dQuote("clmm"), dQuote("clm"))
  }
  indent <- 4L + 4L + nchar(info$varname) + 7L
  index <- info$index[gsub("M_", "", info$resp_mat)]
  linpred <- if (length(info$lp[[info$resp_mat]]) > 0L) {
    paste_linpred(info$parname, info$parelmts[[info$resp_mat]],
                  matnam = info$resp_mat, index = index, cols = info$lp[[info$resp_mat]],
                  scale_pars = info$scale_pars[[info$resp_mat]])
  }
  else {
    "0"
  }
  linpred_nonprop <- if (!is.null(attr(info$parelmts[[info$resp_mat]],
                                       "nonprop"))) {
    rhs <- cvapply(attr(info$parelmts[[info$resp_mat]], "nonprop"),
                   function(par_elmts) {
                     add_linebreaks(paste_linpred(info$parname, par_elmts,
                                                  matnam = info$resp_mat, index = index, cols = attr(info$lp,
                                                                                                     "nonprop")[[info$resp_mat]], scale_pars = info$scale_pars[[info$resp_mat]]),
                                    indent = indent)
                   })
    paste0("\n\n", paste0(tab(4L), "eta_", info$varname,
                          "_", seq_along(rhs), "[", index, "] <- ",
                          rhs, collapse = "\n"))
  }
  dummies <- if (!is.null(info$dummy_cols)) {
    paste0("\n", paste0(paste_dummies(resp_mat = info$resp_mat,
                                      resp_col = info$resp_col, dummy_cols = info$dummy_cols,
                                      index = index, refs = info$refs), collapse = "\n"),
           "\n")
  }
  paste0("\r", tab(), add_dashes(paste0("# Cumulative logit model for ",
                                        info$varname)), "\n", tab(), "for (", index,
         " in 1:", info$N[gsub("M_", "", info$resp_mat)],
         ") {", "\n", tab(4L), info$resp_mat, "[",
         index, ", ", info$resp_col, "] ~ dcat(p_",
         info$varname, "[", index, ", 1:", info$ncat,
         "])", "\n", tab(4L), "eta_", info$varname,
         "[", index, "] <- ", add_linebreaks(linpred,
                                             indent = indent), linpred_nonprop, "\n\n",
         write_opm(info, index), "\n\n",
         #write_logits(info,index, nonprop = !is.null(linpred_nonprop)), "\n",
         dummies,
         info$trafos, tab(), "}", "\n\n",
         tab(), "# Priors for the model for ", info$varname,
         "\n", if (!is.null(info$lp[[info$resp_mat]])) {
           paste0(tab(), "for (k in ", min(unlist(c(info$parelmts,
                                                    lapply(info$parelmts, attr, "nonprop")))),
                  ":", max(unlist(c(info$parelmts, lapply(info$parelmts,
                                                          attr, "nonprop")))), ") {", "\n",
                  get_priordistr(info$shrinkage, type = "ordinal",
                                 parname = info$parname), tab(), "}",
                  "\n\n")
         }, write_priors_opm(info))
}






