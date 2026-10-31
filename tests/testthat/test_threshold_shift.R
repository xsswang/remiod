library("remiod")

# threshold_shift() makes the thresholds gamma_<y> in the (re-scaled) MCMC
# sample fit the JAGS node eta_<y>, which is on the scale of the scaled
# covariates: gamma + shift must equal the threshold on the scale used in JAGS.

test_that("threshold_shift() is 0 when nothing is scaled", {
  MCMC <- matrix(1:6, nrow = 3, dimnames = list(NULL, c("beta[1]", "beta[2]")))
  coefs <- data.frame(varname = c("x", "z"), coef = c("beta[1]", "beta[2]"))

  expect_equal(remiod:::threshold_shift(MCMC, coefs, scale_pars = NULL),
               rep(0, 3))
})


test_that("threshold_shift() uses only centred, proportional-odds coefficients", {
  MCMC <- cbind("beta[1]" = c(1, 2), "beta[2]" = c(10, 20),
                "beta[3]" = c(100, 200), "beta[4]" = c(5, 5))
  coefs <- data.frame(varname = c("x", "z", "w", "x"),
                      coef = paste0("beta[", 1:4, "]"),
                      outcat = c(NA, NA, NA, "2"))
  # z not centred (NA), w not in scale_pars, beta[4] non-proportional
  scale_pars <- data.frame(center = c(3, NA), scale = c(2, NA),
                           row.names = c("x", "z"))

  shift <- remiod:::threshold_shift(MCMC, coefs, scale_pars)

  if (utils::packageVersion("JointAI") <= "1.1.0") {
    # older JointAI versions do not transform the thresholds back
    expect_equal(shift, c(0, 0))
  } else {
    expect_equal(shift, c(1, 2) * 3)
  }
})


test_that("re-scaled thresholds + threshold_shift() = thresholds in JAGS", {
  # model set-up only (no sampling); age is the only scaled covariate
  data(schizow)
  dat <- schizow
  set.seed(2024)
  dat$age <- round(rnorm(nrow(dat), 40, 10))

  # suppressWarnings(): remiod uses functions that are deprecated in JointAI
  fit <- suppressWarnings(
    remiod(y6 ~ tx + age + y0 + y1 + y3, data = dat, trtvar = "tx",
           algorithm = "jags", n.adapt = 0, n.iter = 0, n.chains = 1,
           scale_vars = "age", warn = FALSE, mess = FALSE)
  )
  obj <- fit$mc.mar

  scale_pars <- do.call(rbind, unname(obj$Mlist$scale_pars))
  coefs <- do.call(rbind, obj$coef_list)
  ordinal <- names(obj$info_list)[
    vapply(obj$info_list, function(x) x$modeltype == "clm", logical(1))]
  gammas <- unlist(lapply(ordinal, function(v) {
    paste0("gamma_", v, "[", seq_len(obj$info_list[[v]]$ncat - 1), "]")
  }))

  # an artificial sample on the scale used in JAGS
  pars <- c(unique(coefs$coef), gammas)
  MCMC_jags <- matrix(rnorm(20 * length(pars)), nrow = 20,
                      dimnames = list(NULL, pars))

  # transform back as remiod does after sampling (JointAI's rescale())
  MCMC <- remiod:::rescale(MCMC_jags, coefs = coefs, scale_pars = scale_pars,
                           info_list = obj$info_list,
                           data_list = obj$data_list,
                           groups = obj$Mlist$groups)

  # the test is only meaningful if age is centred in the model of y6
  expect_true(scale_pars["age", "center"] != 0)

  for (v in ordinal) {
    gamma_v <- grep(paste0("^gamma_", v, "\\["), colnames(MCMC), value = TRUE)
    shift <- remiod:::threshold_shift(MCMC, obj$coef_list[[v]], scale_pars)

    expect_equal(MCMC[, gamma_v] + shift, MCMC_jags[, gamma_v],
                 label = paste("thresholds of", v))
  }
})
